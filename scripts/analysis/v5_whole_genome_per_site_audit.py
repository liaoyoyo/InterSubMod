#!/usr/bin/env python3
"""V5 Audit Suite — Whole-Genome per-site HP family + paired concordance.

Multi-threaded version: uses multiprocessing.Pool to split sites by chromosome
across N workers. Each worker opens its own BAM handle.

Inputs:
  5 BAMs: PA_93 (paired truth) + BL_93 + V5_93 + BL_06 + V5_06
  2 VCFs: ClairS-TO @ 0.93 + ClairS-TO @ 0.6 (whole genome PASS sites)

Outputs:
  TSV: data/wg_per_site_hp.tsv.gz       — long format ~450K rows
  TSV: data/wg_summary.tsv              — aggregate per-sample stats
  TSV: data/wg_paired_concordance.tsv   — site-level BL_93/V5_93/BL_06/V5_06 vs PA_93
  PNG: figures/20_whole_genome/wg_summary.png
  PNG: figures/20_whole_genome/wg_paired_concordance.png

Performance:
  - Per BAM: ~50K sites split into ~24 chromosome chunks, processed in parallel
  - Total: 5 BAMs sequential × ~5-10 min/BAM = ~25-50 min
"""
from __future__ import annotations
from pathlib import Path
import sys
import gzip
import time
from collections import Counter
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

plt.rcParams.update({
    "font.family": ["DejaVu Sans", "Droid Sans Fallback"],
    "font.sans-serif": ["DejaVu Sans", "Droid Sans Fallback"],
    "font.size": 10,
    "axes.unicode_minus": False,
})

N_WORKERS = 24

BAMS = {
    "PA_93": "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam",
    "BL_93": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
    "V5_93": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam",
    "BL_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/baseline_06/tumor_tagged.bam",
    "V5_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v5_06/tumor_tagged.bam",
}
VCFS = {
    "0.93": "/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz",
    "0.6":  "/big8_disk/data/HCC1395/ONT/subsample/t30_n20/ClairS_TO_v0_3_0/snv.vcf.gz",
}

OUT_BASE = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite")
DATA_DIR = OUT_BASE / "data"
FIG_DIR = OUT_BASE / "figures" / "20_whole_genome"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

CHROMS = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY")


def load_pass_sites(vcf: str) -> list:
    """Load PASS SNV sites from VCF."""
    sites = []
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in CHROMS:
                continue
            if f[6] != "PASS":
                continue
            ref, alt = f[3], f[4]
            if len(ref) > 1 or len(alt) > 1:
                continue
            sites.append((f[0], int(f[1]), ref, alt))
    return sites


def hp_family(hp):
    if hp is None:    return "untagged"
    if hp in (1, 11): return "HP1"
    if hp in (2, 21): return "HP2"
    if hp == 33:      return "HP33"
    return f"HP{hp}"


def process_chrom_chunk(args):
    """Worker: process all sites within one chromosome for one BAM."""
    bam_path, sites_in_chrom, sample_name, chrom = args
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom_name, pos, ref, alt in sites_in_chrom:
            counter = Counter()
            try:
                for pile in bam.pileup(chrom_name, pos - 1, pos,
                                        min_base_quality=0, min_mapping_quality=1,
                                        truncate=True):
                    if pile.reference_pos != pos - 1:
                        continue
                    for pread in pile.pileups:
                        if pread.is_del or pread.is_refskip:
                            continue
                        read = pread.alignment
                        if read.is_secondary or read.is_supplementary:
                            continue
                        qpos = pread.query_position
                        if qpos is None:
                            continue
                        if read.query_sequence[qpos].upper() != alt.upper():
                            continue
                        try:
                            hp = read.get_tag("HP")
                        except KeyError:
                            hp = None
                        counter[hp_family(hp)] += 1
            except ValueError:
                pass
            rows.append({
                "sample": sample_name, "chrom": chrom_name, "pos": pos,
                "HP1": counter.get("HP1", 0),
                "HP2": counter.get("HP2", 0),
                "HP33": counter.get("HP33", 0),
                "untagged": counter.get("untagged", 0),
                "total_alt": sum(counter.values()),
            })
    return rows


def per_site_hp_alt_parallel(bam_path: str, sites: list, sample_name: str,
                              n_workers: int = N_WORKERS) -> pd.DataFrame:
    """Multi-process per-site HP analysis: split by chromosome."""
    if not Path(bam_path).exists():
        print(f"  [skip] {sample_name}: missing", file=sys.stderr)
        return pd.DataFrame()

    # Group sites by chromosome
    chrom_sites = {}
    for s in sites:
        chrom_sites.setdefault(s[0], []).append(s)

    # Build args for workers
    args_list = [(bam_path, sites_in_chrom, sample_name, chrom)
                 for chrom, sites_in_chrom in chrom_sites.items()]

    print(f"  [{sample_name}] {len(sites):,} sites split into "
          f"{len(args_list)} chrom chunks, {n_workers} workers ...",
          file=sys.stderr, flush=True)

    t0 = time.time()
    all_rows = []
    with Pool(n_workers) as pool:
        for i, chunk_rows in enumerate(pool.imap_unordered(process_chrom_chunk, args_list)):
            all_rows.extend(chunk_rows)
            if (i + 1) % 6 == 0 or (i + 1) == len(args_list):
                elapsed = time.time() - t0
                print(f"    [{sample_name}] {i+1}/{len(args_list)} chunks done "
                      f"({elapsed:.0f}s elapsed)", file=sys.stderr, flush=True)
    elapsed = time.time() - t0
    print(f"  [{sample_name}] DONE in {elapsed:.0f}s "
          f"({len(all_rows):,} sites)", file=sys.stderr, flush=True)
    return pd.DataFrame(all_rows)


def main():
    t_global = time.time()
    print(f"=== WG per-site audit (parallel, {N_WORKERS} workers) ===",
          file=sys.stderr, flush=True)

    # Load PASS sites once
    sites_93 = load_pass_sites(VCFS["0.93"])
    sites_06 = load_pass_sites(VCFS["0.6"])
    print(f"  PASS sites @ 0.93 WG: {len(sites_93):,}", file=sys.stderr, flush=True)
    print(f"  PASS sites @ 0.6  WG: {len(sites_06):,}", file=sys.stderr, flush=True)

    # Per-sample collection (sequential to avoid IO conflict)
    all_dfs = []
    for sample, bam in BAMS.items():
        sites = sites_93 if sample.endswith("_93") else sites_06
        df = per_site_hp_alt_parallel(bam, sites, sample)
        if not df.empty:
            all_dfs.append(df)
            hp1 = df["HP1"].sum()
            hp2 = df["HP2"].sum()
            hp33 = df["HP33"].sum()
            ratio = hp1 / max(hp2, 1)
            ratio_str = f"{ratio:.2f}:1" if ratio >= 1 else f"1:{1/ratio:.2f}"
            print(f"  >> [{sample}] HP1={hp1:,} HP2={hp2:,} HP33={hp33:,} "
                  f"ratio={ratio_str}", file=sys.stderr, flush=True)

    if not all_dfs:
        print("ERROR: no data", file=sys.stderr); return

    long_df = pd.concat(all_dfs, ignore_index=True)
    out_tsv = DATA_DIR / "wg_per_site_hp.tsv.gz"
    long_df.to_csv(out_tsv, sep="\t", index=False, compression="gzip")
    print(f"\n  saved per-site TSV: {out_tsv} ({len(long_df):,} rows)",
          file=sys.stderr, flush=True)

    # Aggregate summary
    summary = (long_df.groupby("sample")
               .agg(n_sites=("total_alt", "size"),
                    total_alt_reads=("total_alt", "sum"),
                    HP1=("HP1", "sum"),
                    HP2=("HP2", "sum"),
                    HP33=("HP33", "sum"),
                    untagged=("untagged", "sum"))
               .reset_index())
    summary["HP1_pct"] = summary["HP1"] / summary["total_alt_reads"] * 100
    summary["HP2_pct"] = summary["HP2"] / summary["total_alt_reads"] * 100
    summary["HP33_pct"] = summary["HP33"] / summary["total_alt_reads"] * 100
    summary["untagged_pct"] = summary["untagged"] / summary["total_alt_reads"] * 100
    summary["HP1_HP2_ratio"] = summary["HP1"] / summary["HP2"].clip(lower=1)
    summary["sites_with_alt"] = summary["total_alt_reads"]  # alias
    summary.to_csv(DATA_DIR / "wg_summary.tsv", sep="\t", index=False)

    print("\n=== Whole-Genome Aggregate Summary ===", file=sys.stderr)
    print(summary[["sample", "n_sites", "total_alt_reads",
                   "HP1", "HP2", "HP33", "HP1_HP2_ratio"]].to_string(index=False),
          file=sys.stderr, flush=True)

    # ============================================================
    # Site-level paired concordance: BL_93/V5_93/BL_06/V5_06 vs PA_93
    # ============================================================
    if "PA_93" in long_df["sample"].values:
        # Pivot to wide: row = (chrom, pos), columns = sample's HP1/HP2 dominant
        def hp_dominant(row):
            """Return 'HP1', 'HP2', or 'tie' or 'no_data'."""
            if row["total_alt"] == 0:
                return "no_data"
            if row["HP1"] > row["HP2"] and row["HP1"] > row["HP33"]:
                return "HP1"
            if row["HP2"] > row["HP1"] and row["HP2"] > row["HP33"]:
                return "HP2"
            if row["HP33"] > row["HP1"] and row["HP33"] > row["HP2"]:
                return "HP33"
            return "tie"

        long_df["dominant"] = long_df.apply(hp_dominant, axis=1)
        wide = long_df.pivot_table(index=["chrom", "pos"], columns="sample",
                                     values="dominant", aggfunc="first")
        wide = wide.reset_index()

        rows = []
        for cmp in ["BL_93", "V5_93", "BL_06", "V5_06"]:
            if cmp not in wide.columns:
                continue
            valid = wide[(wide["PA_93"].isin(["HP1", "HP2"])) &
                         (wide[cmp].isin(["HP1", "HP2", "HP33"]))]
            n = len(valid)
            match = (valid["PA_93"] == valid[cmp]).sum()
            hp33 = (valid[cmp] == "HP33").sum()
            opposite = ((valid["PA_93"] == "HP1") & (valid[cmp] == "HP2") |
                        (valid["PA_93"] == "HP2") & (valid[cmp] == "HP1")).sum()
            rows.append({
                "comparison": f"{cmp} vs PA_93",
                "n_sites_evaluable": n,
                "match": match,
                "match_pct": match/max(n,1)*100,
                "opposite": opposite,
                "opposite_pct": opposite/max(n,1)*100,
                "hp33_conservative": hp33,
                "hp33_pct": hp33/max(n,1)*100,
            })
        paired_df = pd.DataFrame(rows)
        paired_df.to_csv(DATA_DIR / "wg_paired_concordance.tsv",
                         sep="\t", index=False)
        print("\n=== Whole-Genome Paired Concordance vs PA_93 ===",
              file=sys.stderr, flush=True)
        print(paired_df.to_string(index=False), file=sys.stderr, flush=True)

    # ============================================================
    # Figures
    # ============================================================
    samples_order = [s for s in ["PA_93", "BL_93", "V5_93", "BL_06", "V5_06"]
                      if s in summary["sample"].values]
    sub = summary.set_index("sample").reindex(samples_order)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Panel A: HP1:HP2 ratio
    ax = axes[0]
    colors = {"PA_93": "#7B1FA2", "BL_93": "#C62828", "V5_93": "#1565C0",
              "BL_06": "#FF9800", "V5_06": "#2E7D32"}
    bars = ax.bar(sub.index, sub["HP1_HP2_ratio"],
                  color=[colors.get(s, "#888") for s in sub.index],
                  alpha=0.85, edgecolor="black")
    for bar, (s, row) in zip(bars, sub.iterrows()):
        r = row["HP1_HP2_ratio"]
        ratio_str = f"{r:.2f}:1" if r >= 1 else f"1:{1/r:.2f}"
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.05,
                f"{ratio_str}\nn={int(row['total_alt_reads']):,}",
                ha="center", fontsize=9, fontweight="bold")
    ax.axhline(1.0, ls="--", color="#2E7D32", alpha=0.5, label="Expected = 1")
    ax.set_yscale("log")
    ax.set_ylabel("HP1 : HP2 ratio (log scale)")
    ax.set_title(f"Whole-Genome somatic-ALT-only HP1:HP2 ratio",
                 fontsize=12, fontweight="bold")
    ax.legend()

    # Panel B: HP family % stacked
    ax = axes[1]
    bottom = np.zeros(len(sub))
    palette = {"HP1": "#E91E63", "HP2": "#1976D2",
               "HP33": "#7B1FA2", "untagged": "#9E9E9E"}
    for fam in ["HP1", "HP2", "HP33", "untagged"]:
        vals = sub[f"{fam}_pct"].values
        ax.bar(sub.index, vals, bottom=bottom, label=fam,
               color=palette[fam], alpha=0.85)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 3:
                ax.text(i, b + v/2, f"{v:.0f}%",
                        ha="center", va="center", fontsize=9, fontweight="bold",
                        color="white" if fam in ("HP1", "HP2", "HP33") else "black")
        bottom += vals
    ax.set_ylabel("% of ALT reads")
    ax.set_ylim(0, 110)
    ax.set_title("Whole-Genome HP family composition",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", ncol=4)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "wg_summary.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  saved fig: {FIG_DIR}/wg_summary.png", file=sys.stderr, flush=True)

    print(f"\n=== TOTAL ELAPSED: {time.time()-t_global:.0f}s ===",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
