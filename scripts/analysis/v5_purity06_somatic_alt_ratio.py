#!/usr/bin/env python3
"""V5 Audit Suite — somatic-ALT-only HP ratio (對應 17.3:1 metric).

This is the metric that defines self-phasing strength:
  Among reads supporting somatic ALT alleles at PASS sites,
  count HP1 family vs HP2 family.

Reuses ClairS-TO PASS sites + checks each read's base at the somatic position.
"""
from __future__ import annotations
from pathlib import Path
import sys
import gzip
from collections import Counter

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

BAMS = {
    "BL_93": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
    "V5_93": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam",
    "BL_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/baseline_06/tumor_tagged.bam",
    "V5_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v5_06/tumor_tagged.bam",
}

# Each scenario has its own ClairS-TO VCF (PASS sites differ)
VCFS = {
    "0.93": "/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz",
    "0.6": "/big8_disk/data/HCC1395/ONT/subsample/t30_n20/ClairS_TO_v0_3_0/snv.vcf.gz",
}

OUT_BASE = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite")
OUT_DATA = OUT_BASE / "data"
OUT_FIG = OUT_BASE / "figures" / "09_purity06"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

CHROMS = ("chr19", "chr20", "chr21", "chr22")


def load_pass_sites(vcf: str, chroms: tuple) -> list:
    """Load ClairS-TO PASS sites in given chrs. Returns [(chrom, pos, ref, alt)]."""
    sites = []
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in chroms:
                continue
            if f[6] != "PASS":
                continue
            sites.append((f[0], int(f[1]), f[3], f[4]))
    return sites


def hp_family(hp):
    if hp is None:
        return "untagged"
    if hp in (1, 11):
        return "HP1"
    if hp in (2, 21):
        return "HP2"
    if hp == 33:
        return "HP33"
    return f"HP{hp}"


def count_alt_reads_hp(bam_path: str, sites: list) -> Counter:
    """For each PASS site, find reads where read base == ALT, then count HP family."""
    counter = Counter()
    n_sites_done = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, pos, ref, alt in sites:
            if len(alt) > 1 or len(ref) > 1:
                continue  # skip indel
            try:
                for pile in bam.pileup(chrom, pos - 1, pos,
                                        min_base_quality=0,
                                        min_mapping_quality=1,
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
                        base = read.query_sequence[qpos].upper()
                        if base != alt.upper():
                            continue
                        hp = None
                        try:
                            hp = read.get_tag("HP")
                        except KeyError:
                            pass
                        counter[hp_family(hp)] += 1
            except ValueError:
                pass
            n_sites_done += 1
            if n_sites_done % 1000 == 0:
                print(f"    sites processed: {n_sites_done}/{len(sites)}", file=sys.stderr)
    return counter


def main():
    print("=== somatic-ALT-only HP ratio analysis ===", file=sys.stderr)

    # Load PASS sites for each scenario
    sites_93 = load_pass_sites(VCFS["0.93"], CHROMS)
    sites_06 = load_pass_sites(VCFS["0.6"], CHROMS)
    print(f"  PASS sites @ 0.93 chr19-22: {len(sites_93):,}", file=sys.stderr)
    print(f"  PASS sites @ 0.6  chr19-22: {len(sites_06):,}", file=sys.stderr)

    rows = []
    for label, bam_path in BAMS.items():
        if not Path(bam_path).exists():
            print(f"  [{label}] missing, skip", file=sys.stderr)
            continue
        purity = "0.93" if label.endswith("_93") else "0.6"
        sites = sites_93 if purity == "0.93" else sites_06
        print(f"\n  [{label}] processing {len(sites):,} PASS sites ...", file=sys.stderr)
        c = count_alt_reads_hp(bam_path, sites)
        total = sum(c.values())
        hp1 = c.get("HP1", 0)
        hp2 = c.get("HP2", 0)
        ratio = hp1 / max(hp2, 1)
        ratio_str = f"{ratio:.2f}:1" if ratio >= 1 else f"1:{1/ratio:.2f}"
        rows.append({
            "label": label,
            "purity": purity,
            "n_sites": len(sites),
            "total_alt_reads": total,
            "HP1_family": hp1,
            "HP2_family": hp2,
            "HP33": c.get("HP33", 0),
            "untagged_or_HP0": c.get("untagged", 0) + c.get("HP0", 0),
            "ratio_hp1_hp2": ratio,
            "ratio_str": ratio_str,
        })
        print(f"    [{label}] total_alt={total:,} HP1={hp1:,} HP2={hp2:,} HP33={c.get('HP33', 0):,} "
              f"-> ratio = {ratio_str}", file=sys.stderr)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DATA / "purity06_somatic_alt_hp_ratio.tsv", sep="\t", index=False)
    print(f"\n  saved TSV: {OUT_DATA / 'purity06_somatic_alt_hp_ratio.tsv'}", file=sys.stderr)
    print("\n=== Summary ===", file=sys.stderr)
    print(df[["label", "purity", "n_sites", "total_alt_reads",
              "HP1_family", "HP2_family", "ratio_str"]].to_string(index=False), file=sys.stderr)

    # ===== Figure: somatic-ALT-only HP ratio =====
    if not df.empty:
        fig, ax = plt.subplots(figsize=(11, 6.5))
        order = [l for l in ["BL_93", "V5_93", "BL_06", "V5_06"] if l in df["label"].values]
        df_p = df.set_index("label").reindex(order)
        colors = {"BL_93": "#C62828", "V5_93": "#1565C0",
                  "BL_06": "#FF9800", "V5_06": "#2E7D32"}
        ratios = df_p["ratio_hp1_hp2"].values
        # 對 ratio < 1 做 negative log scale 風格
        display_vals = [r if r >= 1 else 1 / r for r in ratios]
        signs = ["HP1 dominant" if r >= 1 else "HP2 dominant" for r in ratios]
        bars = ax.bar(df_p.index, display_vals,
                      color=[colors.get(l, "#888") for l in df_p.index],
                      alpha=0.85, edgecolor="black", linewidth=1.2)
        for bar, (label, row) in zip(bars, df_p.iterrows()):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() * 1.05,
                    f"{row['ratio_str']}\n"
                    f"HP1={row['HP1_family']:,}\n"
                    f"HP2={row['HP2_family']:,}\n"
                    f"({row['n_sites']:,} sites,\n {row['total_alt_reads']:,} ALT reads)",
                    ha="center", fontsize=9, fontweight="bold")
        ax.axhline(1.0, ls="--", color="#2E7D32", alpha=0.5, label="Expected ratio = 1")
        ax.set_ylabel("|ratio| (max(HP1/HP2, HP2/HP1)) [log scale]")
        ax.set_yscale("log")
        ax.set_title("Figure C — Somatic-ALT-only HP1:HP2 ratio (chr19-chr22)\n"
                     "（17.3:1 metric 對應的真實 self-phasing 強度）",
                     fontsize=12, fontweight="bold")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(OUT_FIG / "figC_somatic_alt_ratio.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved figC: {OUT_FIG / 'figC_somatic_alt_ratio.png'}", file=sys.stderr)


if __name__ == "__main__":
    main()
