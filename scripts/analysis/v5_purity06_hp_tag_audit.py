#!/usr/bin/env python3
"""V5 Audit Suite — 0.6 purity HP tag analysis.

Compare HP1:HP2 self-phasing strength + paired concordance between:
  - baseline @ 0.93 purity (existing)
  - V5 @ 0.93 purity (existing)
  - baseline @ 0.6 purity (mixed BAM, NEW)
  - V5 @ 0.6 purity (mixed BAM, NEW)

Hypothesis to verify:
  H1: baseline @ 0.6 self-phasing strength weakens (17.3:1 → 3-5:1 expected)
  H2: V5 @ 0.6 ratio remains ≈ 1:1 (PON-only logic purity-independent)
  H3: V5 vs baseline gap decreases at 0.6 (smaller bias to fix)

Outputs
-------
- TSV: data/purity06_hp_summary.tsv     — aggregate HP family counts × 4 scenarios
- TSV: data/purity06_per_site.tsv       — 15 sites × 4 scenarios
- PNG: figures/09_purity06/figA_aggregate_ratio.png
- PNG: figures/09_purity06/figB_per_site_comparison.png
- PNG: figures/09_purity06/figC_self_phasing_decay.png
"""
from __future__ import annotations
from pathlib import Path
import sys
from collections import Counter, defaultdict

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

# Configuration
BAMS = {
    "BL_93": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam",
    "V5_93": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v5_somatic_fallback/tumor_tagged.bam",
    "BL_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/baseline_06/tumor_tagged.bam",
    "V5_06": "/big7_disk/liaoyoyo2001/longphase-to-mod/output/purity_06_simulation/v5_06/tumor_tagged.bam",
}

# 15 audit sites (reuse existing audit suite sites)
SITES = [
    ("TP01", "chr6", 145444893),
    ("TP02", "chr4", 70548355),
    ("TP03", "chr5", 153209947),
    ("TP04", "chr16", 35118902),
    ("TP05", "chr7", 109185781),
    ("FPA1", "chr8", 93565727),
    ("FPA2", "chr9", 137953060),
    ("FPB1", "chr7", 52087777),
    ("FPB2", "chr9", 75383880),
    ("V5max1", "chr19", 4639528),
    ("V5max2", "chr19", 2235521),
    ("V5max3", "chr19", 7405500),
    ("SP1", "chr19", 17565944),
    ("SP2", "chr19", 12452332),
    ("SP3", "chr19", 12467180),
]

OUT_BASE = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite")
OUT_FIG_DIR = OUT_BASE / "figures" / "09_purity06"
OUT_DATA_DIR = OUT_BASE / "data"
OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)


def hp_family(hp: int | None) -> str:
    """Map HP value to family: HP1 family = {1, 11}, HP2 family = {2, 21}, else."""
    if hp is None:
        return "untagged"
    if hp in (1, 11):
        return "HP1"
    if hp in (2, 21):
        return "HP2"
    if hp == 33:
        return "HP33"
    if hp == 0:
        return "HP0"
    return f"other:{hp}"


def fetch_hp_at_site(bam_path: str, chrom: str, pos: int, window: int = 0) -> Counter:
    """Fetch HP tags for primary reads covering position pos.

    Args:
      window: ±bp around pos. 0 = exact position only.
    """
    counts = Counter()
    if not Path(bam_path).exists():
        return counts
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        try:
            for read in bam.fetch(chrom, max(0, pos - window - 1), pos + window):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                if pos - 1 < read.reference_start or pos > read.reference_end:
                    continue
                hp = None
                try:
                    hp = read.get_tag("HP")
                except KeyError:
                    pass
                counts[hp_family(hp)] += 1
        except ValueError:
            pass  # chrom not in index
    return counts


def aggregate_hp_genome(bam_path: str, sample_chrs: tuple = ("chr19", "chr20", "chr21", "chr22")) -> Counter:
    """Aggregate HP family counts across selected chrs for speed."""
    counts = Counter()
    if not Path(bam_path).exists():
        return counts
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom in sample_chrs:
            try:
                for read in bam.fetch(chrom):
                    if read.is_secondary or read.is_supplementary or read.is_unmapped:
                        continue
                    hp = None
                    try:
                        hp = read.get_tag("HP")
                    except KeyError:
                        pass
                    counts[hp_family(hp)] += 1
            except ValueError:
                continue
    return counts


def main():
    print("=== Phase D: 0.6 purity HP tag audit ===", file=sys.stderr)

    # Verify all BAMs exist
    for label, path in BAMS.items():
        exists = Path(path).exists()
        print(f"  [{label}] {'✓' if exists else '✗ MISSING'}: {path}", file=sys.stderr)
        if not exists and label.endswith("_06"):
            print(f"    -> {label} not yet produced; skip in analysis", file=sys.stderr)

    # ===== Aggregate HP family counts (chr19-22 sample) =====
    print("\n=== Aggregate HP family (chr19-chr22) ===", file=sys.stderr)
    agg_rows = []
    for label, path in BAMS.items():
        if not Path(path).exists():
            continue
        counts = aggregate_hp_genome(path)
        total = sum(counts.values())
        print(f"  [{label}] total={total:,}", file=sys.stderr)
        for fam, n in counts.items():
            agg_rows.append({"label": label, "family": fam, "n": n,
                             "pct": n / max(total, 1) * 100})
    agg_df = pd.DataFrame(agg_rows)
    if not agg_df.empty:
        # Compute HP1:HP2 ratio per label
        ratios = []
        for label in agg_df["label"].unique():
            sub = agg_df[agg_df["label"] == label]
            hp1 = sub[sub["family"] == "HP1"]["n"].sum()
            hp2 = sub[sub["family"] == "HP2"]["n"].sum()
            ratio = hp1 / max(hp2, 1)
            ratios.append({"label": label, "hp1": hp1, "hp2": hp2,
                           "ratio_hp1_hp2": ratio,
                           "ratio_str": f"{ratio:.2f}:1" if ratio >= 1 else f"1:{1/ratio:.2f}"})
        ratio_df = pd.DataFrame(ratios)
        ratio_df.to_csv(OUT_DATA_DIR / "purity06_hp_ratio_summary.tsv",
                        sep="\t", index=False)
        print("\n=== HP1:HP2 ratio summary ===", file=sys.stderr)
        print(ratio_df.to_string(index=False), file=sys.stderr)
    agg_df.to_csv(OUT_DATA_DIR / "purity06_hp_aggregate.tsv", sep="\t", index=False)

    # ===== Per-site analysis (15 audit sites) =====
    print("\n=== Per-site HP analysis ===", file=sys.stderr)
    site_rows = []
    for label, path in BAMS.items():
        if not Path(path).exists():
            continue
        for site_name, chrom, pos in SITES:
            counts = fetch_hp_at_site(path, chrom, pos)
            total = sum(counts.values())
            site_rows.append({
                "label": label,
                "site": site_name,
                "chrom": chrom,
                "pos": pos,
                "n_reads": total,
                "HP1": counts.get("HP1", 0),
                "HP2": counts.get("HP2", 0),
                "HP33": counts.get("HP33", 0),
                "HP0": counts.get("HP0", 0),
                "untagged": counts.get("untagged", 0),
                "ratio_hp1_hp2": counts.get("HP1", 0) / max(counts.get("HP2", 1), 1),
            })
    site_df = pd.DataFrame(site_rows)
    site_df.to_csv(OUT_DATA_DIR / "purity06_per_site.tsv", sep="\t", index=False)
    print(f"  wrote per-site TSV: {len(site_df)} rows", file=sys.stderr)

    # ===== Figures =====
    if not agg_df.empty and len(ratio_df) >= 2:
        # ---- Figure A: HP1:HP2 ratio bar chart (4 scenarios) ----
        fig, ax = plt.subplots(figsize=(11, 6))
        order = [l for l in ["BL_93", "V5_93", "BL_06", "V5_06"] if l in ratio_df["label"].values]
        df_plot = ratio_df.set_index("label").reindex(order)
        colors = {"BL_93": "#C62828", "V5_93": "#1565C0",
                  "BL_06": "#FF9800", "V5_06": "#2E7D32"}
        bars = ax.bar(df_plot.index, df_plot["ratio_hp1_hp2"],
                      color=[colors.get(l, "#888") for l in df_plot.index],
                      alpha=0.85, edgecolor="black", linewidth=1.2)
        for bar, (label, row) in zip(bars, df_plot.iterrows()):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{row['ratio_str']}\n(HP1={row['hp1']:,}, HP2={row['hp2']:,})",
                    ha="center", fontsize=9, fontweight="bold")
        ax.axhline(1.0, ls="--", color="#2E7D32", alpha=0.5, label="Expected ratio = 1")
        ax.set_ylabel("HP1 : HP2 ratio (chr19-chr22)")
        ax.set_title("Figure A — Self-phasing strength (HP1:HP2) across 4 scenarios",
                     fontsize=13, fontweight="bold")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(OUT_FIG_DIR / "figA_hp_ratio_4scenarios.png",
                    dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved figA: {OUT_FIG_DIR / 'figA_hp_ratio_4scenarios.png'}", file=sys.stderr)

    if not site_df.empty:
        # ---- Figure B: Per-site HP1:HP2 ratio comparison ----
        fig, ax = plt.subplots(figsize=(15, 7))
        labels = sorted(site_df["label"].unique(),
                        key=lambda x: ["BL_93", "V5_93", "BL_06", "V5_06"].index(x)
                        if x in ["BL_93", "V5_93", "BL_06", "V5_06"] else 99)
        site_order = [s[0] for s in SITES]
        n_lab = len(labels)
        width = 0.8 / n_lab
        x = np.arange(len(site_order))
        for i, lbl in enumerate(labels):
            sub = site_df[site_df["label"] == lbl].set_index("site").reindex(site_order)
            color = {"BL_93": "#C62828", "V5_93": "#1565C0",
                     "BL_06": "#FF9800", "V5_06": "#2E7D32"}.get(lbl, "#888")
            ax.bar(x + (i - n_lab / 2 + 0.5) * width,
                   sub["ratio_hp1_hp2"].values,
                   width, label=lbl, color=color, alpha=0.85)
        ax.axhline(1.0, ls="--", color="#2E7D32", alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(site_order, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("HP1 : HP2 ratio")
        ax.set_yscale("log")
        ax.set_title("Figure B — Per-site HP1:HP2 ratio (15 audit sites × 4 scenarios)",
                     fontsize=13, fontweight="bold")
        ax.legend(loc="upper right", ncol=4)
        fig.tight_layout()
        fig.savefig(OUT_FIG_DIR / "figB_per_site_ratio.png",
                    dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved figB: {OUT_FIG_DIR / 'figB_per_site_ratio.png'}", file=sys.stderr)

    print("\n=== Phase D complete ===", file=sys.stderr)


if __name__ == "__main__":
    main()
