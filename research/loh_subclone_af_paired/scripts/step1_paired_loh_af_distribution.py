#!/usr/bin/env python3
"""
Step 1: Paired Mode LOH AF Distribution Baseline
==================================================
Parallel to TO step1 — establish AF distribution in LOH regions for paired mode.
LOH regions annotated using TO-derived LOH.bed files (LOH is a genomic property).

Input: all_region_rows.tsv.gz (master dataset, paired mode)
Output:
  - p01: AF distribution LOH vs Non-LOH (per-sample)
  - p02: Intermediate AF proportion across samples
  - p03: Fine-grain LOH AF distribution
  - p04: AF vs Coverage_Multiple scatter
  - p05: AF class by CN tier
  - step1_paired_af_class_statistics.tsv
  - step1_paired_loh_af_by_cn_tier.tsv
"""

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from utils import (
    MASTER, FIGDIR, DATADIR, SAMPLE_ORDER, SAMPLE_LABELS, SAMPLE_SHORT,
    AF_COLORS, classify_af, cn_tier, load_loh_beds, annotate_loh_bed,
    load_paired_data
)
import warnings
warnings.filterwarnings("ignore")

FIGDIR.mkdir(parents=True, exist_ok=True)
DATADIR.mkdir(parents=True, exist_ok=True)


def fig01_af_distribution(df):
    """p01: Per-sample AF histograms, LOH vs Non-LOH."""
    print("\n=== Figure p01: AF Distribution LOH vs Non-LOH ===")

    fig, axes = plt.subplots(2, 4, figsize=(24, 10))
    axes_flat = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        ds = df[df["sample"] == sample]

        for loh_status, color, label in [(True, "#E91E63", "LOH"), (False, "#90A4AE", "Non-LOH")]:
            sub = ds[ds["loh_bed_hit"] == loh_status]
            if len(sub) == 0:
                continue
            # TP only for main plot
            sub_tp = sub[sub["truth_label"] == "TP"]
            if len(sub_tp) > 0:
                ax.hist(sub_tp["caller_af"], bins=50, range=(0, 1), alpha=0.6,
                        color=color, label=f"{label} TP (n={len(sub_tp):,})", density=True)

        ax.set_xlabel("Caller AF")
        ax.set_ylabel("Density")
        ax.set_title(f"{SAMPLE_LABELS.get(sample, sample)}", fontsize=11)
        ax.legend(fontsize=7)

    axes_flat[-1].axis("off")
    fig.suptitle("Step 1: Paired Mode — AF Distribution in LOH vs Non-LOH Regions (TP only)\n"
                 "LOH defined by TO-derived LOH.bed",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p01_paired_af_distribution_loh_vs_nonloh.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig02_intermediate_proportion(df):
    """p02: Intermediate AF proportion across samples."""
    print("\n=== Figure p02: Intermediate AF Proportion ===")

    records = []
    for sample in SAMPLE_ORDER:
        ds = df[(df["sample"] == sample) & (df["truth_label"] == "TP")]
        for loh_status in [True, False]:
            sub = ds[ds["loh_bed_hit"] == loh_status]
            if len(sub) == 0:
                continue
            n_inter = (sub["af_class"] == "Intermediate").sum()
            records.append({
                "sample": SAMPLE_SHORT[sample],
                "loh": "LOH" if loh_status else "Non-LOH",
                "n_total": len(sub),
                "n_intermediate": n_inter,
                "pct_intermediate": 100 * n_inter / len(sub),
            })

    rec_df = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(12, 6))
    samples = [SAMPLE_SHORT[s] for s in SAMPLE_ORDER]
    x = np.arange(len(samples))
    width = 0.35

    loh_data = rec_df[rec_df["loh"] == "LOH"]
    nonloh_data = rec_df[rec_df["loh"] == "Non-LOH"]

    if len(loh_data) > 0:
        ax.bar(x[:len(loh_data)] - width/2, loh_data["pct_intermediate"].values, width,
               label="LOH TP", color="#E91E63", alpha=0.8)
    if len(nonloh_data) > 0:
        ax.bar(x[:len(nonloh_data)] + width/2, nonloh_data["pct_intermediate"].values, width,
               label="Non-LOH TP", color="#90A4AE", alpha=0.8)

    ax.set_xlabel("Sample")
    ax.set_ylabel("% Intermediate AF")
    ax.set_title("Paired Mode: Intermediate AF Proportion in LOH vs Non-LOH (TP)", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, fontsize=9)
    ax.legend()

    plt.tight_layout()
    out = FIGDIR / "p02_paired_intermediate_af_proportion.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return rec_df


def fig03_fine_grain(df):
    """p03: Fine-grain LOH AF distribution."""
    print("\n=== Figure p03: Fine-grain LOH AF ===")

    df_loh_tp = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP")]

    fig, axes = plt.subplots(2, 4, figsize=(24, 10))
    axes_flat = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        sub = df_loh_tp[df_loh_tp["sample"] == sample]
        if len(sub) > 0:
            ax.hist(sub["caller_af"], bins=100, range=(0, 1), color="#E91E63", alpha=0.7)
            n_inter = (sub["af_class"] == "Intermediate").sum()
            pct = 100 * n_inter / len(sub)
            ax.set_title(f"{SAMPLE_SHORT[sample]}\nn={len(sub):,}, {pct:.1f}% intermediate", fontsize=10)
        ax.set_xlabel("Caller AF")
        ax.set_ylabel("Count")

    axes_flat[-1].axis("off")
    fig.suptitle("Paired Mode: Fine-Grain AF Distribution in LOH Regions (TP only)", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p03_paired_loh_af_fine_grain.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig04_af_vs_cm(df):
    """p04: AF vs Coverage_Multiple in LOH."""
    print("\n=== Figure p04: AF vs Coverage_Multiple ===")

    df_loh_tp = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP")]

    fig, axes = plt.subplots(2, 4, figsize=(24, 10))
    axes_flat = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        sub = df_loh_tp[df_loh_tp["sample"] == sample]
        if len(sub) > 0:
            sub_plot = sub.sample(min(5000, len(sub)), random_state=42)
            colors = sub_plot["af_class"].map(AF_COLORS)
            ax.scatter(sub_plot["Coverage_Multiple"], sub_plot["caller_af"],
                       c=colors, alpha=0.3, s=5)
            ax.axhline(0.1, color="gray", ls="--", alpha=0.3)
            ax.axhline(0.9, color="gray", ls="--", alpha=0.3)
            ax.axvline(0.75, color="gray", ls="--", alpha=0.3)
            ax.axvline(1.25, color="gray", ls="--", alpha=0.3)
        ax.set_xlabel("Coverage_Multiple")
        ax.set_ylabel("Caller AF")
        ax.set_title(SAMPLE_SHORT[sample], fontsize=10)
        ax.set_xlim(0, 3)

    axes_flat[-1].axis("off")
    fig.suptitle("Paired Mode: AF vs Coverage_Multiple in LOH (TP)\n"
                 "Blue=Extreme, Green=Near-half, Orange=Intermediate", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p04_paired_loh_af_vs_coverage_multiple.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig05_af_by_cn_tier(df):
    """p05: AF class distribution by CN tier."""
    print("\n=== Figure p05: AF by CN Tier ===")

    df_loh_tp = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP")]
    cn_order = ["CN1", "CN2", "CN3", "CN4+"]
    af_order = ["Extreme", "Near-half", "Intermediate"]

    records = []
    for sample in SAMPLE_ORDER:
        for cn in cn_order:
            sub = df_loh_tp[(df_loh_tp["sample"] == sample) & (df_loh_tp["cn_tier"] == cn)]
            if len(sub) == 0:
                continue
            n_inter = (sub["af_class"] == "Intermediate").sum()
            records.append({
                "sample": SAMPLE_SHORT[sample],
                "cn_tier": cn,
                "truth": "TP",
                "n_total": len(sub),
                "n_intermediate": n_inter,
                "pct_intermediate": 100 * n_inter / len(sub),
                "median_af": sub["caller_af"].median(),
            })

    cn_df = pd.DataFrame(records)
    cn_df.to_csv(DATADIR / "step1_paired_loh_af_by_cn_tier.tsv", sep="\t", index=False)
    print(f"  Saved: {DATADIR / 'step1_paired_loh_af_by_cn_tier.tsv'}")

    # Aggregated plot
    fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
    for idx, cn in enumerate(cn_order):
        ax = axes[idx]
        sub = df_loh_tp[df_loh_tp["cn_tier"] == cn]
        if len(sub) > 0:
            for j, afc in enumerate(af_order):
                n = (sub["af_class"] == afc).sum()
                pct = 100 * n / len(sub)
                ax.bar(j, pct, color=AF_COLORS[afc], alpha=0.8)
                ax.text(j, pct + 1, f"{pct:.1f}%", ha="center", fontsize=9)

        ax.set_xticks(range(3))
        ax.set_xticklabels(af_order, fontsize=8, rotation=15)
        ax.set_title(f"{cn}\n(n={len(sub):,})", fontsize=11, fontweight="bold")
        ax.set_ylabel("% of variants" if idx == 0 else "")

    fig.suptitle("Paired Mode: AF Class Distribution by CN Tier (LOH TP, all samples)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIGDIR / "p05_paired_loh_af_by_cn_tier.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return cn_df


def compute_statistics(df):
    """Compute overall AF class statistics."""
    print("\n=== Computing AF class statistics ===")

    records = []
    for sample in SAMPLE_ORDER:
        ds = df[df["sample"] == sample]
        for loh_status in [True, False]:
            sub_loh = ds[ds["loh_bed_hit"] == loh_status]
            for truth in ["TP", "FP"]:
                sub = sub_loh[sub_loh["truth_label"] == truth]
                if len(sub) == 0:
                    continue
                records.append({
                    "sample": sample,
                    "loh": "LOH" if loh_status else "Non-LOH",
                    "truth": truth,
                    "n_total": len(sub),
                    "n_intermediate": (sub["af_class"] == "Intermediate").sum(),
                    "n_extreme": (sub["af_class"] == "Extreme").sum(),
                    "n_near_half": (sub["af_class"] == "Near-half").sum(),
                    "pct_intermediate": 100 * (sub["af_class"] == "Intermediate").sum() / len(sub),
                    "pct_extreme": 100 * (sub["af_class"] == "Extreme").sum() / len(sub),
                    "pct_near_half": 100 * (sub["af_class"] == "Near-half").sum() / len(sub),
                    "median_af": sub["caller_af"].median(),
                    "mean_af": sub["caller_af"].mean(),
                })

    stats_df = pd.DataFrame(records)
    stats_df.to_csv(DATADIR / "step1_paired_af_class_statistics.tsv", sep="\t", index=False)
    print(f"  Saved: {DATADIR / 'step1_paired_af_class_statistics.tsv'}")

    # Print summary
    loh_tp = stats_df[(stats_df["loh"] == "LOH") & (stats_df["truth"] == "TP")]
    if len(loh_tp) > 0:
        total_n = loh_tp["n_total"].sum()
        total_inter = loh_tp["n_intermediate"].sum()
        print(f"\n  LOH TP total: {total_n:,}, intermediate: {total_inter:,} ({100*total_inter/total_n:.1f}%)")

    return stats_df


def main():
    print("=" * 70)
    print("Step 1: Paired Mode LOH AF Distribution")
    print("=" * 70)

    # Load data
    _, df_paired = load_paired_data()

    # Load LOH.bed and annotate
    print("\nLoading LOH.bed files...")
    loh_beds = load_loh_beds()
    print(f"\nAnnotating {len(df_paired):,} paired variants with LOH.bed...")
    df_paired["loh_bed_hit"] = annotate_loh_bed(df_paired, loh_beds)
    n_hit = df_paired["loh_bed_hit"].sum()
    print(f"  LOH.bed hit: {n_hit:,} ({100*n_hit/len(df_paired):.1f}%)")

    # Generate figures and statistics
    stats = compute_statistics(df_paired)
    fig01_af_distribution(df_paired)
    prop_df = fig02_intermediate_proportion(df_paired)
    fig03_fine_grain(df_paired)
    fig04_af_vs_cm(df_paired)
    cn_df = fig05_af_by_cn_tier(df_paired)

    print("\n" + "=" * 70)
    print("Step 1 complete!")
    print(f"  Figures: {FIGDIR}/p01-p05_*.png")
    print(f"  Data: {DATADIR}/step1_*.tsv")
    print("=" * 70)


if __name__ == "__main__":
    main()
