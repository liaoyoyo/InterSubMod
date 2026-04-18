#!/usr/bin/env python3
"""O13: Cross-Region Methylation Correlation — Quick Validation.

Tests whether nearby SNVs share reads and show correlated methylation patterns.
Core question: Do TP-TP pairs show different cross-region methylation correlation
than FP-FP pairs?

Hypothesis:
  - Somatic (TP) variants from the same subclone → correlated methylation across loci
  - Germline (FP) variants → independent methylation (each driven by own mQTL)

Method:
  1. Find SNV pairs within 50kb on same chromosome (same sample)
  2. For each pair, find shared reads (same read_name in both regions)
  3. For shared reads, compute per-read mean methylation in each region
  4. Correlate these per-read means across the two regions
  5. Compare: TP-TP pairs vs FP-FP pairs vs TP-FP pairs
"""

from __future__ import annotations

import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from sklearn.metrics import roc_auc_score

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure, add_caption,
    ensure_dir, OUTPUT_ROOT, SAMPLE_ORDER, COLOR_TP, COLOR_FP,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value")

# ===== Config =====
OBS_ID = "O13"
RUN_DATE = datetime.now().strftime("%Y%m%d")
WORKSPACE_NAME = f"{RUN_DATE}_O13_cross_region_correlation"
OUT_DIR = OUTPUT_ROOT / WORKSPACE_NAME

MAX_DISTANCE = 50000  # max distance between SNV pairs (bp)
MIN_SHARED_READS = 5  # minimum shared reads to compute correlation
MAX_PAIRS_PER_SAMPLE = 300  # per truth_pair_type (TP-TP, FP-FP, TP-FP)
SEED = 42


def _find_region_dir(region_root: str, chrom: str, pos: int) -> Optional[Path]:
    """Locate actual region directory."""
    variant_dir = Path(region_root) / chrom / f"{chrom}_{pos}"
    if not variant_dir.exists():
        return None
    for child in variant_dir.iterdir():
        if child.is_dir() and (child / "reads" / "reads.tsv").exists():
            return child
    return None


def load_region_reads(region_dir: Path) -> Optional[pd.DataFrame]:
    """Load reads.tsv and methylation.csv, merge by read_id index."""
    reads_path = region_dir / "reads" / "reads.tsv"
    methyl_path = region_dir / "methylation" / "methylation.csv"
    if not reads_path.exists() or not methyl_path.exists():
        return None
    try:
        reads = pd.read_csv(reads_path, sep="\t")
        methyl = pd.read_csv(methyl_path)
        cpg_cols = [c for c in methyl.columns if c != "read_id"]
        if len(cpg_cols) < 5 or len(methyl) < 5:
            return None
        # Per-read mean methylation
        methyl["methyl_mean"] = methyl[cpg_cols].astype(float).mean(axis=1, skipna=True)
        methyl["methyl_std"] = methyl[cpg_cols].astype(float).std(axis=1, skipna=True)
        methyl["n_cpg_obs"] = methyl[cpg_cols].astype(float).notna().sum(axis=1)
        # Merge
        merged = reads.merge(methyl[["read_id", "methyl_mean", "methyl_std", "n_cpg_obs"]],
                             on="read_id", how="inner")
        return merged
    except Exception:
        return None


def find_nearby_pairs(sample_df: pd.DataFrame, max_dist: int) -> pd.DataFrame:
    """Find pairs of SNVs within max_dist on the same chromosome."""
    pairs = []
    for chrom in sample_df["Chr"].unique():
        chrom_df = sample_df[sample_df["Chr"] == chrom].sort_values("Pos")
        positions = chrom_df["Pos"].values
        labels = chrom_df["truth_label"].values
        region_ids = chrom_df["RegionID"].values
        roots = chrom_df["source_region_root"].values

        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                dist = positions[j] - positions[i]
                if dist > max_dist:
                    break
                pairs.append({
                    "chr": chrom,
                    "pos_A": int(positions[i]),
                    "pos_B": int(positions[j]),
                    "distance": int(dist),
                    "label_A": labels[i],
                    "label_B": labels[j],
                    "region_id_A": region_ids[i],
                    "region_id_B": region_ids[j],
                    "root_A": roots[i],
                    "root_B": roots[j],
                    "pair_type": f"{labels[i]}-{labels[j]}",
                })
    return pd.DataFrame(pairs)


def compute_pair_correlation(
    root_A: str, chrom: str, pos_A: int,
    root_B: str, pos_B: int,
) -> Optional[Dict]:
    """For a pair of regions, find shared reads and compute correlation."""
    dir_A = _find_region_dir(root_A, chrom, pos_A)
    dir_B = _find_region_dir(root_B, chrom, pos_B)
    if dir_A is None or dir_B is None:
        return None

    reads_A = load_region_reads(dir_A)
    reads_B = load_region_reads(dir_B)
    if reads_A is None or reads_B is None:
        return None

    # Find shared reads by read_name
    shared = reads_A[["read_name", "methyl_mean", "hp", "alt_support"]].merge(
        reads_B[["read_name", "methyl_mean", "hp", "alt_support"]],
        on="read_name", suffixes=("_A", "_B"),
    )

    n_shared = len(shared)
    n_reads_A = len(reads_A)
    n_reads_B = len(reads_B)

    if n_shared < MIN_SHARED_READS:
        return {"n_shared": n_shared, "n_reads_A": n_reads_A, "n_reads_B": n_reads_B,
                "correlation": np.nan, "p_value": np.nan, "overlap_rate": 0.0}

    # Drop NaN
    valid = shared.dropna(subset=["methyl_mean_A", "methyl_mean_B"])
    if len(valid) < MIN_SHARED_READS:
        return {"n_shared": n_shared, "n_reads_A": n_reads_A, "n_reads_B": n_reads_B,
                "correlation": np.nan, "p_value": np.nan, "overlap_rate": n_shared / min(n_reads_A, n_reads_B)}

    # Pearson correlation of per-read mean methylation
    r, p = sp_stats.pearsonr(valid["methyl_mean_A"], valid["methyl_mean_B"])

    # HP concordance: same HP in both regions?
    hp_match = (valid["hp_A"] == valid["hp_B"]).mean() if len(valid) > 0 else np.nan

    # ALT support concordance
    alt_match = (valid["alt_support_A"] == valid["alt_support_B"]).mean() if len(valid) > 0 else np.nan

    return {
        "n_shared": n_shared,
        "n_reads_A": n_reads_A,
        "n_reads_B": n_reads_B,
        "n_valid": len(valid),
        "overlap_rate": n_shared / min(n_reads_A, n_reads_B),
        "correlation": float(r),
        "p_value": float(p),
        "hp_concordance": float(hp_match),
        "alt_concordance": float(alt_match),
        "methyl_mean_A": float(valid["methyl_mean_A"].mean()),
        "methyl_mean_B": float(valid["methyl_mean_B"].mean()),
    }


def main():
    print(f"{'=' * 70}")
    print(f"O13: Cross-Region Methylation Correlation — Quick Validation")
    print(f"Date: {RUN_DATE}")
    print(f"Output: {OUT_DIR}")
    print(f"{'=' * 70}")

    ensure_dir(OUT_DIR)

    # Load master dataset — TO mode only
    print("\nLoading master dataset...")
    df = load_master_dataset()
    to_df = df[df["mode"] == "to"].copy()
    print(f"TO regions: {len(to_df):,}")

    all_results = []

    for sample in SAMPLE_ORDER:
        s = to_df[to_df["sample"] == sample]
        if len(s) == 0:
            continue

        print(f"\n--- {sample} (n={len(s)}) ---")

        # Find nearby pairs
        pairs = find_nearby_pairs(s, MAX_DISTANCE)
        if len(pairs) == 0:
            print(f"  No nearby pairs found")
            continue

        # Count pair types
        pair_counts = pairs["pair_type"].value_counts()
        print(f"  Nearby pairs (<{MAX_DISTANCE/1000:.0f}kb): {len(pairs)}")
        for pt, cnt in pair_counts.items():
            print(f"    {pt}: {cnt}")

        # Sample balanced pairs per type
        sampled = []
        for pt in ["TP-TP", "FP-FP", "TP-FP", "FP-TP"]:
            sub = pairs[pairs["pair_type"] == pt]
            if len(sub) > MAX_PAIRS_PER_SAMPLE:
                sub = sub.sample(n=MAX_PAIRS_PER_SAMPLE, random_state=SEED)
            sampled.append(sub)
        sampled = pd.concat(sampled, ignore_index=True)
        print(f"  Sampled: {len(sampled)} pairs")

        # Compute correlations
        n_success = 0
        n_fail = 0
        for _, row in sampled.iterrows():
            result = compute_pair_correlation(
                str(row["root_A"]), str(row["chr"]), int(row["pos_A"]),
                str(row["root_B"]), int(row["pos_B"]),
            )
            if result is None:
                n_fail += 1
                continue

            result["sample"] = sample
            result["chr"] = row["chr"]
            result["pos_A"] = row["pos_A"]
            result["pos_B"] = row["pos_B"]
            result["distance"] = row["distance"]
            result["pair_type"] = row["pair_type"]
            result["label_A"] = row["label_A"]
            result["label_B"] = row["label_B"]
            all_results.append(result)
            n_success += 1

        print(f"  Computed: {n_success} success, {n_fail} failed")

    # Compile results
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_DIR / "cross_region_pairs.tsv", sep="\t", index=False)
    print(f"\nTotal pairs: {len(results_df)}")

    if len(results_df) == 0:
        print("No pairs computed. Exiting.")
        return

    # Normalize pair_type: TP-FP and FP-TP → Mixed
    results_df["pair_class"] = results_df["pair_type"].replace({"FP-TP": "Mixed", "TP-FP": "Mixed"})
    results_df.loc[results_df["pair_type"] == "TP-TP", "pair_class"] = "TP-TP"
    results_df.loc[results_df["pair_type"] == "FP-FP", "pair_class"] = "FP-FP"

    # ===== Analysis =====
    print("\n===== Cross-Region Correlation Summary =====")

    # 1. Overall overlap rate
    valid = results_df[results_df["n_shared"] >= MIN_SHARED_READS]
    print(f"\nPairs with >= {MIN_SHARED_READS} shared reads: {len(valid)} / {len(results_df)}")
    print(f"Mean overlap_rate: {results_df['overlap_rate'].mean():.3f}")
    print(f"Median shared reads: {results_df['n_shared'].median():.0f}")

    # 2. Correlation by pair_class
    print(f"\nCorrelation by pair class:")
    for pc in ["TP-TP", "FP-FP", "Mixed"]:
        sub = valid[valid["pair_class"] == pc]
        if len(sub) > 0:
            corrs = sub["correlation"].dropna()
            print(f"  {pc}: n={len(corrs)}, median_r={corrs.median():.4f}, "
                  f"mean_r={corrs.mean():.4f}, std={corrs.std():.4f}")

    # 3. Per-sample breakdown
    print(f"\nPer-sample correlation (TP-TP vs FP-FP):")
    summary_rows = []
    for sample in SAMPLE_ORDER:
        sv = valid[valid["sample"] == sample]
        if len(sv) == 0:
            continue
        for pc in ["TP-TP", "FP-FP", "Mixed"]:
            sub = sv[sv["pair_class"] == pc]
            corrs = sub["correlation"].dropna()
            if len(corrs) >= 5:
                summary_rows.append({
                    "sample": sample,
                    "pair_class": pc,
                    "n_pairs": len(corrs),
                    "median_r": float(corrs.median()),
                    "mean_r": float(corrs.mean()),
                    "std_r": float(corrs.std()),
                    "median_overlap": float(sub["overlap_rate"].median()),
                    "median_shared": float(sub["n_shared"].median()),
                    "median_hp_concordance": float(sub["hp_concordance"].dropna().median()) if sub["hp_concordance"].notna().any() else np.nan,
                })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "per_sample_summary.tsv", sep="\t", index=False)
    print(summary_df.to_string(index=False))

    # 4. Statistical test: TP-TP vs FP-FP correlation difference
    print(f"\n===== Statistical Tests =====")
    tp_corrs = valid[valid["pair_class"] == "TP-TP"]["correlation"].dropna()
    fp_corrs = valid[valid["pair_class"] == "FP-FP"]["correlation"].dropna()
    if len(tp_corrs) >= 10 and len(fp_corrs) >= 10:
        u_stat, u_p = sp_stats.mannwhitneyu(tp_corrs, fp_corrs, alternative="two-sided")
        d = (tp_corrs.mean() - fp_corrs.mean()) / np.sqrt((tp_corrs.var() + fp_corrs.var()) / 2)
        print(f"  TP-TP vs FP-FP correlation:")
        print(f"    TP-TP: n={len(tp_corrs)}, median={tp_corrs.median():.4f}, mean={tp_corrs.mean():.4f}")
        print(f"    FP-FP: n={len(fp_corrs)}, median={fp_corrs.median():.4f}, mean={fp_corrs.mean():.4f}")
        print(f"    Mann-Whitney p={u_p:.4e}, Cohen's d={d:.4f}")

    # 5. Distance effect
    print(f"\n===== Distance Effect =====")
    for dist_bin in [(0, 10000, "<10kb"), (10000, 20000, "10-20kb"), (20000, 50000, "20-50kb")]:
        lo, hi, label = dist_bin
        sub = valid[(valid["distance"] >= lo) & (valid["distance"] < hi)]
        if len(sub) >= 10:
            print(f"  {label}: n={len(sub)}, median_r={sub['correlation'].median():.4f}, "
                  f"median_shared={sub['n_shared'].median():.0f}")

    # ===== Figures =====
    print(f"\n===== Generating Figures =====")
    setup_plot_style()

    # Fig 1: Correlation distribution by pair class
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel A: Histogram
    ax = axes[0]
    for pc, color in [("TP-TP", COLOR_TP), ("FP-FP", COLOR_FP), ("Mixed", "#9E9E9E")]:
        sub = valid[valid["pair_class"] == pc]["correlation"].dropna()
        if len(sub) > 0:
            ax.hist(sub, bins=30, alpha=0.5, color=color, label=f"{pc} (n={len(sub)})", density=True)
    ax.set_xlabel("Pearson r (cross-region methylation)")
    ax.set_ylabel("Density")
    ax.set_title("A: Cross-Region Correlation Distribution")
    ax.legend(fontsize=8)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)

    # Panel B: Correlation vs distance
    ax = axes[1]
    for pc, color, marker in [("TP-TP", COLOR_TP, "o"), ("FP-FP", COLOR_FP, "s")]:
        sub = valid[valid["pair_class"] == pc]
        if len(sub) > 0:
            ax.scatter(sub["distance"] / 1000, sub["correlation"],
                       alpha=0.15, s=8, c=color, marker=marker, label=pc)
    ax.set_xlabel("Distance (kb)")
    ax.set_ylabel("Pearson r")
    ax.set_title("B: Correlation vs Distance")
    ax.legend(fontsize=8, markerscale=3)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)

    # Panel C: Per-sample median correlation
    ax = axes[2]
    if len(summary_df) > 0:
        for pc, color in [("TP-TP", COLOR_TP), ("FP-FP", COLOR_FP)]:
            sub = summary_df[summary_df["pair_class"] == pc]
            if len(sub) > 0:
                ax.barh(
                    [f"{r['sample']}" for _, r in sub.iterrows()],
                    sub["median_r"],
                    color=color, alpha=0.7, label=pc, height=0.35,
                )
        ax.set_xlabel("Median Pearson r")
        ax.set_title("C: Per-Sample Median Correlation")
        ax.legend(fontsize=8)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)

    fig.suptitle("O13: Cross-Region Methylation Correlation (Shared Reads)", fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig01_cross_region_correlation.png")

    # Fig 2: Shared read overlap statistics
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.hist(results_df["n_shared"], bins=50, color="#2196F3", alpha=0.7)
    ax.set_xlabel("Number of shared reads")
    ax.set_ylabel("Count")
    ax.set_title("A: Shared Read Count Distribution")
    ax.axvline(MIN_SHARED_READS, color="red", linestyle=":", label=f"min={MIN_SHARED_READS}")
    ax.legend()

    ax = axes[1]
    ax.hist(results_df["overlap_rate"], bins=50, color="#4CAF50", alpha=0.7)
    ax.set_xlabel("Overlap rate (shared / min(n_A, n_B))")
    ax.set_ylabel("Count")
    ax.set_title("B: Read Overlap Rate")

    ax = axes[2]
    for pc, color in [("TP-TP", COLOR_TP), ("FP-FP", COLOR_FP)]:
        sub = valid[valid["pair_class"] == pc]
        if len(sub) > 0:
            ax.scatter(sub["n_shared"], sub["correlation"],
                       alpha=0.15, s=8, c=color, label=pc)
    ax.set_xlabel("Number of shared reads")
    ax.set_ylabel("Pearson r")
    ax.set_title("C: Correlation vs Shared Read Count")
    ax.legend(fontsize=8, markerscale=3)

    fig.suptitle("O13: Shared Read Statistics", fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig02_shared_read_stats.png")

    # Fig 3: HP and ALT concordance
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    for pc, color in [("TP-TP", COLOR_TP), ("FP-FP", COLOR_FP), ("Mixed", "#9E9E9E")]:
        sub = valid[valid["pair_class"] == pc]["hp_concordance"].dropna()
        if len(sub) > 0:
            ax.hist(sub, bins=20, alpha=0.5, color=color, label=f"{pc} (n={len(sub)})", density=True)
    ax.set_xlabel("HP Concordance (fraction same HP in both regions)")
    ax.set_ylabel("Density")
    ax.set_title("A: HP Tag Concordance Across Regions")
    ax.legend(fontsize=8)

    ax = axes[1]
    for pc, color in [("TP-TP", COLOR_TP), ("FP-FP", COLOR_FP), ("Mixed", "#9E9E9E")]:
        sub = valid[valid["pair_class"] == pc]["alt_concordance"].dropna()
        if len(sub) > 0:
            ax.hist(sub, bins=20, alpha=0.5, color=color, label=f"{pc} (n={len(sub)})", density=True)
    ax.set_xlabel("ALT Support Concordance")
    ax.set_ylabel("Density")
    ax.set_title("B: ALT/REF Concordance Across Regions")
    ax.legend(fontsize=8)

    fig.suptitle("O13: Cross-Region Read Label Concordance", fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "fig03_label_concordance.png")

    # ===== Quick verdict =====
    print(f"\n{'=' * 70}")
    if len(tp_corrs) >= 10 and len(fp_corrs) >= 10:
        delta = tp_corrs.median() - fp_corrs.median()
        if abs(delta) > 0.05 and u_p < 0.01:
            print(f"O13 QUICK VERDICT: SIGNAL DETECTED (delta_r={delta:.4f}, p={u_p:.4e})")
        else:
            print(f"O13 QUICK VERDICT: NO SIGNAL (delta_r={delta:.4f}, p={u_p:.4e})")
    else:
        print(f"O13 QUICK VERDICT: INSUFFICIENT DATA")
    print(f"Output: {OUT_DIR}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
