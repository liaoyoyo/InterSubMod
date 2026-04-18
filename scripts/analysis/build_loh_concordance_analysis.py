#!/usr/bin/env python3
"""LOH Definition Concordance Analysis: HP Imbalance vs LOH.bed

Step 1.1: HP Imbalance (Potential_LOH) vs LOH.bed (to_loh_bed_hit) concordance
Step 1.2: Four-quadrant TP/FP distribution analysis

Data: Master dataset (748K rows), TO mode only, LongPhase-TO baseline.
Terminology: ISM Potential_LOH = "HP Imbalance" (HP_Ratio < 0.1 or > 0.9).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure,
    compute_auc, compute_effect_size, compute_enrichment,
    bootstrap_ci, write_round_context, write_feature_statistics,
    ensure_dir, safe_div, encode_truth_binary,
    mannwhitney_test, chi2_test, format_p, format_ci,
    SAMPLE_ORDER, MODE_ORDER, OUTPUT_ROOT,
    TRUTH_PALETTE, MODE_PALETTE,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from sklearn.metrics import cohen_kappa_score, roc_curve, auc as sklearn_auc
from scipy import stats as sp_stats

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260406_loh_concordance")

# Quadrant colors
QUAD_PALETTE = {
    "Q1_both": "#D32F2F",
    "Q2_ism_only": "#FF9800",
    "Q3_bed_only": "#1976D2",
    "Q4_neither": "#4CAF50",
}
QUAD_ORDER = ["Q1_both", "Q2_ism_only", "Q3_bed_only", "Q4_neither"]
QUAD_LABELS = {
    "Q1_both": "Q1: Both LOH",
    "Q2_ism_only": "Q2: ISM-only",
    "Q3_bed_only": "Q3: LOH.bed-only",
    "Q4_neither": "Q4: Neither",
}


def assign_quadrants(df: pd.DataFrame) -> pd.DataFrame:
    """Assign LOH quadrants based on HP Imbalance and LOH.bed."""
    mask_ism = df["Potential_LOH"].astype(str).str.lower().isin(["true", "1", "yes"])
    mask_bed = df["to_loh_bed_hit"].astype(str).str.lower().isin(["true", "1", "yes"])
    df["loh_quadrant"] = "Q4_neither"
    df.loc[mask_ism & mask_bed, "loh_quadrant"] = "Q1_both"
    df.loc[mask_ism & ~mask_bed, "loh_quadrant"] = "Q2_ism_only"
    df.loc[~mask_ism & mask_bed, "loh_quadrant"] = "Q3_bed_only"
    df["hp_imbalance"] = mask_ism
    df["loh_bed"] = mask_bed
    return df


# ===== Step 1.1: Concordance Analysis =====

def compute_concordance_stats(df_to: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample concordance statistics."""
    rows = []
    for sample in SAMPLE_ORDER:
        ds = df_to[df_to["sample"] == sample]
        if len(ds) == 0:
            continue
        ism = ds["hp_imbalance"].values
        bed = ds["loh_bed"].values

        # Confusion matrix counts
        tp_tp = int(np.sum(ism & bed))       # Q1
        tp_fn = int(np.sum(ism & ~bed))      # Q2
        fn_tp = int(np.sum(~ism & bed))      # Q3
        tn_tn = int(np.sum(~ism & ~bed))     # Q4
        total = len(ds)

        # Kappa
        try:
            kappa = cohen_kappa_score(bed.astype(int), ism.astype(int))
        except Exception:
            kappa = float("nan")

        # Jaccard
        union = tp_tp + tp_fn + fn_tp
        jaccard = safe_div(tp_tp, union) if union > 0 else float("nan")

        # Sensitivity (ISM detecting LOH.bed as reference)
        sensitivity = safe_div(tp_tp, tp_tp + fn_tp)
        specificity = safe_div(tn_tn, tn_tn + tp_fn)

        # Rates
        ism_rate = safe_div(tp_tp + tp_fn, total)
        bed_rate = safe_div(tp_tp + fn_tp, total)

        rows.append({
            "sample": sample,
            "n_total": total,
            "Q1_both": tp_tp,
            "Q2_ism_only": tp_fn,
            "Q3_bed_only": fn_tp,
            "Q4_neither": tn_tn,
            "kappa": round(kappa, 4),
            "jaccard": round(jaccard, 4),
            "sensitivity": round(sensitivity, 4),
            "specificity": round(specificity, 4),
            "hp_imbalance_rate": round(ism_rate, 4),
            "loh_bed_rate": round(bed_rate, 4),
        })
    return pd.DataFrame(rows)


def fig01_confusion_matrices(df_to: pd.DataFrame) -> None:
    """F1: 2x2 confusion matrix heatmap per sample (TO only)."""
    setup_plot_style()
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        ds = df_to[df_to["sample"] == sample]
        ism = ds["hp_imbalance"].values
        bed = ds["loh_bed"].values

        # Build 2x2: rows = ISM (True/False), cols = LOH.bed (True/False)
        cm = np.array([
            [int(np.sum(ism & bed)), int(np.sum(ism & ~bed))],
            [int(np.sum(~ism & bed)), int(np.sum(~ism & ~bed))],
        ])
        cm_pct = cm / cm.sum() * 100

        sns.heatmap(
            cm, annot=True, fmt="d", cmap="YlOrRd", ax=ax,
            xticklabels=["LOH.bed+", "LOH.bed-"],
            yticklabels=["HP-Imb+", "HP-Imb-"],
            cbar=False,
        )
        # Add percentage annotations
        for ri in range(2):
            for ci in range(2):
                ax.text(ci + 0.5, ri + 0.75, f"({cm_pct[ri, ci]:.1f}%)",
                        ha="center", va="center", fontsize=8, color="gray")

        ax.set_title(f"{sample}\n(n={len(ds):,})", fontsize=10)

    # Hide unused subplot
    axes[-1].set_visible(False)

    fig.suptitle("HP Imbalance vs LOH.bed Confusion Matrix (TO mode)", fontsize=14, y=1.02)
    save_figure(fig, OUT_DIR / "01_confusion_matrix.png")


def fig02_hp_ratio_distribution(df_to: pd.DataFrame) -> None:
    """F2: HP_Ratio distribution colored by LOH.bed status."""
    setup_plot_style()
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        ds = df_to[df_to["sample"] == sample]

        for label, color, mask_col in [
            ("In LOH.bed", "#D32F2F", True),
            ("Not in LOH.bed", "#1976D2", False),
        ]:
            subset = ds[ds["loh_bed"] == mask_col]["HP_Ratio"].dropna()
            if len(subset) > 0:
                ax.hist(subset, bins=50, alpha=0.6, label=f"{label} (n={len(subset):,})",
                        color=color, density=True)

        # Mark ISM threshold
        ax.axvline(0.1, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.axvline(0.9, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.set_title(sample, fontsize=10)
        ax.set_xlabel("HP_Ratio")
        ax.set_ylabel("Density")
        ax.legend(fontsize=7)

    axes[-1].set_visible(False)
    fig.suptitle("HP_Ratio Distribution by LOH.bed Status (TO mode)\nDashed lines: HP Imbalance threshold (0.1, 0.9)",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "02_hp_ratio_distribution.png")


def fig03_kappa_barchart(stats_df: pd.DataFrame) -> None:
    """F3: Cohen's kappa bar chart across samples."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(10, 5))

    bars = ax.bar(stats_df["sample"], stats_df["kappa"],
                  color=["#1976D2" if k > 0.4 else "#FF9800" if k > 0.2 else "#D32F2F"
                         for k in stats_df["kappa"]])

    for bar, val in zip(bars, stats_df["kappa"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{val:.3f}", ha="center", va="bottom", fontsize=9)

    ax.set_ylabel("Cohen's Kappa")
    ax.set_title("HP Imbalance vs LOH.bed Agreement (Cohen's Kappa, TO mode)")
    ax.axhline(0.4, color="gray", linestyle=":", alpha=0.5, label="Moderate (0.4)")
    ax.axhline(0.6, color="gray", linestyle="--", alpha=0.5, label="Substantial (0.6)")
    ax.legend(fontsize=8)
    ax.set_ylim(0, max(stats_df["kappa"].max() * 1.2, 0.8))
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "03_kappa_barchart.png")


def fig04_hp_ratio_roc(df_to: pd.DataFrame) -> None:
    """F4: HP_Ratio ROC curve (predicting LOH.bed) per sample."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(8, 8))

    colors = plt.cm.Set2(np.linspace(0, 1, len(SAMPLE_ORDER)))

    for i, sample in enumerate(SAMPLE_ORDER):
        ds = df_to[df_to["sample"] == sample]
        y_true = ds["loh_bed"].astype(int).values
        # HP imbalance = extreme HP_Ratio, so use distance from 0.5
        y_score = (ds["HP_Ratio"] - 0.5).abs().values

        mask = np.isfinite(y_score)
        y_true, y_score = y_true[mask], y_score[mask]

        if len(np.unique(y_true)) < 2:
            continue

        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = sklearn_auc(fpr, tpr)
        ax.plot(fpr, tpr, color=colors[i], label=f"{sample} (AUC={roc_auc:.3f})", linewidth=1.5)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="Random")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("HP_Ratio as LOH.bed Predictor (TO mode)\n|HP_Ratio - 0.5| as score")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_aspect("equal")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "04_hp_ratio_roc.png")


# ===== Step 1.2: Four-Quadrant TP/FP Analysis =====

def compute_quadrant_stats(df_to: pd.DataFrame) -> pd.DataFrame:
    """Compute TP/FP statistics per quadrant per sample."""
    rows = []
    for sample in SAMPLE_ORDER:
        ds = df_to[df_to["sample"] == sample]
        tp_total = int((ds["truth_label"] == "TP").sum())
        fp_total = int((ds["truth_label"] == "FP").sum())

        for quad in QUAD_ORDER:
            qd = ds[ds["loh_quadrant"] == quad]
            tp_count = int((qd["truth_label"] == "TP").sum())
            fp_count = int((qd["truth_label"] == "FP").sum())
            total_q = len(qd)
            fp_rate = safe_div(fp_count, total_q)

            enr = compute_enrichment(tp_count, fp_count, tp_total, fp_total)

            rows.append({
                "sample": sample,
                "quadrant": quad,
                "n_total": total_q,
                "tp_count": tp_count,
                "fp_count": fp_count,
                "fp_rate": round(fp_rate, 4),
                "fp_enrichment": round(enr["enrichment"], 3),
                "fisher_p": enr["p_value"],
                "tp_pct_of_total": round(safe_div(tp_count, tp_total) * 100, 2),
                "fp_pct_of_total": round(safe_div(fp_count, fp_total) * 100, 2),
            })
    return pd.DataFrame(rows)


def fig05_quadrant_stacked_bar(df_to: pd.DataFrame) -> None:
    """F5: Stacked bar of TP/FP per quadrant, per sample."""
    setup_plot_style()
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        ds = df_to[df_to["sample"] == sample]

        tp_counts = []
        fp_counts = []
        for quad in QUAD_ORDER:
            qd = ds[ds["loh_quadrant"] == quad]
            tp_counts.append(int((qd["truth_label"] == "TP").sum()))
            fp_counts.append(int((qd["truth_label"] == "FP").sum()))

        x = np.arange(len(QUAD_ORDER))
        width = 0.6
        ax.bar(x, tp_counts, width, label="TP", color=COLOR_TP, alpha=0.85)
        ax.bar(x, fp_counts, width, bottom=tp_counts, label="FP", color=COLOR_FP, alpha=0.85)

        # Add count labels
        for j in range(len(QUAD_ORDER)):
            total = tp_counts[j] + fp_counts[j]
            if total > 0:
                ax.text(j, total + max(tp_counts) * 0.02,
                        f"{total:,}", ha="center", va="bottom", fontsize=7)

        ax.set_xticks(x)
        ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4"], fontsize=8)
        ax.set_title(f"{sample}", fontsize=10)
        ax.set_ylabel("Count")
        if i == 0:
            ax.legend(fontsize=8)

    axes[-1].set_visible(False)
    fig.suptitle("TP/FP Distribution by LOH Quadrant (TO mode)\n"
                 "Q1=Both LOH, Q2=ISM-only, Q3=LOH.bed-only, Q4=Neither",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "05_quadrant_tp_fp_stacked.png")


def fig06_fp_rate_heatmap(quad_stats: pd.DataFrame) -> None:
    """F6: FP rate heatmap (sample x quadrant)."""
    setup_plot_style()

    pivot = quad_stats.pivot(index="sample", columns="quadrant", values="fp_rate")
    pivot = pivot.reindex(index=SAMPLE_ORDER, columns=QUAD_ORDER)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="YlOrRd", ax=ax,
                xticklabels=[QUAD_LABELS[q] for q in QUAD_ORDER],
                linewidths=0.5, vmin=0, vmax=1)
    ax.set_title("FP Rate by Sample x LOH Quadrant (TO mode)", fontsize=12)
    ax.set_ylabel("Sample")
    ax.set_xlabel("LOH Quadrant")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "06_fp_rate_heatmap.png")


def fig07_fp_enrichment(quad_stats: pd.DataFrame) -> None:
    """F7: FP enrichment ratio bar chart."""
    setup_plot_style()
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        ss = quad_stats[quad_stats["sample"] == sample]

        colors = [QUAD_PALETTE[q] for q in QUAD_ORDER]
        bars = ax.bar(range(4), ss.set_index("quadrant").reindex(QUAD_ORDER)["fp_enrichment"],
                      color=colors)

        ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4"], fontsize=8)
        ax.set_title(sample, fontsize=10)
        ax.set_ylabel("FP Enrichment")

        # Annotate values
        for bar, val in zip(bars, ss.set_index("quadrant").reindex(QUAD_ORDER)["fp_enrichment"]):
            if np.isfinite(val):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                        f"{val:.2f}", ha="center", va="bottom", fontsize=8)

    axes[-1].set_visible(False)
    fig.suptitle("FP Enrichment Ratio by LOH Quadrant (TO mode)\n"
                 "Enrichment > 1 = FP-enriched; dashed line = neutral",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "07_fp_enrichment.png")


# ===== Main =====

def main():
    print("=" * 70)
    print("LOH Definition Concordance Analysis")
    print("Step 1.1: HP Imbalance vs LOH.bed concordance")
    print("Step 1.2: Four-quadrant TP/FP distribution")
    print("=" * 70)

    # Load data
    df = load_master_dataset()
    df_to = df[df["mode"] == "to"].copy()
    print(f"\nTO mode: {len(df_to):,} rows")

    # Assign quadrants
    df_to = assign_quadrants(df_to)

    # Verify HP Imbalance rate (expected ~60% from O3)
    imb_rate = df_to["hp_imbalance"].mean()
    print(f"\n[Verification] HP Imbalance rate (TO): {imb_rate:.4f} ({imb_rate*100:.1f}%)")
    print(f"[Verification] Expected: ~60% (from O3)")

    bed_rate = df_to["loh_bed"].mean()
    print(f"[Verification] LOH.bed rate (TO): {bed_rate:.4f} ({bed_rate*100:.1f}%)")

    # Quadrant distribution
    print("\n[Quadrant distribution (TO)]")
    for q in QUAD_ORDER:
        n = (df_to["loh_quadrant"] == q).sum()
        pct = n / len(df_to) * 100
        print(f"  {QUAD_LABELS[q]}: {n:,} ({pct:.1f}%)")

    # ===== Step 1.1 =====
    print("\n" + "=" * 50)
    print("Step 1.1: Concordance Analysis")
    print("=" * 50)

    concordance_stats = compute_concordance_stats(df_to)
    print("\nPer-sample concordance:")
    print(concordance_stats[["sample", "kappa", "jaccard", "sensitivity", "specificity",
                             "hp_imbalance_rate", "loh_bed_rate"]].to_string(index=False))

    # Save statistics
    concordance_stats.to_csv(OUT_DIR / "concordance_statistics.tsv", sep="\t", index=False)

    # Generate figures
    fig01_confusion_matrices(df_to)
    fig02_hp_ratio_distribution(df_to)
    fig03_kappa_barchart(concordance_stats)
    fig04_hp_ratio_roc(df_to)

    # ===== Step 1.2 =====
    print("\n" + "=" * 50)
    print("Step 1.2: Four-Quadrant TP/FP Analysis")
    print("=" * 50)

    quad_stats = compute_quadrant_stats(df_to)

    # Summary across all samples
    print("\nQuadrant TP/FP summary (all samples combined):")
    for q in QUAD_ORDER:
        qs = quad_stats[quad_stats["quadrant"] == q]
        tp_sum = qs["tp_count"].sum()
        fp_sum = qs["fp_count"].sum()
        total = tp_sum + fp_sum
        fp_r = safe_div(fp_sum, total)
        print(f"  {QUAD_LABELS[q]}: TP={tp_sum:,}, FP={fp_sum:,}, "
              f"FP_rate={fp_r:.3f}, median_enrichment={qs['fp_enrichment'].median():.2f}")

    # Verify: quadrant totals = sample totals
    for sample in SAMPLE_ORDER:
        ds = df_to[df_to["sample"] == sample]
        qs = quad_stats[quad_stats["sample"] == sample]
        assert qs["n_total"].sum() == len(ds), f"Quadrant total mismatch for {sample}!"
    print("\n[Verification] Quadrant totals match sample totals: PASS")

    # Save statistics
    quad_stats.to_csv(OUT_DIR / "quadrant_statistics.tsv", sep="\t", index=False)

    fig05_quadrant_stacked_bar(df_to)
    fig06_fp_rate_heatmap(quad_stats)
    fig07_fp_enrichment(quad_stats)

    # Feature statistics
    write_feature_statistics(
        df_to,
        ["HP_Ratio", "Potential_LOH", "to_loh_bed_hit", "HP1FamilyN", "HP2FamilyN",
         "NumReads", "Coverage_Multiple", "caller_af", "truth_label"],
        OUT_DIR / "feature_statistics.tsv",
    )

    # Round context
    write_round_context(
        OUT_DIR,
        observation_id="loh_concordance",
        title="LOH Definition Concordance: HP Imbalance vs LOH.bed",
        description=(
            "Systematic comparison of ISM HP Imbalance (HP_Ratio threshold) "
            "with LongPhase-TO LOH.bed. Includes 2x2 concordance analysis, "
            "Cohen's kappa, ROC curves, and four-quadrant TP/FP distribution."
        ),
        script_path="scripts/analysis/build_loh_concordance_analysis.py",
        data_sources=[{
            "name": "master_dataset",
            "path": str(load_master_dataset.__defaults__),
            "version": "2026-03-27 HP fix rebuild",
            "longphase_version": "baseline",
        }],
        row_count=len(df_to),
        col_count=len(df_to.columns),
        extra={
            "n_figures": 7,
            "n_tables": 3,
            "terminology": "ISM Potential_LOH = HP Imbalance (HP_Ratio < 0.1 or > 0.9)",
        },
    )

    print("\n" + "=" * 70)
    print(f"DONE. Output: {OUT_DIR}")
    print(f"Figures: 7, Tables: 3")
    print("=" * 70)


if __name__ == "__main__":
    main()
