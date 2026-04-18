#!/usr/bin/env python3
"""G4: Strand Bias & Read Geometry — strand count analysis for germline FP identification.

Uses the newly extracted FAU/FCU/FGU/FTU/RAU/RCU/RGU/RTU fields from G1.
Hypothesis: germline FP are REAL variants → they should have NORMAL strand bias,
making SB a poor discriminator (expected negative result).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    OUTPUT_ROOT, SAMPLE_ORDER, setup_plot_style, save_figure, add_caption,
    compute_auc, ensure_dir,
    COLOR_TP, COLOR_FP,
)

WORKSPACE = OUTPUT_ROOT / "20260401_germline_fp_identification"
G1_DIR = WORKSPACE / "G1_vcf_features"
G4_DIR = ensure_dir(WORKSPACE / "G4_strand_bias")
FIG_DIR = ensure_dir(G4_DIR / "figures")
DATA_DIR = ensure_dir(G4_DIR / "data")


def _print(msg: str):
    print(msg, flush=True)


def load_pass_features() -> pd.DataFrame:
    path = G1_DIR / "all_samples_merged.tsv.gz"
    _print(f"[G4] Loading from {path.name}")
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    for col in ["sb_clairs", "sb_ratio", "sb_fisher_p", "forward_alt", "reverse_alt",
                "forward_ref", "reverse_ref", "af"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    _print(f"  -> {len(df):,} rows")
    return df


def fig01_sb_distribution(df: pd.DataFrame):
    """SB value distribution for PASS TP vs FP."""
    _print("\n[G4] Fig01: SB distribution")
    pass_df = df[df["filter_is_pass"] == 1]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for ax_idx, (col, title) in enumerate([
        ("sb_clairs", "ClairS-TO SB (Fisher p-value)"),
        ("sb_ratio", "SB Ratio (forward_alt / total_alt)"),
        ("sb_fisher_p", "Fisher Exact p-value (2x2)"),
    ]):
        ax = axes[ax_idx]
        for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            vals = pass_df[pass_df["truth_label"] == label][col].dropna()
            if len(vals) > 10:
                vals.plot.kde(ax=ax, label=f"{label} (n={len(vals):,})",
                             color=color, bw_method=0.05, linewidth=1.5)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(col)
        ax.legend(fontsize=8)

    fig.suptitle("G4-Fig01: Strand Bias Distribution (PASS TP vs FP)",
                 fontsize=14, fontweight="bold")
    add_caption(fig, "Germline FP = real variants → expected to have normal SB (overlap with TP)")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G4_fig01_sb_distribution.png")


def fig02_fisher_pvalue(df: pd.DataFrame):
    """Fisher SB p-value distribution."""
    _print("\n[G4] Fig02: Fisher SB p-value")
    pass_df = df[df["filter_is_pass"] == 1]

    fig, ax = plt.subplots(figsize=(10, 5))
    for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
        vals = pass_df[pass_df["truth_label"] == label]["sb_fisher_p"].dropna()
        ax.hist(vals, bins=50, alpha=0.5, color=color, label=f"{label} (n={len(vals):,})",
                density=True, edgecolor="white")

    ax.set_xlabel("Fisher Exact Test p-value")
    ax.set_ylabel("Density")
    ax.set_title("G4-Fig02: Fisher SB p-value Distribution", fontsize=13, fontweight="bold")
    ax.legend()
    add_caption(fig, "Uniform distribution = no strand bias | PASS variants, 7 samples")
    save_figure(fig, FIG_DIR / "G4_fig02_fisher_pvalue.png")


def fig03_sb_ratio(df: pd.DataFrame):
    """SB ratio distribution."""
    _print("\n[G4] Fig03: SB ratio")
    pass_df = df[df["filter_is_pass"] == 1]

    fig, ax = plt.subplots(figsize=(10, 5))
    for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
        vals = pass_df[pass_df["truth_label"] == label]["sb_ratio"].dropna()
        if len(vals) > 10:
            vals.plot.kde(ax=ax, label=f"{label} (n={len(vals):,})",
                         color=color, bw_method=0.05, linewidth=2)

    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5, label="Perfect balance")
    ax.set_xlim(0, 1)
    ax.set_xlabel("SB Ratio (forward_alt / total_alt)")
    ax.set_ylabel("Density")
    ax.set_title("G4-Fig03: Strand Bias Ratio Distribution", fontsize=13, fontweight="bold")
    ax.legend()
    save_figure(fig, FIG_DIR / "G4_fig03_sb_ratio.png")


def fig04_strand_scatter(df: pd.DataFrame):
    """Forward vs reverse alt count scatter."""
    _print("\n[G4] Fig04: Strand scatter")
    pass_df = df[df["filter_is_pass"] == 1]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax_idx, (label, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        ax = axes[ax_idx]
        sub = pass_df[pass_df["truth_label"] == label]
        ax.scatter(sub["forward_alt"], sub["reverse_alt"],
                  alpha=0.05, s=2, c=color, rasterized=True)
        max_val = max(sub["forward_alt"].quantile(0.99), sub["reverse_alt"].quantile(0.99))
        ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="y=x")
        ax.set_xlabel("Forward ALT count")
        ax.set_ylabel("Reverse ALT count")
        ax.set_title(f"{label} (n={len(sub):,})")
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        ax.legend()

    fig.suptitle("G4-Fig04: Forward vs Reverse ALT Read Count",
                 fontsize=14, fontweight="bold")
    add_caption(fig, "Points near y=x → balanced strand support")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_figure(fig, FIG_DIR / "G4_fig04_strand_scatter.png")


def fig05_sb_threshold_sweep(df: pd.DataFrame):
    """SB filter threshold sweep."""
    _print("\n[G4] Fig05: SB threshold sweep")
    from sklearn.metrics import precision_recall_curve

    pass_df = df[df["filter_is_pass"] == 1].copy()
    y_true = (pass_df["truth_label"] == "FP").astype(int).values

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax_idx, (col, label) in enumerate([
        ("sb_clairs", "ClairS SB"),
        ("sb_ratio", "SB Ratio"),
    ]):
        ax = axes[ax_idx]
        score = pass_df[col].values
        mask = np.isfinite(score)
        if mask.sum() < 100:
            continue

        # For SB ratio: distance from 0.5 as "artifact score"
        if col == "sb_ratio":
            s = np.abs(score[mask] - 0.5) * 2  # 0 at 0.5, 1 at extremes
        else:
            s = 1 - score[mask]  # Lower p-value = more biased

        precision, recall, _ = precision_recall_curve(y_true[mask], s)
        auc = compute_auc(y_true[mask], s)
        ax.plot(recall, precision, linewidth=2, label=f"AUC={auc:.3f}")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_title(f"{label}")
        ax.legend()

    fig.suptitle("G4-Fig05: Strand Bias Filter Threshold Sweep",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G4_fig05_sb_threshold_sweep.png")


def fig06_per_sample(df: pd.DataFrame):
    """Per-sample SB comparison."""
    _print("\n[G4] Fig06: Per-sample SB")
    pass_df = df[df["filter_is_pass"] == 1]

    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    axes = axes.flatten()

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = pass_df[pass_df["sample"] == sample]
        for label, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            vals = sub[sub["truth_label"] == label]["sb_ratio"].dropna()
            if len(vals) > 10:
                vals.plot.kde(ax=ax, label=f"{label}", color=color, bw_method=0.05)
        ax.set_xlim(0, 1)
        ax.set_title(sample, fontsize=10)
        ax.legend(fontsize=7)

    if len(SAMPLE_ORDER) < 8:
        axes[-1].set_visible(False)

    fig.suptitle("G4-Fig06: Per-Sample SB Ratio Distribution",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, FIG_DIR / "G4_fig06_per_sample_sb.png")


def fig07_sb_roc(df: pd.DataFrame):
    """ROC for various SB metrics."""
    _print("\n[G4] Fig07: SB ROC")
    from sklearn.metrics import roc_curve

    pass_df = df[df["filter_is_pass"] == 1].copy()
    y_true = (pass_df["truth_label"] == "FP").astype(int).values

    scores = {
        "SB ratio deviation": np.abs(pass_df["sb_ratio"].values - 0.5) * 2,
        "1 - ClairS SB": 1 - pass_df["sb_clairs"].values,
        "1 - Fisher p": 1 - pass_df["sb_fisher_p"].values,
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    for name, score in scores.items():
        mask = np.isfinite(score)
        auc = compute_auc(y_true[mask], score[mask])
        fpr, tpr, _ = roc_curve(y_true[mask], score[mask])
        ax.plot(fpr, tpr, label=f"{name} (AUC={auc:.3f})", linewidth=1.5)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.set_title("G4-Fig07: Strand Bias ROC", fontsize=13, fontweight="bold")
    ax.legend()
    add_caption(fig, "PASS variants, 7 samples | Expected: AUC ≈ 0.50 (no discrimination)")
    save_figure(fig, FIG_DIR / "G4_fig07_sb_roc.png")


def fig08_sb_by_fp_subtype(df: pd.DataFrame):
    """SB by FP subtype (AF-based germline vs artifact)."""
    _print("\n[G4] Fig08: SB by FP subtype")
    pass_df = df[(df["filter_is_pass"] == 1) & (df["truth_label"] == "FP")].copy()

    # Classify FP subtypes by AF
    pass_df["fp_subtype"] = "unknown"
    pass_df.loc[pass_df["af"] >= 0.8, "fp_subtype"] = "likely_hom_germline"
    pass_df.loc[(pass_df["af"] >= 0.35) & (pass_df["af"] <= 0.65), "fp_subtype"] = "likely_het_germline"
    pass_df.loc[pass_df["af"] < 0.2, "fp_subtype"] = "likely_artifact"
    pass_df.loc[(pass_df["af"] >= 0.2) & (pass_df["af"] < 0.35), "fp_subtype"] = "ambiguous"

    fig, ax = plt.subplots(figsize=(10, 6))
    subtypes = ["likely_hom_germline", "likely_het_germline", "ambiguous", "likely_artifact"]
    colors = ["#E53935", "#FF7043", "#FFA726", "#66BB6A"]

    data_for_box = []
    labels_for_box = []
    for st in subtypes:
        vals = pass_df[pass_df["fp_subtype"] == st]["sb_ratio"].dropna()
        data_for_box.append(vals.values)
        labels_for_box.append(f"{st}\n(n={len(vals):,})")

    # Also add TP for reference
    tp_sb = df[(df["filter_is_pass"] == 1) & (df["truth_label"] == "TP")]["sb_ratio"].dropna()
    data_for_box.append(tp_sb.values)
    labels_for_box.append(f"TP (ref)\n(n={len(tp_sb):,})")
    colors.append(COLOR_TP)

    bp = ax.boxplot(data_for_box, labels=labels_for_box, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_ylabel("SB Ratio")
    ax.set_title("G4-Fig08: SB Ratio by FP Subtype (AF-based classification)",
                 fontsize=13, fontweight="bold")
    add_caption(fig, "Germline FP should have SB ratio near 0.5 (same as TP)")
    save_figure(fig, FIG_DIR / "G4_fig08_sb_by_fp_subtype.png")


def main():
    setup_plot_style()
    _print("=" * 70)
    _print("G4: Strand Bias & Read Geometry")
    _print("=" * 70)

    df = load_pass_features()

    fig01_sb_distribution(df)
    fig02_fisher_pvalue(df)
    fig03_sb_ratio(df)
    fig04_strand_scatter(df)
    fig05_sb_threshold_sweep(df)
    fig06_per_sample(df)
    fig07_sb_roc(df)
    fig08_sb_by_fp_subtype(df)

    _print(f"\n[G4] DONE. Output: {G4_DIR}")


if __name__ == "__main__":
    main()
