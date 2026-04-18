#!/usr/bin/env python3
"""
O08: Sample-Level Heterogeneity Observation
============================================

Systematic observation of feature variability across 7 samples x 2 modes,
identifying which features differ most across samples and whether
cross-sample modelling is viable or needs per-sample calibration.

Part of the O1-O10 systematic observation framework for Phase 1A ML baseline.

Generated: 2026-03-31
"""

from __future__ import annotations

import json
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import sys
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure, add_caption,
    compute_auc, compute_effect_size, mannwhitney_test, chi2_test,
    bootstrap_ci, compute_enrichment, format_p, format_ci,
    encode_truth_binary, ensure_dir, safe_div,
    write_round_context, write_feature_statistics, markdown_table,
    OUTPUT_ROOT, SAMPLE_ORDER, MODE_ORDER,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO, TRUTH_PALETTE, MODE_PALETTE,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OBSERVATION_ID = "O08"
OBSERVATION_TITLE = "Sample-Level Heterogeneity Observation"
SCRIPT_PATH = "scripts/analysis/build_observation_O08_sample_heterogeneity.py"

OUT_DIR = OUTPUT_ROOT / "20260401_O08_sample_heterogeneity"
FIG_DIR = OUT_DIR
DATA_DIR = OUT_DIR / "data"

# Core numeric features to analyse
CORE_NUMERIC_FEATURES = [
    "caller_af", "caller_gq", "caller_dp",
    "NumReads", "NumCpGs",
    "CramersV", "HeuristicScore",
    "PairwiseMeanDist", "PairwiseMedianDist",
    "HPMergedDelta", "AlleleDelta",
    "GlobalP", "Quality_Score",
    "HP_Ratio", "Coverage_Multiple",
    "caller_ad_ref", "caller_ad_alt",
    "HPFineF", "LabelHPPermanovaF", "LabelAllelePermanovaF",
    "UnassignedAffinity", "Stability",
    "hp_assign_rate", "hp0_ratio", "hp3_ratio",
    "caller_naf", "caller_ndp",
]

# Features used for PCA
PCA_FEATURES = [
    "caller_af", "caller_gq", "caller_dp",
    "NumReads", "NumCpGs",
    "CramersV", "HeuristicScore",
    "PairwiseMeanDist", "HPMergedDelta", "AlleleDelta",
    "Quality_Score", "HP_Ratio", "Coverage_Multiple",
]

# Features used for per-sample AUC heatmap
AUC_FEATURES = [
    "caller_gq", "caller_af", "caller_dp",
    "CramersV", "HeuristicScore", "Quality_Score",
    "PairwiseMeanDist", "HPMergedDelta", "AlleleDelta",
    "HP_Ratio", "Coverage_Multiple",
    "NumReads", "NumCpGs",
    "GlobalP", "Stability",
]

# Cancer type mapping
CANCER_TYPE_MAP = {
    "HCC1395": "TNBC",
    "HCC1395_DORADO": "TNBC",
    "HCC1937": "TNBC",
    "H1437": "Lung",
    "H2009": "Lung",
    "COLO829": "Melanoma",
    "HCC1954": "Breast",
}

SAMPLE_PALETTE = dict(zip(SAMPLE_ORDER, sns.color_palette("Set2", len(SAMPLE_ORDER))))

MODE_LABELS = {"paired": "Paired", "to": "Tumor-Only (TO)"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def manual_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Fallback AUC via trapezoidal ROC when sklearn unavailable."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    mask = np.isfinite(y_score) & np.isfinite(y_true)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 10:
        return float("nan")
    order = np.argsort(-y_score)
    y_true_sorted = y_true[order]
    npos = y_true_sorted.sum()
    nneg = len(y_true_sorted) - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    tpr_prev, fpr_prev = 0.0, 0.0
    auc_val = 0.0
    tp_cum, fp_cum = 0.0, 0.0
    for i in range(len(y_true_sorted)):
        if y_true_sorted[i] == 1:
            tp_cum += 1
        else:
            fp_cum += 1
        tpr = tp_cum / npos
        fpr = fp_cum / nneg
        auc_val += (fpr - fpr_prev) * (tpr + tpr_prev) / 2
        tpr_prev, fpr_prev = tpr, fpr
    return float(auc_val)


def safe_auc(y_true, y_score):
    """Try sklearn AUC first, fall back to manual."""
    try:
        return compute_auc(y_true, y_score)
    except Exception:
        return manual_auc(y_true, y_score)


def compute_cv(series: pd.Series) -> float:
    """Coefficient of Variation (std / |mean|). Returns NaN if mean ~ 0."""
    s = series.dropna()
    if len(s) < 2:
        return float("nan")
    m = s.mean()
    if abs(m) < 1e-10:
        return float("nan")
    return float(s.std() / abs(m))


def pca_via_svd(X: np.ndarray, n_components: int = 2):
    """PCA via numpy SVD. Returns (scores, explained_variance_ratio, loadings)."""
    # Standardize
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[sd < 1e-10] = 1.0
    Z = (X - mu) / sd
    # Replace NaN with 0 (column mean after standardization)
    Z = np.nan_to_num(Z, nan=0.0)
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    scores = U[:, :n_components] * S[:n_components]
    explained = (S ** 2) / (S ** 2).sum()
    return scores, explained[:n_components], Vt[:n_components, :]


# ---------------------------------------------------------------------------
# Figure generators
# ---------------------------------------------------------------------------


def fig01_feature_cv_across_samples(df: pd.DataFrame) -> Path:
    """CV of each feature across 7 samples (using sample-level medians)."""
    print("[O08] Fig 01: Feature CV across samples ...")

    cv_records = []
    for feat in CORE_NUMERIC_FEATURES:
        if feat not in df.columns:
            continue
        # Compute sample-level median for each sample
        sample_medians = df.groupby("sample")[feat].median()
        cv_val = compute_cv(sample_medians)
        cv_records.append({"feature": feat, "cv": cv_val})

    cv_df = pd.DataFrame(cv_records).dropna(subset=["cv"]).sort_values("cv", ascending=False).head(20)

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = sns.color_palette("viridis", len(cv_df))
    bars = ax.barh(range(len(cv_df)), cv_df["cv"].values, color=colors)
    ax.set_yticks(range(len(cv_df)))
    ax.set_yticklabels(cv_df["feature"].values)
    ax.set_xlabel("Coefficient of Variation (CV)")
    ax.set_title("Fig 01: Feature CV Across Samples (Top 20)\n(based on sample-level medians)")
    ax.invert_yaxis()
    for i, v in enumerate(cv_df["cv"].values):
        ax.text(v + 0.01, i, f"{v:.3f}", va="center", fontsize=8)

    add_caption(fig, f"Source: master dataset ({len(df):,} rows), {len(CORE_NUMERIC_FEATURES)} features")
    fig.tight_layout()

    # Save data
    cv_df.to_csv(DATA_DIR / "fig01_feature_cv.tsv", sep="\t", index=False)
    return save_figure(fig, FIG_DIR / "fig01_feature_cv_across_samples.png")


def fig02_sample_specific_tp_fp_rate(df: pd.DataFrame) -> Path:
    """TP/FP rate by sample x mode grouped bar chart."""
    print("[O08] Fig 02: Sample-specific TP/FP rate ...")

    records = []
    for samp in SAMPLE_ORDER:
        for mode in MODE_ORDER:
            sub = df[(df["sample"] == samp) & (df["mode"] == mode)]
            n_tp = (sub["truth_label"] == "TP").sum()
            n_fp = (sub["truth_label"] == "FP").sum()
            total = n_tp + n_fp
            records.append({
                "sample": samp, "mode": MODE_LABELS[mode],
                "TP_count": n_tp, "FP_count": n_fp, "total": total,
                "TP_rate": safe_div(n_tp, total),
                "FP_rate": safe_div(n_fp, total),
            })
    rate_df = pd.DataFrame(records)
    rate_df.to_csv(DATA_DIR / "fig02_tp_fp_rates.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    x = np.arange(len(SAMPLE_ORDER))
    width = 0.35

    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub = rate_df[rate_df["mode"] == MODE_LABELS[mode]]
        ax.bar(x - width/2, sub["TP_rate"].values, width, label="TP rate", color=COLOR_TP, alpha=0.8)
        ax.bar(x + width/2, sub["FP_rate"].values, width, label="FP rate", color=COLOR_FP, alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLE_ORDER, rotation=45, ha="right")
        ax.set_ylabel("Rate")
        ax.set_title(f"{MODE_LABELS[mode]}")
        ax.legend()
        ax.set_ylim(0, 1.05)
        # Annotate FP rates
        for i, row in enumerate(sub.itertuples()):
            ax.text(i + width/2, row.FP_rate + 0.02, f"{row.FP_rate:.2%}",
                    ha="center", va="bottom", fontsize=7, color=COLOR_FP)

    fig.suptitle("Fig 02: TP/FP Rate by Sample x Mode", fontsize=14, y=1.02)
    add_caption(fig, f"Source: master dataset ({len(df):,} rows)")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig02_sample_specific_tp_fp_rate.png")


def fig03_sample_clustering_by_features(df: pd.DataFrame) -> Path:
    """PCA on sample-level feature medians (14 points: 7 samples x 2 modes)."""
    print("[O08] Fig 03: Sample clustering by features (PCA) ...")

    # Compute sample x mode median for PCA features
    available_pca = [f for f in PCA_FEATURES if f in df.columns]
    medians = df.groupby(["sample", "mode"])[available_pca].median().reset_index()

    X = medians[available_pca].values
    scores, explained, loadings = pca_via_svd(X, n_components=2)

    medians["PC1"] = scores[:, 0]
    medians["PC2"] = scores[:, 1]
    medians.to_csv(DATA_DIR / "fig03_pca_scores.tsv", sep="\t", index=False)

    # Save loadings
    loadings_df = pd.DataFrame(loadings.T, index=available_pca, columns=["PC1", "PC2"])
    loadings_df.to_csv(DATA_DIR / "fig03_pca_loadings.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(10, 8))
    mode_markers = {"paired": "o", "to": "^"}
    for _, row in medians.iterrows():
        samp = row["sample"]
        mode = row["mode"]
        ax.scatter(row["PC1"], row["PC2"],
                   color=SAMPLE_PALETTE.get(samp, "#999"),
                   marker=mode_markers[mode],
                   s=150, edgecolors="black", linewidths=0.5, zorder=5)
        ax.annotate(f"{samp}\n({mode})", (row["PC1"], row["PC2"]),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 8), textcoords="offset points")

    ax.set_xlabel(f"PC1 ({explained[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1]*100:.1f}% variance)")
    ax.set_title("Fig 03: Sample Clustering by Feature Medians (PCA)\n"
                 f"(o = Paired, ^ = TO, {len(available_pca)} features)")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")

    # Legend for samples
    handles = [mpatches.Patch(color=SAMPLE_PALETTE[s], label=s) for s in SAMPLE_ORDER]
    ax.legend(handles=handles, loc="best", fontsize=7, ncol=2)

    add_caption(fig, f"PCA on {len(available_pca)} feature medians, "
                     f"{len(medians)} sample-mode combinations")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig03_sample_clustering_by_features.png")


def fig04_cancer_type_effect(df: pd.DataFrame) -> Path:
    """Feature distributions grouped by cancer type."""
    print("[O08] Fig 04: Cancer type effect ...")

    df_ct = df.copy()
    df_ct["cancer_type"] = df_ct["sample"].map(CANCER_TYPE_MAP)

    # Select key features
    key_feats = ["caller_af", "caller_gq", "CramersV", "Quality_Score",
                 "HPMergedDelta", "AlleleDelta"]
    available = [f for f in key_feats if f in df_ct.columns]
    ct_order = ["TNBC", "Lung", "Melanoma", "Breast"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    kw_records = []
    for idx, feat in enumerate(available):
        ax = axes[idx]
        data_groups = [df_ct.loc[df_ct["cancer_type"] == ct, feat].dropna().values
                       for ct in ct_order]
        # Box plot
        bp = ax.boxplot(data_groups, labels=ct_order, patch_artist=True,
                        showfliers=False, widths=0.6)
        colors_ct = sns.color_palette("Set3", 4)
        for patch, color in zip(bp["boxes"], colors_ct):
            patch.set_facecolor(color)
        ax.set_title(feat, fontsize=11)
        ax.tick_params(axis="x", rotation=30)

        # Kruskal-Wallis test across cancer types
        try:
            H_stat, p_val = sp_stats.kruskal(*[g for g in data_groups if len(g) > 0])
            kw_records.append({"feature": feat, "H_stat": H_stat, "p_value": p_val,
                               "n_groups": len([g for g in data_groups if len(g) > 0])})
            ax.text(0.02, 0.98, f"KW p={format_p(p_val)}", transform=ax.transAxes,
                    fontsize=7, va="top", ha="left",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        except Exception:
            kw_records.append({"feature": feat, "H_stat": np.nan, "p_value": np.nan, "n_groups": 0})

    # Hide unused axes
    for idx in range(len(available), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Fig 04: Feature Distribution by Cancer Type\n"
                 "(TNBC | Lung | Melanoma | Breast)", fontsize=14, y=1.02)
    add_caption(fig, f"Source: master dataset, {len(df):,} rows. Kruskal-Wallis H test shown.")
    fig.tight_layout()

    pd.DataFrame(kw_records).to_csv(DATA_DIR / "fig04_kruskal_wallis.tsv", sep="\t", index=False)
    return save_figure(fig, FIG_DIR / "fig04_cancer_type_effect.png")


def fig05_platform_effect(df: pd.DataFrame) -> Path:
    """HCC1395 (5kHz) vs HCC1395_DORADO platform comparison."""
    print("[O08] Fig 05: Platform effect (5kHz vs DORADO) ...")

    df5k = df[df["sample"] == "HCC1395"].copy()
    df_do = df[df["sample"] == "HCC1395_DORADO"].copy()
    df5k["platform_label"] = "HCC1395 (5kHz)"
    df_do["platform_label"] = "HCC1395_DORADO"
    df_comb = pd.concat([df5k, df_do], ignore_index=True)

    comp_feats = ["caller_af", "caller_gq", "CramersV", "Quality_Score",
                  "HPMergedDelta", "AlleleDelta", "PairwiseMeanDist", "NumReads"]
    available = [f for f in comp_feats if f in df_comb.columns]

    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    axes = axes.flatten()

    mw_records = []
    for idx, feat in enumerate(available[:8]):
        ax = axes[idx]
        sns.violinplot(x="platform_label", y=feat, data=df_comb,
                       ax=ax, cut=0, inner="quartile",
                       palette={"HCC1395 (5kHz)": "#66c2a5", "HCC1395_DORADO": "#fc8d62"})
        ax.set_xlabel("")
        ax.set_title(feat, fontsize=10)
        ax.tick_params(axis="x", rotation=15)

        # Mann-Whitney test
        g1 = df5k[feat].dropna().values
        g2 = df_do[feat].dropna().values
        mw = mannwhitney_test(g1, g2)
        mw_records.append({"feature": feat, **mw})
        ax.text(0.02, 0.98, f"MW p={format_p(mw['p_value'])}\nr={mw['effect_size_r']:.3f}",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    for idx in range(len(available), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Fig 05: Platform Effect -- HCC1395 (5kHz) vs HCC1395_DORADO\n"
                 "(Same cell line, different basecalling)", fontsize=13, y=1.02)
    add_caption(fig, f"5kHz n={len(df5k):,}, DORADO n={len(df_do):,}")
    fig.tight_layout()

    pd.DataFrame(mw_records).to_csv(DATA_DIR / "fig05_platform_mannwhitney.tsv", sep="\t", index=False)
    return save_figure(fig, FIG_DIR / "fig05_platform_effect_5khz_vs_dorado.png")


def fig06_f1_baseline_by_sample(df: pd.DataFrame) -> Path:
    """Precision (and count-based stats) by sample x mode."""
    print("[O08] Fig 06: F1 baseline by sample ...")

    records = []
    for samp in SAMPLE_ORDER:
        for mode in MODE_ORDER:
            sub = df[(df["sample"] == samp) & (df["mode"] == mode)]
            tp_n = (sub["truth_label"] == "TP").sum()
            fp_n = (sub["truth_label"] == "FP").sum()
            total = tp_n + fp_n
            precision = safe_div(tp_n, total)
            # Recall cannot be computed without FN; report as NaN
            records.append({
                "sample": samp, "mode": mode,
                "TP": tp_n, "FP": fp_n, "total": total,
                "precision": precision,
                "recall": float("nan"),  # FN unknown
                "f1": float("nan"),      # Cannot compute without recall
            })
    metrics_df = pd.DataFrame(records)
    metrics_df.to_csv(DATA_DIR / "fig06_precision_metrics.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    x = np.arange(len(SAMPLE_ORDER))
    width = 0.35

    # Left panel: Precision
    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub = metrics_df[metrics_df["mode"] == mode]
        bars = ax.bar(x, sub["precision"].values, width=0.6,
                      color=[COLOR_PAIRED if mode == "paired" else COLOR_TO] * len(x),
                      alpha=0.8, edgecolor="black", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLE_ORDER, rotation=45, ha="right")
        ax.set_ylabel("Precision (TP / (TP+FP))")
        ax.set_title(f"{MODE_LABELS[mode]}")
        ax.set_ylim(0, 1.05)
        for i, bar in enumerate(bars):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f"{sub.iloc[i]['precision']:.3f}\n(n={sub.iloc[i]['total']:,})",
                    ha="center", va="bottom", fontsize=7)

    fig.suptitle("Fig 06: Precision by Sample x Mode\n"
                 "(Recall/F1 not computable: FN count unavailable)", fontsize=13, y=1.02)
    add_caption(fig, f"Source: master dataset ({len(df):,} rows)")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig06_f1_baseline_by_sample.png")


def fig07_feature_stability_heatmap(df: pd.DataFrame) -> Path:
    """Per-sample per-mode AUC heatmap for key features."""
    print("[O08] Fig 07: Feature stability heatmap ...")

    y_bin = encode_truth_binary(df)
    available_feats = [f for f in AUC_FEATURES if f in df.columns]

    auc_records = []
    for samp in SAMPLE_ORDER:
        for mode in MODE_ORDER:
            mask = (df["sample"] == samp) & (df["mode"] == mode)
            sub_y = y_bin[mask]
            for feat in available_feats:
                sub_x = df.loc[mask, feat]
                auc_val = safe_auc(sub_y.values, sub_x.values)
                auc_records.append({
                    "sample": samp, "mode": mode, "feature": feat, "auc": auc_val
                })
    auc_df = pd.DataFrame(auc_records)
    auc_df.to_csv(DATA_DIR / "fig07_per_sample_auc.tsv", sep="\t", index=False)

    # Create two heatmaps side by side
    fig, axes = plt.subplots(1, 2, figsize=(18, 9))

    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub = auc_df[auc_df["mode"] == mode]
        pivot = sub.pivot(index="feature", columns="sample", values="auc")
        pivot = pivot.reindex(columns=SAMPLE_ORDER, index=available_feats)

        # Diverging colormap centered at 0.5
        sns.heatmap(pivot, annot=True, fmt=".3f", cmap="RdYlGn",
                    center=0.5, vmin=0.3, vmax=0.9,
                    ax=ax, linewidths=0.5, linecolor="white",
                    cbar_kws={"label": "AUC"})
        ax.set_title(f"{MODE_LABELS[mode]}", fontsize=12)
        ax.set_xlabel("")
        ax.set_ylabel("")

    fig.suptitle("Fig 07: Feature AUC Stability Heatmap (per-sample per-mode)\n"
                 "Green=discriminative for TP, Red=inverted or poor", fontsize=13, y=1.02)
    add_caption(fig, f"AUC computed on {len(available_feats)} features x "
                     f"{len(SAMPLE_ORDER)} samples x {len(MODE_ORDER)} modes")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig07_feature_stability_heatmap.png")


def fig08_h2009_anomaly(df: pd.DataFrame) -> Path:
    """H2009 LOH anomaly deep look."""
    print("[O08] Fig 08: H2009 anomaly deep look ...")

    h2009 = df[df["sample"] == "H2009"].copy()
    others = df[df["sample"] != "H2009"].copy()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: LOH rate by mode for H2009 vs others
    ax = axes[0, 0]
    loh_rates = []
    for label, sub in [("H2009", h2009), ("Others", others)]:
        for mode in MODE_ORDER:
            m = sub[sub["mode"] == mode]
            n_loh = m["Potential_LOH"].sum() if "Potential_LOH" in m.columns else 0
            loh_rates.append({"group": label, "mode": mode,
                              "LOH_rate": safe_div(n_loh, len(m)),
                              "n_loh": int(n_loh), "n_total": len(m)})
    lr_df = pd.DataFrame(loh_rates)
    lr_df.to_csv(DATA_DIR / "fig08_h2009_loh_rates.tsv", sep="\t", index=False)

    x = np.arange(2)
    w = 0.35
    for i, grp in enumerate(["H2009", "Others"]):
        sub = lr_df[lr_df["group"] == grp]
        ax.bar(x + i * w, sub["LOH_rate"].values, w,
               label=grp, color=["#e74c3c", "#3498db"][i], alpha=0.8)
    ax.set_xticks(x + w/2)
    ax.set_xticklabels(["Paired", "TO"])
    ax.set_ylabel("LOH Rate")
    ax.set_title("A: LOH Rate (H2009 vs Others)")
    ax.legend()

    # Panel B: FP composition in LOH vs non-LOH for H2009
    ax = axes[0, 1]
    for mode in MODE_ORDER:
        sub = h2009[h2009["mode"] == mode]
        loh_mask = sub["Potential_LOH"] == True
        n_fp_loh = ((sub["truth_label"] == "FP") & loh_mask).sum()
        n_fp_noloh = ((sub["truth_label"] == "FP") & ~loh_mask).sum()
        n_tp_loh = ((sub["truth_label"] == "TP") & loh_mask).sum()
        n_tp_noloh = ((sub["truth_label"] == "TP") & ~loh_mask).sum()
        categories = ["FP+LOH", "FP+noLOH", "TP+LOH", "TP+noLOH"]
        vals = [n_fp_loh, n_fp_noloh, n_tp_loh, n_tp_noloh]
        offset = 0 if mode == "paired" else 0.4
        bx = np.arange(4) + offset
        ax.bar(bx, vals, 0.35, label=MODE_LABELS[mode],
               color=COLOR_PAIRED if mode == "paired" else COLOR_TO, alpha=0.7)
        for j, v in enumerate(vals):
            ax.text(bx[j], v + 50, str(v), ha="center", fontsize=7)

    ax.set_xticks(np.arange(4) + 0.2)
    ax.set_xticklabels(["FP+LOH", "FP+noLOH", "TP+LOH", "TP+noLOH"], fontsize=8)
    ax.set_title("B: H2009 TP/FP x LOH Composition")
    ax.legend()

    # Panel C: H2009 Quality_Score distribution vs others
    ax = axes[1, 0]
    for label, sub, color in [("H2009", h2009, "#e74c3c"), ("Others", others, "#3498db")]:
        vals = sub["Quality_Score"].dropna()
        ax.hist(vals, bins=50, alpha=0.5, label=f"{label} (n={len(vals):,})",
                color=color, density=True)
    ax.set_xlabel("Quality_Score")
    ax.set_ylabel("Density")
    ax.set_title("C: Quality Score Distribution")
    ax.legend()

    # Panel D: H2009 caller_af distribution
    ax = axes[1, 1]
    for label, sub, color in [("H2009", h2009, "#e74c3c"), ("Others", others, "#3498db")]:
        vals = sub["caller_af"].dropna()
        ax.hist(vals, bins=50, alpha=0.5, label=f"{label} (n={len(vals):,})",
                color=color, density=True)
    ax.set_xlabel("Caller AF")
    ax.set_ylabel("Density")
    ax.set_title("D: Caller AF Distribution")
    ax.legend()

    fig.suptitle("Fig 08: H2009 Anomaly Deep Look\n"
                 "(H2009 = 36.2% of dataset, 68% FP in LOH area)", fontsize=14, y=1.02)
    add_caption(fig, f"H2009 n={len(h2009):,}, Others n={len(others):,}")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig08_h2009_anomaly_deep_look.png")


def fig09_hcc1395_chr8_hotspot(df: pd.DataFrame) -> Path:
    """HCC1395 chr8 LOH+ASM hotspot analysis."""
    print("[O08] Fig 09: HCC1395 chr8 hotspot ...")

    hcc = df[df["sample"] == "HCC1395"].copy()
    chr_col = "Chr"

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: LOH rate per chromosome for HCC1395
    ax = axes[0, 0]
    chr_loh = hcc.groupby(chr_col).agg(
        n_total=("truth_label", "size"),
        n_loh=("Potential_LOH", "sum"),
    ).reset_index()
    chr_loh["loh_rate"] = chr_loh["n_loh"] / chr_loh["n_total"]
    # Order chromosomes naturally
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_loh_sorted = chr_loh.set_index(chr_col).reindex(
        [c for c in chrom_order if c in chr_loh[chr_col].values]).reset_index()

    colors_chr = ["#e74c3c" if c == "chr8" else "#3498db" for c in chr_loh_sorted[chr_col]]
    ax.bar(range(len(chr_loh_sorted)), chr_loh_sorted["loh_rate"], color=colors_chr)
    ax.set_xticks(range(len(chr_loh_sorted)))
    ax.set_xticklabels(chr_loh_sorted[chr_col].str.replace("chr", ""), rotation=45, fontsize=7)
    ax.set_ylabel("LOH Rate")
    ax.set_title("A: HCC1395 LOH Rate per Chromosome (chr8 = red)")

    # Panel B: chr8 vs non-chr8 TP/FP enrichment
    ax = axes[0, 1]
    is_chr8 = hcc[chr_col] == "chr8"
    for mode in MODE_ORDER:
        m = hcc[hcc["mode"] == mode]
        chr8_m = m[m[chr_col] == "chr8"]
        other_m = m[m[chr_col] != "chr8"]

        tp_chr8 = (chr8_m["truth_label"] == "TP").sum()
        fp_chr8 = (chr8_m["truth_label"] == "FP").sum()
        tp_other = (other_m["truth_label"] == "TP").sum()
        fp_other = (other_m["truth_label"] == "FP").sum()

        fp_rate_chr8 = safe_div(fp_chr8, fp_chr8 + tp_chr8)
        fp_rate_other = safe_div(fp_other, fp_other + tp_other)

        offset = 0 if mode == "paired" else 0.35
        xp = np.array([0, 1]) + offset
        ax.bar(xp, [fp_rate_chr8, fp_rate_other], 0.3,
               label=f"{MODE_LABELS[mode]}",
               color=COLOR_PAIRED if mode == "paired" else COLOR_TO, alpha=0.7)
        ax.text(xp[0], fp_rate_chr8 + 0.01, f"{fp_rate_chr8:.3f}", ha="center", fontsize=8)
        ax.text(xp[1], fp_rate_other + 0.01, f"{fp_rate_other:.3f}", ha="center", fontsize=8)

    ax.set_xticks([0.175, 1.175])
    ax.set_xticklabels(["chr8", "Other chr"])
    ax.set_ylabel("FP Rate")
    ax.set_title("B: HCC1395 FP Rate -- chr8 vs Others")
    ax.legend()

    # Panel C: HPMergedSig rate per chromosome
    ax = axes[1, 0]
    chr_hpsig = hcc.groupby(chr_col)["HPMergedSig"].mean().reset_index()
    chr_hpsig_sorted = chr_hpsig.set_index(chr_col).reindex(
        [c for c in chrom_order if c in chr_hpsig[chr_col].values]).reset_index()
    colors_chr2 = ["#e74c3c" if c == "chr8" else "#2ecc71" for c in chr_hpsig_sorted[chr_col]]
    ax.bar(range(len(chr_hpsig_sorted)), chr_hpsig_sorted["HPMergedSig"], color=colors_chr2)
    ax.set_xticks(range(len(chr_hpsig_sorted)))
    ax.set_xticklabels(chr_hpsig_sorted[chr_col].str.replace("chr", ""), rotation=45, fontsize=7)
    ax.set_ylabel("HPMergedSig Rate")
    ax.set_title("C: HCC1395 HPMergedSig Rate per Chromosome")

    # Panel D: chr8 CramersV vs other chromosomes
    ax = axes[1, 1]
    chr8_cv = hcc.loc[is_chr8, "CramersV"].dropna()
    other_cv = hcc.loc[~is_chr8, "CramersV"].dropna()
    ax.hist(chr8_cv, bins=30, alpha=0.6, label=f"chr8 (n={len(chr8_cv):,})",
            color="#e74c3c", density=True)
    ax.hist(other_cv, bins=30, alpha=0.6, label=f"Other (n={len(other_cv):,})",
            color="#3498db", density=True)
    mw = mannwhitney_test(chr8_cv.values, other_cv.values)
    ax.text(0.98, 0.98, f"MW p={format_p(mw['p_value'])}\nr={mw['effect_size_r']:.3f}",
            transform=ax.transAxes, fontsize=8, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    ax.set_xlabel("CramersV")
    ax.set_ylabel("Density")
    ax.set_title("D: CramersV Distribution -- chr8 vs Others")
    ax.legend()

    fig.suptitle("Fig 09: HCC1395 chr8 LOH+ASM Hotspot\n"
                 "(Known FP enrichment hotspot)", fontsize=14, y=1.02)
    add_caption(fig, f"HCC1395 total n={len(hcc):,}")
    fig.tight_layout()

    # Save chr8 stats
    chr_loh.to_csv(DATA_DIR / "fig09_chr_loh_rates.tsv", sep="\t", index=False)
    return save_figure(fig, FIG_DIR / "fig09_hcc1395_chr8_hotspot.png")


def fig10_sample_purity_effect(df: pd.DataFrame) -> Path:
    """Caller AF median as purity proxy vs precision."""
    print("[O08] Fig 10: Sample purity effect ...")

    records = []
    for samp in SAMPLE_ORDER:
        for mode in MODE_ORDER:
            sub = df[(df["sample"] == samp) & (df["mode"] == mode)]
            af_med = sub["caller_af"].median()
            tp_n = (sub["truth_label"] == "TP").sum()
            fp_n = (sub["truth_label"] == "FP").sum()
            precision = safe_div(tp_n, tp_n + fp_n)
            records.append({
                "sample": samp, "mode": mode,
                "af_median": af_med, "precision": precision,
                "n": len(sub), "TP": tp_n, "FP": fp_n,
            })
    purity_df = pd.DataFrame(records)
    purity_df.to_csv(DATA_DIR / "fig10_purity_proxy.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub = purity_df[purity_df["mode"] == mode]
        for _, row in sub.iterrows():
            ax.scatter(row["af_median"], row["precision"],
                       color=SAMPLE_PALETTE.get(row["sample"], "#999"),
                       s=120, edgecolors="black", linewidths=0.5, zorder=5)
            ax.annotate(row["sample"], (row["af_median"], row["precision"]),
                        fontsize=7, xytext=(5, 5), textcoords="offset points")
        ax.set_xlabel("Median Caller AF (purity proxy)")
        ax.set_ylabel("Precision")
        ax.set_title(f"{MODE_LABELS[mode]}")

        # Correlation
        if len(sub) >= 3:
            r, p = sp_stats.pearsonr(sub["af_median"], sub["precision"])
            ax.text(0.02, 0.02, f"Pearson r={r:.3f}, p={format_p(p)}",
                    transform=ax.transAxes, fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    fig.suptitle("Fig 10: Sample Purity Effect\n"
                 "(Median caller_af as proxy for tumor purity)", fontsize=13, y=1.02)
    handles = [mpatches.Patch(color=SAMPLE_PALETTE[s], label=s) for s in SAMPLE_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=7,
               bbox_to_anchor=(0.5, -0.02))
    add_caption(fig, f"Source: master dataset ({len(df):,} rows)")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig10_sample_purity_effect.png")


def fig11_colo829_low_qs(df: pd.DataFrame) -> Path:
    """COLO829 QS anomaly investigation."""
    print("[O08] Fig 11: COLO829 low QS investigation ...")

    colo = df[df["sample"] == "COLO829"].copy()
    others = df[df["sample"] != "COLO829"].copy()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: QS distribution COLO829 vs others
    ax = axes[0, 0]
    for label, sub, color in [("COLO829", colo, "#e74c3c"), ("Others", others, "#3498db")]:
        vals = sub["Quality_Score"].dropna()
        ax.hist(vals, bins=50, alpha=0.5, label=f"{label} (n={len(vals):,})",
                color=color, density=True)
    ax.set_xlabel("Quality_Score")
    ax.set_ylabel("Density")
    ax.set_title("A: Quality Score Distribution")
    ax.legend()
    mw = mannwhitney_test(colo["Quality_Score"].dropna().values,
                          others["Quality_Score"].dropna().values)
    ax.text(0.98, 0.98, f"MW p={format_p(mw['p_value'])}\nr={mw['effect_size_r']:.3f}",
            transform=ax.transAxes, fontsize=8, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Panel B: QS breakdown by tier (COLO829 vs others)
    ax = axes[0, 1]
    qt_col = "Quality_Tier"
    tier_order = ["Low", "Medium", "High"]
    for i, (label, sub, color) in enumerate([("COLO829", colo, "#e74c3c"),
                                              ("Others", others, "#3498db")]):
        tier_counts = sub[qt_col].value_counts(normalize=True)
        vals = [tier_counts.get(t, 0) for t in tier_order]
        offset = i * 0.35
        ax.bar(np.arange(len(tier_order)) + offset, vals, 0.3,
               label=label, color=color, alpha=0.7)
        for j, v in enumerate(vals):
            ax.text(j + offset, v + 0.01, f"{v:.2%}", ha="center", fontsize=7)
    ax.set_xticks(np.arange(len(tier_order)) + 0.175)
    ax.set_xticklabels(tier_order)
    ax.set_ylabel("Proportion")
    ax.set_title("B: Quality Tier Distribution")
    ax.legend()

    # Panel C: QS component -- LOH rate and Coverage for COLO829
    ax = axes[1, 0]
    features_comp = ["Potential_LOH", "Coverage_Multiple"]
    for feat in features_comp:
        if feat == "Potential_LOH":
            colo_val = colo[feat].mean()
            other_val = others[feat].mean()
            ax.bar([f"COLO829\n{feat}", f"Others\n{feat}"],
                   [colo_val, other_val],
                   color=["#e74c3c", "#3498db"], alpha=0.7)
            ax.text(0, colo_val + 0.01, f"{colo_val:.3f}", ha="center", fontsize=8)
            ax.text(1, other_val + 0.01, f"{other_val:.3f}", ha="center", fontsize=8)
    ax.set_ylabel("Rate / Mean")
    ax.set_title("C: LOH Rate -- COLO829 vs Others")

    # Panel D: COLO829 caller_gq distribution
    ax = axes[1, 1]
    for label, sub, color in [("COLO829", colo, "#e74c3c"), ("Others", others, "#3498db")]:
        vals = sub["caller_gq"].dropna()
        ax.hist(vals, bins=50, alpha=0.5, label=f"{label} (n={len(vals):,})",
                color=color, density=True)
    ax.set_xlabel("Caller GQ")
    ax.set_ylabel("Density")
    ax.set_title("D: Caller GQ Distribution")
    ax.legend()
    mw2 = mannwhitney_test(colo["caller_gq"].dropna().values,
                           others["caller_gq"].dropna().values)
    ax.text(0.98, 0.98, f"MW p={format_p(mw2['p_value'])}\nr={mw2['effect_size_r']:.3f}",
            transform=ax.transAxes, fontsize=8, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    fig.suptitle("Fig 11: COLO829 Low QS Investigation\n"
                 "(Why COLO829 has anomalously low Quality Scores)", fontsize=14, y=1.02)
    add_caption(fig, f"COLO829 n={len(colo):,}, Others n={len(others):,}")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig11_colo829_low_qs_investigation.png")


def fig12_cross_sample_feature_importance(df: pd.DataFrame) -> Path:
    """Top-5 AUC features per sample as facet heatmap."""
    print("[O08] Fig 12: Cross-sample feature importance ...")

    y_bin = encode_truth_binary(df)
    available_feats = [f for f in AUC_FEATURES if f in df.columns]

    # Compute AUC for each sample x feature (both modes combined)
    auc_records = []
    for samp in SAMPLE_ORDER:
        mask = df["sample"] == samp
        sub_y = y_bin[mask]
        for feat in available_feats:
            sub_x = df.loc[mask, feat]
            auc_val = safe_auc(sub_y.values, sub_x.values)
            auc_records.append({"sample": samp, "feature": feat, "auc": auc_val})
    auc_df = pd.DataFrame(auc_records)
    auc_df.to_csv(DATA_DIR / "fig12_cross_sample_auc.tsv", sep="\t", index=False)

    # Find top-5 features per sample
    top5_records = []
    for samp in SAMPLE_ORDER:
        sub = auc_df[auc_df["sample"] == samp].dropna(subset=["auc"])
        # Use max(auc, 1-auc) to capture inverted features
        sub = sub.copy()
        sub["auc_abs"] = sub["auc"].apply(lambda x: max(x, 1 - x) if not np.isnan(x) else 0)
        top5 = sub.nlargest(5, "auc_abs")
        for rank, (_, row) in enumerate(top5.iterrows(), 1):
            top5_records.append({
                "sample": samp, "rank": rank,
                "feature": row["feature"], "auc": row["auc"],
                "auc_abs": row["auc_abs"],
            })
    top5_df = pd.DataFrame(top5_records)
    top5_df.to_csv(DATA_DIR / "fig12_top5_per_sample.tsv", sep="\t", index=False)

    # Build heatmap: union of all top-5 features
    all_top5_feats = top5_df["feature"].unique().tolist()
    pivot = auc_df[auc_df["feature"].isin(all_top5_feats)].pivot(
        index="feature", columns="sample", values="auc")
    pivot = pivot.reindex(columns=SAMPLE_ORDER)

    fig, ax = plt.subplots(figsize=(12, max(6, len(all_top5_feats) * 0.5)))
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="RdYlGn",
                center=0.5, vmin=0.3, vmax=0.9,
                ax=ax, linewidths=0.5, linecolor="white",
                cbar_kws={"label": "AUC"})
    ax.set_title("Fig 12: Cross-Sample Feature Importance\n"
                 "(Union of top-5 AUC features per sample, both modes)", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Mark top-5 cells with asterisks
    for _, row in top5_df.iterrows():
        feat_idx = list(pivot.index).index(row["feature"]) if row["feature"] in pivot.index else -1
        samp_idx = SAMPLE_ORDER.index(row["sample"]) if row["sample"] in SAMPLE_ORDER else -1
        if feat_idx >= 0 and samp_idx >= 0:
            ax.text(samp_idx + 0.5, feat_idx + 0.85, "*",
                    ha="center", va="center", fontsize=10, color="black", fontweight="bold")

    add_caption(fig, f"* = top-5 feature for that sample. "
                     f"AUC computed on TP(1) vs FP(0), both modes combined.")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "fig12_cross_sample_feature_importance.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    start_time = datetime.now()
    setup_plot_style()

    ensure_dir(OUT_DIR)
    ensure_dir(FIG_DIR)
    ensure_dir(DATA_DIR)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print("=" * 70)
    print(f"[O08] {OBSERVATION_TITLE}")
    print(f"[O08] Output: {OUT_DIR}")
    print("=" * 70)

    df = load_master_dataset()
    print(f"[O08] Loaded {len(df):,} rows x {len(df.columns)} cols")

    # Validate expected columns
    required = ["truth_label", "sample", "mode", "caller_af", "Quality_Score"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns: {missing}")

    # ------------------------------------------------------------------
    # Write round context
    # ------------------------------------------------------------------
    write_round_context(
        out_dir=OUT_DIR,
        observation_id=OBSERVATION_ID,
        title=OBSERVATION_TITLE,
        description="Identifies feature variability across 7 samples x 2 modes, "
                    "cancer type effects, platform effects, and sample-specific anomalies.",
        script_path=SCRIPT_PATH,
        data_sources=[{"name": "master_dataset", "path": str(load_master_dataset.__defaults__)}],
        row_count=len(df),
        col_count=len(df.columns),
    )

    # ------------------------------------------------------------------
    # Write feature statistics
    # ------------------------------------------------------------------
    write_feature_statistics(
        df, CORE_NUMERIC_FEATURES,
        DATA_DIR / "feature_descriptive_stats.tsv"
    )

    # ------------------------------------------------------------------
    # Generate all 12 figures
    # ------------------------------------------------------------------
    figures = []

    figures.append(fig01_feature_cv_across_samples(df))
    figures.append(fig02_sample_specific_tp_fp_rate(df))
    figures.append(fig03_sample_clustering_by_features(df))
    figures.append(fig04_cancer_type_effect(df))
    figures.append(fig05_platform_effect(df))
    figures.append(fig06_f1_baseline_by_sample(df))
    figures.append(fig07_feature_stability_heatmap(df))
    figures.append(fig08_h2009_anomaly(df))
    figures.append(fig09_hcc1395_chr8_hotspot(df))
    figures.append(fig10_sample_purity_effect(df))
    figures.append(fig11_colo829_low_qs(df))
    figures.append(fig12_cross_sample_feature_importance(df))

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = (datetime.now() - start_time).total_seconds()
    print("\n" + "=" * 70)
    print(f"[O08] COMPLETE in {elapsed:.1f}s")
    print(f"[O08] Figures generated: {len(figures)}")
    for fp in figures:
        print(f"  - {fp.name}")
    print(f"[O08] Output: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
