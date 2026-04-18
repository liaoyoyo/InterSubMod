#!/usr/bin/env python3
"""
O05: Methylation Feature TP/FP Discriminability Observation
============================================================

Systematic observation of ISM-derived methylation features and their
independent discriminative power for TP vs FP, separated by paired vs TO modes.

Quantifies: distributions, AUC, effect sizes, correlations with caller features,
LOH-stratified behaviour, and PassedGating utility.

Part of the O1-O10 systematic observation framework for Phase 1A ML baseline.

Generated: 2026-03-31
"""

from __future__ import annotations

import json
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

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
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OBSERVATION_ID = "O05"
OBSERVATION_TITLE = "Methylation Feature TP/FP Discriminability Observation"
SCRIPT_PATH = "scripts/analysis/build_observation_O05_methylation_features.py"

OUT_DIR = OUTPUT_ROOT / "20260401_O05_methylation_features"
FIG_DIR = OUT_DIR
DATA_DIR = OUT_DIR / "data"

# ISM methylation features to observe
METHYL_FEATURES = [
    "PairwiseMeanDist", "PairwiseMedianDist",
    "AlleleDelta", "CramersV", "GlobalP",
    "PassedGating", "HeuristicScore",
    "HPMergedDelta", "HPMergedSig",
]

# Continuous methylation features (for violin, AUC, effect size)
METHYL_CONTINUOUS = [
    "PairwiseMeanDist", "PairwiseMedianDist",
    "AlleleDelta", "CramersV", "GlobalP",
    "HeuristicScore", "HPMergedDelta",
]

# Binary methylation features
METHYL_BINARY = ["PassedGating", "HPMergedSig"]

# Caller features for cross-correlation
CALLER_FEATURES = ["caller_af", "caller_gq", "caller_dp", "Quality_Score"]

# All columns we need from the master dataset
USECOLS = [
    "truth_label", "sample", "mode", "variant_key", "af_bin",
    "Potential_LOH",
    # Methylation features
    "PairwiseMeanDist", "PairwiseMedianDist",
    "AlleleDelta", "CramersV", "GlobalP",
    "PassedGating", "HeuristicScore",
    "HPMergedDelta", "HPMergedSig",
    # Caller features (for correlation)
    "caller_af", "caller_gq", "caller_dp", "Quality_Score",
]

MODE_LABELS = {"paired": "Paired", "to": "Tumor-Only (TO)"}
AF_BIN_ORDER = [
    "0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
    "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0",
]
LOH_LABELS = {True: "LOH", False: "non-LOH"}

CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz "
    "(20260327_loh_round1_cross_sample_audit) | "
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _auc_safe(y_true, y_score):
    """AUC with fallback."""
    try:
        return compute_auc(y_true, y_score)
    except Exception:
        return float("nan")


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
    print(f"={'=' * 70}")
    print(f"  {OBSERVATION_ID}: {OBSERVATION_TITLE}")
    print(f"={'=' * 70}")

    df_raw = load_master_dataset(usecols=USECOLS)
    total_rows = len(df_raw)

    # Filter to TP/FP only
    df = df_raw[df_raw["truth_label"].isin(["TP", "FP"])].copy()
    print(f"[O05] Filtered to TP/FP: {len(df):,} rows (dropped {total_rows - len(df):,})")

    # Binary truth label
    df["truth_bin"] = encode_truth_binary(df)

    # Ensure boolean columns
    for bcol in ["PassedGating", "HPMergedSig", "Potential_LOH"]:
        if bcol in df.columns:
            df[bcol] = df[bcol].astype(str).str.lower().isin(["true", "1", "yes"])

    # Mode x Truth counts
    mode_truth_counts = df.groupby(["mode", "truth_label"]).size().unstack(fill_value=0)
    print(f"\n[O05] Mode x Truth distribution:")
    print(mode_truth_counts.to_string())

    # ------------------------------------------------------------------
    # Data Output: source_summary.tsv
    # ------------------------------------------------------------------
    summary_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for tl in ["TP", "FP"]:
            dt = dm[dm["truth_label"] == tl]
            row = {
                "mode": mode,
                "truth_label": tl,
                "n_regions": len(dt),
                "n_samples": dt["sample"].nunique(),
            }
            for feat in METHYL_CONTINUOUS:
                row[f"{feat}_median"] = dt[feat].median() if feat in dt.columns else np.nan
            for feat in METHYL_BINARY:
                if feat in dt.columns:
                    row[f"{feat}_rate"] = dt[feat].mean()
            summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(DATA_DIR / "source_summary.tsv", sep="\t", index=False)
    print(f"[O05] Wrote source_summary.tsv ({len(summary_df)} rows)")

    # ------------------------------------------------------------------
    # Data Output: feature_statistics.tsv
    # ------------------------------------------------------------------
    write_feature_statistics(df, METHYL_FEATURES, DATA_DIR / "feature_statistics.tsv")

    # ------------------------------------------------------------------
    # Statistical Tests: Mann-Whitney U + Cohen's d per feature per mode
    # ------------------------------------------------------------------
    print("\n[O05] Computing statistical tests ...")
    effect_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        tp_mask = dm["truth_label"] == "TP"
        fp_mask = dm["truth_label"] == "FP"
        for feat in METHYL_CONTINUOUS:
            if feat not in dm.columns:
                continue
            tp_vals = dm.loc[tp_mask, feat].dropna().values
            fp_vals = dm.loc[fp_mask, feat].dropna().values

            mwu = mannwhitney_test(tp_vals, fp_vals)
            es = compute_effect_size(tp_vals, fp_vals)

            effect_rows.append({
                "mode": mode,
                "feature": feat,
                "tp_n": len(tp_vals),
                "fp_n": len(fp_vals),
                "tp_median": float(np.median(tp_vals)) if len(tp_vals) > 0 else np.nan,
                "fp_median": float(np.median(fp_vals)) if len(fp_vals) > 0 else np.nan,
                "tp_mean": float(np.mean(tp_vals)) if len(tp_vals) > 0 else np.nan,
                "fp_mean": float(np.mean(fp_vals)) if len(fp_vals) > 0 else np.nan,
                "mwu_statistic": mwu["statistic"],
                "mwu_p_value": mwu["p_value"],
                "mwu_effect_r": mwu["effect_size_r"],
                "cohens_d": es["d"],
                "cohens_d_ci_lo": es["ci_lo"],
                "cohens_d_ci_hi": es["ci_hi"],
            })
    effect_df = pd.DataFrame(effect_rows)
    effect_df.to_csv(DATA_DIR / "methyl_effect_size.tsv", sep="\t", index=False)
    print(f"[O05] Wrote methyl_effect_size.tsv ({len(effect_df)} rows)")

    # ------------------------------------------------------------------
    # AUC by mode (continuous features)
    # ------------------------------------------------------------------
    print("\n[O05] Computing AUC by mode ...")
    auc_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y_true = dm["truth_bin"].values
        for feat in METHYL_CONTINUOUS:
            if feat not in dm.columns:
                continue
            y_score = dm[feat].values
            auc_val = _auc_safe(y_true, y_score)
            auc_rows.append({
                "mode": mode,
                "feature": feat,
                "auc": round(auc_val, 4) if not np.isnan(auc_val) else np.nan,
                "n_total": len(dm),
                "n_valid": int(np.sum(np.isfinite(y_score) & np.isfinite(y_true))),
            })
        # For binary features, compute AUC too
        for feat in METHYL_BINARY:
            if feat not in dm.columns:
                continue
            y_score = dm[feat].astype(float).values
            auc_val = _auc_safe(y_true, y_score)
            auc_rows.append({
                "mode": mode,
                "feature": feat,
                "auc": round(auc_val, 4) if not np.isnan(auc_val) else np.nan,
                "n_total": len(dm),
                "n_valid": int(np.sum(np.isfinite(y_score) & np.isfinite(y_true))),
            })
    auc_df = pd.DataFrame(auc_rows)
    auc_df.to_csv(DATA_DIR / "methyl_auc_by_mode.tsv", sep="\t", index=False)
    print(f"[O05] Wrote methyl_auc_by_mode.tsv ({len(auc_df)} rows)")

    # ------------------------------------------------------------------
    # Chi2 tests for binary features
    # ------------------------------------------------------------------
    print("\n[O05] Computing chi2 tests for binary features ...")
    chi2_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for feat in METHYL_BINARY:
            if feat not in dm.columns:
                continue
            ct = pd.crosstab(dm["truth_label"], dm[feat])
            if ct.shape == (2, 2):
                result = chi2_test(ct.values)
                chi2_rows.append({
                    "mode": mode,
                    "feature": feat,
                    "chi2": result["chi2"],
                    "p_value": result["p_value"],
                    "cramers_v": result["cramers_v"],
                    "n_total": len(dm),
                    "true_tp_rate": dm.loc[dm["truth_label"] == "TP", feat].mean(),
                    "true_fp_rate": dm.loc[dm["truth_label"] == "FP", feat].mean(),
                })
    chi2_df = pd.DataFrame(chi2_rows)
    chi2_df.to_csv(DATA_DIR / "methyl_binary_chi2.tsv", sep="\t", index=False)
    print(f"[O05] Wrote methyl_binary_chi2.tsv ({len(chi2_df)} rows)")

    # ------------------------------------------------------------------
    # PassedGating enrichment by sample x truth
    # ------------------------------------------------------------------
    print("\n[O05] Computing PassedGating rates by sample x truth ...")
    gating_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for samp in SAMPLE_ORDER:
            ds = dm[dm["sample"] == samp]
            if len(ds) == 0:
                continue
            for tl in ["TP", "FP"]:
                dt = ds[ds["truth_label"] == tl]
                if len(dt) == 0:
                    continue
                gating_rows.append({
                    "mode": mode,
                    "sample": samp,
                    "truth_label": tl,
                    "n": len(dt),
                    "passed_gating_rate": dt["PassedGating"].mean(),
                    "passed_gating_n": int(dt["PassedGating"].sum()),
                })
    gating_df = pd.DataFrame(gating_rows)
    gating_df.to_csv(DATA_DIR / "passed_gating_by_sample_truth.tsv", sep="\t", index=False)
    print(f"[O05] Wrote passed_gating_by_sample_truth.tsv ({len(gating_df)} rows)")

    # ------------------------------------------------------------------
    # LOH-stratified statistics
    # ------------------------------------------------------------------
    print("\n[O05] Computing LOH-stratified effect sizes ...")
    loh_effect_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for loh_val in [True, False]:
            dl = dm[dm["Potential_LOH"] == loh_val]
            tp_mask = dl["truth_label"] == "TP"
            fp_mask = dl["truth_label"] == "FP"
            n_tp = tp_mask.sum()
            n_fp = fp_mask.sum()
            loh_label = LOH_LABELS[loh_val]
            for feat in METHYL_CONTINUOUS:
                if feat not in dl.columns:
                    continue
                tp_vals = dl.loc[tp_mask, feat].dropna().values
                fp_vals = dl.loc[fp_mask, feat].dropna().values
                es = compute_effect_size(tp_vals, fp_vals)
                auc_val = _auc_safe(
                    dl.loc[tp_mask | fp_mask, "truth_bin"].values,
                    dl.loc[tp_mask | fp_mask, feat].values,
                )
                loh_effect_rows.append({
                    "mode": mode,
                    "loh_status": loh_label,
                    "feature": feat,
                    "n_tp": n_tp,
                    "n_fp": n_fp,
                    "auc": round(auc_val, 4) if not np.isnan(auc_val) else np.nan,
                    "cohens_d": es["d"],
                    "cohens_d_ci_lo": es["ci_lo"],
                    "cohens_d_ci_hi": es["ci_hi"],
                })
    loh_effect_df = pd.DataFrame(loh_effect_rows)
    loh_effect_df.to_csv(DATA_DIR / "methyl_loh_stratified_effects.tsv", sep="\t", index=False)
    print(f"[O05] Wrote methyl_loh_stratified_effects.tsv ({len(loh_effect_df)} rows)")

    # ------------------------------------------------------------------
    # Correlation matrix (methyl x caller, Spearman)
    # ------------------------------------------------------------------
    print("\n[O05] Computing Spearman correlation matrices ...")
    corr_features = [f for f in METHYL_CONTINUOUS + CALLER_FEATURES if f in df.columns]
    corr_data = {}
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        sub = dm[corr_features].copy()
        # Replace inf with nan
        sub = sub.replace([np.inf, -np.inf], np.nan)
        corr_mat = sub.corr(method="spearman")
        corr_mat.to_csv(DATA_DIR / f"spearman_corr_{mode}.tsv", sep="\t")
        corr_data[mode] = corr_mat
    print(f"[O05] Wrote spearman_corr_paired.tsv and spearman_corr_to.tsv")

    # ==================================================================
    # FIGURES
    # ==================================================================

    caption_base = (
        f"O05 Methylation Features | N={len(df):,} (TP+FP) | "
        f"Source: all_region_rows.tsv.gz | {datetime.now().strftime('%Y-%m-%d')}"
    )

    # ------------------------------------------------------------------
    # Fig01: Methylation feature distributions by truth (8 violins)
    # ------------------------------------------------------------------
    print("\n[O05] Generating fig01 ...")
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes_flat = axes.flatten()
    plot_features = METHYL_CONTINUOUS + ["PassedGating"]  # 8 panels
    for i, feat in enumerate(plot_features):
        ax = axes_flat[i]
        if feat not in df.columns:
            ax.set_title(f"{feat}\n(not found)")
            ax.axis("off")
            continue

        if feat == "PassedGating":
            # Bar plot for binary feature
            rates = df.groupby("truth_label")[feat].mean()
            bars = ax.bar(["TP", "FP"], [rates.get("TP", 0), rates.get("FP", 0)],
                          color=[COLOR_TP, COLOR_FP], alpha=0.8, edgecolor="white")
            ax.set_ylabel("PassedGating Rate")
            ax.set_title("PassedGating", fontsize=11, fontweight="bold")
            for bar, val in zip(bars, [rates.get("TP", 0), rates.get("FP", 0)]):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                        f"{val:.3f}", ha="center", fontsize=8)
        else:
            # Violin plot
            plot_df = df[["truth_label", feat]].dropna()
            if len(plot_df) > 200_000:
                # Subsample for violin performance
                plot_df = plot_df.sample(200_000, random_state=42)

            parts = ax.violinplot(
                [plot_df.loc[plot_df["truth_label"] == "TP", feat].values,
                 plot_df.loc[plot_df["truth_label"] == "FP", feat].values],
                positions=[0, 1], showmedians=True, showextrema=False,
            )
            for pc, color in zip(parts["bodies"], [COLOR_TP, COLOR_FP]):
                pc.set_facecolor(color)
                pc.set_alpha(0.6)
            parts["cmedians"].set_color("black")
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["TP", "FP"])
            ax.set_title(feat, fontsize=11, fontweight="bold")

            # Annotate medians
            tp_med = df.loc[df["truth_label"] == "TP", feat].median()
            fp_med = df.loc[df["truth_label"] == "FP", feat].median()
            ax.text(0.97, 0.95,
                    f"TP med={tp_med:.4f}\nFP med={fp_med:.4f}",
                    transform=ax.transAxes, fontsize=7, va="top", ha="right",
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

    fig.suptitle("Fig01: Methylation Feature Distributions by Truth Label (TP vs FP)",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig01_methyl_feature_distributions_by_truth.png")

    # ------------------------------------------------------------------
    # Fig02: AUC ranking (paired vs TO, grouped bar)
    # ------------------------------------------------------------------
    print("[O05] Generating fig02 ...")
    all_features_for_auc = METHYL_CONTINUOUS + METHYL_BINARY
    auc_pivot = auc_df.pivot(index="feature", columns="mode", values="auc")
    # Reorder by paired AUC descending
    if "paired" in auc_pivot.columns:
        auc_pivot = auc_pivot.sort_values("paired", ascending=False)

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(auc_pivot))
    width = 0.35
    bars_paired = ax.bar(x - width / 2, auc_pivot.get("paired", []), width,
                         label="Paired", color=COLOR_PAIRED, alpha=0.85, edgecolor="white")
    bars_to = ax.bar(x + width / 2, auc_pivot.get("to", []), width,
                     label="TO", color=COLOR_TO, alpha=0.85, edgecolor="white")
    ax.axhline(y=0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7, label="Random (0.5)")
    ax.set_xticks(x)
    ax.set_xticklabels(auc_pivot.index, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("AUC (TP=1 vs FP=0)")
    ax.set_ylim(0.35, 0.75)
    ax.legend(fontsize=9)
    # Annotate values
    for bar in bars_paired:
        h = bar.get_height()
        if not np.isnan(h):
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.005,
                    f"{h:.3f}", ha="center", fontsize=7, color=COLOR_PAIRED)
    for bar in bars_to:
        h = bar.get_height()
        if not np.isnan(h):
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.005,
                    f"{h:.3f}", ha="center", fontsize=7, color=COLOR_TO)

    fig.suptitle("Fig02: Methylation Feature AUC Ranking (Paired vs TO)",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig02_methyl_feature_auc_ranking.png")

    # ------------------------------------------------------------------
    # Fig03: PairwiseMedianDist by sample x mode (facet violin)
    # ------------------------------------------------------------------
    print("[O05] Generating fig03 ...")
    feat03 = "PairwiseMedianDist"
    samples_present = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    n_samp = len(samples_present)

    fig, axes = plt.subplots(2, 1, figsize=(max(14, n_samp * 2), 10), sharex=True)
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[(df["mode"] == mode) & (df["sample"].isin(samples_present))]
        if len(dm) == 0:
            ax.set_title(f"{MODE_LABELS[mode]} -- No Data")
            continue

        plot_data = []
        positions_tp = []
        positions_fp = []
        for j, samp in enumerate(samples_present):
            ds = dm[dm["sample"] == samp]
            tp_v = ds.loc[ds["truth_label"] == "TP", feat03].dropna().values
            fp_v = ds.loc[ds["truth_label"] == "FP", feat03].dropna().values
            if len(tp_v) > 0:
                parts = ax.violinplot([tp_v], positions=[j * 2 - 0.3],
                                      showmedians=True, showextrema=False, widths=0.5)
                for pc in parts["bodies"]:
                    pc.set_facecolor(COLOR_TP)
                    pc.set_alpha(0.6)
                parts["cmedians"].set_color("black")
            if len(fp_v) > 0:
                parts = ax.violinplot([fp_v], positions=[j * 2 + 0.3],
                                      showmedians=True, showextrema=False, widths=0.5)
                for pc in parts["bodies"]:
                    pc.set_facecolor(COLOR_FP)
                    pc.set_alpha(0.6)
                parts["cmedians"].set_color("black")

        ax.set_xticks([j * 2 for j in range(n_samp)])
        ax.set_xticklabels(samples_present, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel(feat03)
        ax.set_title(f"{MODE_LABELS[mode]} -- {feat03} by Sample (TP=blue, FP=red)", fontsize=12)
        ax.legend(
            handles=[mpatches.Patch(color=COLOR_TP, alpha=0.6, label="TP"),
                     mpatches.Patch(color=COLOR_FP, alpha=0.6, label="FP")],
            fontsize=9, loc="upper right",
        )

    fig.suptitle("Fig03: PairwiseMedianDist Distribution by Sample and Mode",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig03_pairwise_dist_by_sample_mode.png")

    # ------------------------------------------------------------------
    # Fig04: AlleleDelta by AF bin (TP vs FP, faceted by mode)
    # ------------------------------------------------------------------
    print("[O05] Generating fig04 ...")
    feat04 = "AlleleDelta"
    af_bins_present = [b for b in AF_BIN_ORDER if b in df["af_bin"].unique()]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]

        medians_tp = []
        medians_fp = []
        ci_tp_lo, ci_tp_hi = [], []
        ci_fp_lo, ci_fp_hi = [], []
        for afb in af_bins_present:
            dsub = dm[dm["af_bin"] == afb]
            tp_v = dsub.loc[dsub["truth_label"] == "TP", feat04].dropna().values
            fp_v = dsub.loc[dsub["truth_label"] == "FP", feat04].dropna().values

            tp_med = float(np.median(tp_v)) if len(tp_v) > 0 else np.nan
            fp_med = float(np.median(fp_v)) if len(fp_v) > 0 else np.nan
            medians_tp.append(tp_med)
            medians_fp.append(fp_med)

            # Bootstrap CI on median
            if len(tp_v) >= 10:
                _, lo, hi = bootstrap_ci(np.median, tp_v, n_boot=500)
                ci_tp_lo.append(lo)
                ci_tp_hi.append(hi)
            else:
                ci_tp_lo.append(np.nan)
                ci_tp_hi.append(np.nan)

            if len(fp_v) >= 10:
                _, lo, hi = bootstrap_ci(np.median, fp_v, n_boot=500)
                ci_fp_lo.append(lo)
                ci_fp_hi.append(hi)
            else:
                ci_fp_lo.append(np.nan)
                ci_fp_hi.append(np.nan)

        x = np.arange(len(af_bins_present))
        medians_tp = np.array(medians_tp)
        medians_fp = np.array(medians_fp)

        ax.plot(x, medians_tp, "o-", color=COLOR_TP, label="TP median", linewidth=2, markersize=5)
        ax.fill_between(x, ci_tp_lo, ci_tp_hi, color=COLOR_TP, alpha=0.15)
        ax.plot(x, medians_fp, "s--", color=COLOR_FP, label="FP median", linewidth=2, markersize=5)
        ax.fill_between(x, ci_fp_lo, ci_fp_hi, color=COLOR_FP, alpha=0.15)

        ax.set_xticks(x)
        ax.set_xticklabels(af_bins_present, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel(f"{feat04} (median)")
        ax.set_title(f"{MODE_LABELS[mode]} -- {feat04} by AF Bin", fontsize=12)
        ax.legend(fontsize=9)
        ax.axhline(y=0, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)

    fig.suptitle("Fig04: AlleleDelta by AF Bin (TP vs FP, with 95% CI)",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig04_allele_delta_by_af_bin.png")

    # ------------------------------------------------------------------
    # Fig05: CramersV vs GlobalP scatter (truth-coloured, mode facets)
    # ------------------------------------------------------------------
    print("[O05] Generating fig05 ...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]
        # Subsample for scatter readability
        if len(dm) > 30_000:
            dm_plot = dm.sample(30_000, random_state=42)
        else:
            dm_plot = dm

        for tl, color in [("FP", COLOR_FP), ("TP", COLOR_TP)]:
            dsub = dm_plot[dm_plot["truth_label"] == tl]
            ax.scatter(dsub["CramersV"], dsub["GlobalP"],
                       c=color, alpha=0.15, s=5, label=f"{tl} (n={len(dsub):,})",
                       rasterized=True)

        ax.set_xlabel("CramersV")
        ax.set_ylabel("GlobalP")
        ax.set_title(f"{MODE_LABELS[mode]}", fontsize=12)
        ax.legend(fontsize=9, markerscale=3)
        ax.set_xlim(-0.02, min(1.0, dm["CramersV"].quantile(0.99) * 1.3))
        ax.set_ylim(-0.02, 1.05)

    fig.suptitle("Fig05: CramersV vs GlobalP Scatter (Truth-Coloured, Paired vs TO)",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig05_cramersv_vs_globalp_scatter.png")

    # ------------------------------------------------------------------
    # Fig06: PassedGating rate heatmap (sample x truth, mode facets)
    # ------------------------------------------------------------------
    print("[O05] Generating fig06 ...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm_gating = gating_df[gating_df["mode"] == mode].copy()
        if len(dm_gating) == 0:
            ax.set_title(f"{MODE_LABELS[mode]} -- No Data")
            ax.axis("off")
            continue

        pivot = dm_gating.pivot(index="sample", columns="truth_label", values="passed_gating_rate")
        # Reorder
        pivot = pivot.reindex(index=[s for s in SAMPLE_ORDER if s in pivot.index])
        if "TP" not in pivot.columns:
            pivot["TP"] = np.nan
        if "FP" not in pivot.columns:
            pivot["FP"] = np.nan
        pivot = pivot[["TP", "FP"]]

        sns.heatmap(pivot, annot=True, fmt=".3f", cmap="YlOrRd", ax=ax,
                    vmin=0, vmax=1, linewidths=0.5, cbar_kws={"shrink": 0.8})
        ax.set_title(f"{MODE_LABELS[mode]} -- PassedGating Rate", fontsize=12)
        ax.set_ylabel("Sample")
        ax.set_xlabel("Truth Label")

    fig.suptitle("Fig06: PassedGating Rate by Sample and Truth Label",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig06_passed_gating_rate_by_sample_truth.png")

    # ------------------------------------------------------------------
    # Fig07: HPMergedSig vs truth cross-table (paired vs TO)
    # ------------------------------------------------------------------
    print("[O05] Generating fig07 ...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]
        ct = pd.crosstab(dm["truth_label"], dm["HPMergedSig"], normalize="index")
        ct = ct.reindex(index=["TP", "FP"])
        ct.columns = [f"Sig={c}" for c in ct.columns]
        ct.plot(kind="bar", ax=ax, color=["#90A4AE", "#FF7043"], alpha=0.8,
                edgecolor="white", width=0.6)
        ax.set_title(f"{MODE_LABELS[mode]} -- HPMergedSig Rate", fontsize=12)
        ax.set_ylabel("Proportion")
        ax.set_xlabel("")
        ax.set_xticklabels(["TP", "FP"], rotation=0)
        ax.legend(fontsize=9)

        # Annotate raw counts
        ct_raw = pd.crosstab(dm["truth_label"], dm["HPMergedSig"])
        for j, tl in enumerate(["TP", "FP"]):
            if tl in ct_raw.index:
                true_n = ct_raw.loc[tl, True] if True in ct_raw.columns else 0
                false_n = ct_raw.loc[tl, False] if False in ct_raw.columns else 0
                ax.text(j, 0.02, f"True={true_n:,}\nFalse={false_n:,}",
                        ha="center", fontsize=7, va="bottom",
                        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

    fig.suptitle("Fig07: HPMergedSig vs Truth Label (Paired vs TO)",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig07_hp_merged_sig_vs_truth.png")

    # ------------------------------------------------------------------
    # Fig08: Methylation feature correlation with caller features
    # ------------------------------------------------------------------
    print("[O05] Generating fig08 ...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        corr_mat = corr_data[mode]
        # Select methyl rows x caller columns + truth_bin
        methyl_cols_avail = [c for c in METHYL_CONTINUOUS if c in corr_mat.index]
        caller_cols_avail = [c for c in CALLER_FEATURES if c in corr_mat.columns]
        all_cols = methyl_cols_avail + caller_cols_avail
        sub_corr = corr_mat.loc[all_cols, all_cols]

        mask = np.triu(np.ones_like(sub_corr, dtype=bool), k=1)
        sns.heatmap(sub_corr, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r",
                    center=0, vmin=-1, vmax=1, ax=ax,
                    linewidths=0.5, cbar_kws={"shrink": 0.8},
                    annot_kws={"fontsize": 7})
        ax.set_title(f"{MODE_LABELS[mode]} -- Spearman Correlation", fontsize=12)
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.tick_params(axis="y", rotation=0, labelsize=8)

    fig.suptitle("Fig08: Methylation x Caller Feature Spearman Correlation Matrix",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig08_methyl_feature_correlation_with_caller.png")

    # ------------------------------------------------------------------
    # Fig09: Methylation features by LOH status (violin, paired vs TO)
    # ------------------------------------------------------------------
    print("[O05] Generating fig09 ...")
    loh_feats_plot = ["CramersV", "HeuristicScore", "AlleleDelta", "HPMergedDelta"]
    fig, axes = plt.subplots(2, len(loh_feats_plot), figsize=(20, 10))
    for i, mode in enumerate(MODE_ORDER):
        dm = df[df["mode"] == mode]
        for j, feat in enumerate(loh_feats_plot):
            ax = axes[i, j]
            if feat not in dm.columns:
                ax.set_title(f"{feat} (N/A)")
                ax.axis("off")
                continue

            groups = []
            labels = []
            colors = []
            for loh_val in [False, True]:
                for tl in ["TP", "FP"]:
                    vals = dm.loc[
                        (dm["Potential_LOH"] == loh_val) & (dm["truth_label"] == tl), feat
                    ].dropna().values
                    if len(vals) > 50000:
                        vals = np.random.RandomState(42).choice(vals, 50000, replace=False)
                    groups.append(vals)
                    labels.append(f"{LOH_LABELS[loh_val]}\n{tl}")
                    colors.append(COLOR_TP if tl == "TP" else COLOR_FP)

            positions = np.arange(len(groups))
            for k, (vals, col) in enumerate(zip(groups, colors)):
                if len(vals) < 5:
                    continue
                parts = ax.violinplot([vals], positions=[k],
                                      showmedians=True, showextrema=False, widths=0.7)
                for pc in parts["bodies"]:
                    pc.set_facecolor(col)
                    pc.set_alpha(0.5)
                parts["cmedians"].set_color("black")

            ax.set_xticks(positions)
            ax.set_xticklabels(labels, fontsize=7)
            title_str = feat if i == 0 else ""
            ax.set_title(title_str, fontsize=10, fontweight="bold")
            if j == 0:
                ax.set_ylabel(f"{MODE_LABELS[mode]}", fontsize=10)

    fig.suptitle("Fig09: Methylation Features by LOH Status (Paired vs TO, TP vs FP)",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    # Add legend
    fig.legend(
        handles=[mpatches.Patch(color=COLOR_TP, alpha=0.5, label="TP"),
                 mpatches.Patch(color=COLOR_FP, alpha=0.5, label="FP")],
        loc="upper right", fontsize=9,
    )
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig09_methyl_feature_by_loh_status.png")

    # ------------------------------------------------------------------
    # Fig10: Cohen's d forest plot (all methyl features, paired vs TO)
    # ------------------------------------------------------------------
    print("[O05] Generating fig10 ...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        edf = effect_df[effect_df["mode"] == mode].copy()
        edf = edf.sort_values("cohens_d", ascending=True)

        y_pos = np.arange(len(edf))
        ax.errorbar(edf["cohens_d"], y_pos,
                    xerr=[edf["cohens_d"] - edf["cohens_d_ci_lo"],
                          edf["cohens_d_ci_hi"] - edf["cohens_d"]],
                    fmt="o", color=COLOR_PAIRED if mode == "paired" else COLOR_TO,
                    capsize=4, markersize=6, linewidth=1.5)
        ax.axvline(x=0, color="gray", linestyle="--", linewidth=1, alpha=0.7)
        ax.axvline(x=0.2, color="green", linestyle=":", linewidth=0.8, alpha=0.4)
        ax.axvline(x=-0.2, color="green", linestyle=":", linewidth=0.8, alpha=0.4)
        ax.axvline(x=0.5, color="orange", linestyle=":", linewidth=0.8, alpha=0.4)
        ax.axvline(x=-0.5, color="orange", linestyle=":", linewidth=0.8, alpha=0.4)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(edf["feature"], fontsize=9)
        ax.set_xlabel("Cohen's d (TP - FP)")
        ax.set_title(f"{MODE_LABELS[mode]}", fontsize=12)

        # Annotate each point
        for j, (_, row) in enumerate(edf.iterrows()):
            ax.text(row["cohens_d"] + 0.02, j,
                    f"d={row['cohens_d']:.3f}",
                    fontsize=7, va="center")

    fig.suptitle("Fig10: Cohen's d Effect Size Forest Plot (Methylation Features, TP vs FP)",
                 fontsize=14, y=1.02)
    add_caption(fig, caption_base + " | Green dotted: |d|=0.2 (small), Orange: |d|=0.5 (medium)")
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig10_methyl_effect_size_forest_plot.png")

    # ------------------------------------------------------------------
    # Round context
    # ------------------------------------------------------------------
    write_round_context(
        OUT_DIR,
        observation_id=OBSERVATION_ID,
        title=OBSERVATION_TITLE,
        description=(
            "Quantifies ISM-derived methylation features' independent "
            "discriminative power for TP vs FP, stratified by paired vs TO mode, "
            "sample, AF bin, and LOH status."
        ),
        script_path=SCRIPT_PATH,
        data_sources=[{
            "file": str(USECOLS),
            "description": "Master dataset columns used",
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "n_figures": 10,
            "methyl_features": METHYL_FEATURES,
            "caller_features": CALLER_FEATURES,
        },
    )

    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 70}")
    print(f"  {OBSERVATION_ID} complete in {elapsed.total_seconds():.1f} seconds")
    print(f"  Output: {OUT_DIR}")
    print(f"={'=' * 70}")


if __name__ == "__main__":
    main()
