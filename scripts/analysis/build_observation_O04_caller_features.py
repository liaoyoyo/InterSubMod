#!/usr/bin/env python3
"""
O04: Caller Feature TP/FP Distribution Observation
====================================================

Systematic observation of caller-level features (AF, GQ, DP, AD, FILTER)
and their discriminative power for TP vs FP, separated by paired vs TO modes.

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
    compute_auc, compute_effect_size, mannwhitney_test,
    encode_truth_binary, write_round_context, write_feature_statistics,
    ensure_dir, safe_div, format_p, format_ci,
    SAMPLE_ORDER, MODE_ORDER, MASTER_DATASET_PATH, OUTPUT_ROOT,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
    TRUTH_PALETTE, MODE_PALETTE,
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

OBSERVATION_ID = "O04"
OBSERVATION_TITLE = "Caller Feature TP/FP Distribution Observation"
SCRIPT_PATH = "scripts/analysis/build_observation_O04_caller_features.py"

OUT_DIR = OUTPUT_ROOT / "20260401_O04_caller_features"
FIG_DIR = OUT_DIR
DATA_DIR = OUT_DIR / "data"

CALLER_NUMERIC_FEATURES = ["caller_af", "caller_gq", "caller_dp", "caller_ad_ref", "caller_ad_alt"]
USECOLS = [
    "truth_label", "sample", "mode", "variant_key",
    "caller_af", "caller_gq", "caller_dp",
    "caller_ad_ref", "caller_ad_alt", "caller_filter", "af_bin",
    "num_reads", "num_cpgs", "quality_score",
]

MODE_LABELS = {"paired": "Paired", "to": "Tumor-Only (TO)"}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compute_ad_ratio(df: pd.DataFrame) -> pd.Series:
    """Compute ad_alt / (ad_ref + ad_alt), returning NaN when denominator is 0."""
    denom = df["caller_ad_ref"] + df["caller_ad_alt"]
    return df["caller_ad_alt"] / denom.replace(0, np.nan)


def manual_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Fallback AUC computation via trapezoidal rule on ROC."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    mask = np.isfinite(y_score) & np.isfinite(y_true)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 10:
        return float("nan")
    # Sort by score descending
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
    print(f"[O04] Filtered to TP/FP: {len(df):,} rows (dropped {total_rows - len(df):,})")

    # Derived column: ad_ratio
    df["ad_ratio"] = compute_ad_ratio(df)

    # Binary truth label
    df["truth_bin"] = encode_truth_binary(df)

    # Basic stats
    mode_truth_counts = df.groupby(["mode", "truth_label"]).size().unstack(fill_value=0)
    print(f"\n[O04] Mode x Truth distribution:")
    print(mode_truth_counts.to_string())

    # ------------------------------------------------------------------
    # Data Output: source_summary.tsv
    # ------------------------------------------------------------------
    summary_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for tl in ["TP", "FP"]:
            dt = dm[dm["truth_label"] == tl]
            summary_rows.append({
                "mode": mode,
                "truth_label": tl,
                "n_regions": len(dt),
                "n_samples": dt["sample"].nunique(),
                "caller_af_median": dt["caller_af"].median(),
                "caller_gq_median": dt["caller_gq"].median(),
                "caller_dp_median": dt["caller_dp"].median(),
                "ad_ratio_median": dt["ad_ratio"].median(),
                "caller_filter_pass_pct": (dt["caller_filter"] == "PASS").mean() * 100,
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(DATA_DIR / "source_summary.tsv", sep="\t", index=False)
    print(f"[O04] Wrote source_summary.tsv ({len(summary_df)} rows)")

    # ------------------------------------------------------------------
    # Data Output: feature_statistics.tsv
    # ------------------------------------------------------------------
    feature_cols = CALLER_NUMERIC_FEATURES + ["ad_ratio"]
    write_feature_statistics(df, feature_cols, DATA_DIR / "feature_statistics.tsv")

    # ------------------------------------------------------------------
    # Statistical Tests: Mann-Whitney U + Cohen's d per feature per mode
    # ------------------------------------------------------------------
    print("\n[O04] Computing statistical tests ...")
    effect_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        tp_mask = dm["truth_label"] == "TP"
        fp_mask = dm["truth_label"] == "FP"
        for feat in feature_cols:
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
                "mwu_statistic": mwu["statistic"],
                "mwu_p_value": mwu["p_value"],
                "mwu_effect_r": mwu["effect_size_r"],
                "cohens_d": es["d"],
                "cohens_d_ci_lo": es["ci_lo"],
                "cohens_d_ci_hi": es["ci_hi"],
            })
    effect_df = pd.DataFrame(effect_rows)
    effect_df.to_csv(DATA_DIR / "caller_effect_size.tsv", sep="\t", index=False)
    print(f"[O04] Wrote caller_effect_size.tsv ({len(effect_df)} rows)")

    # ------------------------------------------------------------------
    # AUC by mode
    # ------------------------------------------------------------------
    print("\n[O04] Computing AUC by mode ...")
    auc_rows = []
    auc_features = feature_cols  # includes ad_ratio
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        y_true = dm["truth_bin"].values
        for feat in auc_features:
            y_score = dm[feat].values
            try:
                auc_val = compute_auc(y_true, y_score)
            except Exception:
                auc_val = manual_auc(y_true, y_score)
            auc_rows.append({
                "mode": mode,
                "feature": feat,
                "auc": round(auc_val, 4) if not np.isnan(auc_val) else np.nan,
                "n_total": len(dm),
                "n_valid": int(np.sum(np.isfinite(y_score) & np.isfinite(y_true))),
            })
    auc_df = pd.DataFrame(auc_rows)
    auc_df.to_csv(DATA_DIR / "caller_auc_by_mode.tsv", sep="\t", index=False)
    print(f"[O04] Wrote caller_auc_by_mode.tsv ({len(auc_df)} rows)")

    # ------------------------------------------------------------------
    # AF bin TP/FP rates
    # ------------------------------------------------------------------
    print("\n[O04] Computing AF bin TP/FP rates ...")
    afbin_rows = []
    for mode in MODE_ORDER:
        dm = df[df["mode"] == mode]
        for afb in sorted(dm["af_bin"].dropna().unique()):
            dsub = dm[dm["af_bin"] == afb]
            n_tp = (dsub["truth_label"] == "TP").sum()
            n_fp = (dsub["truth_label"] == "FP").sum()
            n_tot = n_tp + n_fp
            afbin_rows.append({
                "mode": mode,
                "af_bin": afb,
                "n_total": n_tot,
                "n_tp": n_tp,
                "n_fp": n_fp,
                "tp_rate": safe_div(n_tp, n_tot),
                "fp_rate": safe_div(n_fp, n_tot),
            })
    afbin_df = pd.DataFrame(afbin_rows)
    afbin_df.to_csv(DATA_DIR / "af_bin_tp_fp_rates.tsv", sep="\t", index=False)
    print(f"[O04] Wrote af_bin_tp_fp_rates.tsv ({len(afbin_df)} rows)")

    # ==================================================================
    # FIGURES
    # ==================================================================

    caption_base = (
        f"O04 Caller Features | N={len(df):,} regions (TP+FP) | "
        f"Source: all_region_rows.tsv.gz | {datetime.now().strftime('%Y-%m-%d')}"
    )

    # ------------------------------------------------------------------
    # Fig 01: AF distribution by truth (paired vs TO rows)
    # ------------------------------------------------------------------
    print("\n[O04] Generating fig01 ...")
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]
        bins = np.linspace(0, 1, 51)
        for tl, color in TRUTH_PALETTE.items():
            vals = dm.loc[dm["truth_label"] == tl, "caller_af"].dropna()
            ax.hist(vals, bins=bins, alpha=0.55, color=color, label=f"{tl} (n={len(vals):,})",
                    density=True, edgecolor="white", linewidth=0.3)
        ax.set_title(f"{MODE_LABELS[mode]} -- Caller AF Distribution", fontsize=12)
        ax.set_ylabel("Density")
        ax.legend(fontsize=9)
    axes[-1].set_xlabel("Caller AF (Allele Frequency)")
    fig.suptitle("Fig01: Caller AF Distribution by Truth Label", fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig01_caller_af_distribution_by_truth.png")

    # ------------------------------------------------------------------
    # Fig 02: GQ distribution by truth (paired vs TO)
    # ------------------------------------------------------------------
    print("[O04] Generating fig02 ...")
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]
        for tl, color in TRUTH_PALETTE.items():
            vals = dm.loc[dm["truth_label"] == tl, "caller_gq"].dropna()
            ax.hist(vals, bins=50, alpha=0.55, color=color, label=f"{tl} (n={len(vals):,})",
                    density=True, edgecolor="white", linewidth=0.3)
        ax.set_title(f"{MODE_LABELS[mode]} -- Caller GQ Distribution", fontsize=12)
        ax.set_ylabel("Density")
        ax.legend(fontsize=9)
    axes[-1].set_xlabel("Caller GQ (Genotype Quality)")
    fig.suptitle("Fig02: Caller GQ Distribution by Truth Label", fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig02_caller_gq_distribution_by_truth.png")

    # ------------------------------------------------------------------
    # Fig 03: DP distribution by truth (log scale, paired vs TO)
    # ------------------------------------------------------------------
    print("[O04] Generating fig03 ...")
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]
        for tl, color in TRUTH_PALETTE.items():
            vals = dm.loc[dm["truth_label"] == tl, "caller_dp"].dropna()
            vals_pos = vals[vals > 0]
            bins_log = np.logspace(np.log10(max(vals_pos.min(), 1)), np.log10(vals_pos.max()), 50)
            ax.hist(vals_pos, bins=bins_log, alpha=0.55, color=color,
                    label=f"{tl} (n={len(vals_pos):,})", density=True,
                    edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.set_title(f"{MODE_LABELS[mode]} -- Caller DP Distribution (log scale)", fontsize=12)
        ax.set_ylabel("Density")
        ax.legend(fontsize=9)
    axes[-1].set_xlabel("Caller DP (Read Depth)")
    fig.suptitle("Fig03: Caller DP Distribution by Truth Label", fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig03_caller_dp_distribution_by_truth.png")

    # ------------------------------------------------------------------
    # Fig 04: AUC bar chart by mode for each caller feature
    # ------------------------------------------------------------------
    print("[O04] Generating fig04 ...")
    fig, ax = plt.subplots(figsize=(10, 5))

    feature_labels = {
        "caller_af": "AF", "caller_gq": "GQ", "caller_dp": "DP",
        "caller_ad_ref": "AD_ref", "caller_ad_alt": "AD_alt", "ad_ratio": "AD_ratio",
    }
    auc_plot = auc_df.copy()
    auc_plot["feat_label"] = auc_plot["feature"].map(feature_labels)
    auc_plot = auc_plot.dropna(subset=["feat_label"])

    x = np.arange(len(feature_labels))
    width = 0.35
    for j, mode in enumerate(MODE_ORDER):
        sub = auc_plot[auc_plot["mode"] == mode].set_index("feat_label")
        vals = [sub.loc[fl, "auc"] if fl in sub.index else np.nan
                for fl in feature_labels.values()]
        bars = ax.bar(x + j * width - width / 2, vals, width,
                      label=MODE_LABELS[mode],
                      color=MODE_PALETTE[mode], alpha=0.8, edgecolor="white")
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                        f"{val:.3f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(list(feature_labels.values()), fontsize=10)
    ax.set_ylabel("AUC (TP=1 vs FP=0)")
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.legend(fontsize=10)
    ax.set_title("Fig04: Caller Feature AUC for TP/FP Discrimination (Paired vs TO)", fontsize=13)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig04_caller_feature_auc_by_mode.png")

    # ------------------------------------------------------------------
    # Fig 05: AF vs GQ scatter by truth (paired and TO subplots)
    # ------------------------------------------------------------------
    print("[O04] Generating fig05 ...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]
        # Downsample for readability if very large
        max_pts = 30000
        if len(dm) > max_pts:
            dm_plot = dm.sample(max_pts, random_state=42)
        else:
            dm_plot = dm
        for tl, color in [("FP", COLOR_FP), ("TP", COLOR_TP)]:
            sub = dm_plot[dm_plot["truth_label"] == tl]
            ax.scatter(sub["caller_af"], sub["caller_gq"],
                       alpha=0.05, s=6, c=color, label=tl, rasterized=True)
        ax.set_xlabel("Caller AF")
        ax.set_ylabel("Caller GQ")
        ax.set_title(f"{MODE_LABELS[mode]}", fontsize=12)
        # Custom legend with opaque markers
        handles = [
            mpatches.Patch(color=COLOR_TP, label=f"TP (n={len(dm[dm['truth_label']=='TP']):,})"),
            mpatches.Patch(color=COLOR_FP, label=f"FP (n={len(dm[dm['truth_label']=='FP']):,})"),
        ]
        ax.legend(handles=handles, fontsize=9, loc="lower right")
    fig.suptitle("Fig05: AF vs GQ Scatter by Truth Label", fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig05_af_vs_gq_scatter_by_truth.png")

    # ------------------------------------------------------------------
    # Fig 06: AF bin FP rate (paired vs TO)
    # ------------------------------------------------------------------
    print("[O04] Generating fig06 ...")
    fig, ax = plt.subplots(figsize=(12, 5))

    # Sort af_bins by numeric order
    afbin_plot = afbin_df.copy()
    # Parse bin start for sorting
    def bin_sort_key(b):
        try:
            return float(str(b).split("-")[0])
        except Exception:
            return 999.0
    afbin_plot["sort_key"] = afbin_plot["af_bin"].apply(bin_sort_key)
    afbin_plot = afbin_plot.sort_values("sort_key")

    all_bins = afbin_plot["af_bin"].unique()
    x = np.arange(len(all_bins))
    width = 0.35

    for j, mode in enumerate(MODE_ORDER):
        sub = afbin_plot[afbin_plot["mode"] == mode]
        # Align to all_bins
        fp_rates = []
        for b in all_bins:
            row = sub[sub["af_bin"] == b]
            fp_rates.append(row["fp_rate"].values[0] if len(row) > 0 else 0)
        ax.bar(x + j * width - width / 2, fp_rates, width,
               label=MODE_LABELS[mode], color=MODE_PALETTE[mode], alpha=0.8, edgecolor="white")

    ax.set_xticks(x)
    ax.set_xticklabels(all_bins, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("FP Rate (FP / Total)")
    ax.set_xlabel("AF Bin")
    ax.legend(fontsize=10)
    ax.set_title("Fig06: FP Rate by AF Bin (Paired vs TO)", fontsize=13)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig06_af_bin_tp_fp_rate.png")

    # ------------------------------------------------------------------
    # Fig 07: AD ratio distribution (TP/FP x paired/TO)
    # ------------------------------------------------------------------
    print("[O04] Generating fig07 ...")
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for i, mode in enumerate(MODE_ORDER):
        ax = axes[i]
        dm = df[df["mode"] == mode]
        bins = np.linspace(0, 1, 51)
        for tl, color in TRUTH_PALETTE.items():
            vals = dm.loc[dm["truth_label"] == tl, "ad_ratio"].dropna()
            ax.hist(vals, bins=bins, alpha=0.55, color=color,
                    label=f"{tl} (n={len(vals):,})", density=True,
                    edgecolor="white", linewidth=0.3)
        ax.set_title(f"{MODE_LABELS[mode]} -- AD Ratio Distribution", fontsize=12)
        ax.set_ylabel("Density")
        ax.legend(fontsize=9)
    axes[-1].set_xlabel("AD Ratio = AD_alt / (AD_ref + AD_alt)")
    fig.suptitle("Fig07: Caller AD Ratio Distribution by Truth Label", fontsize=14, y=1.02)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig07_caller_ad_ratio_distribution.png")

    # ------------------------------------------------------------------
    # Fig 08: Per-sample caller feature medians (TP vs FP), 7 rows x 2 cols
    # ------------------------------------------------------------------
    print("[O04] Generating fig08 ...")
    fig, axes = plt.subplots(7, 2, figsize=(12, 22))

    for row_idx, sample in enumerate(SAMPLE_ORDER):
        ds = df[df["sample"] == sample]
        for col_idx, (feat, feat_label) in enumerate([
            ("caller_af", "Caller AF (median)"),
            ("caller_gq", "Caller GQ (median)"),
        ]):
            ax = axes[row_idx, col_idx]
            bar_data = []
            for mode in MODE_ORDER:
                dm = ds[ds["mode"] == mode]
                for tl in ["TP", "FP"]:
                    dt = dm[dm["truth_label"] == tl]
                    med = dt[feat].median() if len(dt) > 0 else np.nan
                    bar_data.append({
                        "mode": mode, "truth": tl, "median": med,
                        "n": len(dt),
                    })
            bar_df = pd.DataFrame(bar_data)

            x_positions = np.arange(len(MODE_ORDER))
            w = 0.35
            for k, tl in enumerate(["TP", "FP"]):
                sub = bar_df[bar_df["truth"] == tl]
                vals = sub["median"].values
                counts = sub["n"].values
                color = COLOR_TP if tl == "TP" else COLOR_FP
                bars = ax.bar(x_positions + k * w - w / 2, vals, w,
                              color=color, alpha=0.8, label=tl, edgecolor="white")
                for bar, val, cnt in zip(bars, vals, counts):
                    if not np.isnan(val):
                        ax.text(bar.get_x() + bar.get_width() / 2,
                                bar.get_height() + 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0] + 1),
                                f"n={cnt}", ha="center", va="bottom", fontsize=6)

            ax.set_xticks(x_positions)
            ax.set_xticklabels([MODE_LABELS[m] for m in MODE_ORDER], fontsize=8)
            ax.set_ylabel(feat_label, fontsize=8)
            ax.set_title(f"{sample}", fontsize=9, fontweight="bold")
            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=7, loc="upper right")

    fig.suptitle("Fig08: Per-Sample Caller Feature Medians (TP vs FP)", fontsize=14, y=1.01)
    add_caption(fig, caption_base)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "fig08_caller_feature_per_sample.png")

    # ------------------------------------------------------------------
    # round_context.json
    # ------------------------------------------------------------------
    write_round_context(
        out_dir=OUT_DIR,
        observation_id=OBSERVATION_ID,
        title=OBSERVATION_TITLE,
        description=(
            "Systematic observation of caller-level features (AF, GQ, DP, AD, FILTER) "
            "and their discriminative power for TP vs FP classification, separated by "
            "paired and tumor-only (TO) calling modes."
        ),
        script_path=SCRIPT_PATH,
        data_sources=[{
            "name": "all_region_rows.tsv.gz",
            "path": str(MASTER_DATASET_PATH),
            "rows": total_rows,
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "figures_generated": 8,
            "features_analysed": feature_cols,
            "elapsed_seconds": (datetime.now() - start_time).total_seconds(),
        },
    )

    # ------------------------------------------------------------------
    # Summary print
    # ------------------------------------------------------------------
    elapsed = (datetime.now() - start_time).total_seconds()
    print(f"\n{'=' * 70}")
    print(f"  {OBSERVATION_ID} COMPLETE")
    print(f"  Output: {OUT_DIR}")
    print(f"  Figures: 8 | Data files: 5 | Elapsed: {elapsed:.1f}s")
    print(f"{'=' * 70}")

    # Return key data for report generation
    return {
        "df": df,
        "effect_df": effect_df,
        "auc_df": auc_df,
        "afbin_df": afbin_df,
        "summary_df": summary_df,
    }


if __name__ == "__main__":
    main()
