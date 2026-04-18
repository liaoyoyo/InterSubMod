#!/usr/bin/env python3
"""
O02: Quality Score Component Decomposition Observation
=======================================================

Systematic decomposition of Quality Score (QS) into its constituent penalties
and bonuses, analysing why QS fails in TO mode (AUC=0.497 ~ random) and which
individual components contribute to or harm discriminative power.

Part of the O1-O10 systematic observation framework for Phase 1A ML baseline.

Generated: 2026-04-01
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
    encode_truth_binary, reconstruct_qs_components,
    ensure_dir, safe_div, write_round_context, write_feature_statistics,
    markdown_table, OUTPUT_ROOT, SAMPLE_ORDER, MODE_ORDER, MASTER_DATASET_PATH,
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
from sklearn.metrics import roc_curve, confusion_matrix

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OBSERVATION_ID = "O02"
OBSERVATION_TITLE = "Quality Score Component Decomposition Observation"
SCRIPT_PATH = "scripts/analysis/build_observation_O02_quality_score_decompose.py"

OUT_DIR = OUTPUT_ROOT / "20260401_O02_quality_score_decompose"
FIG_DIR = OUT_DIR
DATA_DIR = OUT_DIR / "data"

# Columns needed from master dataset
USECOLS = [
    "truth_label", "sample", "mode", "variant_key",
    "Quality_Score", "Quality_Tier",
    "NumReads", "NumCpGs", "Coverage_Multiple", "Potential_LOH",
    "HPMergedSig", "AlleleSig", "CramersV",
    "num_reads", "num_cpgs",
]

# Component column names (from reconstruct_qs_components)
PENALTY_COLS = [
    "qs_penalty_low_reads",
    "qs_penalty_low_cpgs",
    "qs_penalty_cnv_loss",
    "qs_penalty_high_copy",
    "qs_penalty_loh",
]

BONUS_COLS = [
    "qs_bonus_dual_sig",
    "qs_bonus_strong_effect",
]

ALL_COMPONENT_COLS = PENALTY_COLS + BONUS_COLS

# Human-readable labels for components
COMPONENT_LABELS = {
    "qs_penalty_low_reads": "Low Reads\n(<30: -20, <50: -10)",
    "qs_penalty_low_cpgs": "Low CpGs\n(<20: -15, <30: -10)",
    "qs_penalty_cnv_loss": "CNV Loss\n(<0.5: -30, <0.8: -15)",
    "qs_penalty_high_copy": "High Copy\n(>2.0: -20)",
    "qs_penalty_loh": "LOH Penalty\n(LOH=True: -25)",
    "qs_bonus_dual_sig": "Dual Sig Bonus\n(HP+Allele: +10)",
    "qs_bonus_strong_effect": "Strong Effect Bonus\n(CV>=0.3: +15, >=0.2: +5)",
}

# Short labels for bar charts
COMPONENT_SHORT = {
    "qs_penalty_low_reads": "Low Reads",
    "qs_penalty_low_cpgs": "Low CpGs",
    "qs_penalty_cnv_loss": "CNV Loss",
    "qs_penalty_high_copy": "High Copy",
    "qs_penalty_loh": "LOH",
    "qs_bonus_dual_sig": "Dual Sig",
    "qs_bonus_strong_effect": "Strong Effect",
}


# ---------------------------------------------------------------------------
# Data Loading & Preparation
# ---------------------------------------------------------------------------

def load_and_prepare() -> pd.DataFrame:
    """Load master dataset and reconstruct QS components."""
    print(f"\n{'='*70}")
    print(f"  O02: Quality Score Component Decomposition")
    print(f"{'='*70}\n")

    df = load_master_dataset(usecols=USECOLS)

    # Filter to rows with truth labels
    df = df[df["truth_label"].isin(["TP", "FP"])].copy()
    print(f"[O02] After truth filter: {len(df):,} rows")

    # Encode truth binary
    df["truth_binary"] = encode_truth_binary(df)

    # Reconstruct QS components
    qs_df = reconstruct_qs_components(df)
    for col in qs_df.columns:
        df[col] = qs_df[col]

    # Verify reconstruction matches original
    qs_orig = df["Quality_Score"].astype(float)
    qs_recon = df["qs_reconstructed"]
    valid = qs_orig.notna() & qs_recon.notna()
    corr = qs_orig[valid].corr(qs_recon[valid])
    mae = (qs_orig[valid] - qs_recon[valid]).abs().mean()
    print(f"[O02] QS reconstruction: corr={corr:.4f}, MAE={mae:.4f}")

    # Mode counts
    for mode in MODE_ORDER:
        sub = df[df["mode"] == mode]
        tp_n = (sub["truth_label"] == "TP").sum()
        fp_n = (sub["truth_label"] == "FP").sum()
        print(f"[O02]   {mode}: {len(sub):,} rows (TP={tp_n:,}, FP={fp_n:,})")

    return df


# ---------------------------------------------------------------------------
# Fig01: QS Component Waterfall - Paired
# ---------------------------------------------------------------------------

def fig01_waterfall_paired(df: pd.DataFrame) -> Dict[str, Any]:
    """Waterfall chart of mean QS component contributions for paired mode."""
    paired = df[df["mode"] == "paired"]
    stats = {}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for idx, (truth, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        ax = axes[idx]
        sub = paired[paired["truth_label"] == truth]

        # Start from base=100
        base = 100.0
        labels = ["Base"]
        values = [base]
        bottoms = [0]
        colors = ["#E0E0E0"]

        running = base
        for comp in ALL_COMPONENT_COLS:
            mean_val = sub[comp].mean()
            labels.append(COMPONENT_SHORT[comp])
            values.append(mean_val)
            if mean_val >= 0:
                bottoms.append(running)
                colors.append("#4CAF50")  # green for bonus
            else:
                bottoms.append(running + mean_val)
                colors.append("#F44336")  # red for penalty
            running += mean_val

        # Final bar
        labels.append("Final QS")
        values.append(running)
        bottoms.append(0)
        colors.append(color)

        # Plot
        x_pos = np.arange(len(labels))
        bar_vals = [abs(v) for v in values]
        bars = ax.bar(x_pos, bar_vals, bottom=bottoms, color=colors,
                      edgecolor="white", linewidth=0.8, width=0.65)

        # Annotate values
        for i, (v, b) in enumerate(zip(values, bottoms)):
            if i == 0:  # base
                ax.text(i, b + abs(v) / 2, f"{v:.0f}", ha="center", va="center",
                        fontsize=8, fontweight="bold")
            elif i == len(labels) - 1:  # final
                ax.text(i, v / 2, f"{v:.1f}", ha="center", va="center",
                        fontsize=8, fontweight="bold", color="white")
            else:
                y_pos = b + abs(v) / 2 if v < 0 else b + v / 2
                sign = "+" if v > 0 else ""
                ax.text(i, y_pos, f"{sign}{v:.1f}", ha="center", va="center",
                        fontsize=7, fontweight="bold")

        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7.5)
        ax.set_title(f"Paired {truth} (n={len(sub):,})", fontsize=12, fontweight="bold")
        ax.set_ylabel("Quality Score" if idx == 0 else "")
        ax.set_ylim(0, 115)
        ax.axhline(50, color="gray", ls="--", alpha=0.3)

        stats[f"paired_{truth}_final_qs_mean"] = running
        stats[f"paired_{truth}_n"] = len(sub)
        for comp in ALL_COMPONENT_COLS:
            stats[f"paired_{truth}_{comp}_mean"] = float(sub[comp].mean())

    fig.suptitle("Fig01: QS Component Waterfall — Paired Mode", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0.05, 1, 0.93])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | Paired mode only | Base=100, penalties subtract, bonuses add | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig01_qs_component_waterfall_paired.png")

    return stats


# ---------------------------------------------------------------------------
# Fig02: QS Component Waterfall - TO
# ---------------------------------------------------------------------------

def fig02_waterfall_to(df: pd.DataFrame) -> Dict[str, Any]:
    """Waterfall chart of mean QS component contributions for TO mode."""
    to_df = df[df["mode"] == "to"]
    stats = {}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for idx, (truth, color) in enumerate([("TP", COLOR_TP), ("FP", COLOR_FP)]):
        ax = axes[idx]
        sub = to_df[to_df["truth_label"] == truth]

        base = 100.0
        labels = ["Base"]
        values = [base]
        bottoms = [0]
        colors = ["#E0E0E0"]

        running = base
        for comp in ALL_COMPONENT_COLS:
            mean_val = sub[comp].mean()
            labels.append(COMPONENT_SHORT[comp])
            values.append(mean_val)
            if mean_val >= 0:
                bottoms.append(running)
                colors.append("#4CAF50")
            else:
                bottoms.append(running + mean_val)
                colors.append("#F44336")
            running += mean_val

        labels.append("Final QS")
        values.append(running)
        bottoms.append(0)
        colors.append(color)

        x_pos = np.arange(len(labels))
        bar_vals = [abs(v) for v in values]
        bars = ax.bar(x_pos, bar_vals, bottom=bottoms, color=colors,
                      edgecolor="white", linewidth=0.8, width=0.65)

        for i, (v, b) in enumerate(zip(values, bottoms)):
            if i == 0:
                ax.text(i, b + abs(v) / 2, f"{v:.0f}", ha="center", va="center",
                        fontsize=8, fontweight="bold")
            elif i == len(labels) - 1:
                ax.text(i, v / 2, f"{v:.1f}", ha="center", va="center",
                        fontsize=8, fontweight="bold", color="white")
            else:
                y_pos = b + abs(v) / 2 if v < 0 else b + v / 2
                sign = "+" if v > 0 else ""
                ax.text(i, y_pos, f"{sign}{v:.1f}", ha="center", va="center",
                        fontsize=7, fontweight="bold")

        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7.5)
        ax.set_title(f"TO {truth} (n={len(sub):,})", fontsize=12, fontweight="bold")
        ax.set_ylabel("Quality Score" if idx == 0 else "")
        ax.set_ylim(0, 115)
        ax.axhline(50, color="gray", ls="--", alpha=0.3)

        stats[f"to_{truth}_final_qs_mean"] = running
        stats[f"to_{truth}_n"] = len(sub)
        for comp in ALL_COMPONENT_COLS:
            stats[f"to_{truth}_{comp}_mean"] = float(sub[comp].mean())

    fig.suptitle("Fig02: QS Component Waterfall — TO Mode", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0.05, 1, 0.93])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | TO mode only | Base=100, penalties subtract, bonuses add | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig02_qs_component_waterfall_to.png")

    return stats


# ---------------------------------------------------------------------------
# Fig03: LOH Penalty Trigger Rate
# ---------------------------------------------------------------------------

def fig03_loh_trigger_rate(df: pd.DataFrame) -> Dict[str, Any]:
    """LOH penalty trigger rate across TP/FP x paired/TO."""
    stats = {}

    fig, ax = plt.subplots(figsize=(10, 6))

    groups = []
    for mode in MODE_ORDER:
        for truth in ["TP", "FP"]:
            sub = df[(df["mode"] == mode) & (df["truth_label"] == truth)]
            triggered = (sub["qs_penalty_loh"] < 0).sum()
            rate = safe_div(triggered, len(sub)) * 100
            groups.append({
                "mode": mode, "truth": truth, "n": len(sub),
                "triggered": triggered, "rate": rate,
                "label": f"{mode.upper()} {truth}",
            })
            stats[f"loh_trigger_{mode}_{truth}_rate"] = rate
            stats[f"loh_trigger_{mode}_{truth}_n"] = triggered
            stats[f"loh_trigger_{mode}_{truth}_total"] = len(sub)

    labels = [g["label"] for g in groups]
    rates = [g["rate"] for g in groups]
    colors_list = [
        COLOR_TP if g["truth"] == "TP" else COLOR_FP for g in groups
    ]
    hatches = ["" if g["mode"] == "paired" else "///" for g in groups]

    bars = ax.bar(labels, rates, color=colors_list, edgecolor="black", linewidth=0.8)
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)

    for bar, g in zip(bars, groups):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 1,
                f"{g['rate']:.1f}%\n({g['triggered']:,}/{g['n']:,})",
                ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_ylabel("LOH Penalty Trigger Rate (%)", fontsize=11)
    ax.set_title("Fig03: LOH Penalty (-25) Trigger Rate by Truth x Mode", fontsize=13, fontweight="bold")
    ax.set_ylim(0, max(rates) * 1.3)

    # Chi-squared test: mode x LOH trigger
    for mode in MODE_ORDER:
        sub = df[df["mode"] == mode]
        tp_trig = ((sub["truth_label"] == "TP") & (sub["qs_penalty_loh"] < 0)).sum()
        fp_trig = ((sub["truth_label"] == "FP") & (sub["qs_penalty_loh"] < 0)).sum()
        tp_not = ((sub["truth_label"] == "TP") & (sub["qs_penalty_loh"] >= 0)).sum()
        fp_not = ((sub["truth_label"] == "FP") & (sub["qs_penalty_loh"] >= 0)).sum()
        ct = np.array([[tp_trig, tp_not], [fp_trig, fp_not]])
        chi2_res = chi2_test(ct)
        stats[f"loh_chi2_{mode}_chi2"] = chi2_res["chi2"]
        stats[f"loh_chi2_{mode}_p"] = chi2_res["p_value"]
        stats[f"loh_chi2_{mode}_cramers_v"] = chi2_res["cramers_v"]

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_TP, label="TP"),
        mpatches.Patch(facecolor=COLOR_FP, label="FP"),
        mpatches.Patch(facecolor="white", edgecolor="black", hatch="///", label="TO mode"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9)

    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | LOH penalty = -25 when Potential_LOH=True | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig03_qs_loh_penalty_trigger_rate.png")

    return stats


# ---------------------------------------------------------------------------
# Fig04: ROC Comparison With/Without LOH Penalty
# ---------------------------------------------------------------------------

def fig04_roc_without_loh(df: pd.DataFrame) -> Dict[str, Any]:
    """ROC curves comparing QS with and without LOH penalty."""
    stats = {}

    # Compute QS without LOH penalty
    df["qs_no_loh"] = (df["qs_reconstructed"] - df["qs_penalty_loh"]).clip(0, 100)
    # Also compute QS with only LOH (all other components zeroed out)
    df["qs_only_loh"] = (100.0 + df["qs_penalty_loh"]).clip(0, 100)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub = df[df["mode"] == mode].copy()
        y_true = sub["truth_binary"].values

        # Original QS
        auc_orig = compute_auc(y_true, sub["qs_reconstructed"].values)
        fpr_o, tpr_o, _ = roc_curve(y_true, sub["qs_reconstructed"].values)
        ax.plot(fpr_o, tpr_o, color=COLOR_PAIRED if mode == "paired" else COLOR_TO,
                lw=2, label=f"Full QS (AUC={auc_orig:.3f})")

        # QS without LOH
        auc_no_loh = compute_auc(y_true, sub["qs_no_loh"].values)
        fpr_nl, tpr_nl, _ = roc_curve(y_true, sub["qs_no_loh"].values)
        ax.plot(fpr_nl, tpr_nl, color="#FF9800", lw=2, ls="--",
                label=f"QS w/o LOH penalty (AUC={auc_no_loh:.3f})")

        # Diagonal
        ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5)
        ax.set_xlabel("False Positive Rate", fontsize=10)
        ax.set_ylabel("True Positive Rate" if idx == 0 else "", fontsize=10)
        ax.set_title(f"{mode.upper()} Mode", fontsize=12, fontweight="bold")
        ax.legend(loc="lower right", fontsize=9)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)

        delta_auc = auc_no_loh - auc_orig
        ax.text(0.05, 0.85, f"$\\Delta$AUC = {delta_auc:+.3f}",
                transform=ax.transAxes, fontsize=10, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

        stats[f"auc_{mode}_full_qs"] = auc_orig
        stats[f"auc_{mode}_no_loh"] = auc_no_loh
        stats[f"auc_{mode}_delta_no_loh"] = delta_auc

    fig.suptitle("Fig04: ROC Comparison — QS With vs Without LOH Penalty", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0.05, 1, 0.93])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | LOH penalty = -25 removed | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig04_qs_without_loh_penalty_roc.png")

    return stats


# ---------------------------------------------------------------------------
# Fig05: Per-Component AUC Ranking
# ---------------------------------------------------------------------------

def fig05_component_auc_ranking(df: pd.DataFrame) -> Dict[str, Any]:
    """AUC of each individual QS component for TP/FP discrimination."""
    stats = {}

    results = []
    for mode in MODE_ORDER:
        sub = df[df["mode"] == mode]
        y_true = sub["truth_binary"].values
        for comp in ALL_COMPONENT_COLS:
            scores = sub[comp].values
            # For penalties (negative values), higher penalty -> lower score -> should predict FP
            # Use raw value: higher value = more TP-like
            auc = compute_auc(y_true, scores)
            results.append({
                "mode": mode,
                "component": comp,
                "short_name": COMPONENT_SHORT[comp],
                "auc": auc,
            })
            stats[f"auc_{mode}_{comp}"] = auc

        # Also add the full QS and QS-without-LOH
        auc_full = compute_auc(y_true, sub["qs_reconstructed"].values)
        results.append({"mode": mode, "component": "qs_reconstructed",
                        "short_name": "Full QS", "auc": auc_full})
        stats[f"auc_{mode}_full_qs"] = auc_full

    results_df = pd.DataFrame(results)

    fig, ax = plt.subplots(figsize=(14, 7))

    # Grouped bar chart
    components_order = [COMPONENT_SHORT[c] for c in ALL_COMPONENT_COLS] + ["Full QS"]
    x = np.arange(len(components_order))
    width = 0.35

    paired_aucs = []
    to_aucs = []
    for name in components_order:
        p_row = results_df[(results_df["mode"] == "paired") & (results_df["short_name"] == name)]
        t_row = results_df[(results_df["mode"] == "to") & (results_df["short_name"] == name)]
        paired_aucs.append(p_row["auc"].values[0] if len(p_row) else 0.5)
        to_aucs.append(t_row["auc"].values[0] if len(t_row) else 0.5)

    bars_p = ax.bar(x - width / 2, paired_aucs, width, label="Paired", color=COLOR_PAIRED,
                    edgecolor="white", linewidth=0.8)
    bars_t = ax.bar(x + width / 2, to_aucs, width, label="TO", color=COLOR_TO,
                    edgecolor="white", linewidth=0.8)

    # Annotate
    for bar, val in zip(bars_p, paired_aucs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f"{val:.3f}", ha="center", va="bottom", fontsize=7, fontweight="bold",
                color=COLOR_PAIRED)
    for bar, val in zip(bars_t, to_aucs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f"{val:.3f}", ha="center", va="bottom", fontsize=7, fontweight="bold",
                color=COLOR_TO)

    ax.axhline(0.5, color="red", ls="--", lw=1, alpha=0.6, label="Random (AUC=0.5)")
    ax.set_xticks(x)
    ax.set_xticklabels(components_order, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("AUC (TP vs FP)", fontsize=11)
    ax.set_title("Fig05: Individual Component AUC — Paired vs TO", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.set_ylim(0.3, 0.85)

    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | Each component's raw value as standalone predictor | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig05_qs_component_auc_ranking.png")

    # Save data
    results_df.to_csv(DATA_DIR / "component_auc_by_mode.tsv", sep="\t", index=False)

    return stats


# ---------------------------------------------------------------------------
# Fig06: QS Distribution Violin by Truth x Mode
# ---------------------------------------------------------------------------

def fig06_qs_violin(df: pd.DataFrame) -> Dict[str, Any]:
    """Violin plot of QS distribution across TP/FP x paired/TO."""
    stats = {}

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)

    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub = df[df["mode"] == mode].copy()

        parts = ax.violinplot(
            [sub[sub["truth_label"] == t]["qs_reconstructed"].dropna().values for t in ["TP", "FP"]],
            positions=[0, 1],
            showmeans=True,
            showmedians=True,
            showextrema=False,
        )

        # Color the violins
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(COLOR_TP if i == 0 else COLOR_FP)
            pc.set_alpha(0.6)
        parts["cmeans"].set_color("black")
        parts["cmedians"].set_color("white")

        # Add box plot overlay
        bp = ax.boxplot(
            [sub[sub["truth_label"] == t]["qs_reconstructed"].dropna().values for t in ["TP", "FP"]],
            positions=[0, 1],
            widths=0.15,
            patch_artist=False,
            showfliers=False,
            medianprops=dict(color="black", lw=1.5),
            boxprops=dict(color="black"),
            whiskerprops=dict(color="black"),
            capprops=dict(color="black"),
        )

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["TP", "FP"], fontsize=11)
        mode_label = mode.upper()

        # Stats annotations
        for j, truth in enumerate(["TP", "FP"]):
            vals = sub[sub["truth_label"] == truth]["qs_reconstructed"].dropna()
            ax.text(j, 105, f"n={len(vals):,}\nmed={vals.median():.1f}\nmean={vals.mean():.1f}",
                    ha="center", va="bottom", fontsize=8, fontweight="bold")
            stats[f"qs_{mode}_{truth}_mean"] = float(vals.mean())
            stats[f"qs_{mode}_{truth}_median"] = float(vals.median())
            stats[f"qs_{mode}_{truth}_std"] = float(vals.std())

        # Mann-Whitney
        tp_vals = sub[sub["truth_label"] == "TP"]["qs_reconstructed"].dropna().values
        fp_vals = sub[sub["truth_label"] == "FP"]["qs_reconstructed"].dropna().values
        mw = mannwhitney_test(tp_vals, fp_vals)
        es = compute_effect_size(tp_vals, fp_vals)

        ax.text(0.5, 0.02, f"MWU p={format_p(mw['p_value'])}\nCohen's d={es['d']:.3f} {format_ci(es['ci_lo'], es['ci_hi'])}",
                transform=ax.transAxes, ha="center", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

        ax.set_title(f"{mode_label} Mode", fontsize=12, fontweight="bold")
        ax.set_ylabel("Quality Score (Reconstructed)" if idx == 0 else "")
        ax.set_ylim(-5, 125)

        stats[f"qs_{mode}_mw_p"] = mw["p_value"]
        stats[f"qs_{mode}_cohens_d"] = es["d"]
        stats[f"qs_{mode}_cohens_d_ci_lo"] = es["ci_lo"]
        stats[f"qs_{mode}_cohens_d_ci_hi"] = es["ci_hi"]

    fig.suptitle("Fig06: QS Distribution by Truth Label x Mode", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0.05, 1, 0.93])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | Violin + box plot | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig06_qs_final_distribution_by_truth_mode.png")

    return stats


# ---------------------------------------------------------------------------
# Fig07: QS Tier Confusion Matrix
# ---------------------------------------------------------------------------

def fig07_tier_confusion(df: pd.DataFrame) -> Dict[str, Any]:
    """Confusion matrix of QS Tier vs truth label, per mode."""
    stats = {}

    # Define tier thresholds (from InterSubMod code)
    tier_order = ["High", "Medium", "Low"]

    def assign_tier(qs: float) -> str:
        if qs >= 70:
            return "High"
        elif qs >= 40:
            return "Medium"
        else:
            return "Low"

    df["qs_tier_recon"] = df["qs_reconstructed"].apply(assign_tier)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, mode in enumerate(MODE_ORDER):
        ax = axes[idx]
        sub = df[df["mode"] == mode]

        # Build confusion matrix: rows = tier, cols = truth
        cm_data = []
        for tier in tier_order:
            tier_sub = sub[sub["qs_tier_recon"] == tier]
            tp_count = (tier_sub["truth_label"] == "TP").sum()
            fp_count = (tier_sub["truth_label"] == "FP").sum()
            total = tp_count + fp_count
            tp_rate = safe_div(tp_count, total) * 100
            fp_rate = safe_div(fp_count, total) * 100
            cm_data.append({
                "tier": tier, "tp": tp_count, "fp": fp_count, "total": total,
                "tp_rate": tp_rate, "fp_rate": fp_rate,
            })
            stats[f"tier_{mode}_{tier}_tp"] = tp_count
            stats[f"tier_{mode}_{tier}_fp"] = fp_count
            stats[f"tier_{mode}_{tier}_tp_rate"] = tp_rate

        # Plot as heatmap
        cm = np.array([[r["tp"], r["fp"]] for r in cm_data])
        # Normalize by row
        cm_norm = cm / cm.sum(axis=1, keepdims=True) * 100

        im = ax.imshow(cm_norm, cmap="RdYlBu", aspect="auto", vmin=0, vmax=100)

        # Annotate
        for i in range(len(tier_order)):
            for j in range(2):
                text_color = "white" if cm_norm[i, j] > 70 or cm_norm[i, j] < 30 else "black"
                ax.text(j, i, f"{cm_norm[i, j]:.1f}%\n({cm[i, j]:,})",
                        ha="center", va="center", fontsize=9, fontweight="bold",
                        color=text_color)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["TP", "FP"], fontsize=10)
        ax.set_yticks(range(len(tier_order)))
        ax.set_yticklabels([f"{t} (QS {'>=70' if t=='High' else '40-69' if t=='Medium' else '<40'})" for t in tier_order], fontsize=9)
        ax.set_xlabel("Truth Label", fontsize=10)
        ax.set_ylabel("QS Tier" if idx == 0 else "")
        ax.set_title(f"{mode.upper()} Mode", fontsize=12, fontweight="bold")

        # Chi-squared
        chi2_res = chi2_test(cm)
        stats[f"tier_chi2_{mode}_chi2"] = chi2_res["chi2"]
        stats[f"tier_chi2_{mode}_p"] = chi2_res["p_value"]
        stats[f"tier_chi2_{mode}_cramers_v"] = chi2_res["cramers_v"]

        ax.text(0.5, -0.18, f"Chi2={chi2_res['chi2']:.1f}, p={format_p(chi2_res['p_value'])}, V={chi2_res['cramers_v']:.3f}",
                transform=ax.transAxes, ha="center", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    fig.suptitle("Fig07: QS Tier vs Truth Label Confusion Matrix", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0.08, 1, 0.93])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | Tier: High(>=70), Medium(40-69), Low(<40) | Row-normalized % | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig07_qs_tier_confusion_matrix.png")

    return stats


# ---------------------------------------------------------------------------
# Fig08: Leave-One-Out AUC Sensitivity
# ---------------------------------------------------------------------------

def fig08_loo_sensitivity(df: pd.DataFrame) -> Dict[str, Any]:
    """AUC change when each component is individually removed."""
    stats = {}

    results = []
    for mode in MODE_ORDER:
        sub = df[df["mode"] == mode]
        y_true = sub["truth_binary"].values

        # Baseline: full QS
        auc_baseline = compute_auc(y_true, sub["qs_reconstructed"].values)
        stats[f"loo_{mode}_baseline_auc"] = auc_baseline

        for comp in ALL_COMPONENT_COLS:
            # Remove this component
            qs_without = (sub["qs_reconstructed"] - sub[comp]).clip(0, 100).values
            auc_without = compute_auc(y_true, qs_without)
            delta = auc_without - auc_baseline
            results.append({
                "mode": mode,
                "removed": comp,
                "short_name": COMPONENT_SHORT[comp],
                "auc_baseline": auc_baseline,
                "auc_without": auc_without,
                "delta_auc": delta,
            })
            stats[f"loo_{mode}_{comp}_delta"] = delta
            stats[f"loo_{mode}_{comp}_auc_without"] = auc_without

    results_df = pd.DataFrame(results)

    fig, ax = plt.subplots(figsize=(14, 7))

    component_names = [COMPONENT_SHORT[c] for c in ALL_COMPONENT_COLS]
    x = np.arange(len(component_names))
    width = 0.35

    paired_deltas = []
    to_deltas = []
    for name in component_names:
        p_row = results_df[(results_df["mode"] == "paired") & (results_df["short_name"] == name)]
        t_row = results_df[(results_df["mode"] == "to") & (results_df["short_name"] == name)]
        paired_deltas.append(p_row["delta_auc"].values[0] if len(p_row) else 0)
        to_deltas.append(t_row["delta_auc"].values[0] if len(t_row) else 0)

    bars_p = ax.bar(x - width / 2, paired_deltas, width, label="Paired", color=COLOR_PAIRED,
                    edgecolor="white", linewidth=0.8)
    bars_t = ax.bar(x + width / 2, to_deltas, width, label="TO", color=COLOR_TO,
                    edgecolor="white", linewidth=0.8)

    # Annotate
    for bar, val in zip(bars_p, paired_deltas):
        y_offset = 0.002 if val >= 0 else -0.002
        va = "bottom" if val >= 0 else "top"
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + y_offset,
                f"{val:+.4f}", ha="center", va=va, fontsize=7, fontweight="bold",
                color=COLOR_PAIRED)
    for bar, val in zip(bars_t, to_deltas):
        y_offset = 0.002 if val >= 0 else -0.002
        va = "bottom" if val >= 0 else "top"
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + y_offset,
                f"{val:+.4f}", ha="center", va=va, fontsize=7, fontweight="bold",
                color=COLOR_TO)

    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(component_names, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("$\\Delta$AUC (After Removing Component)", fontsize=11)
    ax.set_title("Fig08: QS AUC Sensitivity — Leave-One-Component-Out Analysis", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)

    # Add baseline AUC annotation
    paired_base = stats["loo_paired_baseline_auc"]
    to_base = stats["loo_to_baseline_auc"]
    ax.text(0.02, 0.95, f"Baseline AUC: Paired={paired_base:.3f}, TO={to_base:.3f}",
            transform=ax.transAxes, fontsize=9, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    add_caption(fig, f"Source: all_region_rows.tsv.gz | Positive delta = removing component improves AUC | {OBSERVATION_ID}")
    save_figure(fig, FIG_DIR / "fig08_qs_sensitivity_to_each_component.png")

    # Save data
    results_df.to_csv(DATA_DIR / "loo_auc_sensitivity.tsv", sep="\t", index=False)

    return stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ensure_dir(OUT_DIR)
    ensure_dir(DATA_DIR)
    setup_plot_style()

    # Load data
    df = load_and_prepare()

    # Collect all stats
    all_stats = {}

    # Generate figures
    print("\n[O02] Generating Fig01: Waterfall Paired ...")
    all_stats.update(fig01_waterfall_paired(df))

    print("[O02] Generating Fig02: Waterfall TO ...")
    all_stats.update(fig02_waterfall_to(df))

    print("[O02] Generating Fig03: LOH Trigger Rate ...")
    all_stats.update(fig03_loh_trigger_rate(df))

    print("[O02] Generating Fig04: ROC Without LOH ...")
    all_stats.update(fig04_roc_without_loh(df))

    print("[O02] Generating Fig05: Component AUC Ranking ...")
    all_stats.update(fig05_component_auc_ranking(df))

    print("[O02] Generating Fig06: QS Violin ...")
    all_stats.update(fig06_qs_violin(df))

    print("[O02] Generating Fig07: Tier Confusion ...")
    all_stats.update(fig07_tier_confusion(df))

    print("[O02] Generating Fig08: LOO Sensitivity ...")
    all_stats.update(fig08_loo_sensitivity(df))

    # Write feature statistics
    qs_cols = ALL_COMPONENT_COLS + ["qs_reconstructed"]
    if "Quality_Score" in df.columns:
        qs_cols.append("Quality_Score")
    write_feature_statistics(df, qs_cols, DATA_DIR / "qs_component_statistics.tsv")

    # Write source summary
    source_summary = []
    for mode in MODE_ORDER:
        for truth in ["TP", "FP"]:
            sub = df[(df["mode"] == mode) & (df["truth_label"] == truth)]
            source_summary.append({
                "mode": mode, "truth": truth, "n": len(sub),
                "qs_mean": sub["qs_reconstructed"].mean(),
                "qs_median": sub["qs_reconstructed"].median(),
                "qs_std": sub["qs_reconstructed"].std(),
                "loh_trigger_rate": (sub["qs_penalty_loh"] < 0).mean() * 100,
            })
    pd.DataFrame(source_summary).to_csv(DATA_DIR / "source_summary.tsv", sep="\t", index=False)

    # Save all stats
    with open(DATA_DIR / "all_stats.json", "w") as f:
        json.dump({k: float(v) if isinstance(v, (np.floating, float)) else v
                    for k, v in all_stats.items()}, f, indent=2, default=str)

    # Write round context
    write_round_context(
        out_dir=OUT_DIR,
        observation_id=OBSERVATION_ID,
        title=OBSERVATION_TITLE,
        description="Decomposition of Quality Score into penalty/bonus components; "
                    "analysis of why QS fails in TO mode (AUC~0.5); "
                    "leave-one-out sensitivity analysis of each component.",
        script_path=SCRIPT_PATH,
        data_sources=[{
            "name": "all_region_rows.tsv.gz",
            "path": str(MASTER_DATASET_PATH),
            "description": "748K-row master dataset with QS, LOH, and component inputs",
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "qs_components": ALL_COMPONENT_COLS,
            "key_finding_paired_auc": all_stats.get("auc_paired_full_qs"),
            "key_finding_to_auc": all_stats.get("auc_to_full_qs"),
            "key_finding_to_no_loh_auc": all_stats.get("auc_to_no_loh"),
        },
    )

    print(f"\n{'='*70}")
    print(f"  O02 Complete! Output: {OUT_DIR}")
    print(f"{'='*70}\n")

    return all_stats


if __name__ == "__main__":
    all_stats = main()
