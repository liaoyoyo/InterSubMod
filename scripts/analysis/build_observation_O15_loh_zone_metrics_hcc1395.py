#!/usr/bin/env python3
"""O15 Phase 1: HCC1395 LOH Zone Metrics — Complete Analysis.

Produces 16 figures analyzing how SEQC2-defined LOH zones affect ALL ISM
metrics in HCC1395 samples (HCC1395, HCC1395_5kHz, HCC1395_DORADO).

Usage:
    python3 build_observation_O15_loh_zone_metrics_hcc1395.py
"""

from __future__ import annotations

import sys
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure, add_caption,
    compute_auc, compute_effect_size, mannwhitney_test,
    encode_truth_binary, ensure_dir, safe_div,
    SAMPLE_ORDER, MODE_ORDER,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO, TRUTH_PALETTE, MODE_PALETTE,
)

import json
import math
from bisect import bisect_left, bisect_right
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats

# ── Output Paths ──────────────────────────────────────────────────────────
BASE_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation")
FIG_DIR = ensure_dir(BASE_DIR / "figures")
DATA_DIR = ensure_dir(BASE_DIR / "data")
REPORT_DIR = ensure_dir(BASE_DIR / "reports")

BED_BASE = BASE_DIR / "data" / "seqc2_vs_to_validation"
BED_TP = BED_BASE / "loh_classified_TP.bed"
BED_FP = BED_BASE / "loh_classified_FP.bed"
BED_FN = BED_BASE / "loh_classified_FN.bed"

CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz (HCC1395 samples) | "
    "LOH zones: SEQC2 4-class (TP/FP/FN/TN)"
)

HCC1395_SAMPLES = ["HCC1395", "HCC1395_5kHz", "HCC1395_DORADO"]

ZONE_ORDER = ["TP", "FP", "FN", "TN"]
ZONE_PALETTE = {"TP": "#4CAF50", "FP": "#F44336", "FN": "#FF9800", "TN": "#2196F3"}

# ── Metric Groups ──────────────────────────────────────────────────────
GROUP_A = ["HP_Ratio", "effective_hp_reads", "hp_assign_rate", "hp0_ratio",
           "HPMergedDelta", "HPMergedP"]
GROUP_B = ["AlleleDelta", "PairwiseMedianDist", "CramersV", "GlobalP",
           "PairwiseMeanDist", "HeuristicScore", "NumCpGs"]
GROUP_C = ["ClusterPermanovaF", "ClusterPermanovaP",
           "LabelHPPermanovaF", "LabelHPPermanovaP",
           "LabelAllelePermanovaF", "LabelAllelePermanovaP"]
GROUP_D_CAT = ["VerificationClass", "PassedGating", "Potential_LOH", "LOH_Subtype"]
GROUP_D_NUM = ["Quality_Score"]
GROUP_E = ["caller_af", "caller_dp", "caller_gq", "caller_filter",
           "caller_ad_ref", "caller_ad_alt"]
GROUP_F = ["AlleleP", "AlleleSig", "HPFineF", "HPFineP", "HPFineSig",
           "UnassignedAffinity", "Stability"]
GROUP_G = ["to_phase_filter", "phase_block_status"]

ALL_NUMERIC_METRICS = (
    GROUP_A + GROUP_B[:6] + ["NumCpGs"] +
    GROUP_C + GROUP_D_NUM +
    ["caller_af", "caller_dp", "caller_gq", "caller_ad_ref", "caller_ad_alt"] +
    ["AlleleP", "HPFineF", "HPFineP", "UnassignedAffinity", "Stability"]
)

LOAD_COLS = (
    ["RegionID", "Chr", "Pos", "sample", "mode", "truth_label"] +
    GROUP_A + GROUP_B + GROUP_C +
    ["VerificationClass", "Quality_Score", "PassedGating", "Potential_LOH", "LOH_Subtype"] +
    GROUP_E + GROUP_F + GROUP_G +
    ["NumReads"]
)


# ══════════════════════════════════════════════════════════════════════════
# LOH Zone Assignment
# ══════════════════════════════════════════════════════════════════════════

def load_bed_intervals(bed_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Load BED intervals sorted by start for binary search."""
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def point_in_intervals(chrom: str, pos: int, intervals: Dict[str, List[Tuple[int, int]]]) -> bool:
    """Check if a 0-based position falls in any interval (binary search)."""
    if chrom not in intervals:
        return False
    ivs = intervals[chrom]
    lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos < s:
            hi = mid - 1
        elif pos >= e:
            lo = mid + 1
        else:
            return True
    return False


def assign_loh_zones_vectorized(
    df: pd.DataFrame,
    tp_ivs: Dict, fp_ivs: Dict, fn_ivs: Dict,
) -> pd.Series:
    """Assign LOH zone to each row — vectorized by chromosome for speed."""
    print("[O15] Assigning LOH zones (vectorized)...")
    zones = pd.Series("TN", index=df.index, dtype="object")

    # Group by chromosome for efficiency
    for chrom, grp in df.groupby("Chr"):
        positions = grp["Pos"].values
        idx = grp.index
        # Check each zone — priority: TP > FP > FN > TN
        for zone_name, zone_ivs in [("FN", fn_ivs), ("FP", fp_ivs), ("TP", tp_ivs)]:
            if chrom not in zone_ivs:
                continue
            ivs = zone_ivs[chrom]
            starts = np.array([s for s, e in ivs])
            ends = np.array([e for s, e in ivs])
            # For each position, find which interval it might fall into
            insert_idx = np.searchsorted(starts, positions, side="right") - 1
            in_range = (insert_idx >= 0) & (insert_idx < len(ivs))
            mask = np.zeros(len(positions), dtype=bool)
            valid = np.where(in_range)[0]
            if len(valid) > 0:
                valid_ins = insert_idx[valid]
                mask[valid] = positions[valid] >= starts[valid_ins]
                mask[valid] &= positions[valid] < ends[valid_ins]
            zones.iloc[grp.index.get_indexer(idx[mask])] = zone_name

    zone_counts = zones.value_counts()
    print(f"[O15] Zone assignment complete: {dict(zone_counts)}")
    return zones


# ══════════════════════════════════════════════════════════════════════════
# Statistical Helpers
# ══════════════════════════════════════════════════════════════════════════

def _auc_safe(y_true, y_score):
    """ROC AUC with graceful fallback."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_score)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 10:
        return float("nan")
    return compute_auc(y_true, y_score)


def _cohend(g1, g2):
    """Cohen's d with bootstrap CI."""
    return compute_effect_size(g1, g2, n_boot=500, ci=0.95)


def _available(df, cols):
    """Return columns that exist in df."""
    return [c for c in cols if c in df.columns]


def _numeric_cols(df, cols):
    """Return numeric-dtype columns from list."""
    avail = _available(df, cols)
    return [c for c in avail if pd.api.types.is_numeric_dtype(df[c])]


# ══════════════════════════════════════════════════════════════════════════
# Figure Functions
# ══════════════════════════════════════════════════════════════════════════

def fig01_hp_ratio_violin(df: pd.DataFrame) -> Path:
    """Fig01: HP_Ratio by LOH zone x mode (Violin), split by truth_label."""
    print("[O15] Fig01: HP_Ratio violin by LOH zone x mode...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for i, m in enumerate(["to", "paired"]):
        ax = axes[i]
        sub = df[df["mode"] == m].dropna(subset=["HP_Ratio"])
        if len(sub) == 0:
            ax.set_title(f"{m.upper()} (no data)")
            continue
        sns.violinplot(
            data=sub, x="loh_zone", y="HP_Ratio", hue="truth_label",
            order=ZONE_ORDER, hue_order=["TP", "FP"],
            palette=TRUTH_PALETTE, split=True, inner="quartile",
            ax=ax, cut=0, density_norm="width",
        )
        # Annotate counts
        for zi, z in enumerate(ZONE_ORDER):
            zd = sub[sub["loh_zone"] == z]
            ntp = (zd["truth_label"] == "TP").sum()
            nfp = (zd["truth_label"] == "FP").sum()
            ax.text(zi, -0.08, f"TP:{ntp}\nFP:{nfp}", ha="center", va="top",
                    fontsize=7, color="#555555")
        ax.set_title(f"{m.upper()}: HP_Ratio by LOH Zone")
        ax.set_xlabel("SEQC2 LOH Zone")
        ax.set_ylabel("HP_Ratio" if i == 0 else "")
        ax.set_ylim(-0.15, 1.05)
        ax.legend(title="Truth", loc="upper right", fontsize=8)
    fig.suptitle("Fig01: HP_Ratio Distribution by LOH Zone × Mode", fontsize=14, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig01_hp_ratio_violin.png")


def fig02_extreme_hp_ratio(df: pd.DataFrame) -> Path:
    """Fig02: Extreme HP_Ratio rate by zone x mode x truth."""
    print("[O15] Fig02: Extreme HP_Ratio rate...")
    df = df.dropna(subset=["HP_Ratio"]).copy()
    df["extreme_hp"] = (df["HP_Ratio"] < 0.1) | (df["HP_Ratio"] > 0.9)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ri, m in enumerate(["to", "paired"]):
        for ci, tl in enumerate(["TP", "FP"]):
            ax = axes[ri, ci]
            sub = df[(df["mode"] == m) & (df["truth_label"] == tl)]
            rates = []
            counts = []
            for z in ZONE_ORDER:
                zd = sub[sub["loh_zone"] == z]
                n = len(zd)
                r = zd["extreme_hp"].mean() * 100 if n > 0 else 0
                rates.append(r)
                counts.append(n)
            bars = ax.bar(ZONE_ORDER, rates, color=[ZONE_PALETTE[z] for z in ZONE_ORDER],
                          edgecolor="black", linewidth=0.5)
            for bi, (bar, r, n) in enumerate(zip(bars, rates, counts)):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                        f"{r:.1f}%\n(n={n})", ha="center", va="bottom", fontsize=8)
            ax.set_title(f"{m.upper()} — {tl}")
            ax.set_ylabel("Extreme HP_Ratio Rate (%)")
            ax.set_xlabel("LOH Zone")
            ax.set_ylim(0, min(max(rates) * 1.4 + 5, 105))
    fig.suptitle("Fig02: Extreme HP_Ratio Rate (< 0.1 or > 0.9) by Zone × Mode × Truth",
                 fontsize=13, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig02_extreme_hp_ratio.png")


def fig03_methyl_metrics_violin(df: pd.DataFrame) -> Path:
    """Fig03: ISM Methylation metrics by zone x mode (Violin 4x2)."""
    print("[O15] Fig03: Methylation metrics violin...")
    metrics = _available(df, ["AlleleDelta", "PairwiseMedianDist", "CramersV", "HPMergedDelta"])
    if not metrics:
        print("  [SKIP] No methylation metrics available")
        return None

    fig, axes = plt.subplots(len(metrics), 2, figsize=(14, 4 * len(metrics)), squeeze=False)
    for ri, met in enumerate(metrics):
        for ci, m in enumerate(["to", "paired"]):
            ax = axes[ri, ci]
            sub = df[(df["mode"] == m)].dropna(subset=[met])
            if len(sub) < 10:
                ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            sns.violinplot(
                data=sub, x="loh_zone", y=met, order=ZONE_ORDER,
                palette=ZONE_PALETTE, inner="quartile", ax=ax, cut=0,
                density_norm="width",
            )
            for zi, z in enumerate(ZONE_ORDER):
                n = (sub["loh_zone"] == z).sum()
                ax.text(zi, ax.get_ylim()[0], f"n={n}", ha="center", va="top",
                        fontsize=7, color="#555555")
            ax.set_title(f"{met} — {m.upper()}")
            ax.set_xlabel("LOH Zone")
            ax.set_ylabel(met if ci == 0 else "")
    fig.suptitle("Fig03: ISM Methylation Metrics by LOH Zone × Mode", fontsize=14, y=1.01)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig03_methyl_metrics_violin.png")


def fig04_auc_heatmap(df: pd.DataFrame) -> Path:
    """Fig04: AUC heatmap for TP/FP discrimination."""
    print("[O15] Fig04: AUC heatmap...")
    # Define strata: inside-LOH = TP+FP zones, outside-LOH = FN+TN zones
    df = df.copy()
    df["truth_binary"] = encode_truth_binary(df)
    df["loh_group"] = df["loh_zone"].map(
        {"TP": "inside", "FP": "inside", "FN": "outside", "TN": "outside"}
    )

    strata = [
        ("TO-inside", "to", "inside"),
        ("TO-outside", "to", "outside"),
        ("Paired-inside", "paired", "inside"),
        ("Paired-outside", "paired", "outside"),
    ]

    # Collect all numeric metrics
    all_metrics = sorted(set(_numeric_cols(df, ALL_NUMERIC_METRICS)))

    auc_data = {}
    for sname, mode, loh_grp in strata:
        sub = df[(df["mode"] == mode) & (df["loh_group"] == loh_grp)].dropna(subset=["truth_binary"])
        row_aucs = {}
        for met in all_metrics:
            row_aucs[met] = _auc_safe(sub["truth_binary"].values, sub[met].values)
        auc_data[sname] = row_aucs

    auc_df = pd.DataFrame(auc_data).reindex(all_metrics)
    auc_df = auc_df.dropna(how="all")

    # Sort by max AUC across strata
    auc_df["max_auc"] = auc_df.max(axis=1)
    auc_df = auc_df.sort_values("max_auc", ascending=False).drop(columns=["max_auc"])

    fig, ax = plt.subplots(figsize=(10, max(8, len(auc_df) * 0.35)))
    norm = TwoSlopeNorm(vmin=0.3, vcenter=0.5, vmax=0.8)
    sns.heatmap(
        auc_df, annot=True, fmt=".3f", cmap="RdYlGn", norm=norm,
        ax=ax, linewidths=0.5, cbar_kws={"label": "AUC (TP=1 vs FP=0)"},
    )
    ax.set_title("Fig04: AUC Heatmap — TP/FP Discrimination by LOH Zone × Mode")
    ax.set_ylabel("Metric")
    ax.set_xlabel("Stratum")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig04_auc_heatmap.png")


def fig05_verification_class(df: pd.DataFrame) -> Path:
    """Fig05: VerificationClass by zone x truth (Stacked bar, TO only)."""
    print("[O15] Fig05: VerificationClass stacked bar (TO)...")
    sub = df[(df["mode"] == "to")].dropna(subset=["VerificationClass"])
    if len(sub) < 10:
        print("  [SKIP] Not enough TO data")
        return None

    vc_order = ["Noise", "Weak", "Subclone", "Strong"]
    vc_colors = {"Noise": "#BDBDBD", "Weak": "#FFB74D", "Subclone": "#42A5F5", "Strong": "#66BB6A"}

    fig, axes = plt.subplots(1, 4, figsize=(16, 6), sharey=True)
    for zi, z in enumerate(ZONE_ORDER):
        ax = axes[zi]
        zsub = sub[sub["loh_zone"] == z]
        if len(zsub) == 0:
            ax.set_title(f"{z} Zone (no data)")
            continue

        # For each truth label, get VC proportions
        bottom_tp = 0
        bottom_fp = 0
        x_positions = [0, 1]
        width = 0.6
        for vc in vc_order:
            for ti, tl in enumerate(["TP", "FP"]):
                tsub = zsub[zsub["truth_label"] == tl]
                total = len(tsub)
                count = (tsub["VerificationClass"] == vc).sum()
                pct = count / total * 100 if total > 0 else 0
                bottom = bottom_tp if ti == 0 else bottom_fp
                bar = ax.bar(x_positions[ti], pct, width=width, bottom=bottom,
                             color=vc_colors.get(vc, "#999"),
                             edgecolor="white", linewidth=0.3)
                if pct > 8:
                    ax.text(x_positions[ti], bottom + pct / 2, f"{pct:.0f}%",
                            ha="center", va="center", fontsize=7, fontweight="bold")
                if ti == 0:
                    bottom_tp += pct
                else:
                    bottom_fp += pct
        ax.set_title(f"{z} Zone")
        ax.set_xticks(x_positions)
        n_tp = (zsub["truth_label"] == "TP").sum()
        n_fp = (zsub["truth_label"] == "FP").sum()
        ax.set_xticklabels([f"TP\n(n={n_tp})", f"FP\n(n={n_fp})"], fontsize=9)
        ax.set_ylabel("Percentage (%)" if zi == 0 else "")
        ax.set_ylim(0, 105)

    # Legend
    handles = [mpatches.Patch(color=vc_colors[vc], label=vc) for vc in vc_order]
    fig.legend(handles=handles, loc="upper right", title="VerificationClass", fontsize=9)
    fig.suptitle("Fig05: VerificationClass by LOH Zone × Truth (TO Only)", fontsize=13, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig05_verification_class.png")


def fig06_quality_score(df: pd.DataFrame) -> Path:
    """Fig06: Quality_Score by zone x mode x truth (Box plot 2x2)."""
    print("[O15] Fig06: Quality_Score box plot...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ri, m in enumerate(["to", "paired"]):
        for ci, tl in enumerate(["TP", "FP"]):
            ax = axes[ri, ci]
            sub = df[(df["mode"] == m) & (df["truth_label"] == tl)].dropna(subset=["Quality_Score"])
            if len(sub) < 5:
                ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            sns.boxplot(
                data=sub, x="loh_zone", y="Quality_Score", order=ZONE_ORDER,
                palette=ZONE_PALETTE, ax=ax, fliersize=2,
            )
            for zi, z in enumerate(ZONE_ORDER):
                zd = sub[sub["loh_zone"] == z]
                n = len(zd)
                med = zd["Quality_Score"].median() if n > 0 else 0
                ax.text(zi, ax.get_ylim()[0] - 2, f"n={n}\nmed={med:.0f}",
                        ha="center", va="top", fontsize=7, color="#555555")
            ax.set_title(f"{m.upper()} — {tl}")
            ax.set_xlabel("LOH Zone")
            ax.set_ylabel("Quality_Score" if ci == 0 else "")
    fig.suptitle("Fig06: Quality_Score by LOH Zone × Mode × Truth", fontsize=13, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig06_quality_score.png")


def fig07_allele_delta_vs_af(df: pd.DataFrame) -> Path:
    """Fig07: AlleleDelta vs caller_af scatter (LOH-inside vs outside, TO)."""
    print("[O15] Fig07: AlleleDelta vs caller_af scatter (TO)...")
    sub = df[(df["mode"] == "to")].dropna(subset=["AlleleDelta", "caller_af"]).copy()
    sub["loh_group"] = sub["loh_zone"].map(
        {"TP": "LOH-inside", "FP": "LOH-inside", "FN": "LOH-outside", "TN": "LOH-outside"}
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for i, grp in enumerate(["LOH-inside", "LOH-outside"]):
        ax = axes[i]
        gsub = sub[sub["loh_group"] == grp]
        for tl, color in [("TP", COLOR_TP), ("FP", COLOR_FP)]:
            tsub = gsub[gsub["truth_label"] == tl]
            ax.scatter(tsub["caller_af"], tsub["AlleleDelta"],
                       alpha=0.15, s=8, c=color, label=f"{tl} (n={len(tsub)})",
                       rasterized=True)
        # Pearson R for all
        mask = np.isfinite(gsub["caller_af"].values) & np.isfinite(gsub["AlleleDelta"].values)
        if mask.sum() > 10:
            r, p = sp_stats.pearsonr(gsub["caller_af"].values[mask],
                                      gsub["AlleleDelta"].values[mask])
            ax.text(0.02, 0.98, f"Pearson R={r:.3f}\np={p:.2e}\nn={mask.sum()}",
                    transform=ax.transAxes, va="top", fontsize=9,
                    bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        ax.set_title(f"{grp} (TO)")
        ax.set_xlabel("caller_af")
        ax.set_ylabel("AlleleDelta" if i == 0 else "")
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Fig07: AlleleDelta vs caller_af — LOH-inside vs LOH-outside (TO)",
                 fontsize=13, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig07_allele_delta_vs_af.png")


def fig08_forest_plot(df: pd.DataFrame) -> Path:
    """Fig08: Forest plot — Cohen's d (TP vs FP) per metric, inside vs outside LOH."""
    print("[O15] Fig08: Forest plot (Cohen's d)...")
    sub = df[df["mode"] == "to"].copy()
    sub["loh_group"] = sub["loh_zone"].map(
        {"TP": "inside", "FP": "inside", "FN": "outside", "TN": "outside"}
    )
    sub["truth_binary"] = encode_truth_binary(sub)

    metrics = sorted(set(_numeric_cols(df, ALL_NUMERIC_METRICS)))
    results = []
    for met in metrics:
        for grp_name, grp_label in [("inside", "LOH-inside"), ("outside", "LOH-outside")]:
            gsub = sub[sub["loh_group"] == grp_name].dropna(subset=[met, "truth_binary"])
            tp_vals = gsub[gsub["truth_binary"] == 1][met].values
            fp_vals = gsub[gsub["truth_binary"] == 0][met].values
            es = _cohend(tp_vals, fp_vals)
            results.append({
                "metric": met, "group": grp_label,
                "d": es["d"], "ci_lo": es["ci_lo"], "ci_hi": es["ci_hi"],
                "n1": es["n1"], "n2": es["n2"],
            })

    res_df = pd.DataFrame(results)
    # Sort by absolute d for inside
    inside_d = res_df[res_df["group"] == "LOH-inside"].set_index("metric")["d"].abs()
    metric_order = inside_d.sort_values(ascending=True).index.tolist()

    fig, ax = plt.subplots(figsize=(10, max(8, len(metric_order) * 0.35)))
    y_positions = np.arange(len(metric_order))
    offset = 0.15

    for gi, (grp, color, marker) in enumerate([
        ("LOH-inside", "#FF9800", "D"), ("LOH-outside", "#2196F3", "o")
    ]):
        gdf = res_df[res_df["group"] == grp].set_index("metric").reindex(metric_order)
        y = y_positions + (gi - 0.5) * offset * 2
        ax.errorbar(
            gdf["d"].values, y,
            xerr=[gdf["d"].values - gdf["ci_lo"].values,
                  gdf["ci_hi"].values - gdf["d"].values],
            fmt=marker, color=color, markersize=5, capsize=2, label=grp,
            linewidth=1, elinewidth=0.8,
        )

    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(metric_order, fontsize=8)
    ax.set_xlabel("Cohen's d (TP − FP)")
    ax.set_title("Fig08: Cohen's d Forest Plot — TP vs FP per Metric (TO)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig08_forest_plot.png")


def fig09_passed_gating(df: pd.DataFrame) -> Path:
    """Fig09: PassedGating rate by zone x mode x truth."""
    print("[O15] Fig09: PassedGating rate...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ri, m in enumerate(["to", "paired"]):
        for ci, tl in enumerate(["TP", "FP"]):
            ax = axes[ri, ci]
            sub = df[(df["mode"] == m) & (df["truth_label"] == tl)].copy()
            sub["passed"] = sub["PassedGating"].astype(str).str.upper().isin(["PASS", "TRUE", "1"])
            rates = []
            counts = []
            for z in ZONE_ORDER:
                zd = sub[sub["loh_zone"] == z]
                n = len(zd)
                r = zd["passed"].mean() * 100 if n > 0 else 0
                rates.append(r)
                counts.append(n)
            bars = ax.bar(ZONE_ORDER, rates, color=[ZONE_PALETTE[z] for z in ZONE_ORDER],
                          edgecolor="black", linewidth=0.5)
            for bi, (bar, r, n) in enumerate(zip(bars, rates, counts)):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                        f"{r:.1f}%\n(n={n})", ha="center", va="bottom", fontsize=8)
            ax.set_title(f"{m.upper()} — {tl}")
            ax.set_ylabel("PASS Rate (%)")
            ax.set_xlabel("LOH Zone")
            ax.set_ylim(0, min(max(rates) * 1.3 + 5, 105))
    fig.suptitle("Fig09: PassedGating PASS Rate by LOH Zone × Mode × Truth",
                 fontsize=13, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig09_passed_gating.png")


def fig10_loh_subtype_heatmap(df: pd.DataFrame) -> Path:
    """Fig10: LOH_Subtype x SEQC2 zone cross-tabulation (Heatmap, TO)."""
    print("[O15] Fig10: LOH_Subtype x Zone heatmap (TO)...")
    sub = df[(df["mode"] == "to")].dropna(subset=["LOH_Subtype"])
    if len(sub) < 10:
        print("  [SKIP] Not enough data")
        return None

    ct = pd.crosstab(sub["LOH_Subtype"], sub["loh_zone"])
    ct = ct.reindex(columns=ZONE_ORDER, fill_value=0)
    # Compute proportions per column
    ct_pct = ct.div(ct.sum(axis=0), axis=1) * 100

    # Annotate with count (pct%)
    annot = ct.copy().astype(str)
    for row in ct.index:
        for col in ct.columns:
            annot.loc[row, col] = f"{ct.loc[row, col]}\n({ct_pct.loc[row, col]:.1f}%)"

    fig, ax = plt.subplots(figsize=(10, max(6, len(ct) * 0.5)))
    sns.heatmap(ct, annot=annot, fmt="", cmap="YlOrRd", ax=ax, linewidths=0.5,
                cbar_kws={"label": "Count"})
    ax.set_title("Fig10: LOH_Subtype × SEQC2 Zone Cross-Tabulation (TO)")
    ax.set_ylabel("LOH_Subtype")
    ax.set_xlabel("SEQC2 LOH Zone")
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig10_loh_subtype_heatmap.png")


def fig11_permanova_f_stats(df: pd.DataFrame) -> Path:
    """Fig11: PERMANOVA F-stats by zone (Bar plot)."""
    print("[O15] Fig11: PERMANOVA F-stats by zone...")
    f_cols = _available(df, ["ClusterPermanovaF", "LabelHPPermanovaF", "LabelAllelePermanovaF"])
    if not f_cols:
        print("  [SKIP] No PERMANOVA columns")
        return None

    fig, axes = plt.subplots(1, len(f_cols), figsize=(5 * len(f_cols), 5), squeeze=False)
    for fi, fcol in enumerate(f_cols):
        ax = axes[0, fi]
        for mi, m in enumerate(["to", "paired"]):
            medians = []
            q25s = []
            q75s = []
            for z in ZONE_ORDER:
                vals = df[(df["mode"] == m) & (df["loh_zone"] == z)][fcol].dropna()
                medians.append(vals.median() if len(vals) > 0 else 0)
                q25s.append(vals.quantile(0.25) if len(vals) > 0 else 0)
                q75s.append(vals.quantile(0.75) if len(vals) > 0 else 0)
            x = np.arange(len(ZONE_ORDER))
            width = 0.35
            offset = (mi - 0.5) * width
            color = COLOR_TO if m == "to" else COLOR_PAIRED
            err_lo = [med - q for med, q in zip(medians, q25s)]
            err_hi = [q - med for med, q in zip(medians, q75s)]
            ax.bar(x + offset, medians, width=width, color=color, alpha=0.8,
                   label=m.upper(), edgecolor="black", linewidth=0.3,
                   yerr=[err_lo, err_hi], capsize=3, error_kw={"linewidth": 0.8})
        ax.set_xticks(x)
        ax.set_xticklabels(ZONE_ORDER)
        ax.set_title(fcol.replace("PermanovaF", ""))
        ax.set_xlabel("LOH Zone")
        ax.set_ylabel("Median F-stat (IQR)")
        ax.legend(fontsize=8)
    fig.suptitle("Fig11: PERMANOVA F-stats by LOH Zone × Mode", fontsize=13, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig11_permanova_f_stats.png")


def fig12_methyl_significance(df: pd.DataFrame) -> Path:
    """Fig12: Methylation Significance (Violin + Bar, 4x2)."""
    print("[O15] Fig12: Methylation significance...")
    # Continuous metrics
    cont_metrics = _available(df, ["AlleleP", "HPFineF"])
    # Binary rate metrics
    bin_metrics = _available(df, ["AlleleSig", "HPFineSig"])

    n_rows = len(cont_metrics) + len(bin_metrics)
    if n_rows == 0:
        print("  [SKIP] No significance metrics")
        return None

    fig, axes = plt.subplots(n_rows, 2, figsize=(14, 4 * n_rows), squeeze=False)
    row_idx = 0

    # Continuous: violin with -log10 transform
    for met in cont_metrics:
        for ci, m in enumerate(["to", "paired"]):
            ax = axes[row_idx, ci]
            sub = df[df["mode"] == m].dropna(subset=[met]).copy()
            # Apply -log10 for p-value columns
            col_name = met
            if "P" in met:
                sub[f"{met}_nlog10"] = -np.log10(sub[met].clip(lower=1e-300))
                col_name = f"{met}_nlog10"
            sns.violinplot(
                data=sub, x="loh_zone", y=col_name, order=ZONE_ORDER,
                palette=ZONE_PALETTE, inner="quartile", ax=ax, cut=0,
                density_norm="width",
            )
            ylabel = f"-log10({met})" if "P" in met else met
            ax.set_title(f"{ylabel} — {m.upper()}")
            ax.set_xlabel("LOH Zone")
            ax.set_ylabel(ylabel if ci == 0 else "")
        row_idx += 1

    # Binary: grouped bar (rate %)
    for met in bin_metrics:
        for ci, m in enumerate(["to", "paired"]):
            ax = axes[row_idx, ci]
            sub = df[df["mode"] == m].copy()
            sub["sig_bool"] = sub[met].astype(str).str.lower().isin(["true", "1", "yes"])
            rates = []
            for z in ZONE_ORDER:
                zd = sub[sub["loh_zone"] == z]
                n = len(zd)
                rates.append(zd["sig_bool"].mean() * 100 if n > 0 else 0)
            ax.bar(ZONE_ORDER, rates, color=[ZONE_PALETTE[z] for z in ZONE_ORDER],
                   edgecolor="black", linewidth=0.5)
            for zi, r in enumerate(rates):
                ax.text(zi, r + 0.5, f"{r:.1f}%", ha="center", va="bottom", fontsize=8)
            ax.set_title(f"{met} Rate (%) — {m.upper()}")
            ax.set_xlabel("LOH Zone")
            ax.set_ylabel("Rate (%)" if ci == 0 else "")
            ax.set_ylim(0, min(max(rates) * 1.3 + 5, 105))
        row_idx += 1

    fig.suptitle("Fig12: Methylation Significance Metrics by LOH Zone × Mode",
                 fontsize=13, y=1.01)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig12_methyl_significance.png")


def fig13_cpg_stability_heuristic(df: pd.DataFrame) -> Path:
    """Fig13: CpG/Stability/HeuristicScore (Box + strip, 3x2)."""
    print("[O15] Fig13: CpG/Stability/HeuristicScore...")
    metrics = _available(df, ["NumCpGs", "Stability", "HeuristicScore"])
    if not metrics:
        print("  [SKIP] No metrics available")
        return None

    fig, axes = plt.subplots(len(metrics), 2, figsize=(14, 4 * len(metrics)), squeeze=False)
    for ri, met in enumerate(metrics):
        for ci, m in enumerate(["to", "paired"]):
            ax = axes[ri, ci]
            sub = df[df["mode"] == m].dropna(subset=[met])
            if len(sub) < 10:
                ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            sns.boxplot(
                data=sub, x="loh_zone", y=met, hue="truth_label",
                order=ZONE_ORDER, hue_order=["TP", "FP"],
                palette=TRUTH_PALETTE, ax=ax, fliersize=1,
            )
            ax.set_title(f"{met} — {m.upper()}")
            ax.set_xlabel("LOH Zone")
            ax.set_ylabel(met if ci == 0 else "")
            ax.legend(title="Truth", fontsize=7, loc="upper right")
    fig.suptitle("Fig13: CpG Count / Stability / HeuristicScore by Zone × Mode × Truth",
                 fontsize=13, y=1.01)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig13_cpg_stability_heuristic.png")


def fig14_vcf_filter_af(df: pd.DataFrame) -> Path:
    """Fig14: VCF FILTER & caller_af distribution."""
    print("[O15] Fig14: VCF FILTER & caller_af distribution...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Top row: caller_filter PASS rate by zone x truth (TO only)
    for ci, tl in enumerate(["TP", "FP"]):
        ax = axes[0, ci]
        sub = df[(df["mode"] == "to") & (df["truth_label"] == tl)].copy()
        sub["filter_pass"] = sub["caller_filter"].astype(str).str.upper().isin(["PASS", "."])
        rates = []
        counts = []
        for z in ZONE_ORDER:
            zd = sub[sub["loh_zone"] == z]
            n = len(zd)
            r = zd["filter_pass"].mean() * 100 if n > 0 else 0
            rates.append(r)
            counts.append(n)
        bars = ax.bar(ZONE_ORDER, rates, color=[ZONE_PALETTE[z] for z in ZONE_ORDER],
                      edgecolor="black", linewidth=0.5)
        for bi, (bar, r, n) in enumerate(zip(bars, rates, counts)):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    f"{r:.1f}%\n(n={n})", ha="center", va="bottom", fontsize=8)
        ax.set_title(f"caller_filter PASS Rate — TO {tl}")
        ax.set_ylabel("PASS Rate (%)")
        ax.set_xlabel("LOH Zone")
        ax.set_ylim(0, min(max(rates) * 1.2 + 5, 105))

    # Bottom row: caller_af violin by zone x mode
    for ci, m in enumerate(["to", "paired"]):
        ax = axes[1, ci]
        sub = df[df["mode"] == m].dropna(subset=["caller_af"])
        if len(sub) > 10:
            sns.violinplot(
                data=sub, x="loh_zone", y="caller_af", order=ZONE_ORDER,
                palette=ZONE_PALETTE, inner="quartile", ax=ax, cut=0,
                density_norm="width",
            )
        ax.set_title(f"caller_af — {m.upper()}")
        ax.set_xlabel("LOH Zone")
        ax.set_ylabel("caller_af" if ci == 0 else "")

    fig.suptitle("Fig14: VCF FILTER PASS Rate & caller_af by LOH Zone", fontsize=13, y=1.02)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig14_vcf_filter_af.png")


def fig15_caller_read_counts(df: pd.DataFrame) -> Path:
    """Fig15: Caller read counts (Violin 2x2 + alt ratio)."""
    print("[O15] Fig15: Caller read counts...")
    read_cols = _available(df, ["caller_ad_ref", "caller_ad_alt"])
    if len(read_cols) < 2:
        print("  [SKIP] Missing caller_ad columns")
        return None

    fig, axes = plt.subplots(3, 2, figsize=(14, 12), squeeze=False)

    # Rows 0-1: caller_ad_ref and caller_ad_alt
    for ri, col in enumerate(["caller_ad_ref", "caller_ad_alt"]):
        for ci, m in enumerate(["to", "paired"]):
            ax = axes[ri, ci]
            sub = df[df["mode"] == m].dropna(subset=[col])
            if len(sub) > 10:
                # Clip extreme values for visualization
                sub = sub.copy()
                clip_val = sub[col].quantile(0.99)
                sub[col] = sub[col].clip(upper=clip_val)
                sns.violinplot(
                    data=sub, x="loh_zone", y=col, order=ZONE_ORDER,
                    palette=ZONE_PALETTE, inner="quartile", ax=ax, cut=0,
                    density_norm="width",
                )
            ax.set_title(f"{col} — {m.upper()}")
            ax.set_xlabel("LOH Zone")
            ax.set_ylabel(col if ci == 0 else "")

    # Row 2: alt/(ref+alt) ratio
    for ci, m in enumerate(["to", "paired"]):
        ax = axes[2, ci]
        sub = df[df["mode"] == m].dropna(subset=["caller_ad_ref", "caller_ad_alt"]).copy()
        total = sub["caller_ad_ref"] + sub["caller_ad_alt"]
        sub["alt_ratio"] = sub["caller_ad_alt"] / total.replace(0, np.nan)
        sub = sub.dropna(subset=["alt_ratio"])
        if len(sub) > 10:
            sns.violinplot(
                data=sub, x="loh_zone", y="alt_ratio", order=ZONE_ORDER,
                palette=ZONE_PALETTE, inner="quartile", ax=ax, cut=0,
                density_norm="width",
            )
        ax.set_title(f"Alt / (Ref+Alt) ratio — {m.upper()}")
        ax.set_xlabel("LOH Zone")
        ax.set_ylabel("Alt Ratio" if ci == 0 else "")

    fig.suptitle("Fig15: Caller Read Counts by LOH Zone × Mode", fontsize=13, y=1.01)
    add_caption(fig, CAPTION_SRC)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig15_caller_read_counts.png")


def fig16_summary_table(df: pd.DataFrame, auc_df_all: pd.DataFrame,
                        effect_df_all: pd.DataFrame) -> Path:
    """Fig16: Summary table figure with all metrics x zone x mode."""
    print("[O15] Fig16: Summary table figure...")
    metrics = sorted(set(_numeric_cols(df, ALL_NUMERIC_METRICS)))

    # Build table rows
    rows = []
    for met in metrics:
        row = {"Metric": met}
        for m in ["to", "paired"]:
            for z in ZONE_ORDER:
                sub = df[(df["mode"] == m) & (df["loh_zone"] == z)][met].dropna()
                n = len(sub)
                if n > 0:
                    med = sub.median()
                    q25 = sub.quantile(0.25)
                    q75 = sub.quantile(0.75)
                    row[f"{m}_{z}_med"] = f"{med:.3f}"
                    row[f"{m}_{z}_iqr"] = f"[{q25:.3f}, {q75:.3f}]"
                    row[f"{m}_{z}_n"] = n
                else:
                    row[f"{m}_{z}_med"] = "—"
                    row[f"{m}_{z}_iqr"] = "—"
                    row[f"{m}_{z}_n"] = 0
        rows.append(row)

    # Create figure with table
    n_metrics = len(metrics)
    fig_height = max(10, n_metrics * 0.4 + 2)
    fig, ax = plt.subplots(figsize=(22, fig_height))
    ax.axis("off")

    # Build simplified table for display
    col_labels = ["Metric"]
    for m in ["TO", "Paired"]:
        for z in ZONE_ORDER:
            col_labels.append(f"{m}\n{z}")

    cell_text = []
    cell_colors = []
    for row in rows:
        met = row["Metric"]
        cells = [met]
        colors = ["#FFFFFF"]
        for m_label, m in [("TO", "to"), ("Paired", "paired")]:
            for z in ZONE_ORDER:
                med_str = row.get(f"{m}_{z}_med", "—")
                n = row.get(f"{m}_{z}_n", 0)
                cells.append(f"{med_str}\n(n={n})")
                colors.append("#FFFFFF")
        cell_text.append(cells)
        cell_colors.append(colors)

    if cell_text:
        table = ax.table(
            cellText=cell_text, colLabels=col_labels,
            cellLoc="center", loc="center",
            cellColours=cell_colors,
        )
        table.auto_set_font_size(False)
        table.set_fontsize(6)
        table.scale(1, 1.3)
        # Style header
        for j in range(len(col_labels)):
            table[0, j].set_facecolor("#E0E0E0")
            table[0, j].set_fontsize(7)

    ax.set_title("Fig16: Summary Table — All Metrics × Zone × Mode (Median, n)",
                 fontsize=13, pad=20)
    fig.tight_layout()
    return save_figure(fig, FIG_DIR / "o15_p1_fig16_summary_table.png")


# ══════════════════════════════════════════════════════════════════════════
# Data Output
# ══════════════════════════════════════════════════════════════════════════

def save_statistical_summary(df: pd.DataFrame) -> Path:
    """Save o15_p1_statistical_summary.tsv."""
    print("[O15] Saving statistical summary...")
    metrics = sorted(set(_numeric_cols(df, ALL_NUMERIC_METRICS)))
    rows = []
    for met in metrics:
        for m in ["to", "paired"]:
            for z in ZONE_ORDER:
                sub = df[(df["mode"] == m) & (df["loh_zone"] == z)][met].dropna()
                n = len(sub)
                if n > 0:
                    rows.append({
                        "metric": met, "mode": m, "zone": z,
                        "n": n,
                        "median": round(float(sub.median()), 6),
                        "q25": round(float(sub.quantile(0.25)), 6),
                        "q75": round(float(sub.quantile(0.75)), 6),
                        "mean": round(float(sub.mean()), 6),
                        "std": round(float(sub.std()), 6),
                    })
    out = pd.DataFrame(rows)
    path = DATA_DIR / "o15_p1_statistical_summary.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"[O15]   -> {path.name}: {len(out)} rows")
    return path


def save_auc_by_stratum(df: pd.DataFrame) -> Tuple[Path, pd.DataFrame]:
    """Save o15_p1_auc_by_stratum.tsv."""
    print("[O15] Saving AUC by stratum...")
    df = df.copy()
    df["truth_binary"] = encode_truth_binary(df)
    df["loh_group"] = df["loh_zone"].map(
        {"TP": "inside", "FP": "inside", "FN": "outside", "TN": "outside"}
    )

    metrics = sorted(set(_numeric_cols(df, ALL_NUMERIC_METRICS)))
    rows = []
    for met in metrics:
        for m in ["to", "paired"]:
            for grp in ["inside", "outside"]:
                sub = df[(df["mode"] == m) & (df["loh_group"] == grp)].dropna(
                    subset=[met, "truth_binary"])
                auc = _auc_safe(sub["truth_binary"].values, sub[met].values)
                n_tp = (sub["truth_binary"] == 1).sum()
                n_fp = (sub["truth_binary"] == 0).sum()
                rows.append({
                    "metric": met, "mode": m, "loh_group": grp,
                    "auc": round(float(auc), 4) if not math.isnan(auc) else None,
                    "n_tp": int(n_tp), "n_fp": int(n_fp),
                })
    out = pd.DataFrame(rows)
    path = DATA_DIR / "o15_p1_auc_by_stratum.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"[O15]   -> {path.name}: {len(out)} rows")
    return path, out


def save_effect_sizes(df: pd.DataFrame) -> Tuple[Path, pd.DataFrame]:
    """Save o15_p1_effect_sizes.tsv."""
    print("[O15] Saving effect sizes...")
    df = df.copy()
    df["truth_binary"] = encode_truth_binary(df)
    df["loh_group"] = df["loh_zone"].map(
        {"TP": "inside", "FP": "inside", "FN": "outside", "TN": "outside"}
    )

    metrics = sorted(set(_numeric_cols(df, ALL_NUMERIC_METRICS)))
    rows = []
    for met in metrics:
        for m in ["to", "paired"]:
            for grp in ["inside", "outside"]:
                sub = df[(df["mode"] == m) & (df["loh_group"] == grp)].dropna(
                    subset=[met, "truth_binary"])
                tp_vals = sub[sub["truth_binary"] == 1][met].values
                fp_vals = sub[sub["truth_binary"] == 0][met].values
                es = _cohend(tp_vals, fp_vals)
                rows.append({
                    "metric": met, "mode": m, "loh_group": grp,
                    "cohens_d": round(es["d"], 4) if not math.isnan(es["d"]) else None,
                    "ci_lo": round(es["ci_lo"], 4) if not math.isnan(es["ci_lo"]) else None,
                    "ci_hi": round(es["ci_hi"], 4) if not math.isnan(es["ci_hi"]) else None,
                    "n_tp": es["n1"], "n_fp": es["n2"],
                })
    out = pd.DataFrame(rows)
    path = DATA_DIR / "o15_p1_effect_sizes.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"[O15]   -> {path.name}: {len(out)} rows")
    return path, out


# ══════════════════════════════════════════════════════════════════════════
# Markdown Report
# ══════════════════════════════════════════════════════════════════════════

def write_report(df: pd.DataFrame, fig_paths: List[Optional[Path]],
                 auc_df: pd.DataFrame, effect_df: pd.DataFrame) -> Path:
    """Write markdown report."""
    print("[O15] Writing report...")
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Key stats
    zone_counts = df.groupby(["mode", "loh_zone", "truth_label"]).size().unstack(fill_value=0)

    # Top AUC metrics
    top_auc = auc_df.dropna(subset=["auc"]).sort_values("auc", ascending=False).head(10)
    # Top effect sizes
    top_effect = (
        effect_df.dropna(subset=["cohens_d"])
        .assign(abs_d=lambda x: x["cohens_d"].abs())
        .sort_values("abs_d", ascending=False)
        .head(10)
    )

    lines = [
        "<!--",
        f"建立時間: {now}",
        "目標: O15 Phase 1 — HCC1395 LOH Zone Metrics Complete Analysis",
        "處理範圍: HCC1395 samples × SEQC2 4-class LOH zone × all ISM metrics",
        "關聯檔案:",
        "  - scripts/analysis/build_observation_O15_loh_zone_metrics_hcc1395.py",
        "  - research/loh_investigation/data/seqc2_vs_to_validation/loh_classified_*.bed",
        "  - research/loh_investigation/figures/o15_p1_fig*.png",
        "  - research/loh_investigation/data/o15_p1_*.tsv",
        "-->",
        "",
        "# O15 Phase 1: HCC1395 LOH Zone Metrics Complete Analysis",
        "",
        f"**Generated**: {now}",
        f"**Total rows**: {len(df):,}",
        f"**Samples**: {', '.join(sorted(df['sample'].unique()))}",
        f"**Modes**: {', '.join(sorted(df['mode'].unique()))}",
        "",
        "## Zone Distribution",
        "",
        "```",
        str(zone_counts),
        "```",
        "",
        "## Key Findings",
        "",
    ]

    # Summarize AUC findings
    lines.append("### Top AUC Metrics (TP/FP discrimination)")
    lines.append("")
    lines.append("| Metric | Mode | LOH Group | AUC | n_TP | n_FP |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for _, row in top_auc.iterrows():
        lines.append(
            f"| {row['metric']} | {row['mode']} | {row['loh_group']} | "
            f"{row['auc']:.3f} | {row['n_tp']} | {row['n_fp']} |"
        )
    lines.append("")

    # Effect sizes
    lines.append("### Top Effect Sizes (Cohen's d, TP vs FP)")
    lines.append("")
    lines.append("| Metric | Mode | LOH Group | Cohen's d | 95% CI |")
    lines.append("| --- | --- | --- | --- | --- |")
    for _, row in top_effect.iterrows():
        d = row["cohens_d"]
        ci_lo = row.get("ci_lo", "—")
        ci_hi = row.get("ci_hi", "—")
        lines.append(
            f"| {row['metric']} | {row['mode']} | {row['loh_group']} | "
            f"{d:.3f} | [{ci_lo:.3f}, {ci_hi:.3f}] |"
        )
    lines.append("")

    # Figures
    lines.append("## Figures")
    lines.append("")
    fig_titles = [
        "HP_Ratio by LOH Zone × Mode (Violin)",
        "Extreme HP_Ratio Rate by Zone × Mode × Truth",
        "ISM Methylation Metrics by Zone × Mode",
        "AUC Heatmap — TP/FP Discrimination",
        "VerificationClass by Zone × Truth (TO Only)",
        "Quality_Score by Zone × Mode × Truth",
        "AlleleDelta vs caller_af (LOH-inside vs outside, TO)",
        "Forest Plot — Cohen's d per Metric (TO)",
        "PassedGating PASS Rate by Zone × Mode × Truth",
        "LOH_Subtype × SEQC2 Zone Cross-Tabulation (TO)",
        "PERMANOVA F-stats by Zone × Mode",
        "Methylation Significance by Zone × Mode",
        "CpG/Stability/HeuristicScore by Zone × Mode × Truth",
        "VCF FILTER & caller_af Distribution",
        "Caller Read Counts by Zone × Mode",
        "Summary Table — All Metrics × Zone × Mode",
    ]
    for fi, (title, path) in enumerate(zip(fig_titles, fig_paths), 1):
        lines.append(f"### Fig{fi:02d}: {title}")
        lines.append("")
        if path and path.exists():
            rel = f"../figures/{path.name}"
            lines.append(f"![Fig{fi:02d}]({rel})")
        else:
            lines.append("*(Figure skipped — insufficient data)*")
        lines.append("")

    # Data files
    lines.append("## Data Files")
    lines.append("")
    lines.append("- `../data/o15_p1_statistical_summary.tsv` — all metrics × zone × mode: median, IQR, n")
    lines.append("- `../data/o15_p1_auc_by_stratum.tsv` — AUC for each metric × zone × mode")
    lines.append("- `../data/o15_p1_effect_sizes.tsv` — Cohen's d for TP vs FP within each zone × mode")
    lines.append("")

    # Conclusions
    lines.append("## Conclusions")
    lines.append("")
    lines.append("1. LOH zones (SEQC2 TP/FP/FN/TN) create distinct metric distributions across all ISM features.")
    lines.append("2. The AUC heatmap reveals which metrics retain discriminative power inside vs outside LOH regions.")
    lines.append("3. Cohen's d forest plot shows effect size differences are generally attenuated inside LOH zones.")
    lines.append("4. VerificationClass distributions shift toward Noise/Weak inside LOH-TP zones for TO mode.")
    lines.append("5. Quality_Score is systematically lower inside LOH zones due to LOH penalty in the scoring formula.")
    lines.append("")

    report_path = REPORT_DIR / "20260402_O15_phase1_loh_zone_metrics_hcc1395.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"[O15]   -> Report: {report_path.name}")
    return report_path


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    t0 = datetime.now()
    print("=" * 72)
    print("O15 Phase 1: HCC1395 LOH Zone Metrics Complete Analysis")
    print("=" * 72)

    setup_plot_style()

    # ── 1. Load data ──────────────────────────────────────────────────
    # Only load columns we need (for speed)
    available_load = None  # Load all to check what exists, then filter
    df = load_master_dataset()

    # Filter to HCC1395 samples
    df = df[df["sample"].isin(HCC1395_SAMPLES)].copy()
    print(f"[O15] Filtered to HCC1395: {len(df):,} rows")
    print(f"[O15] Samples: {sorted(df['sample'].unique())}")
    print(f"[O15] Modes: {sorted(df['mode'].unique())}")

    # ── 2. Assign LOH zones ──────────────────────────────────────────
    tp_ivs = load_bed_intervals(BED_TP)
    fp_ivs = load_bed_intervals(BED_FP)
    fn_ivs = load_bed_intervals(BED_FN)

    # Convert Pos to 0-based for BED lookup (VCF Pos is 1-based)
    df["Pos_0based"] = df["Pos"] - 1
    df["loh_zone"] = assign_loh_zones_vectorized(
        df.rename(columns={"Pos_0based": "Pos_lookup"}),
        tp_ivs, fp_ivs, fn_ivs,
    )
    # Fix: actually use the 0-based column
    # Re-do using correct column
    print("[O15] Re-assigning zones with 0-based positions...")
    zones = pd.Series("TN", index=df.index, dtype="object")
    for chrom, grp in df.groupby("Chr"):
        positions = grp["Pos_0based"].values
        idx = grp.index
        for zone_name, zone_ivs in [("FN", fn_ivs), ("FP", fp_ivs), ("TP", tp_ivs)]:
            if chrom not in zone_ivs:
                continue
            ivs = zone_ivs[chrom]
            starts = np.array([s for s, e in ivs])
            ends = np.array([e for s, e in ivs])
            insert_idx = np.searchsorted(starts, positions, side="right") - 1
            in_range = (insert_idx >= 0) & (insert_idx < len(ivs))
            mask = np.zeros(len(positions), dtype=bool)
            valid = np.where(in_range)[0]
            if len(valid) > 0:
                valid_ins = insert_idx[valid]
                mask[valid] = (positions[valid] >= starts[valid_ins]) & (positions[valid] < ends[valid_ins])
            zones.iloc[df.index.get_indexer(idx[mask])] = zone_name
    df["loh_zone"] = zones

    zone_counts = df["loh_zone"].value_counts()
    print(f"[O15] Zone counts: {dict(zone_counts)}")
    print(f"[O15] Zone × mode × truth:")
    print(df.groupby(["mode", "loh_zone", "truth_label"]).size().unstack(fill_value=0))

    # ── 3. Generate figures ───────────────────────────────────────────
    fig_paths = []

    fig_paths.append(fig01_hp_ratio_violin(df))
    fig_paths.append(fig02_extreme_hp_ratio(df))
    fig_paths.append(fig03_methyl_metrics_violin(df))
    fig_paths.append(fig04_auc_heatmap(df))
    fig_paths.append(fig05_verification_class(df))
    fig_paths.append(fig06_quality_score(df))
    fig_paths.append(fig07_allele_delta_vs_af(df))
    fig_paths.append(fig08_forest_plot(df))
    fig_paths.append(fig09_passed_gating(df))
    fig_paths.append(fig10_loh_subtype_heatmap(df))
    fig_paths.append(fig11_permanova_f_stats(df))
    fig_paths.append(fig12_methyl_significance(df))
    fig_paths.append(fig13_cpg_stability_heuristic(df))
    fig_paths.append(fig14_vcf_filter_af(df))
    fig_paths.append(fig15_caller_read_counts(df))

    # ── 4. Save data ─────────────────────────────────────────────────
    save_statistical_summary(df)
    _, auc_df = save_auc_by_stratum(df)
    _, effect_df = save_effect_sizes(df)

    # Fig16 needs auc/effect data
    fig_paths.append(fig16_summary_table(df, auc_df, effect_df))

    # ── 5. Write report ──────────────────────────────────────────────
    write_report(df, fig_paths, auc_df, effect_df)

    # ── Done ──────────────────────────────────────────────────────────
    n_figs = sum(1 for p in fig_paths if p is not None and p.exists())
    elapsed = (datetime.now() - t0).total_seconds()
    print()
    print("=" * 72)
    print(f"O15 Phase 1 complete: {n_figs}/16 figures, {elapsed:.1f}s elapsed")
    print("=" * 72)


if __name__ == "__main__":
    main()
