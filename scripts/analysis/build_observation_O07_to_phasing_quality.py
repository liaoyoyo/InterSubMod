#!/usr/bin/env python3
"""O7: TO Phasing Quality Observation.

Systematically observes tumour-only (TO) phasing quality post HP-fix,
covering HP tag composition, phase assignment rates, HP0/HP3 fractions,
concordance with paired mode, tier distributions, chromosome-level patterns,
and HP imbalance vs LOH.

Usage:
    python3 build_observation_O07_to_phasing_quality.py
"""

from __future__ import annotations

import sys
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

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

import json
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats
from pathlib import Path

# ── Constants ────────────────────────────────────────────────────────────
OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260401_O07_to_phasing_quality")
CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz "
    "(20260327_loh_round1_cross_sample_audit, post-HP-fix) | "
)

HP_PHASING_COLS = [
    "HP1FamilyN", "HP2FamilyN", "NHP0", "NHP3", "NumReads",
    "hp_assign_rate", "hp0_ratio", "hp3_ratio",
    "HP_Ratio", "hp_ratio_core", "effective_hp_reads",
    "Quality_Tier", "Potential_LOH", "core_loh_like",
    "to_phase_ps", "to_phase_gt", "phase_block_status", "to_loh_bed_hit",
    "mode", "truth_label", "sample", "Chr", "variant_key",
    "Coverage_Multiple", "CramersV",
]


# ══════════════════════════════════════════════════════════════════════════
# Helper utilities
# ══════════════════════════════════════════════════════════════════════════

def _available_cols(df, cols):
    """Return list of columns that exist in df."""
    return [c for c in cols if c in df.columns]


def _auc_safe(y_true, y_score):
    """ROC AUC with graceful fallback."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_score)
    y_true, y_score = y_true[mask], y_score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 10:
        return float("nan")
    return compute_auc(y_true, y_score)


def _chrom_sort_key(c):
    """Sort chromosomes: chr1..chr22, chrX, chrY, chrM."""
    c = str(c).replace("chr", "")
    if c.isdigit():
        return (0, int(c))
    elif c == "X":
        return (1, 0)
    elif c == "Y":
        return (1, 1)
    else:
        return (2, 0)


# ══════════════════════════════════════════════════════════════════════════
# Figure 01: TO HP Tag Composition (stacked bar per sample)
# ══════════════════════════════════════════════════════════════════════════

def fig01_to_hp_tag_composition(df, out):
    """Stacked bar: HP1, HP2, HP0, HP3 fractions per sample (TO only)."""
    to = df[df["mode"] == "to"].copy()
    to["HP0_frac"] = to["NHP0"] / to["NumReads"]
    to["HP3_frac"] = to["NHP3"] / to["NumReads"]
    to["HP1_frac"] = to["HP1FamilyN"] / to["NumReads"]
    to["HP2_frac"] = to["HP2FamilyN"] / to["NumReads"]

    sample_order = [s for s in SAMPLE_ORDER if s in to["sample"].unique()]

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Left: stacked bar means
    ax = axes[0]
    means = to.groupby("sample")[["HP1_frac", "HP2_frac", "HP0_frac", "HP3_frac"]].mean()
    means = means.reindex(sample_order)
    colors = ["#1976D2", "#D32F2F", "#9E9E9E", "#FF9800"]
    labels = ["HP1", "HP2", "HP0 (unphased)", "HP3 (ambiguous)"]
    bottoms = np.zeros(len(sample_order))
    x = np.arange(len(sample_order))
    for i, (col, color, label) in enumerate(zip(
        ["HP1_frac", "HP2_frac", "HP0_frac", "HP3_frac"], colors, labels
    )):
        vals = means[col].values
        ax.bar(x, vals, bottom=bottoms, color=color, label=label, edgecolor="white", linewidth=0.5)
        # Annotate percentages
        for j, v in enumerate(vals):
            if v > 0.03:
                ax.text(j, bottoms[j] + v / 2, f"{v:.1%}", ha="center", va="center", fontsize=7,
                        color="white" if i < 2 else "black")
        bottoms += vals
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=35, ha="right")
    ax.set_ylabel("Fraction of reads")
    ax.set_title("Mean HP Tag Composition per Sample (TO only)", fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_ylim(0, 1.05)

    # Right: by truth label
    ax = axes[1]
    for truth_idx, truth in enumerate(["TP", "FP"]):
        sub = to[to["truth_label"] == truth]
        means_t = sub.groupby("sample")[["HP1_frac", "HP2_frac", "HP0_frac", "HP3_frac"]].mean()
        means_t = means_t.reindex(sample_order)
        width = 0.35
        offset = (truth_idx - 0.5) * width
        bottoms_t = np.zeros(len(sample_order))
        for i, (col, color) in enumerate(zip(
            ["HP1_frac", "HP2_frac", "HP0_frac", "HP3_frac"], colors
        )):
            vals = means_t[col].fillna(0).values
            ax.bar(x + offset, vals, width=width, bottom=bottoms_t, color=color,
                   edgecolor="white", linewidth=0.3, alpha=0.85 if truth == "TP" else 0.55)
            bottoms_t += vals
    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=c, label=l) for c, l in zip(colors, labels)
    ] + [
        Patch(facecolor="gray", alpha=0.85, label="TP (left bar)"),
        Patch(facecolor="gray", alpha=0.55, label="FP (right bar)"),
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=35, ha="right")
    ax.set_ylabel("Fraction of reads")
    ax.set_title("HP Composition by Truth Label (TO)", fontweight="bold")
    ax.legend(handles=legend_elements, loc="upper right", fontsize=7, ncol=2)
    ax.set_ylim(0, 1.05)

    fig.suptitle("Fig 01: TO HP Tag Composition", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig01")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig01_to_hp_tag_composition.png")


# ══════════════════════════════════════════════════════════════════════════
# Figure 02: hp_assign_rate Distribution (TO vs Paired)
# ══════════════════════════════════════════════════════════════════════════

def fig02_to_hp_assign_rate_distribution(df, out):
    """hp_assign_rate histogram and violin: TO vs paired."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: histogram overlay
    ax = axes[0]
    for mode, color, label in [("paired", COLOR_PAIRED, "Paired"), ("to", COLOR_TO, "TO")]:
        vals = df.loc[df["mode"] == mode, "hp_assign_rate"].dropna()
        ax.hist(vals, bins=50, alpha=0.55, color=color, label=f"{label} (n={len(vals):,})",
                edgecolor="none")
    ax.set_xlabel("hp_assign_rate")
    ax.set_ylabel("Count")
    ax.set_title("hp_assign_rate Distribution", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel 2: violin by mode x truth
    ax = axes[1]
    plot_df = df[["hp_assign_rate", "mode", "truth_label"]].dropna()
    plot_df["mode"] = pd.Categorical(plot_df["mode"], categories=MODE_ORDER, ordered=True)
    sns.violinplot(
        data=plot_df, x="mode", y="hp_assign_rate", hue="truth_label",
        palette=TRUTH_PALETTE, ax=ax, split=True, inner="quartile", cut=0, linewidth=0.8,
    )
    ax.set_title("hp_assign_rate by Mode x Truth", fontweight="bold")
    ax.legend(title="Truth", fontsize=8)

    # Panel 3: summary stats
    ax = axes[2]
    ax.axis("off")
    rows_tbl = []
    for mode in MODE_ORDER:
        for truth in ["TP", "FP"]:
            mask = (df["mode"] == mode) & (df["truth_label"] == truth)
            vals = df.loc[mask, "hp_assign_rate"].dropna()
            rows_tbl.append({
                "Mode": mode, "Truth": truth,
                "N": f"{len(vals):,}",
                "Mean": f"{vals.mean():.4f}",
                "Median": f"{vals.median():.4f}",
                "P25": f"{vals.quantile(0.25):.4f}",
                "P75": f"{vals.quantile(0.75):.4f}",
            })
    tbl = ax.table(
        cellText=[[r[k] for k in ["Mode", "Truth", "N", "Mean", "Median", "P25", "P75"]] for r in rows_tbl],
        colLabels=["Mode", "Truth", "N", "Mean", "Median", "P25", "P75"],
        loc="center", cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.4)
    ax.set_title("Summary Statistics", fontweight="bold")

    fig.suptitle("Fig 02: hp_assign_rate Distribution (TO vs Paired)", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig02")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig02_to_hp_assign_rate_distribution.png")

    # Return stats for report
    stats = {}
    for mode in MODE_ORDER:
        vals = df.loc[df["mode"] == mode, "hp_assign_rate"].dropna()
        stats[mode] = {"mean": vals.mean(), "median": vals.median(), "n": len(vals)}
    return stats


# ══════════════════════════════════════════════════════════════════════════
# Figure 03: Phase Block Size Distribution
# ══════════════════════════════════════════════════════════════════════════

def fig03_to_phase_block_size_distribution(df, out):
    """Phase block (PS) distribution for TO: block size proxy via grouped counts."""
    to = df[df["mode"] == "to"].copy()

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: phase_block_status pie
    ax = axes[0]
    status_counts = to["phase_block_status"].value_counts()
    colors_pie = ["#4CAF50", "#FF5722", "#9E9E9E"][:len(status_counts)]
    ax.pie(status_counts.values, labels=[f"{k}\n({v:,})" for k, v in status_counts.items()],
           autopct="%1.1f%%", colors=colors_pie, startangle=90)
    ax.set_title("Phase Block Status (TO)", fontweight="bold")

    # Panel 2: effective_hp_reads by phase_block_status
    ax = axes[1]
    has_ps = to[to["phase_block_status"] == "to_ps_available"]
    no_ps = to[to["phase_block_status"] == "to_ps_missing"]
    for sub, label, color in [
        (has_ps, f"PS available (n={len(has_ps):,})", "#4CAF50"),
        (no_ps, f"PS missing (n={len(no_ps):,})", "#FF5722"),
    ]:
        vals = sub["effective_hp_reads"].dropna().clip(lower=1)
        if len(vals) > 0:
            ax.hist(np.log10(vals), bins=60, alpha=0.55, color=color, label=label, edgecolor="none")
    ax.set_xlabel("log10(effective_hp_reads)")
    ax.set_ylabel("Count")
    ax.set_title("Effective HP Reads by PS Status", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel 3: variants per PS block (top PS blocks)
    ax = axes[2]
    ps_valid = to.loc[to["to_phase_ps"].notna(), "to_phase_ps"]
    if len(ps_valid) > 0:
        block_sizes = ps_valid.value_counts()
        ax.hist(block_sizes.values, bins=np.arange(0.5, min(block_sizes.max() + 1.5, 50.5), 1),
                color="#7B1FA2", edgecolor="white", alpha=0.8)
        ax.set_xlabel("Variants per phase block")
        ax.set_ylabel("Number of phase blocks")
        ax.set_title(f"Phase Block Size (n_blocks={len(block_sizes):,})", fontweight="bold")
        # Add stats annotation
        ax.axvline(block_sizes.median(), color="red", ls="--", lw=1.5, label=f"Median={block_sizes.median():.0f}")
        ax.axvline(block_sizes.mean(), color="blue", ls="--", lw=1.5, label=f"Mean={block_sizes.mean():.1f}")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "No PS data available", transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
        ax.set_title("Phase Block Size (no data)", fontweight="bold")

    fig.suptitle("Fig 03: TO Phase Block Distribution", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig03")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig03_to_phase_block_size_distribution.png")


# ══════════════════════════════════════════════════════════════════════════
# Figure 04: HP0 Fraction vs Truth (TO)
# ══════════════════════════════════════════════════════════════════════════

def fig04_to_hp0_fraction_vs_truth(df, out):
    """HP0 fraction (NHP0/NumReads) vs truth label, TO only."""
    to = df[df["mode"] == "to"].copy()
    to["hp0_frac"] = to["NHP0"] / to["NumReads"]

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: violin by truth
    ax = axes[0]
    plot_df = to[["hp0_frac", "truth_label"]].dropna()
    sns.violinplot(
        data=plot_df, x="truth_label", y="hp0_frac",
        palette=TRUTH_PALETTE, ax=ax, inner="quartile", cut=0, linewidth=0.8,
        order=["TP", "FP"],
    )
    ax.set_title("HP0 Fraction by Truth (TO)", fontweight="bold")
    ax.set_ylabel("NHP0 / NumReads")

    # Panel 2: by sample
    ax = axes[1]
    sample_order = [s for s in SAMPLE_ORDER if s in to["sample"].unique()]
    plot_df2 = to[["hp0_frac", "truth_label", "sample"]].dropna()
    plot_df2["sample"] = pd.Categorical(plot_df2["sample"], categories=sample_order, ordered=True)
    sns.boxplot(
        data=plot_df2, x="sample", y="hp0_frac", hue="truth_label",
        palette=TRUTH_PALETTE, ax=ax, showfliers=False, linewidth=0.8,
    )
    ax.set_title("HP0 Fraction by Sample x Truth (TO)", fontweight="bold")
    ax.tick_params(axis="x", rotation=35)
    ax.set_ylabel("NHP0 / NumReads")
    ax.legend(title="Truth", fontsize=7)

    # Panel 3: cumulative distribution
    ax = axes[2]
    for truth, color in TRUTH_PALETTE.items():
        vals = to.loc[to["truth_label"] == truth, "hp0_frac"].dropna().sort_values()
        if len(vals) > 0:
            y = np.arange(1, len(vals) + 1) / len(vals)
            ax.plot(vals, y, color=color, label=f"{truth} (n={len(vals):,})", linewidth=1.5)
    ax.set_xlabel("HP0 Fraction")
    ax.set_ylabel("Cumulative proportion")
    ax.set_title("CDF of HP0 Fraction (TO)", fontweight="bold")
    ax.legend(fontsize=8)

    fig.suptitle("Fig 04: HP0 Fraction vs Truth (TO Only)", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig04")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig04_to_hp0_fraction_vs_truth.png")

    # Stats
    tp_vals = to.loc[to["truth_label"] == "TP", "hp0_frac"].dropna().values
    fp_vals = to.loc[to["truth_label"] == "FP", "hp0_frac"].dropna().values
    mw = mannwhitney_test(tp_vals, fp_vals)
    es = compute_effect_size(tp_vals, fp_vals)
    return {"mannwhitney": mw, "effect_size": es,
            "tp_mean": float(np.nanmean(tp_vals)), "fp_mean": float(np.nanmean(fp_vals))}


# ══════════════════════════════════════════════════════════════════════════
# Figure 05: TO/Paired HP Ratio Concordance (same variant_key)
# ══════════════════════════════════════════════════════════════════════════

def fig05_to_paired_hp_concordance(df, out):
    """HP_Ratio concordance: same variant in TO vs paired mode."""
    paired = df[df["mode"] == "paired"][["variant_key", "HP_Ratio", "truth_label", "sample"]].copy()
    paired = paired.rename(columns={"HP_Ratio": "HP_Ratio_paired"})
    to = df[df["mode"] == "to"][["variant_key", "HP_Ratio", "truth_label", "sample"]].copy()
    to = to.rename(columns={"HP_Ratio": "HP_Ratio_to"})

    merged = paired.merge(to, on=["variant_key", "sample"], suffixes=("_paired", "_to"))
    # Use truth from paired as reference
    merged["truth"] = merged.get("truth_label_paired", merged.get("truth_label"))
    if "truth" not in merged.columns or merged["truth"].isna().all():
        merged["truth"] = merged["truth_label"]

    print(f"[O7 Fig05] Matched variant_keys: {len(merged):,}")

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    # Panel 1: scatter
    ax = axes[0]
    subsample = merged.sample(n=min(50000, len(merged)), random_state=42)
    for truth, color in TRUTH_PALETTE.items():
        mask = subsample["truth"] == truth
        ax.scatter(
            subsample.loc[mask, "HP_Ratio_paired"],
            subsample.loc[mask, "HP_Ratio_to"],
            alpha=0.15, s=3, color=color, label=truth, rasterized=True,
        )
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.5, label="y=x")
    ax.set_xlabel("HP_Ratio (Paired)")
    ax.set_ylabel("HP_Ratio (TO)")
    ax.set_title("HP_Ratio Concordance", fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)

    # Compute correlation
    valid = merged[["HP_Ratio_paired", "HP_Ratio_to"]].dropna()
    if len(valid) > 10:
        r, p_corr = sp_stats.pearsonr(valid["HP_Ratio_paired"], valid["HP_Ratio_to"])
        rho, p_spearman = sp_stats.spearmanr(valid["HP_Ratio_paired"], valid["HP_Ratio_to"])
        ax.text(0.05, 0.95, f"Pearson r={r:.3f}\nSpearman rho={rho:.3f}\nn={len(valid):,}",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    else:
        r, rho = float("nan"), float("nan")

    # Panel 2: difference distribution
    ax = axes[1]
    merged["hp_diff"] = merged["HP_Ratio_to"] - merged["HP_Ratio_paired"]
    for truth, color in TRUTH_PALETTE.items():
        vals = merged.loc[merged["truth"] == truth, "hp_diff"].dropna()
        if len(vals) > 0:
            ax.hist(vals, bins=80, alpha=0.55, color=color,
                    label=f"{truth} (n={len(vals):,}, mean={vals.mean():.4f})", edgecolor="none")
    ax.axvline(0, color="black", ls="--", lw=0.8)
    ax.set_xlabel("HP_Ratio(TO) - HP_Ratio(Paired)")
    ax.set_ylabel("Count")
    ax.set_title("HP_Ratio Difference Distribution", fontweight="bold")
    ax.legend(fontsize=7)

    # Panel 3: per-sample correlation
    ax = axes[2]
    sample_order = [s for s in SAMPLE_ORDER if s in merged["sample"].unique()]
    corr_rows = []
    for samp in sample_order:
        sub = merged[merged["sample"] == samp][["HP_Ratio_paired", "HP_Ratio_to"]].dropna()
        if len(sub) > 10:
            r_s, _ = sp_stats.pearsonr(sub["HP_Ratio_paired"], sub["HP_Ratio_to"])
            corr_rows.append({"sample": samp, "pearson_r": r_s, "n": len(sub)})
    if corr_rows:
        corr_df = pd.DataFrame(corr_rows)
        bars = ax.barh(corr_df["sample"], corr_df["pearson_r"], color=COLOR_TO, edgecolor="white")
        for bar, row in zip(bars, corr_rows):
            ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                    f"r={row['pearson_r']:.3f} (n={row['n']:,})", va="center", fontsize=8)
        ax.set_xlabel("Pearson r")
        ax.set_title("Per-Sample HP_Ratio Correlation (TO vs Paired)", fontweight="bold")
        ax.set_xlim(0, 1.15)
    else:
        ax.text(0.5, 0.5, "Insufficient matched data", transform=ax.transAxes,
                ha="center", va="center")

    fig.suptitle("Fig 05: TO/Paired HP_Ratio Concordance", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig05")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig05_to_paired_hp_concordance.png")

    return {
        "n_matched": len(merged),
        "pearson_r": float(r) if not np.isnan(r) else None,
        "spearman_rho": float(rho) if not np.isnan(rho) else None,
        "mean_diff": float(merged["hp_diff"].mean()) if len(merged) > 0 else None,
        "per_sample_corr": corr_rows,
    }


# ══════════════════════════════════════════════════════════════════════════
# Figure 06: Quality_Tier Distribution Post Fix (TO only)
# ══════════════════════════════════════════════════════════════════════════

def fig06_to_tier_distribution_post_fix(df, out):
    """Quality_Tier stacked bar per sample for TO, compared to paired."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    tier_order = ["High", "Medium", "Low"]
    tier_colors = {"High": "#4CAF50", "Medium": "#FFC107", "Low": "#F44336"}

    for ax_idx, (mode, title_mode) in enumerate([("to", "TO"), ("paired", "Paired")]):
        ax = axes[ax_idx]
        sub = df[df["mode"] == mode]
        sample_order = [s for s in SAMPLE_ORDER if s in sub["sample"].unique()]
        ct = pd.crosstab(sub["sample"], sub["Quality_Tier"], normalize="index")
        ct = ct.reindex(index=sample_order, columns=[t for t in tier_order if t in ct.columns], fill_value=0)

        x = np.arange(len(sample_order))
        bottoms = np.zeros(len(sample_order))
        for tier in [t for t in tier_order if t in ct.columns]:
            vals = ct[tier].values
            ax.bar(x, vals, bottom=bottoms, color=tier_colors.get(tier, "#999"),
                   label=tier, edgecolor="white", linewidth=0.5)
            for j, v in enumerate(vals):
                if v > 0.04:
                    ax.text(j, bottoms[j] + v / 2, f"{v:.0%}", ha="center", va="center", fontsize=7)
            bottoms += vals
        ax.set_xticks(x)
        ax.set_xticklabels(sample_order, rotation=35, ha="right")
        ax.set_ylabel("Proportion")
        ax.set_title(f"Quality_Tier Distribution ({title_mode})", fontweight="bold")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_ylim(0, 1.05)

    fig.suptitle("Fig 06: Quality_Tier Distribution Post HP-Fix", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig06")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig06_to_tier_distribution_post_fix.png")

    # Stats: tier proportions
    to_tiers = df[df["mode"] == "to"]["Quality_Tier"].value_counts(normalize=True).to_dict()
    paired_tiers = df[df["mode"] == "paired"]["Quality_Tier"].value_counts(normalize=True).to_dict()
    return {"to_tiers": to_tiers, "paired_tiers": paired_tiers}


# ══════════════════════════════════════════════════════════════════════════
# Figure 07: Phase Quality by Chromosome (TO only)
# ══════════════════════════════════════════════════════════════════════════

def fig07_to_phase_quality_by_chromosome(df, out):
    """Box plot of hp_assign_rate per chromosome (TO only)."""
    to = df[df["mode"] == "to"].copy()

    # Sort chromosomes
    chroms = sorted(to["Chr"].dropna().unique(), key=_chrom_sort_key)
    # Keep autosomes + X
    chroms = [c for c in chroms if c.replace("chr", "").isdigit() or c.replace("chr", "") == "X"]

    fig, axes = plt.subplots(2, 1, figsize=(20, 12))

    # Panel 1: hp_assign_rate by chrom
    ax = axes[0]
    plot_data = to[to["Chr"].isin(chroms)].copy()
    plot_data["Chr"] = pd.Categorical(plot_data["Chr"], categories=chroms, ordered=True)
    sns.boxplot(
        data=plot_data, x="Chr", y="hp_assign_rate",
        color=COLOR_TO, ax=ax, showfliers=False, linewidth=0.6,
    )
    ax.set_title("hp_assign_rate by Chromosome (TO)", fontweight="bold")
    ax.set_ylabel("hp_assign_rate")
    ax.tick_params(axis="x", rotation=45)

    # Add median line
    overall_median = to["hp_assign_rate"].median()
    ax.axhline(overall_median, color="red", ls="--", lw=1, label=f"Overall median={overall_median:.4f}")
    ax.legend(fontsize=8)

    # Panel 2: HP_Ratio by chrom
    ax = axes[1]
    sns.boxplot(
        data=plot_data, x="Chr", y="HP_Ratio",
        color=COLOR_TO, ax=ax, showfliers=False, linewidth=0.6,
    )
    ax.set_title("HP_Ratio by Chromosome (TO)", fontweight="bold")
    ax.set_ylabel("HP_Ratio")
    ax.tick_params(axis="x", rotation=45)
    ax.axhline(0.5, color="gray", ls="--", lw=0.8, alpha=0.5, label="HP_Ratio=0.5 (balanced)")
    ax.legend(fontsize=8)

    fig.suptitle("Fig 07: TO Phase Quality by Chromosome", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig07")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig07_to_phase_quality_by_chromosome.png")

    # Per-chrom stats
    chrom_stats = to.groupby("Chr").agg(
        n=("hp_assign_rate", "count"),
        hp_assign_mean=("hp_assign_rate", "mean"),
        hp_assign_median=("hp_assign_rate", "median"),
        hp_ratio_mean=("HP_Ratio", "mean"),
    ).reindex(chroms)
    return chrom_stats


# ══════════════════════════════════════════════════════════════════════════
# Figure 08: HP Imbalance vs LOH Determination
# ══════════════════════════════════════════════════════════════════════════

def fig08_to_hp_imbalance_vs_loh(df, out):
    """HP imbalance (|HP_Ratio - 0.5|) vs LOH + truth for TO."""
    to = df[df["mode"] == "to"].copy()
    to["hp_imbalance"] = (to["HP_Ratio"] - 0.5).abs()

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel 1: imbalance by LOH status
    ax = axes[0, 0]
    to["LOH_str"] = to["Potential_LOH"].map({True: "LOH", False: "nonLOH"})
    plot_df = to[["hp_imbalance", "LOH_str"]].dropna()
    sns.violinplot(
        data=plot_df, x="LOH_str", y="hp_imbalance",
        palette={"LOH": "#E65100", "nonLOH": "#1565C0"},
        ax=ax, inner="quartile", cut=0, linewidth=0.8,
    )
    ax.set_title("HP Imbalance by LOH Status (TO)", fontweight="bold")
    ax.set_ylabel("|HP_Ratio - 0.5|")

    # Panel 2: imbalance by truth
    ax = axes[0, 1]
    plot_df2 = to[["hp_imbalance", "truth_label"]].dropna()
    sns.violinplot(
        data=plot_df2, x="truth_label", y="hp_imbalance",
        palette=TRUTH_PALETTE, ax=ax, inner="quartile", cut=0, linewidth=0.8,
        order=["TP", "FP"],
    )
    ax.set_title("HP Imbalance by Truth (TO)", fontweight="bold")
    ax.set_ylabel("|HP_Ratio - 0.5|")

    # Panel 3: 2D density — imbalance vs truth, split by LOH
    ax = axes[1, 0]
    for loh_val, ls, alpha in [(True, "-", 0.8), (False, "--", 0.5)]:
        for truth, color in TRUTH_PALETTE.items():
            mask = (to["Potential_LOH"] == loh_val) & (to["truth_label"] == truth)
            vals = to.loc[mask, "hp_imbalance"].dropna()
            if len(vals) > 50:
                label = f"{'LOH' if loh_val else 'nonLOH'}/{truth} (n={len(vals):,})"
                sns.kdeplot(vals, ax=ax, color=color, linestyle=ls, label=label, linewidth=1.5)
    ax.set_xlabel("|HP_Ratio - 0.5|")
    ax.set_title("HP Imbalance Density: LOH x Truth (TO)", fontweight="bold")
    ax.legend(fontsize=7, loc="upper right")

    # Panel 4: imbalance AUC by LOH group
    ax = axes[1, 1]
    ax.axis("off")
    auc_rows = []
    y_bin = encode_truth_binary(to)
    for loh_val, loh_label in [(True, "LOH"), (False, "nonLOH"), (None, "All")]:
        if loh_val is not None:
            mask = to["Potential_LOH"] == loh_val
        else:
            mask = pd.Series(True, index=to.index)
        sub_y = y_bin[mask].dropna()
        sub_x = to.loc[mask, "hp_imbalance"]
        valid = sub_y.notna() & sub_x.notna()
        if valid.sum() > 10:
            auc_val = _auc_safe(sub_y[valid].values, sub_x[valid].values)
            n_tp = int((sub_y[valid] == 1).sum())
            n_fp = int((sub_y[valid] == 0).sum())
        else:
            auc_val = float("nan")
            n_tp = n_fp = 0
        auc_rows.append({
            "Group": loh_label,
            "AUC": f"{auc_val:.3f}" if not np.isnan(auc_val) else "N/A",
            "N_TP": f"{n_tp:,}",
            "N_FP": f"{n_fp:,}",
        })
    tbl = ax.table(
        cellText=[[r[k] for k in ["Group", "AUC", "N_TP", "N_FP"]] for r in auc_rows],
        colLabels=["LOH Group", "AUC (imbalance)", "N_TP", "N_FP"],
        loc="center", cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.2, 1.6)
    ax.set_title("HP Imbalance AUC for TP/FP Discrimination", fontweight="bold", pad=20)

    fig.suptitle("Fig 08: HP Imbalance vs LOH (TO Only)", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O7 Fig08")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig08_to_hp_imbalance_vs_loh.png")

    # Stats
    loh_vals = to.loc[to["Potential_LOH"] == True, "hp_imbalance"].dropna().values
    non_loh_vals = to.loc[to["Potential_LOH"] == False, "hp_imbalance"].dropna().values
    mw = mannwhitney_test(loh_vals, non_loh_vals)
    es = compute_effect_size(loh_vals, non_loh_vals)
    return {"mannwhitney_loh_vs_nonloh": mw, "effect_size_loh_vs_nonloh": es,
            "loh_mean_imbalance": float(np.nanmean(loh_vals)),
            "nonloh_mean_imbalance": float(np.nanmean(non_loh_vals)),
            "auc_rows": auc_rows}


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    setup_plot_style()
    out = OUT_DIR
    print(f"[O7] Output directory: {out}")

    # Load data
    df = load_master_dataset()
    print(f"[O7] Dataset: {len(df):,} rows, {len(df.columns)} columns")
    print(f"[O7] Columns available: {sorted(df.columns.tolist())[:20]}...")
    print(f"[O7] Mode distribution: {df['mode'].value_counts().to_dict()}")

    # Write round context
    write_round_context(
        out_dir=out,
        observation_id="O07",
        title="TO Phasing Quality Observation",
        description=(
            "Systematic observation of tumour-only (TO) phasing quality post HP-fix. "
            "Covers HP tag composition, assign rate, phase block structure, HP0 fraction, "
            "TO/paired concordance, tier distribution, chromosome patterns, and HP imbalance vs LOH."
        ),
        script_path="scripts/analysis/build_observation_O07_to_phasing_quality.py",
        data_sources=[{
            "name": "master_dataset",
            "path": str(load_master_dataset.__defaults__),
            "description": "748K-row master dataset post HP-fix",
        }],
        row_count=len(df),
        col_count=len(df.columns),
    )

    # Write feature statistics
    available = _available_cols(df, HP_PHASING_COLS)
    write_feature_statistics(df[df["mode"] == "to"], available, out / "feature_statistics_to.tsv")
    write_feature_statistics(df, available, out / "feature_statistics_all.tsv")

    # ── Figures ──────────────────────────────────────────────────────────
    results = {}

    print("\n[O7] === Fig 01: TO HP Tag Composition ===")
    fig01_to_hp_tag_composition(df, out)

    print("\n[O7] === Fig 02: hp_assign_rate Distribution ===")
    results["fig02"] = fig02_to_hp_assign_rate_distribution(df, out)

    print("\n[O7] === Fig 03: Phase Block Size Distribution ===")
    fig03_to_phase_block_size_distribution(df, out)

    print("\n[O7] === Fig 04: HP0 Fraction vs Truth ===")
    results["fig04"] = fig04_to_hp0_fraction_vs_truth(df, out)

    print("\n[O7] === Fig 05: TO/Paired HP Concordance ===")
    results["fig05"] = fig05_to_paired_hp_concordance(df, out)

    print("\n[O7] === Fig 06: Tier Distribution Post Fix ===")
    results["fig06"] = fig06_to_tier_distribution_post_fix(df, out)

    print("\n[O7] === Fig 07: Phase Quality by Chromosome ===")
    chrom_stats = fig07_to_phase_quality_by_chromosome(df, out)

    print("\n[O7] === Fig 08: HP Imbalance vs LOH ===")
    results["fig08"] = fig08_to_hp_imbalance_vs_loh(df, out)

    # ── Save results JSON ────────────────────────────────────────────────
    # Serialize for JSON
    def _serialize(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pd.DataFrame):
            return obj.to_dict(orient="index")
        if isinstance(obj, pd.Series):
            return obj.to_dict()
        return str(obj)

    results_path = out / "results_summary.json"
    with results_path.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=_serialize)
    print(f"\n[O7] Results saved to {results_path}")

    print(f"\n[O7] === All 8 figures generated in {out} ===")
    return results


if __name__ == "__main__":
    results = main()
