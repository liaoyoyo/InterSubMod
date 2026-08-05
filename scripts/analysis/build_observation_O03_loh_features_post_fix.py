#!/usr/bin/env python3
"""O3: LOH Feature Global Observation (Post HP Fix).

Re-establishes LOH feature landscape on the post-HP-fix master dataset.
Generates 10 figures + data files covering HP ratio distributions, LOH
enrichment patterns, chromosome hotspots, concordance, and feature importance.

Usage:
    python3 build_observation_O03_loh_features_post_fix.py
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

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import select_current_view, select_loh_legacy

# ── Constants ────────────────────────────────────────────────────────────
OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260401_O03_loh_features_post_fix")
CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz "
    "(20260327_loh_round1_cross_sample_audit, post-HP-fix) | "
)

LOH_FEATURE_COLS = [
    "HP_Ratio", "hp_ratio_core", "effective_hp_reads",
    "core_loh_like", "Potential_LOH", "LOH_Subtype_LegacyVC",
    "HP1FamilyN", "HP2FamilyN", "NHP0", "NHP3",
    "hp0_ratio", "hp3_ratio", "hp_assign_rate",
    "tool_core_loh_match", "tool_potential_loh",
    "to_loh_bed_hit", "HPMergedSig",
]

LOH_NUMERIC_FEATURES = [
    "HP_Ratio", "hp_ratio_core", "effective_hp_reads",
    "HP1FamilyN", "HP2FamilyN", "NHP0", "NHP3",
    "hp0_ratio", "hp3_ratio", "hp_assign_rate",
]


# ══════════════════════════════════════════════════════════════════════════
# Helper utilities
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


def _available_cols(df, cols):
    """Return list of columns that exist in df."""
    return [c for c in cols if c in df.columns]


def attach_verification_loh_views(df: pd.DataFrame) -> pd.DataFrame:
    """Attach guarded C2 current and canonical legacy-derived LOH views."""
    current = select_current_view(df)
    loh = select_loh_legacy(df)
    prepared = df.copy()
    prepared["_verification_class_current"] = current.values
    prepared["_loh_subtype_legacy"] = loh.values
    prepared.attrs["schema_metadata"] = {
        "verification_selection_field": current.field,
        "verification_schema_status": current.schema_status,
        "unknown_current_class_count": sum(current.unknown_counts.values()),
        "unknown_current_class_values": current.unknown_counts,
        "loh_selection_field": loh.field,
        "loh_schema_status": loh.schema_status,
    }
    return prepared


def write_schema_provenance(df: pd.DataFrame, out: Path) -> None:
    metadata = df.attrs["schema_metadata"]
    pd.DataFrame(
        [
            {
                **{key: value for key, value in metadata.items() if key != "unknown_current_class_values"},
                "unknown_current_class_values": "; ".join(
                    f"{key}:{value}"
                    for key, value in sorted(metadata["unknown_current_class_values"].items())
                ),
            }
        ]
    ).to_csv(out / "schema_provenance.tsv", sep="\t", index=False)


# ══════════════════════════════════════════════════════════════════════════
# Figure 01: HP Ratio Distribution by Sample/Mode
# ══════════════════════════════════════════════════════════════════════════

def fig01_hp_ratio_distribution_by_sample_mode(df, out):
    """HP_Ratio violin/box plots across sample x mode, colored by truth."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Left panel: by sample
    ax = axes[0]
    plot_df = df[["HP_Ratio", "sample", "truth_label"]].dropna()
    # Order samples
    sample_order = [s for s in SAMPLE_ORDER if s in plot_df["sample"].unique()]
    plot_df["sample"] = pd.Categorical(plot_df["sample"], categories=sample_order, ordered=True)
    sns.boxplot(
        data=plot_df, x="sample", y="HP_Ratio", hue="truth_label",
        palette=TRUTH_PALETTE, ax=ax, showfliers=False, linewidth=0.8,
    )
    ax.set_title("HP_Ratio by Sample & Truth Label", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("HP_Ratio")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="Truth", loc="upper right", fontsize=8)

    # Right panel: by mode
    ax = axes[1]
    plot_df2 = df[["HP_Ratio", "mode", "truth_label"]].dropna()
    plot_df2["mode"] = pd.Categorical(plot_df2["mode"], categories=MODE_ORDER, ordered=True)
    sns.violinplot(
        data=plot_df2, x="mode", y="HP_Ratio", hue="truth_label",
        palette=TRUTH_PALETTE, ax=ax, split=True, inner="quartile",
        linewidth=0.8, cut=0,
    )
    ax.set_title("HP_Ratio by Mode & Truth Label", fontweight="bold")
    ax.set_xlabel("Mode")
    ax.set_ylabel("HP_Ratio")
    ax.legend(title="Truth", loc="upper right", fontsize=8)

    fig.suptitle("Fig 01: HP Ratio Distribution by Sample & Mode", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig01")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig01_hp_ratio_distribution_by_sample_mode.png")


# ══════════════════════════════════════════════════════════════════════════
# Figure 02: Effective HP Reads Distribution
# ══════════════════════════════════════════════════════════════════════════

def fig02_effective_hp_reads_distribution(df, out):
    """Effective HP reads distribution: paired vs TO, log scale."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: overall histogram by mode
    ax = axes[0]
    for mode, color in zip(MODE_ORDER, [COLOR_PAIRED, COLOR_TO]):
        vals = df.loc[df["mode"] == mode, "effective_hp_reads"].dropna()
        vals_log = np.log10(vals.clip(lower=1))
        ax.hist(vals_log, bins=80, alpha=0.55, label=f"{mode} (n={len(vals):,})",
                color=color, edgecolor="none")
    ax.set_xlabel("log10(effective_hp_reads)")
    ax.set_ylabel("Count")
    ax.set_title("Effective HP Reads (log10)", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel 2: by mode and truth
    ax = axes[1]
    for mode, ls in zip(MODE_ORDER, ["-", "--"]):
        for truth, color in TRUTH_PALETTE.items():
            mask = (df["mode"] == mode) & (df["truth_label"] == truth)
            vals = df.loc[mask, "effective_hp_reads"].dropna()
            if len(vals) < 10:
                continue
            vals_log = np.log10(vals.clip(lower=1))
            ax.hist(vals_log, bins=60, alpha=0.35, label=f"{mode}/{truth}",
                    color=color, edgecolor="none",
                    linestyle=ls)
    ax.set_xlabel("log10(effective_hp_reads)")
    ax.set_ylabel("Count")
    ax.set_title("By Mode x Truth Label", fontweight="bold")
    ax.legend(fontsize=7)

    # Panel 3: summary stats table
    ax = axes[2]
    ax.axis("off")
    stats_rows = []
    for mode in MODE_ORDER:
        for truth in ["TP", "FP"]:
            mask = (df["mode"] == mode) & (df["truth_label"] == truth)
            vals = df.loc[mask, "effective_hp_reads"].dropna()
            stats_rows.append({
                "Mode": mode, "Truth": truth,
                "N": f"{len(vals):,}",
                "Mean": f"{vals.mean():.1f}",
                "Median": f"{vals.median():.0f}",
                "P25": f"{vals.quantile(0.25):.0f}",
                "P75": f"{vals.quantile(0.75):.0f}",
            })
    tbl = ax.table(
        cellText=[[r[k] for k in ["Mode", "Truth", "N", "Mean", "Median", "P25", "P75"]] for r in stats_rows],
        colLabels=["Mode", "Truth", "N", "Mean", "Median", "P25", "P75"],
        loc="center", cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.4)
    ax.set_title("Summary Statistics", fontweight="bold")

    fig.suptitle("Fig 02: Effective HP Reads Distribution (Paired vs TO)", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig02")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig02_effective_hp_reads_distribution.png")


# ══════════════════════════════════════════════════════════════════════════
# Figure 03: core_loh_like Rate by Tier
# ══════════════════════════════════════════════════════════════════════════

def fig03_core_loh_like_rate_by_tier(df, out):
    """LOH rate by Quality_Tier, mode, and truth_label."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    tier_order = ["High", "Medium", "Low"]
    loh_col = "Potential_LOH"  # same as core_loh_like in post-fix data

    # Panel 1: LOH rate by tier x mode
    ax = axes[0]
    rates = df.groupby(["Quality_Tier", "mode"])[loh_col].mean().unstack("mode")
    rates = rates.reindex(tier_order)
    rates.plot(kind="bar", ax=ax, color=[COLOR_PAIRED, COLOR_TO], edgecolor="white", width=0.7)
    ax.set_title("LOH Rate by Tier & Mode", fontweight="bold")
    ax.set_ylabel("LOH Rate")
    ax.set_xlabel("Quality Tier")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title="Mode", fontsize=8)
    for container in ax.containers:
        ax.bar_label(container, fmt="%.3f", fontsize=7, padding=2)

    # Panel 2: LOH rate by tier x truth
    ax = axes[1]
    rates2 = df.groupby(["Quality_Tier", "truth_label"])[loh_col].mean().unstack("truth_label")
    rates2 = rates2.reindex(tier_order)
    rates2.plot(kind="bar", ax=ax, color=[COLOR_FP, COLOR_TP], edgecolor="white", width=0.7)
    ax.set_title("LOH Rate by Tier & Truth", fontweight="bold")
    ax.set_ylabel("LOH Rate")
    ax.set_xlabel("Quality Tier")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title="Truth", fontsize=8)
    for container in ax.containers:
        ax.bar_label(container, fmt="%.3f", fontsize=7, padding=2)

    # Panel 3: LOH subtype breakdown by tier
    ax = axes[2]
    sub_df = df[["Quality_Tier", "_loh_subtype_legacy"]].copy()
    sub_cts = sub_df.groupby(["Quality_Tier", "_loh_subtype_legacy"]).size().unstack("_loh_subtype_legacy", fill_value=0)
    sub_cts = sub_cts.reindex(tier_order)
    # Normalize to proportions
    sub_pcts = sub_cts.div(sub_cts.sum(axis=1), axis=0)
    subtype_order = ["None", "LOH_Noise", "LOH_Weak", "LOH_Strong", "LOH_Subclone"]
    sub_pcts = sub_pcts.reindex(columns=[c for c in subtype_order if c in sub_pcts.columns])
    sub_pcts.plot(kind="bar", stacked=True, ax=ax, edgecolor="white", width=0.7,
                  colormap="Set2")
    ax.set_title("LOH Subtype Proportion by Tier", fontweight="bold")
    ax.set_ylabel("Proportion")
    ax.set_xlabel("Quality Tier")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title="LOH Subtype", fontsize=7, loc="upper right")

    fig.suptitle("Fig 03: LOH Rate by Quality Tier", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig03")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig03_core_loh_like_rate_by_tier.png")


# ══════════════════════════════════════════════════════════════════════════
# Figure 04: LOH Enrichment Heatmap (sample x mode)
# ══════════════════════════════════════════════════════════════════════════

def fig04_loh_enrichment_heatmap(df, out):
    """LOH enrichment (FP_rate / TP_rate) heatmap across sample x mode."""
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))

    sample_order = [s for s in SAMPLE_ORDER if s in df["sample"].unique()]
    loh_col = "Potential_LOH"

    # Panel 1: LOH rate heatmap (TP)
    ax = axes[0]
    tp_rates = df[df["truth_label"] == "TP"].groupby(["sample", "mode"])[loh_col].mean().unstack("mode")
    tp_rates = tp_rates.reindex(sample_order)
    sns.heatmap(tp_rates, annot=True, fmt=".3f", cmap="Blues", ax=ax,
                vmin=0, vmax=0.8, linewidths=0.5, cbar_kws={"label": "LOH Rate"})
    ax.set_title("TP LOH Rate", fontweight="bold")
    ax.set_ylabel("Sample")
    ax.set_xlabel("Mode")

    # Panel 2: LOH rate heatmap (FP)
    ax = axes[1]
    fp_rates = df[df["truth_label"] == "FP"].groupby(["sample", "mode"])[loh_col].mean().unstack("mode")
    fp_rates = fp_rates.reindex(sample_order)
    sns.heatmap(fp_rates, annot=True, fmt=".3f", cmap="Reds", ax=ax,
                vmin=0, vmax=0.8, linewidths=0.5, cbar_kws={"label": "LOH Rate"})
    ax.set_title("FP LOH Rate", fontweight="bold")
    ax.set_ylabel("Sample")
    ax.set_xlabel("Mode")

    # Panel 3: Enrichment ratio (FP/TP)
    ax = axes[2]
    enrichment = fp_rates / tp_rates.replace(0, np.nan)
    sns.heatmap(enrichment, annot=True, fmt=".2f", cmap="RdYlGn_r", ax=ax,
                center=1.0, vmin=0.3, vmax=3.0, linewidths=0.5,
                cbar_kws={"label": "Enrichment (FP/TP)"})
    ax.set_title("LOH Enrichment (FP rate / TP rate)", fontweight="bold")
    ax.set_ylabel("Sample")
    ax.set_xlabel("Mode")

    fig.suptitle("Fig 04: LOH Enrichment Heatmap (Sample x Mode)", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig04 | enrichment > 1 = FP-enriched")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig04_loh_enrichment_heatmap.png")

    # Save enrichment table
    enr_out = enrichment.copy()
    enr_out.to_csv(out / "fig04_loh_enrichment_table.tsv", sep="\t", float_format="%.4f")
    return enrichment


# ══════════════════════════════════════════════════════════════════════════
# Figure 05: HP0/HP3 Ratio vs LOH
# ══════════════════════════════════════════════════════════════════════════

def fig05_hp0_hp3_ratio_scatter(df, out):
    """HP0 ratio vs HP3 ratio scatter, colored by LOH status."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))

    plot_df = df[["hp0_ratio", "hp3_ratio", "Potential_LOH", "truth_label", "mode"]].dropna()

    # Panel 1: scatter colored by LOH
    ax = axes[0]
    for loh_val, color, label in [(False, "#66BB6A", "Non-LOH"), (True, "#EF5350", "LOH")]:
        mask = plot_df["Potential_LOH"] == loh_val
        sub = plot_df[mask].sample(min(5000, mask.sum()), random_state=42)
        ax.scatter(sub["hp0_ratio"], sub["hp3_ratio"], alpha=0.15, s=6,
                   color=color, label=f"{label} (n={mask.sum():,})", edgecolors="none")
    ax.set_xlabel("HP0 Ratio (unphased)")
    ax.set_ylabel("HP3 Ratio (multi-phase)")
    ax.set_title("HP0 vs HP3 Ratio by LOH Status", fontweight="bold")
    ax.legend(fontsize=8, markerscale=3)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 0.52)

    # Panel 2: hp0_ratio distribution by LOH x truth
    ax = axes[1]
    for loh_val, ls in [(False, "-"), (True, "--")]:
        for truth, color in TRUTH_PALETTE.items():
            mask = (plot_df["Potential_LOH"] == loh_val) & (plot_df["truth_label"] == truth)
            vals = plot_df.loc[mask, "hp0_ratio"]
            if len(vals) < 10:
                continue
            loh_str = "LOH" if loh_val else "noLOH"
            ax.hist(vals, bins=50, alpha=0.35, color=color,
                    label=f"{loh_str}/{truth} (n={len(vals):,})", edgecolor="none")
    ax.set_xlabel("HP0 Ratio")
    ax.set_ylabel("Count")
    ax.set_title("HP0 Ratio by LOH x Truth", fontweight="bold")
    ax.legend(fontsize=7)

    # Panel 3: hp_assign_rate distribution by LOH
    ax = axes[2]
    for loh_val, color, label in [(False, "#66BB6A", "Non-LOH"), (True, "#EF5350", "LOH")]:
        mask = df["Potential_LOH"] == loh_val
        vals = df.loc[mask, "hp_assign_rate"].dropna()
        ax.hist(vals, bins=50, alpha=0.5, color=color,
                label=f"{label} (n={len(vals):,})", edgecolor="none")
    ax.set_xlabel("HP Assign Rate")
    ax.set_ylabel("Count")
    ax.set_title("HP Assignment Rate by LOH Status", fontweight="bold")
    ax.legend(fontsize=8)

    fig.suptitle("Fig 05: HP0/HP3 Ratio & LOH Relationship", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig05")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig05_hp0_hp3_ratio_scatter.png")


# ══════════════════════════════════════════════════════════════════════════
# Figure 06: TO vs Paired LOH Concordance
# ══════════════════════════════════════════════════════════════════════════

def fig06_to_vs_paired_loh_concordance(df, out):
    """LOH concordance between TO and paired mode at matched variant_keys."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    # Build matched pairs via variant_key
    paired_df = df[df["mode"] == "paired"][["variant_key", "sample", "Potential_LOH", "truth_label", "_loh_subtype_legacy"]].copy()
    paired_df = paired_df.rename(columns={
        "Potential_LOH": "LOH_paired",
        "truth_label": "truth_paired",
        "_loh_subtype_legacy": "subtype_paired",
    })
    to_df = df[df["mode"] == "to"][["variant_key", "sample", "Potential_LOH", "truth_label", "_loh_subtype_legacy"]].copy()
    to_df = to_df.rename(columns={
        "Potential_LOH": "LOH_to",
        "truth_label": "truth_to",
        "_loh_subtype_legacy": "subtype_to",
    })

    matched = paired_df.merge(to_df, on=["variant_key", "sample"], how="inner")
    print(f"[O3 Fig06] Matched variant_keys: {len(matched):,}")

    if len(matched) < 100:
        # Fallback: rate-level comparison
        for ax in axes:
            ax.text(0.5, 0.5, "Insufficient matched pairs", ha="center", va="center")
            ax.axis("off")
        fig.suptitle("Fig 06: TO vs Paired LOH Concordance (insufficient data)", fontsize=14)
        save_figure(fig, out / "fig06_to_vs_paired_loh_concordance.png")
        return

    # Panel 1: 2x2 concordance matrix
    ax = axes[0]
    conc = pd.crosstab(matched["LOH_paired"], matched["LOH_to"],
                       rownames=["Paired LOH"], colnames=["TO LOH"])
    sns.heatmap(conc, annot=True, fmt=",", cmap="YlOrBr", ax=ax,
                linewidths=0.5, cbar_kws={"label": "Count"})
    ax.set_title("LOH Concordance Matrix", fontweight="bold")
    # Compute agreement rate
    agree = ((matched["LOH_paired"] == matched["LOH_to"]).sum())
    agree_rate = agree / len(matched)
    ax.text(0.5, -0.15, f"Agreement: {agree_rate:.1%} ({agree:,}/{len(matched):,})",
            transform=ax.transAxes, ha="center", fontsize=9)

    # Panel 2: Concordance by sample
    ax = axes[1]
    sample_conc = matched.groupby("sample").apply(
        lambda g: (g["LOH_paired"] == g["LOH_to"]).mean()
    ).reindex([s for s in SAMPLE_ORDER if s in matched["sample"].unique()])
    sample_conc.plot(kind="barh", ax=ax, color="#5C6BC0", edgecolor="white")
    ax.set_xlabel("Agreement Rate")
    ax.set_title("LOH Agreement by Sample", fontweight="bold")
    ax.set_xlim(0, 1)
    for i, v in enumerate(sample_conc.values):
        ax.text(v + 0.01, i, f"{v:.3f}", va="center", fontsize=8)

    # Panel 3: Discordant cases analysis
    ax = axes[2]
    discord = matched[matched["LOH_paired"] != matched["LOH_to"]]
    discord_type = pd.DataFrame({
        "Paired=LOH, TO=nonLOH": [((discord["LOH_paired"] == True) & (discord["LOH_to"] == False)).sum()],
        "Paired=nonLOH, TO=LOH": [((discord["LOH_paired"] == False) & (discord["LOH_to"] == True)).sum()],
    }).T
    discord_type.columns = ["Count"]
    discord_type["Pct"] = discord_type["Count"] / len(matched) * 100

    # By truth label in discordant set
    discord_truth = discord.groupby(["truth_paired", "truth_to"]).size().reset_index(name="count")

    ax.axis("off")
    text_lines = [
        f"Total matched: {len(matched):,}",
        f"Concordant: {agree:,} ({agree_rate:.1%})",
        f"Discordant: {len(discord):,} ({1 - agree_rate:.1%})",
        "",
        "Discordant breakdown:",
    ]
    for idx, row in discord_type.iterrows():
        text_lines.append(f"  {idx}: {int(row['Count']):,} ({row['Pct']:.1f}%)")
    text_lines.append("")
    text_lines.append("Discordant truth-label pairs:")
    for _, row in discord_truth.iterrows():
        text_lines.append(f"  paired={row['truth_paired']}, to={row['truth_to']}: {row['count']:,}")

    ax.text(0.05, 0.95, "\n".join(text_lines), transform=ax.transAxes,
            fontsize=9, va="top", ha="left", family="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#F5F5F5", alpha=0.8))
    ax.set_title("Discordance Summary", fontweight="bold")

    fig.suptitle("Fig 06: TO vs Paired LOH Concordance (Matched Variants)", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig06")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig06_to_vs_paired_loh_concordance.png")

    # Save concordance data
    conc.to_csv(out / "fig06_concordance_matrix.tsv", sep="\t")
    return matched


# ══════════════════════════════════════════════════════════════════════════
# Figure 07: LOH by Chromosome
# ══════════════════════════════════════════════════════════════════════════

def fig07_loh_by_chromosome(df, out):
    """LOH distribution across chromosomes, identifying hotspots."""
    fig, axes = plt.subplots(2, 2, figsize=(20, 14))

    chrom_col = "Chr"
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chrom_order = [c for c in chrom_order if c in df[chrom_col].unique()]

    # Panel 1: LOH count by chromosome
    ax = axes[0, 0]
    loh_df = df[df["Potential_LOH"] == True]
    chrom_cts = loh_df[chrom_col].value_counts().reindex(chrom_order, fill_value=0)
    total_cts = df[chrom_col].value_counts().reindex(chrom_order, fill_value=0)
    ax.bar(range(len(chrom_order)), chrom_cts.values, color="#EF5350", edgecolor="white", alpha=0.7,
           label="LOH")
    ax.bar(range(len(chrom_order)), total_cts.values - chrom_cts.values, bottom=chrom_cts.values,
           color="#90CAF9", edgecolor="white", alpha=0.5, label="Non-LOH")
    ax.set_xticks(range(len(chrom_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order], fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Count")
    ax.set_title("LOH Count by Chromosome", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel 2: LOH rate by chromosome
    ax = axes[0, 1]
    loh_rate = (chrom_cts / total_cts).fillna(0)
    colors = ["#EF5350" if r > loh_rate.mean() + 2 * loh_rate.std() else "#78909C" for r in loh_rate.values]
    ax.bar(range(len(chrom_order)), loh_rate.values, color=colors, edgecolor="white")
    ax.axhline(loh_rate.mean(), color="black", linestyle="--", linewidth=0.8, label=f"Mean: {loh_rate.mean():.3f}")
    ax.axhline(loh_rate.mean() + 2 * loh_rate.std(), color="red", linestyle=":", linewidth=0.8,
               label=f"+2SD: {loh_rate.mean() + 2 * loh_rate.std():.3f}")
    ax.set_xticks(range(len(chrom_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order], fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("LOH Rate")
    ax.set_title("LOH Rate by Chromosome (red = >2SD)", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel 3: LOH rate by chromosome x mode
    ax = axes[1, 0]
    loh_rate_by_mode = df.groupby([chrom_col, "mode"])["Potential_LOH"].mean().unstack("mode")
    loh_rate_by_mode = loh_rate_by_mode.reindex(chrom_order)
    x = np.arange(len(chrom_order))
    width = 0.35
    if "paired" in loh_rate_by_mode.columns:
        ax.bar(x - width / 2, loh_rate_by_mode["paired"].values, width, color=COLOR_PAIRED,
               label="paired", edgecolor="white", alpha=0.8)
    if "to" in loh_rate_by_mode.columns:
        ax.bar(x + width / 2, loh_rate_by_mode["to"].values, width, color=COLOR_TO,
               label="to", edgecolor="white", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order], fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("LOH Rate")
    ax.set_title("LOH Rate by Chromosome & Mode", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel 4: LOH rate by chromosome x sample (heatmap)
    ax = axes[1, 1]
    loh_rate_samp = df.groupby([chrom_col, "sample"])["Potential_LOH"].mean().unstack("sample")
    sample_order = [s for s in SAMPLE_ORDER if s in loh_rate_samp.columns]
    loh_rate_samp = loh_rate_samp.reindex(index=chrom_order, columns=sample_order)
    sns.heatmap(loh_rate_samp, cmap="YlOrRd", ax=ax, linewidths=0.3,
                cbar_kws={"label": "LOH Rate"}, vmin=0, vmax=0.8)
    ax.set_title("LOH Rate: Chromosome x Sample", fontweight="bold")
    ax.set_ylabel("Chromosome")
    ax.set_xlabel("Sample")

    fig.suptitle("Fig 07: LOH Distribution by Chromosome", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig07")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig07_loh_by_chromosome.png")

    # Save chromosome LOH rates
    loh_rate.to_frame("loh_rate").to_csv(out / "fig07_chrom_loh_rates.tsv", sep="\t", float_format="%.4f")
    return loh_rate


# ══════════════════════════════════════════════════════════════════════════
# Figure 08: LOH vs AF Scatter
# ══════════════════════════════════════════════════════════════════════════

def fig08_loh_vs_af_scatter(df, out):
    """LOH status vs caller AF scatter/density plot."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))

    plot_df = df[["caller_af", "Potential_LOH", "truth_label", "mode"]].dropna()

    # Panel 1: AF distribution by LOH status
    ax = axes[0]
    for loh_val, color, label in [(False, "#66BB6A", "Non-LOH"), (True, "#EF5350", "LOH")]:
        vals = plot_df.loc[plot_df["Potential_LOH"] == loh_val, "caller_af"]
        ax.hist(vals, bins=80, alpha=0.5, color=color, label=f"{label} (n={len(vals):,})",
                edgecolor="none", density=True)
    ax.set_xlabel("Caller AF")
    ax.set_ylabel("Density")
    ax.set_title("AF Distribution by LOH Status", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel 2: LOH rate vs AF bins
    ax = axes[1]
    af_bins = pd.cut(plot_df["caller_af"], bins=np.arange(0, 1.05, 0.05))
    loh_by_af = plot_df.groupby(af_bins)["Potential_LOH"].agg(["mean", "count"])
    # Only plot bins with >100 data points
    loh_by_af = loh_by_af[loh_by_af["count"] >= 100]
    bin_centers = [interval.mid for interval in loh_by_af.index]
    ax.plot(bin_centers, loh_by_af["mean"], "o-", color="#EF5350", linewidth=1.5, markersize=4)
    ax.fill_between(bin_centers, loh_by_af["mean"], alpha=0.15, color="#EF5350")
    ax.set_xlabel("Caller AF (bin center)")
    ax.set_ylabel("LOH Rate")
    ax.set_title("LOH Rate vs Caller AF", fontweight="bold")
    ax.set_xlim(0, 1.0)

    # Add count annotation as secondary y-axis
    ax2 = ax.twinx()
    ax2.bar(bin_centers, loh_by_af["count"], width=0.04, alpha=0.15, color="#78909C")
    ax2.set_ylabel("Count (bars)", color="#78909C")
    ax2.tick_params(axis="y", labelcolor="#78909C")

    # Panel 3: AF >= 0.9 (high-AF LOH) focus
    ax = axes[2]
    high_af = plot_df[plot_df["caller_af"] >= 0.5]
    for truth, color in TRUTH_PALETTE.items():
        for loh_val, ls, marker in [(True, "-", "o"), (False, "--", "s")]:
            mask = (high_af["truth_label"] == truth) & (high_af["Potential_LOH"] == loh_val)
            vals = high_af.loc[mask, "caller_af"]
            loh_str = "LOH" if loh_val else "noLOH"
            ax.hist(vals, bins=40, alpha=0.3, color=color,
                    label=f"{truth}/{loh_str} (n={len(vals):,})", edgecolor="none")
    ax.set_xlabel("Caller AF")
    ax.set_ylabel("Count")
    ax.set_title("High-AF Region (AF >= 0.5) by Truth x LOH", fontweight="bold")
    ax.legend(fontsize=7)

    fig.suptitle("Fig 08: LOH Status vs Caller AF", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig08")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig08_loh_vs_af_scatter.png")


# ══════════════════════════════════════════════════════════════════════════
# Figure 09: Effective HP vs VerificationClass
# ══════════════════════════════════════════════════════════════════════════

def fig09_effective_hp_vs_verificationclass(df, out):
    """Effective HP reads bins x VerificationClass count heatmap."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Define effective_hp_reads bins
    hp_bins = [0, 10, 30, 50, 80, 120, 200, 5000]
    hp_labels = ["0-10", "10-30", "30-50", "50-80", "80-120", "120-200", "200+"]
    df_plot = df[["effective_hp_reads", "_verification_class_current", "mode"]].dropna().copy()
    df_plot["hp_bin"] = pd.cut(df_plot["effective_hp_reads"], bins=hp_bins, labels=hp_labels, right=False)

    vc_order = df_plot["_verification_class_current"].value_counts().index.tolist()

    # Panel 1: Heatmap (all data)
    ax = axes[0]
    ct = pd.crosstab(df_plot["hp_bin"], df_plot["_verification_class_current"])
    ct = ct.reindex(columns=[c for c in vc_order if c in ct.columns])
    sns.heatmap(ct, annot=True, fmt=",", cmap="YlGnBu", ax=ax, linewidths=0.3,
                cbar_kws={"label": "Count"})
    ax.set_title("Effective HP Reads x VerificationClass (All)", fontweight="bold")
    ax.set_ylabel("Effective HP Reads Bin")
    ax.set_xlabel("VerificationClass")

    # Panel 2: Normalized by row (proportion)
    ax = axes[1]
    ct_pct = ct.div(ct.sum(axis=1), axis=0)
    sns.heatmap(ct_pct, annot=True, fmt=".2f", cmap="YlGnBu", ax=ax, linewidths=0.3,
                vmin=0, vmax=1, cbar_kws={"label": "Row Proportion"})
    ax.set_title("Row-Normalized (Proportion)", fontweight="bold")
    ax.set_ylabel("Effective HP Reads Bin")
    ax.set_xlabel("VerificationClass")

    fig.suptitle("Fig 09: Effective HP Reads vs VerificationClass", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig09")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig09_effective_hp_vs_verificationclass.png")

    # Save crosstab
    ct.to_csv(out / "fig09_hp_bin_vs_vc_crosstab.tsv", sep="\t")


# ══════════════════════════════════════════════════════════════════════════
# Figure 10: LOH Feature Importance for TP/FP (AUC Ranking)
# ══════════════════════════════════════════════════════════════════════════

def fig10_loh_feature_importance_for_tp_fp(df, out):
    """AUC ranking of LOH-related features for TP/FP discrimination."""
    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    y_bin = encode_truth_binary(df)

    # Compute AUC for each numeric LOH feature
    all_features = LOH_NUMERIC_FEATURES + [
        "Potential_LOH", "core_loh_like", "HPMergedSig",
        "tool_core_loh_match", "tool_potential_loh", "to_loh_bed_hit",
    ]
    avail_features = _available_cols(df, all_features)

    # Overall AUC
    auc_results = []
    for feat in avail_features:
        vals = df[feat].copy()
        if vals.dtype == bool:
            vals = vals.astype(float)
        elif vals.dtype == object:
            continue
        auc_val = _auc_safe(y_bin, vals)
        # Also compute by mode
        auc_paired = _auc_safe(
            y_bin[df["mode"] == "paired"],
            vals[df["mode"] == "paired"],
        )
        auc_to = _auc_safe(
            y_bin[df["mode"] == "to"],
            vals[df["mode"] == "to"],
        )
        auc_results.append({
            "feature": feat,
            "AUC_all": auc_val,
            "AUC_paired": auc_paired,
            "AUC_to": auc_to,
            "direction": "higher=TP" if auc_val >= 0.5 else "higher=FP",
            "AUC_effective": max(auc_val, 1 - auc_val) if not np.isnan(auc_val) else float("nan"),
        })

    auc_df = pd.DataFrame(auc_results).sort_values("AUC_effective", ascending=False)

    # Panel 1: AUC bar chart (overall)
    ax = axes[0]
    colors = ["#EF5350" if d == "higher=FP" else COLOR_TP for d in auc_df["direction"]]
    y_pos = range(len(auc_df))
    ax.barh(y_pos, auc_df["AUC_effective"].values, color=colors, edgecolor="white")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(auc_df["feature"].values, fontsize=8)
    ax.set_xlabel("AUC (effective = max(AUC, 1-AUC))")
    ax.set_title("LOH Feature AUC for TP/FP", fontweight="bold")
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlim(0.45, 0.75)
    for i, (auc_val, dirn) in enumerate(zip(auc_df["AUC_effective"], auc_df["direction"])):
        if not np.isnan(auc_val):
            ax.text(auc_val + 0.002, i, f"{auc_val:.3f} ({dirn})", va="center", fontsize=7)
    ax.invert_yaxis()

    # Panel 2: AUC by mode comparison
    ax = axes[1]
    x = np.arange(len(auc_df))
    width = 0.35
    ax.barh(x - width / 2, auc_df["AUC_paired"].values, width, color=COLOR_PAIRED,
            label="paired", edgecolor="white")
    ax.barh(x + width / 2, auc_df["AUC_to"].values, width, color=COLOR_TO,
            label="TO", edgecolor="white")
    ax.set_yticks(x)
    ax.set_yticklabels(auc_df["feature"].values, fontsize=8)
    ax.set_xlabel("AUC (raw)")
    ax.set_title("AUC by Mode", fontweight="bold")
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.legend(fontsize=8)
    ax.invert_yaxis()

    # Panel 3: Effect size (Cohen's d) for top features
    ax = axes[2]
    effect_results = []
    for feat in auc_df["feature"].values[:10]:
        vals = df[feat].copy()
        if vals.dtype == bool:
            vals = vals.astype(float)
        elif vals.dtype == object:
            continue
        tp_vals = vals[df["truth_label"] == "TP"].dropna().values
        fp_vals = vals[df["truth_label"] == "FP"].dropna().values
        es = compute_effect_size(tp_vals, fp_vals)
        effect_results.append({"feature": feat, **es})

    es_df = pd.DataFrame(effect_results)
    if len(es_df) > 0:
        colors_es = [COLOR_TP if d >= 0 else COLOR_FP for d in es_df["d"]]
        ax.barh(range(len(es_df)), es_df["d"].values, color=colors_es, edgecolor="white")
        ax.errorbar(es_df["d"].values, range(len(es_df)),
                    xerr=[es_df["d"].values - es_df["ci_lo"].values,
                          es_df["ci_hi"].values - es_df["d"].values],
                    fmt="none", ecolor="gray", capsize=2, linewidth=0.8)
        ax.set_yticks(range(len(es_df)))
        ax.set_yticklabels(es_df["feature"].values, fontsize=8)
        ax.set_xlabel("Cohen's d (TP - FP)")
        ax.set_title("Effect Size (Cohen's d) with 95% CI", fontweight="bold")
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)
        ax.invert_yaxis()

    fig.suptitle("Fig 10: LOH Feature Importance for TP/FP Discrimination", fontsize=14, fontweight="bold")
    add_caption(fig, CAPTION_SRC + "O3 Fig10 | Blue=higher in TP, Red=higher in FP")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    save_figure(fig, out / "fig10_loh_feature_importance_for_tp_fp.png")

    # Save AUC table
    auc_df.to_csv(out / "fig10_loh_feature_auc_table.tsv", sep="\t", index=False, float_format="%.4f")
    return auc_df


# ══════════════════════════════════════════════════════════════════════════
# Statistical summary computation
# ══════════════════════════════════════════════════════════════════════════

def compute_statistical_summary(df, out):
    """Compute comprehensive statistical summaries for LOH features."""
    results = {}

    y_bin = encode_truth_binary(df)

    # 1. LOH rates by mode x truth
    print("\n[O3] Computing LOH rate statistics...")
    for mode in MODE_ORDER:
        for truth in ["TP", "FP"]:
            mask = (df["mode"] == mode) & (df["truth_label"] == truth)
            n = mask.sum()
            loh_n = df.loc[mask, "Potential_LOH"].sum()
            rate = safe_div(loh_n, n)
            key = f"loh_rate_{mode}_{truth}"
            results[key] = {"n": int(n), "loh_n": int(loh_n), "rate": float(rate)}

    # 2. Enrichment by mode
    print("[O3] Computing enrichment by mode...")
    for mode in MODE_ORDER:
        mode_df = df[df["mode"] == mode]
        tp_n = (mode_df["truth_label"] == "TP").sum()
        fp_n = (mode_df["truth_label"] == "FP").sum()
        tp_loh = ((mode_df["truth_label"] == "TP") & (mode_df["Potential_LOH"])).sum()
        fp_loh = ((mode_df["truth_label"] == "FP") & (mode_df["Potential_LOH"])).sum()
        enr = compute_enrichment(int(tp_loh), int(fp_loh), int(tp_n), int(fp_n))
        results[f"enrichment_{mode}"] = enr

    # 3. Mann-Whitney tests for HP features
    print("[O3] Computing Mann-Whitney tests...")
    for feat in LOH_NUMERIC_FEATURES:
        if feat not in df.columns:
            continue
        tp_vals = df.loc[df["truth_label"] == "TP", feat].dropna().values
        fp_vals = df.loc[df["truth_label"] == "FP", feat].dropna().values
        mw = mannwhitney_test(tp_vals, fp_vals)
        es = compute_effect_size(tp_vals, fp_vals)
        results[f"mw_{feat}"] = {**mw, "cohens_d": es["d"], "d_ci": format_ci(es["ci_lo"], es["ci_hi"])}

    # 4. Chi-squared for LOH x truth
    print("[O3] Computing chi-squared tests...")
    ct = pd.crosstab(df["Potential_LOH"], df["truth_label"])
    chi2_res = chi2_test(ct.values)
    results["chi2_loh_truth"] = chi2_res

    # By mode
    for mode in MODE_ORDER:
        mode_df = df[df["mode"] == mode]
        ct_m = pd.crosstab(mode_df["Potential_LOH"], mode_df["truth_label"])
        chi2_m = chi2_test(ct_m.values)
        results[f"chi2_loh_truth_{mode}"] = chi2_m

    # 5. LOH Subtype x Truth chi-squared
    ct_sub = pd.crosstab(df["_loh_subtype_legacy"], df["truth_label"])
    chi2_sub = chi2_test(ct_sub.values)
    results["chi2_loh_subtype_truth"] = chi2_sub

    # Save results
    with open(out / "statistical_summary.json", "w") as f:
        json.dump(results, f, indent=2, default=str, ensure_ascii=False)
    print(f"[O3] Statistical summary saved: {out / 'statistical_summary.json'}")

    return results


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    setup_plot_style()

    print("=" * 70)
    print("O3: LOH Feature Global Observation (Post HP Fix)")
    print("=" * 70)

    # Load data
    df = attach_verification_loh_views(load_master_dataset())
    write_schema_provenance(df, OUT_DIR)
    print(f"\nDataset shape: {df.shape}")
    print(f"Columns: {len(df.columns)}")

    # Verify LOH columns exist
    print("\n[O3] Checking LOH feature columns...")
    avail = _available_cols(df, LOH_FEATURE_COLS)
    missing = [c for c in LOH_FEATURE_COLS if c not in df.columns]
    print(f"  Available: {len(avail)} / {len(LOH_FEATURE_COLS)}")
    if missing:
        print(f"  Missing: {missing}")

    # Write round context
    write_round_context(
        out_dir=OUT_DIR,
        observation_id="O3",
        title="LOH Feature Global Observation (Post HP Fix)",
        description=(
            "Comprehensive LOH feature landscape analysis on post-HP-fix master dataset. "
            "Covers HP ratio distributions, LOH enrichment, chromosome hotspots, "
            "TO/paired concordance, and feature importance for TP/FP discrimination."
        ),
        script_path="scripts/analysis/build_observation_O03_loh_features_post_fix.py",
        data_sources=[{
            "name": "master_dataset",
            "path": str(load_master_dataset.__defaults__),
            "type": "tsv.gz",
        }],
        row_count=len(df),
        col_count=len(df.columns),
        extra={
            "loh_features_available": avail,
            "loh_features_missing": missing,
            "loh_rate_overall": float(df["Potential_LOH"].mean()),
            "post_hp_fix": True,
        },
    )

    # Write feature statistics
    write_feature_statistics(
        df, avail,
        out_path=OUT_DIR / "loh_feature_statistics.tsv",
    )

    # Generate all 10 figures
    print("\n" + "=" * 70)
    print("Generating figures...")
    print("=" * 70)

    print("\n[Fig01] HP Ratio Distribution by Sample & Mode")
    fig01_hp_ratio_distribution_by_sample_mode(df, OUT_DIR)

    print("\n[Fig02] Effective HP Reads Distribution")
    fig02_effective_hp_reads_distribution(df, OUT_DIR)

    print("\n[Fig03] LOH Rate by Quality Tier")
    fig03_core_loh_like_rate_by_tier(df, OUT_DIR)

    print("\n[Fig04] LOH Enrichment Heatmap")
    enrichment = fig04_loh_enrichment_heatmap(df, OUT_DIR)

    print("\n[Fig05] HP0/HP3 Ratio vs LOH")
    fig05_hp0_hp3_ratio_scatter(df, OUT_DIR)

    print("\n[Fig06] TO vs Paired LOH Concordance")
    matched = fig06_to_vs_paired_loh_concordance(df, OUT_DIR)

    print("\n[Fig07] LOH by Chromosome")
    chrom_loh = fig07_loh_by_chromosome(df, OUT_DIR)

    print("\n[Fig08] LOH vs AF Scatter")
    fig08_loh_vs_af_scatter(df, OUT_DIR)

    print("\n[Fig09] Effective HP vs VerificationClass")
    fig09_effective_hp_vs_verificationclass(df, OUT_DIR)

    print("\n[Fig10] LOH Feature Importance for TP/FP")
    auc_df = fig10_loh_feature_importance_for_tp_fp(df, OUT_DIR)

    # Compute statistical summary
    print("\n" + "=" * 70)
    print("Computing statistical summary...")
    print("=" * 70)
    stats = compute_statistical_summary(df, OUT_DIR)

    print("\n" + "=" * 70)
    print("O3 Observation Complete!")
    print(f"Output directory: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
