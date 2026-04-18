#!/usr/bin/env python3
"""LOH Copy Number Signal Detection Analysis.

Investigates the relationship between LOH and copy number alterations
(copy-loss, copy-neutral, copy-gain) across 7 cancer cell lines.

Key questions:
  1. Is LOH predominantly copy-neutral (UPD) or copy-loss (deletion)?
  2. Do specific samples/chromosomes show clear CN abnormalities?
  3. Does copy number influence TO-mode LOH calls?
  4. Is copy-gain LOH real or artifact?

Generates 5 figures + 1 statistics table.

Usage:
    python3 build_loh_copy_number_analysis.py
"""

from __future__ import annotations

import sys
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis")
from observation_common import (
    load_master_dataset,
    setup_plot_style,
    save_figure,
    ensure_dir,
    SAMPLE_ORDER,
    MODE_ORDER,
    COLOR_TP,
    COLOR_FP,
    COLOR_PAIRED,
    COLOR_TO,
    TRUTH_PALETTE,
    MODE_PALETTE,
)

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/"
    "2026/04/20260401_LOH_weekly_report_draft/images"
)

CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz "
    "(20260327_loh_round1_cross_sample_audit) | "
)

# Chromosome order
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# CN category definitions (relative to sample median)
CN_CATEGORIES = {
    "Deep deletion": (0.0, 0.3),
    "Copy loss": (0.3, 0.7),
    "Copy neutral": (0.7, 1.3),
    "Copy gain": (1.3, 2.0),
    "High amplification": (2.0, float("inf")),
}

CN_ORDER = ["Deep deletion", "Copy loss", "Copy neutral", "Copy gain", "High amplification"]
CN_COLORS = {
    "Deep deletion": "#B71C1C",
    "Copy loss": "#EF5350",
    "Copy neutral": "#78909C",
    "Copy gain": "#42A5F5",
    "High amplification": "#1565C0",
}


# ── Data loading ──────────────────────────────────────────────────────────

def load_data() -> pd.DataFrame:
    """Load master dataset with needed columns."""
    usecols = [
        "RegionID", "Chr", "Pos", "NumReads",
        "HP1FamilyN", "HP2FamilyN", "NHP0", "NHP3",
        "HP_Ratio", "core_loh_like", "effective_hp_reads",
        "truth_label", "mode", "sample",
        "caller_af", "caller_dp",
    ]
    df = load_master_dataset(usecols=usecols)

    # Filter to known chromosomes
    df = df[df["Chr"].isin(CHR_ORDER)].copy()

    # Ensure numeric
    for col in ["NumReads", "HP1FamilyN", "HP2FamilyN", "NHP0", "NHP3",
                "HP_Ratio", "caller_af", "caller_dp"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Parse core_loh_like as boolean
    df["loh"] = df["core_loh_like"].astype(str).str.lower() == "true"

    print(f"[CN analysis] Loaded {len(df):,} rows after chr filter")
    return df


def compute_sample_medians(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample genome-wide median NumReads."""
    medians = (
        df.groupby(["sample", "mode"])["NumReads"]
        .median()
        .reset_index()
        .rename(columns={"NumReads": "sample_median_depth"})
    )
    return medians


def assign_cn_category(df: pd.DataFrame, medians: pd.DataFrame) -> pd.DataFrame:
    """Assign CN category to each region based on depth ratio."""
    df = df.merge(medians, on=["sample", "mode"], how="left")
    df["depth_ratio"] = df["NumReads"] / df["sample_median_depth"]

    def _classify(ratio):
        if pd.isna(ratio):
            return "Copy neutral"
        for cat, (lo, hi) in CN_CATEGORIES.items():
            if lo <= ratio < hi:
                return cat
        return "High amplification"

    df["cn_category"] = df["depth_ratio"].apply(_classify)
    return df


def compute_chr_medians(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample per-mode per-chromosome median depth."""
    chr_med = (
        df.groupby(["sample", "mode", "Chr"])["NumReads"]
        .median()
        .reset_index()
        .rename(columns={"NumReads": "chr_median_depth"})
    )
    return chr_med


# ══════════════════════════════════════════════════════════════════════════
# Figure 1: Per-Sample Per-Chromosome Median Read Depth Heatmap
# ══════════════════════════════════════════════════════════════════════════

def fig01_chr_depth_heatmap(df: pd.DataFrame, medians: pd.DataFrame, out: Path):
    """Heatmap of per-sample per-chromosome median depth (Paired / TO)."""
    chr_med = compute_chr_medians(df)
    chr_med = chr_med.merge(medians, on=["sample", "mode"], how="left")
    chr_med["depth_deviation"] = (
        (chr_med["chr_median_depth"] - chr_med["sample_median_depth"])
        / chr_med["sample_median_depth"]
    )

    fig, axes = plt.subplots(1, 2, figsize=(26, 7), sharey=True)

    for ax_i, mode in enumerate(MODE_ORDER):
        ax = axes[ax_i]
        sub = chr_med[chr_med["mode"] == mode]

        # Pivot: rows = samples, cols = chromosomes
        # Use deviation for color, annotate with raw median
        pivot_dev = sub.pivot_table(
            index="sample", columns="Chr",
            values="depth_deviation", aggfunc="first",
        )
        pivot_raw = sub.pivot_table(
            index="sample", columns="Chr",
            values="chr_median_depth", aggfunc="first",
        )

        # Reorder
        pivot_dev = pivot_dev.reindex(
            index=[s for s in SAMPLE_ORDER if s in pivot_dev.index],
            columns=[c for c in CHR_ORDER if c in pivot_dev.columns],
        )
        pivot_raw = pivot_raw.reindex(
            index=pivot_dev.index, columns=pivot_dev.columns,
        )

        # Custom diverging colormap centered at 0
        vmax = max(abs(pivot_dev.min().min()), abs(pivot_dev.max().max()), 0.5)
        vmax = min(vmax, 1.5)  # cap for readability

        sns.heatmap(
            pivot_dev, ax=ax, cmap="RdBu_r", center=0,
            vmin=-vmax, vmax=vmax,
            annot=pivot_raw.round(0).astype(int) if pivot_raw.notna().all().all() else False,
            fmt="d", linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Depth deviation from sample median", "shrink": 0.8},
        )

        # Mark extreme cells
        for i_row in range(pivot_dev.shape[0]):
            for j_col in range(pivot_dev.shape[1]):
                val = pivot_dev.iloc[i_row, j_col]
                if pd.notna(val) and (val > 0.5 or val < -0.3):
                    ax.add_patch(plt.Rectangle(
                        (j_col, i_row), 1, 1,
                        fill=False, edgecolor="black", linewidth=2.5,
                    ))

        ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Sample" if ax_i == 0 else "")
        ax.tick_params(axis="x", rotation=45)

    fig.suptitle(
        "Figure 1: Per-Sample Per-Chromosome Median Read Depth\n"
        "(Color = deviation from sample genome-wide median; Black box = |dev| > 0.5 or < -0.3)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    fig.text(
        0.5, -0.04,
        CAPTION_SRC + "Annotated values = raw median read count per region",
        ha="center", fontsize=7, color="#999999",
    )
    fig.tight_layout()
    save_figure(fig, out / "cn_chr_depth_heatmap.png", dpi=180)

    return chr_med


# ══════════════════════════════════════════════════════════════════════════
# Figure 2: LOH Rate vs Chromosome Depth Deviation
# ══════════════════════════════════════════════════════════════════════════

def fig02_loh_rate_vs_depth(df: pd.DataFrame, medians: pd.DataFrame, out: Path):
    """Scatter: chromosome depth deviation vs LOH rate per sample x chr."""
    # Compute per sample x mode x chr stats
    grp = df.groupby(["sample", "mode", "Chr"]).agg(
        n_total=("loh", "count"),
        n_loh=("loh", "sum"),
        median_depth=("NumReads", "median"),
    ).reset_index()

    grp = grp.merge(medians, on=["sample", "mode"], how="left")
    grp["depth_deviation"] = (grp["median_depth"] - grp["sample_median_depth"]) / grp["sample_median_depth"]
    grp["loh_rate"] = grp["n_loh"] / grp["n_total"]

    # Filter small bins
    grp = grp[grp["n_total"] >= 20].copy()

    sample_palette = {s: sns.color_palette("tab10", n_colors=7)[i]
                      for i, s in enumerate(SAMPLE_ORDER)}

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)

    corr_results = []

    for ax_i, mode in enumerate(MODE_ORDER):
        ax = axes[ax_i]
        sub = grp[grp["mode"] == mode]

        for sample in SAMPLE_ORDER:
            ss = sub[sub["sample"] == sample]
            if len(ss) == 0:
                continue
            ax.scatter(
                ss["depth_deviation"], ss["loh_rate"],
                color=sample_palette[sample], label=sample,
                s=50, alpha=0.7, edgecolors="white", linewidth=0.5,
            )

        # Overall correlation
        if len(sub) >= 5:
            rho, p = sp_stats.spearmanr(sub["depth_deviation"], sub["loh_rate"])
            ax.text(
                0.02, 0.98,
                f"Overall Spearman rho={rho:.3f}, p={p:.2e} (n={len(sub)})",
                transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7),
            )
            corr_results.append({"mode": mode, "scope": "overall", "rho": rho, "p": p, "n": len(sub)})

            # Per-sample correlations
            for sample in SAMPLE_ORDER:
                ss = sub[sub["sample"] == sample]
                if len(ss) >= 5:
                    r, pv = sp_stats.spearmanr(ss["depth_deviation"], ss["loh_rate"])
                    corr_results.append({"mode": mode, "scope": sample, "rho": r, "p": pv, "n": len(ss)})

        ax.axhline(y=0, color="gray", linewidth=0.5, linestyle="--")
        ax.axvline(x=0, color="gray", linewidth=0.5, linestyle="--")
        ax.set_xlabel("Chromosome depth deviation\n(chr_median - sample_median) / sample_median")
        ax.set_ylabel("LOH rate (fraction of regions)" if ax_i == 0 else "")
        ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
        if ax_i == 0:
            ax.legend(fontsize=8, loc="lower right", ncol=2)

    fig.suptitle(
        "Figure 2: LOH Rate vs Chromosome Depth Deviation\n"
        "(Each point = one sample x chromosome; negative correlation = copy-loss LOH)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    fig.text(0.5, -0.04, CAPTION_SRC + "Filtered: chr bins with >= 20 regions",
             ha="center", fontsize=7, color="#999999")
    fig.tight_layout()
    save_figure(fig, out / "cn_loh_rate_vs_depth.png", dpi=180)

    return pd.DataFrame(corr_results)


# ══════════════════════════════════════════════════════════════════════════
# Figure 3: CN Category vs LOH (Stacked bar, Paired/TO x TP/FP)
# ══════════════════════════════════════════════════════════════════════════

def fig03_cn_category_loh_bar(df: pd.DataFrame, out: Path):
    """Stacked bar: LOH/non-LOH proportions within each CN category."""
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    panel_keys = [
        ("paired", "TP"), ("paired", "FP"),
        ("to", "TP"), ("to", "FP"),
    ]

    for idx, (mode, truth) in enumerate(panel_keys):
        ax = axes[idx // 2][idx % 2]
        sub = df[(df["mode"] == mode) & (df["truth_label"] == truth)]

        # Count LOH/non-LOH per CN category
        ct = sub.groupby(["cn_category", "loh"]).size().unstack(fill_value=0)
        ct = ct.reindex(index=CN_ORDER).fillna(0)

        # Normalize to proportions
        ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

        # Plot stacked bar
        if True in ct_pct.columns and False in ct_pct.columns:
            bars_noloh = ax.bar(
                range(len(CN_ORDER)), ct_pct[False],
                color="#B0BEC5", label="Non-LOH", edgecolor="white",
            )
            bars_loh = ax.bar(
                range(len(CN_ORDER)), ct_pct[True],
                bottom=ct_pct[False],
                color="#FF7043", label="LOH", edgecolor="white",
            )
        elif True in ct_pct.columns:
            bars_loh = ax.bar(
                range(len(CN_ORDER)), ct_pct[True],
                color="#FF7043", label="LOH", edgecolor="white",
            )
        else:
            bars_noloh = ax.bar(
                range(len(CN_ORDER)), ct_pct.get(False, pd.Series(0, index=CN_ORDER)),
                color="#B0BEC5", label="Non-LOH", edgecolor="white",
            )

        # Add count annotations
        for i, cat in enumerate(CN_ORDER):
            total = ct.loc[cat].sum() if cat in ct.index else 0
            loh_n = ct.loc[cat].get(True, 0) if cat in ct.index else 0
            loh_pct = (loh_n / total * 100) if total > 0 else 0
            ax.text(
                i, 102, f"n={int(total)}\n({loh_pct:.1f}% LOH)",
                ha="center", fontsize=7, va="bottom",
            )

        ax.set_xticks(range(len(CN_ORDER)))
        ax.set_xticklabels([c.replace(" ", "\n") for c in CN_ORDER], fontsize=8)
        ax.set_ylim(0, 125)
        ax.set_ylabel("Proportion (%)")
        ax.set_title(f"{mode.upper()} - {truth}", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=0))

    fig.suptitle(
        "Figure 3: LOH Proportion by Copy Number Category\n"
        "(CN classified by region depth vs sample median)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    fig.text(0.5, -0.02, CAPTION_SRC + "CN categories: <0.3x=Deep del, 0.3-0.7x=Loss, 0.7-1.3x=Neutral, 1.3-2.0x=Gain, >2.0x=High amp",
             ha="center", fontsize=7, color="#999999")
    fig.tight_layout()
    save_figure(fig, out / "cn_category_loh_bar.png", dpi=180)


# ══════════════════════════════════════════════════════════════════════════
# Figure 4: HP Balance by CN Category
# ══════════════════════════════════════════════════════════════════════════

def fig04_hp_balance_by_cn(df: pd.DataFrame, out: Path):
    """Box plot: HP balance ratio by CN category x LOH status."""
    # Compute HP balance ratio
    hp1 = df["HP1FamilyN"].fillna(0)
    hp2 = df["HP2FamilyN"].fillna(0)
    hp_max = np.maximum(hp1, hp2)
    hp_min = np.minimum(hp1, hp2)
    df = df.copy()
    df["hp_balance"] = np.where(hp_max > 0, hp_min / hp_max, np.nan)

    # Filter: need at least some HP reads
    plot_df = df[(df["hp_balance"].notna()) & (df["HP1FamilyN"] + df["HP2FamilyN"] > 2)].copy()
    plot_df["loh_label"] = plot_df["loh"].map({True: "LOH", False: "Non-LOH"})

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)

    for ax_i, mode in enumerate(MODE_ORDER):
        ax = axes[ax_i]
        sub = plot_df[plot_df["mode"] == mode]

        sns.boxplot(
            data=sub, x="cn_category", y="hp_balance",
            hue="loh_label", order=CN_ORDER,
            palette={"LOH": "#FF7043", "Non-LOH": "#78909C"},
            showfliers=False, ax=ax,
        )

        # Add sample sizes
        for i, cat in enumerate(CN_ORDER):
            for j, loh_val in enumerate([False, True]):
                ss = sub[(sub["cn_category"] == cat) & (sub["loh"] == loh_val)]
                offset = -0.2 + j * 0.4
                ax.text(
                    i + offset, -0.08,
                    f"n={len(ss)}", ha="center", fontsize=6, color="gray",
                )

        ax.set_xticklabels([c.replace(" ", "\n") for c in CN_ORDER], fontsize=8)
        ax.set_xlabel("Copy Number Category")
        ax.set_ylabel("HP Balance Ratio\nmin(HP1,HP2)/max(HP1,HP2)" if ax_i == 0 else "")
        ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
        ax.set_ylim(-0.12, 1.05)
        ax.axhline(y=0.5, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)

    fig.suptitle(
        "Figure 4: HP Balance Ratio by Copy Number Category and LOH Status\n"
        "(Lower balance = more extreme HP skew; copy-gain LOH may show extreme imbalance)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    fig.text(0.5, -0.04, CAPTION_SRC + "HP balance = min(HP1,HP2)/max(HP1,HP2); fliers hidden",
             ha="center", fontsize=7, color="#999999")
    fig.tight_layout()
    save_figure(fig, out / "cn_hp_balance_by_category.png", dpi=180)


# ══════════════════════════════════════════════════════════════════════════
# Figure 5: LOH Composition Per Sample (by CN category)
# ══════════════════════════════════════════════════════════════════════════

def fig05_loh_composition_per_sample(df: pd.DataFrame, out: Path):
    """Stacked bar: CN composition of LOH regions per sample (Paired / TO)."""
    loh_df = df[df["loh"]].copy()

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)

    for ax_i, mode in enumerate(MODE_ORDER):
        ax = axes[ax_i]
        sub = loh_df[loh_df["mode"] == mode]

        # Count per sample x CN category
        ct = sub.groupby(["sample", "cn_category"]).size().unstack(fill_value=0)
        ct = ct.reindex(columns=CN_ORDER, fill_value=0)
        ct = ct.reindex(index=[s for s in SAMPLE_ORDER if s in ct.index])

        # Normalize to proportions
        ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

        # Stacked bar
        bottom = np.zeros(len(ct_pct))
        for cat in CN_ORDER:
            if cat in ct_pct.columns:
                vals = ct_pct[cat].values
                ax.bar(
                    range(len(ct_pct)), vals, bottom=bottom,
                    color=CN_COLORS[cat], label=cat, edgecolor="white", linewidth=0.5,
                )
                # Annotate if > 5%
                for i, (v, b) in enumerate(zip(vals, bottom)):
                    if v > 5:
                        raw_n = ct[cat].iloc[i]
                        ax.text(i, b + v / 2, f"{v:.0f}%\n(n={raw_n})",
                                ha="center", va="center", fontsize=6, color="white",
                                fontweight="bold")
                bottom += vals

        ax.set_xticks(range(len(ct_pct)))
        ax.set_xticklabels(ct_pct.index, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Proportion of LOH regions (%)" if ax_i == 0 else "")
        ax.set_title(f"{mode.upper()} mode", fontsize=13, fontweight="bold")
        ax.set_ylim(0, 110)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=0))
        if ax_i == 1:
            ax.legend(fontsize=8, loc="upper right", title="CN Category")

    fig.suptitle(
        "Figure 5: Copy Number Composition of LOH Regions Per Sample\n"
        "(High Copy-neutral = likely UPD/phasing artifact; high Copy-loss = true deletion)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    fig.text(0.5, -0.04, CAPTION_SRC + "Only LOH-flagged regions included",
             ha="center", fontsize=7, color="#999999")
    fig.tight_layout()
    save_figure(fig, out / "cn_loh_composition_per_sample.png", dpi=180)


# ══════════════════════════════════════════════════════════════════════════
# Statistics table
# ══════════════════════════════════════════════════════════════════════════

def build_stats_table(
    df: pd.DataFrame,
    medians: pd.DataFrame,
    corr_df: pd.DataFrame,
    out: Path,
):
    """Produce comprehensive statistics TSV."""
    rows = []

    # Section 1: Per-sample genome-wide median depth
    rows.append({"section": "genome_median_depth", "metric": "---", "value": "---"})
    for _, r in medians.iterrows():
        rows.append({
            "section": "genome_median_depth",
            "sample": r["sample"],
            "mode": r["mode"],
            "metric": "median_NumReads",
            "value": f"{r['sample_median_depth']:.1f}",
        })

    # Section 2: LOH CN composition (overall)
    rows.append({"section": "loh_cn_composition", "metric": "---", "value": "---"})
    for mode in MODE_ORDER:
        loh_sub = df[(df["loh"]) & (df["mode"] == mode)]
        total = len(loh_sub)
        for cat in CN_ORDER:
            n = (loh_sub["cn_category"] == cat).sum()
            pct = n / total * 100 if total > 0 else 0
            rows.append({
                "section": "loh_cn_composition",
                "mode": mode,
                "metric": cat,
                "value": f"{n} ({pct:.1f}%)",
            })

    # Section 3: LOH CN composition per sample
    rows.append({"section": "loh_cn_per_sample", "metric": "---", "value": "---"})
    for mode in MODE_ORDER:
        for sample in SAMPLE_ORDER:
            loh_sub = df[(df["loh"]) & (df["mode"] == mode) & (df["sample"] == sample)]
            total = len(loh_sub)
            for cat in CN_ORDER:
                n = (loh_sub["cn_category"] == cat).sum()
                pct = n / total * 100 if total > 0 else 0
                rows.append({
                    "section": "loh_cn_per_sample",
                    "sample": sample,
                    "mode": mode,
                    "metric": cat,
                    "value": f"{n} ({pct:.1f}%)",
                })

    # Section 4: Spearman correlations
    rows.append({"section": "spearman_depth_vs_loh", "metric": "---", "value": "---"})
    for _, r in corr_df.iterrows():
        rows.append({
            "section": "spearman_depth_vs_loh",
            "mode": r["mode"],
            "sample": r.get("scope", "overall"),
            "metric": f"rho={r['rho']:.3f}",
            "value": f"p={r['p']:.2e}, n={r['n']}",
        })

    # Section 5: CN category x TP/FP x LOH cross-tab
    rows.append({"section": "cn_tp_fp_loh_crosstab", "metric": "---", "value": "---"})
    for mode in MODE_ORDER:
        for truth in ["TP", "FP"]:
            sub = df[(df["mode"] == mode) & (df["truth_label"] == truth)]
            total = len(sub)
            for cat in CN_ORDER:
                cat_sub = sub[sub["cn_category"] == cat]
                n = len(cat_sub)
                n_loh = cat_sub["loh"].sum()
                loh_rate = n_loh / n * 100 if n > 0 else 0
                rows.append({
                    "section": "cn_tp_fp_loh_crosstab",
                    "mode": mode,
                    "truth_label": truth,
                    "metric": cat,
                    "value": f"n={n}, LOH={n_loh} ({loh_rate:.1f}%)",
                })

    # Section 6: Copy-gain LOH analysis
    rows.append({"section": "copy_gain_loh_analysis", "metric": "---", "value": "---"})
    for mode in MODE_ORDER:
        loh_sub = df[(df["loh"]) & (df["mode"] == mode)]
        gain_loh = loh_sub[loh_sub["cn_category"].isin(["Copy gain", "High amplification"])]
        n_gain = len(gain_loh)
        n_total_loh = len(loh_sub)
        pct = n_gain / n_total_loh * 100 if n_total_loh > 0 else 0
        rows.append({
            "section": "copy_gain_loh_analysis",
            "mode": mode,
            "metric": "gain_or_highamp_LOH",
            "value": f"{n_gain} / {n_total_loh} ({pct:.1f}%)",
        })
        # TP vs FP breakdown in copy-gain LOH
        for truth in ["TP", "FP"]:
            gain_truth = gain_loh[gain_loh["truth_label"] == truth]
            rows.append({
                "section": "copy_gain_loh_analysis",
                "mode": mode,
                "truth_label": truth,
                "metric": "gain_LOH_count",
                "value": f"{len(gain_truth)}",
            })

    stats_df = pd.DataFrame(rows)
    stats_path = out / "cn_analysis_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)
    print(f"[CN analysis] Saved stats: {stats_path}")
    return stats_df


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    setup_plot_style()
    ensure_dir(OUT_DIR)

    print("=" * 70)
    print("LOH Copy Number Signal Detection Analysis")
    print("=" * 70)

    # Load data
    df = load_data()
    medians = compute_sample_medians(df)
    df = assign_cn_category(df, medians)

    print(f"\n[CN analysis] Sample medians:")
    print(medians.to_string(index=False))

    print(f"\n[CN analysis] CN category distribution:")
    print(df["cn_category"].value_counts().reindex(CN_ORDER))

    # Generate figures
    print("\n--- Figure 1: Chromosome depth heatmap ---")
    chr_med = fig01_chr_depth_heatmap(df, medians, OUT_DIR)

    print("\n--- Figure 2: LOH rate vs depth deviation ---")
    corr_df = fig02_loh_rate_vs_depth(df, medians, OUT_DIR)

    print("\n--- Figure 3: CN category LOH bar ---")
    fig03_cn_category_loh_bar(df, OUT_DIR)

    print("\n--- Figure 4: HP balance by CN ---")
    fig04_hp_balance_by_cn(df, OUT_DIR)

    print("\n--- Figure 5: LOH composition per sample ---")
    fig05_loh_composition_per_sample(df, OUT_DIR)

    # Statistics table
    print("\n--- Statistics table ---")
    stats = build_stats_table(df, medians, corr_df, OUT_DIR)

    # Print summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Output directory: {OUT_DIR}")
    print(f"Figures: 5")
    print(f"Stats table: cn_analysis_stats.tsv")

    # Quick summary for conclusions
    print("\n--- Quick Summary ---")
    for mode in MODE_ORDER:
        loh_sub = df[(df["loh"]) & (df["mode"] == mode)]
        cn_neutral = (loh_sub["cn_category"] == "Copy neutral").sum()
        cn_loss = (loh_sub["cn_category"].isin(["Copy loss", "Deep deletion"])).sum()
        cn_gain = (loh_sub["cn_category"].isin(["Copy gain", "High amplification"])).sum()
        total = len(loh_sub)
        print(f"\n  {mode.upper()} LOH composition (n={total}):")
        print(f"    Copy-neutral: {cn_neutral} ({cn_neutral/total*100:.1f}%)")
        print(f"    Copy-loss/deletion: {cn_loss} ({cn_loss/total*100:.1f}%)")
        print(f"    Copy-gain/amplification: {cn_gain} ({cn_gain/total*100:.1f}%)")

    return df, medians, corr_df, stats


if __name__ == "__main__":
    main()
