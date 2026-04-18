#!/usr/bin/env python3
"""LOH Region Read Depth Analysis -- Is Read Depth Really Halved? Copy Gain?

Investigates whether LOH regions show the expected read depth reduction
(~0.5x of non-LOH) or unexpectedly high depth (suggesting copy number gain),
and whether TP vs FP differ systematically in LOH read depth.

Generates:
  - Figure 1: LOH vs non-LOH read depth distribution (violin + box)
  - Figure 2: Per-sample LOH read depth ratio (bar chart)
  - Figure 3: HP1 vs HP2 balance scatter (LOH vs non-LOH)
  - Figure 4: Per-chromosome read depth heatmap by LOH status
  - Statistics table: loh_read_depth_stats.tsv

Usage:
    python3 build_loh_read_depth_analysis.py
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
    add_caption,
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
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sp_stats
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/"
    "2026/04/20260401_LOH_weekly_report_draft/images"
)

CAPTION_SRC = (
    "Source: all_region_rows.tsv.gz "
    "(20260327_loh_round1_cross_sample_audit) | "
)

# Chromosome display order
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ── Data Loading ─────────────────────────────────────────────────────────────
def load_data() -> pd.DataFrame:
    """Load master dataset with only required columns."""
    usecols = [
        "NumReads", "HP1FamilyN", "HP2FamilyN", "NHP0", "NHP3",
        "HP_Ratio", "Potential_LOH", "core_loh_like",
        "effective_hp_reads", "truth_label", "mode", "sample", "Chr",
    ]
    df = load_master_dataset(usecols=usecols)

    # Convert boolean-like columns
    for col in ["Potential_LOH", "core_loh_like"]:
        if df[col].dtype == object:
            df[col] = df[col].map({"True": True, "False": False})
        df[col] = df[col].astype(bool)

    # Derived columns
    df["loh_status"] = df["core_loh_like"].map({True: "LOH", False: "non-LOH"})
    hp_sum = df["HP1FamilyN"] + df["HP2FamilyN"]
    df["hp_imbalance"] = np.where(
        hp_sum > 0,
        np.abs(df["HP1FamilyN"] - df["HP2FamilyN"]) / hp_sum,
        np.nan,
    )
    df["minor_hp_frac"] = np.where(
        hp_sum > 0,
        np.minimum(df["HP1FamilyN"], df["HP2FamilyN"]) / hp_sum,
        np.nan,
    )
    print(f"[load_data] LOH rows: {df['core_loh_like'].sum():,}, "
          f"non-LOH: {(~df['core_loh_like']).sum():,}")
    return df


# ── Figure 1: LOH vs non-LOH Read Depth Distribution ────────────────────────
def figure1_depth_distribution(df: pd.DataFrame) -> Path:
    """Violin + strip plot of NumReads by LOH status, faceted by mode."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax_idx, mode_val in enumerate(MODE_ORDER):
        ax = axes[ax_idx]
        sub = df[df["mode"] == mode_val].copy()

        # Clip for visualization (keep originals for stats annotation)
        clip_upper = sub["NumReads"].quantile(0.99)

        # Violin plot (LOH status)
        sns.violinplot(
            data=sub, x="loh_status", y="NumReads",
            hue="truth_label", hue_order=["TP", "FP"],
            palette=TRUTH_PALETTE, split=True,
            inner="quartile", linewidth=0.8, cut=0,
            ax=ax, order=["non-LOH", "LOH"],
        )
        ax.set_ylim(0, clip_upper * 1.05)
        ax.set_title(f"{mode_val.upper()}", fontsize=13, fontweight="bold")
        ax.set_xlabel("LOH Status (core_loh_like)")
        if ax_idx == 0:
            ax.set_ylabel("NumReads")
        else:
            ax.set_ylabel("")

        # Annotate medians
        for loh_val_idx, loh_val in enumerate(["non-LOH", "LOH"]):
            for tl in ["TP", "FP"]:
                vals = sub[(sub["loh_status"] == loh_val) & (sub["truth_label"] == tl)]["NumReads"]
                if len(vals) > 0:
                    med = vals.median()
                    color = COLOR_TP if tl == "TP" else COLOR_FP
                    offset = -0.12 if tl == "TP" else 0.12
                    ax.annotate(
                        f"med={med:.0f}",
                        xy=(loh_val_idx + offset, min(med, clip_upper * 0.95)),
                        fontsize=7, color=color, fontweight="bold",
                        ha="center", va="bottom",
                    )

        # Legend only on first panel
        if ax_idx > 0:
            ax.get_legend().remove()
        else:
            ax.legend(title="Truth", loc="upper right", fontsize=8)

    fig.suptitle(
        "Figure 1: Read Depth (NumReads) in LOH vs non-LOH Regions",
        fontsize=14, fontweight="bold", y=1.02,
    )
    add_caption(fig, CAPTION_SRC + "core_loh_like definition: HP_Ratio<0.1 or >0.9, effective_hp_reads>=30")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "loh_read_depth_distribution.png")


# ── Figure 2: Per-Sample LOH Depth Ratio ─────────────────────────────────────
def figure2_depth_ratio(df: pd.DataFrame) -> tuple[Path, pd.DataFrame]:
    """Per-sample median(NumReads|LOH) / median(NumReads|non-LOH) bar chart."""
    setup_plot_style()

    # Compute ratios per sample x mode x truth_label
    records = []
    for sample in SAMPLE_ORDER:
        for mode_val in MODE_ORDER:
            for tl in ["TP", "FP"]:
                sub = df[(df["sample"] == sample) & (df["mode"] == mode_val) & (df["truth_label"] == tl)]
                loh_reads = sub[sub["core_loh_like"]]["NumReads"]
                nonloh_reads = sub[~sub["core_loh_like"]]["NumReads"]
                if len(loh_reads) >= 5 and len(nonloh_reads) >= 5:
                    ratio = loh_reads.median() / nonloh_reads.median()
                    # Mann-Whitney U test
                    u_stat, u_p = sp_stats.mannwhitneyu(
                        loh_reads, nonloh_reads, alternative="two-sided"
                    )
                else:
                    ratio = np.nan
                    u_p = np.nan
                records.append({
                    "sample": sample, "mode": mode_val, "truth_label": tl,
                    "ratio": ratio,
                    "n_loh": len(loh_reads), "n_nonloh": len(nonloh_reads),
                    "med_loh": loh_reads.median() if len(loh_reads) > 0 else np.nan,
                    "med_nonloh": nonloh_reads.median() if len(nonloh_reads) > 0 else np.nan,
                    "u_pvalue": u_p,
                })
    ratio_df = pd.DataFrame(records)

    # Plot: 2x2 grid (mode x truth_label)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharey=True)
    combos = [("paired", "TP"), ("paired", "FP"), ("to", "TP"), ("to", "FP")]

    for ax_idx, (mode_val, tl) in enumerate(combos):
        ax = axes[ax_idx // 2, ax_idx % 2]
        sub = ratio_df[(ratio_df["mode"] == mode_val) & (ratio_df["truth_label"] == tl)].copy()
        sub = sub.dropna(subset=["ratio"])

        if len(sub) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"{mode_val.upper()} - {tl}")
            continue

        # Sort by ratio
        sub = sub.sort_values("ratio", ascending=True)
        colors = [COLOR_TP if tl == "TP" else COLOR_FP] * len(sub)

        bars = ax.barh(sub["sample"], sub["ratio"], color=colors, alpha=0.8, edgecolor="white")

        # Reference lines
        ax.axvline(0.5, color="red", linestyle="--", linewidth=1.2, alpha=0.7, label="0.5x (expected LOH)")
        ax.axvline(1.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.5, label="1.0x (no change)")

        # Annotate values
        for bar, ratio_val in zip(bars, sub["ratio"]):
            ax.text(
                bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2,
                f"{ratio_val:.2f}", va="center", fontsize=8, fontweight="bold",
            )

        ax.set_title(f"{mode_val.upper()} - {tl}", fontsize=12, fontweight="bold")
        ax.set_xlabel("Median Read Depth Ratio (LOH / non-LOH)")
        ax.legend(loc="lower right", fontsize=7)
        ax.set_xlim(0, max(sub["ratio"].max() * 1.2, 1.5))

    fig.suptitle(
        "Figure 2: Per-Sample LOH Read Depth Ratio\n"
        "median(NumReads|LOH) / median(NumReads|non-LOH)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    add_caption(fig, CAPTION_SRC + "Red dashed=0.5x (expected if LOH loses half reads), gray dashed=1.0x")
    fig.tight_layout()
    path = save_figure(fig, OUT_DIR / "loh_read_depth_ratio_per_sample.png")
    return path, ratio_df


# ── Figure 3: HP1 vs HP2 Balance Scatter ─────────────────────────────────────
def figure3_hp_balance(df: pd.DataFrame) -> Path:
    """Scatter plot of HP1FamilyN vs HP2FamilyN, colored by LOH status."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    loh_palette = {"LOH": "#E53935", "non-LOH": "#1E88E5"}

    for ax_idx, mode_val in enumerate(MODE_ORDER):
        ax = axes[ax_idx]
        sub = df[df["mode"] == mode_val].copy()

        # Subsample for plotting performance (max 20K per group)
        rng = np.random.RandomState(42)
        groups = []
        for loh_val in ["non-LOH", "LOH"]:
            g = sub[sub["loh_status"] == loh_val]
            if len(g) > 20000:
                g = g.sample(20000, random_state=rng)
            groups.append(g)
        plot_df = pd.concat(groups)

        # Plot non-LOH first (background), then LOH on top
        for loh_val, alpha_val, marker_size in [("non-LOH", 0.15, 3), ("LOH", 0.35, 5)]:
            mask = plot_df["loh_status"] == loh_val
            ax.scatter(
                plot_df.loc[mask, "HP1FamilyN"],
                plot_df.loc[mask, "HP2FamilyN"],
                c=loh_palette[loh_val], s=marker_size, alpha=alpha_val,
                label=loh_val, rasterized=True,
            )

        # Diagonal line (balanced HP)
        max_val = max(sub["HP1FamilyN"].quantile(0.99), sub["HP2FamilyN"].quantile(0.99))
        ax.plot([0, max_val], [0, max_val], "k--", linewidth=0.8, alpha=0.4, label="HP1=HP2")

        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        ax.set_xlabel("HP1FamilyN")
        ax.set_ylabel("HP2FamilyN")
        ax.set_title(f"{mode_val.upper()}", fontsize=13, fontweight="bold")
        ax.legend(loc="upper right", fontsize=8, markerscale=3)
        ax.set_aspect("equal")

    fig.suptitle(
        "Figure 3: HP1 vs HP2 Read Count Balance\n"
        "LOH regions should cluster near axes (one HP ~0)",
        fontsize=14, fontweight="bold", y=1.04,
    )
    add_caption(fig, CAPTION_SRC + "LOH=red (core_loh_like), non-LOH=blue; diagonal=balanced HP")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "loh_hp_balance_scatter.png")


# ── Figure 4: Per-Chromosome Depth Heatmap ───────────────────────────────────
def figure4_chr_heatmap(df: pd.DataFrame) -> Path:
    """Heatmap: per-sample x per-chromosome median read depth, LOH vs non-LOH."""
    setup_plot_style()

    # Filter to standard chromosomes
    df_chr = df[df["Chr"].isin(CHR_ORDER)].copy()

    # We'll make a heatmap of ratio: median(LOH) / median(non-LOH) per sample x chr
    # Plus separate panels showing absolute values

    fig, axes = plt.subplots(2, 1, figsize=(18, 12), gridspec_kw={"height_ratios": [1, 1]})

    for ax_idx, mode_val in enumerate(MODE_ORDER):
        ax = axes[ax_idx]
        sub = df_chr[df_chr["mode"] == mode_val]

        # Compute median NumReads per sample x chr x loh_status
        pivot_data = {}
        for sample in SAMPLE_ORDER:
            row = {}
            for chrom in CHR_ORDER:
                cell = sub[(sub["sample"] == sample) & (sub["Chr"] == chrom)]
                loh_med = cell[cell["core_loh_like"]]["NumReads"].median()
                nonloh_med = cell[~cell["core_loh_like"]]["NumReads"].median()
                if pd.notna(loh_med) and pd.notna(nonloh_med) and nonloh_med > 0:
                    row[chrom] = loh_med / nonloh_med
                else:
                    row[chrom] = np.nan
            pivot_data[sample] = row

        heatmap_df = pd.DataFrame(pivot_data).T
        # Reorder columns
        heatmap_df = heatmap_df[[c for c in CHR_ORDER if c in heatmap_df.columns]]

        sns.heatmap(
            heatmap_df, ax=ax, cmap="RdYlBu_r", center=1.0,
            vmin=0.3, vmax=1.7,
            annot=True, fmt=".2f", annot_kws={"fontsize": 6},
            linewidths=0.3, linecolor="white",
            cbar_kws={"label": "Depth Ratio (LOH/non-LOH)", "shrink": 0.8},
        )
        ax.set_title(f"{mode_val.upper()} -- Median Read Depth Ratio (LOH / non-LOH) per Chromosome",
                      fontsize=12, fontweight="bold")
        ax.set_ylabel("Sample")
        ax.set_xlabel("Chromosome")

    fig.suptitle(
        "Figure 4: Per-Chromosome LOH Read Depth Ratio\n"
        "Values <0.5 suggest pure LOH loss; >1.0 may indicate copy number gain",
        fontsize=14, fontweight="bold", y=1.02,
    )
    add_caption(fig, CAPTION_SRC + "Ratio = median(NumReads|LOH) / median(NumReads|non-LOH); gray=insufficient data")
    fig.tight_layout()
    return save_figure(fig, OUT_DIR / "loh_chr_depth_heatmap.png")


# ── Statistics Table ─────────────────────────────────────────────────────────
def compute_stats_table(df: pd.DataFrame) -> pd.DataFrame:
    """Compute summary statistics per sample x mode x truth_label x LOH_status."""
    groups = df.groupby(["sample", "mode", "truth_label", "loh_status"])

    records = []
    for (sample, mode_val, tl, loh_val), g in groups:
        reads = g["NumReads"]
        imbal = g["hp_imbalance"]
        records.append({
            "sample": sample,
            "mode": mode_val,
            "truth_label": tl,
            "loh_status": loh_val,
            "count": len(g),
            "median_NumReads": reads.median(),
            "mean_NumReads": reads.mean(),
            "std_NumReads": reads.std(),
            "q25_NumReads": reads.quantile(0.25),
            "q75_NumReads": reads.quantile(0.75),
            "median_hp_imbalance": imbal.median(),
            "mean_hp_imbalance": imbal.mean(),
            "median_HP1FamilyN": g["HP1FamilyN"].median(),
            "median_HP2FamilyN": g["HP2FamilyN"].median(),
            "median_NHP0": g["NHP0"].median(),
        })

    stats_df = pd.DataFrame(records)
    # Sort
    stats_df["sample"] = pd.Categorical(stats_df["sample"], categories=SAMPLE_ORDER, ordered=True)
    stats_df["mode"] = pd.Categorical(stats_df["mode"], categories=MODE_ORDER, ordered=True)
    stats_df = stats_df.sort_values(["sample", "mode", "truth_label", "loh_status"]).reset_index(drop=True)
    return stats_df


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    ensure_dir(OUT_DIR)
    df = load_data()

    print("\n=== Figure 1: LOH vs non-LOH Read Depth Distribution ===")
    p1 = figure1_depth_distribution(df)
    print(f"  -> {p1}")

    print("\n=== Figure 2: Per-Sample LOH Read Depth Ratio ===")
    p2, ratio_df = figure2_depth_ratio(df)
    print(f"  -> {p2}")
    print("\n  Ratio summary:")
    print(ratio_df.to_string(index=False))

    print("\n=== Figure 3: HP1 vs HP2 Balance Scatter ===")
    p3 = figure3_hp_balance(df)
    print(f"  -> {p3}")

    print("\n=== Figure 4: Per-Chromosome Depth Heatmap ===")
    p4 = figure4_chr_heatmap(df)
    print(f"  -> {p4}")

    print("\n=== Statistics Table ===")
    stats_df = compute_stats_table(df)
    stats_path = OUT_DIR / "loh_read_depth_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False, float_format="%.2f")
    print(f"  -> {stats_path}")

    # Print key summary findings
    print("\n" + "=" * 70)
    print("KEY FINDINGS SUMMARY")
    print("=" * 70)

    # Overall LOH vs non-LOH depth
    for mode_val in MODE_ORDER:
        for tl in ["TP", "FP"]:
            sub = df[(df["mode"] == mode_val) & (df["truth_label"] == tl)]
            loh = sub[sub["core_loh_like"]]["NumReads"]
            nonloh = sub[~sub["core_loh_like"]]["NumReads"]
            if len(loh) > 0 and len(nonloh) > 0:
                ratio = loh.median() / nonloh.median()
                u_stat, u_p = sp_stats.mannwhitneyu(loh, nonloh, alternative="two-sided")
                print(f"  {mode_val.upper()}/{tl}: "
                      f"LOH median={loh.median():.1f}, non-LOH median={nonloh.median():.1f}, "
                      f"ratio={ratio:.3f}, Mann-Whitney p={u_p:.2e}")

    # HP imbalance comparison
    print("\n  HP Imbalance (|HP1-HP2|/(HP1+HP2)):")
    for mode_val in MODE_ORDER:
        sub = df[df["mode"] == mode_val]
        loh_imbal = sub[sub["core_loh_like"]]["hp_imbalance"].dropna()
        nonloh_imbal = sub[~sub["core_loh_like"]]["hp_imbalance"].dropna()
        if len(loh_imbal) > 0 and len(nonloh_imbal) > 0:
            print(f"    {mode_val.upper()}: LOH median={loh_imbal.median():.3f}, "
                  f"non-LOH median={nonloh_imbal.median():.3f}")

    # Per-sample ratio summary
    print("\n  Per-sample ratio range:")
    for mode_val in MODE_ORDER:
        sub = ratio_df[ratio_df["mode"] == mode_val].dropna(subset=["ratio"])
        if len(sub) > 0:
            print(f"    {mode_val.upper()}: min={sub['ratio'].min():.3f}, "
                  f"max={sub['ratio'].max():.3f}, median={sub['ratio'].median():.3f}")

    print("\n  Done. All outputs in:", OUT_DIR)


if __name__ == "__main__":
    main()
