#!/usr/bin/env python3
"""
Concordance Quadrant x Read Depth Analysis
-------------------------------------------
Investigate whether TO-only LOH is driven by read depth artefacts
(low effective_hp_reads causing random HP ratio extremes) or by
self-phasing circular dependency in LongPhase-TO.

Input : same_locus_compare.tsv (288,609 matched loci, 'both' overlap class)
Output: 5 figures + 1 stats TSV
"""

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

# ── paths ──────────────────────────────────────────────────────────────
INPUT_TSV = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260327_loh_round1_cross_sample_audit/same_locus_compare.tsv"
)
OUT_DIR = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/"
    "2026/04/20260401_LOH_weekly_report_draft/images"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── palette & order ───────────────────────────────────────────────────
QUAD_ORDER = ["both_LOH", "TO_only_LOH", "Paired_only_LOH", "neither_LOH"]
QUAD_PALETTE = {
    "both_LOH":        "#4C72B0",
    "TO_only_LOH":     "#DD8452",
    "Paired_only_LOH": "#55A868",
    "neither_LOH":     "#C44E52",
}

SAMPLE_ORDER = [
    "COLO829", "HCC1395", "HCC1395_DORADO", "HCC1937",
    "HCC1954", "H1437", "H2009",
]


def load_data() -> pd.DataFrame:
    """Load same_locus_compare.tsv, keep 'both' overlap class, derive quadrant."""
    print("[INFO] Loading data ...")
    cols_needed = [
        "sample", "locus_overlap_class",
        "paired_to_concordance",
        "core_loh_like_paired", "core_loh_like_to",
        "hp_ratio_core_paired", "hp_ratio_core_to",
        "effective_hp_reads_paired", "effective_hp_reads_to",
        "truth_label_paired", "truth_label_to",
        "caller_dp_paired", "caller_dp_to",
    ]
    df = pd.read_csv(INPUT_TSV, sep="\t", usecols=cols_needed, low_memory=False)
    df = df[df["locus_overlap_class"] == "both"].copy()
    print(f"[INFO] Loaded {len(df):,} matched loci (overlap=both)")

    # Derive quadrant from core_loh_like columns (independently computed)
    def _quadrant(row):
        p = str(row["core_loh_like_paired"]).strip().lower() == "true"
        t = str(row["core_loh_like_to"]).strip().lower() == "true"
        if p and t:
            return "both_LOH"
        elif (not p) and t:
            return "TO_only_LOH"
        elif p and (not t):
            return "Paired_only_LOH"
        else:
            return "neither_LOH"

    df["quadrant"] = df.apply(_quadrant, axis=1)
    print("[INFO] Quadrant counts:")
    for q in QUAD_ORDER:
        n = (df["quadrant"] == q).sum()
        print(f"       {q:20s} = {n:>8,} ({100*n/len(df):.1f}%)")

    # Convert numeric columns
    for c in ["hp_ratio_core_paired", "hp_ratio_core_to",
              "effective_hp_reads_paired", "effective_hp_reads_to",
              "caller_dp_paired", "caller_dp_to"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Verify truth labels
    tp_p = (df["truth_label_paired"] == "TP").sum()
    tp_t = (df["truth_label_to"] == "TP").sum()
    fp_p = (df["truth_label_paired"] == "FP").sum()
    fp_t = (df["truth_label_to"] == "FP").sum()
    print(f"[INFO] truth_label_paired: TP={tp_p:,}, FP={fp_p:,}, other={len(df)-tp_p-fp_p:,}")
    print(f"[INFO] truth_label_to:     TP={tp_t:,}, FP={fp_t:,}, other={len(df)-tp_t-fp_t:,}")

    return df


# ── Figure 1: Quadrant Read Depth Distribution ───────────────────────
def fig1_quadrant_depth(df: pd.DataFrame):
    """Box + strip plot of effective_hp_reads by quadrant, split by Paired/TO."""
    print("[FIG1] Quadrant Read Depth Distribution ...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)

    for ax, mode in zip(axes, ["paired", "to"]):
        col = f"effective_hp_reads_{mode}"
        data = df[["quadrant", col]].dropna()

        sns.boxplot(
            data=data, x="quadrant", y=col, order=QUAD_ORDER,
            palette=QUAD_PALETTE, showfliers=False, width=0.6, ax=ax,
        )
        # add swarm-like jitter for small groups, or just medians annotated
        medians = data.groupby("quadrant")[col].median()
        for i, q in enumerate(QUAD_ORDER):
            med = medians.get(q, np.nan)
            n = (data["quadrant"] == q).sum()
            ax.text(
                i, med, f"  med={med:.0f}\n  n={n:,}",
                fontsize=8, va="bottom", ha="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8),
            )

        ax.set_title(f"effective_hp_reads ({mode.upper()})", fontsize=13, fontweight="bold")
        ax.set_xlabel("Concordance Quadrant", fontsize=11)
        ax.set_ylabel("effective_hp_reads" if mode == "paired" else "", fontsize=11)
        ax.set_ylim(0, min(300, data[col].quantile(0.99) * 1.2))
        ax.tick_params(axis="x", rotation=25)

    fig.suptitle(
        "Figure 1: Effective HP Reads by Concordance Quadrant",
        fontsize=14, fontweight="bold", y=1.01,
    )
    fig.tight_layout()
    out = OUT_DIR / "concordance_quadrant_depth.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"       Saved -> {out}")


# ── Figure 2: HP Ratio Scatter (Paired vs TO) ────────────────────────
def fig2_hp_ratio_scatter(df: pd.DataFrame):
    """Scatter of hp_ratio_core_paired vs hp_ratio_core_to, coloured by quadrant."""
    print("[FIG2] HP Ratio Scatter ...")
    data = df[["quadrant", "hp_ratio_core_paired", "hp_ratio_core_to"]].dropna()

    fig, ax = plt.subplots(figsize=(10, 9))

    # Draw each quadrant separately (back-to-front by size for visibility)
    for q in ["neither_LOH", "both_LOH", "Paired_only_LOH", "TO_only_LOH"]:
        sub = data[data["quadrant"] == q]
        ax.scatter(
            sub["hp_ratio_core_paired"], sub["hp_ratio_core_to"],
            c=QUAD_PALETTE[q], label=f"{q} (n={len(sub):,})",
            s=3, alpha=0.15, edgecolors="none", rasterized=True,
        )

    # density contours for TO_only_LOH
    to_only = data[data["quadrant"] == "TO_only_LOH"]
    if len(to_only) > 100:
        try:
            sns.kdeplot(
                x=to_only["hp_ratio_core_paired"],
                y=to_only["hp_ratio_core_to"],
                levels=5, color=QUAD_PALETTE["TO_only_LOH"],
                linewidths=1.2, ax=ax,
            )
        except Exception:
            pass

    ax.set_xlabel("HP Ratio (Paired)", fontsize=12)
    ax.set_ylabel("HP Ratio (TO)", fontsize=12)
    ax.set_title(
        "Figure 2: HP Ratio Paired vs TO, Coloured by Concordance Quadrant",
        fontsize=13, fontweight="bold",
    )
    ax.legend(markerscale=4, fontsize=10, loc="upper left")

    # reference lines for LOH threshold (0.1 / 0.9)
    for v in [0.1, 0.9]:
        ax.axhline(v, ls="--", lw=0.7, color="gray", alpha=0.5)
        ax.axvline(v, ls="--", lw=0.7, color="gray", alpha=0.5)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()
    out = OUT_DIR / "concordance_hp_ratio_scatter.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"       Saved -> {out}")


# ── Figure 3: Depth vs HP Ratio 2D Hexbin ────────────────────────────
def fig3_depth_vs_hpratio(df: pd.DataFrame):
    """Hexbin of effective_hp_reads vs hp_ratio for TO and Paired, by quadrant group."""
    print("[FIG3] Depth vs HP Ratio Hexbin ...")
    groups = ["TO_only_LOH", "both_LOH", "neither_LOH"]

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    for row_idx, mode in enumerate(["to", "paired"]):
        col_depth = f"effective_hp_reads_{mode}"
        col_ratio = f"hp_ratio_core_{mode}"

        for col_idx, grp in enumerate(groups):
            ax = axes[row_idx, col_idx]
            sub = df[df["quadrant"] == grp][[col_depth, col_ratio]].dropna()
            sub = sub[sub[col_depth] <= 200]  # clip for visibility

            if len(sub) < 10:
                ax.text(0.5, 0.5, f"n={len(sub)}", transform=ax.transAxes, ha="center")
                continue

            hb = ax.hexbin(
                sub[col_depth], sub[col_ratio],
                gridsize=50, cmap="YlOrRd", mincnt=1,
                extent=[0, 200, 0, 1],
            )
            fig.colorbar(hb, ax=ax, shrink=0.7, label="count")

            ax.set_title(
                f"{mode.upper()} | {grp}\n(n={len(sub):,})",
                fontsize=11, fontweight="bold",
            )
            ax.set_xlabel(f"effective_hp_reads_{mode}")
            ax.set_ylabel(f"hp_ratio_core_{mode}" if col_idx == 0 else "")

    fig.suptitle(
        "Figure 3: Read Depth vs HP Ratio (Top=TO, Bottom=Paired)",
        fontsize=14, fontweight="bold", y=1.01,
    )
    fig.tight_layout()
    out = OUT_DIR / "concordance_depth_vs_hpratio.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"       Saved -> {out}")


# ── Figure 4: Per-Sample Quadrant Depth ──────────────────────────────
def fig4_per_sample_depth(df: pd.DataFrame):
    """Grouped bar chart: median effective_hp_reads per sample per quadrant."""
    print("[FIG4] Per-Sample Quadrant Depth ...")

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)

    for ax, mode in zip(axes, ["paired", "to"]):
        col = f"effective_hp_reads_{mode}"
        tbl = (
            df.groupby(["sample", "quadrant"])[col]
            .median()
            .unstack("quadrant")
            .reindex(index=SAMPLE_ORDER, columns=QUAD_ORDER)
        )

        tbl.plot.bar(
            ax=ax, color=[QUAD_PALETTE[q] for q in QUAD_ORDER],
            width=0.78, edgecolor="white", linewidth=0.5,
        )
        ax.set_title(f"Median effective_hp_reads ({mode.upper()})", fontsize=13, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Median effective_hp_reads" if mode == "paired" else "")
        ax.tick_params(axis="x", rotation=30)
        ax.legend(fontsize=8, ncol=2, loc="upper right")

        # annotate bar values
        for container in ax.containers:
            ax.bar_label(container, fmt="%.0f", fontsize=6, padding=2, rotation=90)

    fig.suptitle(
        "Figure 4: Per-Sample Median Effective HP Reads by Quadrant",
        fontsize=14, fontweight="bold", y=1.01,
    )
    fig.tight_layout()
    out = OUT_DIR / "concordance_per_sample_depth.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"       Saved -> {out}")


# ── Figure 5: Depth Difference (TO - Paired) by Quadrant ─────────────
def fig5_depth_diff(df: pd.DataFrame):
    """Violin plot of (effective_hp_reads_to - effective_hp_reads_paired) per quadrant."""
    print("[FIG5] Depth Difference (TO - Paired) ...")
    data = df[["quadrant", "effective_hp_reads_to", "effective_hp_reads_paired"]].dropna()
    data["depth_diff"] = data["effective_hp_reads_to"] - data["effective_hp_reads_paired"]

    # Pre-compute ylim for proper annotation
    ymin = data["depth_diff"].quantile(0.01) * 1.2
    ymax = data["depth_diff"].quantile(0.99) * 1.2

    fig, ax = plt.subplots(figsize=(12, 7))

    sns.violinplot(
        data=data, x="quadrant", y="depth_diff",
        order=QUAD_ORDER, palette=QUAD_PALETTE,
        cut=0, inner="box", scale="width", ax=ax,
    )
    ax.axhline(0, ls="--", lw=1, color="black", alpha=0.6)
    ax.set_ylim(ymin, ymax)

    # Annotate medians and means within the visible ylim
    for i, q in enumerate(QUAD_ORDER):
        sub = data[data["quadrant"] == q]["depth_diff"]
        med = sub.median()
        mean = sub.mean()
        n = len(sub)
        ax.text(
            i, ymax * 0.82,
            f"med={med:.1f}\nmean={mean:.1f}\nn={n:,}",
            ha="center", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8),
        )

    ax.set_title(
        "Figure 5: Read Depth Difference (TO - Paired) by Quadrant",
        fontsize=13, fontweight="bold",
    )
    ax.set_xlabel("Concordance Quadrant", fontsize=11)
    ax.set_ylabel("effective_hp_reads_to - effective_hp_reads_paired", fontsize=11)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    out = OUT_DIR / "concordance_depth_diff.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"       Saved -> {out}")


# ── Statistical Tests ─────────────────────────────────────────────────
def run_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Wilcoxon rank-sum tests, Cohen's d, low-depth proportions."""
    print("[STATS] Running statistical tests ...")
    results = []

    # Helper: Cohen's d
    def cohens_d(a, b):
        na, nb = len(a), len(b)
        if na < 2 or nb < 2:
            return np.nan
        pooled_std = np.sqrt(((na - 1) * a.std()**2 + (nb - 1) * b.std()**2) / (na + nb - 2))
        if pooled_std == 0:
            return 0.0
        return (a.mean() - b.mean()) / pooled_std

    # Extract groups
    grps = {}
    for q in QUAD_ORDER:
        sub = df[df["quadrant"] == q]
        grps[q] = {
            "eff_to": sub["effective_hp_reads_to"].dropna().values,
            "eff_paired": sub["effective_hp_reads_paired"].dropna().values,
            "hp_ratio_to": sub["hp_ratio_core_to"].dropna().values,
            "hp_ratio_paired": sub["hp_ratio_core_paired"].dropna().values,
        }

    # ── Test 1: TO_only_LOH vs both_LOH on effective_hp_reads_to
    a = grps["TO_only_LOH"]["eff_to"]
    b = grps["both_LOH"]["eff_to"]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    results.append({
        "test": "TO_only vs both_LOH | eff_hp_reads_TO",
        "group_A": "TO_only_LOH",
        "group_B": "both_LOH",
        "metric": "effective_hp_reads_to",
        "n_A": len(a), "n_B": len(b),
        "median_A": np.median(a), "median_B": np.median(b),
        "mean_A": np.mean(a), "mean_B": np.mean(b),
        "U_statistic": U, "p_value": p, "cohens_d": d,
    })

    # ── Test 2: TO_only_LOH vs neither on effective_hp_reads_to
    a = grps["TO_only_LOH"]["eff_to"]
    b = grps["neither_LOH"]["eff_to"]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    results.append({
        "test": "TO_only vs neither | eff_hp_reads_TO",
        "group_A": "TO_only_LOH",
        "group_B": "neither_LOH",
        "metric": "effective_hp_reads_to",
        "n_A": len(a), "n_B": len(b),
        "median_A": np.median(a), "median_B": np.median(b),
        "mean_A": np.mean(a), "mean_B": np.mean(b),
        "U_statistic": U, "p_value": p, "cohens_d": d,
    })

    # ── Test 3: TO_only_LOH vs both_LOH on effective_hp_reads_paired
    a = grps["TO_only_LOH"]["eff_paired"]
    b = grps["both_LOH"]["eff_paired"]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    results.append({
        "test": "TO_only vs both_LOH | eff_hp_reads_Paired",
        "group_A": "TO_only_LOH",
        "group_B": "both_LOH",
        "metric": "effective_hp_reads_paired",
        "n_A": len(a), "n_B": len(b),
        "median_A": np.median(a), "median_B": np.median(b),
        "mean_A": np.mean(a), "mean_B": np.mean(b),
        "U_statistic": U, "p_value": p, "cohens_d": d,
    })

    # ── Test 4: TO_only_LOH vs neither on effective_hp_reads_paired
    a = grps["TO_only_LOH"]["eff_paired"]
    b = grps["neither_LOH"]["eff_paired"]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    results.append({
        "test": "TO_only vs neither | eff_hp_reads_Paired",
        "group_A": "TO_only_LOH",
        "group_B": "neither_LOH",
        "metric": "effective_hp_reads_paired",
        "n_A": len(a), "n_B": len(b),
        "median_A": np.median(a), "median_B": np.median(b),
        "mean_A": np.mean(a), "mean_B": np.mean(b),
        "U_statistic": U, "p_value": p, "cohens_d": d,
    })

    # ── Test 5: HP ratio TO — TO_only vs both_LOH
    a = grps["TO_only_LOH"]["hp_ratio_to"]
    b = grps["both_LOH"]["hp_ratio_to"]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    results.append({
        "test": "TO_only vs both_LOH | hp_ratio_TO",
        "group_A": "TO_only_LOH",
        "group_B": "both_LOH",
        "metric": "hp_ratio_core_to",
        "n_A": len(a), "n_B": len(b),
        "median_A": np.median(a), "median_B": np.median(b),
        "mean_A": np.mean(a), "mean_B": np.mean(b),
        "U_statistic": U, "p_value": p, "cohens_d": d,
    })

    # ── Test 6: HP ratio Paired — TO_only vs both_LOH
    a = grps["TO_only_LOH"]["hp_ratio_paired"]
    b = grps["both_LOH"]["hp_ratio_paired"]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    results.append({
        "test": "TO_only vs both_LOH | hp_ratio_Paired",
        "group_A": "TO_only_LOH",
        "group_B": "both_LOH",
        "metric": "hp_ratio_core_paired",
        "n_A": len(a), "n_B": len(b),
        "median_A": np.median(a), "median_B": np.median(b),
        "mean_A": np.mean(a), "mean_B": np.mean(b),
        "U_statistic": U, "p_value": p, "cohens_d": d,
    })

    stats_df = pd.DataFrame(results)

    # ── Low-depth proportion by quadrant ──────────────────────────────
    print("\n[STATS] Low-depth proportion (effective_hp_reads <= 30):")
    lowdepth_rows = []
    for q in QUAD_ORDER:
        sub = df[df["quadrant"] == q]
        for mode in ["paired", "to"]:
            col = f"effective_hp_reads_{mode}"
            valid = sub[col].dropna()
            low = (valid <= 30).sum()
            pct = 100 * low / len(valid) if len(valid) > 0 else 0
            print(f"       {q:20s} | {mode:6s} | low<=30: {low:>7,}/{len(valid):>7,} ({pct:.1f}%)")
            lowdepth_rows.append({
                "test": f"low_depth_pct | {q} | {mode}",
                "group_A": q, "group_B": "-",
                "metric": f"effective_hp_reads_{mode} <= 30",
                "n_A": int(low), "n_B": int(len(valid)),
                "median_A": pct, "median_B": np.nan,
                "mean_A": np.nan, "mean_B": np.nan,
                "U_statistic": np.nan, "p_value": np.nan, "cohens_d": np.nan,
            })
    stats_df = pd.concat([stats_df, pd.DataFrame(lowdepth_rows)], ignore_index=True)

    # ── Additional: depth quantile summary per quadrant ───────────────
    print("\n[STATS] Depth quantiles per quadrant:")
    for q in QUAD_ORDER:
        sub = df[df["quadrant"] == q]
        for mode in ["paired", "to"]:
            col = f"effective_hp_reads_{mode}"
            valid = sub[col].dropna()
            if len(valid) == 0:
                continue
            qs = valid.quantile([0.05, 0.25, 0.50, 0.75, 0.95]).values
            print(f"       {q:20s} | {mode:6s} | p5={qs[0]:.0f}  p25={qs[1]:.0f}  "
                  f"p50={qs[2]:.0f}  p75={qs[3]:.0f}  p95={qs[4]:.0f}")

    # ── HP ratio summary for TO_only_LOH (Paired side) ───────────────
    print("\n[STATS] HP ratio distribution for TO_only_LOH on Paired side:")
    to_only_paired_ratio = grps["TO_only_LOH"]["hp_ratio_paired"]
    mid_range = ((to_only_paired_ratio >= 0.2) & (to_only_paired_ratio <= 0.8)).sum()
    extreme = ((to_only_paired_ratio < 0.1) | (to_only_paired_ratio > 0.9)).sum()
    total_v = len(to_only_paired_ratio)
    print(f"       0.2 <= hp_ratio_paired <= 0.8 : {mid_range:,} / {total_v:,} ({100*mid_range/total_v:.1f}%)")
    print(f"       hp_ratio_paired < 0.1 or > 0.9: {extreme:,} / {total_v:,} ({100*extreme/total_v:.1f}%)")

    print("\n[STATS] HP ratio distribution for TO_only_LOH on TO side:")
    to_only_to_ratio = grps["TO_only_LOH"]["hp_ratio_to"]
    extreme_to = ((to_only_to_ratio < 0.1) | (to_only_to_ratio > 0.9)).sum()
    total_t = len(to_only_to_ratio)
    print(f"       hp_ratio_to < 0.1 or > 0.9: {extreme_to:,} / {total_t:,} ({100*extreme_to/total_t:.1f}%)")

    out = OUT_DIR / "concordance_depth_stats.tsv"
    stats_df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    print(f"\n[STATS] Saved -> {out}")

    return stats_df


# ── Conclusion printer ────────────────────────────────────────────────
def print_conclusion(df: pd.DataFrame, stats_df: pd.DataFrame):
    """Print structured conclusion based on analysis results."""
    print("\n" + "=" * 80)
    print("CONCLUSION SUMMARY")
    print("=" * 80)

    # Q1: Is TO-only LOH driven by low read depth?
    # Compare median effective_hp_reads_to for TO_only vs both_LOH
    row_t1 = stats_df[stats_df["test"].str.contains("TO_only vs both_LOH.*eff_hp_reads_TO")].iloc[0]
    row_t2 = stats_df[stats_df["test"].str.contains("TO_only vs neither.*eff_hp_reads_TO")].iloc[0]
    d_vs_both = row_t1["cohens_d"]
    d_vs_neither = row_t2["cohens_d"]

    print(f"\nQ1: Is TO-only LOH driven by low read depth?")
    print(f"    TO_only median eff_hp_reads_TO = {row_t1['median_A']:.0f}")
    print(f"    both_LOH median                = {row_t1['median_B']:.0f}")
    print(f"    neither  median                = {row_t2['median_B']:.0f}")
    print(f"    Cohen's d (vs both_LOH)   = {d_vs_both:.3f}")
    print(f"    Cohen's d (vs neither)    = {d_vs_neither:.3f}")
    if abs(d_vs_both) < 0.2 and abs(d_vs_neither) < 0.2:
        print("    --> CONCLUSION: Read depth difference is NEGLIGIBLE (|d| < 0.2).")
        print("        TO-only LOH is NOT primarily caused by low read depth.")
    elif abs(d_vs_both) < 0.5:
        print("    --> CONCLUSION: Read depth difference is SMALL (0.2 <= |d| < 0.5).")
        print("        Read depth plays a minor role but is not the dominant factor.")
    else:
        print("    --> CONCLUSION: Read depth difference is MODERATE-LARGE (|d| >= 0.5).")
        print("        Read depth is likely a significant contributor to TO-only LOH.")

    # Q2: Self-phasing bias evidence
    row_hp = stats_df[stats_df["test"].str.contains("TO_only vs both_LOH.*hp_ratio_Paired")].iloc[0]
    print(f"\nQ2: Self-phasing bias evidence?")
    print(f"    TO_only_LOH: Paired hp_ratio median = {row_hp['median_A']:.3f}")
    print(f"    both_LOH:    Paired hp_ratio median = {row_hp['median_B']:.3f}")
    print(f"    Cohen's d (hp_ratio_paired) = {row_hp['cohens_d']:.3f}")

    to_only = df[df["quadrant"] == "TO_only_LOH"]
    mid = ((to_only["hp_ratio_core_paired"] >= 0.2) & (to_only["hp_ratio_core_paired"] <= 0.8))
    n_valid = to_only["hp_ratio_core_paired"].notna().sum()
    print(f"    TO_only loci with 0.2 <= paired_hp_ratio <= 0.8: {mid.sum():,}/{n_valid:,} ({100*mid.sum()/n_valid:.1f}%)")
    if mid.sum() / n_valid > 0.5:
        print("    --> Most TO-only LOH loci show BALANCED Paired HP ratio,")
        print("        meaning the locus is genuinely not LOH under germline phasing,")
        print("        but TO self-phasing pushes HP ratio to extremes.")
        print("        This strongly supports SELF-PHASING BIAS as the primary cause.")

    # Q3: Most anomalous samples
    print(f"\nQ3: Per-sample TO-only LOH proportion:")
    for s in SAMPLE_ORDER:
        sub = df[df["sample"] == s]
        n_total = len(sub)
        n_to_only = (sub["quadrant"] == "TO_only_LOH").sum()
        if n_total > 0:
            print(f"    {s:20s}: {n_to_only:>6,} / {n_total:>6,} ({100*n_to_only/n_total:.1f}%)")

    # Q4: Depth vs HP ratio extreme quantification
    print(f"\nQ4: Read depth vs HP ratio extreme relationship:")
    for mode in ["to", "paired"]:
        col_d = f"effective_hp_reads_{mode}"
        col_r = f"hp_ratio_core_{mode}"
        valid = df[[col_d, col_r, "quadrant"]].dropna()
        # For TO: at low depth (<=30), what fraction has extreme HP ratio?
        low = valid[valid[col_d] <= 30]
        high = valid[valid[col_d] > 30]
        ext_low = ((low[col_r] < 0.1) | (low[col_r] > 0.9)).sum()
        ext_high = ((high[col_r] < 0.1) | (high[col_r] > 0.9)).sum()
        pct_low = 100 * ext_low / len(low) if len(low) > 0 else 0
        pct_high = 100 * ext_high / len(high) if len(high) > 0 else 0
        print(f"    {mode.upper():6s}: depth<=30 extreme ratio = {pct_low:.1f}%  |  depth>30 extreme ratio = {pct_high:.1f}%")

    print("\n" + "=" * 80)


# ── Main ──────────────────────────────────────────────────────────────
def main():
    df = load_data()
    fig1_quadrant_depth(df)
    fig2_hp_ratio_scatter(df)
    fig3_depth_vs_hpratio(df)
    fig4_per_sample_depth(df)
    fig5_depth_diff(df)
    stats_df = run_stats(df)
    print_conclusion(df, stats_df)
    print("\n[DONE] All outputs saved to:", OUT_DIR)


if __name__ == "__main__":
    main()
