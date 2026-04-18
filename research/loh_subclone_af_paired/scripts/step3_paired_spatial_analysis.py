#!/usr/bin/env python3
"""
Step 3: Paired Mode LOH Segment Spatial Analysis
==================================================
Parallel to TO step3 — segment-level AF consistency and NGroups correlation.

Core principle (TITAN/Battenberg):
  - Multiple adjacent variants with consistent intermediate AF → segmental subclonal event
  - AF-SD within segment correlates with NGroups if subclonal structure is real
  - Uniform (low AF-SD) vs Mixed (high AF-SD) segments differ in methylation groups

Input: all_region_rows.tsv.gz (master dataset, paired mode) + LOH.bed
Output:
  - p11: Segment AF-SD vs NGroups scatter + binned means
  - p12: Uniform vs Mixed segment comparison
  - p13: Per-sample segment-level Spearman correlation
  - step3_paired_segment_statistics.tsv
  - step3_paired_per_sample_segment_correlation.tsv
"""

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats
from utils import (
    MASTER, FIGDIR, DATADIR, SAMPLE_ORDER, SAMPLE_SHORT,
    AF_COLORS, cn_tier, load_loh_beds, load_paired_data
)
import warnings
warnings.filterwarnings("ignore")

FIGDIR.mkdir(parents=True, exist_ok=True)
DATADIR.mkdir(parents=True, exist_ok=True)


def annotate_loh_and_assign_segments(df, beds):
    """Annotate variants with LOH.bed hit AND assign segment_id simultaneously.

    More efficient than two-pass: annotate + assign in a single loop.
    """
    print("\nAnnotating LOH.bed and assigning segments...")

    df["loh_bed_hit"] = False
    df["segment_id"] = None

    for sample, bed in beds.items():
        mask = df["sample"] == sample
        if mask.sum() == 0:
            continue

        n_assigned = 0
        for chrom in df.loc[mask, "Chr"].unique():
            chr_mask = mask & (df["Chr"] == chrom)
            if chr_mask.sum() == 0:
                continue

            bed_chr = bed[bed["chr"] == chrom].sort_values("start")
            if len(bed_chr) == 0:
                continue

            starts = bed_chr["start"].values
            ends = bed_chr["end"].values
            seg_ids = bed_chr["segment_id"].values
            positions = df.loc[chr_mask, "Pos"].values

            idx = np.searchsorted(starts, positions, side="right") - 1
            valid = (idx >= 0) & (idx < len(starts))
            in_segment = valid & (positions < ends[np.clip(idx, 0, len(ends) - 1)])

            chr_indices = df.index[chr_mask]
            df.loc[chr_indices[in_segment], "loh_bed_hit"] = True

            # Assign segment IDs for variants inside LOH segments
            for k in np.where(in_segment)[0]:
                df.loc[chr_indices[k], "segment_id"] = seg_ids[idx[k]]
                n_assigned += 1

        print(f"  {sample}: {n_assigned} variants assigned to LOH segments")

    n_hit = df["loh_bed_hit"].sum()
    n_seg = df["segment_id"].notna().sum()
    print(f"  Total LOH hits: {n_hit:,}, with segment: {n_seg:,}")

    return df


def compute_segment_stats(df):
    """Compute per-segment statistics from assigned data."""
    print("\nComputing per-segment statistics...")

    df_valid = df.dropna(subset=["segment_id"])
    df_tp = df_valid[df_valid["truth_label"] == "TP"]

    records = []
    for seg_id, group in df_tp.groupby("segment_id"):
        n = len(group)
        if n < 2:  # Need >=2 variants for SD
            continue

        sample = group["sample"].iloc[0]

        # AF statistics
        af_mean = group["caller_af"].mean()
        af_sd = group["caller_af"].std()
        af_range = group["caller_af"].max() - group["caller_af"].min()

        # Intermediate AF proportion
        n_inter = ((group["caller_af"] >= 0.1) & (group["caller_af"] <= 0.4) |
                   (group["caller_af"] >= 0.6) & (group["caller_af"] <= 0.9)).sum()
        pct_inter = 100 * n_inter / n

        # NGroups
        mean_ngroups = group["HPFineNGroups"].mean()
        pct_ngroups_ge2 = 100 * (group["HPFineNGroups"] >= 2).sum() / n

        # Methylation features
        mean_allele_delta = group["AlleleDelta"].mean()
        mean_hpfinef = group["HPFineF"].mean()

        # Coverage (CN proxy)
        mean_cm = group["Coverage_Multiple"].mean()
        dominant_cn = cn_tier(mean_cm)

        records.append({
            "segment_id": seg_id,
            "sample": sample,
            "n_variants": n,
            "af_mean": af_mean,
            "af_sd": af_sd,
            "af_range": af_range,
            "pct_intermediate": pct_inter,
            "mean_ngroups": mean_ngroups,
            "pct_ngroups_ge2": pct_ngroups_ge2,
            "mean_allele_delta": mean_allele_delta,
            "mean_hpfinef": mean_hpfinef,
            "mean_cm": mean_cm,
            "dominant_cn": dominant_cn,
            "mean_numreads": group["NumReads"].mean(),
        })

    seg_df = pd.DataFrame(records)
    seg_df.to_csv(DATADIR / "step3_paired_segment_statistics.tsv", sep="\t", index=False)
    print(f"  {len(seg_df)} segments with >=2 TP variants")
    print(f"  Saved: {DATADIR / 'step3_paired_segment_statistics.tsv'}")

    return seg_df


def fig11_segment_af_sd_vs_ngroups(seg_df):
    """p11: Segment AF-SD vs mean NGroups, colored by CN tier."""
    print("\n=== Figure p11: Segment AF-SD vs NGroups ===")

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    cn_colors = {"CN1": "#E53935", "CN2": "#1E88E5", "CN3": "#43A047", "CN4+": "#FB8C00"}

    # Left: scatter by CN tier
    ax = axes[0]
    for cn, color in cn_colors.items():
        sub = seg_df[seg_df["dominant_cn"] == cn]
        if len(sub) > 0:
            ax.scatter(sub["af_sd"], sub["mean_ngroups"], alpha=0.3,
                       s=sub["n_variants"] * 2,
                       c=color, label=f"{cn} (n={len(sub)})", edgecolors="none")

    ax.set_xlabel("Within-segment AF Standard Deviation")
    ax.set_ylabel("Mean HPFineNGroups")
    ax.set_title("Segment AF-SD vs NGroups (size proportional to n_variants)", fontsize=11, fontweight="bold")
    ax.legend()
    ax.set_xlim(-0.02, seg_df["af_sd"].quantile(0.99))

    # Right: binned means for CN1
    ax2 = axes[1]
    cn1 = seg_df[seg_df["dominant_cn"] == "CN1"].copy()
    if len(cn1) >= 10:
        cn1["af_sd_bin"] = pd.qcut(cn1["af_sd"], q=5, duplicates="drop")
        bin_stats = cn1.groupby("af_sd_bin", observed=True).agg(
            mean_ng=("mean_ngroups", "mean"),
            sem_ng=("mean_ngroups", "sem"),
            n=("mean_ngroups", "count"),
            mean_sd=("af_sd", "mean"),
        ).reset_index()

        ax2.errorbar(bin_stats["mean_sd"], bin_stats["mean_ng"],
                     yerr=bin_stats["sem_ng"],
                     fmt="o-", color="#E53935", capsize=5, markersize=8, linewidth=2)
        for _, row in bin_stats.iterrows():
            ax2.annotate(f"n={row['n']:.0f}", (row["mean_sd"], row["mean_ng"]),
                         textcoords="offset points", xytext=(5, 10), fontsize=8)

        rho, pval = scipy_stats.spearmanr(cn1["af_sd"], cn1["mean_ngroups"])
        ax2.text(0.05, 0.95, f"Spearman rho={rho:.3f}\np={pval:.2e}",
                 transform=ax2.transAxes, fontsize=11, va="top",
                 bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.9))

    ax2.set_xlabel("Within-segment AF Standard Deviation")
    ax2.set_ylabel("Mean HPFineNGroups")
    ax2.set_title("Deletion LOH (CN~1) - Binned AF-SD vs NGroups", fontsize=11, fontweight="bold")

    fig.suptitle("Paired Mode Step 3: Within-Segment AF Variability vs Methylation Groups",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "p11_paired_segment_af_sd_vs_ngroups.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig12_segment_uniform_vs_mixed(seg_df):
    """p12: Uniform vs Mixed segment comparison."""
    print("\n=== Figure p12: Uniform vs Mixed segments ===")

    cn1 = seg_df[seg_df["dominant_cn"] == "CN1"].copy()
    if len(cn1) < 10:
        print("  Insufficient CN1 segments, skipping")
        return

    cn1["segment_type"] = np.where(
        cn1["pct_intermediate"] <= 10, "Uniform (<=10% inter)",
        np.where(cn1["pct_intermediate"] >= 50, "Mixed (>=50% inter)", "Middle")
    )

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    # Row 1: NGroups, AlleleDelta, HPFineF
    for j, (metric, label) in enumerate([
        ("mean_ngroups", "Mean NGroups"),
        ("mean_allele_delta", "Mean AlleleDelta"),
        ("mean_hpfinef", "Mean HPFineF"),
    ]):
        ax = axes[0, j]
        for stype, color in [("Uniform (<=10% inter)", "#1976D2"),
                              ("Mixed (>=50% inter)", "#FF9800")]:
            vals = cn1[cn1["segment_type"] == stype][metric].dropna()
            if len(vals) > 0:
                ax.hist(vals, bins=30, alpha=0.5, density=True,
                        label=f"{stype}\n(n={len(vals)}, mu={vals.mean():.3f})",
                        color=color, edgecolor="none")

        ax.set_xlabel(label)
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)

        uniform = cn1[cn1["segment_type"] == "Uniform (<=10% inter)"][metric].dropna()
        mixed = cn1[cn1["segment_type"] == "Mixed (>=50% inter)"][metric].dropna()
        if len(uniform) >= 5 and len(mixed) >= 5:
            _, pval = scipy_stats.mannwhitneyu(mixed, uniform, alternative="greater")
            d = (mixed.mean() - uniform.mean()) / np.sqrt((mixed.var() + uniform.var()) / 2 + 1e-10)
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
            ax.set_title(f"{label}\nd={d:.3f} {sig}", fontsize=11, fontweight="bold")

    # Row 2: per-sample counts, AF-SD distribution, n_variants
    ax = axes[1, 0]
    per_sample = cn1.groupby(["sample", "segment_type"]).size().unstack(fill_value=0)
    cols_to_plot = [c for c in ["Uniform (<=10% inter)", "Mixed (>=50% inter)"] if c in per_sample.columns]
    if len(cols_to_plot) > 0:
        per_sample[cols_to_plot].plot(
            kind="bar", ax=ax, color=["#1976D2", "#FF9800"][:len(cols_to_plot)], alpha=0.8)
    ax.set_xlabel("Sample")
    ax.set_ylabel("# Segments")
    ax.set_title("Segment Count by Type per Sample", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)

    ax2 = axes[1, 1]
    for stype, color in [("Uniform (<=10% inter)", "#1976D2"),
                          ("Mixed (>=50% inter)", "#FF9800")]:
        vals = cn1[cn1["segment_type"] == stype]["af_sd"].dropna()
        if len(vals) > 0:
            ax2.hist(vals, bins=30, alpha=0.5, density=True,
                     label=f"{stype} (n={len(vals)})", color=color, edgecolor="none")
    ax2.set_xlabel("AF Standard Deviation")
    ax2.set_ylabel("Density")
    ax2.set_title("AF Variability Distribution", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=8)

    ax3 = axes[1, 2]
    for stype, color in [("Uniform (<=10% inter)", "#1976D2"),
                          ("Mixed (>=50% inter)", "#FF9800")]:
        vals = cn1[cn1["segment_type"] == stype]["n_variants"].dropna()
        if len(vals) > 0:
            ax3.hist(vals, bins=30, alpha=0.5, density=True,
                     label=f"{stype} (n={len(vals)})", color=color, edgecolor="none")
    ax3.set_xlabel("# Variants in Segment")
    ax3.set_ylabel("Density")
    ax3.set_title("Variant Count per Segment", fontsize=11, fontweight="bold")
    ax3.legend(fontsize=8)

    fig.suptitle("Paired Mode: Deletion LOH (CN~1) - Uniform vs Mixed Segments\n"
                 "Uniform = <=10% intermediate AF; Mixed = >=50% intermediate AF",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p12_paired_segment_uniform_vs_mixed.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig13_per_sample_segment_consistency(seg_df):
    """p13: Per-sample Spearman correlation of AF-SD vs NGroups."""
    print("\n=== Figure p13: Per-sample segment consistency ===")

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes_flat = axes.flatten()

    consistency_records = []

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes_flat[i]
        sub = seg_df[(seg_df["sample"] == sample) & (seg_df["dominant_cn"] == "CN1")]

        if len(sub) < 5:
            ax.text(0.5, 0.5, f"{SAMPLE_SHORT[sample]}\nInsufficient data (n={len(sub)})",
                    transform=ax.transAxes, ha="center", va="center")
            ax.set_title(SAMPLE_SHORT[sample])
            consistency_records.append({
                "sample": sample,
                "n_segments": len(sub),
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
                "direction": "insufficient",
            })
            continue

        ax.scatter(sub["af_sd"], sub["mean_ngroups"], alpha=0.4,
                   s=sub["n_variants"] * 3, c="#E53935", edgecolors="none")

        rho, pval = scipy_stats.spearmanr(sub["af_sd"], sub["mean_ngroups"])
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

        ax.text(0.05, 0.95, f"rho={rho:.3f} {sig}\nn={len(sub)}",
                transform=ax.transAxes, fontsize=10, va="top",
                bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.9))

        ax.set_xlabel("AF SD")
        ax.set_ylabel("Mean NGroups")
        ax.set_title(SAMPLE_SHORT[sample], fontsize=12, fontweight="bold")

        consistency_records.append({
            "sample": sample,
            "n_segments": len(sub),
            "spearman_rho": rho,
            "spearman_p": pval,
            "direction": "positive" if rho > 0 else "negative",
        })

    axes_flat[-1].axis("off")

    fig.suptitle("Paired Mode: Per-Sample Segment AF-SD vs NGroups (Deletion LOH CN~1)\n"
                 "Each dot = one LOH segment; size proportional to variant count",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "p13_paired_per_sample_segment_consistency.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    cons_df = pd.DataFrame(consistency_records)
    cons_df.to_csv(DATADIR / "step3_paired_per_sample_segment_correlation.tsv",
                   sep="\t", index=False)
    print(f"  Saved: {DATADIR / 'step3_paired_per_sample_segment_correlation.tsv'}")

    return cons_df


def print_step3_summary(seg_df, cons_df):
    """Print Step 3 summary."""
    print("\n" + "=" * 80)
    print("STEP 3 SUMMARY: Paired Mode Spatial Analysis")
    print("=" * 80)

    print(f"\n  Total segments with >=2 TP: {len(seg_df)}")
    for cn in ["CN1", "CN2", "CN3", "CN4+"]:
        n = len(seg_df[seg_df["dominant_cn"] == cn])
        print(f"    {cn}: {n} segments")

    cn1 = seg_df[seg_df["dominant_cn"] == "CN1"]
    if len(cn1) >= 5:
        rho, pval = scipy_stats.spearmanr(cn1["af_sd"], cn1["mean_ngroups"])
        print(f"\n  CN1 overall: AF-SD vs NGroups Spearman rho={rho:.3f}, p={pval:.2e}")

        uniform = cn1[cn1["pct_intermediate"] <= 10]
        mixed = cn1[cn1["pct_intermediate"] >= 50]
        if len(uniform) >= 5 and len(mixed) >= 5:
            print(f"\n  Uniform segments (<=10% intermediate): n={len(uniform)}")
            print(f"    Mean NGroups: {uniform['mean_ngroups'].mean():.3f}")
            print(f"    Mean AlleleDelta: {uniform['mean_allele_delta'].mean():.4f}")
            print(f"    Mean AF-SD: {uniform['af_sd'].mean():.4f}")
            print(f"  Mixed segments (>=50% intermediate): n={len(mixed)}")
            print(f"    Mean NGroups: {mixed['mean_ngroups'].mean():.3f}")
            print(f"    Mean AlleleDelta: {mixed['mean_allele_delta'].mean():.4f}")
            print(f"    Mean AF-SD: {mixed['af_sd'].mean():.4f}")

    if cons_df is not None and len(cons_df) > 0:
        valid = cons_df.dropna(subset=["spearman_rho"])
        print("\n  Per-sample segment-level Spearman rho (AF-SD vs NGroups):")
        n_pos = 0
        n_sig = 0
        for _, row in valid.iterrows():
            sig = "***" if row["spearman_p"] < 0.001 else "**" if row["spearman_p"] < 0.01 else "*" if row["spearman_p"] < 0.05 else "ns"
            print(f"    {row['sample']}: rho={row['spearman_rho']:.3f} {sig} (n={row['n_segments']})")
            if row["spearman_rho"] > 0:
                n_pos += 1
            if row["spearman_p"] < 0.05:
                n_sig += 1
        print(f"  Direction consistency: {n_pos}/{len(valid)} positive")
        print(f"  Significant: {n_sig}/{len(valid)}")

        # H3p verdict
        if n_pos >= 5:
            print(f"\n  >>> H3p POSITIVE: {n_pos}/{len(valid)} samples show positive rho")
        else:
            print(f"\n  >>> H3p NEGATIVE: only {n_pos}/{len(valid)} positive (need >=5)")


def main():
    print("=" * 70)
    print("Step 3: Paired Mode LOH Segment Spatial Analysis")
    print("=" * 70)

    # Load paired data
    _, df_paired = load_paired_data()

    # Load LOH.bed
    print("\nLoading LOH.bed files...")
    beds = load_loh_beds()

    # Annotate and assign segments in one pass
    df_paired = annotate_loh_and_assign_segments(df_paired, beds)

    # Compute segment stats
    seg_df = compute_segment_stats(df_paired)

    # Generate figures
    fig11_segment_af_sd_vs_ngroups(seg_df)
    fig12_segment_uniform_vs_mixed(seg_df)
    cons_df = fig13_per_sample_segment_consistency(seg_df)

    # Summary
    print_step3_summary(seg_df, cons_df)

    print("\n" + "=" * 70)
    print("Step 3 complete!")
    print(f"  Figures: {FIGDIR}/p11-p13_*.png")
    print(f"  Data: {DATADIR}/step3_*.tsv")
    print("=" * 70)


if __name__ == "__main__":
    main()
