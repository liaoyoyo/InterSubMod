#!/usr/bin/env python3
"""
Step 3: LOH Segment 內多位點空間分析
======================================
研究計劃：LOH 區域 Intermediate AF 與亞克隆偵測
目的：在 LOH segments 內，相鄰 variants 的 AF 一致性是否支持 subclonal LOH？
      Low AF-SD segments vs High AF-SD segments 的 NGroups 和 methylation 差異。

核心原理（TITAN/Battenberg）：
  - 多個相鄰 variants 一致 intermediate AF → segmental subclonal event
  - AF 在 segment 內一致 = purity/CN effect（非 subclone）
  - AF 在 segment 內變異大 = 可能混合事件（boundary or subclone transition）

輸入：all_region_rows.tsv.gz + LOH.bed files
輸出：
  - Per-segment statistics
  - Segment-level correlation plots
  - AF consistency analysis
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── Paths ──
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUTDIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_subclone_af")
FIGDIR = OUTDIR / "figures"
DATADIR = OUTDIR / "data"

LOH_BEDS = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "COLO829": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260317_colo829_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed",
    "H1437": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "H2009": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1937": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1954": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
}

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_SHORT = {
    "HCC1395": "HCC1395", "HCC1395_DORADO": "HCC1395D",
    "COLO829": "COLO829", "H1437": "H1437", "H2009": "H2009",
    "HCC1937": "HCC1937", "HCC1954": "HCC1954",
}


def cn_tier(cm):
    if pd.isna(cm):
        return "Unknown"
    if cm < 0.75:
        return "CN1"
    elif cm < 1.25:
        return "CN2"
    elif cm < 1.75:
        return "CN3"
    else:
        return "CN4+"


def load_loh_beds():
    """Load LOH.bed files into dict of DataFrames."""
    beds = {}
    for sample, path in LOH_BEDS.items():
        try:
            df = pd.read_csv(path, sep="\t", header=None,
                             names=["chr", "start", "end", "label", "het_rate",
                                    "dot", "thick_start", "thick_end", "rgb"])
            df["sample"] = sample
            df["segment_id"] = [f"{sample}_{i}" for i in range(len(df))]
            df["segment_size"] = df["end"] - df["start"]
            beds[sample] = df
            print(f"  {sample}: {len(df)} segments, median size={df['segment_size'].median():,.0f}bp")
        except Exception as e:
            print(f"  {sample}: ERROR loading LOH.bed: {e}")
    return beds


def load_data():
    """Load master dataset, TO mode, LOH regions only."""
    print("Loading master dataset...")
    cols = [
        "RegionID", "Chr", "Pos", "sample", "mode", "truth_label", "caller_af",
        "to_loh_bed_hit", "HPFineNGroups", "HPFineF",
        "Coverage_Multiple", "NumReads", "HP_Ratio",
        "AlleleDelta", "Quality_Score",
        "HP1FamilyN", "HP2FamilyN",
    ]
    df = pd.read_csv(MASTER, sep="\t", usecols=cols, low_memory=False)
    df_to = df[df["mode"] == "to"].copy()
    df_to["loh"] = df_to["to_loh_bed_hit"].astype(str).str.lower() == "true"
    df_to["cn_tier"] = df_to["Coverage_Multiple"].apply(cn_tier)
    return df_to


def assign_segments(df, beds):
    """Assign each variant to its LOH segment using vectorized merge_asof."""
    print("\nAssigning variants to LOH segments...")

    all_assignments = []
    for sample in SAMPLE_ORDER:
        if sample not in beds:
            continue

        bed = beds[sample]
        ds = df[(df["sample"] == sample) & (df["loh"] == True)].copy()

        if len(ds) == 0:
            continue

        # Vectorized: per-chromosome merge_asof (find segment containing each variant)
        ds["segment_id"] = None
        for chrom in ds["Chr"].unique():
            ds_chr = ds[ds["Chr"] == chrom].sort_values("Pos")
            bed_chr = bed[bed["chr"] == chrom].sort_values("start")

            if len(bed_chr) == 0:
                continue

            # For each variant, find the last segment whose start <= Pos
            idx = np.searchsorted(bed_chr["start"].values, ds_chr["Pos"].values, side="right") - 1
            valid = (idx >= 0) & (idx < len(bed_chr))

            for k, (ds_idx, pos) in enumerate(zip(ds_chr.index, ds_chr["Pos"].values)):
                if valid[k]:
                    bed_row = bed_chr.iloc[idx[k]]
                    if pos < bed_row["end"]:
                        ds.loc[ds_idx, "segment_id"] = bed_row["segment_id"]

        all_assignments.append(ds)
        n_assigned = ds["segment_id"].notna().sum()
        print(f"  {sample}: {n_assigned}/{len(ds)} variants assigned to segments")

    return pd.concat(all_assignments, ignore_index=True)


def compute_segment_stats(df_assigned, beds):
    """Compute per-segment statistics."""
    print("\nComputing per-segment statistics...")

    df_valid = df_assigned.dropna(subset=["segment_id"])

    # Only TP for biological analysis
    df_tp = df_valid[df_valid["truth_label"] == "TP"]

    records = []
    for seg_id, group in df_tp.groupby("segment_id"):
        n = len(group)
        if n < 2:  # Need at least 2 variants for SD
            continue

        sample = seg_id.rsplit("_", 1)[0]

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

        # Methylation
        mean_allele_delta = group["AlleleDelta"].mean()
        mean_hpfinef = group["HPFineF"].mean()

        # Coverage (CN proxy)
        mean_cm = group["Coverage_Multiple"].mean()
        dominant_cn = cn_tier(mean_cm)

        # Segment info from bed
        bed_match = None
        for s in SAMPLE_ORDER:
            if s in beds and seg_id in beds[s]["segment_id"].values:
                bed_match = beds[s][beds[s]["segment_id"] == seg_id].iloc[0]
                break

        seg_size = bed_match["segment_size"] if bed_match is not None else np.nan

        records.append({
            "segment_id": seg_id,
            "sample": sample,
            "n_variants": n,
            "segment_size_bp": seg_size,
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
    seg_df.to_csv(DATADIR / "step3_segment_statistics.tsv", sep="\t", index=False)
    print(f"  {len(seg_df)} segments with ≥2 TP variants")
    print(f"  Stats saved: {DATADIR / 'step3_segment_statistics.tsv'}")

    return seg_df


def fig11_segment_af_sd_vs_ngroups(seg_df):
    """Figure 11: Segment AF-SD vs mean NGroups, colored by CN tier."""
    print("\n=== Figure 11: Segment AF-SD vs NGroups ===")

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    cn_colors = {"CN1": "#E53935", "CN2": "#1E88E5", "CN3": "#43A047", "CN4+": "#FB8C00"}

    # Left: scatter
    ax = axes[0]
    for cn, color in cn_colors.items():
        sub = seg_df[seg_df["dominant_cn"] == cn]
        if len(sub) > 0:
            ax.scatter(sub["af_sd"], sub["mean_ngroups"], alpha=0.3, s=sub["n_variants"] * 2,
                       c=color, label=f"{cn} (n={len(sub)})", edgecolors="none")

    ax.set_xlabel("Within-segment AF Standard Deviation")
    ax.set_ylabel("Mean HPFineNGroups")
    ax.set_title("Segment AF-SD vs NGroups (size ∝ n_variants)", fontsize=12, fontweight="bold")
    ax.legend()
    ax.set_xlim(-0.02, seg_df["af_sd"].quantile(0.99))

    # Right: binned means with error bars (CN1 only)
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

        ax2.errorbar(bin_stats["mean_sd"], bin_stats["mean_ng"], yerr=bin_stats["sem_ng"],
                     fmt="o-", color="#E53935", capsize=5, markersize=8, linewidth=2)
        for _, row in bin_stats.iterrows():
            ax2.annotate(f"n={row['n']:.0f}", (row["mean_sd"], row["mean_ng"]),
                        textcoords="offset points", xytext=(5, 10), fontsize=8)

        # Spearman correlation
        rho, pval = scipy_stats.spearmanr(cn1["af_sd"], cn1["mean_ngroups"])
        ax2.text(0.05, 0.95, f"Spearman ρ={rho:.3f}\np={pval:.2e}",
                 transform=ax2.transAxes, fontsize=11, va="top",
                 bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.9))

    ax2.set_xlabel("Within-segment AF Standard Deviation")
    ax2.set_ylabel("Mean HPFineNGroups")
    ax2.set_title("Deletion LOH (CN≈1) — Binned AF-SD vs NGroups", fontsize=12, fontweight="bold")

    fig.suptitle("Step 3: Spatial Analysis — Within-Segment AF Variability vs Methylation Groups",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "11_segment_af_sd_vs_ngroups.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig12_segment_intermediate_proportion(seg_df):
    """Figure 12: Segments with high vs low intermediate AF proportion."""
    print("\n=== Figure 12: High vs Low intermediate AF segments ===")

    cn1 = seg_df[seg_df["dominant_cn"] == "CN1"].copy()
    if len(cn1) < 10:
        print("  Insufficient CN1 segments, skipping")
        return

    # Split into uniform (low intermediate) vs mixed (high intermediate)
    median_inter = cn1["pct_intermediate"].median()
    cn1["segment_type"] = np.where(cn1["pct_intermediate"] <= 10, "Uniform (≤10% inter)",
                                    np.where(cn1["pct_intermediate"] >= 50, "Mixed (≥50% inter)", "Middle"))

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    # Row 1: NGroups, AlleleDelta, HPFineF distributions
    for j, (metric, label) in enumerate([
        ("mean_ngroups", "Mean NGroups"),
        ("mean_allele_delta", "Mean AlleleDelta"),
        ("mean_hpfinef", "Mean HPFineF"),
    ]):
        ax = axes[0, j]
        for stype, color in [("Uniform (≤10% inter)", "#1976D2"), ("Mixed (≥50% inter)", "#FF9800")]:
            vals = cn1[cn1["segment_type"] == stype][metric].dropna()
            if len(vals) > 0:
                ax.hist(vals, bins=30, alpha=0.5, density=True, label=f"{stype}\n(n={len(vals)}, μ={vals.mean():.3f})",
                        color=color, edgecolor="none")

        ax.set_xlabel(label)
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)

        # Test
        uniform = cn1[cn1["segment_type"] == "Uniform (≤10% inter)"][metric].dropna()
        mixed = cn1[cn1["segment_type"] == "Mixed (≥50% inter)"][metric].dropna()
        if len(uniform) >= 5 and len(mixed) >= 5:
            _, pval = scipy_stats.mannwhitneyu(mixed, uniform, alternative="greater")
            d = (mixed.mean() - uniform.mean()) / np.sqrt((mixed.var() + uniform.var()) / 2 + 1e-10)
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
            ax.set_title(f"{label}\nd={d:.3f} {sig}", fontsize=11, fontweight="bold")

    # Row 2: Per-sample segment counts, AF-SD, n_variants
    ax = axes[1, 0]
    per_sample = cn1.groupby(["sample", "segment_type"]).size().unstack(fill_value=0)
    if "Uniform (≤10% inter)" in per_sample.columns and "Mixed (≥50% inter)" in per_sample.columns:
        per_sample[["Uniform (≤10% inter)", "Mixed (≥50% inter)"]].plot(
            kind="bar", ax=ax, color=["#1976D2", "#FF9800"], alpha=0.8)
    ax.set_xlabel("Sample")
    ax.set_ylabel("# Segments")
    ax.set_title("Segment Count by Type per Sample", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)

    ax2 = axes[1, 1]
    for stype, color in [("Uniform (≤10% inter)", "#1976D2"), ("Mixed (≥50% inter)", "#FF9800")]:
        vals = cn1[cn1["segment_type"] == stype]["af_sd"].dropna()
        if len(vals) > 0:
            ax2.hist(vals, bins=30, alpha=0.5, density=True, label=f"{stype} (n={len(vals)})",
                     color=color, edgecolor="none")
    ax2.set_xlabel("AF Standard Deviation")
    ax2.set_ylabel("Density")
    ax2.set_title("AF Variability Distribution", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=8)

    ax3 = axes[1, 2]
    for stype, color in [("Uniform (≤10% inter)", "#1976D2"), ("Mixed (≥50% inter)", "#FF9800")]:
        vals = cn1[cn1["segment_type"] == stype]["n_variants"].dropna()
        if len(vals) > 0:
            ax3.hist(vals, bins=30, alpha=0.5, density=True, label=f"{stype} (n={len(vals)})",
                     color=color, edgecolor="none")
    ax3.set_xlabel("# Variants in Segment")
    ax3.set_ylabel("Density")
    ax3.set_title("Variant Count per Segment", fontsize=11, fontweight="bold")
    ax3.legend(fontsize=8)

    fig.suptitle("Step 3: Deletion LOH (CN≈1) — Uniform vs Mixed Segments\n"
                 "Uniform = ≤10% intermediate AF; Mixed = ≥50% intermediate AF",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "12_segment_uniform_vs_mixed.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def fig13_per_sample_segment_consistency(seg_df):
    """Figure 13: Per-sample consistency of segment-level AF-SD vs NGroups correlation."""
    print("\n=== Figure 13: Per-sample segment consistency ===")

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    consistency_records = []

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        sub = seg_df[(seg_df["sample"] == sample) & (seg_df["dominant_cn"] == "CN1")]

        if len(sub) < 5:
            ax.text(0.5, 0.5, f"{SAMPLE_SHORT[sample]}\nInsufficient data (n={len(sub)})",
                    transform=ax.transAxes, ha="center", va="center")
            ax.set_title(SAMPLE_SHORT[sample])
            continue

        ax.scatter(sub["af_sd"], sub["mean_ngroups"], alpha=0.4, s=sub["n_variants"] * 3,
                   c="#E53935", edgecolors="none")

        rho, pval = scipy_stats.spearmanr(sub["af_sd"], sub["mean_ngroups"])
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

        ax.text(0.05, 0.95, f"ρ={rho:.3f} {sig}\nn={len(sub)}",
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

    axes[-1].axis("off")

    fig.suptitle("Step 3: Per-Sample Segment AF-SD vs NGroups (Deletion LOH CN≈1)\n"
                 "Each dot = one LOH segment; size ∝ variant count",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "13_per_sample_segment_consistency.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    cons_df = pd.DataFrame(consistency_records)
    cons_df.to_csv(DATADIR / "step3_per_sample_segment_correlation.tsv", sep="\t", index=False)
    print(f"  Consistency saved")

    return cons_df


def print_step3_summary(seg_df, cons_df):
    """Print Step 3 summary."""
    print("\n" + "=" * 80)
    print("STEP 3 SUMMARY: Spatial Analysis — Within-Segment AF Consistency")
    print("=" * 80)

    print(f"\n  Total segments with ≥2 TP: {len(seg_df)}")
    for cn in ["CN1", "CN2", "CN3", "CN4+"]:
        n = len(seg_df[seg_df["dominant_cn"] == cn])
        print(f"    {cn}: {n} segments")

    cn1 = seg_df[seg_df["dominant_cn"] == "CN1"]
    if len(cn1) >= 5:
        rho, pval = scipy_stats.spearmanr(cn1["af_sd"], cn1["mean_ngroups"])
        print(f"\n  CN1 overall: AF-SD vs NGroups Spearman ρ={rho:.3f}, p={pval:.2e}")

        # Uniform vs Mixed
        uniform = cn1[cn1["pct_intermediate"] <= 10]
        mixed = cn1[cn1["pct_intermediate"] >= 50]
        if len(uniform) >= 5 and len(mixed) >= 5:
            print(f"\n  Uniform segments (≤10% intermediate): n={len(uniform)}")
            print(f"    Mean NGroups: {uniform['mean_ngroups'].mean():.3f}")
            print(f"    Mean AlleleDelta: {uniform['mean_allele_delta'].mean():.4f}")
            print(f"    Mean AF-SD: {uniform['af_sd'].mean():.4f}")
            print(f"  Mixed segments (≥50% intermediate): n={len(mixed)}")
            print(f"    Mean NGroups: {mixed['mean_ngroups'].mean():.3f}")
            print(f"    Mean AlleleDelta: {mixed['mean_allele_delta'].mean():.4f}")
            print(f"    Mean AF-SD: {mixed['af_sd'].mean():.4f}")

    if cons_df is not None and len(cons_df) > 0:
        print("\n  Per-sample segment-level Spearman ρ (AF-SD vs NGroups):")
        n_pos = 0
        n_sig = 0
        for _, row in cons_df.iterrows():
            sig = "***" if row["spearman_p"] < 0.001 else "**" if row["spearman_p"] < 0.01 else "*" if row["spearman_p"] < 0.05 else "ns"
            print(f"    {row['sample']}: ρ={row['spearman_rho']:.3f} {sig} (n={row['n_segments']})")
            if row["spearman_rho"] > 0:
                n_pos += 1
            if row["spearman_p"] < 0.05:
                n_sig += 1
        print(f"  Direction consistency: {n_pos}/{len(cons_df)} positive")
        print(f"  Significant: {n_sig}/{len(cons_df)}")


def main():
    # Load LOH.bed files
    print("Loading LOH.bed files...")
    beds = load_loh_beds()

    # Load master data
    df = load_data()

    # Assign to segments
    df_assigned = assign_segments(df, beds)

    # Compute segment stats
    seg_df = compute_segment_stats(df_assigned, beds)

    # Generate figures
    fig11_segment_af_sd_vs_ngroups(seg_df)
    fig12_segment_intermediate_proportion(seg_df)
    cons_df = fig13_per_sample_segment_consistency(seg_df)

    # Summary
    print_step3_summary(seg_df, cons_df)


if __name__ == "__main__":
    main()
