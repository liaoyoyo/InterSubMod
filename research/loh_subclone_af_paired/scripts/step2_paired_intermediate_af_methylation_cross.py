#!/usr/bin/env python3
"""
Step 2: Paired Mode — Intermediate AF × NGroups × Methylation Cross-Analysis
=============================================================================
Parallel to TO step2. Tests H1p and H4p.

Core analysis:
  1. NGroups by AF class × CN tier (p06)
  2. Per-sample Mann-Whitney: Intermediate > Extreme NGroups in CN1 (p07)
  3. NumReads-controlled analysis (p08)
  4. Methylation features Cohen's d (p09)
  5. Per-sample consistency (p10)
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
    FIGDIR, DATADIR, SAMPLE_ORDER, SAMPLE_SHORT, AF_COLORS,
    NR_BINS, NR_LABELS, classify_af, cn_tier, load_loh_beds,
    annotate_loh_bed, load_paired_data
)
import warnings
warnings.filterwarnings("ignore")


def fig06_ngroups_by_af_class_cn_tier(df):
    """p06: NGroups distribution by AF class in each CN tier."""
    print("\n=== Figure p06: NGroups by AF class × CN tier ===")

    df_loh_tp = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP")]
    cn_order = ["CN1", "CN2", "CN3", "CN4+"]
    af_classes = ["Extreme", "Near-half", "Intermediate"]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    records = []

    for ax_idx, cn_val in enumerate(cn_order):
        ax = axes.flatten()[ax_idx]
        sub = df_loh_tp[df_loh_tp["cn_tier"] == cn_val]

        ngroups_range = [0, 1, 2, 3, 4]
        x = np.arange(len(ngroups_range))
        width = 0.25

        for j, afc in enumerate(af_classes):
            sub_af = sub[sub["af_class"] == afc]
            if len(sub_af) == 0:
                continue
            counts = sub_af["HPFineNGroups"].value_counts().reindex(ngroups_range, fill_value=0)
            pcts = 100 * counts / counts.sum()
            ax.bar(x + (j - 1) * width, pcts.values, width,
                   label=f"{afc} (n={len(sub_af):,})", color=AF_COLORS[afc], alpha=0.8)

            records.append({
                "cn_tier": cn_val, "af_class": afc, "n": len(sub_af),
                "mean_ngroups": sub_af["HPFineNGroups"].mean(),
                "median_ngroups": sub_af["HPFineNGroups"].median(),
                "pct_ngroups_ge2": 100 * (sub_af["HPFineNGroups"] >= 2).sum() / len(sub_af),
                "mean_numreads": sub_af["NumReads"].mean(),
            })

        ax.set_xlabel("HPFineNGroups")
        ax.set_ylabel("% of variants")
        ax.set_title(f"{cn_val} (LOH TP)", fontsize=12, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(ngroups_range)
        ax.legend(fontsize=9)

    fig.suptitle("Step 2 (Paired): HPFineNGroups by AF Class in LOH TP\n"
                 "H1p: Intermediate AF → more NGroups?", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "p06_paired_ngroups_by_af_class_cn_tier.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    ngroups_stats = pd.DataFrame(records)
    ngroups_stats.to_csv(DATADIR / "step2_paired_ngroups_by_af_class.tsv", sep="\t", index=False)
    return ngroups_stats


def fig07_ngroups_per_sample(df):
    """p07: Per-sample Mann-Whitney test, CN1 LOH TP."""
    print("\n=== Figure p07: Per-sample NGroups CN1 ===")

    df_cn1_tp = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP") & (df["cn_tier"] == "CN1")]

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Left: mean NGroups per sample
    ax = axes[0]
    bar_records = []
    for sample in SAMPLE_ORDER:
        ds = df_cn1_tp[df_cn1_tp["sample"] == sample]
        for afc in ["Extreme", "Intermediate"]:
            sub = ds[ds["af_class"] == afc]
            if len(sub) < 5:
                continue
            bar_records.append({
                "sample": SAMPLE_SHORT[sample], "af_class": afc,
                "n": len(sub), "mean_ngroups": sub["HPFineNGroups"].mean(),
                "sem_ngroups": sub["HPFineNGroups"].sem(),
            })

    bar_df = pd.DataFrame(bar_records)
    if len(bar_df) > 0:
        samples_list = bar_df["sample"].unique()
        x = np.arange(len(samples_list))
        width = 0.35
        ext = bar_df[bar_df["af_class"] == "Extreme"]
        inter = bar_df[bar_df["af_class"] == "Intermediate"]
        if len(ext) > 0:
            ax.bar(x[:len(ext)] - width/2, ext["mean_ngroups"].values, width,
                   yerr=ext["sem_ngroups"].values, label="Extreme AF",
                   color="#1976D2", alpha=0.8, capsize=3)
        if len(inter) > 0:
            ax.bar(x[:len(inter)] + width/2, inter["mean_ngroups"].values, width,
                   yerr=inter["sem_ngroups"].values, label="Intermediate AF",
                   color="#FF9800", alpha=0.8, capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels(samples_list, fontsize=9)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Mean HPFineNGroups")
    ax.set_title("Paired: Deletion LOH (CN≈1) TP — Mean NGroups", fontsize=12, fontweight="bold")
    ax.legend()

    # Right: MW test per sample
    ax2 = axes[1]
    test_records = []
    for sample in SAMPLE_ORDER:
        ds = df_cn1_tp[df_cn1_tp["sample"] == sample]
        extreme = ds[ds["af_class"] == "Extreme"]["HPFineNGroups"].dropna()
        intermediate = ds[ds["af_class"] == "Intermediate"]["HPFineNGroups"].dropna()

        if len(extreme) >= 5 and len(intermediate) >= 5:
            stat, pval = scipy_stats.mannwhitneyu(intermediate, extreme, alternative="greater")
            n1, n2 = len(intermediate), len(extreme)
            r = 1 - 2 * stat / (n1 * n2)
            test_records.append({
                "sample": SAMPLE_SHORT[sample],
                "n_extreme": len(extreme), "n_intermediate": len(intermediate),
                "mean_ng_extreme": extreme.mean(), "mean_ng_intermediate": intermediate.mean(),
                "mw_U": stat, "mw_p": pval, "rank_biserial_r": r,
                "mean_nr_extreme": ds[ds["af_class"] == "Extreme"]["NumReads"].mean(),
                "mean_nr_intermediate": ds[ds["af_class"] == "Intermediate"]["NumReads"].mean(),
            })

    test_df = pd.DataFrame(test_records)
    if len(test_df) > 0:
        colors = ["#4CAF50" if p < 0.05 else "#F44336" for p in test_df["mw_p"]]
        ax2.barh(range(len(test_df)), -np.log10(test_df["mw_p"].clip(1e-300)), color=colors, alpha=0.8)
        ax2.axvline(-np.log10(0.05), color="red", ls="--", alpha=0.7, label="p=0.05")
        ax2.set_yticks(range(len(test_df)))
        ax2.set_yticklabels(test_df["sample"])
        ax2.set_xlabel("-log10(p-value)")
        ax2.set_title("Mann-Whitney: Intermediate > Extreme NGroups?", fontsize=12, fontweight="bold")
        ax2.legend()
        for idx, row in test_df.iterrows():
            ax2.text(-np.log10(row["mw_p"]) + 0.3, idx,
                     f"r={row['rank_biserial_r']:.3f}", fontsize=8, va="center")

    fig.suptitle("Step 2 (Paired): NGroups in Deletion LOH (CN≈1)", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "p07_paired_ngroups_per_sample_cn1.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    test_df.to_csv(DATADIR / "step2_paired_ngroups_mw_test_cn1.tsv", sep="\t", index=False)
    return test_df


def fig08_numreads_controlled(df):
    """p08: NR-bin controlled NGroups analysis."""
    print("\n=== Figure p08: NumReads-controlled analysis ===")

    df_cn1_tp = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP") & (df["cn_tier"] == "CN1")].copy()
    df_cn1_tp = df_cn1_tp.dropna(subset=["HPFineNGroups", "NumReads"])

    fig, axes = plt.subplots(1, len(NR_BINS), figsize=(22, 6), sharey=True)
    all_records = []

    for ax_idx, (nr_lo, nr_hi) in enumerate(NR_BINS):
        ax = axes[ax_idx]
        sub = df_cn1_tp[(df_cn1_tp["NumReads"] >= nr_lo) & (df_cn1_tp["NumReads"] < nr_hi)]

        extreme = sub[sub["af_class"] == "Extreme"]
        intermediate = sub[sub["af_class"] == "Intermediate"]

        ngroups_range = [0, 1, 2, 3, 4]
        x = np.arange(len(ngroups_range))
        width = 0.35

        if len(extreme) > 0:
            counts_e = extreme["HPFineNGroups"].value_counts().reindex(ngroups_range, fill_value=0)
            pcts_e = 100 * counts_e / counts_e.sum()
            ax.bar(x - width/2, pcts_e.values, width, label=f"Extreme (n={len(extreme):,})",
                   color="#1976D2", alpha=0.8)

        if len(intermediate) > 0:
            counts_i = intermediate["HPFineNGroups"].value_counts().reindex(ngroups_range, fill_value=0)
            pcts_i = 100 * counts_i / counts_i.sum()
            ax.bar(x + width/2, pcts_i.values, width, label=f"Intermediate (n={len(intermediate):,})",
                   color="#FF9800", alpha=0.8)

        ax.set_xlabel("HPFineNGroups")
        ax.set_ylabel("% of variants" if ax_idx == 0 else "")
        ax.set_title(f"NR={NR_LABELS[ax_idx]}", fontsize=11, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(ngroups_range)
        ax.legend(fontsize=7)

        # MW test
        if len(extreme) >= 5 and len(intermediate) >= 5:
            stat, pval = scipy_stats.mannwhitneyu(
                intermediate["HPFineNGroups"], extreme["HPFineNGroups"], alternative="greater")
            n1, n2 = len(intermediate), len(extreme)
            r = 1 - 2 * stat / (n1 * n2)
            all_records.append({
                "nr_bin": NR_LABELS[ax_idx],
                "n_extreme": len(extreme), "n_intermediate": len(intermediate),
                "mean_ng_extreme": extreme["HPFineNGroups"].mean(),
                "mean_ng_intermediate": intermediate["HPFineNGroups"].mean(),
                "mw_p": pval, "rank_biserial_r": r,
            })

    fig.suptitle("Step 2 (Paired): NGroups by AF Class — NumReads Controlled (CN1 LOH TP)\n"
                 "H4p: Effect persists after NR control?", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIGDIR / "p08_paired_ngroups_numreads_controlled.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    nr_df = pd.DataFrame(all_records)
    nr_df.to_csv(DATADIR / "step2_paired_ngroups_numreads_controlled.tsv", sep="\t", index=False)
    return nr_df


def fig09_methylation_features(df):
    """p09: Methylation features comparison."""
    print("\n=== Figure p09: Methylation features ===")

    features = ["HPFineF", "AlleleDelta", "PairwiseMeanDist", "hp_imbalance", "CramersV", "Quality_Score"]
    cn_tiers = ["CN1", "CN2"]

    records = []
    for feat in features:
        for cn_val in cn_tiers:
            sub = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP") & (df["cn_tier"] == cn_val)]
            extreme = sub[sub["af_class"] == "Extreme"][feat].dropna()
            intermediate = sub[sub["af_class"] == "Intermediate"][feat].dropna()

            if len(extreme) >= 5 and len(intermediate) >= 5:
                pooled_std = np.sqrt((extreme.var() + intermediate.var()) / 2)
                d = (intermediate.mean() - extreme.mean()) / pooled_std if pooled_std > 0 else 0
                _, pval = scipy_stats.mannwhitneyu(intermediate, extreme, alternative="two-sided")
                records.append({
                    "feature": feat, "cn_tier": cn_val,
                    "n_extreme": len(extreme), "n_intermediate": len(intermediate),
                    "mean_extreme": extreme.mean(), "mean_intermediate": intermediate.mean(),
                    "cohen_d": d, "mw_p": pval,
                })

    feat_df = pd.DataFrame(records)
    feat_df.to_csv(DATADIR / "step2_paired_methylation_features_comparison.tsv", sep="\t", index=False)
    print(f"  Saved: {DATADIR / 'step2_paired_methylation_features_comparison.tsv'}")

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for idx, feat in enumerate(features):
        ax = axes.flatten()[idx]
        feat_data = feat_df[feat_df["feature"] == feat]
        if len(feat_data) > 0:
            x = np.arange(len(feat_data))
            ax.bar(x, feat_data["cohen_d"].values, color=["#E91E63", "#9C27B0"][:len(feat_data)], alpha=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(feat_data["cn_tier"].values)
            ax.axhline(0, color="gray", ls="-", alpha=0.3)
            for i, row in feat_data.iterrows():
                ax.text(x[list(feat_data.index).index(i)], row["cohen_d"] + 0.02,
                        f"d={row['cohen_d']:.3f}", ha="center", fontsize=8)
        ax.set_title(feat, fontsize=11, fontweight="bold")
        ax.set_ylabel("Cohen's d")

    fig.suptitle("Step 2 (Paired): Methylation Feature Effect Sizes (Intermediate vs Extreme AF)\n"
                 "Positive d = Intermediate > Extreme", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIGDIR / "p09_paired_methylation_features_by_af_class.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return feat_df


def fig10_per_sample_consistency(df):
    """p10: Per-sample consistency."""
    print("\n=== Figure p10: Per-sample consistency ===")

    df_cn1_tp = df[(df["loh_bed_hit"] == True) & (df["truth_label"] == "TP") & (df["cn_tier"] == "CN1")]

    records = []
    for sample in SAMPLE_ORDER:
        ds = df_cn1_tp[df_cn1_tp["sample"] == sample]
        extreme = ds[ds["af_class"] == "Extreme"]
        intermediate = ds[ds["af_class"] == "Intermediate"]

        if len(extreme) >= 5 and len(intermediate) >= 5:
            stat, pval = scipy_stats.mannwhitneyu(
                intermediate["HPFineNGroups"], extreme["HPFineNGroups"], alternative="greater")
            n1, n2 = len(intermediate), len(extreme)
            r = 1 - 2 * stat / (n1 * n2)
            records.append({
                "sample": SAMPLE_SHORT[sample],
                "n_extreme": len(extreme), "n_intermediate": len(intermediate),
                "mean_ng_extreme": extreme["HPFineNGroups"].mean(),
                "mean_ng_intermediate": intermediate["HPFineNGroups"].mean(),
                "delta_ng": intermediate["HPFineNGroups"].mean() - extreme["HPFineNGroups"].mean(),
                "mw_p": pval, "rank_biserial_r": r,
                "direction": "Inter > Ext" if intermediate["HPFineNGroups"].mean() > extreme["HPFineNGroups"].mean() else "Ext > Inter",
                "mean_nr_extreme": extreme["NumReads"].mean(),
                "mean_nr_intermediate": intermediate["NumReads"].mean(),
            })

    cons_df = pd.DataFrame(records)
    cons_df.to_csv(DATADIR / "step2_paired_per_sample_consistency.tsv", sep="\t", index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    if len(cons_df) > 0:
        x = np.arange(len(cons_df))
        colors = ["#4CAF50" if d > 0 else "#F44336" for d in cons_df["delta_ng"]]
        ax.bar(x, cons_df["delta_ng"], color=colors, alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(cons_df["sample"], fontsize=9)
        ax.axhline(0, color="gray", ls="-", alpha=0.5)
        for i, row in cons_df.iterrows():
            ax.text(x[i], row["delta_ng"] + 0.02,
                    f"r={row['rank_biserial_r']:.2f}\np={row['mw_p']:.1e}",
                    ha="center", fontsize=7)

    ax.set_xlabel("Sample")
    ax.set_ylabel("Δ NGroups (Intermediate - Extreme)")
    ax.set_title("Paired Mode: Per-Sample NGroups Delta (CN1 LOH TP)\n"
                 f"Direction consistent: {sum(1 for d in cons_df['delta_ng'] if d > 0)}/{len(cons_df)}",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = FIGDIR / "p10_paired_ngroups_per_sample_consistency.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    return cons_df


def main():
    print("=" * 70)
    print("Step 2: Paired Mode — Intermediate AF × NGroups Cross-Analysis")
    print("=" * 70)

    _, df_paired = load_paired_data()

    print("\nLoading LOH.bed files...")
    loh_beds = load_loh_beds()
    print(f"\nAnnotating paired variants with LOH.bed...")
    df_paired["loh_bed_hit"] = annotate_loh_bed(df_paired, loh_beds)
    n_hit = df_paired["loh_bed_hit"].sum()
    print(f"  LOH.bed hit: {n_hit:,} ({100*n_hit/len(df_paired):.1f}%)")

    # Run analyses
    ng_stats = fig06_ngroups_by_af_class_cn_tier(df_paired)
    print("\n  NGroups stats (CN1):")
    cn1 = ng_stats[ng_stats["cn_tier"] == "CN1"]
    for _, row in cn1.iterrows():
        print(f"    {row['af_class']}: mean={row['mean_ngroups']:.3f}, n={row['n']:,}")

    test_df = fig07_ngroups_per_sample(df_paired)
    print("\n  Per-sample MW test:")
    for _, row in test_df.iterrows():
        sig = "***" if row["mw_p"] < 0.001 else ("**" if row["mw_p"] < 0.01 else ("*" if row["mw_p"] < 0.05 else "ns"))
        print(f"    {row['sample']}: ΔNG={row['mean_ng_intermediate']-row['mean_ng_extreme']:.3f}, "
              f"r={row['rank_biserial_r']:.3f}, p={row['mw_p']:.2e} {sig}")

    nr_df = fig08_numreads_controlled(df_paired)
    print("\n  NR-controlled:")
    for _, row in nr_df.iterrows():
        print(f"    NR {row['nr_bin']}: r={row['rank_biserial_r']:.3f}, p={row['mw_p']:.2e}")

    feat_df = fig09_methylation_features(df_paired)
    cons_df = fig10_per_sample_consistency(df_paired)

    n_positive = sum(1 for d in cons_df["delta_ng"] if d > 0)
    n_sig = sum(1 for p in cons_df["mw_p"] if p < 0.05)
    print(f"\n  H1p result: {n_positive}/{len(cons_df)} positive direction, {n_sig}/{len(cons_df)} significant")

    print("\n" + "=" * 70)
    print("Step 2 complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
