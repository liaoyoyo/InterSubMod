#!/usr/bin/env python3
"""
Step 2: Intermediate AF × HPFineNGroups × Methylation Cross-Analysis
=====================================================================
研究計劃：LOH 區域 Intermediate AF 與亞克隆偵測
目的：在 LOH 區域中，intermediate AF variants 是否有更高的
      HPFineNGroups 和不同的甲基化模式？

核心假說：
  H1: Intermediate AF LOH variants → 更多 NGroups (subclone diversity)
  H4: HPFineNGroups 反映亞克隆複雜度

輸入：all_region_rows.tsv.gz (master dataset)
輸出：
  - 圖表：NGroups 分布 × AF class × CN tier
  - 圖表：Methylation features × AF class
  - 統計檢驗結果
  - NumReads-controlled analysis
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
FIGDIR.mkdir(parents=True, exist_ok=True)
DATADIR.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_SHORT = {
    "HCC1395": "HCC1395", "HCC1395_DORADO": "HCC1395D",
    "COLO829": "COLO829", "H1437": "H1437", "H2009": "H2009",
    "HCC1937": "HCC1937", "HCC1954": "HCC1954",
}


def classify_af(af):
    if af < 0.1 or af > 0.9:
        return "Extreme"
    elif 0.4 <= af <= 0.6:
        return "Near-half"
    else:
        return "Intermediate"


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


def load_data():
    """Load master dataset, filter TO mode."""
    print("Loading master dataset...")
    cols = [
        "sample", "mode", "truth_label", "caller_af",
        "to_loh_bed_hit", "HPFineNGroups", "HPFineF", "HPFineP",
        "Coverage_Multiple", "NumReads", "HP_Ratio",
        "AlleleDelta", "Potential_LOH", "LOH_Subtype",
        "Quality_Score", "PairwiseMeanDist", "PairwiseMedianDist",
        "CramersV", "ClusterPermanovaF",
        "HP1FamilyN", "HP2FamilyN",
    ]
    df = pd.read_csv(MASTER, sep="\t", usecols=cols, low_memory=False)
    df_to = df[df["mode"] == "to"].copy()
    df_to["loh"] = df_to["to_loh_bed_hit"].astype(str).str.lower() == "true"
    df_to["af_class"] = df_to["caller_af"].apply(classify_af)
    df_to["cn_tier"] = df_to["Coverage_Multiple"].apply(cn_tier)

    # Derived features
    df_to["hp_imbalance"] = (df_to["HP_Ratio"] - 0.5).abs()
    total_hp = df_to["HP1FamilyN"].fillna(0) + df_to["HP2FamilyN"].fillna(0)
    df_to["effective_hp_reads"] = total_hp

    print(f"  TO mode rows: {len(df_to):,}")
    return df_to


def fig6_ngroups_by_af_class_cn_tier(df):
    """Figure 6: HPFineNGroups distribution by AF class in each CN tier (LOH only, TP only)."""
    print("\n=== Figure 6: NGroups by AF class × CN tier ===")

    df_loh_tp = df[(df["loh"] == True) & (df["truth_label"] == "TP")].copy()
    cn_order = ["CN1", "CN2", "CN3", "CN4+"]
    af_classes = ["Extreme", "Near-half", "Intermediate"]
    colors = {"Extreme": "#1976D2", "Near-half": "#4CAF50", "Intermediate": "#FF9800"}

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    records = []
    for ax_idx, cn in enumerate(cn_order):
        ax = axes.flatten()[ax_idx]
        sub = df_loh_tp[df_loh_tp["cn_tier"] == cn]

        ngroups_range = [0, 1, 2, 3, 4]
        x = np.arange(len(ngroups_range))
        width = 0.25

        for j, afc in enumerate(af_classes):
            sub_af = sub[sub["af_class"] == afc]
            if len(sub_af) == 0:
                continue

            # NGroups distribution
            counts = sub_af["HPFineNGroups"].value_counts().reindex(ngroups_range, fill_value=0)
            pcts = 100 * counts / counts.sum()

            ax.bar(x + (j - 1) * width, pcts.values, width, label=f"{afc} (n={len(sub_af):,})",
                   color=colors[afc], alpha=0.8)

            # Statistics
            records.append({
                "cn_tier": cn,
                "af_class": afc,
                "n": len(sub_af),
                "mean_ngroups": sub_af["HPFineNGroups"].mean(),
                "median_ngroups": sub_af["HPFineNGroups"].median(),
                "pct_ngroups_ge2": 100 * (sub_af["HPFineNGroups"] >= 2).sum() / len(sub_af),
                "mean_numreads": sub_af["NumReads"].mean(),
            })

        ax.set_xlabel("HPFineNGroups")
        ax.set_ylabel("% of variants")
        ax.set_title(f"{cn} (LOH TP)", fontsize=12, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(ngroups_range)
        ax.legend(fontsize=9)

    fig.suptitle("Step 2: HPFineNGroups Distribution by AF Class in LOH TP\n"
                 "Hypothesis: Intermediate AF → more NGroups (subclone diversity)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "06_ngroups_by_af_class_cn_tier.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    ngroups_stats = pd.DataFrame(records)
    ngroups_stats.to_csv(DATADIR / "step2_ngroups_by_af_class.tsv", sep="\t", index=False)
    print(f"  Stats saved: {DATADIR / 'step2_ngroups_by_af_class.tsv'}")

    return ngroups_stats


def fig7_ngroups_per_sample(df):
    """Figure 7: Per-sample NGroups mean in deletion LOH (CN≈1), Extreme vs Intermediate."""
    print("\n=== Figure 7: Per-sample NGroups in deletion LOH ===")

    df_cn1_tp = df[(df["loh"] == True) & (df["truth_label"] == "TP") & (df["cn_tier"] == "CN1")].copy()

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Left: mean NGroups per sample
    ax = axes[0]
    records = []
    for sample in SAMPLE_ORDER:
        ds = df_cn1_tp[df_cn1_tp["sample"] == sample]
        for afc in ["Extreme", "Intermediate"]:
            sub = ds[ds["af_class"] == afc]
            if len(sub) < 5:
                continue
            records.append({
                "sample": SAMPLE_SHORT[sample],
                "af_class": afc,
                "n": len(sub),
                "mean_ngroups": sub["HPFineNGroups"].mean(),
                "sem_ngroups": sub["HPFineNGroups"].sem(),
                "mean_numreads": sub["NumReads"].mean(),
            })

    rec_df = pd.DataFrame(records)
    samples_list = rec_df["sample"].unique()
    x = np.arange(len(samples_list))
    width = 0.35

    ext = rec_df[rec_df["af_class"] == "Extreme"]
    inter = rec_df[rec_df["af_class"] == "Intermediate"]

    ax.bar(x - width/2, ext["mean_ngroups"].values, width, yerr=ext["sem_ngroups"].values,
           label="Extreme AF", color="#1976D2", alpha=0.8, capsize=3)
    ax.bar(x + width/2, inter["mean_ngroups"].values, width, yerr=inter["sem_ngroups"].values,
           label="Intermediate AF", color="#FF9800", alpha=0.8, capsize=3)

    ax.set_xlabel("Sample")
    ax.set_ylabel("Mean HPFineNGroups")
    ax.set_title("Deletion LOH (CN≈1) TP: Mean NGroups", fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(samples_list, fontsize=9)
    ax.legend()

    # Right: Mann-Whitney test per sample
    ax2 = axes[1]
    test_records = []
    for sample in SAMPLE_ORDER:
        ds = df_cn1_tp[df_cn1_tp["sample"] == sample]
        extreme = ds[ds["af_class"] == "Extreme"]["HPFineNGroups"].dropna()
        intermediate = ds[ds["af_class"] == "Intermediate"]["HPFineNGroups"].dropna()

        if len(extreme) >= 5 and len(intermediate) >= 5:
            stat, pval = scipy_stats.mannwhitneyu(intermediate, extreme, alternative="greater")
            # Effect size: rank-biserial correlation
            n1, n2 = len(intermediate), len(extreme)
            r = 1 - 2 * stat / (n1 * n2)

            # Also test with NumReads-matched subsample
            # Match on NumReads quartile
            extreme_nr_q = extreme.index.map(lambda idx: ds.loc[idx, "NumReads"])
            intermediate_nr_q = intermediate.index.map(lambda idx: ds.loc[idx, "NumReads"])

            test_records.append({
                "sample": SAMPLE_SHORT[sample],
                "n_extreme": len(extreme),
                "n_intermediate": len(intermediate),
                "mean_ng_extreme": extreme.mean(),
                "mean_ng_intermediate": intermediate.mean(),
                "mw_U": stat,
                "mw_p": pval,
                "rank_biserial_r": r,
                "mean_nr_extreme": ds[ds["af_class"] == "Extreme"]["NumReads"].mean(),
                "mean_nr_intermediate": ds[ds["af_class"] == "Intermediate"]["NumReads"].mean(),
            })

    test_df = pd.DataFrame(test_records)

    # Bar chart of p-values (log scale)
    colors = ["#4CAF50" if p < 0.05 else "#F44336" for p in test_df["mw_p"]]
    ax2.barh(range(len(test_df)), -np.log10(test_df["mw_p"].clip(1e-300)), color=colors, alpha=0.8)
    ax2.axvline(-np.log10(0.05), color="red", ls="--", alpha=0.7, label="p=0.05")
    ax2.set_yticks(range(len(test_df)))
    ax2.set_yticklabels(test_df["sample"])
    ax2.set_xlabel("-log10(p-value)")
    ax2.set_title("Mann-Whitney: Intermediate > Extreme NGroups?", fontsize=12, fontweight="bold")
    ax2.legend()

    # Add effect sizes as text
    for idx, row in test_df.iterrows():
        ax2.text(-np.log10(row["mw_p"]) + 0.3, idx,
                 f"r={row['rank_biserial_r']:.3f}", fontsize=8, va="center")

    fig.suptitle("Step 2: HPFineNGroups in Deletion LOH (CN≈1) — Extreme vs Intermediate AF",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGDIR / "07_ngroups_per_sample_cn1.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    test_df.to_csv(DATADIR / "step2_ngroups_mw_test_cn1.tsv", sep="\t", index=False)
    print(f"  Test results saved: {DATADIR / 'step2_ngroups_mw_test_cn1.tsv'}")

    return test_df


def fig8_numreads_controlled(df):
    """Figure 8: NGroups comparison after controlling for NumReads (the known confound)."""
    print("\n=== Figure 8: NumReads-controlled NGroups analysis ===")

    df_cn1_tp = df[(df["loh"] == True) & (df["truth_label"] == "TP") & (df["cn_tier"] == "CN1")].copy()
    df_cn1_tp = df_cn1_tp.dropna(subset=["HPFineNGroups", "NumReads"])

    # Bin by NumReads
    nr_bins = [(10, 30), (30, 50), (50, 80), (80, 150), (150, 500)]
    nr_labels = ["10-30", "30-50", "50-80", "80-150", "150+"]

    fig, axes = plt.subplots(1, len(nr_bins), figsize=(22, 6), sharey=True)

    all_records = []
    for ax_idx, (nr_lo, nr_hi) in enumerate(nr_bins):
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

        # Statistical test
        if len(extreme) >= 5 and len(intermediate) >= 5:
            stat, pval = scipy_stats.mannwhitneyu(
                intermediate["HPFineNGroups"], extreme["HPFineNGroups"], alternative="greater"
            )
            n1, n2 = len(intermediate), len(extreme)
            r = 1 - 2 * stat / (n1 * n2)
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
            ax.set_title(f"NR={nr_labels[ax_idx]}\np={pval:.2e} {sig}\nr={r:.3f}",
                         fontsize=10, fontweight="bold")
            all_records.append({
                "nr_bin": nr_labels[ax_idx],
                "n_extreme": len(extreme),
                "n_intermediate": len(intermediate),
                "mean_ng_extreme": extreme["HPFineNGroups"].mean(),
                "mean_ng_intermediate": intermediate["HPFineNGroups"].mean(),
                "mw_p": pval,
                "rank_biserial_r": r,
            })
        else:
            ax.set_title(f"NR={nr_labels[ax_idx]}\n(insufficient data)", fontsize=10)

        ax.set_xlabel("NGroups")
        ax.set_xticks(x)
        ax.set_xticklabels(ngroups_range)
        ax.legend(fontsize=8)
        if ax_idx == 0:
            ax.set_ylabel("% of variants")

    fig.suptitle("Step 2: NGroups in Deletion LOH (CN≈1) — Controlled by NumReads\n"
                 "Test: Do intermediate AF variants have more NGroups after controlling NR confound?",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIGDIR / "08_ngroups_numreads_controlled.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    ctrl_df = pd.DataFrame(all_records)
    ctrl_df.to_csv(DATADIR / "step2_ngroups_numreads_controlled.tsv", sep="\t", index=False)
    print(f"  Controlled stats saved")

    return ctrl_df


def fig9_methylation_features(df):
    """Figure 9: Methylation feature comparison between AF classes in LOH."""
    print("\n=== Figure 9: Methylation features × AF class ===")

    df_loh_tp = df[(df["loh"] == True) & (df["truth_label"] == "TP")].copy()

    # Features to compare
    features = [
        ("HPFineF", "HP Fine F-statistic"),
        ("AlleleDelta", "Allele Delta"),
        ("PairwiseMeanDist", "Pairwise Mean Distance"),
        ("hp_imbalance", "HP Imbalance"),
        ("CramersV", "Cramér's V"),
        ("Quality_Score", "Quality Score"),
    ]

    cn_tiers = ["CN1", "CN2"]
    af_classes = ["Extreme", "Intermediate"]

    fig, axes = plt.subplots(len(features), len(cn_tiers), figsize=(14, 4 * len(features)))

    records = []
    for i, (feat, feat_label) in enumerate(features):
        for j, cn in enumerate(cn_tiers):
            ax = axes[i, j]
            sub = df_loh_tp[df_loh_tp["cn_tier"] == cn]

            for afc, color in zip(af_classes, ["#1976D2", "#FF9800"]):
                vals = sub[sub["af_class"] == afc][feat].dropna()
                if len(vals) > 0:
                    ax.hist(vals, bins=50, alpha=0.5, density=True, label=f"{afc} (n={len(vals):,})",
                            color=color, edgecolor="none")

            # Stats test
            extreme_vals = sub[sub["af_class"] == "Extreme"][feat].dropna()
            inter_vals = sub[sub["af_class"] == "Intermediate"][feat].dropna()

            if len(extreme_vals) >= 10 and len(inter_vals) >= 10:
                stat, pval = scipy_stats.mannwhitneyu(inter_vals, extreme_vals, alternative="two-sided")
                d = (inter_vals.mean() - extreme_vals.mean()) / np.sqrt(
                    (inter_vals.var() + extreme_vals.var()) / 2 + 1e-10)
                sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

                records.append({
                    "feature": feat,
                    "cn_tier": cn,
                    "n_extreme": len(extreme_vals),
                    "n_intermediate": len(inter_vals),
                    "mean_extreme": extreme_vals.mean(),
                    "mean_intermediate": inter_vals.mean(),
                    "cohen_d": d,
                    "mw_p": pval,
                })

                ax.text(0.98, 0.95, f"d={d:.3f} {sig}", transform=ax.transAxes,
                        ha="right", va="top", fontsize=9,
                        bbox=dict(boxstyle="round", fc="white", alpha=0.8))

            ax.set_title(f"{feat_label} — {cn}", fontsize=10)
            ax.legend(fontsize=8)
            if j == 0:
                ax.set_ylabel("Density")

    fig.suptitle("Step 2: Methylation Features — Extreme vs Intermediate AF in LOH TP",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = FIGDIR / "09_methylation_features_by_af_class.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    feat_df = pd.DataFrame(records)
    feat_df.to_csv(DATADIR / "step2_methylation_features_comparison.tsv", sep="\t", index=False)
    print(f"  Feature comparison saved")

    return feat_df


def fig10_per_sample_consistency(df):
    """Figure 10: Cross-sample consistency of NGroups difference (CN1 deletion LOH)."""
    print("\n=== Figure 10: Cross-sample consistency ===")

    df_cn1_tp = df[(df["loh"] == True) & (df["truth_label"] == "TP") & (df["cn_tier"] == "CN1")].copy()

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    ngroups_range = [0, 1, 2, 3, 4]
    consistency_records = []

    for i, sample in enumerate(SAMPLE_ORDER):
        ax = axes[i]
        ds = df_cn1_tp[df_cn1_tp["sample"] == sample]

        extreme = ds[ds["af_class"] == "Extreme"]
        intermediate = ds[ds["af_class"] == "Intermediate"]

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

        # Test
        if len(extreme) >= 5 and len(intermediate) >= 5:
            stat, pval = scipy_stats.mannwhitneyu(
                intermediate["HPFineNGroups"], extreme["HPFineNGroups"], alternative="greater"
            )
            n1, n2 = len(intermediate), len(extreme)
            r = 1 - 2 * stat / (n1 * n2)
            direction = "Inter > Ext" if intermediate["HPFineNGroups"].mean() > extreme["HPFineNGroups"].mean() else "Ext ≥ Inter"
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

            consistency_records.append({
                "sample": sample,
                "n_extreme": len(extreme),
                "n_intermediate": len(intermediate),
                "mean_ng_extreme": extreme["HPFineNGroups"].mean(),
                "mean_ng_intermediate": intermediate["HPFineNGroups"].mean(),
                "delta_ng": intermediate["HPFineNGroups"].mean() - extreme["HPFineNGroups"].mean(),
                "mw_p": pval,
                "rank_biserial_r": r,
                "direction": direction,
                "mean_nr_extreme": extreme["NumReads"].mean(),
                "mean_nr_intermediate": intermediate["NumReads"].mean(),
            })

            ax.set_title(f"{SAMPLE_SHORT[sample]}\n{direction} {sig} (r={r:.3f})",
                         fontsize=10, fontweight="bold")
        else:
            ax.set_title(f"{SAMPLE_SHORT[sample]}\n(insufficient data)", fontsize=10)

        ax.set_xlabel("NGroups")
        ax.set_ylabel("%")
        ax.set_xticks(x)
        ax.set_xticklabels(ngroups_range)
        ax.legend(fontsize=8)

    axes[-1].axis("off")

    fig.suptitle("Step 2: NGroups in Deletion LOH (CN≈1) — Per-Sample Consistency\n"
                 "H1: Intermediate AF variants have more methylation groups",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGDIR / "10_ngroups_per_sample_consistency.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    cons_df = pd.DataFrame(consistency_records)
    cons_df.to_csv(DATADIR / "step2_per_sample_consistency.tsv", sep="\t", index=False)
    print(f"  Consistency results saved")

    return cons_df


def print_step2_summary(ngroups_stats, test_df, ctrl_df, feat_df, cons_df):
    """Print comprehensive Step 2 summary."""
    print("\n" + "=" * 80)
    print("STEP 2 SUMMARY: Intermediate AF × HPFineNGroups × Methylation")
    print("=" * 80)

    print("\n--- NGroups by AF class (all samples, LOH TP) ---")
    for _, row in ngroups_stats.iterrows():
        print(f"  {row['cn_tier']} × {row['af_class']}: n={row['n']:,}, "
              f"mean NGroups={row['mean_ngroups']:.3f}, "
              f"NGroups≥2: {row['pct_ngroups_ge2']:.1f}%, "
              f"mean NR={row['mean_numreads']:.0f}")

    if test_df is not None and len(test_df) > 0:
        print("\n--- Per-sample Mann-Whitney: Intermediate > Extreme NGroups (CN≈1) ---")
        n_positive = 0
        for _, row in test_df.iterrows():
            sig = "***" if row["mw_p"] < 0.001 else "**" if row["mw_p"] < 0.01 else "*" if row["mw_p"] < 0.05 else "ns"
            direction = "✓" if row["mean_ng_intermediate"] > row["mean_ng_extreme"] else "✗"
            print(f"  {row['sample']}: Intermediate={row['mean_ng_intermediate']:.3f} vs "
                  f"Extreme={row['mean_ng_extreme']:.3f}, p={row['mw_p']:.2e} {sig}, "
                  f"r={row['rank_biserial_r']:.3f} {direction}"
                  f"  [NR: Inter={row['mean_nr_intermediate']:.0f}, Ext={row['mean_nr_extreme']:.0f}]")
            if row["mean_ng_intermediate"] > row["mean_ng_extreme"]:
                n_positive += 1
        print(f"  Direction consistency: {n_positive}/{len(test_df)} samples positive")

    if ctrl_df is not None and len(ctrl_df) > 0:
        print("\n--- NumReads-controlled NGroups test (CN≈1, all samples) ---")
        for _, row in ctrl_df.iterrows():
            sig = "***" if row["mw_p"] < 0.001 else "**" if row["mw_p"] < 0.01 else "*" if row["mw_p"] < 0.05 else "ns"
            print(f"  NR={row['nr_bin']}: Intermediate={row['mean_ng_intermediate']:.3f} vs "
                  f"Extreme={row['mean_ng_extreme']:.3f}, p={row['mw_p']:.2e} {sig}, r={row['rank_biserial_r']:.3f}")

    if feat_df is not None and len(feat_df) > 0:
        print("\n--- Methylation features: Intermediate vs Extreme (top effects) ---")
        feat_df_sorted = feat_df.reindex(feat_df["cohen_d"].abs().sort_values(ascending=False).index)
        for _, row in feat_df_sorted.head(10).iterrows():
            sig = "***" if row["mw_p"] < 0.001 else "**" if row["mw_p"] < 0.01 else "*" if row["mw_p"] < 0.05 else "ns"
            print(f"  {row['feature']}×{row['cn_tier']}: d={row['cohen_d']:.3f} {sig}, "
                  f"Inter={row['mean_intermediate']:.4f} vs Ext={row['mean_extreme']:.4f}")

    if cons_df is not None and len(cons_df) > 0:
        print("\n--- Cross-sample consistency (CN≈1 deletion LOH) ---")
        n_pos = (cons_df["delta_ng"] > 0).sum()
        n_sig = (cons_df["mw_p"] < 0.05).sum()
        print(f"  Direction: {n_pos}/{len(cons_df)} samples have Intermediate > Extreme NGroups")
        print(f"  Significant (p<0.05): {n_sig}/{len(cons_df)} samples")
        print(f"  NumReads confound check:")
        for _, row in cons_df.iterrows():
            nr_ratio = row["mean_nr_intermediate"] / row["mean_nr_extreme"] if row["mean_nr_extreme"] > 0 else float("inf")
            print(f"    {row['sample']}: NR ratio (Inter/Ext) = {nr_ratio:.2f}")


def main():
    df = load_data()

    ngroups_stats = fig6_ngroups_by_af_class_cn_tier(df)
    test_df = fig7_ngroups_per_sample(df)
    ctrl_df = fig8_numreads_controlled(df)
    feat_df = fig9_methylation_features(df)
    cons_df = fig10_per_sample_consistency(df)

    print_step2_summary(ngroups_stats, test_df, ctrl_df, feat_df, cons_df)


if __name__ == "__main__":
    main()
