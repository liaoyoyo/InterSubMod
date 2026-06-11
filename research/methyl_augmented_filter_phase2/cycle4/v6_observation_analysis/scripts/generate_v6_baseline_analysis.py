#!/usr/bin/env python3
"""V6 vs baseline 觀察驅動分析 — for cycle 5+ feature engineering hypothesis.

User 2026-05-20 要求 (after LOSO confirms HCC1395 AF overfit):
  "先清楚合適的觀察數據與機制才能找到更有價值的特徵與用算法分析出更有效的特徵"

非 fit-and-test 框架。純觀察 5 樣本 V6 features distribution × TP/FP +
HCC1395 V6 vs baseline NG distribution shift + chr-level + LOH interaction
→ 為 cycle 5+ 提供 hypothesis input。

Outputs (8 figures + 3 tables in cycle4/v6_observation_analysis/):
  figures/01_v6_vs_baseline_NG_dist.png   HCC1395 V6 vs baseline NG dist shift
  figures/02_auc_heatmap.png              5 samples × 10 features AUC
  figures/03_cohen_d_heatmap.png          5 samples × 10 features Cohen's d
  figures/04_non_methyl_dist.png          4 non-methyl features × 5 samples TP/FP KDE
  figures/05_methyl_dist.png              5 methyl features × 5 samples TP/FP KDE
  figures/06_chr8_hotspot.png             HCC1395 V6 marker rate per-chr
  figures/07_loh_interaction.png          LOH=Inner vs Outer/NA AUC delta per feature
  figures/08_5goal_map.png                Feature × 5-goal applicability matrix

  data/feature_auc_cohen.tsv              全 sample × feature AUC + Cohen's d
  data/hcc1395_v6_baseline_ng_shift.tsv   NG distribution shift
  data/loh_stratified_auc.tsv             LOH-stratified AUC per feature
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from sklearn.metrics import roc_auc_score

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CYCLE1_MASTER = REPO / "research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv"
FOUR_WAY = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv"
CYCLE2_DATA = REPO / "research/methyl_augmented_filter_phase2/cycle2/data"

OUT = REPO / "research/methyl_augmented_filter_phase2/cycle4/v6_observation_analysis"
FIG = OUT / "figures"
DATA = OUT / "data"
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]
SAMPLE_COLORS = {"HCC1395": "#D97757", "HCC1937": "#15803D", "HCC1954": "#A16207",
                 "H1437": "#1E3A8A", "H2009": "#6B7280"}

NON_METHYL = {
    "caller_af": "Caller AF (BAM-invariant)",
    "loh_inner_flag": "LOH Inner flag (BAM-invariant)",
    "Coverage_Multiple_imp": "Coverage Multiple (BAM-invariant)",
    "V6_off_NG": "V6 NG count (V6-specific)",
}
METHYL = {
    "V6_off_meth_HPMergedDelta": "HP Merged β Delta",
    "V6_off_meth_HPFineF": "HP Fine F-stat",
    "V6_off_meth_NME_imbalance": "NME imbalance",
    "V6_off_meth_Epipoly_Delta": "Epipoly Delta",
    "V6_off_meth_ClusterPermanovaF": "Cluster PERMANOVA F",
}
ALL_FEATURES = {**NON_METHYL, **METHYL}


def load_sample(sample: str) -> pd.DataFrame:
    if sample == "HCC1395":
        df = pd.read_csv(CYCLE1_MASTER, sep="\t", low_memory=False)
    else:
        df = pd.read_csv(CYCLE2_DATA / f"{sample}_master_augmented.tsv", sep="\t", low_memory=False)
    if "y" not in df.columns and "label" in df.columns:
        df["y"] = (df["label"] == "TP").astype(int)
    if "loh_inner_flag" not in df.columns and "loh_side" in df.columns:
        df["loh_inner_flag"] = (df["loh_side"] == "Inner").astype(int)
    if "Coverage_Multiple_imp" not in df.columns and "Coverage_Multiple" in df.columns:
        df["Coverage_Multiple_imp"] = df["Coverage_Multiple"].fillna(df["Coverage_Multiple"].median())
    return df


# ─────────── Figure 1: HCC1395 V6 vs baseline NG distribution shift ─────────
def fig_01_ng_shift():
    df = pd.read_csv(FOUR_WAY, sep="\t", low_memory=False)
    df["y"] = (df["label"] == "TP").astype(int)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, label_val, label_name in zip(axes, [1, 0], ["TP", "FP"]):
        sub = df[df["y"] == label_val]
        ng_bins = [0, 1, 2, 3, 4, 5, 100]
        ng_labels = ["0", "1", "2", "3", "4", "5+"]
        baseline_ng = pd.cut(sub["baseline_off_NG"], bins=ng_bins, labels=ng_labels,
                              right=False, include_lowest=True).value_counts(normalize=True).sort_index()
        v6_ng = pd.cut(sub["V6_off_NG"], bins=ng_bins, labels=ng_labels,
                        right=False, include_lowest=True).value_counts(normalize=True).sort_index()

        x = np.arange(len(ng_labels))
        ax.bar(x - 0.2, baseline_ng.reindex(ng_labels).fillna(0).values, 0.4,
               label="Baseline", color="#6B7280", edgecolor="black", linewidth=0.4)
        ax.bar(x + 0.2, v6_ng.reindex(ng_labels).fillna(0).values, 0.4,
               label="V6", color="#D97757", edgecolor="black", linewidth=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(ng_labels)
        ax.set_xlabel("NG group count")
        ax.set_ylabel(f"Fraction of {label_name} regions")
        ax.set_title(f"HCC1395 {label_name}: V6 vs baseline NG distribution\n"
                     f"n_{label_name} = {len(sub):,}", fontsize=11, fontweight="bold")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Figure 1 — V6 vs baseline NG distribution shift (HCC1395)\n"
                 "V6 顯著增加 NG≥3 比例（marker候選增加 +52.4% 全局）",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "01_v6_vs_baseline_NG_dist.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ 01_v6_vs_baseline_NG_dist.png")


# ─────────── Figure 2+3: AUC + Cohen's d heatmap ────────────────────────────
def compute_auc_cohen_per_sample(sample_dfs: dict) -> pd.DataFrame:
    rows = []
    for sample, df in sample_dfs.items():
        for feat, _ in ALL_FEATURES.items():
            if feat not in df.columns:
                rows.append({"sample": sample, "feature": feat,
                              "n_TP": int((df["y"]==1).sum()),
                              "n_FP": int((df["y"]==0).sum()),
                              "auc": np.nan, "cohen_d": np.nan, "ks_d": np.nan})
                continue
            sub = df[[feat, "y"]].dropna()
            if len(sub) < 10 or sub["y"].nunique() < 2:
                rows.append({"sample": sample, "feature": feat,
                              "n_TP": int((df["y"]==1).sum()),
                              "n_FP": int((df["y"]==0).sum()),
                              "auc": np.nan, "cohen_d": np.nan, "ks_d": np.nan})
                continue
            try:
                auc = roc_auc_score(sub["y"], sub[feat])
            except Exception:
                auc = np.nan
            tp_vals = sub.loc[sub["y"]==1, feat].values
            fp_vals = sub.loc[sub["y"]==0, feat].values
            if len(tp_vals) > 1 and len(fp_vals) > 1:
                pooled_std = np.sqrt((np.var(tp_vals, ddof=1) + np.var(fp_vals, ddof=1)) / 2)
                cohen_d = (np.mean(tp_vals) - np.mean(fp_vals)) / pooled_std if pooled_std > 0 else 0
                ks_d, _ = ks_2samp(tp_vals, fp_vals)
            else:
                cohen_d = np.nan
                ks_d = np.nan
            rows.append({"sample": sample, "feature": feat,
                          "n_TP": int((df["y"]==1).sum()),
                          "n_FP": int((df["y"]==0).sum()),
                          "auc": auc, "cohen_d": cohen_d, "ks_d": ks_d})
    return pd.DataFrame(rows)


def fig_02_03_heatmaps(stats_df: pd.DataFrame):
    feat_order = list(ALL_FEATURES.keys())
    for metric, fname, vmin, vmax, cmap, title_suffix in [
        ("auc", "02_auc_heatmap.png", 0.4, 0.7, "RdYlGn", "AUC (0.5=random)"),
        ("cohen_d", "03_cohen_d_heatmap.png", -1.5, 1.5, "RdBu_r", "Cohen's d (effect size)"),
    ]:
        mat = np.zeros((len(SAMPLE_ORDER), len(feat_order)))
        for i, s in enumerate(SAMPLE_ORDER):
            for j, f in enumerate(feat_order):
                row = stats_df[(stats_df["sample"]==s) & (stats_df["feature"]==f)]
                mat[i, j] = float(row[metric].iloc[0]) if len(row) else np.nan

        fig, ax = plt.subplots(figsize=(11, 5))
        im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        for i in range(len(SAMPLE_ORDER)):
            for j in range(len(feat_order)):
                v = mat[i, j]
                if not np.isnan(v):
                    color = "white" if (metric == "auc" and (v < 0.45 or v > 0.65)) or \
                            (metric == "cohen_d" and abs(v) > 0.8) else "black"
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9,
                            color=color, fontweight="bold")
        ax.set_xticks(range(len(feat_order)))
        ax.set_xticklabels([ALL_FEATURES[f] for f in feat_order], rotation=30, ha="right", fontsize=9)
        ax.set_yticks(range(len(SAMPLE_ORDER)))
        ax.set_yticklabels(SAMPLE_ORDER)
        ax.set_title(f"{title_suffix} — TP vs FP separability per feature × sample",
                     fontsize=12, fontweight="bold")
        cbar = fig.colorbar(im, ax=ax, label=metric)
        fig.tight_layout()
        fig.savefig(FIG / fname, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ {fname}")


# ─────────── Figure 4+5: Distribution overlays ──────────────────────────────
def fig_04_05_distributions(sample_dfs: dict):
    for feats, fname, suptitle, ncols in [
        (NON_METHYL, "04_non_methyl_dist.png", "4 Non-Methyl Features × 5 Samples TP/FP Distribution", 2),
        (METHYL, "05_methyl_dist.png", "5 ISM Methylation Features × 5 Samples TP/FP Distribution", 3),
    ]:
        nrows = (len(feats) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
        if nrows == 1:
            axes = axes.reshape(1, -1)

        for idx, (feat, feat_name) in enumerate(feats.items()):
            ax = axes[idx // ncols, idx % ncols]
            for s in SAMPLE_ORDER:
                df = sample_dfs[s]
                if feat not in df.columns:
                    continue
                tp_vals = df.loc[df["y"]==1, feat].dropna().values
                fp_vals = df.loc[df["y"]==0, feat].dropna().values
                if len(tp_vals) > 10 and len(fp_vals) > 10:
                    # Determine bins
                    all_vals = np.concatenate([tp_vals, fp_vals])
                    bins = np.linspace(np.percentile(all_vals, 1),
                                        np.percentile(all_vals, 99), 30)
                    ax.hist(tp_vals, bins=bins, alpha=0.4, density=True,
                            color=SAMPLE_COLORS[s], label=f"{s} TP", histtype="step", linewidth=1.8)
                    ax.hist(fp_vals, bins=bins, alpha=0.6, density=True,
                            color=SAMPLE_COLORS[s], label=f"{s} FP", histtype="step",
                            linewidth=1.0, linestyle="--")
            ax.set_xlabel(feat_name, fontsize=10)
            ax.set_ylabel("density")
            ax.set_title(feat_name, fontsize=10, fontweight="bold")
            if idx == 0:
                ax.legend(fontsize=6.5, ncol=2, loc="upper right")
            ax.grid(alpha=0.3)

        # Hide empty axes
        for idx in range(len(feats), nrows * ncols):
            axes[idx // ncols, idx % ncols].set_visible(False)

        fig.suptitle(f"Figure — {suptitle}\n"
                     "實線 = TP, 虛線 = FP；5 samples color-coded",
                     fontsize=12, fontweight="bold", y=1.00)
        fig.tight_layout()
        fig.savefig(FIG / fname, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ {fname}")


# ─────────── Figure 6: chr-level marker rate ─────────────────────────────────
def fig_06_chr8_hotspot():
    df = pd.read_csv(FOUR_WAY, sep="\t", low_memory=False)
    df["y"] = (df["label"] == "TP").astype(int)
    # Marker = NG >= 3
    df["baseline_marker"] = (df["baseline_off_NG"] >= 3).astype(int)
    df["v6_marker"] = (df["V6_off_NG"] >= 3).astype(int)

    chr_order = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    rows = []
    for ch in chr_order:
        sub = df[df["chr"] == ch]
        if len(sub) == 0:
            continue
        n_tp = int((sub["y"]==1).sum())
        n_fp = int((sub["y"]==0).sum())
        # Marker count among TP / FP
        n_tp_marker_v6 = int(((sub["y"]==1) & (sub["v6_marker"]==1)).sum())
        n_fp_marker_v6 = int(((sub["y"]==0) & (sub["v6_marker"]==1)).sum())
        n_tp_marker_b = int(((sub["y"]==1) & (sub["baseline_marker"]==1)).sum())
        n_fp_marker_b = int(((sub["y"]==0) & (sub["baseline_marker"]==1)).sum())
        rows.append({
            "chr": ch, "n_TP": n_tp, "n_FP": n_fp,
            "v6_marker_rate_TP": n_tp_marker_v6/n_tp if n_tp else 0,
            "v6_marker_rate_FP": n_fp_marker_v6/n_fp if n_fp else 0,
            "b_marker_rate_TP": n_tp_marker_b/n_tp if n_tp else 0,
            "b_marker_rate_FP": n_fp_marker_b/n_fp if n_fp else 0,
        })
    chrdf = pd.DataFrame(rows)

    fig, axes = plt.subplots(2, 1, figsize=(14, 9))

    # Top: TP marker rate
    x = np.arange(len(chrdf))
    axes[0].bar(x - 0.2, chrdf["b_marker_rate_TP"], 0.4, label="Baseline", color="#6B7280",
                 edgecolor="black", linewidth=0.4)
    axes[0].bar(x + 0.2, chrdf["v6_marker_rate_TP"], 0.4, label="V6", color="#D97757",
                 edgecolor="black", linewidth=0.4)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(chrdf["chr"], rotation=45, fontsize=8)
    axes[0].set_ylabel("TP Marker rate (NG≥3)", fontsize=10)
    axes[0].set_title("HCC1395 per-chr TP marker rate (V6 vs Baseline)", fontsize=11, fontweight="bold")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.3)

    # Bottom: FP marker rate (idea: chr8 hotspot)
    axes[1].bar(x - 0.2, chrdf["b_marker_rate_FP"], 0.4, label="Baseline", color="#6B7280",
                 edgecolor="black", linewidth=0.4)
    axes[1].bar(x + 0.2, chrdf["v6_marker_rate_FP"], 0.4, label="V6", color="#C2410C",
                 edgecolor="black", linewidth=0.4)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(chrdf["chr"], rotation=45, fontsize=8)
    axes[1].set_ylabel("FP Marker rate (NG≥3)", fontsize=10)
    axes[1].set_title("HCC1395 per-chr FP marker rate (V6 vs Baseline) — chr8 hotspot 視覺化",
                       fontsize=11, fontweight="bold")
    axes[1].legend()
    axes[1].grid(axis="y", alpha=0.3)

    fig.suptitle("Figure 6 — chr-level V6 vs Baseline marker rate (HCC1395)\n"
                 "比較 V6 改進在哪些 chr; chr8 特殊熱點 in FP marker",
                 fontsize=12, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(FIG / "06_chr8_hotspot.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ 06_chr8_hotspot.png")


# ─────────── Figure 7: LOH × feature interaction ─────────────────────────────
def fig_07_loh_interaction(sample_dfs: dict):
    feat_order = list(ALL_FEATURES.keys())
    rows = []
    for sample, df in sample_dfs.items():
        for feat in feat_order:
            if feat not in df.columns or "loh_inner_flag" not in df.columns:
                continue
            inner = df[df["loh_inner_flag"] == 1]
            outer = df[df["loh_inner_flag"] == 0]
            for zone, sub in [("Inner", inner), ("Outer/NA", outer)]:
                sub2 = sub[[feat, "y"]].dropna()
                if len(sub2) < 10 or sub2["y"].nunique() < 2:
                    rows.append({"sample": sample, "feature": feat, "zone": zone, "auc": np.nan})
                    continue
                try:
                    auc = roc_auc_score(sub2["y"], sub2[feat])
                except Exception:
                    auc = np.nan
                rows.append({"sample": sample, "feature": feat, "zone": zone, "auc": auc})
    lohdf = pd.DataFrame(rows)
    lohdf.to_csv(DATA / "loh_stratified_auc.tsv", sep="\t", index=False)

    # Compute AUC delta (Inner - Outer/NA) per sample × feature
    fig, ax = plt.subplots(figsize=(11, 5))
    mat = np.zeros((len(SAMPLE_ORDER), len(feat_order)))
    for i, s in enumerate(SAMPLE_ORDER):
        for j, f in enumerate(feat_order):
            inner = lohdf[(lohdf["sample"]==s) & (lohdf["feature"]==f) & (lohdf["zone"]=="Inner")]
            outer = lohdf[(lohdf["sample"]==s) & (lohdf["feature"]==f) & (lohdf["zone"]=="Outer/NA")]
            if len(inner) and len(outer):
                a_in = inner["auc"].iloc[0]
                a_out = outer["auc"].iloc[0]
                if not (np.isnan(a_in) or np.isnan(a_out)):
                    mat[i, j] = a_in - a_out
                else:
                    mat[i, j] = np.nan
            else:
                mat[i, j] = np.nan

    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.15, vmax=0.15, aspect="auto")
    for i in range(len(SAMPLE_ORDER)):
        for j in range(len(feat_order)):
            v = mat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.3f}", ha="center", va="center", fontsize=9, fontweight="bold",
                        color="white" if abs(v) > 0.08 else "black")
    ax.set_xticks(range(len(feat_order)))
    ax.set_xticklabels([ALL_FEATURES[f] for f in feat_order], rotation=30, ha="right", fontsize=9)
    ax.set_yticks(range(len(SAMPLE_ORDER)))
    ax.set_yticklabels(SAMPLE_ORDER)
    ax.set_title("Figure 7 — AUC delta: Inner vs Outer/NA per feature × sample\n"
                 "正值 = feature 在 LOH Inner 更有區分力; 揭示 LOH × feature interaction",
                 fontsize=12, fontweight="bold")
    fig.colorbar(im, ax=ax, label="AUC(Inner) − AUC(Outer/NA)")
    fig.tight_layout()
    fig.savefig(FIG / "07_loh_interaction.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ 07_loh_interaction.png")


# ─────────── Figure 8: 5-goal map ────────────────────────────────────────────
def fig_08_5goal_map():
    feat_order = list(ALL_FEATURES.keys())
    goals = ["G1 Per-CpG×HP×ALT", "G2 HPFineN marker", "G3 Two-hit order",
              "G4 TO normal model", "G5 LR filter"]
    # Manually map applicability (1=high, 0.5=medium, 0=none)
    applicability = {
        "caller_af":                {"G1": 0.5, "G2": 0,   "G3": 0,   "G4": 1,   "G5": 1},
        "loh_inner_flag":           {"G1": 0,   "G2": 0.5, "G3": 1,   "G4": 0,   "G5": 1},
        "Coverage_Multiple_imp":    {"G1": 0,   "G2": 0,   "G3": 1,   "G4": 0.5, "G5": 1},
        "V6_off_NG":                {"G1": 0.5, "G2": 1,   "G3": 0.5, "G4": 0,   "G5": 1},
        "V6_off_meth_HPMergedDelta": {"G1": 1,   "G2": 0.5, "G3": 0.5, "G4": 0,   "G5": 0.5},
        "V6_off_meth_HPFineF":       {"G1": 1,   "G2": 0.5, "G3": 0.5, "G4": 0,   "G5": 0.5},
        "V6_off_meth_NME_imbalance": {"G1": 0.5, "G2": 0,   "G3": 0.5, "G4": 0,   "G5": 0.5},
        "V6_off_meth_Epipoly_Delta": {"G1": 0.5, "G2": 0,   "G3": 0.5, "G4": 0,   "G5": 0.5},
        "V6_off_meth_ClusterPermanovaF": {"G1": 1, "G2": 0.5, "G3": 0.5, "G4": 0,   "G5": 0.5},
    }
    mat = np.zeros((len(feat_order), len(goals)))
    for i, f in enumerate(feat_order):
        for j, g in enumerate(goals):
            g_key = g.split()[0]
            mat[i, j] = applicability[f].get(g_key, 0)

    fig, ax = plt.subplots(figsize=(9, 6))
    im = ax.imshow(mat, cmap="YlGn", vmin=0, vmax=1, aspect="auto")
    for i in range(len(feat_order)):
        for j in range(len(goals)):
            v = mat[i, j]
            label = "★★" if v == 1 else ("★" if v == 0.5 else "")
            if label:
                ax.text(j, i, label, ha="center", va="center", fontsize=13,
                        color="white" if v == 1 else "black", fontweight="bold")
    ax.set_xticks(range(len(goals)))
    ax.set_xticklabels(goals, rotation=15, fontsize=9)
    ax.set_yticks(range(len(feat_order)))
    ax.set_yticklabels([ALL_FEATURES[f] for f in feat_order], fontsize=9)
    ax.set_title("Figure 8 — Feature × 5 終極目標 Applicability Map\n"
                 "★★=高度相關 / ★=中度相關 / 空白=低/無；以 InterSubMod 5 目標為軸",
                 fontsize=11, fontweight="bold")
    fig.colorbar(im, ax=ax, label="Applicability")
    fig.tight_layout()
    fig.savefig(FIG / "08_5goal_map.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ 08_5goal_map.png")


if __name__ == "__main__":
    print("Loading 5 sample master TSVs ...")
    sample_dfs = {s: load_sample(s) for s in SAMPLE_ORDER}
    for s, df in sample_dfs.items():
        print(f"  {s}: rows={len(df):,} TP={int((df['y']==1).sum()):,} FP={int((df['y']==0).sum()):,}")

    print("\nGenerating 8 V6 observation figures ...")
    fig_01_ng_shift()

    print("Computing AUC + Cohen's d ...")
    stats_df = compute_auc_cohen_per_sample(sample_dfs)
    stats_df.to_csv(DATA / "feature_auc_cohen.tsv", sep="\t", index=False)
    fig_02_03_heatmaps(stats_df)

    fig_04_05_distributions(sample_dfs)
    fig_06_chr8_hotspot()
    fig_07_loh_interaction(sample_dfs)
    fig_08_5goal_map()

    print(f"\n✓ 8 figures saved to: {FIG}")
    print(f"✓ Data files saved to: {DATA}")
