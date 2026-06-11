#!/usr/bin/env python3
"""Cycle 4 Stage 1d — Gap matrix + 12 methyl features distribution boxplot.

Reads inventory TSVs + cycle 1+2 master TSVs.

Outputs:
    cycle4/figures/stage1_gap_matrix.png       3-panel: tried vs not-tried heatmap
    cycle4/figures/stage1_features_boxplot.png 12 features × 5 samples distribution
"""
from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
CYCLE4 = REPO / "research/methyl_augmented_filter_phase2/cycle4"
FIG = CYCLE4 / "figures"
FIG.mkdir(parents=True, exist_ok=True)

FEAT_INV = CYCLE4 / "stage1_mapping/feature_inventory.tsv"
ALGO_INV = CYCLE4 / "stage1_mapping/algorithm_inventory.tsv"
OBS_INV = CYCLE4 / "stage1_mapping/observation_inventory.tsv"

CYCLE1_MASTER = REPO / "research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv"
CYCLE2_DATA = REPO / "research/methyl_augmented_filter_phase2/cycle2/data"

SAMPLE_ORDER = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]

# All 12 methyl features in master TSV (V6_off prefix)
METHYL_FEATURES = [
    ("HPMergedDelta", "F_A1"),
    ("HPFineF", "F_A2"),
    ("NME_imbalance", "F_A3"),
    ("Epipoly_Delta", "F_A4"),
    ("ClusterPermanovaF", "F_A5"),
    ("HPMergedP", "F_B1"),
    ("HPFineSig", "F_B2"),
    ("HPFineNGroups", "F_B3"),
    ("AlleleP", "F_B4"),
    ("AlleleDelta", "F_B5"),
    ("Entropy_Imbalance", "F_B6"),
    ("Epipoly_HP1", "F_B7_HP1"),
]


def fig_gap_matrix():
    """3-panel: feature × used / algorithm × tried / observation × done."""
    feat = pd.read_csv(FEAT_INV, sep="\t")
    algo = pd.read_csv(ALGO_INV, sep="\t")
    obs = pd.read_csv(OBS_INV, sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(16, 6), gridspec_kw={"width_ratios": [1, 1.2, 1.3]})

    # Panel 1: Features
    ax = axes[0]
    feat["used_score"] = feat["used_in_cycle1_LR"].map({"YES": 2, "NO": 0}) + \
                         feat["used_in_cycle3_ablation"].map({"YES": 1, "NO": 0})
    feat = feat.sort_values("used_score", ascending=False)
    colors = ["#15803D" if s == 3 else ("#A16207" if s == 2 else "#6B7280")
              for s in feat["used_score"]]
    ax.barh(feat["feature_id"] + " " + feat["feature_name"].str.replace("V6_off_meth_", ""),
            feat["used_score"], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_xlim(0, 3.5)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["unused", "ablation\nonly", "cycle 1\nonly", "both"], fontsize=9)
    ax.set_title("12 Methyl Features — usage status", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    ax.tick_params(axis="y", labelsize=8.5)
    ax.grid(axis="x", alpha=0.3)

    # Panel 2: Algorithms
    ax = axes[1]
    algo["status_num"] = algo["status"].map({"tried": 1, "tried (sanity)": 1, "not tried": 0})
    colors = ["#15803D" if s == 1 else "#C2410C" for s in algo["status_num"]]
    sizes = pd.to_numeric(algo["cost_est_days"], errors="coerce").fillna(0.3)
    ax.barh(algo["algo_id"] + " " + algo["algorithm"],
            sizes,
            color=colors, edgecolor="black", linewidth=0.4)
    ax.set_xlabel("estimated cost (days)\ngreen = tried, red = not tried", fontsize=9)
    ax.set_title("Algorithms — tried vs not tried", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    ax.tick_params(axis="y", labelsize=8.5)
    ax.grid(axis="x", alpha=0.3)

    # Panel 3: Observations
    ax = axes[2]
    obs["status_num"] = obs["status"].map({"tried": 1, "not done": 0, "not done (reserve)": 0, "not done (low ROI)": 0})
    colors = []
    for _, r in obs.iterrows():
        if r["status_num"] == 1:
            colors.append("#15803D")
        elif "reserve" in str(r["status"]) or "low ROI" in str(r["status"]):
            colors.append("#A16207")
        else:
            colors.append("#C2410C")
    sizes = pd.to_numeric(obs["cost_est_days"], errors="coerce").fillna(0.3)
    ax.barh(obs["obs_id"] + " " + obs["observation"].str[:45],
            sizes, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_xlabel("estimated cost (days)\ngreen = tried, orange = reserve/low ROI, red = priority not done", fontsize=8.5)
    ax.set_title("Observations — done vs not done", fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    ax.tick_params(axis="y", labelsize=7.5)
    ax.grid(axis="x", alpha=0.3)

    fig.suptitle("Cycle 4 Stage 1 — Gap Matrix (feature × algorithm × observation)\n"
                 "Reveals tried vs not-tried 三維空間 for cycle 3 vestigial 結論 surrounding context",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "stage1_gap_matrix.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {FIG}/stage1_gap_matrix.png")


def fig_features_boxplot():
    """12 methyl features × 5 samples distribution boxplot."""
    # Load all 5 sample master TSVs
    samples_data = {}
    samples_data["HCC1395"] = pd.read_csv(CYCLE1_MASTER, sep="\t", low_memory=False)
    for s in ["H1437", "H2009", "HCC1954", "HCC1937"]:
        samples_data[s] = pd.read_csv(CYCLE2_DATA / f"{s}_master_augmented.tsv", sep="\t", low_memory=False)

    # Collect feature values per sample
    fig, axes = plt.subplots(3, 4, figsize=(16, 11))
    for idx, (feat_name, feat_id) in enumerate(METHYL_FEATURES):
        ax = axes[idx // 4, idx % 4]
        col = f"V6_off_meth_{feat_name}"
        plot_data = []
        plot_labels = []
        for s in SAMPLE_ORDER:
            df = samples_data[s]
            if col in df.columns:
                vals = df[col].dropna().values
                plot_data.append(vals)
                plot_labels.append(f"{s}\n(n={len(vals):,})")
            else:
                plot_data.append([])
                plot_labels.append(f"{s}\n(N/A)")
        bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True,
                         showfliers=False, widths=0.6)
        colors = ["#D97757", "#15803D", "#A16207", "#1E3A8A", "#6B7280"]
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        ax.set_title(f"{feat_id}: {feat_name}", fontsize=9.5, fontweight="bold")
        ax.tick_params(axis="x", labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Cycle 4 Stage 1 — 12 Methylation Features Distribution × 5 Samples\n"
                 "(Used in cycle 1 LR: F_A1-F_A5; Unused: F_B1-F_B7)",
                 fontsize=12, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(FIG / "stage1_features_boxplot.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {FIG}/stage1_features_boxplot.png")


if __name__ == "__main__":
    print("Generating Stage 1 figures …")
    fig_gap_matrix()
    fig_features_boxplot()
