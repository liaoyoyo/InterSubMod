#!/usr/bin/env python3
"""Cycle 4 Stage 2 — 4 priority figures for HTML.

Outputs:
    cycle4/figures/stage2_roi_bubble.png        ROI bubble chart (cost vs reward, size = challenge power)
    cycle4/figures/stage2_priority_ranking.png  13 candidates ranked horizontal bar
    cycle4/figures/stage2_decision_tree.png     Trial A/B/C decision branches
    cycle4/figures/stage2_h_coverage_matrix.png H_C4 × candidate coverage heatmap
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

CYCLE4 = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/methyl_augmented_filter_phase2/cycle4")
FIG = CYCLE4 / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# 13 candidates: (id, name, challenge, cost_days, reward, trial, H)
CANDIDATES = [
    ("O3",  "Interaction LR + LRT",            5, 1.0, 4, "A", "H_C4_INT"),
    ("F1",  "Methyl entropy per region",       4, 0.5, 4, "A", "H_C4_NEW_FEAT"),
    ("F2",  "Per-CpG × HP χ² Fisher",          4, 0.5, 3, "A", "H_C4_NEW_FEAT"),
    ("O1",  "Methyl entropy histogram + AUC",  3, 0.5, 3, "A", "H_C4_NEW_FEAT"),
    ("F4",  "Methyl variance per region",      3, 0.5, 2, "A", "H_C4_NEW_FEAT"),
    ("A_U1", "Random Forest",                  4, 1.5, 3, "B", "H_C4_NL"),
    ("A_U2", "XGBoost",                        4, 1.5, 3, "B", "H_C4_NL"),
    ("O4",  "Per-zone heterogeneous LR",       4, 2.0, 3, "C", "H_C4_ZONE"),
    ("F3",  "CpG density per region",          2, 0.25, 2, "A", "H_C4_NEW_FEAT"),
    ("F5",  "HP1/HP2 methyl correlation",      2, 0.5, 2, "A", "H_C4_NEW_FEAT"),
    ("O2",  "Per-CpG × HP × ALT χ²",           3, 1.0, 2, "A", "H_C4_NEW_FEAT"),
    ("A_U3", "Polynomial / Spline LR",         3, 1.5, 2, "B", "H_C4_NL"),
    ("A_U8", "Elastic Net",                    2, 1.0, 2, "—", "—"),
]

TRIAL_COLOR = {"A": "#15803D", "B": "#1E3A8A", "C": "#A16207", "—": "#6B7280"}


def fig_roi_bubble():
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for cid, name, ch, cost, reward, trial, h in CANDIDATES:
        col = TRIAL_COLOR[trial]
        size = ch * 80
        ax.scatter(cost, reward, s=size, c=col, alpha=0.75, edgecolors="black", linewidth=1, zorder=3)
        ax.annotate(f"{cid}\n{name[:18]}", (cost, reward),
                    xytext=(7, 5), textcoords="offset points", fontsize=8.5)

    ax.set_xlabel("Estimated cost (days)", fontsize=11)
    ax.set_ylabel("Predicted reward score (1-5)", fontsize=11)
    ax.set_xlim(-0.1, 2.4)
    ax.set_ylim(0.5, 5)
    ax.grid(alpha=0.3)
    handles = [mpatches.Patch(color=TRIAL_COLOR[t], label=f"Trial {t}") for t in ["A", "B", "C", "—"]]
    handles.append(mpatches.Patch(color="white", label="(bubble size ∝ challenge power 1-5)"))
    ax.legend(handles=handles, loc="lower right", fontsize=9)
    ax.set_title("Cycle 4 Stage 2 — 13 Candidates ROI Bubble\n"
                 "Lower-left + larger = higher ROI; Trial A green ROI 最高",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "stage2_roi_bubble.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {FIG}/stage2_roi_bubble.png")


def fig_priority_ranking():
    # Compute ROI score
    rows = []
    for cid, name, ch, cost, reward, trial, h in CANDIDATES:
        roi = (ch * reward) / cost if cost > 0 else 0
        rows.append({"id": cid, "name": name, "trial": trial, "roi": roi,
                     "cost": cost, "challenge": ch, "reward": reward, "H": h})
    df = pd.DataFrame(rows).sort_values("roi", ascending=True)

    fig, ax = plt.subplots(figsize=(11, 7.5))
    colors = [TRIAL_COLOR[t] for t in df["trial"]]
    ax.barh(df["id"] + " " + df["name"], df["roi"], color=colors, edgecolor="black", linewidth=0.4)
    for i, (_, r) in enumerate(df.iterrows()):
        ax.text(r["roi"] + 0.5, i, f"ROI={r['roi']:.1f}  ch={r['challenge']} cost={r['cost']:.2f}d rwd={r['reward']}",
                ha="left", va="center", fontsize=8)

    ax.set_xlabel("ROI score = challenge × reward / cost", fontsize=11)
    handles = [mpatches.Patch(color=TRIAL_COLOR[t], label=f"Trial {t}") for t in ["A", "B", "C", "—"]]
    ax.legend(handles=handles, loc="lower right", fontsize=9)
    ax.set_title("Cycle 4 Stage 2 — 13 Candidates Ranked by ROI\n"
                 "Top of list = highest ROI; F1/O3/F2/O1/F3 should go into Trial A",
                 fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "stage2_priority_ranking.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {FIG}/stage2_priority_ranking.png")


def fig_decision_tree():
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")

    # Trial A box
    ax.add_patch(mpatches.FancyBboxPatch((4.5, 6.2), 3, 1.2, boxstyle="round,pad=0.1",
                                          fc="#DCFCE7", ec="#15803D", linewidth=2))
    ax.text(6, 7.0, "Trial A (~1.5 day)", ha="center", fontsize=11, fontweight="bold", color="#15803D")
    ax.text(6, 6.55, "Interaction LR + 5 Python new features", ha="center", fontsize=9)

    # H_C4_INT branch
    ax.annotate("", xy=(3.5, 4.7), xytext=(5, 6.2),
                arrowprops=dict(arrowstyle="->", color="#1E3A8A", lw=1.5))
    ax.add_patch(mpatches.FancyBboxPatch((1.5, 4.0), 3.5, 0.7, boxstyle="round,pad=0.1",
                                          fc="#EFF6FF", ec="#1E3A8A", linewidth=1.5))
    ax.text(3.25, 4.35, "H_C4_INT PASS + H_C4_NEW_FEAT PASS", ha="center", fontsize=9, fontweight="bold", color="#1E3A8A")

    ax.annotate("", xy=(8.5, 4.7), xytext=(7, 6.2),
                arrowprops=dict(arrowstyle="->", color="#A16207", lw=1.5))
    ax.add_patch(mpatches.FancyBboxPatch((6.5, 4.0), 3.5, 0.7, boxstyle="round,pad=0.1",
                                          fc="#FEF7E6", ec="#A16207", linewidth=1.5))
    ax.text(8.25, 4.35, "1 of 2 PASS (marginal)", ha="center", fontsize=9, fontweight="bold", color="#A16207")

    # All FAIL
    ax.annotate("", xy=(11, 4.7), xytext=(7.5, 6.2),
                arrowprops=dict(arrowstyle="->", color="#C2410C", lw=1.5))
    ax.add_patch(mpatches.FancyBboxPatch((10, 4.0), 1.8, 0.7, boxstyle="round,pad=0.1",
                                          fc="#FEF2F0", ec="#C2410C", linewidth=1.5))
    ax.text(10.9, 4.35, "全 FAIL", ha="center", fontsize=9, fontweight="bold", color="#C2410C")

    # Trial B
    ax.annotate("", xy=(3.25, 2.7), xytext=(3.25, 4.0),
                arrowprops=dict(arrowstyle="->", color="#1E3A8A", lw=1.5))
    ax.add_patch(mpatches.FancyBboxPatch((1, 2.0), 4.5, 0.7, boxstyle="round,pad=0.1",
                                          fc="#DBEAFE", ec="#1E3A8A", linewidth=2))
    ax.text(3.25, 2.35, "Trial B RF/XGBoost (~2 day)", ha="center", fontsize=10, fontweight="bold", color="#1E3A8A")

    # Trial C
    ax.annotate("", xy=(3.25, 1.0), xytext=(3.25, 2.0),
                arrowprops=dict(arrowstyle="->", color="#A16207", lw=1.5))
    ax.add_patch(mpatches.FancyBboxPatch((1, 0.3), 4.5, 0.7, boxstyle="round,pad=0.1",
                                          fc="#FEF7E6", ec="#A16207", linewidth=2))
    ax.text(3.25, 0.65, "Trial C Per-zone LR (~2 day)", ha="center", fontsize=10, fontweight="bold", color="#A16207")

    # Marginal branch action
    ax.annotate("", xy=(8.25, 2.7), xytext=(8.25, 4.0),
                arrowprops=dict(arrowstyle="->", color="#A16207", lw=1.2, linestyle="dashed"))
    ax.add_patch(mpatches.FancyBboxPatch((6.5, 2.0), 3.5, 0.7, boxstyle="round,pad=0.1",
                                          fc="#FEF7E6", ec="#A16207", linewidth=1, linestyle="dashed"))
    ax.text(8.25, 2.35, "Trial B/C optional", ha="center", fontsize=9, color="#A16207")

    # All FAIL action
    ax.annotate("", xy=(10.9, 2.7), xytext=(10.9, 4.0),
                arrowprops=dict(arrowstyle="->", color="#C2410C", lw=1.5))
    ax.add_patch(mpatches.FancyBboxPatch((9.7, 1.8), 2.4, 0.9, boxstyle="round,pad=0.1",
                                          fc="#FEF2F0", ec="#C2410C", linewidth=2))
    ax.text(10.9, 2.45, "Cycle 4 NEGATIVE\nfinal; archive ", ha="center", fontsize=9, fontweight="bold", color="#C2410C")
    ax.text(10.9, 2.15, "methyl軌", ha="center", fontsize=9, color="#C2410C")

    ax.text(6, 7.7, "Cycle 4 Stage 2 — Trial A→B→C Serial Decision Tree (user 5/20 選擇)",
            ha="center", fontsize=13, fontweight="bold")

    fig.savefig(FIG / "stage2_decision_tree.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {FIG}/stage2_decision_tree.png")


def fig_h_coverage_matrix():
    """H_C4 × candidate coverage heatmap."""
    h_list = ["H_C4_INT", "H_C4_NL", "H_C4_ZONE", "H_C4_NEW_FEAT"]
    cand_list = ["O3", "F1", "F2", "O1", "F4", "A_U1", "A_U2", "O4", "F3", "F5", "O2"]

    mat = np.zeros((len(cand_list), len(h_list)))
    for i, cid in enumerate(cand_list):
        for c in CANDIDATES:
            if c[0] == cid:
                h = c[6]
                if h in h_list:
                    j = h_list.index(h)
                    mat[i, j] = c[2]  # challenge power
                break

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat, cmap="YlGn", aspect="auto", vmin=0, vmax=5)
    for i in range(len(cand_list)):
        for j in range(len(h_list)):
            v = mat[i, j]
            if v > 0:
                ax.text(j, i, f"{int(v)}", ha="center", va="center",
                        fontsize=11, fontweight="bold",
                        color="white" if v > 3 else "black")

    ax.set_xticks(range(len(h_list)))
    ax.set_xticklabels(h_list, fontsize=10, rotation=20)
    ax.set_yticks(range(len(cand_list)))
    ax.set_yticklabels(cand_list, fontsize=10)
    ax.set_title("Cycle 4 Stage 2 — H_C4 × Candidate Challenge Coverage\n"
                 "Each cell = challenge power of candidate against the H (0-5)",
                 fontsize=11, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, label="Challenge power")
    fig.tight_layout()
    fig.savefig(FIG / "stage2_h_coverage_matrix.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {FIG}/stage2_h_coverage_matrix.png")


if __name__ == "__main__":
    print("Generating Stage 2 figures …")
    fig_roi_bubble()
    fig_priority_ranking()
    fig_decision_tree()
    fig_h_coverage_matrix()
