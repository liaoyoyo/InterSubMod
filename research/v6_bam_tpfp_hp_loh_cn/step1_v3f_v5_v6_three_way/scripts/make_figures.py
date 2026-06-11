#!/usr/bin/env python3
"""Step 1 stage 8: figures.

Inputs:
- step1_master_three_way.tsv
- step1_trajectory.tsv

Outputs to figures/step1/:
- step1_three_panel_heatmap.png   (V3F / V5 / V6 NG × LOH side TP rate)
- step1_delta_heatmap.png         (3 delta heatmaps: V5-V3F, V6-V5, V6-V3F)
- step1_trajectory_sankey.png     (5-class trajectory bars)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch

# CJK font fallback (per InterSubMod feedback_matplotlib_cjk_font_rule)
plt.rcParams["font.family"] = ["DejaVu Sans", "Droid Sans Fallback", "Noto Sans CJK TC", "Noto Sans CJK SC"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 9

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"
FIG = STEP1 / "figures/step1"
FIG.mkdir(parents=True, exist_ok=True)

BAMS = ["V3F", "V5", "V6"]


def col(bam: str, flag: str, feat: str) -> str:
    return f"{bam}_{flag}_{feat}"


def build_panel(df: pd.DataFrame, bam: str):
    """Return a 2D matrix: rows=NG bin, cols=LOH side; values=TP rate."""
    sub = df.copy()
    sub["NG"] = pd.to_numeric(sub[col(bam, "off", "NG")], errors="coerce")
    # Derive loh_side_final from master.tsv `loh_side` column; UNKNOWN if missing.
    if "loh_side" in sub.columns:
        sub["loh_side_final"] = sub["loh_side"].fillna("UNKNOWN").replace("", "UNKNOWN")
    else:
        sub["loh_side_final"] = "UNKNOWN"
    sub = sub[sub["NG"].notna() & sub["loh_side_final"].notna()]
    sub["NG_bin"] = sub["NG"].clip(upper=5).astype(int)  # cap at 5+
    # only consider master-joined rows so loh_side is meaningful
    sub = sub[sub["master_join_ok"] == 1]
    sides = ["Inner", "Outer", "UNKNOWN"]
    sub["loh_side_disp"] = sub["loh_side_final"].apply(
        lambda x: x if x in sides else "UNKNOWN"
    )

    mat = pd.DataFrame(index=range(0, 6), columns=sides, dtype=float)
    n_mat = pd.DataFrame(index=range(0, 6), columns=sides, dtype=int)
    for ng in range(0, 6):
        for side in sides:
            cell = sub[(sub["NG_bin"] == ng) & (sub["loh_side_disp"] == side)]
            n_total = len(cell)
            n_tp = (cell["label"] == "TP").sum()
            mat.loc[ng, side] = n_tp / n_total if n_total else np.nan
            n_mat.loc[ng, side] = n_total
    return mat, n_mat


def plot_three_panel(df: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    mats = []
    for ax, bam in zip(axes, BAMS):
        mat, n_mat = build_panel(df, bam)
        mats.append((bam, mat, n_mat))
    vmin = min(m[1].values[np.isfinite(m[1].values)].min() for m in mats if np.any(np.isfinite(m[1].values)))
    vmax = max(m[1].values[np.isfinite(m[1].values)].max() for m in mats if np.any(np.isfinite(m[1].values)))
    for ax, (bam, mat, n_mat) in zip(axes, mats):
        im = ax.imshow(mat.values, vmin=vmin, vmax=vmax, cmap="RdYlGn", aspect="auto")
        ax.set_title(f"{bam} BAM — TP rate by NG × LOH side")
        ax.set_xticks(range(len(mat.columns)))
        ax.set_xticklabels(mat.columns)
        ax.set_yticks(range(len(mat.index)))
        ax.set_yticklabels([f"NG={i}" + ("+" if i == 5 else "") for i in mat.index])
        ax.set_xlabel("LOH side")
        if ax is axes[0]:
            ax.set_ylabel("HPFineNGroups (NG_off)")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.iloc[i, j]
                n = n_mat.iloc[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f"{v:.2f}\n(n={n})", ha="center", va="center",
                            color="white" if v < (vmin + vmax) / 2 else "black", fontsize=7)
                else:
                    ax.text(j, i, "—", ha="center", va="center", color="gray", fontsize=8)
    fig.suptitle("Step 1: V3F / V5 / V6 three-panel — TP rate by NG × LOH side (HCC1395 paired-pileup, master-joined)",
                 fontsize=11, fontweight="bold")
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("TP rate", fontsize=9)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def plot_delta(df: pd.DataFrame, out: Path) -> None:
    """Plot 3 delta heatmaps (NG × LOH side, TP rate delta)."""
    panels = {bam: build_panel(df, bam)[0] for bam in BAMS}
    deltas = [
        ("V5 − V3F", panels["V5"] - panels["V3F"]),
        ("V6 − V5", panels["V6"] - panels["V5"]),
        ("V6 − V3F", panels["V6"] - panels["V3F"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    # Use symmetric scale
    all_vals = np.concatenate([d[1].values.flatten() for d in deltas])
    all_vals = all_vals[np.isfinite(all_vals)]
    vmax = max(abs(all_vals.min()), abs(all_vals.max())) if len(all_vals) else 0.1
    for ax, (title, mat) in zip(axes, deltas):
        im = ax.imshow(mat.values, vmin=-vmax, vmax=vmax, cmap="RdBu_r", aspect="auto")
        ax.set_title(f"Δ TP rate ({title})")
        ax.set_xticks(range(len(mat.columns)))
        ax.set_xticklabels(mat.columns)
        ax.set_yticks(range(len(mat.index)))
        ax.set_yticklabels([f"NG={i}" + ("+" if i == 5 else "") for i in mat.index])
        ax.set_xlabel("LOH side")
        if ax is axes[0]:
            ax.set_ylabel("HPFineNGroups (NG_off)")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.iloc[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f"{v:+.3f}", ha="center", va="center",
                            color="black", fontsize=8)
                else:
                    ax.text(j, i, "—", ha="center", va="center", color="gray", fontsize=8)
    fig.suptitle("Step 1: Δ TP rate heatmaps — V3F → V5 → V6 (HCC1395 paired-pileup)",
                 fontsize=11, fontweight="bold")
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("Δ TP rate", fontsize=9)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def plot_trajectory_sankey(traj: pd.DataFrame, out: Path) -> None:
    """Plot a Sankey-like flow: source = label (TP/FP), middle = trajectory class.
    Uses stacked bar approximation since true matplotlib Sankey is verbose.
    """
    counts = traj.groupby(["label", "class"]).size().unstack(fill_value=0)
    # Order classes
    class_order = ["A_both_improve", "B_only_V5_improve", "C_only_V6_improve",
                   "D_no_change", "E_reverse_or_decrease", "MISSING"]
    cols_present = [c for c in class_order if c in counts.columns]
    counts = counts[cols_present]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: stacked bar
    ax = axes[0]
    colors = {"A_both_improve": "#1b9e77", "B_only_V5_improve": "#7570b3",
              "C_only_V6_improve": "#d95f02", "D_no_change": "#999999",
              "E_reverse_or_decrease": "#e7298a", "MISSING": "#666666"}
    bottoms = np.zeros(len(counts.index))
    for cls in counts.columns:
        vals = counts[cls].values
        ax.bar(range(len(counts.index)), vals, bottom=bottoms,
               color=colors.get(cls, "gray"), label=cls.replace("_", " "))
        # annotate
        for i, v in enumerate(vals):
            if v > 0:
                ax.text(i, bottoms[i] + v / 2, f"{v}\n({v/counts.sum(axis=1).iloc[i]:.1%})",
                        ha="center", va="center", fontsize=8,
                        color="white" if cls != "D_no_change" else "black")
        bottoms += vals
    ax.set_xticks(range(len(counts.index)))
    ax.set_xticklabels(counts.index)
    ax.set_ylabel("Region count")
    ax.set_title("Per-region NG trajectory class — V3F→V5→V6")
    ax.legend(loc="upper right", fontsize=8)

    # Right: fraction normalized
    ax = axes[1]
    frac = counts.div(counts.sum(axis=1), axis=0)
    bottoms = np.zeros(len(frac.index))
    for cls in frac.columns:
        vals = frac[cls].values
        ax.bar(range(len(frac.index)), vals, bottom=bottoms,
               color=colors.get(cls, "gray"), label=cls.replace("_", " "))
        for i, v in enumerate(vals):
            if v > 0.02:
                ax.text(i, bottoms[i] + v / 2, f"{v:.1%}",
                        ha="center", va="center", fontsize=8,
                        color="white" if cls != "D_no_change" else "black")
        bottoms += vals
    ax.set_xticks(range(len(frac.index)))
    ax.set_xticklabels(frac.index)
    ax.set_ylabel("Fraction within label")
    ax.set_ylim(0, 1)
    ax.set_title("Per-region NG trajectory — normalized")

    fig.suptitle("Step 1: V3F → V5 → V6 per-region NG trajectory classification",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}", flush=True)


def main() -> int:
    print("[Step 1.8] loading data ...", flush=True)
    df = pd.read_csv(STEP1 / "step1_master_three_way.tsv", sep="\t", low_memory=False)
    traj = pd.read_csv(STEP1 / "step1_trajectory.tsv", sep="\t", low_memory=False)
    print(f"  master rows={len(df)}, trajectory rows={len(traj)}", flush=True)

    plot_three_panel(df, FIG / "step1_three_panel_heatmap.png")
    plot_delta(df, FIG / "step1_delta_heatmap.png")
    plot_trajectory_sankey(traj, FIG / "step1_trajectory_sankey.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
