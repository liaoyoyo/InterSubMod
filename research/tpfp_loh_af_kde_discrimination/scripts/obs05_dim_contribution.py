#!/usr/bin/env python3
"""
obs05 | Single-dimension marginal contribution to TP rate (HCC1395 TO).

Question: if we collapse the 5D cube along each axis in turn, which axis
gives the largest swing in TP rate?  Answers "which dimension matters most?"

Output: 5 sub-panels (one per axis).  Each panel = bar of TP rate per category
with overlay of baseline dashed line; marginals have n annotated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (DATA_DIR, PALETTE, apply_style, ensure_fig_dir, load_master)


def main() -> None:
    apply_style()
    fig_dir = ensure_fig_dir()
    df = load_master()

    to = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")].copy()
    edges = [0.65, 0.99, 1.33, 1.82]
    to["cn_tier_F"] = pd.cut(to["CovM_used"], bins=[-np.inf] + edges + [np.inf],
                              labels=[f"T{i}" for i in range(len(edges) + 1)], include_lowest=True)
    to["nr_band"] = to["NumReads"].apply(
        lambda x: "low" if x < 60 else ("mid" if x < 120 else "high")
    )
    to["LOH_Subtype"] = to["LOH_Subtype"].fillna("None")
    baseline_tp = to["tp_label"].mean()

    dims = [
        ("LOH_Subtype", ["None", "LOH_Noise", "LOH_Weak", "LOH_Subclone", "LOH_Strong"]),
        ("AF_class", ["Extreme", "Intermediate", "Near-half"]),
        ("cn_tier_F", ["T0", "T1", "T2", "T3", "T4"]),
        ("HPFineNGroups", [1, 2, 3, 4]),
        ("nr_band", ["low", "mid", "high"]),
    ]

    fig, axes = plt.subplots(1, 5, figsize=(20, 5))

    for ax, (col, order) in zip(axes, dims):
        g = to.groupby(col, observed=True)
        agg = g.agg(n=("tp_label", "size"), n_tp=("tp_label", "sum"))
        agg["tp_rate"] = agg["n_tp"] / agg["n"]
        agg = agg.reindex(order)
        colors = [PALETTE["green"] if r > baseline_tp + 0.05 else (
                   PALETTE["red"] if r < baseline_tp - 0.05 else PALETTE["grey"])
                   for r in agg["tp_rate"]]
        bars = ax.bar(range(len(order)), agg["tp_rate"].fillna(0), color=colors,
                      edgecolor=PALETTE["dark"], linewidth=0.7)
        ax.axhline(baseline_tp, color=PALETTE["accent"], linestyle="--", linewidth=1.5,
                    label=f"baseline {baseline_tp:.3f}")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([str(x) for x in order], rotation=30, ha="right", fontsize=9)
        ax.set_title(col, fontsize=11, fontweight="bold", color=PALETTE["dark"])
        ax.set_ylabel("TP rate")
        ax.set_ylim(0, 1.0)
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(fontsize=8, loc="upper right")

        # Annotate n under each bar
        for i, (cat, r) in enumerate(agg.iterrows()):
            if pd.isna(r["n"]):
                continue
            n = int(r["n"]); rate = r["tp_rate"]
            ax.text(i, rate + 0.02, f"{rate:.2f}", ha="center", fontsize=8, color=PALETTE["dark"])
            ax.text(i, 0.02, f"n={n:,}", ha="center", fontsize=7, color=PALETTE["dark"])

    fig.suptitle(
        "obs05 | HCC1395 TO — TP rate marginal contribution per dimension (green ≥baseline+0.05, red ≤baseline−0.05)",
        fontsize=13, fontweight="bold", color=PALETTE["dark"], y=1.01,
    )
    fig.tight_layout()
    out = fig_dir / "obs05_dim_contribution.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs05] wrote {out}")


if __name__ == "__main__":
    main()
