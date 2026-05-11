#!/usr/bin/env python3
"""
obs07 | Pareto envelope trajectory — zoom into top cells and mark 3 operating points.

Question: how does the cumulative envelope transition from Top-17 (96.1%)
through Top-28 (90.4%) to Top-40?  Visualize the slope change and the 3
key operating points, overlaid with S3 and S5 for reference.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import DATA_DIR, PALETTE, apply_style, ensure_fig_dir


def main() -> None:
    apply_style()
    fig_dir = ensure_fig_dir()

    env = pd.read_csv(DATA_DIR / "tpfp_5d_cumulative_envelope_HCC1395_TO.tsv", sep="\t")
    env = env.reset_index(drop=True)
    env["k"] = env.index + 1

    # Baseline (HCC1395 TO)
    baseline_tp = 0.7107   # from step6 report
    # S3 / S5 measured
    s3_prec, s3_recall = 0.9550, 0.0128
    s5_prec, s5_recall = 0.9180, 0.0291

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Panel A: purity vs k with operating points
    ax = axes[0]
    ax.plot(env["k"], env["cum_tp_rate"], color=PALETTE["dark"], linewidth=2,
             label="cumulative TP rate (purity)")
    ax.axhline(baseline_tp, color=PALETTE["accent"], linestyle="--", linewidth=1.5,
                label=f"HCC1395 TO baseline {baseline_tp:.3f}")

    # Key operating points
    markers = [(17, "Top-17", PALETTE["green"]),
               (28, "Top-28", PALETTE["blue"]),
               (32, "Top-32 (probe)", PALETTE["accent"])]
    for k, lbl, col in markers:
        if k - 1 < len(env):
            r = env.iloc[k - 1]
            ax.scatter([k], [r["cum_tp_rate"]], s=90, color=col,
                        edgecolor=PALETTE["dark"], zorder=5,
                        label=f"{lbl}: p={r['cum_tp_rate']:.3f}, recall={r['cum_tp_recall']:.2%}")
            ax.axvline(k, color=col, linestyle=":", alpha=0.4, linewidth=1)
    ax.set_xlabel("k = number of top-ranked cells accumulated")
    ax.set_ylabel("cumulative TP rate (purity)")
    ax.set_title("Panel A | Purity decay as top-k cells accumulate", fontsize=11, fontweight="bold",
                  color=PALETTE["dark"], loc="left")
    ax.set_xlim(0, min(60, len(env)))
    ax.set_ylim(0.6, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    # Panel B: operating point trajectory (recall vs purity) with S3/S5 overlaid
    ax = axes[1]
    ax.plot(env["cum_tp_recall"], env["cum_tp_rate"],
             color=PALETTE["dark"], linewidth=2, label="5D cube cumulative envelope")
    for k, lbl, col in markers:
        if k - 1 < len(env):
            r = env.iloc[k - 1]
            ax.scatter([r["cum_tp_recall"]], [r["cum_tp_rate"]], s=110, color=col,
                        edgecolor=PALETTE["dark"], zorder=5,
                        label=f"{lbl}")
            ax.annotate(f"{lbl}\n({r['cum_tp_recall']:.2%}, {r['cum_tp_rate']:.3f})",
                         (r["cum_tp_recall"], r["cum_tp_rate"]),
                         textcoords="offset points", xytext=(6, -18), fontsize=8,
                         color=PALETTE["dark"])
    # Biology scheme overlay
    ax.scatter([s3_recall], [s3_prec], s=160, color=PALETTE["red"], marker="P",
                edgecolor=PALETTE["dark"], zorder=6, label=f"S3 Diploid Het ({s3_recall:.2%}, {s3_prec:.3f})")
    ax.scatter([s5_recall], [s5_prec], s=160, color=PALETTE["accent"], marker="X",
                edgecolor=PALETTE["dark"], zorder=6, label=f"S5 Combo ({s5_recall:.2%}, {s5_prec:.3f})")
    ax.axhline(baseline_tp, color=PALETTE["accent"], linestyle="--", linewidth=1.5, alpha=0.6)
    ax.set_xlabel("cumulative TP recall")
    ax.set_ylabel("cumulative TP rate (purity)")
    ax.set_title("Panel B | Envelope vs biology schemes (S3/S5 Pareto-dominated)",
                  fontsize=11, fontweight="bold", color=PALETTE["dark"], loc="left")
    ax.set_xlim(0, 0.3)
    ax.set_ylim(0.6, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    fig.suptitle("obs07 | Pareto envelope trajectory — Top-17/28/32 operating points",
                  fontsize=14, fontweight="bold", color=PALETTE["dark"], y=1.01)
    fig.tight_layout()
    out = fig_dir / "obs07_pareto_trajectory.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs07] wrote {out}")


if __name__ == "__main__":
    main()
