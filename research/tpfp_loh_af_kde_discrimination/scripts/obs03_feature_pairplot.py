#!/usr/bin/env python3
"""
obs03 | Feature pairplot (AF, CovM_used, NumReads, AlleleDelta) — HCC1395 TO.

Question: are TP/FP structurally separable in 2-D projections of these
4 features?  Focus on HCC1395 TO (11,606 FP, main FP library).

Layout: 4x4 grid.  Diagonals = kernel density; off-diagonals = scatter with
transparency (TP over FP or FP over TP by sample size).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (PALETTE, apply_style, ensure_fig_dir, load_master)


FEATURES = [
    ("AF", (0, 1), "AF"),
    ("CovM_used", (0, 3), "CovM_used"),
    ("NumReads", (0, 300), "NumReads"),
    ("AlleleDelta", (-1, 1), "AlleleDelta"),
]


def main() -> None:
    apply_style()
    fig_dir = ensure_fig_dir()
    df = load_master()

    sub = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")].copy()
    tp = sub[sub["tp_label"] == 1]
    fp = sub[sub["tp_label"] == 0]

    n_feat = len(FEATURES)
    fig, axes = plt.subplots(n_feat, n_feat, figsize=(14, 14))

    for i, (fi, xlim_i, lbl_i) in enumerate(FEATURES):
        for j, (fj, xlim_j, lbl_j) in enumerate(FEATURES):
            ax = axes[i, j]
            if i == j:
                tp_vals = tp[fi].dropna()
                fp_vals = fp[fi].dropna()
                bins = np.linspace(xlim_i[0], xlim_i[1], 50)
                ax.hist(tp_vals.clip(*xlim_i), bins=bins, density=True,
                        alpha=0.7, color=PALETTE["tp"], label=f"TP n={len(tp_vals):,}",
                        edgecolor="white", linewidth=0.2)
                ax.hist(fp_vals.clip(*xlim_i), bins=bins, density=True,
                        alpha=0.55, color=PALETTE["fp"], label=f"FP n={len(fp_vals):,}",
                        edgecolor="white", linewidth=0.2)
                if i == 0:
                    ax.legend(fontsize=8, loc="upper right")
                ax.set_xlim(*xlim_i)
                ax.grid(True, alpha=0.25)
            else:
                # Scatter: FP then TP (so TP dots on top)
                ax.scatter(fp[fj].clip(*xlim_j), fp[fi].clip(*xlim_i),
                           s=3, alpha=0.18, color=PALETTE["fp"], label="FP")
                # thin TP to 5k for perf
                tp_s = tp.sample(min(5000, len(tp)), random_state=0) if len(tp) > 5000 else tp
                ax.scatter(tp_s[fj].clip(*xlim_j), tp_s[fi].clip(*xlim_i),
                           s=3, alpha=0.22, color=PALETTE["tp"], label="TP")
                ax.set_xlim(*xlim_j)
                ax.set_ylim(*xlim_i)
                ax.grid(True, alpha=0.2)
            if i == n_feat - 1:
                ax.set_xlabel(lbl_j, fontsize=10)
            if j == 0:
                ax.set_ylabel(lbl_i, fontsize=10)

    fig.suptitle("obs03 | Feature pairplot HCC1395 TO (AF × CovM × NumReads × AlleleDelta)",
                  fontsize=14, fontweight="bold", color=PALETTE["dark"], y=0.995)
    fig.tight_layout()
    out = fig_dir / "obs03_feature_pairplot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs03] wrote {out}")


if __name__ == "__main__":
    main()
