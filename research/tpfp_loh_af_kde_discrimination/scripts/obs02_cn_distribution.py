#!/usr/bin/env python3
"""
obs02 | Per-sample CovM_used distribution, TP vs FP overlay.

Question: do FP and TP differ in CN profile?  The 4-boundary strategy F
(SEQC2-grounded) partitions CovM into T0-T4.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (PALETTE, apply_style, ensure_fig_dir, load_master)


CN_EDGES = [0.65, 0.99, 1.33, 1.82]


def main() -> None:
    apply_style()
    fig_dir = ensure_fig_dir()
    df = load_master()

    pf_samples = ["HCC1395", "HCC1937", "H2009", "COLO829"]

    fig, axes = plt.subplots(nrows=len(pf_samples), ncols=2,
                              figsize=(13, 3 * len(pf_samples)), sharex=True)

    edges = np.linspace(0, 3.0, 61)

    for r, smpl in enumerate(pf_samples):
        sub = df[(df["sample"] == smpl) & (df["mode"] == "paired_full")]
        tp = sub[sub["tp_label"] == 1]["CovM_used"].dropna()
        fp = sub[sub["tp_label"] == 0]["CovM_used"].dropna()
        ax = axes[r, 0]
        ax.hist(tp.clip(0, 3), bins=edges, alpha=0.72, density=True, color=PALETTE["tp"],
                label=f"TP n={len(tp):,}", edgecolor="white", linewidth=0.3)
        ax.hist(fp.clip(0, 3), bins=edges, alpha=0.62, density=True, color=PALETTE["fp"],
                label=f"FP n={len(fp):,}", edgecolor="white", linewidth=0.3)
        for e in CN_EDGES:
            ax.axvline(e, color=PALETTE["dark"], linestyle=":", linewidth=0.9, alpha=0.7)
        ax.set_title(f"{smpl} | paired_full", fontsize=11, loc="left", color=PALETTE["dark"])
        ax.set_ylabel("density")
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(True, alpha=0.25)

        ax = axes[r, 1]
        if smpl == "HCC1395":
            sub_to = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")]
            tp_to = sub_to[sub_to["tp_label"] == 1]["CovM_used"].dropna()
            fp_to = sub_to[sub_to["tp_label"] == 0]["CovM_used"].dropna()
            ax.hist(tp_to.clip(0, 3), bins=edges, alpha=0.72, density=True, color=PALETTE["tp"],
                    label=f"TP n={len(tp_to):,}", edgecolor="white", linewidth=0.3)
            ax.hist(fp_to.clip(0, 3), bins=edges, alpha=0.62, density=True, color=PALETTE["fp"],
                    label=f"FP n={len(fp_to):,}", edgecolor="white", linewidth=0.3)
            for e in CN_EDGES:
                ax.axvline(e, color=PALETTE["dark"], linestyle=":", linewidth=0.9, alpha=0.7)
            ax.set_title(f"{smpl} | HCC1395 TO mode", fontsize=11, loc="left", color=PALETTE["dark"])
            ax.legend(fontsize=9, loc="upper right")
            ax.grid(True, alpha=0.25)
        else:
            ax.axis("off")
            ax.text(0.5, 0.5, "(TO mode only HCC1395)", ha="center", va="center",
                    color=PALETTE["grey"], fontsize=10, transform=ax.transAxes)

    for ax in axes[-1]:
        ax.set_xlabel("CovM_used (coverage multiple)")

    fig.suptitle("obs02 | CovM distribution overlay (TP vs FP, per sample)",
                  fontsize=14, fontweight="bold", color=PALETTE["dark"], y=0.995)
    fig.text(0.5, 0.005,
              "Dotted vertical lines = F-strategy SEQC2-grounded CN boundaries [0.65, 0.99, 1.33, 1.82]",
              ha="center", color=PALETTE["grey"], fontsize=9)
    fig.tight_layout(rect=[0, 0.015, 1, 0.985])
    out = fig_dir / "obs02_cn_distribution_per_sample.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs02] wrote {out}")


if __name__ == "__main__":
    main()
