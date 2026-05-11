#!/usr/bin/env python3
"""
obs01 | Per-sample AF distribution, TP vs FP overlay.

Question: In which AF ranges do FPs pile up vs TPs? Does the answer change
between paired_full and HCC1395 TO?

Layout: 2 columns (paired_full | HCC1395 TO) x 4 rows (HCC1395, HCC1937,
H2009, COLO829).  Samples with <100 FP are skipped (low power).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (PALETTE, apply_style, ensure_fig_dir, load_master)


def main() -> None:
    apply_style()
    fig_dir = ensure_fig_dir()
    df = load_master()

    # samples with enough FP to plot AF shape (paired_full)
    pf_samples = ["HCC1395", "HCC1937", "H2009", "COLO829"]
    to_samples = ["HCC1395"]

    fig, axes = plt.subplots(nrows=max(len(pf_samples), 1), ncols=2,
                              figsize=(13, 3 * len(pf_samples)), sharex=True)

    edges = np.linspace(0, 1.0, 41)

    for r, smpl in enumerate(pf_samples):
        # paired_full
        sub = df[(df["sample"] == smpl) & (df["mode"] == "paired_full")]
        tp = sub[sub["tp_label"] == 1]["AF"].dropna()
        fp = sub[sub["tp_label"] == 0]["AF"].dropna()
        ax = axes[r, 0]
        ax.hist(tp, bins=edges, alpha=0.72, density=True, color=PALETTE["tp"],
                label=f"TP n={len(tp):,}", edgecolor="white", linewidth=0.3)
        ax.hist(fp, bins=edges, alpha=0.62, density=True, color=PALETTE["fp"],
                label=f"FP n={len(fp):,}", edgecolor="white", linewidth=0.3)
        ax.set_title(f"{smpl} | paired_full", fontsize=11, loc="left", color=PALETTE["dark"])
        ax.set_ylabel("density")
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(True, alpha=0.25)
        ax.axvspan(0.4, 0.6, color=PALETTE["blue"], alpha=0.08)

        # HCC1395 TO mode only populated for r=0 (HCC1395)
        ax = axes[r, 1]
        if smpl == "HCC1395":
            sub_to = df[(df["sample"] == "HCC1395") & (df["mode"] == "to_pileup")]
            tp_to = sub_to[sub_to["tp_label"] == 1]["AF"].dropna()
            fp_to = sub_to[sub_to["tp_label"] == 0]["AF"].dropna()
            ax.hist(tp_to, bins=edges, alpha=0.72, density=True, color=PALETTE["tp"],
                    label=f"TP n={len(tp_to):,}", edgecolor="white", linewidth=0.3)
            ax.hist(fp_to, bins=edges, alpha=0.62, density=True, color=PALETTE["fp"],
                    label=f"FP n={len(fp_to):,}", edgecolor="white", linewidth=0.3)
            ax.set_title(f"{smpl} | HCC1395 TO mode (main FP library)", fontsize=11,
                         loc="left", color=PALETTE["dark"])
            ax.legend(fontsize=9, loc="upper right")
            ax.grid(True, alpha=0.25)
            ax.axvspan(0.4, 0.6, color=PALETTE["blue"], alpha=0.08)
        else:
            ax.axis("off")
            ax.text(0.5, 0.5, "(TO mode only HCC1395)", ha="center", va="center",
                    color=PALETTE["grey"], fontsize=10, transform=ax.transAxes)

    for ax in axes[-1]:
        ax.set_xlabel("AF (allele frequency)")

    fig.suptitle("obs01 | AF distribution overlay (TP vs FP, per sample)",
                  fontsize=14, fontweight="bold", color=PALETTE["dark"], y=0.995)
    fig.text(0.5, 0.005, "Blue band = Near-half window [0.4, 0.6]",
              ha="center", color=PALETTE["grey"], fontsize=9)
    fig.tight_layout(rect=[0, 0.015, 1, 0.985])
    out = fig_dir / "obs01_af_distribution_per_sample.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs01] wrote {out}")


if __name__ == "__main__":
    main()
