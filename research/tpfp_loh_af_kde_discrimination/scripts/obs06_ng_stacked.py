#!/usr/bin/env python3
"""
obs06 | Per-sample NG distribution stacked (TP + FP segments).

Question: does NG tile density differ between samples and modes?  This
helps validate whether NG-heavy cells are a consistent biology signal or
sample-dependent artefact.

Layout: 7 samples paired_full + HCC1395 TO = 8 bars, each stacked into NG=1/2/3/4
sub-bars (color = NG) with TP/FP split overlaid via hatch for FP.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (PALETTE, SAMPLE_ORDER, apply_style, ensure_fig_dir,
                         load_master)


NG_COLORS = {1: "#D55E00", 2: "#5B8DB7", 3: "#009E73", 4: "#5A3A7F"}


def build_stack(df: pd.DataFrame):
    g = df.groupby(["HPFineNGroups", "tp_label"]).size().unstack(fill_value=0)
    g = g.reindex([1, 2, 3, 4])
    tp_counts = g.get(1, pd.Series(0, index=[1, 2, 3, 4])).fillna(0)
    fp_counts = g.get(0, pd.Series(0, index=[1, 2, 3, 4])).fillna(0)
    return tp_counts.astype(int), fp_counts.astype(int)


def main() -> None:
    apply_style()
    fig_dir = ensure_fig_dir()
    df = load_master()

    groups = [(s, "paired_full") for s in SAMPLE_ORDER] + [("HCC1395", "to_pileup")]
    labels = [f"{s}\n{m}" for s, m in groups]

    fig, ax = plt.subplots(figsize=(14, 7))
    xs = np.arange(len(groups))
    bar_w = 0.35

    max_total = 0

    for i, (s, m) in enumerate(groups):
        sub = df[(df["sample"] == s) & (df["mode"] == m)]
        tp_c, fp_c = build_stack(sub)
        # TP stack
        bottom = 0
        for ng in [1, 2, 3, 4]:
            v = tp_c.get(ng, 0)
            ax.bar(i - bar_w / 2, v, width=bar_w, bottom=bottom,
                   color=NG_COLORS[ng], edgecolor=PALETTE["dark"], linewidth=0.3,
                   label=f"TP NG={ng}" if i == 0 else None)
            bottom += v
        # FP stack (striped)
        bottom2 = 0
        for ng in [1, 2, 3, 4]:
            v = fp_c.get(ng, 0)
            ax.bar(i + bar_w / 2, v, width=bar_w, bottom=bottom2,
                   color=NG_COLORS[ng], edgecolor=PALETTE["dark"], linewidth=0.3,
                   hatch="///", alpha=0.85,
                   label=f"FP NG={ng}" if i == 0 else None)
            bottom2 += v
        total = bottom + bottom2
        max_total = max(max_total, total)
        ax.text(i, total + 500, f"n={total:,}", ha="center", fontsize=8, color=PALETTE["dark"])

    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("region count")
    ax.set_title("obs06 | HPFineNGroups (NG) distribution per sample/mode — left bar TP, right bar FP (hatched)",
                  fontsize=12, fontweight="bold", color=PALETTE["dark"])
    ax.grid(True, axis="y", alpha=0.25)

    # Build dedupe legend
    handles, labs = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labs):
        if l not in seen:
            seen[l] = h
    ax.legend(seen.values(), seen.keys(), ncol=4, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.18), framealpha=0.9)

    fig.tight_layout()
    out = fig_dir / "obs06_ng_by_sample_stacked.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[obs06] wrote {out}")


if __name__ == "__main__":
    main()
