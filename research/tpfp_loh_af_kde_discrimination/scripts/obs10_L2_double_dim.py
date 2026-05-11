#!/usr/bin/env python3
"""
obs10 | L2 — Double-dimension per-sample heatmap (V2 2026-04-22).

For each of the 10 dim pairs, produce TWO figures:
  - `L2_{x}_x_{y}_TO.png`     — 6 TO samples in 2x3 grid
  - `L2_{x}_x_{y}_paired.png` — 7 paired_full samples in 3x3 grid

Each figure:
  - Research-question suptitle
  - Per-panel metadata badge (tier, FP, baseline, n)
  - Computed takeaway caption box (highlights best cell + consensus)

Writes to ``figures/new/L2/``.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (
    DIM_ORDERS, DIM_SHORT_LABEL, FIG_NEW_DIR, PALETTE,
    annotate_power_tier, apply_style, build_5d_key, compute_tp_rate_by_key,
    load_master, mode_display_name, mode_subset_grid, research_suptitle,
    sample_label, takeaway_caption,
)

L2_DIR = FIG_NEW_DIR / "L2"

DIM_PAIRS = list(combinations(
    ["LOH_Subtype", "AF_class", "cn_tier_F", "HPFineNGroups", "nr_band"], 2,
))


def _pair_question(dim_x: str, dim_y: str, mode: str) -> str:
    mode_label = mode_display_name(mode)
    return (f"{DIM_SHORT_LABEL[dim_x]} × {DIM_SHORT_LABEL[dim_y]} 聯合切分在 {mode_label} "
            f"是否呈現跨樣本一致的 TP 熱區？")


def _pair_takeaway(dim_x: str, dim_y: str, mode: str, df5: pd.DataFrame) -> tuple[str, str]:
    """Compute takeaway: how many samples have a cell with TP>0.90 at which (x,y) slot?"""
    sm_list = [sm for sm in df5[df5["mode"] == mode][["sample", "mode"]].drop_duplicates().values.tolist()]
    x_order = DIM_ORDERS[dim_x]
    y_order = DIM_ORDERS[dim_y]

    # Count per-cell: how many samples have (n>=20 AND tp_rate>0.90)
    cell_hits: dict[tuple, int] = {}
    for sample, m in sm_list:
        sub = df5[(df5["sample"] == sample) & (df5["mode"] == m)]
        agg = compute_tp_rate_by_key(sub, [dim_x, dim_y], min_n=20)
        for _, row in agg.iterrows():
            key = (row[dim_x], row[dim_y])
            if pd.notna(row["tp_rate"]) and row["tp_rate"] > 0.90:
                cell_hits[key] = cell_hits.get(key, 0) + 1

    if not cell_hits:
        return (f"○ Takeaway：{mode_display_name(mode)} 下 (dim_x, dim_y) 無任何 cell TP>0.90。", "grey")

    best_key, best_hits = max(cell_hits.items(), key=lambda x: x[1])
    total_samples = len(sm_list)
    n_high_cells = sum(1 for v in cell_hits.values() if v >= max(2, total_samples // 2))

    if mode == "to_pileup":
        if best_hits >= 5:
            colour = "green"
            msg = (f"✓ Takeaway：({DIM_SHORT_LABEL[dim_x]}={best_key[0]}, "
                   f"{DIM_SHORT_LABEL[dim_y]}={best_key[1]}) 在 {best_hits}/{total_samples} TO 樣本達 TP>0.90；"
                   f"{n_high_cells} cells 達 ≥{max(2, total_samples//2)} 樣本一致。跨樣本共識強。")
        elif best_hits >= 3:
            colour = "blue"
            msg = (f"◐ Takeaway：最佳 cell ({best_key[0]}, {best_key[1]}) 在 {best_hits}/{total_samples} "
                   f"TO 樣本 TP>0.90。跨樣本共識中等；{total_samples - best_hits} 樣本例外。")
        else:
            colour = "red"
            msg = (f"⚠ Takeaway：無 cell 在 ≥3/{total_samples} TO 樣本達 TP>0.90。"
                   f"此 2D 切分**不具跨樣本共識**，不建議作為 filter 軸。")
    else:
        # paired baseline already ≥ 0.94 → nearly all cells high, not informative
        colour = "grey"
        msg = (f"○ Takeaway：paired_full baseline 已 0.94-0.99 飽和；"
               f"{n_high_cells} cells TP>0.90 符合預期但無區分力。L2 在 paired mode 視覺貢獻低，詳見 §7.5 consistency。")
    return (msg, colour)


def plot_pair_mode(dim_x: str, dim_y: str, mode: str, df5: pd.DataFrame) -> Path:
    x_order = DIM_ORDERS[dim_x]
    y_order = DIM_ORDERS[dim_y]
    fig, axes, sm_list = mode_subset_grid(mode)
    research_suptitle(
        fig,
        _pair_question(dim_x, dim_y, mode),
        context="cell color = TP rate (RdYlGn) · cell text = TP% / n · grey cell = n<20 · 紅框 = FP<100",
        y=0.975, fontsize=18,
    )

    norm = Normalize(vmin=0.3, vmax=1.0)
    cmap = plt.get_cmap("RdYlGn")
    cmap.set_bad(color=PALETTE["grey"], alpha=0.4)
    last_im = None

    for ax, (sample, m) in zip(axes[:len(sm_list)], sm_list):
        sub = df5[(df5["sample"] == sample) & (df5["mode"] == m)]
        if len(sub) == 0:
            ax.axis("off")
            continue
        baseline = sub["tp_label"].mean()
        n_total = len(sub)
        agg = compute_tp_rate_by_key(sub, [dim_x, dim_y], min_n=20)
        pivot_rate = agg.pivot(index=dim_y, columns=dim_x, values="tp_rate") \
            .reindex(index=y_order, columns=x_order)
        pivot_n = agg.pivot(index=dim_y, columns=dim_x, values="n") \
            .reindex(index=y_order, columns=x_order).fillna(0).astype(int)

        data = pivot_rate.values.astype(float)
        masked = np.ma.masked_invalid(data)
        im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="auto")
        last_im = im

        ax.set_xticks(range(len(x_order)))
        ax.set_xticklabels([str(x) for x in x_order], rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(y_order)))
        ax.set_yticklabels([str(y) for y in y_order], fontsize=8)
        ax.set_title(sample_label(sample, m), fontsize=11, fontweight="bold",
                     color=PALETTE["dark"])
        ax.set_xlabel(DIM_SHORT_LABEL[dim_x], fontsize=9)
        ax.set_ylabel(DIM_SHORT_LABEL[dim_y], fontsize=9)

        for yi in range(len(y_order)):
            for xi in range(len(x_order)):
                n = int(pivot_n.iat[yi, xi])
                rate = pivot_rate.iat[yi, xi]
                if n == 0:
                    continue
                if np.isnan(rate):
                    ax.text(xi, yi, f"n={n}\n<20", ha="center", va="center",
                            fontsize=6, color=PALETTE["dark"])
                else:
                    colour = "black" if rate > 0.55 else "white"
                    ax.text(xi, yi, f"{rate*100:.0f}%\nn={n}",
                            ha="center", va="center", fontsize=6.5, color=colour)

        annotate_power_tier(ax, sample, m, baseline_tp=baseline, n_total=n_total)

    # shared colorbar
    cbar_ax = fig.add_axes([0.93, 0.12, 0.012, 0.78])
    if last_im is not None:
        cbar = fig.colorbar(last_im, cax=cbar_ax)
        cbar.set_label("TP rate", fontsize=9)

    takeaway_text, colour = _pair_takeaway(dim_x, dim_y, mode, df5)
    takeaway_caption(fig, takeaway_text, color=colour, y=0.018, fontsize=12)

    fig.tight_layout(rect=[0, 0.065, 0.92, 0.93])
    mode_tag = "TO" if mode == "to_pileup" else "paired"
    out = L2_DIR / f"L2_{DIM_SHORT_LABEL[dim_x].lower()}_x_{DIM_SHORT_LABEL[dim_y].lower()}_{mode_tag}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    apply_style()
    L2_DIR.mkdir(parents=True, exist_ok=True)
    df = load_master()
    df5 = build_5d_key(df)
    df5 = df5[df5["HPFineNGroups"].isin(DIM_ORDERS["HPFineNGroups"])]
    for dx, dy in DIM_PAIRS:
        for mode in ("to_pileup", "paired_full"):
            out = plot_pair_mode(dx, dy, mode, df5)
            print(f"[obs10] wrote {out}")


if __name__ == "__main__":
    main()
