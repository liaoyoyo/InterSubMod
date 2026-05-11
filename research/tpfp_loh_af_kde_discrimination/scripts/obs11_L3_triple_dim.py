#!/usr/bin/env python3
"""
obs11 | L3 — Triple-dimension per-sample heatmap (4 figures).

For each of 4 high-value trios, produce one figure with 8 outer panels
(sample-mode) x inner col facets (the 3rd dim).  Inside each inner facet is
a heatmap of (dim_row x dim_col) coloured by TP rate.

Trios (following plan §2.3):
  - LOH x AF x CN      -> ``L3_loh_af_cn.png``  (fig_v2_1 per-sample extension)
  - LOH x AF x NG      -> ``L3_loh_af_ng.png``
  - AF x CN x NR       -> ``L3_af_cn_nr.png``
  - LOH x CN x NG      -> ``L3_loh_cn_ng.png``

Layout: figure = (n_samples=8) outer rows x (inner_facets) cols.  We use the
transposed convention used by ``fig_v2_1``: sample = outer row, inner
dimension = outer col within the row.  A 2x4 macro grid of outer panels,
each containing a strip of ``len(z_order)`` heatmaps.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (
    DIM_ORDERS, DIM_SHORT_LABEL, FIG_NEW_DIR, PALETTE, SAMPLE_FP_META,
    SAMPLE_MODE_ORDER, apply_style, build_5d_key, compute_tp_rate_by_key,
    load_master, research_suptitle, sample_label, takeaway_caption,
)

L3_DIR = FIG_NEW_DIR / "L3"

TRIOS: list[tuple[str, str, str]] = [
    ("LOH_Subtype", "AF_class", "cn_tier_F"),   # core — fig_v2_1 per-sample
    ("LOH_Subtype", "AF_class", "HPFineNGroups"),
    ("AF_class", "cn_tier_F", "nr_band"),
    ("LOH_Subtype", "cn_tier_F", "HPFineNGroups"),
]

# Research questions per trio (2026-04-22 V2 update)
TRIO_QUESTIONS = {
    ("LOH_Subtype", "AF_class", "cn_tier_F"):
        "fig_v2_1 (LOH × AF × CN) 的高 TP 熱區是否在全 13 sample-mode 重現？",
    ("LOH_Subtype", "AF_class", "HPFineNGroups"):
        "LOH × AF × NG 聯合下，NG=3 的 HCC1395_TO 高 TP 訊號是否為 per-sample 特例？",
    ("AF_class", "cn_tier_F", "nr_band"):
        "NR × CN × AF 聯合下，high NR + Near-half + CN≈2 是否跨樣本穩定？",
    ("LOH_Subtype", "cn_tier_F", "HPFineNGroups"):
        "LOH × CN × NG 聯合下，(LOH_Strong, T2, NG=3) 是否為跨樣本可複製熱區？",
}


def _trio_takeaway(dim_row: str, dim_col: str, dim_outer: str, df5: pd.DataFrame) -> tuple[str, str]:
    """Compute takeaway: how many samples have (row, col, outer) with TP>0.85?"""
    cell_hits: dict[tuple, int] = {}
    total_samples = 0
    for sample, mode in SAMPLE_MODE_ORDER:
        sub = df5[(df5["sample"] == sample) & (df5["mode"] == mode)]
        if len(sub) < 50:
            continue
        total_samples += 1
        agg = compute_tp_rate_by_key(sub, [dim_row, dim_col, dim_outer], min_n=20)
        for _, row in agg.iterrows():
            if pd.notna(row["tp_rate"]) and row["tp_rate"] > 0.85:
                key = (row[dim_row], row[dim_col], row[dim_outer])
                cell_hits[key] = cell_hits.get(key, 0) + 1

    if not cell_hits:
        return ("○ Takeaway：3D 組合下無任何 cell 達 TP>0.85 且跨樣本可複製。", "grey")

    # Top-3
    top3 = sorted(cell_hits.items(), key=lambda x: x[1], reverse=True)[:3]
    top_str = "; ".join(f"({k[0]}, {k[1]}, {k[2]})={v}/{total_samples}" for k, v in top3)
    best_hits = top3[0][1]

    if best_hits >= total_samples * 0.75:
        colour = "green"
        msg = f"✓ Takeaway：Top 3D cells 最高 {best_hits}/{total_samples} 樣本達 TP>0.85。Top-3: {top_str}。跨樣本共識強。"
    elif best_hits >= total_samples * 0.5:
        colour = "blue"
        msg = f"◐ Takeaway：Top cell 在 {best_hits}/{total_samples} 樣本一致。Top-3: {top_str}。中等共識。"
    else:
        colour = "red"
        msg = f"⚠ Takeaway：Top cell 僅 {best_hits}/{total_samples} 樣本一致。Top-3: {top_str}。3D 組合跨樣本泛化性弱。"
    return (msg, colour)


def plot_trio(dim_row: str, dim_col: str, dim_outer: str, df5: pd.DataFrame) -> Path:
    row_order = DIM_ORDERS[dim_row]
    col_order = DIM_ORDERS[dim_col]
    outer_order = DIM_ORDERS[dim_outer]

    n_outer = len(outer_order)
    n_samples = len(SAMPLE_MODE_ORDER)

    # Dynamic macro layout based on n_samples (2026-04-22 update: now 13 combos).
    # For ≤8 samples: 4×2; 9-12: 6×2; 13-15: 5×3; ≥16: 4×4.
    if n_samples <= 8:
        macro_rows, macro_cols = 4, 2
    elif n_samples <= 12:
        macro_rows, macro_cols = 6, 2
    elif n_samples <= 15:
        macro_rows, macro_cols = 5, 3
    else:
        macro_rows, macro_cols = 4, 4
    # Larger panels for readability (2026-04-22 V2b fix)
    fig_w = max(18.0, 2.8 * n_outer * macro_cols + 2.5)
    fig_h = 4.6 * macro_rows + 2.0
    fig = plt.figure(figsize=(fig_w, fig_h))

    outer_gs = gridspec.GridSpec(macro_rows, macro_cols, figure=fig,
                                 hspace=0.80, wspace=0.20,
                                 left=0.06, right=0.93, top=0.91, bottom=0.08)

    norm = Normalize(vmin=0.3, vmax=1.0)
    cmap = plt.get_cmap("RdYlGn")
    cmap.set_bad(color=PALETTE["grey"], alpha=0.4)
    last_im = None

    question = TRIO_QUESTIONS.get(
        (dim_row, dim_col, dim_outer),
        f"{DIM_SHORT_LABEL[dim_row]} × {DIM_SHORT_LABEL[dim_col]} × {DIM_SHORT_LABEL[dim_outer]} per-sample",
    )
    research_suptitle(
        fig, question,
        context=f"inner facet = {DIM_SHORT_LABEL[dim_outer]} slice · cell color = TP rate · grey = n<20 · 紅框 = FP<100",
        y=0.975, fontsize=18,
    )

    for idx, (sample, mode) in enumerate(SAMPLE_MODE_ORDER):
        macro_r = idx // macro_cols
        macro_c = idx % macro_cols
        sub = df5[(df5["sample"] == sample) & (df5["mode"] == mode)]
        if len(sub) == 0:
            continue

        inner_gs = gridspec.GridSpecFromSubplotSpec(
            1, n_outer, subplot_spec=outer_gs[macro_r, macro_c], wspace=0.08,
        )

        baseline = sub["tp_label"].mean()
        n_total = len(sub)
        agg = compute_tp_rate_by_key(sub, [dim_row, dim_col, dim_outer], min_n=20)

        meta = SAMPLE_FP_META.get((sample, mode), {})
        panel_title = (f"{sample_label(sample, mode)}  {meta.get('tier','')} "
                       f"FP={meta.get('fp','?'):,}  base={baseline:.2f}  n={n_total:,}")
        if meta.get("flag"):
            panel_title += f"  ({meta['flag']})"

        for j, z_val in enumerate(outer_order):
            ax = fig.add_subplot(inner_gs[0, j])
            slab = agg[agg[dim_outer] == z_val]
            pivot_rate = slab.pivot(index=dim_row, columns=dim_col, values="tp_rate") \
                .reindex(index=row_order, columns=col_order)
            pivot_n = slab.pivot(index=dim_row, columns=dim_col, values="n") \
                .reindex(index=row_order, columns=col_order).fillna(0).astype(int)

            data = pivot_rate.values.astype(float)
            im = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, norm=norm, aspect="auto")
            last_im = im

            ax.set_xticks(range(len(col_order)))
            ax.set_yticks(range(len(row_order)))
            ax.set_xticklabels([str(x) for x in col_order], rotation=45, ha="right", fontsize=8)
            if j == 0:
                ax.set_yticklabels([str(y) for y in row_order], fontsize=8)
                ax.set_ylabel(DIM_SHORT_LABEL[dim_row], fontsize=9)
            else:
                ax.set_yticklabels([])
            ax.set_title(f"{DIM_SHORT_LABEL[dim_outer]}={z_val}", fontsize=9,
                         color=PALETTE["dark"])

            for yi in range(len(row_order)):
                for xi in range(len(col_order)):
                    n = int(pivot_n.iat[yi, xi])
                    rate = pivot_rate.iat[yi, xi]
                    if n == 0:
                        continue
                    if np.isnan(rate):
                        ax.text(xi, yi, f"{n}", ha="center", va="center",
                                fontsize=7, color=PALETTE["dark"])
                    else:
                        colour = "black" if rate > 0.55 else "white"
                        ax.text(xi, yi, f"{rate*100:.0f}\n{n}",
                                ha="center", va="center", fontsize=7, color=colour)

            # FP<100 red frame
            if meta.get("fp", 0) < 100:
                for spine in ax.spines.values():
                    spine.set_edgecolor(PALETTE["red"])
                    spine.set_linewidth(1.2)

            # Outer panel title only on middle inner facet
            if j == n_outer // 2:
                ax.text(0.5, 1.35, panel_title, transform=ax.transAxes,
                        ha="center", va="bottom", fontsize=11, fontweight="bold",
                        color=PALETTE["dark"])

    # shared colorbar
    if last_im is not None:
        cbar_ax = fig.add_axes([0.945, 0.10, 0.008, 0.82])
        cbar = fig.colorbar(last_im, cax=cbar_ax)
        cbar.set_label("TP rate", fontsize=9)

    # takeaway caption at the bottom
    takeaway_text, colour = _trio_takeaway(dim_row, dim_col, dim_outer, df5)
    takeaway_caption(fig, takeaway_text, color=colour, y=0.005, fontsize=10)

    out = L3_DIR / f"L3_{DIM_SHORT_LABEL[dim_row].lower()}_{DIM_SHORT_LABEL[dim_col].lower()}_{DIM_SHORT_LABEL[dim_outer].lower()}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    apply_style()
    L3_DIR.mkdir(parents=True, exist_ok=True)
    df = load_master()
    df5 = build_5d_key(df)
    df5 = df5[df5["HPFineNGroups"].isin(DIM_ORDERS["HPFineNGroups"])]
    for trio in TRIOS:
        out = plot_trio(*trio, df5)
        print(f"[obs11] wrote {out}")


if __name__ == "__main__":
    main()
