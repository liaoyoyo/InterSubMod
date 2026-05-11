#!/usr/bin/env python3
"""
obs09 | L1 — Single-dimension per-sample marginal TP rate (V2 2026-04-22).

For each of the 5 dimensions (LOH / AF / CN / NG / NR), produce TWO figures:
  - `L1_{dim}_TO.png`     — 6 TO samples in a 2x3 grid (baseline 25-91%, strong gradient)
  - `L1_{dim}_paired.png` — 7 paired_full samples in a 3x3 grid (baseline 94-99%, saturated)

Each figure:
  - Research-question suptitle
  - Per-panel metadata badge (tier, FP, baseline, n)
  - Computed takeaway caption box at the bottom

Writes to ``figures/new/L1/``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _obs_common import (
    DIM_ORDERS, DIM_SHORT_LABEL, FIG_NEW_DIR, PALETTE,
    annotate_power_tier, apply_style, build_5d_key, compute_tp_rate_by_key,
    load_master, mode_display_name, mode_subset_grid, research_suptitle,
    sample_label, takeaway_caption,
)

L1_DIR = FIG_NEW_DIR / "L1"

# Research questions per (dim, mode). Filled at file-level for reuse.
QUESTIONS = {
    ("LOH_Subtype", "to_pileup"):
        "LOH_Subtype 對 TP rate 的邊際影響在 6 TO 樣本間是否方向一致？",
    ("LOH_Subtype", "paired_full"):
        "LOH_Subtype 在 paired_full 7 樣本的 TP rate 是否因 baseline 飽和而失去區分力？",
    ("AF_class", "to_pileup"):
        "AF_class (Extreme / Intermediate / Near-half) 在 6 TO 樣本的 TP rate 方向是否一致？",
    ("AF_class", "paired_full"):
        "AF_class 在 paired_full 樣本的 TP rate 有否可辨識差異？",
    ("cn_tier_F", "to_pileup"):
        "CN tier (T0-T4) 對 6 TO 樣本 TP rate 的影響是否跨樣本一致？",
    ("cn_tier_F", "paired_full"):
        "CN tier 在 paired_full 樣本間是否呈現一致 TP 梯度？",
    ("HPFineNGroups", "to_pileup"):
        "HPFineNGroups (NG) 梯度在 6 TO 樣本是否跨樣本穩健？",
    ("HPFineNGroups", "paired_full"):
        "HPFineNGroups 在 paired_full 是否因飽和失去區分？",
    ("nr_band", "to_pileup"):
        "NumReads band (low/mid/high) 對 6 TO 樣本 TP rate 的影響方向？",
    ("nr_band", "paired_full"):
        "NumReads band 在 paired_full 樣本的微弱梯度觀察",
}


def _compute_takeaway(dim: str, mode: str, df5: pd.DataFrame, order: list) -> tuple[str, str]:
    """Return (takeaway_text, colour_tag) for this (dim, mode) combination."""
    if mode == "to_pileup":
        samples = 6
    else:
        samples = 7

    # Compute direction: how many samples have last-category TP > first-category TP
    directions = {"up": 0, "down": 0, "flat": 0}
    gaps = []
    outliers = []
    for sample, m in [(s, mode) for s in df5[df5["mode"] == mode]["sample"].unique()]:
        sub = df5[(df5["sample"] == sample) & (df5["mode"] == m)]
        if len(sub) < 20:
            continue
        agg = compute_tp_rate_by_key(sub, [dim], min_n=20).set_index(dim).reindex(order)
        rates = agg["tp_rate"].dropna()
        if len(rates) < 2:
            continue
        first, last = rates.iloc[0], rates.iloc[-1]
        gap = last - first
        gaps.append((sample, gap))
        if gap > 0.05:
            directions["up"] += 1
        elif gap < -0.05:
            directions["down"] += 1
        else:
            directions["flat"] += 1

    up, down, flat = directions["up"], directions["down"], directions["flat"]
    total = up + down + flat
    if total == 0:
        return (f"n 不足：{mode_display_name(mode)} 下 {dim} 無可解讀訊號", "grey")

    dominant = max(directions, key=directions.get)
    gap_vals = [abs(g) for _, g in gaps]
    gap_range = f"{min(gap_vals):.2f}-{max(gap_vals):.2f}" if gap_vals else "n/a"

    if dominant == "up" and up >= total * 0.75:
        # Find outliers (wrong direction or flat)
        outl = [s for s, g in gaps if g < 0.05]
        colour = "green"
        msg = (f"✓ Takeaway：{up}/{total} 樣本呈現 {DIM_SHORT_LABEL[dim]} 正向梯度（gap {gap_range}）。"
               f"跨樣本一致性{'強' if up == total else '中-強'}")
        if outl:
            msg += f"；例外：{', '.join(outl)}"
        return (msg + "。", colour)
    elif dominant == "down" and down >= total * 0.75:
        colour = "blue"
        return (f"↓ Takeaway：{down}/{total} 樣本呈現 {DIM_SHORT_LABEL[dim]} 反向梯度（gap {gap_range}）。"
                f"高值類別反而 TP 較低，需檢查 confound。", colour)
    elif flat >= total * 0.75:
        colour = "grey"
        if mode == "paired_full":
            return (f"○ Takeaway：{flat}/{total} 樣本在 {DIM_SHORT_LABEL[dim]} 軸上差異 <0.05。"
                    f"paired_full baseline 已 ≥0.94，視覺上飽和無梯度。", colour)
        return (f"○ Takeaway：{flat}/{total} 樣本差異 <0.05；{DIM_SHORT_LABEL[dim]} 邊際效果不顯著。", colour)
    else:
        colour = "red"
        return (f"⚠ Takeaway：方向不一致（↑{up} / ↓{down} / ○{flat}）。"
                f"{DIM_SHORT_LABEL[dim]} 的影響存在樣本異質性，不可視為通用 filter。", colour)


def plot_dim_mode(dim: str, mode: str, df5: pd.DataFrame) -> Path:
    order = DIM_ORDERS[dim]
    fig, axes, sm_list = mode_subset_grid(mode)
    question = QUESTIONS.get((dim, mode), f"{DIM_SHORT_LABEL[dim]} × {mode_display_name(mode)}")
    research_suptitle(fig, question,
                      context="bars = TP rate · red dashed = sample baseline · error bars = Wilson 95% CI · 紅框 = FP<100 低信度",
                      y=0.975, fontsize=18)

    for ax, (sample, m) in zip(axes[:len(sm_list)], sm_list):
        sub = df5[(df5["sample"] == sample) & (df5["mode"] == m)]
        if len(sub) == 0:
            ax.axis("off")
            continue
        baseline = sub["tp_label"].mean()
        n_total = len(sub)
        agg = compute_tp_rate_by_key(sub, [dim], min_n=20)
        agg = agg.set_index(dim).reindex(order).reset_index()

        xs = np.arange(len(order))
        rates = agg["tp_rate"].values
        lo = agg["wilson_lo"].values
        hi = agg["wilson_hi"].values
        err_lower = np.where(np.isnan(rates), 0, np.clip(rates - lo, 0, None))
        err_upper = np.where(np.isnan(rates), 0, np.clip(hi - rates, 0, None))

        colors = []
        for r in rates:
            if np.isnan(r):
                colors.append(PALETTE["grey"])
            elif r > baseline + 0.10:
                colors.append(PALETTE["green"])
            elif r < baseline - 0.10:
                colors.append(PALETTE["red"])
            else:
                colors.append(PALETTE["blue"])

        ax.bar(xs, np.nan_to_num(rates), color=colors,
               edgecolor=PALETTE["dark"], linewidth=0.6)
        valid = ~np.isnan(rates)
        if valid.any():
            ax.errorbar(xs[valid], rates[valid],
                        yerr=[err_lower[valid], err_upper[valid]],
                        fmt="none", ecolor=PALETTE["dark"], capsize=3, linewidth=0.8)
        ax.axhline(baseline, color=PALETTE["accent"], linestyle="--", linewidth=1.3)

        ax.set_xticks(xs)
        ax.set_xticklabels([str(x) for x in order], rotation=20, ha="right", fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("TP rate", fontsize=9)
        ax.set_title(sample_label(sample, m), fontsize=11, fontweight="bold",
                     color=PALETTE["dark"])
        ax.grid(True, axis="y", alpha=0.25)

        for i, row in agg.iterrows():
            n = int(row["n"]) if pd.notna(row["n"]) else 0
            rate = row["tp_rate"]
            if np.isnan(rate):
                label = f"n={n}\n(n<20)" if n > 0 else "n=0"
                ax.text(i, 0.03, label, ha="center", fontsize=6.5, color=PALETTE["grey"])
            else:
                ax.text(i, min(rate + (err_upper[i] if valid[i] else 0) + 0.02, 1.01),
                        f"{rate:.2f}", ha="center", fontsize=7, color=PALETTE["dark"])
                ax.text(i, 0.03, f"n={n:,}", ha="center", fontsize=6.5, color=PALETTE["dark"])

        annotate_power_tier(ax, sample, m, baseline_tp=baseline, n_total=n_total)

    takeaway_text, colour = _compute_takeaway(dim, mode, df5, order)
    takeaway_caption(fig, takeaway_text, color=colour, y=0.018, fontsize=12)

    fig.tight_layout(rect=[0, 0.065, 1, 0.93])
    mode_tag = "TO" if mode == "to_pileup" else "paired"
    out = L1_DIR / f"L1_{DIM_SHORT_LABEL[dim].lower()}_{mode_tag}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    apply_style()
    L1_DIR.mkdir(parents=True, exist_ok=True)
    df = load_master()
    df5 = build_5d_key(df)
    df5 = df5[df5["HPFineNGroups"].isin(DIM_ORDERS["HPFineNGroups"])]

    for dim in ("LOH_Subtype", "AF_class", "cn_tier_F", "HPFineNGroups", "nr_band"):
        for mode in ("to_pileup", "paired_full"):
            out = plot_dim_mode(dim, mode, df5)
            print(f"[obs09] wrote {out}")


if __name__ == "__main__":
    main()
