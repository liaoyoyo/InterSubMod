#!/usr/bin/env python3
"""S06 — 三條路 codeflow figure.

Three vertical columns:
    Path 1 baseline (no flag) — gray
    Path 2 V5 only (somaticCalling reruns) — green
    Path 3 second round (no rerun, just re-phase) — red

Each column shows the code flow as stacked boxes with arrows.
Bottom: outcome row showing HP1:HP2 ratio.

Output: fig_s06_three_paths_codeflow.png (1600 x 900 px @ 200 dpi)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from matplotlib import font_manager

CJK_FONT = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
font_manager.fontManager.addfont(CJK_FONT)
plt.rcParams["font.family"] = ["Droid Sans Fallback", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).parent / "fig_s06_three_paths_codeflow.png"

# 統一 4 步：layout 對齊三路 outcome strip / ratio 位置
# 改用「肯定描述」+「形狀+顏色雙編碼」（path number 1/2/3 為形狀區分，避免紅綠 colorblind 風險）
paths = [
    {
        "number": "1",
        "title": "路 1 baseline",
        "color": "#7F8C8D",
        "steps": [
            "原始 phasing pass",
            "PON 不參與重 phase",
            "Pass 2 不啟動",
            "self-phasing 殘留",
        ],
        "outcome": "self-phasing 17.3:1 殘留",
        "ratio": "1.328",
        "ratio_color": "#7F8C8D",
        "verdict": "(基準)",
    },
    {
        "number": "2",
        "title": "路 2 V5-only (反轉)",
        "color": "#1976D2",  # blue 取代 green，避開紅綠並列
        "steps": [
            "V5 flag PON-only",
            "重算 somatic + 重分 hap",
            "(somaticCalling 重跑)",
            "self-phasing 反轉",
        ],
        "outcome": "舊 V5 純路 2",
        "ratio": "0.735",
        "ratio_color": "#1976D2",
        "verdict": "[OK] 反轉",
    },
    {
        "number": "3",
        "title": "路 3 second-round (抵消)",
        "color": "#C0392B",
        "steps": [
            "highPurity > 0.9 觸發",
            "只重分 hap (不重算)",
            "(無 somaticCalling)",
            "前一輪反轉被覆寫",
        ],
        "outcome": "self-phasing 偏移殘留",
        "ratio": "1.400",
        "ratio_color": "#C0392B",
        "verdict": "[X] 抵消",
    },
]


def draw_step_box(ax, x, y, w, h, color, text, *, bold=False):
    box = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.04",
        facecolor="#FAFAFA", edgecolor=color, linewidth=1.5,
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=11, fontweight="bold" if bold else "normal",
            color="#1F1F1F")


def main():
    fig, ax = plt.subplots(figsize=(16, 9.5), dpi=100)
    ax.set_xlim(-0.2, 3.2)
    ax.set_ylim(-3.5, 3.3)

    # Top title
    ax.text(1.5, 3.0,
            "三條路 audit：路 3 second round 抵消路 2 反轉效果",
            ha="center", va="center", fontsize=18, fontweight="bold",
            color="#C0392B")

    col_x = [0.5, 1.5, 2.5]
    col_w = 0.85
    step_h = 0.45

    for i, p in enumerate(paths):
        x = col_x[i]
        # Column title (large header bar with path number circle for shape encoding)
        head = mpatches.FancyBboxPatch(
            (x - col_w / 2, 2.15), col_w, 0.55,
            boxstyle="round,pad=0.04",
            facecolor=p["color"], edgecolor=p["color"], linewidth=2,
        )
        ax.add_patch(head)
        # Number circle (shape encoding for colorblind safety)
        num_circ = mpatches.Circle(
            (x - col_w / 2 + 0.15, 2.42), 0.13,
            facecolor="white", edgecolor=p["color"], linewidth=2.5, zorder=5,
        )
        ax.add_patch(num_circ)
        ax.text(x - col_w / 2 + 0.15, 2.42, p["number"], ha="center", va="center",
                fontsize=14, fontweight="bold", color=p["color"], zorder=6)
        ax.text(x + 0.05, 2.42, p["title"], ha="center", va="center",
                fontsize=14, fontweight="bold", color="white")

        # Steps stacked vertically
        n = len(p["steps"])
        step_y_start = 1.6
        step_gap = 0.65
        for j, step_text in enumerate(p["steps"]):
            y = step_y_start - j * step_gap
            draw_step_box(ax, x, y, col_w, step_h, p["color"], step_text)
            # Arrow to next
            if j < n - 1:
                arrow = FancyArrowPatch(
                    (x, y - step_h / 2), (x, y - step_gap + step_h / 2),
                    arrowstyle="-|>", mutation_scale=15,
                    color=p["color"], linewidth=1.5,
                )
                ax.add_patch(arrow)

        # Outcome strip — fixed y so all 3 columns align
        out_y = -1.30
        out_box = mpatches.FancyBboxPatch(
            (x - col_w / 2, out_y - 0.30), col_w, 0.55,
            boxstyle="round,pad=0.04",
            facecolor="#F5F5F5", edgecolor=p["color"], linewidth=2,
        )
        ax.add_patch(out_box)
        ax.text(x, out_y - 0.04, p["outcome"], ha="center", va="center",
                fontsize=10, color="#1F1F1F")

        # Ratio big number — fixed y, well below outcome
        ratio_y = -2.30
        ax.text(x, ratio_y + 0.45, "HP1:HP2",
                ha="center", va="center", fontsize=11, color="#7F8C8D")
        ax.text(x, ratio_y + 0.05, p["ratio"],
                ha="center", va="center", fontsize=30, fontweight="bold",
                color=p["ratio_color"])
        ax.text(x, ratio_y - 0.45, p["verdict"],
                ha="center", va="center", fontsize=13, fontweight="bold",
                color=p["ratio_color"])

    # Bottom verdict strip — pushed down below ratio numbers
    verd_box = mpatches.FancyBboxPatch(
        (-0.15, -3.4), 3.3, 0.55,
        boxstyle="round,pad=0.05",
        facecolor="#FFF3E0", edgecolor="#E67E22", linewidth=2,
    )
    ax.add_patch(verd_box)
    ax.text(1.5, -3.13,
            "NEW V5 (路 2+3) HP1:HP2=1.400 抵消 ｜ NEW V5 noPath3 (強制跳過路 3) HP1:HP2=1.127 復現舊 V5 反轉",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color="#1F1F1F")

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"[OK] {OUT}")


if __name__ == "__main__":
    main()
