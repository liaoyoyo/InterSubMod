#!/usr/bin/env python3
"""S05 — ploidy bug → Pass 2 從未觸發 因果鏈 figure.

5-stage horizontal causal flow:
    4/12 BAM → ploidy=0 → purity=0 → highPurity=false → Pass 1 only

Output: fig_s05_ploidy_causal_chain.png (1600 x 700 px @ 200 dpi)
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

OUT = Path(__file__).parent / "fig_s05_ploidy_causal_chain.png"

stages = [
    {"label": "4/12 BAM", "body": "pononly_v2b\ntumor_tagged.bam", "color": "#1F4E79"},
    {"label": "ploidy=0", "body": "ploidy bug\n計算短路", "color": "#C0392B"},
    {"label": "purity=0", "body": "三次多項式\n→ 0", "color": "#C0392B"},
    {"label": "highPurity\n=false", "body": "Pass 2 second\nround 跳過", "color": "#C0392B"},
    {"label": "Pass 1 only\nresult", "body": "= 4/29 PI 報告\n所有 V5 數值", "color": "#E67E22"},
]


def draw_box(ax, x, y, w, h, color, label, body):
    box = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.05",
        facecolor="white", edgecolor=color, linewidth=2.5,
    )
    ax.add_patch(box)
    # Label (top portion, large bold)
    ax.text(x, y + h / 4, label, ha="center", va="center",
            fontsize=14, fontweight="bold", color=color)
    # Body (bottom portion, smaller)
    ax.text(x, y - h / 4, body, ha="center", va="center",
            fontsize=10, color="#1F1F1F")


def main():
    fig, ax = plt.subplots(figsize=(16, 6.5), dpi=100)
    ax.set_xlim(-0.5, 5)
    ax.set_ylim(-2.0, 2.5)

    # Top title
    ax.text(2.25, 2.0,
            "ploidy bug → purity=0 → highPurity=false → Pass 2 從未觸發",
            ha="center", va="center", fontsize=18, fontweight="bold",
            color="#C0392B")

    # Boxes
    box_w = 0.82
    box_h = 1.0
    y_box = 0.7
    for i, s in enumerate(stages):
        x = i * 0.95 + 0.05
        draw_box(ax, x, y_box, box_w, box_h, s["color"], s["label"], s["body"])
        # Arrow to next
        if i < len(stages) - 1:
            arrow = FancyArrowPatch(
                (x + box_w / 2, y_box), (x + 0.95 - box_w / 2, y_box),
                arrowstyle="-|>", mutation_scale=20,
                color="#7F8C8D", linewidth=2.0,
            )
            ax.add_patch(arrow)

    # Solution box (4/30 d0bcd8c)
    sol_box = mpatches.FancyBboxPatch(
        (-0.4, -1.3), 4.9, 0.8,
        boxstyle="round,pad=0.05",
        facecolor="#E8F5E9", edgecolor="#27AE60", linewidth=2.0,
    )
    ax.add_patch(sol_box)
    ax.text(2.05, -0.95, "[已解決] 4/30 d0bcd8c — ploidy 計算分支修補",
            ha="center", va="center", fontsize=13, fontweight="bold",
            color="#1B5E20")
    ax.text(2.05, -1.15,
            "Pass 2 真實觸發 → 但 4/29 PI 報告 BAM 仍是 4/12 舊版 → 數值不變 → 必須重跑",
            ha="center", va="center", fontsize=11, color="#1F1F1F")

    # P0 caveat strip
    p0_box = mpatches.FancyBboxPatch(
        (-0.4, -1.95), 4.9, 0.5,
        boxstyle="round,pad=0.05",
        facecolor="#FFEBEE", edgecolor="#C0392B", linewidth=2.0,
    )
    ax.add_patch(p0_box)
    ax.text(2.05, -1.7,
            "[!] P0 — Pass 2 重驗（~21 hr 平行 × 7 樣本 + ~4 hr 分析 = 25 hr）",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color="#B71C1C")

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"[OK] {OUT}")


if __name__ == "__main__":
    main()
