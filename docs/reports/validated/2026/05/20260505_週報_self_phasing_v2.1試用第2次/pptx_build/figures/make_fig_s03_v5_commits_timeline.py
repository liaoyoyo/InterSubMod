#!/usr/bin/env python3
"""S03 — V5 5 commits timeline figure.

Renders a horizontal timeline:
    early phase (3 commits) → 4/30 critical phase (2 commits)
    last 2 commits highlighted in red (ploidy fix + threshold cherry-pick)

Output: fig_s03_v5_5commits_timeline.png (1600 x 700 px @ 200 dpi)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import font_manager

# CJK font registration
CJK_FONT = "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"
font_manager.fontManager.addfont(CJK_FONT)
plt.rcParams["font.family"] = ["Droid Sans Fallback", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).parent / "fig_s03_v5_5commits_timeline.png"

# Commits — 簡化早期 3 commit 為「Phase 0/1/2」標籤，4/30 critical 保留 SHA
commits = [
    {"sha": "Phase 0", "label": "PON-only flag", "date": "early", "kind": "early"},
    {"sha": "Phase 1", "label": "tag layer fix", "date": "early", "kind": "early"},
    {"sha": "Phase 2", "label": "somatic fallback", "date": "early", "kind": "early"},
    {"sha": "d0bcd8c", "label": "ploidy bug fix\n→ Pass 2 觸發", "date": "2026-04-30", "kind": "critical"},
    {"sha": "938f0df", "label": "threshold 0.95→0.9\ncherry-pick", "date": "2026-04-30", "kind": "critical"},
]

EARLY = "#1F4E79"
CRITICAL = "#C0392B"
LIGHT_BG = "#FAFAFA"


def main():
    fig, ax = plt.subplots(figsize=(16, 6.5), dpi=100)
    ax.set_xlim(-0.5, len(commits) - 0.5)
    ax.set_ylim(-1.0, 2.5)

    # Background
    fig.patch.set_facecolor("white")

    # Timeline horizontal axis
    ax.axhline(y=0.4, xmin=0.05, xmax=0.95, color="#BDC3C7", linewidth=2.5, zorder=1)

    # Phase shading
    early_rect = mpatches.FancyBboxPatch(
        (-0.4, -0.7), 3.0, 2.8, boxstyle="round,pad=0.05",
        facecolor="#ECF0F1", edgecolor=EARLY, linewidth=1.5, alpha=0.5, zorder=0,
    )
    ax.add_patch(early_rect)
    ax.text(1.0, 1.95, "Phase 0-2 早期 V5 修補（3 commits）",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color=EARLY, zorder=5)

    critical_rect = mpatches.FancyBboxPatch(
        (2.6, -0.7), 2.0, 2.8, boxstyle="round,pad=0.05",
        facecolor="#FFEBEE", edgecolor=CRITICAL, linewidth=1.5, alpha=0.6, zorder=0,
    )
    ax.add_patch(critical_rect)
    ax.text(3.5, 1.95, "* 4/30 Critical 修補（2 commits）",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color=CRITICAL, zorder=5)

    # Nodes
    for i, c in enumerate(commits):
        color = CRITICAL if c["kind"] == "critical" else EARLY
        # Circle node on the timeline
        circ = mpatches.Circle((i, 0.4), 0.18, facecolor=color, edgecolor="white",
                                linewidth=2.5, zorder=4)
        ax.add_patch(circ)
        # SHA above the node (large)
        ax.text(i, 1.1, c["sha"], ha="center", va="center",
                fontsize=15, fontweight="bold", color=color, zorder=5,
                family="monospace")
        # Date below SHA (small)
        ax.text(i, 0.85, c["date"], ha="center", va="center",
                fontsize=9, color="#7F8C8D", zorder=5)
        # Label below the node
        ax.text(i, -0.2, c["label"], ha="center", va="top",
                fontsize=11, color="#1F1F1F", zorder=5)

    # Verdict strip below
    ax.text(2.0, -0.75,
            "[已解決] git diff --stat HEAD = 0 byte（working tree clean） · 舊 caveat R1 已隨 938f0df 解決",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color="#27AE60",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#E8F5E9",
                      edgecolor="#27AE60", linewidth=1.5))

    # Title
    ax.set_title("V5 修補鏈 = 5 commits（含 4/30 ploidy fix + threshold cherry-pick）",
                 fontsize=18, fontweight="bold", color="#1F4E79", pad=15)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"[OK] {OUT}")


if __name__ == "__main__":
    main()
