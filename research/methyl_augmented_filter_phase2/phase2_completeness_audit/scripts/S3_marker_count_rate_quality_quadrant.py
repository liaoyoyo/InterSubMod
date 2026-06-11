#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S3_marker_count_rate_quality_quadrant.py

Build a single self-explanatory combined figure for slide S3 of
presenter_BCG_6slide.html. 3 panels:

  Panel A (~30%): NG / Marker / Rate definition diagram (cluster cartoon).
  Panel B (~40%): 4-quadrant scatter of (Marker count, Marker rate)
                  showing baseline -> V6 movement with quadrant tints.
  Panel C (~30%): Mechanism credibility (why it is trustworthy)
                  + counterexample (transparency).

Aspect: 16:9 (1920x1080 at dpi=150).
Read-only. No external TSV dependency: numbers are inlined from
20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01 (4-way audit 110/110 PASS).

Output PNG:
  research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/
      S3_marker_count_rate_quality_quadrant.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import font_manager, rcParams
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ---------------------------------------------------------------------------
# CJK font injection (Droid Sans Fallback) -- graceful fallback if absent.
# ---------------------------------------------------------------------------
_CJK_CANDIDATES = [
    "Droid Sans Fallback",
    "Noto Sans CJK TC",
    "Noto Sans CJK SC",
    "WenQuanYi Zen Hei",
    "DejaVu Sans",
]
_available = {f.name for f in font_manager.fontManager.ttflist}
_picked = next((c for c in _CJK_CANDIDATES if c in _available), "DejaVu Sans")
rcParams["font.family"] = [_picked, "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False

# ---------------------------------------------------------------------------
# Colour palette (aligned with presenter_BCG_6slide.html CSS tokens).
# ---------------------------------------------------------------------------
COL_BASELINE = "#DC2626"   # red (baseline)
COL_V6 = "#2563EB"         # blue (V6 fix)
COL_HP33 = "#9467BD"       # purple (HP=33 ambig)
COL_HP1 = "#F59E0B"        # amber (HP=1 cluster)
COL_HP2 = "#10B981"        # green (HP=2 cluster)
COL_TEXT = "#111827"
COL_MUTED = "#4B5563"
COL_FAINT = "#9CA3AF"
COL_GOOD_BG = "#ECFDF5"     # light green (good quadrant)
COL_BAD_BG = "#FEF2F2"      # light red  (bad quadrant)
COL_NEUTRAL_BG = "#F9FAFB"
COL_GOOD_BORDER = "#10B981"
COL_BAD_BORDER = "#EF4444"

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
OUT_PATH = OUT_DIR / "S3_marker_count_rate_quality_quadrant.png"
COPY_PATH = (
    Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/"
         "20260522_PI_V6_signoff_4goal/assets/figures/target2/"
         "S3_marker_count_rate_quality_quadrant.png")
)

# ---------------------------------------------------------------------------
# Hard-coded canonical data (20260519 4-way audit 110/110 PASS).
# ---------------------------------------------------------------------------
BASELINE_COUNT = 15738
BASELINE_RATE = 0.8967
V6_COUNT = 23980
V6_RATE = 0.9093
D_COUNT_PCT = (V6_COUNT - BASELINE_COUNT) / BASELINE_COUNT * 100.0  # +52.4 %
D_RATE_PP = (V6_RATE - BASELINE_RATE) * 100.0                       # +1.26 pp


# ===========================================================================
# Panel A — NG / Marker / Rate cluster cartoon + formula box
# ===========================================================================
def _draw_read_row(ax, y, n_reads, color, hp_label, *, bar_w=0.085, gap=0.012,
                   x0=0.10):
    """Render a horizontal row of square 'read' tiles to mimic a cluster."""
    for i in range(n_reads):
        x = x0 + i * (bar_w + gap)
        rect = FancyBboxPatch(
            (x, y), bar_w, 0.045,
            boxstyle="round,pad=0.002,rounding_size=0.005",
            linewidth=0.6, edgecolor=color, facecolor=color, alpha=0.85,
        )
        ax.add_patch(rect)
    # HP label on the right
    ax.text(
        x0 + n_reads * (bar_w + gap) + 0.02, y + 0.0225,
        hp_label, ha="left", va="center",
        fontsize=9.5, color=color, fontweight="bold",
    )


def _draw_cluster_box(ax, y0, h, color, label, label_color=None):
    """Bracket-style box around a cluster row."""
    rect = patches.Rectangle(
        (0.07, y0), 0.86, h,
        linewidth=1.0, edgecolor=color, facecolor="none",
        linestyle="--", alpha=0.55,
    )
    ax.add_patch(rect)
    ax.text(
        0.075, y0 + h - 0.012, label,
        ha="left", va="top",
        fontsize=8.5, color=label_color or color, style="italic",
    )


def panel_a_definition(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Title
    ax.text(
        0.5, 0.99, "Panel A · NG / Marker / Rate 定義",
        ha="center", va="top", fontsize=12.5, fontweight="bold", color=COL_TEXT,
    )
    ax.text(
        0.5, 0.940,
        "ISM region 內 read 在 phasing-space\n的聚類群數 (NGroups = NG)",
        ha="center", va="top", fontsize=9.5, color=COL_MUTED,
    )

    # Three clusters (HP=1, HP=2, HP=33)
    # Cluster 1: HP=1, 8 reads
    _draw_cluster_box(ax, 0.755, 0.085, COL_HP1, "group 1 · HP=1 cluster")
    _draw_read_row(ax, 0.770, 8, COL_HP1, "HP=1")

    # Cluster 2: HP=2, 4 reads
    _draw_cluster_box(ax, 0.635, 0.085, COL_HP2, "group 2 · HP=2 cluster")
    _draw_read_row(ax, 0.650, 4, COL_HP2, "HP=2")

    # Cluster 3: HP=33, 3 reads (V6 fix — ambig exposed)
    _draw_cluster_box(ax, 0.515, 0.085, COL_HP33,
                      "group 3 · HP=33 (V6 fix 暴露 ambig reads)")
    _draw_read_row(ax, 0.530, 3, COL_HP33, "HP=33")

    # NG = 3 callout under clusters
    ng_box = FancyBboxPatch(
        (0.10, 0.395), 0.80, 0.072,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        linewidth=1.2, edgecolor=COL_V6, facecolor="#EFF6FF",
    )
    ax.add_patch(ng_box)
    ax.text(
        0.5, 0.431,
        "NG = 3   (3 個 phasing clusters)",
        ha="center", va="center", fontsize=11,
        fontweight="bold", color=COL_V6,
    )

    # Formula box
    formula_box = FancyBboxPatch(
        (0.07, 0.075), 0.86, 0.290,
        boxstyle="round,pad=0.014,rounding_size=0.012",
        linewidth=1.0, edgecolor=COL_FAINT, facecolor="#FFFFFF",
    )
    ax.add_patch(formula_box)
    ax.text(
        0.5, 0.335, "公式定義",
        ha="center", va="center", fontsize=10.5,
        fontweight="bold", color=COL_TEXT,
    )
    ax.text(
        0.10, 0.270,
        "Marker count = Σ region with NG ≥ 3",
        ha="left", va="center", fontsize=10, color=COL_TEXT,
    )
    ax.text(
        0.10, 0.225,
        "      (region 觸發 somatic signature)",
        ha="left", va="center", fontsize=8.8, color=COL_MUTED, style="italic",
    )
    ax.text(
        0.10, 0.165,
        "Marker rate = n_TP / (n_TP + n_FP)",
        ha="left", va="center", fontsize=10, color=COL_TEXT,
    )
    ax.text(
        0.10, 0.120,
        "      within NG ≥ 3 region (TP purity)",
        ha="left", va="center", fontsize=8.8, color=COL_MUTED, style="italic",
    )


# ===========================================================================
# Panel B — 4-quadrant scatter (Marker count × Marker rate)
# ===========================================================================
def panel_b_quadrant(ax) -> None:
    # Axis ranges chosen so the V6 point and baseline both visible with margin.
    x_min, x_max = 12500.0, 27500.0
    y_min, y_max = 0.878, 0.922
    x_mid = (BASELINE_COUNT + V6_COUNT) / 2.0  # 19,859
    y_mid = (BASELINE_RATE + V6_RATE) / 2.0    # 0.9030

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # ---- Quadrant tints ----
    # Top-right (good) — V6 落點
    ax.add_patch(patches.Rectangle(
        (x_mid, y_mid), x_max - x_mid, y_max - y_mid,
        facecolor=COL_GOOD_BG, edgecolor="none", zorder=0,
    ))
    # Top-left (overly strict) — count down, rate up
    ax.add_patch(patches.Rectangle(
        (x_min, y_mid), x_mid - x_min, y_max - y_mid,
        facecolor=COL_NEUTRAL_BG, edgecolor="none", zorder=0,
    ))
    # Bottom-right (noise) — count up, rate down
    ax.add_patch(patches.Rectangle(
        (x_mid, y_min), x_max - x_mid, y_mid - y_min,
        facecolor=COL_BAD_BG, edgecolor="none", zorder=0,
    ))
    # Bottom-left (full failure)
    ax.add_patch(patches.Rectangle(
        (x_min, y_min), x_mid - x_min, y_mid - y_min,
        facecolor=COL_NEUTRAL_BG, edgecolor="none", zorder=0,
    ))

    # Quadrant dividers
    ax.axvline(x_mid, color=COL_FAINT, linestyle=":", linewidth=0.9, zorder=1)
    ax.axhline(y_mid, color=COL_FAINT, linestyle=":", linewidth=0.9, zorder=1)

    # ---- Quadrant labels (corners) — keep small to avoid clutter ----
    ax.text(
        x_min + 250, y_max - 0.0018,
        "左上 · 量↓ 質↑  ✗ 過嚴", ha="left", va="top",
        fontsize=8.5, color=COL_MUTED,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor=COL_FAINT, linewidth=0.6, alpha=0.9),
        zorder=3,
    )
    ax.text(
        x_max - 250, y_min + 0.0018,
        "右下 · 量↑ 質↓  ✗ 雜訊化", ha="right", va="bottom",
        fontsize=8.5, color=COL_BAD_BORDER,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor=COL_BAD_BORDER, linewidth=0.7, alpha=0.95),
        zorder=3,
    )
    ax.text(
        x_min + 250, y_min + 0.0018,
        "左下 · 量↓ 質↓  ✗ 全敗", ha="left", va="bottom",
        fontsize=8.5, color=COL_MUTED,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor=COL_FAINT, linewidth=0.6, alpha=0.9),
        zorder=3,
    )
    # Top-right corner highlighted (good quadrant where V6 lives)
    ax.text(
        x_max - 250, y_max - 0.0018,
        "右上 · 量↑ 質↑  ✓ 量質同進", ha="right", va="top",
        fontsize=9, color=COL_GOOD_BORDER, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white",
                  edgecolor=COL_GOOD_BORDER, linewidth=0.9, alpha=0.97),
        zorder=3,
    )

    # ---- baseline and V6 points ----
    ax.scatter(
        [BASELINE_COUNT], [BASELINE_RATE],
        s=260, marker="o", color=COL_BASELINE, edgecolor="white",
        linewidth=2.0, zorder=5, label="baseline",
    )
    ax.scatter(
        [V6_COUNT], [V6_RATE],
        s=520, marker="*", color=COL_V6, edgecolor="white",
        linewidth=2.0, zorder=5, label="V6",
    )

    # Arrow baseline -> V6
    arrow = FancyArrowPatch(
        (BASELINE_COUNT + 280, BASELINE_RATE + 0.0006),
        (V6_COUNT - 380, V6_RATE - 0.0006),
        arrowstyle="-|>", mutation_scale=22,
        color=COL_V6, linewidth=2.2, zorder=4,
        connectionstyle="arc3,rad=-0.08",
    )
    ax.add_patch(arrow)

    # Delta label near arrow midpoint
    mid_x = (BASELINE_COUNT + V6_COUNT) / 2.0
    mid_y = (BASELINE_RATE + V6_RATE) / 2.0
    ax.text(
        mid_x, mid_y + 0.0042,
        f"+{D_COUNT_PCT:.1f}%  count\n+{D_RATE_PP:.2f}pp  rate",
        ha="center", va="bottom",
        fontsize=10.5, color=COL_V6, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor=COL_V6, linewidth=1.0, alpha=0.97),
        zorder=6,
    )

    # baseline point annotation
    ax.annotate(
        f"baseline\n({BASELINE_COUNT:,}, {BASELINE_RATE:.4f})",
        xy=(BASELINE_COUNT, BASELINE_RATE),
        xytext=(BASELINE_COUNT - 100, BASELINE_RATE - 0.0095),
        fontsize=9.5, color=COL_BASELINE, fontweight="bold",
        ha="center", va="top",
        arrowprops=dict(arrowstyle="-", color=COL_BASELINE,
                        linewidth=0.8, alpha=0.7),
        zorder=6,
    )

    # V6 point annotation (star) — place below-right so it doesn't collide
    # with the +52.4% label above.
    ax.annotate(
        f"★ V6  ({V6_COUNT:,}, {V6_RATE:.4f})\nnet TP:FP = 14:1",
        xy=(V6_COUNT, V6_RATE),
        xytext=(V6_COUNT + 50, V6_RATE - 0.0070),
        fontsize=9.5, color=COL_V6, fontweight="bold",
        ha="center", va="top",
        arrowprops=dict(arrowstyle="-", color=COL_V6,
                        linewidth=0.8, alpha=0.7),
        zorder=6,
    )

    # "量質同進 ✓" badge — place above the V6 star
    ax.text(
        V6_COUNT, V6_RATE + 0.0085,
        "量質同進 ✓", ha="center", va="bottom",
        fontsize=10, color=COL_GOOD_BORDER, fontweight="bold", zorder=7,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#D1FAE5",
                  edgecolor=COL_GOOD_BORDER, linewidth=1.0, alpha=0.97),
    )

    # ---- Axis cosmetics ----
    ax.set_xlabel("Marker count  (ISM region 標為 NG≥3 的數量)",
                  fontsize=10.5, color=COL_TEXT)
    ax.set_ylabel("Marker rate  (NG≥3 region 中 TP 占比 / TP purity)",
                  fontsize=10.5, color=COL_TEXT)
    ax.set_title("Panel B · Marker count × Marker rate 4-象限  (HCC1395 4-way ISM)",
                 fontsize=12.5, fontweight="bold", color=COL_TEXT, pad=8)

    # Grid
    ax.grid(True, linestyle=":", linewidth=0.5, color=COL_FAINT, alpha=0.6,
            zorder=0)
    # X tick formatting (use thousands separator)
    ax.set_xticks([13000, 15000, 17000, 19000, 21000, 23000, 25000, 27000])
    ax.set_xticklabels([f"{int(t / 1000)}K" for t in ax.get_xticks()],
                       fontsize=9.5)
    ax.set_yticks([0.88, 0.89, 0.90, 0.91, 0.92])
    ax.tick_params(axis="y", labelsize=9.5)

    # Frame
    for spine in ax.spines.values():
        spine.set_color(COL_MUTED)
        spine.set_linewidth(0.8)


# ===========================================================================
# Panel C — credibility + counterexample
# ===========================================================================
def panel_c_credibility(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.5, 0.99, "Panel C · 機制可信度 + Counterexample",
        ha="center", va="top", fontsize=12.5, fontweight="bold", color=COL_TEXT,
    )

    # ---- Credibility box (top half) ----
    cred_box = FancyBboxPatch(
        (0.04, 0.50), 0.92, 0.42,
        boxstyle="round,pad=0.014,rounding_size=0.014",
        linewidth=1.2, edgecolor=COL_GOOD_BORDER, facecolor=COL_GOOD_BG,
    )
    ax.add_patch(cred_box)
    ax.text(
        0.07, 0.885, "為何可信 (機制鏈)",
        ha="left", va="top",
        fontsize=11, fontweight="bold", color=COL_GOOD_BORDER,
    )

    cred_lines = [
        "✓ baseline priority bug → ambig reads",
        "    被強制標 HP=1 (read-level 偏移)",
        "✓ V6 fix → HP=33 暴露 ambig reads",
        "    (Layer 1.5 ambig firing)",
        "✓ cluster 多樣化 → region NG≥3 ↑",
        "    → Marker count 上升 +52.4%",
        "✓ 新進 NG≥3 region 多為真實 subclone",
        "    → Marker rate 同步上升 +1.26pp",
    ]
    y_top = 0.835
    line_step = 0.040
    for i, line in enumerate(cred_lines):
        ax.text(
            0.075, y_top - i * line_step, line,
            ha="left", va="top",
            fontsize=9.2, color=COL_TEXT,
        )

    # ---- Counterexample box (bottom half) ----
    counter_box = FancyBboxPatch(
        (0.04, 0.045), 0.92, 0.405,
        boxstyle="round,pad=0.014,rounding_size=0.014",
        linewidth=1.2, edgecolor=COL_BAD_BORDER, facecolor=COL_BAD_BG,
    )
    ax.add_patch(counter_box)
    ax.text(
        0.07, 0.415, "Counterexample (transparency)",
        ha="left", va="top",
        fontsize=11, fontweight="bold", color=COL_BAD_BORDER,
    )

    counter_lines = [
        "⚠ chr12  ΔF1 = −0.0158",
        "⚠ AF<0.1 stratum  ΔF1 = −0.0304",
        "⚠ 350 lost TP 全落在 LOH 邊界",
        "—",
        "→ V6 是 net 14:1 trade-off",
        "    (非 universally better)",
    ]
    y_top = 0.360
    line_step = 0.048
    for i, line in enumerate(counter_lines):
        weight = "bold" if line.startswith("→") else "normal"
        color = COL_BAD_BORDER if line.startswith("→") else COL_TEXT
        ax.text(
            0.075, y_top - i * line_step, line,
            ha="left", va="top",
            fontsize=10 if line.startswith("→") else 9.2,
            color=color, fontweight=weight,
        )


# ===========================================================================
# Top-level figure assembly
# ===========================================================================
def build_figure() -> plt.Figure:
    fig = plt.figure(figsize=(12.8, 7.2), dpi=150, facecolor="white")
    # 1920 x 1080 at dpi=150 -> 12.8 x 7.2 inches.

    # Title strip
    fig.text(
        0.5, 0.972,
        "S3 · V6 ISM Region 層  —  Marker count × rate 量質 2-axis  +  機制可信度  +  counterexample",
        ha="center", va="top",
        fontsize=14.5, fontweight="bold", color=COL_TEXT,
    )
    fig.text(
        0.5, 0.935,
        "HCC1395 4-way ISM · baseline → V6 · +52.4% count · +1.26pp rate · net TP:FP = 14:1 · 110/110 PASS",
        ha="center", va="top",
        fontsize=10, color=COL_MUTED,
    )

    # Source footer
    fig.text(
        0.5, 0.012,
        "source: master TSV 112 cols × 35,332 regions (HCC1395 4-way ISM) · "
        "20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01 audit 110/110 PASS",
        ha="center", va="bottom",
        fontsize=8.5, color=COL_MUTED, style="italic",
    )

    # Layout: 3 columns (30 / 40 / 30) within central band.
    gs = gridspec.GridSpec(
        1, 3,
        figure=fig,
        left=0.030, right=0.990,
        bottom=0.072, top=0.855,
        wspace=0.22,
        width_ratios=[0.30, 0.40, 0.30],
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    panel_a_definition(ax_a)
    panel_b_quadrant(ax_b)
    panel_c_credibility(ax_c)

    return fig


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build_figure()
    fig.savefig(OUT_PATH, dpi=150, bbox_inches=None, facecolor="white")
    print(f"[S3] wrote {OUT_PATH}  ({OUT_PATH.stat().st_size / 1024:.1f} KiB)")

    # Copy to presentation assets dir
    COPY_PATH.parent.mkdir(parents=True, exist_ok=True)
    COPY_PATH.write_bytes(OUT_PATH.read_bytes())
    print(f"[S3] copied  {COPY_PATH}")

    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
