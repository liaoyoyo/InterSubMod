#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S2_priority_bug_mechanism_plus_ALTonly.py

Build a single combined figure for slide S2 of presenter_BCG_6slide.html:
- Panel A (35%): priority bug mechanism diagram (baseline vs V6)
- Panel B (65%): grouped bar of ALT-only HP1:HP2 ratio across baseline + 5 V6 samples

Aspect: 16:9 (1920x1080 at dpi=150)
Read-only on A2 TSV; outputs PNG only.

Data source:
  research/methyl_augmented_filter_phase2/phase2_completeness_audit/A2_ALT_only_HP_ratio_5sample.tsv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib import font_manager, rcParams

# ---------------------------------------------------------------------------
# CJK font injection (Droid Sans Fallback) — graceful fallback if absent.
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
# Colour palette (matches presenter_BCG_6slide.html CSS tokens).
# ---------------------------------------------------------------------------
COL_BASELINE = "#DC2626"   # red (priority bug)
COL_V6 = "#2563EB"         # blue (fix)
COL_AMBIG = "#9467BD"      # purple (HP=33)
COL_REF = "#6B7280"        # gray dashed (1:1 reference)
COL_TEXT = "#111827"
COL_MUTED = "#4B5563"
COL_BG_BASELINE = "#FEE2E2"
COL_BG_V6 = "#DBEAFE"
COL_BG_AMBIG = "#EDE9FE"

ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = ROOT / "A2_ALT_only_HP_ratio_5sample.tsv"
OUT_DIR = ROOT / "figures"
OUT_PATH = OUT_DIR / "S2_priority_bug_mechanism_plus_ALTonly.png"
COPY_PATH = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/"
    "20260522_PI_V6_signoff_4goal/assets/figures/target1/"
    "S2_priority_bug_mechanism_plus_ALTonly.png"
)


def load_ratios() -> list[tuple[str, str, float]]:
    """Return list of (sample, bam_version, ALT_only_ratio) from A2 TSV."""
    rows: list[tuple[str, str, float]] = []
    with TSV_PATH.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            try:
                ratio = float(r["ALT_only_ratio"])
            except (KeyError, ValueError):
                continue
            rows.append((r["sample"], r["bam_version"], ratio))
    return rows


# ---------------------------------------------------------------------------
# Panel A — mechanism diagram (matplotlib patches).
# ---------------------------------------------------------------------------

def _box(ax, xy, w, h, text, fc, ec, fontsize=9, fontweight="normal", textcolor=None):
    box = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.05",
        linewidth=1.2,
        facecolor=fc,
        edgecolor=ec,
    )
    ax.add_patch(box)
    cx = xy[0] + w / 2
    cy = xy[1] + h / 2
    ax.text(
        cx,
        cy,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=fontweight,
        color=textcolor or COL_TEXT,
        wrap=True,
    )
    return (cx, cy)


def _arrow(ax, start, end, color=COL_MUTED, label=None, lw=1.3, label_offset=(0, 0)):
    arr = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=lw,
        color=color,
    )
    ax.add_patch(arr)
    if label:
        mx = (start[0] + end[0]) / 2 + label_offset[0]
        my = (start[1] + end[1]) / 2 + label_offset[1]
        ax.text(mx, my, label, ha="center", va="center", fontsize=8, color=COL_MUTED, style="italic")


def draw_mechanism(ax):
    """Draw the priority-bug mechanism diagram on a 0-10 / 0-100 logical canvas.

    Layout map (y-axis, 0=bottom):
        100  panel title
        92   baseline section header
        78-86 src box A  (baseline)
        70-77 votes box  (arrow above)
        58-67 HP=1 bug box
        51   separator
        47   V6 section header
        33-41 src box B
        20-30 HP=33 ambig box (longer arrow above)
        6-13  contrast banner
    """
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 100)
    ax.set_axis_off()

    # Panel title
    ax.text(
        5.0,
        97,
        "Panel A · Priority bug 機制",
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
        color=COL_TEXT,
    )

    # ------ Upper half: baseline ------
    ax.text(0.3, 91, "baseline (bug active)", fontsize=10, fontweight="bold", color=COL_BASELINE)

    # Source box
    _box(ax, (0.4, 80), 9.2, 7.5, "ALT read in somatic-only region", "#F3F4F6", "#9CA3AF", fontsize=9)
    # arrow down + label to the side
    _arrow(ax, (5.0, 80), (5.0, 75), color=COL_MUTED, lw=1.4)
    ax.text(5.0, 77.5, "無 germline het 可投票", ha="center", va="center", fontsize=8.5, color=COL_MUTED, style="italic",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none"))
    # votes box
    _box(ax, (0.4, 67), 9.2, 7.5, "votes = [ ]  (空票)", "#F3F4F6", "#9CA3AF", fontsize=9)
    _arrow(ax, (5.0, 67), (5.0, 62), color=COL_MUTED, lw=1.4)
    ax.text(5.0, 64.5, "hardcode default → HP=1", ha="center", va="center", fontsize=8.5, color=COL_MUTED, style="italic",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none"))
    # bug box
    _box(
        ax,
        (0.4, 53),
        9.2,
        8.5,
        "HP = 1  (bug)\n4.41 HP1-leaning fingerprint",
        COL_BG_BASELINE,
        COL_BASELINE,
        fontsize=9.5,
        fontweight="bold",
        textcolor=COL_BASELINE,
    )

    # Separator
    ax.plot([0.3, 9.7], [49, 49], color="#D1D5DB", linewidth=0.8, linestyle=(0, (4, 3)))

    # ------ Lower half: V6 ------
    ax.text(0.3, 46, "V6 fix (Layer 1.5 ambig)", fontsize=10, fontweight="bold", color=COL_V6)

    _box(ax, (0.4, 33), 9.2, 7.5, "ALT read in somatic-only region", "#F3F4F6", "#9CA3AF", fontsize=9)
    _arrow(ax, (5.0, 33), (5.0, 25), color=COL_MUTED, lw=1.4)
    ax.text(5.0, 28.5, "V6 偵測 germline-absent\n→ 改判 ambig", ha="center", va="center",
            fontsize=8.5, color=COL_MUTED, style="italic",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none"))
    _box(
        ax,
        (0.4, 15),
        9.2,
        9.0,
        "HP = 33  (ambig)\n0.43 HP2-leaning balanced",
        COL_BG_AMBIG,
        COL_AMBIG,
        fontsize=9.5,
        fontweight="bold",
        textcolor=COL_AMBIG,
    )

    # Contrast annotation
    ax.text(
        5.0,
        7,
        "baseline 4.41 → V6 0.43  ·  ~10× HP2 reverse",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color=COL_TEXT,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFF7ED", edgecolor="#F59E0B", linewidth=1.0),
    )


# ---------------------------------------------------------------------------
# Panel B — grouped bar (baseline vs 5 V6 samples).
# ---------------------------------------------------------------------------

def draw_ratio_bar(ax, rows):
    # Build canonical x order: HCC1395_baseline, then 5 V6.
    order = [
        ("HCC1395", "baseline"),
        ("HCC1395", "V6"),
        ("HCC1937", "V6"),
        ("HCC1954", "V6"),
        ("H1437", "V6"),
        ("H2009", "V6"),
    ]
    lut = {(s, v): r for (s, v, r) in rows}
    ratios = [lut.get(k, float("nan")) for k in order]
    labels = [
        "HCC1395\nbaseline",
        "HCC1395\nV6",
        "HCC1937\nV6",
        "HCC1954\nV6",
        "H1437\nV6",
        "H2009\nV6",
    ]
    colors = [COL_BASELINE] + [COL_V6] * 5

    x = list(range(len(order)))
    bars = ax.bar(
        x,
        ratios,
        color=colors,
        edgecolor="#1F2937",
        linewidth=0.7,
        width=0.62,
        zorder=3,
    )

    # 1:1 reference
    ax.axhline(y=1.0, color=COL_REF, linestyle="--", linewidth=1.1, zorder=2)
    ax.text(
        len(order) - 0.5,
        1.05,
        "1:1 balanced reference",
        ha="right",
        va="bottom",
        fontsize=9,
        color=COL_REF,
        style="italic",
    )

    # Value labels on each bar
    for i, (rect, v) in enumerate(zip(bars, ratios)):
        if v != v:  # nan
            continue
        ax.text(
            rect.get_x() + rect.get_width() / 2,
            v + 0.08,
            f"{v:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
            color=COL_BASELINE if i == 0 else COL_V6,
        )

    # Set y-lim early so annotations sit inside frame
    y_top = max(5.6, max(ratios) + 1.2)
    ax.set_ylim(0, y_top)

    # Highlight + annotation: baseline priority bug fingerprint (label to the RIGHT of bar)
    ax.annotate(
        "priority bug\nfingerprint",
        xy=(0.32, ratios[0] - 0.3),
        xytext=(0.95, ratios[0] - 0.6),
        fontsize=9.5,
        color=COL_BASELINE,
        fontweight="bold",
        ha="left",
        va="center",
        arrowprops=dict(
            arrowstyle="->",
            color=COL_BASELINE,
            lw=1.2,
            connectionstyle="arc3,rad=0.22",
        ),
    )

    # Sweep annotation: 10x HP2 reverse (baseline -> HCC1395 V6), positioned above baseline bar
    y_arc = ratios[0] + 0.8
    ax.annotate(
        "",
        xy=(1, ratios[1] + 0.25),
        xytext=(0, ratios[0] + 0.25),
        arrowprops=dict(
            arrowstyle="->",
            color="#F59E0B",
            lw=1.8,
            connectionstyle="arc3,rad=-0.42",
        ),
    )
    ax.text(
        0.5,
        y_arc,
        "~10× HP2 reverse",
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
        color="#B45309",
    )

    # Axes
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("ALT-only HP1 : HP2 ratio", fontsize=11)
    ax.tick_params(axis="y", labelsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.45, zorder=1)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#9CA3AF")
    ax.spines["bottom"].set_color("#9CA3AF")

    # Title for panel
    ax.set_title(
        "Panel B · ALT-only HP1:HP2 ratio 跨 5 V6 樣本 fingerprint reverse",
        fontsize=12,
        fontweight="bold",
        loc="left",
        pad=8,
        color=COL_TEXT,
    )

    # Footer note inside axes (top-right, well inside frame)
    ax.text(
        len(order) - 1.0,
        y_top * 0.93,
        "5/5 V6 樣本 ratio 全 < 1.0 ✓",
        ha="right",
        va="top",
        fontsize=10,
        fontweight="bold",
        color=COL_V6,
        bbox=dict(boxstyle="round,pad=0.3", facecolor=COL_BG_V6, edgecolor=COL_V6, linewidth=1.0),
    )


# ---------------------------------------------------------------------------
# Main composition
# ---------------------------------------------------------------------------

def main() -> int:
    if not TSV_PATH.exists():
        print(f"[FATAL] TSV not found: {TSV_PATH}", file=sys.stderr)
        return 2
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    COPY_PATH.parent.mkdir(parents=True, exist_ok=True)

    rows = load_ratios()
    if not rows:
        print("[FATAL] No data rows parsed from TSV", file=sys.stderr)
        return 2

    fig = plt.figure(figsize=(12.8, 7.2), dpi=150)  # 1920 x 1080
    gs = gridspec.GridSpec(
        nrows=1,
        ncols=2,
        width_ratios=[0.35, 0.65],
        left=0.045,
        right=0.985,
        top=0.80,
        bottom=0.10,
        wspace=0.16,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    draw_mechanism(ax_a)
    draw_ratio_bar(ax_b, rows)

    fig.text(
        0.045,
        0.95,
        "ALT-only HP1:HP2 ratio — priority bug indicator",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
        color=COL_TEXT,
    )
    fig.text(
        0.045,
        0.895,
        "baseline 4.41 (bug active) → V6 0.43–0.79 全 < 1.0  ·  跨 5 樣本 fingerprint reverse 5/5",
        ha="left",
        va="top",
        fontsize=12,
        color=COL_MUTED,
    )

    # Provenance footer
    fig.text(
        0.045,
        0.018,
        "source: A2_ALT_only_HP_ratio_5sample.tsv  ·  scope: chr8+chr19 (H2009 chr19 only)  ·  V6 = germline-absent revert",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#6B7280",
    )

    fig.savefig(OUT_PATH, dpi=150, facecolor="white")
    fig.savefig(COPY_PATH, dpi=150, facecolor="white")
    plt.close(fig)

    print(f"[OK] wrote {OUT_PATH}")
    print(f"[OK] copied to {COPY_PATH}")
    print(f"[INFO] font used: {_picked}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
