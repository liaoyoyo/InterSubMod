"""
gen_v4_evolution_figures.py — HKU v4 演化樹 figure regenerator (5/23 user revisions)

Produces 3 v2 PNGs into ../figures/ (replacing the original F2/F3/F4 v1 conceptually,
but the v1 files are PRESERVED as historical reference):

    F2v2_6state_with_tree_inset.png        — 3x3 R/A/X matrix + 4 inset evolution trees (A/B/C/D)
    F3v2_x_imputation_flow.png             — 4-panel X imputation flow with 3-step icons
    F4v2_evolution_tree_4scenarios.png     — 2x2 grid of A/B/C/D evolution trees

Design contract (per user 5/23 evening):
  - Linear depth naming: HP1-1-1 = child-of HP1-1 (linear evolution, +1 mutation per step)
  - HP1-1, HP1-2 are first-generation children of HP1 (same depth, different mutation event)
  - NO sibling numbering convention; depth notation only

Scenario definitions:
  A: {RR, AR, AA}           Linear: RR(HP1) -> AR(HP1-1) -> AA(HP1-1-1)
  B: {RR, AR, RA}           Branching: HP1 -> {HP1-1 (AR), HP1-2 (RA)}
  C: {RR, AR, RA, AA}       Branching + linear; AA ambiguous: HP1-1-1 or HP1-2-1
  D: HP1 + HP2 simultaneous Two parallel lineages with full A/C structure each

Design constraints:
  - matplotlib + seaborn (re-use _common.setup_cjk_font + same style as gen_concept_figures.py)
  - CJK font: Droid Sans Fallback (per memory feedback_matplotlib_cjk_font_rule.md)
  - DPI=150 PNG; 16:9 aspect for PPT/HTML embed
  - Independent legends/titles/source notes per figure
  - ASCII fallback for box-drawing chars (no U+2500 — per A3 prior matplotlib bug)

Repro:
    cd /big7_disk/liaoyoyo2001/InterSubMod
    python research/hku_collaboration/scripts/gen_v4_evolution_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Project root for plot_setup
PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.lib.plot_setup import setup_plot_style  # noqa: E402

import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle  # noqa: E402
import numpy as np  # noqa: E402
import seaborn as sns  # noqa: E402

# CJK font chain (verified: DejaVu Sans + Droid Sans Fallback covers TC/SC/JP)
setup_plot_style(base_size=11, dpi=150)
_cjk_chain = ["DejaVu Sans", "Droid Sans Fallback"]
matplotlib.rcParams["font.family"] = _cjk_chain
matplotlib.rcParams["font.sans-serif"] = _cjk_chain

FIG_DIR = PROJECT_ROOT / "research" / "hku_collaboration" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# --- Color palette (aligned with gen_concept_figures.py for visual continuity) ---
COLOR_REF = "#4DABF7"        # blue (Ref allele)
COLOR_ALT = "#FF6B6B"        # red (Alt allele)
COLOR_MISSING = "#CED4DA"    # gray (X missing)

# HP lineage colors (warm = HP1 lineage; cool = HP2 lineage)
COLOR_HP1_ROOT = "#FFE3E3"      # very light red (HP1 root)
COLOR_HP1_GEN1 = "#FFA8A8"      # light red (HP1-N first gen)
COLOR_HP1_GEN2 = "#FF6B6B"      # red (HP1-N-M second gen)
COLOR_HP2_ROOT = "#D0EBFF"      # very light blue
COLOR_HP2_GEN1 = "#74C0FC"      # light blue
COLOR_HP2_GEN2 = "#4DABF7"      # blue
COLOR_AMBIG = "#9775FA"         # purple (ambiguous AA)
COLOR_EDGE = "#343A40"

SOURCE_NOTE = "Source: HKU collaboration v4 evolution figures · gen_v4_evolution_figures.py"


# ============================================================
# Helper: draw an evolution tree node (rounded rect with label)
# ============================================================
def draw_tree_node(ax, x, y, label, sublabel="", facecolor=COLOR_HP1_GEN1,
                   edgecolor=COLOR_EDGE, width=1.2, height=0.55,
                   fontsize=9, sub_fontsize=7, text_color="#212529"):
    """Draw a rounded-rectangle tree node with HP code label + optional sublabel."""
    box = FancyBboxPatch(
        (x - width / 2, y - height / 2), width, height,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        linewidth=1.3, edgecolor=edgecolor, facecolor=facecolor, alpha=0.95,
    )
    ax.add_patch(box)
    if sublabel:
        ax.text(x, y + 0.10, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color=text_color)
        ax.text(x, y - 0.15, sublabel, ha="center", va="center",
                fontsize=sub_fontsize, color="#495057", style="italic")
    else:
        ax.text(x, y, label, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color=text_color)


def draw_tree_edge(ax, x1, y1, x2, y2, label="", color="#868E96", linewidth=1.4,
                   fontsize=7, ambiguous=False):
    """Draw a parent->child edge with optional mutation event label."""
    style = "->" if not ambiguous else "->"
    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle=style, mutation_scale=11,
        color=color, linewidth=linewidth,
        linestyle="--" if ambiguous else "-",
    )
    ax.add_patch(arrow)
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mx + 0.05, my, label, ha="left", va="center",
                fontsize=fontsize, color="#495057", style="italic",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="#DEE2E6", linewidth=0.5, alpha=0.85))


# ============================================================
# Scenario tree drawing primitives (reused across F2 inset + F4)
# ============================================================
def draw_scenario_A(ax, x0=0, y0=0, scale=1.0, compact=False):
    """A: RR(HP1) -> AR(HP1-1) -> AA(HP1-1-1) — linear evolution."""
    fs = 9 if not compact else 7
    sfs = 7 if not compact else 5.5
    w = 1.3 * scale
    h = 0.55 * scale
    # Three nodes vertically
    nodes = [
        (x0 + 1.5 * scale, y0 + 2.6 * scale, "HP1", "RR (germline)", COLOR_HP1_ROOT),
        (x0 + 1.5 * scale, y0 + 1.5 * scale, "HP1-1", "AR (1st somatic)", COLOR_HP1_GEN1),
        (x0 + 1.5 * scale, y0 + 0.4 * scale, "HP1-1-1", "AA (2nd somatic)", COLOR_HP1_GEN2),
    ]
    for (x, y, lab, sub, color) in nodes:
        draw_tree_node(ax, x, y, lab, sub, facecolor=color,
                       width=w, height=h, fontsize=fs, sub_fontsize=sfs)
    # Edges
    draw_tree_edge(ax, nodes[0][0], nodes[0][1] - h / 2,
                   nodes[1][0], nodes[1][1] + h / 2,
                   label="+AR site" if not compact else "",
                   fontsize=sfs)
    draw_tree_edge(ax, nodes[1][0], nodes[1][1] - h / 2,
                   nodes[2][0], nodes[2][1] + h / 2,
                   label="+RA site" if not compact else "",
                   fontsize=sfs)


def draw_scenario_B(ax, x0=0, y0=0, scale=1.0, compact=False):
    """B: RR(HP1) branches to HP1-1 (AR) and HP1-2 (RA) — branching 1st gen."""
    fs = 9 if not compact else 7
    sfs = 7 if not compact else 5.5
    w = 1.3 * scale
    h = 0.55 * scale
    # Root + 2 children
    root = (x0 + 1.5 * scale, y0 + 2.6 * scale, "HP1", "RR (germline)", COLOR_HP1_ROOT)
    left = (x0 + 0.7 * scale, y0 + 1.2 * scale, "HP1-1", "AR", COLOR_HP1_GEN1)
    right = (x0 + 2.3 * scale, y0 + 1.2 * scale, "HP1-2", "RA", COLOR_HP1_GEN1)
    for (x, y, lab, sub, color) in [root, left, right]:
        draw_tree_node(ax, x, y, lab, sub, facecolor=color,
                       width=w, height=h, fontsize=fs, sub_fontsize=sfs)
    draw_tree_edge(ax, root[0], root[1] - h / 2,
                   left[0], left[1] + h / 2,
                   label="+AR site" if not compact else "", fontsize=sfs)
    draw_tree_edge(ax, root[0], root[1] - h / 2,
                   right[0], right[1] + h / 2,
                   label="+RA site" if not compact else "", fontsize=sfs)


def draw_scenario_C(ax, x0=0, y0=0, scale=1.0, compact=False):
    """C: RR -> {AR, RA} -> {AA, AA} — branching + linear, AA ambiguous."""
    fs = 9 if not compact else 7
    sfs = 7 if not compact else 5.5
    w = 1.3 * scale
    h = 0.55 * scale
    root = (x0 + 1.5 * scale, y0 + 3.0 * scale, "HP1", "RR", COLOR_HP1_ROOT)
    left = (x0 + 0.6 * scale, y0 + 1.8 * scale, "HP1-1", "AR", COLOR_HP1_GEN1)
    right = (x0 + 2.4 * scale, y0 + 1.8 * scale, "HP1-2", "RA", COLOR_HP1_GEN1)
    aa_l = (x0 + 0.6 * scale, y0 + 0.5 * scale, "HP1-1-1", "AA", COLOR_HP1_GEN2)
    aa_r = (x0 + 2.4 * scale, y0 + 0.5 * scale, "HP1-2-1", "AA", COLOR_HP1_GEN2)
    for (x, y, lab, sub, color) in [root, left, right, aa_l, aa_r]:
        draw_tree_node(ax, x, y, lab, sub, facecolor=color,
                       width=w, height=h, fontsize=fs, sub_fontsize=sfs)
    draw_tree_edge(ax, root[0], root[1] - h / 2, left[0], left[1] + h / 2,
                   label="+AR" if not compact else "", fontsize=sfs)
    draw_tree_edge(ax, root[0], root[1] - h / 2, right[0], right[1] + h / 2,
                   label="+RA" if not compact else "", fontsize=sfs)
    draw_tree_edge(ax, left[0], left[1] - h / 2, aa_l[0], aa_l[1] + h / 2,
                   label="+RA" if not compact else "", fontsize=sfs)
    draw_tree_edge(ax, right[0], right[1] - h / 2, aa_r[0], aa_r[1] + h / 2,
                   label="+AR" if not compact else "", fontsize=sfs)
    # Ambiguity bracket
    if not compact:
        ax.annotate(
            "", xy=(aa_r[0] - w / 2 - 0.08, aa_r[1]),
            xytext=(aa_l[0] + w / 2 + 0.08, aa_l[1]),
            arrowprops=dict(arrowstyle="<->", color=COLOR_AMBIG, lw=1.5, linestyle=":"),
        )
        ax.text((aa_l[0] + aa_r[0]) / 2, aa_l[1] - 0.45 * scale,
                "AA ambiguous\n(HP1-1-1 or HP1-2-1)",
                ha="center", va="top", fontsize=sfs, color=COLOR_AMBIG,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#F3F0FF",
                          edgecolor=COLOR_AMBIG, linewidth=1.0))


def draw_scenario_D(ax, x0=0, y0=0, scale=1.0, compact=False):
    """D: HP1 + HP2 parallel lineages, each full A/C structure."""
    fs = 8 if not compact else 6.5
    sfs = 6 if not compact else 5.0
    w = 1.1 * scale
    h = 0.48 * scale

    # ------ HP1 lineage (left half, warm colors) ------
    hp1_root = (x0 + 1.2 * scale, y0 + 3.5 * scale, "HP1", "RR", COLOR_HP1_ROOT)
    hp1_l = (x0 + 0.45 * scale, y0 + 2.4 * scale, "HP1-1", "AR", COLOR_HP1_GEN1)
    hp1_r = (x0 + 1.95 * scale, y0 + 2.4 * scale, "HP1-2", "RA", COLOR_HP1_GEN1)
    hp1_aa_l = (x0 + 0.45 * scale, y0 + 1.2 * scale, "HP1-1-1", "AA", COLOR_HP1_GEN2)
    hp1_aa_r = (x0 + 1.95 * scale, y0 + 1.2 * scale, "HP1-2-1", "AA", COLOR_HP1_GEN2)

    # ------ HP2 lineage (right half, cool colors) ------
    hp2_root = (x0 + 4.0 * scale, y0 + 3.5 * scale, "HP2", "RR", COLOR_HP2_ROOT)
    hp2_l = (x0 + 3.25 * scale, y0 + 2.4 * scale, "HP2-1", "AR", COLOR_HP2_GEN1)
    hp2_r = (x0 + 4.75 * scale, y0 + 2.4 * scale, "HP2-2", "RA", COLOR_HP2_GEN1)
    hp2_aa_l = (x0 + 3.25 * scale, y0 + 1.2 * scale, "HP2-1-1", "AA", COLOR_HP2_GEN2)
    hp2_aa_r = (x0 + 4.75 * scale, y0 + 1.2 * scale, "HP2-2-1", "AA", COLOR_HP2_GEN2)

    all_nodes = [hp1_root, hp1_l, hp1_r, hp1_aa_l, hp1_aa_r,
                 hp2_root, hp2_l, hp2_r, hp2_aa_l, hp2_aa_r]
    for (x, y, lab, sub, color) in all_nodes:
        draw_tree_node(ax, x, y, lab, sub, facecolor=color,
                       width=w, height=h, fontsize=fs, sub_fontsize=sfs)

    # Edges HP1
    for parent, child, lab in [
        (hp1_root, hp1_l, "+AR"), (hp1_root, hp1_r, "+RA"),
        (hp1_l, hp1_aa_l, "+RA"), (hp1_r, hp1_aa_r, "+AR"),
    ]:
        draw_tree_edge(ax, parent[0], parent[1] - h / 2,
                       child[0], child[1] + h / 2,
                       label=lab if not compact else "",
                       fontsize=sfs, color="#FA5252")
    # Edges HP2
    for parent, child, lab in [
        (hp2_root, hp2_l, "+AR"), (hp2_root, hp2_r, "+RA"),
        (hp2_l, hp2_aa_l, "+RA"), (hp2_r, hp2_aa_r, "+AR"),
    ]:
        draw_tree_edge(ax, parent[0], parent[1] - h / 2,
                       child[0], child[1] + h / 2,
                       label=lab if not compact else "",
                       fontsize=sfs, color="#339AF0")

    # Germline het root annotation (parent of both HP1/HP2 roots)
    if not compact:
        germ_x = x0 + 2.6 * scale
        germ_y = y0 + 4.6 * scale
        draw_tree_node(ax, germ_x, germ_y, "germline het",
                       "Ref / Alt (normal diploid)",
                       facecolor="#FFF9DB", width=2.0 * scale, height=0.55 * scale,
                       fontsize=fs, sub_fontsize=sfs)
        # Edges to HP1 root + HP2 root
        for child in [hp1_root, hp2_root]:
            draw_tree_edge(ax, germ_x, germ_y - 0.55 * scale / 2,
                           child[0], child[1] + h / 2,
                           label="", fontsize=sfs, color="#868E96")

    # Separator between two lineages
    if not compact:
        ax.axvline(x0 + 2.6 * scale, y0 + 0.6 * scale, y0 + 3.4 * scale,
                   color="#DEE2E6", linewidth=0.8, linestyle=":")


# ============================================================
# F2v2: 3x3 R/A/X matrix + 4 inset evolution trees
# ============================================================
def fig_F2v2():
    fig = plt.figure(figsize=(16, 9), dpi=150)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.1], height_ratios=[1.0, 1.0],
                          wspace=0.18, hspace=0.30,
                          left=0.05, right=0.97, top=0.90, bottom=0.06)

    # Left column: 3x3 matrix spanning both rows
    ax_mat = fig.add_subplot(gs[:, 0])

    # Build 3x3 intensity (Alt signal)
    states = ["X", "Ref", "Alt"]
    alt_signal = np.array([
        [0.25, 0.0, 0.5],
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 1.0],
    ])
    labels_state = np.array([
        ["XX", "XR", "XA"],
        ["RX", "RR", "RA"],
        ["AX", "AR", "AA"],
    ])
    # v2 renames: linear depth notation
    labels_rule = np.array([
        ["need imputation\n(X,X)", "need imputation\n(X,R)", "need imputation\n(X,A)"],
        ["need imputation\n(R,X)", "HP1\n(germline root)", "HP1-1\n(1st somatic)"],
        ["need imputation\n(A,X)", "HP1-1\n(1st somatic)", "HP1-1-1\n(2nd somatic)"],
    ])

    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    sns.heatmap(
        alt_signal, annot=False, cmap=cmap, vmin=0, vmax=1,
        cbar_kws={"label": "Alt signal strength (Ref→Alt)", "shrink": 0.65, "pad": 0.02},
        linewidths=1.5, linecolor="white",
        xticklabels=[f"2nd site: {s}" for s in states],
        yticklabels=[f"1st site: {s}" for s in states],
        ax=ax_mat, square=True,
    )
    for i in range(3):
        for j in range(3):
            ax_mat.text(j + 0.5, i + 0.28, labels_state[i, j],
                        ha="center", va="center",
                        fontsize=15, fontweight="bold",
                        color="white" if alt_signal[i, j] > 0.6 or alt_signal[i, j] < 0.05 else "#212529")
            ax_mat.text(j + 0.5, i + 0.70, labels_rule[i, j],
                        ha="center", va="center",
                        fontsize=7.5,
                        color="white" if alt_signal[i, j] > 0.6 else "#495057",
                        style="italic")

    ax_mat.set_title("3 x 3 R/A/X matrix (diagonal = homotype; off-diag = mixed)",
                     fontsize=11, color="#495057", pad=10)
    ax_mat.set_xlabel("")
    ax_mat.set_ylabel("")
    ax_mat.tick_params(axis="both", labelsize=10)

    # Right side: 4 inset trees (A, B, C, D)
    ax_A = fig.add_subplot(gs[0, 1])
    ax_B = fig.add_subplot(gs[1, 1])

    # Split right column further into 2x2 manually using ax_A/ax_B and child positions
    # Simpler: use sub-gridspec for 2x2 inset grid
    inset_gs = gs[:, 1].subgridspec(2, 2, wspace=0.18, hspace=0.30)
    ax_A.remove()
    ax_B.remove()
    ax_sA = fig.add_subplot(inset_gs[0, 0])
    ax_sB = fig.add_subplot(inset_gs[0, 1])
    ax_sC = fig.add_subplot(inset_gs[1, 0])
    ax_sD = fig.add_subplot(inset_gs[1, 1])

    for ax, label, drawer, xlim, ylim in [
        (ax_sA, "A: linear {RR, AR, AA}", draw_scenario_A, (0, 3.0), (0, 3.5)),
        (ax_sB, "B: branching {RR, AR, RA}", draw_scenario_B, (0, 3.0), (0, 3.5)),
        (ax_sC, "C: full {RR, AR, RA, AA}", draw_scenario_C, (0, 3.0), (0, 3.8)),
        (ax_sD, "D: HP1 + HP2 parallel", draw_scenario_D, (0, 5.5), (0, 5.5)),
    ]:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("auto")
        ax.axis("off")
        scale = 1.0 if drawer is not draw_scenario_D else 1.0
        drawer(ax, x0=0, y0=0, scale=scale, compact=True)
        ax.set_title(label, fontsize=9.5, fontweight="bold",
                     color="#212529", pad=4)

    # Title + subtitle
    fig.suptitle("F2v2 · Two-locus R/A/X 6-state matrix + 4 evolution scenarios (inset)",
                 fontsize=15, fontweight="bold", y=0.98)
    fig.text(0.5, 0.93,
             "v2 linear depth naming: HP1-1-1 = child-of HP1-1 (no sibling numbering)",
             ha="center", fontsize=10, color="#495057", style="italic")
    fig.text(0.99, 0.005, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
             ha="right", va="bottom")

    out = FIG_DIR / "F2v2_6state_with_tree_inset.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F2v2] {out}")
    return out


# ============================================================
# F3v2: X imputation flow with 3-step icons
# ============================================================
def fig_F3v2():
    fig, axes = plt.subplots(4, 1, figsize=(16, 11), dpi=150)
    fig.subplots_adjust(hspace=0.55, top=0.93, bottom=0.06, left=0.04, right=0.97)

    panels = [
        # (start, end_options, X_position, before_colors, after_top, after_bot)
        ("RX", ["RR", "RA"], 1, [COLOR_REF, COLOR_MISSING],
         [COLOR_REF, COLOR_REF], [COLOR_REF, COLOR_ALT]),
        ("AX", ["AA", "AR"], 1, [COLOR_ALT, COLOR_MISSING],
         [COLOR_ALT, COLOR_ALT], [COLOR_ALT, COLOR_REF]),
        ("XR", ["RR", "AR"], 0, [COLOR_MISSING, COLOR_REF],
         [COLOR_REF, COLOR_REF], [COLOR_ALT, COLOR_REF]),
        ("XA", ["RA", "AA"], 0, [COLOR_MISSING, COLOR_ALT],
         [COLOR_REF, COLOR_ALT], [COLOR_ALT, COLOR_ALT]),
    ]

    for ax, (start, end_opts, x_pos, before, after_top, after_bot) in zip(axes, panels):
        ax.set_xlim(0, 16)
        ax.set_ylim(0, 4)
        ax.axis("off")

        # Panel title bar
        ax.add_patch(Rectangle((0, 3.5), 16, 0.45, facecolor="#E7F5FF",
                               edgecolor="#74C0FC", linewidth=1.0))
        ax.text(0.3, 3.72, f"{start}  ->  {{ {end_opts[0]}, {end_opts[1]} }}",
                fontsize=14, fontweight="bold", color="#1864AB",
                va="center", ha="left")
        ax.text(15.7, 3.72,
                f"X at {'1st' if x_pos == 0 else '2nd'} site -> impute via 3-step procedure",
                fontsize=9.5, color="#495057", style="italic",
                va="center", ha="right")

        # --- LEFT: original state ---
        ax.text(1.3, 2.9, "Original (X)", fontsize=10,
                fontweight="bold", color="#495057", ha="center")
        for i, c in enumerate(before):
            circ = Circle((0.8 + i * 1.0, 1.9), 0.25,
                          color=c, ec=COLOR_EDGE, linewidth=1.5)
            ax.add_patch(circ)
            lbl = "X" if c == COLOR_MISSING else ("R" if c == COLOR_REF else "A")
            ax.text(0.8 + i * 1.0, 1.9, lbl, ha="center", va="center",
                    fontsize=10, fontweight="bold", color="white")
        ax.plot([0.3, 2.3], [1.9, 1.9], color="#495057", linewidth=2, zorder=0)
        ax.text(1.3, 1.0, "longphase\ncannot resolve",
                ha="center", va="top", fontsize=7.5, color="#868E96")

        # --- MIDDLE: 3-step icons ---
        # Arrow into middle
        ax.add_patch(FancyArrowPatch((2.7, 1.9), (3.7, 1.9),
                                     arrowstyle="->", mutation_scale=14,
                                     color="#868E96", linewidth=1.5))

        step_x = [4.5, 7.5, 10.5]
        step_titles = [
            "(1) same-PS reads",
            "(2) methyl. cluster",
            "(3) flanking SNV vote",
        ]
        step_desc = [
            "Compare reads sharing\nthe same phaseset (PS)",
            "Hierarchical cluster on\nread x CpG matrix",
            "Vote from nearby SNV\nallele pair (R/A)",
        ]

        for k, (sx, st, sd) in enumerate(zip(step_x, step_titles, step_desc)):
            # Step box
            ax.add_patch(FancyBboxPatch(
                (sx - 1.2, 0.4), 2.4, 2.6,
                boxstyle="round,pad=0.04,rounding_size=0.10",
                linewidth=1.2, edgecolor="#74C0FC", facecolor="#F8F9FA",
                alpha=0.95,
            ))
            ax.text(sx, 2.75, st, ha="center", va="center",
                    fontsize=9, fontweight="bold", color="#1864AB")

            # Mini icon per step
            if k == 0:
                # Step 1: 3 short read segments stacked w/ same color bar
                for r in range(3):
                    ry = 1.95 - r * 0.25
                    ax.plot([sx - 0.7, sx + 0.7], [ry, ry],
                            color="#5C7CFA", linewidth=2.4, zorder=0)
                    # PS tag dot at left
                    ax.add_patch(Circle((sx - 0.85, ry), 0.06,
                                        color="#5C7CFA", ec="none"))
                ax.text(sx, 1.05, "PS=12345", ha="center", fontsize=7,
                        color="#5C7CFA", fontweight="bold")
            elif k == 1:
                # Step 2: mini methyl heatmap
                rng = np.random.default_rng(7)
                heat = rng.random((3, 5))
                for i in range(3):
                    for j in range(5):
                        ax.add_patch(Rectangle(
                            (sx - 0.6 + j * 0.22, 2.05 - i * 0.22),
                            0.20, 0.20,
                            facecolor=plt.cm.RdYlBu_r(heat[i, j]),
                            edgecolor="#DEE2E6", linewidth=0.4,
                        ))
                ax.text(sx, 1.05, "read x CpG", ha="center",
                        fontsize=7, color="#495057")
            else:
                # Step 3: two flanking SNV pairs (R-?-R and A-?-A)
                # Top read RXR
                for i, c in enumerate([COLOR_REF, COLOR_MISSING, COLOR_REF]):
                    ax.add_patch(Circle((sx - 0.5 + i * 0.45, 2.05),
                                        0.13, color=c, ec=COLOR_EDGE, lw=0.8))
                ax.plot([sx - 0.65, sx + 0.45], [2.05, 2.05],
                        color="#495057", lw=1.2, zorder=0)
                # Bottom read AXA
                for i, c in enumerate([COLOR_ALT, COLOR_MISSING, COLOR_ALT]):
                    ax.add_patch(Circle((sx - 0.5 + i * 0.45, 1.55),
                                        0.13, color=c, ec=COLOR_EDGE, lw=0.8))
                ax.plot([sx - 0.65, sx + 0.45], [1.55, 1.55],
                        color="#495057", lw=1.2, zorder=0)
                ax.text(sx, 1.05, "flanking pair vote",
                        ha="center", fontsize=7, color="#495057")

            # Description below icon
            ax.text(sx, 0.55, sd, ha="center", va="center",
                    fontsize=7.5, color="#495057")

            # Arrow chaining steps
            if k < 2:
                ax.add_patch(FancyArrowPatch(
                    (step_x[k] + 1.3, 1.9), (step_x[k + 1] - 1.3, 1.9),
                    arrowstyle="->", mutation_scale=12,
                    color="#74C0FC", linewidth=1.2,
                ))

        # Arrow out of last step into result column
        ax.add_patch(FancyArrowPatch((11.8, 1.9), (12.8, 1.9),
                                     arrowstyle="->", mutation_scale=14,
                                     color="#868E96", linewidth=1.5))

        # --- RIGHT: post-imputation options (top + bottom) ---
        ax.text(14.4, 2.9, "Imputed outcomes", fontsize=10,
                fontweight="bold", color="#495057", ha="center")
        # Top option
        for i, c in enumerate(after_top):
            circ = Circle((13.6 + i * 0.7, 2.35), 0.22,
                          color=c, ec=COLOR_EDGE, linewidth=1.3)
            ax.add_patch(circ)
            lbl = "R" if c == COLOR_REF else "A"
            ax.text(13.6 + i * 0.7, 2.35, lbl, ha="center", va="center",
                    fontsize=9, fontweight="bold", color="white")
        ax.plot([13.2, 14.5], [2.35, 2.35], color="#495057", lw=1.8, zorder=0)
        ax.text(15.3, 2.35, f"-> {end_opts[0]}", ha="left", va="center",
                fontsize=10, color="#1864AB", fontweight="bold")
        # Bottom option
        for i, c in enumerate(after_bot):
            circ = Circle((13.6 + i * 0.7, 1.45), 0.22,
                          color=c, ec=COLOR_EDGE, linewidth=1.3)
            ax.add_patch(circ)
            lbl = "R" if c == COLOR_REF else "A"
            ax.text(13.6 + i * 0.7, 1.45, lbl, ha="center", va="center",
                    fontsize=9, fontweight="bold", color="white")
        ax.plot([13.2, 14.5], [1.45, 1.45], color="#495057", lw=1.8, zorder=0)
        ax.text(15.3, 1.45, f"-> {end_opts[1]}", ha="left", va="center",
                fontsize=10, color="#C92A2A", fontweight="bold")

    # Legend
    legend_elems = [
        mpatches.Patch(facecolor=COLOR_REF, edgecolor=COLOR_EDGE, label="Ref allele"),
        mpatches.Patch(facecolor=COLOR_ALT, edgecolor=COLOR_EDGE, label="Alt allele"),
        mpatches.Patch(facecolor=COLOR_MISSING, edgecolor=COLOR_EDGE,
                       label="X (missing / longphase unresolved)"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", fontsize=10,
               frameon=True, ncol=3, bbox_to_anchor=(0.5, -0.005))

    fig.suptitle("F3v2 · X imputation: 3-step ISM procedure across 4 X patterns",
                 fontsize=15, fontweight="bold", y=0.98)
    fig.text(0.5, 0.955,
             "(1) same-PS reads -> (2) methylation cluster -> (3) flanking SNV pair vote",
             ha="center", fontsize=10, color="#495057", style="italic")
    fig.text(0.99, 0.005, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
             ha="right", va="bottom")

    out = FIG_DIR / "F3v2_x_imputation_flow.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F3v2] {out}")
    return out


# ============================================================
# F4v2: 2x2 grid of A/B/C/D evolution trees
# ============================================================
def fig_F4v2():
    fig = plt.figure(figsize=(16, 10), dpi=150)
    gs = fig.add_gridspec(2, 2, wspace=0.10, hspace=0.30,
                          left=0.04, right=0.97, top=0.91, bottom=0.06)

    panel_specs = [
        ("A", "Scenario A: linear evolution {RR, AR, AA}",
         "HP1 -> HP1-1 -> HP1-1-1 (each step = +1 somatic mutation)",
         draw_scenario_A, 1.3, (0, 4), (0, 4)),
        ("B", "Scenario B: branching at 1st generation {RR, AR, RA}",
         "HP1 -> {HP1-1 (AR), HP1-2 (RA)} (two distinct 1st-gen mutation events)",
         draw_scenario_B, 1.3, (0, 4), (0, 4)),
        ("C", "Scenario C: branching + linear {RR, AR, RA, AA}",
         "AA reachable from either HP1-1-1 or HP1-2-1 (ambiguous - resolve via methyl. cluster or flanking SNV)",
         draw_scenario_C, 1.3, (0, 4), (0, 4.5)),
        ("D", "Scenario D: HP1 + HP2 parallel lineages (germline het root)",
         "Two haplotypes each carry full A/C-style somatic substructure",
         draw_scenario_D, 1.4, (0, 8), (0, 6)),
    ]

    for idx, (tag, title, subtitle, drawer, scale, xlim, ylim) in enumerate(panel_specs):
        r, c = divmod(idx, 2)
        ax = fig.add_subplot(gs[r, c])
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("auto")
        ax.axis("off")

        # Panel title bar
        ax.add_patch(Rectangle((xlim[0], ylim[1] - 0.45), xlim[1] - xlim[0], 0.45,
                               facecolor="#FFF9DB", edgecolor="#FFD43B",
                               linewidth=1.0, alpha=0.7))
        ax.text(xlim[0] + 0.15, ylim[1] - 0.22, f"({tag}) {title}",
                fontsize=11, fontweight="bold", color="#5C3D00",
                va="center", ha="left")

        # Tree
        drawer(ax, x0=0, y0=0, scale=scale, compact=False)

        # Subtitle under panel
        ax.text((xlim[0] + xlim[1]) / 2, 0.05, subtitle,
                ha="center", va="bottom", fontsize=8.5, color="#495057",
                style="italic")

    # Legend (HP1 lineage warm, HP2 lineage cool, ambiguous purple, germline root yellow)
    legend_elems = [
        mpatches.Patch(facecolor=COLOR_HP1_ROOT, edgecolor=COLOR_EDGE, label="HP1 root (germline)"),
        mpatches.Patch(facecolor=COLOR_HP1_GEN1, edgecolor=COLOR_EDGE, label="HP1 1st-gen (+1 somatic)"),
        mpatches.Patch(facecolor=COLOR_HP1_GEN2, edgecolor=COLOR_EDGE, label="HP1 2nd-gen (+2 somatic)"),
        mpatches.Patch(facecolor=COLOR_HP2_ROOT, edgecolor=COLOR_EDGE, label="HP2 root"),
        mpatches.Patch(facecolor=COLOR_HP2_GEN1, edgecolor=COLOR_EDGE, label="HP2 1st-gen"),
        mpatches.Patch(facecolor=COLOR_HP2_GEN2, edgecolor=COLOR_EDGE, label="HP2 2nd-gen"),
        mpatches.Patch(facecolor="#FFF9DB", edgecolor="#FFD43B", label="germline het root"),
        mpatches.Patch(facecolor="#F3F0FF", edgecolor=COLOR_AMBIG, label="ambiguous AA"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", fontsize=9,
               frameon=True, ncol=4, bbox_to_anchor=(0.5, -0.005))

    fig.suptitle("F4v2 · 4-scenario HP evolution trees (linear depth naming)",
                 fontsize=15, fontweight="bold", y=0.985)
    fig.text(0.5, 0.945,
             "Naming convention: HP1-1-1 = child-of HP1-1 (depth = number of somatic mutations from germline root)",
             ha="center", fontsize=10, color="#495057", style="italic")
    fig.text(0.99, 0.005, SOURCE_NOTE, fontsize=7, color="#ADB5BD",
             ha="right", va="bottom")

    out = FIG_DIR / "F4v2_evolution_tree_4scenarios.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[F4v2] {out}")
    return out


# ============================================================
# Main
# ============================================================
def main():
    outputs = []
    outputs.append(fig_F2v2())
    outputs.append(fig_F3v2())
    outputs.append(fig_F4v2())
    print("\n[DONE] generated:")
    for o in outputs:
        print(f"  - {o}")


if __name__ == "__main__":
    main()
