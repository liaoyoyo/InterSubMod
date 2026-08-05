#!/usr/bin/env python3
"""Generate figures for Self-Phasing PI report (2026-04-22).

Inputs:
  - /tmp/ism_hp_fix_phase1/merged_off.csv (40,115 sites, flag off)
  - /tmp/ism_hp_fix_phase1/merged_on.csv  (40,115 sites, flag on)

Outputs:
  - docs/reports/pi_reports/2026/04/figures/fig{1..10}_*.png

Figures:
  1. Paired vs TO pipeline comparison (conceptual)
  3. AF=0.3 phasing graph walkthrough (conceptual)
  4. Quantitative evidence summary (bar chart)
  5. Feature impact 3-tier matrix
  6. Fix option decision tree
  7. HPFineNGroups distribution collapse (stacked bar)
  8. Phase 1 AUC off vs on (grouped bar)
  9. HP_Ratio violin plot (off vs on, LOH/Non-LOH)
  10. Next-step decision tree

(Figure 2 reuses existing landscape doc figure.)
"""
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from sklearn.metrics import roc_auc_score

OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/pi_reports/2026/04/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_OFF = Path("/tmp/ism_hp_fix_phase1/merged_off.csv")
DATA_ON = Path("/tmp/ism_hp_fix_phase1/merged_on.csv")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "figure.dpi": 120,
})

# Color palette
C_PAIRED = "#2E7D32"  # green
C_TO = "#C62828"      # red
C_FIX = "#1565C0"     # blue
C_WARN = "#F9A825"    # amber
C_NEUTRAL = "#757575"


# -----------------------------------------------------------------------------
# Helper: rounded box with text
# -----------------------------------------------------------------------------
def box(ax, x, y, w, h, text, facecolor="#E8F5E9", edgecolor="#1B5E20",
        fontsize=9, fontweight="normal"):
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                          facecolor=facecolor, edgecolor=edgecolor, linewidth=1.3)
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, fontweight=fontweight, wrap=True)


def arrow(ax, x1, y1, x2, y2, color="#424242", style="->"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=16, color=color, linewidth=1.6))


# =============================================================================
# Figure 1: Paired vs TO pipeline comparison
# =============================================================================
def fig1_pipeline_comparison():
    fig, (axp, axt) = plt.subplots(1, 2, figsize=(14, 7))
    for ax in (axp, axt):
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis("off")

    # Paired panel
    axp.set_title("Paired Mode ✅ (Normal + Tumor)", fontsize=13, color=C_PAIRED,
                  fontweight="bold")
    axp.add_patch(Rectangle((0, 0), 10, 10, facecolor="#F1F8E9",
                            edgecolor=C_PAIRED, linewidth=2))

    boxes_p = [
        (0.5, 8.5, 9, 1.0, "Input: Normal BAM + Tumor BAM", "#FFFFFF"),
        (0.5, 7.0, 9, 1.0, "Step 1: Call germline SNPs from Normal\n(millions of high-confidence het sites)", "#E8F5E9"),
        (0.5, 5.5, 9, 1.0, "Step 2: Build phasing scaffold from germline SNPs only", "#E8F5E9"),
        (0.5, 4.0, 9, 1.0, "Step 3: Somatic variants passively attach to scaffold", "#E8F5E9"),
        (0.5, 2.5, 9, 1.0, "Step 4: Haplotag reads → HP:i:1 / HP:i:2", "#E8F5E9"),
        (0.5, 0.8, 9, 1.2, "Result: HP_Ratio reflects TRUE haplotype\n(balanced ~0.5; LOH → skewed by real CNV)",
         "#C8E6C9"),
    ]
    for x, y, w, h, t, fc in boxes_p:
        box(axp, x, y, w, h, t, facecolor=fc, edgecolor=C_PAIRED)
    for i in range(5):
        y_from = 8.5 - i * 1.5
        y_to = y_from - 0.5
        arrow(axp, 5.0, y_from, 5.0, y_to, color=C_PAIRED)

    # TO panel
    axt.set_title("TO Mode ❌ (Tumor-Only, self-phasing)", fontsize=13, color=C_TO,
                  fontweight="bold")
    axt.add_patch(Rectangle((0, 0), 10, 10, facecolor="#FFEBEE",
                            edgecolor=C_TO, linewidth=2))

    boxes_t = [
        (0.5, 8.5, 9, 1.0, "Input: Tumor BAM only (no matched Normal)", "#FFFFFF"),
        (0.5, 7.0, 9, 1.0, "Step 1: Use PON to guess germline\n(+ unlabelled somatic mixed in)", "#FFEBEE"),
        (0.5, 5.5, 9, 1.0, "Step 2: Germline + Somatic jointly build scaffold\n⚠ Somatic acts as phasing anchor", "#FFCDD2"),
        (0.5, 4.0, 9, 1.0, "Step 3: Somatic ALT reads link via shared variants\n→ forced into ONE haplotype", "#FFCDD2"),
        (0.5, 2.5, 9, 1.0, "Step 4: Haplotag → HP:i:11 / 21 / 33\n(somatic self-phase blocks)", "#FFCDD2"),
        (0.5, 0.8, 9, 1.2, "Result: 94.6% somatic reads → HP1\nHP_Ratio ARTIFICIALLY skewed", "#EF9A9A"),
    ]
    for x, y, w, h, t, fc in boxes_t:
        box(axt, x, y, w, h, t, facecolor=fc, edgecolor=C_TO)
    for i in range(5):
        y_from = 8.5 - i * 1.5
        y_to = y_from - 0.5
        arrow(axt, 5.0, y_from, 5.0, y_to, color=C_TO)

    fig.suptitle("Figure 1 — Paired vs TO Processing Pipeline",
                 fontsize=15, fontweight="bold", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig1_pipeline_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 3: AF=0.3 phasing graph walkthrough
# =============================================================================
def fig3_af03_walkthrough():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    # --- Left: Paired mode (no self-phasing) ---
    ax1.set_title("Paired: AF=0.3 somatic variant splits randomly",
                  fontsize=12, color=C_PAIRED, fontweight="bold")
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.axis("off")

    # 10 reads: 3 ALT (circles), 7 REF
    # ALT reads go randomly to HP1/HP2 based on germline SNPs around
    ax1.text(5, 9.3, "10 reads total: 3 ALT + 7 REF", ha="center", fontsize=11)
    hp_line1_y, hp_line2_y = 7.0, 2.5
    ax1.plot([0.5, 9.5], [hp_line1_y, hp_line1_y], "-", color="#66BB6A", lw=1.5)
    ax1.plot([0.5, 9.5], [hp_line2_y, hp_line2_y], "-", color="#66BB6A", lw=1.5)
    ax1.text(0.2, hp_line1_y, "HP1", fontsize=11, fontweight="bold", color=C_PAIRED, ha="right")
    ax1.text(0.2, hp_line2_y, "HP2", fontsize=11, fontweight="bold", color=C_PAIRED, ha="right")

    # ALT reads (red) split 1-2 and 2-1
    alt_positions = [(1.5, hp_line1_y), (4.5, hp_line2_y), (7.5, hp_line1_y)]
    for x, y in alt_positions:
        ax1.plot(x, y, "o", markersize=18, color=C_TO, markeredgecolor="black", zorder=5)
    # REF reads
    ref_positions = [(2.5, hp_line1_y), (3.5, hp_line2_y), (5.5, hp_line1_y),
                     (6.5, hp_line2_y), (8.5, hp_line2_y), (6.0, hp_line1_y),
                     (2.5, hp_line2_y)]
    for x, y in ref_positions:
        ax1.plot(x, y, "o", markersize=14, color="#90CAF9", markeredgecolor="black", zorder=4)

    ax1.text(5, 1.0, "HP1: 2 ALT / 2 REF | HP2: 1 ALT / 4 REF",
             ha="center", fontsize=10, fontweight="bold")
    ax1.text(5, 0.3, "HP_Ratio ≈ 0.5 (balanced, as expected)",
             ha="center", fontsize=10, color=C_PAIRED, fontweight="bold")

    # --- Right: TO mode (self-phasing collapses ALT to one side) ---
    ax2.set_title("TO: ALT reads link via shared variants → all on HP1",
                  fontsize=12, color=C_TO, fontweight="bold")
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis("off")
    ax2.text(5, 9.3, "Phasing graph edges connect ALT reads at positions X, Y, Z",
             ha="center", fontsize=11)

    # Phasing graph: 3 nodes linked as triangle
    node_xy = [(3.0, 7.0), (7.0, 7.0), (5.0, 5.0)]
    node_labels = ["ALT @ X", "ALT @ Y", "ALT @ Z"]
    for (x, y), lbl in zip(node_xy, node_labels):
        ax2.add_patch(plt.Circle((x, y), 0.5, color=C_TO, zorder=5))
        ax2.text(x, y, lbl, ha="center", va="center", fontsize=8,
                 color="white", fontweight="bold", zorder=6)
    # Triangle edges with evidence labels
    edge_labels = [("via Read A", (5.0, 7.3)), ("via Read B", (4.0, 5.8)),
                   ("via Read C", (6.0, 5.8))]
    for (xa, ya), (xb, yb) in [(node_xy[0], node_xy[1]), (node_xy[0], node_xy[2]),
                                (node_xy[1], node_xy[2])]:
        ax2.plot([xa, xb], [ya, yb], "-", color="#555", lw=2, zorder=3)
    for lbl, (x, y) in edge_labels:
        ax2.text(x, y, lbl, ha="center", fontsize=8, style="italic",
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none"))

    # Result arrow
    ax2.annotate("", xy=(5.0, 2.5), xytext=(5.0, 4.3),
                 arrowprops=dict(arrowstyle="->", color=C_TO, lw=2))
    ax2.text(5.0, 3.5, "self-phased", ha="center", fontsize=9, color=C_TO,
             fontweight="bold", bbox=dict(boxstyle="round,pad=0.2", fc="#FFF3E0", ec=C_TO))

    ax2.text(5, 1.8, "All ALT reads merged into ONE phase block",
             ha="center", fontsize=10, fontweight="bold")
    ax2.text(5, 1.0, "Haplotag assigns all ALT → HP:i:11/21 (self-phasing tag)",
             ha="center", fontsize=10)
    ax2.text(5, 0.3, "HP_Ratio → 0.94 (FALSE LOH)",
             ha="center", fontsize=10, color=C_TO, fontweight="bold")

    fig.suptitle("Figure 3 — AF=0.3 Somatic Variant: Why Self-Phasing Creates False LOH",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_af03_walkthrough.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 4: Quantitative evidence summary
# =============================================================================
def fig4_evidence_summary():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel A: 4 key evidence metrics as horizontal bar chart
    ax = axes[0]
    metrics = ["Somatic HP1:HP2 bias\n(lower = better)",
               "Somatic reads to HP1 (%)",
               "ISM HP_Ratio LOH artifact (%)",
               "Paired-TO HP_Ratio\ncorrelation (r)"]
    values = [17.3, 94.6, 62.0, 0.001]
    ideal = [1.0, 50.0, 0.0, 1.0]
    colors_bar = [C_TO, C_TO, C_TO, C_TO]

    y_pos = np.arange(len(metrics))
    ax.barh(y_pos, values, color=colors_bar, alpha=0.85, edgecolor="black")
    for i, (v, ideal_v) in enumerate(zip(values, ideal)):
        ax.axvline(ideal_v, ymin=(i) / len(metrics), ymax=(i + 1) / len(metrics),
                   color=C_PAIRED, lw=2, linestyle="--")
        ax.text(v, i, f"  {v}", va="center", fontsize=10, fontweight="bold")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(metrics, fontsize=9)
    ax.set_xlabel("Observed value (dashed green = expected under Paired)")
    ax.set_title("A. Four key pieces of evidence for self-phasing (HCC1395 TO)",
                 fontsize=11, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)
    ax.invert_yaxis()

    # Panel B: Cross-sample direction consistency (7/7)
    ax2 = axes[1]
    samples = ["HCC1395", "HCC1395_DORADO", "HCC1954", "HCC1937",
               "H1437", "H2009", "COLO829"]
    direction = [1] * 7  # all confirm
    colors_s = [C_TO if d == 1 else C_PAIRED for d in direction]
    ax2.bar(range(len(samples)), [1] * 7, color=colors_s, alpha=0.85, edgecolor="black")
    ax2.set_xticks(range(len(samples)))
    ax2.set_xticklabels(samples, rotation=30, ha="right", fontsize=9)
    ax2.set_ylim(0, 1.3)
    ax2.set_yticks([])
    ax2.set_title("B. Self-phasing observed in 7/7 samples (CV-2 pass)",
                  fontsize=11, fontweight="bold")
    for i in range(7):
        ax2.text(i, 0.5, "✓", ha="center", va="center", fontsize=20,
                 color="white", fontweight="bold")
    ax2.text(3, 1.2, "Cross-sample direction: CONSISTENT (excludes sample-specific artifact)",
             ha="center", fontsize=10, style="italic", color=C_NEUTRAL)
    for spine in ["top", "right", "left"]:
        ax2.spines[spine].set_visible(False)

    fig.suptitle("Figure 4 — Quantitative Evidence (stability level 4/5)",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig4_evidence_summary.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 5: Feature impact 3-tier matrix
# =============================================================================
def fig5_impact_matrix():
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")

    tiers = [
        ("🔴 Severe impact (unreliable)", "#FFCDD2", "#B71C1C", 1,
         ["HP_Ratio → false LOH",
          "Potential_LOH → 62% artifact",
          "HPMergedDelta / Sig → direction flip",
          "hp_assign_rate → biased high",
          "effective_hp_reads → off",
          "HPFineNGroups≥3 → artifact of somatic tags"]),
        ("🟡 Moderate impact (indirect)", "#FFF9C4", "#F57F17", 2,
         ["Quality_Score → AUC 0.497",
          "GlobalP → HP component noisy",
          "CramersV → HP component noisy",
          "Historical legacy VerificationClass (v1) → label_sig contaminated"]),
        ("🟢 No impact (conclusions solid)", "#C8E6C9", "#1B5E20", 3,
         ["PairwiseMean / MedianDist",
          "AlleleDelta / AlleleP",
          "Caller features (AF/GQ/DP/SB)",
          "Methylation matrix (raw)",
          "CpG coords, region_methyl_mean"]),
    ]

    col_w = 3.8
    for i, (title, fc, ec, _col, items) in enumerate(tiers):
        x = 0.3 + i * col_w
        # Header
        box(ax, x, 8.6, col_w - 0.3, 1.0, title, facecolor=fc, edgecolor=ec,
            fontsize=11, fontweight="bold")
        # Items
        for j, item in enumerate(items):
            y_item = 7.8 - j * 1.05
            box(ax, x + 0.15, y_item, col_w - 0.6, 0.9, item, facecolor="white",
                edgecolor=ec, fontsize=9)

    # Counts footer
    box(ax, 0.3, 0.3, 3.5, 0.9, "B class: ~29 features (38%)\nNEEDS RE-TEST",
        facecolor="#FFCDD2", edgecolor="#B71C1C", fontweight="bold", fontsize=10)
    box(ax, 4.1, 0.3, 3.5, 0.9, "C class: ~14 features (7%)\nmostly removed",
        facecolor="#FFF9C4", edgecolor="#F57F17", fontweight="bold", fontsize=10)
    box(ax, 7.9, 0.3, 3.8, 0.9, "A class: ~42 features (55%)\nSOLID (HP-free)",
        facecolor="#C8E6C9", edgecolor="#1B5E20", fontweight="bold", fontsize=10)

    ax.set_title("Figure 5 — Self-Phasing Impact Classification of ISM Features",
                 fontsize=14, fontweight="bold", pad=15)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig5_impact_matrix.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 6: Fix option decision tree
# =============================================================================
def fig6_decision_tree():
    fig, ax = plt.subplots(figsize=(13, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Root
    box(ax, 4.0, 8.5, 4.0, 1.0,
        "Problem: Self-phasing\npollutes ISM HP features",
        facecolor="#FFEBEE", edgecolor=C_TO, fontweight="bold", fontsize=11)

    # Level 1: two approaches
    box(ax, 0.5, 6.3, 4.5, 1.3,
        "Long-term: Fix LongPhase-TO\nhaplotag module\n(upstream C++, 3rd-party)",
        facecolor="#ECEFF1", edgecolor=C_NEUTRAL, fontsize=10)
    box(ax, 7.0, 6.3, 4.5, 1.3,
        "Short-term: ISM-side filter\n(this work)",
        facecolor="#E3F2FD", edgecolor=C_FIX, fontsize=10, fontweight="bold")

    # Lines from root
    arrow(ax, 5.0, 8.4, 2.75, 7.7)
    arrow(ax, 7.0, 8.4, 9.25, 7.7)

    # Long-term evaluation
    box(ax, 0.5, 4.5, 4.5, 1.3,
        "Risk: HIGH (external repo, unclear\ntimeline, deep C++ changes)\n→ DEFERRED",
        facecolor="#FFF3E0", edgecolor=C_WARN, fontsize=9)
    arrow(ax, 2.75, 6.2, 2.75, 5.85)

    # Short-term split: Option A vs B
    box(ax, 5.5, 4.5, 2.7, 1.3,
        "Option A:\nReadParser\nupstream filter\n(1 switch)",
        facecolor="#C8E6C9", edgecolor=C_PAIRED, fontsize=9, fontweight="bold")
    box(ax, 8.8, 4.5, 2.7, 1.3,
        "Option B:\nDownstream\ndispersed filter\n(4+ files)",
        facecolor="#FFEBEE", edgecolor=C_TO, fontsize=9)
    arrow(ax, 8.0, 6.2, 6.85, 5.85)
    arrow(ax, 10.0, 6.2, 10.15, 5.85)

    # Decision criteria table
    ax.text(6.0, 3.8, "Comparison:", fontsize=10, fontweight="bold")
    criteria = [
        ("Change points", "1", "4+"),
        ("Reversibility (flag off)", "full", "per-module"),
        ("Audit capability", "preserved", "native"),
        ("Regression risk", "LOW", "MEDIUM"),
    ]
    for i, (c, a, b) in enumerate(criteria):
        ax.text(5.7, 3.3 - i * 0.35, c, fontsize=8)
        ax.text(7.3, 3.3 - i * 0.35, a, fontsize=8, color=C_PAIRED, fontweight="bold")
        ax.text(9.8, 3.3 - i * 0.35, b, fontsize=8, color=C_TO)

    # Final decision
    box(ax, 3.0, 0.5, 6.0, 1.0,
        "✅ DECISION: Option A\n(`--germline-hp-only` flag, default off, audit-preserved)",
        facecolor="#C8E6C9", edgecolor=C_PAIRED, fontweight="bold", fontsize=11)
    arrow(ax, 6.85, 4.4, 6.0, 1.55, color=C_PAIRED)

    ax.set_title("Figure 6 — Fix-Option Decision Tree", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig6_fix_decision_tree.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 7: HPFineNGroups distribution collapse (data-driven)
# =============================================================================
def fig7_ngroups_collapse():
    df_off = pd.read_csv(DATA_OFF, usecols=["HPFineNGroups", "label"])
    df_on = pd.read_csv(DATA_ON, usecols=["HPFineNGroups", "label"])

    # Smoke (chr19 615 sites) — from report text
    smoke = {"off": {0: 2, 1: 3, 2: 110, 3: 394, 4: 106},
             "on": {0: 2, 1: 10, 2: 603, 3: 0, 4: 0}}

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: chr19 smoke
    ax = axes[0]
    ng_vals = [0, 1, 2, 3, 4]
    off_vals = [smoke["off"].get(k, 0) for k in ng_vals]
    on_vals = [smoke["on"].get(k, 0) for k in ng_vals]
    x = np.arange(len(ng_vals))
    w = 0.35
    ax.bar(x - w / 2, off_vals, w, label="flag off", color=C_NEUTRAL, edgecolor="black")
    ax.bar(x + w / 2, on_vals, w, label="flag on", color=C_FIX, edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels([f"NG={n}" for n in ng_vals])
    ax.set_ylabel("Count of regions")
    ax.set_title("A. Phase 0 smoke (chr19, 615 sites)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    # Panel B: Phase 1 TP
    ax = axes[1]
    tp_off = df_off[df_off.label == 1].HPFineNGroups.value_counts().sort_index()
    tp_on = df_on[df_on.label == 1].HPFineNGroups.value_counts().sort_index()
    tp_off = tp_off.reindex(ng_vals, fill_value=0)
    tp_on = tp_on.reindex(ng_vals, fill_value=0)
    ax.bar(x - w / 2, tp_off, w, label="flag off", color=C_NEUTRAL, edgecolor="black")
    ax.bar(x + w / 2, tp_on, w, label="flag on", color=C_FIX, edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels([f"NG={n}" for n in ng_vals])
    ax.set_title("B. Phase 1 HCC1395 TO — TP (28,509 sites)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    # Annotation
    ax.annotate(f"NG≥3: {int(tp_off[3] + tp_off[4])} → 0\n(100% disappear)",
                xy=(3.2, tp_off[3]), xytext=(2.5, tp_off[2] * 0.7),
                fontsize=10, fontweight="bold", color=C_TO,
                arrowprops=dict(arrowstyle="->", color=C_TO, lw=1.5))

    # Panel C: Phase 1 FP
    ax = axes[2]
    fp_off = df_off[df_off.label == 0].HPFineNGroups.value_counts().sort_index()
    fp_on = df_on[df_on.label == 0].HPFineNGroups.value_counts().sort_index()
    fp_off = fp_off.reindex(ng_vals, fill_value=0)
    fp_on = fp_on.reindex(ng_vals, fill_value=0)
    ax.bar(x - w / 2, fp_off, w, label="flag off", color=C_NEUTRAL, edgecolor="black")
    ax.bar(x + w / 2, fp_on, w, label="flag on", color=C_FIX, edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels([f"NG={n}" for n in ng_vals])
    ax.set_title("C. Phase 1 HCC1395 TO — FP (11,606 sites)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Figure 7 — HPFineNGroups Distribution Collapses When Somatic HP Tags Removed",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig7_ngroups_collapse.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 8: Phase 1 AUC off vs on comparison
# =============================================================================
def fig8_auc_comparison():
    df_off = pd.read_csv(DATA_OFF)
    df_on = pd.read_csv(DATA_ON)
    # Assert row alignment by RegionID
    assert df_off.RegionID.equals(df_on.RegionID)
    y = df_off.label.astype(int).values  # 1 = TP, 0 = FP

    features = [
        "AlleleDelta", "HP_Ratio", "HPFineNGroups", "HeuristicScore",
        "Quality_Score", "Coverage_Multiple", "CramersV", "CramersV_HPFine",
        "HPFineF", "HP1FamilyN", "HP2FamilyN", "NHP3", "NHP0", "NumReads",
        "Fisher_Frac_Sig" if "Fisher_Frac_Sig" in df_off.columns else "NHP3",
        "HPMergedDelta", "HPFine_NGroups_CF",
    ]
    features = [f for f in features if f in df_off.columns]

    # Compute AUC (best direction)
    def best_auc(y_true, scores):
        s = pd.Series(scores)
        s = s.fillna(s.median())
        if s.nunique() < 2:
            return 0.5
        try:
            a = roc_auc_score(y_true, s)
        except Exception:
            return np.nan
        return max(a, 1 - a)

    aucs_off, aucs_on = [], []
    for f in features:
        aucs_off.append(best_auc(y, df_off[f].values))
        aucs_on.append(best_auc(y, df_on[f].values))

    # Sort by off AUC descending
    order = np.argsort(aucs_off)[::-1]
    features_s = [features[i] for i in order]
    aucs_off_s = [aucs_off[i] for i in order]
    aucs_on_s = [aucs_on[i] for i in order]

    fig, ax = plt.subplots(figsize=(14, 7))
    x = np.arange(len(features_s))
    w = 0.38
    ax.bar(x - w / 2, aucs_off_s, w, label="flag off (baseline)", color=C_NEUTRAL,
           edgecolor="black")
    ax.bar(x + w / 2, aucs_on_s, w, label="flag on (germline-hp-only)", color=C_FIX,
           edgecolor="black")
    ax.axhline(0.5, color="black", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(0.55, color=C_WARN, ls="--", lw=1.2, alpha=0.7,
               label="plan Quality_Score target = 0.55")
    ax.set_xticks(x)
    ax.set_xticklabels(features_s, rotation=40, ha="right", fontsize=9)
    ax.set_ylabel("AUC (best direction)")
    ax.set_ylim(0.48, 0.68)
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    # Delta annotations
    for i, (a_off, a_on) in enumerate(zip(aucs_off_s, aucs_on_s)):
        delta = a_on - a_off
        color = C_PAIRED if delta >= 0.02 else (C_TO if delta <= -0.02 else C_NEUTRAL)
        ax.text(i, max(a_off, a_on) + 0.005, f"{delta:+.3f}",
                ha="center", fontsize=7, color=color, fontweight="bold")

    ax.set_title("Figure 8 — Phase 1 AUC: flag off vs flag on (HCC1395 TO, 40,115 sites)\n"
                 "No feature achieves Δ ≥ +0.02; 4 HP-derived features drop ≤ -0.025",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig8_phase1_auc.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    return list(zip(features_s, aucs_off_s, aucs_on_s))


# =============================================================================
# Figure 9: HP_Ratio violin plot
# =============================================================================
def fig9_hp_ratio_violin():
    df_off = pd.read_csv(DATA_OFF, usecols=["HP_Ratio", "Potential_LOH", "label"])
    df_on = pd.read_csv(DATA_ON, usecols=["HP_Ratio", "Potential_LOH", "label"])
    df_off["flag"] = "off"
    df_on["flag"] = "on"
    d = pd.concat([df_off, df_on], ignore_index=True)
    d["loh"] = d.Potential_LOH.map({True: "LOH", False: "Non-LOH"})
    d["tp_fp"] = d.label.map({1: "TP", 0: "FP"})
    d["group"] = d["tp_fp"] + " / " + d["loh"]

    fig, ax = plt.subplots(figsize=(13, 6.5))
    groups = ["TP / Non-LOH", "TP / LOH", "FP / Non-LOH", "FP / LOH"]
    positions = []
    data_lists = []
    labels = []
    colors = []
    for i, g in enumerate(groups):
        for j, flag in enumerate(["off", "on"]):
            sel = d[(d.group == g) & (d.flag == flag)]
            arr = sel.HP_Ratio.dropna().values
            if arr.size < 2 or np.allclose(arr.min(), arr.max()):
                continue  # violin needs non-degenerate data
            positions.append(i * 3 + j)
            data_lists.append(arr)
            labels.append(f"{g}\n(flag {flag})")
            colors.append(C_NEUTRAL if flag == "off" else C_FIX)

    parts = ax.violinplot(data_lists, positions=positions, showmedians=True, widths=0.8)
    for pc, color in zip(parts["bodies"], colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.6)
        pc.set_edgecolor("black")
    for part_name in ["cmins", "cmaxes", "cbars", "cmedians"]:
        if part_name in parts:
            parts[part_name].set_color("black")

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("HP_Ratio (min allele frac on dominant HP)")
    ax.axhline(0.5, color="black", ls="--", lw=0.8, alpha=0.5, label="balanced = 0.5")
    ax.axhline(0.1, color=C_TO, ls=":", lw=1, alpha=0.5, label="LOH threshold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(-0.05, 1.05)

    # Annotation: median shift on Non-LOH
    ax.annotate("Non-LOH median\n0.554 → 0.531\n(Δ = -0.023)",
                xy=(0.5, 0.54), xytext=(0.5, 0.78), fontsize=9, ha="center",
                color=C_FIX, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="#E3F2FD", ec=C_FIX),
                arrowprops=dict(arrowstyle="->", color=C_FIX))

    ax.set_title("Figure 9 — HP_Ratio distribution: flag off vs on\n"
                 "Non-LOH shifts toward 0.5 (expected direction); LOH barely moves (already extreme)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig9_hp_ratio_violin.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 10: Next-step decision tree
# =============================================================================
def fig10_next_step_tree():
    fig, ax = plt.subplots(figsize=(13, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Root
    box(ax, 4.0, 8.7, 4.0, 1.0,
        "Phase 1 result: filter NEGATIVE\nNo feature gains AUC ≥ +0.02",
        facecolor="#FFF3E0", edgecolor=C_WARN, fontweight="bold", fontsize=11)

    # 4 options
    options = [
        (0.1, 5.8, 2.7, 1.4, "Phase 2:\n7 samples x 2 modes\nfull rerun", "#FFCDD2",
         C_TO, "STOP"),
        (3.1, 5.8, 2.7, 1.4, "P2: rebuild\nwithin_dom_alt_frac\nPython pipeline", "#FFF9C4",
         C_WARN, "SKIP"),
        (6.1, 5.8, 2.7, 1.4, "P3: Change flag\ndefault to ON\nin future release", "#FFF9C4",
         C_WARN, "STATUS QUO"),
        (9.1, 5.8, 2.7, 1.4, "P4: Master dataset\nx 2 flags x 2 modes\n(if publishing marker)", "#E3F2FD",
         C_FIX, "DEFER"),
    ]
    for x, y, w, h, t, fc, ec, verdict in options:
        box(ax, x, y, w, h, t, facecolor=fc, edgecolor=ec, fontsize=9)
        # Verdict pill below each option (filled with ec color, white text)
        rect = FancyBboxPatch((x + 0.3, y - 1.2), w - 0.6, 0.9,
                              boxstyle="round,pad=0.02", facecolor=ec, edgecolor=ec)
        ax.add_patch(rect)
        ax.text(x + w / 2, y - 1.2 + 0.45, verdict, ha="center", va="center",
                fontsize=10, fontweight="bold", color="white")

    # Connect root to options
    for x, y, w, _, _, _, _, _ in options:
        arrow(ax, 6.0, 8.6, x + w / 2, y + 1.4)

    # Reasoning
    reasons = [
        (0.1, 2.8, 2.7, "Phase 1 gate failed\ndefinitively.\nNo reversal potential."),
        (3.1, 2.8, 2.7, "Downstream rebuild\n~1 day. Expected\ngain < +0.02."),
        (6.1, 2.8, 2.7, "Would immediately\nbreak HPFineNGroups\nsubclone applications."),
        (9.1, 2.8, 2.7, "Required only if\npublishing marker\npaper. Not now."),
    ]
    for x, y, w, r in reasons:
        box(ax, x, y, w, 1.4, r, facecolor="white", edgecolor=C_NEUTRAL, fontsize=8)

    # Recommended direction
    box(ax, 2.5, 0.3, 7.0, 1.4,
        "Recommended: return to higher-ROI tracks\n(CovM baseline blocker / HCC1395 chr8 hotspot / Phase 2 characterization)",
        facecolor="#C8E6C9", edgecolor=C_PAIRED, fontweight="bold", fontsize=10)

    ax.set_title("Figure 10 — Next-Step Decision Tree After Phase 1 Gate",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig10_next_steps.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 11 — Three-track × two-flag impact matrix
# =============================================================================
def fig11_three_track_matrix():
    """2D heatmap showing flag effect across three pipeline tracks × analysis dimensions.

    Rows: seven representative analysis dimensions
    Cols: three tracks × two flag states = 6 columns
    Cells colored by impact severity (unaffected / no-op / mechanism-only / AUC-affected)
    """
    fig, ax = plt.subplots(figsize=(14, 7))

    # Analysis dimensions (rows)
    rows = [
        "Baseline F1 (TP/FP discrimination)",
        "HP_Ratio distribution",
        "HPFineNGroups bucket counts",
        "Methylation features (PairwiseMean/Allele*)",
        "LOH.bed region-level LOH",
        "Downstream phasing (VCF/BAM output)",
        "ISM TSV audit columns (NHP_Somatic*)",
    ]
    # Columns: track × flag
    cols = [
        ("paired_full", "off"), ("paired_full", "on"),
        ("paired_pileup", "off"), ("paired_pileup", "on"),
        ("to_pileup", "off"), ("to_pileup", "on"),
    ]

    # Impact codes: 0=baseline, 1=no-op, 2=expected shift (intended), 3=unaffected
    # Severity colors: green=unaffected, grey=no-op, blue=intended shift
    codes = np.array([
        # paired_full off/on, paired_pileup off/on, to_pileup off/on
        [0, 1, 0, 1, 0, 1],   # Baseline F1: paired no-op; TO no-AUC-delta (已 Phase1 驗證)
        [0, 1, 0, 1, 0, 2],   # HP_Ratio: paired 不變; TO Non-LOH 朝 0.5 (intended)
        [0, 1, 0, 1, 0, 2],   # NGroups: paired 不變; TO NG≥3 歸零 (intended reveals mechanism)
        [3, 3, 3, 3, 3, 3],   # Methylation: 完全不依賴 HP
        [3, 3, 3, 3, 3, 3],   # LOH.bed: Jaccard=1.0 已證
        [3, 3, 3, 3, 3, 3],   # Downstream phasing: flag 在 ReadParser 下游
        [0, 1, 0, 1, 0, 1],   # Audit: 永遠反映 raw
    ])

    color_map = {
        0: "#E0E0E0",  # baseline (off, unchanged reference)
        1: "#C8E6C9",  # no-op (flag has no effect)
        2: "#90CAF9",  # intended shift (TO with flag on)
        3: "#81C784",  # unaffected (structurally independent)
    }
    label_map = {
        0: "baseline",
        1: "no-op",
        2: "intended\nshift",
        3: "unaffected",
    }

    for i in range(codes.shape[0]):
        for j in range(codes.shape[1]):
            code = codes[i, j]
            ax.add_patch(Rectangle((j, len(rows) - 1 - i), 1, 1,
                                   facecolor=color_map[code],
                                   edgecolor="black", linewidth=0.6))
            ax.text(j + 0.5, len(rows) - 1 - i + 0.5, label_map[code],
                    ha="center", va="center", fontsize=8)

    # Track group headers (top)
    track_headers = [("paired_full", 0, 2), ("paired_pileup", 2, 2), ("to_pileup", 4, 2)]
    track_colors = {"paired_full": "#E8F5E9", "paired_pileup": "#FFF9C4", "to_pileup": "#FFEBEE"}
    for name, x0, w in track_headers:
        ax.add_patch(Rectangle((x0, len(rows)), w, 0.6,
                               facecolor=track_colors[name], edgecolor="black"))
        ax.text(x0 + w / 2, len(rows) + 0.3, name,
                ha="center", va="center", fontsize=11, fontweight="bold")

    # Flag sub-headers
    for j, (_, flag) in enumerate(cols):
        ax.text(j + 0.5, len(rows) + 0.8, f"flag {flag}",
                ha="center", va="center", fontsize=9, style="italic")

    # Row labels
    for i, row in enumerate(rows):
        ax.text(-0.1, len(rows) - 1 - i + 0.5, row,
                ha="right", va="center", fontsize=9)

    # Legend
    legend_items = [
        mpatches.Patch(facecolor="#E0E0E0", edgecolor="black", label="baseline (unchanged)"),
        mpatches.Patch(facecolor="#C8E6C9", edgecolor="black", label="no-op (no target tag)"),
        mpatches.Patch(facecolor="#90CAF9", edgecolor="black", label="intended shift"),
        mpatches.Patch(facecolor="#81C784", edgecolor="black", label="structurally unaffected"),
    ]
    ax.legend(handles=legend_items, loc="lower center",
              bbox_to_anchor=(0.5, -0.18), ncol=4, fontsize=9, frameon=False)

    ax.set_xlim(-7, len(cols) + 0.5)
    ax.set_ylim(-0.5, len(rows) + 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Figure 11 — Three Pipeline Tracks × Two Flag States: Impact Matrix",
                 fontsize=13, fontweight="bold", pad=20)

    # Footer annotation
    ax.text(len(cols) / 2, -0.9,
            "Key insight: paired tracks are no-op; only to_pileup gets intended shifts. "
            "Methylation & LOH.bed structurally unaffected.",
            ha="center", va="top", fontsize=9, style="italic", color="#424242")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig11_three_track_matrix.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 12 — Adversarial review decision tree (7 challenges × 3 perspectives)
# =============================================================================
def fig12_adversarial_tree():
    """Layered tree: 3 perspectives → 7 challenges → verdicts."""
    fig, ax = plt.subplots(figsize=(16, 10))

    # Perspectives at top
    perspectives = [
        ("Code-level\n(C1-C3)", 4, 9, "#1565C0"),
        ("Biology\n(B1-B3)", 9, 9, "#2E7D32"),
        ("Paper methodology\n(M1-M3)", 14, 9, "#6A1B9A"),
    ]
    for name, x, y, color in perspectives:
        box(ax, x - 1.3, y - 0.5, 2.6, 1.0, name,
            facecolor="white", edgecolor=color, fontweight="bold", fontsize=11)

    # Challenges (middle row)
    challenges = [
        ("C1", 2.0, 6.5, "PON quality?",
         "Jaccard=1.0 PON-only", "Structural, not\nPON-quality issue"),
        ("C2", 4.0, 6.5, "Break downstream?",
         "Phase 0 conservation", "Mechanism\nclean"),
        ("C3", 6.0, 6.5, "Why keep flag?",
         "NG mechanism\nreveal", "Research-level\nnull test tool"),
        ("B1", 8.0, 6.5, "Tumor biology?",
         "ITH clonal\nexpansion", "Consistent with\nknown biology"),
        ("B2", 10.0, 6.5, "Biological meaning\nof HP:i:11/21/33?",
         "GT=0|0 parsing\nedge case", "None; algorithmic\nartifact"),
        ("B3", 12.0, 6.5, "Lose subclone info?",
         "NG is bucket count,\nnot methylation", "Flag is negative\ncontrol (pivotal)"),
        ("M1", 13.5, 6.5, "LongPhase design?",
         "All phasers\ngermline-focused", "TO extension\nknown side-effect"),
        ("M2", 15.0, 6.5, "Others solve?",
         "ClairS-TO also\nno anchor filter", "ISM novelty"),
        ("M3", 16.5, 6.5, "Fix upstream?",
         "3-6mo vs 1-2d\ntrade-off", "Downstream\npragmatic"),
    ]

    for cid, x, y, question, evidence, verdict in challenges:
        # Question box
        q_color = "#FFF3E0" if cid.startswith("C") else ("#F1F8E9" if cid.startswith("B") else "#F3E5F5")
        q_edge = "#1565C0" if cid.startswith("C") else ("#2E7D32" if cid.startswith("B") else "#6A1B9A")
        box(ax, x - 0.7, y - 0.3, 1.4, 0.6, f"{cid}: {question}",
            facecolor=q_color, edgecolor=q_edge, fontsize=7, fontweight="bold")
        # Evidence box
        ax.annotate("", xy=(x, y - 0.8), xytext=(x, y - 0.3),
                    arrowprops=dict(arrowstyle="->", color="grey"))
        box(ax, x - 0.7, y - 1.6, 1.4, 0.7, evidence,
            facecolor="white", edgecolor="grey", fontsize=6.5)
        # Verdict box
        ax.annotate("", xy=(x, y - 2.1), xytext=(x, y - 1.6),
                    arrowprops=dict(arrowstyle="->", color="#1B5E20"))
        box(ax, x - 0.7, y - 2.9, 1.4, 0.7, verdict,
            facecolor="#C8E6C9", edgecolor="#1B5E20", fontsize=6.5, fontweight="bold")

    # Lines from perspectives to their challenges
    def draw_line(x0, y0, x1, y1, color):
        ax.plot([x0, x1], [y0, y1], "-", color=color, linewidth=1.2, alpha=0.5)
    # Code perspective to C1-C3
    for _, x, y, *_ in challenges[:3]:
        draw_line(perspectives[0][1], perspectives[0][2] - 0.5, x, y + 0.3, "#1565C0")
    # Biology perspective to B1-B3
    for _, x, y, *_ in challenges[3:6]:
        draw_line(perspectives[1][1], perspectives[1][2] - 0.5, x, y + 0.3, "#2E7D32")
    # Methodology perspective to M1-M3
    for _, x, y, *_ in challenges[6:9]:
        draw_line(perspectives[2][1], perspectives[2][2] - 0.5, x, y + 0.3, "#6A1B9A")

    # Bottom convergence
    box(ax, 6, 1.5, 6, 1.3,
        "CONVERGENT VERDICT:\n"
        "Self-phasing is real, fix is correct,\n"
        "three pipeline tracks unaffected",
        facecolor="#1B5E20", edgecolor="black", fontsize=11, fontweight="bold")
    # White override for dark bg text
    for txt in ax.texts:
        pass  # keep as is; FancyBboxPatch + text above
    # Add explicit white text
    ax.text(9, 2.15,
            "CONVERGENT VERDICT:\n"
            "Self-phasing is real, fix is correct,\n"
            "three pipeline tracks unaffected",
            ha="center", va="center", fontsize=11, fontweight="bold", color="white")

    # Arrows from verdicts to convergence
    for cid, x, y, *_ in challenges:
        ax.annotate("", xy=(9, 3), xytext=(x, y - 2.9),
                    arrowprops=dict(arrowstyle="->", color="grey", alpha=0.25))

    ax.set_xlim(0, 18)
    ax.set_ylim(0.5, 10.5)
    ax.axis("off")
    ax.set_title("Figure 12 — Adversarial Review Decision Tree: 3 Perspectives × 9 Challenges",
                 fontsize=13, fontweight="bold", pad=10)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig12_adversarial_tree.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 13 — Full-case storyboard: Problem → Fix → Result
# =============================================================================
def fig13_case_storyboard():
    """Three-panel storyboard tracing ONE concrete case (AF=0.3 somatic variant)
    through: (A) Problem state (flag=off, self-phasing), (B) Fix mechanism
    (demotion logic), (C) Result state (flag=on, artifact resolved).

    This figure consolidates the before/after narrative that otherwise would
    require jumping between Fig 3 / Fig 7 / Fig 8 / Fig 9.
    """
    fig = plt.figure(figsize=(18, 10))

    # Layout: 3 panels side by side + top title + bottom summary
    gs = fig.add_gridspec(2, 3, height_ratios=[5, 1], hspace=0.3, wspace=0.25)

    # ---- Panel A: Problem state ----
    axA = fig.add_subplot(gs[0, 0])
    axA.set_facecolor("#FFEBEE")
    axA.set_xlim(0, 10)
    axA.set_ylim(0, 10)
    axA.set_title("A. PROBLEM (flag=off)\nSelf-phasing active", fontsize=12,
                  fontweight="bold", color="#C62828")

    # Draw 10 reads for AF=0.3 variant
    # 3 ALT reads clustered HP1, 7 REF reads split
    reads_A = [
        # (x, y, is_alt, assigned_hp, color)
        (1, 8, True, "HP1-1", "#D32F2F"),   # ALT
        (2.5, 8, True, "HP1-1", "#D32F2F"), # ALT
        (4, 8, True, "HP1-1", "#D32F2F"),   # ALT (all forced HP1-1)
        (1, 6, False, "HP1", "#1976D2"),    # REF
        (2.5, 6, False, "HP1", "#1976D2"),  # REF
        (4, 6, False, "HP1", "#1976D2"),    # REF
        (5.5, 6, False, "HP1", "#1976D2"),  # REF
        (1, 4, False, "HP2", "#388E3C"),    # REF
        (2.5, 4, False, "HP2", "#388E3C"),  # REF
        (4, 4, False, "HP2", "#388E3C"),    # REF
    ]
    for x, y, is_alt, hp, color in reads_A:
        axA.add_patch(Rectangle((x, y - 0.3), 1.3, 0.6, facecolor=color,
                                edgecolor="black", linewidth=0.8))
        symbol = "ALT" if is_alt else "REF"
        axA.text(x + 0.65, y, f"{symbol}\n{hp}", ha="center", va="center",
                 fontsize=7, color="white", fontweight="bold")

    # HP_Ratio and NG annotations
    axA.text(7.5, 8, "3 ALT\n→ HP1-1\n(forced)", ha="center", va="center",
             fontsize=9, fontweight="bold", color="#C62828",
             bbox=dict(facecolor="white", edgecolor="#C62828", boxstyle="round,pad=0.3"))
    axA.text(7.5, 5, "HP buckets seen:\n{HP1, HP1-1, HP2}\nNG=3", ha="center", va="center",
             fontsize=9, fontweight="bold",
             bbox=dict(facecolor="white", edgecolor="grey", boxstyle="round,pad=0.3"))
    axA.text(5, 1.5, "HP_Ratio = min(HP1+HP1-1, HP2) / total\n"
                     "= min(7, 3) / 10 = 0.30\n"
                     "⚠ skewed toward HP1 (artifact)",
             ha="center", va="center", fontsize=9, color="#C62828", fontweight="bold",
             bbox=dict(facecolor="white", edgecolor="#C62828", boxstyle="round,pad=0.3"))
    axA.axis("off")

    # ---- Panel B: Fix mechanism ----
    axB = fig.add_subplot(gs[0, 1])
    axB.set_facecolor("#E3F2FD")
    axB.set_xlim(0, 10)
    axB.set_ylim(0, 10)
    axB.set_title("B. FIX (ReadParser demotion logic)\n--germline-hp-only=true",
                  fontsize=12, fontweight="bold", color="#1565C0")

    # Code-like block
    code_text = (
        "if (config.germline_hp_only &&\n"
        "    (hp_raw == \"1-1\" ||\n"
        "     hp_raw == \"2-1\" ||\n"
        "     hp_raw == \"3\")) {\n"
        "    info.hp_tag = \"0\";  ← demote\n"
        "} else {\n"
        "    info.hp_tag = hp_raw;\n"
        "}\n"
        "\n"
        "// raw preserved:\n"
        "info.hp_tag_raw = hp_raw;"
    )
    axB.text(0.5, 8.5, code_text, ha="left", va="top", fontsize=9,
             family="monospace",
             bbox=dict(facecolor="white", edgecolor="#1565C0", boxstyle="round,pad=0.5"))

    # Arrow showing demotion
    axB.annotate("", xy=(5, 3.2), xytext=(5, 4.5),
                 arrowprops=dict(arrowstyle="->", color="#1565C0", lw=2))

    axB.text(5, 2.3, "HP:i:11/21/33 → \"0\"\n(audit via NHP_Somatic11/21/33)\n"
                     "Phase 0 mass conservation PASS",
             ha="center", va="center", fontsize=9, color="#1565C0", fontweight="bold",
             bbox=dict(facecolor="white", edgecolor="#1565C0", boxstyle="round,pad=0.3"))
    axB.axis("off")

    # ---- Panel C: Result state ----
    axC = fig.add_subplot(gs[0, 2])
    axC.set_facecolor("#E8F5E9")
    axC.set_xlim(0, 10)
    axC.set_ylim(0, 10)
    axC.set_title("C. RESULT (flag=on)\nArtifact removed",
                  fontsize=12, fontweight="bold", color="#2E7D32")

    # After demotion: ALT reads go to HP=0 (unphased)
    reads_C = [
        (1, 8, True, "0", "#757575"),       # ALT → unphased
        (2.5, 8, True, "0", "#757575"),     # ALT → unphased
        (4, 8, True, "0", "#757575"),       # ALT → unphased
        (1, 6, False, "HP1", "#1976D2"),
        (2.5, 6, False, "HP1", "#1976D2"),
        (4, 6, False, "HP1", "#1976D2"),
        (5.5, 6, False, "HP1", "#1976D2"),
        (1, 4, False, "HP2", "#388E3C"),
        (2.5, 4, False, "HP2", "#388E3C"),
        (4, 4, False, "HP2", "#388E3C"),
    ]
    for x, y, is_alt, hp, color in reads_C:
        axC.add_patch(Rectangle((x, y - 0.3), 1.3, 0.6, facecolor=color,
                                edgecolor="black", linewidth=0.8))
        symbol = "ALT" if is_alt else "REF"
        axC.text(x + 0.65, y, f"{symbol}\n{hp}", ha="center", va="center",
                 fontsize=7, color="white", fontweight="bold")

    axC.text(7.5, 8, "3 ALT\n→ HP=0\n(unphased)", ha="center", va="center",
             fontsize=9, fontweight="bold", color="#2E7D32",
             bbox=dict(facecolor="white", edgecolor="#2E7D32", boxstyle="round,pad=0.3"))
    axC.text(7.5, 5, "HP buckets seen:\n{HP1, HP2}\nNG=2", ha="center", va="center",
             fontsize=9, fontweight="bold",
             bbox=dict(facecolor="white", edgecolor="grey", boxstyle="round,pad=0.3"))
    axC.text(5, 1.5, "HP_Ratio = min(HP1, HP2) / (HP1+HP2)\n"
                     "= min(4, 3) / 7 = 0.43\n"
                     "✓ closer to balanced (0.5)",
             ha="center", va="center", fontsize=9, color="#2E7D32", fontweight="bold",
             bbox=dict(facecolor="white", edgecolor="#2E7D32", boxstyle="round,pad=0.3"))
    axC.axis("off")

    # ---- Bottom summary row ----
    ax_sum = fig.add_subplot(gs[1, :])
    ax_sum.set_xlim(0, 10)
    ax_sum.set_ylim(0, 1)
    ax_sum.axis("off")

    summary = (
        "Before → After summary for this case (AF=0.3 somatic variant, 10 reads):  "
        "HP_Ratio 0.30 → 0.43 (closer to expected 0.5)  |  "
        "NGroups 3 → 2 (artifact HP1-1 bucket removed)  |  "
        "Raw HP tags preserved in hp_tag_raw for audit  |  "
        "Phase 0 mass conservation: NHP0 gain = Σ(NHP_Somatic11/21/33) ✓"
    )
    ax_sum.text(5, 0.5, summary, ha="center", va="center", fontsize=10,
                fontweight="bold",
                bbox=dict(facecolor="#FFF9C4", edgecolor="#F57F17",
                          boxstyle="round,pad=0.5"))

    fig.suptitle("Figure 13 — End-to-End Case Storyboard: One AF=0.3 Somatic Variant "
                 "Through Problem → Fix → Result",
                 fontsize=14, fontweight="bold", y=0.98)
    fig.savefig(OUT_DIR / "fig13_case_storyboard.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 14 — PON-only phasing upstream-fix storyboard (PPT-ready, 4 panels)
# =============================================================================
def fig14_pon_only_storyboard():
    """Four-panel storyboard explaining PON-only phasing (upstream fix).

    Panels (PPT-optimised):
      A. Problem    — self-phasing mechanism in TO mode
      B. PON-fix    — convertNonGermlineToSomatic() upstream logic
      C. VCF result — Jaccard=1.0, bias eliminated, N50 +99.7%
      D. Haplotag caveat + Combined solution (PON-only + ISM --germline-hp-only)

    Optimised for 16:9 PPT slide: low information density per panel,
    large fonts, headline-style verdicts.
    """
    fig = plt.figure(figsize=(20, 11))
    gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.18)

    # -------------------- Panel A: Problem --------------------
    axA = fig.add_subplot(gs[0, 0])
    axA.set_facecolor("#FFEBEE")
    axA.set_xlim(0, 10)
    axA.set_ylim(0, 10)
    axA.set_title("A. PROBLEM: Self-phasing in TO mode",
                  fontsize=15, fontweight="bold", color="#C62828")

    # Left: variants entering graph
    axA.text(1.5, 9, "TO-mode VCF input:", fontsize=11, fontweight="bold")
    var_boxes = [
        (0.3, 7.5, "germline", "#81C784"),
        (2.0, 7.5, "germline", "#81C784"),
        (3.7, 7.5, "somatic", "#EF5350"),
        (0.3, 6.3, "somatic", "#EF5350"),
        (2.0, 6.3, "unknown", "#BDBDBD"),
        (3.7, 6.3, "somatic", "#EF5350"),
    ]
    for x, y, label, color in var_boxes:
        axA.add_patch(Rectangle((x, y), 1.4, 0.6, facecolor=color,
                                edgecolor="black", linewidth=0.8))
        axA.text(x + 0.7, y + 0.3, label, ha="center", va="center",
                 fontsize=8, color="white", fontweight="bold")

    # Arrow to phasing graph
    axA.annotate("", xy=(5.5, 6.8), xytext=(5.1, 6.8),
                 arrowprops=dict(arrowstyle="->", color="black", lw=2))

    # Right: "all as anchors" graph with strong somatic-somatic edges
    axA.text(7.5, 9, "All as anchors →", fontsize=11, fontweight="bold")
    # Draw somatic cluster with heavy edges
    somatic_pos = [(6.5, 7.5), (8.5, 7.5), (7.5, 5.5)]
    for x, y in somatic_pos:
        axA.add_patch(mpatches.Circle((x, y), 0.35, facecolor="#EF5350",
                                      edgecolor="black", linewidth=0.8))
        axA.text(x, y, "S", ha="center", va="center", fontsize=10,
                 color="white", fontweight="bold")
    # Heavy red edges
    for i in range(len(somatic_pos)):
        for j in range(i + 1, len(somatic_pos)):
            x1, y1 = somatic_pos[i]
            x2, y2 = somatic_pos[j]
            axA.plot([x1, x2], [y1, y2], "-", color="#C62828", lw=3.5, alpha=0.8)

    axA.text(7.5, 4.2, "Somatic edges\nDOMINATE graph",
             ha="center", va="center", fontsize=10, color="#C62828",
             fontweight="bold",
             bbox=dict(facecolor="white", edgecolor="#C62828",
                       boxstyle="round,pad=0.3"))

    # Bottom verdict
    axA.text(5, 2, "★ Consequence ★\n"
                   "94.6% somatic reads forced → HP1  (bias 17.3 : 1)\n"
                   "Cross-mode HP_Ratio correlation  r = 0.001\n"
                   "62% TO-mode ISM LOH = artifact",
             ha="center", va="center", fontsize=11, fontweight="bold",
             color="#B71C1C",
             bbox=dict(facecolor="white", edgecolor="#B71C1C",
                       boxstyle="round,pad=0.5", linewidth=2))
    axA.axis("off")

    # -------------------- Panel B: PON-only fix --------------------
    axB = fig.add_subplot(gs[0, 1])
    axB.set_facecolor("#E3F2FD")
    axB.set_xlim(0, 10)
    axB.set_ylim(0, 10)
    axB.set_title("B. FIX: --pon-only-phasing (upstream)",
                  fontsize=15, fontweight="bold", color="#1565C0")

    axB.text(1.5, 9, "Upstream filter before phasing graph:",
             fontsize=11, fontweight="bold")

    # Code-like explanation
    code_B = (
        "LongPhase-TO / PhasingProcess.cpp:\n"
        "\n"
        "  convertNonGermlineToSomatic()\n"
        "    ▸ PON-confirmed germline ── primary anchor\n"
        "    ▸ All other variants ─── reduced edge weight\n"
        "\n"
        "CLI: --pon-only-phasing"
    )
    axB.text(0.5, 7.2, code_B, fontsize=10, family="monospace",
             verticalalignment="top",
             bbox=dict(facecolor="white", edgecolor="#1565C0",
                       boxstyle="round,pad=0.5"))

    # Arrow showing transformation
    axB.annotate("", xy=(5, 3.8), xytext=(5, 4.6),
                 arrowprops=dict(arrowstyle="->", color="#1565C0", lw=3))

    # Bottom result (graph with only germline anchors)
    germline_pos = [(3, 2.5), (5, 2.5), (7, 2.5)]
    for x, y in germline_pos:
        axB.add_patch(mpatches.Circle((x, y), 0.35, facecolor="#81C784",
                                      edgecolor="black", linewidth=0.8))
        axB.text(x, y, "G", ha="center", va="center", fontsize=10,
                 color="white", fontweight="bold")
    # Light green edges
    for i in range(len(germline_pos)):
        for j in range(i + 1, len(germline_pos)):
            x1, y1 = germline_pos[i]
            x2, y2 = germline_pos[j]
            axB.plot([x1, x2], [y1, y2], "-", color="#2E7D32", lw=2)

    axB.text(5, 0.8, "Germline anchors only → clean scaffold",
             ha="center", va="center", fontsize=11, fontweight="bold",
             color="#1B5E20",
             bbox=dict(facecolor="white", edgecolor="#1B5E20",
                       boxstyle="round,pad=0.3"))
    axB.axis("off")

    # -------------------- Panel C: VCF layer results --------------------
    axC = fig.add_subplot(gs[1, 0])
    axC.set_facecolor("#E8F5E9")
    axC.set_xlim(0, 10)
    axC.set_ylim(0, 10)
    axC.set_title("C. RESULT (VCF layer): Clean, faster, higher quality",
                  fontsize=15, fontweight="bold", color="#2E7D32")

    # Results table
    results = [
        ("LOH.bed Jaccard",          "—",          "1.0000",       "✓ identical"),
        ("Somatic HP1:HP2 bias",     "17.3 : 1",    "eliminated",   "✓ self-phasing resolved"),
        ("Phase block N50",          "4,061",       "8,109",        "✓ +99.7%"),
        ("Phased rate",              "54.9%",       "78.5%",        "✓ +23.6 pp"),
        ("Runtime",                  "2,693 s",     "1,976 s",      "✓ 1.36× faster"),
    ]
    col_x = [0.3, 3.3, 5.3, 7.3]
    headers = ["Metric", "Baseline", "PON-only", "Verdict"]
    # header row
    for cx, h in zip(col_x, headers):
        axC.text(cx, 8.7, h, fontsize=11, fontweight="bold", color="#1B5E20")
    axC.plot([0.2, 9.8], [8.4, 8.4], "-", color="#1B5E20", lw=1.5)
    # data rows
    for i, row in enumerate(results):
        y = 7.5 - i * 1.2
        for cx, val in zip(col_x, row):
            axC.text(cx, y, val, fontsize=10, verticalalignment="center")
        axC.plot([0.2, 9.8], [y - 0.5, y - 0.5], ":", color="grey", lw=0.5)

    axC.text(5, 0.7, "★ VCF layer: PERFECTLY FIXED ★",
             ha="center", va="center", fontsize=12, fontweight="bold",
             color="#1B5E20",
             bbox=dict(facecolor="white", edgecolor="#1B5E20",
                       boxstyle="round,pad=0.5", linewidth=2))
    axC.axis("off")

    # -------------------- Panel D: Haplotag caveat + Combined solution --------------------
    axD = fig.add_subplot(gs[1, 1])
    axD.set_facecolor("#FFF9C4")
    axD.set_xlim(0, 10)
    axD.set_ylim(0, 10)
    axD.set_title("D. CAVEAT & COMBINED SOLUTION",
                  fontsize=15, fontweight="bold", color="#F57C00")

    # Top: Haplotag caveat
    axD.text(0.5, 9.3, "⚠ Haplotag layer caveat (new artifact)",
             fontsize=12, fontweight="bold", color="#E65100")

    caveat = (
        "After PON-only phasing, somatic sites become GT=0|0.\n"
        "LongPhase haplotag sets refHaplotype=UNDEFINED →\n"
        "all reads tagged as a single HP (usually HP:i:21).\n"
        "\n"
        "ISM HP_Ratio TP median:\n"
        "  Paired    = 0.5000   ← expected baseline\n"
        "  Baseline  = 0.8358   ← self-phasing bias\n"
        "  PON-only  = 0.0000   ← new artifact (opposite extreme)\n"
        "\n"
        "ISM-only LOH excess:  15.4% → 54.8%   ⚠"
    )
    axD.text(0.5, 8.7, caveat, fontsize=9.5, family="monospace",
             verticalalignment="top",
             bbox=dict(facecolor="white", edgecolor="#E65100",
                       boxstyle="round,pad=0.4"))

    # Bottom: Combined solution
    axD.text(0.5, 3.5, "✓ Combined solution (production)",
             fontsize=12, fontweight="bold", color="#1B5E20")

    combined = (
        "  Upstream:  --pon-only-phasing   (LongPhase-TO)\n"
        "             └─► fixes VCF layer\n"
        "\n"
        "  Downstream: --germline-hp-only  (ISM ReadParser)\n"
        "             └─► bypasses haplotag artifact"
    )
    axD.text(0.5, 2.9, combined, fontsize=10, family="monospace",
             verticalalignment="top",
             bbox=dict(facecolor="white", edgecolor="#1B5E20",
                       boxstyle="round,pad=0.4"))
    axD.axis("off")

    # Overall title
    fig.suptitle("Figure 14 — PON-only Phasing (Upstream Fix): "
                 "Problem → Fix → VCF Result → Haplotag Caveat → Combined",
                 fontsize=16, fontweight="bold", y=0.995)

    fig.savefig(OUT_DIR / "fig14_pon_only_storyboard.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 15 — Evidence chain dependency matrix (methodology)
# =============================================================================
def fig15_evidence_dependency_matrix():
    """Dependency matrix answering: 'Does our self-phasing judgment rely on
    Paired HP1/HP2 balance as ground truth?'

    6 evidence lines × 4 dependency dimensions:
      - Paired as reference (does this evidence use Paired HP_Ratio as truth?)
      - HP1/HP2 balance metric (does it use the ratio directly?)
      - Biological null hypothesis (does it rely on 50:50 expectation?)
      - Algorithm/theory (does it derive from phasing graph math?)

    Color coding:
      green = independent from dimension, red = dependent, amber = partial
    """
    fig, ax = plt.subplots(figsize=(16, 9))

    evidence_lines = [
        ("E1: Somatic HP1:HP2 bias = 17.3:1",       "genome-wide read-level count"),
        ("E2: Cross-mode HP_Ratio  r = 0.001",      "Paired × TO same-site pairs"),
        ("E3: 62% TO LOH = artifact",               "TO LOH - (Paired non-LOH intersect)"),
        ("E4: 7/7 samples direction consistency",   "CV-2 within-TO replication"),
        ("E5: PON-only phasing experiment",         "perturbation: remove somatic, observe"),
        ("E6: Phasing graph edge weight derivation","phasing math: w = Σ co-occurrence"),
    ]
    dimensions = [
        ("Paired\nas reference",      "paired"),
        ("HP1/HP2 balance\n(HP_Ratio)", "hp_ratio"),
        ("Biological null\n(50:50 expectation)", "bio_null"),
        ("Algorithm /\ntheory", "theory"),
    ]
    # Dependency codes: 0=independent (green), 1=partial (amber), 2=dependent (red)
    # Rows: E1..E6; Columns: paired, hp_ratio, bio_null, theory
    dep = np.array([
        # paired, hp_ratio, bio_null, theory
        [0, 0, 2, 1],   # E1: genome-wide HP1/HP2 read count; null is 50:50 biology
        [2, 2, 0, 0],   # E2: cross-mode; uses Paired as reference and HP_Ratio metric
        [2, 2, 0, 0],   # E3: 62% artifact; uses Paired HP_Ratio 0.4-0.6 as baseline
        [0, 1, 0, 0],   # E4: within-TO CV-2; independent of Paired
        [0, 1, 1, 2],   # E5: PON-only perturbation; mainly theory-driven
        [0, 0, 0, 2],   # E6: pure math derivation
    ])

    color_map = {0: "#81C784", 1: "#FFD54F", 2: "#E57373"}
    label_map = {0: "INDEPENDENT", 1: "PARTIAL", 2: "DEPENDENT"}

    n_rows = len(evidence_lines)
    n_cols = len(dimensions)

    # Draw cells
    for i in range(n_rows):
        for j in range(n_cols):
            code = dep[i, j]
            ax.add_patch(Rectangle((j + 4, n_rows - 1 - i), 1, 1,
                                   facecolor=color_map[code],
                                   edgecolor="black", linewidth=0.8))
            ax.text(j + 4.5, n_rows - 1 - i + 0.5, label_map[code],
                    ha="center", va="center", fontsize=8.5,
                    fontweight="bold", color="#263238")

    # Column headers (dimensions)
    for j, (name, _) in enumerate(dimensions):
        ax.text(j + 4.5, n_rows + 0.5, name,
                ha="center", va="center", fontsize=10,
                fontweight="bold", color="#263238")

    # Row labels (evidence)
    for i, (title, method) in enumerate(evidence_lines):
        ax.text(3.9, n_rows - 1 - i + 0.65, title,
                ha="right", va="center", fontsize=10.5, fontweight="bold")
        ax.text(3.9, n_rows - 1 - i + 0.25, f"({method})",
                ha="right", va="center", fontsize=8.5, style="italic", color="#616161")

    # Right-side verdict column
    verdicts = [
        "Biology-null based",
        "Paired-triangulated",
        "Paired-triangulated",
        "TO-internal",
        "Experimental",
        "Theoretical",
    ]
    for i, v in enumerate(verdicts):
        vcolor = ("#E57373" if "Paired" in v
                  else "#81C784" if ("Theoretical" in v or "TO-internal" in v or "Experimental" in v)
                  else "#FFD54F")
        ax.add_patch(Rectangle((n_cols + 4, n_rows - 1 - i), 2.5, 1,
                               facecolor=vcolor, edgecolor="black", linewidth=0.8))
        ax.text(n_cols + 4 + 1.25, n_rows - 1 - i + 0.5, v,
                ha="center", va="center", fontsize=9, fontweight="bold")

    ax.text(n_cols + 4 + 1.25, n_rows + 0.5, "Evidence\nlineage",
            ha="center", va="center", fontsize=10,
            fontweight="bold", color="#263238")

    # Summary box
    ax.text(6.75, -2,
            "KEY INSIGHT\n"
            "─────────────────────────────────────────────────\n"
            "4/6 evidence lines are INDEPENDENT of Paired-as-reference:\n"
            "  • E1 (17.3:1 bias)  → biological null 50:50 expectation\n"
            "  • E4 (7/7 consistency) → purely TO-internal replication\n"
            "  • E5 (PON-only experiment) → perturbation/causal\n"
            "  • E6 (phasing graph math) → theoretical derivation\n"
            "\n"
            "Only E2 and E3 use Paired HP_Ratio as ground truth.\n"
            "➜ Even if Paired had calibration issues, the core conclusion\n"
            "  (self-phasing exists in TO mode) remains supported by 4 of 6 lines.",
            ha="center", va="center", fontsize=10, family="monospace",
            bbox=dict(facecolor="#FFFDE7", edgecolor="#F57F17",
                      boxstyle="round,pad=0.8", linewidth=2))

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor="#81C784", edgecolor="black",
                       label="INDEPENDENT of this dimension"),
        mpatches.Patch(facecolor="#FFD54F", edgecolor="black",
                       label="PARTIAL (uses the metric but not as sole basis)"),
        mpatches.Patch(facecolor="#E57373", edgecolor="black",
                       label="DEPENDENT on this dimension"),
    ]
    ax.legend(handles=legend_patches, loc="lower left",
              bbox_to_anchor=(0, -0.4), ncol=1, fontsize=9, frameon=True)

    ax.set_xlim(-4, n_cols + 8)
    ax.set_ylim(-4.5, n_rows + 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Figure 15 — Evidence Chain Dependency Matrix: "
                 "Is Self-Phasing Conclusion Reliant on Paired HP1/HP2 Balance?",
                 fontsize=13, fontweight="bold", pad=15)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig15_evidence_dependency_matrix.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 16 — V5 three-layer getVote() logic + code excerpt
# =============================================================================
def fig16_v5_threelayer_logic():
    """Visualise V5 three-layer logic + minimal code excerpt + scenario table.

    Demonstrates the key change vs V3-Fixed: addition of Layer 1.5
    (somatic fallback) when germline votes = 0.
    """
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 2], hspace=0.3, wspace=0.2)

    # ---- Top-left: Three-layer flow diagram ----
    axL = fig.add_subplot(gs[0, 0])
    axL.set_xlim(0, 10)
    axL.set_ylim(0, 10)
    axL.set_title("V5 getVote(): Three-Layer Determination Logic",
                  fontsize=13, fontweight="bold", color="#1565C0")

    # Input box
    box(axL, 3, 8.5, 4, 1.0,
        "INPUT: countMap[HAPLOTYPE1/2/1_1/2_1/3]",
        facecolor="#E3F2FD", edgecolor="#1565C0", fontsize=9.5, fontweight="bold")

    # Layer 1
    axL.annotate("", xy=(5, 7.7), xytext=(5, 8.5),
                 arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
    box(axL, 1.5, 6.5, 7, 1.2,
        "Layer 1 — GERMLINE PRIORITY\n"
        "if (germlineHP1 > 0 || germlineHP2 > 0)\n"
        "→ min/max from germline votes; germlineResult = 1 or 2",
        facecolor="#C8E6C9", edgecolor="#2E7D32", fontsize=9, fontweight="bold")

    # Layer 1.5 (V5 NEW)
    axL.annotate("", xy=(5, 5.7), xytext=(5, 6.5),
                 arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
    box(axL, 1.5, 4.5, 7, 1.2,
        "★ Layer 1.5 — SOMATIC FALLBACK (V5 NEW) ★\n"
        "else if (somaticHP1_1 > 0 || somaticHP2_1 > 0)\n"
        "→ min/max from somatic directional votes; germlineResult = 1 or 2",
        facecolor="#FFE0B2", edgecolor="#E65100", fontsize=9, fontweight="bold")

    # Layer 2
    axL.annotate("", xy=(5, 3.7), xytext=(5, 4.5),
                 arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
    box(axL, 1.5, 2.5, 7, 1.2,
        "Layer 2 — SOMATIC ANNOTATION (encoding)\n"
        "germlineResult = 1 + somatic → HP:i:11\n"
        "germlineResult = 2 + somatic → HP:i:21\n"
        "germlineResult = 0 + somatic-only → HP:i:33 (true ambiguous)",
        facecolor="#F8BBD0", edgecolor="#AD1457", fontsize=9, fontweight="bold")

    # Confidence threshold
    axL.annotate("", xy=(5, 1.7), xytext=(5, 2.5),
                 arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
    box(axL, 1.5, 0.5, 7, 1.0,
        "judgeHaplotype() threshold: max/(max+min) ≥ 0.6\n"
        "→ split votes intercepted as HP:i:33",
        facecolor="#E0E0E0", edgecolor="#424242", fontsize=8.5)

    axL.axis("off")

    # ---- Top-right: Code excerpt ----
    axR = fig.add_subplot(gs[0, 1])
    axR.set_xlim(0, 10)
    axR.set_ylim(0, 10)
    axR.set_title("HaplotagProcess.cpp:512-563 (V5)",
                  fontsize=13, fontweight="bold", color="#1565C0")

    code_excerpt = (
        "// Layer 1: Germline priority\n"
        "if (germlineHP1 > 0 || germlineHP2 > 0) {\n"
        "    if (germlineHP1 >= germlineHP2) {\n"
        "        min = germlineHP2; max = germlineHP1;\n"
        "        germlineResult = 1;\n"
        "    } else { ... germlineResult = 2; }\n"
        "}\n"
        "// Layer 1.5: Somatic fallback (V5 NEW)\n"
        "else if (somaticHP1 > 0 || somaticHP2 > 0) {\n"
        "    if (somaticHP1 >= somaticHP2) {\n"
        "        min = somaticHP2; max = somaticHP1;\n"
        "        germlineResult = 1;\n"
        "    } else { ... germlineResult = 2; }\n"
        "}\n"
        "// Layer 2: Tag encoding\n"
        "if (somaticTotal > 0) {\n"
        "    if (germlineResult == 1)       hpResult = 11;\n"
        "    else if (germlineResult == 2)  hpResult = 21;\n"
        "    else                           hpResult = 33;\n"
        "}"
    )
    axR.text(0.3, 9.7, code_excerpt, fontsize=9, family="monospace",
             verticalalignment="top",
             bbox=dict(facecolor="#FAFAFA", edgecolor="#1565C0",
                       boxstyle="round,pad=0.5", linewidth=1.5))
    axR.axis("off")

    # ---- Bottom: Scenario table ----
    axB = fig.add_subplot(gs[1, :])
    axB.set_xlim(0, 12)
    axB.set_ylim(0, 8)
    axB.set_title("Five Scenarios: V3-Fixed vs V5 Behavior — only Scenario C differs",
                  fontsize=12, fontweight="bold")

    headers = ["Scenario", "Germline votes", "Somatic votes",
               "V3-Fixed result", "V5 result", "Difference"]
    rows = [
        ("A: germline majority",  "HP1=20, HP2=3",  "HP1_1=2",                  "11", "11", "same"),
        ("B: germline minority",  "HP1=20, HP2=3",  "HP2_1=5",                  "11", "11", "same"),
        ("★ C: germline absent",  "HP1=0, HP2=0",   "HP2_1=5, HP1_1=1",         "33", "21", "V5 recovers direction"),
        ("D: pure HP3",           "HP1=0, HP2=0",   "HP3=3",                    "33", "33", "same"),
        ("E: no somatic",         "HP1=0, HP2=0",   "—",                        "0",  "0",  "same"),
    ]
    col_x = [0.2, 2.6, 5.0, 7.4, 8.6, 9.8]
    # Header
    for cx, h in zip(col_x, headers):
        axB.text(cx, 7.3, h, fontsize=10, fontweight="bold")
    axB.plot([0.1, 11.9], [7.0, 7.0], "-", color="black", lw=1)

    for i, row in enumerate(rows):
        y = 6.0 - i * 1.1
        is_critical = row[0].startswith("★")
        bg = "#FFF8E1" if is_critical else "white"
        if is_critical:
            axB.add_patch(Rectangle((0.05, y - 0.45), 11.9, 1.0,
                                    facecolor=bg, edgecolor="#E65100",
                                    linewidth=1.5, zorder=0))
        for cx, val in zip(col_x, row):
            weight = "bold" if is_critical else "normal"
            axB.text(cx, y, val, fontsize=9, fontweight=weight,
                     verticalalignment="center")

    axB.text(6, 0.4,
             "Scenario C is the V5 KEY UPGRADE: previously discarded somatic directional "
             "evidence (HP1_1/HP2_1) is now recovered as fallback when germline = 0.",
             ha="center", va="center", fontsize=10, fontweight="bold",
             color="#E65100",
             bbox=dict(facecolor="#FFF3E0", edgecolor="#E65100",
                       boxstyle="round,pad=0.4"))
    axB.axis("off")

    fig.suptitle("Figure 16 — V5 Somatic Fallback: Three-Layer getVote() Logic + "
                 "Five-Scenario Behavior Comparison",
                 fontsize=14, fontweight="bold", y=0.995)

    fig.savefig(OUT_DIR / "fig16_v5_threelayer_logic.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 17 — HP tag distribution across 5 versions (stacked bar + delta)
# =============================================================================
def fig17_hp_tag_5versions():
    """Two-panel comparison:
    Left: chr1:1-10M sample (showing per-version HP tag composition)
    Right: full-genome V3-F vs V5 (showing HP:i:33 collapse)
    """
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={"width_ratios": [3, 2]})

    # ---- Left: chr1:1-10M sample, 4 versions ----
    versions_L = ["Baseline", "V2b\nPON-Only", "V3-Fixed", "V5\n(current)"]
    # Data from session report Section 7.3
    germ1   = [26126, 26206, 26206, 26206]
    germ2   = [13598, 13682, 13682, 13682]
    hp11    = [1357,  1410,  1158,  1359]
    hp21    = [1256,  1203,  1070,  1243]
    hp33    = [0,     0,     385,   11]

    x = np.arange(len(versions_L))
    width = 0.6

    axes[0].bar(x, germ1, width, label="HP:i:1 (germline HP1)", color="#2E7D32", alpha=0.85)
    axes[0].bar(x, germ2, width, bottom=germ1, label="HP:i:2 (germline HP2)", color="#66BB6A", alpha=0.85)
    bottom_3 = np.array(germ1) + np.array(germ2)
    axes[0].bar(x, hp11, width, bottom=bottom_3, label="HP:i:11 (somatic HP1)", color="#1565C0", alpha=0.85)
    bottom_4 = bottom_3 + np.array(hp11)
    axes[0].bar(x, hp21, width, bottom=bottom_4, label="HP:i:21 (somatic HP2)", color="#42A5F5", alpha=0.85)
    bottom_5 = bottom_4 + np.array(hp21)
    axes[0].bar(x, hp33, width, bottom=bottom_5, label="HP:i:33 (ambiguous)", color="#E65100", alpha=0.95)

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(versions_L, fontsize=10)
    axes[0].set_ylabel("Read count", fontsize=11)
    axes[0].set_title("A. chr1:1-10M sample — 4 versions\n"
                      "(Baseline & V2b: HP:i:33 = 0 due to enum/integer bug)",
                      fontsize=11, fontweight="bold")
    axes[0].legend(loc="upper right", fontsize=8.5, framealpha=0.95)
    # Annotate HP:i:33 emergence
    axes[0].annotate("V3-F: HP:i:33\nfinally appears\n(but over-marked)",
                     xy=(2, bottom_5[2] + hp33[2] / 2),
                     xytext=(2.5, 38000),
                     fontsize=8, color="#E65100", fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="#E65100"))
    axes[0].annotate("V5: HP:i:33\ncollapses to 11\n(directional recovered)",
                     xy=(3, bottom_5[3] + hp33[3] / 2),
                     xytext=(3, 38000),
                     fontsize=8, color="#1565C0", fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="#1565C0"))

    # ---- Right: Full-genome V3-F vs V5 (only somatic tags) ----
    versions_R = ["V3-Fixed", "V5"]
    full_hp11  = [584117,  666997]
    full_hp21  = [547118,  593720]
    full_hp33  = [239679,  110197]

    x2 = np.arange(len(versions_R))
    width2 = 0.5
    axes[1].bar(x2, full_hp11, width2, label="HP:i:11", color="#1565C0", alpha=0.85)
    axes[1].bar(x2, full_hp21, width2, bottom=full_hp11, label="HP:i:21", color="#42A5F5", alpha=0.85)
    bot = np.array(full_hp11) + np.array(full_hp21)
    axes[1].bar(x2, full_hp33, width2, bottom=bot, label="HP:i:33 (amb)", color="#E65100", alpha=0.95)

    # Delta annotations
    deltas = ["V3-F", "V5\nΔHP:i:11=+82,880\nΔHP:i:21=+46,602\nΔHP:i:33=-129,482 (-54%)"]
    for i, d in enumerate(deltas):
        total = full_hp11[i] + full_hp21[i] + full_hp33[i]
        axes[1].text(i, total + 60000, d, ha="center", fontsize=9,
                     fontweight="bold" if i == 1 else "normal",
                     color="#1565C0" if i == 1 else "black")

    axes[1].set_xticks(x2)
    axes[1].set_xticklabels(versions_R, fontsize=10)
    axes[1].set_ylabel("Read count (full genome, somatic tags)", fontsize=11)
    axes[1].set_title("B. Full-genome — V3-F vs V5 somatic tags\n"
                      "129K reads recovered from HP:i:33 → HP:i:11/21",
                      fontsize=11, fontweight="bold")
    axes[1].legend(loc="upper left", fontsize=9, framealpha=0.95)
    axes[1].set_ylim(0, max(np.array(full_hp11) + np.array(full_hp21) + np.array(full_hp33)) * 1.30)

    fig.suptitle("Figure 17 — HP Tag Distribution Across Versions: "
                 "Bug Discovery (V3-F) → Information Recovery (V5)",
                 fontsize=14, fontweight="bold", y=0.99)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig17_hp_tag_5versions.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 18 — Concordance + AMB% + F1 evolution across versions
# =============================================================================
def fig18_concordance_amb_f1():
    """Three-panel comparison of quality metrics across versions.

    Panel A: Clean-PS Concordance (Baseline / V5 / Paired ground truth)
    Panel B: AMB% closeness to Paired
    Panel C: SEQC2 F1 (essentially flat — confirms F1 doesn't measure tag quality)
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel A: Clean-PS Somatic + Germline Concordance
    versions = ["Baseline", "V5", "Paired\n(ground truth)"]
    somatic_acc = [82.2, 90.5, 100.0]   # % from session report §9.5
    germline_acc = [87.2, 91.7, 100.0]
    paired_target = 100

    x = np.arange(len(versions))
    width = 0.35
    bars_s = axes[0].bar(x - width/2, somatic_acc, width,
                         label="Somatic", color="#1565C0", alpha=0.85)
    bars_g = axes[0].bar(x + width/2, germline_acc, width,
                         label="Germline", color="#2E7D32", alpha=0.85)
    axes[0].axhline(100, ls="--", color="grey", alpha=0.5)

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(versions, fontsize=10)
    axes[0].set_ylabel("Clean-PS Concordance (%)", fontsize=11)
    axes[0].set_title("A. Concordance vs Paired ground truth\n"
                      "(higher = closer to truth)",
                      fontsize=11, fontweight="bold")
    axes[0].legend(loc="lower right", fontsize=9)
    axes[0].set_ylim(75, 105)

    for bars in [bars_s, bars_g]:
        for bar in bars:
            h = bar.get_height()
            axes[0].text(bar.get_x() + bar.get_width()/2, h + 0.5,
                         f"{h:.1f}%", ha="center", fontsize=8.5)

    # Annotate gap
    axes[0].annotate("Baseline gap = 17.8pp\n(somatic)", xy=(0, 82.2),
                     xytext=(-0.1, 76),
                     fontsize=8, color="#C62828",
                     arrowprops=dict(arrowstyle="->", color="#C62828"))
    axes[0].annotate("V5 gap = 9.5pp", xy=(1, 90.5),
                     xytext=(1.1, 95),
                     fontsize=8, color="#1565C0", fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="#1565C0"))

    # Panel B: AMB% closeness to Paired
    versions_B = ["Baseline", "V3-Fixed", "V5", "Paired\n(target)"]
    amb_global = [1.3, 17.5, 8.0, 5.4]   # global AMB%
    amb_clean  = [None, None, 4.7, 5.4]   # clean-PS AMB%

    x_B = np.arange(len(versions_B))
    bars_amb = axes[1].bar(x_B, amb_global,
                           color=["#C62828", "#FFA726", "#1565C0", "#2E7D32"],
                           alpha=0.85)
    axes[1].axhline(5.4, ls="--", color="#2E7D32", alpha=0.6,
                    label="Paired target (5.4%)")
    axes[1].set_xticks(x_B)
    axes[1].set_xticklabels(versions_B, fontsize=10)
    axes[1].set_ylabel("AMB% (HP:i:33 / somatic total)", fontsize=11)
    axes[1].set_title("B. AMB% Honesty\n"
                      "(closer to Paired = more truthful)",
                      fontsize=11, fontweight="bold")
    axes[1].legend(loc="upper right", fontsize=9)

    for bar, v in zip(bars_amb, amb_global):
        h = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2, h + 0.4,
                     f"{v}%", ha="center", fontsize=9, fontweight="bold")

    # Annotate Baseline = false-low
    axes[1].annotate("Baseline 1.3%:\nFALSE-LOW\n(bug forces direction)",
                     xy=(0, 1.3), xytext=(0.3, 12),
                     fontsize=8.5, color="#C62828", fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="#C62828"))
    axes[1].annotate("V3-F 17.5%:\nover-conservative", xy=(1, 17.5),
                     xytext=(1.3, 16),
                     fontsize=8.5, color="#FFA726", fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="#FFA726"))
    axes[1].annotate("V5 8.0%:\nbalanced", xy=(2, 8.0),
                     xytext=(2.4, 13),
                     fontsize=8.5, color="#1565C0", fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="#1565C0"))

    # Panel C: SEQC2 F1 evolution
    versions_C = ["Baseline", "V2b", "V3-Fixed", "V5"]
    f1_seqc2 = [0.7157, 0.7155, 0.7154, 0.7154]
    f1_ism   = [0.0154, 0.0177, 0.0124, 0.0125]

    x_C = np.arange(len(versions_C))
    color_seqc2 = ["#FFB74D"] * 4
    color_seqc2[3] = "#1565C0"  # highlight V5

    axes[2].bar(x_C, f1_seqc2, color=color_seqc2, alpha=0.85)
    axes[2].set_xticks(x_C)
    axes[2].set_xticklabels(versions_C, fontsize=10)
    axes[2].set_ylabel("SEQC2 F1", fontsize=11)
    axes[2].set_title("C. SEQC2 F1 (calling pipeline)\n"
                      "Flat range 0.7154-0.7157 — F1 doesn't measure tag quality",
                      fontsize=11, fontweight="bold")
    axes[2].set_ylim(0.710, 0.720)

    for i, v in enumerate(f1_seqc2):
        axes[2].text(i, v + 0.0003, f"{v:.4f}", ha="center", fontsize=9,
                     fontweight="bold" if i == 3 else "normal")

    # Bottom annotation strip
    axes[2].text(1.5, 0.7115,
                 "ΔF1 (V5 vs Baseline) = -0.0003 (noise)\n"
                 "All quality improvements invisible in F1",
                 ha="center", va="center", fontsize=8.5, style="italic",
                 color="#424242",
                 bbox=dict(facecolor="#FFFDE7", edgecolor="#F57F17",
                           boxstyle="round,pad=0.3"))

    fig.suptitle("Figure 18 — Quality Metrics Evolution: Tag-quality (A,B) Improves while F1 (C) is Flat",
                 fontsize=14, fontweight="bold", y=1.01)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig18_concordance_amb_f1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Figure 19 — Conclusion impact matrix (V5 effect on existing research)
# =============================================================================
def fig19_conclusion_impact_matrix():
    """Heatmap mapping 10 existing conclusions to V5 impact tiers."""
    fig, ax = plt.subplots(figsize=(16, 10))

    conclusions = [
        ("CL-001",  "Self-phasing 17.3:1 bias",                    "Tier 1"),
        ("CL-009",  "62% TO LOH = artifact",                       "Tier 1"),
        ("CL-010",  "HPFineNGroups 89.1% TP rate (subclone marker)","Tier 1"),
        ("CL-011",  "LOH-constrained phasing ΔNG = +0.385",         "Tier 2"),
        ("CL-006",  "LOH.bed AUC = 0.70",                            "Tier 2"),
        ("CL-007",  "Germline ASM >> somatic (3-6×)",               "Tier 3"),
        ("CL-025",  "LOH.bed Jaccard = 1.0 (PON-only)",             "Tier 3"),
        ("CL-023",  "Coverage_Multiple confound",                   "Tier 3"),
        ("E6",      "Phasing graph math (theory)",                  "Tier 3"),
        ("E4",      "7/7 cross-sample direction",                   "Tier 3"),
    ]

    dims = [
        ("HP-tag\ndependent",      "hp_dep"),
        ("Mechanism\nrobust",       "mech"),
        ("Quantitative\nshift on V5", "shift"),
        ("Re-eval\npriority",       "prio"),
    ]

    # Score matrix: 0=no/low, 1=partial/medium, 2=high
    # Tier 1: hp_dep=2, mech=1 (mostly), shift=2, prio=2
    # Tier 2: hp_dep=1-2, mech=2, shift=1, prio=1
    # Tier 3: hp_dep=0, mech=2, shift=0, prio=0
    scores = np.array([
        [2, 1, 2, 2],   # CL-001
        [2, 1, 2, 2],   # CL-009
        [2, 1, 2, 2],   # CL-010 highest priority
        [1, 2, 1, 1],   # CL-011
        [1, 2, 0, 1],   # CL-006
        [0, 2, 0, 0],   # CL-007
        [0, 2, 0, 0],   # CL-025
        [0, 2, 0, 0],   # CL-023
        [0, 2, 0, 0],   # E6
        [0, 2, 0, 0],   # E4
    ])

    # Color: rows: Tier 1 red-tinged; Tier 2 amber; Tier 3 green
    tier_colors = {"Tier 1": "#FFCDD2", "Tier 2": "#FFE0B2", "Tier 3": "#C8E6C9"}
    cell_colors = {0: "#E0E0E0", 1: "#FFD54F", 2: "#EF5350"}
    # For mechanism robust, invert: 2=green=robust, 1=partial, 0=concern
    mech_colors = {0: "#EF5350", 1: "#FFD54F", 2: "#81C784"}

    n_rows = len(conclusions)
    n_cols = len(dims)

    # Tier indicator column (left side)
    for i, (cid, desc, tier) in enumerate(conclusions):
        ax.add_patch(Rectangle((-0.5, n_rows - 1 - i), 0.5, 1,
                               facecolor=tier_colors[tier],
                               edgecolor="black", linewidth=0.6))
        ax.text(-0.25, n_rows - 1 - i + 0.5, tier.replace("Tier ", "T"),
                ha="center", va="center", fontsize=9, fontweight="bold")

    # Cells
    for i in range(n_rows):
        for j in range(n_cols):
            v = scores[i, j]
            # Use mech_colors for column 1 ("mech robust"), else cell_colors
            color = mech_colors[v] if j == 1 else cell_colors[v]
            label = ["LOW", "MED", "HIGH"][v] if j != 1 else \
                    ["WEAK", "PARTIAL", "ROBUST"][v]
            ax.add_patch(Rectangle((j, n_rows - 1 - i), 1, 1,
                                   facecolor=color, edgecolor="black", linewidth=0.6))
            ax.text(j + 0.5, n_rows - 1 - i + 0.5, label,
                    ha="center", va="center", fontsize=8.5, fontweight="bold")

    # Column headers
    for j, (name, _) in enumerate(dims):
        ax.text(j + 0.5, n_rows + 0.4, name,
                ha="center", va="center", fontsize=9.5, fontweight="bold")
    ax.text(-0.25, n_rows + 0.4, "Tier",
            ha="center", va="center", fontsize=9.5, fontweight="bold")

    # Row labels (right side)
    for i, (cid, desc, tier) in enumerate(conclusions):
        ax.text(n_cols + 0.2, n_rows - 1 - i + 0.65,
                f"{cid}", fontsize=10, fontweight="bold")
        ax.text(n_cols + 0.2, n_rows - 1 - i + 0.25,
                desc, fontsize=8.5, color="#424242")

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor="#FFCDD2", edgecolor="black",
                       label="Tier 1: needs re-evaluation under V5"),
        mpatches.Patch(facecolor="#FFE0B2", edgecolor="black",
                       label="Tier 2: possibly shifted, mechanism robust"),
        mpatches.Patch(facecolor="#C8E6C9", edgecolor="black",
                       label="Tier 3: unaffected"),
    ]
    ax.legend(handles=legend_patches, loc="lower center",
              bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=10,
              frameon=True)

    # Summary
    ax.text(n_cols / 2 - 0.25, -1.5,
            "Conclusion stability (V5 impact):  3/10 Tier 1 (re-eval)  |  "
            "2/10 Tier 2 (monitor)  |  5/10 Tier 3 (unaffected)\n"
            "All Tier 1 conclusions retain mechanism but require quantitative re-computation under V5.",
            ha="center", va="center", fontsize=10, fontweight="bold",
            bbox=dict(facecolor="#FFFDE7", edgecolor="#F57F17",
                      boxstyle="round,pad=0.5", linewidth=1.5))

    ax.set_xlim(-1, n_cols + 8)
    ax.set_ylim(-2.5, n_rows + 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Figure 19 — V5 Impact on 10 Existing Conclusions: "
                 "Tier-1 needs re-eval, Tier-3 unaffected",
                 fontsize=13, fontweight="bold", pad=15)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig19_conclusion_impact_matrix.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main():
    print("[1/11] Pipeline comparison...")
    fig1_pipeline_comparison()
    print("[2/11] AF=0.3 walkthrough...")
    fig3_af03_walkthrough()
    print("[3/11] Evidence summary...")
    fig4_evidence_summary()
    print("[4/11] Impact matrix...")
    fig5_impact_matrix()
    print("[5/11] Fix decision tree...")
    fig6_decision_tree()
    print("[6/11] NGroups collapse (reads Phase 1 CSVs)...")
    fig7_ngroups_collapse()
    print("[7/11] Phase 1 AUC comparison...")
    auc_data = fig8_auc_comparison()
    print("[8/11] HP_Ratio violin...")
    fig9_hp_ratio_violin()
    print("[9/11] Next-step tree...")
    fig10_next_step_tree()
    print("[10/12] Three-track × two-flag matrix (NEW)...")
    fig11_three_track_matrix()
    print("[11/12] Adversarial review tree (NEW)...")
    fig12_adversarial_tree()
    print("[12/12] End-to-end case storyboard (NEW)...")
    fig13_case_storyboard()
    print("\nAll figures written to:", OUT_DIR)
    print("\nAUC summary (flag off → flag on):")
    for name, a_off, a_on in auc_data:
        print(f"  {name:25s}  {a_off:.4f} -> {a_on:.4f}  Δ={a_on - a_off:+.4f}")


if __name__ == "__main__":
    main()
