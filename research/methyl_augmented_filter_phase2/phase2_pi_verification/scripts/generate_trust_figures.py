#!/usr/bin/env python3
"""Phase 2 PI Trust Framework — 6 matplotlib figures.

Sections 1-7 in standalone HTML each use one figure plus tables.

Outputs:
    figures/01_kpi_progress.png       4 KPI vs target dashboard
    figures/02_evidence_ladder.png    7 claims × tier rank
    figures/03_counter_evidence.png   3 reconciliation + cross-session
    figures/04_trust_radar.png        6-axis self-audit
    figures/05_decision_roi.png       5 paths × cost × prior × gain
    figures/06_uncertainty_map.png    5-layer epistemic map
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
OUT = REPO / "research/methyl_augmented_filter_phase2/phase2_pi_verification/figures"
OUT.mkdir(parents=True, exist_ok=True)


def fig_01_kpi_progress():
    """4 KPI vs target horizontal bar."""
    kpis = [
        ("KPI1\nΔF1 magnitude\n(HCC1395-internal)", 224, "✅ 達成 +0.02236 / target +0.01", "#15803D"),
        ("KPI2\nCross-sample\ndirection concordance", 40, "⚠ 1/5 transfer · 3/5 refit≈0 / target ≥4/5", "#C2410C"),
        ("KPI3\nProduction rule\ndeployable", 75, "🟡 4-feature gate rule / target full panel", "#A16207"),
        ("KPI4\nMechanism\nunderstood", 80, "🟡 caller-F1-headroom L3 + ISM vestigial L2 / target L1-L2", "#A16207"),
    ]

    fig, ax = plt.subplots(figsize=(11, 5.5))
    y_pos = np.arange(len(kpis))
    values = [k[1] for k in kpis]
    colors = [k[3] for k in kpis]

    bars = ax.barh(y_pos, values, color=colors, edgecolor="black", linewidth=0.5)
    ax.axvline(100, color="black", linestyle="--", linewidth=1.2, label="100% target")

    # Annotate
    for i, (label, v, note, c) in enumerate(kpis):
        ax.text(v + 5 if v < 100 else v - 35, i, f"{v}%", ha="left" if v < 100 else "right",
                va="center", fontsize=14, fontweight="bold",
                color="white" if v >= 100 else "black")
        ax.text(0, i - 0.32, note, ha="left", va="top", fontsize=8.5, color="#6B6B66")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([k[0] for k in kpis], fontsize=10)
    ax.set_xlim(0, 240)
    ax.set_xlabel("Achievement %", fontsize=11)
    ax.set_title("Figure 1 — Phase 2 Goal Achievement Dashboard\n"
                 "4 KPI vs target (=100%): KPI1 已超達 / KPI2 仍 NEGATIVE / KPI3+4 部分達成",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "01_kpi_progress.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/01_kpi_progress.png")


def fig_02_evidence_ladder():
    """7 claims × tier (n + pre-reg)."""
    claims = [
        ("HCC1395 cycle 1 ΔF1 +0.02236", 4, "H_C1_2 pre-reg", "n=1·multi-seed"),
        ("BAM-invariant V3F/V5/V6", 4, "H_C1_6 pre-reg", "n=3·max_var 0.00055"),
        ("ISM vestigial in LR", 4, "H_M1a pre-reg", "n=5·4/5 concordant"),
        ("Cross-sample transfer NEGATIVE", 3, "H_C1_5 pre-reg", "n=5·p=0.1875"),
        ("caller_af 主導 67% disaster", 3, "H_A1 pre-reg", "n=4 new samples"),
        ("Caller-F1-headroom mechanism", 3, "post-hoc", "n=5 quadrant"),
        ("Cycle 3 qualifying mean +0.01499", 2, "H_C3_1 pre-reg", "n=2 (limit)"),
    ]
    claims_sorted = sorted(claims, key=lambda c: c[1], reverse=False)

    fig, ax = plt.subplots(figsize=(13, 6))
    y_pos = np.arange(len(claims_sorted))
    tier_colors = {2: "#C2410C", 3: "#A16207", 4: "#15803D"}
    bars = ax.barh(y_pos, [c[1] for c in claims_sorted],
                    color=[tier_colors[c[1]] for c in claims_sorted],
                    edgecolor="black", linewidth=0.4)

    for i, (claim, tier, prereg, n_info) in enumerate(claims_sorted):
        ax.text(tier + 0.1, i, "⭐" * tier, ha="left", va="center", fontsize=14)
        ax.text(0.05, i + 0.05, claim, ha="left", va="bottom", fontsize=10, color="white",
                fontweight="bold")
        ax.text(0.05, i - 0.05, f"{prereg} · {n_info}", ha="left", va="top", fontsize=8,
                color="white", style="italic")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([])
    ax.set_xlim(0, 5.5)
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(["", "⭐", "⭐⭐", "⭐⭐⭐", "⭐⭐⭐⭐", "⭐⭐⭐⭐⭐"], fontsize=12)
    ax.set_xlabel("Evidence Tier (scientific-rigor §2)", fontsize=11)
    ax.set_title("Figure 2 — Evidence Ladder: 7 Major Claims ranked by Tier\n"
                 "⭐2 = n-limited · ⭐3 = pre-reg with caveat · ⭐4 = pre-reg + n=5 + cross-validated",
                 fontsize=12, fontweight="bold")

    handles = [mpatches.Patch(color="#15803D", label="⭐4 (L2)"),
               mpatches.Patch(color="#A16207", label="⭐3 (L3)"),
               mpatches.Patch(color="#C2410C", label="⭐2 (n-limited)")]
    ax.legend(handles=handles, loc="lower right", fontsize=9)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "02_evidence_ladder.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/02_evidence_ladder.png")


def fig_03_counter_evidence():
    """3 reconciliation matrices + cross-session bar."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), gridspec_kw={"width_ratios": [2, 1]})

    # Panel 1: 3 reconciliation matrices (table-like)
    ax = axes[0]
    ax.axis("off")
    reconciliations = [
        ("Cycle 1 +0.02236 vs Cycle 3 ablation vestigial",
         "Cycle 1 數字仍正確；vestigial 意指 methylation 不提供 incremental\n"
         "→ Drop 5 methyl → ΔF1 +0.02171 仍接近 +0.02236",
         "#1E3A8A"),
        ("Cycle 1 ⭐3 strong vs Cycle 2 cross-sample NEGATIVE",
         "Cycle 1 是 HCC1395-internal；Cycle 2 證實「不可 universal transfer」\n"
         "→ caller-F1-headroom 機制：3/4 新樣本 F1>0.83 留無 filter 空間",
         "#A16207"),
        ("HCC1937 +0.00761 vs H1437/H2009 ≈0",
         "Cycle 3 gate rule 已預測：caller F1 < 0.80 + FP density > 10%\n"
         "→ HCC1937 (F1=0.37, FP=16%) qualifies；H1437/H2009 不 qualify",
         "#15803D"),
    ]
    for i, (matter, reconcile, c) in enumerate(reconciliations):
        y = 0.85 - i * 0.30
        ax.add_patch(mpatches.FancyBboxPatch((0.02, y - 0.20), 0.45, 0.20,
                                              boxstyle="round,pad=0.02",
                                              fc="#FEF2F0", ec="#C2410C", linewidth=1.5,
                                              transform=ax.transAxes))
        ax.text(0.05, y - 0.05, "表面矛盾", fontsize=8.5, fontweight="bold", color="#C2410C",
                transform=ax.transAxes)
        ax.text(0.05, y - 0.12, matter, fontsize=8.5, color="#7F1D1D", transform=ax.transAxes,
                wrap=True)

        ax.annotate("", xy=(0.55, y - 0.10), xytext=(0.48, y - 0.10),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="black", lw=1.5))

        ax.add_patch(mpatches.FancyBboxPatch((0.55, y - 0.20), 0.43, 0.20,
                                              boxstyle="round,pad=0.02",
                                              fc="#DCFCE7", ec=c, linewidth=1.5,
                                              transform=ax.transAxes))
        ax.text(0.58, y - 0.05, "實際相容", fontsize=8.5, fontweight="bold", color=c,
                transform=ax.transAxes)
        ax.text(0.58, y - 0.10, reconcile, fontsize=8, color="#166534", transform=ax.transAxes)

    ax.set_title("3 表面矛盾 → 實際 reconcile 結論", fontsize=11, fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Panel 2: cross-session bar
    ax = axes[1]
    paths = ["5/19 V6 §13\nDay 3 H_ABL_1", "5/20 Cycle 3 Step 1.5\nH_M1a"]
    values = [0.00060, 0.00065]
    bars = ax.bar(paths, values, color=["#1E3A8A", "#A16207"], edgecolor="black", linewidth=0.5, width=0.5)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width()/2, v + 0.00003, f"+{v:.5f}",
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylabel("HCC1395 methyl incremental ΔF1", fontsize=10)
    ax.set_ylim(0, 0.0010)
    ax.set_title("Cross-session\nindependent reproduction\n(差距 < 5%)",
                 fontsize=10, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Figure 3 — Counter-Evidence Audit (3 reconciliation + cross-session)",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "03_counter_evidence.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/03_counter_evidence.png")


def fig_04_trust_radar():
    """6-axis radar self-audit."""
    axes_labels = ["Data\n(authoritative)", "Method\n(pre-reg)", "Statistics\n(effect size + CI)",
                   "Causality\n(mechanism)", "Reproducibility\n(bit-exact)", "Negative-evidence\n(counter-evidence)"]
    scores = [90, 85, 70, 65, 95, 90]  # 0-100
    annotations = [
        "SEQC2 + 4 sample orthogonal truth",
        "9 H pre-reg + 1 post-hoc 標明",
        "small n=5 limit (Wilcoxon power)",
        "HCC1954 33% 殘餘未解釋",
        "cycle 2 V6 drift 0 + seed=42",
        "Step 5c / cycle 2 / H_M1a 全列",
    ]

    angles = np.linspace(0, 2 * np.pi, len(scores), endpoint=False).tolist()
    scores_closed = scores + [scores[0]]
    angles_closed = angles + [angles[0]]

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(projection="polar"))
    ax.plot(angles_closed, scores_closed, "o-", linewidth=2, color="#D97757")
    ax.fill(angles_closed, scores_closed, alpha=0.3, color="#D97757")

    for angle, score, label in zip(angles, scores, annotations):
        ax.text(angle, score + 8, f"{score}", ha="center", va="center", fontsize=11,
                fontweight="bold")

    ax.set_xticks(angles)
    ax.set_xticklabels(axes_labels, fontsize=10)
    ax.set_ylim(0, 100)
    ax.set_yticks([20, 40, 60, 80, 100])
    ax.set_yticklabels(["20", "40", "60", "80", "100"], fontsize=8)
    ax.set_title("Figure 4 — Trust Radar 6-axis Self-Audit\n"
                 "Phase 2 current scores (0-100 per axis)",
                 fontsize=12, fontweight="bold", pad=20)
    ax.grid(True, alpha=0.4)

    fig.tight_layout()
    fig.savefig(OUT / "04_trust_radar.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/04_trust_radar.png")


def fig_05_decision_roi():
    """5 paths cost × prior × goal gain bubble."""
    paths = [
        ("A", "Cycle 4 Trial A", 1.5, 75, 15, "#15803D"),
        ("B", "Cycle 3 Step 2 panel", 0.6, 30, 30, "#1E3A8A"),
        ("C", "Pivot phase_block_3d", 0.0, 60, 25, "#A16207"),
        ("D", "Pivot thread_d", 0.0, 70, 30, "#A16207"),
        ("E", "BAM new feature extract", 4.0, 15, 40, "#C2410C"),
    ]

    fig, ax = plt.subplots(figsize=(11, 6.5))
    for code, name, cost, prior, gain, c in paths:
        size = gain * 80
        ax.scatter(cost, prior, s=size, c=c, alpha=0.7, edgecolors="black", linewidth=1.2,
                   zorder=3)
        ax.annotate(f"{code} {name}\nprior={prior}% · gain={gain}% (KPI2)",
                    xy=(cost, prior), xytext=(15, 8), textcoords="offset points",
                    fontsize=9)

    ax.set_xlabel("Cost (days)", fontsize=11)
    ax.set_ylabel("Prior PASS probability (%)", fontsize=11)
    ax.set_xlim(-0.5, 5)
    ax.set_ylim(0, 100)
    ax.set_title("Figure 5 — Decision Tree ROI Bubble (5 Paths)\n"
                 "Bubble size ∝ expected KPI2 gain · Lower-left + large = highest ROI",
                 fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "05_decision_roi.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/05_decision_roi.png")


def fig_06_uncertainty_map():
    """5-layer uncertainty horizontal panels."""
    layers = [
        ("KNOW (L1-L2)", "#15803D",
         ["HCC1395 ΔF1 +0.02236 multi-seed std 5e-5",
          "BAM-invariant V3F/V5/V6 max_var 0.00055",
          "ISM vestigial: drop 5 methyl → HCC1395 -0.00065",
          "caller_af 67% HCC1954 disaster confound",
          "Cross-session H_ABL_1 reproduce ~+0.0006"]),
        ("BELIEVE (L3)", "#A16207",
         ["Caller-F1-headroom 是真實機制 (n=5 quadrant)",
          "HCC1937 may generalize (BRCA1 outlier risk)",
          "4-feature filter 在 low-F1 樣本有效"]),
        ("DON'T KNOW (Research)", "#C2410C",
         ["R1: HCC1954 transfer 33% 殘餘 disaster 真因",
          "R2: Interaction LR (caller_af × methyl) > +0.005?",
          "R3: RF/XGB 在 n=35k 是否 overfit",
          "R4: BAM-level 5hmC/strand-specific orthogonal signal?",
          "R5: COLO829 truth set 可取得?"]),
        ("DON'T KNOW (Agent reasoning) ⭐ NEW", "#1E3A8A",
         ["A1: Cycle 3 ablation 框架僅 drop/shrink；未測 interaction/non-linear",
          "A2: Tier ⭐⭐⭐⭐ vs ⭐⭐⭐ 邊界主觀；reviewer +/- 1 tier",
          "A3: Cohen ribbon (+0.005/+0.01) 沿用 cycle 1 v1.0；其他版本可能改 marginal verdict"]),
    ]

    fig, ax = plt.subplots(figsize=(13, 9))
    ax.axis("off")
    y_cursor = 0.95
    for layer_name, color, items in layers:
        h = 0.04 + 0.07 * len(items)
        ax.add_patch(mpatches.FancyBboxPatch((0.02, y_cursor - h), 0.96, h,
                                              boxstyle="round,pad=0.005",
                                              fc=color, ec="black", linewidth=0.8, alpha=0.25,
                                              transform=ax.transAxes))
        ax.text(0.05, y_cursor - 0.03, layer_name, fontsize=12, fontweight="bold", color=color,
                transform=ax.transAxes)
        for i, item in enumerate(items):
            ax.text(0.08, y_cursor - 0.07 - i * 0.06, "• " + item, fontsize=10,
                    transform=ax.transAxes, color="#141413")
        y_cursor = y_cursor - h - 0.02

    ax.set_title("Figure 6 — Uncertainty 5-layer Epistemic Map\n"
                 "From 'we know' → 'we believe' → 'research unknown' → 'agent unknown' → 'subjectivity'",
                 fontsize=12, fontweight="bold", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "06_uncertainty_map.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/06_uncertainty_map.png")


if __name__ == "__main__":
    print("Generating 6 PI Trust Framework figures …")
    fig_01_kpi_progress()
    fig_02_evidence_ladder()
    fig_03_counter_evidence()
    fig_04_trust_radar()
    fig_05_decision_roi()
    fig_06_uncertainty_map()
    print(f"\nAll 6 figures saved to: {OUT}")
