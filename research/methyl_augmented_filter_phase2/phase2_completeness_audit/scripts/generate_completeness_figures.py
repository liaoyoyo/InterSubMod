#!/usr/bin/env python3
"""Phase 2 Engineering Science Completeness — 2 figures.

Figure 1: Validation pipeline flowchart (7 stages, methods, verdicts)
Figure 2: Completeness gap heatmap (14 hypotheses × 6 audit dimensions)
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures")
OUT.mkdir(parents=True, exist_ok=True)


# ─────────── Figure 1: Validation pipeline flowchart ──────────────────────
def fig_pipeline():
    fig, ax = plt.subplots(figsize=(15, 9))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 11)
    ax.axis("off")

    # 7 stage boxes: (x, y, w, h, label, method, verdict, color)
    stages = [
        (0.3, 8.5, 2.5, 1.5,
         "Stage 0: v1.0 Pilot",
         "Cell-gated filter\n(HCC1395, paired-pileup)",
         "ΔF1 +0.00242 marginal\n(< Cohen 0.005)", "#FCD34D"),
        (3.3, 8.5, 2.5, 1.5,
         "Stage 1: Cycle 1\nGlobal LR",
         "10-feat L2 LR\n5-fold OOF (row-level)",
         "ΔF1 +0.02236\nHCC1395 in-distribution", "#DCFCE7"),
        (6.3, 8.5, 2.5, 1.5,
         "Stage 2: Cycle 2\nCross-sample",
         "Transfer cycle1 coef\n→ 4 new samples",
         "DIRECTION_NEGATIVE\n4/5 ΔF1 < 0", "#FEF2F0"),
        (9.3, 8.5, 2.5, 1.5,
         "Stage 3: Cycle 3\nGate Rule",
         "F1<0.80 + FP>10%\nper-sample re-fit",
         "Qualifying mean\n+0.01499 (n=2)", "#FEF7E6"),
        (12.3, 8.5, 2.5, 1.5,
         "Stage 4: Cycle 3.5\nAblation",
         "Drop methyl block\n× 5 samples refit",
         "ISM vestigial\n+0.00065 incremental", "#FEF2F0"),
        (3.3, 5.0, 2.5, 1.5,
         "Stage 5: LOSO\n(triggered by user)",
         "True sample-level CV\n10-feat",
         "5/5 ΔF1 ≈ 0\nHCC1395 LOSO -0.00012", "#FEF2F0"),
        (6.3, 5.0, 2.5, 1.5,
         "Stage 6: V6 Observation",
         "AUC/Cohen per-feat\n× 5 samples (no LR)",
         "caller_af 0.20-0.92\nLOH/HPFineF 4-5/5 pos", "#DCFCE7"),
        (9.3, 5.0, 2.5, 1.5,
         "Stage 7: H_NEW_2\n2-feat LOSO",
         "loh + HPFineF\n4-sample train",
         "0/5 above +0.002\nFAIL", "#FEF2F0"),
        (12.3, 5.0, 2.5, 1.5,
         "Stage 8: H_NEW_4\n9-feat drop caller_af",
         "Sanity check\nHCC1395 hold-out",
         "+0.00699 unexpected\n(post-hoc only)", "#FEF7E6"),
    ]
    for x, y, w, h, label, method, verdict, c in stages:
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.05",
            fc=c, ec="black", linewidth=1.2))
        ax.text(x + w/2, y + h - 0.25, label, ha="center", va="top",
                fontsize=9.5, fontweight="bold")
        ax.text(x + w/2, y + h/2 - 0.1, method, ha="center", va="center",
                fontsize=8, color="#444")
        ax.text(x + w/2, y + 0.25, verdict, ha="center", va="bottom",
                fontsize=8, fontweight="bold", color="#7F1D1D" if "FAIL" in verdict or "NEGATIVE" in verdict or "vestigial" in verdict or "≈ 0" in verdict else "#166534" if "+0.02" in verdict or "PASS" in verdict else "#78350F")

    # Arrows
    def arr(x1, y1, x2, y2, label=""):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color="#333", lw=1.3))
        if label:
            ax.text((x1+x2)/2, (y1+y2)/2 + 0.12, label, ha="center", fontsize=7.5, style="italic", color="#666")

    # Top row sequence
    arr(2.85, 9.25, 3.3, 9.25, "pivot")
    arr(5.85, 9.25, 6.3, 9.25, "cross-sample")
    arr(8.85, 9.25, 9.3, 9.25, "gate redesign")
    arr(11.85, 9.25, 12.3, 9.25, "ablation")

    # Cycle3.5 → LOSO (user triggers)
    ax.annotate("", xy=(4.55, 6.5), xytext=(13.55, 8.5),
                arrowprops=dict(arrowstyle="->", color="#C2410C", lw=2,
                                connectionstyle="arc3,rad=0.3"))
    ax.text(7.5, 7.7, "User 質疑\ntriggered LOSO", ha="center", fontsize=9, fontweight="bold",
            color="#C2410C", style="italic")

    # LOSO row
    arr(5.85, 5.75, 6.3, 5.75, "observation-driven")
    arr(8.85, 5.75, 9.3, 5.75, "H_NEW_2 follows obs")
    arr(11.85, 5.75, 12.3, 5.75, "H_NEW_4 sanity")

    # Final verdict box
    ax.add_patch(mpatches.FancyBboxPatch(
        (4.5, 1.8), 7, 1.8, boxstyle="round,pad=0.1",
        fc="#FEF2F0", ec="#C2410C", linewidth=2.5))
    ax.text(8, 3.2, "Final Phase 2 Verdict (2026-05-20)", ha="center", va="center",
            fontsize=12, fontweight="bold", color="#7F1D1D")
    ax.text(8, 2.6, "LR-based universal filter FAILED at sample level", ha="center", va="center",
            fontsize=10.5, fontweight="bold")
    ax.text(8, 2.15, "ISM characterization 保留 (v0.3 ⭐3); HCC1395 +0.02236 是 AF over-fit",
            ha="center", va="center", fontsize=9.5)

    arr(5, 5, 7, 3.6)
    arr(13, 5, 10, 3.6)

    ax.text(8, 10.5, "Phase 2 Validation Pipeline — 9 Stages over 4 Days (2026-05-18 → 2026-05-20)",
            ha="center", fontsize=14, fontweight="bold")

    fig.savefig(OUT / "01_pipeline_flowchart.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/01_pipeline_flowchart.png")


# ─────────── Figure 2: Completeness gap heatmap ──────────────────────────
def fig_gap_matrix():
    # 14 hypotheses × 6 audit dimensions
    # Dim: D(Data) M(Method) S(Stats) C(Causality) R(Reproducibility) N(Negative-evidence)
    # Score: 1=verified, 0.5=partial, 0=missing/violated

    hypotheses = [
        ("H_C1_2 ΔF1 > +0.00242", [1, 1, 1, 0.5, 1, 1]),
        ("H_C1_3 ΔF1 ≥ +0.01", [1, 1, 1, 0.5, 1, 1]),
        ("H_C1_5 cross-sample", [1, 1, 0.5, 0.5, 1, 1]),  # n=5 limited
        ("H_C1_6 BAM-invariant", [1, 1, 1, 1, 1, 1]),
        ("H_C3_1 qualify mean +0.01", [0.5, 1, 0.5, 0.5, 1, 1]),  # n=2 limit
        ("H_C3_3 high-F1 skip", [1, 1, 0.5, 0.5, 1, 1]),
        ("H_C3_2 panel n≥4", [0, 0, 0, 0, 0, 0]),  # PENDING
        ("H_M1a ISM vestigial", [1, 1, 1, 0.5, 1, 1]),
        ("H_A1 caller_af confound", [1, 1, 1, 1, 1, 1]),
        ("LOSO sample-level", [1, 1, 1, 1, 1, 1]),  # NEW core finding
        ("V6 observation per-feat AUC", [1, 1, 1, 0.5, 1, 1]),
        ("H_NEW_2 2-feat LOSO", [1, 1, 1, 0.5, 1, 1]),  # FAIL
        ("H_NEW_4 drop caller_af", [1, 0.5, 0.5, 0.5, 1, 1]),  # pre-reg violated
        ("Caller-F1-headroom mech", [0.5, 0.5, 0.5, 0.5, 1, 1]),  # L3 post-hoc
    ]
    dims = ["D Data", "M Method", "S Stats", "C Causality", "R Reproducibility", "N Negative-evidence"]

    mat = np.array([h[1] for h in hypotheses])
    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    for i in range(len(hypotheses)):
        for j in range(len(dims)):
            v = mat[i, j]
            mark = "✓" if v == 1 else ("△" if v == 0.5 else "✗")
            color = "white" if v < 0.5 else ("black" if v == 0.5 else "darkgreen")
            ax.text(j, i, mark, ha="center", va="center", fontsize=14, fontweight="bold", color=color)
    ax.set_xticks(range(len(dims)))
    ax.set_xticklabels(dims, fontsize=10, rotation=15)
    ax.set_yticks(range(len(hypotheses)))
    ax.set_yticklabels([h[0] for h in hypotheses], fontsize=10)
    ax.set_title("Phase 2 Completeness Audit — 14 Hypotheses × 6 Verification Dimensions\n"
                 "✓ = verified · △ = partial · ✗ = missing/violated/pending",
                 fontsize=12, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, label="Verification Status (0=missing → 1=verified)")
    fig.tight_layout()
    fig.savefig(OUT / "02_completeness_gap_matrix.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/02_completeness_gap_matrix.png")


if __name__ == "__main__":
    print("Generating Phase 2 completeness audit figures ...")
    fig_pipeline()
    fig_gap_matrix()
