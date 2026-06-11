#!/usr/bin/env python3
"""Generate 6 review figures for PI HTML standalone.

Uses cycle 2 + cycle 3 + ablation data TSVs (real numbers, not synthetic).

Outputs (all to cycle3/ablation/figures/review/):
    01_ablation_heatmap.png       4 variants × 5 samples refit ΔF1 heatmap
    02_vestigial_delta.png        Per-sample Δ(full − no-methyl)
    03_transfer_mitigation.png    Caller_af shrinkage disaster mitigation
    04_refit_coef.png             Refit coef magnitude per variant
    05_phase2_timeline.png        v1.0 → cycle 1 → cycle 2 → cycle 3 → step 1.5 timeline
    06_mechanism_quadrant.png     caller F1 × FP density 2D with gate boundaries
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# Ensure CJK font rendering (per known-pitfalls feedback_matplotlib_cjk_font_rule.md)
matplotlib.rcParams["font.sans-serif"] = ["Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
ABLATION = REPO / "research/methyl_augmented_filter_phase2/cycle3/ablation"
CYCLE2_DATA = REPO / "research/methyl_augmented_filter_phase2/cycle2/data/cycle2_cross_sample_delta_f1.tsv"
ABL_DATA = ABLATION / "data/cycle3_step1_5_min_ablation.tsv"
GATE_JSON = REPO / "research/methyl_augmented_filter_phase2/cycle3/cycle3_caller_f1_gate.json"
SUMMARY = ABLATION / "intermediate/cycle3_step1_5_summary.json"

OUT = ABLATION / "figures/review"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]
VARIANT_ORDER = ["full", "no-methyl", "no-caller-af", "no-both"]
VARIANT_LABEL = {
    "full": "Full\n(10 feat)",
    "no-methyl": "No methyl\n(5 feat)",
    "no-caller-af": "No caller_af\n(9 feat)",
    "no-both": "No both\n(4 feat)",
}

# ─────────── Figure 1: ablation heatmap ─────────────────────────────
def fig_ablation_heatmap():
    df = pd.read_csv(ABL_DATA, sep="\t")
    refit = df[df["mode"] == "refit"]
    mat = np.zeros((len(SAMPLE_ORDER), len(VARIANT_ORDER)))
    for i, s in enumerate(SAMPLE_ORDER):
        for j, v in enumerate(VARIANT_ORDER):
            mat[i, j] = float(
                refit[(refit["sample"] == s) & (refit["variant"] == v)]["delta_F1"].iloc[0]
            )

    fig, ax = plt.subplots(figsize=(8.5, 5))
    vmax = max(abs(mat.min()), mat.max(), 0.025)
    im = ax.imshow(mat, cmap="RdYlGn", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(len(SAMPLE_ORDER)):
        for j in range(len(VARIANT_ORDER)):
            val = mat[i, j]
            txt_color = "white" if abs(val) > vmax * 0.5 else "black"
            ax.text(j, i, f"{val:+.5f}", ha="center", va="center",
                    fontsize=11, color=txt_color, fontweight="bold")

    ax.set_xticks(range(len(VARIANT_ORDER)))
    ax.set_xticklabels([VARIANT_LABEL[v] for v in VARIANT_ORDER], fontsize=10)
    ax.set_yticks(range(len(SAMPLE_ORDER)))
    ax.set_yticklabels(SAMPLE_ORDER, fontsize=11)
    ax.set_title("Figure 1 — Refit OOF ΔF1 across 4 variants × 5 samples\n"
                 "(Cycle 3 Step 1.5 ablation; refit mode)",
                 fontsize=12, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, label="ΔF1")
    cbar.ax.set_ylabel("ΔF1 (F1_post − caller_F1)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "01_ablation_heatmap.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/01_ablation_heatmap.png")


# ─────────── Figure 2: vestigial delta ──────────────────────────────
def fig_vestigial_delta():
    df = pd.read_csv(ABL_DATA, sep="\t")
    refit = df[df["mode"] == "refit"]
    full_vals = {s: float(refit[(refit["sample"] == s) & (refit["variant"] == "full")]["delta_F1"].iloc[0])
                 for s in SAMPLE_ORDER}
    nomethyl_vals = {s: float(refit[(refit["sample"] == s) & (refit["variant"] == "no-methyl")]["delta_F1"].iloc[0])
                     for s in SAMPLE_ORDER}
    delta = {s: full_vals[s] - nomethyl_vals[s] for s in SAMPLE_ORDER}

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#15803D" if delta[s] > 0 else ("#9CA3AF" if abs(delta[s]) < 1e-6 else "#C2410C")
              for s in SAMPLE_ORDER]
    bars = ax.bar(SAMPLE_ORDER, [delta[s] for s in SAMPLE_ORDER],
                   color=colors, edgecolor="black", linewidth=0.5)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axhline(0.003, color="#15803D", linestyle="--", linewidth=1, label="Cohen +0.003 threshold (H_M1a PASS)")
    ax.axhline(0.001, color="#A16207", linestyle=":", linewidth=1, label="FAIL boundary +0.001")
    for s, d in zip(SAMPLE_ORDER, [delta[s] for s in SAMPLE_ORDER]):
        ax.text(s, d + (0.00015 if d >= 0 else -0.00015), f"{d:+.5f}",
                ha="center", va="bottom" if d >= 0 else "top", fontsize=10, fontweight="bold")

    ax.set_ylabel("Δ(full − no-methyl) ΔF1", fontsize=11)
    ax.set_xlabel("Sample")
    ax.set_title("Figure 2 — Per-sample ISM incremental contribution\n"
                 "Δ(full ΔF1 − no-methyl ΔF1) = methylation block's incremental value\n"
                 "Positive = ISM adds value · Negative = ISM hurts · 0 = vestigial",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "02_vestigial_delta.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/02_vestigial_delta.png")


# ─────────── Figure 3: transfer mitigation ──────────────────────────
def fig_transfer_mitigation():
    df = pd.read_csv(ABL_DATA, sep="\t")
    transfer = df[df["mode"] == "transfer_fixed"]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    width = 0.20
    x = np.arange(len(SAMPLE_ORDER))
    colors = {"full": "#7F1D1D", "no-methyl": "#C2410C",
              "no-caller-af": "#A16207", "no-both": "#15803D"}
    for i, v in enumerate(VARIANT_ORDER):
        vals = [float(transfer[(transfer["sample"] == s) & (transfer["variant"] == v)]["delta_F1"].iloc[0])
                for s in SAMPLE_ORDER]
        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, vals, width, color=colors[v],
                      label=VARIANT_LABEL[v].replace("\n", " "), edgecolor="black", linewidth=0.4)
        for j, val in enumerate(vals):
            ax.text(x[j] + offset, val + (0.005 if val >= 0 else -0.005),
                    f"{val:+.3f}", ha="center", va="bottom" if val >= 0 else "top",
                    fontsize=7.5, fontweight="bold" if val < -0.05 else "normal")

    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLE_ORDER, fontsize=11)
    ax.set_ylabel("ΔF1 (transfer τ=0.39 fixed)", fontsize=11)
    ax.set_title("Figure 3 — Transfer mode ΔF1 with coef shrinkage\n"
                 "(Apply cycle 1 filter with target coefs forced to 0)\n"
                 "Caller_af shrinkage (gold) mitigates HCC1954 disaster from −0.38 to −0.13",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "03_transfer_mitigation.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/03_transfer_mitigation.png")


# ─────────── Figure 4: refit coef ──────────────────────────────────
def fig_refit_coef():
    df = pd.read_csv(ABL_DATA, sep="\t")
    refit = df[df["mode"] == "refit"]
    # Extract refit_coef_json per (sample, variant)
    rows = []
    for s in SAMPLE_ORDER:
        for v in VARIANT_ORDER:
            row = refit[(refit["sample"] == s) & (refit["variant"] == v)].iloc[0]
            try:
                coefs = json.loads(row["refit_coef_json"])
            except (json.JSONDecodeError, TypeError):
                continue
            for feat, coef in coefs.items():
                rows.append({"sample": s, "variant": v, "feature": feat, "coef": coef})
    cdf = pd.DataFrame(rows)

    # Focus on full variant — show coef magnitude across 5 samples
    full = cdf[cdf["variant"] == "full"]
    pivot = full.pivot(index="feature", columns="sample", values="coef")
    # Reorder samples
    pivot = pivot[SAMPLE_ORDER]
    # Sort features by abs mean
    pivot["_abs_mean"] = pivot.abs().mean(axis=1)
    pivot = pivot.sort_values("_abs_mean", ascending=False).drop(columns="_abs_mean")

    fig, ax = plt.subplots(figsize=(9.5, 6))
    vmax = max(abs(pivot.values.min()), pivot.values.max())
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            v = pivot.values[i, j]
            txt_color = "white" if abs(v) > vmax * 0.5 else "black"
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                    fontsize=9, color=txt_color, fontweight="bold")

    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, fontsize=11)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=10)
    ax.set_title("Figure 4 — Refit LR coefficients (full 10-feature model)\n"
                 "per-sample magnitude (sorted by |mean| across 5 samples)\n"
                 "caller_af top across all samples; methyl features bottom rank",
                 fontsize=11, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, label="Standardized LR coef")
    fig.tight_layout()
    fig.savefig(OUT / "04_refit_coef.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/04_refit_coef.png")


# ─────────── Figure 5: timeline ─────────────────────────────────────
def fig_timeline():
    fig, ax = plt.subplots(figsize=(12, 5))
    milestones = [
        ("v1.0 cell-gated\npilot", "2026-05-18 上", 0.00242, "marginal", "#FCD34D"),
        ("Cycle 1 global LR\n(HCC1395)", "5-18 晚", 0.02236, "⭐3 strong", "#15803D"),
        ("Cycle 2 transfer\n4 new samples mean", "5-19 凌", -0.094, "NEGATIVE", "#C2410C"),
        ("Cycle 2 re-fit\n4 new samples mean", "5-19 凌", 0.006, "MIXED", "#A16207"),
        ("Cycle 3 gate rule\nqualifying mean", "5-19 凌", 0.01499, "PASS gate", "#15803D"),
        ("Cycle 3 ablation\nno-methyl HCC1395", "5-20", 0.02171, "ISM vestigial", "#1E3A8A"),
        ("Cycle 3 ablation\nno-methyl 5-sample mean", "5-20", 0.00630, "略勝 full", "#1E3A8A"),
    ]

    x = list(range(len(milestones)))
    y = [m[2] for m in milestones]
    colors = [m[4] for m in milestones]
    bars = ax.bar(x, y, color=colors, edgecolor="black", linewidth=0.5, width=0.6)

    for i, (label, date, val, verdict, col) in enumerate(milestones):
        ax.text(i, val + (0.001 if val >= 0 else -0.005),
                f"{val:+.5f}\n({verdict})",
                ha="center", va="bottom" if val >= 0 else "top",
                fontsize=8.5, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"{m[0]}\n{m[1]}" for m in milestones], fontsize=9)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axhline(0.01, color="#15803D", linestyle="--", linewidth=1, label="Cohen +0.01 threshold")
    ax.axhline(0.003, color="#A16207", linestyle=":", linewidth=1, label="Cohen +0.003 ribbon")
    ax.set_ylabel("ΔF1", fontsize=11)
    ax.set_title("Figure 5 — Phase 2 ΔF1 evolution timeline (2 days, 4 days)\n"
                 "v1.0 marginal → cycle 1 internal valid → cycle 2 cross-sample fail → cycle 3 gate rule → ablation reframe",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "05_phase2_timeline.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/05_phase2_timeline.png")


# ─────────── Figure 6: mechanism quadrant ───────────────────────────
def fig_mechanism_quadrant():
    # 5-sample data
    cycle2 = pd.read_csv(CYCLE2_DATA, sep="\t")
    refit = cycle2[cycle2["mode"] == "refit"].set_index("sample")
    samples = {}
    for s in SAMPLE_ORDER:
        r = refit.loc[s]
        tp = float(r["tp_total_used"])
        fp = float(r["fp_total_used"])
        samples[s] = {
            "caller_F1": float(r["caller_F1"]),
            "fp_density": fp / (tp + fp) * 100,
            "delta_F1": float(r["delta_F1"]),
        }

    fig, ax = plt.subplots(figsize=(10, 7))
    # Gate boundaries: caller F1 < 0.80, FP density > 10%
    ax.add_patch(mpatches.Rectangle((0.30, 10), 0.50, 12,
                                     facecolor="#DCFCE7", edgecolor="#15803D",
                                     linewidth=1.5, linestyle="--", alpha=0.5))
    ax.text(0.55, 20, "Qualifying quadrant\ncaller F1 < 0.80 AND FP density > 10%",
            ha="center", va="center", fontsize=11, fontweight="bold", color="#15803D")

    for s, d in samples.items():
        if d["caller_F1"] < 0.80 and d["fp_density"] > 10:
            color, marker = "#15803D", "o"
        elif d["delta_F1"] < -0.10:
            color, marker = "#C2410C", "X"
        else:
            color, marker = "#6B7280", "s"
        ax.scatter(d["caller_F1"], d["fp_density"], s=300, c=color,
                   marker=marker, edgecolors="white", linewidth=2, zorder=3)
        ax.annotate(f"{s}\nΔF1={d['delta_F1']:+.4f}",
                    xy=(d["caller_F1"], d["fp_density"]),
                    xytext=(8, 8), textcoords="offset points",
                    fontsize=10, fontweight="bold")

    ax.axvline(0.80, color="#15803D", linestyle=":", linewidth=1)
    ax.axhline(10, color="#15803D", linestyle=":", linewidth=1)
    ax.set_xlim(0.30, 1.0)
    ax.set_ylim(0, 22)
    ax.set_xlabel("Caller F1 baseline", fontsize=12)
    ax.set_ylabel("FP density (%)", fontsize=12)
    ax.set_title("Figure 6 — Caller-F1-headroom mechanism (5-sample scatter)\n"
                 "Filter only effective in qualifying quadrant (low caller F1 + high FP density)\n"
                 "HCC1395 + HCC1937 qualify; 3 high-F1 samples skip filter",
                 fontsize=11, fontweight="bold")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "06_mechanism_quadrant.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {OUT}/06_mechanism_quadrant.png")


if __name__ == "__main__":
    print("Generating review figures …")
    fig_ablation_heatmap()
    fig_vestigial_delta()
    fig_transfer_mitigation()
    fig_refit_coef()
    fig_timeline()
    fig_mechanism_quadrant()
    print(f"\nAll 6 figures saved to: {OUT}")
