#!/usr/bin/env python3
"""F1 analysis figures — 4 plots for baseline vs V6 F1 layer.

F6: F1 vs NG threshold (4 BAMs as lines + caller-alone dashed reference)
F7: Precision-Recall curve (4 BAMs scatter at T=1-5 + caller-alone star)
F8: F1 stratified by LOH zone (grouped bars 4 BAMs × 2 zones × T=2,3)
F9: ΔF1 (V6 - baseline) sweep across thresholds
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = ['DejaVu Sans', 'Droid Sans Fallback', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

STEP1 = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way")
FIGURES = STEP1.parent / "figures" / "baseline_vs_v6"
FIGURES.mkdir(parents=True, exist_ok=True)

BAMS = ["baseline", "V3F", "V5", "V6"]
COLORS = {"baseline": "#d62728", "V3F": "#ff7f0e", "V5": "#2ca02c", "V6": "#1f77b4"}
MARKERS = {"baseline": "o", "V3F": "s", "V5": "^", "V6": "D"}
CALLER_ALONE_F1 = 0.7166


def load_thresholds():
    rows = []
    with (STEP1 / "step1_f1_analysis_thresholds.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "BAM": r["BAM"],
                "T": int(r["NG_threshold"]),
                "TP": int(r["TP"]),
                "FP": int(r["FP"]),
                "FN": int(r["FN"]),
                "P": float(r["precision"]),
                "R": float(r["recall"]),
                "F1": float(r["F1"]),
                "delta_caller": float(r["delta_F1_vs_caller_alone"]),
            })
    return rows


def load_per_zone():
    rows = []
    with (STEP1 / "step1_f1_analysis_per_zone.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "BAM": r["BAM"],
                "T": int(r["NG_threshold"]),
                "zone": r["LOH_zone"],
                "F1": float(r["F1"]),
                "P": float(r["precision"]),
                "R": float(r["recall"]),
                "F1_caller_zone": float(r["F1_caller_zone"]),
            })
    return rows


def f6_f1_vs_threshold(rows):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    thresholds = sorted(set(r["T"] for r in rows))
    for bam in BAMS:
        f1s = [next(r["F1"] for r in rows if r["BAM"] == bam and r["T"] == t) for t in thresholds]
        ax.plot(thresholds, f1s, marker=MARKERS[bam], color=COLORS[bam],
                label=f"{bam} (best F1={max(f1s):.4f} @ T={thresholds[np.argmax(f1s)]})",
                linewidth=2.2, markersize=10)
    # caller-alone dashed reference
    ax.axhline(CALLER_ALONE_F1, color='#222', linestyle='--', linewidth=2, alpha=0.7,
               label=f"caller-alone F1 = {CALLER_ALONE_F1:.4f} (no filter)")
    # Phase 2 Cycle 1 LR reference
    ax.axhline(CALLER_ALONE_F1 + 0.02236, color='#15803D', linestyle=':', linewidth=1.8, alpha=0.7,
               label=f"Phase 2 Cycle 1 multi-axis LR (V6 only) = {CALLER_ALONE_F1+0.02236:.4f}")
    ax.set_xlabel("NG_off threshold T (keep regions with NG_off ≥ T)", fontsize=11)
    ax.set_ylabel("F1 on pileup subset (TP=30,490 / FP=4,842 / FN_caller=19,288)", fontsize=11)
    ax.set_title("F6 — F1 vs ISM NG hard threshold filter, 4 BAMs comparison\n"
                 "Hard NG threshold barely matches caller-alone at T=1; multi-axis LR (Cycle 1) needed to beat caller", fontsize=12)
    ax.set_xticks(thresholds)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='center right', fontsize=9)
    ax.set_ylim([0, 0.78])
    plt.tight_layout()
    out = FIGURES / "F6_f1_vs_ng_threshold.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def f7_pr_curve(rows):
    fig, ax = plt.subplots(figsize=(11, 7))
    for bam in BAMS:
        ps = []
        rs = []
        for r in sorted([r for r in rows if r["BAM"] == bam], key=lambda x: x["T"]):
            ps.append(r["P"])
            rs.append(r["R"])
        ax.plot(rs, ps, color=COLORS[bam], marker=MARKERS[bam], markersize=12,
                linewidth=2, label=bam, alpha=0.85)
        # annotate T values
        for i, r in enumerate(sorted([r for r in rows if r["BAM"] == bam], key=lambda x: x["T"])):
            ax.annotate(f"T={r['T']}", (rs[i], ps[i]),
                        xytext=(5, 5), textcoords='offset points', fontsize=8, color=COLORS[bam])
    # caller-alone star
    p_caller = 30490 / (30490 + 4842)
    r_caller = 30490 / (30490 + 19288)
    ax.scatter([r_caller], [p_caller], color='#222', marker='*', s=400, zorder=10,
               label=f"caller-alone (P={p_caller:.3f}, R={r_caller:.3f}, F1={CALLER_ALONE_F1:.4f})")
    # F1 iso-contour lines
    for f1_iso in [0.3, 0.5, 0.7, 0.72]:
        r_iso = np.linspace(0.01, 1, 100)
        p_iso = f1_iso * r_iso / (2 * r_iso - f1_iso)
        valid = (p_iso > 0) & (p_iso < 1)
        ax.plot(r_iso[valid], p_iso[valid], color='gray', linewidth=0.6, linestyle=':', alpha=0.5)
        # label at mid
        mid_idx = len(r_iso) // 2
        if valid[mid_idx]:
            ax.annotate(f"F1={f1_iso}", (r_iso[mid_idx], p_iso[mid_idx]),
                        fontsize=7, color='gray', alpha=0.7)
    ax.set_xlabel("Recall", fontsize=11)
    ax.set_ylabel("Precision", fontsize=11)
    ax.set_title("F7 — Precision-Recall curve, 4 BAMs at NG thresholds T=1-5\n"
                 "★ caller-alone reference; iso-F1 contours dashed gray", fontsize=12)
    ax.set_xlim([0, 0.8])
    ax.set_ylim([0.5, 1.0])
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower left', fontsize=10)
    plt.tight_layout()
    out = FIGURES / "F7_precision_recall_curve.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def f8_zone_f1(rows):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, zone in zip(axes, ["LOH-positive", "LOH-NA"]):
        zone_rows = [r for r in rows if r["zone"] == zone]
        thresholds = sorted(set(r["T"] for r in zone_rows))
        x = np.arange(len(BAMS))
        width = 0.35

        f1s_t2 = [next((r["F1"] for r in zone_rows if r["BAM"] == b and r["T"] == 2), 0) for b in BAMS]
        f1s_t3 = [next((r["F1"] for r in zone_rows if r["BAM"] == b and r["T"] == 3), 0) for b in BAMS]
        f1_caller_zone = zone_rows[0]["F1_caller_zone"] if zone_rows else 0

        ax.bar(x - width/2, f1s_t2, width, label='NG≥2 filter', color=[COLORS[b] for b in BAMS], alpha=0.85)
        ax.bar(x + width/2, f1s_t3, width, label='NG≥3 filter',
               color=[COLORS[b] for b in BAMS], alpha=0.85, hatch='//')

        ax.axhline(f1_caller_zone, color='#222', linestyle='--', linewidth=2,
                   label=f"caller-alone zone F1 = {f1_caller_zone:.4f}")

        ax.set_xticks(x)
        ax.set_xticklabels(BAMS)
        ax.set_ylabel("F1 (zone-stratified)")
        ax.set_title(f"{zone} zone")
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, axis='y', alpha=0.3)
        for i, (t2, t3) in enumerate(zip(f1s_t2, f1s_t3)):
            ax.text(i - width/2, t2 + 0.005, f"{t2:.3f}", ha='center', fontsize=8)
            ax.text(i + width/2, t3 + 0.005, f"{t3:.3f}", ha='center', fontsize=8)

    fig.suptitle("F8 — F1 stratified by LOH zone (HCC1395 full genome)\n"
                 "LOH-positive zone: 1,567 truth positives (small) — V6 ≈ V3F; LOH-NA: 48,211 truth positives — V6 dominates",
                 fontsize=12)
    plt.tight_layout()
    out = FIGURES / "F8_f1_by_loh_zone.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def f9_delta_f1(rows):
    fig, ax = plt.subplots(figsize=(11, 6))
    thresholds = sorted(set(r["T"] for r in rows))
    deltas_v6_base = []
    deltas_v6_v3f = []
    deltas_v6_v5 = []
    for T in thresholds:
        f1_b = next(r["F1"] for r in rows if r["BAM"] == "baseline" and r["T"] == T)
        f1_v3f = next(r["F1"] for r in rows if r["BAM"] == "V3F" and r["T"] == T)
        f1_v5 = next(r["F1"] for r in rows if r["BAM"] == "V5" and r["T"] == T)
        f1_v6 = next(r["F1"] for r in rows if r["BAM"] == "V6" and r["T"] == T)
        deltas_v6_base.append(f1_v6 - f1_b)
        deltas_v6_v3f.append(f1_v6 - f1_v3f)
        deltas_v6_v5.append(f1_v6 - f1_v5)

    x = np.arange(len(thresholds))
    width = 0.25

    ax.bar(x - width, deltas_v6_base, width, label='V6 - baseline', color='#d62728', alpha=0.85)
    ax.bar(x, deltas_v6_v3f, width, label='V6 - V3F', color='#ff7f0e', alpha=0.85)
    ax.bar(x + width, deltas_v6_v5, width, label='V6 - V5', color='#2ca02c', alpha=0.85)

    ax.axhline(0, color='#222', linewidth=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"T={t}" for t in thresholds])
    ax.set_xlabel("NG_off threshold T")
    ax.set_ylabel("Δ F1 (V6 - other BAM)")
    ax.set_title("F9 — V6 F1 advantage over baseline / V3F / V5 across NG thresholds\n"
                 "V6 strictly dominates baseline at all T; V6 ≥ V5 at all T; V6 ≥ V3F at T=1-3, V6 < V3F at T=4-5", fontsize=12)
    ax.grid(True, axis='y', alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    for i, (db, dv3f, dv5) in enumerate(zip(deltas_v6_base, deltas_v6_v3f, deltas_v6_v5)):
        ax.text(i - width, db + 0.005, f"{db:+.3f}", ha='center', fontsize=8)
        ax.text(i, dv3f + 0.005, f"{dv3f:+.3f}", ha='center', fontsize=8)
        ax.text(i + width, dv5 + 0.005, f"{dv5:+.3f}", ha='center', fontsize=8)
    plt.tight_layout()
    out = FIGURES / "F9_delta_f1_v6_vs_others.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def main():
    print("[plot F6-F9] reading data ...")
    rows = load_thresholds()
    zone_rows = load_per_zone()
    print(f"  threshold rows: {len(rows)}, zone rows: {len(zone_rows)}")
    f6_f1_vs_threshold(rows)
    f7_pr_curve(rows)
    f8_zone_f1(zone_rows)
    f9_delta_f1(rows)
    print(f"\n[plot] 4 figures saved to {FIGURES}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
