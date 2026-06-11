#!/usr/bin/env python3
"""Plot baseline vs V6 head-to-head distributions (HCC1395).

Produces 6 figures in figures/baseline_vs_v6/:

F1. global_summary_4panel.png — 4-panel summary bar plots
    (a) HP1:HP2 ratio across 4 BAMs (log scale)
    (b) marker count NG≥3 TP vs FP across 4 BAMs
    (c) hp=33 conservative-ambiguous reads
    (d) marker rate (TP / (TP+FP))

F2. per_chr_hp_ratio.png — per-chromosome HP1:HP2 ratio across 4 BAMs

F3. marker_rate_by_loh_zone.png — marker rate stratified by LOH bucket

F4. hp_family_distribution.png — stacked bar of hp class distribution per BAM

F5. tp_fp_density_off_vs_on.png — NG distribution off vs on, TP vs FP

F6. before_after_v6.png — single-panel narrative: baseline → V6 key shift
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


def load_summary():
    rows = []
    with (STEP1 / "step1_baseline_vs_v6_summary.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    return {r["BAM"]: r for r in rows}


def load_per_chr():
    with (STEP1 / "step1_baseline_vs_v6_per_chr.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_per_zone():
    with (STEP1 / "step1_baseline_vs_v6_per_zone.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f1_global_summary(summary):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # (a) HP1:HP2 ratio
    ax = axes[0, 0]
    ratios = [float(summary[b]["hp1_hp2_ratio"]) for b in BAMS]
    bars = ax.bar(BAMS, ratios, color=[COLORS[b] for b in BAMS])
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='balanced (1:1)')
    ax.set_ylabel("HP1 series : HP2 series ratio")
    ax.set_title("(a) Genome-wide HP1:HP2 ratio (ALT-supporting reads)")
    ax.legend()
    for bar, val in zip(bars, ratios):
        ax.text(bar.get_x() + bar.get_width()/2, val * 1.05, f"{val:.2f}",
                ha='center', fontsize=10)

    # (b) marker NG>=3 counts
    ax = axes[0, 1]
    tp_counts = [int(summary[b]["marker_NGge3_TP"]) for b in BAMS]
    fp_counts = [int(summary[b]["marker_NGge3_FP"]) for b in BAMS]
    x = np.arange(len(BAMS))
    width = 0.35
    ax.bar(x - width/2, tp_counts, width, label='TP markers', color='steelblue')
    ax.bar(x + width/2, fp_counts, width, label='FP markers', color='salmon')
    ax.set_xticks(x)
    ax.set_xticklabels(BAMS)
    ax.set_ylabel("Region count")
    ax.set_title("(b) Marker filter NG≥3 region counts")
    ax.legend()
    for i, (t, f) in enumerate(zip(tp_counts, fp_counts)):
        ax.text(i - width/2, t + max(tp_counts)*0.01, f"{t:,}", ha='center', fontsize=8)
        ax.text(i + width/2, f + max(tp_counts)*0.01, f"{f:,}", ha='center', fontsize=8)

    # (c) hp=3 + hp=33 ambiguous reads (in HCC1395 TO data, all are in hp=3 bucket)
    ax = axes[1, 0]
    hp_amb = [int(summary[b]["hp_ambig_reads"]) for b in BAMS]
    bars = ax.bar(BAMS, hp_amb, color=[COLORS[b] for b in BAMS])
    ax.set_ylabel("hp=3+33 ambiguous reads (TO mode: 100% in hp=3)")
    ax.set_title("(c) Somatic-ambiguous tag reads (V3F-style conservative)")
    for bar, val in zip(bars, hp_amb):
        ax.text(bar.get_x() + bar.get_width()/2, val + max(hp_amb)*0.01,
                f"{val:,}", ha='center', fontsize=9)

    # (d) marker rate
    ax = axes[1, 1]
    rates = [float(summary[b]["marker_rate_off"]) for b in BAMS]
    bars = ax.bar(BAMS, rates, color=[COLORS[b] for b in BAMS])
    ax.axhline(0.85, color='red', linestyle='--', alpha=0.5, label='target ≥ 0.85')
    ax.set_ylabel("Marker rate (TP / (TP+FP))")
    ax.set_title("(d) Marker rate (purity at NG≥3)")
    ax.set_ylim([0.7, 1.0])
    ax.legend()
    for bar, val in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.005,
                f"{val:.4f}", ha='center', fontsize=9)

    fig.suptitle("HCC1395 Full Genome — baseline → V3F → V5 → V6 4-way comparison\n"
                 "TP=30,490 | FP=4,842 | Window=±5kb | ISM flag=off",
                 fontsize=12, y=0.995)
    plt.tight_layout()
    out = FIGURES / "F1_global_summary_4panel.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def f2_per_chr_ratio(per_chr):
    fig, ax = plt.subplots(figsize=(15, 6))
    chroms = [r["chr"] for r in per_chr if r["chr"].startswith("chr")]
    # natural sort
    def chr_key(c):
        s = c.replace("chr", "")
        return (s != "X", s != "Y", int(s) if s.isdigit() else 99)
    chroms = sorted(chroms, key=chr_key)
    x = np.arange(len(chroms))
    width = 0.2

    for i, bam in enumerate(BAMS):
        ratios = []
        for c in chroms:
            r = next((row for row in per_chr if row["chr"] == c), None)
            if r is None:
                ratios.append(0)
                continue
            try:
                ratios.append(float(r[f"{bam}_ratio"]))
            except (ValueError, KeyError):
                ratios.append(0)
        ax.bar(x + (i - 1.5) * width, ratios, width, label=bam, color=COLORS[bam])

    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='balanced')
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel("HP1 series : HP2 series ratio")
    ax.set_title("Per-chromosome HP1:HP2 ratio across 4 BAMs (HCC1395 ALT reads, ISM flag=off)\n"
                 "baseline expected to show systematic HP1 bias; V6 expected near-balance except known cnLOH chrs")
    ax.legend(loc='upper right')
    ax.set_yscale('symlog', linthresh=1)
    plt.tight_layout()
    out = FIGURES / "F2_per_chr_hp_ratio.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def f3_marker_rate_by_zone(per_zone):
    fig, ax = plt.subplots(figsize=(13, 6))
    zones = [r["LOH_zone"] for r in per_zone if r["LOH_zone"] != "NA" and r["LOH_zone"] != "Unknown"]
    zones = sorted(set(zones))
    x = np.arange(len(zones))
    width = 0.2

    for i, bam in enumerate(BAMS):
        rates = []
        for z in zones:
            r = next((row for row in per_zone if row["LOH_zone"] == z), None)
            try:
                rates.append(float(r[f"{bam}_marker_rate"]))
            except (ValueError, TypeError, KeyError):
                rates.append(0)
        ax.bar(x + (i - 1.5) * width, rates, width, label=bam, color=COLORS[bam])

    ax.axhline(0.85, color='red', linestyle='--', alpha=0.5, label='target ≥ 0.85')
    ax.set_xticks(x)
    ax.set_xticklabels(zones, rotation=30, ha='right')
    ax.set_ylabel("Marker rate (TP / (TP+FP))")
    ax.set_title("Marker filter purity (NG≥3) by LOH zone — baseline vs V3F/V5/V6")
    ax.legend()
    plt.tight_layout()
    out = FIGURES / "F3_marker_rate_by_loh_zone.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def f4_hp_distribution(summary):
    # need to read findings.json for hp_family breakdown
    findings = json.load((STEP1 / "step1_baseline_vs_v6_findings.json").open())
    # find HP family from summary_rows.. but summary doesn't have all HP keys
    # Skip if not enough data — derive from totals
    fig, ax = plt.subplots(figsize=(11, 6))
    # Use bars: hp1 / hp2 / hp33 / others stacked
    bottoms = np.zeros(len(BAMS))
    cats = [
        ("hp1 series", "hp1_series", "#4daf4a"),
        ("hp2 series", "hp2_series", "#377eb8"),
        ("hp=33 (conservative)", "hp_ambig_reads", "#ff7f00"),
    ]
    totals = [int(summary[b]["total_reads"]) for b in BAMS]
    for label, key, color in cats:
        vals = [int(summary[b][key]) for b in BAMS]
        ax.bar(BAMS, vals, bottom=bottoms, label=label, color=color)
        bottoms += np.array(vals)

    # remaining = total - sum
    remain = np.array(totals) - bottoms
    ax.bar(BAMS, remain, bottom=bottoms, label="hp=0/other", color="lightgray")

    ax.set_ylabel("ALT-supporting read count")
    ax.set_title("HP class distribution across 4 BAMs (HCC1395 full genome, flag=off)")
    ax.legend()
    plt.tight_layout()
    out = FIGURES / "F4_hp_family_distribution.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def f5_narrative_before_after(summary):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    base_ratio = float(summary["baseline"]["hp1_hp2_ratio"])
    v6_ratio = float(summary["V6"]["hp1_hp2_ratio"])
    ax.bar(["baseline", "V6"], [base_ratio, v6_ratio],
           color=["#d62728", "#1f77b4"])
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    ax.set_yscale('log')
    ax.set_ylabel("HP1:HP2 ratio (log scale)")
    ax.set_title("HP1:HP2 ratio\n(closer to 1 = better)")
    for i, val in enumerate([base_ratio, v6_ratio]):
        ax.text(i, val * 1.1, f"{val:.2f}", ha='center', fontsize=11)

    ax = axes[1]
    base_rate = float(summary["baseline"]["marker_rate_off"])
    v6_rate = float(summary["V6"]["marker_rate_off"])
    ax.bar(["baseline", "V6"], [base_rate, v6_rate],
           color=["#d62728", "#1f77b4"])
    ax.axhline(0.85, color='red', linestyle='--', alpha=0.5, label='target ≥ 0.85')
    ax.set_ylim([0.7, 1.0])
    ax.set_ylabel("Marker rate at NG≥3")
    ax.set_title("Marker purity\n(higher = better)")
    for i, val in enumerate([base_rate, v6_rate]):
        ax.text(i, val + 0.005, f"{val:.4f}", ha='center', fontsize=11)
    ax.legend()

    ax = axes[2]
    base_amb = int(summary["baseline"]["hp_ambig_reads"])
    v6_amb = int(summary["V6"]["hp_ambig_reads"])
    ax.bar(["baseline", "V6"], [base_amb, v6_amb],
           color=["#d62728", "#1f77b4"])
    ax.set_ylabel("hp=3 (somatic-ambiguous) reads")
    ax.set_title("Conservative ambiguous tagging\n(V6 expected ≫ baseline)")
    for i, val in enumerate([base_amb, v6_amb]):
        ax.text(i, val + max(base_amb, v6_amb)*0.02, f"{val:,}", ha='center', fontsize=11)

    fig.suptitle("baseline → V6: 3 key shifts (HCC1395)", fontsize=13, y=1.02)
    plt.tight_layout()
    out = FIGURES / "F5_narrative_before_after.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"  → {out}")


def main() -> int:
    print(f"[plot] reading from {STEP1}")
    summary = load_summary()
    per_chr = load_per_chr()
    per_zone = load_per_zone()
    print(f"  4 BAMs loaded: {list(summary.keys())}")
    print(f"  {len(per_chr)} chromosomes")
    print(f"  {len(per_zone)} LOH zones")

    f1_global_summary(summary)
    f2_per_chr_ratio(per_chr)
    f3_marker_rate_by_zone(per_zone)
    f4_hp_distribution(summary)
    f5_narrative_before_after(summary)

    print(f"\n[plot] all figures saved to {FIGURES}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
