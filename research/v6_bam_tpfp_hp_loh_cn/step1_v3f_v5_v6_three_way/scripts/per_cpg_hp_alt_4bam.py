#!/usr/bin/env python3
"""Goal 1 — per-CpG × HP × ALT/REF association strength across 4 BAMs.

Fallback design (per plan R1 — significance_summary aggregator too complex):
  Use reads.tsv hp + alt_support columns directly.
  Per ISM run (16 total = 4 BAM × 2 flag × 2 label), walk all 35,332 region dirs,
  accumulate global HP × ALT × TP/FP × BAM contingency.

Output: per-BAM HP × ALT contingency + Cramer's V + Odds Ratio (HP1 vs HP2 for ALT vs REF).

This quantifies Goal 1 foundation: "is HP × ALT association inflated by priority bug?"
  - baseline: priority bug → HP1 過度與 ALT 關聯 (spurious)
  - V3F/V6: balanced
  - V5: partial revert
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PHASEC = REPO / "research/paired_priority_bug_audit/phaseC_genome_three_way"
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"
FIG_DIR = REPO / "research/v6_bam_tpfp_hp_loh_cn/figures/baseline_vs_v6"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# HP collapse to 2 buckets for HP1 series vs HP2 series
HP1_SERIES = {"1", "1-1", "11"}
HP2_SERIES = {"2", "2-1", "21"}


def count_reads_in_file(reads_tsv):
    """Return (hp1_alt, hp1_ref, hp2_alt, hp2_ref, other_alt, other_ref)."""
    counts = {"hp1_alt": 0, "hp1_ref": 0, "hp2_alt": 0, "hp2_ref": 0,
              "other_alt": 0, "other_ref": 0}
    try:
        with open(reads_tsv, "rb") as fh:
            header = fh.readline().decode().rstrip().split("\t")
            try:
                hp_i = header.index("hp")
                alt_i = header.index("alt_support")
            except ValueError:
                return counts
            for line in fh:
                cols = line.decode().rstrip().split("\t")
                if len(cols) <= max(hp_i, alt_i):
                    continue
                hp = cols[hp_i]
                alt = cols[alt_i]
                if hp in HP1_SERIES:
                    key = "hp1_alt" if alt == "ALT" else "hp1_ref"
                elif hp in HP2_SERIES:
                    key = "hp2_alt" if alt == "ALT" else "hp2_ref"
                else:
                    key = "other_alt" if alt == "ALT" else "other_ref"
                counts[key] += 1
    except Exception:
        pass
    return counts


def find_reads_tsv(run_dir):
    files = []
    for root, _, fn in os.walk(run_dir):
        if "reads.tsv" in fn:
            files.append(os.path.join(root, "reads.tsv"))
    return files


def aggregate_run(run_dir, n_workers=16):
    """Sum HP × ALT counts across all regions in this run."""
    files = find_reads_tsv(run_dir)
    total = Counter()
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for c in ex.map(count_reads_in_file, files, chunksize=200):
            for k, v in c.items():
                total[k] += v
    return dict(total), len(files)


def cramers_v(table_2x2):
    """Cramer's V for 2x2 contingency. table = [[a, b], [c, d]]."""
    a, b = table_2x2[0]
    c, d = table_2x2[1]
    n = a + b + c + d
    if n == 0:
        return 0.0
    chi2_num = n * (a * d - b * c) ** 2
    chi2_denom = (a + b) * (c + d) * (a + c) * (b + d)
    if chi2_denom == 0:
        return 0.0
    chi2 = chi2_num / chi2_denom
    # phi for 2x2 = sqrt(chi2 / n) = Cramer's V
    return float(np.sqrt(chi2 / n))


def odds_ratio(table_2x2):
    """OR = (a * d) / (b * c), HP1-ALT odds vs HP2-ALT odds."""
    a, b = table_2x2[0]
    c, d = table_2x2[1]
    if b == 0 or c == 0:
        return float("inf") if a * d > 0 else 0.0
    return (a * d) / (b * c)


def main():
    print("Aggregating reads HP × ALT across 4 BAM × off × {tp, fp} (8 runs)", flush=True)
    print(f"  Reading from {PHASEC}", flush=True)

    bams = ["baseline", "V3F", "V5", "V6"]
    flags = ["off"]  # restrict to flag=off for foundational HP × ALT analysis
    labels = ["tp", "fp"]

    results = {b: {} for b in bams}
    summary_rows = []

    for bam in bams:
        for flag in flags:
            for label in labels:
                run_dir = PHASEC / f"{bam}_{flag}_{label}"
                if not run_dir.exists():
                    print(f"  SKIP {bam}_{flag}_{label}: missing", flush=True)
                    continue
                t0 = time.time()
                counts, n_regions = aggregate_run(run_dir, n_workers=16)
                elapsed = time.time() - t0
                results[bam][label] = counts
                total_reads = sum(counts.values())
                hp1_alt = counts.get("hp1_alt", 0)
                hp1_ref = counts.get("hp1_ref", 0)
                hp2_alt = counts.get("hp2_alt", 0)
                hp2_ref = counts.get("hp2_ref", 0)

                print(f"  {bam}_{flag}_{label}: {n_regions} regions, {total_reads} reads, "
                      f"HP1[ALT={hp1_alt} REF={hp1_ref}] HP2[ALT={hp2_alt} REF={hp2_ref}]  "
                      f"({elapsed:.1f}s)", flush=True)

                summary_rows.append({
                    "BAM": bam, "flag": flag, "label": label,
                    "n_regions": n_regions, "total_reads": total_reads,
                    "HP1_ALT": hp1_alt, "HP1_REF": hp1_ref,
                    "HP2_ALT": hp2_alt, "HP2_REF": hp2_ref,
                    "other_ALT": counts.get("other_alt", 0),
                    "other_REF": counts.get("other_ref", 0),
                })

    # Compute per-BAM stats: HP × ALT (within TP regions, FP regions, combined)
    print(f"\n=== Per-BAM HP1 × HP2 vs ALT × REF (TP + FP regions combined) ===", flush=True)
    bam_stats = []
    for bam in bams:
        if "tp" not in results[bam] or "fp" not in results[bam]:
            continue
        tp_c = results[bam]["tp"]
        fp_c = results[bam]["fp"]
        # Combined (TP+FP)
        hp1_alt = tp_c["hp1_alt"] + fp_c["hp1_alt"]
        hp1_ref = tp_c["hp1_ref"] + fp_c["hp1_ref"]
        hp2_alt = tp_c["hp2_alt"] + fp_c["hp2_alt"]
        hp2_ref = tp_c["hp2_ref"] + fp_c["hp2_ref"]
        table_combined = [[hp1_alt, hp1_ref], [hp2_alt, hp2_ref]]
        cv = cramers_v(table_combined)
        or_val = odds_ratio(table_combined)
        # HP-ratio for ALT-only and REF-only
        hp1_hp2_alt = hp1_alt / max(hp2_alt, 1)
        hp1_hp2_ref = hp1_ref / max(hp2_ref, 1)
        # Imbalance: difference of ALT-only ratio from REF-only ratio (should be ~0 if no priority bug)
        ratio_imbalance = abs(hp1_hp2_alt - hp1_hp2_ref) / max((hp1_hp2_alt + hp1_hp2_ref) / 2, 1e-9)
        bam_stats.append({
            "BAM": bam,
            "HP1_ALT": hp1_alt, "HP1_REF": hp1_ref,
            "HP2_ALT": hp2_alt, "HP2_REF": hp2_ref,
            "cramers_V": cv,
            "odds_ratio_HP1ALT_vs_HP2ALT_REF": or_val,
            "HP1_HP2_ratio_ALT_only": hp1_hp2_alt,
            "HP1_HP2_ratio_REF_only": hp1_hp2_ref,
            "imbalance_ALT_vs_REF": ratio_imbalance,
        })
        print(f"  {bam}: Cramer's V={cv:.4f}, OR={or_val:.3f}, "
              f"HP1:HP2 ALT={hp1_hp2_alt:.3f}, HP1:HP2 REF={hp1_hp2_ref:.3f}, "
              f"imbalance={ratio_imbalance:.3f}", flush=True)

    # Save TSVs
    raw_path = STEP1 / "step1_per_cpg_hp_alt_raw.tsv"
    with raw_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\n→ {raw_path}")

    stats_path = STEP1 / "step1_per_cpg_hp_alt_summary.tsv"
    with stats_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bam_stats[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(bam_stats)
    print(f"→ {stats_path}")

    # F12 figure
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['DejaVu Sans', 'Droid Sans Fallback']
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    bam_names = [s["BAM"] for s in bam_stats]
    colors = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"]

    # (a) Cramer's V
    ax = axes[0]
    cv_vals = [s["cramers_V"] for s in bam_stats]
    bars = ax.bar(bam_names, cv_vals, color=colors, alpha=0.85)
    ax.set_ylabel("Cramer's V (HP × ALT association)")
    ax.set_title("(a) HP × ALT 關聯強度\n(0 = independent / 1 = perfect association)")
    ax.grid(True, axis='y', alpha=0.3)
    for bar, v in zip(bars, cv_vals):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.001, f"{v:.4f}", ha='center', fontsize=9)

    # (b) HP1:HP2 ratio ALT-only vs REF-only
    ax = axes[1]
    x = np.arange(len(bam_names))
    width = 0.35
    alt_ratios = [s["HP1_HP2_ratio_ALT_only"] for s in bam_stats]
    ref_ratios = [s["HP1_HP2_ratio_REF_only"] for s in bam_stats]
    ax.bar(x - width/2, alt_ratios, width, label='ALT-only', color='#1f77b4', alpha=0.85)
    ax.bar(x + width/2, ref_ratios, width, label='REF-only', color='#aaaaaa', alpha=0.85)
    ax.axhline(1.0, color='#222', linestyle='--', linewidth=1, label='balanced 1:1')
    ax.set_xticks(x)
    ax.set_xticklabels(bam_names)
    ax.set_ylabel("HP1 : HP2 ratio")
    ax.set_title("(b) HP1:HP2 ratio split by ALT vs REF\n(priority bug = ALT 顯著偏離 REF)")
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_yscale('log')
    for i, (a, r) in enumerate(zip(alt_ratios, ref_ratios)):
        ax.text(i - width/2, a * 1.1, f"{a:.2f}", ha='center', fontsize=8)
        ax.text(i + width/2, r * 1.1, f"{r:.2f}", ha='center', fontsize=8)

    # (c) Imbalance (ALT vs REF difference)
    ax = axes[2]
    imb = [s["imbalance_ALT_vs_REF"] for s in bam_stats]
    bars = ax.bar(bam_names, imb, color=colors, alpha=0.85)
    ax.set_ylabel("|HP1:HP2(ALT) - HP1:HP2(REF)| / mean")
    ax.set_title("(c) Priority bug imbalance proxy\n(理想 = 0：ALT/REF 同一 HP1:HP2 ratio)")
    ax.grid(True, axis='y', alpha=0.3)
    for bar, v in zip(bars, imb):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.02, f"{v:.3f}", ha='center', fontsize=9)

    fig.suptitle("F12 — Goal 1 foundation: HP × ALT/REF association 跨 4 BAM (HCC1395 全基因組, flag=off)",
                 fontsize=12)
    plt.tight_layout()
    out_fig = FIG_DIR / "F12_per_cpg_hp_alt_4bam.png"
    plt.savefig(out_fig, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"→ {out_fig}")

    # Findings JSON
    findings = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "scope": "HCC1395 paired-pileup, 35,332 regions × 4 BAM × flag=off × {TP, FP}",
        "methodology": "Aggregated reads.tsv hp × alt_support counts; HP1 series = {1,1-1,11}, HP2 = {2,2-1,21}",
        "bam_stats": bam_stats,
        "verdict": {
            "expected_baseline_highest_cv": cv_vals[0] == max(cv_vals),
            "v3f_v6_balanced": bam_stats[1]["imbalance_ALT_vs_REF"] < 0.5 and bam_stats[3]["imbalance_ALT_vs_REF"] < 0.5,
        },
    }
    findings_path = STEP1 / "step1_per_cpg_hp_alt_findings.json"
    with findings_path.open("w") as fh:
        json.dump(findings, fh, indent=2)
    print(f"→ {findings_path}")


if __name__ == "__main__":
    main()
