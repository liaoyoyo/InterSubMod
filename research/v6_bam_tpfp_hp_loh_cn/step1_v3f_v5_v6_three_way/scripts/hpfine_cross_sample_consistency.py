#!/usr/bin/env python3
"""Goal 2 — HPFineNGroups cross-sample consistency (V6 across 5 samples).

Aggregate NG distribution per sample (HCC1395 + H1437 + H2009 + HCC1954 + HCC1937)
using V6 BAM × flag=off × {TP, FP} reads.tsv data.

For each sample × {TP, FP}, count #{regions with NG_off = k} for k in {0,1,2,3,4,5,6+}.

Output: per-sample NG distribution + cross-sample Spearman correlation matrix.

Caveat: only V6 ISM data exists for 4 non-HCC1395 samples (Phase D). V3F/V5/baseline
not run on these 4 samples (~50 min ISM × 4 samples × 3 BAMs ≈ 10 hr if added — out of scope).
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
from scipy import stats

REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PHASED_DIR = REPO / "research/paired_priority_bug_audit/phaseD_v6_5sample"
PHASEC_DIR = REPO / "research/paired_priority_bug_audit/phaseC_genome_three_way"  # HCC1395 V6
STEP1 = REPO / "research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way"
FIG_DIR = REPO / "research/v6_bam_tpfp_hp_loh_cn/figures/baseline_vs_v6"
FIG_DIR.mkdir(parents=True, exist_ok=True)

NG_VALID = {"1", "2", "1-1", "2-1", "3", "11", "21", "33"}


def count_ng_in_file(reads_tsv):
    """Return NG count = #{distinct hp values from NG_VALID with count >= 1}."""
    hp_set = set()
    try:
        with open(reads_tsv, "rb") as fh:
            header = fh.readline().decode().rstrip().split("\t")
            try:
                hp_i = header.index("hp")
            except ValueError:
                return 0
            for line in fh:
                cols = line.decode().rstrip().split("\t")
                if len(cols) > hp_i:
                    val = cols[hp_i]
                    if val in NG_VALID:
                        hp_set.add(val)
    except Exception:
        return 0
    return len(hp_set)


def find_reads_tsv(run_dir):
    out = []
    for root, _, files in os.walk(run_dir):
        if "reads.tsv" in files:
            out.append(os.path.join(root, "reads.tsv"))
    return out


def aggregate_sample(run_dir, n_workers=16):
    """Return Counter of NG values across all regions in this run."""
    files = find_reads_tsv(run_dir)
    counts = Counter()
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for ng in ex.map(count_ng_in_file, files, chunksize=200):
            counts[ng] += 1
    return dict(counts), len(files)


def main():
    print("Aggregating V6 NG distribution per sample × {TP, FP}", flush=True)

    # 5 samples × 2 labels
    samples = [
        ("HCC1395", PHASEC_DIR / "V6_off_tp", PHASEC_DIR / "V6_off_fp"),
        ("H1437", PHASED_DIR / "H1437/off_tp", PHASED_DIR / "H1437/off_fp"),
        ("H2009", PHASED_DIR / "H2009/off_tp", PHASED_DIR / "H2009/off_fp"),
        ("HCC1954", PHASED_DIR / "HCC1954/off_tp", PHASED_DIR / "HCC1954/off_fp"),
        ("HCC1937", PHASED_DIR / "HCC1937/off_tp", PHASED_DIR / "HCC1937/off_fp"),
    ]

    all_data = {}  # sample → {label: NG_counter}
    summary_rows = []
    for s, tp_dir, fp_dir in samples:
        if not tp_dir.exists() or not fp_dir.exists():
            print(f"  SKIP {s}: missing dirs", flush=True)
            continue
        t0 = time.time()
        tp_counts, n_tp = aggregate_sample(tp_dir, n_workers=16)
        fp_counts, n_fp = aggregate_sample(fp_dir, n_workers=16)
        elapsed = time.time() - t0
        all_data[s] = {"tp": tp_counts, "fp": fp_counts, "n_tp": n_tp, "n_fp": n_fp}
        print(f"  {s}: TP {n_tp} regions, FP {n_fp} regions  ({elapsed:.1f}s)", flush=True)
        # TP counts (NG=k → count)
        for ng_k in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
            summary_rows.append({
                "sample": s, "label": "TP", "NG_value": ng_k,
                "region_count": tp_counts.get(ng_k, 0),
                "total_regions": n_tp,
                "fraction": tp_counts.get(ng_k, 0) / max(n_tp, 1),
            })
            summary_rows.append({
                "sample": s, "label": "FP", "NG_value": ng_k,
                "region_count": fp_counts.get(ng_k, 0),
                "total_regions": n_fp,
                "fraction": fp_counts.get(ng_k, 0) / max(n_fp, 1),
            })

    raw_path = STEP1 / "step1_hpfine_cross_sample_raw.tsv"
    with raw_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "label", "NG_value",
                                            "region_count", "total_regions", "fraction"],
                            delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\n→ {raw_path}", flush=True)

    # Cross-sample Spearman on TP NG distribution (NG vector per sample)
    print(f"\n=== Cross-sample TP NG distribution Spearman ρ ===", flush=True)
    sample_names = list(all_data.keys())
    ng_vectors = {}
    for s in sample_names:
        tp = all_data[s]["tp"]
        # Vector: fractions at NG=0..8
        vec = [tp.get(k, 0) / max(all_data[s]["n_tp"], 1) for k in range(9)]
        ng_vectors[s] = vec

    rho_matrix = []
    for s1 in sample_names:
        row = []
        for s2 in sample_names:
            rho, _ = stats.spearmanr(ng_vectors[s1], ng_vectors[s2])
            row.append(rho)
        rho_matrix.append(row)
        print(f"  {s1}: " + ", ".join(f"{rho_matrix[-1][j]:.4f}" for j in range(len(sample_names))))

    # Wilcoxon: NG=3 vs NG=2 ratio (canonical filter zone)
    print(f"\n=== NG distribution focus: NG=3 vs NG=2 ratio (TP) ===", flush=True)
    ratio_rows = []
    for s in sample_names:
        tp = all_data[s]["tp"]
        ng2 = tp.get(2, 0)
        ng3 = tp.get(3, 0)
        ratio = ng3 / max(ng2, 1)
        ratio_rows.append({"sample": s, "TP_NG2": ng2, "TP_NG3": ng3, "NG3_NG2_ratio": ratio})
        print(f"  {s}: NG=2={ng2}, NG=3={ng3}, ratio={ratio:.3f}", flush=True)

    # Summary
    summary_path = STEP1 / "step1_hpfine_cross_sample_summary.tsv"
    with summary_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ratio_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(ratio_rows)
    print(f"→ {summary_path}", flush=True)

    # F13 figure — NG distribution per sample + Spearman heatmap
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['DejaVu Sans', 'Droid Sans Fallback']
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # (a) NG distribution per sample (TP)
    ax = axes[0]
    colors_5 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    x = list(range(9))
    width = 0.15
    for i, s in enumerate(sample_names):
        offsets = x_vals = [k + (i - 2) * width for k in x]
        ax.bar(offsets, ng_vectors[s], width, label=s, color=colors_5[i], alpha=0.85)
    ax.set_xlabel("NG_off value")
    ax.set_ylabel("Region fraction (TP regions)")
    ax.set_title("F13a — V6 NG distribution per sample (TP regions)")
    ax.set_xticks(x)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)

    # (b) Spearman heatmap
    ax = axes[1]
    rho_arr = np.array(rho_matrix)
    im = ax.imshow(rho_arr, cmap='RdYlBu_r', vmin=0.0, vmax=1.0, aspect='auto')
    ax.set_xticks(range(len(sample_names)))
    ax.set_yticks(range(len(sample_names)))
    ax.set_xticklabels(sample_names, rotation=45, ha='right')
    ax.set_yticklabels(sample_names)
    ax.set_title("F13b — Cross-sample Spearman ρ on NG distribution\n(higher = more consistent ISM behavior)")
    for i in range(len(sample_names)):
        for j in range(len(sample_names)):
            ax.text(j, i, f"{rho_arr[i,j]:.3f}", ha='center', va='center',
                    color='black' if abs(rho_arr[i,j]) < 0.5 else 'white', fontsize=9)
    plt.colorbar(im, ax=ax, label='Spearman ρ')

    fig.suptitle("F13 — Goal 2: V6 HPFineN cross-sample consistency (HCC1395 + Phase D 4 samples)",
                 fontsize=12)
    plt.tight_layout()
    out_fig = FIG_DIR / "F13_hpfine_cross_sample.png"
    plt.savefig(out_fig, dpi=130, bbox_inches='tight')
    plt.close()
    print(f"→ {out_fig}")

    # Findings JSON
    avg_offdiag_rho = float(np.mean([rho_arr[i, j] for i in range(len(sample_names))
                                      for j in range(len(sample_names)) if i != j]))
    findings = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "samples": sample_names,
        "ng_distributions_TP_fractions": ng_vectors,
        "spearman_matrix": rho_arr.tolist(),
        "avg_offdiagonal_rho": avg_offdiag_rho,
        "ng3_ng2_ratios": ratio_rows,
        "verdict": (
            "PASS — Spearman ρ avg > 0.7 across samples" if avg_offdiag_rho > 0.7
            else "PARTIAL — avg < 0.7, some samples diverge"
        ),
        "caveat": "V3F/V5/baseline 在 4 個非 HCC1395 樣本未跑 (Phase D V6 only)；計算成本 ~10hr 補完成本未在 scope",
    }
    findings_path = STEP1 / "step1_hpfine_cross_sample_findings.json"
    with findings_path.open("w") as fh:
        json.dump(findings, fh, indent=2, default=str)
    print(f"→ {findings_path}")
    print(f"\n  Avg off-diagonal Spearman ρ = {avg_offdiag_rho:.4f}")
    print(f"  Verdict: {findings['verdict']}")


if __name__ == "__main__":
    main()
