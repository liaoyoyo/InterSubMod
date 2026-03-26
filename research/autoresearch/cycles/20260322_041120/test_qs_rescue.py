#!/usr/bin/env python3
"""
Cycle 20260322_041120: Quality_Score rescue rules analysis
Tests QS >= 50 alone (no GQ gate) vs various combinations
Dataset: HCC1395 5kHz TO (pilot)
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import numpy as np

# Baseline metrics (from dataset_overview.tsv)
BASELINE_TP = 28396
BASELINE_FP = 11843
TRUTH_TOTAL = 39447  # TP + FN = 28396 + 11051
BASELINE_F1 = 0.712697

DATA_PATH = Path("/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
                 "20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv")

DORADO_DATA_PATH = Path("/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
                         "20260311_gq_methylation_rescue_matrix_with_dorado_to/"
                         "gq_plus_methylation_rule_sweep.tsv")


def compute_f1(tp: int, fp: int, truth_total: int):
    fn = truth_total - tp
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / truth_total if truth_total > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}


def simulate_rescue(tp_df, fp_df, tp_mask, fp_mask, label: str, baseline_tp=BASELINE_TP,
                    baseline_fp=BASELINE_FP, truth_total=TRUTH_TOTAL):
    rescued_tp = tp_mask.sum()
    rescued_fp = fp_mask.sum()
    new_tp = baseline_tp + rescued_tp
    new_fp = baseline_fp + rescued_fp
    metrics = compute_f1(new_tp, new_fp, truth_total)
    delta = metrics["f1"] - BASELINE_F1
    return {
        "label": label,
        "rescued_tp": int(rescued_tp),
        "rescued_fp": int(rescued_fp),
        "delta_f1": delta,
        "result_f1": metrics["f1"],
        "precision": metrics["precision"],
        "recall": metrics["recall"],
    }


def main():
    print("=" * 70)
    print("Cycle 20260322_041120: QS Rescue Rule Analysis")
    print(f"Dataset: HCC1395 5kHz TO (pilot)")
    print(f"Baseline: F1={BASELINE_F1:.6f}, TP={BASELINE_TP}, FP={BASELINE_FP}")
    print("=" * 70)

    df = pd.read_csv(DATA_PATH, sep="\t")
    df = df.reset_index(drop=True)
    tp = df[df["truth_status"] == "TP"].copy()
    fp = df[df["truth_status"] == "FP"].copy()
    print(f"\nCandidate pool: TP={len(tp)}, FP={len(fp)} (ISM-analyzed rescue candidates)")
    print(f"Note: These are currently FN (TP) or filtered-out FP in the baseline")

    print("\n--- QS-based Rescue Rules (Pilot: 5kHz TO) ---")
    rules = [
        ("QS>=40 alone", tp["Quality_Score"] >= 40, fp["Quality_Score"] >= 40),
        ("QS>=50 alone [H004]", tp["Quality_Score"] >= 50, fp["Quality_Score"] >= 50),
        ("QS>=60 alone", tp["Quality_Score"] >= 60, fp["Quality_Score"] >= 60),
        ("QS>=70 alone", tp["Quality_Score"] >= 70, fp["Quality_Score"] >= 70),
        ("QS>=75 alone", tp["Quality_Score"] >= 75, fp["Quality_Score"] >= 75),
        ("GQ>=10 alone", tp["gq"] >= 10, fp["gq"] >= 10),
        ("GQ>=10+QS>=60 [prev best]", (tp["gq"] >= 10) & (tp["Quality_Score"] >= 60),
                                       (fp["gq"] >= 10) & (fp["Quality_Score"] >= 60)),
        ("QS>=50+PMD>=0.15", (tp["Quality_Score"] >= 50) & (tp["PairwiseMedianDist"] >= 0.15),
                              (fp["Quality_Score"] >= 50) & (fp["PairwiseMedianDist"] >= 0.15)),
        ("QS>=50+PMD>=0.20", (tp["Quality_Score"] >= 50) & (tp["PairwiseMedianDist"] >= 0.20),
                              (fp["Quality_Score"] >= 50) & (fp["PairwiseMedianDist"] >= 0.20)),
        ("hp_assign<0.5+QS>=60 EXCLUDE",  # filter: don't rescue if hp_assign low
         (tp["Quality_Score"] >= 60) & (tp["hp_assign_rate"] >= 0.5),
         (fp["Quality_Score"] >= 60) & (fp["hp_assign_rate"] >= 0.5)),
    ]

    results = []
    for label, tp_mask, fp_mask in rules:
        r = simulate_rescue(tp, fp, tp_mask, fp_mask, label)
        results.append(r)
        marker = " <<<" if r["delta_f1"] > 0.007 else ("" if r["delta_f1"] <= 0 else "")
        print(f"  {label:<40s}: TP={r['rescued_tp']:3d}, FP={r['rescued_fp']:3d}, "
              f"delta={r['delta_f1']:+.6f}, F1={r['result_f1']:.6f}{marker}")

    best = max(results, key=lambda x: x["delta_f1"])
    print(f"\n[BEST] {best['label']}: delta_F1={best['delta_f1']:+.6f}, F1={best['result_f1']:.6f}")

    # Analyze QS thresholds more finely
    print("\n--- Fine-grained QS sweep (no GQ gate) ---")
    for qs in range(40, 85, 5):
        tp_cnt = (tp["Quality_Score"] >= qs).sum()
        fp_cnt = (fp["Quality_Score"] >= qs).sum()
        new_tp = BASELINE_TP + tp_cnt
        new_fp = BASELINE_FP + fp_cnt
        prec = new_tp / (new_tp + new_fp)
        rec = new_tp / TRUTH_TOTAL
        f1 = 2 * prec * rec / (prec + rec)
        delta = f1 - BASELINE_F1
        print(f"  QS>={qs}: rescued TP={tp_cnt}, FP={fp_cnt}, delta={delta:+.6f}, F1={f1:.6f}")

    # H009 observation: CramersV > 0.3 precision in TO candidates
    print("\n--- H009 Observational: CramersV thresholds in TO candidates ---")
    for cv in [0.1, 0.2, 0.3, 0.4, 0.5]:
        cv_subset = df[df["CramersV"] >= cv]
        if len(cv_subset) > 0:
            prec = (cv_subset["truth_status"] == "TP").sum() / len(cv_subset)
            tp_cnt = (cv_subset["truth_status"] == "TP").sum()
            fp_cnt = (cv_subset["truth_status"] == "FP").sum()
            print(f"  CramersV>={cv}: N={len(cv_subset)}, TP={tp_cnt}, FP={fp_cnt}, precision={prec*100:.1f}%")

    # DORADO TO comparison (from existing sweep)
    print("\n--- DORADO TO comparison (from existing sweep) ---")
    try:
        dorado_df = pd.read_csv(DORADO_DATA_PATH, sep="\t")
        dorado_to = dorado_df[dorado_df["dataset_id"] == "hcc1395_dorado_to"]
        best_dorado = dorado_to.sort_values("delta_f1_vs_baseline", ascending=False).iloc[0]
        print(f"  DORADO TO best rule: {best_dorado['rule_id']}")
        print(f"  delta_F1 vs baseline: {best_dorado['delta_f1_vs_baseline']:+.6f}")
        print(f"  TP rescued: {best_dorado['tp_rescued']}, FP reintroduced: {best_dorado['fp_reintroduced']}")
        print(f"  DORADO TO baseline F1: {best_dorado['f1'] - best_dorado['delta_f1_vs_baseline']:.6f}")
    except Exception as e:
        print(f"  [DORADO data unavailable: {e}]")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Best rule for 5kHz TO: {best['label']}")
    print(f"  Predicted delta_F1: {best['delta_f1']:+.6f}")
    print(f"  Predicted F1: {best['result_f1']:.6f}")
    print(f"  TP rescued: {best['rescued_tp']}, FP reintroduced: {best['rescued_fp']}")
    print(f"  Precision after rescue: {(BASELINE_TP + best['rescued_tp']) / (BASELINE_TP + best['rescued_tp'] + BASELINE_FP + best['rescued_fp']):.4f}")

    # Write result summary
    result = {
        "cycle_id": "20260322_041120",
        "hypothesis": "QS>=50 as standalone rescue rule (no GQ gate)",
        "dataset": "HCC1395_5kHz_TO",
        "baseline_f1": BASELINE_F1,
        "best_rule": best["label"],
        "best_delta_f1": best["delta_f1"],
        "best_result_f1": best["result_f1"],
        "best_rescued_tp": best["rescued_tp"],
        "best_rescued_fp": best["rescued_fp"],
        "all_rules": results,
    }
    out_path = Path(__file__).parent / "result_summary.json"
    with open(out_path, "w") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
