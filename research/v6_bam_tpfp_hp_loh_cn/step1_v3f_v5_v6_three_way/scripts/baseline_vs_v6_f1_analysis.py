#!/usr/bin/env python3
"""baseline vs V6 F1 analysis — can ISM filter beat caller-alone F1?

Reads step1_master_four_way.tsv (35,332 caller PASS calls).
Treats master's TP/FP labels (vs SEQC2 truth) as ground truth.
Computes F1 of "keep NG_off >= T" filter strategy at multiple T thresholds for each BAM.

KEY METRICS:
1. F1 vs caller-alone (no filter applied — all 35,332 PASS calls kept)
   - caller-alone F1 on pileup subset = 2 × 30,490 / (2 × 30,490 + 4,842 + 0) = 0.9265
2. F1 of ISM marker filter at NG_off ∈ {1, 2, 3, 4, 5}
3. Precision-Recall curve per BAM
4. F1 stratified by LOH zone (LOH-positive vs LOH-NA)
5. Δ F1 (V6 - baseline) at each threshold

Outputs:
- step1_f1_analysis_thresholds.tsv — F1/P/R per (BAM, NG threshold)
- step1_f1_analysis_per_zone.tsv — F1 per (BAM, NG threshold, LOH zone)
- step1_f1_analysis_summary.json — best thresholds + caller-headroom
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

STEP1 = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way")
MASTER_PATH = STEP1 / "step1_master_four_way.tsv"


def safe_int(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return 0


def f1_pr(tp, fp, fn):
    p = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    r = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0.0
    return f1, p, r


def main():
    if not MASTER_PATH.exists():
        print(f"ERROR: {MASTER_PATH} not found", file=sys.stderr)
        return 1

    bams = ["baseline", "V3F", "V5", "V6"]
    thresholds = [1, 2, 3, 4, 5]

    # FN_CALLER_PILEUP = truth positives in pileup subset NOT called by caller (caller missed)
    # Source: InterSubMod/research/methyl_augmented_filter_phase2/cycle1/scripts/final_filter_and_verdict.py:44
    # Canonical caller-alone F1 on pileup subset = 2*30490 / (2*30490 + 4842 + 19288) = 0.7166 (matches V6 doc §8.6)
    FN_CALLER_PILEUP = 19288

    # Counters
    # at_ng[bam][T][label] = count of regions where NG_off >= T
    at_ng = {b: {T: {"TP": 0, "FP": 0} for T in thresholds} for b in bams}
    total_label = {"TP": 0, "FP": 0}

    # zone-stratified
    # at_ng_zone[bam][T][zone][label]
    at_ng_zone = {b: {T: defaultdict(lambda: {"TP": 0, "FP": 0}) for T in thresholds} for b in bams}
    total_zone = defaultdict(lambda: {"TP": 0, "FP": 0})

    with MASTER_PATH.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            label = row["label"]
            loh_annot = row.get("LOH_Bed_Annotation", "") or ""
            zone = "LOH-positive" if loh_annot else "LOH-NA"

            total_label[label] += 1
            total_zone[zone][label] += 1

            for b in bams:
                ng = safe_int(row.get(f"{b}_off_NG", "0"))
                for T in thresholds:
                    if ng >= T:
                        at_ng[b][T][label] += 1
                        at_ng_zone[b][T][zone][label] += 1

    TP_TOTAL = total_label["TP"]
    FP_TOTAL = total_label["FP"]
    print(f"Master subset: TP_total={TP_TOTAL} | FP_total={FP_TOTAL}")

    # Caller-alone F1 on pileup subset (no filter applied, but FN_CALLER_PILEUP counted)
    f1_caller_alone, p_caller, r_caller = f1_pr(TP_TOTAL, FP_TOTAL, FN_CALLER_PILEUP)
    print(f"Caller-alone F1 (pileup subset, no filter): {f1_caller_alone:.4f} (P={p_caller:.4f}, R={r_caller:.4f})")
    print(f"  Canonical caller F1 vs SEQC2 truth = 0.7166 (V6 doc §8.6)")
    sanity = abs(f1_caller_alone - 0.7166) < 0.001
    print(f"  Sanity: {'PASS ✓' if sanity else 'FAIL ✗ (formula bug?)'}")
    print()

    # ─── F1 sweep across NG thresholds ───
    rows = []
    print("=== F1 sweep: ISM NG_off >= T filter ===")
    print(f"{'BAM':<10} {'T':>3} {'TP':>7} {'FP':>6} {'FN':>7} {'P':>7} {'R':>7} {'F1':>7} {'ΔF1(vs caller-alone)':>20}")
    for T in thresholds:
        for b in bams:
            tp = at_ng[b][T]["TP"]
            fp = at_ng[b][T]["FP"]
            fn = FN_CALLER_PILEUP + (TP_TOTAL - tp)  # caller-missed + filter-rejected caller TPs
            f1, p, r = f1_pr(tp, fp, fn)
            delta = f1 - f1_caller_alone
            print(f"{b:<10} {T:>3} {tp:>7} {fp:>6} {fn:>7} {p:>7.4f} {r:>7.4f} {f1:>7.4f} {delta:>+20.4f}")
            rows.append({
                "BAM": b, "NG_threshold": T,
                "TP": tp, "FP": fp, "FN": fn,
                "precision": f"{p:.4f}", "recall": f"{r:.4f}", "F1": f"{f1:.4f}",
                "delta_F1_vs_caller_alone": f"{delta:+.4f}",
            })

    tsv_path = STEP1 / "step1_f1_analysis_thresholds.tsv"
    with tsv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\n  → {tsv_path}")

    # ─── Per-zone F1 ───
    zone_rows = []
    print("\n=== F1 stratified by LOH zone (NG_off >= 3) ===")
    print(f"{'BAM':<10} {'Zone':<14} {'TP':>6} {'FP':>5} {'FN':>6} {'P':>7} {'R':>7} {'F1':>7}")
    # Apportion FN_CALLER_PILEUP across zones proportionally to TP_zone
    fn_caller_zone = {}
    for zone in ["LOH-positive", "LOH-NA"]:
        TP_zone = total_zone[zone]["TP"]
        fn_caller_zone[zone] = round(FN_CALLER_PILEUP * TP_zone / TP_TOTAL)
    print(f"  FN_CALLER_PILEUP apportioned: LOH-positive={fn_caller_zone['LOH-positive']} / LOH-NA={fn_caller_zone['LOH-NA']} (sum={sum(fn_caller_zone.values())})")

    for zone in ["LOH-positive", "LOH-NA"]:
        TP_zone = total_zone[zone]["TP"]
        FP_zone = total_zone[zone]["FP"]
        f1_zone_caller, _, _ = f1_pr(TP_zone, FP_zone, fn_caller_zone[zone])
        print(f"  Zone {zone}: TP_total={TP_zone}, FP_total={FP_zone}, FN_caller_zone={fn_caller_zone[zone]}, caller-alone F1={f1_zone_caller:.4f}")
        for T in [2, 3]:
            for b in bams:
                tp = at_ng_zone[b][T][zone]["TP"]
                fp = at_ng_zone[b][T][zone]["FP"]
                fn = fn_caller_zone[zone] + (TP_zone - tp)
                f1, p, r = f1_pr(tp, fp, fn)
                print(f"{b:<10} {zone:<14} {tp:>6} {fp:>5} {fn:>6} {p:>7.4f} {r:>7.4f} {f1:>7.4f} (T={T})")
                zone_rows.append({
                    "BAM": b, "NG_threshold": T, "LOH_zone": zone,
                    "TP": tp, "FP": fp, "FN": fn,
                    "precision": f"{p:.4f}", "recall": f"{r:.4f}", "F1": f"{f1:.4f}",
                    "F1_caller_zone": f"{f1_zone_caller:.4f}",
                    "delta_F1_zone": f"{f1 - f1_zone_caller:+.4f}",
                })

    zone_tsv = STEP1 / "step1_f1_analysis_per_zone.tsv"
    with zone_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(zone_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(zone_rows)
    print(f"\n  → {zone_tsv}")

    # ─── Δ F1 V6 - baseline ───
    print("\n=== Δ F1 (V6 - baseline) sweep ===")
    delta_v6_base = []
    for T in thresholds:
        tp_v6, fp_v6 = at_ng["V6"][T]["TP"], at_ng["V6"][T]["FP"]
        tp_b, fp_b = at_ng["baseline"][T]["TP"], at_ng["baseline"][T]["FP"]
        f1_v6, _, _ = f1_pr(tp_v6, fp_v6, FN_CALLER_PILEUP + TP_TOTAL - tp_v6)
        f1_b, _, _ = f1_pr(tp_b, fp_b, FN_CALLER_PILEUP + TP_TOTAL - tp_b)
        delta = f1_v6 - f1_b
        delta_v6_base.append({"NG_threshold": T, "f1_v6": f1_v6, "f1_baseline": f1_b, "delta": delta})
        print(f"  T={T}: F1_v6={f1_v6:.4f} - F1_baseline={f1_b:.4f} = Δ {delta:+.4f}")

    # ─── Best F1 per BAM ───
    print("\n=== Best F1 per BAM (across thresholds) ===")
    best_per_bam = {}
    for b in bams:
        best_T = None
        best_f1 = -1
        for T in thresholds:
            tp = at_ng[b][T]["TP"]
            fp = at_ng[b][T]["FP"]
            fn = FN_CALLER_PILEUP + (TP_TOTAL - tp)
            f1, _, _ = f1_pr(tp, fp, fn)
            if f1 > best_f1:
                best_f1 = f1
                best_T = T
        best_per_bam[b] = {"best_NG_threshold": best_T, "best_F1": best_f1}
        print(f"  {b}: best F1 = {best_f1:.4f} at T={best_T}")
    print(f"  caller-alone F1 = {f1_caller_alone:.4f} (no filter)")

    # ─── Findings JSON ───
    findings = {
        "data_source": str(MASTER_PATH),
        "tp_total_pileup_subset": TP_TOTAL,
        "fp_total_pileup_subset": FP_TOTAL,
        "caller_alone_F1_pileup_subset": f1_caller_alone,
        "canonical_caller_F1_vs_seqc2_truth": 0.7166,
        "caller_alone_F1_at_purity_0_6": 0.6273,
        "best_per_bam": best_per_bam,
        "delta_v6_minus_baseline_by_threshold": delta_v6_base,
        "verdict": {
            "ISM_NG_filter_beats_caller_alone_on_pileup_subset": any(
                best_per_bam[b]["best_F1"] > f1_caller_alone for b in bams
            ),
            "V6_dominates_baseline_at_all_thresholds": all(
                d["delta"] > 0 for d in delta_v6_base
            ),
            "note": (
                "F1 on pileup subset (master labels) is high because subset excludes caller FN. "
                "Canonical caller F1 vs SEQC2 truth = 0.7166 (includes 10,938 FN). "
                "To improve canonical F1, would need smarter filter (multi-axis LR like Phase 2 Cycle 1 +0.02236)."
            ),
        },
    }
    findings_path = STEP1 / "step1_f1_analysis_summary.json"
    with findings_path.open("w") as fh:
        json.dump(findings, fh, indent=2)
    print(f"\n  → {findings_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
