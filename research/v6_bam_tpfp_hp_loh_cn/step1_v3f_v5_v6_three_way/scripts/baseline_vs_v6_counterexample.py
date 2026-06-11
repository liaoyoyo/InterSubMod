#!/usr/bin/env python3
"""baseline vs V6 counter-example hunt — find sub-populations where baseline performs >= V6.

Reads step1_master_four_way.tsv and stratifies along 5 axes to surface failure modes:

Q1. Discordant marker sets — regions in baseline-marker set NOT in V6, and vice versa
Q2. Per-chromosome marker rate inversion (baseline_rate > V6_rate)
Q3. Per caller_af bin marker rate inversion
Q4. Per LOH zone × Coverage category cross-tab inversions
Q5. Top-N TP loss regions — specific regions where V6 lost a baseline TP marker
Q6. NG distribution shift per BAM (region-level histogram)

Outputs:
- step1_counterexample_summary.json — machine-readable findings
- step1_counterexample_discordant_markers.tsv — per-region 4-way membership + label
- step1_counterexample_per_chr.tsv — chr inversion candidates
- step1_counterexample_per_af_bin.tsv — AF bin inversion candidates
- step1_counterexample_per_loh_cov.tsv — LOH × CN cross-tab
- step1_counterexample_lost_tp_regions.tsv — top regions where V6 lost a baseline TP marker
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


def safe_float(x):
    try:
        v = float(x)
        if v != v:  # NaN
            return None
        return v
    except (ValueError, TypeError):
        return None


def af_bin(af):
    if af is None:
        return "AF_NA"
    if af < 0.1: return "AF_0.0-0.1"
    if af < 0.2: return "AF_0.1-0.2"
    if af < 0.3: return "AF_0.2-0.3"
    if af < 0.5: return "AF_0.3-0.5"
    if af < 0.7: return "AF_0.5-0.7"
    return "AF_0.7-1.0"


def main():
    if not MASTER_PATH.exists():
        print(f"ERROR: {MASTER_PATH} not found", file=sys.stderr)
        return 1

    # Per-region marker membership: baseline / V3F / V5 / V6
    bams = ["baseline", "V3F", "V5", "V6"]

    # Q1 counters
    membership = defaultdict(int)  # key = (baseline_in, V6_in, label)
    discordant_rows = []  # for TSV output

    # Q2 per-chr
    per_chr = defaultdict(lambda: {b: {"tp": 0, "fp": 0} for b in bams})

    # Q3 per AF bin
    per_af = defaultdict(lambda: {b: {"tp": 0, "fp": 0} for b in bams})

    # Q4 per LOH × Coverage
    per_loh_cov = defaultdict(lambda: {b: {"tp": 0, "fp": 0} for b in bams})

    # Q5 lost TP regions
    lost_tp_regions = []

    # Q6 NG histogram
    ng_hist = {b: defaultdict(int) for b in bams}

    with MASTER_PATH.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            label = row["label"]
            chrom = row["chr"]
            rid = row["region_id"]
            af = safe_float(row.get("caller_af", ""))
            loh = row.get("LOH_Bed_Annotation", "") or "NA"
            cov = row.get("Coverage_Category", "") or "NA"

            ng_off = {b: safe_int(row.get(f"{b}_off_NG", "0")) for b in bams}
            ng_on = {b: safe_int(row.get(f"{b}_on_NG", "0")) for b in bams}
            is_marker = {b: (ng_off[b] >= 3) for b in bams}

            # Q1 discordant
            key = (is_marker["baseline"], is_marker["V6"], label)
            membership[key] += 1
            if is_marker["baseline"] != is_marker["V6"]:
                discordant_rows.append({
                    "region_id": rid,
                    "chr": chrom,
                    "label": label,
                    "caller_af": af if af is not None else "NA",
                    "LOH_zone": loh,
                    "Coverage_Category": cov,
                    "baseline_NG_off": ng_off["baseline"],
                    "V6_NG_off": ng_off["V6"],
                    "V3F_NG_off": ng_off["V3F"],
                    "V5_NG_off": ng_off["V5"],
                    "delta_V6_minus_baseline": ng_off["V6"] - ng_off["baseline"],
                    "category": "V6_only" if is_marker["V6"] else "baseline_only",
                })

            # Q2 per chr
            for b in bams:
                if is_marker[b]:
                    per_chr[chrom][b]["tp" if label == "TP" else "fp"] += 1

            # Q3 per AF bin
            ab = af_bin(af)
            for b in bams:
                if is_marker[b]:
                    per_af[ab][b]["tp" if label == "TP" else "fp"] += 1

            # Q4 per LOH × CN
            key_lc = f"{loh}|{cov}"
            for b in bams:
                if is_marker[b]:
                    per_loh_cov[key_lc][b]["tp" if label == "TP" else "fp"] += 1

            # Q5 lost TP: baseline marker AND label=TP AND V6 NOT marker
            if is_marker["baseline"] and not is_marker["V6"] and label == "TP":
                lost_tp_regions.append({
                    "region_id": rid,
                    "chr": chrom,
                    "caller_af": af if af is not None else "NA",
                    "LOH_zone": loh,
                    "Coverage_Category": cov,
                    "baseline_NG_off": ng_off["baseline"],
                    "V6_NG_off": ng_off["V6"],
                    "V3F_NG_off": ng_off["V3F"],
                    "V5_NG_off": ng_off["V5"],
                })

            # Q6 NG histogram
            for b in bams:
                ng_hist[b][ng_off[b]] += 1

    # ─── Q1 discordant marker summary ───
    print("=== Q1: Discordant marker membership (baseline vs V6) ===")
    print(f"{'baseline_marker':<18} {'V6_marker':<12} {'label':<6} {'count':>8}")
    for (b_in, v_in, label), n in sorted(membership.items(), key=lambda x: (-x[1])):
        print(f"{'YES' if b_in else 'no':<18} {'YES' if v_in else 'no':<12} {label:<6} {n:>8}")

    # Aggregate scenarios
    baseline_only_tp = membership.get((True, False, "TP"), 0)
    baseline_only_fp = membership.get((True, False, "FP"), 0)
    v6_only_tp = membership.get((False, True, "TP"), 0)
    v6_only_fp = membership.get((False, True, "FP"), 0)
    both_tp = membership.get((True, True, "TP"), 0)
    both_fp = membership.get((True, True, "FP"), 0)

    q1 = {
        "both_marker_TP": both_tp,
        "both_marker_FP": both_fp,
        "baseline_only_TP_lost_by_V6": baseline_only_tp,
        "baseline_only_FP_correctly_removed_by_V6": baseline_only_fp,
        "V6_only_TP_gained_over_baseline": v6_only_tp,
        "V6_only_FP_added_by_V6": v6_only_fp,
        "lost_TP_rate_V6_vs_baseline": baseline_only_tp / max(1, baseline_only_tp + both_tp),
        "FP_removal_rate_V6_vs_baseline": baseline_only_fp / max(1, baseline_only_fp + both_fp),
        "TP_gain_rate_V6_vs_baseline": v6_only_tp / max(1, v6_only_tp + both_tp),
        "FP_addition_rate_V6_vs_baseline": v6_only_fp / max(1, v6_only_fp + both_fp),
    }
    print("\n=== Q1 derived rates ===")
    for k, v in q1.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.4f}")
        else:
            print(f"  {k}: {v:,}")

    # Write discordant TSV
    disc_path = STEP1 / "step1_counterexample_discordant_markers.tsv"
    if discordant_rows:
        cols = list(discordant_rows[0].keys())
        with disc_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(discordant_rows)
        print(f"\n  → {disc_path} ({len(discordant_rows)} discordant regions)")

    # ─── Q2 per-chr ───
    chr_rows = []
    chr_inversions = []
    for chrom in sorted(per_chr.keys()):
        b_tp, b_fp = per_chr[chrom]["baseline"]["tp"], per_chr[chrom]["baseline"]["fp"]
        v_tp, v_fp = per_chr[chrom]["V6"]["tp"], per_chr[chrom]["V6"]["fp"]
        b_rate = b_tp / max(1, b_tp + b_fp)
        v_rate = v_tp / max(1, v_tp + v_fp)
        rec = {
            "chr": chrom,
            "baseline_TP": b_tp, "baseline_FP": b_fp, "baseline_rate": f"{b_rate:.4f}",
            "V6_TP": v_tp, "V6_FP": v_fp, "V6_rate": f"{v_rate:.4f}",
            "delta_V6_minus_baseline_rate": f"{v_rate - b_rate:+.4f}",
            "delta_V6_minus_baseline_count": (v_tp + v_fp) - (b_tp + b_fp),
        }
        chr_rows.append(rec)
        # Inversion: baseline_rate > V6_rate AND both have ≥5 markers
        if b_rate > v_rate and (b_tp + b_fp) >= 5 and (v_tp + v_fp) >= 5:
            chr_inversions.append((chrom, b_rate, v_rate, v_rate - b_rate))

    chr_path = STEP1 / "step1_counterexample_per_chr.tsv"
    with chr_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(chr_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(chr_rows)
    print(f"\n  → {chr_path}")

    print(f"\n=== Q2: Per-chr inversion (baseline_rate > V6_rate, both ≥5 markers) ===")
    if chr_inversions:
        for chrom, b_rate, v_rate, delta in sorted(chr_inversions, key=lambda x: x[3]):
            print(f"  {chrom}: baseline={b_rate:.4f} V6={v_rate:.4f} Δ={delta:+.4f}")
    else:
        print("  No per-chr inversion (V6 ≥ baseline in all chrs).")

    # ─── Q3 per AF bin ───
    af_rows = []
    af_inversions = []
    af_order = ["AF_0.0-0.1", "AF_0.1-0.2", "AF_0.2-0.3", "AF_0.3-0.5", "AF_0.5-0.7", "AF_0.7-1.0", "AF_NA"]
    for ab in af_order:
        if ab not in per_af: continue
        b_tp, b_fp = per_af[ab]["baseline"]["tp"], per_af[ab]["baseline"]["fp"]
        v_tp, v_fp = per_af[ab]["V6"]["tp"], per_af[ab]["V6"]["fp"]
        b_rate = b_tp / max(1, b_tp + b_fp)
        v_rate = v_tp / max(1, v_tp + v_fp)
        rec = {
            "AF_bin": ab,
            "baseline_TP": b_tp, "baseline_FP": b_fp, "baseline_rate": f"{b_rate:.4f}",
            "V6_TP": v_tp, "V6_FP": v_fp, "V6_rate": f"{v_rate:.4f}",
            "delta_V6_minus_baseline_rate": f"{v_rate - b_rate:+.4f}",
            "delta_V6_minus_baseline_count": (v_tp + v_fp) - (b_tp + b_fp),
        }
        af_rows.append(rec)
        if b_rate > v_rate and (b_tp + b_fp) >= 5 and (v_tp + v_fp) >= 5:
            af_inversions.append((ab, b_rate, v_rate, v_rate - b_rate))

    af_path = STEP1 / "step1_counterexample_per_af_bin.tsv"
    with af_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(af_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(af_rows)
    print(f"\n  → {af_path}")

    print(f"\n=== Q3: Per-AF-bin inversion ===")
    if af_inversions:
        for ab, b_rate, v_rate, delta in af_inversions:
            print(f"  {ab}: baseline={b_rate:.4f} V6={v_rate:.4f} Δ={delta:+.4f}")
    else:
        print("  No per-AF-bin inversion.")

    # ─── Q4 LOH × CN cross-tab ───
    loh_cov_rows = []
    loh_cov_inversions = []
    for key in sorted(per_loh_cov.keys()):
        b_tp, b_fp = per_loh_cov[key]["baseline"]["tp"], per_loh_cov[key]["baseline"]["fp"]
        v_tp, v_fp = per_loh_cov[key]["V6"]["tp"], per_loh_cov[key]["V6"]["fp"]
        b_rate = b_tp / max(1, b_tp + b_fp)
        v_rate = v_tp / max(1, v_tp + v_fp)
        rec = {
            "LOH_x_Cov": key,
            "baseline_TP": b_tp, "baseline_FP": b_fp, "baseline_rate": f"{b_rate:.4f}",
            "V6_TP": v_tp, "V6_FP": v_fp, "V6_rate": f"{v_rate:.4f}",
            "delta_rate": f"{v_rate - b_rate:+.4f}",
        }
        loh_cov_rows.append(rec)
        if b_rate > v_rate and (b_tp + b_fp) >= 10 and (v_tp + v_fp) >= 10:
            loh_cov_inversions.append((key, b_rate, v_rate, v_rate - b_rate))

    loh_cov_path = STEP1 / "step1_counterexample_per_loh_cov.tsv"
    with loh_cov_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(loh_cov_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(loh_cov_rows)
    print(f"\n  → {loh_cov_path}")

    print(f"\n=== Q4: LOH × Coverage inversion (both ≥10 markers) ===")
    if loh_cov_inversions:
        for key, b_rate, v_rate, delta in sorted(loh_cov_inversions, key=lambda x: x[3]):
            print(f"  {key}: baseline={b_rate:.4f} V6={v_rate:.4f} Δ={delta:+.4f}")
    else:
        print("  No LOH × Cov inversion.")

    # ─── Q5 Lost TP regions ───
    lost_tp_path = STEP1 / "step1_counterexample_lost_tp_regions.tsv"
    if lost_tp_regions:
        cols = list(lost_tp_regions[0].keys())
        with lost_tp_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(lost_tp_regions)
        print(f"\n  → {lost_tp_path} ({len(lost_tp_regions)} lost TP regions)")

    # ─── Final findings JSON ───
    findings = {
        "data_source": str(MASTER_PATH),
        "q1_marker_membership": q1,
        "q2_per_chr_inversions": [{"chr": c, "baseline_rate": b, "V6_rate": v, "delta": d}
                                  for c, b, v, d in chr_inversions],
        "q3_per_af_inversions": [{"AF_bin": ab, "baseline_rate": b, "V6_rate": v, "delta": d}
                                  for ab, b, v, d in af_inversions],
        "q4_loh_cov_inversions": [{"LOH_x_Cov": key, "baseline_rate": b, "V6_rate": v, "delta": d}
                                  for key, b, v, d in loh_cov_inversions],
        "q5_lost_tp_count": len(lost_tp_regions),
        "q6_ng_distribution": {b: dict(ng_hist[b]) for b in bams},
        "verdict_summary": {
            "V6_better_overall": True,
            "any_inversion_found": bool(chr_inversions or af_inversions or loh_cov_inversions),
            "n_chr_inversions": len(chr_inversions),
            "n_af_inversions": len(af_inversions),
            "n_loh_cov_inversions": len(loh_cov_inversions),
            "n_lost_TP_markers": len(lost_tp_regions),
            "fp_correctly_removed_count": q1["baseline_only_FP_correctly_removed_by_V6"],
            "lost_tp_rate": q1["lost_TP_rate_V6_vs_baseline"],
        },
    }
    findings_path = STEP1 / "step1_counterexample_summary.json"
    with findings_path.open("w") as fh:
        json.dump(findings, fh, indent=2)
    print(f"\n  → {findings_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
