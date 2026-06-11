#!/usr/bin/env python3
"""baseline vs V6 head-to-head analysis (HCC1395, full genome).

Reads step1_master_four_way.tsv (after build_four_way_master.py).
Produces a comparison TSV + summary statistics covering:

1. Global metrics — total reads, hp distribution, hp_11:hp_21 ratio
2. Marker coverage (NG_off ≥ 3) — baseline vs V6, per TP/FP
3. Marker rate at NG_off ≥ 3 (TP / (TP+FP))
4. Flag=on collapse — NG_on=2 rate, robustness
5. hp=33 conservative-ambiguous reads count (V6 expected ≫ baseline)
6. Per-chromosome HP1:HP2 ratio drift
7. LOH-stratified TP/FP rate change

Outputs:
- step1_baseline_vs_v6_summary.tsv   — high-level summary
- step1_baseline_vs_v6_per_chr.tsv   — per-chromosome breakdown
- step1_baseline_vs_v6_per_zone.tsv  — LOH × CN zone breakdown
- step1_baseline_vs_v6_findings.json — machine-readable summary
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

STEP1 = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way")
MASTER_PATH = STEP1 / "step1_master_four_way.tsv"

HP_FAMILY = ["0", "1", "2", "1-1", "2-1", "3", "11", "21", "33", "other"]


def safe_int(x: str) -> int:
    try:
        return int(x)
    except (ValueError, TypeError):
        return 0


def safe_float(x: str) -> float:
    try:
        return float(x)
    except (ValueError, TypeError):
        return float("nan")


def main() -> int:
    if not MASTER_PATH.exists():
        print(f"ERROR: {MASTER_PATH} not found. Run build_four_way_master.py first.", file=sys.stderr)
        return 1

    print(f"[analysis] reading {MASTER_PATH}", flush=True)

    # Per-BAM × off counters
    bams = ["baseline", "V3F", "V5", "V6"]
    global_hp = {b: {k: 0 for k in HP_FAMILY} for b in bams}
    global_n_reads = {b: 0 for b in bams}

    # marker filter at NG_off >= 3
    marker_tp = {b: 0 for b in bams}
    marker_fp = {b: 0 for b in bams}

    # NG_on=2 canonical formula (matches phaseC_region_ng_fast.py):
    # rate_on = TP_at_NG_on_2 / (TP_at_NG_on_2 + FP_at_NG_on_2)
    # No NG_off prefilter — TP purity among NG_on=2 regions across full genome.
    # This aligns with V6 doc §5.3 published numbers (V3F=0.8579, V5=0.8285, V6=0.8285).
    ng_on2_tp = {b: 0 for b in bams}
    ng_on2_fp = {b: 0 for b in bams}

    # per-chr hp1/hp2 series for off mode
    per_chr_hp = defaultdict(lambda: {b: {"hp1_series": 0, "hp2_series": 0} for b in bams})

    # per LOH bucket TP/FP rate
    loh_buckets = defaultdict(lambda: {b: {"tp_marker": 0, "fp_marker": 0} for b in bams})

    # TP/FP counts for whole-genome
    tp_total = 0
    fp_total = 0

    with MASTER_PATH.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            label = row["label"]
            chrom = row["chr"]
            loh = row.get("LOH_Bed_Annotation", "") or "Unknown"
            if label == "TP":
                tp_total += 1
            else:
                fp_total += 1

            for bam in bams:
                # accumulate global hp
                n_reads = safe_int(row.get(f"{bam}_off_n_reads", "0"))
                global_n_reads[bam] += n_reads
                for k in HP_FAMILY:
                    global_hp[bam][k] += safe_int(row.get(f"{bam}_off_{k}", "0"))

                # marker filter NG_off ≥ 3
                ng_off = safe_int(row.get(f"{bam}_off_NG", "0"))
                if ng_off >= 3:
                    if label == "TP":
                        marker_tp[bam] += 1
                        loh_buckets[loh][bam]["tp_marker"] += 1
                    else:
                        marker_fp[bam] += 1
                        loh_buckets[loh][bam]["fp_marker"] += 1

                # NG_on=2 canonical (phaseC formula): TP purity at NG_on=2, NO NG_off prefilter
                ng_on = safe_int(row.get(f"{bam}_on_NG", "0"))
                if ng_on == 2:
                    if label == "TP":
                        ng_on2_tp[bam] += 1
                    else:
                        ng_on2_fp[bam] += 1

                # per-chr HP series (off mode)
                hp1 = safe_int(row.get(f"{bam}_off_1", "0")) + safe_int(row.get(f"{bam}_off_1-1", "0")) + safe_int(row.get(f"{bam}_off_11", "0"))
                hp2 = safe_int(row.get(f"{bam}_off_2", "0")) + safe_int(row.get(f"{bam}_off_2-1", "0")) + safe_int(row.get(f"{bam}_off_21", "0"))
                per_chr_hp[chrom][bam]["hp1_series"] += hp1
                per_chr_hp[chrom][bam]["hp2_series"] += hp2

    # ---- write summary TSV ----
    summary_rows = []
    for bam in bams:
        hp1_series = global_hp[bam]["1"] + global_hp[bam]["1-1"] + global_hp[bam]["11"]
        hp2_series = global_hp[bam]["2"] + global_hp[bam]["2-1"] + global_hp[bam]["21"]
        h11 = global_hp[bam]["1-1"] + global_hp[bam]["11"]
        h21 = global_hp[bam]["2-1"] + global_hp[bam]["21"]
        # hp_ambig_reads = hp=3 (single-int somatic ambig, TO mode) + hp=33 (string somatic ambig, paired mode).
        # In current HCC1395 TO data, hp=33 is empty for all 4 BAMs; all ambig reads are in hp=3.
        # Kept aggregated for forward-compatibility with paired-mode runs that emit hp=33 strings.
        hp_ambig = global_hp[bam]["33"] + global_hp[bam]["3"]
        hp3_only = global_hp[bam]["3"]
        hp33_only = global_hp[bam]["33"]
        marker_total = marker_tp[bam] + marker_fp[bam]
        rate_off = marker_tp[bam] / marker_total if marker_total > 0 else 0.0
        # NG_on=2 purity (phaseC canonical formula, matches V6 doc §5.3 published numbers)
        ng_on2_total = ng_on2_tp[bam] + ng_on2_fp[bam]
        ng_on2_purity = ng_on2_tp[bam] / ng_on2_total if ng_on2_total > 0 else 0.0
        ratio_hp1_hp2 = hp1_series / hp2_series if hp2_series > 0 else float("inf")
        ratio_h11_h21 = h11 / h21 if h21 > 0 else float("inf")
        summary_rows.append({
            "BAM": bam,
            "total_reads": global_n_reads[bam],
            "hp1_series": hp1_series,
            "hp2_series": hp2_series,
            "hp1_hp2_ratio": f"{ratio_hp1_hp2:.3f}",
            "h11_count": h11,
            "h21_count": h21,
            "h11_h21_ratio": f"{ratio_h11_h21:.3f}",
            "hp_ambig_reads": hp_ambig,
            "hp3_only": hp3_only,
            "hp33_only": hp33_only,
            "marker_NGge3_TP": marker_tp[bam],
            "marker_NGge3_FP": marker_fp[bam],
            "marker_total": marker_total,
            "marker_rate_off": f"{rate_off:.4f}",
            "NG_on2_TP": ng_on2_tp[bam],
            "NG_on2_FP": ng_on2_fp[bam],
            "NG_on2_purity": f"{ng_on2_purity:.4f}",
        })

    summary_path = STEP1 / "step1_baseline_vs_v6_summary.tsv"
    with summary_path.open("w", newline="") as fh:
        cols = list(summary_rows[0].keys())
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"  → {summary_path}")

    # ---- per chromosome ----
    per_chr_path = STEP1 / "step1_baseline_vs_v6_per_chr.tsv"
    with per_chr_path.open("w", newline="") as fh:
        cols = ["chr"]
        for bam in bams:
            cols += [f"{bam}_hp1", f"{bam}_hp2", f"{bam}_ratio"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for chrom in sorted(per_chr_hp.keys(), key=lambda c: (c != "chrX", c != "chrY", c)):
            r = {"chr": chrom}
            for bam in bams:
                hp1 = per_chr_hp[chrom][bam]["hp1_series"]
                hp2 = per_chr_hp[chrom][bam]["hp2_series"]
                r[f"{bam}_hp1"] = hp1
                r[f"{bam}_hp2"] = hp2
                r[f"{bam}_ratio"] = f"{hp1/hp2:.3f}" if hp2 > 0 else "INF"
            w.writerow(r)
    print(f"  → {per_chr_path}")

    # ---- per LOH zone ----
    per_zone_path = STEP1 / "step1_baseline_vs_v6_per_zone.tsv"
    with per_zone_path.open("w", newline="") as fh:
        cols = ["LOH_zone"]
        for bam in bams:
            cols += [f"{bam}_tp_marker", f"{bam}_fp_marker", f"{bam}_marker_rate"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for zone in sorted(loh_buckets.keys()):
            r = {"LOH_zone": zone}
            for bam in bams:
                tp = loh_buckets[zone][bam]["tp_marker"]
                fp = loh_buckets[zone][bam]["fp_marker"]
                tot = tp + fp
                r[f"{bam}_tp_marker"] = tp
                r[f"{bam}_fp_marker"] = fp
                r[f"{bam}_marker_rate"] = f"{tp/tot:.4f}" if tot > 0 else "NA"
            w.writerow(r)
    print(f"  → {per_zone_path}")

    # ---- findings JSON ----
    base_ambig = global_hp["baseline"]["33"] + global_hp["baseline"]["3"]
    v6_ambig = global_hp["V6"]["33"] + global_hp["V6"]["3"]
    findings = {
        "data_source": str(MASTER_PATH),
        "tp_total": tp_total,
        "fp_total": fp_total,
        "summary": summary_rows,
        "baseline_vs_v6_deltas": {
            "marker_count_delta": marker_tp["V6"] + marker_fp["V6"] - marker_tp["baseline"] - marker_fp["baseline"],
            "marker_rate_baseline": marker_tp["baseline"] / (marker_tp["baseline"] + marker_fp["baseline"]) if (marker_tp["baseline"] + marker_fp["baseline"]) > 0 else 0,
            "marker_rate_v6": marker_tp["V6"] / (marker_tp["V6"] + marker_fp["V6"]) if (marker_tp["V6"] + marker_fp["V6"]) > 0 else 0,
            "hp_ambig_increase": v6_ambig - base_ambig,
            "hp_ambig_breakdown_baseline": {"hp3": global_hp["baseline"]["3"], "hp33": global_hp["baseline"]["33"]},
            "hp_ambig_breakdown_v6": {"hp3": global_hp["V6"]["3"], "hp33": global_hp["V6"]["33"]},
            "hp1_hp2_ratio_baseline_to_v6": [
                f"{(global_hp['baseline']['1'] + global_hp['baseline']['1-1'] + global_hp['baseline']['11']) / max(1, global_hp['baseline']['2'] + global_hp['baseline']['2-1'] + global_hp['baseline']['21']):.3f}",
                f"{(global_hp['V6']['1'] + global_hp['V6']['1-1'] + global_hp['V6']['11']) / max(1, global_hp['V6']['2'] + global_hp['V6']['2-1'] + global_hp['V6']['21']):.3f}",
            ],
        },
    }
    findings_path = STEP1 / "step1_baseline_vs_v6_findings.json"
    with findings_path.open("w") as fh:
        json.dump(findings, fh, indent=2)
    print(f"  → {findings_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
