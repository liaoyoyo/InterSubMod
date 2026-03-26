#!/usr/bin/env python3
"""03_filter_analysis.py - Filter condition analysis and F1 score calculation.

Reads InterSubMod significance_summary.csv for TP and FP sites,
optionally enriches with VCF features (QUAL, VAF),
applies methylation-based filtering conditions,
and computes baseline/filtered Precision, Recall, F1.

Usage:
    python3 03_filter_analysis.py \
        --tp-summary intersubmod_tp/significance_summary.csv \
        --fp-summary intersubmod_fp/significance_summary.csv \
        --tp-count 30476 --fp-count 4822 --truth-total 39447 \
        --output-dir /path/to/output \
        [--somatic-vcf /path/to/somatic.vcf.gz] \
        [--sample HCC1395] [--mode s-pure] [--purity pure]
"""

import argparse
import csv
import gzip
import json
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Methylation filter analysis and F1 calculation")
    parser.add_argument("--tp-summary", required=True, help="Path to TP significance_summary.csv")
    parser.add_argument("--fp-summary", required=True, help="Path to FP significance_summary.csv")
    parser.add_argument("--tp-count", type=int, required=True, help="Total TP variant count")
    parser.add_argument("--fp-count", type=int, required=True, help="Total FP variant count")
    parser.add_argument("--truth-total", type=int, required=True, help="Total truth set variant count")
    parser.add_argument("--output-dir", required=True, help="Output directory for metrics.json")
    parser.add_argument("--somatic-vcf", default=None, help="Original somatic VCF for QUAL/VAF extraction (legacy)")
    parser.add_argument("--tp-vcf-file", default=None, help="TP VCF file for QUAL/VAF extraction (preferred)")
    parser.add_argument("--fp-vcf-file", default=None, help="FP VCF file for QUAL/VAF extraction (preferred)")
    parser.add_argument("--sample", default="unknown", help="Sample name")
    parser.add_argument("--mode", default="s-pure", help="Verification mode")
    parser.add_argument("--purity", default="pure", help="Purity level")
    return parser.parse_args()


def parse_purity_value(raw_value):
    try:
        return float(raw_value)
    except (TypeError, ValueError):
        return None


def purity_bin_from_value(purity_value):
    if purity_value is None:
        return "unknown"
    if purity_value < 40.0:
        return "lt40"
    if purity_value < 60.0:
        return "40to60"
    return "ge60"


def determine_rule_config(raw_purity):
    purity_value = parse_purity_value(raw_purity)
    purity_bin = purity_bin_from_value(purity_value)

    if purity_bin == "lt40":
        return {
            "rule_id": "purity_lt40_annotation_only_v1",
            "purity_bin": purity_bin,
            "apply_filter": False,
            "allele_delta_min": None,
            # allele_delta_only_min disabled: cross-sample AUROC inconsistent
            # (HCC1395-5kHz=0.763, HCC1395-DORADO=0.412, H2009=0.379)
            "allele_delta_only_min": None,
            "vaf_max": None,
            "require_cv_support": False,
            "cv_support_max": 0.03,
        }

    if purity_bin == "40to60":
        return {
            "rule_id": "purity_40to60_conservative_ad_vaf_cv_v1",
            "purity_bin": purity_bin,
            "apply_filter": True,
            "allele_delta_min": 0.20,
            # allele_delta_only_min disabled: cross-sample AUROC inconsistent
            "allele_delta_only_min": None,
            "vaf_max": 0.10,
            "require_cv_support": True,
            "cv_support_max": 0.03,
        }

    # ge60 / unknown: VCF-based AD+VAF rule retained; standalone AlleleDelta
    # filter disabled pending cross-sample validation.
    # Evidence: AUROC=0.763 (HCC1395-5kHz) but 0.412 (HCC1395-DORADO),
    # 0.545 (COLO829), 0.379 (H2009) — not cross-sample consistent.
    return {
        "rule_id": "purity_ge60_ad_vaf_core_v1" if purity_bin == "ge60" else "purity_unknown_ad_vaf_core_v1",
        "purity_bin": purity_bin,
        "apply_filter": True,
        "allele_delta_min": 0.15,
        "allele_delta_only_min": None,  # disabled — see evidence above
        "experimental_hcc1395_5khz_only": True,  # annotation flag
        "vaf_max": 0.15,
        "require_cv_support": False,
        "cv_support_max": 0.03,
    }


def read_summary_csv(filepath):
    """Read significance_summary.csv and return list of dicts."""
    rows = []
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def parse_vcf_features(vcf_path):
    """Extract QUAL and VAF from somatic VCF. Returns dict keyed by (chr, pos)."""
    features = {}
    if vcf_path is None:
        return features

    open_func = gzip.open if vcf_path.endswith(".gz") else open
    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            qual = float(fields[5]) if fields[5] != "." else 0.0

            # Try to extract VAF from FORMAT/SAMPLE fields
            vaf = None
            fmt_fields = fields[8].split(":")
            sample_fields = fields[9].split(":")
            fmt_dict = dict(zip(fmt_fields, sample_fields))

            if "VAF" in fmt_dict:
                try:
                    vaf = float(fmt_dict["VAF"])
                except (ValueError, TypeError):
                    pass
            elif "AF" in fmt_dict:
                try:
                    vaf = float(fmt_dict["AF"])
                except (ValueError, TypeError):
                    pass
            elif "AD" in fmt_dict:
                try:
                    ad_values = fmt_dict["AD"].split(",")
                    if len(ad_values) >= 2:
                        ref_count = int(ad_values[0])
                        alt_count = int(ad_values[1])
                        total = ref_count + alt_count
                        vaf = alt_count / total if total > 0 else None
                except (ValueError, TypeError):
                    pass

            features[(chrom, pos)] = {"qual": qual, "vaf": vaf}

    return features


def should_filter_variant(row, vcf_features=None, rule_config=None):
    """Apply methylation-based filter conditions.

    Purity-aware default:
        <40% purity   -> annotation only, do not hard filter
        40-60% purity -> conservative AD+VAF core with CV support
        >=60% purity  -> AD+VAF core, CV as auxiliary signal

    Additional methylation-only rule (no VCF required):
        AlleleDelta > 0.25 indicates strong haplotype-specific methylation
        asymmetry consistent with germline ASM — applied after VCF core rule.

    Returns True if the variant should be REMOVED (filtered out).
    """
    if rule_config is None:
        rule_config = determine_rule_config(None)

    try:
        allele_delta = abs(float(row.get("AlleleDelta", 0) or 0))
        cramers_v = float(row.get("CramersV", 0) or 0)
    except (ValueError, TypeError):
        return False

    if not rule_config["apply_filter"]:
        return False

    chrom = row.get("Chr", "")
    pos = int(row.get("Pos", 0))

    # Get VCF features if available
    qual = None
    vaf = None
    if vcf_features and (chrom, pos) in vcf_features:
        feat = vcf_features[(chrom, pos)]
        qual = feat.get("qual")
        vaf = feat.get("vaf")

    # Primary rule: AlleleDelta + VAF (requires VCF)
    if vaf is not None:
        core_trigger = (
            allele_delta > rule_config["allele_delta_min"]
            and vaf < rule_config["vaf_max"]
        )
        if core_trigger:
            if rule_config["require_cv_support"] and cramers_v >= rule_config["cv_support_max"]:
                return False
            return True

    # Methylation-only fallback: disabled — cross-sample AUROC inconsistent.
    # AlleleDelta AUROC(FP-discriminating): HCC1395-5kHz=0.763 but
    # HCC1395-DORADO=0.412, COLO829=0.545, H2009=0.379.
    # Not safe to apply as a cross-sample hard filter.
    # Feature retained in output (AlleleDelta column) for future analysis.
    allele_delta_only_threshold = rule_config.get("allele_delta_only_min", None)
    if allele_delta_only_threshold is not None and allele_delta > allele_delta_only_threshold:
        return True

    return False


def compute_metrics(tp, fp, truth_total):
    """Compute precision, recall, and F1 score."""
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / truth_total if truth_total > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"tp": tp, "fp": fp, "precision": round(precision, 4), "recall": round(recall, 4), "f1": round(f1, 4)}


def main():
    args = parse_args()
    rule_config = determine_rule_config(args.purity)

    # Read significance summaries
    print(f"[Filter Analysis] Reading TP summary: {args.tp_summary}")
    tp_rows = read_summary_csv(args.tp_summary)
    print(f"[Filter Analysis] Reading FP summary: {args.fp_summary}")
    fp_rows = read_summary_csv(args.fp_summary)

    print(f"[Filter Analysis] TP regions: {len(tp_rows)}, FP regions: {len(fp_rows)}")
    print(f"[Filter Analysis] Rule: {rule_config['rule_id']} | purity_bin={rule_config['purity_bin']}")

    # Load VCF features
    # Prefer separate TP/FP VCF files (which have LongPhase-S probability-scale QUAL)
    # Fall back to single somatic VCF if provided
    tp_vcf_features = None
    fp_vcf_features = None

    if args.tp_vcf_file or args.fp_vcf_file:
        if args.tp_vcf_file:
            print(f"[Filter Analysis] Loading TP VCF features: {args.tp_vcf_file}")
            tp_vcf_features = parse_vcf_features(args.tp_vcf_file)
            print(f"[Filter Analysis] Loaded TP features for {len(tp_vcf_features)} variants")
        if args.fp_vcf_file:
            print(f"[Filter Analysis] Loading FP VCF features: {args.fp_vcf_file}")
            fp_vcf_features = parse_vcf_features(args.fp_vcf_file)
            print(f"[Filter Analysis] Loaded FP features for {len(fp_vcf_features)} variants")
    elif args.somatic_vcf:
        print(f"[Filter Analysis] Loading VCF features: {args.somatic_vcf}")
        tp_vcf_features = parse_vcf_features(args.somatic_vcf)
        fp_vcf_features = tp_vcf_features
        print(f"[Filter Analysis] Loaded features for {len(tp_vcf_features)} variants")

    # Baseline metrics
    baseline = compute_metrics(args.tp_count, args.fp_count, args.truth_total)
    print(f"[Filter Analysis] Baseline: P={baseline['precision']}, R={baseline['recall']}, F1={baseline['f1']}")

    # Apply filter to TP variants (count how many would be incorrectly removed)
    tp_filtered_count = 0
    for row in tp_rows:
        if should_filter_variant(row, tp_vcf_features, rule_config):
            tp_filtered_count += 1

    # Apply filter to FP variants (count how many would be correctly removed)
    fp_filtered_count = 0
    for row in fp_rows:
        if should_filter_variant(row, fp_vcf_features, rule_config):
            fp_filtered_count += 1

    # Filtered metrics
    filtered_tp = args.tp_count - tp_filtered_count
    filtered_fp = args.fp_count - fp_filtered_count
    filtered = compute_metrics(filtered_tp, filtered_fp, args.truth_total)

    print(f"[Filter Analysis] Filtered: P={filtered['precision']}, R={filtered['recall']}, F1={filtered['f1']}")
    print(f"[Filter Analysis] TP removed: {tp_filtered_count}, FP removed: {fp_filtered_count}")

    # Build output
    result = {
        "sample": args.sample,
        "mode": args.mode,
        "purity": args.purity,
        "purity_bin": rule_config["purity_bin"],
        "truth_total": args.truth_total,
        "tp_regions_analyzed": len(tp_rows),
        "fp_regions_analyzed": len(fp_rows),
        "rule": rule_config,
        "baseline": baseline,
        "filtered": filtered,
        "improvement": {
            "f1_delta": round(filtered["f1"] - baseline["f1"], 4),
            "fp_removed": fp_filtered_count,
            "tp_removed": tp_filtered_count,
            "fp_removed_rate": round(fp_filtered_count / args.fp_count, 6) if args.fp_count else 0.0,
            "tp_removed_rate": round(tp_filtered_count / args.tp_count, 6) if args.tp_count else 0.0,
            "trigger_count": fp_filtered_count + tp_filtered_count,
            "trigger_precision_fp": (
                round(fp_filtered_count / (fp_filtered_count + tp_filtered_count), 6)
                if (fp_filtered_count + tp_filtered_count) > 0
                else 0.0
            ),
        },
    }

    # Write output
    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, "metrics.json")
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"[Filter Analysis] Metrics written to: {output_path}")

    # Print summary
    print()
    print("=" * 60)
    print(f"  Sample: {args.sample} | Mode: {args.mode} | Purity: {args.purity}")
    print(f"  Truth total: {args.truth_total}")
    print(f"  Baseline:  TP={baseline['tp']}, FP={baseline['fp']}, "
          f"P={baseline['precision']}, R={baseline['recall']}, F1={baseline['f1']}")
    print(f"  Filtered:  TP={filtered['tp']}, FP={filtered['fp']}, "
          f"P={filtered['precision']}, R={filtered['recall']}, F1={filtered['f1']}")
    print(f"  F1 delta:  {result['improvement']['f1_delta']:+.4f}")
    print(f"  FP removed: {fp_filtered_count}, TP removed: {tp_filtered_count}")
    print("=" * 60)


if __name__ == "__main__":
    main()
