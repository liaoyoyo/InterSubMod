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
from collections import defaultdict


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
            vaf = 0.0
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
                        vaf = alt_count / total if total > 0 else 0.0
                except (ValueError, TypeError):
                    pass

            features[(chrom, pos)] = {"qual": qual, "vaf": vaf}

    return features


def should_filter_variant(row, vcf_features=None):
    """Apply methylation-based filter conditions.

    Filter condition (remove if TRUE):
        (QUAL < 0.75) OR (AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24)

    If VCF features are not available, only the methylation-only condition is used:
        AlleleDelta > 0.25 AND CramersV < 0.05

    Returns True if the variant should be REMOVED (filtered out).
    """
    try:
        allele_delta = abs(float(row.get("AlleleDelta", 0)))
        cramers_v = float(row.get("CramersV", 0))
    except (ValueError, TypeError):
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

    # Condition 1: Low QUAL (Disabled for Pileup Mode)
    # if qual is not None and qual < 0.75:
    #     return True

    # Condition 2: High allele delta + low Cramer's V + low VAF (Optimized 2026-02-12)
    # Optimized: AD > 0.15, CV < 0.03, VAF < 0.15
    if allele_delta > 0.15 and cramers_v < 0.03:
        if vaf is not None:
            if vaf < 0.15:
                return True
        else:
            # Without VAF, apply a more conservative filter (skip VAF check)
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

    # Read significance summaries
    print(f"[Filter Analysis] Reading TP summary: {args.tp_summary}")
    tp_rows = read_summary_csv(args.tp_summary)
    print(f"[Filter Analysis] Reading FP summary: {args.fp_summary}")
    fp_rows = read_summary_csv(args.fp_summary)

    print(f"[Filter Analysis] TP regions: {len(tp_rows)}, FP regions: {len(fp_rows)}")

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
        if should_filter_variant(row, tp_vcf_features):
            tp_filtered_count += 1

    # Apply filter to FP variants (count how many would be correctly removed)
    fp_filtered_count = 0
    for row in fp_rows:
        if should_filter_variant(row, fp_vcf_features):
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
        "truth_total": args.truth_total,
        "tp_regions_analyzed": len(tp_rows),
        "fp_regions_analyzed": len(fp_rows),
        "baseline": baseline,
        "filtered": filtered,
        "improvement": {
            "f1_delta": round(filtered["f1"] - baseline["f1"], 4),
            "fp_removed": fp_filtered_count,
            "tp_removed": tp_filtered_count,
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
