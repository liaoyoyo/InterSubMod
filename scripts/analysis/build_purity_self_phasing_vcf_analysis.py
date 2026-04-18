#!/usr/bin/env python3
"""
Purity-dependent self-phasing analysis from existing ClairS-TO phased VCFs.

Analyzes how self-phasing circular dependency changes with tumor purity using
the phased VCF outputs that already exist for HCC1395 subsamples.

Key VCF fields used:
  - GT: phased genotype (0|1 / 1|0)
  - PS: phase set ID (scaffold membership)
  - AF: allele frequency
  - FILTER: PASS / NonSomatic / etc.
  - INFO H flag: "Variant found only in one haplotype" (= LOH-like proxy)
  - INFO Verdict_Somatic / Verdict_Germline flags

Output: TSV tables + summary to workspace directory.
"""
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# --- Configuration ---
SEQC2_TRUTH = Path("/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz")
SUBSAMPLE_BASE = Path("/big8_disk/data/HCC1395/ONT/subsample")
OUTPUT_DIR = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260402_phasing_causal_chain")

PURITY_LEVELS = [
    ("t10_n40", 0.20),
    ("t20_n30", 0.40),
    ("t30_n20", 0.60),
    ("t40_n10", 0.80),
    ("t50_n00", 1.00),
]

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]

# --- Helper functions ---

def load_seqc2_truth_set():
    """Load SEQC2 truth set SNV positions into a set of (chrom, pos)."""
    truth = set()
    with gzip.open(SEQC2_TRUTH, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos = fields[0], int(fields[1])
            truth.add((chrom, pos))
    print(f"  Loaded {len(truth)} SEQC2 truth set positions")
    return truth


def parse_phased_vcf(vcf_path):
    """Parse a single phased VCF, return list of variant dicts."""
    variants = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            filt = fields[6]
            info = fields[7]
            fmt_keys = fields[8].split(":")
            fmt_vals = fields[9].split(":")

            fmt = dict(zip(fmt_keys, fmt_vals))

            # Parse INFO flags
            info_fields = info.split(";")
            info_set = set()
            for f_item in info_fields:
                if "=" not in f_item:
                    info_set.add(f_item)

            gt = fmt.get("GT", ".")
            is_phased = "|" in gt
            ps = fmt.get("PS", ".")
            has_ps = ps != "." and ps != "0"

            try:
                af = float(fmt.get("AF", "0"))
            except ValueError:
                af = 0.0
            try:
                dp = int(fmt.get("DP", "0"))
            except ValueError:
                dp = 0

            variants.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "filter": filt,
                "gt": gt,
                "is_phased": is_phased,
                "has_ps": has_ps,
                "ps": ps,
                "af": af,
                "dp": dp,
                "H_flag": "H" in info_set,
                "verdict_somatic": "Verdict_Somatic" in info_set,
                "verdict_germline": "Verdict_Germline" in info_set,
                "verdict_subclonal": "Verdict_SubclonalSomatic" in info_set,
            })
    return variants


def analyze_purity_level(subsample_dir, purity_name, purity_value, truth_set):
    """Analyze all chromosomes for one purity level."""
    phased_dir = subsample_dir / "ClairS_TO_v0_3_0" / "tmp" / "phasing_output" / "phased_vcf_output"

    all_variants = []
    for chrom in CHROMOSOMES:
        vcf_path = phased_dir / f"tumor_phased_{chrom}.vcf.gz"
        if not vcf_path.exists():
            continue
        variants = parse_phased_vcf(vcf_path)
        all_variants.extend(variants)

    if not all_variants:
        print(f"  WARNING: No variants found for {purity_name}")
        return None

    df = pd.DataFrame(all_variants)

    # Classify: is this a SEQC2 truth TP?
    df["is_seqc2_tp"] = df.apply(lambda r: (r["chrom"], r["pos"]) in truth_set, axis=1)

    # Classify by FILTER
    df["is_pass"] = df["filter"] == "PASS"
    df["is_nonsomatic"] = df["filter"].str.contains("NonSomatic", na=False)

    # PASS + SEQC2 TP = true positive somatic in TO pipeline
    df["is_tp_pass"] = df["is_pass"] & df["is_seqc2_tp"]
    # PASS + not SEQC2 TP = FP somatic call
    df["is_fp_pass"] = df["is_pass"] & ~df["is_seqc2_tp"]

    # --- Compute metrics ---
    n_total = len(df)
    n_pass = df["is_pass"].sum()
    n_nonsomatic = df["is_nonsomatic"].sum()
    n_tp_pass = df["is_tp_pass"].sum()
    n_fp_pass = df["is_fp_pass"].sum()

    # Self-phasing: TP somatic variants with PS assignment in phased VCF
    tp_in_scaffold = df[df["is_seqc2_tp"] & df["has_ps"]]
    n_tp_in_scaffold = len(tp_in_scaffold)
    n_tp_total_in_vcf = df["is_seqc2_tp"].sum()
    tp_scaffold_rate = n_tp_in_scaffold / max(n_tp_total_in_vcf, 1)

    # H flag analysis: proxy for LOH-like (variant on single haplotype only)
    pass_df = df[df["is_pass"]].copy()
    tp_pass_h = pass_df[pass_df["is_tp_pass"]]["H_flag"].mean() if n_tp_pass > 0 else 0
    fp_pass_h = pass_df[pass_df["is_fp_pass"]]["H_flag"].mean() if n_fp_pass > 0 else 0

    # H flag for TP somatic in scaffold (self-phased TPs)
    tp_scaffold_h = tp_in_scaffold["H_flag"].mean() if len(tp_in_scaffold) > 0 else 0

    # NonSomatic H flag rate
    nonsom_df = df[df["is_nonsomatic"]]
    nonsom_h = nonsom_df["H_flag"].mean() if len(nonsom_df) > 0 else 0

    # AF distribution for TP PASS variants
    tp_af = pass_df[pass_df["is_tp_pass"]]["af"]
    tp_af_median = tp_af.median() if len(tp_af) > 0 else 0
    tp_af_mean = tp_af.mean() if len(tp_af) > 0 else 0

    # Phase block analysis
    ps_variants = df[df["has_ps"]]
    if len(ps_variants) > 0:
        block_sizes = ps_variants.groupby(["chrom", "ps"]).size()
        n_blocks = len(block_sizes)
        block_median = block_sizes.median()
        block_max = block_sizes.max()
        # Count blocks containing at least one SEQC2 TP
        blocks_with_tp = ps_variants[ps_variants["is_seqc2_tp"]].groupby(["chrom", "ps"]).size()
        n_blocks_with_tp = len(blocks_with_tp)
        block_contamination_rate = n_blocks_with_tp / max(n_blocks, 1)
    else:
        n_blocks = 0
        block_median = 0
        block_max = 0
        n_blocks_with_tp = 0
        block_contamination_rate = 0

    # Verdict analysis
    n_verdict_somatic = df["verdict_somatic"].sum()
    n_verdict_germline = df["verdict_germline"].sum()

    # AF-stratified H flag rate for TP PASS
    af_bins = [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 0.8), (0.8, 1.01)]
    af_strata = []
    for lo, hi in af_bins:
        subset = pass_df[(pass_df["is_tp_pass"]) & (pass_df["af"] >= lo) & (pass_df["af"] < hi)]
        af_strata.append({
            "purity": purity_value,
            "purity_name": purity_name,
            "af_bin": f"{lo:.1f}-{hi:.1f}",
            "n": len(subset),
            "h_flag_rate": subset["H_flag"].mean() if len(subset) > 0 else np.nan,
            "af_median": subset["af"].median() if len(subset) > 0 else np.nan,
        })

    result = {
        "purity": purity_value,
        "purity_name": purity_name,
        "n_total_variants": n_total,
        "n_pass": n_pass,
        "n_nonsomatic": n_nonsomatic,
        "n_seqc2_tp_in_vcf": n_tp_total_in_vcf,
        "n_tp_pass": n_tp_pass,
        "n_fp_pass": n_fp_pass,
        "tp_precision": n_tp_pass / max(n_pass, 1),
        "n_tp_in_scaffold": n_tp_in_scaffold,
        "tp_scaffold_rate": tp_scaffold_rate,
        "tp_pass_H_rate": tp_pass_h,
        "fp_pass_H_rate": fp_pass_h,
        "tp_scaffold_H_rate": tp_scaffold_h,
        "nonsomatic_H_rate": nonsom_h,
        "H_enrichment_TP_vs_FP": tp_pass_h / max(fp_pass_h, 1e-6),
        "tp_af_median": tp_af_median,
        "tp_af_mean": tp_af_mean,
        "expected_het_vaf": purity_value / 2,
        "n_phase_blocks": n_blocks,
        "block_median_size": block_median,
        "block_max_size": block_max,
        "n_blocks_with_tp": n_blocks_with_tp,
        "block_contamination_rate": block_contamination_rate,
        "n_verdict_somatic": n_verdict_somatic,
        "n_verdict_germline": n_verdict_germline,
    }

    return result, af_strata


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading SEQC2 truth set...")
    truth_set = load_seqc2_truth_set()

    results = []
    all_af_strata = []

    for purity_name, purity_value in PURITY_LEVELS:
        print(f"\nAnalyzing {purity_name} (purity={purity_value:.0%})...")
        subsample_dir = SUBSAMPLE_BASE / purity_name

        result, af_strata = analyze_purity_level(subsample_dir, purity_name, purity_value, truth_set)
        if result:
            results.append(result)
            all_af_strata.extend(af_strata)
            # Print key metrics
            print(f"  PASS variants: {result['n_pass']:,}")
            print(f"  TP PASS: {result['n_tp_pass']:,}, FP PASS: {result['n_fp_pass']:,}")
            print(f"  TP in scaffold: {result['n_tp_in_scaffold']:,} ({result['tp_scaffold_rate']:.1%})")
            print(f"  TP H-flag rate: {result['tp_pass_H_rate']:.3f}, FP H-flag rate: {result['fp_pass_H_rate']:.3f}")
            print(f"  TP AF median: {result['tp_af_median']:.3f} (expected: {result['expected_het_vaf']:.3f})")
            print(f"  Phase blocks: {result['n_phase_blocks']:,}, contaminated: {result['n_blocks_with_tp']:,} ({result['block_contamination_rate']:.1%})")

    # Save results
    df_results = pd.DataFrame(results)
    out_path = OUTPUT_DIR / "purity_self_phasing_vcf_summary.tsv"
    df_results.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved summary: {out_path}")

    df_af = pd.DataFrame(all_af_strata)
    out_af = OUTPUT_DIR / "purity_self_phasing_af_stratified.tsv"
    df_af.to_csv(out_af, sep="\t", index=False)
    print(f"Saved AF-stratified: {out_af}")

    # Print comparison table
    print("\n" + "=" * 100)
    print("PURITY-DEPENDENT SELF-PHASING SUMMARY")
    print("=" * 100)
    cols = ["purity", "n_tp_pass", "n_fp_pass", "tp_scaffold_rate",
            "tp_pass_H_rate", "fp_pass_H_rate", "H_enrichment_TP_vs_FP",
            "tp_af_median", "block_contamination_rate"]
    print(df_results[cols].to_string(index=False, float_format="%.4f"))

    # Key hypothesis test
    print("\n" + "-" * 80)
    print("KEY HYPOTHESIS: Self-phasing effect weakens with lower purity")
    print("-" * 80)
    for _, row in df_results.iterrows():
        purity = row["purity"]
        h_tp = row["tp_pass_H_rate"]
        h_fp = row["fp_pass_H_rate"]
        scaffold = row["tp_scaffold_rate"]
        contam = row["block_contamination_rate"]
        print(f"  Purity {purity:.0%}: TP H-flag={h_tp:.3f}, FP H-flag={h_fp:.3f}, "
              f"TP-in-scaffold={scaffold:.1%}, block-contam={contam:.1%}")

    if len(df_results) >= 2:
        h_100 = df_results[df_results["purity"] == 1.0]["tp_pass_H_rate"].values[0]
        h_20 = df_results[df_results["purity"] == 0.2]["tp_pass_H_rate"].values[0]
        print(f"\n  H-flag TP rate: 100% purity = {h_100:.3f}, 20% purity = {h_20:.3f}")
        print(f"  Ratio: {h_20 / max(h_100, 1e-6):.2f}x")
        if h_20 < h_100:
            print("  => CONFIRMED: Self-phasing H-flag effect weakens at lower purity")
        else:
            print("  => UNEXPECTED: H-flag rate does not decrease with lower purity")


if __name__ == "__main__":
    main()
