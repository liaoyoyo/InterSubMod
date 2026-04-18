#!/usr/bin/env python3
"""
Purity-dependent LOH proxy analysis: AF-based HP imbalance from phased VCFs.

Uses AF field in ClairS-TO phased VCF as LOH-like proxy:
  AF > 0.9 → extreme allele imbalance → LOH-like (analogous to HP_Ratio > 0.9)

Key questions:
  1. Does TP LOH-like rate decrease with purity? (self-phasing weakens)
  2. Does the TP>FP LOH enrichment pattern reverse at lower purity?
  3. At what purity does the effect disappear?
"""
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd

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
LOH_THRESH = 0.9  # AF > this = LOH-like


def load_truth():
    truth = set()
    with gzip.open(SEQC2_TRUTH, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 3)
            truth.add((fields[0], int(fields[1])))
    return truth


def analyze_purity(subsample_dir, purity_name, purity_value, truth_set):
    phased_dir = subsample_dir / "ClairS_TO_v0_3_0" / "tmp" / "phasing_output" / "phased_vcf_output"

    records = []
    for chrom in CHROMOSOMES:
        vcf_path = phased_dir / f"tumor_phased_{chrom}.vcf.gz"
        if not vcf_path.exists():
            continue
        with gzip.open(vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 10:
                    continue

                chrom_v = fields[0]
                pos = int(fields[1])
                filt = fields[6]

                # Only PASS variants (somatic calls)
                if filt != "PASS":
                    continue

                fmt_keys = fields[8].split(":")
                fmt_vals = fields[9].strip().split(":")
                fmt = dict(zip(fmt_keys, fmt_vals))

                try:
                    af = float(fmt.get("AF", "0"))
                except ValueError:
                    continue
                try:
                    dp = int(fmt.get("DP", "0"))
                except ValueError:
                    dp = 0

                gt = fmt.get("GT", ".")
                ps = fmt.get("PS", ".")
                has_ps = "|" in gt and ps not in (".", "0")
                is_tp = (chrom_v, pos) in truth_set

                records.append({
                    "af": af,
                    "dp": dp,
                    "has_ps": has_ps,
                    "is_tp": is_tp,
                    "loh_like": af > LOH_THRESH,
                    "loh_like_08": af > 0.8,
                })

    df = pd.DataFrame(records)
    n = len(df)
    tp = df[df["is_tp"]]
    fp = df[~df["is_tp"]]

    # LOH-like rates
    tp_loh_rate = tp["loh_like"].mean() if len(tp) > 0 else 0
    fp_loh_rate = fp["loh_like"].mean() if len(fp) > 0 else 0
    tp_loh08_rate = tp["loh_like_08"].mean() if len(tp) > 0 else 0
    fp_loh08_rate = fp["loh_like_08"].mean() if len(fp) > 0 else 0

    # LOH enrichment (TP/FP style)
    enrichment_09 = fp_loh_rate / max(tp_loh_rate, 1e-6) if tp_loh_rate > 0 else float("nan")
    enrichment_08 = fp_loh08_rate / max(tp_loh08_rate, 1e-6) if tp_loh08_rate > 0 else float("nan")

    # Scaffold-specific
    tp_scaffold = tp[tp["has_ps"]]
    tp_no_scaffold = tp[~tp["has_ps"]]
    tp_scaffold_loh = tp_scaffold["loh_like"].mean() if len(tp_scaffold) > 0 else 0
    tp_noscaffold_loh = tp_no_scaffold["loh_like"].mean() if len(tp_no_scaffold) > 0 else 0

    # AF distribution
    tp_af_median = tp["af"].median() if len(tp) > 0 else 0
    fp_af_median = fp["af"].median() if len(fp) > 0 else 0

    # DP-stratified LOH (for depth >=30, analogous to Tier A+)
    tp_highDP = tp[tp["dp"] >= 30]
    fp_highDP = fp[fp["dp"] >= 30]
    tp_highDP_loh = tp_highDP["loh_like"].mean() if len(tp_highDP) > 0 else 0
    fp_highDP_loh = fp_highDP["loh_like"].mean() if len(fp_highDP) > 0 else 0

    return {
        "purity": purity_value,
        "purity_name": purity_name,
        "n_pass": n,
        "n_tp": len(tp),
        "n_fp": len(fp),
        "tp_af_median": tp_af_median,
        "fp_af_median": fp_af_median,
        "tp_loh09_rate": tp_loh_rate,
        "fp_loh09_rate": fp_loh_rate,
        "loh09_enrichment_FP_over_TP": enrichment_09,
        "tp_loh08_rate": tp_loh08_rate,
        "fp_loh08_rate": fp_loh08_rate,
        "loh08_enrichment_FP_over_TP": enrichment_08,
        "tp_scaffold_loh09": tp_scaffold_loh,
        "tp_noscaffold_loh09": tp_noscaffold_loh,
        "n_tp_scaffold": len(tp_scaffold),
        "n_tp_noscaffold": len(tp_no_scaffold),
        "tp_highDP_loh09": tp_highDP_loh,
        "fp_highDP_loh09": fp_highDP_loh,
        "n_tp_highDP": len(tp_highDP),
        "n_fp_highDP": len(fp_highDP),
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    truth = load_truth()
    print(f"SEQC2 truth: {len(truth)} positions\n")

    results = []
    for pname, pval in PURITY_LEVELS:
        print(f"Analyzing {pname} (purity={pval:.0%})...")
        r = analyze_purity(SUBSAMPLE_BASE / pname, pname, pval, truth)
        results.append(r)
        print(f"  TP={r['n_tp']:,}, FP={r['n_fp']:,}")
        print(f"  TP AF median={r['tp_af_median']:.3f}, FP AF median={r['fp_af_median']:.3f}")
        print(f"  TP LOH@0.9={r['tp_loh09_rate']:.4f}, FP LOH@0.9={r['fp_loh09_rate']:.4f}, "
              f"enrichment(FP/TP)={r['loh09_enrichment_FP_over_TP']:.3f}")
        print(f"  TP(scaffold) LOH@0.9={r['tp_scaffold_loh09']:.4f}, "
              f"TP(no-scaffold) LOH@0.9={r['tp_noscaffold_loh09']:.4f}")
        print(f"  TP(DP>=30) LOH@0.9={r['tp_highDP_loh09']:.4f}, "
              f"FP(DP>=30) LOH@0.9={r['fp_highDP_loh09']:.4f}")

    df = pd.DataFrame(results)
    out = OUTPUT_DIR / "purity_loh_proxy_by_af_threshold.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nSaved: {out}")

    # Print key comparison
    print("\n" + "=" * 90)
    print("LOH ENRICHMENT DIRECTION vs PURITY (core hypothesis)")
    print("=" * 90)
    print(f"{'Purity':>8} {'TP LOH%':>10} {'FP LOH%':>10} {'FP/TP':>8} {'Direction':>20} {'TP AF':>8}")
    print("-" * 90)
    for _, r in df.iterrows():
        tp_l = r["tp_loh09_rate"]
        fp_l = r["fp_loh09_rate"]
        enr = r["loh09_enrichment_FP_over_TP"]
        direction = "TP-enriched" if enr < 1 else ("FP-enriched" if enr > 1 else "Neutral")
        print(f"{r['purity']:>7.0%} {tp_l:>10.4f} {fp_l:>10.4f} {enr:>8.3f} {direction:>20} {r['tp_af_median']:>8.3f}")

    print("\n" + "-" * 90)
    print("INTERPRETATION:")
    enrichments = df["loh09_enrichment_FP_over_TP"].values
    purities = df["purity"].values
    tp_enriched = [p for p, e in zip(purities, enrichments) if e < 1]
    fp_enriched = [p for p, e in zip(purities, enrichments) if e > 1]
    if tp_enriched:
        print(f"  TP-enriched (LOH marks TP): purity = {[f'{p:.0%}' for p in tp_enriched]}")
    if fp_enriched:
        print(f"  FP-enriched (LOH marks FP): purity = {[f'{p:.0%}' for p in fp_enriched]}")

    # Check the self-phasing scaffold effect
    print("\n" + "=" * 90)
    print("SELF-PHASING SCAFFOLD EFFECT: TP in scaffold vs not in scaffold")
    print("=" * 90)
    print(f"{'Purity':>8} {'TP scaffold LOH%':>18} {'TP no-scaffold LOH%':>22} {'Ratio':>8}")
    print("-" * 90)
    for _, r in df.iterrows():
        s = r["tp_scaffold_loh09"]
        ns = r["tp_noscaffold_loh09"]
        ratio = s / max(ns, 1e-6) if ns > 0 else float("nan")
        print(f"{r['purity']:>7.0%} {s:>18.4f} {ns:>22.4f} {ratio:>8.2f}x")


if __name__ == "__main__":
    main()
