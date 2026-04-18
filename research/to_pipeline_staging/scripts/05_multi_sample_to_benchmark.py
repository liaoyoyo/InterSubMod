#!/usr/bin/env python3
"""
Multi-Sample ClairS-TO Benchmarking
====================================
Benchmark all 7 samples' ClairS-TO VCFs against truth sets,
extract VCF features, and compare TP/FP characteristics.

For HCC1395: uses canonical pipeline data (already benchmarked).
For others: uses bcftools isec to benchmark TO VCFs.

Output:
  data/multi_sample_to_summary.json — per-sample F1 + feature summary
  data/multi_sample_to_features.csv — all samples' PASS variant features
"""

import os, sys, gzip, json, subprocess, tempfile
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
os.makedirs(OUTDIR, exist_ok=True)

# ═══════════════════════════════════════
# Sample Definitions
# ═══════════════════════════════════════

SAMPLES = {
    'HCC1395': {
        'to_vcf': '/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/step01_clairs_to/snv.vcf.gz',
        'truth_vcf': '/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz',
        'truth_bed': '/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed',
        'truth_total': 39447,
        'cnv_bed': '/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed',
        'source': 'canonical TO pipeline (ONT_5kHz BAM, liaoyoyo2001)',
    },
    'HCC1395_DORADO': {
        'to_vcf': '/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/snv.vcf.gz',
        'truth_vcf': '/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz',
        'truth_bed': '/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed',
        'truth_total': 39447,
        'cnv_bed': '/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed',
        'source': 'zhenyu112',
    },
    'COLO829': {
        'to_vcf': '/big8_disk/data/COLO829/ONT_PAO/ClairS_TO_v0_3_0/snv.vcf.gz',
        'truth_vcf': '/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz',
        'truth_bed': '',
        'truth_total': 41427,
        'cnv_bed': '',
        'source': 'zhenyu112',
    },
    'H1437': {
        'to_vcf': '/big8_disk/data/H1437/ONT/ClairS_TO_v0_3_0/snv.vcf.gz',
        'truth_vcf': '/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark_somatic-only.vcf.gz',
        'truth_bed': '/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark.bed',
        'truth_total': 90129,
        'cnv_bed': '',
        'source': 'zhenyu112',
    },
    'H2009': {
        'to_vcf': '/big8_disk/data/H2009/ONT/ClairS_TO_v0_3_0/snv.vcf.gz',
        'truth_vcf': '/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark_somatic-only.vcf.gz',
        'truth_bed': '/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark.bed',
        'truth_total': 168529,
        'cnv_bed': '',
        'source': 'zhenyu112',
    },
    'HCC1937': {
        'to_vcf': '/big8_disk/data/HCC1937/ONT/ClairS_TO_v0_3_0/snv.vcf.gz',
        'truth_vcf': '/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark_somatic-only.vcf.gz',
        'truth_bed': '/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark.bed',
        'truth_total': 60691,
        'cnv_bed': '',
        'source': 'zhenyu112',
    },
    'HCC1954': {
        'to_vcf': '/big8_disk/data/HCC1954/ONT/ClairS_TO_v0_3_0/snv.vcf.gz',
        'truth_vcf': '/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark_somatic-only.vcf.gz',
        'truth_bed': '/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark.bed',
        'truth_total': 26567,
        'cnv_bed': '',
        'source': 'zhenyu112',
    },
}


def benchmark_to_vcf(sample_name, config, workdir):
    """Benchmark a TO VCF against truth using bcftools isec."""
    to_vcf = config['to_vcf']
    truth_vcf = config['truth_vcf']
    truth_bed = config['truth_bed']
    truth_total = config['truth_total']

    print(f"\n{'='*60}")
    print(f"Benchmarking: {sample_name}")
    print(f"  TO VCF: {to_vcf}")
    print(f"  Truth: {truth_vcf}")
    print(f"  BED: {truth_bed or 'NONE'}")
    print(f"{'='*60}")

    # Step 1: Extract PASS from TO VCF
    pass_vcf = os.path.join(workdir, f"{sample_name}_pass.vcf.gz")
    cmd = f"bcftools view -f PASS -Oz -o {pass_vcf} {to_vcf} && bcftools index {pass_vcf}"
    subprocess.run(cmd, shell=True, check=True, capture_output=True)
    pass_count = int(subprocess.run(
        f"bcftools view {pass_vcf} | grep -v '^#' | wc -l",
        shell=True, capture_output=True, text=True).stdout.strip())
    print(f"  PASS variants: {pass_count:,}")

    # Step 2: Scope to BED if available
    if truth_bed and os.path.exists(truth_bed):
        scoped_pass = os.path.join(workdir, f"{sample_name}_pass_scoped.vcf.gz")
        scoped_truth = os.path.join(workdir, f"{sample_name}_truth_scoped.vcf.gz")
        subprocess.run(
            f"bcftools view -R {truth_bed} -Oz -o {scoped_pass} {pass_vcf} && bcftools index {scoped_pass}",
            shell=True, check=True, capture_output=True)
        subprocess.run(
            f"bcftools view -R {truth_bed} -Oz -o {scoped_truth} {truth_vcf} && bcftools index {scoped_truth}",
            shell=True, check=True, capture_output=True)
        call_vcf = scoped_pass
        ref_vcf = scoped_truth
    else:
        call_vcf = pass_vcf
        # Need to index truth
        truth_idx = os.path.join(workdir, f"{sample_name}_truth.vcf.gz")
        subprocess.run(
            f"cp {truth_vcf} {truth_idx} && bcftools index {truth_idx}",
            shell=True, check=True, capture_output=True)
        ref_vcf = truth_idx

    # Step 3: bcftools isec
    isec_dir = os.path.join(workdir, f"{sample_name}_isec")
    os.makedirs(isec_dir, exist_ok=True)
    cmd_isec = f"bcftools isec -p {isec_dir} {call_vcf} {ref_vcf}"
    subprocess.run(cmd_isec, shell=True, check=True, capture_output=True)

    # 0000.vcf = in call but not truth (FP)
    # 0001.vcf = in truth but not call (FN)
    # 0002.vcf = in both from call (TP)
    # 0003.vcf = in both from truth
    fp_count = int(subprocess.run(
        f"grep -v '^#' {isec_dir}/0000.vcf | wc -l",
        shell=True, capture_output=True, text=True).stdout.strip())
    fn_count = int(subprocess.run(
        f"grep -v '^#' {isec_dir}/0001.vcf | wc -l",
        shell=True, capture_output=True, text=True).stdout.strip())
    tp_count = int(subprocess.run(
        f"grep -v '^#' {isec_dir}/0002.vcf | wc -l",
        shell=True, capture_output=True, text=True).stdout.strip())

    prec = tp_count / max(tp_count + fp_count, 1)
    rec = tp_count / max(tp_count + fn_count, 1)
    f1 = 2 * prec * rec / max(prec + rec, 1e-10)

    print(f"  TP={tp_count:,}  FP={fp_count:,}  FN={fn_count:,}")
    print(f"  Precision={prec:.4f}  Recall={rec:.4f}  F1={f1:.4f}")

    # Step 4: Extract VCF features from TP and FP
    tp_features = extract_vcf_features(f"{isec_dir}/0002.vcf", sample_name, is_tp=True)
    fp_features = extract_vcf_features(f"{isec_dir}/0000.vcf", sample_name, is_tp=False)

    result = {
        'sample': sample_name,
        'source': config['source'],
        'total_pass': pass_count,
        'tp': tp_count, 'fp': fp_count, 'fn': fn_count,
        'precision': round(prec, 4), 'recall': round(rec, 4), 'f1': round(f1, 4),
        'truth_total': truth_total,
    }

    return result, tp_features + fp_features


def extract_vcf_features(vcf_path, sample_name, is_tp):
    """Extract key VCF features from TP or FP VCF."""
    records = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            qual = float(fields[5]) if fields[5] != '.' else 0.0
            filt = fields[6]

            # Parse INFO
            info = {}
            for item in fields[7].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info[k] = v
                else:
                    info[item] = True

            # Parse FORMAT
            fmt_keys = fields[8].split(':')
            fmt_vals = fields[9].split(':')
            fmt = dict(zip(fmt_keys, fmt_vals))

            try:
                af = float(fmt.get('AF', 0))
            except (ValueError, TypeError):
                af = 0.0
            try:
                dp = int(fmt.get('DP', 0))
            except (ValueError, TypeError):
                dp = 0
            try:
                gq = float(fmt.get('GQ', 0))
            except (ValueError, TypeError):
                gq = 0.0

            gt = fmt.get('GT', '.')
            sb_val = 1.0
            if 'SB' in info and isinstance(info['SB'], str):
                try:
                    sb_val = float(info['SB'])
                except:
                    pass

            records.append({
                'sample': sample_name,
                'chrom': chrom, 'pos': pos,
                'is_tp': 1 if is_tp else 0,
                'qual': qual, 'gq': gq, 'dp': dp, 'af': af,
                'sb': sb_val,
                'phased': 1 if '|' in gt else 0,
                'haplotype_flag': 1 if 'H' in info else 0,
            })
    return records


def compute_auc(labels, scores):
    if len(set(labels)) < 2 or len(set(scores)) <= 1:
        return 0.5
    pairs = sorted(zip(scores, labels), reverse=True)
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    tp, fp = 0, 0
    auc = 0.0
    prev_score = None
    for score, label in pairs:
        if score != prev_score:
            prev_score = score
        if label == 1:
            tp += 1
        else:
            fp += 1
            auc += tp
    return auc / (n_pos * n_neg)


# ═══════════════════════════════════════
# Main
# ═══════════════════════════════════════

if __name__ == '__main__':
    workdir = tempfile.mkdtemp(prefix='to_benchmark_')
    print(f"Working directory: {workdir}")

    all_results = {}
    all_features = []

    for sample_name, config in SAMPLES.items():
        if not os.path.exists(config['to_vcf']):
            print(f"[SKIP] {sample_name}: TO VCF not found")
            continue

        result, features = benchmark_to_vcf(sample_name, config, workdir)
        all_results[sample_name] = result
        all_features.extend(features)

    # ═══════════════════════════════════════
    # Summary Table
    # ═══════════════════════════════════════
    print("\n" + "=" * 90)
    print("CROSS-SAMPLE SUMMARY: ClairS-TO PASS Benchmark")
    print("=" * 90)
    print(f"{'Sample':<18} {'Source':<15} {'PASS':>8} {'TP':>8} {'FP':>8} {'FN':>8} {'Prec':>8} {'Rec':>8} {'F1':>8}")
    print("-" * 90)
    for name, r in all_results.items():
        print(f"{name:<18} {r['source']:<15} {r['total_pass']:>8} {r['tp']:>8} {r['fp']:>8} {r['fn']:>8} {r['precision']:>8.4f} {r['recall']:>8.4f} {r['f1']:>8.4f}")

    # ═══════════════════════════════════════
    # Per-Sample Feature Statistics
    # ═══════════════════════════════════════
    df = pd.DataFrame(all_features)
    print(f"\nTotal feature records: {len(df)}")

    print("\n" + "=" * 90)
    print("PER-SAMPLE FEATURE DISCRIMINATION (AUC for TP vs FP)")
    print("=" * 90)
    print(f"{'Sample':<18} {'AF':>8} {'QUAL':>8} {'GQ':>8} {'DP':>8} {'SB':>8} {'Phased':>8} {'H-flag':>8}")
    print("-" * 90)

    feature_cols = ['af', 'qual', 'gq', 'dp', 'sb', 'phased', 'haplotype_flag']
    per_sample_auc = {}
    for sample_name in all_results:
        sub = df[df['sample'] == sample_name]
        labels = sub['is_tp'].values
        aucs = {}
        for feat in feature_cols:
            vals = sub[feat].values.astype(float)
            valid = ~np.isnan(vals)
            if valid.sum() < 50:
                aucs[feat] = np.nan
                continue
            aucs[feat] = compute_auc(labels[valid].tolist(), vals[valid].tolist())
        per_sample_auc[sample_name] = aucs
        print(f"{sample_name:<18} {aucs['af']:>8.4f} {aucs['qual']:>8.4f} {aucs['gq']:>8.4f} {aucs['dp']:>8.4f} {aucs['sb']:>8.4f} {aucs['phased']:>8.4f} {aucs['haplotype_flag']:>8.4f}")

    # AF statistics per sample
    print("\n" + "=" * 90)
    print("PER-SAMPLE AF DISTRIBUTION")
    print("=" * 90)
    print(f"{'Sample':<18} {'TP AF mean':>12} {'FP AF mean':>12} {'TP AF med':>12} {'FP AF med':>12} {'Δ(FP-TP)':>10}")
    print("-" * 90)
    for sample_name in all_results:
        sub = df[df['sample'] == sample_name]
        tp_af = sub[sub['is_tp']==1]['af']
        fp_af = sub[sub['is_tp']==0]['af']
        delta = fp_af.mean() - tp_af.mean()
        print(f"{sample_name:<18} {tp_af.mean():>12.4f} {fp_af.mean():>12.4f} {tp_af.median():>12.4f} {fp_af.median():>12.4f} {delta:>10.4f}")

    # Phasing rate per sample
    print("\n" + "=" * 90)
    print("PER-SAMPLE PHASING RATE")
    print("=" * 90)
    print(f"{'Sample':<18} {'TP phased%':>12} {'FP phased%':>12} {'Δ(FP-TP)':>10}")
    print("-" * 90)
    for sample_name in all_results:
        sub = df[df['sample'] == sample_name]
        tp_ph = sub[sub['is_tp']==1]['phased'].mean()
        fp_ph = sub[sub['is_tp']==0]['phased'].mean()
        delta = fp_ph - tp_ph
        print(f"{sample_name:<18} {100*tp_ph:>12.1f} {100*fp_ph:>12.1f} {delta:>10.4f}")

    # ═══════════════════════════════════════
    # Save Outputs
    # ═══════════════════════════════════════
    print("\n" + "=" * 60)
    print("Saving Outputs")
    print("=" * 60)

    # Save summary JSON
    summary = {
        'analysis': 'ClairS-TO multi-sample benchmark',
        'samples': all_results,
        'per_sample_auc': {k: {f: round(v, 4) for f, v in aucs.items()} for k, aucs in per_sample_auc.items()},
    }
    summary_path = os.path.join(OUTDIR, "multi_sample_to_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  Saved: {summary_path}")

    # Save features CSV
    features_path = os.path.join(OUTDIR, "multi_sample_to_features.csv")
    df.to_csv(features_path, index=False)
    print(f"  Saved: {features_path} ({len(df)} rows)")

    print(f"\n[DONE] Multi-sample benchmark complete.")
    print(f"  Working directory: {workdir}")
