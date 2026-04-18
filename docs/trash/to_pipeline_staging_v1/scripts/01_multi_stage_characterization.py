#!/usr/bin/env python3
"""
Multi-Stage TO Pipeline TP/FP Characterization
================================================
Analyzes TP/FP at each stage of the ClairS-TO → LongPhase-TO → ISM pipeline,
using VCF features, LOH/CNV annotation, and ISM methylation features.

Stages:
  Stage 1 (ClairS-TO):  All called variants with VCF features → label TP/FP
  Stage 2 (LongPhase):  Add phasing features (PS, phased/unphased)
  Stage 3 (ISM):        Add ISM methylation features from significance_summary.csv

Output:
  - Per-stage F1 metrics
  - Per-feature AUC at each stage
  - Multi-modal characterization plots

Usage:
  python3 01_multi_stage_characterization.py
"""

import os, sys, gzip, csv, json
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

# ── Data paths (HCC1395) ──
CLAIRS_TO_VCF = "/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz"
LONGPHASE_DIR = "/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/tmp/phasing_output/phased_vcf_output"
TRUTH_VCF = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
HC_BED = "/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
CNV_LOH_BED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"

# ISM outputs (paired mode as reference, TO three_way_comparison)
ISM_TO_SIG = "/big7_disk/liaoyoyo2001/big7_disk_output/big8_output_archive/three_way_comparison/tumor_only_full/significance_summary.csv"
ISM_PAIRED_TP = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_pileup/20260314_HCC1395_paired_pileup_pileup_complete_matrix/intersubmod_tp/significance_summary.csv"
ISM_PAIRED_FP = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_pileup/20260314_HCC1395_paired_pileup_pileup_complete_matrix/intersubmod_fp/significance_summary.csv"

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")

# ── Helper: parse VCF ──
def parse_vcf_gz(path, pass_only=False):
    """Parse a VCF.gz, return list of dicts with extracted features."""
    records = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                continue
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
            qual = float(fields[5]) if fields[5] != '.' else 0.0
            filt = fields[6]

            if pass_only and filt != 'PASS':
                continue

            # Parse INFO
            info_str = fields[7]
            info = {}
            for item in info_str.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info[k] = v
                else:
                    info[item] = True

            # Parse FORMAT
            fmt_keys = fields[8].split(':') if len(fields) > 8 else []
            fmt_vals = fields[9].split(':') if len(fields) > 9 else []
            fmt = dict(zip(fmt_keys, fmt_vals))

            rec = {
                'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                'qual': qual, 'filter': filt,
                # FORMAT fields
                'gt': fmt.get('GT', '.'),
                'gq': float(fmt.get('GQ', 0)),
                'dp': int(fmt.get('DP', 0)),
                'af': float(fmt.get('AF', 0)),
                'ps': fmt.get('PS', '.'),
                # INFO flags
                'pon_1': 'PoN_1' in info,
                'pon_2': 'PoN_2' in info,
                'pon_3': 'PoN_3' in info,
                'pon_4': 'PoN_4' in info,
                'verdict_germline': 'Verdict_Germline' in info,
                'verdict_somatic': 'Verdict_Somatic' in info,
                'verdict_subclonal': 'Verdict_SubclonalSomatic' in info,
                'haplotype_flag': 'H' in info,
                'sb': float(info.get('SB', 1.0)) if 'SB' in info and info['SB'] is not True else 1.0,
                # Strand counts
                'fau': int(info.get('FAU', 0)) if 'FAU' in info and info['FAU'] is not True else 0,
                'fcu': int(info.get('FCU', 0)) if 'FCU' in info and info['FCU'] is not True else 0,
                'fgu': int(info.get('FGU', 0)) if 'FGU' in info and info['FGU'] is not True else 0,
                'ftu': int(info.get('FTU', 0)) if 'FTU' in info and info['FTU'] is not True else 0,
                'rau': int(info.get('RAU', 0)) if 'RAU' in info and info['RAU'] is not True else 0,
                'rcu': int(info.get('RCU', 0)) if 'RCU' in info and info['RCU'] is not True else 0,
                'rgu': int(info.get('RGU', 0)) if 'RGU' in info and info['RGU'] is not True else 0,
                'rtu': int(info.get('RTU', 0)) if 'RTU' in info and info['RTU'] is not True else 0,
            }
            # Derived features
            fwd = rec['fau'] + rec['fcu'] + rec['fgu'] + rec['ftu']
            rev = rec['rau'] + rec['rcu'] + rec['rgu'] + rec['rtu']
            rec['strand_ratio'] = fwd / max(fwd + rev, 1)
            rec['pon_count'] = sum([rec['pon_1'], rec['pon_2'], rec['pon_3'], rec['pon_4']])
            rec['phased'] = '|' in rec['gt']
            records.append(rec)
    return records

def load_truth_set(path):
    """Load truth VCF as set of (chrom, pos, ref, alt)."""
    truth = set()
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            truth.add((fields[0], int(fields[1]), fields[3], fields[4]))
    return truth

def load_hc_regions(path):
    """Load high-confidence BED regions as list of (chrom, start, end)."""
    regions = []
    with open(path) as f:
        for line in f:
            fields = line.strip().split('\t')
            regions.append((fields[0], int(fields[1]), int(fields[2])))
    return regions

def load_cnv_loh_bed(path):
    """Load CNV/LOH BED as list of (chrom, start, end, type)."""
    regions = []
    with open(path) as f:
        for line in f:
            fields = line.strip().split('\t')
            regions.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    return regions

def annotate_cnv_loh(records, cnv_loh_regions):
    """Annotate each record with CNV/LOH status using interval search."""
    # Build interval tree per chrom
    from collections import defaultdict
    tree = defaultdict(list)
    for chrom, start, end, rtype in cnv_loh_regions:
        tree[chrom].append((start, end, rtype))
    # Sort for binary search
    for chrom in tree:
        tree[chrom].sort()

    for rec in records:
        pos = rec['pos']
        rec['cnv_type'] = 'neutral'
        rec['in_loh'] = False
        rec['in_gain'] = False
        rec['in_loss'] = False
        for start, end, rtype in tree.get(rec['chrom'], []):
            if start > pos:
                break
            if start <= pos <= end:
                rec['cnv_type'] = rtype
                if rtype == 'loh':
                    rec['in_loh'] = True
                elif rtype == 'gain':
                    rec['in_gain'] = True
                elif rtype == 'loss':
                    rec['in_loss'] = True
    return records

def load_ism_significance(path):
    """Load ISM significance_summary.csv, return dict keyed by (chr, pos)."""
    ism = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row['Chr'], int(row['Pos']))
            ism[key] = row
    return ism

def compute_f1(tp, fp, fn):
    prec = tp / max(tp + fp, 1)
    rec = tp / max(tp + fn, 1)
    f1 = 2 * prec * rec / max(prec + rec, 1e-10)
    return {'tp': tp, 'fp': fp, 'fn': fn, 'precision': prec, 'recall': rec, 'f1': f1}

def compute_auc(labels, scores):
    """Simple AUC computation without sklearn."""
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

# ═══════════════════════════════════════════════
# PART A: Load Data
# ═══════════════════════════════════════════════
print("=" * 70)
print("PART A: Loading data sources")
print("=" * 70)

# A1: Truth set
print("\n[A1] Loading truth set...", flush=True)
truth_set = load_truth_set(TRUTH_VCF)
print(f"  Truth somatic SNVs: {len(truth_set)}")

# A2: CNV/LOH BED
print("\n[A2] Loading CNV/LOH BED...", flush=True)
cnv_loh = load_cnv_loh_bed(CNV_LOH_BED)
n_loh = sum(1 for r in cnv_loh if r[3] == 'loh')
n_gain = sum(1 for r in cnv_loh if r[3] == 'gain')
n_loss = sum(1 for r in cnv_loh if r[3] == 'loss')
print(f"  CNV/LOH regions: {len(cnv_loh)} (gain={n_gain}, loss={n_loss}, loh={n_loh})")

# A3: ClairS-TO final VCF (all variants)
print("\n[A3] Loading ClairS-TO final VCF (all variants)...", flush=True)
all_variants = parse_vcf_gz(CLAIRS_TO_VCF)
print(f"  Total variants: {len(all_variants)}")

# Label TP/FP
for rec in all_variants:
    key = (rec['chrom'], rec['pos'], rec['ref'], rec['alt'])
    rec['is_truth'] = key in truth_set

# Filter categories
filter_counts = Counter(rec['filter'] for rec in all_variants)
pass_variants = [r for r in all_variants if r['filter'] == 'PASS']
nonsomatic = [r for r in all_variants if 'NonSomatic' in r['filter']]
lowqual = [r for r in all_variants if r['filter'] != 'PASS' and 'NonSomatic' not in r['filter']]

print(f"  PASS: {len(pass_variants)}")
print(f"  NonSomatic: {len(nonsomatic)}")
print(f"  LowQual/Other: {len(lowqual)}")

# Annotate CNV/LOH
print("\n[A4] Annotating CNV/LOH...", flush=True)
all_variants = annotate_cnv_loh(all_variants, cnv_loh)

# A5: Load ISM significance summary (TO mode)
print("\n[A5] Loading ISM significance summary (TO mode)...", flush=True)
ism_to = load_ism_significance(ISM_TO_SIG)
print(f"  ISM TO regions: {len(ism_to)}")

# A6: Load LongPhase phased data (merge all chr)
print("\n[A6] Loading LongPhase phased VCF data...", flush=True)
phased_info = {}  # (chrom, pos) -> {ps, phased}
import glob
phased_files = sorted(glob.glob(os.path.join(LONGPHASE_DIR, "tumor_phased_chr*.vcf.gz")))
for pf in phased_files:
    with gzip.open(pf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            fmt_keys = fields[8].split(':')
            fmt_vals = fields[9].split(':')
            fmt = dict(zip(fmt_keys, fmt_vals))
            gt = fmt.get('GT', '.')
            ps = fmt.get('PS', '.')
            phased_info[(chrom, pos)] = {
                'lp_phased': '|' in gt,
                'lp_ps': ps,
                'lp_gt': gt
            }
print(f"  LongPhase phased records: {len(phased_info)}")

# ═══════════════════════════════════════════════
# PART B: Stage-by-Stage F1 Analysis
# ═══════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART B: Stage-by-Stage F1 Analysis")
print("=" * 70)

truth_total = len(truth_set)

# Stage 1: ClairS-TO PASS
pass_tp = sum(1 for r in pass_variants if r['is_truth'])
pass_fp = sum(1 for r in pass_variants if not r['is_truth'])
pass_fn = truth_total - pass_tp
stage1 = compute_f1(pass_tp, pass_fp, pass_fn)
print(f"\n[Stage 1] ClairS-TO PASS variants:")
print(f"  TP={stage1['tp']}, FP={stage1['fp']}, FN={stage1['fn']}")
print(f"  Precision={stage1['precision']:.4f}, Recall={stage1['recall']:.4f}, F1={stage1['f1']:.4f}")

# Also check: how many truth variants are in NonSomatic?
truth_in_nonsomatic = sum(1 for r in nonsomatic if r['is_truth'])
truth_in_lowqual = sum(1 for r in lowqual if r['is_truth'])
print(f"\n  Truth variants lost to NonSomatic filter: {truth_in_nonsomatic}")
print(f"  Truth variants lost to LowQual filter: {truth_in_lowqual}")
print(f"  Truth variants not called at all: {pass_fn - truth_in_nonsomatic - truth_in_lowqual}")

# Stage 2: LongPhase adds phasing info - same variant set, different features
# Check how many PASS variants are phased
pass_phased = sum(1 for r in pass_variants if (r['chrom'], r['pos']) in phased_info and phased_info[(r['chrom'], r['pos'])]['lp_phased'])
pass_unphased = len(pass_variants) - pass_phased
print(f"\n[Stage 2] LongPhase-TO Phasing:")
print(f"  PASS variants phased: {pass_phased} ({100*pass_phased/max(len(pass_variants),1):.1f}%)")
print(f"  PASS variants unphased: {pass_unphased} ({100*pass_unphased/max(len(pass_variants),1):.1f}%)")

# Check if phasing correlates with TP/FP
tp_phased = sum(1 for r in pass_variants if r['is_truth'] and (r['chrom'], r['pos']) in phased_info and phased_info[(r['chrom'], r['pos'])]['lp_phased'])
fp_phased = sum(1 for r in pass_variants if not r['is_truth'] and (r['chrom'], r['pos']) in phased_info and phased_info[(r['chrom'], r['pos'])]['lp_phased'])
tp_total = sum(1 for r in pass_variants if r['is_truth'])
fp_total = sum(1 for r in pass_variants if not r['is_truth'])
print(f"  TP phased: {tp_phased}/{tp_total} ({100*tp_phased/max(tp_total,1):.1f}%)")
print(f"  FP phased: {fp_phased}/{fp_total} ({100*fp_phased/max(fp_total,1):.1f}%)")

# Stage 3: ISM adds methylation features
# Match PASS variants to ISM regions
pass_with_ism = 0
for r in pass_variants:
    key = (r['chrom'], r['pos'])
    if key in ism_to:
        r['has_ism'] = True
        r['ism'] = ism_to[key]
        pass_with_ism += 1
    else:
        r['has_ism'] = False

ism_tp = sum(1 for r in pass_variants if r['has_ism'] and r['is_truth'])
ism_fp = sum(1 for r in pass_variants if r['has_ism'] and not r['is_truth'])
ism_fn = truth_total - ism_tp
stage3 = compute_f1(ism_tp, ism_fp, ism_fn)

# ISM SuggestFilter
ism_filtered_tp = 0
ism_filtered_fp = 0
for r in pass_variants:
    if r['has_ism']:
        sf = r['ism'].get('SuggestFilter', 'false')
        if sf == 'true':
            if r['is_truth']:
                ism_filtered_tp += 1
            else:
                ism_filtered_fp += 1

stage3_filtered = compute_f1(ism_tp - ism_filtered_tp, ism_fp - ism_filtered_fp, truth_total - (ism_tp - ism_filtered_tp))

print(f"\n[Stage 3] ISM Analysis:")
print(f"  PASS variants with ISM data: {pass_with_ism}/{len(pass_variants)}")
print(f"  ISM TP={ism_tp}, FP={ism_fp}")
print(f"  ISM SuggestFilter removes: TP={ism_filtered_tp}, FP={ism_filtered_fp}")
print(f"  After ISM filter: TP={stage3_filtered['tp']}, FP={stage3_filtered['fp']}")
print(f"  Precision={stage3_filtered['precision']:.4f}, Recall={stage3_filtered['recall']:.4f}, F1={stage3_filtered['f1']:.4f}")

# Summary table
print(f"\n{'='*70}")
print(f"Stage-by-Stage Summary:")
print(f"{'Stage':<30} {'TP':>7} {'FP':>7} {'FN':>7} {'Prec':>8} {'Rec':>8} {'F1':>8}")
print(f"{'-'*70}")
print(f"{'1. ClairS-TO PASS':<30} {stage1['tp']:>7} {stage1['fp']:>7} {stage1['fn']:>7} {stage1['precision']:>8.4f} {stage1['recall']:>8.4f} {stage1['f1']:>8.4f}")
print(f"{'2. + LongPhase (same set)':<30} {stage1['tp']:>7} {stage1['fp']:>7} {stage1['fn']:>7} {stage1['precision']:>8.4f} {stage1['recall']:>8.4f} {stage1['f1']:>8.4f}")
print(f"{'3. + ISM (before filter)':<30} {stage3['tp']:>7} {stage3['fp']:>7} {stage3['fn']:>7} {stage3['precision']:>8.4f} {stage3['recall']:>8.4f} {stage3['f1']:>8.4f}")
print(f"{'3b. + ISM (after filter)':<30} {stage3_filtered['tp']:>7} {stage3_filtered['fp']:>7} {stage3_filtered['fn']:>7} {stage3_filtered['precision']:>8.4f} {stage3_filtered['recall']:>8.4f} {stage3_filtered['f1']:>8.4f}")

# ═══════════════════════════════════════════════
# PART C: Per-Feature AUC at Each Stage
# ═══════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART C: Per-Feature TP vs FP Characterization")
print("=" * 70)

# Build DataFrame for PASS variants
rows = []
for r in pass_variants:
    row = {
        'chrom': r['chrom'], 'pos': r['pos'],
        'is_tp': 1 if r['is_truth'] else 0,
        # VCF features (Stage 1)
        'qual': r['qual'],
        'gq': r['gq'],
        'dp': r['dp'],
        'af': r['af'],
        'sb': r['sb'],
        'strand_ratio': r['strand_ratio'],
        'pon_count': r['pon_count'],
        'verdict_germline': int(r['verdict_germline']),
        'verdict_somatic': int(r['verdict_somatic']),
        'verdict_subclonal': int(r['verdict_subclonal']),
        'haplotype_flag': int(r['haplotype_flag']),
        # CNV/LOH (annotation)
        'in_loh': int(r['in_loh']),
        'in_gain': int(r['in_gain']),
        'in_loss': int(r['in_loss']),
        'cnv_type': r['cnv_type'],
    }
    # Phasing features (Stage 2)
    pi = phased_info.get((r['chrom'], r['pos']), {})
    row['lp_phased'] = int(pi.get('lp_phased', False))
    # Phase set block size will be computed later

    # ISM features (Stage 3)
    if r['has_ism']:
        ism = r['ism']
        for col in ['NumReads', 'NumCpGs', 'CramersV', 'HeuristicScore',
                     'HPMergedDelta', 'AlleleDelta', 'UnassignedAffinity',
                     'HP_Ratio', 'Quality_Score', 'HPFineNGroups',
                     'HP1FamilyN', 'HP2FamilyN', 'Coverage_Multiple']:
            try:
                row[f'ism_{col}'] = float(ism.get(col, np.nan))
            except (ValueError, TypeError):
                row[f'ism_{col}'] = np.nan
        for col in ['HPMergedSig', 'AlleleSig', 'PassedGating', 'Significant',
                     'Potential_LOH', 'SuggestFilter']:
            val = ism.get(col, '')
            row[f'ism_{col}'] = 1 if val.lower() == 'true' else 0
        row[f'ism_VerificationClass'] = ism.get('VerificationClass', '')
        row[f'ism_Quality_Tier'] = ism.get('Quality_Tier', '')
        row[f'ism_DominantLabel'] = ism.get('DominantLabel', '')
        row['has_ism'] = 1
    else:
        row['has_ism'] = 0

    rows.append(row)

df = pd.DataFrame(rows)
print(f"\nTotal PASS variants in DataFrame: {len(df)}")
print(f"  TP: {df['is_tp'].sum()}, FP: {(1-df['is_tp']).sum()}")

# Compute Phase Set block sizes
print("\n[C1] Computing phase set block sizes...", flush=True)
ps_blocks = defaultdict(list)
for (chrom, pos), info in phased_info.items():
    if info['lp_phased'] and info['lp_ps'] != '.':
        ps_blocks[(chrom, info['lp_ps'])].append(pos)

ps_block_size = {}
for (chrom, ps), positions in ps_blocks.items():
    block_size = len(positions)
    for p in positions:
        ps_block_size[(chrom, p)] = block_size

df['ps_block_size'] = df.apply(lambda r: ps_block_size.get((r['chrom'], r['pos']), 0), axis=1)
print(f"  Phase blocks: {len(ps_blocks)}")
print(f"  Median block size: {np.median([len(v) for v in ps_blocks.values()]):.0f}")

# ── Feature AUC computation ──
print("\n[C2] Per-Feature AUC (TP=positive):", flush=True)

# Stage 1 features (VCF)
vcf_features = ['qual', 'gq', 'dp', 'af', 'sb', 'strand_ratio', 'pon_count',
                 'verdict_germline', 'verdict_somatic', 'verdict_subclonal', 'haplotype_flag']

# Stage 2 features (Phasing + CNV/LOH annotation)
phasing_features = ['lp_phased', 'ps_block_size', 'in_loh', 'in_gain', 'in_loss']

# Stage 3 features (ISM)
ism_features = ['ism_NumReads', 'ism_NumCpGs', 'ism_CramersV', 'ism_HeuristicScore',
                'ism_HPMergedDelta', 'ism_AlleleDelta', 'ism_UnassignedAffinity',
                'ism_HP_Ratio', 'ism_Quality_Score', 'ism_HPFineNGroups',
                'ism_HP1FamilyN', 'ism_HP2FamilyN', 'ism_Coverage_Multiple',
                'ism_HPMergedSig', 'ism_AlleleSig', 'ism_PassedGating',
                'ism_Significant', 'ism_Potential_LOH', 'ism_SuggestFilter']

all_features = vcf_features + phasing_features + ism_features

auc_results = {}
labels_all = df['is_tp'].values

print(f"\n{'Feature':<30} {'AUC':>7} {'TP mean':>10} {'FP mean':>10} {'N valid':>8} {'Stage':>8}")
print(f"{'-'*73}")

for feat in all_features:
    if feat not in df.columns:
        continue
    vals = df[feat].values
    valid = ~np.isnan(vals.astype(float))
    if valid.sum() < 100:
        continue
    v = vals[valid].astype(float)
    l = labels_all[valid]

    auc = compute_auc(l.tolist(), v.tolist())
    tp_mean = v[l == 1].mean() if (l == 1).sum() > 0 else np.nan
    fp_mean = v[l == 0].mean() if (l == 0).sum() > 0 else np.nan

    stage = "S1-VCF" if feat in vcf_features else ("S2-Phase" if feat in phasing_features else "S3-ISM")
    auc_results[feat] = {'auc': auc, 'tp_mean': tp_mean, 'fp_mean': fp_mean, 'n': valid.sum(), 'stage': stage}
    print(f"  {feat:<28} {auc:>7.4f} {tp_mean:>10.4f} {fp_mean:>10.4f} {valid.sum():>8d} {stage:>8}")

# ═══════════════════════════════════════════════
# PART D: Stratified Analysis
# ═══════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART D: Stratified Analysis")
print("=" * 70)

# D1: TP/FP by FILTER category (all variants, not just PASS)
print("\n[D1] TP/FP by FILTER category (all 3.18M variants):")
filter_groups = defaultdict(lambda: {'tp': 0, 'fp': 0})
for r in all_variants:
    filt = r['filter']
    # Simplify filter
    if filt == 'PASS':
        fgroup = 'PASS'
    elif 'NonSomatic' in filt and 'LowQual' in filt:
        fgroup = 'LowQual+NonSomatic'
    elif 'NonSomatic' in filt:
        fgroup = 'NonSomatic'
    else:
        fgroup = 'LowQual/Other'

    if r['is_truth']:
        filter_groups[fgroup]['tp'] += 1
    else:
        filter_groups[fgroup]['fp'] += 1

print(f"  {'Filter':<25} {'TP':>8} {'FP':>10} {'TP rate':>10}")
for fg in ['PASS', 'NonSomatic', 'LowQual+NonSomatic', 'LowQual/Other']:
    d = filter_groups[fg]
    total = d['tp'] + d['fp']
    rate = d['tp'] / max(total, 1)
    print(f"  {fg:<25} {d['tp']:>8} {d['fp']:>10} {rate:>10.6f}")

# D2: TP/FP by LOH/CNV status (PASS variants)
print("\n[D2] PASS TP/FP by CNV/LOH status:")
for cnv_cat in ['neutral', 'gain', 'loss', 'loh']:
    sub = df[df['cnv_type'] == cnv_cat]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {cnv_cat:<10} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# D3: TP/FP by Verdict (PASS variants)
print("\n[D3] PASS TP/FP by ClairS-TO Verdict:")
for v in ['Somatic', 'SubclonalSomatic', 'Germline', 'No Verdict']:
    if v == 'Somatic':
        sub = df[df['verdict_somatic'] == 1]
    elif v == 'SubclonalSomatic':
        sub = df[df['verdict_subclonal'] == 1]
    elif v == 'Germline':
        sub = df[df['verdict_germline'] == 1]
    else:
        sub = df[(df['verdict_somatic']==0) & (df['verdict_subclonal']==0) & (df['verdict_germline']==0)]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {v:<22} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# D4: TP/FP by phasing status (PASS variants)
print("\n[D4] PASS TP/FP by phasing status:")
for ph_status, ph_label in [(1, 'Phased'), (0, 'Unphased')]:
    sub = df[df['lp_phased'] == ph_status]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {ph_label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# D5: ISM VerificationClass distribution
print("\n[D5] ISM VerificationClass (PASS variants with ISM):")
ism_sub = df[df['has_ism'] == 1]
for vc in ['Strong', 'Weak', 'Noise', 'Subclone', 'Insufficient', '']:
    sub = ism_sub[ism_sub['ism_VerificationClass'] == vc]
    if len(sub) == 0:
        continue
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    label = vc if vc else 'Empty'
    print(f"  {label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# D6: ISM Quality_Tier distribution
print("\n[D6] ISM Quality_Tier (PASS variants with ISM):")
for qt in ['HIGH', 'MEDIUM', 'LOW', '']:
    sub = ism_sub[ism_sub['ism_Quality_Tier'] == qt]
    if len(sub) == 0:
        continue
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    label = qt if qt else 'Empty'
    print(f"  {label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# ═══════════════════════════════════════════════
# PART E: Save DataFrame and AUC results
# ═══════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART E: Saving outputs")
print("=" * 70)

# Save PASS variants DataFrame
out_csv = os.path.join(OUTDIR, "hcc1395_pass_multimodal.csv")
df.to_csv(out_csv, index=False)
print(f"  Saved: {out_csv} ({len(df)} rows)")

# Save AUC results
auc_df = pd.DataFrame([
    {'feature': k, **v} for k, v in auc_results.items()
]).sort_values('auc', ascending=False)
auc_csv = os.path.join(OUTDIR, "hcc1395_feature_auc_by_stage.csv")
auc_df.to_csv(auc_csv, index=False)
print(f"  Saved: {auc_csv} ({len(auc_df)} features)")

# Save stage summary
stages_json = {
    'sample': 'HCC1395',
    'truth_total': truth_total,
    'stages': {
        'stage1_clairs_to_pass': stage1,
        'stage2_longphase': stage1,  # same F1, adds features
        'stage3_ism_before_filter': {
            'tp': stage3['tp'], 'fp': stage3['fp'], 'fn': stage3['fn'],
            'precision': stage3['precision'], 'recall': stage3['recall'], 'f1': stage3['f1']
        },
        'stage3b_ism_after_filter': {
            'tp': stage3_filtered['tp'], 'fp': stage3_filtered['fp'], 'fn': stage3_filtered['fn'],
            'precision': stage3_filtered['precision'], 'recall': stage3_filtered['recall'], 'f1': stage3_filtered['f1']
        }
    },
    'filter_breakdown': dict(filter_groups),
    'truth_in_nonsomatic': truth_in_nonsomatic,
    'truth_in_lowqual': truth_in_lowqual,
}
stages_json_path = os.path.join(OUTDIR, "hcc1395_stage_metrics.json")
with open(stages_json_path, 'w') as f:
    json.dump(stages_json, f, indent=2, default=str)
print(f"  Saved: {stages_json_path}")

print("\n[DONE] Part A-E complete. Run 02_generate_plots.py for figures.")
