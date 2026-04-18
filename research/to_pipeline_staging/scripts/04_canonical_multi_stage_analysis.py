#!/usr/bin/env python3
"""
Canonical Multi-Stage TO Pipeline TP/FP Characterization (v2)
=============================================================
CORRECTED version using canonical pipeline data instead of three_way_comparison.

Data sources (HCC1395 canonical TO run):
  - ClairS-TO VCF: archive/step01_clairs_to/snv.vcf.gz (88MB, ONT_5kHz BAM, Mar 2026)
  - Benchmark TP/FP: archive/step02_benchmark_clairs_to/{tp,fp}.vcf
  - LongPhase-TO phased VCF: archive/step03_longphase_to/tumor_phased.vcf.gz
  - LongPhase-TO LOH.bed: archive/step03_longphase_to/tumor_phased_LOH.bed
  - ISM TP: archive/step05_intersubmod/intersubmod_tp/significance_summary.csv
  - ISM FP: archive/step05_intersubmod/intersubmod_fp/significance_summary.csv

Key correction from v1:
  - v1 used zhenyu112 VCF (47,798 PASS, Sep 2025) + three_way ISM (112K rows)
  - v2 uses canonical VCF (48,085 PASS, Mar 2026) + canonical ISM (40,213 rows)
  - v1 had data alignment error causing inflated ISM implicit filter effect
  - v2 uses pre-benchmarked TP/FP VCFs → no BED filtering ambiguity

Output:
  data/hcc1395_canonical_multimodal.csv     — all features per variant
  data/hcc1395_canonical_feature_auc.csv    — per-feature AUC
  data/hcc1395_canonical_stage_metrics.json  — per-stage F1
"""

import os, sys, gzip, csv, json
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

# ═══════════════════════════════════════
# Data Paths
# ═══════════════════════════════════════
ARCHIVE = "/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1"

# Step 01: ClairS-TO output
CLAIRS_TO_VCF = f"{ARCHIVE}/step01_clairs_to/snv.vcf.gz"

# Step 02: Benchmark (pre-classified TP/FP within truth BED)
TP_VCF = f"{ARCHIVE}/step02_benchmark_clairs_to/tp.vcf"
FP_VCF = f"{ARCHIVE}/step02_benchmark_clairs_to/fp.vcf"
VARIANT_COUNTS = f"{ARCHIVE}/step02_benchmark_clairs_to/variant_counts.txt"

# Step 03: LongPhase-TO phased VCF + LOH/GE beds
PHASED_VCF = f"{ARCHIVE}/step03_longphase_to/tumor_phased.vcf.gz"
LOH_BED = f"{ARCHIVE}/step03_longphase_to/tumor_phased_LOH.bed"
GE_BED = f"{ARCHIVE}/step03_longphase_to/tumor_phased_GE.bed"

# Step 05: ISM significance_summary (separate TP/FP)
ISM_TP_SIG = f"{ARCHIVE}/step05_intersubmod/intersubmod_tp/significance_summary.csv"
ISM_FP_SIG = f"{ARCHIVE}/step05_intersubmod/intersubmod_fp/significance_summary.csv"

# External annotation
CNV_LOH_BED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
TRUTH_TOTAL = 39447

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

# ═══════════════════════════════════════
# Helper Functions
# ═══════════════════════════════════════

def parse_vcf_records(path, is_gzipped=False):
    """Parse VCF (gzipped or plain), return list of dicts."""
    records = []
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    with opener(path, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            chrom, pos, _, ref, alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
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

            rec = {
                'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                'qual': qual, 'filter': filt,
                'gt': fmt.get('GT', '.'),
                'gq': float(fmt.get('GQ', 0)),
                'dp': int(fmt.get('DP', 0)),
                'af': float(fmt.get('AF', 0)),
                'ps': fmt.get('PS', '.'),
                'gt2': fmt.get('GT2', '.'),
                'ps2': fmt.get('PS2', '.'),
                # INFO: strand bias
                'sb': float(info.get('SB', 1.0)) if 'SB' in info and info['SB'] is not True else 1.0,
                'haplotype_flag': 'H' in info,
                # INFO: strand counts
                'fau': int(info.get('FAU', 0)) if isinstance(info.get('FAU'), str) else 0,
                'fcu': int(info.get('FCU', 0)) if isinstance(info.get('FCU'), str) else 0,
                'fgu': int(info.get('FGU', 0)) if isinstance(info.get('FGU'), str) else 0,
                'ftu': int(info.get('FTU', 0)) if isinstance(info.get('FTU'), str) else 0,
                'rau': int(info.get('RAU', 0)) if isinstance(info.get('RAU'), str) else 0,
                'rcu': int(info.get('RCU', 0)) if isinstance(info.get('RCU'), str) else 0,
                'rgu': int(info.get('RGU', 0)) if isinstance(info.get('RGU'), str) else 0,
                'rtu': int(info.get('RTU', 0)) if isinstance(info.get('RTU'), str) else 0,
                # INFO: verdict
                'verdict_germline': 'Verdict_Germline' in info,
                'verdict_somatic': 'Verdict_Somatic' in info,
                'verdict_subclonal': 'Verdict_SubclonalSomatic' in info,
                # INFO: PoN
                'pon_1': 'PoN_1' in info,
                'pon_2': 'PoN_2' in info,
                'pon_3': 'PoN_3' in info,
                'pon_4': 'PoN_4' in info,
            }
            fwd = rec['fau'] + rec['fcu'] + rec['fgu'] + rec['ftu']
            rev = rec['rau'] + rec['rcu'] + rec['rgu'] + rec['rtu']
            rec['strand_ratio'] = fwd / max(fwd + rev, 1)
            rec['pon_count'] = sum([rec['pon_1'], rec['pon_2'], rec['pon_3'], rec['pon_4']])
            rec['phased'] = '|' in rec['gt']
            records.append(rec)
    return records

def load_bed_intervals(path):
    """Load BED as dict of chrom -> sorted list of (start, end, label)."""
    from collections import defaultdict
    tree = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            label = fields[3] if len(fields) > 3 else '.'
            tree[chrom].append((start, end, label))
    for chrom in tree:
        tree[chrom].sort()
    return dict(tree)

def annotate_bed(records, bed_tree, prefix):
    """Annotate records with BED region membership."""
    for rec in records:
        pos = rec['pos']
        rec[f'{prefix}_hit'] = False
        rec[f'{prefix}_label'] = 'none'
        for start, end, label in bed_tree.get(rec['chrom'], []):
            if start > pos:
                break
            if start <= pos <= end:
                rec[f'{prefix}_hit'] = True
                rec[f'{prefix}_label'] = label
    return records

def load_ism_significance(path):
    """Load ISM significance_summary.csv → dict keyed by (chr, pos)."""
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
    return {'tp': int(tp), 'fp': int(fp), 'fn': int(fn),
            'precision': round(prec, 6), 'recall': round(rec, 6), 'f1': round(f1, 6)}

def compute_auc(labels, scores):
    """ROC AUC without sklearn."""
    if len(set(labels)) < 2:
        return 0.5
    # Check for constant scores
    if len(set(scores)) <= 1:
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
# PART A: Load All Data Sources
# ═══════════════════════════════════════
print("=" * 70)
print("PART A: Loading Canonical Data Sources")
print("=" * 70)
print(f"  Archive: {ARCHIVE}")

# A1: Verify data files
for label, path in [
    ("ClairS-TO VCF", CLAIRS_TO_VCF),
    ("TP VCF", TP_VCF),
    ("FP VCF", FP_VCF),
    ("Phased VCF", PHASED_VCF),
    ("LOH BED", LOH_BED),
    ("ISM TP sig", ISM_TP_SIG),
    ("ISM FP sig", ISM_FP_SIG),
    ("CNV/LOH BED", CNV_LOH_BED),
]:
    exists = os.path.exists(path)
    print(f"  {label}: {'OK' if exists else 'MISSING'} — {path}")
    if not exists:
        sys.exit(f"[FATAL] Missing: {path}")

# A2: Load benchmark counts
print("\n[A2] Benchmark variant counts:")
counts = {}
with open(VARIANT_COUNTS) as f:
    for line in f:
        k, v = line.strip().split('=')
        counts[k] = int(v)
print(f"  PASS_TOTAL (in BED) = {counts['PASS_TOTAL']}")
print(f"  TP = {counts['TP_COUNT']}, FP = {counts['FP_COUNT']}")
print(f"  FN = {counts['FN_COUNT']}, TRUTH_TOTAL = {counts['TRUTH_TOTAL']}")

# A3: Load ClairS-TO FILTER breakdown
print("\n[A3] ClairS-TO VCF FILTER breakdown (full 3.2M variants):", flush=True)
filter_counts = Counter()
total_vcf = 0
with gzip.open(CLAIRS_TO_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        total_vcf += 1
        filt = line.split('\t')[6]
        # Simplify
        if filt == 'PASS':
            filter_counts['PASS'] += 1
        elif 'NonSomatic' in filt and 'LowQual' in filt:
            filter_counts['LowQual+NonSomatic'] += 1
        elif 'NonSomatic' in filt:
            filter_counts['NonSomatic'] += 1
        elif 'NoAncestry' in filt:
            filter_counts['LowQual+NoAncestry'] += 1
        elif 'MultiHap' in filt:
            filter_counts['LowQual+MultiHap'] += 1
        else:
            filter_counts['LowQual/Other'] += 1

print(f"  Total variants in VCF: {total_vcf:,}")
for filt, count in sorted(filter_counts.items(), key=lambda x: -x[1]):
    print(f"    {filt:<25} {count:>10,} ({100*count/total_vcf:.2f}%)")

# A4: Parse TP and FP VCFs
print("\n[A4] Parsing benchmark TP/FP VCFs...", flush=True)
tp_records = parse_vcf_records(TP_VCF, is_gzipped=False)
fp_records = parse_vcf_records(FP_VCF, is_gzipped=False)
print(f"  TP records: {len(tp_records)}")
print(f"  FP records: {len(fp_records)}")
assert len(tp_records) == counts['TP_COUNT'], f"TP mismatch: {len(tp_records)} vs {counts['TP_COUNT']}"
assert len(fp_records) == counts['FP_COUNT'], f"FP mismatch: {len(fp_records)} vs {counts['FP_COUNT']}"

# Label TP/FP
for r in tp_records:
    r['is_tp'] = True
for r in fp_records:
    r['is_tp'] = False
all_pass_in_bed = tp_records + fp_records

# A5: Load LongPhase-TO phasing info
print("\n[A5] Loading LongPhase-TO phased VCF...", flush=True)
phased_info = {}  # (chrom, pos) -> {phased, ps, gt2, ps2}
with gzip.open(PHASED_VCF, 'rt') as f:
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
        gt2 = fmt.get('GT2', '.')
        ps2 = fmt.get('PS2', '.')
        phased_info[(chrom, pos)] = {
            'lp_phased': '|' in gt,
            'lp_ps': ps,
            'lp_gt': gt,
            'lp_gt2': gt2,
            'lp_ps2': ps2,
        }
print(f"  Phased VCF records loaded: {len(phased_info):,}")
n_phased = sum(1 for v in phased_info.values() if v['lp_phased'])
print(f"  Phased variants: {n_phased:,} ({100*n_phased/len(phased_info):.1f}%)")

# A6: Load LOH/GE BED
print("\n[A6] Loading LongPhase-TO LOH/GE BEDs...", flush=True)
loh_tree = load_bed_intervals(LOH_BED)
ge_tree = load_bed_intervals(GE_BED)
loh_regions = sum(len(v) for v in loh_tree.values())
ge_regions = sum(len(v) for v in ge_tree.values())
print(f"  LOH regions: {loh_regions}")
print(f"  GE regions: {ge_regions}")

# A7: Load CNV/LOH external annotation
print("\n[A7] Loading external CNV/LOH BED...", flush=True)
cnv_tree = load_bed_intervals(CNV_LOH_BED)
cnv_total = sum(len(v) for v in cnv_tree.values())
cnv_types = Counter()
for intervals in cnv_tree.values():
    for _, _, label in intervals:
        cnv_types[label] += 1
print(f"  CNV regions: {cnv_total} ({dict(cnv_types)})")

# A8: Load ISM significance summaries
print("\n[A8] Loading ISM significance summaries (canonical TP/FP)...", flush=True)
ism_tp = load_ism_significance(ISM_TP_SIG)
ism_fp = load_ism_significance(ISM_FP_SIG)
print(f"  ISM TP regions: {len(ism_tp)}")
print(f"  ISM FP regions: {len(ism_fp)}")

# Merge into single dict with labels
ism_all = {}
for key, val in ism_tp.items():
    val['_ism_label'] = 'tp'
    ism_all[key] = val
for key, val in ism_fp.items():
    val['_ism_label'] = 'fp'
    ism_all[key] = val
print(f"  ISM total: {len(ism_all)}")

# ═══════════════════════════════════════
# PART B: Build Feature DataFrame
# ═══════════════════════════════════════
print("\n" + "=" * 70)
print("PART B: Building Feature DataFrame")
print("=" * 70)

# Annotate all records
annotate_bed(all_pass_in_bed, loh_tree, 'lp_loh')
annotate_bed(all_pass_in_bed, cnv_tree, 'cnv')

rows = []
for r in all_pass_in_bed:
    row = {
        'chrom': r['chrom'], 'pos': r['pos'], 'ref': r['ref'], 'alt': r['alt'],
        'is_tp': 1 if r['is_tp'] else 0,
        # VCF features (Stage 1: ClairS-TO)
        'qual': r['qual'], 'gq': r['gq'], 'dp': r['dp'], 'af': r['af'],
        'sb': r['sb'], 'strand_ratio': r['strand_ratio'],
        'pon_count': r['pon_count'],
        'verdict_germline': int(r['verdict_germline']),
        'verdict_somatic': int(r['verdict_somatic']),
        'verdict_subclonal': int(r['verdict_subclonal']),
        'haplotype_flag': int(r['haplotype_flag']),
        # CNV/LOH external annotation
        'cnv_type': r['cnv_label'],
        'in_loh': int(r['cnv_label'] == 'loh'),
        'in_gain': int(r['cnv_label'] == 'gain'),
        'in_loss': int(r['cnv_label'] == 'loss'),
        # LongPhase-TO LOH
        'lp_loh': int(r['lp_loh_hit']),
    }

    # Phasing features (Stage 2: LongPhase-TO)
    pi = phased_info.get((r['chrom'], r['pos']), {})
    row['lp_phased'] = int(pi.get('lp_phased', False))
    row['lp_ps'] = pi.get('lp_ps', '.')

    # ISM features (Stage 3: InterSubMod)
    key = (r['chrom'], r['pos'])
    ism_data = ism_all.get(key)
    if ism_data:
        row['has_ism'] = 1
        # Numeric ISM features
        for col in ['NumReads', 'NumCpGs', 'CramersV', 'HeuristicScore',
                     'HPMergedDelta', 'AlleleDelta', 'UnassignedAffinity',
                     'HP_Ratio', 'Quality_Score', 'HPFineNGroups',
                     'HP1FamilyN', 'HP2FamilyN', 'Coverage_Multiple',
                     'PairwiseMeanDist', 'PairwiseMedianDist']:
            try:
                row[f'ism_{col}'] = float(ism_data.get(col, np.nan))
            except (ValueError, TypeError):
                row[f'ism_{col}'] = np.nan
        # Boolean ISM features
        for col in ['HPMergedSig', 'AlleleSig', 'PassedGating', 'Significant',
                     'Potential_LOH', 'SuggestFilter']:
            val = ism_data.get(col, '')
            row[f'ism_{col}'] = 1 if str(val).lower() == 'true' else 0
        # Categorical ISM features
        row['ism_VerificationClass'] = ism_data.get('VerificationClass', '')
        row['ism_Quality_Tier'] = ism_data.get('Quality_Tier', '')
        row['ism_DominantLabel'] = ism_data.get('DominantLabel', '')
        row['ism_Coverage_Category'] = ism_data.get('Coverage_Category', '')
        row['ism_LOH_Subtype'] = ism_data.get('LOH_Subtype', '')
    else:
        row['has_ism'] = 0

    rows.append(row)

df = pd.DataFrame(rows)
print(f"  Total records: {len(df)}")
print(f"  TP: {df['is_tp'].sum()}, FP: {(1-df['is_tp']).sum():.0f}")
print(f"  With ISM: {df['has_ism'].sum()} ({100*df['has_ism'].mean():.2f}%)")

# Phase set block sizes
print("\n[B1] Computing phase set block sizes...", flush=True)
ps_blocks = defaultdict(list)
for (chrom, pos), info in phased_info.items():
    if info['lp_phased'] and info['lp_ps'] != '.':
        ps_blocks[(chrom, info['lp_ps'])].append(pos)

ps_block_size = {}
for (chrom, ps), positions in ps_blocks.items():
    sz = len(positions)
    for p in positions:
        ps_block_size[(chrom, p)] = sz

df['ps_block_size'] = df.apply(lambda r: ps_block_size.get((r['chrom'], r['pos']), 0), axis=1)
print(f"  Phase blocks: {len(ps_blocks):,}")
if ps_blocks:
    sizes = [len(v) for v in ps_blocks.values()]
    print(f"  Block size: median={np.median(sizes):.0f}, mean={np.mean(sizes):.1f}, max={max(sizes)}")

# ═══════════════════════════════════════
# PART C: Stage-by-Stage Metrics
# ═══════════════════════════════════════
print("\n" + "=" * 70)
print("PART C: Stage-by-Stage Metrics (Canonical)")
print("=" * 70)

# Stage 1: ClairS-TO PASS (within truth BED)
stage1 = compute_f1(counts['TP_COUNT'], counts['FP_COUNT'], counts['FN_COUNT'])
print(f"\n[Stage 1] ClairS-TO PASS (in truth BED):")
print(f"  Total PASS in VCF: {filter_counts['PASS']:,}")
print(f"  PASS in truth BED: {counts['PASS_TOTAL']:,}")
print(f"  TP={stage1['tp']:,}  FP={stage1['fp']:,}  FN={stage1['fn']:,}")
print(f"  Precision={stage1['precision']:.4f}  Recall={stage1['recall']:.4f}  F1={stage1['f1']:.4f}")

# Stage 2: LongPhase-TO (same variant set, adds phasing)
tp_phased = df[df['is_tp']==1]['lp_phased'].sum()
fp_phased = df[df['is_tp']==0]['lp_phased'].sum()
tp_total = df['is_tp'].sum()
fp_total = len(df) - tp_total
print(f"\n[Stage 2] LongPhase-TO Phasing (no FILTER change):")
print(f"  TP phased: {tp_phased}/{tp_total} ({100*tp_phased/tp_total:.1f}%)")
print(f"  FP phased: {fp_phased}/{fp_total} ({100*fp_phased/fp_total:.1f}%)")
print(f"  → LongPhase does NOT change PASS/FILTER → same TP/FP/F1 as Stage 1")
print(f"  → All phasing-related FILTERs (NoAncestry, MultiHap) apply ONLY to LowQual variants")

# Stage 3: ISM processing
ism_tp_count = df[(df['is_tp']==1) & (df['has_ism']==1)].shape[0]
ism_fp_count = df[(df['is_tp']==0) & (df['has_ism']==1)].shape[0]
lost_tp = counts['TP_COUNT'] - ism_tp_count
lost_fp = counts['FP_COUNT'] - ism_fp_count
stage3 = compute_f1(ism_tp_count, ism_fp_count, TRUTH_TOTAL - ism_tp_count)

# ISM SuggestFilter
sf_tp = df[(df['is_tp']==1) & (df.get('ism_SuggestFilter', 0)==1)].shape[0] if 'ism_SuggestFilter' in df.columns else 0
sf_fp = df[(df['is_tp']==0) & (df.get('ism_SuggestFilter', 0)==1)].shape[0] if 'ism_SuggestFilter' in df.columns else 0
stage3f = compute_f1(ism_tp_count - sf_tp, ism_fp_count - sf_fp, TRUTH_TOTAL - (ism_tp_count - sf_tp))

print(f"\n[Stage 3] ISM Processing:")
print(f"  ISM processed TP: {ism_tp_count}/{counts['TP_COUNT']} (lost {lost_tp}, {100*lost_tp/counts['TP_COUNT']:.3f}%)")
print(f"  ISM processed FP: {ism_fp_count}/{counts['FP_COUNT']} (lost {lost_fp}, {100*lost_fp/counts['FP_COUNT']:.3f}%)")
print(f"  ISM implicit filter: removes {lost_tp+lost_fp} variants ({lost_tp} TP + {lost_fp} FP)")
if lost_fp > 0:
    print(f"  Selectivity: {lost_fp/max(lost_tp,1):.1f}× (FP removed per TP lost)")
print(f"\n  Before SuggestFilter: {stage3}")
print(f"  SuggestFilter removes: {sf_tp} TP + {sf_fp} FP = {sf_tp+sf_fp}")
print(f"  After SuggestFilter:  {stage3f}")

# Summary table
print(f"\n{'='*80}")
print(f"{'Stage':<35} {'TP':>7} {'FP':>7} {'FN':>7} {'Prec':>8} {'Rec':>8} {'F1':>8}")
print(f"{'-'*80}")
print(f"{'1. ClairS-TO PASS (in BED)':<35} {stage1['tp']:>7} {stage1['fp']:>7} {stage1['fn']:>7} {stage1['precision']:>8.4f} {stage1['recall']:>8.4f} {stage1['f1']:>8.4f}")
print(f"{'2. + LongPhase-TO (no change)':<35} {stage1['tp']:>7} {stage1['fp']:>7} {stage1['fn']:>7} {stage1['precision']:>8.4f} {stage1['recall']:>8.4f} {stage1['f1']:>8.4f}")
print(f"{'3. + ISM implicit filter':<35} {stage3['tp']:>7} {stage3['fp']:>7} {stage3['fn']:>7} {stage3['precision']:>8.4f} {stage3['recall']:>8.4f} {stage3['f1']:>8.4f}")
print(f"{'3b. + ISM SuggestFilter':<35} {stage3f['tp']:>7} {stage3f['fp']:>7} {stage3f['fn']:>7} {stage3f['precision']:>8.4f} {stage3f['recall']:>8.4f} {stage3f['f1']:>8.4f}")

# ═══════════════════════════════════════
# PART D: Per-Feature AUC
# ═══════════════════════════════════════
print("\n" + "=" * 70)
print("PART D: Per-Feature AUC (TP vs FP discrimination)")
print("=" * 70)

vcf_features = ['qual', 'gq', 'dp', 'af', 'sb', 'strand_ratio', 'pon_count',
                 'verdict_germline', 'verdict_somatic', 'verdict_subclonal', 'haplotype_flag']
phasing_features = ['lp_phased', 'ps_block_size', 'in_loh', 'in_gain', 'in_loss', 'lp_loh']
ism_features = [c for c in df.columns if c.startswith('ism_') and df[c].dtype in ['float64', 'int64', 'int32', 'float32']]

all_features = vcf_features + phasing_features + ism_features

auc_results = {}
labels_all = df['is_tp'].values

print(f"\n{'Feature':<30} {'AUC':>7} {'TP mean':>10} {'FP mean':>10} {'N':>8} {'Stage':>8}")
print(f"{'-'*73}")

for feat in all_features:
    if feat not in df.columns:
        continue
    vals = pd.to_numeric(df[feat], errors='coerce').values
    valid = ~np.isnan(vals)
    if valid.sum() < 100:
        continue
    v = vals[valid].astype(float)
    l = labels_all[valid]
    if len(set(l)) < 2:
        continue

    auc = compute_auc(l.tolist(), v.tolist())
    tp_mean = v[l == 1].mean() if (l == 1).sum() > 0 else np.nan
    fp_mean = v[l == 0].mean() if (l == 0).sum() > 0 else np.nan

    stage = "S1-VCF" if feat in vcf_features else ("S2-Phase" if feat in phasing_features else "S3-ISM")
    auc_results[feat] = {'auc': round(auc, 4), 'tp_mean': round(float(tp_mean), 4),
                          'fp_mean': round(float(fp_mean), 4), 'n': int(valid.sum()), 'stage': stage}
    print(f"  {feat:<28} {auc:>7.4f} {tp_mean:>10.4f} {fp_mean:>10.4f} {valid.sum():>8d} {stage:>8}")

# ═══════════════════════════════════════
# PART E: Stratified Analysis
# ═══════════════════════════════════════
print("\n" + "=" * 70)
print("PART E: Stratified Analysis")
print("=" * 70)

# E1: TP/FP by CNV type
print("\n[E1] TP/FP by CNV type:")
for cnv_cat in ['neutral', 'gain', 'loss', 'loh', 'none']:
    sub = df[df['cnv_type'] == cnv_cat]
    if len(sub) == 0:
        continue
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {cnv_cat:<10} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f} ({len(sub)} total)")

# E2: TP/FP by LongPhase-TO LOH
print("\n[E2] TP/FP by LongPhase-TO LOH status:")
for loh_val, label in [(1, 'In LOH'), (0, 'Not in LOH')]:
    sub = df[df['lp_loh'] == loh_val]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# E3: TP/FP by phasing
print("\n[E3] TP/FP by LongPhase-TO phasing:")
for ph_val, label in [(1, 'Phased'), (0, 'Unphased')]:
    sub = df[df['lp_phased'] == ph_val]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# E4: TP/FP by Verdict
print("\n[E4] TP/FP by ClairS-TO Verdict:")
for v_name, v_col in [('Somatic', 'verdict_somatic'), ('SubclonalSomatic', 'verdict_subclonal'),
                       ('Germline', 'verdict_germline')]:
    sub = df[df[v_col] == 1]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {v_name:<22} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")
# No verdict
sub = df[(df['verdict_somatic']==0) & (df['verdict_subclonal']==0) & (df['verdict_germline']==0)]
tp = sub['is_tp'].sum()
fp = len(sub) - tp
print(f"  {'No Verdict':<22} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# E5: ISM VerificationClass (ISM subset only)
print("\n[E5] ISM VerificationClass (ISM-processed variants only):")
ism_sub = df[df['has_ism'] == 1]
for vc in ['Strong', 'Weak', 'Noise', 'Subclone', 'Insufficient', '']:
    sub = ism_sub[ism_sub['ism_VerificationClass'] == vc]
    if len(sub) == 0:
        continue
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    label = vc if vc else '(empty)'
    print(f"  {label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# E6: ISM Quality_Tier
print("\n[E6] ISM Quality_Tier:")
for qt in ['High', 'Medium', 'Low', '']:
    sub = ism_sub[ism_sub['ism_Quality_Tier'] == qt]
    if len(sub) == 0:
        continue
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    label = qt if qt else '(empty)'
    print(f"  {label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# E7: ISM Coverage_Category
print("\n[E7] ISM Coverage_Category:")
if 'ism_Coverage_Category' in ism_sub.columns:
    for cc in ['Normal', 'Low', 'Elevated', 'CNV_Loss', 'CNV_Gain', 'High_Copy']:
        sub = ism_sub[ism_sub['ism_Coverage_Category'] == cc]
        if len(sub) == 0:
            continue
        tp = sub['is_tp'].sum()
        fp = len(sub) - tp
        print(f"  {cc:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# E8: ISM Potential_LOH
print("\n[E8] ISM Potential_LOH:")
for loh_val, label in [(1, 'LOH'), (0, 'Not LOH')]:
    sub = ism_sub[ism_sub['ism_Potential_LOH'] == loh_val]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    print(f"  {label:<15} TP={tp:>6} FP={fp:>6} TP_rate={tp/max(tp+fp,1):.4f}")

# ═══════════════════════════════════════
# PART F: ISM Feature AUC on ISM Subset
# ═══════════════════════════════════════
print("\n" + "=" * 70)
print("PART F: ISM Feature AUC (ISM-processed subset only)")
print("=" * 70)

ism_sub = df[df['has_ism'] == 1].copy()
ism_labels = ism_sub['is_tp'].values

print(f"\n  ISM subset: {len(ism_sub)} variants (TP={ism_sub['is_tp'].sum()}, FP={len(ism_sub)-ism_sub['is_tp'].sum()})")
print(f"\n{'Feature':<30} {'AUC (full)':>10} {'AUC (ISM)':>10} {'TP mean':>10} {'FP mean':>10}")
print(f"{'-'*70}")

ism_auc_results = {}
for feat in ism_features + vcf_features + phasing_features:
    if feat not in ism_sub.columns:
        continue
    vals_ism = pd.to_numeric(ism_sub[feat], errors='coerce').values
    valid = ~np.isnan(vals_ism)
    if valid.sum() < 100:
        continue
    v = vals_ism[valid].astype(float)
    l = ism_labels[valid]
    if len(set(l)) < 2:
        continue

    auc_ism = compute_auc(l.tolist(), v.tolist())
    auc_full = auc_results.get(feat, {}).get('auc', np.nan)
    tp_mean = v[l == 1].mean()
    fp_mean = v[l == 0].mean()

    ism_auc_results[feat] = {'auc_full': auc_full, 'auc_ism': round(auc_ism, 4)}
    print(f"  {feat:<28} {auc_full:>10.4f} {auc_ism:>10.4f} {tp_mean:>10.4f} {fp_mean:>10.4f}")

# ═══════════════════════════════════════
# PART G: Save Outputs
# ═══════════════════════════════════════
print("\n" + "=" * 70)
print("PART G: Saving Outputs")
print("=" * 70)

# G1: Full DataFrame
out_csv = os.path.join(OUTDIR, "hcc1395_canonical_multimodal.csv")
df.to_csv(out_csv, index=False)
print(f"  Saved: {out_csv} ({len(df)} rows)")

# G2: AUC results
auc_df = pd.DataFrame([
    {'feature': k, **v} for k, v in auc_results.items()
]).sort_values('auc', ascending=False)
auc_csv = os.path.join(OUTDIR, "hcc1395_canonical_feature_auc.csv")
auc_df.to_csv(auc_csv, index=False)
print(f"  Saved: {auc_csv} ({len(auc_df)} features)")

# G3: Stage metrics JSON
stage_metrics = {
    'sample': 'HCC1395',
    'vcf_source': 'canonical TO pipeline (ONT_5kHz BAM, liaoyoyo2001, Mar 2026)',
    'vcf_path': CLAIRS_TO_VCF,
    'truth_total': TRUTH_TOTAL,
    'total_vcf_variants': total_vcf,
    'total_pass': filter_counts['PASS'],
    'pass_in_bed': counts['PASS_TOTAL'],
    'stages': {
        'stage1_clairs_to_pass': stage1,
        'stage2_longphase_to': stage1,
        'stage3_ism_implicit': stage3,
        'stage3b_ism_suggest_filter': stage3f,
    },
    'ism_implicit_filter': {
        'tp_lost': int(lost_tp),
        'fp_lost': int(lost_fp),
        'total_lost': int(lost_tp + lost_fp),
        'tp_loss_rate': round(lost_tp / counts['TP_COUNT'], 6),
        'fp_loss_rate': round(lost_fp / counts['FP_COUNT'], 6),
    },
    'filter_breakdown': dict(filter_counts),
    'data_provenance': {
        'canonical_archive': ARCHIVE,
        'previous_wrong_vcf': '/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz (zhenyu112, Sep 2025)',
        'correction_note': 'v1 used zhenyu112 VCF (ONT BAM, no MM/ML) + three_way_comparison ISM; v2 uses canonical VCF (ONT_5kHz BAM, with MM/ML) + canonical ISM'
    }
}
stage_json_path = os.path.join(OUTDIR, "hcc1395_canonical_stage_metrics.json")
with open(stage_json_path, 'w') as f:
    json.dump(stage_metrics, f, indent=2, default=str)
print(f"  Saved: {stage_json_path}")

print("\n[DONE] Canonical analysis complete.")
print(f"  Next: Run 05_canonical_plots.py for figures.")
