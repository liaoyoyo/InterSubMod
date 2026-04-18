#!/usr/bin/env python3
"""
Investigation 1: GC-Content Correction of Coverage_Multiple
Using TO pipeline data (ClairS-TO + LongPhase-TO + ISM)
With SEQC2 ground truth CN for validation (HCC1395 only)

Data source:
  ClairS-TO VCF: step01_clairs_to/snv.vcf.gz
  LongPhase-TO phased VCF: step03_longphase_to/tumor_phased.vcf.gz
  ISM TP: step05_intersubmod/intersubmod_tp/significance_summary.csv
  ISM FP: step05_intersubmod/intersubmod_fp/significance_summary.csv
  F1 baseline: 0.7127, filtered: 0.7130
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pysam
from scipy.stats import spearmanr, pearsonr
from scipy.ndimage import gaussian_filter1d
from sklearn.metrics import roc_auc_score
import subprocess
import warnings
warnings.filterwarnings('ignore')

# === Paths ===
TO_BASE = "/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1"
REF_FASTA = "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
SEQC2_GAIN = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_gain_cn.bed"
SEQC2_LOSS = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnv_loss_cn.bed"
FIG_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/figures"
DATA_DIR = "/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data"

# === Load TO ISM data ===
tp = pd.read_csv(f'{TO_BASE}/step05_intersubmod/intersubmod_tp/significance_summary.csv')
fp = pd.read_csv(f'{TO_BASE}/step05_intersubmod/intersubmod_fp/significance_summary.csv')
tp['Label'] = 'TP'
fp['Label'] = 'FP'
all_data = pd.concat([tp, fp], ignore_index=True)
print(f'TO pipeline: TP={len(tp)}, FP={len(fp)}, Total={len(all_data)}')
print(f'Pipeline: ClairS-TO → LongPhase-TO → ISM')
print(f'F1: baseline=0.7127, filtered=0.7130')

# === Step 1: Annotate with SEQC2 CN ===
print('\n=== Annotating with SEQC2 CN ground truth ===')

# Build CN lookup from SEQC2 BED files
cn_intervals = []
for bed_path, cn_type in [(SEQC2_GAIN, 'gain'), (SEQC2_LOSS, 'loss')]:
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            cn = float(parts[3])
            cn_intervals.append((chrom, start, end, cn))

print(f'SEQC2 CNV intervals: {len(cn_intervals)}')

# Annotate each region
def get_seqc2_cn(chrom, pos):
    for c, s, e, cn in cn_intervals:
        if c == chrom and s <= pos < e:
            return cn
    return 2.0  # Default diploid

all_data['SEQC2_CN'] = all_data.apply(lambda r: get_seqc2_cn(r['Chr'], r['Pos']), axis=1)
print(f'SEQC2_CN distribution:')
print(all_data['SEQC2_CN'].value_counts().sort_index().head(15))

# === Step 2: Compute GC content per region ===
print('\n=== Computing GC content per region ===')

ref = pysam.FastaFile(REF_FASTA)
window = 1000  # ±1000bp around SNV

gc_contents = []
for _, row in all_data.iterrows():
    chrom = row['Chr']
    pos = int(row['Pos'])
    start = max(0, pos - window)
    end = pos + window
    try:
        seq = ref.fetch(chrom, start, end).upper()
        gc = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else np.nan
    except:
        gc = np.nan
    gc_contents.append(gc)

all_data['GC_Content'] = gc_contents
ref.close()

valid_gc = all_data['GC_Content'].notna()
print(f'GC computed: {valid_gc.sum()}/{len(all_data)} valid')
print(f'GC range: {all_data["GC_Content"].min():.3f} - {all_data["GC_Content"].max():.3f}')
print(f'GC mean: {all_data["GC_Content"].mean():.3f}')

# === Step 3: GC-Coverage Curve ===
print('\n=== Building GC-Coverage Curve ===')

# Bin GC content and compute median NumReads per bin
gc_bins = np.arange(0.20, 0.75, 0.02)
gc_midpoints = (gc_bins[:-1] + gc_bins[1:]) / 2
gc_medians = []
gc_ns = []

for i in range(len(gc_bins)-1):
    mask = (all_data['GC_Content'] >= gc_bins[i]) & (all_data['GC_Content'] < gc_bins[i+1])
    n = mask.sum()
    gc_ns.append(n)
    if n >= 50:
        gc_medians.append(all_data.loc[mask, 'NumReads'].median())
    else:
        gc_medians.append(np.nan)

gc_medians = np.array(gc_medians, dtype=float)
valid_bins = ~np.isnan(gc_medians)

# Fit polynomial (degree 2) for GC correction
from numpy.polynomial import polynomial as P
gc_fit = P.polyfit(gc_midpoints[valid_bins], gc_medians[valid_bins], deg=2,
                    w=np.sqrt(np.array(gc_ns)[valid_bins]))
gc_expected = P.polyval(all_data['GC_Content'].values, gc_fit)

# GC-corrected NumReads
global_median = all_data['NumReads'].median()
gc_correction_factor = gc_expected / global_median
gc_correction_factor = np.clip(gc_correction_factor, 0.5, 2.0)  # Sanity clamp
all_data['NumReads_GCcorr'] = all_data['NumReads'] / gc_correction_factor

# KDE diploid estimation on GC-corrected NumReads
nr_valid = all_data.loc[all_data['NumReads_GCcorr'] > 5, 'NumReads_GCcorr']
hist, edges = np.histogram(nr_valid, bins=np.arange(0, nr_valid.max()+2, 2))
smoothed = gaussian_filter1d(hist.astype(float), sigma=3)
mode_bin = np.argmax(smoothed)
diploid_gc = (edges[mode_bin] + edges[mode_bin+1]) / 2

# Also original diploid
nr_orig = all_data.loc[all_data['NumReads'] > 5, 'NumReads']
hist_o, edges_o = np.histogram(nr_orig, bins=np.arange(0, nr_orig.max()+2, 2))
smoothed_o = gaussian_filter1d(hist_o.astype(float), sigma=3)
mode_bin_o = np.argmax(smoothed_o)
diploid_orig = (edges_o[mode_bin_o] + edges_o[mode_bin_o+1]) / 2

print(f'Diploid estimate (original): {diploid_orig:.0f}x')
print(f'Diploid estimate (GC-corrected): {diploid_gc:.0f}x')
print(f'SEQC2 ground truth diploid: 44x')

# Compute Coverage_Multiple variants
all_data['CovMult_75x'] = all_data['NumReads'] / 75.0
all_data['CovMult_KDE'] = all_data['NumReads'] / diploid_orig
all_data['CovMult_GCcorr'] = all_data['NumReads_GCcorr'] / diploid_gc

# === Step 4: Validation against SEQC2 CN ===
print('\n=== Pearson r with SEQC2_CN ===')
cn = all_data['SEQC2_CN']

for label, col in [('75x hardcoded', 'CovMult_75x'),
                    ('KDE auto', 'CovMult_KDE'),
                    ('GC-corrected', 'CovMult_GCcorr')]:
    valid = cn.notna() & all_data[col].notna()
    r, p = pearsonr(all_data.loc[valid, col], cn[valid])
    print(f'  {label:<20}: r = {r:.4f}')

# Category accuracy
def ism_category(cm):
    if cm < 0.5: return 'CNV_Loss'
    elif cm < 0.8: return 'Low'
    elif cm <= 1.2: return 'Normal'
    elif cm <= 1.5: return 'Elevated'
    elif cm <= 2.0: return 'CNV_Gain'
    else: return 'High_Copy'

def true_cn_category(c):
    if c <= 0.5: return 'CNV_Loss'
    elif c <= 1.0: return 'Low'
    elif c <= 2.5: return 'Normal'
    elif c <= 3.5: return 'Elevated'
    elif c <= 4.5: return 'CNV_Gain'
    else: return 'High_Copy'

true_cats = cn.apply(true_cn_category)
print(f'\n=== CN Zone Classification Accuracy ===')
for label, col in [('75x', 'CovMult_75x'), ('KDE', 'CovMult_KDE'), ('GC-corr', 'CovMult_GCcorr')]:
    pred = all_data[col].apply(ism_category)
    acc = (pred == true_cats).mean()
    print(f'  {label:<12}: {acc:.1%}')

# Per-CN median Coverage_Multiple
print(f'\n=== Coverage_Multiple Median by True CN ===')
print(f'{"CN":>5} {"n":>6} {"75x":>8} {"KDE":>8} {"GC-corr":>8} {"Expected":>8}')
for c in [1, 2, 3, 4, 5, 6, 8]:
    mask = cn == c
    n = mask.sum()
    if n < 10: continue
    expected = c / 2.0
    print(f'{c:>5} {n:>6} {all_data.loc[mask, "CovMult_75x"].median():>8.3f} '
          f'{all_data.loc[mask, "CovMult_KDE"].median():>8.3f} '
          f'{all_data.loc[mask, "CovMult_GCcorr"].median():>8.3f} '
          f'{expected:>8.2f}')

# === Step 5: TP/FP AUC comparison ===
print(f'\n=== TP/FP Discrimination AUC ===')
labels = (all_data['Label'] == 'TP').astype(int)
for feat_label, col in [('Quality_Score', 'Quality_Score'), ('CovMult_75x', 'CovMult_75x'),
                          ('CovMult_KDE', 'CovMult_KDE'), ('CovMult_GCcorr', 'CovMult_GCcorr'),
                          ('NumReads', 'NumReads')]:
    valid = all_data[col].notna()
    if valid.sum() > 0:
        auc = roc_auc_score(labels[valid], all_data.loc[valid, col])
        print(f'  {feat_label:<20}: AUC = {auc:.4f}')

# GC deviation analysis
print(f'\n=== GC Correction Impact ===')
gc_dev = (all_data['NumReads_GCcorr'] - all_data['NumReads']).abs() / all_data['NumReads']
print(f'Mean relative change from GC correction: {gc_dev.mean():.3%}')
print(f'Max relative change: {gc_dev.max():.3%}')
print(f'Regions changed > 10%: {(gc_dev > 0.1).sum()} ({(gc_dev > 0.1).mean():.1%})')

# === Figure 32: GC-Coverage Curve ===
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# GC-Coverage scatter
ax = axes[0]
sample = all_data.sample(min(5000, len(all_data)), random_state=42)
ax.scatter(sample['GC_Content'], sample['NumReads'], alpha=0.1, s=1, color='steelblue')
gc_plot = np.linspace(0.25, 0.70, 100)
gc_fit_vals = P.polyval(gc_plot, gc_fit)
ax.plot(gc_plot, gc_fit_vals, 'r-', linewidth=2, label='Polynomial fit (deg=2)')
ax.axhline(global_median, color='green', linestyle='--', alpha=0.5, label=f'Global median={global_median:.0f}')
ax.set_xlabel('GC Content')
ax.set_ylabel('NumReads')
ax.set_title('GC-Coverage Relationship', fontweight='bold')
ax.legend(fontsize=9)
ax.set_ylim(0, 250)
ax.grid(True, alpha=0.3)

# GC correction factor
ax = axes[1]
gc_range = np.linspace(0.25, 0.70, 100)
gc_factors = P.polyval(gc_range, gc_fit) / global_median
ax.plot(gc_range, gc_factors, 'b-', linewidth=2)
ax.axhline(1.0, color='green', linestyle='--', alpha=0.5)
ax.set_xlabel('GC Content')
ax.set_ylabel('Correction Factor')
ax.set_title('GC Correction Factor', fontweight='bold')
ax.grid(True, alpha=0.3)
ax.set_ylim(0.5, 1.5)

# Before/after Coverage_Multiple vs SEQC2_CN
ax = axes[2]
valid = cn.notna()
r_before = pearsonr(all_data.loc[valid, 'CovMult_KDE'], cn[valid])[0]
r_after = pearsonr(all_data.loc[valid, 'CovMult_GCcorr'], cn[valid])[0]
sample_v = all_data[valid].sample(min(3000, valid.sum()), random_state=42)
ax.scatter(sample_v['CovMult_KDE'], sample_v['SEQC2_CN'], alpha=0.1, s=1, color='blue', label=f'KDE (r={r_before:.4f})')
ax.scatter(sample_v['CovMult_GCcorr'], sample_v['SEQC2_CN']+0.1, alpha=0.1, s=1, color='red', label=f'GC-corr (r={r_after:.4f})')
ax.plot([0,4],[0,8], 'k--', alpha=0.3)
ax.set_xlabel('Coverage_Multiple')
ax.set_ylabel('SEQC2 True CN')
ax.set_title('CN Estimation Accuracy', fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

fig.suptitle('Investigation 1: GC-Content Correction (TO Pipeline, HCC1395)\n'
             f'ClairS-TO → LongPhase-TO → ISM | TP={len(tp)}, FP={len(fp)}, F1=0.7127',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/32_gc_correction_to.png', dpi=150, bbox_inches='tight')
print(f'\nSaved: 32_gc_correction_to.png')

# === Figure 33: Coverage_Multiple accuracy comparison ===
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
cn_vals = [2, 3, 4, 5, 6]

for idx, (label, col) in enumerate([('75x Hardcoded', 'CovMult_75x'),
                                      ('KDE Auto', 'CovMult_KDE'),
                                      ('GC-Corrected', 'CovMult_GCcorr')]):
    ax = axes[idx]
    data = [all_data.loc[cn==c, col].values for c in cn_vals if (cn==c).sum() >= 10]
    labels_cn = [f'CN={c}\n(n={(cn==c).sum()})' for c in cn_vals if (cn==c).sum() >= 10]

    bp = ax.boxplot(data, labels=labels_cn, patch_artist=True, showfliers=False)
    colors = ['#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # Expected values
    for i, c in enumerate([c for c in cn_vals if (cn==c).sum() >= 10]):
        ax.plot(i+1, c/2.0, 'r*', markersize=12, zorder=5)

    r_val = pearsonr(all_data.loc[cn.notna(), col], cn[cn.notna()])[0]
    acc = (all_data[col].apply(ism_category) == true_cats).mean()
    ax.set_ylabel('Coverage_Multiple')
    ax.set_title(f'{label}\nr={r_val:.4f}, acc={acc:.1%}', fontweight='bold')
    ax.set_ylim(0, 5)
    ax.grid(True, alpha=0.3, axis='y')

fig.suptitle('Coverage_Multiple by True CN: 75x vs KDE vs GC-Corrected\n'
             '(Red stars = expected values | TO Pipeline)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/33_covm_accuracy_comparison_to.png', dpi=150, bbox_inches='tight')
print(f'Saved: 33_covm_accuracy_comparison_to.png')

# Save annotated TO data
out_path = f'{DATA_DIR}/annotated_hcc1395_to_cnv.tsv'
all_data.to_csv(out_path, sep='\t', index=False)
print(f'Saved: {out_path}')

# === SUMMARY ===
r_75 = pearsonr(all_data.loc[cn.notna(), 'CovMult_75x'], cn[cn.notna()])[0]
r_kde = pearsonr(all_data.loc[cn.notna(), 'CovMult_KDE'], cn[cn.notna()])[0]
r_gc = pearsonr(all_data.loc[cn.notna(), 'CovMult_GCcorr'], cn[cn.notna()])[0]
acc_75 = (all_data['CovMult_75x'].apply(ism_category) == true_cats).mean()
acc_kde = (all_data['CovMult_KDE'].apply(ism_category) == true_cats).mean()
acc_gc = (all_data['CovMult_GCcorr'].apply(ism_category) == true_cats).mean()

print(f'\n{"="*70}')
print(f'SUMMARY (TO Pipeline: ClairS-TO → LongPhase-TO → ISM)')
print(f'{"="*70}')
print(f'{"Method":<20} {"Pearson r":>10} {"CN Acc":>8} {"Delta r":>10}')
print(f'{"-"*50}')
print(f'{"75x hardcoded":<20} {r_75:>10.4f} {acc_75:>7.1%}')
print(f'{"KDE auto":<20} {r_kde:>10.4f} {acc_kde:>7.1%} {r_kde-r_75:>+10.4f}')
print(f'{"GC-corrected":<20} {r_gc:>10.4f} {acc_gc:>7.1%} {r_gc-r_75:>+10.4f}')
print(f'\nGo/No-Go threshold: delta-r >= 0.03')
print(f'GC correction delta-r = {r_gc-r_kde:+.4f}')
if abs(r_gc - r_kde) >= 0.03:
    print(f'→ GO: GC correction improves CN estimation')
elif abs(r_gc - r_kde) >= 0.01:
    print(f'→ MARGINAL: Small improvement, may not justify C++ implementation')
else:
    print(f'→ NO-GO: GC correction does not improve CN estimation')
