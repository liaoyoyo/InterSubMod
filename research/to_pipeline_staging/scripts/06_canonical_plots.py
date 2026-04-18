#!/usr/bin/env python3
"""
Canonical TO Pipeline Visualization (v2)
=========================================
Generate figures from canonical HCC1395 TO pipeline analysis.
Uses data from 04_canonical_multi_stage_analysis.py output.
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
os.makedirs(FIGDIR, exist_ok=True)

# Load data
df = pd.read_csv(os.path.join(DATADIR, "hcc1395_canonical_multimodal.csv"))
auc_df = pd.read_csv(os.path.join(DATADIR, "hcc1395_canonical_feature_auc.csv"))
with open(os.path.join(DATADIR, "hcc1395_canonical_stage_metrics.json")) as f:
    metrics = json.load(f)

print(f"Loaded: {len(df)} variants, {len(auc_df)} features")

tp = df[df['is_tp']==1]
fp = df[df['is_tp']==0]

# ═══════════════════════════════════════
# Figure C1: Pipeline Waterfall (Canonical)
# ═══════════════════════════════════════
print("[Fig C1] Pipeline waterfall...", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Panel A: Full FILTER breakdown
stages = metrics['stages']
s1 = stages['stage1_clairs_to_pass']
s3 = stages['stage3_ism_implicit']
s3b = stages['stage3b_ism_suggest_filter']

waterfall_data = [
    ('ClairS-TO\nAll 3.2M', metrics['total_vcf_variants'], 0, '#bdc3c7'),
    ('NonSomatic\nfilter', metrics['filter_breakdown']['NonSomatic'] + metrics['filter_breakdown'].get('LowQual+NonSomatic', 0), 0, '#e74c3c'),
    ('LowQual\nfilters', sum(v for k, v in metrics['filter_breakdown'].items() if k not in ['PASS', 'NonSomatic', 'LowQual+NonSomatic']), 0, '#e67e22'),
    ('PASS', metrics['total_pass'], 0, '#3498db'),
    ('In truth BED', s1['tp'] + s1['fp'], 0, '#2ecc71'),
]

ax = axes[0]
labels = [d[0] for d in waterfall_data]
values = [d[1] for d in waterfall_data]
colors = [d[3] for d in waterfall_data]
bars = ax.bar(range(len(values)), values, color=colors, edgecolor='black', linewidth=0.5)
for i, (label, val) in enumerate(zip(labels, values)):
    if val > 100000:
        ax.text(i, val + max(values)*0.02, f'{val:,}', ha='center', fontsize=8, fontweight='bold')
    else:
        ax.text(i, val + max(values)*0.02, f'{val:,}', ha='center', fontsize=8)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel('Variant Count')
ax.set_yscale('log')
ax.set_ylim(1e3, 1e7)
ax.set_title('A. ClairS-TO Variant Filtering (Log Scale)')

# Panel B: TP/FP pipeline stages
ax = axes[1]
stage_labels = ['ClairS-TO\nPASS\n(in BED)', 'LongPhase-TO\n(unchanged)', 'ISM\nimplicit', 'ISM\nSuggestFilter']
tp_vals = [s1['tp'], s1['tp'], s3['tp'], s3b['tp']]
fp_vals = [s1['fp'], s1['fp'], s3['fp'], s3b['fp']]
f1_vals = [s1['f1'], s1['f1'], s3['f1'], s3b['f1']]

x = np.arange(len(stage_labels))
w = 0.35
bars_tp = ax.bar(x - w/2, tp_vals, w, color='#2ecc71', alpha=0.8, label='TP', edgecolor='black', linewidth=0.3)
bars_fp = ax.bar(x + w/2, fp_vals, w, color='#e74c3c', alpha=0.8, label='FP', edgecolor='black', linewidth=0.3)

ax2 = ax.twinx()
ax2.plot(x, f1_vals, 'ko-', markersize=8, linewidth=2, label='F1')
for i, f1 in enumerate(f1_vals):
    ax2.annotate(f'F1={f1:.4f}', xy=(i, f1), xytext=(0, 10),
                 textcoords='offset points', ha='center', fontsize=8, fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(stage_labels, fontsize=8)
ax.set_ylabel('Count')
ax2.set_ylabel('F1 Score')
ax2.set_ylim(0.70, 0.72)
ax.legend(loc='upper left', fontsize=8)
ax2.legend(loc='upper right', fontsize=8)
ax.set_title('B. TP/FP at Each Pipeline Stage (Canonical)')

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "C01_canonical_pipeline_waterfall.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: C01_canonical_pipeline_waterfall.png")

# ═══════════════════════════════════════
# Figure C2: Feature AUC Ranking
# ═══════════════════════════════════════
print("[Fig C2] Feature AUC ranking...", flush=True)

fig, ax = plt.subplots(figsize=(10, 10))

# Sort by AUC (use max(auc, 1-auc) for inverted features)
auc_df_plot = auc_df.copy()
auc_df_plot['auc_effective'] = auc_df_plot['auc'].apply(lambda x: max(x, 1-x))
auc_df_plot = auc_df_plot.sort_values('auc_effective', ascending=True)

# Filter out features with AUC=0.5 (no discrimination)
auc_plot = auc_df_plot[auc_df_plot['auc_effective'] > 0.51]

colors_map = {'S1-VCF': '#3498db', 'S2-Phase': '#e67e22', 'S3-ISM': '#2ecc71'}
colors = [colors_map.get(s, '#bdc3c7') for s in auc_plot['stage']]
y = range(len(auc_plot))
bars = ax.barh(y, auc_plot['auc_effective'].values, color=colors, alpha=0.8, edgecolor='black', linewidth=0.3)

ax.set_yticks(y)
ax.set_yticklabels(auc_plot['feature'].values, fontsize=8)
ax.axvline(0.5, color='gray', linestyle='--', alpha=0.5, label='Random (0.5)')
ax.axvline(0.6, color='orange', linestyle='--', alpha=0.3, label='Weak (0.6)')
ax.set_xlabel('AUC (effective: max(AUC, 1-AUC))')
ax.set_title('Feature AUC for TP vs FP Discrimination\n(HCC1395 Canonical TO, 40,239 PASS variants)')
ax.set_xlim(0.5, 1.0)

# Add value labels
for i, (auc, feat) in enumerate(zip(auc_plot['auc_effective'].values, auc_plot['feature'].values)):
    ax.text(auc + 0.005, i, f'{auc:.3f}', va='center', fontsize=7)

# Legend for stages
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#3498db', label='VCF (Stage 1)'),
                   Patch(facecolor='#e67e22', label='Phase/CNV (Stage 2)'),
                   Patch(facecolor='#2ecc71', label='ISM (Stage 3)')]
ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "C02_canonical_feature_auc_ranking.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: C02_canonical_feature_auc_ranking.png")

# ═══════════════════════════════════════
# Figure C3: VCF Feature Distributions
# ═══════════════════════════════════════
print("[Fig C3] VCF feature distributions...", flush=True)

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# AF distribution
ax = axes[0, 0]
ax.hist(tp['af'], bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={tp["af"].mean():.3f})')
ax.hist(fp['af'], bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={fp["af"].mean():.3f})')
ax.set_xlabel('AF (Allele Frequency)')
ax.set_title('A. AF Distribution (FP > TP)')
ax.legend(fontsize=8)

# DP distribution
ax = axes[0, 1]
ax.hist(tp['dp'].clip(0, 300), bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={tp["dp"].mean():.0f})')
ax.hist(fp['dp'].clip(0, 300), bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={fp["dp"].mean():.0f})')
ax.set_xlabel('DP (Depth)')
ax.set_title('B. Depth Distribution (nearly identical)')
ax.legend(fontsize=8)

# QUAL distribution
ax = axes[0, 2]
ax.hist(tp['qual'], bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={tp["qual"].mean():.1f})')
ax.hist(fp['qual'], bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={fp["qual"].mean():.1f})')
ax.set_xlabel('QUAL')
ax.set_title('C. QUAL Distribution (nearly identical)')
ax.legend(fontsize=8)

# GQ distribution
ax = axes[1, 0]
ax.hist(tp['gq'], bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={tp["gq"].mean():.1f})')
ax.hist(fp['gq'], bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={fp["gq"].mean():.1f})')
ax.set_xlabel('GQ (Genotype Quality)')
ax.set_title('D. GQ Distribution')
ax.legend(fontsize=8)

# SB distribution
ax = axes[1, 1]
ax.hist(tp['sb'].clip(0, 2), bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={tp["sb"].mean():.3f})')
ax.hist(fp['sb'].clip(0, 2), bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={fp["sb"].mean():.3f})')
ax.set_xlabel('SB (Strand Bias)')
ax.set_title('E. Strand Bias (nearly identical)')
ax.legend(fontsize=8)

# Phasing status
ax = axes[1, 2]
tp_phased_rate = tp['lp_phased'].mean()
fp_phased_rate = fp['lp_phased'].mean()
x_ph = [0, 1]
ax.bar([0-0.15, 1-0.15], [1-tp_phased_rate, tp_phased_rate], 0.3, color='#2ecc71', label='TP')
ax.bar([0+0.15, 1+0.15], [1-fp_phased_rate, fp_phased_rate], 0.3, color='#e74c3c', label='FP')
ax.set_xticks([0, 1])
ax.set_xticklabels(['Unphased', 'Phased'])
ax.set_ylabel('Rate')
ax.set_title(f'F. LongPhase-TO Phasing\n(TP={100*tp_phased_rate:.1f}%, FP={100*fp_phased_rate:.1f}%)')
ax.legend(fontsize=8)

plt.suptitle('VCF Feature Distributions: TP vs FP (HCC1395 Canonical TO)', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "C03_canonical_vcf_features.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: C03_canonical_vcf_features.png")

# ═══════════════════════════════════════
# Figure C4: CNV/LOH Stratification
# ═══════════════════════════════════════
print("[Fig C4] CNV/LOH stratification...", flush=True)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: TP rate by CNV type
ax = axes[0]
cnv_cats = ['gain', 'loh', 'loss', 'none']
tp_rates = []
counts_cnv = []
for cat in cnv_cats:
    sub = df[df['cnv_type'] == cat]
    if len(sub) > 0:
        tp_rates.append(sub['is_tp'].mean())
        counts_cnv.append(len(sub))
    else:
        tp_rates.append(0)
        counts_cnv.append(0)

colors_cnv = ['#2ecc71', '#3498db', '#e74c3c', '#bdc3c7']
bars = ax.bar(cnv_cats, tp_rates, color=colors_cnv, edgecolor='black', linewidth=0.5)
ax.axhline(df['is_tp'].mean(), color='red', linestyle='--', alpha=0.5, label=f'Overall={df["is_tp"].mean():.3f}')
for bar, rate, n in zip(bars, tp_rates, counts_cnv):
    ax.text(bar.get_x()+bar.get_width()/2, rate+0.01, f'{rate:.3f}\n(n={n:,})', ha='center', fontsize=8)
ax.set_ylabel('TP Rate')
ax.set_title('A. TP Rate by CNV Type')
ax.legend(fontsize=8)

# Panel B: AF by CNV × TP/FP
ax = axes[1]
cnv_cats_b = ['gain', 'loh', 'none']
x = np.arange(len(cnv_cats_b))
tp_af = [tp[tp['cnv_type']==c]['af'].mean() for c in cnv_cats_b]
fp_af = [fp[fp['cnv_type']==c]['af'].mean() for c in cnv_cats_b]
ax.bar(x-0.15, tp_af, 0.3, color='#2ecc71', label='TP')
ax.bar(x+0.15, fp_af, 0.3, color='#e74c3c', label='FP')
ax.set_xticks(x)
ax.set_xticklabels(cnv_cats_b)
ax.set_ylabel('Mean AF')
ax.set_title('B. Mean AF by CNV × TP/FP')
ax.legend(fontsize=8)

# Panel C: LongPhase LOH TP rate
ax = axes[2]
loh_cats = ['In LOH', 'Not in LOH']
loh_tp_rates = []
for loh_val in [1, 0]:
    sub = df[df['lp_loh'] == loh_val]
    loh_tp_rates.append(sub['is_tp'].mean())
bars = ax.bar(loh_cats, loh_tp_rates, color=['#3498db', '#bdc3c7'], edgecolor='black', linewidth=0.5)
ax.axhline(df['is_tp'].mean(), color='red', linestyle='--', alpha=0.5, label=f'Overall={df["is_tp"].mean():.3f}')
for bar, rate in zip(bars, loh_tp_rates):
    ax.text(bar.get_x()+bar.get_width()/2, rate+0.01, f'{rate:.3f}', ha='center', fontsize=9, fontweight='bold')
ax.set_ylabel('TP Rate')
ax.set_title('C. TP Rate by LongPhase-TO LOH')
ax.legend(fontsize=8)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "C04_canonical_cnv_stratification.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: C04_canonical_cnv_stratification.png")

# ═══════════════════════════════════════
# Figure C5: ISM Features (Canonical)
# ═══════════════════════════════════════
print("[Fig C5] ISM feature distributions...", flush=True)

ism_sub = df[df['has_ism']==1].copy()
ism_tp = ism_sub[ism_sub['is_tp']==1]
ism_fp = ism_sub[ism_sub['is_tp']==0]

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Panel A: VerificationClass
ax = axes[0, 0]
vc_cats = ['Strong', 'Weak', 'Noise', 'Subclone']
vc_tp_rates = []
vc_counts = []
for vc in vc_cats:
    sub = ism_sub[ism_sub['ism_VerificationClass'] == vc]
    if len(sub) > 0:
        vc_tp_rates.append(sub['is_tp'].mean())
        vc_counts.append(len(sub))
    else:
        vc_tp_rates.append(0)
        vc_counts.append(0)
vc_colors = ['#27ae60', '#f39c12', '#e74c3c', '#9b59b6']
bars = ax.bar(vc_cats, vc_tp_rates, color=vc_colors, edgecolor='black', linewidth=0.5)
for bar, rate, n in zip(bars, vc_tp_rates, vc_counts):
    ax.text(bar.get_x()+bar.get_width()/2, rate+0.01, f'{rate:.3f}\n(n={n:,})', ha='center', fontsize=7)
ax.axhline(ism_sub['is_tp'].mean(), color='red', linestyle='--', alpha=0.5)
ax.set_ylabel('TP Rate')
ax.set_title('A. TP Rate by VerificationClass')

# Panel B: Quality_Score
ax = axes[0, 1]
qs_tp = ism_tp['ism_Quality_Score'].dropna()
qs_fp = ism_fp['ism_Quality_Score'].dropna()
ax.hist(qs_tp, bins=30, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={qs_tp.mean():.1f})')
ax.hist(qs_fp, bins=30, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={qs_fp.mean():.1f})')
ax.set_xlabel('Quality_Score')
ax.set_title('B. Quality_Score Distribution')
ax.legend(fontsize=8)

# Panel C: HPFineNGroups
ax = axes[0, 2]
hpng_tp = ism_tp['ism_HPFineNGroups'].dropna()
hpng_fp = ism_fp['ism_HPFineNGroups'].dropna()
bins_hpng = np.arange(0, 10)
ax.hist(hpng_tp.clip(0, 9), bins=bins_hpng, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={hpng_tp.mean():.2f})')
ax.hist(hpng_fp.clip(0, 9), bins=bins_hpng, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={hpng_fp.mean():.2f})')
ax.set_xlabel('HPFineNGroups')
ax.set_title('C. HPFineNGroups (AUC=0.826)')
ax.legend(fontsize=8)

# Panel D: AlleleDelta
ax = axes[1, 0]
ad_tp = ism_tp['ism_AlleleDelta'].dropna()
ad_fp = ism_fp['ism_AlleleDelta'].dropna()
ax.hist(ad_tp, bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={ad_tp.mean():.4f})')
ax.hist(ad_fp, bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={ad_fp.mean():.4f})')
ax.set_xlabel('AlleleDelta')
ax.set_xlim(0, 0.3)
ax.set_title('D. AlleleDelta (AUC=0.569)')
ax.legend(fontsize=8)

# Panel E: Coverage_Category TP rate
ax = axes[1, 1]
if 'ism_Coverage_Category' in ism_sub.columns:
    cc_cats = ['Normal', 'Low', 'Elevated', 'CNV_Loss', 'CNV_Gain', 'High_Copy']
    cc_rates = []
    cc_ns = []
    for cc in cc_cats:
        sub = ism_sub[ism_sub['ism_Coverage_Category'] == cc]
        if len(sub) > 0:
            cc_rates.append(sub['is_tp'].mean())
            cc_ns.append(len(sub))
        else:
            cc_rates.append(0)
            cc_ns.append(0)
    cc_colors = ['#3498db', '#95a5a6', '#e67e22', '#e74c3c', '#2ecc71', '#9b59b6']
    bars = ax.bar(range(len(cc_cats)), cc_rates, color=cc_colors, edgecolor='black', linewidth=0.3)
    ax.set_xticks(range(len(cc_cats)))
    ax.set_xticklabels(cc_cats, fontsize=7, rotation=30)
    ax.axhline(ism_sub['is_tp'].mean(), color='red', linestyle='--', alpha=0.5)
    for bar, rate, n in zip(bars, cc_rates, cc_ns):
        ax.text(bar.get_x()+bar.get_width()/2, rate+0.005, f'{rate:.3f}\n({n:,})', ha='center', fontsize=6)
    ax.set_ylabel('TP Rate')
    ax.set_title('E. TP Rate by Coverage_Category')

# Panel F: HP_Ratio
ax = axes[1, 2]
hr_tp = ism_tp['ism_HP_Ratio'].dropna()
hr_fp = ism_fp['ism_HP_Ratio'].dropna()
ax.hist(hr_tp, bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={hr_tp.mean():.3f})')
ax.hist(hr_fp, bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={hr_fp.mean():.3f})')
ax.set_xlabel('HP_Ratio')
ax.set_title('F. HP_Ratio (AUC=0.602)')
ax.legend(fontsize=8)

plt.suptitle('ISM Feature Analysis (HCC1395 Canonical TO, 40,213 ISM-processed)', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "C05_canonical_ism_features.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: C05_canonical_ism_features.png")

# ═══════════════════════════════════════
# Figure C6: ISM Implicit Filter Detail
# ═══════════════════════════════════════
print("[Fig C6] ISM implicit filter detail...", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: Implicit filter comparison
ax = axes[0]
implicit = metrics['ism_implicit_filter']
categories = ['TP lost', 'FP lost']
values = [implicit['tp_lost'], implicit['fp_lost']]
colors_if = ['#2ecc71', '#e74c3c']
bars = ax.bar(categories, values, color=colors_if, edgecolor='black', linewidth=0.5)
for bar, val in zip(bars, values):
    ax.text(bar.get_x()+bar.get_width()/2, val+0.5, str(val), ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('Count')
ax.set_title(f'A. ISM Implicit Filter: Only {sum(values)} variants removed\n(TP loss={implicit["tp_loss_rate"]*100:.3f}%, FP loss={implicit["fp_loss_rate"]*100:.3f}%)')

# Panel B: SuggestFilter analysis
ax = axes[1]
sf_tp = df[(df['is_tp']==1) & (df.get('ism_SuggestFilter', 0)==1)].shape[0] if 'ism_SuggestFilter' in df.columns else 0
sf_fp = df[(df['is_tp']==0) & (df.get('ism_SuggestFilter', 0)==1)].shape[0] if 'ism_SuggestFilter' in df.columns else 0
categories_sf = ['TP removed', 'FP removed']
values_sf = [sf_tp, sf_fp]
bars = ax.bar(categories_sf, values_sf, color=colors_if, edgecolor='black', linewidth=0.5)
for bar, val in zip(bars, values_sf):
    ax.text(bar.get_x()+bar.get_width()/2, val+1, str(val), ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('Count')
f1_before = metrics['stages']['stage3_ism_implicit']['f1']
f1_after = metrics['stages']['stage3b_ism_suggest_filter']['f1']
ax.set_title(f'B. ISM SuggestFilter: {sf_tp+sf_fp} variants removed\n(F1: {f1_before:.4f} → {f1_after:.4f}, Δ={f1_after-f1_before:+.4f})')

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "C06_canonical_ism_filter_detail.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: C06_canonical_ism_filter_detail.png")

print("\n[DONE] All canonical plots generated.")
