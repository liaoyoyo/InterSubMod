#!/usr/bin/env python3
"""
Generate visualization plots for TO pipeline multi-stage characterization.
Reads from data/ outputs of 01_multi_stage_characterization.py.
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
os.makedirs(FIGDIR, exist_ok=True)

# Load data
df = pd.read_csv(os.path.join(DATADIR, "hcc1395_pass_multimodal.csv"))
auc_df = pd.read_csv(os.path.join(DATADIR, "hcc1395_feature_auc_by_stage.csv"))
with open(os.path.join(DATADIR, "hcc1395_stage_metrics.json")) as f:
    stages = json.load(f)

# Fix AUC for constant features (all same value → should be 0.5)
for idx, row in auc_df.iterrows():
    feat = row['feature']
    if feat in df.columns:
        vals = df[feat].dropna()
        if vals.nunique() <= 1:
            auc_df.loc[idx, 'auc'] = 0.5
        # Also check if TP mean == FP mean within tolerance
        elif abs(row['tp_mean'] - row['fp_mean']) < 1e-6:
            auc_df.loc[idx, 'auc'] = 0.5

print("Corrected AUC values:")
for _, r in auc_df.sort_values('auc', ascending=False).iterrows():
    if r['auc'] != 0.5:
        print(f"  {r['feature']:<30} AUC={r['auc']:.4f}  stage={r['stage']}")

# ═══════════════════════════════════════
# Figure 1: Stage-by-Stage F1 Pipeline
# ═══════════════════════════════════════
print("\n[Fig 1] Stage-by-stage F1 pipeline...", flush=True)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: TP/FP counts at each stage
stage_names = ['ClairS-TO\nPASS', 'ISM\n(before filter)', 'ISM\n(after filter)']
s = stages['stages']
tps = [s['stage1_clairs_to_pass']['tp'], s['stage3_ism_before_filter']['tp'], s['stage3b_ism_after_filter']['tp']]
fps = [s['stage1_clairs_to_pass']['fp'], s['stage3_ism_before_filter']['fp'], s['stage3b_ism_after_filter']['fp']]
fns = [s['stage1_clairs_to_pass']['fn'], s['stage3_ism_before_filter']['fn'], s['stage3b_ism_after_filter']['fn']]

x = np.arange(len(stage_names))
w = 0.25
axes[0].bar(x - w, tps, w, color='#2ecc71', label='TP')
axes[0].bar(x, fps, w, color='#e74c3c', label='FP')
axes[0].bar(x + w, fns, w, color='#95a5a6', label='FN')
axes[0].set_xticks(x)
axes[0].set_xticklabels(stage_names, fontsize=9)
axes[0].set_ylabel('Count')
axes[0].set_title('A. TP/FP/FN by Pipeline Stage')
axes[0].legend(fontsize=8)
for i, (tp, fp, fn) in enumerate(zip(tps, fps, fns)):
    axes[0].text(i - w, tp + 300, f'{tp}', ha='center', va='bottom', fontsize=7)
    axes[0].text(i, fp + 300, f'{fp}', ha='center', va='bottom', fontsize=7, color='red')

# Panel B: Precision / Recall / F1
precs = [s['stage1_clairs_to_pass']['precision'], s['stage3_ism_before_filter']['precision'], s['stage3b_ism_after_filter']['precision']]
recs = [s['stage1_clairs_to_pass']['recall'], s['stage3_ism_before_filter']['recall'], s['stage3b_ism_after_filter']['recall']]
f1s = [s['stage1_clairs_to_pass']['f1'], s['stage3_ism_before_filter']['f1'], s['stage3b_ism_after_filter']['f1']]

axes[1].plot(x, precs, 'o-', color='#3498db', label='Precision', linewidth=2, markersize=8)
axes[1].plot(x, recs, 's-', color='#e67e22', label='Recall', linewidth=2, markersize=8)
axes[1].plot(x, f1s, 'D-', color='#9b59b6', label='F1', linewidth=2, markersize=8)
axes[1].set_xticks(x)
axes[1].set_xticklabels(stage_names, fontsize=9)
axes[1].set_ylabel('Score')
axes[1].set_ylim(0.5, 0.9)
axes[1].set_title('B. Precision / Recall / F1')
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.3)
for i in range(len(stage_names)):
    axes[1].annotate(f'{f1s[i]:.3f}', (x[i], f1s[i]), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=8, color='#9b59b6')

# Panel C: FILTER breakdown (all 3.18M variants)
fb = stages['filter_breakdown']
categories = ['PASS', 'NonSomatic', 'LowQual+NonSomatic', 'LowQual/Other']
tp_counts = [fb[c]['tp'] for c in categories]
fp_counts = [fb[c]['fp'] for c in categories]

x2 = np.arange(len(categories))
axes[2].barh(x2, fp_counts, 0.4, color='#e74c3c', alpha=0.7, label='FP')
axes[2].barh(x2, tp_counts, 0.4, color='#2ecc71', alpha=0.9, label='TP (truth)')
axes[2].set_yticks(x2)
axes[2].set_yticklabels(categories, fontsize=9)
axes[2].set_xlabel('Count (log scale)')
axes[2].set_xscale('log')
axes[2].set_title('C. All 3.18M Variants by FILTER')
axes[2].legend(fontsize=8)
for i, (tp, fp) in enumerate(zip(tp_counts, fp_counts)):
    if tp > 0:
        axes[2].text(max(tp, fp) * 1.5, i, f'TP={tp}', va='center', fontsize=7)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "01_stage_pipeline_f1.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 01_stage_pipeline_f1.png")

# ═══════════════════════════════════════
# Figure 2: Feature AUC by Stage
# ═══════════════════════════════════════
print("\n[Fig 2] Feature AUC by stage...", flush=True)

# Filter out constant features (AUC=0.5)
auc_real = auc_df[auc_df['auc'] != 0.5].copy()
auc_real['abs_auc'] = abs(auc_real['auc'] - 0.5) + 0.5
auc_real = auc_real.sort_values('abs_auc', ascending=True)

fig, ax = plt.subplots(figsize=(10, max(6, len(auc_real) * 0.35)))
colors = {'S1-VCF': '#3498db', 'S2-Phase': '#e67e22', 'S3-ISM': '#2ecc71'}
bars = ax.barh(range(len(auc_real)), auc_real['auc'].values, color=[colors.get(s, '#95a5a6') for s in auc_real['stage'].values], alpha=0.8)
ax.set_yticks(range(len(auc_real)))
ax.set_yticklabels(auc_real['feature'].values, fontsize=8)
ax.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('AUC (TP=positive)')
ax.set_title('Feature Discriminability: TP vs FP in ClairS-TO PASS Variants\n(Features with AUC=0.50 constant values excluded)')
ax.legend(handles=[Patch(color=c, label=l) for l, c in colors.items()], fontsize=8, loc='lower right')
for i, (auc_val, feat) in enumerate(zip(auc_real['auc'].values, auc_real['feature'].values)):
    ax.text(auc_val + 0.005, i, f'{auc_val:.3f}', va='center', fontsize=7)
ax.set_xlim(0.35, max(auc_real['auc'].values) + 0.06)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "02_feature_auc_by_stage.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 02_feature_auc_by_stage.png")

# ═══════════════════════════════════════
# Figure 3: CNV/LOH Stratification
# ═══════════════════════════════════════
print("\n[Fig 3] CNV/LOH stratification...", flush=True)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: TP rate by CNV type
cnv_types = ['neutral', 'gain', 'loss', 'loh']
tp_rates = []
tp_ns = []
fp_ns = []
for ct in cnv_types:
    sub = df[df['cnv_type'] == ct]
    tp = sub['is_tp'].sum()
    fp = len(sub) - tp
    tp_rates.append(tp / max(tp + fp, 1))
    tp_ns.append(tp)
    fp_ns.append(fp)

bars_a = axes[0].bar(cnv_types, tp_rates, color=['#95a5a6', '#2ecc71', '#e74c3c', '#3498db'])
axes[0].set_ylabel('TP Rate (precision)')
axes[0].set_title('A. TP Rate by CNV/LOH Status')
axes[0].axhline(0.5964, color='red', linestyle='--', alpha=0.5, label=f'Overall={0.5964:.3f}')
axes[0].legend(fontsize=8)
for i, (rate, tp, fp) in enumerate(zip(tp_rates, tp_ns, fp_ns)):
    axes[0].text(i, rate + 0.01, f'{rate:.3f}\n(n={tp+fp})', ha='center', va='bottom', fontsize=8)

# Panel B: TP/FP stacked bar by CNV type
tp_arr = np.array(tp_ns)
fp_arr = np.array(fp_ns)
x3 = np.arange(len(cnv_types))
axes[1].bar(x3, tp_arr, 0.5, color='#2ecc71', label='TP')
axes[1].bar(x3, fp_arr, 0.5, bottom=tp_arr, color='#e74c3c', label='FP')
axes[1].set_xticks(x3)
axes[1].set_xticklabels(cnv_types)
axes[1].set_ylabel('Count')
axes[1].set_title('B. TP/FP Distribution by CNV/LOH')
axes[1].legend(fontsize=8)

# Panel C: AF distribution by TP/FP and CNV status
for ct, color, marker in [('neutral', '#95a5a6', 'o'), ('gain', '#2ecc71', 's'), ('loh', '#3498db', 'D')]:
    sub = df[df['cnv_type'] == ct]
    tp_af = sub[sub['is_tp'] == 1]['af'].values
    fp_af = sub[sub['is_tp'] == 0]['af'].values
    if len(tp_af) > 100 and len(fp_af) > 100:
        # Box positions
        positions_tp = cnv_types.index(ct) * 2
        positions_fp = cnv_types.index(ct) * 2 + 0.7
        bp_tp = axes[2].boxplot([tp_af], positions=[positions_tp], widths=0.5, patch_artist=True,
                                boxprops=dict(facecolor='#2ecc71', alpha=0.5),
                                medianprops=dict(color='black'), showfliers=False)
        bp_fp = axes[2].boxplot([fp_af], positions=[positions_fp], widths=0.5, patch_artist=True,
                                boxprops=dict(facecolor='#e74c3c', alpha=0.5),
                                medianprops=dict(color='black'), showfliers=False)

axes[2].set_xticks([0.35, 2.35, 6.35])
axes[2].set_xticklabels(['neutral', 'gain', 'loh'])
axes[2].set_ylabel('Allele Frequency (AF)')
axes[2].set_title('C. AF Distribution by CNV × TP/FP')
axes[2].legend(handles=[Patch(color='#2ecc71', alpha=0.5, label='TP'),
                        Patch(color='#e74c3c', alpha=0.5, label='FP')], fontsize=8)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "03_cnv_loh_stratification.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 03_cnv_loh_stratification.png")

# ═══════════════════════════════════════
# Figure 4: VCF Feature Distributions (TP vs FP)
# ═══════════════════════════════════════
print("\n[Fig 4] VCF feature distributions...", flush=True)

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
tp = df[df['is_tp'] == 1]
fp = df[df['is_tp'] == 0]

# AF
axes[0, 0].hist(tp['af'].values, bins=50, alpha=0.6, color='#2ecc71', label=f'TP (n={len(tp)})', density=True)
axes[0, 0].hist(fp['af'].values, bins=50, alpha=0.6, color='#e74c3c', label=f'FP (n={len(fp)})', density=True)
axes[0, 0].set_xlabel('AF')
axes[0, 0].set_title('A. Allele Frequency')
axes[0, 0].legend(fontsize=8)

# QUAL
axes[0, 1].hist(tp['qual'].clip(0, 50).values, bins=50, alpha=0.6, color='#2ecc71', density=True)
axes[0, 1].hist(fp['qual'].clip(0, 50).values, bins=50, alpha=0.6, color='#e74c3c', density=True)
axes[0, 1].set_xlabel('QUAL (clipped to 50)')
axes[0, 1].set_title('B. Variant Quality')

# DP
axes[0, 2].hist(tp['dp'].clip(0, 300).values, bins=50, alpha=0.6, color='#2ecc71', density=True)
axes[0, 2].hist(fp['dp'].clip(0, 300).values, bins=50, alpha=0.6, color='#e74c3c', density=True)
axes[0, 2].set_xlabel('DP (clipped to 300)')
axes[0, 2].set_title('C. Read Depth')

# GQ
axes[1, 0].hist(tp['gq'].clip(0, 50).values, bins=50, alpha=0.6, color='#2ecc71', density=True)
axes[1, 0].hist(fp['gq'].clip(0, 50).values, bins=50, alpha=0.6, color='#e74c3c', density=True)
axes[1, 0].set_xlabel('GQ (clipped to 50)')
axes[1, 0].set_title('D. Genotype Quality')

# SB (strand bias)
sb_tp = tp['sb'].values
sb_fp = fp['sb'].values
axes[1, 1].hist(sb_tp[sb_tp < 1], bins=50, alpha=0.6, color='#2ecc71', density=True)
axes[1, 1].hist(sb_fp[sb_fp < 1], bins=50, alpha=0.6, color='#e74c3c', density=True)
axes[1, 1].set_xlabel('SB (p-value, < 1)')
axes[1, 1].set_title('E. Strand Bias')

# Phasing status by TP/FP
phased_rates = {
    'TP\nPhased': len(tp[tp['lp_phased']==1]) / max(len(tp), 1),
    'TP\nUnphased': len(tp[tp['lp_phased']==0]) / max(len(tp), 1),
    'FP\nPhased': len(fp[fp['lp_phased']==1]) / max(len(fp), 1),
    'FP\nUnphased': len(fp[fp['lp_phased']==0]) / max(len(fp), 1),
}
bars = axes[1, 2].bar(phased_rates.keys(), phased_rates.values(),
                       color=['#2ecc71', '#27ae60', '#e74c3c', '#c0392b'])
axes[1, 2].set_ylabel('Fraction')
axes[1, 2].set_title('F. Phasing Status')
for bar, val in zip(bars, phased_rates.values()):
    axes[1, 2].text(bar.get_x() + bar.get_width()/2, val + 0.01, f'{val:.3f}', ha='center', fontsize=8)

plt.suptitle('ClairS-TO PASS Variants: VCF Feature Distributions (HCC1395)', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "04_vcf_feature_distributions.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 04_vcf_feature_distributions.png")

# ═══════════════════════════════════════
# Figure 5: ISM Feature Distributions (TP vs FP)
# ═══════════════════════════════════════
print("\n[Fig 5] ISM feature distributions...", flush=True)

ism_df = df[df['has_ism'] == 1].copy()
ism_tp = ism_df[ism_df['is_tp'] == 1]
ism_fp = ism_df[ism_df['is_tp'] == 0]

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# AlleleDelta
ad_tp = ism_tp['ism_AlleleDelta'].dropna().values
ad_fp = ism_fp['ism_AlleleDelta'].dropna().values
axes[0, 0].hist(ad_tp, bins=50, alpha=0.6, color='#2ecc71', label=f'TP (n={len(ad_tp)})', density=True)
axes[0, 0].hist(ad_fp, bins=50, alpha=0.6, color='#e74c3c', label=f'FP (n={len(ad_fp)})', density=True)
axes[0, 0].set_xlabel('AlleleDelta')
axes[0, 0].set_title('A. ISM AlleleDelta')
axes[0, 0].legend(fontsize=8)

# AlleleSig rate
sig_rates = {
    'TP': ism_tp['ism_AlleleSig'].mean(),
    'FP': ism_fp['ism_AlleleSig'].mean(),
}
bars = axes[0, 1].bar(sig_rates.keys(), sig_rates.values(), color=['#2ecc71', '#e74c3c'])
axes[0, 1].set_ylabel('AlleleSig Rate')
axes[0, 1].set_title('B. AlleleSig Significance Rate')
for bar, val in zip(bars, sig_rates.values()):
    axes[0, 1].text(bar.get_x() + bar.get_width()/2, val + 0.01, f'{val:.3f}', ha='center', fontsize=9)

# Quality_Score
qs_tp = ism_tp['ism_Quality_Score'].dropna().values
qs_fp = ism_fp['ism_Quality_Score'].dropna().values
axes[0, 2].hist(qs_tp[qs_tp > 0], bins=50, alpha=0.6, color='#2ecc71', density=True)
axes[0, 2].hist(qs_fp[qs_fp > 0], bins=50, alpha=0.6, color='#e74c3c', density=True)
axes[0, 2].set_xlabel('Quality_Score')
axes[0, 2].set_title('C. ISM Quality Score')

# NumReads
nr_tp = ism_tp['ism_NumReads'].dropna().clip(0, 500).values
nr_fp = ism_fp['ism_NumReads'].dropna().clip(0, 500).values
axes[1, 0].hist(nr_tp, bins=50, alpha=0.6, color='#2ecc71', density=True)
axes[1, 0].hist(nr_fp, bins=50, alpha=0.6, color='#e74c3c', density=True)
axes[1, 0].set_xlabel('NumReads (clipped to 500)')
axes[1, 0].set_title('D. ISM NumReads')

# CramersV
cv_tp = ism_tp['ism_CramersV'].dropna().values
cv_fp = ism_fp['ism_CramersV'].dropna().values
axes[1, 1].hist(cv_tp[cv_tp > 0], bins=50, alpha=0.6, color='#2ecc71', density=True)
axes[1, 1].hist(cv_fp[cv_fp > 0], bins=50, alpha=0.6, color='#e74c3c', density=True)
axes[1, 1].set_xlabel("Cramer's V (> 0)")
axes[1, 1].set_title("E. ISM Cramer's V")

# VerificationClass breakdown
vc_cats = ['Strong', 'Weak', 'Noise', 'Subclone']
vc_tp_counts = [len(ism_tp[ism_tp['ism_VerificationClass'] == vc]) for vc in vc_cats]
vc_fp_counts = [len(ism_fp[ism_fp['ism_VerificationClass'] == vc]) for vc in vc_cats]
x4 = np.arange(len(vc_cats))
axes[1, 2].bar(x4 - 0.2, vc_tp_counts, 0.35, color='#2ecc71', label='TP')
axes[1, 2].bar(x4 + 0.2, vc_fp_counts, 0.35, color='#e74c3c', label='FP')
axes[1, 2].set_xticks(x4)
axes[1, 2].set_xticklabels(vc_cats)
axes[1, 2].set_ylabel('Count')
axes[1, 2].set_title('F. ISM VerificationClass')
axes[1, 2].legend(fontsize=8)
# Add TP rates
for i, (tpc, fpc) in enumerate(zip(vc_tp_counts, vc_fp_counts)):
    rate = tpc / max(tpc + fpc, 1)
    axes[1, 2].text(i, max(tpc, fpc) + 100, f'TP%={rate:.1%}', ha='center', fontsize=7)

plt.suptitle('ISM Features: TP vs FP (HCC1395, TO Mode, variants with ISM data)', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "05_ism_feature_distributions.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 05_ism_feature_distributions.png")

# ═══════════════════════════════════════
# Figure 6: Multi-Modal Feature Combination
# ═══════════════════════════════════════
print("\n[Fig 6] Multi-modal feature summary...", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Best features per stage
non_trivial = auc_df[(auc_df['auc'] != 0.5) & (auc_df['auc'] < 0.99)].copy()
non_trivial['auc_abs'] = abs(non_trivial['auc'] - 0.5)
top_per_stage = non_trivial.groupby('stage').apply(lambda g: g.nlargest(5, 'auc_abs')).reset_index(drop=True)
top_per_stage = top_per_stage.sort_values(['stage', 'auc_abs'], ascending=[True, True])

colors_stage = {'S1-VCF': '#3498db', 'S2-Phase': '#e67e22', 'S3-ISM': '#2ecc71'}
y_pos = range(len(top_per_stage))
axes[0].barh(y_pos, top_per_stage['auc'].values,
             color=[colors_stage.get(s, '#999') for s in top_per_stage['stage'].values], alpha=0.8)
axes[0].set_yticks(list(y_pos))
axes[0].set_yticklabels(top_per_stage['feature'].values, fontsize=8)
axes[0].axvline(0.5, color='gray', linestyle='--', alpha=0.5)
axes[0].set_xlabel('AUC')
axes[0].set_title('A. Top 5 Features per Stage\n(excluding constant features)')
axes[0].legend(handles=[Patch(color=c, label=l) for l, c in colors_stage.items()], fontsize=8, loc='lower right')
for i, auc_val in enumerate(top_per_stage['auc'].values):
    axes[0].text(auc_val + 0.005, i, f'{auc_val:.3f}', va='center', fontsize=7)

# Panel B: Key insight - neutral vs CNV regions
regions = ['Neutral\n(no CNV)', 'Gain', 'Loss', 'LOH']
tp_rates_cnv = []
for ct in ['neutral', 'gain', 'loss', 'loh']:
    sub = df[df['cnv_type'] == ct]
    tp_rates_cnv.append(sub['is_tp'].mean())

bar_colors = ['#bdc3c7', '#2ecc71', '#e74c3c', '#3498db']
bars = axes[1].bar(regions, tp_rates_cnv, color=bar_colors, edgecolor='black', linewidth=0.5)
axes[1].axhline(df['is_tp'].mean(), color='red', linestyle='--', alpha=0.5, label=f'Overall TP rate={df["is_tp"].mean():.3f}')
axes[1].set_ylabel('TP Rate')
axes[1].set_title('B. CNV Status: Strongest Single Predictor\n(HCC1395 PASS variants)')
axes[1].legend(fontsize=8)
for bar, rate in zip(bars, tp_rates_cnv):
    axes[1].text(bar.get_x() + bar.get_width()/2, rate + 0.01, f'{rate:.3f}', ha='center', fontsize=10, fontweight='bold')

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "06_multimodal_summary.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 06_multimodal_summary.png")

# ═══════════════════════════════════════
# Figure 7: AF × CNV × TP/FP Scatter
# ═══════════════════════════════════════
print("\n[Fig 7] AF × CNV scatter...", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: AF vs DP colored by TP/FP
sample_n = min(10000, len(df))
sample_idx = np.random.RandomState(42).choice(len(df), sample_n, replace=False)
df_s = df.iloc[sample_idx]

tp_s = df_s[df_s['is_tp'] == 1]
fp_s = df_s[df_s['is_tp'] == 0]

axes[0].scatter(fp_s['af'], fp_s['dp'].clip(0, 300), alpha=0.15, s=3, color='#e74c3c', label='FP')
axes[0].scatter(tp_s['af'], tp_s['dp'].clip(0, 300), alpha=0.15, s=3, color='#2ecc71', label='TP')
axes[0].set_xlabel('AF')
axes[0].set_ylabel('DP (clipped to 300)')
axes[0].set_title('A. AF vs DP (10K sample)')
axes[0].legend(fontsize=8, markerscale=5)

# Panel B: AF distribution by CNV type × TP/FP
data_bp = []
labels_bp = []
positions_bp = []
colors_bp = []
pos = 0
for ct in ['neutral', 'gain', 'loh']:
    sub = df[df['cnv_type'] == ct]
    tp_af_sub = sub[sub['is_tp'] == 1]['af'].values
    fp_af_sub = sub[sub['is_tp'] == 0]['af'].values
    if len(tp_af_sub) > 50:
        data_bp.append(tp_af_sub)
        positions_bp.append(pos)
        colors_bp.append('#2ecc71')
        labels_bp.append(f'{ct}\nTP')
        pos += 1
    if len(fp_af_sub) > 50:
        data_bp.append(fp_af_sub)
        positions_bp.append(pos)
        colors_bp.append('#e74c3c')
        labels_bp.append(f'{ct}\nFP')
        pos += 1
    pos += 0.5

bp = axes[1].boxplot(data_bp, positions=positions_bp, widths=0.6, patch_artist=True, showfliers=False)
for patch, c in zip(bp['boxes'], colors_bp):
    patch.set_facecolor(c)
    patch.set_alpha(0.6)
axes[1].set_xticks(positions_bp)
axes[1].set_xticklabels(labels_bp, fontsize=8)
axes[1].set_ylabel('AF')
axes[1].set_title('B. AF by CNV Type × TP/FP')

plt.suptitle('HCC1395 ClairS-TO PASS: AF and Depth Patterns', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "07_af_cnv_scatter.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 07_af_cnv_scatter.png")

print("\n[DONE] All 7 figures generated.")
