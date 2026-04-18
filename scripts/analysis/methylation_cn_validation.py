#!/usr/bin/env python3
"""
Investigation 2A + 2B + 3A: Methylation Feature-CN Correlation & Validation
Using SEQC2 ground truth CN for HCC1395 (31,640 regions)

Questions:
  2A: Do ISM methylation features correlate with true CN? (residualized against NumReads)
  2B: Can HPFineNGroups predict CN >= 4? (sensitivity + specificity)
  3A: Does clustering quality degrade in high-CN regions?
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr, pearsonr, wilcoxon, mannwhitneyu
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

# === Data ===
DATA_PATH = '/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv'
FIG_DIR = '/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/figures'
DATA_DIR = '/big7_disk/liaoyoyo2001/InterSubMod/research/seqc2_cnv_stratification/data'

df = pd.read_csv(DATA_PATH, sep='\t')
print(f'Loaded {len(df)} regions, TP={sum(df["Label"]=="TP")}, FP={sum(df["Label"]=="FP")}')

# Work with TP only for CN-methylation analysis (avoid FP confound)
tp = df[df['Label'] == 'TP'].copy()
# Also keep full for reference
# Filter to regions with valid SEQC2_CN
tp = tp[tp['SEQC2_CN'].notna() & (tp['SEQC2_CN'] > 0)].copy()
print(f'TP with valid SEQC2_CN: {len(tp)}')
print(f'SEQC2_CN distribution:\n{tp["SEQC2_CN"].value_counts().sort_index()}')

# ============================================================
# Investigation 2A: Feature-CN Correlation Matrix
# ============================================================
print('\n' + '='*70)
print('Investigation 2A: ISM Feature vs SEQC2_CN Correlation')
print('='*70)

features = [
    ('PairwiseMeanDist', 'HP-free'),
    ('PairwiseMedianDist', 'HP-free'),
    ('ClusterPermanovaF', 'HP-free'),
    ('Stability', 'HP-free'),
    ('CramersV', 'HP-free'),
    ('HPFineNGroups', 'HP-dependent'),
    ('HPMergedDelta', 'HP-dependent'),
    ('HPFineF', 'HP-dependent'),
    ('AlleleDelta', 'HP-dependent'),
    ('LabelHPPermanovaF', 'HP-dependent'),
    ('LabelAllelePermanovaF', 'HP-dependent'),
    ('HP_Ratio', 'HP-dependent'),
    ('HP1FamilyN', 'Structural'),
    ('HP2FamilyN', 'Structural'),
    ('NumReads', 'Structural'),
    ('NumCpGs', 'Structural'),
    ('Coverage_Multiple', 'Coverage'),
    ('Quality_Score', 'Composite'),
]

results_2a = []
cn = tp['SEQC2_CN'].values
nr = tp['NumReads'].values.reshape(-1, 1)

print(f'\n{"Feature":<30} {"Type":<15} {"Raw r":>8} {"Raw p":>10} {"Resid r":>8} {"Resid p":>10} {"Verdict":>10}')
print('-' * 100)

for feat_name, feat_type in features:
    if feat_name not in tp.columns:
        continue
    vals = tp[feat_name].values.astype(float)
    valid = ~np.isnan(vals) & ~np.isnan(cn)
    if valid.sum() < 100:
        continue

    # Raw correlation
    r_raw, p_raw = spearmanr(vals[valid], cn[valid])

    # Residualized: regress feature on NumReads, use residuals
    lr = LinearRegression()
    lr.fit(nr[valid], vals[valid])
    residuals = vals[valid] - lr.predict(nr[valid])
    r_resid, p_resid = spearmanr(residuals, cn[valid])

    verdict = ''
    if abs(r_resid) >= 0.25:
        verdict = 'SIGNAL'
    elif abs(r_resid) >= 0.15:
        verdict = 'weak'
    else:
        verdict = 'confound'

    print(f'{feat_name:<30} {feat_type:<15} {r_raw:>+8.3f} {p_raw:>10.2e} {r_resid:>+8.3f} {p_resid:>10.2e} {verdict:>10}')
    results_2a.append({
        'Feature': feat_name, 'Type': feat_type,
        'Raw_Spearman': r_raw, 'Raw_p': p_raw,
        'Residualized_Spearman': r_resid, 'Residualized_p': p_resid,
        'n': valid.sum(), 'Verdict': verdict
    })

# Save results
results_2a_df = pd.DataFrame(results_2a)
results_2a_df.to_csv(f'{DATA_DIR}/methylation_cn_correlation.tsv', sep='\t', index=False)

# === Figure 28: Feature-CN Correlation Heatmap ===
fig, ax = plt.subplots(figsize=(10, 8))
plot_df = results_2a_df.sort_values('Raw_Spearman', ascending=False)
x = np.arange(len(plot_df))
width = 0.35

bars1 = ax.barh(x + width/2, plot_df['Raw_Spearman'], width, label='Raw Spearman', color='steelblue', alpha=0.7)
bars2 = ax.barh(x - width/2, plot_df['Residualized_Spearman'], width, label='Residualized (NumReads removed)', color='coral', alpha=0.7)

ax.set_yticks(x)
ax.set_yticklabels([f"{r['Feature']} ({r['Type']})" for _, r in plot_df.iterrows()], fontsize=9)
ax.axvline(0.25, color='green', linestyle='--', alpha=0.5, label='|r|=0.25 threshold')
ax.axvline(-0.25, color='green', linestyle='--', alpha=0.5)
ax.axvline(0.15, color='orange', linestyle=':', alpha=0.5, label='|r|=0.15 threshold')
ax.axvline(-0.15, color='orange', linestyle=':', alpha=0.5)
ax.set_xlabel('Spearman Correlation with SEQC2_CN')
ax.set_title('Investigation 2A: ISM Feature vs True CN Correlation\n(HCC1395 TP only, n={})'.format(len(tp)), fontweight='bold')
ax.legend(loc='lower right', fontsize=9)
ax.grid(True, alpha=0.3, axis='x')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/28_feature_cn_correlation_heatmap.png', dpi=150, bbox_inches='tight')
print(f'\nSaved: 28_feature_cn_correlation_heatmap.png')


# ============================================================
# Investigation 2B: Cluster Count vs CN Agreement
# ============================================================
print('\n' + '='*70)
print('Investigation 2B: HPFineNGroups vs SEQC2_CN Agreement')
print('='*70)

# 2B.1: HPFineNGroups distribution by CN
print(f'\n{"CN":>5} {"n":>6} {"NGroups_mean":>12} {"NGroups=2":>10} {"NGroups=3":>10} {"NGroups=4":>10}')
print('-' * 60)
cn_groups_data = []
for c in sorted(tp['SEQC2_CN'].unique()):
    mask = tp['SEQC2_CN'] == c
    n = mask.sum()
    if n < 10: continue
    ng = tp.loc[mask, 'HPFineNGroups']
    ng_valid = ng.dropna()
    if len(ng_valid) == 0: continue
    print(f'{c:>5.1f} {n:>6} {ng_valid.mean():>12.2f} {(ng_valid==2).mean():>10.1%} {(ng_valid==3).mean():>10.1%} {(ng_valid==4).mean():>10.1%}')
    cn_groups_data.append({'CN': c, 'n': n, 'mean_NGroups': ng_valid.mean(),
                           'pct_2': (ng_valid==2).mean(), 'pct_3': (ng_valid==3).mean(), 'pct_4': (ng_valid==4).mean()})

# 2B.2: Residualized HPFineNGroups vs CN
tp_valid = tp.dropna(subset=['HPFineNGroups', 'SEQC2_CN']).copy()
lr = LinearRegression()
lr.fit(tp_valid[['NumReads']], tp_valid['HPFineNGroups'])
tp_valid['HPFineNGroups_resid'] = tp_valid['HPFineNGroups'] - lr.predict(tp_valid[['NumReads']])

r_raw, _ = spearmanr(tp_valid['HPFineNGroups'], tp_valid['SEQC2_CN'])
r_resid, _ = spearmanr(tp_valid['HPFineNGroups_resid'], tp_valid['SEQC2_CN'])
print(f'\nHPFineNGroups vs SEQC2_CN:')
print(f'  Raw Spearman:         r = {r_raw:.3f}')
print(f'  Residualized (NR):    r = {r_resid:.3f}')
print(f'  Confound fraction:    {1 - abs(r_resid)/abs(r_raw):.1%} of signal is NumReads confound')

# 2B.3: Sensitivity/Specificity for CN >= 4 prediction
cn_high = (tp_valid['SEQC2_CN'] >= 4).astype(int)
print(f'\n=== HPFineNGroups as CN >= 4 predictor ===')
print(f'CN >= 4: {cn_high.sum()} ({cn_high.mean():.1%}), CN < 4: {(1-cn_high).sum()} ({(1-cn_high).mean():.1%})')

for threshold in [2, 3, 4]:
    pred = (tp_valid['HPFineNGroups'] >= threshold).astype(int)
    tp_count = ((pred == 1) & (cn_high == 1)).sum()
    fp_count = ((pred == 1) & (cn_high == 0)).sum()
    fn_count = ((pred == 0) & (cn_high == 1)).sum()
    tn_count = ((pred == 0) & (cn_high == 0)).sum()
    sens = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    spec = tn_count / (tn_count + fp_count) if (tn_count + fp_count) > 0 else 0
    print(f'  NGroups >= {threshold}: sensitivity={sens:.3f}, specificity={spec:.3f}, '
          f'TP={tp_count}, FP={fp_count}, FN={fn_count}, TN={tn_count}')

# Residualized version
for threshold_pct in [50, 60, 70]:
    threshold_val = np.percentile(tp_valid['HPFineNGroups_resid'], threshold_pct)
    pred = (tp_valid['HPFineNGroups_resid'] >= threshold_val).astype(int)
    tp_count = ((pred == 1) & (cn_high == 1)).sum()
    fp_count = ((pred == 1) & (cn_high == 0)).sum()
    fn_count = ((pred == 0) & (cn_high == 1)).sum()
    tn_count = ((pred == 0) & (cn_high == 0)).sum()
    sens = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
    spec = tn_count / (tn_count + fp_count) if (tn_count + fp_count) > 0 else 0
    print(f'  Residualized >= P{threshold_pct}: sensitivity={sens:.3f}, specificity={spec:.3f}')

# AUC for CN >= 4 prediction
auc_raw = roc_auc_score(cn_high, tp_valid['HPFineNGroups'])
auc_resid = roc_auc_score(cn_high, tp_valid['HPFineNGroups_resid'])
auc_numreads = roc_auc_score(cn_high, tp_valid['NumReads'])
auc_covm = roc_auc_score(cn_high, tp_valid['Coverage_Multiple'])
print(f'\n  AUC for CN >= 4 prediction:')
print(f'    HPFineNGroups (raw):    {auc_raw:.3f}')
print(f'    HPFineNGroups (resid):  {auc_resid:.3f}')
print(f'    NumReads:               {auc_numreads:.3f}')
print(f'    Coverage_Multiple:      {auc_covm:.3f}')

# === Figure 29: HPFineNGroups vs CN box plot ===
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
cn_values = [2, 3, 4, 5, 6, 7, 8]
cn_data_raw = [tp_valid.loc[tp_valid['SEQC2_CN']==c, 'HPFineNGroups'].values for c in cn_values if (tp_valid['SEQC2_CN']==c).sum() >= 10]
cn_labels = [f'CN={c}\n(n={(tp_valid["SEQC2_CN"]==c).sum()})' for c in cn_values if (tp_valid['SEQC2_CN']==c).sum() >= 10]
cn_data_resid = [tp_valid.loc[tp_valid['SEQC2_CN']==c, 'HPFineNGroups_resid'].values for c in cn_values if (tp_valid['SEQC2_CN']==c).sum() >= 10]

bp1 = axes[0].boxplot(cn_data_raw, labels=cn_labels, patch_artist=True)
for patch in bp1['boxes']:
    patch.set_facecolor('steelblue')
    patch.set_alpha(0.6)
axes[0].set_ylabel('HPFineNGroups')
axes[0].set_title(f'Raw (Spearman r={r_raw:.3f})', fontsize=12, fontweight='bold')
axes[0].grid(True, alpha=0.3, axis='y')

bp2 = axes[1].boxplot(cn_data_resid, labels=cn_labels, patch_artist=True)
for patch in bp2['boxes']:
    patch.set_facecolor('coral')
    patch.set_alpha(0.6)
axes[1].set_ylabel('HPFineNGroups (NumReads residualized)')
axes[1].set_title(f'Residualized (Spearman r={r_resid:.3f})', fontsize=12, fontweight='bold')
axes[1].grid(True, alpha=0.3, axis='y')

fig.suptitle('Investigation 2B: HPFineNGroups vs True CN\n(HCC1395 TP, SEQC2 ground truth)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/29_hpfinengroups_vs_cn_boxplot.png', dpi=150, bbox_inches='tight')
print(f'\nSaved: 29_hpfinengroups_vs_cn_boxplot.png')

# === Figure 30: Sensitivity/Specificity curve ===
fig, ax = plt.subplots(figsize=(8, 6))
# ROC-like curve for HPFineNGroups raw and residualized predicting CN >= 4
from sklearn.metrics import roc_curve
fpr_raw, tpr_raw, _ = roc_curve(cn_high, tp_valid['HPFineNGroups'])
fpr_resid, tpr_resid, _ = roc_curve(cn_high, tp_valid['HPFineNGroups_resid'])
fpr_nr, tpr_nr, _ = roc_curve(cn_high, tp_valid['NumReads'])

ax.plot(fpr_raw, tpr_raw, 'b-', linewidth=2, label=f'HPFineNGroups raw (AUC={auc_raw:.3f})')
ax.plot(fpr_resid, tpr_resid, 'r--', linewidth=2, label=f'HPFineNGroups resid (AUC={auc_resid:.3f})')
ax.plot(fpr_nr, tpr_nr, 'g:', linewidth=2, label=f'NumReads (AUC={auc_numreads:.3f})')
ax.plot([0,1], [0,1], 'k--', alpha=0.3)
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate (Sensitivity)')
ax.set_title('Investigation 2B: ROC for CN >= 4 Prediction\n(HCC1395 TP)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/30_cn4_prediction_roc.png', dpi=150, bbox_inches='tight')
print(f'Saved: 30_cn4_prediction_roc.png')


# ============================================================
# Investigation 3A: Clustering Quality by CN Zone
# ============================================================
print('\n' + '='*70)
print('Investigation 3A: Clustering Quality (ClusterPermanovaF) by CN Zone')
print('='*70)

# ClusterPermanovaF by CN
tp_cluster = tp[tp['ClusterPermanovaValid'].astype(str).isin(['True', 'true', '1', '1.0'])].copy()
print(f'Regions with valid ClusterPermanovaF: {len(tp_cluster)}')

cn_zones = [(2, 'CN=2'), (3, 'CN=3'), (4, 'CN=4'), (5, 'CN=5+')]
print(f'\n{"Zone":<10} {"n":>6} {"F mean":>8} {"F median":>10} {"F std":>8}')
print('-' * 50)
zone_data = []
for cn_val, label in cn_zones:
    if cn_val < 5:
        mask = tp_cluster['SEQC2_CN'] == cn_val
    else:
        mask = tp_cluster['SEQC2_CN'] >= cn_val
    n = mask.sum()
    if n < 10: continue
    f_vals = tp_cluster.loc[mask, 'ClusterPermanovaF'].astype(float)
    print(f'{label:<10} {n:>6} {f_vals.mean():>8.2f} {f_vals.median():>10.2f} {f_vals.std():>8.2f}')
    zone_data.append({'Zone': label, 'n': n, 'F_mean': f_vals.mean(), 'F_median': f_vals.median(), 'values': f_vals.values})

# Statistical tests
cn2_f = tp_cluster.loc[tp_cluster['SEQC2_CN']==2, 'ClusterPermanovaF'].astype(float).dropna()
cn4_f = tp_cluster.loc[tp_cluster['SEQC2_CN']==4, 'ClusterPermanovaF'].astype(float).dropna()
cn5p_f = tp_cluster.loc[tp_cluster['SEQC2_CN']>=5, 'ClusterPermanovaF'].astype(float).dropna()

if len(cn2_f) > 10 and len(cn4_f) > 10:
    stat, pval = mannwhitneyu(cn2_f, cn4_f, alternative='two-sided')
    d_cn2_cn4 = (cn2_f.mean() - cn4_f.mean()) / np.sqrt((cn2_f.var() + cn4_f.var()) / 2)
    print(f'\nCN=2 vs CN=4: Mann-Whitney U p={pval:.2e}, Cohen\'s d={d_cn2_cn4:.3f}')

if len(cn2_f) > 10 and len(cn5p_f) > 10:
    stat, pval = mannwhitneyu(cn2_f, cn5p_f, alternative='two-sided')
    d_cn2_cn5 = (cn2_f.mean() - cn5p_f.mean()) / np.sqrt((cn2_f.var() + cn5p_f.var()) / 2)
    print(f'CN=2 vs CN>=5: Mann-Whitney U p={pval:.2e}, Cohen\'s d={d_cn2_cn5:.3f}')

# 3A also: HP_Ratio by CN (expect LOH regions to show extreme HP_Ratio)
print(f'\n=== HP_Ratio Distribution by CN ===')
print(f'{"CN":>5} {"n":>6} {"HP_Ratio mean":>14} {"HP_Ratio std":>13} {"LOH frac":>10}')
print('-' * 55)
for c in [1, 2, 3, 4, 5, 6]:
    if c < 5:
        mask = tp['SEQC2_CN'] == c
    else:
        mask = tp['SEQC2_CN'] >= c
    n = mask.sum()
    if n < 10: continue
    hr = tp.loc[mask, 'HP_Ratio'].astype(float)
    loh = tp.loc[mask, 'Potential_LOH'].astype(str).isin(['True', 'true', '1', '1.0']).mean()
    label = f'CN={c}' if c < 5 else f'CN>={c}'
    print(f'{label:>5} {n:>6} {hr.mean():>14.3f} {hr.std():>13.3f} {loh:>10.1%}')

# === Figure 31: ClusterPermanovaF by CN zone ===
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Violin plot
valid_zones = [zd for zd in zone_data if len(zd['values']) >= 10]
vp = axes[0].violinplot([zd['values'] for zd in valid_zones], showmedians=True, showextrema=False)
for i, pc in enumerate(vp['bodies']):
    pc.set_alpha(0.6)
axes[0].set_xticks(range(1, len(valid_zones)+1))
axes[0].set_xticklabels([zd['Zone'] for zd in valid_zones])
axes[0].set_ylabel('ClusterPermanovaF')
axes[0].set_title('Clustering Quality by CN Zone', fontsize=12, fontweight='bold')
axes[0].grid(True, alpha=0.3, axis='y')

# HP_Ratio distribution by CN
cn_for_hr = [2, 3, 4]
hr_data = []
for c in cn_for_hr:
    mask = tp['SEQC2_CN'] == c
    hr_vals = tp.loc[mask, 'HP_Ratio'].astype(float).dropna().values
    if len(hr_vals) >= 10:
        hr_data.append(hr_vals)
cn5p_hr = tp.loc[tp['SEQC2_CN']>=5, 'HP_Ratio'].astype(float).dropna().values
if len(cn5p_hr) >= 10:
    hr_data.append(cn5p_hr)
    cn_for_hr.append('5+')

vp2 = axes[1].violinplot(hr_data, showmedians=True, showextrema=False)
for pc in vp2['bodies']:
    pc.set_alpha(0.6)
axes[1].set_xticks(range(1, len(hr_data)+1))
axes[1].set_xticklabels([f'CN={c}' for c in cn_for_hr])
axes[1].set_ylabel('HP_Ratio (HP1 / total)')
axes[1].set_title('HP Assignment Balance by CN Zone', fontsize=12, fontweight='bold')
axes[1].axhline(0.5, color='green', linestyle='--', alpha=0.5, label='Balanced (0.5)')
axes[1].legend()
axes[1].grid(True, alpha=0.3, axis='y')

fig.suptitle('Investigation 3A: Read Assignment Quality by CN Zone\n(HCC1395 TP, SEQC2 ground truth)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/31_clustering_quality_by_cn.png', dpi=150, bbox_inches='tight')
print(f'\nSaved: 31_clustering_quality_by_cn.png')

# ============================================================
# Final Summary
# ============================================================
print('\n' + '='*70)
print('SUMMARY & GO/NO-GO DECISIONS')
print('='*70)

# 2A verdict
hp_free_feats = results_2a_df[results_2a_df['Type'] == 'HP-free']
max_resid_hp_free = hp_free_feats['Residualized_Spearman'].abs().max()
best_hp_free = hp_free_feats.loc[hp_free_feats['Residualized_Spearman'].abs().idxmax()]
print(f'\n2A: Best HP-free residualized |r| = {max_resid_hp_free:.3f} ({best_hp_free["Feature"]})')
if max_resid_hp_free >= 0.25:
    print(f'    → SIGNAL: HP-free methylation features have genuine CN correlation')
elif max_resid_hp_free >= 0.15:
    print(f'    → WEAK: Some signal exists but below actionable threshold')
else:
    print(f'    → CONFOUND: HP-free methylation features are CN-blind after NumReads correction')

# 2B verdict
print(f'\n2B: HPFineNGroups CN prediction:')
print(f'    Raw Spearman r = {r_raw:.3f}, Residualized r = {r_resid:.3f}')
print(f'    AUC for CN>=4: raw={auc_raw:.3f}, resid={auc_resid:.3f}, NumReads={auc_numreads:.3f}')
if auc_resid >= 0.60 and auc_resid > auc_numreads * 0.95:
    print(f'    → SIGNAL: HPFineNGroups has independent CN prediction value')
elif auc_resid >= 0.55:
    print(f'    → WEAK: Some independent signal but largely NumReads confound')
else:
    print(f'    → CONFOUND: HPFineNGroups-CN correlation is primarily NumReads artifact')

# 3A verdict
if len(cn2_f) > 10 and len(cn4_f) > 10:
    print(f'\n3A: Clustering quality degradation:')
    print(f'    CN=2 vs CN=4 Cohen\'s d = {d_cn2_cn4:.3f}')
    if abs(d_cn2_cn4) > 0.3:
        print(f'    → SIGNAL: Clustering quality changes with CN (d > 0.3)')
    elif abs(d_cn2_cn4) > 0.1:
        print(f'    → WEAK: Small effect, may not warrant CN-aware parameters')
    else:
        print(f'    → NONE: Clustering quality is CN-independent')

print(f'\n{"="*70}')
print(f'All figures saved to: {FIG_DIR}/')
print(f'Data saved to: {DATA_DIR}/methylation_cn_correlation.tsv')
