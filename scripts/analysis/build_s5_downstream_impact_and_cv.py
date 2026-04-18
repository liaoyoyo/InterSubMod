#!/usr/bin/env python3
"""
S5: Downstream Impact Mediation Test + Cross-Validation
========================================================
Agent D — 驗證 LOH 差異如何傳播到 ISM 下游指標，以及 cross-validation checks。

S5.1: VerificationClass Concordance Mediation Test
S5.2: QualityScore LOH Mediation
S5.3: Methylation Feature Direction Mediation
CV-1: Self-phasing effect size predicts LOH over-call
CV-2: 7-sample direction consistency
CV-5: Self-phasing LOH disappearance per sample
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, mannwhitneyu, spearmanr

warnings.filterwarnings('ignore')

# ── paths ──────────────────────────────────────────────────────────────────
MASTER = '/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz'
OUTDIR = '/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260402_phasing_causal_chain/'
os.makedirs(OUTDIR, exist_ok=True)

# ── helpers ────────────────────────────────────────────────────────────────
def compute_auc(tp_vals, fp_vals):
    """Mann-Whitney U based AUC. Returns NaN if either group < 5."""
    tp_vals = tp_vals.dropna()
    fp_vals = fp_vals.dropna()
    if len(tp_vals) < 5 or len(fp_vals) < 5:
        return np.nan
    u, _ = mannwhitneyu(tp_vals, fp_vals, alternative='two-sided')
    return u / (len(tp_vals) * len(fp_vals))


def cramers_v(ct):
    """Cramér's V from a contingency table (DataFrame or ndarray)."""
    if isinstance(ct, pd.DataFrame):
        ct = ct.values
    if ct.sum() == 0:
        return np.nan
    chi2, p, dof, _ = chi2_contingency(ct)
    n = ct.sum()
    k = min(ct.shape) - 1
    return np.sqrt(chi2 / (n * k)) if k > 0 else 0.0


# ── load master dataset ──────────────────────────────────────────────────
print('[INFO] Loading master dataset …')
df = pd.read_csv(MASTER, sep='\t', low_memory=False)
print(f'  rows = {len(df):,}   cols = {df.shape[1]}')

# Normalise column access (lowercase aliases already exist in dataset)
# Use lowercase versions: verification_class, quality_score, core_loh_like, etc.
# The dataset has both CamelCase originals and lowercase duplicates.
# We'll use whichever is available.

vc_col = 'verification_class' if 'verification_class' in df.columns else 'VerificationClass'
qs_col = 'quality_score'      if 'quality_score'      in df.columns else 'Quality_Score'
loh_col = 'core_loh_like'
mode_col = 'mode'
truth_col = 'truth_label'
vk_col = 'variant_key'
sample_col = 'sample'
hp_ratio_col = 'HP_Ratio'  # CamelCase original

# Ensure booleans
if df[loh_col].dtype == object:
    df[loh_col] = df[loh_col].map({'True': True, 'False': False, True: True, False: False})
df[loh_col] = df[loh_col].astype(bool)

# Ensure numeric
for c in [qs_col, hp_ratio_col, 'CramersV', 'AlleleDelta', 'HPMergedDelta',
          'PairwiseMedianDist', 'HeuristicScore']:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')

print(f'  modes: {df[mode_col].value_counts().to_dict()}')
print(f'  truth: {df[truth_col].value_counts().to_dict()}')
print(f'  core_loh_like: {df[loh_col].value_counts().to_dict()}')

# ══════════════════════════════════════════════════════════════════════════
# Build same-locus pairs  (variant_key present in both paired and to)
# ══════════════════════════════════════════════════════════════════════════
print('\n[INFO] Building same-locus pairs …')

paired_df = df[df[mode_col] == 'paired'].copy()
to_df     = df[df[mode_col] == 'to'].copy()

common_vk = set(paired_df[vk_col].unique()) & set(to_df[vk_col].unique())
print(f'  variant_keys: paired={paired_df[vk_col].nunique():,}  to={to_df[vk_col].nunique():,}  common={len(common_vk):,}')

paired_sub = paired_df[paired_df[vk_col].isin(common_vk)].copy()
to_sub     = to_df[to_df[vk_col].isin(common_vk)].copy()

# Some variant_keys may have duplicates within the same mode (multiple datasets).
# Keep the row with highest source_priority, then first.
def dedup(sub, mode_label):
    if 'source_priority' in sub.columns:
        sub = sub.sort_values('source_priority', ascending=False)
    return sub.drop_duplicates(subset=[vk_col], keep='first')

paired_sub = dedup(paired_sub, 'paired')
to_sub     = dedup(to_sub, 'to')

pairs = paired_sub.set_index(vk_col).join(
    to_sub.set_index(vk_col),
    lsuffix='_P', rsuffix='_T',
    how='inner'
)
print(f'  same-locus pairs: {len(pairs):,}')

# LOH status categories
pairs['loh_cat'] = 'neither_LOH'
both_loh = pairs[f'{loh_col}_P'] & pairs[f'{loh_col}_T']
to_only  = (~pairs[f'{loh_col}_P']) & pairs[f'{loh_col}_T']
p_only   = pairs[f'{loh_col}_P'] & (~pairs[f'{loh_col}_T'])
pairs.loc[both_loh, 'loh_cat'] = 'both_LOH'
pairs.loc[to_only,  'loh_cat'] = 'TO_only_LOH'
pairs.loc[p_only,   'loh_cat'] = 'paired_only_LOH'

print(f'  loh_cat distribution:\n{pairs["loh_cat"].value_counts().to_string()}')

# ══════════════════════════════════════════════════════════════════════════
# S5.1  VerificationClass Concordance Mediation
# ══════════════════════════════════════════════════════════════════════════
print('\n' + '='*72)
print('S5.1  VerificationClass Concordance Mediation Test')
print('='*72)

vc_p = f'{vc_col}_P'
vc_t = f'{vc_col}_T'

# Overall cross-tab
ct_all = pd.crosstab(pairs[vc_p], pairs[vc_t])
cv_all = cramers_v(ct_all)
print(f'\n[Overall] Cramér\'s V = {cv_all:.4f}   n = {len(pairs):,}')
print(ct_all.to_string())

# Per LOH category
s51_rows = []
for cat in ['neither_LOH', 'TO_only_LOH', 'both_LOH', 'paired_only_LOH']:
    sub = pairs[pairs['loh_cat'] == cat]
    if len(sub) < 10:
        continue
    ct = pd.crosstab(sub[vc_p], sub[vc_t])
    cv = cramers_v(ct)
    # Concordance rate = fraction where VC matches
    concordance = (sub[vc_p] == sub[vc_t]).mean()
    s51_rows.append({
        'loh_category': cat,
        'n_pairs': len(sub),
        'cramers_v': round(cv, 4),
        'concordance_rate': round(concordance, 4),
    })
    print(f'\n[{cat}] Cramér\'s V = {cv:.4f}  concordance = {concordance:.4f}  n = {len(sub):,}')
    print(ct.to_string())

s51_rows.append({
    'loh_category': 'OVERALL',
    'n_pairs': len(pairs),
    'cramers_v': round(cv_all, 4),
    'concordance_rate': round((pairs[vc_p] == pairs[vc_t]).mean(), 4),
})

s51_df = pd.DataFrame(s51_rows)
s51_path = os.path.join(OUTDIR, 's5_1_verification_class_concordance_by_loh.tsv')
s51_df.to_csv(s51_path, sep='\t', index=False)
print(f'\n  → saved {s51_path}')

# Interpretation
neither = s51_df[s51_df['loh_category'] == 'neither_LOH']
overall = s51_df[s51_df['loh_category'] == 'OVERALL']
if len(neither) and len(overall):
    delta_conc = neither.iloc[0]['concordance_rate'] - overall.iloc[0]['concordance_rate']
    print(f'\n  Interpretation: neither_LOH concordance − overall = {delta_conc:+.4f}')
    if delta_conc > 0.02:
        print('  → LOH IS a mediator of VC discordance (concordance improves when LOH removed)')
    else:
        print('  → LOH alone does NOT fully mediate VC discordance')


# ══════════════════════════════════════════════════════════════════════════
# S5.2  QualityScore LOH Mediation
# ══════════════════════════════════════════════════════════════════════════
print('\n' + '='*72)
print('S5.2  QualityScore LOH Mediation')
print('='*72)

qs_p = f'{qs_col}_P'
qs_t = f'{qs_col}_T'

s52_rows = []

# ---- Spearman correlation paired QS vs TO QS by LOH ----
for cat in ['OVERALL', 'neither_LOH', 'TO_only_LOH', 'both_LOH']:
    sub = pairs if cat == 'OVERALL' else pairs[pairs['loh_cat'] == cat]
    valid = sub[[qs_p, qs_t]].dropna()
    if len(valid) < 10:
        continue
    rho, p = spearmanr(valid[qs_p], valid[qs_t])
    row = {'metric': 'QS_spearman_PairedVsTO', 'loh_category': cat,
           'n': len(valid), 'value': round(rho, 4), 'p_value': round(p, 6)}
    s52_rows.append(row)
    print(f'  [{cat}] Spearman(QS_paired, QS_to) = {rho:.4f}  p={p:.2e}  n={len(valid):,}')

# ---- AUC (TP vs FP) by mode × LOH ----
print('\n  AUC (TP vs FP) by mode × LOH status:')
for mode_label, mode_sub in [('paired', paired_df), ('to', to_df)]:
    for loh_val, loh_label in [(None, 'ALL'), (False, 'non_LOH'), (True, 'LOH')]:
        if loh_val is not None:
            sub = mode_sub[mode_sub[loh_col] == loh_val]
        else:
            sub = mode_sub
        tp_qs = sub.loc[sub[truth_col] == 'TP', qs_col]
        fp_qs = sub.loc[sub[truth_col] == 'FP', qs_col]
        auc = compute_auc(tp_qs, fp_qs)
        row = {'metric': f'QS_AUC_{mode_label}', 'loh_category': loh_label,
               'n': len(sub), 'value': round(auc, 4) if not np.isnan(auc) else np.nan,
               'p_value': np.nan}
        s52_rows.append(row)
        print(f'    {mode_label:>6s} | {loh_label:>7s} | AUC = {auc:.4f}  (TP={len(tp_qs):,}  FP={len(fp_qs):,})')

s52_df = pd.DataFrame(s52_rows)
s52_path = os.path.join(OUTDIR, 's5_2_qs_mediation_by_loh.tsv')
s52_df.to_csv(s52_path, sep='\t', index=False)
print(f'\n  → saved {s52_path}')

# Interpretation
to_all_auc = s52_df[(s52_df['metric'] == 'QS_AUC_to') & (s52_df['loh_category'] == 'ALL')]['value']
to_non_auc = s52_df[(s52_df['metric'] == 'QS_AUC_to') & (s52_df['loh_category'] == 'non_LOH')]['value']
if len(to_all_auc) and len(to_non_auc):
    delta = to_non_auc.iloc[0] - to_all_auc.iloc[0]
    print(f'\n  TO QS AUC: non_LOH={to_non_auc.iloc[0]:.4f}  ALL={to_all_auc.iloc[0]:.4f}  delta={delta:+.4f}')
    if delta > 0.02:
        print('  → LOH mediates QS failure (non-LOH QS performs better)')
    else:
        print('  → LOH alone does NOT fully explain QS failure in TO mode')


# ══════════════════════════════════════════════════════════════════════════
# S5.3  Methylation Feature Direction Mediation
# ══════════════════════════════════════════════════════════════════════════
print('\n' + '='*72)
print('S5.3  Methylation Feature Direction Mediation')
print('='*72)

features = ['CramersV', 'AlleleDelta', 'HPMergedDelta', 'PairwiseMedianDist', 'HeuristicScore']
# Check which features exist
features = [f for f in features if f in df.columns]
print(f'  Features available: {features}')

s53_rows = []

for feat in features:
    for mode_label, mode_sub in [('paired', paired_df), ('to', to_df)]:
        for loh_val, loh_label in [(None, 'ALL'), (False, 'non_LOH'), (True, 'LOH')]:
            if loh_val is not None:
                sub = mode_sub[mode_sub[loh_col] == loh_val]
            else:
                sub = mode_sub
            tp_v = sub.loc[sub[truth_col] == 'TP', feat]
            fp_v = sub.loc[sub[truth_col] == 'FP', feat]
            auc = compute_auc(tp_v, fp_v)
            tp_mean = tp_v.mean() if len(tp_v) > 0 else np.nan
            fp_mean = fp_v.mean() if len(fp_v) > 0 else np.nan
            direction = 'TP>FP' if (tp_mean > fp_mean) else ('FP>TP' if fp_mean > tp_mean else 'equal')
            s53_rows.append({
                'feature': feat,
                'mode': mode_label,
                'loh_category': loh_label,
                'auc': round(auc, 4) if not np.isnan(auc) else np.nan,
                'tp_mean': round(tp_mean, 6) if not np.isnan(tp_mean) else np.nan,
                'fp_mean': round(fp_mean, 6) if not np.isnan(fp_mean) else np.nan,
                'direction': direction,
                'n_tp': len(tp_v.dropna()),
                'n_fp': len(fp_v.dropna()),
            })

s53_df = pd.DataFrame(s53_rows)
s53_path = os.path.join(OUTDIR, 's5_3_methylation_feature_mediation.tsv')
s53_df.to_csv(s53_path, sep='\t', index=False)
print(f'\n  → saved {s53_path}')

# Check direction reversal
print('\n  Direction analysis (TP>FP expected in paired):')
reversal_count = 0
mediated_count = 0
for feat in features:
    fd = s53_df[s53_df['feature'] == feat]
    p_all = fd[(fd['mode'] == 'paired') & (fd['loh_category'] == 'ALL')].iloc[0]
    t_all = fd[(fd['mode'] == 'to')     & (fd['loh_category'] == 'ALL')].iloc[0]
    t_non = fd[(fd['mode'] == 'to')     & (fd['loh_category'] == 'non_LOH')].iloc[0]
    p_non = fd[(fd['mode'] == 'paired') & (fd['loh_category'] == 'non_LOH')].iloc[0]

    reversed_all  = (p_all['direction'] != t_all['direction'])
    reversed_non  = (p_non['direction'] != t_non['direction'])

    status = ''
    if reversed_all and not reversed_non:
        status = 'LOH MEDIATES reversal'
        mediated_count += 1
    elif reversed_all and reversed_non:
        status = 'STILL reversed (extra confound beyond LOH)'
        reversal_count += 1
    elif not reversed_all:
        status = 'No reversal (same direction)'

    print(f'    {feat:>22s}:  P_all={p_all["auc"]:.3f}({p_all["direction"]})  '
          f'T_all={t_all["auc"]:.3f}({t_all["direction"]})  '
          f'T_nonLOH={t_non["auc"]:.3f}({t_non["direction"]})  '
          f'→ {status}')

# ── Plot S5.3 ──
fig, axes = plt.subplots(1, len(features), figsize=(4*len(features), 5), sharey=False)
if len(features) == 1:
    axes = [axes]
for ax, feat in zip(axes, features):
    fd = s53_df[s53_df['feature'] == feat].copy()
    fd['label'] = fd['mode'] + '\n' + fd['loh_category']
    bars = ax.bar(range(len(fd)), fd['auc'].values, color=[
        '#1f77b4' if m == 'paired' else '#ff7f0e' for m in fd['mode']
    ], alpha=0.8)
    ax.set_xticks(range(len(fd)))
    ax.set_xticklabels(fd['label'].values, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('AUC (TP vs FP)')
    ax.set_title(feat, fontsize=10)
    ax.axhline(0.5, color='grey', ls='--', lw=0.8)
    ax.set_ylim(0.3, 0.9)

plt.suptitle('S5.3: Feature AUC by Mode × LOH Status', fontsize=13, y=1.02)
plt.tight_layout()
fig_path = os.path.join(OUTDIR, 's5_3_feature_auc_by_loh.png')
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.close()
print(f'\n  → saved {fig_path}')


# ══════════════════════════════════════════════════════════════════════════
# CV: Cross-Validation Checks
# ══════════════════════════════════════════════════════════════════════════
print('\n' + '='*72)
print('CV: Cross-Validation Checks')
print('='*72)

cv_rows = []

# ── CV-1: Self-phasing effect size predicts LOH over-call ──
print('\n--- CV-1: HP imbalance ↔ self-phasing LOH rate (per-sample Spearman) ---')

# Per-sample TO TP mean |HP_Ratio - 0.5| (imbalance)
to_tp = to_df[to_df[truth_col] == 'TP'].copy()
to_tp['hp_imbalance'] = (to_tp[hp_ratio_col] - 0.5).abs()
imbalance_by_sample = to_tp.groupby(sample_col)['hp_imbalance'].mean()

# Self-phasing LOH rate per sample:
# = (TO LOH sites where Paired is NOT LOH) / (total TO LOH sites), from same-locus pairs
to_loh_pairs = pairs[pairs[f'{loh_col}_T'] == True]
selfphasing_rates = {}
for samp in pairs[f'{sample_col}_P'].unique():
    # Filter pairs for this sample (use _P or _T, they should be the same sample)
    sp = to_loh_pairs[to_loh_pairs[f'{sample_col}_P'] == samp]
    if len(sp) == 0:
        continue
    # Self-phasing LOH = TO is LOH but Paired is NOT
    sp_loh = sp[~sp[f'{loh_col}_P']].shape[0]
    selfphasing_rates[samp] = sp_loh / len(sp)

sp_series = pd.Series(selfphasing_rates, name='selfphasing_loh_rate')

# Merge
cv1_df = pd.DataFrame({
    'hp_imbalance': imbalance_by_sample,
    'selfphasing_loh_rate': sp_series
}).dropna()

if len(cv1_df) >= 4:
    rho, p = spearmanr(cv1_df['hp_imbalance'], cv1_df['selfphasing_loh_rate'])
    print(f'  Spearman r = {rho:.4f}  p = {p:.4f}  n_samples = {len(cv1_df)}')
    print(cv1_df.to_string())
    cv_rows.append({'check': 'CV-1', 'metric': 'spearman_r_imbalance_vs_selfphasing',
                    'value': round(rho, 4), 'target': '>0.6',
                    'pass': 'YES' if rho > 0.6 else 'NO', 'n': len(cv1_df)})
else:
    print(f'  Insufficient samples: {len(cv1_df)}')
    cv_rows.append({'check': 'CV-1', 'metric': 'spearman_r_imbalance_vs_selfphasing',
                    'value': np.nan, 'target': '>0.6', 'pass': 'INSUFFICIENT', 'n': len(cv1_df)})

# ── CV-2: 7-sample direction consistency ──
print('\n--- CV-2: 7-sample direction consistency ---')

samples = sorted(df[sample_col].unique())
cv2_checks = {'TO_TP_imbalance_gt_Paired': [],
              'TO_LOH_rate_gt_Paired': [],
              'selfphasing_LOH_exists': []}

for samp in samples:
    # Check 1: TO TP imbalance > Paired TP imbalance
    p_tp = paired_df[(paired_df[sample_col] == samp) & (paired_df[truth_col] == 'TP')]
    t_tp = to_df[(to_df[sample_col] == samp) & (to_df[truth_col] == 'TP')]
    p_imb = (p_tp[hp_ratio_col] - 0.5).abs().mean() if len(p_tp) > 0 else np.nan
    t_imb = (t_tp[hp_ratio_col] - 0.5).abs().mean() if len(t_tp) > 0 else np.nan
    check1 = t_imb > p_imb if not (np.isnan(t_imb) or np.isnan(p_imb)) else None
    cv2_checks['TO_TP_imbalance_gt_Paired'].append(check1)

    # Check 2: TO LOH rate > Paired LOH rate
    p_all = paired_df[paired_df[sample_col] == samp]
    t_all = to_df[to_df[sample_col] == samp]
    p_loh_rate = p_all[loh_col].mean() if len(p_all) > 0 else np.nan
    t_loh_rate = t_all[loh_col].mean() if len(t_all) > 0 else np.nan
    check2 = t_loh_rate > p_loh_rate if not (np.isnan(t_loh_rate) or np.isnan(p_loh_rate)) else None
    cv2_checks['TO_LOH_rate_gt_Paired'].append(check2)

    # Check 3: Self-phasing LOH exists
    check3 = samp in selfphasing_rates and selfphasing_rates[samp] > 0
    cv2_checks['selfphasing_LOH_exists'].append(check3)

    print(f'  {samp:>20s}: imb_TO>P={str(check1):>5s}  LOH_TO>P={str(check2):>5s}  '
          f'sp_LOH={str(check3):>5s}  '
          f'(P_imb={p_imb:.4f} T_imb={t_imb:.4f} P_LOH={p_loh_rate:.3f} T_LOH={t_loh_rate:.3f})')

for check_name, results in cv2_checks.items():
    valid = [r for r in results if r is not None]
    yes_count = sum(1 for r in valid if r)
    consistency = yes_count / len(valid) if valid else 0
    cv_rows.append({'check': 'CV-2', 'metric': check_name,
                    'value': round(consistency, 4),
                    'target': '7/7 = 1.0',
                    'pass': 'YES' if consistency == 1.0 else f'{yes_count}/{len(valid)}',
                    'n': len(valid)})
    print(f'  {check_name}: {yes_count}/{len(valid)} = {consistency:.2%}')


# ── CV-5: Per-sample self-phasing LOH disappearance rate ──
print('\n--- CV-5: Self-phasing LOH disappearance rate per sample ---')
# Self-phasing LOH disappearance = among sites that are TO LOH + Paired NOT LOH,
# how many is the rate? Already computed in selfphasing_rates.
# But the user wants "disappearance rate" = the fraction of TO LOH that disappears
# when self-phasing is removed = same as selfphasing_rates.
cv5_all_pass = True
for samp in sorted(selfphasing_rates.keys()):
    rate = selfphasing_rates[samp]
    passed = rate > 0.30
    if not passed:
        cv5_all_pass = False
    print(f'  {samp:>20s}: self-phasing LOH disappearance = {rate:.3f}  '
          f'{"PASS" if passed else "FAIL"} (target > 0.30)')
    cv_rows.append({'check': 'CV-5', 'metric': f'selfphasing_disappear_{samp}',
                    'value': round(rate, 4), 'target': '>0.30',
                    'pass': 'YES' if passed else 'NO',
                    'n': int(to_loh_pairs[to_loh_pairs[f'{sample_col}_P'] == samp].shape[0])})

# Overall
overall_rate = np.mean(list(selfphasing_rates.values()))
print(f'  {"MEAN":>20s}: {overall_rate:.3f}')
cv_rows.append({'check': 'CV-5', 'metric': 'selfphasing_disappear_MEAN',
                'value': round(overall_rate, 4), 'target': '>0.30',
                'pass': 'YES' if cv5_all_pass else 'PARTIAL', 'n': len(selfphasing_rates)})


# ── Save CV summary ──
cv_df = pd.DataFrame(cv_rows)
cv_path = os.path.join(OUTDIR, 'cv_validation_summary.tsv')
cv_df.to_csv(cv_path, sep='\t', index=False)
print(f'\n  → saved {cv_path}')


# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print('\n' + '='*72)
print('SUMMARY')
print('='*72)

print('\n── S5.1: VerificationClass Concordance Mediation ──')
for _, r in s51_df.iterrows():
    print(f'  {r["loh_category"]:>20s}: Cramér\'s V = {r["cramers_v"]:.4f}  concordance = {r["concordance_rate"]:.4f}  n = {r["n_pairs"]:,}')
neither_row = s51_df[s51_df['loh_category'] == 'neither_LOH']
overall_row = s51_df[s51_df['loh_category'] == 'OVERALL']
if len(neither_row) and len(overall_row):
    d = neither_row.iloc[0]['concordance_rate'] - overall_row.iloc[0]['concordance_rate']
    print(f'  → neither_LOH − OVERALL concordance delta = {d:+.4f}')
    print(f'  → Conclusion: LOH {"IS" if d > 0.02 else "is NOT fully"} a mediator of VC discordance')

print('\n── S5.2: QualityScore LOH Mediation ──')
for _, r in s52_df.iterrows():
    print(f'  {r["metric"]:>30s} | {r["loh_category"]:>10s}: value = {r["value"]:.4f}  n = {r["n"]:,}')
to_all_row = s52_df[(s52_df['metric'] == 'QS_AUC_to') & (s52_df['loh_category'] == 'ALL')]
to_non_row = s52_df[(s52_df['metric'] == 'QS_AUC_to') & (s52_df['loh_category'] == 'non_LOH')]
if len(to_all_row) and len(to_non_row):
    d = to_non_row.iloc[0]['value'] - to_all_row.iloc[0]['value']
    print(f'  → TO QS AUC: non_LOH − ALL = {d:+.4f}')
    print(f'  → Conclusion: LOH {"mediates" if d > 0.02 else "does NOT fully explain"} QS failure in TO mode')

print('\n── S5.3: Feature Direction Mediation ──')
for feat in features:
    fd = s53_df[s53_df['feature'] == feat]
    p_all = fd[(fd['mode'] == 'paired') & (fd['loh_category'] == 'ALL')].iloc[0]
    t_all = fd[(fd['mode'] == 'to')     & (fd['loh_category'] == 'ALL')].iloc[0]
    t_non = fd[(fd['mode'] == 'to')     & (fd['loh_category'] == 'non_LOH')].iloc[0]
    rev_all = p_all['direction'] != t_all['direction']
    rev_non = fd[(fd['mode'] == 'paired') & (fd['loh_category'] == 'non_LOH')].iloc[0]['direction'] != t_non['direction']
    tag = ''
    if rev_all and not rev_non:
        tag = '→ LOH mediates'
    elif rev_all and rev_non:
        tag = '→ Extra confound'
    else:
        tag = '→ No reversal'
    print(f'  {feat:>22s}: P={p_all["auc"]:.3f}({p_all["direction"]})  '
          f'T={t_all["auc"]:.3f}({t_all["direction"]})  '
          f'T_nonLOH={t_non["auc"]:.3f}({t_non["direction"]})  {tag}')

print('\n── CV Validation Summary ──')
for _, r in cv_df.iterrows():
    v = r['value']
    v_str = f'{v:.4f}' if not (isinstance(v, float) and np.isnan(v)) else 'NaN'
    print(f'  [{r["check"]}] {r["metric"]:>45s}: {v_str}  target={r["target"]}  pass={r["pass"]}')

print('\n── Output Files ──')
for f in sorted(os.listdir(OUTDIR)):
    fpath = os.path.join(OUTDIR, f)
    if os.path.isfile(fpath):
        sz = os.path.getsize(fpath)
        print(f'  {f}  ({sz:,} bytes)')

print('\n[DONE] S5 Downstream Impact Mediation + CV completed.')
