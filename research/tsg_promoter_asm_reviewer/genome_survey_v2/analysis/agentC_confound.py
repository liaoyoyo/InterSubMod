#!/usr/bin/env python3
"""
Agent C -- Adversarial confound test: permutation null + ALLELE-axis directional bias.

Goal: try to REFUTE Agent B's strong-ASM / clear-consistent loci as confound artifacts
rather than true somatic ASM.

Focus:
  (1) Permutation null -- if we shuffle germ/som labels, what does the strong-ASM
      count distribution look like? Does observed 172 strong-ASM exceed the null?
  (2) ALLELE-axis directional bias -- if ALT-allele reads have systematic mapping/length
      bias making ALT-end beta artificially low, hypo should systematically dominate hyper
      on the ALLELE axis. Test sign asymmetry, and contrast with HP axis (which does NOT
      involve an allele switch -> if bias were allele-specific, HP axis should be cleaner).
  (3) ALLELE vs HP magnitude/coverage confound for the clear-consistent loci.

read-only: only reads asm_dualaxis_tp.tsv / asm_dualaxis_fp.tsv, writes to analysis/.
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2'
tp = pd.read_csv(f'{BASE}/asm_dualaxis_tp.tsv', sep='\t')
fp = pd.read_csv(f'{BASE}/asm_dualaxis_fp.tsv', sep='\t')

rng = np.random.default_rng(20260530)

# Agent B thresholds (recompute Bonferroni on valid-p)
N_valid = tp['wilcoxon_p'].notna().sum()
ALPHA_BONF = 0.05 / N_valid
DELTA_THR = 0.10
print(f'[setup] N_valid_p={N_valid}  Bonferroni alpha={ALPHA_BONF:.3e}  delta_thr={DELTA_THR}')

def strong_mask(df):
    return (df['wilcoxon_p'] < ALPHA_BONF) & (df['max_abs_delta'].abs() >= DELTA_THR)

# We use max_abs_delta as the effect (Agent B used |Δβ|>=0.1; for sign we use mean_delta)
obs_strong = (tp['wilcoxon_p'] < ALPHA_BONF) & (tp['mean_delta'].abs() >= DELTA_THR)
n_obs_strong = int(obs_strong.sum())
print(f'[obs] strong-ASM (Bonf AND |mean_delta|>=0.1) = {n_obs_strong}')

# =====================================================================
# TEST 1: Direction asymmetry on ALLELE vs HP axis
# If ALT reads carry a systematic hypo bias, ALLELE-axis hypo% >> 50% (and >> HP-axis hypo%)
# =====================================================================
print('\n===== TEST 1: ALLELE-axis directional (hypo-bias) confound =====')
def dir_break(df, label):
    nz = df[df['mean_delta'] != 0]
    hypo = (nz['mean_delta'] < 0).sum(); hyper = (nz['mean_delta'] > 0).sum()
    tot = hypo + hyper
    # binomial test vs 50/50
    bt = stats.binomtest(hypo, tot, 0.5)
    print(f'  [{label}] n={tot}  hypo={hypo} ({100*hypo/tot:.2f}%)  hyper={hyper} ({100*hyper/tot:.2f}%)  '
          f'signed_mean_delta={nz["mean_delta"].mean():+.4f}  binom_p(vs50/50)={bt.pvalue:.3g}')
    return hypo, hyper, tot, nz['mean_delta'].mean()

allele = tp[tp['axis_type']=='ALLELE']
hp = tp[tp['axis_type']=='HP']
dir_break(tp, 'ALL records')
a_hypo,a_hyper,a_tot,a_mean = dir_break(allele, 'ALLELE-axis (ALT_vs_REF)')
h_hypo,h_hyper,h_tot,h_mean = dir_break(hp, 'HP-axis (HP*_vs_*-1)')

# Strong-only direction by axis
print('  --- strong-ASM only (Bonf AND |delta|>=0.1) by axis ---')
for label, sub in [('ALLELE', allele), ('HP', hp)]:
    s = sub[(sub['wilcoxon_p']<ALPHA_BONF) & (sub['mean_delta'].abs()>=DELTA_THR)]
    if len(s):
        hypo=(s['mean_delta']<0).sum(); hyper=(s['mean_delta']>0).sum()
        bt=stats.binomtest(hypo,len(s),0.5)
        print(f'    [{label}] strong n={len(s)} hypo={hypo}({100*hypo/len(s):.1f}%) hyper={hyper}({100*hyper/len(s):.1f}%) binom_p={bt.pvalue:.3g}')

# Cross-axis: is ALLELE more hypo-skewed than HP? (2x2 hypo/hyper x axis)
cont = np.array([[a_hypo, a_hyper],[h_hypo, h_hyper]])
chi2, pchi, _, _ = stats.chi2_contingency(cont)
print(f'  [ALLELE vs HP hypo/hyper 2x2] chi2={chi2:.3f} p={pchi:.4g}  '
      f'(ALLELE hypo%={100*a_hypo/a_tot:.2f} vs HP hypo%={100*h_hypo/h_tot:.2f})')
print(f'  => If ALT-read hypo bias real: expect ALLELE hypo% >> HP hypo% AND ALLELE signed_mean << 0')

# =====================================================================
# TEST 2: Permutation null -- shuffle germ/som within each record
# The per-record stat is mean_delta = mean_som_beta - mean_germ_beta.
# Record-level permutation: under H0 (no ASM), germ and som labels are exchangeable,
# so sign of delta is random -> a label flip gives -delta with prob 0.5.
# We simulate the null strong-ASM COUNT by randomly flipping signs of mean_delta and
# asking how many records would still pass |delta|>=0.1. BUT magnitude is preserved under
# sign flip, so this tests SIGN structure not magnitude. The sharper null: the p-value
# field is derived from within-locus read-level wilcoxon; we cannot re-permute reads from
# aggregates. So we do TWO complementary nulls:
#   (2a) Sign-flip null: tests whether the DIRECTION concordance / asymmetry is real.
#   (2b) p-value uniformity null: under global H0, wilcoxon_p ~ Uniform(0,1).
#        Count expected p<ALPHA_BONF by chance and compare to observed Bonferroni passers.
# =====================================================================
print('\n===== TEST 2a: Sign-flip permutation (direction structure) =====')
# Observed signed mean across ALLELE axis
obs_allele_signed = allele['mean_delta'].mean()
B = 20000
flip_means = np.empty(B)
ad = allele['mean_delta'].values
for b in range(B):
    signs = rng.choice([-1,1], size=len(ad))
    flip_means[b] = (ad*signs).mean()
p_signflip = (np.abs(flip_means) >= abs(obs_allele_signed)).mean()
print(f'  ALLELE signed_mean_delta obs={obs_allele_signed:+.5f}  '
      f'sign-flip null mean={flip_means.mean():+.5f} sd={flip_means.std():.5f}  '
      f'p(|null|>=|obs|)={p_signflip:.4g}')
print(f'  => Large p means observed signed bias is NOT distinguishable from random sign -> '
      f'no systematic directional (hypo) artifact at the population level.')

print('\n===== TEST 2b: p-value uniformity / Bonferroni null =====')
valid_p = tp['wilcoxon_p'].dropna().values
n_bonf_obs = (valid_p < ALPHA_BONF).sum()
exp_bonf_null = N_valid * ALPHA_BONF  # expected under uniform
print(f'  observed Bonferroni passers={n_bonf_obs}  expected under Uniform null={exp_bonf_null:.3f}')
print(f'  enrichment over null = {n_bonf_obs/max(exp_bonf_null,1e-9):.1f}x')
# Poisson tail prob of seeing >= n_bonf_obs under null
from scipy.stats import poisson
p_poisson = poisson.sf(n_bonf_obs-1, exp_bonf_null)
print(f'  Poisson P(>= {n_bonf_obs} | lambda={exp_bonf_null:.3f}) = {p_poisson:.3g}')
# strong-ASM observed vs null (Bonf passers that ALSO randomly have |delta|>=0.1)
frac_big_delta = (tp['mean_delta'].abs() >= DELTA_THR).mean()
exp_strong_null = exp_bonf_null * frac_big_delta
print(f'  frac records |delta|>=0.1 (global) = {frac_big_delta:.4f}')
print(f'  expected strong-ASM under null ~ {exp_strong_null:.3f}  vs observed {n_obs_strong}  '
      f'-> enrichment {n_obs_strong/max(exp_strong_null,1e-9):.1f}x')

# =====================================================================
# TEST 3: ALLELE-axis hypo NOT explained by germ-beta ceiling (mathematical bound)
# If germ_beta is near 1.0, delta can only go down (hypo) -> structural asymmetry source.
# Conversely if germ_beta near 0, only hyper possible. Check whether observed asymmetry is
# explained by germ_beta distribution (a confound), independent of any ALT-read bias.
# =====================================================================
print('\n===== TEST 3: germ-beta ceiling/floor structural asymmetry =====')
def ceiling_diag(df, label):
    g = df['mean_germ_beta']
    # potential downside = g (room to go hypo), potential upside = 1-g (room to go hyper)
    print(f'  [{label}] mean germ_beta={g.mean():.3f} median={g.median():.3f}  '
          f'frac germ_beta>0.9={ (g>0.9).mean():.3f}  frac<0.1={ (g<0.1).mean():.3f}')
    # mechanistic prediction: high germ_beta records should be hypo-skewed regardless of allele
    hi = df[df['mean_germ_beta']>0.9]; lo = df[df['mean_germ_beta']<0.1]
    for n2,s2 in [('germ_beta>0.9',hi),('germ_beta<0.1',lo)]:
        nz=s2[s2['mean_delta']!=0]
        if len(nz):
            print(f'     {n2}: n={len(nz)} hypo%={100*(nz.mean_delta<0).mean():.1f} hyper%={100*(nz.mean_delta>0).mean():.1f}')
ceiling_diag(allele,'ALLELE')
ceiling_diag(hp,'HP')

# =====================================================================
# TEST 4: Re-examine the 15 clear-consistent loci under permutation/bias lens
# For each: does its strong signal depend on a single axis only? Is direction consistent
# with germ-beta ceiling (mechanistic alt explanation) or genuinely bidirectional?
# =====================================================================
print('\n===== TEST 4: clear-consistent loci robustness audit =====')
cc_loci = [
 ('chr17',9577080),('chr9',42376881),('chr1',55277717),('chr22',34684238),
 ('chr17',21303695),('chr12',31601630),('chr8',132411612),('chr17',39675811),
 ('chr8',142166567),('chr22',37696857),('chr16',17774746),('chr7',1169953),
 ('chr17',79992447),('chr2',23386385),('chr10',131311892)]
rows=[]
for ch,pos in cc_loci:
    sub = tp[(tp['chrom']==ch)&(tp['somatic_pos']==pos)]
    al = sub[sub['axis_type']=='ALLELE']
    hpx = sub[sub['axis_type']=='HP']
    # pick maxabs HP sub-axis
    hp_use = hpx.loc[hpx['mean_delta'].abs().idxmax()] if len(hpx) else None
    al_use = al.loc[al['mean_delta'].abs().idxmax()] if len(al) else None
    d_hp = hp_use['mean_delta'] if hp_use is not None else np.nan
    d_al = al_use['mean_delta'] if al_use is not None else np.nan
    g_al = al_use['mean_germ_beta'] if al_use is not None else np.nan
    # ceiling artifact flag: hypo AND germ_beta very high, OR hyper AND germ_beta very low
    ceil_flag=''
    if al_use is not None:
        if d_al<0 and g_al>0.9: ceil_flag='HYPO@high-germ(ceiling-suspect)'
        elif d_al>0 and g_al<0.1: ceil_flag='HYPER@low-germ(floor-suspect)'
        else: ceil_flag='mid-range-OK'
    same_sign = (np.sign(d_hp)==np.sign(d_al)) if (not np.isnan(d_hp) and not np.isnan(d_al)) else None
    rows.append(dict(locus=f'{ch}:{pos}', d_hp=round(d_hp,3) if not np.isnan(d_hp) else None,
                     d_allele=round(d_al,3) if not np.isnan(d_al) else None,
                     germ_beta_allele=round(g_al,3) if not np.isnan(g_al) else None,
                     same_sign=same_sign, ceiling_flag=ceil_flag,
                     n_cpg_allele=int(al_use['n_paired_cpg']) if al_use is not None else None))
cc = pd.DataFrame(rows)
pd.set_option('display.width',200, 'display.max_columns',20)
print(cc.to_string(index=False))
n_ceil_suspect = cc['ceiling_flag'].str.contains('suspect').sum()
n_same = cc['same_sign'].sum()
print(f'\n  clear-consistent: {n_same}/{len(cc)} dual-axis same-sign survive; '
      f'{n_ceil_suspect}/{len(cc)} flagged ceiling/floor structural-suspect')

# =====================================================================
# TEST 5: KEY adversarial -- is ALLELE-axis hypo enriched at HIGH coverage (a read-bias signature)?
# Mapping/length bias would scale with read depth differences. Check correlation of |delta|
# and direction with n_paired_cpg and with germ/som beta asymmetry.
# =====================================================================
print('\n===== TEST 5: ALLELE hypo vs coverage (read-bias signature) =====')
al = allele.copy()
al['is_hypo'] = (al['mean_delta']<0).astype(int)
# does hypo probability rise with coverage? logistic-ish via correlation
r_cov = stats.pointbiserialr(al['is_hypo'], al['n_paired_cpg'])
print(f'  point-biserial(is_hypo, n_paired_cpg) r={r_cov.correlation:+.4f} p={r_cov.pvalue:.3g}')
# split by coverage quartile
al['cov_q'] = pd.qcut(al['n_paired_cpg'], 4, labels=['Q1','Q2','Q3','Q4'], duplicates='drop')
print(al.groupby('cov_q')['mean_delta'].agg(['mean','median','count']))
print(al.groupby('cov_q')['is_hypo'].mean())
print('  => If read-bias real, hypo fraction should rise monotonically with coverage. Flat = no bias.')

cc.to_csv(f'{BASE}/analysis/agentC_clear_consistent_audit.csv', index=False)
print('\n[written] agentC_clear_consistent_audit.csv')
