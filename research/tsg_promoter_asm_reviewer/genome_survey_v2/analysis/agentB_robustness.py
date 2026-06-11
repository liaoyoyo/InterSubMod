#!/usr/bin/env python3
"""
Agent B - Data validation: within-sample robustness analysis of dual-axis ASM
A-characterization. HCC1395 paired_full single sample.

Inputs (read-only):
  ../asm_dualaxis_tp.tsv  (51,171 dual-axis records, 30,511 unique variants)
  ../asm_dualaxis_fp.tsv  (5,149 records, caller false-positive)

Columns: chrom, somatic_pos, axis, axis_type, loh_status, n_paired_cpg,
         mean_germ_beta, mean_som_beta, mean_delta(=som-germ), max_abs_delta, wilcoxon_p

Outputs to ./ :
  - fig1_delta_distribution.png  (hypo/hyper delta dist, TP)
  - fig2_dualaxis_concordance.png (HP-axis vs ALLELE-axis scatter)
  - fig3_tp_vs_fp_strongASM.png  (TP vs FP strong-ASM proportion)
  - prints all numbers with provenance
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2"
OUT = BASE + "/analysis"

DELTA_THRESH = 0.10  # effect size |mean_delta| >= 0.1

def load(path):
    df = pd.read_csv(path, sep="\t")
    # ensure numeric
    for c in ["n_paired_cpg","mean_germ_beta","mean_som_beta","mean_delta","max_abs_delta","wilcoxon_p"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

tp = load(BASE + "/asm_dualaxis_tp.tsv")
fp = load(BASE + "/asm_dualaxis_fp.tsv")

print("="*80)
print("AGENT B - WITHIN-SAMPLE ROBUSTNESS (A-characterization)")
print("source TP: asm_dualaxis_tp.tsv  | source FP: asm_dualaxis_fp.tsv")
print("="*80)
print(f"TP records: {len(tp)}  | FP records: {len(fp)}")
print(f"TP unique variants (chrom,somatic_pos): {tp[['chrom','somatic_pos']].drop_duplicates().shape[0]}")
print(f"FP unique variants (chrom,somatic_pos): {fp[['chrom','somatic_pos']].drop_duplicates().shape[0]}")
print("\naxis breakdown TP:\n", tp['axis'].value_counts().to_string())
print("\naxis_type x loh_status TP:\n", tp.groupby(['axis_type','loh_status']).size().to_string())
print("\nwilcoxon_p NA count TP:", tp['wilcoxon_p'].isna().sum(), " FP:", fp['wilcoxon_p'].isna().sum())

# ============================================================
# 1. MULTIPLE TESTING (col: wilcoxon_p, mean_delta) -- on TP records
# ============================================================
print("\n" + "="*80)
print("[1] MULTIPLE TESTING  (TP, col wilcoxon_p + mean_delta)")
print("="*80)
# Drop NA p for testing-count denominators (record NA separately)
tp_p = tp['wilcoxon_p'].dropna()
N = len(tp_p)
N_all = len(tp)
n_na_p = tp['wilcoxon_p'].isna().sum()
print(f"N records with valid wilcoxon_p: {N}  (of {N_all}; {n_na_p} NA)")

p05 = (tp_p < 0.05).sum()
bonf_alpha = 0.05 / N
bonf_pass = (tp_p < bonf_alpha).sum()
print(f"p<0.05 uncorrected: {p05}  ({100*p05/N:.2f}% of valid-p)")
print(f"Bonferroni alpha = 0.05/{N} = {bonf_alpha:.3e}")
print(f"Bonferroni pass (p < {bonf_alpha:.3e}): {bonf_pass}  ({100*bonf_pass/N:.2f}%)")

# BH-FDR 5%
def bh_fdr(pvals, q=0.05):
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * (np.arange(1, n+1) / n)
    passed = ranked <= thresh
    if not passed.any():
        return np.zeros(n, dtype=bool), 0
    kmax = np.max(np.where(passed)[0])
    crit = ranked[kmax]
    rej = p <= crit
    return rej, rej.sum()

p_arr = tp_p.values
rej, fdr_pass = bh_fdr(p_arr, 0.05)
print(f"BH-FDR 5% pass: {fdr_pass}  ({100*fdr_pass/N:.2f}%)")

# strong-ASM = Bonferroni AND |mean_delta| >= 0.1  (record level)
tp_valid = tp.dropna(subset=['wilcoxon_p']).copy()
tp_valid['bonf'] = tp_valid['wilcoxon_p'] < bonf_alpha
tp_valid['delta_ok'] = tp_valid['mean_delta'].abs() >= DELTA_THRESH
tp_valid['strong'] = tp_valid['bonf'] & tp_valid['delta_ok']
strong_n = tp_valid['strong'].sum()
print(f"strong-ASM (Bonferroni AND |delta|>={DELTA_THRESH}): {strong_n}  ({100*strong_n/N:.2f}% of valid-p)")
print(f"  of which |delta|>={DELTA_THRESH} alone: {tp_valid['delta_ok'].sum()} ({100*tp_valid['delta_ok'].sum()/N:.2f}%)")

# ============================================================
# 2. DUAL-AXIS CONCORDANCE  -- key consistency evidence
# ============================================================
print("\n" + "="*80)
print("[2] DUAL-AXIS CONCORDANCE  (col axis, mean_delta; per unique variant)")
print("="*80)
# For each unique variant: gather HP-axis delta (HP1_vs_HP1-1 or HP2_vs_HP2-1)
#                          and ALLELE-axis delta (ALT_vs_REF)
hp_mask = tp['axis'].isin(['HP1_vs_HP1-1','HP2_vs_HP2-1'])
allele_mask = tp['axis'] == 'ALT_vs_REF'
print(f"HP-axis records: {hp_mask.sum()}  | ALLELE-axis records: {allele_mask.sum()}")

# HP-axis: a variant may have HP1 and/or HP2 sub-axis. Aggregate per variant
# by taking the delta with max |delta| (dominant signal), but also compute
# mean to be robust. Use mean of HP deltas per variant.
hp = tp[hp_mask].groupby(['chrom','somatic_pos']).agg(
    hp_delta_mean=('mean_delta','mean'),
    hp_delta_maxabs=('mean_delta', lambda s: s.loc[s.abs().idxmax()]),
    hp_p_min=('wilcoxon_p','min'),
    hp_n_axes=('mean_delta','size'),
    hp_loh=('loh_status','first'),
).reset_index()

allele = tp[allele_mask].groupby(['chrom','somatic_pos']).agg(
    al_delta=('mean_delta','mean'),
    al_p=('wilcoxon_p','min'),
    al_loh=('loh_status','first'),
).reset_index()

merged = pd.merge(hp, allele, on=['chrom','somatic_pos'], how='inner')
print(f"Unique variants with BOTH HP-axis AND ALLELE-axis: {len(merged)}")

# concordance = same sign of delta (use HP mean-delta vs allele delta)
def sign_concord(a, b):
    # exclude exact zeros
    m = (a != 0) & (b != 0)
    same = (np.sign(a[m]) == np.sign(b[m]))
    return same.sum(), m.sum()

# Use maxabs HP delta as the dominant HP signal
same_all, valid_all = sign_concord(merged['hp_delta_maxabs'].values, merged['al_delta'].values)
print(f"\nALL dual-axis variants (n={valid_all} nonzero both):")
print(f"  same-sign (concordant) = {same_all}  ({100*same_all/valid_all:.2f}%)")
# Also report using HP mean delta
same_mean, valid_mean = sign_concord(merged['hp_delta_mean'].values, merged['al_delta'].values)
print(f"  [using HP mean-delta]: {same_mean}/{valid_mean} = {100*same_mean/valid_mean:.2f}%")

# Pearson / Spearman of the two axes
m_nonzero = (merged['hp_delta_maxabs']!=0)&(merged['al_delta']!=0)
r_p, p_p = stats.pearsonr(merged.loc[m_nonzero,'hp_delta_maxabs'], merged.loc[m_nonzero,'al_delta'])
r_s, p_s = stats.spearmanr(merged.loc[m_nonzero,'hp_delta_maxabs'], merged.loc[m_nonzero,'al_delta'])
print(f"  Pearson r (HP maxabs vs ALLELE) = {r_p:.3f} (p={p_p:.2e})")
print(f"  Spearman rho = {r_s:.3f} (p={p_s:.2e})")

# concordance within strong-ASM subset:
# strong = both axes pass Bonferroni AND |delta|>=0.1
merged['hp_strong'] = (merged['hp_p_min'] < bonf_alpha) & (merged['hp_delta_maxabs'].abs() >= DELTA_THRESH)
merged['al_strong'] = (merged['al_p'] < bonf_alpha) & (merged['al_delta'].abs() >= DELTA_THRESH)
strong_both = merged[merged['hp_strong'] & merged['al_strong']].copy()
print(f"\nStrong-ASM on BOTH axes (Bonf AND |delta|>={DELTA_THRESH}): n={len(strong_both)}")
if len(strong_both) > 0:
    s2, v2 = sign_concord(strong_both['hp_delta_maxabs'].values, strong_both['al_delta'].values)
    print(f"  same-sign = {s2}/{v2} = {100*s2/v2:.2f}%")

# concordance within "either-strong" (at least one axis strong) -- broader
either_strong = merged[merged['hp_strong'] | merged['al_strong']].copy()
print(f"\nStrong-ASM on AT LEAST ONE axis: n={len(either_strong)}")
if len(either_strong) > 0:
    s3, v3 = sign_concord(either_strong['hp_delta_maxabs'].values, either_strong['al_delta'].values)
    print(f"  same-sign = {s3}/{v3} = {100*s3/v3:.2f}%")

# Note: HP-axis only valid in nonLOH per methods. Restrict to nonLOH HP variants.
merged_nonloh = merged[merged['hp_loh']=='nonLOH'].copy()
if len(merged_nonloh) > 0:
    s4, v4 = sign_concord(merged_nonloh['hp_delta_maxabs'].values, merged_nonloh['al_delta'].values)
    print(f"\nDual-axis restricted to nonLOH HP variants (HP-axis valid): n={v4}")
    print(f"  same-sign = {s4}/{v4} = {100*s4/v4:.2f}%")

# ============================================================
# 3. LOH vs nonLOH  (col loh_status, axis_type, mean_delta)
# ============================================================
print("\n" + "="*80)
print("[3] LOH vs nonLOH  (col loh_status, axis_type)")
print("="*80)
# ALLELE-axis records are the comparable basis (valid for both LOH & nonLOH)
al = tp[tp['axis']=='ALT_vs_REF'].copy()
al = al.dropna(subset=['wilcoxon_p'])
al['bonf'] = al['wilcoxon_p'] < bonf_alpha
al['delta_ok'] = al['mean_delta'].abs() >= DELTA_THRESH
al['strong'] = al['bonf'] & al['delta_ok']
for grp in ['LOH','nonLOH']:
    sub = al[al['loh_status']==grp]
    if len(sub)==0:
        continue
    print(f"\nALLELE-axis {grp}: n={len(sub)}")
    print(f"  mean |delta| = {sub['mean_delta'].abs().mean():.4f}  median |delta| = {sub['mean_delta'].abs().median():.4f}")
    print(f"  mean delta (signed) = {sub['mean_delta'].mean():.4f}")
    print(f"  hypo (delta<0): {(sub['mean_delta']<0).sum()} ({100*(sub['mean_delta']<0).mean():.1f}%)  hyper: {(sub['mean_delta']>0).sum()} ({100*(sub['mean_delta']>0).mean():.1f}%)")
    print(f"  strong-ASM: {sub['strong'].sum()} ({100*sub['strong'].mean():.2f}%)")
# Mann-Whitney on |delta| LOH vs nonLOH
loh_d = al[al['loh_status']=='LOH']['mean_delta'].abs().dropna()
nonloh_d = al[al['loh_status']=='nonLOH']['mean_delta'].abs().dropna()
if len(loh_d)>0 and len(nonloh_d)>0:
    u, pu = stats.mannwhitneyu(loh_d, nonloh_d, alternative='two-sided')
    print(f"\nMann-Whitney |delta| LOH vs nonLOH (ALLELE-axis): U={u:.0f}, p={pu:.2e}")

# ============================================================
# 4. DIRECTION BREAKDOWN within strong-ASM (record level)
# ============================================================
print("\n" + "="*80)
print("[4] DIRECTION BREAKDOWN in strong-ASM (col mean_delta)")
print("="*80)
ss = tp_valid[tp_valid['strong']]
hypo = (ss['mean_delta']<0).sum()
hyper = (ss['mean_delta']>0).sum()
print(f"strong-ASM records: {len(ss)}")
print(f"  hypo (delta<0): {hypo} ({100*hypo/len(ss):.1f}%)")
print(f"  hyper (delta>0): {hyper} ({100*hyper/len(ss):.1f}%)")
# binomial test vs 50/50
bt = stats.binomtest(hypo, len(ss), 0.5) if hasattr(stats,'binomtest') else None
if bt is not None:
    print(f"  binomial test (hypo vs 50/50): p={bt.pvalue:.2e}")
# breakdown by axis_type
print("  by axis_type:")
for at in ss['axis_type'].unique():
    sub = ss[ss['axis_type']==at]
    print(f"    {at}: n={len(sub)}, hypo={100*(sub['mean_delta']<0).mean():.1f}%, hyper={100*(sub['mean_delta']>0).mean():.1f}%")

# ============================================================
# 5. TP vs FP strong-ASM (B-confirm; expect NEGATIVE)
# ============================================================
print("\n" + "="*80)
print("[5] TP vs FP strong-ASM  (B-confirm; expect NEGATIVE = no discrimination)")
print("="*80)
# Use SAME threshold. Recompute Bonferroni per-set? Use TP-derived bonf_alpha
# for both AND also report FP-own-N bonf for fairness.
fp_valid = fp.dropna(subset=['wilcoxon_p']).copy()
N_fp = len(fp_valid)
fp_bonf_alpha = 0.05 / N_fp
print(f"FP valid-p records: {N_fp}")

# Apply TP-derived alpha to both (common threshold), record-level
def strong_prop(df, alpha):
    b = df['wilcoxon_p'] < alpha
    d = df['mean_delta'].abs() >= DELTA_THRESH
    s = (b & d)
    return s.sum(), len(df), s

print(f"\n-- Common threshold (TP-Bonferroni alpha={bonf_alpha:.3e}) + |delta|>={DELTA_THRESH} --")
tp_s, tp_n, _ = strong_prop(tp_valid, bonf_alpha)
fp_s, fp_n, _ = strong_prop(fp_valid, bonf_alpha)
print(f"  TP strong: {tp_s}/{tp_n} = {100*tp_s/tp_n:.2f}%")
print(f"  FP strong: {fp_s}/{fp_n} = {100*fp_s/fp_n:.2f}%")
# Fisher exact 2x2
table = [[tp_s, tp_n - tp_s],[fp_s, fp_n - fp_s]]
odds, pf = stats.fisher_exact(table)
chi2, pc, dof, exp = stats.chi2_contingency(table)
print(f"  Fisher exact OR={odds:.3f}, p={pf:.3e}")
print(f"  Chi2={chi2:.2f}, p={pc:.3e}")

# Also each-set own-Bonferroni (fair within-set FDR control)
print(f"\n-- Each-set own Bonferroni (TP a={bonf_alpha:.2e}, FP a={fp_bonf_alpha:.2e}) + |delta|>={DELTA_THRESH} --")
tp_s2, _, _ = strong_prop(tp_valid, bonf_alpha)
fp_s2, _, _ = strong_prop(fp_valid, fp_bonf_alpha)
print(f"  TP strong: {tp_s2}/{tp_n} = {100*tp_s2/tp_n:.2f}%")
print(f"  FP strong: {fp_s2}/{fp_n} = {100*fp_s2/fp_n:.2f}%")
table2 = [[tp_s2, tp_n - tp_s2],[fp_s2, fp_n - fp_s2]]
odds2, pf2 = stats.fisher_exact(table2)
print(f"  Fisher exact OR={odds2:.3f}, p={pf2:.3e}")

# Also compare on ALLELE-axis only (cleanest, both sets dominated by ALT_vs_REF)
print(f"\n-- ALLELE-axis only, common threshold --")
tp_al = tp_valid[tp_valid['axis']=='ALT_vs_REF']
fp_al = fp_valid[fp_valid['axis']=='ALT_vs_REF']
tp_als,_,_ = strong_prop(tp_al, bonf_alpha)
fp_als,_,_ = strong_prop(fp_al, bonf_alpha)
print(f"  TP ALLELE strong: {tp_als}/{len(tp_al)} = {100*tp_als/len(tp_al):.2f}%")
print(f"  FP ALLELE strong: {fp_als}/{len(fp_al)} = {100*fp_als/len(fp_al):.2f}%")
tabl3 = [[tp_als, len(tp_al)-tp_als],[fp_als, len(fp_al)-fp_als]]
odds3, pf3 = stats.fisher_exact(tabl3)
print(f"  Fisher exact OR={odds3:.3f}, p={pf3:.3e}")
# mean |delta| compare
print(f"  TP ALLELE mean|delta|={tp_al['mean_delta'].abs().mean():.4f}  FP ALLELE mean|delta|={fp_al['mean_delta'].abs().mean():.4f}")

# ============================================================
# 6. CLEAR & CONSISTENT loci table
#    = Bonferroni (both axes) + |delta|>=0.1 (both) + dual-axis same sign
# ============================================================
print("\n" + "="*80)
print("[6] CLEAR & CONSISTENT loci (Bonf both axes + |delta|>=0.1 both + same sign)")
print("="*80)
cc = strong_both.copy()
# same sign
cc = cc[np.sign(cc['hp_delta_maxabs']) == np.sign(cc['al_delta'])]
# also require nonLOH (HP-axis valid)
cc_nonloh = cc[cc['hp_loh']=='nonLOH'].copy()
print(f"clear&consistent candidates (strong both axes + same sign): {len(cc)}")
print(f"  of which HP-axis-valid (nonLOH): {len(cc_nonloh)}")
# rank by combined effect (mean of |hp| and |al| delta) and significance
cc['combined_abs_delta'] = (cc['hp_delta_maxabs'].abs() + cc['al_delta'].abs())/2
cc['min_p'] = cc[['hp_p_min','al_p']].max(axis=1)  # the larger (weaker) p of the two
cc_sorted = cc.sort_values(['combined_abs_delta'], ascending=False)
top = cc_sorted.head(15)
print("\nTOP clear&consistent loci:")
print(top[['chrom','somatic_pos','hp_loh','hp_delta_maxabs','al_delta','hp_p_min','al_p','combined_abs_delta']].to_string(index=False))

# save full clear&consistent table
cc_sorted.to_csv(OUT+"/clear_consistent_loci.csv", index=False)
print(f"\n[saved] {OUT}/clear_consistent_loci.csv  ({len(cc_sorted)} loci)")

# ============================================================
# FIGURES
# ============================================================
plt.rcParams['font.size'] = 10

# Fig 1: delta distribution hypo/hyper (TP, all valid records)
fig, axes = plt.subplots(1,2, figsize=(12,4.5))
ax = axes[0]
d = tp_valid['mean_delta'].dropna()
ax.hist(d[d<0], bins=60, color='#2c7fb8', alpha=0.8, label=f'hypo (n={(d<0).sum()})')
ax.hist(d[d>0], bins=60, color='#d95f0e', alpha=0.8, label=f'hyper (n={(d>0).sum()})')
ax.axvline(0, color='k', lw=0.8)
ax.axvline(-DELTA_THRESH, color='gray', ls='--', lw=0.8)
ax.axvline(DELTA_THRESH, color='gray', ls='--', lw=0.8)
ax.set_xlabel('mean_delta (som - germ beta)')
ax.set_ylabel('count')
ax.set_title('TP all valid records: delta distribution')
ax.legend()
# strong-ASM only
ax = axes[1]
ds = ss['mean_delta'].dropna()
ax.hist(ds[ds<0], bins=40, color='#2c7fb8', alpha=0.8, label=f'hypo (n={(ds<0).sum()})')
ax.hist(ds[ds>0], bins=40, color='#d95f0e', alpha=0.8, label=f'hyper (n={(ds>0).sum()})')
ax.axvline(0, color='k', lw=0.8)
ax.set_xlabel('mean_delta')
ax.set_ylabel('count')
ax.set_title(f'strong-ASM only (Bonf & |d|>={DELTA_THRESH}): n={len(ss)}')
ax.legend()
plt.tight_layout()
plt.savefig(OUT+"/fig1_delta_distribution.png", dpi=130)
plt.close()
print(f"[fig] {OUT}/fig1_delta_distribution.png")

# Fig 2: dual-axis concordance scatter
fig, ax = plt.subplots(figsize=(6.5,6))
mm = merged[(merged['hp_delta_maxabs']!=0)&(merged['al_delta']!=0)]
# color by strong
col = np.where(mm['hp_strong']&mm['al_strong'], '#cb181d',
       np.where(mm['hp_strong']|mm['al_strong'], '#fb6a4a', '#bdbdbd'))
ax.scatter(mm['hp_delta_maxabs'], mm['al_delta'], s=10, c=col, alpha=0.5, edgecolors='none')
ax.axhline(0, color='k', lw=0.6); ax.axvline(0, color='k', lw=0.6)
lim = 1.05
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.plot([-1,1],[-1,1], color='gray', ls=':', lw=0.8)
ax.set_xlabel('HP-axis mean_delta (maxabs sub-axis)')
ax.set_ylabel('ALLELE-axis mean_delta (ALT_vs_REF)')
ax.set_title(f'Dual-axis concordance  n={len(mm)}  same-sign={100*same_all/valid_all:.1f}%\n'
             f'Pearson r={r_p:.2f}  (red=strong both, orange=strong one)')
# shade concordant quadrants
ax.axhspan(0,lim, xmin=0.5, xmax=1.0, color='green', alpha=0.04)
ax.axhspan(-lim,0, xmin=0.0, xmax=0.5, color='green', alpha=0.04)
plt.tight_layout()
plt.savefig(OUT+"/fig2_dualaxis_concordance.png", dpi=130)
plt.close()
print(f"[fig] {OUT}/fig2_dualaxis_concordance.png")

# Fig 3: TP vs FP strong-ASM proportion (common threshold, overall + ALLELE-only)
fig, ax = plt.subplots(figsize=(6.5,5))
labels = ['Overall\n(common Bonf)', 'ALLELE-axis only']
tp_props = [100*tp_s/tp_n, 100*tp_als/len(tp_al)]
fp_props = [100*fp_s/fp_n, 100*fp_als/len(fp_al)]
x = np.arange(len(labels)); w=0.35
b1 = ax.bar(x-w/2, tp_props, w, label=f'TP (n={tp_n})', color='#3182bd')
b2 = ax.bar(x+w/2, fp_props, w, label=f'FP (n={fp_n})', color='#e6550d')
ax.bar_label(b1, fmt='%.1f%%', padding=2, fontsize=9)
ax.bar_label(b2, fmt='%.1f%%', padding=2, fontsize=9)
ax.set_ylabel('strong-ASM proportion (%)')
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_title(f'TP vs FP strong-ASM  (Fisher OR={odds:.2f}, p={pf:.1e})\n'
             'Expected NEGATIVE: ASM common in FP too')
ax.legend()
plt.tight_layout()
plt.savefig(OUT+"/fig3_tp_vs_fp_strongASM.png", dpi=130)
plt.close()
print(f"[fig] {OUT}/fig3_tp_vs_fp_strongASM.png")

# ============================================================
# Dump a machine-readable summary
# ============================================================
import json
summary = {
    "n_tp_records": int(N_all), "n_tp_valid_p": int(N), "n_tp_na_p": int(n_na_p),
    "p05_uncorrected": int(p05), "bonf_alpha": float(bonf_alpha),
    "bonf_pass": int(bonf_pass), "fdr_pass": int(fdr_pass),
    "strong_asm": int(strong_n), "strong_asm_pct": float(100*strong_n/N),
    "dual_axis_n": int(valid_all), "dual_axis_concord_pct": float(100*same_all/valid_all),
    "dual_axis_pearson_r": float(r_p),
    "strong_both_n": int(len(strong_both)),
    "strong_both_concord_pct": float(100*s2/v2) if len(strong_both)>0 else None,
    "direction_hypo": int(hypo), "direction_hyper": int(hyper),
    "direction_hypo_pct": float(100*hypo/len(ss)),
    "tp_strong_pct": float(100*tp_s/tp_n), "fp_strong_pct": float(100*fp_s/fp_n),
    "tpfp_fisher_or": float(odds), "tpfp_fisher_p": float(pf),
    "tp_allele_strong_pct": float(100*tp_als/len(tp_al)), "fp_allele_strong_pct": float(100*fp_als/len(fp_al)),
    "tpfp_allele_fisher_or": float(odds3), "tpfp_allele_fisher_p": float(pf3),
    "clear_consistent_n": int(len(cc_sorted)),
}
with open(OUT+"/agentB_summary.json","w") as f:
    json.dump(summary, f, indent=2)
print(f"\n[saved] {OUT}/agentB_summary.json")
print("\nDONE.")
