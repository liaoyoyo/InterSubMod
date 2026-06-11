#!/usr/bin/env python3
"""
Agent C — Adversarial confound check: coverage / CpG-density (n_paired_cpg).

Hypothesis to REFUTE Agent B:
  strong-ASM loci are an artifact of high n_paired_cpg.
  Wilcoxon p is inflated by large n (more CpG = more "data points" = smaller p),
  so "Bonferroni-passing strong-ASM" may just be the high-n_cpg tail, NOT real ASM.

Tests:
  T1. Is n_paired_cpg of strong-ASM shifted vs all TP? (Mann-Whitney + medians)
  T2. Does Wilcoxon p depend on n_paired_cpg (Spearman p vs n)? Does |Δβ| depend on n?
      -- decouple "p driven by n" (confound) from "p driven by effect size" (real).
  T3. Within n_cpg-bin re-test: inside each n_cpg stratum, is |Δβ| of strong-ASM
      still elevated vs same-bin background? If strong-ASM = high-n artifact, the
      |Δβ| separation should COLLAPSE after binning on n.
  T4. n-matched permutation/resampling: draw background records matched on n_cpg
      distribution to strong-ASM; compare |Δβ|. If effect persists at matched n,
      confound refuted.
  T5. Effect-size-only screen (ignore p entirely): if we threshold ONLY on |Δβ|>=0.1
      (no Wilcoxon), is THAT set still n-cpg biased? Tests whether the |Δβ| criterion
      itself smuggles in coverage.
"""
import numpy as np
import pandas as pd
from scipy import stats

TP = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/asm_dualaxis_tp.tsv"

df = pd.read_csv(TP, sep="\t")
n0 = len(df)
df = df.dropna(subset=["wilcoxon_p", "mean_delta", "n_paired_cpg"]).copy()
df["abs_delta"] = df["mean_delta"].abs()
N = len(df)
print(f"[load] {n0} rows -> {N} with valid p+delta+n_cpg ({n0-N} dropped)")

# Bonferroni alpha per Agent B (n_total valid = 51091 -> alpha = 0.05/51091)
ALPHA = 0.05 / 51091
print(f"[alpha] Bonferroni = {ALPHA:.3e}")

df["bonf"] = df["wilcoxon_p"] < ALPHA
df["strong"] = df["bonf"] & (df["abs_delta"] >= 0.1)
strong = df[df["strong"]]
allset = df
print(f"[strong] n_strong = {len(strong)} ; n_all = {N}")

def q(s):
    return dict(n=int(len(s)), mean=float(s.mean()), median=float(s.median()),
                p25=float(s.quantile(.25)), p75=float(s.quantile(.75)),
                p90=float(s.quantile(.90)), max=float(s.max()))

print("\n=== T1: n_paired_cpg distribution strong vs all ===")
print("strong n_cpg:", q(strong["n_paired_cpg"]))
print("all    n_cpg:", q(allset["n_paired_cpg"]))
u, pu = stats.mannwhitneyu(strong["n_paired_cpg"], allset["n_paired_cpg"], alternative="two-sided")
print(f"Mann-Whitney n_cpg strong vs all: U={u:.3e} p={pu:.3e}")
# also strong vs non-strong (cleaner)
nonstrong = df[~df["strong"]]
u2, pu2 = stats.mannwhitneyu(strong["n_paired_cpg"], nonstrong["n_paired_cpg"], alternative="two-sided")
print(f"Mann-Whitney n_cpg strong vs NON-strong: U={u2:.3e} p={pu2:.3e}")
med_ratio = strong["n_paired_cpg"].median() / allset["n_paired_cpg"].median()
print(f"median n_cpg ratio strong/all = {med_ratio:.2f}x")

print("\n=== T2: does Wilcoxon p depend on n_cpg vs on effect size? ===")
# use -log10(p) so larger = more significant
df["neglogp"] = -np.log10(df["wilcoxon_p"].clip(lower=1e-300))
rho_pn, p_pn = stats.spearmanr(df["neglogp"], df["n_paired_cpg"])
rho_pd, p_pd = stats.spearmanr(df["neglogp"], df["abs_delta"])
rho_dn, p_dn = stats.spearmanr(df["abs_delta"], df["n_paired_cpg"])
print(f"Spearman(-log10p , n_cpg)    = {rho_pn:+.3f}  p={p_pn:.2e}   <- if strongly + => p inflated by n (confound)")
print(f"Spearman(-log10p , |delta|)  = {rho_pd:+.3f}  p={p_pd:.2e}   <- if strongly + => p driven by real effect")
print(f"Spearman(|delta| , n_cpg)    = {rho_dn:+.3f}  p={p_dn:.2e}   <- KEY: does effect size itself rise with n?")

print("\n=== T3: WITHIN n_cpg-bin re-test of |delta| separation ===")
# Define n_cpg bins; within each bin compare strong vs non-strong |delta|,
# and (key) compare whether strong-ASM |delta| exceeds same-bin background.
bins = [0,3,5,8,12,20,40,10**9]
labels = ["<=3","4-5","6-8","9-12","13-20","21-40",">40"]
df["ncbin"] = pd.cut(df["n_paired_cpg"], bins=bins, labels=labels, right=True, include_lowest=True)
# re-derive strong slice AFTER ncbin column added (avoid stale-copy KeyError)
strong = df[df["strong"]]
nonstrong = df[~df["strong"]]
rows = []
for lab in labels:
    sub = df[df["ncbin"] == lab]
    if len(sub) == 0:
        continue
    s = sub[sub["strong"]]
    b = sub[~sub["strong"]]
    bg_mean = float(b["abs_delta"].mean()) if len(b) else np.nan
    st_mean = float(s["abs_delta"].mean()) if len(s) else np.nan
    # fraction reaching |delta|>=0.1 within this bin (effect-size only, no p)
    frac_big = float((sub["abs_delta"] >= 0.1).mean())
    rows.append((lab, len(sub), len(s),
                 round(bg_mean,4), round(st_mean,4),
                 round(st_mean - bg_mean,4) if (len(s) and len(b)) else np.nan,
                 round(frac_big,4)))
binsdf = pd.DataFrame(rows, columns=["ncbin","n_total","n_strong","bg_meanAbsD","strong_meanAbsD","strong_minus_bg","frac_absD>=0.1"])
print(binsdf.to_string(index=False))
# Critical: does strong-ASM still exceed bin background EVERYWHERE? (real effect)
# Or does the gap vanish in high-n bins (artifact)?
present = binsdf.dropna(subset=["strong_minus_bg"])
n_bins_pos = int((present["strong_minus_bg"] > 0).sum())
print(f"\nbins where strong |delta| > bin-background |delta|: {n_bins_pos}/{len(present)}")

print("\n=== T3b: strong-ASM count concentrated in high-n bins? ===")
# If strong-ASM is purely high-n artifact, ALL strong should sit in top n bins.
print("strong-ASM count per n_cpg bin:")
print(strong["ncbin"].value_counts().reindex(labels).to_string())
print("\nstrong-ASM as fraction of records WITHIN each bin (strong-rate per bin):")
rate = df.groupby("ncbin", observed=True)["strong"].mean().reindex(labels)
print((rate*100).round(3).to_string())

print("\n=== T4: n_cpg-matched resampling — is |delta| of strong elevated at MATCHED n? ===")
# For each strong locus, sample K background (non-strong) records with SAME n_cpg.
rng = np.random.default_rng(42)
K = 50
bg = df[~df["strong"]]
bg_by_n = {n: g["abs_delta"].values for n, g in bg.groupby("n_paired_cpg")}
matched_means = []
strong_means = []
unmatched = 0
for _, r in strong.iterrows():
    n = r["n_paired_cpg"]
    pool = bg_by_n.get(n)
    if pool is None or len(pool) == 0:
        unmatched += 1
        continue
    draw = rng.choice(pool, size=min(K, len(pool)), replace=len(pool) < K)
    matched_means.append(draw.mean())
    strong_means.append(r["abs_delta"])
strong_means = np.array(strong_means); matched_means = np.array(matched_means)
print(f"matched {len(strong_means)} strong loci (unmatched={unmatched})")
print(f"strong |delta| mean = {strong_means.mean():.4f} ; n-matched background |delta| mean = {matched_means.mean():.4f}")
print(f"ratio strong / matched-bg = {strong_means.mean()/matched_means.mean():.2f}x")
# paired test: each strong vs its own matched-n background mean
w, pw = stats.wilcoxon(strong_means, matched_means)
print(f"paired Wilcoxon strong vs matched-n background: W={w:.1f} p={pw:.3e}")
print(f"strong > its matched-n bg in {int((strong_means>matched_means).sum())}/{len(strong_means)} loci")

print("\n=== T5: does the |delta|>=0.1 criterion ALONE smuggle in n_cpg bias? ===")
big = df[df["abs_delta"] >= 0.1]
small = df[df["abs_delta"] < 0.1]
print("|delta|>=0.1 set n_cpg:", q(big["n_paired_cpg"]))
print("|delta|<0.1  set n_cpg:", q(small["n_paired_cpg"]))
u5, pu5 = stats.mannwhitneyu(big["n_paired_cpg"], small["n_paired_cpg"], alternative="two-sided")
print(f"Mann-Whitney n_cpg (|delta|>=0.1 vs <0.1): U={u5:.3e} p={pu5:.3e}")
print(f"median n_cpg: big={big['n_paired_cpg'].median()} small={small['n_paired_cpg'].median()}")
# DIRECTION CHECK: low-n inflates |mean_delta| (fewer CpG averaged => noisier mean => larger |delta|)
# vs high-n inflates p. These are OPPOSITE biases.

print("\n=== T6: are Agent B's named clear-consistent loci high-n outliers? ===")
named = [("chr17",9577080),("chr9",42376881),("chr1",55277717),("chr22",34684238),
         ("chr17",21303695),("chr12",31601630),("chr8",132411612),("chr17",39675811),
         ("chr8",142166567),("chr22",37696857),("chr16",17774746),("chr7",1169953),
         ("chr17",79992447),("chr2",23386385),("chr10",131311892)]
med_all = df["n_paired_cpg"].median()
for c,p in named:
    sub = df[(df["chrom"]==c) & (df["somatic_pos"]==p)]
    if len(sub)==0:
        print(f"  {c}:{p}  NOT FOUND"); continue
    ncs = sorted(sub["n_paired_cpg"].unique())
    dmax = sub["abs_delta"].max()
    print(f"  {c}:{p:>10}  n_cpg={ncs}  maxAbsD={dmax:.3f}  (cohort median n_cpg={med_all})")

print("\nDONE")
