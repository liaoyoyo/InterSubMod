#!/usr/bin/env python3
"""Agent C follow-up: decompose the two opposite n_cpg biases.

KEY QUESTION: Bonferroni-ONLY (no |delta| gate) = pure high-n artifact?
If yes, then the |delta|>=0.1 co-criterion is what rescues strong-ASM from
being a coverage artifact. Quantify how much the |delta| gate matters.
"""
import numpy as np, pandas as pd
from scipy import stats

TP = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/asm_dualaxis_tp.tsv"
df = pd.read_csv(TP, sep="\t").dropna(subset=["wilcoxon_p","mean_delta","n_paired_cpg"]).copy()
df["abs_delta"] = df["mean_delta"].abs()
ALPHA = 0.05/51091
df["bonf"] = df["wilcoxon_p"] < ALPHA

bonf_only = df[df["bonf"]]                                  # Agent B's 313
strong    = df[df["bonf"] & (df["abs_delta"]>=0.1)]         # Agent B's 172
bonf_subthresh = df[df["bonf"] & (df["abs_delta"]<0.1)]     # passed p but small effect

print("=== Bonferroni-ONLY (313) vs +|delta|>=0.1 gate (172) ===")
for name,s in [("ALL TP",df),("Bonf-only(313)",bonf_only),
               ("Bonf+|d|>=0.1 STRONG(172)",strong),
               ("Bonf but |d|<0.1 (rejected by gate)",bonf_subthresh)]:
    nc=s["n_paired_cpg"]; ad=s["abs_delta"]
    print(f"{name:38s} n={len(s):5d}  med_n_cpg={nc.median():6.1f}  mean_n_cpg={nc.mean():6.1f}  med|d|={ad.median():.3f}  mean|d|={ad.mean():.3f}")

# How n-biased is each set vs all? (median n_cpg ratio)
med_all = df["n_paired_cpg"].median()
print(f"\nmedian n_cpg ratio vs ALL ({med_all}):")
print(f"  Bonf-only      = {bonf_only['n_paired_cpg'].median()/med_all:.2f}x")
print(f"  Strong(+gate)  = {strong['n_paired_cpg'].median()/med_all:.2f}x")
print(f"  gate-rejected  = {bonf_subthresh['n_paired_cpg'].median()/med_all:.2f}x")

# The gate-rejected set = the part the |delta| gate REMOVES. Is it MORE n-biased?
# If gate-rejected has even higher n_cpg, the |delta| gate is actively stripping
# the pure-coverage artifact (small effect made significant only by huge n).
print(f"\ngate-rejected (Bonf p-pass but |d|<0.1) median n_cpg = {bonf_subthresh['n_paired_cpg'].median()} "
      f"vs strong = {strong['n_paired_cpg'].median()}")
u,p = stats.mannwhitneyu(bonf_subthresh['n_paired_cpg'], strong['n_paired_cpg'], alternative='greater')
print(f"  MWU gate-rejected n_cpg > strong n_cpg (one-sided): U={u:.3e} p={p:.3e}")

# What is the smallest |delta| that still passed Bonferroni, and at what n?
print("\n=== smallest-effect Bonferroni passers (the coverage-driven tail) ===")
tail = bonf_only.nsmallest(10, "abs_delta")[["chrom","somatic_pos","n_paired_cpg","abs_delta","wilcoxon_p"]]
print(tail.to_string(index=False))

# Robustness: re-rank strong-ASM by EFFECT SIZE alone (no p). Top effect-size loci n_cpg?
print("\n=== top-15 by |delta| ALONE (ignore p): are these the high-n ones? ===")
top_d = df.nlargest(15,"abs_delta")[["chrom","somatic_pos","n_paired_cpg","abs_delta","wilcoxon_p"]]
print(top_d.to_string(index=False))
print(f"\ntop-15-by-|delta| median n_cpg = {top_d['n_paired_cpg'].median()} (cohort {med_all}) "
      f"-- LOW n => effect-size tail is the OPPOSITE (noisy small-n), not coverage")
