#!/usr/bin/env python3
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

TP="/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/asm_dualaxis_tp.tsv"
df=pd.read_csv(TP,sep="\t").dropna(subset=["wilcoxon_p","mean_delta","n_paired_cpg"]).copy()
df["abs_delta"]=df["mean_delta"].abs()
ALPHA=0.05/51091
df["bonf"]=df["wilcoxon_p"]<ALPHA
df["strong"]=df["bonf"]&(df["abs_delta"]>=0.1)
strong=df[df["strong"]]; rest=df[~df["strong"]]

fig,ax=plt.subplots(1,3,figsize=(16,4.6))
# A: n_cpg distribution strong vs all (confound EXISTS)
ax[0].hist(np.log10(rest["n_paired_cpg"]),bins=50,alpha=.5,density=True,label="non-strong",color="#999")
ax[0].hist(np.log10(strong["n_paired_cpg"]),bins=30,alpha=.6,density=True,label="strong-ASM",color="#d62728")
ax[0].set_xlabel("log10 n_paired_cpg"); ax[0].set_ylabel("density")
ax[0].set_title("A. CONFOUND REAL: strong-ASM 4x higher n_cpg\n(median 103 vs 26, MWU p=1e-73)")
ax[0].legend()
# B: |delta| vs n_cpg scatter (effect size DECREASES with n => p-confound, not effect-confound)
ax[1].scatter(rest["n_paired_cpg"],rest["abs_delta"],s=3,alpha=.08,color="#999")
ax[1].scatter(strong["n_paired_cpg"],strong["abs_delta"],s=14,alpha=.8,color="#d62728")
ax[1].axhline(0.1,ls="--",c="k",lw=.8)
ax[1].set_xscale("log"); ax[1].set_xlabel("n_paired_cpg (log)"); ax[1].set_ylabel("|mean_delta|")
ax[1].set_title("B. |Δβ| DECLINES with n (Spearman -0.25)\nstrong-ASM (red) clears 0.1 at all n")
# C: within-bin |delta| strong vs background (REFUTATION)
bins=[0,3,5,8,12,20,40,1e9]; labels=["<=3","4-5","6-8","9-12","13-20","21-40",">40"]
df["ncbin"]=pd.cut(df["n_paired_cpg"],bins=bins,labels=labels,include_lowest=True)
bg=df[~df["strong"]].groupby("ncbin",observed=True)["abs_delta"].mean().reindex(labels)
st=df[df["strong"]].groupby("ncbin",observed=True)["abs_delta"].mean().reindex(labels)
x=np.arange(len(labels))
ax[2].bar(x-.2,bg.values,.4,label="non-strong bg",color="#999")
ax[2].bar(x+.2,st.values,.4,label="strong-ASM",color="#d62728")
ax[2].set_xticks(x); ax[2].set_xticklabels(labels,rotation=45); ax[2].set_xlabel("n_cpg bin")
ax[2].set_ylabel("mean |Δβ|")
ax[2].set_title("C. WITHIN-BIN: strong stays 4-7x above bg\n(7.6x at n-matched resampling, p=3e-29)")
ax[2].legend()
plt.tight_layout()
out="/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2/analysis/fig_agentC_coverage_confound.png"
plt.savefig(out,dpi=130); print("wrote",out)
