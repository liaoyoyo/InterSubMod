#!/usr/bin/env python3
"""全基因組 TP vs FP 五類對比圖(特異性)。純讀 cluster_redesign_wg_summary.json。"""
import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]=["Droid Sans Fallback","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
s=json.load(open(f"{A}/cluster_redesign_wg_summary.json"))
tp,fp=s["fine_TP"],s["fine_FP"]; ntp,nfp=s["by_set"]["TP"],s["by_set"]["FP"]
order=["CONFIRMED","NEAR_CONFIRMED","REAL_NOVEL","REAL_DIFFUSE","NO_CLEAR_SPLIT"]
lab=["CONFIRMED\n(對齊germline)","NEAR_CONF","REAL_NOVEL\n(subclone候選)","REAL_DIFFUSE","NO_CLEAR\n(真無結構)"]
tpct=[100*tp.get(k,0)/ntp for k in order]; fpct=[100*fp.get(k,0)/nfp for k in order]
en=[(t/f if f>0 else 0) for t,f in zip(tpct,fpct)]
fig,(ax,ax2)=plt.subplots(1,2,figsize=(13,5),gridspec_kw={"width_ratios":[1.5,1]})
x=np.arange(len(order)); w=0.38
ax.bar(x-w/2,tpct,w,label=f"TP somatic (n={ntp})",color="#0d9488",edgecolor="#222")
ax.bar(x+w/2,fpct,w,label=f"FP artifact (n={nfp})",color="#db2777",edgecolor="#222")
for i,(t,f) in enumerate(zip(tpct,fpct)): ax.text(i,max(t,f)+0.8,f"{en[i]:.2f}×",ha="center",fontsize=9,fontweight="bold",
    color="#0d9488" if en[i]>1.3 else ("#db2777" if en[i]<0.8 else "#777"))
ax.set_xticks(x); ax.set_xticklabels(lab,fontsize=8.5); ax.set_ylabel("% of loci in set"); ax.legend(fontsize=9)
ax.set_title("全基因組五類 TP vs FP(數字=TP/FP 富集倍率)",fontsize=11); ax.grid(axis="y",alpha=0.3)
# enrichment 條
cols=["#0d9488" if e>1.3 else "#db2777" if e<0.8 else "#9ca3af" for e in en]
ax2.barh(x,en,color=cols,edgecolor="#222"); ax2.axvline(1,ls="--",c="#333")
ax2.set_yticks(x); ax2.set_yticklabels([k.replace("_SPLIT","") for k in order],fontsize=8.5); ax2.invert_yaxis()
ax2.set_xlabel("TP/FP 富集倍率(>1=TP特異)"); ax2.set_title("特異性:>1 somatic-相關 / <1 非特異",fontsize=10)
for i,e in enumerate(en): ax2.text(e+0.05,i,f"{e:.2f}×",va="center",fontsize=8.5)
fig.suptitle("🔴 CONFIRMED(germline-cis)TP 富集 3.3× = somatic-相關;REAL_NOVEL(subclone候選)0.65×(FP更多)= 非somatic特異",fontsize=10.5,y=1.0)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(f"{A}/figs_dashboard/wg_tpfp_contrast.png",dpi=120,bbox_inches="tight")
print("WROTE wg_tpfp_contrast.png")
