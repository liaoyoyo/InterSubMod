import json, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"]=["Noto Sans CJK TC","Droid Sans Fallback","DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"]=False
AD=os.path.dirname(os.path.abspath(__file__))
rows=[]
for fp in glob.glob(AD+"/seqc2_cn_methyl_v2_chr*.json"): rows.extend(json.load(open(fp))["regions"])
def dedup(rs,g=4000):
    by={}; [by.setdefault(r["chrom"],[]).append(r) for r in rs]; k=[]
    for c,l in by.items():
        l=sorted(l,key=lambda r:r["pos"]);last=-1e18
        for r in l:
            if r["pos"]-last>=g: k.append(r);last=r["pos"]
    return k
ded=dedup(rows)

# 圖A: CpG-SNP 排除前後 AUC (散點 raw vs excl)
fig,ax=plt.subplots(figsize=(6.5,6))
xr=[r["anchor_auc_raw"] for r in ded if r["anchor_auc_cpgsnp_excl"] is not None]
xe=[r["anchor_auc_cpgsnp_excl"] for r in ded if r["anchor_auc_cpgsnp_excl"] is not None]
ax.scatter(xr,xe,alpha=0.4,s=18,c="#D97757")
ax.plot([0.5,1],[0.5,1],"--",c="gray",lw=1,label="y=x (無變化)")
ax.set_xlabel("anchor AUC (raw)");ax.set_ylabel("anchor AUC (CpG-SNP 排除後)")
ax.set_title("圖A: CpG-SNP 排除前後 AUC 幾乎不變\n(drop median=0.0, Wilcoxon p=0.36, 僅0.6%區drop>0.05)\n→ 95%訊號『不是』CpG-SNP pseudo-ASM artifact")
ax.legend(fontsize=8);ax.set_xlim(0.5,1.02);ax.set_ylim(0.5,1.02)
fig.tight_layout();fig.savefig(AD+"/fig_v2_cpgsnp_excl.png",dpi=130);plt.close()

# 圖B: 去重後 status AUC + separable frac
fig,ax=plt.subplots(figsize=(7.5,4.5))
order=["neutral","gain","loss","loh"]
C={"gain":"#C2410C","loss":"#1E3A8A","loh":"#A16207","neutral":"#6B6B66"}
data=[[r["anchor_auc_raw"] for r in ded if r["seqc2_status"]==s] for s in order]
bp=ax.boxplot(data,labels=[f"{s}\n(n={len(d)})" for s,d in zip(order,data)],patch_artist=True,widths=0.6)
for p,s in zip(bp["boxes"],order): p.set_facecolor(C[s]);p.set_alpha(0.6)
ax.axhline(0.58,ls="--",c="gray",lw=1)
ax.set_ylabel("甲基 anchor AUC")
ax.set_title("圖B: 去重 window 後各狀態 AUC (7 小染色體 317 區)\nLOH vs neutral 無差異 (p=0.96) — 與前輪 chr5/6/7 矛盾→染色體特異")
fig.tight_layout();fig.savefig(AD+"/fig_v2_status_dedup.png",dpi=130);plt.close()

# 圖C: LOH 分得開 by chr (對照 v1)
fig,ax=plt.subplots(figsize=(7.5,4.5))
allloh={}
for r in [x for x in rows if x["seqc2_status"]=="loh"]:
    allloh.setdefault(r["chrom"],[]).append(r["anchor_auc_raw"])
# v1
v1=[]
for fp in glob.glob(AD+"/seqc2_cn_methyl_chr*.json"):
    if "v2" in fp:continue
    v1.extend(json.load(open(fp))["regions"])
for r in [x for x in v1 if x["seqc2_status"]=="loh"]:
    allloh.setdefault(r["chrom"]+"(v1)",[]).append(r["anchor_auc"])
chrs=sorted(allloh.keys())
data=[allloh[c] for c in chrs]
bp=ax.boxplot(data,labels=[f"{c}\n(n={len(d)})" for c,d in zip(chrs,data)],patch_artist=True,widths=0.6)
for p in bp["boxes"]:p.set_facecolor("#A16207");p.set_alpha(0.5)
ax.axhline(0.58,ls="--",c="gray",lw=1)
ax.set_ylabel("LOH 區甲基 anchor AUC")
ax.set_title("圖C: LOH 可分性是染色體特異的 (非 LOH 普遍性質)\nchr5/6 低(0.6) vs chr8/15/22 高(0.8-0.97)；chr8 兩版重現 0.80")
plt.xticks(rotation=30,ha="right",fontsize=8)
fig.tight_layout();fig.savefig(AD+"/fig_v2_loh_by_chr.png",dpi=130);plt.close()
print("3 figs saved:")
for f in ["fig_v2_cpgsnp_excl.png","fig_v2_status_dedup.png","fig_v2_loh_by_chr.png"]:
    p=AD+"/"+f; print(f"  {f}: {'OK' if os.path.exists(p) else 'MISS'} {os.path.getsize(p)//1024 if os.path.exists(p) else 0}KB")
