import json,glob,os,numpy as np
import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"]=["Noto Sans CJK TC","Droid Sans Fallback","DejaVu Sans"];matplotlib.rcParams["axes.unicode_minus"]=False
AD=os.path.dirname(os.path.abspath(__file__))
sites=[]
for fp in glob.glob(AD+"/loh_tumor_vaf_chr*.json"):
    d=json.load(open(fp))
    if d.get('n_loh_het',0)>0: sites.extend(d['sites'])
tv=[s['tumor_vaf'] for s in sites]; nv=[s['normal_vaf'] for s in sites if s['normal_vaf']]
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,4.6))
# 1 tumor VAF 直方圖 (U形=cnLOH)
ax1.hist(tv,bins=40,range=(0,1),color="#C2410C",alpha=0.75,label=f"tumor (median {np.median(tv):.3f})")
ax1.hist(nv,bins=40,range=(0,1),color="#2563EB",alpha=0.5,label=f"normal germline (median {np.median(nv):.3f})")
ax1.axvline(0.5,ls="--",c="gray",lw=1)
ax1.set_xlabel("VAF (該 het 位點 alt 比例)");ax1.set_ylabel("位點數")
ax1.set_title("圖1: SEQC2 LOH 區 het 位點 VAF (2693 位點)\ntumor U形(→0/1, 98.6% cnLOH) vs normal 集中0.5(雜合)")
ax1.legend(fontsize=9)
# 2 normal vs tumor 散點
ax2.scatter(nv,[s['tumor_vaf'] for s in sites if s['normal_vaf']],s=6,alpha=0.3,c="#A16207")
ax2.axhline(0.5,ls="--",c="gray",lw=0.8);ax2.axvline(0.5,ls="--",c="gray",lw=0.8)
ax2.set_xlabel("normal germline VAF (~0.5 雜合)");ax2.set_ylabel("tumor VAF")
ax2.set_title("圖2: normal(雜合) → tumor(失雜合) 映射\nLOH = germline 雜合但 tumor→0/1，符合定義")
fig.suptitle("LOH 區 tumor VAF 系統分布：98.6% 真 cnLOH（更正 V7 單點外推）",fontsize=12,y=1.04)
fig.tight_layout();fig.savefig(AD+"/fig_loh_tumor_vaf.png",dpi=130,bbox_inches="tight");plt.close()
print("saved fig_loh_tumor_vaf.png  (tumor median %.3f, normal median %.3f)"%(np.median(tv),np.median(nv)))
