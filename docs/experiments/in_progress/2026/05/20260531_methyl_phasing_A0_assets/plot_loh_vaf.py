import json,os,numpy as np
import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
matplotlib.rcParams["font.sans-serif"]=["Noto Sans CJK TC","Droid Sans Fallback","DejaVu Sans"];matplotlib.rcParams["axes.unicode_minus"]=False
AD=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(AD+"/loh_characterize.json"))
R=d["regions"]
sep=[r for r in R if r["separable"]];notg=[r for r in R if not r["separable"]]
fig,axes=plt.subplots(1,3,figsize=(13,4.2))
# 1 het VAF
for ax,key,ttl,ref in [(axes[0],"het_vaf_median","LOH區 germline het VAF\n(真cnLOH該偏離0.5)",0.5),
                        (axes[1],"hp_balance","HP1/HP2 read 平衡度\n(真LOH該極不平衡→0)",None),
                        (axes[2],"loh_boundary_dist","到LOH邊界距離(bp)\n(邊界效應該小)",None)]:
    sv=[r[key] for r in sep if r[key] is not None];nv=[r[key] for r in notg if r[key] is not None]
    ax.boxplot([nv,sv],labels=[f"分不開\n(n={len(nv)})",f"分得開\n(n={len(sv)})"],patch_artist=True,widths=0.6)
    if ref:ax.axhline(ref,ls="--",c="gray",lw=1)
    ax.set_title(ttl,fontsize=9)
fig.suptitle("LOH 定性：分得開 vs 分不開 LOH 的判別指標（43 區）— 兩者 VAF 都~0.5(都雜合)",fontsize=11,y=1.04)
fig.tight_layout();fig.savefig(AD+"/fig_loh_characterize.png",dpi=130,bbox_inches="tight");plt.close()
print("saved fig_loh_characterize.png")
