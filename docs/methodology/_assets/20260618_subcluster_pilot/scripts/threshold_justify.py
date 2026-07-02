#!/usr/bin/env python3
"""為何 sil>=0.4: best_sil 分布圖 + 門檻敏感度(0.3..0.5)。+ 沒訊號拆解已在 nosignal_breakdown.json。
產 figs (base64-ready PNG) + sensitivity.json。"""
import json, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS=f"{ASSET}/figs"
SIG="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_wg_bdcprime_verify/significance_summary.csv"
recs=json.load(open(f"{ASSET}/records_wg2.json")); N=len(recs)
NOISE={"Noise_Uniform","Noise_Uncorrelated","Noise_Chaotic"}
MIN_CELL=3; TH_WITHIN=0.5

def has_within(w): return any(w[l]["best_sil"] is not None and w[l]["best_sil"]>=TH_WITHIN for l in ("1-1","2-1"))
def corr(a,th):
    """return (has_struct, one2one, one2many, many2one) at TH_struct=th."""
    hs = a["best_k"] is not None and a["best_k"]>=2 and a["best_sil"] is not None and a["best_sil"]>=th
    if not (hs and a["contingency"]): return hs,False,False,False
    clu=defaultdict(set); lab=defaultdict(set)
    for L,d in a["contingency"].items():
        for cl,c in d.items():
            if c>=MIN_CELL: clu[cl].add(L); lab[L].add(cl)
    o1m=any(len(s)>=2 for s in clu.values()); m1=any(len(s)>=2 for s in lab.values())
    return hs,(not o1m and not m1),o1m,m1

# ---- 1. best_sil distribution (loci with a valid all-read split) ----
sils=[r["all"]["best_sil"] for r in recs if r["all"]["best_sil"] is not None]
nosplit=sum(1 for r in recs if r["all"]["best_sil"] is None)
fig,ax=plt.subplots(figsize=(6.4,3.4))
ax.hist(sils,bins=40,range=(0,1),color="#4A77A8",alpha=0.8,edgecolor="white",linewidth=0.3)
for th,col,lab in [(0.3,"#B8862F","0.3"),(0.4,"#C0432F","0.4 (採用)"),(0.5,"#3F7E5B","0.5")]:
    ax.axvline(th,color=col,lw=1.8,ls="--"); ax.text(th,ax.get_ylim()[1]*0.92,lab,color=col,fontsize=9,ha="center",rotation=0)
ge04=sum(1 for s in sils if s>=0.4)
ax.set_xlabel("all-read best silhouette (有切出≥2群各≥3 的位點)",fontsize=9)
ax.set_ylabel("位點數",fontsize=9)
ax.set_title(f"silhouette 分布: 有切={len(sils)} / 切不出(best_k=None)={nosplit}({round(100*nosplit/N)}%)\n≥0.4: {ge04} 個 = 採用門檻落在右尾",fontsize=9.5)
plt.rcParams["font.sans-serif"]=["Noto Sans CJK TC","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False
fig.tight_layout(); fig.savefig(f"{FIGS}/sil_distribution.png",dpi=100); plt.close(fig)

# ---- 2. threshold sensitivity ----
THS=[0.30,0.35,0.40,0.45,0.50]
sens=[]
for th in THS:
    hs=o2o=o1m=m1=0
    for r in recs:
        a,b,c,d=corr(r["all"],th)
        hs+=a; o2o+=b; o1m+=c; m1+=d
    nosig=sum(1 for r in recs if not (r["all"]["best_k"] is not None and r["all"]["best_k"]>=2 and r["all"]["best_sil"] is not None and r["all"]["best_sil"]>=th) and not has_within(r["within"]))
    sens.append({"th":th,"has_struct":hs,"one2one":o2o,"one2many":o1m,"many2one":m1,"nosig":nosig})
json.dump({"N":N,"rows":sens},open(f"{ASSET}/sensitivity.json","w"),indent=1)

fig,ax=plt.subplots(figsize=(6.4,3.4))
x=[s["th"] for s in sens]
for key,col,lab in [("has_struct","#4A77A8","有結構"),("one2one","#3F7E5B","1對1"),
                    ("many2one","#D97757","多對1"),("one2many","#8a8a8a","1對多")]:
    ax.plot(x,[100*s[key]/N for s in sens],"o-",color=col,label=lab,lw=1.8,ms=5)
ax.axvline(0.4,color="#C0432F",lw=1.2,ls=":")
ax.set_xlabel("結構門檻 TH_struct (silhouette)",fontsize=9); ax.set_ylabel("% of 30490 位點",fontsize=9)
ax.set_title("門檻敏感度: 各類比例 vs 門檻 (曲線平滑→非 knife-edge)\n0.4(紅虛線)是穩定區段, 非陡降點",fontsize=9.5)
ax.legend(fontsize=8,loc="upper right"); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(f"{FIGS}/threshold_sensitivity.png",dpi=100); plt.close(fig)

print("WROTE figs/sil_distribution.png + threshold_sensitivity.png + sensitivity.json")
print(f"\n=== 門檻敏感度 (% of {N}) ===")
print(f"{'TH':>5}{'有結構':>9}{'1對1':>8}{'多對1':>8}{'沒訊號':>9}")
for s in sens:
    print(f"{s['th']:>5}{100*s['has_struct']/N:>8.1f}%{100*s['one2one']/N:>7.1f}%{100*s['many2one']/N:>7.1f}%{100*s['nosig']/N:>8.1f}%")
