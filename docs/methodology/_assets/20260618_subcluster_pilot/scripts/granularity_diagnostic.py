#!/usr/bin/env python3
"""精細度階梯診斷:對指定位點，掃 k=2..6 maxclust，每階看『多出來的群』是否
  (a) balance: 最小群 ≥ grp_min   (b) real: stab_excess(扣 within-1-group null)≥0.10   (c) align: 可靠對齊
→ 找精細度天花板(excess 掉到 null 以下那一階)= 客觀『切到多細才合理』的標準。
產 granularity_curve.png + console 表。import cluster_redesign 重用三閘函式。"""
import os, csv, glob, json, sys
import numpy as np
from scipy.cluster.hierarchy import fcluster
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]=["Droid Sans Fallback","DejaVu Sans"]; plt.rcParams["axes.unicode_minus"]=False

WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; OUT=f"{WT}/output/_kprofile_heatmap"
dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd

def load(key):
    rd=dirmap[key]
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me=CR.loadm(f"{rd}/methylation/methylation.csv"); mi={r:i for i,r in enumerate(mids)}
    it=lambda t:str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and it(reads[r]["is_tumor"]) and reads[r]["hp"] in CR.LABMAP and r in mi]
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=CR.peel(sub)
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]
    P=np.array([Me[mi[r]] for r in kids])
    hp=[CR.LABMAP[reads[r]["hp"]] for r in kids]; alle=[reads[r]["alt_support"] for r in kids]
    return sub,P,hp,alle

KEYS=["chr4_190112507","chr8_138619384","chr8_132099309"]
rng=np.random.default_rng(20260622)
fig,axs=plt.subplots(1,len(KEYS),figsize=(15,4.6))
report={}
for ax,key in zip(axs,KEYS):
    sub,P,hp,alle=load(key); n=sub.shape[0]; gm=CR.grp_min(n); Z,s=CR.linkZ(sub)
    axes={"hp":["HP1" if l in("1","1-1") else "HP2" for l in hp],
          "carrier":["G" if l in("1","2") else "C" for l in hp],
          "allele":[a if a in("REF","ALT") else None for a in alle]}
    ks=list(range(2,7)); excs=[]; rows=[]
    for k in ks:
        cl=fcluster(Z,k,"maxclust"); cnt=Counter(cl); sizes=sorted(cnt.values(),reverse=True)
        minsz=min(sizes); balance=minsz>=gm
        exc=CR.stab_excess(sub,P,k,rng); excs.append(exc if exc is not None else np.nan)
        out_idx=[i for i in range(n) if cnt[cl[i]]<gm]
        bestV=0; beste=0; bestax="-"
        for axn,lab in axes.items():
            V,p,e=CR.align_resolution(cl.tolist(),out_idx,lab)
            if V is not None and V>bestV: bestV,beste,bestax=V,(e or 0),axn
        real=(exc is not None and exc>=0.10)
        align_ok=(bestV>=0.3 and beste>=5)
        rows.append({"k":k,"sizes":sizes,"min":minsz,"balance":balance,"excess":exc,"real":real,
                     "bestaxis":bestax,"V":round(bestV,2),"e":round(beste,1),"align":align_ok})
    report[key]={"n":n,"grp_min":gm,"ladder":rows}
    # curve
    ax.plot(ks,excs,"o-",color="#0d9488",lw=2,label="excess(扣 null)")
    ax.axhline(0.10,ls="--",c="#db2777",lw=1.2,label="real 門檻 0.10")
    for r in rows:
        m="o"; c="#0d9488" if (r["real"] and r["balance"]) else "#db2777"
        ax.annotate(("✓both" if (r["real"] and r["align"] and r["balance"]) else
                     ("real" if r["real"] and r["balance"] else ("min<gm" if not r["balance"] else "弱"))),
                    (r["k"],(r["excess"] or 0)+0.03),ha="center",fontsize=7,color=c)
        ax.annotate(f"min={r['min']}",(r["k"],(r["excess"] or 0)-0.06),ha="center",fontsize=6.5,color="#666")
    ax.set_title(f"{key} n={n} (grp_min={gm})",fontsize=9.5); ax.set_xlabel("k (切幾群)"); ax.set_ylim(-0.15,1.1)
    ax.set_xticks(ks); ax.grid(alpha=0.3); ax.legend(fontsize=7)
axs[0].set_ylabel("stab_excess")
fig.suptitle("精細度階梯:excess 掉到 0.10(null)以下 = 精細度天花板;min<grp_min = 多出的是小群該歸邊緣",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(f"{A}/figs_redesign/granularity_curve.png",dpi=120,bbox_inches="tight")
json.dump(report,open(f"{A}/granularity_diagnostic.json","w"),indent=1,default=str)
for key,R in report.items():
    print(f"=== {key} n={R['n']} grp_min={R['grp_min']} ===")
    for r in R["ladder"]:
        print(f"  k={r['k']} sizes={r['sizes']} min={r['min']} balance={'Y' if r['balance'] else 'N'} "
              f"excess={r['excess']} real={'Y' if r['real'] else 'N'} align={r['bestaxis']}(V={r['V']},e={r['e']},{'ok' if r['align'] else 'no'})")
print("\ncurve → figs_redesign/granularity_curve.png")
