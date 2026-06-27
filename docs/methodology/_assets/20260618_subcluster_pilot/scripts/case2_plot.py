#!/usr/bin/env python3
"""個案 dual-panel: 按 HP family 排序，看甲基/距離是否真按 HP 分開(+ UPGMA k 對齊cut)。"""
import glob, csv, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts")
import ism_heatmap_std as H
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from collections import Counter
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
OUT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_case2_heatmap"
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
import os; FIGS=f"{A}/figs_case2"; os.makedirs(FIGS,exist_ok=True)
def _f(x):
    x=str(x).strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p,cols=False):
    r=open(p).read().strip().split("\n"); ids=[];M=[];c=r[0].split(",")[1:]
    for ln in r[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return (ids,np.array(M),c) if cols else (ids,np.array(M))
def peel(s):
    idx=list(range(s.shape[0]))
    while True:
        S=s[np.ix_(idx,idx)];b=(S<0)|np.isnan(S);np.fill_diagonal(b,False)
        if not b.any(): return idx
        idx.remove(idx[int(np.argmax(b.sum(1)))])
        if len(idx)<MIN_SZ*2: return idx
dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
for key,title_k in (("chr1_210988832",3),("chr1_78328503",2)):
    rd=dirmap.get(key)
    if not rd: continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me,cpgs=loadm(f"{rd}/methylation/methylation.csv",cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=len(kids); Z=linkage(squareform(sub,checks=False),method="average")
    hpf=[0 if LABMAP[reads[r]["hp"]] in("1","1-1") else 1 for r in kids]
    rowm=[np.nanmean(Me[mi[kids[i]]]) if np.isfinite(np.nanmean(Me[mi[kids[i]]])) else 0 for i in range(n)]
    cut=list(fcluster(Z,title_k,'maxclust'))
    order=sorted(range(n),key=lambda i:(hpf[i],rowm[i]))
    kids=[kids[i] for i in order]; sub=sub[np.ix_(order,order)]; meth=np.array([Me[mi[r]] for r in kids])
    cutO=[cut[i] for i in order]
    cl=sorted(set(cutO)); cm={c:H.CLUSTER_COL[i%len(H.CLUSTER_COL)] for i,c in enumerate(cl)}
    sidebars=[(f"UPGMA k={title_k}",[cm[v] for v in cutO])]
    sidebars+=H.sidebar_specs({r:reads[r] for r in kids},kids,cluster_of=None,include_tn=False,include_strand=True)
    title=f"{key.replace('_',':')} n={n} | 按 HP family 排序 | best_k=None(被拒) 但 HP-PERMANOVA 顯著、UPGMA k={title_k} 對齊HP CramerV=1.0(僅1-2 outlier 害它被拒)"
    H.mpl_dual_panel(meth,sub,sidebars,[int(c) for c in cpgs],int(key.split('_')[1]),title,f"{FIGS}/{key}.png",n_cluster=len(cl),dpi=95)
    print("plotted",key)
