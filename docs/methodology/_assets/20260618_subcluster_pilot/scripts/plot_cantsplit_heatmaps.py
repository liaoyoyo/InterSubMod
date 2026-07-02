#!/usr/bin/env python3
"""切不出驗證 dual-panel: 切不出按 HP 排序看『甲基 mean-shift 對齊 HP 但距離無離散塊』vs 切得出有塊。
+ 距離分佈直方圖(off-diag): 切不出單峰(無 gap) vs 切得出雙峰(有 gap)。SoT ism_heatmap_std。"""
import os, glob, csv, json, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts")
import ism_heatmap_std as H
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
OUT=f"{WT}/output/_cantsplit_heatmap"; A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS=f"{A}/figs_cantsplit"; os.makedirs(FIGS,exist_ok=True)
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
items=json.load(open(f"{A}/cantsplit_examples.json"))["items"]; done=[]
distvals={"切不出·有訊號":[],"切不出·真null":[],"切得出·對照":[]}
for L in items:
    key=L["key"]; rd=dirmap.get(key)
    if not rd: continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me,cpgs=loadm(f"{rd}/methylation/methylation.csv",cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    if len(kids)<MIN_SZ*2: continue
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    if len(kp)<MIN_SZ*2: continue
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=len(kids); meth0=np.array([Me[mi[r]] for r in kids])
    hp=["HP1" if LABMAP[reads[r]["hp"]] in("1","1-1") else "HP2" for r in kids]
    rowmean=[np.nanmean(meth0[i]) if np.isfinite(np.nanmean(meth0[i])) else 0 for i in range(n)]
    split = L["group"]=="切得出·對照"
    cut=None
    if split:
        Z=linkage(squareform(sub,checks=False),method="average"); cut=list(fcluster(Z,2,'maxclust'))
        order=sorted(range(n),key=lambda i:(cut[i],hp[i],rowmean[i]))
    else:
        order=sorted(range(n),key=lambda i:(hp[i],rowmean[i]))   # 切不出按 HP 再 mean-β 看 gradient
    kids=[kids[i] for i in order]; sub=sub[np.ix_(order,order)]; meth=meth0[order]
    # 距離分佈(off-diag upper)
    iu=np.triu_indices(n,1); distvals[L["group"]].extend(sub[iu].tolist())
    sidebars=[]
    if split and cut is not None:
        co=[cut[i] for i in order]; cl=sorted(set(co)); cm={c:H.CLUSTER_COL[i%len(H.CLUSTER_COL)] for i,c in enumerate(cl)}
        sidebars=[("k2",[cm[v] for v in co])]
    sidebars+=H.sidebar_specs({r:reads[r] for r in kids},kids,cluster_of=None,include_tn=False,include_strand=True)
    title=f"[{L['group']}] {L['chrom']}:{L['pos']} n={n} | |germΔβ|={L['germ']} HP-PERMANOVA p={L['hpP']}"
    fn=f"{FIGS}/cs_{L['group']}_{key}.png"
    H.mpl_dual_panel(meth,sub,sidebars,[int(c) for c in cpgs],int(L["pos"]),title,fn,
                     n_cluster=(2 if split else 0),dpi=90)
    done.append({"key":key,"group":L["group"],"chrom":L["chrom"],"pos":L["pos"],"n":n,
        "germ":L["germ"],"hpP":L["hpP"],"png":f"figs_cantsplit/cs_{L['group']}_{key}.png"})
# 距離分佈直方圖
bins=np.linspace(0,1.0,21); hist={}
for g,v in distvals.items():
    h,_=np.histogram(np.clip(v,0,1.0),bins=bins,density=True)
    hist[g]={"density":[round(x,3) for x in h.tolist()],"median":round(float(np.median(v)),3) if v else None,
             "n_pairs":len(v),"iqr":[round(float(np.percentile(v,25)),3),round(float(np.percentile(v,75)),3)] if v else None}
hist["bins"]=[round(x,3) for x in bins.tolist()]
json.dump({"n":len(done),"items":done,"distdist":hist},open(f"{A}/cantsplit_heatmap_index.json","w"),indent=1)
from collections import Counter
print(f"plotted {len(done)}/{len(items)} dual-panel; by group:",dict(Counter(d['group'] for d in done)))
print("距離分佈 median(切不出 vs 切得出):", {g:hist[g]['median'] for g in distvals})
