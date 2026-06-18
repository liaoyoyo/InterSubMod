#!/usr/bin/env python3
"""W0/W2/W3: per-case heatmap PNG + enriched metrics for carrier-1-1 sub-cluster candidates.
Output: figs/{pos}_meth.png each + cases.json (locked metrics for workstation spec)."""
import os, csv, glob, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
from collections import Counter

ROOT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_chr1_subcluster_pilot"
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS=f"{ASSET}/figs"; os.makedirs(FIGS,exist_ok=True)
MIN_SZ=3
recs=json.load(open(f"{ASSET}/records.json"))
# candidates: carrier_1-1 sil>=0.4 (the flagged set) + low-sil controls (<=0.25) for calibration
c11=[r for r in recs if r["subset"]=="carrier_1-1" and r["best_sil"] is not None]
cand=sorted([r for r in c11 if r["best_sil"]>=0.4],key=lambda r:-r["best_sil"])
ctrl=sorted([r for r in c11 if r["best_sil"]<=0.25],key=lambda r:r["best_sil"])[:4]
items=cand+ctrl
print(f"candidates sil>=0.4: {len(cand)}; low controls: {len(ctrl)}; total items: {len(items)}")

def _f(x):
    x=x.strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def find_region(pos):
    g=glob.glob(f"{ROOT}/**/chr1_{pos}/*/distance/BERNOULLI/matrix.csv",recursive=True)
    return os.path.dirname(os.path.dirname(os.path.dirname(g[0]))) if g else None
def load(rdir):
    reads=list(csv.DictReader(open(f"{rdir}/reads/reads.tsv"),delimiter="\t"))
    rid2={r["read_id"]:r for r in reads}
    rows=open(f"{rdir}/distance/BERNOULLI/matrix.csv").read().strip().split("\n")
    ids=[]; M=[]
    for ln in rows[1:]:
        p=ln.split(","); ids.append(p[0]); M.append([_f(x) for x in p[1:]])
    ml=open(f"{rdir}/methylation/methylation.csv").read().strip().split("\n")
    meth={}
    for ln in ml[1:]:
        p=ln.split(","); meth[p[0]]=[_f(x) for x in p[1:]]
    return ids,np.array(M),rid2,meth,len(ml[0].split(","))-1
def peel(sub):
    idx=list(range(sub.shape[0]))
    while True:
        S=sub[np.ix_(idx,idx)]; bad=(S<0)|np.isnan(S); np.fill_diagonal(bad,False)
        if not bad.any(): return idx
        idx.remove(idx[int(np.argmax(bad.sum(1)))])
        if len(idx)<MIN_SZ*2: return idx
def subclust(sub):
    np.fill_diagonal(sub,0.0); sub=np.maximum(sub,sub.T)
    Z=linkage(squareform(sub,checks=False),method="average"); best=None
    for k in range(2,min(5,sub.shape[0]//MIN_SZ)+1):
        lab=fcluster(Z,k,criterion="maxclust"); sz=Counter(lab)
        if len(sz)<k or min(sz.values())<MIN_SZ: continue
        s=silhouette_score(sub,lab,metric="precomputed")
        if best is None or s>best[2]: best=(k,lab,s)
    return best

out=[]
for r in items:
    pos=r["pos"]; rdir=find_region(pos)
    if not rdir: continue
    ids,M,rid2,meth,ncpg=load(rdir)
    id2row={rid:i for i,rid in enumerate(ids)}
    c_ids=[rid for rid in ids if rid in rid2 and rid2[rid]["hp"]=="1-1" and rid2[rid]["is_tumor"] in("1","true","True")]
    rows=[id2row[x] for x in c_ids]; sub=M[np.ix_(rows,rows)].copy(); keep=peel(sub)
    kept=[c_ids[i] for i in keep]; bs=subclust(sub[np.ix_(keep,keep)].copy())
    if not bs: continue
    k,labels,sil=bs
    # cluster mean-beta + per-CpG delta + strand purity
    def cmean(cl):
        ids_c=[kept[i] for i in range(len(kept)) if labels[i]==cl]
        prr=[]
        for rid in ids_c:
            v=[x for x in meth.get(rid,[]) if x is not None and not np.isnan(x)]
            if v: prr.append(np.mean(v))
        return np.mean(prr) if prr else np.nan, ids_c
    cls=sorted(set(labels)); m0,ids0=cmean(cls[0]); m1,ids1=cmean(cls[1])
    # per-CpG |delta| between the two largest clusters
    def prof(ids_c):
        P=[]
        for j in range(ncpg):
            vs=[meth[r][j] for r in ids_c if r in meth and j<len(meth[r]) and meth[r][j] is not None and not np.isnan(meth[r][j])]
            P.append(np.mean(vs) if vs else np.nan)
        return P
    p0,p1=prof(ids0),prof(ids1)
    nd=[abs(a-b) for a,b in zip(p0,p1) if not np.isnan(a) and not np.isnan(b)]
    def strand_purity(ids_c):
        st=[rid2[r].get("strand","") for r in ids_c]
        fwd=sum(1 for s in st if s in("+","forward","1")); n=len(st)
        return round(max(fwd,n-fwd)/n,2) if n else 0
    sp0,sp1=strand_purity(ids0),strand_purity(ids1)
    # heatmap
    order=np.argsort(labels,kind="stable"); okept=[kept[i] for i in order]; olab=labels[order]
    mat=np.array([[(meth[rr][j] if (rr in meth and j<len(meth[rr])) else np.nan) for j in range(ncpg)] for rr in okept])
    colok=~np.all(np.isnan(mat),axis=0); mat=mat[:,colok]
    fig,ax=plt.subplots(figsize=(5.2,3.4))
    im=ax.imshow(mat,aspect="auto",cmap="RdYlBu_r",vmin=0,vmax=1,interpolation="nearest")
    for b in np.where(np.diff(olab)!=0)[0]: ax.axhline(b+0.5,color="black",lw=1.4)
    sb=np.zeros((len(okept),2,3)); pal={1:(0.85,0.37,0.34),2:(0.25,0.49,0.36),3:(0.72,0.53,0.18),4:(0.3,0.3,0.6)}
    for i,(rr,l) in enumerate(zip(okept,olab)):
        sb[i,0]=pal.get(int(l),(.5,.5,.5)); st=rid2[rr].get("strand","")
        sb[i,1]=(0.2,0.4,0.8) if st in("+","forward","1") else (0.9,0.6,0.1)
    axsb=ax.inset_axes([-0.13,0,0.1,1],transform=ax.transAxes)
    axsb.imshow(sb,aspect="auto"); axsb.set_xticks([0,1]); axsb.set_xticklabels(["clu","str"],fontsize=7); axsb.set_yticks([])
    ax.set_title(f"chr1:{pos}  k={k} sil={sil:.3f}  n={len(okept)}",fontsize=9)
    ax.set_xlabel(f"CpG ({mat.shape[1]})",fontsize=8); ax.set_yticks([])
    plt.colorbar(im,ax=ax,fraction=0.03,pad=0.02,label="β")
    fig.savefig(f"{FIGS}/{pos}_meth.png",dpi=95,bbox_inches="tight"); plt.close(fig)
    tier="明顯(≥0.5)" if sil>=0.5 else ("邊際(0.4-0.5)" if sil>=0.4 else "低/控制(<0.25)")
    out.append({"pos":pos,"vc":r["vc"],"sil":round(sil,3),"k":int(k),"n_raw":r["n_raw"],"n_complete":r["n_complete"],
        "min_frac":r["min_frac"],"level_dbeta":round(float(m0-m1),3),"cluster_mean0":round(float(m0),3),
        "cluster_mean1":round(float(m1),3),"cpg_max_dbeta":round(float(max(nd)),3) if nd else None,
        "cpg_big":int(sum(1 for d in nd if d>0.3)),"strand_purity_max":round(max(sp0,sp1),2),"tier":tier,
        "png":f"figs/{pos}_meth.png"})
json.dump(out,open(f"{ASSET}/cases.json","w"),ensure_ascii=False,indent=1)
print(f"WROTE {ASSET}/cases.json ({len(out)} cases) + {len(out)} PNGs")
print("tiers:",Counter(c["tier"] for c in out))
