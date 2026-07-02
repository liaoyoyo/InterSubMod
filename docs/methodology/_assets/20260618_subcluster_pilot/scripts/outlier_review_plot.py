#!/usr/bin/env python3
"""挑 outlier-tolerant 新規則代表案 + dual-panel(標 valid群 vs outlier) + index → 供檢視 HTML。
類別: 救回+對齊(好) / 救回+不對齊(precision問題) / chr1 兩案 / 舊已可分(對照)。"""
import glob, csv, os, json
import numpy as np
import sys; sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts")
import ism_heatmap_std as H
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter
MIN_SZ=3; MAXOUT=0.20; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; FIGS=f"{A}/figs_outlier"; os.makedirs(FIGS,exist_ok=True)
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
def cv(g,c):
    G=sorted(set(g));C=sorted(set(c))
    if len(G)<2 or len(C)<2: return None
    t=np.zeros((len(G),len(C)));gi={x:i for i,x in enumerate(G)};ci={x:i for i,x in enumerate(C)}
    for a,b in zip(g,c): t[gi[a],ci[b]]+=1
    t=t[t.sum(1)>0][:,t.sum(0)>0]
    if min(t.shape)<2: return None
    try: chi2,_,_,_=chi2_contingency(t)
    except: return None
    return float(np.sqrt(chi2/(t.sum()*(min(t.shape)-1))))
def old_split(Z,n):
    for k in range(2,min(6,n//MIN_SZ)+1):
        sz=Counter(fcluster(Z,k,'maxclust'))
        if len(sz)>=k and min(sz.values())>=MIN_SZ: return True
    return False
def new_split(sub,Z,n):
    best=None
    for k in range(2,min(8,n)+1):
        cl=fcluster(Z,k,'maxclust');sz=Counter(cl)
        valid=[c for c,s in sz.items() if s>=MIN_SZ]
        if len(valid)<2: continue
        nout=sum(s for c,s in sz.items() if s<MIN_SZ)
        if nout/n>MAXOUT: continue
        mask=[i for i in range(n) if cl[i] in valid]
        if len(set(cl[i] for i in mask))<2: continue
        try: sil=silhouette_score(sub[np.ix_(mask,mask)],[cl[i] for i in mask],metric="precomputed")
        except: continue
        if best is None or sil>best[0]: best=(sil,len(valid),nout,list(cl),valid)
    return best
def align(sub,kids,reads,cl,valid):
    n=len(kids); mask=[i for i in range(n) if cl[i] in valid]
    lab=[LABMAP[reads[kids[i]]["hp"]] for i in mask]
    axes={"hp":["A" if l in("1","1-1") else "B" for l in lab],"carrier":["G" if l in("1","2") else "C" for l in lab],"allele":[reads[kids[i]]["alt_support"] for i in mask]}
    out=None
    for ax,gl in axes.items():
        idx=[j for j in range(len(mask)) if (ax!="allele" or gl[j] in("REF","ALT"))]
        if len(idx)<6: continue
        gg=[gl[j] for j in idx]; cc=[cl[mask[j]] for j in idx]; vr=cv(gg,cc)
        if vr is None or vr<0.3: continue
        nulls=[cv(list(np.random.default_rng(s).permutation(gg)),cc) for s in range(30)]; nulls=[x for x in nulls if x is not None]
        if nulls and vr>np.percentile(nulls,95) and (out is None or vr>out[1]): out=(ax,round(vr,3))
    return out

mats=glob.glob(f"{WT}/output/_thresh_cal_2122/**/distance/BERNOULLI/matrix.csv",recursive=True)
mats+=glob.glob(f"{WT}/output/_case2_heatmap/**/distance/BERNOULLI/matrix.csv",recursive=True)
buckets={"救回+對齊":[],"救回+不對齊":[],"舊已可分":[],"chr1案":[]}
recs=[]
for mp in mats:
    rd=mp.rsplit("/distance/",1)[0]; rt=f"{rd}/reads/reads.tsv"
    if not os.path.exists(rt): continue
    pos=None
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): pos=part; break
    reads={r["read_id"]:r for r in csv.DictReader(open(rt),delimiter="\t")}
    dids,D=loadm(mp); di={r:i for i,r in enumerate(dids)}
    mcs=f"{rd}/methylation/methylation.csv"
    if not os.path.exists(mcs): continue
    mids,Me,cpgs=loadm(mcs,cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    if len(kids)<6: continue
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    if len(kp)<6: continue
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=len(kids); Z=linkage(squareform(sub,checks=False),method="average")
    o=old_split(Z,n); ns=new_split(sub,Z,n)
    is_chr1 = pos in ("chr1_210988832","chr1_78328503")
    cat=None
    if is_chr1: cat="chr1案"
    elif (not o) and ns: cat="救回+對齊" if align(sub,kids,reads,ns[3],ns[4]) else "救回+不對齊"
    elif o and ns: cat="舊已可分"
    if cat and (len(buckets[cat])<(2 if cat=="chr1案" else 4)):
        recs.append((pos,n,o,ns,cat,kids,sub,Me,mi,cpgs,reads,Z))
        buckets[cat].append(pos)
# plot
idx=[]
for pos,n,o,ns,cat,kids,sub,Me,mi,cpgs,reads,Z in recs:
    if ns: sil,nv,nout,cl,valid=ns
    else:  # chr1_78328503 minority-fail: 用 k=2 cut 顯示
        cl=list(fcluster(Z,2,'maxclust')); sz=Counter(cl); valid=[c for c,s in sz.items() if s>=MIN_SZ]; nv=len(valid); nout=sum(s for c,s in sz.items() if s<MIN_SZ)
    al=align(sub,kids,reads,cl,valid) if (valid and len(valid)>=2) else None
    hp=[LABMAP[reads[r]["hp"]] for r in kids]
    isout=[cl[i] not in valid for i in range(n)]
    rowm=[np.nanmean(Me[mi[kids[i]]]) if np.isfinite(np.nanmean(Me[mi[kids[i]]])) else 0 for i in range(n)]
    order=sorted(range(n),key=lambda i:(isout[i],cl[i],rowm[i]))  # valid 群在前,outlier 最後
    K=[kids[i] for i in order]; S=sub[np.ix_(order,order)]; meth=np.array([Me[mi[r]] for r in K])
    clO=[cl[i] for i in order]; outO=[isout[i] for i in order]
    vmap={c:H.CLUSTER_COL[i%len(H.CLUSTER_COL)] for i,c in enumerate(sorted(valid))}
    cutcol=[("#9ca3af" if outO[j] else vmap.get(clO[j],"#9ca3af")) for j in range(n)]  # outlier=灰
    sidebars=[("new-cut(灰=outlier)",cutcol)]
    sidebars+=H.sidebar_specs({r:reads[r] for r in K},K,cluster_of=None,include_tn=False,include_strand=True)
    al_s=f"對齊 {al[0]} V={al[1]}" if al else "不對齊a-priori"
    title=f"[{cat}] {pos.replace('_',':')} n={n} | 舊:{'可分' if o else 'best_k=None'} → 新:valid群={nv}+outlier={nout}({round(100*nout/n)}%) | {al_s}"
    fn=f"{FIGS}/{pos}.png"
    H.mpl_dual_panel(meth,S,sidebars,[int(c) for c in cpgs],int(pos.split('_')[1]),title,fn,n_cluster=len(valid),dpi=92)
    idx.append({"pos":pos,"cat":cat,"n":n,"old":"可分" if o else "best_k=None","n_valid":nv,"n_outlier":nout,
        "outlier_pct":round(100*nout/n,1),"aligned":al_s,"png":f"figs_outlier/{pos}.png"})
json.dump({"n":len(idx),"items":idx},open(f"{A}/outlier_review_index.json","w"),indent=1)
print(f"plotted {len(idx)} cases; cats:",{k:len(v) for k,v in buckets.items()})
