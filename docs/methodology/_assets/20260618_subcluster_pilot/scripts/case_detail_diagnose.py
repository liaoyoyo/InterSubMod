#!/usr/bin/env python3
"""挖 reject/doubt 案的 split 被什麼驅動: strand 偏差? 少數 CpG? per-group 組成。"""
import glob, csv, os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter
MIN_SZ=3; MAXOUT=0.20; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
TARGETS={"chr22_33271491":"reject","chr22_30104872":"reject","chr22_19212572":"doubt(out19%)",
         "chr22_36962377":"doubt(舊可分)","chr22_39580028":"agree(不對齊)","chr22_22816187":"agree(對齊,對照)"}
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
def new_split(sub,Z,n):
    best=None
    for k in range(2,min(8,n)+1):
        cl=fcluster(Z,k,'maxclust');sz=Counter(cl)
        valid=[c for c,s in sz.items() if s>=MIN_SZ]
        if len(valid)<2: continue
        nout=sum(s for c,s in sz.items() if s<MIN_SZ)
        if nout/n>MAXOUT: continue
        mask=[i for i in range(n) if cl[i] in valid]
        try: sil=silhouette_score(sub[np.ix_(mask,mask)],[cl[i] for i in mask],metric="precomputed")
        except: continue
        if best is None or sil>best[0]: best=(sil,len(valid),nout,list(cl),valid)
    return best
allm=glob.glob(f"{WT}/output/_thresh_cal_2122/**/distance/BERNOULLI/matrix.csv",recursive=True)
dirmap={}
for mp in allm:
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
for tgt,verd in TARGETS.items():
    rd=dirmap.get(tgt)
    if not rd: print(tgt,"no matrix"); continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me,cpgs=loadm(f"{rd}/methylation/methylation.csv",cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=len(kids); Z=linkage(squareform(sub,checks=False),method="average"); ns=new_split(sub,Z,n)
    print(f"\n===== {tgt} [{verd}] n={n} CpG={len(cpgs)} =====")
    if not ns: print("  新規則不可分"); continue
    sil,nv,nout,cl,valid=ns
    mask=[i for i in range(n) if cl[i] in valid]
    strand=[reads[kids[i]].get("strand","?") for i in mask]; clm=[cl[i] for i in mask]
    hp=[LABMAP[reads[kids[i]]["hp"]] for i in mask]; alle=[reads[kids[i]]["alt_support"] for i in mask]
    hpf=["HP1" if h in("1","1-1") else "HP2" for h in hp]
    vs=cv(strand,clm); vh=cv(hpf,clm); va=cv([a for a in alle],clm)
    print(f"  新valid群={nv} outlier={nout}({round(100*nout/n)}%) sil={round(sil,3)}")
    print(f"  cut 對齊: STRAND CramerV={round(vs,3) if vs else None}  HP={round(vh,3) if vh else None}  ALLELE={round(va,3) if va else None}")
    for g in sorted(valid):
        gi=[j for j in range(len(mask)) if clm[j]==g]
        sc=Counter(strand[j] for j in gi); hc=Counter(hpf[j] for j in gi)
        mb=np.nanmean([np.nanmean(Me[mi[kids[mask[j]]]]) for j in gi])
        print(f"    群{g}: n={len(gi)} strand={dict(sc)} HP={dict(hc)} mean-β={mb:.3f}")
    flag=[]
    if vs and vs>0.6: flag.append("🔴STRAND驅動(strand偏差假象)")
    if len(cpgs)<5: flag.append("🔴CpG太少")
    if nout/n>0.15: flag.append("⚠outlier比例高")
    print(f"  → {' '.join(flag) if flag else '無明顯技術假象(可能真epiallele子結構)'}")
