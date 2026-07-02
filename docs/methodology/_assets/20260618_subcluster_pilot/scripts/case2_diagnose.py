#!/usr/bin/env python3
"""診斷「HP family 可分但 best_k=None」個案: 追 UPGMA 各 k 切群大小(為何被拒) + HP-supervised 分離度。"""
import glob, csv, sys
import numpy as np
sys.path.insert(0,"/big7_disk/liaoyoyo2001/InterSubMod/scripts")
import ism_heatmap_std as H
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from collections import Counter
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
OUT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_case2_heatmap"
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
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
    chi2,p,_,_=chi2_contingency(t); return float(np.sqrt(chi2/(t.sum()*(min(t.shape)-1)))),float(p)
def permF(D2,lab):
    n=len(lab);gs=set(lab)
    SS_T=np.triu(D2,1).sum()/n; SS_W=0
    for g in gs:
        idx=[i for i in range(n) if lab[i]==g];ng=len(idx)
        if ng<2: continue
        SS_W+=np.triu(D2[np.ix_(idx,idx)],1).sum()/ng
    K=len(gs);
    if n-K<=0 or SS_W<=0: return None
    return ((SS_T-SS_W)/(K-1))/(SS_W/(n-K))

dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd
for key in ("chr1_210988832","chr1_78328503"):
    rd=dirmap.get(key)
    if not rd: print(key,"no matrix"); continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    mids,Me,cpgs=loadm(f"{rd}/methylation/methylation.csv",cols=True); mi={r:i for i,r in enumerate(mids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP and r in mi]
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    kids=[kids[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=len(kids); Z=linkage(squareform(sub,checks=False),method="average")
    hpf=[0 if LABMAP[reads[r]["hp"]] in("1","1-1") else 1 for r in kids]  # HP-family
    print(f"\n===== {key}  n(after peel)={n}  HP1-fam={hpf.count(0)} HP2-fam={hpf.count(1)} =====")
    print("  --- UPGMA 各 k 切群大小(為何 best_k=None) ---")
    for k in range(2,min(6,n//MIN_SZ)+1):
        cl=fcluster(Z,k,'maxclust'); sz=sorted(Counter(cl).values(),reverse=True)
        ok = len(sz)>=k and min(sz)>=MIN_SZ
        al=cv(hpf,list(cl))
        print(f"    k={k}: cluster sizes={sz}  {'✓通過' if ok else '✗拒(有群<3)'}  對齊HP CramerV={round(al[0],3) if al else None}")
    # HP-supervised 分離度
    D2=sub**2
    Fhp=permF(D2,hpf)
    win1=np.triu(D2[np.ix_([i for i in range(n) if hpf[i]==0]*1,[i for i in range(n) if hpf[i]==0])],1)
    idx1=[i for i in range(n) if hpf[i]==0]; idx2=[i for i in range(n) if hpf[i]==1]
    w1=np.triu(D2[np.ix_(idx1,idx1)],1); w1=w1[w1>0]
    w2=np.triu(D2[np.ix_(idx2,idx2)],1); w2=w2[w2>0]
    bt=D2[np.ix_(idx1,idx2)].ravel()
    mb=[np.nanmean([Me[mi[kids[i]]] for i in idx1]),np.nanmean([Me[mi[kids[i]]] for i in idx2])]
    # permutation p for HP PERMANOVA
    rng=np.random.default_rng(1); cnt=1
    for _ in range(199):
        ls=list(hpf); rng.shuffle(ls); Fp=permF(D2,ls)
        if Fp is not None and Fp>=Fhp: cnt+=1
    print(f"  --- HP-supervised 分離度(直接按 HP family 切) ---")
    print(f"    HP1 內距²中位={np.median(w1):.3f}  HP2 內距²中位={np.median(w2):.3f}  HP1-HP2 間距²中位={np.median(bt):.3f}")
    print(f"    HP-PERMANOVA pseudo-F={Fhp:.2f} p={cnt/200:.3f}  | mean-β HP1={mb[0]:.3f} HP2={mb[1]:.3f} Δβ={abs(mb[0]-mb[1]):.3f}")
    print(f"    → {'HP 在距離空間真的分開(間距>內距, PERMANOVA 顯著)=可監督切' if Fhp and cnt/200<0.05 and np.median(bt)>max(np.median(w1),np.median(w2)) else 'HP 距離上沒明顯分開(看起來分但實際未必)'}")
