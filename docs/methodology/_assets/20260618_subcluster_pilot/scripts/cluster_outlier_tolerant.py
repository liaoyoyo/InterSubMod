#!/usr/bin/env python3
"""outlier-tolerant 切群規則驗證(pilot chr21+22): 舊規則「任一群<3 即拒」vs 新規則「≥2 valid 群(≥3)即可分,
小群當 outlier 記錄不為準, outlier 比例 <=20%」。置換 null 驗救回的是真結構(對齊 a-priori)還是雜訊。"""
import glob, csv, os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter
rng=np.random.default_rng(11)
MIN_SZ=3; MAXOUT=0.20; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
OUT=f"{WT}/output/_thresh_cal_2122"; A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
def _f(x):
    x=str(x).strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p):
    r=open(p).read().strip().split("\n"); ids=[];M=[]
    for ln in r[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return ids,np.array(M)
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
def old_bestk(sub,Z,n):
    for k in range(2,min(6,n//MIN_SZ)+1):
        cl=fcluster(Z,k,'maxclust');sz=Counter(cl)
        if len(sz)<k or min(sz.values())<MIN_SZ: continue
        return True
    return False
def new_bestk(sub,Z,n):
    """回 (splittable, n_valid_groups, n_outliers, valid_labels_or_None)"""
    best=None
    for k in range(2,min(8,n)+1):  # 允許更高 k(outlier 消耗名額)
        cl=fcluster(Z,k,'maxclust');sz=Counter(cl)
        valid=[c for c,s in sz.items() if s>=MIN_SZ]
        if len(valid)<2: continue
        nout=sum(s for c,s in sz.items() if s<MIN_SZ)
        if nout/n>MAXOUT: continue
        mask=[i for i in range(n) if cl[i] in valid]
        if len(set(cl[i] for i in mask))<2: continue
        try: sil=silhouette_score(sub[np.ix_(mask,mask)],[cl[i] for i in mask],metric="precomputed")
        except: continue
        if best is None or sil>best[0]: best=(sil,len(valid),nout,cl,mask)
    if best is None: return (False,0,0,None,None)
    return (True,best[1],best[2],best[3],best[4])
def aligned(sub,kids,reads,cl,mask):
    """valid 群是否對齊任一 a-priori 軸(V_real>null_95 & V>=.3)"""
    sm=sub[np.ix_(mask,mask)]; clm=[cl[i] for i in mask]
    lab=[LABMAP[reads[kids[i]]["hp"]] for i in mask]
    axes={"hp":["A" if l in("1","1-1") else "B" for l in lab],
          "carrier":["G" if l in("1","2") else "C" for l in lab],
          "allele":[reads[kids[i]]["alt_support"] for i in mask]}
    for ax,gl in axes.items():
        idx=[j for j in range(len(mask)) if (ax!="allele" or gl[j] in("REF","ALT"))]
        if len(idx)<6: continue
        gg=[gl[j] for j in idx]; cc=[clm[j] for j in idx]
        vr=cv(gg,cc)
        if vr is None or vr<0.3: continue
        nulls=[cv(list(np.random.default_rng(s).permutation(gg)),cc) for s in range(30)]
        nulls=[x for x in nulls if x is not None]
        if nulls and vr>np.percentile(nulls,95): return True,ax,round(vr,3)
    return False,None,None

mats=glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True)
N=0; old_split=0; new_split=0; rescued=0; resc_aligned=0; outfrac=[]; ex=[]
for mp in mats:
    rd=mp.rsplit("/distance/",1)[0]; rt=f"{rd}/reads/reads.tsv"
    if not os.path.exists(rt): continue
    reads={r["read_id"]:r for r in csv.DictReader(open(rt),delimiter="\t")}
    ids,M=loadm(mp); di={r:i for i,r in enumerate(ids)}
    def istum(t): return str(t) in("1","true","True")
    base=[r for r in ids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP]
    if len(base)<6: continue
    sub=M[np.ix_([di[r] for r in base],[di[r] for r in base])]; kp=peel(sub)
    if len(kp)<6: continue
    kids=[base[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=len(kids); Z=linkage(squareform(sub,checks=False),method="average")
    N+=1
    o=old_bestk(sub,Z,n); s,nv,nout,cl,mask=new_bestk(sub,Z,n)
    old_split+=o; new_split+=s
    if (not o) and s:
        rescued+=1; outfrac.append(nout/n)
        al,ax,vr=aligned(sub,kids,reads,cl,mask)
        if al:
            resc_aligned+=1
            if len(ex)<5: ex.append((rd.split('/')[-1],nv,nout,ax,vr))
print(f"=== outlier-tolerant 規則 pilot (chr21+22, N={N} 可聚類位點) ===")
print(f"舊規則 splittable: {old_split} ({100*old_split/N:.1f}%)")
print(f"新規則 splittable: {new_split} ({100*new_split/N:.1f}%)")
print(f"新規則救回(舊None→新可分): {rescued}; 其中對齊 a-priori(真結構,V>null): {resc_aligned} ({100*resc_aligned/rescued:.1f}% of rescued)" if rescued else "無救回")
print(f"  救回的 outlier 比例 中位={np.median(outfrac):.3f} (max 限 {MAXOUT})" if outfrac else "")
print(f"  範例(救回且對齊): {ex}")
print(f"\n→ 救回對齊率 >>5%(null) = 新規則救回的是真結構非雜訊" if rescued and resc_aligned/rescued>0.3 else "")
