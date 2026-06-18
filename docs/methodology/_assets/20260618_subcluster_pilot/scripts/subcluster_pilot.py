#!/usr/bin/env python3
"""描述性 pilot: 在 carrier-only(1-1/2-1) + HP-level 子集上偵測「可能有子分群」。
標準: min_sz>=3, 掃 k>=2 取最佳 silhouette(=群內聚集 vs 群外差異)。
輸出: 覆蓋漏斗 + silhouette 分布 + 各門檻通過數 + best-k 分布 + 與 VC 交叉。
口徑: 偵測「候選 apparent 子分群」, 非驗證子克隆 (含 epiallele/技術假象)。"""
import os, csv, glob, json
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
from collections import Counter, defaultdict

ROOT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_chr1_subcluster_pilot"
SIG="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/output/_chr1_strength_verify/significance_summary.csv"
MIN_SZ=3
NOISE={"Noise_Uniform","Noise_Uncorrelated","Noise_Chaotic"}

# VC lookup by Pos
vc_by_pos={}
for r in csv.DictReader(open(SIG)):
    vc_by_pos[r["Pos"]]=r["VerificationClass"]

def _f(x):
    x=x.strip()
    if x in ("","NA","nan","NaN"): return np.nan
    try: return float(x)
    except: return np.nan

def load_matrix(path):
    rows=open(path).read().strip().split("\n")
    ids=[]; M=[]
    for ln in rows[1:]:
        p=ln.split(","); ids.append(p[0]); M.append([_f(x) for x in p[1:]])
    return ids, np.array(M)

def peel_complete(sub):
    """drop reads with -1/NaN (invalid pair) greedily until complete submatrix."""
    idx=list(range(sub.shape[0]))
    while True:
        S=sub[np.ix_(idx,idx)]
        bad=(S<0)|np.isnan(S)
        np.fill_diagonal(bad,False)
        if not bad.any(): return idx
        cnt=bad.sum(axis=1)
        drop=idx[int(np.argmax(cnt))]
        idx.remove(drop)
        if len(idx)<MIN_SZ*2: return idx

def best_split(sub):
    """UPGMA, scan k=2..min(5,n//3), require all clusters>=MIN_SZ, pick best silhouette.
    returns (best_k, best_sil, min_frac) or None."""
    n=sub.shape[0]
    if n<MIN_SZ*2: return None
    np.fill_diagonal(sub,0.0)
    sub=np.maximum(sub,sub.T)
    try: Z=linkage(squareform(sub,checks=False),method="average")
    except Exception: return None
    best=None
    for k in range(2,min(5,n//MIN_SZ)+1):
        lab=fcluster(Z,k,criterion="maxclust")
        sizes=Counter(lab)
        if len(sizes)<k: continue
        if min(sizes.values())<MIN_SZ: continue
        try: s=silhouette_score(sub,lab,metric="precomputed")
        except Exception: continue
        mf=min(sizes.values())/n
        if best is None or s>best[1]: best=(k,round(float(s),3),round(mf,3))
    return best

# gather region dirs
mats=glob.glob(f"{ROOT}/**/distance/BERNOULLI/matrix.csv",recursive=True)
print(f"region matrices found: {len(mats)}")

records=[]
for mp in mats:
    rdir=os.path.dirname(os.path.dirname(os.path.dirname(mp)))
    rt=os.path.join(rdir,"reads","reads.tsv")
    if not os.path.exists(rt): continue
    # pos from path chr1_POS
    pos=None
    for part in rdir.split("/"):
        if part.startswith("chr1_") and part.count("_")==1:
            pos=part.split("_")[1]; break
    reads=list(csv.DictReader(open(rt),delimiter="\t"))
    rid2hp={r["read_id"]:(r["hp"],r["is_tumor"]) for r in reads}
    ids,M=load_matrix(mp)
    id2row={rid:i for i,rid in enumerate(ids)}
    def subset(pred):
        return [id2row[rid] for rid in ids if rid in rid2hp and pred(*rid2hp[rid]) and rid in id2row]
    def istum(t): return str(t) in ("1","true","True")
    subsets={
        "carrier_1-1": subset(lambda h,t: h=="1-1" and istum(t)),
        "carrier_2-1": subset(lambda h,t: h=="2-1" and istum(t)),
        "HP1_germ+carrier": subset(lambda h,t: h in("1","HP1","1-1") and istum(t)),
        "HP2_germ+carrier": subset(lambda h,t: h in("2","HP2","2-1") and istum(t)),
    }
    vc=vc_by_pos.get(pos,"?")
    for name,rows in subsets.items():
        n_raw=len(rows)
        rec={"pos":pos,"subset":name,"vc":vc,"n_raw":n_raw,"eligible":n_raw>=MIN_SZ*2,
             "n_complete":None,"best_k":None,"best_sil":None,"min_frac":None}
        if n_raw>=MIN_SZ*2:
            sub=M[np.ix_(rows,rows)]
            keep=peel_complete(sub)
            rec["n_complete"]=len(keep)
            if len(keep)>=MIN_SZ*2:
                bs=best_split(sub[np.ix_(keep,keep)].copy())
                if bs: rec["best_k"],rec["best_sil"],rec["min_frac"]=bs
        records.append(rec)

# ===== aggregate =====
OUT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
os.makedirs(OUT,exist_ok=True)
json.dump(records,open(f"{OUT}/records.json","w"))

def report(name):
    sub=[r for r in records if r["subset"]==name]
    tot=len(sub); elig=[r for r in sub if r["eligible"]]
    analyzed=[r for r in sub if r["best_sil"] is not None]
    print(f"\n=== {name} ===")
    print(f"  標籤總數: {tot}")
    print(f"  夠 read 可問(n>=6): {len(elig)} ({round(100*len(elig)/tot,1)}%)")
    print(f"  實際算出 best split: {len(analyzed)} ({round(100*len(analyzed)/tot,1)}% of all; {round(100*len(analyzed)/len(elig),1) if elig else 0}% of eligible)")
    if analyzed:
        sils=sorted(r["best_sil"] for r in analyzed)
        n=len(sils)
        print(f"  silhouette 分布: min={sils[0]} q25={sils[n//4]} med={sils[n//2]} q75={sils[3*n//4]} max={sils[-1]}")
        for th in (0.4,0.5,0.6,0.7):
            p=sum(1 for s in sils if s>=th)
            print(f"    >= {th}: {p} ({round(100*p/len(analyzed),1)}% of analyzed, {round(100*p/tot,1)}% of all 標籤)")
        kd=Counter(r["best_k"] for r in analyzed)
        print(f"  best-k 分布: {dict(sorted(kd.items()))}")
        # cross VC: 候選(sil>=0.5) 落在 noise vs structure?
        cand=[r for r in analyzed if r["best_sil"]>=0.5]
        if cand:
            nv=sum(1 for r in cand if r["vc"] in NOISE); sv=sum(1 for r in cand if r["vc"] not in NOISE and r["vc"]!="?")
            print(f"  候選(sil>=0.5) n={len(cand)}: 落 Noise-VC {nv} / 落 結構-VC {sv}")

for nm in ("carrier_1-1","carrier_2-1","HP1_germ+carrier","HP2_germ+carrier"):
    report(nm)
print(f"\n[records -> {OUT}/records.json]")
