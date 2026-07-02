#!/usr/bin/env python3
"""診斷三個切群失敗模式是否真實存在（用 29 個 example 的既有 binary 距離矩陣，純讀）。
FM1 單離群吃掉群位: maxclust=k 下某群只有 1-2 read → 現規則 min<3 整個 k 被丟 → 失去其餘乾淨多群。
FM2 maxclust 切碎 vs 距離門檻: maxclust=k 強切 k 群 vs 依樹高 gap 自然切 → 後者群數/群大小是否更乾淨。
FM3 多解析度被折成標籤數: k=2 與 k=3 都乾淨對齊但只輸出一個。
輸出 cluster_failure_modes_probe.json + console。"""
import os, csv, glob, json
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from collections import Counter

WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
OUT=f"{WT}/output/_kprofile_heatmap"
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}

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

def natural_cut_by_gap(Z, n):
    """依 merge 高度最大 gap 自然切（distance 門檻）→ 回 (k_natural, t, sizes)。"""
    h=Z[:,2]  # n-1 merge heights, ascending
    if len(h)<2: return 1,0.0,[n]
    # 從頂端找最大相對 gap（合併高度跳最大處 = 自然群界）
    gaps=np.diff(h); gi=int(np.argmax(gaps))
    t=(h[gi]+h[gi+1])/2.0
    cl=fcluster(Z,t,criterion="distance")
    return len(set(cl)),float(t),sorted(Counter(cl).values(),reverse=True)

dirmap={}
for mp in glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True):
    rd=mp.rsplit("/distance/",1)[0]
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): dirmap[part]=rd

items=json.load(open(f"{A}/kprofile_examples.json"))["items"]
out=[]; fm1=[]; fm2=[]; fm3=[]
for L in items:
    key=f"{L['chrom']}_{L['pos']}"; rd=dirmap.get(key)
    if not rd: continue
    reads={r["read_id"]:r for r in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
    dids,D=loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={r:i for i,r in enumerate(dids)}
    def istum(t): return str(t) in("1","true","True")
    kids=[r for r in dids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP]
    if len(kids)<MIN_SZ*2: continue
    sub=D[np.ix_([di[r] for r in kids],[di[r] for r in kids])]; kp=peel(sub)
    if len(kp)<MIN_SZ*2: continue
    sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
    n=sub.shape[0]
    Z=linkage(squareform(sub,checks=False),method="average")
    rec={"key":key,"group":L["group"],"n":n,"maxclust":{},"big_clusters":{},"accepted_ks":[]}
    # FM1/FM3: maxclust 各 k 的群大小 + 現規則接受與否 + 「大群(≥3)數」
    max_big=0; max_big_k=0
    for k in range(2,min(6,n//1)+1):
        if k>=n: break
        cl=fcluster(Z,k,"maxclust"); sz=sorted(Counter(cl).values(),reverse=True)
        big=sum(1 for s in sz if s>=MIN_SZ); small=sum(1 for s in sz if s<MIN_SZ)
        accepted=(len(sz)>=k and min(sz)>=MIN_SZ)
        rec["maxclust"][k]={"sizes":sz,"big":big,"small_singletons":small,"accepted_current":accepted}
        if accepted: rec["accepted_ks"].append(k)
        if big>max_big: max_big,max_big_k=big,k
    rec["max_big_clusters"]=max_big; rec["max_big_at_k"]=max_big_k
    rec["accepted_max_k"]=max(rec["accepted_ks"]) if rec["accepted_ks"] else 1
    # FM2: 自然 gap 切
    kn,t,szn=natural_cut_by_gap(Z,n)
    rec["natural_gap"]={"k":kn,"t":round(t,4),"sizes":szn,"big":sum(1 for s in szn if s>=MIN_SZ)}
    out.append(rec)
    # 失敗模式判定
    # FM1: 存在某 k 有 ≥3 big 群但被現規則丟(因有 singleton)，且現接受的 max_k 較小
    if max_big>=3 and max_big>rec["accepted_max_k"]:
        fm1.append({"key":key,"max_big_clusters":max_big,"at_k":max_big_k,"accepted_max_k":rec["accepted_max_k"],
                    "maxclust_sizes":rec["maxclust"].get(max_big_k,{}).get("sizes")})
    # FM2: 自然 gap 切的 big 群數 > 現接受 max_k 的 big 群數（自然切看到更多乾淨群）
    acc_big=rec["maxclust"].get(rec["accepted_max_k"],{}).get("big",1) if rec["accepted_ks"] else 1
    if rec["natural_gap"]["big"]>acc_big and rec["natural_gap"]["k"]!=rec["accepted_max_k"]:
        fm2.append({"key":key,"natural_k":kn,"natural_big":rec["natural_gap"]["big"],"natural_sizes":szn,
                    "accepted_max_k":rec["accepted_max_k"],"accepted_big":acc_big})
    # FM3: 多個 accepted k（≥2）= 多解析度被折疊
    if len(rec["accepted_ks"])>=2:
        fm3.append({"key":key,"accepted_ks":rec["accepted_ks"],"meaningful_ks":L.get("meaningful_ks")})

summary={"n_loci":len(out),
    "FM1_singleton_eats_slot":{"n":len(fm1),"examples":fm1[:8]},
    "FM2_maxclust_vs_natural_gap":{"n":len(fm2),"examples":fm2[:8]},
    "FM3_multiresolution_collapsed":{"n":len(fm3),"examples":fm3[:8]},
    "loci":out}
json.dump(summary,open(f"{A}/cluster_failure_modes_probe.json","w"),indent=1)
print(f"=== n={len(out)} 代表位點 ===")
print(f"FM1 單離群吃群位（有≥3大群但現規則只接受更少k）: {len(fm1)}")
for e in fm1[:6]: print(f"   {e['key']}: maxclust k={e['at_k']} sizes={e['maxclust_sizes']} → 現只接受到 k={e['accepted_max_k']}")
print(f"FM2 自然gap切比maxclust看到更多乾淨群: {len(fm2)}")
for e in fm2[:6]: print(f"   {e['key']}: 自然 k={e['natural_k']} sizes={e['natural_sizes']} vs 現 max_k={e['accepted_max_k']}(big={e['accepted_big']})")
print(f"FM3 多解析度(≥2 accepted k 被折成1個輸出): {len(fm3)}")
for e in fm3[:6]: print(f"   {e['key']}: accepted_ks={e['accepted_ks']} meaningful_ks={e['meaningful_ks']}")
