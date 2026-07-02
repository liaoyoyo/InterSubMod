#!/usr/bin/env python3
"""門檻校準 + 置換 null: 「切不出/不顯著 是真無法處理還是門檻太嚴」+「位點能分幾群」。
每位點 MIN_SZ∈{2,3} × k=2..5: V_real(對齊 a-priori) vs B=30 打亂標籤 null → V>null_95=真顯著。
① MIN_SZ 3→2 sensitivity: 救回的是否真結構(real>null) vs 噪音(~5% FP by construction)。
② 門檻校準: 各 V 門檻下 真實率 vs null 率分離點。
③ per-locus k 曲線(sil + V_real + null_95 band)。chr21+chr22 pilot。"""
import os, glob, csv, json
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter
rng=np.random.default_rng(42)
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
OUT=f"{WT}/output/_thresh_cal_2122"; A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}; B=30
def _f(x):
    x=str(x).strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def loadm(p):
    r=open(p).read().strip().split("\n"); ids=[];M=[]
    for ln in r[1:]:
        q=ln.split(","); ids.append(q[0]); M.append([_f(x) for x in q[1:]])
    return ids,np.array(M)
def peel(s,mn):
    idx=list(range(s.shape[0]))
    while True:
        S=s[np.ix_(idx,idx)];b=(S<0)|np.isnan(S);np.fill_diagonal(b,False)
        if not b.any(): return idx
        idx.remove(idx[int(np.argmax(b.sum(1)))])
        if len(idx)<mn*2: return idx
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

mats=glob.glob(f"{OUT}/**/distance/BERNOULLI/matrix.csv",recursive=True)
recs=[]; all_pairs=[]   # (V_real, sig_by_null) for 校準
for mp in mats:
    rd=mp.rsplit("/distance/",1)[0]; rt=f"{rd}/reads/reads.tsv"
    if not os.path.exists(rt): continue
    pos=None
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): pos=part; break
    reads={r["read_id"]:r for r in csv.DictReader(open(rt),delimiter="\t")}
    ids,M=loadm(mp); di={r:i for i,r in enumerate(ids)}
    def istum(t): return str(t) in("1","true","True")
    base=[r for r in ids if r in reads and istum(reads[r]["is_tumor"]) and reads[r]["hp"] in LABMAP]
    if len(base)<4: continue
    rec={"pos":pos,"by_min":{}}
    for mn in (2,3):
        sub0=M[np.ix_([di[r] for r in base],[di[r] for r in base])]; kp=peel(sub0,mn)
        if len(kp)<mn*2: rec["by_min"][mn]={"splittable":False,"aligned":False,"curve":[]}; continue
        kids=[base[i] for i in kp]; sub=sub0[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T)
        n=len(kids); Z=linkage(squareform(sub,checks=False),method="average")
        lab=[LABMAP[reads[r]["hp"]] for r in kids]
        axes={"hp":["A" if l in("1","1-1") else "B" for l in lab],
              "carrier":["G" if l in("1","2") else "C" for l in lab],
              "allele":[reads[r]["alt_support"] for r in kids]}
        valid=False; aligned=False; curve=[]
        for k in range(2,min(5,n//mn)+1):
            cl=fcluster(Z,k,'maxclust'); sz=Counter(cl)
            if len(sz)<k or min(sz.values())<mn: continue
            valid=True
            try: sil=round(float(silhouette_score(sub,cl,metric="precomputed")),3)
            except: sil=None
            ent={"k":k,"sil":sil}
            for ax,gl in axes.items():
                m=[j for j in range(n) if (ax!="allele" or gl[j] in("REF","ALT"))]
                if len(m)<4: ent[f"V_{ax}"]=None; ent[f"n95_{ax}"]=None; continue
                gg=[gl[j] for j in m]; cc=[cl[j] for j in m]
                vr=cv(gg,cc)
                if vr is None: ent[f"V_{ax}"]=None; ent[f"n95_{ax}"]=None; continue
                nulls=[]
                for _ in range(B):
                    gs=list(gg); rng.shuffle(gs); vn=cv(gs,cc)
                    if vn is not None: nulls.append(vn)
                n95=float(np.percentile(nulls,95)) if nulls else None
                ent[f"V_{ax}"]=round(vr,3); ent[f"n95_{ax}"]=round(n95,3) if n95 is not None else None
                if n95 is not None:
                    sig=vr>n95; all_pairs.append((vr,sig))
                    if sig and vr>=0.3: aligned=True
            curve.append(ent)
        rec["by_min"][mn]={"splittable":valid,"aligned":aligned,"curve":curve,"n":n}
    recs.append(rec)

# ① MIN_SZ sensitivity
def cnt(mn,key): return sum(1 for r in recs if r["by_min"].get(mn,{}).get(key))
N=len(recs)
notsplit3=[r for r in recs if not r["by_min"].get(3,{}).get("splittable")]
rescued=[r for r in notsplit3 if r["by_min"].get(2,{}).get("splittable")]
rescued_aligned=[r for r in rescued if r["by_min"].get(2,{}).get("aligned")]
# ② 門檻校準: 各 V 門檻 真實率
pairs=np.array(all_pairs) if all_pairs else np.zeros((0,2))
cal=[]
for t in [0.1,0.2,0.3,0.4,0.5,0.6,0.7]:
    sel=pairs[pairs[:,0]>=t] if len(pairs) else pairs
    cal.append({"V_th":t,"n_partitions":int(len(sel)),"frac_real":round(float(sel[:,1].mean()),3) if len(sel) else None})
out={"sample":"chr21+chr22 pilot","N":N,
  "split3":cnt(3,"splittable"),"aligned3":cnt(3,"aligned"),
  "split2":cnt(2,"splittable"),"aligned2":cnt(2,"aligned"),
  "notsplit3":len(notsplit3),"rescued_3to2":len(rescued),"rescued_aligned":len(rescued_aligned),
  "rescued_aligned_pct":round(100*len(rescued_aligned)/len(rescued),1) if rescued else None,
  "null_fp_expect":5.0,"calibration":cal,
  "curves":[{"pos":r["pos"],"by_min":{str(m):r["by_min"][m] for m in r["by_min"]}} for r in recs]}
json.dump(out,open(f"{A}/threshold_calibration.json","w"),indent=1)
print(f"N={N}")
print(f"@MIN_SZ=3: splittable={out['split3']} aligned(real>null)={out['aligned3']}")
print(f"@MIN_SZ=2: splittable={out['split2']} aligned={out['aligned2']}")
print(f"切不出@3={out['notsplit3']} → 放寬到 MIN_SZ=2 救回 splittable={out['rescued_3to2']}, 其中 aligned(真結構)={out['rescued_aligned']} ({out['rescued_aligned_pct']}%)")
print(f"  (null FP by construction ~5%; rescued-aligned 若 >>5% = 真被門檻擋, <=5% = 噪音)")
print("門檻校準(V_th → 真實率 frac_real):")
for c in cal: print(f"  V>={c['V_th']}: {c['n_partitions']} 切法, 真實率={c['frac_real']}")
