#!/usr/bin/env python3
"""Q1: PERMANOVA clean-location(扣 dispersion) vs Fisher+V(cansplit) Venn。
Q2: 實作 Python PERMANOVA(Anderson2001 pseudo-F + 置換), 跑 2組(HP-family) vs 4組(HP-fine) 證明能處理多組。"""
import json, csv, glob, os
import numpy as np
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
SIG=f"{WT}/output/_wg_bdcprime_verify/significance_summary.csv"
recs={f"{r['chrom']}_{r['pos']}":r for r in json.load(open(f"{A}/records_wg2.json"))}
N=len(recs)
def tru(v): return str(v).strip().lower() in("true","1","yes")
def ff(v):
    try: return float(v)
    except: return None

# ===== Q1: clean-location Venn =====
cansplit={k for k,r in recs.items() if r['all']['best_k'] is not None and r['all']['best_k']>=2}
clean=set()  # 任一軸 PERMANOVA sig & 非 dispersion
for row in csv.DictReader(open(SIG)):
    key=f"{row['Chr']}_{row['Pos']}"
    if key not in recs: continue
    hv=tru(row.get('LabelHPPermanovaValid')); hp=ff(row.get('LabelHPPermanovaP')); hw=tru(row.get('LabelHPDispersionWarn'))
    av=tru(row.get('LabelAllelePermanovaValid')); ap=ff(row.get('LabelAllelePermanovaP')); aw=tru(row.get('LabelAlleleDispersionWarn'))
    hp_clean=hv and hp is not None and hp<0.05 and not hw
    al_clean=av and ap is not None and ap<0.05 and not aw
    if hp_clean or al_clean: clean.add(key)
A_=cansplit; C=clean
print("===== Q1: Fisher+V(cansplit) vs PERMANOVA clean-location(扣dispersion) Venn =====")
print(f"  A=Fisher+V cansplit: {len(A_)} ({100*len(A_)/N:.1f}%)")
print(f"  C=PERMANOVA clean-location(任一軸): {len(C)} ({100*len(C)/N:.1f}%)")
print(f"  A∩C: {len(A_&C)} ({100*len(A_&C)/N:.1f}%)")
print(f"  A−C(能切群但clean-location不顯著): {len(A_-C)} ({100*len(A_-C)/N:.1f}%)")
print(f"  C−A(clean-location有但切不出群): {len(C-A_)} ({100*len(C-A_)/N:.1f}%)")
print(f"  Jaccard(A,C)={len(A_&C)/len(A_|C):.3f}")

# ===== Q2: 4-group HP-fine PERMANOVA (Python 實作) =====
LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}
FINE={"1":0,"1-1":1,"2":2,"2-1":3}; FAM={"1":0,"1-1":0,"2":1,"2-1":1}
rng=np.random.default_rng(7)
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
        if len(idx)<4: return idx
def pseudoF(D2,lab):
    n=len(lab); gs=set(lab); K=len(gs)
    if K<2: return None
    SS_T=np.triu(D2,1).sum()/n
    SS_W=0.0
    for g in gs:
        idx=[i for i in range(n) if lab[i]==g]; ng=len(idx)
        if ng<2: continue
        SS_W+=np.triu(D2[np.ix_(idx,idx)],1).sum()/ng
    SS_A=SS_T-SS_W; dfb=K-1; dfw=n-K
    if dfw<=0 or SS_W<=0: return None
    return (SS_A/dfb)/(SS_W/dfw)
def permanova(D2,lab,B=199,ming=3):
    sz={}
    for l in lab: sz[l]=sz.get(l,0)+1
    if sum(1 for v in sz.values() if v>=ming)<2: return None  # 需≥2 組各≥ming
    F=pseudoF(D2,lab)
    if F is None: return None
    cnt=1
    for _ in range(B):
        ls=list(lab); rng.shuffle(ls); Fp=pseudoF(D2,ls)
        if Fp is not None and Fp>=F: cnt+=1
    return F,cnt/(B+1),len([v for v in sz.values() if v>=ming])

mats=glob.glob(f"{WT}/output/_thresh_cal_2122/**/distance/BERNOULLI/matrix.csv",recursive=True)
res={"fam2":{"testable":0,"sig":0},"fine4":{"testable":0,"sig":0,"ngroups":{}}}
ex=[]
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
    kids=[base[i] for i in kp]; sub=sub[np.ix_(kp,kp)].copy(); np.fill_diagonal(sub,0); sub=np.maximum(sub,sub.T); D2=sub**2
    hp=[LABMAP[reads[r]["hp"]] for r in kids]
    fam=[FAM[h] for h in hp]; fine=[FINE[h] for h in hp]
    r2=permanova(D2,fam); r4=permanova(D2,fine)
    if r2: res["fam2"]["testable"]+=1; res["fam2"]["sig"]+= (r2[1]<0.05)
    if r4:
        res["fine4"]["testable"]+=1; res["fine4"]["sig"]+= (r4[1]<0.05)
        res["fine4"]["ngroups"][r4[2]]=res["fine4"]["ngroups"].get(r4[2],0)+1
        if len(ex)<3 and r4[2]>=3: ex.append((rd.split("/")[-1],r4[2],round(r4[0],2),round(r4[1],3)))
print("\n===== Q2: 2組(HP-family) vs 4組(HP-fine) PERMANOVA (chr21+chr22, B=199) =====")
print(f"  2組 HP-family: testable={res['fam2']['testable']}, sig(p<.05)={res['fam2']['sig']}")
print(f"  4組 HP-fine:   testable={res['fine4']['testable']}, sig={res['fine4']['sig']}")
print(f"  4組 testable 中實際有效組數分佈(≥3 reads): {dict(sorted(res['fine4']['ngroups'].items()))}")
print(f"  → 4組 PERMANOVA 能跑(omnibus,任意K); 但需各 sub-tag ≥3 reads, 1-1/2-1 少→多數仍只 2-3 有效組")
print(f"  範例(實際≥3組): {ex}")
json.dump({"Q1":{"cansplit":len(A_),"clean":len(C),"AnC":len(A_&C),"A_only":len(A_-C),"C_only":len(C-A_),"N":N},
           "Q2":res},open(f"{A}/permanova_clean_4group.json","w"),indent=1)
