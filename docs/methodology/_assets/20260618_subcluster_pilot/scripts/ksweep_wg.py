#!/usr/bin/env python3
"""全基因組 k-sweep: 每 region 掃 k=2..min(5,n//3), 每 k 記 silhouette + a-priori 軸(HP/CARRIER/ALLELE)
CramerV + chi2 p (raw, 未 gate/未校正 — 校正在 ksweep_analyze.py 便宜後算)。
disk-safe 逐 chr: 跑 binary → sweep → 刪 tmp → 累積 dump。基於 wg_contingency.py + ksweep_pilot.py。
SoT: records → ksweep_wg_records.json (全 22 chr)。單樣本 HCC1395 tumor reads BERNOULLI ±5000。"""
import os, csv, glob, json, subprocess, shutil, time
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.metrics import silhouette_score
from collections import Counter

WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
BIN=f"{WT}/build/bin/inter_sub_mod"
TUMOR="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF="/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
VCFDIR="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
ASSET=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"
MIN_SZ=3; LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}

def _f(x):
    x=x.strip(); return np.nan if x in("","NA","nan","NaN") else float(x)
def load_matrix(path):
    rows=open(path).read().strip().split("\n"); ids=[]; M=[]
    for ln in rows[1:]:
        p=ln.split(","); ids.append(p[0]); M.append([_f(x) for x in p[1:]])
    return ids,np.array(M)
def peel(sub):
    idx=list(range(sub.shape[0]))
    while True:
        S=sub[np.ix_(idx,idx)]; bad=(S<0)|np.isnan(S); np.fill_diagonal(bad,False)
        if not bad.any(): return idx
        idx.remove(idx[int(np.argmax(bad.sum(1)))])
        if len(idx)<MIN_SZ*2: return idx
def cramerv(groups, clusters):
    g=sorted(set(groups)); c=sorted(set(clusters))
    if len(g)<2 or len(c)<2: return None,None,None
    tab=np.zeros((len(g),len(c)))
    gi={x:i for i,x in enumerate(g)}; ci={x:i for i,x in enumerate(c)}
    for a,b in zip(groups,clusters): tab[gi[a],ci[b]]+=1
    tab=tab[tab.sum(1)>0][:,tab.sum(0)>0]
    if tab.shape[0]<2 or tab.shape[1]<2: return None,None,None
    try: chi2,p,dof,exp=chi2_contingency(tab)
    except Exception: return None,None,None
    nn=tab.sum(); v=float(np.sqrt(chi2/(nn*(min(tab.shape)-1)))) if nn>0 else None
    min_exp=float(exp.min())   # Cochran reliability proxy
    return v,float(p),min_exp

for d in glob.glob(f"{WT}/output/_ksw_tmp_*"): shutil.rmtree(d,ignore_errors=True)
records=[]; t0=time.time()
for c in range(1,23):
    chrom=f"chr{c}"; vcf=f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"
    if not os.path.exists(vcf): continue
    out=f"{WT}/output/_ksw_tmp_{chrom}"; shutil.rmtree(out,ignore_errors=True); os.makedirs(out)
    subprocess.run([BIN,"-t",TUMOR,"-n",NORMAL,"-r",REF,"-v",vcf,"-w","5000","-j","16",
        "--distance-metric","BERNOULLI","--nan-distance-strategy","SKIP","-o",out],
        stdout=open(f"{out}/run.log","w"),stderr=subprocess.STDOUT)
    mats=glob.glob(f"{out}/**/distance/BERNOULLI/matrix.csv",recursive=True); n0=len(records)
    for mp in mats:
        rdir=os.path.dirname(os.path.dirname(os.path.dirname(mp)))
        rt=os.path.join(rdir,"reads","reads.tsv")
        if not os.path.exists(rt): continue
        pos=None
        for part in rdir.split("/"):
            if part.startswith(f"{chrom}_") and part.count("_")==1: pos=part.split("_")[1]; break
        reads={r["read_id"]:r for r in csv.DictReader(open(rt),delimiter="\t")}
        ids,M=load_matrix(mp)
        def istum(t): return str(t) in ("1","true","True")
        rows=[]; lab=[]; alle=[]
        for i,rid in enumerate(ids):
            r=reads.get(rid)
            if r and istum(r["is_tumor"]) and r["hp"] in LABMAP:
                rows.append(i); lab.append(LABMAP[r["hp"]]); alle.append(r["alt_support"])
        if len(rows)<MIN_SZ*2: continue
        sub=M[np.ix_(rows,rows)]; keep=peel(sub)
        if len(keep)<MIN_SZ*2: continue
        lab=[lab[i] for i in keep]; alle=[alle[i] for i in keep]
        sub=sub[np.ix_(keep,keep)].copy(); np.fill_diagonal(sub,0.0); sub=np.maximum(sub,sub.T)
        n=len(keep)
        try: Z=linkage(squareform(sub,checks=False),method="average")
        except Exception: continue
        hp=["HP1" if l in ("1","1-1") else "HP2" for l in lab]
        carrier=["G" if l in ("1","2") else "C" for l in lab]
        per_k=[]
        for k in range(2,min(5,n//MIN_SZ)+1):
            cl=fcluster(Z,k,criterion="maxclust"); sz=Counter(cl)
            if len(sz)<k or min(sz.values())<MIN_SZ: continue
            try: sil=float(silhouette_score(sub,cl,metric="precomputed"))
            except Exception: continue
            vh,ph,eh=cramerv(hp,list(cl)); vc,pc,ec=cramerv(carrier,list(cl)); va,pa,ea=cramerv(alle,list(cl))
            per_k.append({"k":k,"sil":round(sil,3),
                "V_hp":vh,"p_hp":ph,"e_hp":eh,"V_carrier":vc,"p_carrier":pc,"e_carrier":ec,
                "V_allele":va,"p_allele":pa,"e_allele":ea})
        if not per_k: continue
        records.append({"chrom":chrom,"pos":pos,"n":n,"n_k":len(per_k),"per_k":per_k})
    shutil.rmtree(out,ignore_errors=True)
    json.dump(records,open(f"{ASSET}/ksweep_wg_records.json","w"))
    print(f"[{chrom}] regions={len(mats)} +{len(records)-n0} total={len(records)} elapsed={int(time.time()-t0)}s",flush=True)
json.dump(records,open(f"{ASSET}/ksweep_wg_records.json","w"))
print(f"DONE total={len(records)} elapsed={int(time.time()-t0)}s")
