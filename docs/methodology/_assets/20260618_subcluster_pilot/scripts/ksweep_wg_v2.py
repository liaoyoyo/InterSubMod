#!/usr/bin/env python3
"""全基因組 k-sweep dual-mode: 每 region 同時算 (A)tumor-only (B)merged tumor+normal。
一次 binary pass，兩種 read-set 各自分群+三軸 CramerV。輸出兩個 records 檔(同 ksweep_wg_records 格式):
  ksweep_records_tumor.json   (tumor-only, 應重現原 ksweep_wg_records)
  ksweep_records_merged.json  (tumor+normal)
分開輸出便於 ksweep_analyze.py 各跑各的、對比差距。"""
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
    return v,float(p),float(exp.min())

def sweep(M, rows, hp, alle):
    """rows: 該 read-set 在矩陣的 index; hp/alle: 對應標籤。回 per_k 或 None。"""
    if len(rows)<MIN_SZ*2: return None,0
    sub=M[np.ix_(rows,rows)]; keep=peel(sub)
    if len(keep)<MIN_SZ*2: return None,0
    hp=[hp[i] for i in keep]; alle=[alle[i] for i in keep]
    sub=sub[np.ix_(keep,keep)].copy(); np.fill_diagonal(sub,0.0); sub=np.maximum(sub,sub.T)
    n=len(keep)
    try: Z=linkage(squareform(sub,checks=False),method="average")
    except Exception: return None,0
    HP=["HP1" if l in("1","1-1") else "HP2" for l in hp]
    CAR=["G" if l in("1","2") else "C" for l in hp]
    am=[j for j in range(n) if alle[j] in("REF","ALT")]   # allele 只用 REF/ALT(排 UNKNOWN)
    per_k=[]
    for k in range(2,min(5,n//MIN_SZ)+1):
        cl=fcluster(Z,k,criterion="maxclust"); sz=Counter(cl)
        if len(sz)<k or min(sz.values())<MIN_SZ: continue
        try: sil=float(silhouette_score(sub,cl,metric="precomputed"))
        except Exception: continue
        vh,ph,eh=cramerv(HP,list(cl)); vc,pc,ec=cramerv(CAR,list(cl))
        if am: va,pa,ea=cramerv([alle[j] for j in am],[cl[j] for j in am])
        else: va,pa,ea=None,None,None
        per_k.append({"k":k,"sil":round(sil,3),"V_hp":vh,"p_hp":ph,"e_hp":eh,
            "V_carrier":vc,"p_carrier":pc,"e_carrier":ec,"V_allele":va,"p_allele":pa,"e_allele":ea})
    return (per_k if per_k else None),n

for d in glob.glob(f"{WT}/output/_kswv2_tmp_*"): shutil.rmtree(d,ignore_errors=True)
rec_t=[]; rec_m=[]; t0=time.time()
for c in range(1,23):
    chrom=f"chr{c}"; vcf=f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"
    if not os.path.exists(vcf): continue
    out=f"{WT}/output/_kswv2_tmp_{chrom}"; shutil.rmtree(out,ignore_errors=True); os.makedirs(out)
    subprocess.run([BIN,"-t",TUMOR,"-n",NORMAL,"-r",REF,"-v",vcf,"-w","5000","-j","16",
        "--distance-metric","BERNOULLI","--nan-distance-strategy","SKIP","-o",out],
        stdout=open(f"{out}/run.log","w"),stderr=subprocess.STDOUT)
    mats=glob.glob(f"{out}/**/distance/BERNOULLI/matrix.csv",recursive=True)
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
        # read-sets
        rows_t=[]; hp_t=[]; al_t=[]; rows_m=[]; hp_m=[]; al_m=[]
        for i,rid in enumerate(ids):
            r=reads.get(rid)
            if not r or r["hp"] not in LABMAP: continue
            rows_m.append(i); hp_m.append(LABMAP[r["hp"]]); al_m.append(r["alt_support"])
            if istum(r["is_tumor"]):
                rows_t.append(i); hp_t.append(LABMAP[r["hp"]]); al_t.append(r["alt_support"])
        pk_t,nt=sweep(M,rows_t,hp_t,al_t)
        pk_m,nm=sweep(M,rows_m,hp_m,al_m)
        if pk_t: rec_t.append({"chrom":chrom,"pos":pos,"n":nt,"per_k":pk_t})
        if pk_m: rec_m.append({"chrom":chrom,"pos":pos,"n":nm,"per_k":pk_m})
    shutil.rmtree(out,ignore_errors=True)
    json.dump(rec_t,open(f"{ASSET}/ksweep_records_tumor.json","w"))
    json.dump(rec_m,open(f"{ASSET}/ksweep_records_merged.json","w"))
    print(f"[{chrom}] tumor={len(rec_t)} merged={len(rec_m)} elapsed={int(time.time()-t0)}s",flush=True)
json.dump(rec_t,open(f"{ASSET}/ksweep_records_tumor.json","w"))
json.dump(rec_m,open(f"{ASSET}/ksweep_records_merged.json","w"))
print(f"DONE tumor={len(rec_t)} merged={len(rec_m)} elapsed={int(time.time()-t0)}s")
