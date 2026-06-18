#!/usr/bin/env python3
"""全基因組 cluster×label 對應 (1對1/1對多/多對1) — disk-safe 逐chr。
per 位點存: all-read 分群 best_k/sil + cluster×label 列聯表原始counts + within-label 子分群(4標籤)。
分類門檻之後可便宜重算 (records_wg2.json)。"""
import os, csv, glob, json, subprocess, shutil, time
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
from collections import Counter, defaultdict

WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
BIN=f"{WT}/build/bin/inter_sub_mod"
TUMOR="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF="/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
VCFDIR="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
ASSET=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"
MIN_SZ=3
LABMAP={"1":"1","HP1":"1","1-1":"1-1","2":"2","HP2":"2","2-1":"2-1"}

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
def cluster_bestk(sub):
    """UPGMA scan k>=2, min cluster>=MIN_SZ, best silhouette. returns (k,labels_array,sil) or None."""
    n=sub.shape[0]
    if n<MIN_SZ*2: return None
    np.fill_diagonal(sub,0.0); sub=np.maximum(sub,sub.T)
    try: Z=linkage(squareform(sub,checks=False),method="average")
    except Exception: return None
    best=None
    for k in range(2,min(6,n//MIN_SZ)+1):
        lab=fcluster(Z,k,criterion="maxclust"); sz=Counter(lab)
        if len(sz)<k or min(sz.values())<MIN_SZ: continue
        try: s=silhouette_score(sub,lab,metric="precomputed")
        except Exception: continue
        if best is None or s>best[2]: best=(k,lab,float(s))
    return best

# cleanup leftover tmp
for d in glob.glob(f"{WT}/output/_wg_sub_tmp_*")+glob.glob(f"{WT}/output/_wgc_tmp_*"):
    shutil.rmtree(d,ignore_errors=True)

records=[]; t0=time.time()
for c in range(1,23):
    chrom=f"chr{c}"; vcf=f"{VCFDIR}/filtered_snv_tp_{chrom}.vcf.gz"
    if not os.path.exists(vcf): continue
    out=f"{WT}/output/_wgc_tmp_{chrom}"; shutil.rmtree(out,ignore_errors=True); os.makedirs(out)
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
        reads=list(csv.DictReader(open(rt),delimiter="\t"))
        rid2=( {x["read_id"]:x for x in reads} )
        ids,M=load_matrix(mp); id2row={rid:i for i,rid in enumerate(ids)}
        def istum(t): return str(t) in("1","true","True")
        # ----- all-read tumor clustering + contingency -----
        all_rows=[]; all_lab=[]
        for rid in ids:
            if rid in rid2 and istum(rid2[rid]["is_tumor"]) and rid2[rid]["hp"] in LABMAP:
                all_rows.append(id2row[rid]); all_lab.append(LABMAP[rid2[rid]["hp"]])
        allrec={"n":len(all_rows),"n_complete":None,"best_k":None,"best_sil":None,"contingency":None,
                "labels":dict(Counter(all_lab))}
        if len(all_rows)>=MIN_SZ*2:
            sub=M[np.ix_(all_rows,all_rows)]; keep=peel(sub); allrec["n_complete"]=len(keep)
            if len(keep)>=MIN_SZ*2:
                cb=cluster_bestk(sub[np.ix_(keep,keep)].copy())
                if cb:
                    k,labs,sil=cb; allrec["best_k"]=int(k); allrec["best_sil"]=round(sil,3)
                    kept_lab=[all_lab[i] for i in keep]
                    cont=defaultdict(lambda: defaultdict(int))
                    for cl,lb in zip(labs,kept_lab): cont[lb][int(cl)]+=1
                    allrec["contingency"]={lb:dict(d) for lb,d in cont.items()}
        # ----- within-label subcluster (多對1 pieces) -----
        within={}
        for lab in ("1-1","2-1","1","2"):
            rows=[id2row[rid] for rid in ids if rid in rid2 and istum(rid2[rid]["is_tumor"]) and LABMAP.get(rid2[rid]["hp"])==lab]
            w={"n":len(rows),"best_k":None,"best_sil":None,"min_frac":None}
            if len(rows)>=MIN_SZ*2:
                sub=M[np.ix_(rows,rows)]; keep=peel(sub)
                if len(keep)>=MIN_SZ*2:
                    cb=cluster_bestk(sub[np.ix_(keep,keep)].copy())
                    if cb:
                        k,labs,sil=cb; sz=Counter(labs)
                        w["best_k"]=int(k); w["best_sil"]=round(sil,3); w["min_frac"]=round(min(sz.values())/len(keep),3)
            within[lab]=w
        records.append({"chrom":chrom,"pos":pos,"all":allrec,"within":within})
    shutil.rmtree(out,ignore_errors=True)
    json.dump(records,open(f"{ASSET}/records_wg2.json","w"))
    print(f"[{chrom}] regions={len(mats)} +{len(records)-n0} total={len(records)} elapsed={int(time.time()-t0)}s",flush=True)
json.dump(records,open(f"{ASSET}/records_wg2.json","w"))
print(f"DONE total={len(records)} elapsed={int(time.time()-t0)}s")
