#!/usr/bin/env python3
"""全基因組 cluster_redesign — TP+FP, big7 本機, 平行化。
每 chr: binary 跑 TP + FP VCF → 各 locus 矩陣 → analyze_locus(五類) → 收 {set:TP/FP, fine/coarse conf} → 刪暫存。
平行: Pool over loci(單執行緒 BLAS)。null 減量(Rnull=15)提速。輸出 cluster_redesign_wg.json + summary(5類×TP/FP)。
用法: python3 cluster_redesign_wg.py [CHRS]  e.g. "22"(smoke) 或留空=全 1-22。"""
import os, csv, glob, json, sys, subprocess, shutil, time
os.environ.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
import numpy as np
sys.path.insert(0,os.path.dirname(__file__)); import cluster_redesign as CR
from multiprocessing import Pool

WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; BIN=f"{WT}/build/bin/inter_sub_mod"
TUMOR="/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL="/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF="/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
VCFDIR="/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"; os.makedirs("/big7_disk/liaoyoyo2001/tmp",exist_ok=True)
CR.stab_excess.__defaults__=(15,)   # Rnull 25→15 提速(全基因組)
NPROC=24
_raw=sys.argv[1].split(",") if len(sys.argv)>1 else [str(c) for c in range(1,23)]
CHRS=[c if str(c).startswith("chr") else f"chr{c}" for c in _raw]

def proc_locus(arg):
    rd,setlab=arg
    try:
        reads={x["read_id"]:x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"),delimiter="\t")}
        dids,D=CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di={x:i for i,x in enumerate(dids)}
        rows=open(f"{rd}/methylation/methylation.csv").read().strip().split("\n")
        mi={}; M=[]
        for j,ln in enumerate(rows[1:]):
            q=ln.split(","); mi[q[0]]=j; M.append([np.nan if v in("","NA","nan","NaN") else float(v) for v in q[1:]])
        M=np.array(M)
        it=lambda t:str(t) in("1","true","True")
        ids=[x for x in dids if x in reads and it(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
        if len(ids)<CR.MIN_SZ*2: return None
        sub=D[np.ix_([di[x] for x in ids],[di[x] for x in ids])]; kp=CR.peel(sub)
        if len(kp)<CR.MIN_SZ*2: return None
        ids=[ids[i] for i in kp]; sub=sub[np.ix_(kp,kp)]
        P=np.array([M[mi[x]] for x in ids]); hp=[CR.LABMAP[reads[x]["hp"]] for x in ids]; al=[reads[x]["alt_support"] for x in ids]
        rng=np.random.default_rng(20260622)
        a=CR.analyze_locus(sub,P,hp,al,rng)
        pos=None
        for part in rd.split("/"):
            if part.startswith("chr") and part.count("_")==1: pos=part.split("_")[1]; break
        return {"chrom":rd.split("/")[-1].split("_")[0] if False else None,"pos":pos,"set":setlab,"n":len(ids),
                "coarse_k":a["coarse_k"],"coarse_conf":a["coarse_confidence"],"fine_k":a["fine_k"],
                "fine_conf":a["fine_confidence"],"n_confirmed":a["n_confirmed"],"n_novel":a["n_novel"]}
    except Exception as e:
        return {"err":str(e)[:80],"rd":rd,"set":setlab}

def run():
    out=[]; t0=time.time()
    for chrom in CHRS:
        for setlab,pref in (("TP","filtered_snv_tp"),("FP","filtered_snv_fp")):
            vcf=f"{VCFDIR}/{pref}_{chrom}.vcf.gz"
            if not os.path.exists(vcf): continue
            od=f"{WT}/output/_crwg_{chrom}_{setlab}"; shutil.rmtree(od,ignore_errors=True); os.makedirs(od)
            subprocess.run([BIN,"-t",TUMOR,"-n",NORMAL,"-r",REF,"-v",vcf,"-w","5000","-j","16",
                "--distance-metric","BERNOULLI","--nan-distance-strategy","SKIP","-o",od],
                stdout=open(f"{od}/run.log","w"),stderr=subprocess.STDOUT)
            mats=glob.glob(f"{od}/**/distance/BERNOULLI/matrix.csv",recursive=True)
            args=[(m.rsplit("/distance/",1)[0],setlab) for m in mats]
            with Pool(NPROC) as pool:
                res=[r for r in pool.map(proc_locus,args) if r]
            for r in res:
                if r and "chrom" in r: r["chrom"]=chrom
            out+=[r for r in res if r and "err" not in r]
            errs=[r for r in res if r and "err" in r]
            shutil.rmtree(od,ignore_errors=True)
            print(f"[{chrom}/{setlab}] loci={len(mats)} kept={len([r for r in res if 'err' not in r])} err={len(errs)} elapsed={int(time.time()-t0)}s",flush=True)
            json.dump(out,open(f"{A}/cluster_redesign_wg_records.json","w"))
    # summary: 5類 × TP/FP
    from collections import Counter
    def dist(s): return dict(Counter(r["fine_conf"] for r in out if r["set"]==s))
    summ={"n":len(out),"by_set":{s:sum(1 for r in out if r["set"]==s) for s in("TP","FP")},
          "fine_TP":dist("TP"),"fine_FP":dist("FP"),"elapsed_s":int(time.time()-t0),
          "params":{"Rnull":15,"nproc":NPROC,"chrs":CHRS}}
    json.dump(summ,open(f"{A}/cluster_redesign_wg_summary.json","w"),indent=2)
    print("DONE",json.dumps(summ,ensure_ascii=False))

if __name__=="__main__": run()
