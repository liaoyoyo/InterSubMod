#!/usr/bin/env python3
"""smoke test: import decisionflow_wg functions, run process() on real chr21+22 matrices.
驗證 PERMANOVA/PERMDISP 不炸 + record schema 完整。不重跑 binary。"""
import sys, glob, csv, time, os
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
import decisionflow_wg as D
WT=D.WT; LABMAP=D.LABMAP
mats=glob.glob(f"{WT}/output/_thresh_cal_2122/**/distance/BERNOULLI/matrix.csv",recursive=True)
print(f"found {len(mats)} matrices")
n_ok=n_perm=n_split=n_err=0; t0=time.time(); samples=[]
for mp in mats[:120]:
    rd=mp.rsplit("/distance/",1)[0]; rt=f"{rd}/reads/reads.tsv"
    if not os.path.exists(rt): continue
    pos=None
    for part in rd.split("/"):
        if part.count("_")==1 and part.startswith("chr"): chrom,pos=part.split("_"); break
    if pos is None: continue
    reads={r["read_id"]:r for r in csv.DictReader(open(rt),delimiter="\t")}
    ids,M=D.loadm(mp)
    def istum(t): return str(t) in("1","true","True")
    tum=[i for i,rid in enumerate(ids) if rid in reads and istum(reads[rid]["is_tumor"]) and reads[rid]["hp"] in LABMAP]
    mer=[i for i,rid in enumerate(ids) if rid in reads and reads[rid]["hp"] in LABMAP]
    try:
        rt_rec=D.process(M,ids,reads,tum,chrom,pos)
        rm_rec=D.process(M,ids,reads,mer,chrom,pos)
        n_ok+=1
        if rt_rec.get("perm") and any(v for v in rt_rec["perm"].values()): n_perm+=1
        if rt_rec.get("byT",{}).get(4,{}).get("split"): n_split+=1
        if len(samples)<3 and rt_rec.get("n_complete") and rt_rec["n_complete"]>=10:
            samples.append((pos,rt_rec))
    except Exception as e:
        n_err+=1
        if n_err<=3: print(f"  ERR {chrom}:{pos}: {type(e).__name__}: {e}")
print(f"\nprocessed_ok={n_ok} with_perm={n_perm} split@T4={n_split} errors={n_err} elapsed={time.time()-t0:.1f}s")
print("\n=== sample records (n_complete>=10) ===")
import json
for pos,rec in samples:
    print(f"\nchr_:{pos} n={rec['n']} lab_n={rec['lab_n']} n_complete={rec['n_complete']}")
    print("  perm:",json.dumps(rec['perm']))
    print("  byT :",json.dumps(rec['byT']))
