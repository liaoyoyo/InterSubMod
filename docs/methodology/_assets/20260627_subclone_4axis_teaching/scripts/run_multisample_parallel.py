#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
多樣本平行 orchestrator(2026-06-29):充分利用 48 核(原 5-way serial 只用 ~7 核)。
4 階段(per-sample barrier 在 merge/region):
 P1 linkage: 全域池 MAXW 跑所有 (sample×chrom) sm_linkage → part_chrN.json
 P2 merge+pair_lists per sample(快)
 P3 multilocus: 全域池跑所有 (sample×chrom) → ml_part_chrN.json
 P4 region_integration+topology+candidate per sample(快)
per-chrom(22-way)讓密集 chrom(長 pole)立即開跑,不卡在 5-way batch。
"""
import os, sys, json, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed

PY = "/bip7_disk/liaoyoyo2001/miniconda3/bin/python3"
SCD = os.path.dirname(os.path.abspath(__file__))
C = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
OUT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
MAXW = int(os.environ.get("MAXW", "40"))
CHROMS = [f"chr{c}" for c in range(1, 23)]
import glob as _glob

def lp(s):
    g = _glob.glob(f"{C}/{s}/paired_full/*complete_matrix/longphase_s")
    return g[0] if g else None

NB = {"H2009": "/big8_disk/data/H2009/ONT/H2009BL.bam",
      "HCC1395_DORADO": "/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam",
      "HCC1937": "/big8_disk/data/HCC1937/ONT/HCC1937BL.bam",
      "HCC1954": "/big8_disk/data/HCC1954/ONT/HCC1954BL.bam"}
SAMPLES = [s for s in ["H2009", "HCC1395_DORADO", "HCC1937", "HCC1954"] if lp(s) and os.path.exists(NB[s])]

def env_for(s):
    e = dict(os.environ)
    e.update({"SM_TBAM": f"{lp(s)}/{s}_tagged.bam", "SM_NBAM": NB[s], "SM_VD": lp(s),
              "SM_VCF_MODE": "single", "SM_CNBED": "", "SM_WORKDIR": f"{OUT}/{s}", "SM_DATA": f"{OUT}/{s}"})
    return e

def run(script, args, env, log):
    with open(log, "w") as f:
        return subprocess.run([PY, f"{SCD}/{script}"] + args, env=env, stdout=f, stderr=subprocess.STDOUT).returncode

def pooled(jobs, label):
    t0 = time.time(); done = 0; fail = 0
    with ThreadPoolExecutor(max_workers=MAXW) as ex:
        futs = {ex.submit(run, *j[:4]): j[4] for j in jobs}
        for fu in as_completed(futs):
            done += 1; rc = fu.result()
            if rc != 0: fail += 1
            if done % 10 == 0 or done == len(jobs):
                print(f"  [{label}] {done}/{len(jobs)} done ({fail} fail) t={int(time.time()-t0)}s", flush=True)
    return fail

for s in SAMPLES:
    os.makedirs(f"{OUT}/{s}", exist_ok=True)
print(f"=== PARALLEL {len(SAMPLES)} samples × {len(CHROMS)} chroms, MAXW={MAXW} ===", flush=True)

# P1 linkage
jobs = [("sm_linkage_genomewide.py", [c, f"{OUT}/{s}/sm_linkage_gw_part_{c}.json"], env_for(s),
         f"{OUT}/{s}/lk_{c}.log", f"{s}:{c}") for s in SAMPLES for c in CHROMS]
print(f"P1 linkage: {len(jobs)} jobs", flush=True); pooled(jobs, "linkage")

# P2 merge + pair_lists per sample
for s in SAMPLES:
    run("sm_merge_parts.py", [], env_for(s), f"{OUT}/{s}/merge.log")
    run("sm_pair_lists.py", [], env_for(s), f"{OUT}/{s}/pair_lists.log")
    print(f"  P2 {s} merged+pair_lists", flush=True)

# P3 multilocus
jobs = [("sm_multilocus_combinations.py", [c, f"{OUT}/{s}/ml_part_{c}.json"], env_for(s),
         f"{OUT}/{s}/ml_{c}.log", f"{s}:{c}") for s in SAMPLES for c in CHROMS]
print(f"P3 multilocus: {len(jobs)} jobs", flush=True); pooled(jobs, "multilocus")

# P4 region + topology + candidate + gene
for s in SAMPLES:
    e = env_for(s)
    for sc in ("sm_region_integration.py", "topology_analysis.py", "candidate_scoring.py", "region_gene_annotation.py"):
        rc = run(sc, [], e, f"{OUT}/{s}/{sc}.log")
        if rc != 0:
            print(f"  P4 {s} {sc} FAIL(rc={rc})", flush=True); break
    else:
        print(f"  P4 {s} DONE", flush=True)
print(f"=== ALL PARALLEL DONE {SAMPLES} ===", flush=True)
