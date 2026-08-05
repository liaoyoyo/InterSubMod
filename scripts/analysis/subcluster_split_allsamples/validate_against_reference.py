#!/usr/bin/env python3
"""驗證 compute_split_accounting.py 的 Q1 + 分帶邏輯能否重現 HCC1395 pileup-track 參考數字。
讀既有 records_wg2.json + significance_summary.csv,套 driver 同一套邏輯,對比參考表。
參考(用戶貼): A=5997(19.7%) C=8833(29.0%) AnC=1560(5.1%) A-C=4437(14.6%) C-A=7273(23.9%) Jaccard=0.118
split_accounting.json: N=30490 insuff=412 cant=24081 cansplit=5997 clear=2941 mid=765 weak=506 vweak=1249 degen=536
"""
import sys, json, csv
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/subcluster_split_allsamples")
from compute_split_accounting import (cramers_v_and_sig, band_of, tru, ff, MIN_SZ)
from collections import Counter

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
SIG = f"{WT}/output/_wg_bdcprime_verify/significance_summary.csv"

records = json.load(open(f"{A}/records_wg2.json"))
N = len(records)
recs = {f"{r['chrom']}_{r['pos']}": r for r in records}

cansplit = {k for k, r in recs.items() if r["all"]["best_k"] is not None and r["all"]["best_k"] >= 2}
clean = set()
for row in csv.DictReader(open(SIG)):
    key = f"{row['Chr']}_{row['Pos']}"
    if key not in recs:
        continue
    hv = tru(row.get("LabelHPPermanovaValid")); hp = ff(row.get("LabelHPPermanovaP")); hw = tru(row.get("LabelHPDispersionWarn"))
    av = tru(row.get("LabelAllelePermanovaValid")); ap = ff(row.get("LabelAllelePermanovaP")); aw = tru(row.get("LabelAlleleDispersionWarn"))
    if (hv and hp is not None and hp < 0.05 and not hw) or (av and ap is not None and ap < 0.05 and not aw):
        clean.add(key)

A_, C = cansplit, clean
inter = len(A_ & C); uni = len(A_ | C)

def pc(x):
    return round(100 * x / N, 1)

print("===== Q1 Venn (driver 邏輯 vs 參考) =====")
print(f"  N={N}  (參考 30490)")
print(f"  A=cansplit: {len(A_)} ({pc(len(A_))}%)   參考 5997 (19.7%)")
print(f"  C=clean PERMANOVA: {len(C)} ({pc(len(C))}%)   參考 8833 (29.0%)")
print(f"  A∩C: {inter} ({pc(inter)}%)   參考 1560 (5.1%)")
print(f"  A−C: {len(A_-C)} ({pc(len(A_-C))}%)   參考 4437 (14.6%)")
print(f"  C−A: {len(C-A_)} ({pc(len(C-A_))}%)   參考 7273 (23.9%)")
print(f"  Jaccard: {round(inter/uni,3)}   參考 0.118")

# 分帶
insuff = cant = 0; bands = Counter()
for k, r in recs.items():
    a = r["all"]
    if a["n"] < MIN_SZ * 2:
        insuff += 1; continue
    if a["best_k"] is None or a["best_k"] < 2:
        cant += 1; continue
    V, sig, nlab, nclu = cramers_v_and_sig(a["contingency"]) if a["contingency"] else (None, False, 0, 0)
    bands[band_of(V, sig)] += 1

print("\n===== 切割品質分帶 (driver 重建 vs split_accounting.json) =====")
print(f"  insuff(n<6): {insuff}   參考 412")
print(f"  cant(n>=6 但 best_k<2): {cant}   參考 24081")
print(f"  cansplit: {len(A_)}   參考 5997")
print(f"  ├ clear(V>=0.7+sig): {bands['clear']}   參考 2941")
print(f"  ├ mid(0.5-0.7+sig): {bands['mid']}   參考 765")
print(f"  ├ weak(0.3-0.5): {bands['weak']}   參考 506")
print(f"  ├ vweak(<0.3): {bands['vweak']}   參考 1249")
print(f"  └ degen(單標籤/單群): {bands['degen']}   參考 536")
print(f"  分帶加總={bands['clear']+bands['mid']+bands['weak']+bands['vweak']+bands['degen']} (應=cansplit {len(A_)})")
print(f"  總加總={insuff+cant+len(A_)} (應=N {N})")
