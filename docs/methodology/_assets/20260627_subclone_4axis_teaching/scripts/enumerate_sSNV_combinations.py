#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sSNV 關聯結果完整列舉（§13-A 從凍結源 derive）：
  A. 每區 sSNV 數分布（1/2/3.../bucket）
  B. tree_shape（區域級組合）分布 + n_sSNV × shape 交叉
  C. n_populations 分布
  D. 基因型字串組合（2-sSNV / 3-sSNV 的 populations 實際組合 + 各組合區域數）
  E. more-reads：每區 read 數分布
  F. classify() 不對稱容忍量化（從 per-pair lists 的 RR/RA/AR/AA cell）
輸出 grep-able JSON + 印表。
用法：python3 enumerate_sSNV_combinations.py
"""
import json, os, csv
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
RI = json.load(open(os.path.join(DATA, "sm_region_integration.json"), encoding="utf-8"))
regs = RI["regions"]
N = len(regs)

def bucket_n(n):
    if n <= 4: return str(n)
    if n <= 9: return "5-9"
    if n <= 49: return "10-49"
    return "50+"

# A. n_sSNV 分布
by_nsnv = Counter(r["n_sSNV"] for r in regs)
nsnv_bucket = Counter(bucket_n(r["n_sSNV"]) for r in regs)
# B. tree_shape + 交叉
by_shape = Counter(r["tree_shape"] for r in regs)
shape_x_nsnv = defaultdict(Counter)
for r in regs:
    shape_x_nsnv[r["tree_shape"]][bucket_n(r["n_sSNV"])] += 1
# C. n_populations
by_npop = Counter(r["n_populations"] for r in regs)
# D. 基因型組合（populations 的 key 集合）
def combo_sig(r):
    pops = r.get("populations", {}) or {}
    return "|".join(sorted(pops.keys()))
combo2 = Counter(combo_sig(r) for r in regs if r["n_sSNV"] == 2)
combo3 = Counter(combo_sig(r) for r in regs if r["n_sSNV"] == 3)
# E. reads/region（populations 值總和）
def reads_of(r):
    return sum((r.get("populations", {}) or {}).values())
reads = [reads_of(r) for r in regs]
reads_sorted = sorted(reads)
def pct(p): return reads_sorted[min(len(reads_sorted)-1, int(p*len(reads_sorted)))]
reads_summary = {"min": min(reads), "p25": pct(.25), "median": pct(.5), "p75": pct(.75), "max": max(reads),
                 "ge20_regions": sum(1 for x in reads if x >= 20), "ge50_regions": sum(1 for x in reads if x >= 50)}

# F. 不對稱量化（per-pair lists）
def read_pairs(name):
    p = os.path.join(DATA, "lists", name)
    if not os.path.exists(p): return []
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
def ic(d, k): return int(d[k])
asym = {}
co = read_pairs("co_linked_sameHP.tsv") + read_pairs("co_linked_diffHP.tsv")
na = read_pairs("nested_a_in_b_sameHP.tsv") + read_pairs("nested_a_in_b_diffHP.tsv")
nb = read_pairs("nested_b_in_a_sameHP.tsv") + read_pairs("nested_b_in_a_diffHP.tsv")
me = read_pairs("mutual_excl_sameHP.tsv") + read_pairs("mutual_excl_diffHP.tsv")
ind = read_pairs("independent_sameHP.tsv") + read_pairs("independent_diffHP.tsv")
# nested_a_in_b: off-diagonal present = RA(>0), zero = AR. 若 RA==1 → 可能本是 co_linked(被1雜訊降級)
na_ra1 = sum(1 for d in na if ic(d, "RA") == 1)
nb_ar1 = sum(1 for d in nb if ic(d, "AR") == 1)
# independent: 兩 off-diagonal 皆>0；min(RA,AR)==1 → 可能本是 nested(被1雜訊推成 independent)
ind_minoff1 = sum(1 for d in ind if min(ic(d, "RA"), ic(d, "AR")) == 1)
# mutual_excl: aa<2 容忍；AA==1 = 既有不對稱容忍受益數
me_aa1 = sum(1 for d in me if ic(d, "AA") == 1)
nested_tot = len(na) + len(nb)
asym = {
    "co_linked_total": len(co), "nested_total": nested_tot, "independent_total": len(ind),
    "mutual_excl_total": len(me),
    "nested_with_offdiag_eq1(could_be_colinked)": na_ra1 + nb_ar1,
    "nested_offdiag1_pct_of_nested": round(100 * (na_ra1 + nb_ar1) / nested_tot, 2) if nested_tot else None,
    "independent_with_minoffdiag_eq1(could_be_nested)": ind_minoff1,
    "independent_minoff1_pct_of_independent": round(100 * ind_minoff1 / len(ind), 2) if ind else None,
    "mutual_excl_with_AA_eq1(existing_tolerance_benefit)": me_aa1,
    "mutual_excl_AA1_pct": round(100 * me_aa1 / len(me), 2) if me else None,
}

out = {
    "source": "sm_region_integration.json + lists/*.tsv (frozen @ branch feat/summary-nreadsvalid 5308d9e)",
    "total_regions": N,
    "A_by_n_sSNV_exact_le10": {str(k): by_nsnv[k] for k in sorted(by_nsnv) if k <= 10},
    "A_by_n_sSNV_bucket": dict(nsnv_bucket),
    "B_by_tree_shape": dict(by_shape),
    "B_shape_x_nSSNV_bucket": {s: dict(c) for s, c in shape_x_nsnv.items()},
    "C_by_n_populations": {str(k): by_npop[k] for k in sorted(by_npop, key=lambda x: (x is None, x))},
    "D_combos_2sSNV(populations_keyset -> regions)": dict(combo2.most_common()),
    "D_combos_3sSNV_top15": dict(combo3.most_common(15)),
    "D_n_distinct_3sSNV_combos": len(combo3),
    "E_reads_per_region": reads_summary,
    "F_classify_asymmetry": asym,
}
with open(os.path.join(DATA, "sSNV_combination_enumeration.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

# 驗證 + 印表
assert N == 7143
assert sum(by_shape.values()) == N
print(f"=== total regions = {N} ===")
print(f"A. n_sSNV bucket: {dict(nsnv_bucket)}  (exact 2={by_nsnv[2]} 3={by_nsnv[3]} 4={by_nsnv[4]})")
print(f"   n_sSNV==1 regions: {by_nsnv.get(1,0)}  (1-sSNV 不可連鎖)")
print(f"B. tree_shape: {dict(by_shape)}")
print(f"C. n_populations: {dict(sorted(by_npop.items(), key=lambda kv: (kv[0] is None, kv[0])))}")
print(f"D. 2-sSNV 組合(populations keyset): {dict(combo2.most_common(8))}")
print(f"   3-sSNV distinct 組合數={len(combo3)}; top: {dict(combo3.most_common(6))}")
print(f"E. reads/region: {reads_summary}")
print(f"F. 不對稱量化:")
for k, v in asym.items(): print(f"     {k} = {v}")
print("OK wrote sSNV_combination_enumeration.json")
