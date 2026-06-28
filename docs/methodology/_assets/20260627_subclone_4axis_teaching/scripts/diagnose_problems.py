#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
診斷 4 個特殊問題並實測「有沒有方法處理」：
  P1 單讀脆弱性 → 實跑兩種對稱化規則(≥2 / binomial 檢定)重分類,量化真正改變比例
  P2 populations vs tree_shape 不一致 → 量化矛盾率
  P4 no_confirmed 是低 power 還是真沒結構 → 按 read 數/n_sSNV 拆解
用法：python3 diagnose_problems.py
"""
import json, os, csv, math
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

def binom_sf(c, N, eps):
    """P(X>=c | Binomial(N,eps))。c<=0→1。"""
    if c <= 0: return 1.0
    if c > N: return 0.0
    return sum(math.comb(N, k) * eps**k * (1-eps)**(N-k) for k in range(c, N+1))

def classify_orig(rr, ra, ar, aa):
    if aa < 2: return "mutual_excl" if (ra+ar) >= 2 else "sparse"
    if ra == 0 and ar == 0: return "co_linked"
    if ar == 0: return "nested_a_in_b"
    if ra == 0: return "nested_b_in_a"
    return "independent"

def classify_flags(aaP, raP, arP, ra, ar):
    """以 present flags 重分類(對稱)。"""
    if not aaP: return "mutual_excl" if (raP or arP) else "sparse"
    if not raP and not arP: return "co_linked"
    if not arP: return "nested_a_in_b"
    if not raP: return "nested_b_in_a"
    return "independent"

# ---- 載入全部 per-pair ----
pairs = []
for c in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl", "independent"):
    for hp in ("sameHP", "diffHP"):
        p = os.path.join(DATA, "lists", f"{c}_{hp}.tsv")
        if os.path.exists(p):
            pairs += list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
def cell(d): return int(d["RR"]), int(d["RA"]), int(d["AR"]), int(d["AA"])

# ---- P1: 重分類 ----
def reclassify(rule):
    trans = defaultdict(Counter); newc = Counter()
    for d in pairs:
        rr, ra, ar, aa = cell(d); N = int(d["coread"]); old = d["config"]
        if rule == "ge2":
            aaP, raP, arP = aa >= 2, ra >= 2, ar >= 2
        else:  # binom eps
            eps = rule
            aaP = binom_sf(aa, N, eps) < 0.05
            raP = binom_sf(ra, N, eps) < 0.05
            arP = binom_sf(ar, N, eps) < 0.05
        nw = classify_flags(aaP, raP, arP, ra, ar)
        trans[old][nw] += 1; newc[nw] += 1
    return trans, newc

oldc = Counter(d["config"] for d in pairs)
tr2, new2 = reclassify("ge2")
trB, newB = reclassify(0.02)
changed2 = sum(v for o in tr2 for n, v in tr2[o].items() if n != o)
changedB = sum(v for o in trB for n, v in trB[o].items() if n != o)

# ---- P2: populations vs tree_shape 矛盾(2-sSNV) ----
RI = json.load(open(os.path.join(DATA, "sm_region_integration.json"), encoding="utf-8"))
regs = RI["regions"]
def has(pops, ch): return any(ch in g for g in pops) if pops else False
p2 = {"empty_pops_2snv": 0, "co_linked_no_AA": 0, "sibling_no_mutualexcl": 0, "total_2snv": 0}
for r in regs:
    if r["n_sSNV"] != 2: continue
    p2["total_2snv"] += 1
    pops = list((r.get("populations") or {}).keys())
    if not pops: p2["empty_pops_2snv"] += 1; continue
    sh = r["tree_shape"]
    if sh == "co_linked_lineage" and "AA" not in pops: p2["co_linked_no_AA"] += 1
    if sh == "sibling_only" and not ("AR" in pops and "RA" in pops): p2["sibling_no_mutualexcl"] += 1

# ---- P4: no_confirmed 是低 power 還是真沒結構 ----
def reads_of(r): return sum((r.get("populations") or {}).values())
nc = [r for r in regs if r["tree_shape"] == "no_confirmed_structure"]
p4 = {"total": len(nc),
      "reads_lt6": sum(1 for r in nc if reads_of(r) < 6),
      "reads_6to19": sum(1 for r in nc if 6 <= reads_of(r) < 20),
      "reads_ge20": sum(1 for r in nc if reads_of(r) >= 20),
      "n_sSNV_2": sum(1 for r in nc if r["n_sSNV"] == 2),
      "n_sSNV_ge5": sum(1 for r in nc if r["n_sSNV"] >= 5)}

out = {
    "P1_reclassification": {
        "orig_config_counts": dict(oldc),
        "rule_ge2_new_counts": dict(new2), "rule_ge2_n_changed": changed2,
        "rule_ge2_pct_changed": round(100*changed2/len(pairs), 1),
        "rule_binom_eps0.02_new_counts": dict(newB), "rule_binom_n_changed": changedB,
        "rule_binom_pct_changed": round(100*changedB/len(pairs), 1),
        "transition_ge2": {o: dict(c) for o, c in tr2.items()},
        "total_pairs": len(pairs),
    },
    "P2_pop_vs_shape": p2,
    "P4_no_confirmed_breakdown": p4,
}
with open(os.path.join(DATA, "problem_diagnosis.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

print(f"=== P1 單讀脆弱性：實跑重分類({len(pairs):,} pairs) ===")
print(f"  原始 config: {dict(oldc)}")
print(f"  對稱≥2 後  : {dict(new2)}  → 改變 {changed2:,} ({round(100*changed2/len(pairs),1)}%)")
print(f"  binom ε=2% : {dict(newB)}  → 改變 {changedB:,} ({round(100*changedB/len(pairs),1)}%)")
print(f"  主要轉移(≥2):")
for o in tr2:
    ch = {n: v for n, v in tr2[o].items() if n != o}
    if ch: print(f"    {o} → {ch}")
print(f"=== P2 populations vs tree_shape 矛盾(2-sSNV n={p2['total_2snv']}) ===")
print(f"  空 populations: {p2['empty_pops_2snv']} ; co_linked 卻無AA: {p2['co_linked_no_AA']} ; sibling 卻非互斥: {p2['sibling_no_mutualexcl']}")
print(f"=== P4 no_confirmed({p4['total']}) 拆解 ===")
print(f"  reads<6: {p4['reads_lt6']} / 6-19: {p4['reads_6to19']} / ≥20: {p4['reads_ge20']}")
print(f"  n_sSNV==2: {p4['n_sSNV_2']} / ≥5: {p4['n_sSNV_ge5']}")
print("OK wrote problem_diagnosis.json")
