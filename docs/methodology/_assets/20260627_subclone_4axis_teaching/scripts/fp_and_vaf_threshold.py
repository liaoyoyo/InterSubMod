#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
(A) 確認弱判定(單讀)是否 FP-source 較高(用 src_a/src_b,可驗證)。
(B) 實作用戶定義:保留最低 1 條,下限 = max(ONT 噪聲期望 coread×ε, 真實期望 coread×VAF 的判讀)。
    核心可驗證量:
      noise_expect = coread × ε   (ONT 錯誤率)
      real_expect  = coread × min_VAF   (若這個 lineage 真存在於此頻率,該有幾條 off-diag)
      observed     = 決定格 read 數
    判讀:observed ≤ noise_expect → 噪聲;observed ≪ real_expect → lineage 該更多卻沒 → 可疑。
(C) 用 cell> coread×ε 當「present」重分類,報 band(ε=1/2/3%)。
用法：python3 fp_and_vaf_threshold.py
"""
import os, csv, json, statistics as st
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

def load(name):
    p = os.path.join(DATA, "lists", name)
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")) if os.path.exists(p) else []

allp = []
for c in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl", "independent"):
    for hp in ("sameHP", "diffHP"):
        allp += load(f"{c}_{hp}.tsv")

def cells(d): return int(d["RR"]), int(d["RA"]), int(d["AR"]), int(d["AA"])
def is_fp(d): return d.get("src_a") == "FP" or d.get("src_b") == "FP"
def both_fp(d): return d.get("src_a") == "FP" and d.get("src_b") == "FP"
def minvaf(d):
    try: return min(float(d["vaf_a"]), float(d["vaf_b"]))
    except: return None
def deciding(d):
    rr, ra, ar, aa = cells(d); cfg = d["config"]
    if cfg == "nested_a_in_b": return "nested", ra
    if cfg == "nested_b_in_a": return "nested", ar
    if cfg == "independent": return "independent", min(ra, ar)
    return cfg, None

# ---- A. FP 比例 ----
dec = [d for d in allp if deciding(d)[0] in ("nested", "independent")]
weak = [d for d in dec if deciding(d)[1] == 1]
strong = [d for d in dec if (deciding(d)[1] or 0) >= 2]
def fp_rate(rows): return round(100*sum(1 for d in rows if is_fp(d))/len(rows), 1) if rows else None
def bothfp_rate(rows): return round(100*sum(1 for d in rows if both_fp(d))/len(rows), 1) if rows else None
A = {
    "all_decided_pairs": len(dec), "all_any_FP_pct": fp_rate(dec), "all_both_FP_pct": bothfp_rate(dec),
    "weak_any_FP_pct": fp_rate(weak), "weak_both_FP_pct": bothfp_rate(weak),
    "strong_any_FP_pct": fp_rate(strong), "strong_both_FP_pct": bothfp_rate(strong),
    "co_linked_any_FP_pct": fp_rate([d for d in allp if d["config"] == "co_linked"]),
}

# ---- B. noise vs real expectation(weak 單讀) ----
B = {"eps_grid": [0.01, 0.02, 0.03], "weak_n": len(weak)}
real_exp = [int(d["coread"]) * (minvaf(d) or 0) for d in weak]
B["weak_real_expect(coread*minVAF)_median"] = round(st.median(real_exp), 2)
B["weak_observed_is_1_but_real_expect_ge3"] = sum(1 for x in real_exp if x >= 3)
B["weak_observed_is_1_but_real_expect_ge5"] = sum(1 for x in real_exp if x >= 5)
for eps in (0.01, 0.02, 0.03):
    # 單讀(=1) ≤ 噪聲期望 → 判噪聲(被拒)
    rejected = sum(1 for d in weak if 1 <= int(d["coread"]) * eps)
    B[f"weak_rejected_as_noise@eps{eps}"] = rejected
    B[f"weak_rejected_pct@eps{eps}"] = round(100*rejected/len(weak), 1)
    B[f"coread_cutoff@eps{eps}(單讀需coread<此值才保留)"] = round(1/eps, 1)

# ---- C. 用 cell>coread×ε 重分類報 band ----
def classify_flags(aaP, raP, arP):
    if not aaP: return "mutual_excl" if (raP or arP) else "sparse"
    if not raP and not arP: return "co_linked"
    if not arP: return "nested_a_in_b"
    if not raP: return "nested_b_in_a"
    return "independent"
orig = Counter(d["config"] for d in allp)
C = {"orig": dict(orig)}
for eps in (0.01, 0.02, 0.03):
    newc = Counter()
    for d in allp:
        rr, ra, ar, aa = cells(d); N = int(d["coread"]); fl = N*eps
        newc[classify_flags(aa > fl, ra > fl, ar > fl)] += 1
    C[f"eps{eps}"] = dict(newc)

out = {"A_FP_enrichment": A, "B_noise_vs_real_expectation": B, "C_band_by_eps": C}
with open(os.path.join(DATA, "fp_and_vaf_threshold.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

print("=== A. FP 比例(可驗證 src) ===")
for k, v in A.items(): print(f"  {k} = {v}")
print("=== B. 單讀 noise vs real 期望 ===")
print(f"  weak n={B['weak_n']}; real_expect(coread×minVAF) 中位={B['weak_real_expect(coread*minVAF)_median']}")
print(f"  單讀但「若真存在該有≥3條」: {B['weak_observed_is_1_but_real_expect_ge3']} ; ≥5條: {B['weak_observed_is_1_but_real_expect_ge5']}")
for eps in (0.01, 0.02, 0.03):
    print(f"  ε={eps}: 單讀被判噪聲 {B[f'weak_rejected_as_noise@eps{eps}']} ({B[f'weak_rejected_pct@eps{eps}']}%); 單讀保留需 coread<{B[f'coread_cutoff@eps{eps}(單讀需coread<此值才保留)']}")
print("=== C. cell>coread×ε 重分類 band ===")
print(f"  orig: {dict(orig)}")
for eps in (0.01, 0.02, 0.03): print(f"  ε={eps}: {C[f'eps{eps}']}")
print("OK wrote fp_and_vaf_threshold.json")
