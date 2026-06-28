#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
比較「絕對下限 vs VAF 調整 vs ONT-ε」三類門檻,用 FP 比例當可驗證裁判:
  好門檻 = 丟掉的那堆 FP 比例高(抓到噪聲)、留下的那堆 FP 比例低。
  太嚴 = 連 FP-poor(像 TP)的也丟;太鬆 = FP-rich 的還留著。
對「弱判定(單讀)」評各規則的 reject/keep 與其 both-FP 比例,並報整體結構 band。
用法：python3 compare_thresholds.py
"""
import os, csv, json
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
def both_fp(d): return d.get("src_a") == "FP" and d.get("src_b") == "FP"
def minvaf(d):
    try: return min(float(d["vaf_a"]), float(d["vaf_b"]))
    except: return 0.0
def deciding(d):
    rr, ra, ar, aa = cells(d); cfg = d["config"]
    if cfg == "nested_a_in_b": return "nested", ra
    if cfg == "nested_b_in_a": return "nested", ar
    if cfg == "independent": return "independent", min(ra, ar)
    return cfg, None

weak = [d for d in allp if deciding(d)[0] in ("nested", "independent") and deciding(d)[1] == 1]

# 各規則對「單讀」的 reject 判斷
def rules(d):
    N = int(d["coread"]); mv = minvaf(d)
    return {
        "絕對≥2(flat)": True,                         # 單讀一律 reject
        "ε=1%(coread×.01)": 1 <= N*0.01,
        "ε=2%(coread×.02)": 1 <= N*0.02,
        "ε=3%(coread×.03)": 1 <= N*0.03,
        "VAF期望≥3則拒": N*mv >= 3,                  # 若真 lineage 該≥3條卻只1條
        "VAF期望≥5則拒": N*mv >= 5,
        "組合(ε2% 或 VAF≥3)": (1 <= N*0.02) or (N*mv >= 3),
    }

names = list(rules(weak[0]).keys())
calib = {}
for nm in names:
    rej = [d for d in weak if rules(d)[nm]]
    keep = [d for d in weak if not rules(d)[nm]]
    def fp(rows): return round(100*sum(1 for d in rows if both_fp(d))/len(rows), 1) if rows else None
    calib[nm] = {
        "reject_n": len(rej), "reject_pct": round(100*len(rej)/len(weak), 1), "reject_bothFP%": fp(rej),
        "keep_n": len(keep), "keep_bothFP%": fp(keep),
        "separation(rejFP-keepFP)": round((fp(rej) or 0) - (fp(keep) or 0), 1) if keep and rej else None,
    }

# 整體結構 band(cell-present 規則套全 pairs)
def classify_flags(aaP, raP, arP):
    if not aaP: return "mutual_excl" if (raP or arP) else "sparse"
    if not raP and not arP: return "co_linked"
    if not arP: return "nested_a_in_b"
    if not raP: return "nested_b_in_a"
    return "independent"
def band(present):
    cc = Counter()
    for d in allp:
        rr, ra, ar, aa = cells(d); N = int(d["coread"]); mv = minvaf(d)
        cc[classify_flags(present(aa, N, mv), present(ra, N, mv), present(ar, N, mv))] += 1
    # 合併 nested
    return {"co_linked": cc["co_linked"], "nested": cc["nested_a_in_b"]+cc["nested_b_in_a"],
            "independent": cc["independent"], "mutual_excl": cc["mutual_excl"], "sparse": cc.get("sparse", 0)}
bands = {
    "現行(≥1)": band(lambda c, N, mv: c >= 1),
    "絕對≥2": band(lambda c, N, mv: c >= 2),
    "ε=2%": band(lambda c, N, mv: c > N*0.02),
    "組合(ε2% 且 單讀需VAF期望<3)": band(lambda c, N, mv: (c > N*0.02) and not (c == 1 and N*mv >= 3)),
}

out = {"weak_n": len(weak), "calibration_on_weak": calib, "structure_band": bands}
with open(os.path.join(DATA, "threshold_comparison.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

print(f"=== 弱判定 n={len(weak)};各規則 reject/keep 與 FP 裁判 ===")
print(f"{'規則':<22}{'拒n':>7}{'拒%':>6}{'拒FP%':>7}{'留FP%':>7}{'分離度':>7}")
for nm in names:
    c = calib[nm]
    print(f"{nm:<22}{c['reject_n']:>7}{c['reject_pct']:>6}{str(c['reject_bothFP%']):>7}{str(c['keep_bothFP%']):>7}{str(c['separation(rejFP-keepFP)']):>7}")
print(f"\n=== 整體結構 band ===")
print(f"{'規則':<26}{'co_linked':>10}{'nested':>8}{'indep':>7}{'mut_ex':>7}{'sparse':>7}")
for nm, b in bands.items():
    print(f"{nm:<26}{b['co_linked']:>10}{b['nested']:>8}{b['independent']:>7}{b['mutual_excl']:>7}{b['sparse']:>7}")
print("OK wrote threshold_comparison.json")
