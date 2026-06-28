#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
單-sSNV(無法建樹)完整 accounting：全 35,332 sSNV 宇宙分桶 + 「無法處理 vs 其實有資訊」判讀。
回答：整體各情況比例 + 單個位點是否真無法處理(or 有 CCF/Tier-PS 可用)。
來源：per_sSNV_census.tsv(n_partners_le50k/n_realized_links/vaf/somatic/cn) + completeness_ledger。
輸出：single_snv_accounting.json
用法：python3 single_snv_accounting.py
"""
import csv, json, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
rows = list(csv.DictReader(open(os.path.join(DATA, "per_sSNV_census.tsv"), encoding="utf-8"), delimiter="\t"))
def I(x):
    try: return int(x)
    except: return 0
def F(x):
    try: return float(x)
    except: return None
N = len(rows)

buck = Counter(); src_buck = Counter(); up_ccf = Counter(); up_cn = Counter()
for r in rows:
    npart, nlink, v = I(r["n_partners_le50k"]), I(r["n_realized_links"]), F(r["vaf"])
    if nlink > 0: b = "linked"
    elif npart > 0: b = "underpowered"
    else: b = "isolated"
    buck[b] += 1; src_buck[f"{r['src']}|{b}"] += 1
    if b == "underpowered" and v is not None:
        up_ccf["high_ge0.4(clonal)" if v >= 0.4 else ("mid_0.1-0.4" if v >= 0.1 else "low_lt0.1(subclonal/noise)")] += 1
        up_cn[r["cn_state"]] += 1

single = N - buck["linked"]
out = {
    "universe_total": N, "tp": sum(v for k, v in src_buck.items() if k.startswith("TP")),
    "fp": sum(v for k, v in src_buck.items() if k.startswith("FP")),
    "buckets": {
        "linked": {"n": buck["linked"], "pct": round(100*buck["linked"]/N, 1),
                   "status": "可建樹(進拓樸分析;見 topology_workstation)"},
        "underpowered": {"n": buck["underpowered"], "pct": round(100*buck["underpowered"]/N, 1),
                         "status": "有 partner 但無共讀 link → 加深 ONT 覆蓋可救;有 census VAF→CCF 可刻畫",
                         "has_census_vaf": buck["underpowered"], "ccf_tendency": dict(up_ccf), "by_cn": dict(up_cn)},
        "isolated": {"n": buck["isolated"], "pct": round(100*buck["isolated"]/N, 1),
                     "status": "read-span(≤50kb)內無 partner → Tier-R 樹外;此 census 無共讀 VAF",
                     "rescue": "(a) caller(ClairS) VCF 仍有 VAF→CCF 可刻畫 (b) 可能有 same-PS >50kb partner(Tier-PS deferred,未做)→非永久無法處理"},
    },
    "single_sSNV_total": single, "single_pct": round(100*single/N, 1),
    "by_source_bucket": dict(src_buck),
    "verdict": {
        "整體可建樹": f"{buck['linked']} ({round(100*buck['linked']/N,1)}%)",
        "單位點_可救(underpowered,深覆蓋+有CCF)": f"{buck['underpowered']} ({round(100*buck['underpowered']/N,1)}%)",
        "單位點_Tier-R樹外(isolated,有caller VAF+可能Tier-PS)": f"{buck['isolated']} ({round(100*buck['isolated']/N,1)}%)",
        "結論": "單位點非『全無法處理』:underpowered 有 CCF+可深覆蓋救;isolated 有 caller VAF 可刻畫 clonal 譜+可能 Tier-PS 連鎖。真正『拓樸-dead』= isolated 中無 same-PS partner 者(待 Tier-PS 細分量化)。",
    },
}
with open(os.path.join(DATA, "single_snv_accounting.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

print(f"全 sSNV 宇宙 = {N}")
for b in ("linked", "underpowered", "isolated"):
    print(f"  {b}: {buck[b]} ({round(100*buck[b]/N,1)}%) — {out['buckets'][b]['status']}")
print(f"單-sSNV(非linked) = {single} ({round(100*single/N,1)}%)")
print(f"  underpowered CCF: {dict(up_ccf)}")
print(f"by_source: {dict(src_buck)}")
print("OK wrote single_snv_accounting.json")
