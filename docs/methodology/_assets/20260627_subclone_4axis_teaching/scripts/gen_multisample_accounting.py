#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
逐樣本 single_snv_accounting.json 生成器(R1)。
複用 single_snv_accounting.py 的分桶邏輯(linked/underpowered/isolated + tp/fp + CCF tendency),
對 6 multisample 樣本各讀 lists/per_sSNV_census.tsv → 寫 single_snv_accounting.json 進該樣本目錄。
誠實旗標:fp_truth_sparse(無外部 truth set)、cn_all_neutral(census 未併 CN)、derived_from_census(衍生非凍結 pipeline)。
HCC1395(4axis)同邏輯重算→與凍結 json 對照(不覆寫凍結檔,只印驗證)。
用法: python3 gen_multisample_accounting.py
"""
import csv, json, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA4AXIS = os.path.normpath(os.path.join(HERE, "..", "data"))
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"

def I(x):
    try: return int(x)
    except: return 0
def F(x):
    try: return float(x)
    except: return None

def find_census(dr):
    for cand in (os.path.join(dr, "lists", "per_sSNV_census.tsv"),
                 os.path.join(dr, "per_sSNV_census.tsv")):
        if os.path.exists(cand): return cand
    return None

def build_accounting(census_path):
    rows = list(csv.DictReader(open(census_path, encoding="utf-8"), delimiter="\t"))
    N = len(rows)
    buck = Counter(); src_buck = Counter(); up_ccf = Counter(); up_cn = Counter(); cn_all = Counter()
    for r in rows:
        npart, nlink, v = I(r["n_partners_le50k"]), I(r["n_realized_links"]), F(r["vaf"])
        b = "linked" if nlink > 0 else ("underpowered" if npart > 0 else "isolated")
        buck[b] += 1; src_buck[f"{r['src']}|{b}"] += 1
        cn_all[r.get("cn_state", "?")] += 1
        if b == "underpowered" and v is not None:
            up_ccf["high_ge0.4(clonal)" if v >= 0.4 else ("mid_0.1-0.4" if v >= 0.1 else "low_lt0.1(subclonal/noise)")] += 1
            up_cn[r["cn_state"]] += 1
    if not N: return None
    tp = sum(v for k, v in src_buck.items() if k.startswith("TP"))
    fp = sum(v for k, v in src_buck.items() if k.startswith("FP"))
    single = N - buck["linked"]
    cn_states = set(cn_all.keys())
    out = {
        "universe_total": N, "tp": tp, "fp": fp,
        "buckets": {
            "linked": {"n": buck["linked"], "pct": round(100*buck["linked"]/N, 1),
                       "status": "可建樹(進拓樸分析)"},
            "underpowered": {"n": buck["underpowered"], "pct": round(100*buck["underpowered"]/N, 1),
                             "status": "有 partner 但無共讀 link → 加深覆蓋可救;有 census VAF→CCF 可刻畫",
                             "has_census_vaf": buck["underpowered"], "ccf_tendency": dict(up_ccf), "by_cn": dict(up_cn)},
            "isolated": {"n": buck["isolated"], "pct": round(100*buck["isolated"]/N, 1),
                         "status": "read-span(≤50kb)內無 partner → Tier-R 樹外"},
        },
        "single_sSNV_total": single, "single_pct": round(100*single/N, 1),
        "by_source_bucket": dict(src_buck),
        # R1 誠實旗標(衍生 artifact 用)
        "derived_from_census": True,
        "fp_truth_sparse": (fp / N < 0.005),
        "cn_all_neutral": (cn_states == {"neutral"} or cn_states == {"neutral", "?"}),
    }
    return out

samples = []
if os.path.isdir(MSROOT):
    for s in sorted(os.listdir(MSROOT)):
        if os.path.exists(os.path.join(MSROOT, s, "topology_per_region.json")):
            samples.append((s, os.path.join(MSROOT, s)))

print("%-16s %8s %8s %7s %8s %6s %6s  flags" % ("sample", "universe", "linked", "tp", "fp", "under", "iso"))
for name, dr in samples:
    cp = find_census(dr)
    if not cp:
        print("%-16s  NO census (skip)" % name); continue
    acc = build_accounting(cp)
    with open(os.path.join(dr, "single_snv_accounting.json"), "w", encoding="utf-8") as f:
        json.dump(acc, f, ensure_ascii=False, indent=1)
    flags = []
    if acc["fp_truth_sparse"]: flags.append("FP-sparse")
    if acc["cn_all_neutral"]: flags.append("CN-neutral")
    print("%-16s %8d %8d %7d %8d %6d %6d  %s" % (name, acc["universe_total"], acc["buckets"]["linked"]["n"],
          acc["tp"], acc["fp"], acc["buckets"]["underpowered"]["n"], acc["buckets"]["isolated"]["n"], ",".join(flags)))

# HCC1395 4axis 對照(不覆寫凍結檔)
cp = find_census(DATA4AXIS)
if cp:
    acc = build_accounting(cp)
    frozen = json.load(open(os.path.join(DATA4AXIS, "single_snv_accounting.json"), encoding="utf-8"))
    match = (acc["universe_total"] == frozen["universe_total"] and acc["tp"] == frozen["tp"] and
             acc["fp"] == frozen["fp"] and acc["buckets"]["linked"]["n"] == frozen["buckets"]["linked"]["n"])
    print("HCC1395(4axis) 重算 vs 凍結 json: universe %d/%d tp %d/%d linked %d/%d → %s" % (
        acc["universe_total"], frozen["universe_total"], acc["tp"], frozen["tp"],
        acc["buckets"]["linked"]["n"], frozen["buckets"]["linked"]["n"], "MATCH ✓" if match else "MISMATCH ✗"))
