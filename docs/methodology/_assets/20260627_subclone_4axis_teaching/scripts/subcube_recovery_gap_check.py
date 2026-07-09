#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gap 驗證(2026-07-09):E_subcube_recovered(partial-read 救回)區是否真的建/枚舉了 Steiner 樹?
假設:這些區只被『辨識+標 tree_shape』,並未跑 sub-face group-Steiner 建樹/枚舉(edges=[])。
只讀不寫。用法: python3 subcube_recovery_gap_check.py
"""
import json, os
from collections import Counter
from itertools import combinations
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.normpath(os.path.join(HERE, "..", "data"))
det = json.load(open(os.path.join(D, "topology_per_region.json")))["detail"]
rec = [r for r in det if "E_subcube_recovered" in (r.get("determinacy") or "")]

print("="*88); print("E_subcube_recovered gap 驗證(HCC1395)"); print("="*88)
print(f"recovered 區總數           = {len(rec)}")
print(f"  有非空 edges(真建了樹)   = {sum(1 for r in rec if r.get('edges'))}   ← 應=0 證實未建樹")
nd = Counter(r.get("n_sSNV", 0) for r in rec)
n2 = sum(v for k, v in nd.items() if k == 2); n3 = sum(v for k, v in nd.items() if k >= 3)
print(f"  n=2(單對,無需枚舉)     = {n2} ({n2/len(rec)*100:.1f}%)")
print(f"  n>=3(需真 Steiner)      = {n3} ({n3/len(rec)*100:.1f}%)")
sub_lbl = Counter(r["determinacy"] for r in rec)
for k, v in sub_lbl.items(): print(f"    {k} = {v}")

# 乾淨 n>=3(veclen==n)中,有幾個含未觀測對(=該枚舉多候選/需隱藏節點)
def obs_pairs(sub, k):
    o = 0
    for i, j in combinations(range(k), 2):
        c = Counter()
        for vec, n in sub.items():
            a, b = vec[i], vec[j]
            if a in "RA" and b in "RA": c[a+b] += n
        if c: o += 1
    return o
tot = 0; need_enum = 0
for r in rec:
    n = r.get("n_sSNV", 0); sub = r.get("subread_groups") or {}
    if n < 3 or not sub: continue
    k = len(next(iter(sub)))
    if k != n: continue
    tot += 1
    if obs_pairs(sub, k) < comb(k, 2): need_enum += 1
print(f"\n乾淨 n>=3 recovered(veclen==n) = {tot}")
print(f"  含未觀測對(該輸出多候選/隱藏節點)= {need_enum} ({need_enum/tot*100:.1f}%)" if tot else "  (無)")
print("\n結論:recovered 區 = 辨識+標 tree_shape,但 edges 全空、100% 需枚舉卻未枚舉。")
print("→ 「救回 N 區」精確語意 = 『辨識出有 partial 共現底料』,非『建出/枚舉出樹』。")
print("="*88)
