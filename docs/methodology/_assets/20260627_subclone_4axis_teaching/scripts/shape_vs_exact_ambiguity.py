#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VAF 破不掉的歧義:是「同形狀不同順序」還是「真多形狀」? (2026-07-10)
每 ambiguous unit:候選樹的 distinct SHAPE 數 × VAF 是否破(top≥0.6) × CN。
形狀 = single/linear/branched/star(拓撲);精確樹 = 含順序+隱藏節點擺放。
問題:VAF 破不掉的(對稱/CN/bulk),形狀是否仍唯一(只內部順序歧義),還是真多形狀。
"""
import json, os, math, glob, sys
from collections import Counter, defaultdict
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts")
import tree_enumeration_solver as S

GLOB = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/mlhp_part_*.json"
MARGIN = 0.05; TEMP = 0.05; GERM = ("1", "2", "3")

def tree_shape(edges):
    if not edges: return "single"
    ch = defaultdict(list); nodes = set()
    for p, c in edges:
        ch[p].append(c); nodes.add(p); nodes.add(c)
    if len([n for n in nodes if n != "ROOT"]) <= 1: return "single"
    cc = {p: len(cs) for p, cs in ch.items()}
    if any(v >= 2 for p, v in cc.items() if p != "ROOT"): return "branched"
    if cc.get("ROOT", 0) >= 2: return "star"
    return "linear"

def _unlab(s):
    if s == "ROOT": return frozenset()
    s2 = s[2:] if s.startswith("H_") else s
    return frozenset(i for i, c in enumerate(s2) if c == "A")

def score_viol(edges, ccf):
    sc = 0.0
    for (p, c) in edges:
        P, C = _unlab(p), _unlab(c); acq = C - P
        if len(acq) != 1: continue
        j = next(iter(acq))
        for anc in P: sc += (ccf[anc] - ccf[j])
    return sc

def ccf_cc(colcov, positions, k):
    if not colcov or len(positions) < k: return None
    ccf = {}
    for j in range(k):
        v = colcov.get(str(positions[j])) or colcov.get(positions[j])
        if not v: return None
        tot = v[0] + v[1]; ccf[j] = v[1] / tot if tot else 0.0
    return ccf

groups = []
for f in sorted(glob.glob(GLOB)):
    groups += json.load(open(f)).get("groups", [])

# 累計:(broke, n_shapes==1?, cn) → count
tab = Counter()
shape_when_unbroken = Counter()
shape_when_broke_top = Counter()
n_amb = 0
for r in groups:
    if r.get("n_sSNV", 0) < 2: continue
    pbh = r.get("populations_by_hp", {}) or {}; sbh = r.get("subread_groups_by_hp", {}) or {}
    cbh = r.get("col_coverage_by_hp", {}) or {}; positions = r.get("positions", []); cn = r.get("cn")
    for fam in sorted(set(pbh) | set(sbh)):
        if fam not in GERM: continue
        full = pbh.get(fam, {}) or {}; part = list((sbh.get(fam, {}) or {}).keys())
        if not full and not part: continue
        k = len(next(iter(full))) if full else len(part[0])
        if k < 2 or k > 8: continue
        res = S.enumerate_min_trees(full, part, k)
        cls = S.classify(res)
        if res["capped"] or "ambiguous" not in cls or res["n_trees"] < 2 or len(res["trees"]) < 2:
            continue
        n_amb += 1
        shapes = [tree_shape(t["edges"]) for t in res["trees"]]
        n_shapes = len(set(shapes))
        cn_class = "gain" if cn == "gain" else ("clean" if cn in ("neutral", "loh", "loss") else "unknown")
        # VAF 打分
        ccf = ccf_cc(cbh.get(fam, {}), positions, k)
        if ccf is None:
            tot = sum(full.values()) or 1
            ccf = {j: sum(n for g, n in full.items() if g[j] == "A") / tot for j in range(k)}
        scored = [(score_viol(t["edges"], ccf), tree_shape(t["edges"])) for t in res["trees"]]
        mx = max(s for s, _ in scored)
        exps = [math.exp((s - mx) / TEMP) for s, _ in scored]; tot = sum(exps) or 1
        post = [e / tot for e in exps]; top = max(post)
        broke = top >= 0.6
        single_shape = (n_shapes == 1)
        tab[(broke, single_shape, cn_class)] += 1
        if not broke:
            shape_when_unbroken[n_shapes] += 1
        else:
            # 破掉時,top 樹的形狀是否是唯一被選中的形狀
            topshape = scored[post.index(top)][1]
            # top-posterior 諸樹(≥0.6 佔比者)的形狀集合
            shape_when_broke_top[topshape] += 1

print(f"n_ambiguous = {n_amb}")
print("\n=== 全 ambiguous:形狀是否唯一(n_distinct_shapes==1 = 拓撲已定,只內部順序/隱藏擺放歧義) ===")
single_total = sum(v for (b, s, c), v in tab.items() if s)
print(f"  單一形狀(拓撲已定): {single_total} ({100*single_total/n_amb:.1f}%)")
print(f"  多形狀(真的多拓撲): {n_amb-single_total} ({100*(n_amb-single_total)/n_amb:.1f}%)")

print("\n=== VAF 破不掉的(top<0.6:對稱/難) → 形狀分布 ===")
ub = sum(shape_when_unbroken.values())
for ns in sorted(shape_when_unbroken):
    print(f"  n_shapes={ns}: {shape_when_unbroken[ns]} ({100*shape_when_unbroken[ns]/ub:.1f}%)")
ub_single = shape_when_unbroken.get(1, 0)
print(f"  → VAF 破不掉者中,單一形狀(拓撲仍定)= {ub_single} ({100*ub_single/ub:.1f}%);真多形狀 = {ub-ub_single} ({100*(ub-ub_single)/ub:.1f}%)")

print("\n=== 交叉:(VAF破? × 單一形狀? × CN) ===")
for (b, s, c), v in sorted(tab.items()):
    print(f"  破={b!s:5} 單形狀={s!s:5} CN={c:7}: {v}")

print("\n=== CN-gain vs clean 的形狀唯一性(CN 是否製造多形狀?) ===")
for cnc in ["clean", "gain"]:
    n = sum(v for (b, s, c), v in tab.items() if c == cnc)
    ns1 = sum(v for (b, s, c), v in tab.items() if c == cnc and s)
    if n: print(f"  {cnc}: 單一形狀 {ns1}/{n} = {100*ns1/n:.1f}%")
