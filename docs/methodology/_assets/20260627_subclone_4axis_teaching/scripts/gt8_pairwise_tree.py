#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[>8 sSNV 全 pairwise 建樹] 去 MAX_SNV=8 cap — 對 truncated(>8)區用『全部 sSNV 的 pairwise 關係』
建 position-level perfect-phylogeny 樹(不需單一全跨 read):
  co_linked(RA=AR=0)  → union 併節點(同一 lineage 事件)
  nested(小邊<ε)       → 祖→裔 有向邊(ε 去噪)
  4-gamete(雙邊≥ε)     → conflict pair(計數)
  → union-find 群 + 祖裔 DAG + 環偵測 + transitive reduction → 樹。
🔴 CN-clean(loh/neutral/loss)可信;CN-gain 標 multiplicity-caveat(co_linked 可能是擴增假共現非真共event)。
provenance: 節點為 pairwise 組合(B_pairwise)非單分子(A_determined);無 read 跨全部。
來源:topology(truncated 名單)+ sm_linkage pairs+census。輸出:gt8_trees.json。§13 純分析。
用法: SM_DATA=<dir> python3 gt8_pairwise_tree.py  (pairs 從 SM_DATA 或 {M.A})"""
import os, sys, json
from bisect import bisect_left, bisect_right
from itertools import combinations
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(HERE, "..", "data")))
det = json.load(open(os.path.join(DATA, "topology_per_region.json"), encoding="utf-8"))["detail"]
_pp = os.path.join(DATA, "sm_linkage_genomewide.json")
if not os.path.exists(_pp):
    sys.path.insert(0, HERE)
    import sm_linkage_genomewide as M  # noqa
    _pp = os.path.join(M.A, "sm_linkage_genomewide.json")
L = json.load(open(_pp, encoding="utf-8"))
census = L["census"]
plut = {(p["chrom"], min(p["a"], p["b"]), max(p["a"], p["b"])): p for p in L["pairs"]}
som_by_chrom = defaultdict(list)
for k, c in census.items():
    if c.get("somatic") is True:
        cc, pp = k.rsplit(":", 1); som_by_chrom[cc].append(int(pp))
for cc in som_by_chrom:
    som_by_chrom[cc].sort()
EPS = 0.02


def build(chrom, som, vafmap):
    par = {p: p for p in som}

    def find(x):
        r = x
        while par[r] != r:
            r = par[r]
        while par[x] != r:
            par[x], x = r, par[x]
        return r

    anc, n_conf, npow = [], 0, 0
    for a, b in combinations(som, 2):
        pr = plut.get((chrom, a, b))
        if pr is None or not pr.get("powered"):
            continue
        npow += 1
        c = pr["cells"]; ra, ar, aa, cr = c["RA"], c["AR"], c["AA"], pr["coread"]
        if aa < 2:
            continue                    # mutual_excl/sparse → sibling(無祖裔邊)
        thr = max(2, EPS * cr); big, small = max(ra, ar), min(ra, ar)
        if big == 0:
            par[find(a)] = find(b)      # co_linked
        elif small >= thr:
            n_conf += 1                 # 4-gamete
        else:
            anc.append((a, b) if ar >= ra else (b, a))   # ar 大=a 單獨=a 祖
    grp = defaultdict(set)
    for p in som:
        grp[find(p)].add(p)
    gid = {p: find(p) for p in som}
    gedges = set()
    for a, d in anc:
        if gid[a] != gid[d]:
            gedges.add((gid[a], gid[d]))
    # 環偵測 + 去背向邊(DFS spanning) → DAG
    adj = defaultdict(set)
    for u, v in gedges:
        adj[u].add(v)
    color, has_cycle, kept = {}, False, set()
    WHITE, GRAY, BLACK = 0, 1, 2

    def dfs(u):
        nonlocal has_cycle
        color[u] = GRAY
        for v in sorted(adj[u]):
            if color.get(v, WHITE) == GRAY:
                has_cycle = True; continue      # back-edge → 去除
            kept.add((u, v))
            if color.get(v, WHITE) == WHITE:
                dfs(v)
        color[u] = BLACK
    for g in list(grp):
        if color.get(g, WHITE) == WHITE:
            dfs(g)
    # transitive reduction on kept
    radj = defaultdict(set)
    for u, v in kept:
        radj[u].add(v)

    def reach(u, v, skip):
        st = [u]; seen = set()
        while st:
            x = st.pop()
            for y in radj[x]:
                if (x, y) == skip:
                    continue
                if y == v:
                    return True
                if y not in seen:
                    seen.add(y); st.append(y)
        return False
    tree = [[u, v] for (u, v) in kept if not reach(u, v, (u, v))]
    indeg = defaultdict(int)
    for u, v in tree:
        indeg[v] += 1
    roots = [g for g in grp if indeg[g] == 0]
    # depth
    depth = {g: 0 for g in roots}
    order = list(roots)
    ci = defaultdict(list)
    for u, v in tree:
        ci[u].append(v)
    i = 0
    while i < len(order):
        u = order[i]; i += 1
        for v in ci[u]:
            depth[v] = depth.get(u, 0) + 1; order.append(v)
    nodes = [{"gid": str(g), "pos": sorted(grp[g]), "size": len(grp[g])} for g in grp]
    return {"n_positions": len(som), "n_nodes": len(grp), "n_tree_edges": len(tree),
            "n_roots": len(roots), "max_depth": max(depth.values()) if depth else 0,
            "n_conflict_pairs": n_conf, "n_powered_pairs": npow, "has_cycle": has_cycle,
            "conflict_frac": round(n_conf / max(1, npow), 3),
            "nodes": nodes, "edges": [[str(u), str(v)] for u, v in tree],
            "roots": [str(g) for g in roots]}


out, agg = {}, defaultdict(int)
for r in det:
    if not r.get("truncated"):
        continue
    chrom, start, end = r["chrom"], r["start"], r["start"] + r["span"]
    arr = som_by_chrom.get(chrom, [])
    som = [p for p in arr[bisect_left(arr, start):bisect_right(arr, end)] if start <= p <= end]
    if len(som) < 2:
        continue
    t = build(chrom, som, r.get("pos_vaf", {}))
    cn = r.get("cn")
    caveat = "CN-gain:pairwise 一致≠真樹(co_linked 可能 multiplicity 假共現);需 SAVANA CN" if cn == "gain" else None
    conflict = t["conflict_frac"] >= 0.05 or t["has_cycle"]
    t.update({"region": r["region"], "cn": cn, "span": r["span"],
              "cn_class": "clean" if cn in ("loh", "neutral", "loss") else "gain_caveat",
              "verdict": "CONFLICT(CN-multiplicity)" if conflict else "pairwise-tree",
              "cn_caveat": caveat, "vector_cap_nodes": r.get("n_clusters"),
              "provenance": "B_pairwise(組合·非單分子;無 read 跨全部 sSNV)"})
    out[r["region"]] = t
    agg["total"] += 1
    agg[t["verdict"]] += 1
    agg["cn_" + t["cn_class"]] += 1
    if t["verdict"] == "pairwise-tree" and t["cn_class"] == "clean":
        agg["clean_pairwise_tree(可信)"] += 1

json.dump({"params": {"EPS": EPS, "sample": "HCC1395"}, "aggregate": dict(agg), "trees": out},
          open(os.path.join(DATA, "gt8_trees.json"), "w"), ensure_ascii=False)
print("GT8 DONE →", os.path.join(DATA, "gt8_trees.json"))
print(json.dumps(dict(agg), ensure_ascii=False, indent=1))
# 展示:recover 的節點 vs 8-cap
print("=== recover 範例(全 pairwise 節點 vs 8-cap n_clusters)===")
for reg, t in sorted(out.items(), key=lambda x: -x[1]["n_positions"])[:8]:
    print(f"  {reg:<26} {t['n_positions']}sSNV → {t['n_nodes']}節點/{t['max_depth']}深 (8-cap群={t['vector_cap_nodes']}) {t['verdict']} cn={t['cn']}")
