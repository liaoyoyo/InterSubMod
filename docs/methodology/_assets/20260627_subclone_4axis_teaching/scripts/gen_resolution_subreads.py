#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""每區完整解析 + 所有子read + 比例（供工作站顯示）。
對每個『多突變邊』區 + worked examples：
  - 所有子read genotype（含 partial 覆蓋,X=未覆蓋;grouped by 向量→count）
  - full-span 向量(絕對) + partial(部分覆蓋)
  - pairwise 2×2(全對) + relation
  - 每條 ≥2-gained 邊的解析(BLOCK/ORDERED/CONFLICT + ε tier + VAF)
  - provenance: 絕對 direct / 組合 pairwise / 比例
  - 盲點旗標
§13: per_read_calls 真實資料 → 落 JSON。"""
import json
import os
import sys
from bisect import bisect_left, bisect_right
from itertools import combinations
from collections import Counter, defaultdict
import pysam

SCD = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"
sys.path.insert(0, SCD)
import sm_linkage_genomewide as M  # noqa: E402

SC = os.path.dirname(os.path.abspath(__file__))
DATA = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data"
T = json.load(open(f"{DATA}/topology_per_region.json"))["detail"]
L = json.load(open(os.path.join(M.A, "sm_linkage_genomewide.json")))  # canonical HCC1395 pairs
census = L["census"]
plut = {(p["chrom"], min(p["a"], p["b"]), max(p["a"], p["b"])): p for p in L["pairs"]}
V2 = json.load(open(os.path.join(DATA, "resolve_v2_complete.json")))
multi_regions = {r["region"] for r in V2["regions"] if r["resolved_edges"]}
WORKED = ["chr8:58338760-58349709", "chr5:170542421-170563845", "chr17:48357368-48365161"]
targets = sorted(multi_regions | set(WORKED))
by_region = {r["region"]: r for r in T}
som_by_chrom = defaultdict(list)
for key, c in census.items():
    if c.get("somatic") is True:
        cc, pp = key.rsplit(":", 1)
        som_by_chrom[cc].append(int(pp))
for cc in som_by_chrom:
    som_by_chrom[cc].sort()
EPS = 0.02
MINREAD = 3


def positions_for(r):
    chrom, start = r["chrom"], r["start"]
    end = start + r["span"]
    glen = r.get("genotype_len", r["n_sSNV"])
    arr = som_by_chrom.get(chrom, [])
    s = [p for p in arr[bisect_left(arr, start):bisect_right(arr, end)] if start <= p <= end]
    if len(s) == glen:
        return s
    if len(s) > glen == 8:
        best = min(range(len(s) - 8 + 1), key=lambda i: s[i + 7] - s[i])
        return s[best:best + 8]
    return None


def edge_resolution(chrom, gained_pos, vafmap):
    rels = []
    n_block = n_nested = n_conflict = n_nocoread = 0
    tent = False
    order = []
    for a, b in combinations(sorted(gained_pos), 2):
        pr = plut.get((chrom, a, b))
        if pr is None:
            n_nocoread += 1
            rels.append({"pair": [a, b], "rel": "none"})
            continue
        c = pr["cells"]
        ra, ar, aa, cr = c["RA"], c["AR"], c["AA"], pr["coread"]
        thr = max(2, EPS * cr)
        big, small = max(ra, ar), min(ra, ar)
        if aa < 2:
            rel = "mutual_excl"
            n_conflict += 1
        elif big == 0:
            rel = "co_linked"
            n_block += 1
        elif small >= thr:
            rel = "independent_4gamete"
            n_conflict += 1
        else:
            rel = "nested"
            n_nested += 1
            order.append([a, b] if ar >= ra else [b, a])
            if big < thr:
                tent = True
        rels.append({"pair": [a, b], "cells": c, "coread": cr, "rel": rel})
    if n_conflict:
        klass = "CONFLICT"
    elif n_nocoread:
        klass = "AMBIG_NOCOREAD"
    elif n_nested == 0:
        klass = "BLOCK"
    elif n_block == 0:
        klass = "ORDERED"
    else:
        klass = "PARTIAL_ORDER"
    vaf_ok = None
    if order and vafmap:
        ck = [(vafmap.get(f"{chrom}:{e}"), vafmap.get(f"{chrom}:{l}")) for e, l in order]
        ck = [(x, y) for x, y in ck if x is not None and y is not None]
        vaf_ok = all(x >= y - 0.05 for x, y in ck) if ck else None
    return {"klass": klass, "tentative": tent, "order": order, "vaf_ok": vaf_ok, "pair_rels": rels}


tb = pysam.AlignmentFile(M.TBAM, "rb")
out = {}
for reg in targets:
    r = by_region.get(reg)
    if r is None:
        continue
    positions = positions_for(r)
    if positions is None:
        continue
    chrom = r["chrom"]
    som = [(p, census[f"{chrom}:{p}"]["ref"], census[f"{chrom}:{p}"]["alt"], "") for p in positions]
    calls, hp, ps = M.per_read_calls(tb, chrom, som)
    k = len(positions)
    # 所有子read: grouped vector(含 X=未覆蓋)
    subread = Counter()
    full = Counter()
    cover_hist = Counter()
    for rn, c in calls.items():
        vec = "".join("A" if c.get(p) == "ALT" else ("R" if c.get(p) == "REF" else "X") for p in positions)
        nc = k - vec.count("X")
        if nc == 0:
            continue
        subread[vec] += 1
        cover_hist[nc] += 1
        if nc == k:
            full[vec.replace("X", "")] += 1
    n_full = sum(full.values())
    # pairwise 2×2 全對
    pw = []
    for i, j in combinations(range(k), 2):
        pi, pj = positions[i], positions[j]
        t = Counter()
        for rn, c in calls.items():
            if c.get(pi) in ("REF", "ALT") and c.get(pj) in ("REF", "ALT"):
                t[("A" if c[pi] == "ALT" else "R") + ("A" if c[pj] == "ALT" else "R")] += 1
        pw.append({"i": i, "j": j, "pos": [pi, pj],
                   "cells": {cc: t.get(cc, 0) for cc in ("RR", "RA", "AR", "AA")}})
    # 邊解析
    vafmap = r.get("pos_vaf", {})
    glen = r.get("genotype_len", r["n_sSNV"])
    edge_res = []
    for P, C in r.get("edges", []):
        if P == "ROOT":
            continue
        gi = [ii for ii in range(glen) if P[ii] != C[ii]]
        if len(gi) < 2:
            continue
        gp = [positions[ii] for ii in gi]
        res = edge_resolution(chrom, gp, vafmap)
        res.update({"parent": P, "child": C, "gained_pos": gp})
        edge_res.append(res)
    # provenance: 絕對(full-span pop>=MINREAD) vs 組合(pairwise 可推但無全跨)
    abs_pops = {v: n for v, n in full.items() if n >= MINREAD}
    n_abs = len(abs_pops)
    tot_abs_reads = sum(abs_pops.values())
    flags = []
    for e in edge_res:
        if e["klass"] == "BLOCK" and r.get("cn") == "gain":
            flags.append("block@CN-gain(疑multiplicity)")
        if e["tentative"]:
            flags.append("ε-tentative(中間群<2%)")
        if e["klass"] == "CONFLICT" and r.get("cn") == "gain":
            flags.append("conflict@CN-gain")
    out[reg] = {
        "n_sSNV": k, "positions": positions, "cn": r.get("cn"),
        "determinacy": r.get("determinacy"), "topology_type": r.get("topology_type"),
        "node_paths": r.get("node_paths", {}), "pos_vaf": vafmap,
        "n_reads_touch": len(subread), "n_full_span": n_full,
        "cover_hist": dict(sorted(cover_hist.items())),
        "subread_vectors": dict(sorted(subread.items(), key=lambda x: -x[1])),
        "full_span_pops": dict(sorted(full.items(), key=lambda x: -x[1])),
        "pairwise": pw, "edge_resolution": edge_res,
        "provenance": {"n_absolute_pops": n_abs, "absolute_reads": tot_abs_reads,
                       "abs_frac": round(tot_abs_reads / max(1, sum(subread.values())), 3)},
        "blind_flags": sorted(set(flags)),
    }
tb.close()
json.dump({"params": {"EPS": EPS, "MINREAD": MINREAD, "sample": "HCC1395", "n_regions": len(out)},
           "regions": out}, open(os.path.join(DATA, "resolution_subreads.json"), "w"), ensure_ascii=False)
print(f"DONE regions={len(out)} -> resolution_subreads.json")
# 快速比例摘要
kc = Counter()
for reg, d in out.items():
    for e in d["edge_resolution"]:
        kc[e["klass"]] += 1
print("edge klass:", dict(kc))
print("含盲點旗標區:", sum(1 for d in out.values() if d["blind_flags"]))
