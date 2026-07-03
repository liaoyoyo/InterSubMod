#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全 3885 區子read + provenance（per-chrom 平行）。B: 擴到每區都能看子read。
argv[1]=chroms(逗號) argv[2]=out_part。多突變區含完整 edge_resolution+pairwise;單突變區只 subread+provenance。
§13: per_read_calls 真實 → 落 part JSON,由 merge 合。"""
import json
import os
import sys
from bisect import bisect_left, bisect_right
from itertools import combinations
from collections import Counter, defaultdict
import pysam

SCD = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCD)
import sm_linkage_genomewide as M  # noqa: E402

DATA = os.path.normpath(os.path.join(SCD, "..", "data"))
CHROMS = set(sys.argv[1].split(",")) if len(sys.argv) > 1 else None
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(DATA, "subreads_all.json")
EPS, MINREAD = 0.02, 3

T = json.load(open(os.path.join(DATA, "topology_per_region.json")))["detail"]
L = json.load(open(os.path.join(M.A, "sm_linkage_genomewide.json")))
census = L["census"]
plut = {(p["chrom"], min(p["a"], p["b"]), max(p["a"], p["b"])): p for p in L["pairs"]}
som_by_chrom = defaultdict(list)
for key, c in census.items():
    if c.get("somatic") is True:
        cc, pp = key.rsplit(":", 1)
        som_by_chrom[cc].append(int(pp))
for cc in som_by_chrom:
    som_by_chrom[cc].sort()


def positions_for(r, glen):
    arr = som_by_chrom.get(r["chrom"], [])
    start, end = r["start"], r["start"] + r["span"]
    s = [p for p in arr[bisect_left(arr, start):bisect_right(arr, end)] if start <= p <= end]
    if len(s) == glen:
        return s
    if len(s) > glen == 8:
        best = min(range(len(s) - 8 + 1), key=lambda i: s[i + 7] - s[i])
        return s[best:best + 8]
    return None


def edge_res(chrom, gpos, vafmap):
    order, nb, nn, nc, nx, tent = [], 0, 0, 0, 0, False
    for a, b in combinations(sorted(gpos), 2):
        pr = plut.get((chrom, a, b))
        if pr is None:
            nx += 1; continue
        c = pr["cells"]; ra, ar, aa, cr = c["RA"], c["AR"], c["AA"], pr["coread"]
        thr = max(2, EPS * cr); big, small = max(ra, ar), min(ra, ar)
        if aa < 2:
            nc += 1
        elif big == 0:
            nb += 1
        elif small >= thr:
            nc += 1
        else:
            nn += 1; order.append([a, b] if ar >= ra else [b, a])
            if big < thr:
                tent = True
    klass = "CONFLICT" if nc else ("AMBIG_NOCOREAD" if nx else ("BLOCK" if nn == 0 else ("ORDERED" if nb == 0 else "PARTIAL_ORDER")))
    vaf_ok = None
    if order and vafmap:
        ck = [(vafmap.get(f"{chrom}:{e}"), vafmap.get(f"{chrom}:{l}")) for e, l in order]
        ck = [(x, y) for x, y in ck if x is not None and y is not None]
        vaf_ok = all(x >= y - 0.05 for x, y in ck) if ck else None
    return {"klass": klass, "order": order, "tentative": tent, "vaf_ok": vaf_ok}


tb = pysam.AlignmentFile(M.TBAM, "rb")
out = {}
for r in T:
    if CHROMS and r["chrom"] not in CHROMS:
        continue
    glen = r.get("genotype_len", r["n_sSNV"])
    positions = positions_for(r, glen)
    if positions is None:
        continue
    chrom = r["chrom"]
    try:
        som = [(p, census[f"{chrom}:{p}"]["ref"], census[f"{chrom}:{p}"]["alt"], "") for p in positions]
    except KeyError:
        continue
    calls, hp, ps = M.per_read_calls(tb, chrom, som)
    k = len(positions)
    subread, full, cover = Counter(), Counter(), Counter()
    for rn, c in calls.items():
        vec = "".join("A" if c.get(p) == "ALT" else ("R" if c.get(p) == "REF" else "X") for p in positions)
        ncov = k - vec.count("X")
        if ncov == 0:
            continue
        subread[vec] += 1; cover[ncov] += 1
        if ncov == k:
            full[vec.replace("X", "")] += 1
    abs_reads = sum(n for v, n in full.items() if n >= MINREAD)
    tot = sum(subread.values()) or 1
    rec = {"n_sSNV": k, "positions": positions, "cn": r.get("cn"),
           "determinacy": r.get("determinacy"), "n_full_span": sum(full.values()),
           "cover_hist": dict(sorted(cover.items())),
           "subread_vectors": dict(sorted(subread.items(), key=lambda x: -x[1])[:50]),
           "provenance": {"n_absolute_pops": sum(1 for v, n in full.items() if n >= MINREAD),
                          "absolute_reads": abs_reads, "abs_frac": round(abs_reads / tot, 3)}}
    # 多突變區: 加 pairwise + edge_resolution
    vafmap = r.get("pos_vaf", {})
    er = []
    for P, C in r.get("edges", []):
        if P == "ROOT":
            continue
        gi = [i for i in range(glen) if P[i] != C[i]]
        if len(gi) >= 2:
            res = edge_res(chrom, [positions[i] for i in gi], vafmap)
            res.update({"parent": P, "child": C, "gained_pos": [positions[i] for i in gi]})
            er.append(res)
    if er:
        pw = []
        for i, j in combinations(range(k), 2):
            pi, pj = positions[i], positions[j]
            t = Counter()
            for rn, c in calls.items():
                if c.get(pi) in ("REF", "ALT") and c.get(pj) in ("REF", "ALT"):
                    t[("A" if c[pi] == "ALT" else "R") + ("A" if c[pj] == "ALT" else "R")] += 1
            pw.append({"i": i, "j": j, "cells": {cc: t.get(cc, 0) for cc in ("RR", "RA", "AR", "AA")}})
        rec["edge_resolution"] = er
        rec["pairwise"] = pw
        flags = []
        for e in er:
            if e["klass"] == "BLOCK" and r.get("cn") == "gain":
                flags.append("block@CN-gain")
            if e["tentative"]:
                flags.append("ε-tentative")
            if e["klass"] == "CONFLICT" and r.get("cn") == "gain":
                flags.append("conflict@CN-gain")
        rec["blind_flags"] = sorted(set(flags))
    out[r["region"]] = rec
tb.close()
json.dump(out, open(OUT, "w"), ensure_ascii=False)
print(f"DONE chroms={sys.argv[1] if len(sys.argv)>1 else 'ALL'} regions={len(out)} -> {OUT}", flush=True)
