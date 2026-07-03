#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""完整 resolver v2 — 三盲點處理 + VAF 交叉核對 + 樹結構分類統計。
盲點1: truncated → densest-8 對映(復現 sm_multilocus L57-60)後分類
盲點2: CONFLICT/ORDERED 皆做 ε=2% 分層(min(RA,AR)<ε → 衝突降級為 nested;中間群<ε → 定序 tentative)
盲點3: ORDERED 方向用 pos_vaf 交叉核對(較早突變應較高 VAF=較 clonal)
輸出: 邊/區分類 + 樹結構 taxonomy + 盲點量化 + 每區 resolved 結構。§13 落 JSON。"""
import json
import os
import sys
from bisect import bisect_left, bisect_right
from itertools import combinations
from collections import Counter, defaultdict

SC = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SC)
import sm_linkage_genomewide as M  # noqa: E402  (canonical HCC1395 pairs = {M.A}/sm_linkage_genomewide.json)
DATA = os.path.normpath(os.path.join(SC, "..", "data"))
T = json.load(open(os.path.join(DATA, "topology_per_region.json")))
det = T["detail"]
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
EPS = 0.02


def region_positions(chrom, start, end, glen):
    arr = som_by_chrom.get(chrom, [])
    s = [p for p in arr[bisect_left(arr, start):bisect_right(arr, end)] if start <= p <= end]
    if len(s) == glen:
        return s
    if len(s) > glen == 8:  # truncated: densest 8 (min span 連續)
        best = min(range(len(s) - 8 + 1), key=lambda i: s[i + 7] - s[i])
        return s[best:best + 8]
    return None  # 仍不齊


def edge_class(chrom, gained_pos, vafmap):
    """回 dict: klass, tier(solid/tentative), order(list), vaf_ok, rels。含 ε 分層。"""
    order_edges = []
    n_block = n_nested = n_conflict = n_nocoread = 0
    tentative = False
    for a, b in combinations(sorted(gained_pos), 2):
        pr = plut.get((chrom, a, b))
        if pr is None:
            n_nocoread += 1
            continue
        c = pr["cells"]
        ra, ar, aa, cr = c["RA"], c["AR"], c["AA"], pr["coread"]
        thr = max(2, EPS * cr)
        if aa < 2:
            n_conflict += 1  # gained 卻互斥=sibling 混入(異常)
            continue
        big, small = max(ra, ar), min(ra, ar)
        if big == 0:  # RA=AR=0 → co_linked
            n_block += 1
        elif small >= thr:  # 雙邊都夠 → 真 4-gamete
            n_conflict += 1
        else:  # 一邊為主、小邊<ε → nested(ε 去噪);方向由 big 邊定
            n_nested += 1
            # ar 大 → (a=A,b=R)存在 → a 單獨出現 → a 先; ra 大 → b 先
            earlier, later = (a, b) if ar >= ra else (b, a)
            order_edges.append((earlier, later))
            if big < thr:
                tentative = True
    if n_conflict > 0:
        klass = "CONFLICT"
    elif n_nocoread > 0:
        klass = "AMBIG_NOCOREAD"
    elif n_nested == 0:
        klass = "BLOCK"
    elif n_block == 0:
        klass = "ORDERED"
    else:
        klass = "PARTIAL_ORDER"
    # VAF 交叉核對
    vaf_ok = None
    if order_edges and vafmap:
        checks = [(vafmap.get(f"{chrom}:{e}"), vafmap.get(f"{chrom}:{l}")) for e, l in order_edges]
        checks = [(x, y) for x, y in checks if x is not None and y is not None]
        if checks:
            vaf_ok = all(x >= y - 0.05 for x, y in checks)  # 較早≥較晚(容忍 0.05)
    return {"klass": klass, "tentative": tentative, "order": order_edges, "vaf_ok": vaf_ok}


edge_agg = Counter()
region_class = Counter()
tree_taxo = Counter()
blind = Counter()
vaf_check = Counter()
regions_out = []
for r in det:
    chrom, start = r["chrom"], r["start"]
    end = start + r["span"]
    glen = r.get("genotype_len", r["n_sSNV"])
    positions = region_positions(chrom, start, end, glen)
    vafmap = r.get("pos_vaf", {})
    if positions is None:
        region_class["mapping_unresolved"] += 1
        blind["對映不齊(非8截斷/複雜)"] += 1
        continue
    if r.get("truncated"):
        blind["truncated(densest-8已對映)"] += 1
    ec = []
    order_all = []
    tent = conflict = block = ordered = False
    for P, C in r.get("edges", []):
        if P == "ROOT":
            continue
        gi = [i for i in range(glen) if P[i] != C[i]]
        if len(gi) <= 1:
            edge_agg["single_mut"] += 1
            continue
        gp = [positions[i] for i in gi]
        res = edge_class(chrom, gp, vafmap)
        edge_agg[res["klass"]] += 1
        ec.append(res["klass"])
        order_all += res["order"]
        tent |= res["tentative"]
        conflict |= (res["klass"] == "CONFLICT")
        block |= (res["klass"] == "BLOCK")
        ordered |= (res["klass"] in ("ORDERED", "PARTIAL_ORDER"))
        if res["klass"] == "BLOCK" and r.get("cn") == "gain":
            blind["block@CN-gain(疑multiplicity非真共event)"] += 1
        if res["tentative"]:
            blind["ε-tentative(中間群<2%,需放寬層)"] += 1
        if res["vaf_ok"] is False:
            vaf_check["VAF方向不一致(需檢)"] += 1
        elif res["vaf_ok"] is True:
            vaf_check["VAF方向一致(佐證)"] += 1
    # 樹結構 taxonomy
    tt = r.get("topology_type", "?")
    if not ec:
        tree_taxo[f"{tt}|全單突變"] += 1
        region_class["clean_single_mut"] += 1
    else:
        if conflict:
            region_class["CONFLICT"] += 1
            tree_taxo["conflict(非單樹)"] += 1
        elif block and not ordered:
            region_class["BLOCK(共event)"] += 1
            tree_taxo["block-collapsed(多突變一起)"] += 1
        elif ordered and not block:
            region_class["ORDERED(可定序)" + ("|tentative" if tent else "|solid")] += 1
            tree_taxo["ordered-lineage(定序)"] += 1
        else:
            region_class["MIXED(block+定序)"] += 1
            tree_taxo["mixed"] += 1
    regions_out.append({"region": r["region"], "n_sSNV": r["n_sSNV"], "cn": r.get("cn"),
                        "orig_determinacy": r.get("determinacy"), "resolved_edges": ec,
                        "tentative": tent, "n_order_constraints": len(order_all)})

out = {"params": {"EPS": EPS, "sample": "HCC1395", "n_regions": len(det)},
       "edge_classification": dict(edge_agg), "region_classification": dict(region_class),
       "tree_taxonomy": dict(tree_taxo), "blind_spots": dict(blind), "vaf_crosscheck": dict(vaf_check),
       "regions": regions_out}
json.dump(out, open(os.path.join(DATA, "resolve_v2_complete.json"), "w"), ensure_ascii=False)
for title, dd in [("EDGE 分類", edge_agg), ("REGION 分類", region_class),
                  ("樹結構 taxonomy", tree_taxo), ("盲點量化", blind), ("VAF 交叉核對", vaf_check)]:
    print(f"=== {title} ===")
    for k, v in dd.most_common():
        print(f"  {k:<38} {v}")
