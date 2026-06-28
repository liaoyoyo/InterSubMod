#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用「最後 region 關聯結構的差異 + 穩定性」裁判門檻太嚴/太鬆(非外部 FP 標籤)。
做法：把各門檻傳播到 region 層 → 每區的 nested edge 接回 per-pair cell → 決定該 edge 在門檻下是否存活
     → 重新推導 region 樹形 → 比較最終結構分布差異 + 穩定性 + 變動集中在哪(coread/cn)。
門檻只會「移除」單讀 nested edge(變嚴),故 full_tree/linear 可能塌陷;sibling(aa=0)穩定。
用法：python3 region_threshold_impact.py
"""
import os, csv, json, statistics as st
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

# pair lookup: (chrom, lo, hi) -> deciding off-diag, coread, config
lookup = {}
for c in ("nested_a_in_b", "nested_b_in_a", "co_linked", "independent", "mutual_excl"):
    for hp in ("sameHP", "diffHP"):
        p = os.path.join(DATA, "lists", f"{c}_{hp}.tsv")
        if not os.path.exists(p): continue
        for d in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"):
            pa, pb = int(d["pos_a"]), int(d["pos_b"])
            key = (d["chrom"], min(pa, pb), max(pa, pb))
            ra, ar = int(d["RA"]), int(d["AR"]); cfg = d["config"]
            off = ra if cfg == "nested_a_in_b" else (ar if cfg == "nested_b_in_a" else min(ra, ar))
            lookup[key] = {"off": off, "coread": int(d["coread"]), "config": cfg}

RI = json.load(open(os.path.join(DATA, "sm_region_integration.json"), encoding="utf-8"))
regs = RI["regions"]

def parse(p):  # "chr1:1574734" -> ("chr1",1574734)
    ch, po = p.rsplit(":", 1); return ch, int(po)

def edge_survives(region_chrom, e, eps):
    (ca, pa), (cb, pb) = parse(e[0]), parse(e[1])
    key = (ca, min(pa, pb), max(pa, pb))
    rec = lookup.get(key)
    if rec is None: return None  # 查不到
    return rec["off"] > rec["coread"] * eps

def new_shape(r, eps):
    if r["has_cycle"]: return "inconsistent"
    nz = 0; miss = 0
    for e in r.get("nested_edges", []):
        s = edge_survives(r["chrom"], e, eps)
        if s is None: miss += 1; nz += 1   # 查不到→保守當存活(不誤殺)
        elif s: nz += 1
    nsib = len(r.get("sibling_pairs", []))
    if nz > 0 and nsib > 0: return "full_tree"
    if nz > 0: return "linear_nested"
    if nsib > 0: return "sibling_only"
    # 失去所有 nested+sibling → 看是否還有共連
    return "co_linked_lineage" if (r.get("n_colinked_merges", 0) > 0 or (r.get("populations") and any("A" * r["n_sSNV"] == g for g in r["populations"]))) else "no_confirmed_structure"

EPS = [0.0, 0.01, 0.02, 0.03]  # 0=現行
dists = {}
shapes_by_eps = {}
for eps in EPS:
    if eps == 0.0:
        sh = [r["tree_shape"] for r in regs]
    else:
        sh = [new_shape(r, eps) for r in regs]
    shapes_by_eps[eps] = sh
    dists[eps] = dict(Counter(sh))

# 穩定性：相鄰 ε 之間有多少區換 shape
def changed(a, b): return sum(1 for x, y in zip(shapes_by_eps[a], shapes_by_eps[b]) if x != y)
stability = {
    "現行→ε1%": changed(0.0, 0.01), "ε1%→ε2%": changed(0.01, 0.02), "ε2%→ε3%": changed(0.02, 0.03),
    "現行→ε2%": changed(0.0, 0.02),
}

# 變動集中度：現行→ε2% 從 full_tree/linear 塌陷的區,其 coread/cn vs 維持的
def region_med_coread(r):
    cs = [lookup[k]["coread"] for e in r.get("nested_edges", []) for k in [(parse(e[0])[0], min(parse(e[0])[1], parse(e[1])[1]), max(parse(e[0])[1], parse(e[1])[1]))] if k in lookup]
    return st.median(cs) if cs else None
collapsed, kept = [], []
for i, r in enumerate(regs):
    if shapes_by_eps[0.0][i] in ("full_tree", "linear_nested"):
        (collapsed if shapes_by_eps[0.02][i] not in ("full_tree", "linear_nested") else kept).append(r)
def cn_mix(rows): return dict(Counter(r["cn"] for r in rows))
def med_cr(rows):
    v = [x for r in rows for x in [region_med_coread(r)] if x is not None]
    return round(st.median(v), 1) if v else None
concentration = {
    "structured→collapsed_n": len(collapsed), "structured_kept_n": len(kept),
    "collapsed_cn": cn_mix(collapsed), "kept_cn": cn_mix(kept),
    "collapsed_med_edge_coread": med_cr(collapsed), "kept_med_edge_coread": med_cr(kept),
}
miss_edges = sum(1 for r in regs for e in r.get("nested_edges", []) if edge_survives(r["chrom"], e, 0.02) is None)
tot_edges = sum(len(r.get("nested_edges", [])) for r in regs)

# region 級 flat≥2(零單讀)：edge 存活當 off>=2(非 >coread×eps)
def shape_flat2(r):
    if r["has_cycle"]: return "inconsistent"
    nz = 0
    for e in r.get("nested_edges", []):
        (ca, pa), (cb, pb) = parse(e[0]), parse(e[1])
        rec = lookup.get((ca, min(pa, pb), max(pa, pb)))
        if rec is None or rec["off"] >= 2: nz += 1
    nsib = len(r.get("sibling_pairs", []))
    if nz > 0 and nsib > 0: return "full_tree"
    if nz > 0: return "linear_nested"
    if nsib > 0: return "sibling_only"
    return "collapsed"
flat2 = Counter(shape_flat2(r) for r in regs)

out = {"shape_dist_by_eps": {str(k): v for k, v in dists.items()}, "stability_regions_changed": stability,
       "concentration_collapsed_vs_kept": concentration,
       "flat2_zero_single_read_shape": dict(flat2), "flat2_full_tree": flat2["full_tree"],
       "edge_lookup_miss": miss_edges, "edge_total": tot_edges}
with open(os.path.join(DATA, "region_threshold_impact.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

print(f"=== region 樹形分布 vs 門檻(edge lookup miss {miss_edges}/{tot_edges}) ===")
order = ["full_tree", "linear_nested", "sibling_only", "co_linked_lineage", "no_confirmed_structure", "inconsistent"]
print(f"{'shape':<24}" + "".join(f"{('現行' if e==0 else 'ε'+str(int(e*100))+'%'):>9}" for e in EPS))
for s in order:
    print(f"{s:<24}" + "".join(f"{dists[e].get(s,0):>9}" for e in EPS))
print(f"\n=== 穩定性(換 shape 的區數 / 7143) ===")
for k, v in stability.items(): print(f"  {k}: {v} ({round(100*v/7143,1)}%)")
print(f"\n=== 變動集中度(現行→ε2%,從 full_tree/linear 塌陷者) ===")
for k, v in concentration.items(): print(f"  {k} = {v}")
print("OK wrote region_threshold_impact.json")
