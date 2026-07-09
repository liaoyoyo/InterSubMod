#!/usr/bin/env python3
"""
analyze_mixed_shape_readlik.py — 只對 mixed_shape 區(形狀未定)算 read-likelihood,看能否排出形狀(2026-07-09)。
快版:讀 region-view 的 trees(edges) + mlhp 的 populations_by_hp(read 數),不重跑 solver。
score(tree)=Σ_{edge p→x} log(reads(p)+1),reads(p)=父 genotype 觀測 read 數(隱藏/ROOT=0)→ 走觀測高-read 父的樹更可能。
判:top-likelihood 的樹若集中在單一形狀且明顯>其他形狀 → 形狀 read-可定;否則仍平手(all-hidden/genuine)。
用法: python3 analyze_mixed_shape_readlik.py
"""
import os, sys, json, glob, math
from collections import Counter, defaultdict

MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
MARGIN = 0.5   # top 形狀分數 - 次形狀分數 > 此(log 空間)→ 明顯;否則平手(穩健)

def rv_path(s):
    return f"{PILOT}/layered_region_view_HCC1395.json" if s == "HCC1395" else f"{MSROOT}/{s}/layered_region_view_{s}.json"

def mlhp_lookup(s):
    """(chrom,start) -> {fam: {genotype:count}} 從 mlhp populations_by_hp。"""
    wd = PILOT if s == "HCC1395" else f"{MSROOT}/{s}"
    lk = {}
    for f in sorted(glob.glob(f"{wd}/mlhp_part_*.json")):
        for g in json.load(open(f))["groups"]:
            lk[(g["chrom"], g["start"])] = g.get("populations_by_hp", {}) or {}
    return lk

def tree_shape(edges):
    if not edges:
        return "single"
    ch = defaultdict(list); nodes = set()
    for p, c in edges:
        ch[p].append(c); nodes.add(p); nodes.add(c)
    if len([n for n in nodes if n != "ROOT"]) <= 1:
        return "single"
    cc = {p: len(cs) for p, cs in ch.items()}
    if any(v >= 2 for p, v in cc.items() if p != "ROOT"):
        return "branched_internal"
    if cc.get("ROOT", 0) >= 2:
        return "star_root"
    return "linear"

def tree_score(edges, reads):
    return sum(math.log(reads.get(p, 0) + 1) for p, c in edges if p != "ROOT" and not str(p).startswith("H_"))

def run(s):
    lk = mlhp_lookup(s)
    regions = json.load(open(rv_path(s)))["regions"]
    agg = Counter(); ex = []
    for r in regions:
        for L in r["lineages"]:
            if L["family"] not in ("1", "2") or L.get("capped"):
                continue
            trees = L.get("trees") or []
            shapes = {tree_shape(t["edges"]) for t in trees}
            if len(shapes) <= 1:
                continue   # 只看 mixed_shape(形狀已未定)
            agg["mixed_total"] += 1
            reads = (lk.get((r["chrom"], r["start"]), {}) or {}).get(L["family"], {}) or {}
            # 每形狀的最高 score
            shape_best = {}
            for t in trees:
                sh = tree_shape(t["edges"]); sc = tree_score(t["edges"], reads)
                shape_best[sh] = max(shape_best.get(sh, -1), sc)
            ranked = sorted(shape_best.items(), key=lambda x: -x[1])
            top_sh, top_sc = ranked[0]; sec_sc = ranked[1][1] if len(ranked) > 1 else -1
            if top_sc <= 0:
                agg["tied_all_hidden"] += 1          # 所有父隱藏(read 分不出形狀)
            elif top_sc - sec_sc > MARGIN:
                agg["shape_resolved"] += 1
                agg[f"→{top_sh}"] += 1
                if s == "HCC1395" and len(ex) < 5:
                    ex.append((r["chrom"], r["start"], L["family"], top_sh, round(top_sc - sec_sc, 2)))
            else:
                agg["tied_genuine"] += 1             # 多形狀同分
    return agg, ex

print("=" * 92)
print("mixed_shape(形狀未定)區 read-likelihood 能否排出形狀")
print(f"判:top形狀分數-次形狀 > {MARGIN}(log)→ read-可定;≤0→all-hidden(分不出);其間→genuine 平手")
print("=" * 92)
print(f"{'樣本':15}{'mixed總':>9}{'read可定形狀':>14}{'%':>6}{'all-hidden':>12}{'genuine平手':>12}")
print("-" * 92)
tot = Counter()
for s in SAMPLES:
    agg, ex = run(s)
    mt = agg["mixed_total"]; rs = agg["shape_resolved"]; ah = agg["tied_all_hidden"]; tg = agg["tied_genuine"]
    for k, v in agg.items():
        tot[k] += v
    print(f"{s:15}{mt:>9}{rs:>14}{f'{100*rs/mt:.0f}%' if mt else '-':>6}{ah:>12}{tg:>12}")
    if s == "HCC1395":
        byshape = {k[1:]: v for k, v in agg.items() if k.startswith('→')}
        print(f"    HCC1395 read-可定的形狀分布: {byshape}")
        if ex:
            print(f"    範例: " + " · ".join(f"{c}:{st}(fam{f}→{sh},Δ{d})" for c, st, f, sh, d in ex[:4]))
print("-" * 92)
mt = tot["mixed_total"]; rs = tot["shape_resolved"]
print(f"{'全 7 樣本':15}{mt:>9}{rs:>14}{f'{100*rs/mt:.0f}%' if mt else '-':>6}{tot['tied_all_hidden']:>12}{tot['tied_genuine']:>12}")
print(f"\n裁決:read-likelihood 能對 {100*rs/mt:.0f}% 的 mixed_shape 區排出更可能的形狀;"
      f"{100*tot['tied_all_hidden']/mt:.0f}% 純隱藏(分不出);{100*tot['tied_genuine']/mt:.0f}% 多形同分。")
