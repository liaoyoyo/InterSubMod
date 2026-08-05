#!/usr/bin/env python3
"""
analyze_read_likelihood.py — 量化:對 ambiguous 區加 read-likelihood 能排出多少 top 樹 vs 平手(2026-07-09)。
方法(非循環,遺傳):每棵枚舉樹 score = Σ_{edge p→x} log(reads(p)+1),reads(p)=父 genotype 觀測 read 數(隱藏/ROOT=0)。
  → 走「觀測到的高-read 中間節點」的樹 score 高 = 更可能。unique top(有 margin)=read-resolvable;
  平手分兩類:all-hidden(所有父都隱藏 0-read,如 {11}→隱藏 10/01)/genuine-tie(多樹同觀測父分數)。
判準:resolvable = 恰 1 樹 score 為最高 且 最高>0 且 (次高/最高)<TIE_RATIO(穩健,不過度解讀噪聲)。
用法: python3 analyze_read_likelihood.py   (7 樣本;HCC1395 出細節)
"""
import os, sys, json, glob, math
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tree_enumeration_solver as S

C = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
PILOT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
TIE_RATIO = 0.90     # 次高/最高 ≥ 此值 → 視為平手(穩健);<此值 → top 明顯勝出
MIN_TOP = 1.0        # 最高 score 需 > 此(至少一觀測父有 read)才算 resolvable

def mlhp_paths(s):
    wd = PILOT if s == "HCC1395" else f"{MSROOT}/{s}"
    return sorted(glob.glob(f"{wd}/mlhp_part_*.json"))

def tree_score(edges, reads):
    """Σ log(reads(parent)+1);parent=觀測 genotype(非 ROOT/H_)才算。"""
    sc = 0.0
    for p, c in edges:
        if p != "ROOT" and not str(p).startswith("H_"):
            sc += math.log(reads.get(p, 0) + 1)
    return sc

def analyze_unit(full, part, k):
    """回 (class, ntrees, verdict, topgap)。verdict∈{resolvable,tied_all_hidden,tied_genuine,capped,not_ambiguous}。"""
    res = S.enumerate_min_trees(full, part, k)
    if res["capped"]:
        return ("capped", res["n_trees"], "capped", None)
    if res["n_trees"] <= 1:
        return (S.classify(res), res["n_trees"], "not_ambiguous", None)
    # ambiguous: reads = 觀測 genotype→count(full-pop);partial-only 則 full={} → 全隱藏父
    reads = dict(full)
    scores = sorted((tree_score(t["edges"], reads) for t in res["trees"]), reverse=True)
    top = scores[0]; second = scores[1] if len(scores) > 1 else 0.0
    if top <= 0:
        return ("ambiguous", res["n_trees"], "tied_all_hidden", 0.0)     # 所有父隱藏(read 分不出)
    ratio = second / top if top > 0 else 1.0
    if ratio < TIE_RATIO and top > MIN_TOP:
        return ("ambiguous", res["n_trees"], "resolvable", round(1 - ratio, 3))  # top 明顯勝
    return ("ambiguous", res["n_trees"], "tied_genuine", round(ratio, 3))        # 多樹同分

def run_sample(s, detail=False):
    groups = []
    for f in mlhp_paths(s):
        groups += json.load(open(f))["groups"]
    agg = Counter(); det_examples = []
    for g in groups:
        if g.get("n_sSNV", 0) < 2:
            continue
        pbh = g.get("populations_by_hp", {}) or {}; sbh = g.get("subread_groups_by_hp", {}) or {}
        for fam in ("1", "2"):   # germline lineage
            full = pbh.get(fam, {}) or {}; part = list((sbh.get(fam, {}) or {}).keys())
            if not full and not part:
                continue
            k = len(next(iter(full))) if full else len(part[0])
            cls, nt, verdict, gap = analyze_unit(full, part, k)
            if verdict in ("capped", "not_ambiguous"):
                continue
            has_full = len(full) > 0
            key = verdict + ("" if has_full else "_partialonly")
            agg[key] += 1
            if detail and verdict == "resolvable" and len(det_examples) < 5:
                det_examples.append((g["chrom"], g["start"], fam, nt, gap))
    return agg, det_examples

print("=" * 96)
print("read-likelihood 量化:ambiguous 區能否用 read 數排出更可能的樹")
print(f"判準: score=Σlog(父觀測read+1); resolvable=unique top 且 次高/最高<{TIE_RATIO}; 平手分 all-hidden/genuine")
print("=" * 96)
print(f"{'樣本':15}{'ambig總':>8}{'resolvable':>12}{'%':>6}{'tied_all_hidden':>16}{'tied_genuine':>13}{'partial-only 平手':>16}")
print("-" * 96)
tot = Counter()
for s in SAMPLES:
    agg, ex = run_sample(s, detail=(s == "HCC1395"))
    resolv = agg["resolvable"] + agg["resolvable_partialonly"]
    tah = agg["tied_all_hidden"] + agg["tied_all_hidden_partialonly"]
    tg = agg["tied_genuine"]
    po = agg["tied_all_hidden_partialonly"] + agg["tied_genuine_partialonly"] + agg["resolvable_partialonly"]
    ambig = resolv + tah + tg + agg["tied_genuine_partialonly"]
    for k, v in agg.items():
        tot[k] += v
    pctr = 100 * resolv / ambig if ambig else 0
    print(f"{s:15}{ambig:>8}{resolv:>12}{f'{pctr:.0f}%':>6}{tah:>16}{tg:>13}{po:>16}")
    if s == "HCC1395" and ex:
        print(f"    HCC1395 resolvable 範例: " + " · ".join(f"{c}:{st}(fam{fm},{n}樹,gap{gp})" for c, st, fm, n, gp in ex[:4]))
print("-" * 96)
tresolv = tot["resolvable"] + tot["resolvable_partialonly"]
tah = tot["tied_all_hidden"] + tot["tied_all_hidden_partialonly"]
tg = tot["tied_genuine"] + tot["tied_genuine_partialonly"]
tambig = tresolv + tah + tg
print(f"{'全 7 樣本':15}{tambig:>8}{tresolv:>12}{f'{100*tresolv/tambig:.0f}%' if tambig else '-':>6}{tah:>16}{tg:>13}")
print(f"\n裁決:read-likelihood 能對 {100*tresolv/tambig:.0f}% 的 ambiguous 區排出明顯 top 樹;"
      f"{100*tah/tambig:.0f}% 純隱藏祖先(read 分不出);{100*tg/tambig:.0f}% 觀測父同分平手。")
