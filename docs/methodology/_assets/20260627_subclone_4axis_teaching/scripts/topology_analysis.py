#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
拓樸分析（cluster-count-first 算法落到實際資料）：
  每區 → ① n_clusters(observed genotype 向量) → ② HP 定根 → ③ 從向量的突變集 laminar 關係解樹
  → ④ topology_type(linear/branched/star/mixed/incompatible) → ⑤ determinacy + 候選數
理論：perfect-phylogeny → 突變集須巢狀或互斥(laminar);觀察 c 群 → 拓樸限縮到 c-節點樹。
輸出：topology_per_region.json(detail 有 genotype 向量的區 + 全域 stats)；chr17 驗證。
用法：python3 topology_analysis.py
"""
import json, os, csv
from collections import Counter, defaultdict
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

# locus -> HP
locus_hp = {}
for c in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl", "independent"):
    for hp in ("sameHP", "diffHP"):
        p = os.path.join(DATA, "lists", f"{c}_{hp}.tsv")
        if not os.path.exists(p): continue
        for d in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"):
            locus_hp[(d["chrom"], int(d["pos_a"]))] = d["hp_a"]
            locus_hp[(d["chrom"], int(d["pos_b"]))] = d["hp_b"]
hpL = {"1-1": "H1", "2-1": "H2", "3": "H3?"}

def altset(g): return frozenset(i for i, ch in enumerate(g) if ch == "A")

def _laminar(sets):
    for a, b in combinations(sets, 2):
        if a & b and not (a <= b or b <= a): return False
    return True

def solve_topology(pops):
    """pops={geno:count}. 回 (type, edges, nodes, dropped_reads)。
    population 噪聲過濾：貪婪移除最小衝突群直到 laminar(perfect-phylogeny);記 dropped。"""
    alt = {g: altset(g) for g in pops if "A" in g}
    if not alt: return ("germline_only", [], list(pops), 0, 0)
    total = sum(pops.values()) or 1
    work = dict(alt); cnt = {g: pops[g] for g in alt}; dropped = 0
    while len(work) >= 2 and not _laminar(list(work.values())):
        conf = set()
        for g1, g2 in combinations(work, 2):
            a, b = work[g1], work[g2]
            if a & b and not (a <= b or b <= a): conf.add(g1); conf.add(g2)
        if not conf: break
        victim = min(conf, key=lambda g: cnt[g])
        dropped += cnt[victim]; del work[victim]
    if len(work) >= 2 and not _laminar(list(work.values())):
        return ("incompatible", [], list(alt), dropped, 0)
    alt = work  # 過濾後
    # 建樹：每節點父=最大的真子集(否則 germline 根);偵測 ambiguous parentage(缺中間群→順序未定)
    nodes = sorted(alt, key=lambda g: len(alt[g]))
    edges = []; ambig = 0
    for g in nodes:
        cand = [h for h in nodes if alt[h] < alt[g]]
        if cand:
            msz = max(len(alt[h]) for h in cand)
            mcand = [h for h in cand if len(alt[h]) == msz]
            parent = mcand[0]
            # 跳 >1 突變(缺中間群) 或 多個同 size 父 → 內部順序未定
            if (len(alt[g]) - msz) > 1 or len(mcand) > 1: ambig += 1
        else:
            parent = "ROOT"
        edges.append((parent, g))
    # 分型
    depths = {}
    def dep(g):
        if g == "ROOT": return 0
        par = [p for p, c in edges if c == g][0]
        return 1 + dep(par)
    maxd = max(dep(g) for g in nodes)
    roots = [c for p, c in edges if p == "ROOT"]
    has_branch = any(len([c for p, c in edges if p == par]) >= 2 for par in set(p for p, c in edges))
    if len(nodes) == 1: t = "single"
    elif len(roots) >= 2 and not has_branch: t = "star(全姊妹)"
    elif maxd == len(nodes) and not has_branch: t = "linear(全直系)"
    elif has_branch: t = "branched(直系+姊妹)"
    else: t = "mixed"
    return (t, edges, nodes, dropped, ambig)

RI = json.load(open(os.path.join(DATA, "sm_region_integration.json"), encoding="utf-8"))
regs = [r for r in RI["regions"] if r["n_sSNV"] >= 2]

stats = {"n_roots": Counter(), "n_clusters": Counter(), "topology_type": Counter(),
         "determinacy": Counter(), "with_genotype_vectors": 0}
detail = []
for r in regs:
    pops = r.get("populations") or {}
    # HP 根
    posset = set()
    for e in r.get("nested_edges", []): posset.add(e[0]); posset.add(e[1])
    for e in r.get("sibling_pairs", []): posset.add(e[0]); posset.add(e[1])
    hps = Counter(hpL.get(locus_hp.get((p.split(":")[0], int(p.split(":")[1])))) for p in posset)
    germ = {h for h in hps if h in ("H1", "H2")}
    n_roots = len(germ) + (1 if "H3?" in hps and not germ else 0)
    stats["n_roots"][f"{len(germ)}根" + ("+HP3" if "H3?" in hps else "")] += 1
    # cluster + topology（用 genotype 向量）
    alt_vecs = [g for g in pops if "A" in g]
    nclust = len(alt_vecs)
    if pops:
        ttype, edges, nodes, dropped, ambig = solve_topology(pops)
    else:
        ttype, edges, nodes, dropped, ambig = ("no_genotype_vectors", [], [], 0, 0)
    drop_frac = round(dropped / (sum(pops.values()) or 1), 3)
    stats["n_clusters"][nclust] += 1
    stats["topology_type"][ttype] += 1
    # determinacy
    if r["has_cycle"]: det = "incompatible"
    elif len(alt_vecs) >= 2 and sum(pops.values()) >= 6: det = "A_determined(單分子向量)"
    elif r["tree_shape"] in ("full_tree", "linear_nested", "sibling_only"): det = "B_pairwise_structure"
    elif r["tree_shape"] == "no_confirmed_structure": det = "C_underdetermined"
    else: det = "other"
    if ambig>0 and det.startswith('A'): det='A_ambiguous_order(缺中間群)'
    stats["determinacy"][det] += 1
    # detail：有 genotype 向量(可畫樹)的區
    if len(alt_vecs) >= 1 and len(pops) >= 1:
        stats["with_genotype_vectors"] += 1
        detail.append({"region": r["region"], "chrom": r["chrom"], "span": r["span"], "n_sSNV": r["n_sSNV"],
                       "cn": r["cn"], "haplotypes": "".join(sorted(set(h for h in hps if h))) or "?",
                       "n_clusters": nclust, "topology_type": ttype, "determinacy": det,
                       "drop_noise_frac": drop_frac, "ambig_nodes": ambig,
                       "tree_shape": r["tree_shape"], "populations": pops, "edges": edges})

detail.sort(key=lambda d: (-d["n_clusters"], -d["n_sSNV"]))
out = {"stats": {k: dict(v) if isinstance(v, Counter) else v for k, v in stats.items()},
       "n_detail": len(detail), "detail": detail}
with open(os.path.join(DATA, "topology_per_region.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False)

print(f"=== 全域 stats (regions ≥2 sSNV = {len(regs)}) ===")
print(f"n_roots: {dict(stats['n_roots'])}")
print(f"n_clusters: {dict(sorted(stats['n_clusters'].items()))}")
print(f"topology_type: {dict(stats['topology_type'])}")
print(f"determinacy: {dict(stats['determinacy'])}")
print(f"有 genotype 向量可畫樹的區(detail): {len(detail)}")
# chr17 驗證(per-read 向量)
c17 = json.load(open(os.path.join(DATA, "chr17_subclone_data.json"), encoding="utf-8"))
A,B1,B2="48365089","48362515","48365161"
vp=Counter()
for rd in c17["reads"]:
    g=rd["geno"]; s=("A" if g.get(B1)=="ALT" else "R")+("A" if g.get(A)=="ALT" else "R")+("A" if g.get(B2)=="ALT" else "R")
    vp[s]+=1
t,e,n,dr,amb=solve_topology(dict(vp))
print(f"=== chr17 驗證: 向量{dict(vp)} → topology={t}, dropped={dr}, ambig={amb}, edges={e} (應 linear: α祖先→β後代) ===")
print("OK wrote topology_per_region.json")
