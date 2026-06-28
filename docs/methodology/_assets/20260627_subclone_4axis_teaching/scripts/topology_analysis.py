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

# locus -> src(TP/FP, SEQC2 觀察)
locus_src = {}
for c in ("co_linked", "nested_a_in_b", "nested_b_in_a", "mutual_excl", "independent"):
    for hp in ("sameHP", "diffHP"):
        p = os.path.join(DATA, "lists", f"{c}_{hp}.tsv")
        if not os.path.exists(p): continue
        for d in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"):
            locus_src[(d["chrom"], int(d["pos_a"]))] = d.get("src_a", "?")
            locus_src[(d["chrom"], int(d["pos_b"]))] = d.get("src_b", "?")

# hg38 染色體長度(.fai) + centromere 中點(approx) → genome context(telomere/centromere/arm)
CHRLEN = {"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,
"chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,"chr10":133797422,"chr11":135086622,
"chr12":133275309,"chr13":114364328,"chr14":107043718,"chr15":101991189,"chr16":90338345,"chr17":83257441,
"chr18":80373285,"chr19":58617616,"chr20":64444167,"chr21":46709983,"chr22":50818468}
CENTRO = {"chr1":123400000,"chr2":93900000,"chr3":90900000,"chr4":50000000,"chr5":48800000,"chr6":59800000,
"chr7":60100000,"chr8":45200000,"chr9":43000000,"chr10":39800000,"chr11":53400000,"chr12":35500000,
"chr13":17700000,"chr14":17200000,"chr15":19000000,"chr16":36800000,"chr17":25100000,"chr18":18500000,
"chr19":26200000,"chr20":28100000,"chr21":12000000,"chr22":15000000}
def genome_ctx(chrom, start, end, pad=3000000):
    L = CHRLEN.get(chrom); ce = CENTRO.get(chrom)
    if L is None: return "?"
    if start < pad or end > L - pad: return "telomere"
    if ce and not (end < ce - pad or start > ce + pad): return "centromere"
    return "arm"

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
        victim = min(conf, key=lambda g: (cnt[g], g))  # tiebreak: genotype 字串 → 決定性可重現
        dropped += cnt[victim]; del work[victim]
    if len(work) >= 2 and not _laminar(list(work.values())):
        return ("incompatible", [], list(alt), dropped, 0)
    alt = work  # 過濾後
    # 建樹：每節點父=最大的真子集(否則 germline 根);偵測 ambiguous parentage(缺中間群→順序未定)
    nodes = sorted(alt, key=lambda g: (len(alt[g]), g))  # tiebreak: genotype 字串 → 父選擇決定性
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

def dewey_paths(edges, pops, hdigit):
    """每節點(genotype 向量)的 lineage path 標籤 HP{h}-{b1}-{b2}…
    姊妹編號 = **子樹總 read 數遞減**(自己+所有子孫=該 lineage 分支佔比/CCF proxy;
    多者=?-1、少者=?-2);tiebreak: 自己 read 數 → genotype 字串(穩定)。"""
    from collections import defaultdict as dd
    ch = dd(list); haspar = set()
    for p, c in edges:
        ch[p].append(c)
        if p != "ROOT": haspar.add(c)
    cnt = lambda g: pops.get(g, 0)
    _sub = {}
    def subtree(n):  # 自己 + 所有子孫 read 數(樹結構無環,memoized 遞迴)
        if n not in _sub:
            _sub[n] = cnt(n) + sum(subtree(c) for c in ch.get(n, []))
        return _sub[n]
    skey = lambda x: (-subtree(x), -cnt(x), x)  # 子樹總和遞減 → self → 穩定
    paths = {}
    def walk(n, path):
        paths[n] = f"H{hdigit}-" + "-".join(map(str, path))
        for j, k in enumerate(sorted(ch.get(n, []), key=skey), 1):
            if k not in paths: walk(k, path + [j])
    for i, rt in enumerate(sorted(ch.get("ROOT", []), key=skey), 1):
        walk(rt, [i])
    return paths

RI = json.load(open(os.path.join(DATA, "sm_region_integration.json"), encoding="utf-8"))
regs = [r for r in RI["regions"] if r["n_sSNV"] >= 2]

stats = {"n_roots": Counter(), "n_clusters": Counter(), "topology_type": Counter(),
         "determinacy": Counter(), "with_genotype_vectors": 0, "region_coverage": Counter()}
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
    # G1 修(2026-06-29):determinacy 是「單分子建樹」分類,只對有 genotype 向量的區有意義
    # → canonical denominator = with-vector(3885)。無向量區改記 region_coverage(避免 7143/3885 混引)。
    has_vec = len(alt_vecs) >= 1 and len(pops) >= 1
    stats["region_coverage"]["with_genotype_vector" if has_vec else
                              ("germline_only" if ttype == "germline_only" else "no_genotype_vector")] += 1
    # detail：有 genotype 向量(可畫樹)的區
    if has_vec:
        stats["with_genotype_vectors"] += 1
        stats["determinacy"][det] += 1  # G1: 只在此計數 → stats.determinacy == detail 分布(3885)
        tpfp = Counter(locus_src.get((r["chrom"], int(p.split(":")[1])), "?") for p in posset)
        node_hp = {p: hpL.get(locus_hp.get((r["chrom"], int(p.split(":")[1])))) for p in sorted(posset)}
        germ_groot = set(v for v in node_hp.values() if v in ("H1", "H2"))
        germline_reads = pops.get("R" * r["n_sSNV"], 0)
        hap_str = "".join(sorted(set(h for h in hps if h))) or "?"
        hdig = "1" if hap_str == "H1" else ("2" if hap_str == "H2" else "?")
        node_paths = dewey_paths(edges, pops, hdig)
        undefined = (ttype == "incompatible") or (ambig > 0)
        detail.append({"region": r["region"], "chrom": r["chrom"], "start": r["start"], "span": r["span"],
                       "n_sSNV": r["n_sSNV"], "cn": r["cn"],
                       "haplotypes": hap_str, "node_paths": node_paths, "undefined": undefined,
                       "n_roots": len(germ_groot), "germline_reads": germline_reads,
                       "n_clusters": nclust, "topology_type": ttype, "determinacy": det,
                       "has_cycle": bool(r.get("has_cycle")),
                       "cycle_cause": ("CN-gain-multiplicity(同突變多拷貝→pairwise 矛盾)" if r.get("has_cycle") and r["cn"] == "gain"
                                       else ("other-pairwise-cycle" if r.get("has_cycle") else None)),
                       "drop_noise_frac": drop_frac, "ambig_nodes": ambig,
                       "genome_ctx": genome_ctx(r["chrom"], r["start"], r["end"]),
                       "tp": tpfp.get("TP", 0), "fp": tpfp.get("FP", 0),
                       "tree_shape": r["tree_shape"], "populations": pops, "edges": edges,
                       "node_hp": node_hp, "pos_nested": r.get("nested_edges", []),
                       "pos_sibling": r.get("sibling_pairs", []), "pos_vaf": r.get("vaf", {})})

detail.sort(key=lambda d: (-d["n_clusters"], -d["n_sSNV"]))

# chr17 worked example 完整資料(S/r/m 一致標籤)
c17 = json.load(open(os.path.join(DATA, "chr17_subclone_data.json"), encoding="utf-8"))
c17_snvs = sorted(c17["snvs"], key=lambda s: s["pos"])  # 依座標 S1<S2<S3
slabel = {str(s["pos"]): f"S{i+1}" for i, s in enumerate(c17_snvs)}
ALPHA = "48365089"
def sstat(pos):
    st = c17["snv_stat"][str(pos)]; tot = st["tumor_REF"] + st["tumor_ALT"]
    return {"S": slabel[str(pos)], "pos": pos, "change": f'{st["ref"]}>{st["alt"]}',
            "vaf": round(st["tumor_ALT"]/tot, 2) if tot else 0, "hp": "H1",
            "src": locus_src.get(("chr17", pos), "?"),
            "somatic_confirmed": st["somatic_confirmed"],
            "role": "α祖先" if str(pos) == ALPHA else "β後代"}
vp17 = Counter()
ordpos = [str(s["pos"]) for s in c17_snvs]
for rd in c17["reads"]:
    g = rd["geno"]; vp17["".join("A" if g.get(p) == "ALT" else "R" for p in ordpos)] += 1
tt17, ed17, nd17, dr17, am17 = solve_topology(dict(vp17))
tot17 = sum(vp17.values())
chr17_worked = {
    "locus": "chr17:48360161", "genome_ctx": genome_ctx("chr17", 48360161, 48365161),
    "snvs": [sstat(s["pos"]) for s in c17_snvs],
    "populations": [{"vec": k, "reads": v, "pct": round(100*v/tot17, 1),
                     "muts": ",".join(slabel[ordpos[i]] for i, ch in enumerate(k) if ch == "A") or "germline"}
                    for k, v in sorted(vp17.items(), key=lambda x: -x[1])],
    "edges": ed17, "topology_type": tt17, "dropped_noise": dr17,
    "n_cpg": c17["n_cpg"], "n_sig_diff_cpg": c17["n_sig_diff_cpg"],
    "sig_cpg": [{"m": f"m{i+1}", "cpg": x["cpg"], "L1": x["L1"], "L2": x["L2"], "dbeta": x["dbeta"]}
                for i, x in enumerate(c17["sig_diff_cpg"])],
}
out = {"stats": {k: dict(v) if isinstance(v, Counter) else v for k, v in stats.items()},
       "n_detail": len(detail), "detail": detail, "chr17_worked": chr17_worked,
       "chroms": sorted(set(d["chrom"] for d in detail), key=lambda c: int(c[3:]))}
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
