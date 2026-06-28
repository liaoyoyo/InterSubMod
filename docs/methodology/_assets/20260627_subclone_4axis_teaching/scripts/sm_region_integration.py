"""[per-region 整合 — 最終綜合] 單位 = 最大可關聯區域 (≤50kb chain group)。
每區域整合: 樹狀結構(nested 祖先→後代 + sibling 分支 + co_linked 併 lineage) + multi-locus populations
+ CN context + clean 狀態 → 樹形分類。genome-wide 比例 + per-region TSV，供之後驗證/查詢。
依賴: sm_linkage_genomewide.json(pairs+census) + ml_part_*.json(multilocus populations) + CN bed。
"""
import json
import glob
import os
from collections import defaultdict, Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
CNBED = "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed"
LD = f"{A}/lists"
os.makedirs(LD, exist_ok=True)
TIER_R = 50000


def load_cn():
    iv = defaultdict(list)
    for ln in open(CNBED):
        p = ln.split()
        if len(p) < 4 or not p[1].isdigit():
            continue
        iv[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in iv:
        iv[c].sort()
    return iv


def cn_state(cn, chrom, pos):
    for s, e, lab in cn.get(chrom, []):
        if s <= pos <= e:
            return lab
    return "neutral"


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


def transitive_reduction(nodes, edges):
    succ = defaultdict(set)
    for u, v in edges:
        succ[u].add(v)

    def reach(start, avoid):
        seen = set()
        st = [w for w in succ[start] if (start, w) != avoid]
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x)
            st.extend(succ[x])
        return seen
    reduced = set((u, v) for u, v in edges if v not in reach(u, (u, v)))
    indeg = defaultdict(int)
    adj = defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        indeg[v] += 1
    q = [n for n in nodes if indeg[n] == 0]
    seen = 0
    while q:
        x = q.pop()
        seen += 1
        for w in adj[x]:
            indeg[w] -= 1
            if indeg[w] == 0:
                q.append(w)
    return reduced, seen < len(nodes)


def build_tree(chrom, members, pairs_idx, census):
    """members=set of pos. 從 confirmed same-HP pairs 建樹。回 tree dict。"""
    node_uf = UF()
    for m in members:
        node_uf.find(f"{chrom}:{m}")
    rel_pairs = []
    for p in pairs_idx:
        if p["a"] in members and p["b"] in members:
            rel_pairs.append(p)
            if p["rel"] == "co_linked":
                node_uf.union(f"{chrom}:{p['a']}", f"{chrom}:{p['b']}")
    node_of = {f"{chrom}:{m}": node_uf.find(f"{chrom}:{m}") for m in members}
    nodes = set(node_of.values())
    nested = set()
    sibling = set()
    for p in rel_pairs:
        na, nb = node_of[f"{chrom}:{p['a']}"], node_of[f"{chrom}:{p['b']}"]
        if na == nb:
            continue
        if p["rel"] == "nested_a_in_b":
            nested.add((nb, na))   # b 祖先 -> a 後代
        elif p["rel"] == "nested_b_in_a":
            nested.add((na, nb))   # a 祖先 -> b 後代
        elif p["rel"] == "mutual_excl":
            sibling.add(tuple(sorted((na, nb))))
    reduced, cyc = transitive_reduction(nodes, nested)
    vaf = {}
    for m in members:
        v = census[f"{chrom}:{m}"].get("vaf")
        if v is not None:
            n = node_of[f"{chrom}:{m}"]
            vaf[n] = max(vaf.get(n, 0.0), v)
    succ = defaultdict(set)
    for u, v in reduced:
        succ[u].add(v)

    def depth(n, seen):
        return 0 if n in seen else 1 + max([depth(c, seen | {n}) for c in succ[n]], default=0)
    roots = nodes - {v for _, v in reduced}
    md = max([depth(r, set()) for r in roots], default=0)
    return {"n_nodes": len(nodes), "n_nested": len(reduced), "n_sibling": len(sibling),
            "nested_edges": sorted(list(reduced)), "sibling_pairs": sorted(list(sibling)),
            "max_depth": md, "has_cycle": cyc, "vaf": {k: round(v, 2) for k, v in vaf.items()},
            "n_colinked_merges": len(members) - len(nodes)}


def classify_shape(t):
    if t["has_cycle"]:
        return "inconsistent"
    if t["n_nested"] >= 1 and t["n_sibling"] >= 1:
        return "full_tree"          # 分支+深度 (chr17 式)
    if t["n_nested"] >= 1:
        return "linear_nested" if t["max_depth"] >= 2 else "single_nested"
    if t["n_sibling"] >= 1:
        return "sibling_only"       # 平行分支
    if t["n_colinked_merges"] >= 1:
        return "co_linked_lineage"  # 單一 lineage 無分支
    return "no_confirmed_structure"


def main():
    d = json.load(open(f"{A}/sm_linkage_genomewide.json"))
    census = d["census"]
    cn = load_cn()
    # index confirmed same-HP+somatic pairs by chrom
    pairs_by_chrom = defaultdict(list)
    for p in d["pairs"]:
        if p["powered"] and p["somatic_a"] and p["somatic_b"] and p["same_hp"]:
            pairs_by_chrom[p["chrom"]].append(p)
    # multilocus populations: index by chrom -> list of regions
    ml = defaultdict(list)
    for f in glob.glob(f"{A}/ml_part_*.json"):
        for g in json.load(open(f))["groups"]:
            ml[g["chrom"]].append(g)

    def ml_match(chrom, lo, hi):
        for g in ml.get(chrom, []):
            if lo <= g["start"] <= hi or lo <= g["end"] <= hi:
                return g
        return None

    # build regions: per chrom, somatic sSNV chained <=50kb
    regions = []
    shape_cnt = Counter()
    shape_clean = Counter()
    for chrom in [f"chr{c}" for c in range(1, 23)]:
        som = sorted(int(k.split(":")[1]) for k, c in census.items()
                     if c["chrom"] == chrom and c.get("somatic") is True)
        if not som:
            continue
        cur = [som[0]]
        groups = []
        for x in som[1:]:
            if x - cur[-1] > TIER_R:
                groups.append(cur)
                cur = []
            cur.append(x)
        groups.append(cur)
        for g in groups:
            if len(g) < 2:
                continue
            members = set(g)
            tree = build_tree(chrom, members, pairs_by_chrom[chrom], census)
            shape = classify_shape(tree)
            cst = cn_state(cn, chrom, g[len(g) // 2])
            mlrec = ml_match(chrom, g[0], g[-1])
            rec = {"region": f"{chrom}:{g[0]}-{g[-1]}", "chrom": chrom, "start": g[0], "end": g[-1],
                   "span": g[-1] - g[0], "n_sSNV": len(g), "cn": cst,
                   "tree_shape": shape, "n_nodes": tree["n_nodes"], "n_nested": tree["n_nested"],
                   "n_sibling": tree["n_sibling"], "max_depth": tree["max_depth"],
                   "has_cycle": tree["has_cycle"], "n_colinked_merges": tree["n_colinked_merges"],
                   "nested_edges": tree["nested_edges"], "sibling_pairs": tree["sibling_pairs"],
                   "vaf": tree["vaf"],
                   "n_populations": mlrec["n_populations"] if mlrec else None,
                   "populations": mlrec["populations"] if mlrec else None}
            regions.append(rec)
            shape_cnt[shape] += 1
            if cst in ("loh", "neutral"):
                shape_clean[shape] += 1

    n = len(regions)
    agg = {
        "n_regions": n,
        "shape_distribution": dict(shape_cnt),
        "shape_proportion": {k: round(v / n, 3) for k, v in shape_cnt.items()},
        "shape_distribution_clean_loh_neutral": dict(shape_clean),
        "n_branching_full_tree": shape_cnt["full_tree"],
        "n_clean_full_tree_loh_neutral": shape_clean["full_tree"],
        "by_cn": dict(Counter(r["cn"] for r in regions)),
        "by_n_populations": dict(Counter(r["n_populations"] for r in regions if r["n_populations"])),
    }
    out = {"unit": "maximal linkage region (≤50kb somatic-sSNV chain)", "aggregate": agg,
           "regions": sorted(regions, key=lambda r: (-r["n_nested"] - r["n_sibling"], -r["n_sSNV"]))}
    json.dump(out, open(f"{A}/sm_region_integration.json", "w"), ensure_ascii=False)

    # per-region TSV (queryable)
    with open(f"{LD}/regions.tsv", "w") as fh:
        fh.write("region\tchrom\tstart\tend\tspan\tn_sSNV\tcn\ttree_shape\tn_nodes\tn_nested\t"
                 "n_sibling\tmax_depth\thas_cycle\tn_populations\n")
        for r in out["regions"]:
            fh.write(f'{r["region"]}\t{r["chrom"]}\t{r["start"]}\t{r["end"]}\t{r["span"]}\t{r["n_sSNV"]}\t'
                     f'{r["cn"]}\t{r["tree_shape"]}\t{r["n_nodes"]}\t{r["n_nested"]}\t{r["n_sibling"]}\t'
                     f'{r["max_depth"]}\t{r["has_cycle"]}\t{r["n_populations"]}\n')

    print(f"REGIONS n={n}")
    print("shape distribution (all / clean loh+neutral):")
    for k in sorted(shape_cnt, key=lambda x: -shape_cnt[x]):
        print(f"  {k:24s} {shape_cnt[k]:>6} ({100*shape_cnt[k]/n:.1f}%)  clean={shape_clean.get(k,0)}")
    print(f"full_tree (branching+depth): {shape_cnt['full_tree']} all, {shape_clean['full_tree']} clean")
    print(f"-> sm_region_integration.json + lists/regions.tsv")


if __name__ == "__main__":
    main()
