#!/usr/bin/env python3
"""[S3 演化重建 — 從 per-pair 連鎖地圖建 local clonal forest]

讀 sm_linkage_genomewide.json (pairs + census)。只取 confirmed 邊 (powered=coread>=COREAD_MIN & 兩端 somatic)。
同-HP 邊 = subclone 相關；異-HP = allelic(另記)。
- co_linked  → union-find 併同 lineage 節點。
- nested     → 有向邊 祖先→後代。 nested_a_in_b: a⊂b → b 祖先→a；nested_b_in_a: a 祖先→b。
- mutual_excl(same_hp) → sibling 對 (平行分支, 無邊)。
每個弱連通分量(local region):transitive reduction + cycle 偵測 + VAF 梯度檢查(Foltz: 祖先 VAF >= 後代)。
輸出 sm_evolution.json: per-component 樹拓樸 + aggregate。數字落檔, 不手打。
"""
import json
from collections import defaultdict

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
IN = f"{A}/sm_linkage_genomewide.json"
OUT = f"{A}/sm_evolution.json"


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
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def transitive_reduction(nodes, edges):
    """edges = set of (u,v) directed (祖先->後代). 移除可由其他路徑達成的直接邊。回 reduced set + has_cycle。"""
    succ = defaultdict(set)
    for u, v in edges:
        succ[u].add(v)
    # reachability (excluding direct edge) via DFS
    def reach(start, avoid):
        seen = set()
        stack = [w for w in succ[start] if (start, w) != avoid]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(succ[x])
        return seen
    reduced = set()
    for u, v in edges:
        if v not in reach(u, (u, v)):
            reduced.add((u, v))
    # cycle detection (Kahn)
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
    has_cycle = seen < len(nodes)
    return reduced, has_cycle


def main():
    d = json.load(open(IN))
    census = d["census"]
    pairs = d["pairs"]
    cm = d["params"]["COREAD_MIN"]

    # confirmed same-HP edges, grouped by chrom
    by_chrom = defaultdict(list)
    allelic = 0
    for p in pairs:
        if not (p["powered"] and p["somatic_a"] and p["somatic_b"]):
            continue
        if not p["same_hp"]:
            if p["rel"] == "mutual_excl":
                allelic += 1
            continue
        by_chrom[p["chrom"]].append(p)

    components = []
    for chrom, ps in by_chrom.items():
        # union-find: co_linked merges; all edges connect for weak components
        uf = UF()
        for p in ps:
            a = f"{chrom}:{p['a']}"
            b = f"{chrom}:{p['b']}"
            uf.find(a); uf.find(b)
            if p["rel"] == "co_linked":
                uf.union(a, b)
        # build per weakly-connected component (connect any related pair)
        for p in ps:
            uf.union(f"{chrom}:{p['a']}", f"{chrom}:{p['b']}")  # weak connectivity
        # NOTE: weak union above merges co_linked AND links any related nodes into one component;
        # to keep co_linked as same NODE but others as separate nodes, track node-merge separately:
        node_uf = UF()
        for p in ps:
            if p["rel"] == "co_linked":
                node_uf.union(f"{chrom}:{p['a']}", f"{chrom}:{p['b']}")

        comp_members = defaultdict(set)
        for p in ps:
            for pos in (p["a"], p["b"]):
                key = f"{chrom}:{pos}"
                comp_members[uf.find(key)].add(key)

        for croot, members in comp_members.items():
            # nodes = co_linked-merged groups
            node_of = {m: node_uf.find(m) for m in members}
            nodes = set(node_of.values())
            nested_edges = set()
            sibling = set()
            for p in ps:
                ka, kb = f"{chrom}:{p['a']}", f"{chrom}:{p['b']}"
                if ka not in members:
                    continue
                na, nb = node_of[ka], node_of[kb]
                if na == nb:
                    continue
                if p["rel"] == "nested_a_in_b":
                    nested_edges.add((nb, na))  # b 祖先 -> a 後代
                elif p["rel"] == "nested_b_in_a":
                    nested_edges.add((na, nb))  # a 祖先 -> b 後代
                elif p["rel"] == "mutual_excl":
                    sibling.add(tuple(sorted((na, nb))))
            reduced, has_cycle = transitive_reduction(nodes, nested_edges)
            # VAF per node (max vaf among merged sSNV)
            vaf = {}
            for m in members:
                n = node_of[m]
                v = census[m].get("vaf")
                if v is not None:
                    vaf[n] = max(vaf.get(n, 0.0), v)
            # VAF gradient: 祖先 vaf >= 後代 vaf ?
            grad_viol = 0
            for u, v in reduced:
                if u in vaf and v in vaf and vaf[u] + 0.05 < vaf[v]:
                    grad_viol += 1
            # depth
            succ = defaultdict(set)
            for u, v in reduced:
                succ[u].add(v)

            def depth(n, seen):
                if n in seen:
                    return 0
                seen = seen | {n}
                return 1 + max([depth(c, seen) for c in succ[n]], default=0)
            roots = nodes - {v for _, v in reduced}
            maxdepth = max([depth(r, set()) for r in roots], default=0)
            components.append({
                "chrom": chrom, "n_sSNV": len(members), "n_nodes": len(nodes),
                "members": sorted(members), "node_of": node_of,
                "nested_edges": sorted(list(reduced)), "n_nested": len(reduced),
                "sibling_pairs": sorted(list(sibling)), "n_sibling": len(sibling),
                "has_cycle": has_cycle, "vaf_gradient_violations": grad_viol,
                "max_depth": maxdepth, "vaf": vaf,
                "is_tree": (len(reduced) >= 1 and not has_cycle),
            })

    # aggregate
    agg = {
        "n_components": len(components),
        "n_with_nested": sum(1 for c in components if c["n_nested"] >= 1),
        "n_with_sibling": sum(1 for c in components if c["n_sibling"] >= 1),
        "n_multilevel_tree": sum(1 for c in components if c["max_depth"] >= 2),
        "n_sibling_and_nested": sum(1 for c in components if c["n_nested"] >= 1 and c["n_sibling"] >= 1),
        "n_with_cycle": sum(1 for c in components if c["has_cycle"]),
        "n_vaf_violations": sum(c["vaf_gradient_violations"] for c in components),
        "allelic_diffHP_excluded": allelic,
    }
    out = {"params": d["params"], "aggregate": agg,
           "components": sorted(components, key=lambda c: (-c["n_nodes"], -c["max_depth"]))}
    json.dump(out, open(OUT, "w"), ensure_ascii=False)
    print(f"DONE components={agg['n_components']} with_nested={agg['n_with_nested']} "
          f"with_sibling={agg['n_with_sibling']} multilevel={agg['n_multilevel_tree']} "
          f"sib+nested={agg['n_sibling_and_nested']} cycles={agg['n_with_cycle']} -> {OUT}")


if __name__ == "__main__":
    main()
