#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
enumerate_candidate_trees(誠實版 / Part B)。
對 B1 ambiguous(缺中間群)區:枚舉每個「一次獲≥2變異」節點的累積順序排列,materialize 成完整候選樹
(插入虛擬中間節點),用 read-based parsimony 算分 + softmax。因中間群未觀察 → 各候選 read_score 相同
→ softmax 均勻 → 標 equiprobable:true(誠實:無 read 證據可排,非在給答案)。
B2 C_underdetermined:不列舉,記「真欠定·非可列舉」。A1 incompatible:跳過(無單樹)。
🔴 排序絕不用甲基(06-28 cis-control SAME-HP double-dip);只用讀數+parsimony,VAF 僅 tiebreak。
來源:topology_per_region.json(+ backbone_resolution.json 取 determinacy/missing)。輸出:candidate_trees.json。
用法: SM_DATA=<sample_dir> python3 enumerate_candidate_trees.py  (預設 HCC1395 4axis)
"""
import os, json, math
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(HERE, "..", "data")))
det = json.load(open(os.path.join(DATA, "topology_per_region.json"), encoding="utf-8"))["detail"]
DET = lambda r: r.get("determinacy", "?")
CAP = 24  # 每區候選樹上限(≤(k+1)!;k_max≈3→24)

def gained(parent_vec, child_vec):
    """child 比 parent 多獲得的 A 位置(parent_vec=None 代表 germline 全 R)。"""
    return [i for i in range(len(child_vec)) if child_vec[i] == "A" and (parent_vec is None or parent_vec[i] != "A")]

def materialize(edges, ambig_choice, glen):
    """把每個 ambiguous 邊 (p,c) 依選定順序 order 展成 p→虛擬→…→c 的鏈。回 (new_edges, virtual_nodes)。"""
    repl = {(p, c): order for (p, c, order) in ambig_choice}
    new_edges, virtual = [], []
    for p, c in edges:
        if (p, c) in repl:
            order = repl[(p, c)]
            pvec = ("R" * glen) if p == "ROOT" else p
            acc = list(pvec); prev = p
            for pos in order[:-1]:               # 除最後一步外皆為「虛擬中間群」
                acc[pos] = "A"
                vname = "".join(acc)
                new_edges.append([prev, vname]); virtual.append(vname); prev = vname
            new_edges.append([prev, c])          # 最後一步接回真實觀察節點 c
        else:
            new_edges.append([p, c])
    return new_edges, virtual

def read_score(edges, pops):
    """read-based parsimony score:觀察節點 reads 多且淺(近根)= 符合。虛擬節點 reads=0 不貢獻。
    🔴 不用甲基。各候選共享同一組觀察節點 → 通常相同(=equiprobable)。"""
    depth = {}
    ch = {}
    for p, c in edges: ch.setdefault(p, []).append(c)
    # BFS depth from ROOT
    stack = [("ROOT", 0)]
    while stack:
        n, d = stack.pop()
        for c in ch.get(n, []):
            depth[c] = d + 1; stack.append((c, d + 1))
    s = 0.0
    for g, reads in pops.items():
        if reads > 0 and g in depth:
            s += (1.0 / depth[g]) * math.log(reads + 1)
    return round(s, 4)

def enumerate_region(r):
    pops = r["populations"]; edges = r.get("edges", [])
    if not edges or not pops: return None
    glen = len(next(iter(pops)))
    ambig = []
    for p, c in edges:
        g = gained(None if p == "ROOT" else p, c)
        if len(g) >= 2: ambig.append((p, c, g))
    if not ambig: return None
    per_node = [(p, c, (list(permutations(g)) if len(g) <= 4 else [tuple(g)])) for (p, c, g) in ambig]
    true_count = 1
    for (p, c, g) in ambig: true_count *= math.factorial(len(g)) if len(g) <= 4 else 1
    # cartesian product(capped)
    combos = [[]]
    for (p, c, orders) in per_node:
        nxt = []
        for combo in combos:
            for o in orders:
                if len(nxt) < CAP: nxt.append(combo + [(p, c, o)])
        combos = nxt
    cands = []
    for combo in combos:
        ne, vn = materialize(edges, combo, glen)
        cands.append({"edges": ne, "virtual_nodes": vn, "read_score": read_score(ne, pops)})
    # softmax over read_score
    scores = [c["read_score"] for c in cands]
    mx = max(scores) if scores else 0
    exps = [math.exp(s - mx) for s in scores]; tot = sum(exps) or 1
    for c, e in zip(cands, exps): c["softmax_prob"] = round(e / tot, 4)
    cands.sort(key=lambda c: -c["softmax_prob"])
    for i, c in enumerate(cands): c["rank"] = i + 1
    equ = (len(set(scores)) <= 1)  # 所有候選 read_score 相同 = 真等機率
    return {"region": r["region"], "base_determinacy": DET(r), "n_sSNV": r["n_sSNV"],
            "ambig_nodes": r.get("ambig_nodes", 0), "missing_intermediates": sum(len(g) - 1 for (_, _, g) in ambig),
            "n_candidates": len(cands), "true_count": true_count, "capped": true_count > CAP,
            "equiprobable": equ, "candidate_set": cands,
            "honest_note": ("等機率·無 read 證據可排(中間群未觀察)" if equ else "read_score 有差異但弱"),
            "truncated": bool(r.get("truncated")) or r["n_sSNV"] > glen}

b1, b2_nonenum, capped_n, equ_n = [], 0, 0, 0
for r in det:
    d = DET(r)
    if "ambiguous" in d:
        res = enumerate_region(r)
        if res:
            b1.append(res)
            if res["capped"]: capped_n += 1
            if res["equiprobable"]: equ_n += 1
    elif d == "C_underdetermined":
        b2_nonenum += 1

out = {
    "schema_version": "20260630",
    "summary": {
        "B1_enumerated": len(b1),
        "B1_equiprobable": equ_n,
        "B1_capped(>%d)" % CAP: capped_n,
        "candidates_total": sum(x["n_candidates"] for x in b1),
        "B2_underdetermined_nonenumerable": b2_nonenum,
        "note": "B1=缺中間群順序排列(materialize 完整候選樹,虛擬中間節點);因中間群未觀察→等機率,非 ranking。"
                "B2=單一ALT群/缺對→真欠定非可列舉。排序只用讀數+parsimony,🔴絕不用甲基(06-28)。"},
    "candidate_trees": b1,
    "B2_underdetermined": {"n": b2_nonenum, "status": "真欠定·單一ALT群或缺對·需深覆蓋·非可列舉"},
}
with open(os.path.join(DATA, "candidate_trees.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False)
print("ENUMERATE DONE →", os.path.join(DATA, "candidate_trees.json"))
print(json.dumps(out["summary"], ensure_ascii=False, indent=1))
