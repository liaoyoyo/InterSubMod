#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
enumerate_candidate_trees(誠實版 / Part B) — 2026-07-04 gained-pair pairwise 定序修正。

原盲點(用戶 catch):對「缺中間群」jump 邊,只用『全跨 read 的 population(節點)』判 read_score →
各候選相同 → 全 equiprobable。但**沒查那兩個缺失中間群突變之間的 pairwise 2×2**(該 2×2 已含定序資訊)。

修正:枚舉排列『之前』先查 gained 突變間 pairwise(sm_linkage pairs):
  co_linked(RA=AR=0)     → 併成 block group(一起獲得·無中間群·不排序)
  nested(小邊<ε=2%)       → 定序約束(中間群被觀測·方向由大邊定)
  4-gamete(雙邊≥ε)        → CONFLICT(flag)
  無直接 coread           → 才是真等機率 → 保留排列枚舉
→ 用觀察約束 prune 排列集(linear extensions of the partial order);block/ordered → 收斂單一候選、
  equiprobable=false;唯真無 coread 仍 equiprobable。ε 對齊 ONT 錯誤率;方向再用 pos_vaf 交叉核對。
🔴 排序絕不用甲基(06-28 cis-control SAME-HP double-dip)。B2/跨-PS non-identifiability 仍成立(不涉)。

來源:topology_per_region.json + sm_linkage_genomewide.json(pairs+census;SM_DATA 內或 {M.A})。
輸出:candidate_trees.json(backward-compat 欄位 + 每區 resolution + 每邊 edge_resolution)。
用法: SM_DATA=<sample_dir> python3 enumerate_candidate_trees.py
"""
import os, sys, json, math
from itertools import permutations, combinations
from bisect import bisect_left, bisect_right
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SM_DATA", os.path.normpath(os.path.join(HERE, "..", "data")))
det = json.load(open(os.path.join(DATA, "topology_per_region.json"), encoding="utf-8"))["detail"]
DET = lambda r: r.get("determinacy", "?")
CAP = 24        # 每區候選樹上限
EPS = 0.02      # ONT 噪聲底線(小邊<EPS*coread=噪聲)

# --- pairs + census(gained-pair pairwise 定序用) ---
_pp = os.path.join(DATA, "sm_linkage_genomewide.json")
if not os.path.exists(_pp):
    sys.path.insert(0, HERE)
    import sm_linkage_genomewide as M  # noqa: E402
    _pp = os.path.join(M.A, "sm_linkage_genomewide.json")
_HAVE_PAIRS = os.path.exists(_pp)
plut, som_by_chrom = {}, defaultdict(list)
if _HAVE_PAIRS:
    _L = json.load(open(_pp, encoding="utf-8"))
    for p in _L["pairs"]:
        plut[(p["chrom"], min(p["a"], p["b"]), max(p["a"], p["b"]))] = p
    for key, c in _L["census"].items():
        if c.get("somatic") is True:
            cc, pos = key.rsplit(":", 1)
            som_by_chrom[cc].append(int(pos))
    for cc in som_by_chrom:
        som_by_chrom[cc].sort()


def positions_for(r, glen):
    """vector index → sSNV 位點(census somatic in span;truncated 取 densest-8)。"""
    arr = som_by_chrom.get(r["chrom"], [])
    start, end = r["start"], r["start"] + r["span"]
    s = [p for p in arr[bisect_left(arr, start):bisect_right(arr, end)] if start <= p <= end]
    if len(s) == glen:
        return s
    if len(s) > glen == 8:
        best = min(range(len(s) - 8 + 1), key=lambda i: s[i + 7] - s[i])
        return s[best:best + 8]
    return None


def _geno(lbl):  # gap#2: 正規化節點 label(build_hidden_node_tree 的隱藏節點是 'H_<geno>')
    return lbl[2:] if isinstance(lbl, str) and lbl.startswith("H_") else lbl

def gained(parent_vec, child_vec):
    if parent_vec is not None: parent_vec = _geno(parent_vec)
    child_vec = _geno(child_vec)
    return [i for i in range(len(child_vec)) if child_vec[i] == "A" and (parent_vec is None or parent_vec[i] != "A")]

# === gap#2(2026-07-04):parent-choice 拓樸歧義枚舉 + 全區 identifiability ===
def _altset(g): return frozenset(i for i, ch in enumerate(g) if ch == "A")

def parent_choice_enum(r):
    """A_ambiguous 但無 multi-flip 邊(單-flip parent-choice 歧義,舊碼 silently return None 的區):
    對每觀測節點列『所有等大最大觀測真子集父』,笛卡兒積枚舉 → dedup 得所有 minimal 拓樸。k<=8 暴力。"""
    pops = r["populations"]
    genos = sorted([g for g in pops if "A" in g], key=lambda g: (g.count("A"), g))
    if not genos:
        return []
    k = len(genos[0]); obs = [_altset(g) for g in genos]; lab = {_altset(g): g for g in genos}
    def L(a): return "ROOT" if not a else lab.get(a, "H_" + "".join("A" if i in a else "R" for i in range(k)))
    per = []
    for X in obs:
        cands = [a for a in obs if a < X]
        par = [a for a in cands if len(a) == max(len(c) for c in cands)] if cands else [frozenset()]
        per.append((X, par))
    trees = [[]]
    for X, par in per:
        trees = [t + [(L(p), L(X))] for t in trees for p in par][:CAP]
    seen = set(); uniq = []
    for t in trees:
        key = frozenset(t)
        if key not in seen:
            seen.add(key); uniq.append({"edges": t, "virtual_nodes": [], "read_score": read_score(t, pops)})
    return uniq


def resolve_and_order(chrom, gidx, positions, vafmap):
    """gained 突變(vector index)→ block group 併合 + nested 定序約束 → 候選 group-序列(linear extensions)。"""
    par = {i: i for i in gidx}

    def find(x):
        root = x
        while par[root] != root:
            root = par[root]
        while par[x] != root:
            par[x], x = root, par[x]
        return root

    order, conflict, nocoread, tent = [], False, False, False
    for a, b in combinations(gidx, 2):
        pa, pb = positions[a], positions[b]
        lo, hi = min(pa, pb), max(pa, pb)
        pr = plut.get((chrom, lo, hi))
        if pr is None:
            nocoread = True
            continue
        c = pr["cells"]
        ra, ar, aa, cr = c["RA"], c["AR"], c["AA"], pr["coread"]
        thr = max(2, EPS * cr)
        big, small = max(ra, ar), min(ra, ar)
        if aa < 2:
            conflict = True
        elif big == 0:                 # co_linked → 併 group
            par[find(a)] = find(b)
        elif small >= thr:             # 4-gamete
            conflict = True
        else:                          # nested: 小邊<ε → 定序
            lo_idx = a if pa == lo else b
            hi_idx = b if lo_idx == a else a
            earlier = lo_idx if ar >= ra else hi_idx   # ar 大=lo 單獨出現=lo 先
            later = hi_idx if earlier == lo_idx else lo_idx
            order.append((earlier, later))
            if big < thr:
                tent = True
    grp = {}
    for i in gidx:
        grp.setdefault(find(i), set()).add(i)
    gid = {i: find(i) for i in gidx}
    gorder = {(gid[e], gid[l]) for e, l in order if gid[e] != gid[l]}
    keys = list(grp.keys())
    exts = []
    for perm in permutations(keys):
        if all(perm.index(e) < perm.index(l) for e, l in gorder):
            exts.append([grp[k] for k in perm])
            if len(exts) >= CAP:
                break
    if conflict:
        klass = "CONFLICT"
    elif nocoread:
        klass = "AMBIG_NOCOREAD"
    elif len(grp) == 1:
        klass = "BLOCK"
    elif order and len(exts) == 1:
        klass = "ORDERED"
    else:
        klass = "PARTIAL"
    # VAF 交叉核對(較早突變應較高 VAF)
    vaf_ok = None
    if order and vafmap:
        ck = [(vafmap.get(f"{chrom}:{positions[e]}"), vafmap.get(f"{chrom}:{positions[l]}")) for e, l in order]
        ck = [(x, y) for x, y in ck if x is not None and y is not None]
        vaf_ok = all(x >= y - 0.05 for x, y in ck) if ck else None
    return {"klass": klass, "group_orders": exts, "tentative": tent, "vaf_ok": vaf_ok,
            "n_order_constraints": len(order)}


def materialize_groups(edges, combo, glen):
    """combo=[(p,c,group_chain)]。co_linked group 一步加完(無中間虛擬)。"""
    repl = {(p, c): chain for (p, c, chain) in combo}
    new_edges, virtual = [], []
    for p, c in edges:
        if (p, c) in repl:
            chain = repl[(p, c)]
            acc = list("R" * glen if p == "ROOT" else p)
            prev = p
            for i, G in enumerate(chain):
                for idx in G:
                    acc[idx] = "A"
                vname = "".join(acc)
                if i < len(chain) - 1:
                    new_edges.append([prev, vname]); virtual.append(vname); prev = vname
                else:
                    new_edges.append([prev, c])
        else:
            new_edges.append([p, c])
    return new_edges, virtual


def materialize_perm(edges, ambig_choice, glen):
    """fallback(無 pairs/對映不齊):原逐突變排列。"""
    repl = {(p, c): order for (p, c, order) in ambig_choice}
    new_edges, virtual = [], []
    for p, c in edges:
        if (p, c) in repl:
            acc = list("R" * glen if p == "ROOT" else p)
            prev = p
            for pos in repl[(p, c)][:-1]:
                acc[pos] = "A"; vname = "".join(acc)
                new_edges.append([prev, vname]); virtual.append(vname); prev = vname
            new_edges.append([prev, c])
        else:
            new_edges.append([p, c])
    return new_edges, virtual


def read_score(edges, pops):
    """read-based parsimony;虛擬節點 reads=0 不貢獻。🔴 不用甲基。"""
    depth, ch = {}, {}
    for p, c in edges:
        ch.setdefault(p, []).append(c)
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


def _finalize(r, cands, ambig, true_count, resolution, edge_res, glen):
    scores = [c["read_score"] for c in cands]
    mx = max(scores) if scores else 0
    exps = [math.exp(s - mx) for s in scores]; tot = sum(exps) or 1
    for c, e in zip(cands, exps):
        c["softmax_prob"] = round(e / tot, 4)
    cands.sort(key=lambda c: -c["softmax_prob"])
    for i, c in enumerate(cands):
        c["rank"] = i + 1
    equ = len(cands) > 1 and len(set(scores)) <= 1
    note = {"BLOCK": "block 共event(gained 對 co_linked·一起獲得·非缺中間群)",
            "ORDERED": "順序可定(gained 對 nested·中間群被觀測·pairwise 定向)",
            "CONFLICT": "gained 對 4-gamete 衝突(可能 CN-multiplicity)",
            "PARTIAL": "部分定序(混 block+nested)+殘餘等機率",
            "AMBIG_NOCOREAD": "真等機率·gained 對無直接 coread(中間群未觀察)",
            "MIXED": "多邊混合解析"}.get(resolution, "")
    # gap#2: isomorphism dedup(canonical frozenset of labeled edges)→ n_minimal_trees + identifiability
    seen = set(); nmin = 0
    for c in cands:
        key = frozenset((p, cc) for p, cc in c["edges"])
        if key not in seen:
            seen.add(key); nmin += 1
    ident = ("conflict" if resolution == "CONFLICT"
             else "ambiguous" if (nmin > 1 or true_count > 1) else "determined")
    return {"region": r["region"], "base_determinacy": DET(r), "n_sSNV": r["n_sSNV"],
            "ambig_nodes": r.get("ambig_nodes", 0),
            "missing_intermediates": sum(len(g) - 1 for (_, _, g) in ambig),
            "n_candidates": len(cands), "true_count": true_count, "capped": true_count > CAP,
            "n_minimal_trees": nmin, "identifiability": ident,  # gap#2: 拓樸可辨識性(dedup 後)
            "equiprobable": equ, "resolution": resolution, "edge_resolution": edge_res,
            "candidate_set": cands, "honest_note": note,
            "truncated": bool(r.get("truncated")) or r["n_sSNV"] > glen}


def enumerate_region(r):
    pops = r["populations"]; edges = r.get("edges", [])
    if not edges or not pops:
        return None
    glen = len(next(iter(pops)))
    ambig = [(p, c, gained(None if p == "ROOT" else p, c)) for p, c in edges]
    ambig = [(p, c, g) for (p, c, g) in ambig if len(g) >= 2]
    if not ambig:
        return None
    positions = positions_for(r, glen) if _HAVE_PAIRS else None
    if positions is None:
        # fallback: 原排列枚舉(對映不齊/無 pairs)
        per = [(p, c, (list(permutations(g)) if len(g) <= 4 else [tuple(g)])) for (p, c, g) in ambig]
        tc = 1
        for (_, _, g) in ambig:
            tc *= math.factorial(len(g)) if len(g) <= 4 else 1
        combos = [[]]
        for (p, c, orders) in per:
            combos = [combo + [(p, c, o)] for combo in combos for o in orders][:CAP]
        cands = []
        for combo in combos:
            ne, vn = materialize_perm(edges, combo, glen)
            cands.append({"edges": ne, "virtual_nodes": vn, "read_score": read_score(ne, pops)})
        return _finalize(r, cands, ambig, tc, "UNMAPPED_PERM", [], glen)
    vafmap = r.get("pos_vaf", {})
    chrom = r["chrom"]
    edge_infos, edge_res, klasses = [], [], []
    for (p, c, g) in ambig:
        info = resolve_and_order(chrom, g, positions, vafmap)
        edge_infos.append((p, c, info))
        klasses.append(info["klass"])
        edge_res.append({"parent": p, "child": c, "gained_pos": [positions[i] for i in g],
                         "klass": info["klass"], "tentative": info["tentative"], "vaf_ok": info["vaf_ok"]})
    combos = [[]]
    for (p, c, info) in edge_infos:
        combos = [combo + [(p, c, go)] for combo in combos for go in info["group_orders"]][:CAP]
    cands = []
    for combo in combos:
        ne, vn = materialize_groups(edges, combo, glen)
        cands.append({"edges": ne, "virtual_nodes": vn, "read_score": read_score(ne, pops)})
    tc = 1
    for (_, _, info) in edge_infos:
        tc *= max(1, len(info["group_orders"]))
    resolution = klasses[0] if len(set(klasses)) == 1 else "MIXED"
    return _finalize(r, cands, ambig, tc, resolution, edge_res, glen)


b1, b2_nonenum = [], 0
rescnt = defaultdict(int)
ident_table = defaultdict(int)  # gap#2: 全 with-vector 區 identifiability 分佈(不 bloat candidate_set)
for r in det:
    d = DET(r)
    if "ambiguous" in d:
        res = enumerate_region(r)
        if res is None:
            # gap#2 修:舊碼對「無 multi-flip 邊」的 parent-choice/單-flip 歧義 silently return None(25 區憑空消失)。
            #   改枚舉等大觀測父選擇,顯式標 identifiability,不再丟。
            pc = parent_choice_enum(r); nmin = len(pc)
            res = {"region": r["region"], "base_determinacy": d, "n_sSNV": r["n_sSNV"],
                   "ambig_nodes": r.get("ambig_nodes", 0), "missing_intermediates": 0,
                   "n_candidates": len(pc), "true_count": nmin, "capped": nmin > CAP,
                   "n_minimal_trees": nmin, "identifiability": "ambiguous" if nmin > 1 else "determined",
                   "equiprobable": nmin > 1, "resolution": "PARENT_CHOICE", "edge_resolution": [],
                   "candidate_set": pc, "honest_note": "parent-choice 拓樸歧義(節點有多個等大觀測父)",
                   "truncated": bool(r.get("truncated"))}
        b1.append(res); rescnt[res["resolution"]] += 1
        ident_table[res["identifiability"]] += 1
    elif "E_subcube_recovered" in d:
        ident_table["subcube_recovered(gap#1;partial建樹)"] += 1          # gap#1: partial-read 救回(無 full-cov)
    elif "A_determined" in d:
        ident_table["determined"] += 1                                  # gap#2: 蓋章(n_minimal_trees=1)
    elif "recurrence" in d:
        ident_table["recurrence(m通道;非乾淨可辨識)"] += 1
    elif d == "incompatible":
        ident_table["conflict(真artifact/cycle)"] += 1
    elif d == "C_underdetermined":
        b2_nonenum += 1; ident_table["nonenumerable(C_underdetermined)"] += 1
    elif d == "B_pairwise_structure":
        ident_table["nonenumerable(B_pairwise)"] += 1
    else:
        ident_table["other"] += 1

equ_n = sum(1 for x in b1 if x["equiprobable"])
b1_ambiguous = sum(1 for x in b1 if x["identifiability"] == "ambiguous")
b1_determined = sum(1 for x in b1 if x["identifiability"] == "determined")
edge_klass = defaultdict(int)
for x in b1:
    for e in x.get("edge_resolution", []):
        edge_klass[e["klass"]] += 1
out = {
    "schema_version": "20260704",
    "summary": {
        "B1_enumerated": len(b1),
        "B1_ambiguous": b1_ambiguous,          # gap#2: n_minimal_trees>1(真拓樸歧義)
        "B1_determined": b1_determined,         # gap#2: enumerate 後收斂唯一
        "B1_equiprobable": equ_n,
        "B1_resolved_nonequiprobable": len(b1) - equ_n,
        "identifiability_table": dict(ident_table),  # gap#2: 全 with-vector 區 determined/ambiguous/recurrence/conflict/nonenumerable 完整分佈
        "resolution_dist": dict(rescnt),
        "edge_klass_dist": dict(edge_klass),
        "candidates_total": sum(x["n_candidates"] for x in b1),
        "B2_underdetermined_nonenumerable": b2_nonenum,
        "note": "2026-07-04 gained-pair pairwise 定序:co_linked→block group / nested→定序 / 4-gamete→conflict / "
                "無coread→真等機率。用觀察約束 prune 排列。🔴絕不用甲基。B2/跨-PS non-identifiability 仍成立。"},
    "candidate_trees": b1,
    "B2_underdetermined": {"n": b2_nonenum, "status": "真欠定·單一ALT群或缺對·非可列舉"},
}
with open(os.path.join(DATA, "candidate_trees.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False)
print("ENUMERATE DONE →", os.path.join(DATA, "candidate_trees.json"))
print(json.dumps(out["summary"], ensure_ascii=False, indent=1))
