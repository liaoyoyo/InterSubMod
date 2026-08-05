#!/usr/bin/env python3
"""
tree_enumeration_solver.py — 覆蓋全節點 + 枚舉全最小 Model-A group-Steiner solver (2026-07-06)

Formal spec (SoT): InterSubMod/docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md §5-§7.5
複雜度定位 : InterSubMod/docs/method_comparison/20260705_ism_computational_complexity_positioning_01.md §1.3/§2.3

核心（使用者 2026-07-06 定案）:
  §7.5.1 覆蓋公理 : 每觀測(full=singleton group / partial=subcube group)必被最終樹覆蓋;完整觀測其基因型節點必 ∈ V;禁為求乾淨樹丟節點。
  §7.5.3 枚舉紅線 : 枚舉「所有」最小(min hidden-node)相容 unit-flip arborescence,非求單一 min-cost;不用任何通道對樹集排序(甲基事後)。
  §4/§7.5.2 Model A: recurrence 允許 → post-hoc 偵測翻兩次的 bit → 送獨立 m-通道拆 artifact vs 真 recurrence(不在結構目標懲罰)。
  §7.5.2 capped k≤8 : |X|=2^k≤256 → 暴力窮舉 tractable ≠ 解通用 NP-hard directed/group Steiner。
  §7.5.2 IDP      : partial(=「?」entry)補全 = Incomplete Directed Perfect Phylogeny (Pe'er,Pupko,Shamir,Sharan 2004),多完成 non-identifiable → 枚舉。

節點內部表示 = frozenset(bit index)（altset,R=0/A=1；root=frozenset()）。輸出時 label 成基因型字串。

V1-V7 驗證器(§step2)：獨立重算,不信任 enumerate 自報。
"""
import math
from itertools import combinations, product
from collections import Counter

# ---- 觀測解析 ----------------------------------------------------------------

def parse_full(g):
    """完整基因型字串(R/A) → altset frozenset。"""
    return frozenset(i for i, c in enumerate(g) if c == "A")

def parse_partial(s):
    """部分觀測字串(R/A/X) → (ones, zeros)。X=未觀測(「?」)。"""
    ones = frozenset(i for i, c in enumerate(s) if c == "A")
    zeros = frozenset(i for i, c in enumerate(s) if c == "R")
    return ones, zeros

def subcube_members(ones, zeros, universe):
    """partial 的相容子立方體(§3):Y⊇ones、Y∩zeros=∅、Y⊆universe。回 frozenset 節點列表。"""
    unknown = [i for i in universe if i not in ones and i not in zeros]
    mem = []
    for r in range(len(unknown) + 1):
        for extra in combinations(unknown, r):
            mem.append(frozenset(ones) | frozenset(extra))
    return mem

# ---- label ------------------------------------------------------------------

def label(x, full_set, k):
    """frozenset 節點 → 顯示 label。root→'ROOT';完整觀測→基因型字串;其他→'H_...'(隱藏/partial-supported)。"""
    if not x:
        return "ROOT"
    s = "".join("A" if i in x else "R" for i in range(k))
    return s if x in full_set else "H_" + s

# ---- 樹合法性基元 -----------------------------------------------------------

def _unit_preds(x):
    """x 的所有 unit-flip 前驅(移除恰一個 bit)。"""
    return [x - {j} for j in x]

def _is_closed(N):
    """unit-flip-closed:每非根節點至少一 unit-pred ∈ N → 可建 root-connected arborescence。"""
    for x in N:
        if x and not any((x - {j}) in N for j in x):
            return False
    return True

def _covers(N, full, partial_subcubes):
    """覆蓋公理(§7.5.1):full ⊆ N 且每 partial subcube 與 N 相交。"""
    if not all(f in N for f in full):
        return False
    for sub in partial_subcubes:
        if not (N & sub):
            return False
    return True

def _is_compatible(N):
    """Rooted three-gamete 相容性：每對字元若同時出現 (1,1)(1,0)(0,1) 即衝突。
    N 恆含 all-zero ROOT，故這與在 N∪{ROOT} 上做 four-gamete test 完全等價。
    定理:N 相容 ⟺ 存在無-recurrence 生成樹(唯一 perfect phylogeny);不相容 ⟺ 每棵生成樹皆 recurrence。
    🔴 不可用「節點是否有≥2 pred」代替(反例:全≤1pred 仍可 rooted-three-gamete 衝突→唯一樹含 recurrence)。"""
    chars = sorted(set().union(*N)) if N else []
    for i, j in combinations(chars, 2):
        pats = set()
        for x in N:
            pats.add((i in x, j in x))
        if (True, True) in pats and (True, False) in pats and (False, True) in pats:
            return False
    return True

def recurrence_bits(edges):
    """post-hoc(§4):某 bit 在 ≥2 條邊被翻(homoplasy)= recurrence。回翻兩次以上的 bit 集合。"""
    flip = Counter()
    for (p, c) in edges:
        d = c - p                       # unit-flip → 恰一個 bit
        flip[next(iter(d))] += 1
    return frozenset(b for b, n in flip.items() if n >= 2)

# ---- 核心:枚舉所有最小相容樹 -----------------------------------------------

def enumerate_min_trees(full_genos, partial_strs, k, extra_cap=4, per_level_budget=150_000, tree_cap=32):
    """
    full_genos : {基因型字串(R/A): count} 或 [基因型字串]
    partial_strs: [部分觀測字串(R/A/X)]
    回 dict:
      n_hidden      : 最小 extra(hidden/partial-supported)節點數
      trees         : [{edges:[(labP,labC)], recurrence:[bit...], n_hidden:int}]  全部同分最小樹
      n_trees       : exact analytical candidate count
      trees_complete: whether trees contains every candidate (tree_cap <= 0 means unlimited)
      capped        : bool(超 budget → fallback 單樹)
      cap_reason    : str|None
      universe      : sorted bit indices
      full_labels   : [完整觀測 label]
    """
    if isinstance(full_genos, dict):
        full_counts = dict(full_genos)
        full_genos = list(full_genos.keys())
    else:
        full_counts = {g: 1 for g in full_genos}
    full = set(parse_full(g) for g in full_genos)
    partials = [parse_partial(s) for s in partial_strs]
    # universe = 任何觀測(full/partial)曾為 A 的位置;X-only 且從不 A 的位置視為恆 0(不入 universe)
    universe = sorted(set().union(*[set(x) for x in full], *[o for (o, z) in partials]) if (full or partials) else set())
    partial_subcubes = [set(subcube_members(o, z, universe)) for (o, z) in partials]

    ROOT = frozenset()
    base = {ROOT} | full                # 覆蓋公理:root + 完整觀測恆在
    # 候選 extra pool = 任何觀測/partial-member 的子集(潛在祖先),排除 root 與 full
    maximal = set(full)
    for sub in partial_subcubes:
        maximal |= sub
    pool = set()
    for m in maximal:
        ml = sorted(m)
        for r in range(len(ml) + 1):
            for c in combinations(ml, r):
                pool.add(frozenset(c))
    extra_pool = sorted(pool - base, key=lambda s: (len(s), sorted(s)))

    capped = False
    cap_reason = None
    feasible_N = None
    n_hidden = None
    # 逐層(extra 節點數 e 遞增)搜尋最小相容節點集。每層先算 C(pool,e),超 budget → fast-cap(不逐一燒),
    # 保證小 e_min(0-2,可枚舉的有意義區)秒回,只有 dense 大 pool 需要的高 e 層 cap(NP-hard territory,§7.5.2)。
    for e in range(0, min(extra_cap, len(extra_pool)) + 1):
        if math.comb(len(extra_pool), e) > per_level_budget:
            capped = True
            cap_reason = f"level e={e}: C({len(extra_pool)},{e})={math.comb(len(extra_pool),e)} > budget {per_level_budget}(太密)"
            break
        found = []
        for extra in combinations(extra_pool, e):
            N = base | set(extra)
            if _covers(N, full, partial_subcubes) and _is_closed(N):
                found.append(frozenset(N))
        if found:
            feasible_N = found
            n_hidden = e
            break

    if feasible_N is None:
        # 未在預算/上限內找到 → greedy fallback 單樹(仍滿足覆蓋公理),flag capped
        capped = True
        cap_reason = cap_reason or f"no feasible N within extra_cap={extra_cap}(pool={len(extra_pool)})"
        N = _greedy_closure(base, full, partial_subcubes)
        feasible_N = [frozenset(N)]
        n_hidden = len(N) - len(full) - 1

    # 分析式計數(2026-07-06,效能:取代生成所有樹)。定理:單一 feasible N 內
    #   節點 x 的 unit-pred 數 npred(x) → 該 N 的樹數 = Π npred(x);且
    #   節點有 npred≥2 ⟺ 該 N rooted-three-gamete 不相容
    #   (等價於明示 ROOT 後的 four-gamete 不相容) ⟺ 其所有樹皆 recurrence;
    #   npred 全=1 → 相容 → 唯一無-rec 樹。
    #   跨 feasible N 求和(不同節點集→不同樹,不重疊)。has_recfree=∃相容N(有乾淨 perfect-phylogeny 解)。
    n_trees, n_recfree_N, n_recurrent_N = 0, 0, 0
    for N in feasible_N:
        prod = 1
        for x in N:
            if x:
                prod *= sum(1 for j in x if (x - {j}) in N)  # unit-pred 數 ≥1(feasible=closed 保證)
        n_trees += prod
        # recurrence = N rooted-three-gamete 不相容(=root-augmented four-gamete 不相容;
        # ⟺ 所有生成樹皆含 recurrence);相容 → 唯一無-rec 樹(prod=1)
        if _is_compatible(N):
            n_recfree_N += 1
        else:
            n_recurrent_N += 1
    has_recfree = n_recfree_N > 0        # ∃ 無-recurrence 的最小解
    has_recurrence = n_recurrent_N > 0   # ∃ 需 recurrence 的最小解
    # Generate candidates for downstream analysis. A non-positive cap means all
    # candidates; a positive cap is display-only and must never be treated as a
    # complete analytical set.
    store_limit = None if tree_cap is None or tree_cap <= 0 else int(tree_cap)
    trees = []
    for N in feasible_N:
        for edges in _parent_choice_trees(N):
            lab_edges = [(label(p, full, k), label(c, full, k)) for (p, c) in edges]
            trees.append({"edges": lab_edges, "recurrence": sorted(recurrence_bits(edges)),
                          "n_hidden": len(set(N) - ({ROOT} | full))})
            if store_limit is not None and len(trees) >= store_limit:
                break
        if store_limit is not None and len(trees) >= store_limit:
            break

    return {"n_hidden": n_hidden, "trees": trees, "n_trees": n_trees, "n_trees_stored": len(trees),
            "trees_complete": (len(trees) == n_trees), "tree_cap": tree_cap,
            "has_recfree": has_recfree, "has_recurrence": has_recurrence,
            "n_feasible_N": len(feasible_N), "n_recfree_N": n_recfree_N, "n_recurrent_N": n_recurrent_N,
            "capped": capped, "cap_reason": cap_reason,
            "universe": universe, "full_labels": sorted(label(f, full, k) for f in full),
            "_full": full, "_partial_subcubes": partial_subcubes, "_feasible_N": feasible_N}

def _parent_choice_trees(N):
    """給定節點集 N,枚舉每個非根節點獨立選 unit-pred∈N 的所有 arborescence。
    邊皆 popcount 遞增 → 天然無環 + root-connected(每節點有父且 root 為終祖)。"""
    nodes = [x for x in N if x]
    choice_lists = []
    for x in nodes:
        preds = [x - {j} for j in x if (x - {j}) in N]
        choice_lists.append([(p, x) for p in preds])
    for combo in product(*choice_lists):
        yield list(combo)

def _greedy_closure(base, full, partial_subcubes):
    """fallback:覆蓋所有 group 的一個 closed 節點集(每 partial 取最小成本 member,補 unit-pred 鏈)。"""
    N = set(base)
    # 每個未覆蓋 partial 取 bit 最少的 member
    for sub in partial_subcubes:
        if not (N & sub):
            rep = min(sub, key=lambda s: (len(s), sorted(s)))
            N.add(rep)
    # 補 unit-pred 鏈直到 closed(每次沿 min-index 前驅補)
    changed = True
    while changed:
        changed = False
        for x in list(N):
            if x and not any((x - {j}) in N for j in x):
                N.add(x - {min(x)})
                changed = True
    return N

# ---- 可辨識性分類(§6) ------------------------------------------------------

def classify(result):
    """determined / ambiguous / recurrence_required / underdetermined。用分析式 has_recfree/has_recurrence。
    🔴 capped(枚舉未完,§7.5.2 NP-hard-dense)→ 絕不宣稱 determined/枚舉全(V7 會擋);誠實標 (capped;枚舉未完)。
    recurrence_required = 「所有最小解都需 recurrence」(無乾淨 perfect-phylogeny 解);若 ∃ 乾淨解則歸 determined/ambiguous
    (recurrence 只是替代選項非必需)。"""
    n_trees = result.get("n_trees", 0)
    if n_trees == 0:
        return "underdetermined"
    if result.get("capped"):
        return ("recurrence_required(capped;枚舉未完)" if result.get("has_recurrence") and not result.get("has_recfree")
                else "underdetermined(capped;枚舉未完;太密)")
    if result.get("has_recfree"):          # ∃ 乾淨 perfect-phylogeny 解 → 由樹數判 determined/ambiguous
        return "determined" if n_trees == 1 else "ambiguous"
    if result.get("has_recurrence"):       # 所有最小解都需 recurrence
        return "recurrence_required"
    return "underdetermined"

# ================= V1-V7 驗證器(獨立重算,不信任 enumerate 自報) =================

def _edges_to_fs(lab_edges, k):
    """label 邊 → frozenset 邊(還原 bit set),供獨立幾何檢查。"""
    def unlab(s):
        if s == "ROOT":
            return frozenset()
        s2 = s[2:] if s.startswith("H_") else s
        return frozenset(i for i, c in enumerate(s2) if c == "A")
    return [(unlab(p), unlab(c)) for (p, c) in lab_edges]

def V1_valid_arborescence(tree, k):
    """每邊 unit-flip(c=p+1bit)、每非根節點恰一父、全 root-可達、無環。"""
    fe = _edges_to_fs(tree["edges"], k)
    for (p, c) in fe:
        if not (p < c and len(c - p) == 1):
            return False, f"non-unit-flip edge {p}->{c}"
    parent = {}
    for (p, c) in fe:
        if c in parent:
            return False, f"node {c} has >1 parent"
        parent[c] = p
    # root-可達 + 無環(沿 parent 上溯必達 ∅)
    for c in parent:
        seen, cur = set(), c
        while cur:
            if cur in seen:
                return False, f"cycle at {cur}"
            seen.add(cur)
            if cur not in parent:
                return False, f"{cur} not root-connected"
            cur = parent[cur]
    return True, "ok"

def V2_coverage(tree, full, partial_subcubes, k):
    """覆蓋公理(§7.5.1):每完整觀測 ∈ 節點;每 partial subcube ∩ 節點 ≠ ∅。"""
    fe = _edges_to_fs(tree["edges"], k)
    nodes = {frozenset()} | {c for (_, c) in fe} | {p for (p, _) in fe}
    for f in full:
        if f not in nodes:
            return False, f"full obs {f} not covered"
    for sub in partial_subcubes:
        if not (nodes & sub):
            return False, f"partial group {min(sub, key=len) if sub else '?'} not covered"
    return True, "ok"

def V3_consistency(tree, full, k):
    """觀測完整基因型 ⊆ 樹節點(不引入與觀測矛盾的節點:所有節點皆 universe 內合法 altset,恆真;檢 full 全在)。"""
    fe = _edges_to_fs(tree["edges"], k)
    nodes = {frozenset()} | {c for (_, c) in fe} | {p for (p, _) in fe}
    return (full <= nodes), ("ok" if full <= nodes else "observed not subset of nodes")

def V6_model_correct(result):
    """分析旗標一致性:recurrence_required ⟺ 有 recurrent-N 且無 rec-free-N;determined/ambiguous ⟹ 有 rec-free-N。
    另檢:stored 範例樹的 recurrence 標註與 has_recurrence 不矛盾。"""
    cls = classify(result)
    hr = result.get("has_recurrence"); hf = result.get("has_recfree")
    if cls == "recurrence_required":
        if not hr:
            return False, "recurrence_required but has_recurrence=False"
        if hf:
            return False, "recurrence_required but a rec-free minimal solution exists"
    if cls in ("determined", "ambiguous") and not hf:
        return False, f"{cls} but no rec-free minimal solution"
    # stored 範例樹若全 rec,則 has_recurrence 必 True
    if result["trees"] and all(t["recurrence"] for t in result["trees"]) and not hr:
        return False, "all stored trees recurrent but has_recurrence=False"
    return True, "ok"

def V7_no_overclaim(result):
    """determined ⟹ 恰 1 樹 且 無 recurrence 且 有 rec-free 且 not capped。"""
    cls = classify(result)
    if cls == "determined":
        if result.get("n_trees") != 1:
            return False, f"determined but n_trees={result.get('n_trees')}"
        if result.get("has_recurrence"):
            return False, "determined but has_recurrence"
        if not result.get("has_recfree"):
            return False, "determined but no rec-free solution"
        if result.get("capped"):
            return False, "determined but enumeration capped(可能漏樹)"
    return True, "ok"

def V4_minimality_independent(result, full_genos, partial_strs, k):
    """獨立重算最小 hidden:確認 n_hidden-1 層全 infeasible(不信任 enumerate 的早停)。"""
    if result["capped"]:
        return None, "skipped(capped)"
    e_min = result["n_hidden"]
    if e_min == 0:
        return True, "min hidden=0(trivially minimal)"
    # 重建 pool 並窮舉 e_min-1 層,須全 infeasible
    full = set(parse_full(g) for g in (full_genos.keys() if isinstance(full_genos, dict) else full_genos))
    partials = [parse_partial(s) for s in partial_strs]
    universe = sorted(set().union(*[set(x) for x in full], *[o for (o, z) in partials]) if (full or partials) else set())
    subs = [set(subcube_members(o, z, universe)) for (o, z) in partials]
    base = {frozenset()} | full
    maximal = set(full)
    for sub in subs:
        maximal |= sub
    pool = set()
    for m in maximal:
        ml = sorted(m)
        for r in range(len(ml) + 1):
            for c in combinations(ml, r):
                pool.add(frozenset(c))
    extra_pool = sorted(pool - base)
    for extra in combinations(extra_pool, e_min - 1):
        N = base | set(extra)
        if _covers(N, full, subs) and _is_closed(N):
            return False, f"found feasible N with {e_min-1} hidden(non-minimal!)"
    return True, f"confirmed min hidden={e_min}"

def V5_enumerate_all_independent(result, full_genos, partial_strs, k):
    """獨立重算最小層全樹數,與 enumerate 回報比對(no missing / no extra)。capped 則 skip。"""
    if result["capped"]:
        return None, "skipped(capped)"
    full = set(parse_full(g) for g in (full_genos.keys() if isinstance(full_genos, dict) else full_genos))
    partials = [parse_partial(s) for s in partial_strs]
    universe = sorted(set().union(*[set(x) for x in full], *[o for (o, z) in partials]) if (full or partials) else set())
    subs = [set(subcube_members(o, z, universe)) for (o, z) in partials]
    base = {frozenset()} | full
    maximal = set(full)
    for sub in subs:
        maximal |= sub
    pool = set()
    for m in maximal:
        ml = sorted(m)
        for r in range(len(ml) + 1):
            for c in combinations(ml, r):
                pool.add(frozenset(c))
    extra_pool = sorted(pool - base)
    e_min = result["n_hidden"]
    # 獨立重算 feasible_N + 分析式 n_trees(不生成所有樹),比對 main 回報
    n_trees_ind = 0
    for extra in combinations(extra_pool, e_min):
        N = base | set(extra)
        if _covers(N, full, subs) and _is_closed(N):
            prod = 1
            for x in N:
                if x:
                    prod *= sum(1 for j in x if (x - {j}) in N)
            n_trees_ind += prod
    return (n_trees_ind == result["n_trees"],
            f"independent={n_trees_ind} vs reported={result['n_trees']}")

def verify_all(result, full_genos, partial_strs, k, light=False):
    """跑 V1-V7,回 {V1..V7: (pass, msg)} + overall。
    light=True:跳過 V4/V5(獨立重算最小性/完整性,昂貴)→ 只跑 V1/V2/V3/V6/V7(便宜 per-tree)。
    全量/多樣本時用 light=True 全跑 + 抽樣做完整 V4/V5(避免 18k units 各重算兩次的成本)。"""
    full = set(parse_full(g) for g in (full_genos.keys() if isinstance(full_genos, dict) else full_genos))
    partials = [parse_partial(s) for s in partial_strs]
    universe = result["universe"]
    subs = [set(subcube_members(o, z, universe)) for (o, z) in partials]
    out = {}
    # per-tree V1/V2/V3
    for name, fn in [("V1", lambda t: V1_valid_arborescence(t, k)),
                     ("V2", lambda t: V2_coverage(t, full, subs, k)),
                     ("V3", lambda t: V3_consistency(t, full, k))]:
        ok, msg = True, "ok"
        for t in result["trees"]:
            o, m = fn(t)
            if not o:
                ok, msg = False, f"tree fail: {m}"
                break
        out[name] = (ok, msg)
    if light:
        out["V4"] = (None, "skipped(light)")
        out["V5"] = (None, "skipped(light)")
    else:
        out["V4"] = V4_minimality_independent(result, full_genos, partial_strs, k)
        out["V5"] = V5_enumerate_all_independent(result, full_genos, partial_strs, k)
    out["V6"] = V6_model_correct(result)
    out["V7"] = V7_no_overclaim(result)
    out["overall"] = all(v[0] for v in out.values() if v[0] is not None)
    return out

# ================= golden 手算對照(§step3) =================

def _run_golden():
    """手算真值對照。每案 (名, full, partials, k, 期望分類, 期望樹數, 期望 hidden)。"""
    cases = [
        # 1. 兩獨立 sibling → 唯一樹(各只 root 為父)
        ("two_siblings {10,01}", ["AR", "RA"], [], 2, "determined", 1, 0),
        # 2. 單一雙突變,缺中間群 → 順序未定 2 樹
        ("single_double {11}", ["AA"], [], 2, "ambiguous", 2, 1),
        # 3. rooted three-gamete 衝突 {10,01,11};ROOT 補 00 後為 four-gamete。
        #    legacy case 名保留，避免既有 log/snapshot 名稱漂移。
        ("four_gamete {10,01,11}", ["AR", "RA", "AA"], [], 2, "recurrence_required", 2, 0),
        # 4. k=3 同一 rooted three-gamete 條件(bit0,bit1){100,010,110};legacy 名保留
        ("four_gamete_k3 {100,010,110}", ["ARR", "RAR", "AAR"], [], 3, "recurrence_required", 2, 0),
        # 5. 乾淨 nested {100,110} → 唯一 linear
        ("nested {100,110}", ["ARR", "AAR"], [], 3, "determined", 1, 0),
        # 6. 三獨立 sibling {100,010,001} → 唯一星狀
        ("three_siblings", ["ARR", "RAR", "RRA"], [], 3, "determined", 1, 0),
        # 7. partial-only IDP 多完成:AX(bit0=1)+XA(bit1=1) → 3 個 e=2 完成 {10,11}/{01,11}/{10,01}
        ("partial_IDP {AX,XA}", [], ["AX", "XA"], 2, None, None, 2),
        # 8. full+partial:full {100},partial 1X? (bit0=1) 已被 100 覆蓋 → 不新增
        ("full_covers_partial", ["ARR"], ["AR"], 2, "determined", 1, 0),
    ]
    print("=" * 92)
    print("GOLDEN 手算對照 — tree_enumeration_solver")
    print("=" * 92)
    allpass = True
    for name, full, part, k, exp_cls, exp_nt, exp_h in cases:
        res = enumerate_min_trees(full, part, k)
        cls = classify(res)
        ver = verify_all(res, full, part, k)
        ok_cls = (exp_cls is None) or (cls == exp_cls)
        ok_nt = (exp_nt is None) or (res["n_trees"] == exp_nt)
        ok_h = (exp_h is None) or (res["n_hidden"] == exp_h)
        ok_v = ver["overall"]
        good = ok_cls and ok_nt and ok_h and ok_v
        allpass = allpass and good
        print(f"\n[{'PASS' if good else 'FAIL'}] {name}")
        print(f"   class={cls}(exp {exp_cls})  n_trees={res['n_trees']}(exp {exp_nt})  "
              f"n_hidden={res['n_hidden']}(exp {exp_h})  capped={res['capped']}")
        vflags = " ".join(f"{k2}={'✓' if v[0] else ('–' if v[0] is None else '✗')}" for k2, v in ver.items() if k2 != "overall")
        print(f"   verify: {vflags}  overall={'✓' if ok_v else '✗'}")
        for t in res["trees"][:4]:
            print(f"     tree: {t['edges']}  rec={t['recurrence']}")
        if len(res["trees"]) > 4:
            print(f"     ... (+{len(res['trees'])-4} more)")
        for k2, v in ver.items():
            if k2 != "overall" and v[0] is False:
                print(f"   ⚠ {k2}: {v[1]}")
    print("\n" + "=" * 92)
    print(f"GOLDEN 總結: {'ALL PASS ✓' if allpass else 'SOME FAIL ✗'}")
    print("=" * 92)
    return allpass

if __name__ == "__main__":
    import sys
    ok = _run_golden()
    sys.exit(0 if ok else 1)
