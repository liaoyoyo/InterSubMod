#!/usr/bin/env python3
"""
complexity_bench.py — 逐步驟計算量實測 (tree_enumeration_solver)

目的:對 solver 的 7 個階段分別計時 + 記錄組合量,產生可引用的真值表。
方法:忠實鏡像 enumerate_min_trees 內部階段(逐行對照 :118-201),加 timer/counter。
自我驗證:鏡像版的 e_min / n_trees / capped 必須與官方 enumerate_min_trees 完全一致,
         否則該 case 標 MIRROR_MISMATCH(不可採信)。
輸出:complexity_bench.json + complexity_bench.tsv
"""
import json
import math
import random
import sys
import time
from itertools import combinations, product

sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/"
                   "20260627_subclone_4axis_teaching/scripts")
import tree_enumeration_solver as S  # noqa: E402

T = time.perf_counter


def instrumented(full_genos, partial_strs, k, extra_cap=4, budget=150_000,
                 materialize_limit=None, materialize_guard_s=20.0):
    """鏡像 enumerate_min_trees,逐階段計時。materialize_limit=None → 全生成(生產預設)。"""
    m = {}

    # ---- Step 1: 解析觀測 (:118-119) ----
    t0 = T()
    if isinstance(full_genos, dict):
        full_list = list(full_genos.keys())
    else:
        full_list = list(full_genos)
    full = set(S.parse_full(g) for g in full_list)
    partials = [S.parse_partial(s) for s in partial_strs]
    m["t1_parse"] = T() - t0

    # ---- Step 2: universe (:121) ----
    t0 = T()
    universe = sorted(set().union(*[set(x) for x in full], *[o for (o, z) in partials])
                      if (full or partials) else set())
    m["t2_universe"] = T() - t0

    # ---- Step 3: subcube groups (:122) ----
    t0 = T()
    partial_subcubes = [set(S.subcube_members(o, z, universe)) for (o, z) in partials]
    m["t3_subcubes"] = T() - t0
    m["subcube_sizes"] = [len(s) for s in partial_subcubes]
    m["subcube_total_members"] = sum(len(s) for s in partial_subcubes)

    ROOT = frozenset()
    base = {ROOT} | full

    # ---- Step 4: candidate Steiner pool (:127-136) ----
    t0 = T()
    maximal = set(full)
    for sub in partial_subcubes:
        maximal |= sub
    pool = set()
    subset_ops = 0
    for mm in maximal:
        ml = sorted(mm)
        for r in range(len(ml) + 1):
            for c in combinations(ml, r):
                pool.add(frozenset(c))
                subset_ops += 1
    extra_pool = sorted(pool - base, key=lambda s: (len(s), sorted(s)))
    m["t4_pool"] = T() - t0
    m["n_maximal"] = len(maximal)
    m["pool_subset_ops"] = subset_ops          # 冗餘度指標
    m["n_pool"] = len(pool)
    m["P_extra_pool"] = len(extra_pool)

    # ---- Step 5: 逐層最小解搜尋 (:144-157) — 主導階段 ----
    t0 = T()
    capped = False
    cap_reason = None
    feasible_N = None
    n_hidden = None
    level_stats = []
    cand_tested = 0
    for e in range(0, min(extra_cap, len(extra_pool)) + 1):
        c_pe = math.comb(len(extra_pool), e)
        if c_pe > budget:
            capped = True
            cap_reason = f"level e={e}: C({len(extra_pool)},{e})={c_pe} > budget {budget}"
            level_stats.append({"e": e, "C_P_e": c_pe, "tested": 0, "found": 0, "capped": True})
            break
        te = T()
        found = []
        for extra in combinations(extra_pool, e):
            N = base | set(extra)
            cand_tested += 1
            if S._covers(N, full, partial_subcubes) and S._is_closed(N):
                found.append(frozenset(N))
        level_stats.append({"e": e, "C_P_e": c_pe, "tested": c_pe,
                            "found": len(found), "capped": False,
                            "secs": T() - te})
        if found:
            feasible_N = found
            n_hidden = e
            break
    if feasible_N is None:
        capped = True
        cap_reason = cap_reason or f"no feasible N within extra_cap={extra_cap}"
        N = S._greedy_closure(base, full, partial_subcubes)
        feasible_N = [frozenset(N)]
        n_hidden = len(N) - len(full) - 1
    m["t5_level_search"] = T() - t0
    m["level_stats"] = level_stats
    m["candidates_tested"] = cand_tested
    m["e_min"] = n_hidden
    m["capped"] = capped
    m["cap_reason"] = cap_reason
    m["n_feasible_N"] = len(feasible_N)
    m["mean_N_size"] = (sum(len(N) for N in feasible_N) / len(feasible_N)) if feasible_N else 0

    # ---- Step 6: 分析式計數 + 相容性 (:173-187) ----
    t0 = T()
    n_trees, n_recfree, n_rec = 0, 0, 0
    for N in feasible_N:
        prod = 1
        for x in N:
            if x:
                prod *= sum(1 for j in x if (x - {j}) in N)
        n_trees += prod
        if S._is_compatible(N):
            n_recfree += 1
        else:
            n_rec += 1
    m["t6_analytic_count"] = T() - t0
    m["n_trees"] = n_trees
    m["n_recfree_N"] = n_recfree
    m["n_recurrent_N"] = n_rec

    # ---- Step 7: 實體化樹 (:191-201) ----
    t0 = T()
    trees = 0
    truncated = False
    for N in feasible_N:
        for edges in S._parent_choice_trees(N):
            [(S.label(p, full, k), S.label(c, full, k)) for (p, c) in edges]
            S.recurrence_bits(edges)
            trees += 1
            if materialize_limit is not None and trees >= materialize_limit:
                truncated = True
                break
            if (trees & 1023) == 0 and (T() - t0) > materialize_guard_s:
                truncated = True
                break
        if truncated:
            break
    m["t7_materialize"] = T() - t0
    m["trees_materialized"] = trees
    m["materialize_truncated"] = truncated

    m["t_total_no_verify"] = sum(m[key] for key in
                                 ["t1_parse", "t2_universe", "t3_subcubes", "t4_pool",
                                  "t5_level_search", "t6_analytic_count", "t7_materialize"])
    return m


def bench_case(name, full_genos, partial_strs, k, materialize_limit=None,
               run_v45=True, v45_guard_s=30.0):
    rec = {"case": name, "k": k,
           "n_full": len(full_genos), "n_partial": len(partial_strs)}
    m = instrumented(full_genos, partial_strs, k, materialize_limit=materialize_limit)
    rec.update(m)

    # 官方 solver 交叉核對(鏡像忠實度自我驗證)
    t0 = T()
    off = S.enumerate_min_trees(full_genos, partial_strs, k, tree_cap=0)
    rec["t_official_total"] = T() - t0
    rec["official_n_trees"] = off["n_trees"]
    rec["official_e_min"] = off["n_hidden"]
    rec["official_capped"] = bool(off["capped"])
    rec["mirror_ok"] = (off["n_trees"] == m["n_trees"]
                        and off["n_hidden"] == m["e_min"]
                        and bool(off["capped"]) == bool(m["capped"]))

    # V1/V2/V3/V6/V7 (light) vs V4/V5 (完整) 成本分離
    t0 = T()
    S.verify_all(off, full_genos, partial_strs, k, light=True)
    rec["t_verify_light"] = T() - t0

    rec["t_verify_V4"] = None
    rec["t_verify_V5"] = None
    if run_v45 and not off["capped"]:
        est = math.comb(m["P_extra_pool"], max(m["e_min"] - 1, 0)) + \
              math.comb(m["P_extra_pool"], m["e_min"])
        if est <= 3_000_000:
            t0 = T()
            S.V4_minimality_independent(off, full_genos, partial_strs, k)
            rec["t_verify_V4"] = T() - t0
            t0 = T()
            S.V5_enumerate_all_independent(off, full_genos, partial_strs, k)
            rec["t_verify_V5"] = T() - t0
        else:
            rec["v45_skipped_reason"] = f"est {est} combos > 3e6 guard"
    return rec


def rand_fulls(k, f, rng):
    seen, out = set(), []
    while len(out) < f:
        g = "".join(rng.choice("RA") for _ in range(k))
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def window_partials(k, p, width, rng):
    """模擬 read 只穿過連續 width 個位點,其餘 X(未覆蓋)。"""
    out = []
    for _ in range(p):
        st = rng.randrange(0, max(k - width, 0) + 1)
        s = ["X"] * k
        for i in range(st, min(st + width, k)):
            s[i] = rng.choice("RA")
        out.append("".join(s))
    return out


def main():
    rng = random.Random(20260725)
    cases = []

    # A. golden 手算案例(最小規模基準)
    cases.append(("A1_two_siblings", ["AR", "RA"], [], 2))
    cases.append(("A2_single_double", ["AA"], [], 2))
    cases.append(("A3_three_gamete", ["AR", "RA", "AA"], [], 2))
    cases.append(("A4_partial_IDP", [], ["AX", "XA"], 2))

    # B. full-only,k 遞增(觀察 k 對各階段的影響)
    for k in [2, 3, 4, 5, 6, 7, 8]:
        cases.append((f"B_fullonly_k{k}_f3", rand_fulls(k, min(3, 2 ** k), rng), [], k))
    # C. full-only,k=8,full 數遞增
    for f in [1, 2, 4, 8, 16, 32]:
        cases.append((f"C_fullonly_k8_f{f}", rand_fulls(8, f, rng), [], 8))
    # D. partial-heavy(read 窗寬 3),k=8,partial 數遞增 — 最貼近 pure-partial 區
    for p in [1, 2, 4, 8, 16]:
        cases.append((f"D_partial_k8_w3_p{p}", [], window_partials(8, p, 3, rng), 8))
    # E. 混合:少量 full + 多 partial
    for p in [2, 4, 8]:
        cases.append((f"E_mix_k8_f2_p{p}", rand_fulls(8, 2, rng),
                      window_partials(8, p, 4, rng), 8))
    # F. 極端:k=8 單一 all-ALT(最大 hidden 需求)
    cases.append(("F_k8_all_alt", ["AAAAAAAA"], [], 8))
    cases.append(("F_k6_all_alt", ["AAAAAA"], [], 6))
    cases.append(("F_k4_all_alt", ["AAAA"], [], 4))

    out = []
    for name, full, part, k in cases:
        try:
            t0 = T()
            rec = bench_case(name, full, part, k)
            rec["wall_total_s"] = T() - t0
            out.append(rec)
            print(f"[ok] {name:24s} P={rec['P_extra_pool']:4d} e_min={rec['e_min']} "
                  f"n_trees={rec['n_trees']:>10} cand={rec['candidates_tested']:>8} "
                  f"t5={rec['t5_level_search']*1000:8.2f}ms t7={rec['t7_materialize']*1000:8.2f}ms "
                  f"mirror={'OK' if rec['mirror_ok'] else 'MISMATCH'}", flush=True)
        except Exception as ex:  # noqa: BLE001
            print(f"[FAIL] {name}: {type(ex).__name__}: {ex}", flush=True)
            out.append({"case": name, "error": f"{type(ex).__name__}: {ex}"})

    base = ("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
            "efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad/complexity_bench")
    with open(base + ".json", "w") as fh:
        json.dump({"cases": out,
                   "solver": "tree_enumeration_solver.py (2026-07-06)",
                   "params": {"extra_cap": 4, "per_level_budget": 150000,
                              "materialize": "unlimited (production ANALYSIS_TREE_CAP=0)"}},
                  fh, indent=1)

    cols = ["case", "k", "n_full", "n_partial", "subcube_total_members", "n_maximal",
            "pool_subset_ops", "n_pool", "P_extra_pool", "e_min", "candidates_tested",
            "n_feasible_N", "n_trees", "trees_materialized", "capped", "mirror_ok",
            "t3_subcubes", "t4_pool", "t5_level_search", "t6_analytic_count",
            "t7_materialize", "t_verify_light", "t_verify_V4", "t_verify_V5",
            "t_total_no_verify"]
    with open(base + ".tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            if "error" in r:
                continue
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    n_mm = sum(1 for r in out if r.get("mirror_ok") is False)
    print(f"\nwrote {base}.json / .tsv   cases={len(out)}  mirror_mismatch={n_mm}")


if __name__ == "__main__":
    main()
