#!/usr/bin/env python3
"""
build_data.py — 計算量漏斗教學頁的資料層（純計算 + 自我驗證，零手打數字）

設計原則（§13-A 由構造防捏造）：
  * 每個數字都是「純數學重算」或「從既有 verified JSON 讀出」，本檔不硬寫任何 empirical 數字。
  * 每個數字都附 formula（可手算）+ verify（可一行重跑）+ source（檔案路徑）。
  * 最後跑 SELF-CHECK：把重算的純數學值與獨立來源（complexity bench）比對，不符即 exit 1。

輸出: 20260726_complexity_funnel_teaching_01.data.json
用法: python3 build_data.py            # 產生 data.json + 印驗證表
      python3 build_data.py --check    # 只驗證不寫檔
"""
import json
import math
import sys
from itertools import combinations
from pathlib import Path

REPO = Path(__file__).resolve().parents[6]
HERE = Path(__file__).resolve().parent

SRC = {
    "bench": REPO / "docs/reports/in_progress/2026/07/20260725_solver_complexity_stepwise_01"
                    "/20260725_solver_complexity_stepwise_01.data.json",
    "census": REPO / "docs/methodology/_assets/20260627_subclone_4axis_teaching/data/funnel_census_HCC1395.json",
    "layered": REPO / "docs/methodology/_assets/20260618_subcluster_pilot/layered_reconstruction_HCC1395.json",
    "ccf": REPO / "docs/methodology/_assets/20260627_subclone_4axis_teaching/data"
                  "/ccf_tree_weighting_full_observe.json",
    "solver": REPO / "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py",
    "spec": REPO / "docs/methodology/20260704_formal_problem_statement_topology_from_cooccurrence_01.md",
}


def load(key):
    p = SRC[key]
    if not p.exists():
        sys.exit(f"REFUSE: source missing -> {p}")
    return json.loads(p.read_text())


# ---------------------------------------------------------------- 純數學重算

def digits(n: int) -> int:
    """大整數位數（不轉 float，避免溢位）。"""
    return len(str(n))


def sci(n: int, sig=4) -> str:
    """大整數 → 科學記號字串（純整數運算，無 float 溢位）。"""
    s = str(n)
    if len(s) <= sig:
        return s
    return f"{s[0]}.{s[1:sig]}e+{len(s) - 1}"


def hypercube_stats(k: int) -> dict:
    """有向布林超立方體 {0,1}^k（只准 0->1 單位翻轉）的基本量。"""
    v = 2 ** k                      # 頂點
    e = k * 2 ** (k - 1)            # 有向 unit-flip 邊（每邊 = 選一個 0 位翻成 1）
    return {"k": k, "vertices": v, "edges": e,
            "all_subsets": 2 ** v, "all_subsets_digits": digits(2 ** v)}


def arborescence_count_hypercube(k: int) -> int:
    """超立方體上以 0^k 為根的 spanning arborescence 數。

    定理：每個非根頂點 x 獨立選一個 unit-pred（少一個 1 的鄰居），共 popcount(x) 種；
    邊皆 popcount 遞增 → 天然無環且 root-connected。故
        A(k) = prod_{x != 0} popcount(x) = prod_{j=1..k} j ^ C(k,j)
    """
    a = 1
    for j in range(1, k + 1):
        a *= j ** math.comb(k, j)
    return a


def arborescence_count_bruteforce(k: int) -> int:
    """A(k) 的獨立重算（逐頂點數 pred，不用乘冪公式）——交叉驗證用。"""
    a = 1
    for x in range(1, 2 ** k):
        a *= bin(x).count("1")
    return a


def complete_digraph_arborescences(n: int) -> int:
    """對照組：n 個帶標號頂點的完全圖，指定根的 spanning arborescence 數 = n^(n-1)
    （Cayley 公式的有根版；未加任何超立方體/單調限制）。"""
    return n ** (n - 1)


def closed_sets_bruteforce(k: int) -> int:
    """暴力數「合法節點集」：含 root、且每個非根成員至少有一個 unit-pred 也在集合內
    （= solver 的 _is_closed；此條件 ⟺ 存在 root-connected arborescence）。
    只對 k<=4 可行（子集空間 2^(2^k)）。"""
    verts = list(range(2 ** k))
    preds = {x: [x & ~(1 << b) for b in range(k) if x >> b & 1] for x in verts}
    cnt = 0
    nonroot = [x for x in verts if x]
    for mask in range(2 ** len(nonroot)):
        S = {0}
        for i, x in enumerate(nonroot):
            if mask >> i & 1:
                S.add(x)
        if all(any(p in S for p in preds[x]) for x in S if x):
            cnt += 1
    return cnt


def closed_sets_leveldp(k: int) -> int:
    """同上但用分層 DP（狀態=上一層被選中的頂點子集），可算到 k=5。
    正確性由 SELF-CHECK 對 k<=4 與暴力版逐一比對保證。"""
    levels = [[x for x in range(2 ** k) if bin(x).count("1") == lv] for lv in range(k + 1)]
    dp = {frozenset(levels[0]): 1}          # level0 = {root}，必選
    for lv in range(1, k + 1):
        cur, nxt = levels[lv], {}
        for prev, ways in dp.items():
            ok = [y for y in cur if any((y & ~(1 << b)) in prev for b in range(k) if y >> b & 1)]
            for m in range(2 ** len(ok)):
                T = frozenset(ok[i] for i in range(len(ok)) if m >> i & 1)
                nxt[T] = nxt.get(T, 0) + ways
        dp = nxt
    return sum(dp.values())


def cand_curve(max_P: int, extra_cap: int, budget: int):
    """S5 逐層搜尋的候選節點集總數（忠實鏡像 solver:144-157 的 fast-cap 行為）：
        E(P) = min(extra_cap, min{e : C(P,e) > budget} - 1)
        Cand(P) = sum_{e=0..E(P)} C(P,e)
    """
    out = []
    for P in range(max_P + 1):
        E = extra_cap
        for e in range(0, extra_cap + 1):
            if math.comb(P, e) > budget:
                E = e - 1
                break
        E = max(E, 0)
        out.append({"P": P, "E": E, "cand": sum(math.comb(P, e) for e in range(E + 1))})
    return out


# ---------------------------------------------------------------- 組裝

def build():
    bench, census, layered, ccf = load("bench"), load("census"), load("layered"), load("ccf")
    params = bench["_meta"]["params"]
    K_MAX, EXTRA_CAP, BUDGET = params["MAX_SNV"], params["extra_cap"], params["per_level_budget"]

    checks = []

    def ck(name, got, want, note=""):
        ok = (got == want)
        checks.append({"name": name, "got": str(got)[:60], "want": str(want)[:60], "pass": ok, "note": note})
        return ok

    # ---- L2 超立方體 ----
    cube = [hypercube_stats(k) for k in range(1, K_MAX + 1)]
    ck("超立方體頂點 k=8 = 256", cube[-1]["vertices"], 256, "2^8")
    bench_ss = {r["k"]: r["vertices"] for r in bench["state_space"]}
    for c in cube:
        if c["k"] in bench_ss:
            ck(f"頂點數 k={c['k']} 對得上 bench state_space", c["vertices"], bench_ss[c["k"]])

    # ---- L3 演化樹限制 ----
    arb = []
    bench_A = {r["k"]: r["value"] for r in bench["A_max"]}
    for k in range(1, K_MAX + 1):
        a_formula = arborescence_count_hypercube(k)
        a_brute = arborescence_count_bruteforce(k)
        ck(f"A(k={k}) 乘冪公式 == 逐頂點暴力", a_formula, a_brute, "prod j^C(k,j) vs prod popcount(x)")
        if k in bench_A:
            ck(f"A(k={k}) 與 bench A_max 一致", str(a_formula), bench_A[k], "獨立來源交叉驗證")
        n = 2 ** k
        free = complete_digraph_arborescences(n)
        arb.append({
            "k": k, "n_vertices": n,
            "unconstrained": str(free), "unconstrained_digits": digits(free), "unconstrained_sci": sci(free),
            "hypercube": str(a_formula), "hypercube_digits": digits(a_formula), "hypercube_sci": sci(a_formula),
            "cut_orders": digits(free) - digits(a_formula),
        })

    # ---- L4 合法節點集（closed）----
    closed = []
    for k in range(1, 6):
        dp = closed_sets_leveldp(k)
        bf = closed_sets_bruteforce(k) if k <= 4 else None
        if bf is not None:
            ck(f"closed-set 計數 k={k}：分層DP == 暴力", dp, bf, "DP 正確性由暴力背書")
        closed.append({"k": k, "vertices": 2 ** k, "all_subsets": str(2 ** (2 ** k)),
                       "all_subsets_digits": digits(2 ** (2 ** k)),
                       "closed": dp, "bruteforce_agrees": (bf == dp) if bf is not None else None,
                       "frac_ppm": round(dp / 2 ** (2 ** k) * 1e6, 4) if k <= 4 else None})

    # ---- L5 H 選擇 ----
    curve = cand_curve(2 ** K_MAX - 1, EXTRA_CAP, BUDGET)
    peak = max(curve, key=lambda r: r["cand"])
    b5 = bench["step5_bound"]
    ck("Cand 全域最大值", peak["cand"], b5["global_max"], "獨立重算 vs bench")
    ck("Cand 最大值出現的 P", peak["P"], b5["at_P"])
    ck("該點的 E", peak["E"], b5["at_E"])
    max_P_e4 = max(r["P"] for r in curve if r["E"] >= 4)
    max_P_e3 = max(r["P"] for r in curve if r["E"] >= 3)
    ck("允許 e=4 的最大 P", max_P_e4, b5["max_P_allowing_e4"])
    ck("允許 e=3 的最大 P", max_P_e3, b5["max_P_allowing_e3"])
    saw = [{"P": p, **{kk: vv for kk, vv in next(r for r in curve if r["P"] == p).items() if kk != "P"}}
           for p in (44, 45, 46, 96, 97, 98, 255)]

    # ---- L1 分區（empirical，來自 census）----
    f = census["funnel_sSNV"]
    ck("census 六層 funnel 加總 == universe",
       f["L2_out_of_scope_chrXY"] + f["L3_positional_singleton"] + f["L4_cap_excluded_densest8"]
       + f["L5_read_unsupported"] + f["L6_retained_sSNV"], f["L1_all_pass_universe"], "分層互斥且窮盡")
    ck("chr1-22 = universe - chrXY",
       f["L1_all_pass_universe"] - f["L2_out_of_scope_chrXY"], f["autosomal_chr1_22"])

    # 未分區的全域狀態空間（純數學）
    N_auto = f["autosomal_chr1_22"]
    global_digits = int(N_auto * math.log10(2)) + 1

    # 分區後的實際狀態空間總和：sum 2^k_i over solver units
    detail = layered["detail"]
    ks = [u["n_sSNV"] for u in detail]
    sum_states = sum(2 ** x for x in ks)
    ck("layered detail 筆數 == n_detail_units", len(detail), layered["n_detail_units"])
    bench_k_hist = {r["v"]: r["n"] for r in bench["real"]["k"]["hist"]}
    my_k_hist = {}
    for x in ks:
        my_k_hist[x] = my_k_hist.get(x, 0) + 1
    ck("k 直方圖與 bench real.k.hist 一致", my_k_hist, bench_k_hist, "同一份 detail 獨立重數")

    # ---- L6 樹的量（empirical）----
    nt = bench["real"]["n_trees"]
    tv = bench["theory_vs_real"]
    ck("理論 A(8) 位數 == bench 記錄", digits(arborescence_count_hypercube(8)), tv["A_max_k8_digits"])

    # ---- L7 VAF/CCF ----
    cs = ccf["summary"]
    ck("CCF: 破等 + 維持均勻 == ambiguous 總數",
       cs["tie_broken(≥0.6)"] + cs["stay_uniform"], cs["n_ambiguous_units"])
    ck("CCF: broke_frac 重算",
       round(cs["tie_broken(≥0.6)"] / cs["n_ambiguous_units"], 3), round(cs["broke_frac"], 3))
    ck("CCF: winner clean frac 重算",
       round(cs["winner_pigeonhole_clean"] / cs["tie_broken(≥0.6)"], 3), round(cs["winner_clean_frac_of_broke"], 3))
    ck("CCF ambiguous 數 == layered determinacy_lineage ambiguous",
       cs["n_ambiguous_units"], layered["L1_ssnv_algorithm"]["determinacy_lineage"]["ambiguous_structure(多完成/多結構)"])

    all_pass = all(c["pass"] for c in checks)

    data = {
        "_meta": {
            "title": "從全基因組到 Group Steiner Arborescence — 計算量漏斗",
            "date": "2026-07-26",
            "task_type": "F demo / teaching（教學頁；數字全為既有 verified 來源重算，非新分析）",
            "sample_scope": "HCC1395 單樣本；純數學層與樣本無關",
            "generated_by": "build_data.py",
            "params": {"MAX_SNV": K_MAX, "extra_cap": EXTRA_CAP, "per_level_budget": BUDGET},
            "sources": {k: str(v.relative_to(REPO)) for k, v in SRC.items()},
            "self_check_pass": all_pass,
            "n_checks": len(checks),
        },
        "self_check": checks,
        "L0_global": {
            "n_autosomal_sSNV": N_auto,
            "genotype_space_digits": global_digits,
            "formula": r"|X_{global}| = 2^{N},\quad N = 80{,}234",
            "verify": "python3 -c \"import math;print(int(80234*math.log10(2))+1)\"",
        },
        "L1_partition": {
            "funnel": f,
            "n_solver_units": layered["n_detail_units"],
            "sum_states_after_partition": sum_states,
            "sum_states_digits": digits(sum_states),
            "upper_bound_units_times_256": layered["n_detail_units"] * 2 ** K_MAX,
            "k_hist": [{"k": kk, "n": my_k_hist[kk], "pct": round(my_k_hist[kk] / len(ks) * 100, 2),
                        "states_2k": 2 ** kk} for kk in sorted(my_k_hist)],
            "formula": r"\sum_i 2^{k_i}\ \le\ U\cdot 2^{k_{max}}",
            "verify": "python3 -c \"import json;d=json.load(open(SRC));print(sum(2**u['n_sSNV'] for u in d['detail']))\"",
        },
        "L2_hypercube": {"rows": cube,
                         "formula": r"|V|=2^k,\quad |E|=k\,2^{k-1},\quad \#\text{任意節點子集}=2^{2^k}",
                         "verify": "手算：k=8 → 256 頂點、8·128=1024 條有向邊、2^256 個子集"},
        "L3_tree_constraint": {"rows": arb,
                               "formula": r"A(k)=\prod_{x\neq 0}\mathrm{popcount}(x)=\prod_{j=1}^{k} j^{\binom{k}{j}}"
                                          r"\quad\text{vs}\quad n^{\,n-1},\ n=2^k",
                               "verify": "本檔以「乘冪公式」與「逐頂點 popcount 連乘」兩法獨立算並比對（SELF-CHECK）"},
        "L4_closed": {"rows": closed,
                      "formula": r"N\ \text{合法} \iff \forall x\in N\setminus\{0\}\ \exists j\in x:\ x\setminus\{j\}\in N",
                      "verify": "k≤4 暴力枚舉全部 2^(2^k) 子集；k=5 用分層 DP，且 DP 對 k≤4 與暴力逐一吻合"},
        "L5_hidden_choice": {
            "curve": curve,
            "sawtooth": saw,
            "peak": peak,
            "max_P_allowing_e4": max_P_e4,
            "max_P_allowing_e3": max_P_e3,
            "e_min_hist": bench["real"]["e_min"]["hist"],
            "e_min_stats": bench["real"]["e_min"]["stats"],
            "capped": bench["real"]["capped"],
            "pool_hist": bench["real"]["p"]["hist"],
            "pool_stats": bench["real"]["p"]["stats"],
            "stage_share": bench["stage_share"],
            "stage_total_ms": bench["stage_share_total_ms"],
            "walltime": bench["step5_bound"]["worst_walltime_s"],
            "throughput_us": bench["throughput_us_per_candidate"],
            "formula": r"\mathrm{Cand}(P)=\sum_{e=0}^{E(P)}\binom{P}{e},\quad "
                       r"E(P)=\min\!\Big(4,\ \min\{e:\tbinom{P}{e}>150000\}-1\Big)",
            "verify": "python3 -c \"from math import comb;print(sum(comb(45,e) for e in range(5)))\"  # 164221",
        },
        "L6_trees": {
            "n_trees_stats": nt["stats"], "n_trees_total": nt["total"], "n_trees_tail": nt["tail"],
            "theory_vs_real": tv,
            "bench_cases": bench["bench_cases"],
            "formula": r"n_{\text{trees}}=\sum_{N\in\mathcal N_{\min}}\ \prod_{x\in N\setminus\{0\}}"
                       r"\big|\{j\in x:\ x\setminus\{j\}\in N\}\big|",
            "verify": "solver V5_enumerate_all_independent 獨立重算全樹數並與回報值比對（不信任自報）",
        },
        "L7_vaf": {
            "summary": cs,
            "formula": r"P(T\mid \text{reads})\propto \prod \text{read-count}\ \text{(pigeonhole: 祖先 read 數} \ge \text{後裔)}",
            "verify": "破等數 + 維持均勻數 == ambiguous 總數；比例逐項重算（SELF-CHECK 4 項）",
        },
        "L8_group_steiner": {
            "tree_level": census["tree_level"],
            "citations": [
                {"key": "foulds1982steiner", "what": "Steiner problem in phylogeny is NP-complete",
                 "where": "Adv. Appl. Math. 3(1):43–49 (1982)", "doi": "10.1016/S0196-8858(82)80004-3"},
                {"key": "garg2000group", "what": "Group Steiner Tree 多對數近似 O(log^2 n log k)",
                 "where": "J. Algorithms 37(1):66–84 (2000)", "doi": "10.1006/jagm.2000.1096"},
                {"key": "halperin2003polylog", "what": "Group Steiner 不可近似下界 Ω(log^{2-ε} n)",
                 "where": "STOC 2003:585–594", "doi": "10.1145/780542.780628"},
                {"key": "charikar1999approximation", "what": "Directed Steiner 近似 O(k^ε)",
                 "where": "J. Algorithms 33(1):73–91 (1999)", "doi": "10.1006/jagm.1999.1042"},
                {"key": "mahapatra2025parameterized", "what": "超立方體上 Steiner arborescence 的 FPT 演算法",
                 "where": "Acta Inf. 62(1):6 (2025)", "doi": "10.1007/s00236-024-00474-8"},
                {"key": "peer2004incomplete", "what": "Incomplete Directed Perfect Phylogeny（poly→linear）",
                 "where": "SIAM J. Comput. 33(3):590–607 (2004)", "doi": "10.1137/S0097539702406510"},
                {"key": "gusfield1991efficient", "what": "directed perfect phylogeny 多項式建構",
                 "where": "Networks 21(1):19–28 (1991)", "doi": "10.1002/net.3230210104"},
            ],
        },
    }
    return data, checks, all_pass


def main():
    data, checks, ok = build()
    w = max(len(c["name"]) for c in checks) + 2
    print("=" * (w + 46))
    print("SELF-CHECK — 每個純數學重算 vs 獨立來源")
    print("=" * (w + 46))
    for c in checks:
        print(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['name']:<{w}} got={c['got']:<18} want={c['want']}")
    print("-" * (w + 46))
    print(f"  {sum(c['pass'] for c in checks)}/{len(checks)} PASS")
    if not ok:
        sys.exit("SELF-CHECK FAILED — 不寫出 data.json（§13-A refuse-on-mismatch）")
    if "--check" in sys.argv:
        print("  (--check 模式：不寫檔)")
        return
    out = HERE / "20260726_complexity_funnel_teaching_01.data.json"
    out.write_text(json.dumps(data, ensure_ascii=False, indent=1))
    print(f"  wrote {out.name}  ({out.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
