#!/usr/bin/env python3
"""
complexity_theory.py — 理論最大量精確計算(大整數,不估算)

對應 tree_enumeration_solver.py 各階段的 worst-case 上界。
參數同生產:k<=MAX_SNV=8、extra_cap=4、per_level_budget=150000、materialize 無上限。

定理(Step 6/7 的基礎):
  給定 unit-flip-closed 節點集 N (含 ROOT=∅),單調生成樹(arborescence)數
      A(N) = Π_{x∈N, x≠∅} d_N(x),  d_N(x) = #{ j∈x : x∖{j} ∈ N }
  證明:每個非根節點獨立選一個 unit-pred;沿 parent 指標 popcount 嚴格遞減
       → 必無環且終止於唯一無父節點 ∅ → 每組選擇恰對應一棵合法 arborescence,
       不同選擇給不同邊集 → 雙射。∎
  推論:計數是 O(|N|·k) 多項式;列舉是 Θ(A(N)) 指數。
"""
import json
from math import comb

OUT = ("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
       "efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad/complexity_theory.json")

MAXK = 8
EXTRA_CAP = 4
BUDGET = 150_000

res = {"params": {"MAX_SNV": MAXK, "extra_cap": EXTRA_CAP,
                  "per_level_budget": BUDGET}}

# ---------- 1. 狀態空間 ----------
res["state_space"] = [{"k": k, "vertices_2^k": 2 ** k,
                       "max_pool": 2 ** k, "max_P_extra_pool": 2 ** k - 1}
                      for k in range(1, MAXK + 1)]

# ---------- 2. Step 3/4:group 與 pool 上界 ----------
# subcube 成員數 = 2^u (u=未覆蓋位點數) <= 2^k
# pool 去重後 <= 2^k;但建構操作數 = Σ_{m∈maximal} 2^{|m|} <= |maximal|·2^k
res["step34_bounds"] = {
    "max_subcube_members_per_partial": 2 ** MAXK,
    "max_pool_after_dedup": 2 ** MAXK,
    "pool_build_ops_formula": "sum_{m in maximal} 2^{|m|} <= (f + p*2^k) * 2^k",
    "example_k8_p16_f4_ops_upper": (4 + 16 * 2 ** 8) * 2 ** 8,
}

# ---------- 3. Step 5:逐層搜尋的候選數上界 ----------
# 迴圈在「第一個 C(P,e) > BUDGET 的 e」中斷。最壞情況(到最後允許層才找到解或全無解)
# 測試候選總數 = Σ_{e=0}^{E(P)} C(P,e),E(P)=min(extra_cap, 首個超budget的e - 1)
def E_of_P(P):
    for e in range(0, min(EXTRA_CAP, P) + 1):
        if comb(P, e) > BUDGET:
            return e - 1
    return min(EXTRA_CAP, P)

rows = []
for P in range(0, 2 ** MAXK):
    E = E_of_P(P)
    total = sum(comb(P, e) for e in range(0, max(E, -1) + 1)) if E >= 0 else 0
    rows.append({"P": P, "E_max_level": E, "candidates_tested_total": total})
best = max(rows, key=lambda r: r["candidates_tested_total"])
res["step5_candidate_bound"] = {
    "formula": "sum_{e=0}^{E(P)} C(P,e), E(P)=min(4, first e with C(P,e)>150000, minus 1)",
    "global_max_candidates_tested": best["candidates_tested_total"],
    "attained_at_P": best["P"],
    "attained_E": best["E_max_level"],
    "max_P_allowing_e4": max(r["P"] for r in rows if r["E_max_level"] >= 4),
    "max_P_allowing_e3": max(r["P"] for r in rows if r["E_max_level"] >= 3),
    "at_P_255": next(r for r in rows if r["P"] == 255),
    "selected_profile": [r for r in rows if r["P"] in
                         (10, 20, 45, 46, 50, 96, 97, 98, 128, 200, 255)],
}
# 每候選成本:_covers O(f + Σ_i min(|N|,|sub_i|)) + _is_closed O(Σ_{x∈N}|x|) <= O(|N|·k)
res["step5_per_candidate_cost"] = {
    "covers": "O(f + p * min(|N|, 2^k))",
    "is_closed": "O(|N| * k)",
    "max_N_size_under_extra_cap": "1 + f + 4",
}

# ---------- 4. Step 6/7:樹數的精確理論最大 ----------
# A(N) 在 N = 整個超立方體時最大:每個 popcount=j 的節點有 j 個 unit-pred,
# 且 popcount=j 的節點有 C(k,j) 個 → A_max(k) = Π_{j=1}^{k} j^{C(k,j)}
def A_max_full_cube(k):
    v = 1
    for j in range(1, k + 1):
        v *= j ** comb(k, j)
    return v

amax = []
for k in range(1, MAXK + 1):
    a = A_max_full_cube(k)
    amax.append({"k": k, "vertices": 2 ** k,
                 "A_max_full_hypercube": str(a),
                 "digits": len(str(a)),
                 "requires_hidden_nodes_e": 2 ** k - 1})  # 若 f=0
res["step67_absolute_max_trees"] = {
    "theorem": "A(N) = prod_{x in N, x != empty} d_N(x)",
    "max_formula": "A_max(k) = prod_{j=1}^{k} j^C(k,j)  (N = full hypercube)",
    "per_k": amax,
    "note": "此絕對上界要求 N 含全部 2^k 頂點,即 f+e = 2^k-1;"
            " extra_cap=4 下需 f >= 2^k-5 個相異完整觀測,實務不可達",
}

# 受 caps 限制的可達上界:|N| = 1 + f + e, e <= 4, d_N(x) <= min(|x|, k)
# 對固定 |N|=m,A(N) <= (k)^(m-1);更緊:節點 popcount 受 |N| 限制
def A_bound_under_cap(f, k, e=EXTRA_CAP):
    m_nonroot = f + e
    return k ** m_nonroot

res["step67_reachable_bound_under_caps"] = {
    "formula": "A(N) <= k^(f+extra_cap) per feasible N; total <= |feasible_N| * k^(f+4)",
    "examples": [{"k": 8, "f": f, "bound_per_N": str(A_bound_under_cap(f, 8)),
                  "digits": len(str(A_bound_under_cap(f, 8)))}
                 for f in [1, 2, 3, 5, 8, 12]],
    "max_feasible_N": "<= C(P, e_min) <= 150000 (budget)",
}

# ---------- 5. 全域最壞情況組合 ----------
res["global_worst_case"] = {
    "step5_candidates": best["candidates_tested_total"],
    "step6_count_cost": "O(|feasible_N| * |N| * k^2)  — 多項式,即使 A(N) 天文數字",
    "step7_materialize_cost": "Theta(total_A) — 生產 ANALYSIS_TREE_CAP=0 表示無上限",
    "why_step7_is_the_real_risk": "計數多項式但列舉指數;capped 區已被 Step5 擋掉,"
                                  "但未 capped 卻 A(N) 大的區會在 Step7 爆炸",
}

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1, ensure_ascii=False)

print("=== 狀態空間 ===")
for r in res["state_space"]:
    print(f"  k={r['k']}  2^k={r['vertices_2^k']:>4}  max_P={r['max_P_extra_pool']:>4}")
print("\n=== Step5 候選數上界 ===")
print(f"  全域最大候選測試數 = {best['candidates_tested_total']:,} (P={best['P']}, E={best['E_max_level']})")
print(f"  允許到 e=4 的最大 P = {res['step5_candidate_bound']['max_P_allowing_e4']}")
print(f"  允許到 e=3 的最大 P = {res['step5_candidate_bound']['max_P_allowing_e3']}")
print(f"  P=255 時: E={res['step5_candidate_bound']['at_P_255']['E_max_level']}, "
      f"tested={res['step5_candidate_bound']['at_P_255']['candidates_tested_total']:,}")
for r in res["step5_candidate_bound"]["selected_profile"]:
    print(f"    P={r['P']:>3}  E={r['E_max_level']}  tested={r['candidates_tested_total']:>9,}")
print("\n=== 樹數絕對理論最大 A_max(k) = prod j^C(k,j) ===")
for r in amax:
    print(f"  k={r['k']}  2^k={r['vertices']:>4}  A_max={r['A_max_full_hypercube'][:28]}"
          f"{'...' if r['digits'] > 28 else ''}  ({r['digits']} 位數)")
print("\n=== caps 下每個 N 的可達上界 k^(f+4), k=8 ===")
for r in res["step67_reachable_bound_under_caps"]["examples"]:
    print(f"  f={r['f']:>2}  bound={r['bound_per_N']:>22}  ({r['digits']} 位數)")
print(f"\nwrote {OUT}")
