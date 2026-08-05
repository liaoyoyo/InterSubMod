<!--
建立時間：2026-07-16
目標：獨立稽核 Boolean-hypercube rooted group-Steiner 中 partial-read 的數學語意、實際程式行為與可擴充解法
處理範圍：新 symbolic MILP、凍結 canonical v5 solver、單元測試與 RRA+AAA+RAX 小例子；不包含全樣本數據結論
關聯檔案：InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py；canonical v5 source bundle 的 003_tree_enumeration_solver.py
-->

# Partial-read Hypercube group-Steiner 方法獨立稽核

用 SCQA + Step→Verify：先釐清「partial completion」與「候選樹」不是同一層，再用程式碼、測試與反例驗證正確流程。

**結論：partial pattern 有 `u` 個 `X` 時，在完整 `k` 維立方體中概念上對應 `2^u` 個 completion；但正確求解不是把每種 completion 分開跑、看到任一成功就停止。應把每個達 `structural_exact_pattern_minread` 的 distinct partial pattern保留為一個 group constraint `N ∩ G(p) ≠ ∅`，先求全域最少 hidden vertices，再列出所有同分 optimal vertex sets；低於門檻的 informative molecules不形成structural constraint，但仍進read-pattern likelihood。likelihood只可排序不同 vertex sets；相同 vertex set 的不同 parent edges 對目前 snapshot read/VAF 資訊不可辨。**（影響：高；程式行為信心：高；跨樣本生物結論：本文件未評估）

## 1. 稽核範圍與證據版本

| 角色 | 輸入路徑 | SHA-256 | 稽核重點 |
|---|---|---|---|
| 新 symbolic solver | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py` | `9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95` | symbolic mask、exact reductions、MILP、optimal vertex-set enumeration、mixture likelihood |
| 凍結 canonical solver | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/source_bundle/files/imported/003_tree_enumeration_solver.py` | `36727f4e1d8d7ce8abf869606211c93d8c1a0506dd7d142e855863c412ca0d61` | canonical v5 實際 subcube expansion、最小 hidden 搜尋、parent enumeration |
| 單元測試 | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py` | `a86c23e981cf92b2a41c1d86e778d5bba6e6a56d0ba1743180de5534fb3313de` | symbolic/explicit 等價、independent exhaustive oracle、exact reductions、例子、likelihood invariance；2026-07-17 current targeted replay 23/23 PASS |
| current-solver identity-bound pilot | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot_v3_identity_bound/pilot_receipt.json` | `5cab8c37f31908d13b0c921619264598a4e6e3ef29cbcf723d1a6c421b6297c8` | receipt 綁定 current runner／solver SHA；k≤8 exhaustive、legacy cross-check、k9–12 pilot；明示 `PILOT_NOT_FINAL_VALIDATION` |

本稽核服務 G3（read-level 訊息建模）、G4（可重現驗證）與 G5（可外部稽核的方法契約）。

## 2. 問題定義：一條 partial read 是一個 group，不是一棵完成後的樹

### 2.1 狀態、partial pattern 與 completion 數

令區域含 `k` 個 sSNV，mutation-state vertex 為

\[
v\in\{0,1\}^k,\qquad 0=R,\;1=A.
\]

一條 read pattern `p ∈ {R,A,X}^k` 以三個 bitmask 表示：

- `fixed_mask M_p`：`R/A` 已觀測位置；
- `alt_mask A_p`：已觀測為 `A` 的位置；
- `free_mask`：`X` 未觀測位置。

與 `p` 相容的節點群組為

\[
G(p)=\{v\in U:((v\oplus A_p)\;\&\;M_p)=0\}.
\]

若完整 `k` 維 hypercube 都屬於 solver universe，且有 `u` 個 `X`，則

\[
|G(p)|=2^u.
\]

這是**概念 completion 數**。新程式 `SymbolicPattern.from_string()` 只存三個 masks；`compatible()` 用一個 bitwise expression 判斷節點是否屬於 group。每次 `_build_problem()` 重建一個 MILP 時，會對每個 reduction 後仍 active 的 group 呼叫 `_compatible_vertices()`，列出 active universe 內至多 `2^u_eff` 個相容 vertex indices，寫成該次 sparse hit row 的非零係數。也就是說，indices 會在**每次 MILP build × 每個 reduced group** materialize，並非整個 workflow 只產生一次；但它們仍只是同一條限制式的係數，不是 `2^u_eff` 個分開執行的 tree worlds，也不會跨 reads 建立 joint-completion 笛卡兒積。source-bound method contract據此分開宣告`active_compatible_vertex_indices_materialized_for_sparse_rows=true`、`completion_wise_tree_worlds_materialized=false`與`cross_read_cartesian_products_materialized=false`（`hypercube_exact.py:36-91,208-282,285-408`）。

### 2.2 canonical v5 有一個需明示的 universe 假設

凍結 canonical solver 的 `U` 不是必然為全部 `k` loci。它只保留「至少在某個 full/partial pattern 曾明確觀測為 A」的位置；某位置若只出現在 `X`、從未出現 `A`，會被視為固定 0（`003_tree_enumeration_solver.py:120-122`）。因此其有效 completion 數是

\[
2^{|free\_mask\cap U|},
\]

可能小於 `2^u`。這是**現行模型假設**，不是 partial read 的一般數學定理。新 solver 可用 `universe_mode=all_loci` 保留所有維度，或以 `observed_alt` 重現 canonical 行為（`hypercube_exact.py:119-125,185-190,302-307`）。若要對論文宣稱「每個 X 皆有兩種可能」，必須使用 full-universe；若維持 canonical 口徑，應寫成「在 active observed-ALT universe 內的有效可能數」。

### 2.3 Production 中哪些 patterns 會形成 structural constraints

M2 production 先在同一個 `component×HP family×exact PS` unit 內，依**完整 exact R/A/X 字串**計數。只有 `count >= structural_exact_pattern_minread` 的 distinct informative patterns 會進 structural stage：其中無 `X` 者成為 mandatory full states，有 `X` 者各自形成一個 group constraint。`AXR`、`ARX`、`AXX` 即使彼此相容，也不會合併湊門檻。低於門檻的 informative molecules 不會形成結構限制，但仍連同達門檻 molecules 一起進 quality-aware likelihood；所以門檻決定「哪些 patterns 可約束候選結構」，不是丟棄低支持 reads（`build_m2_patterns_and_rank.py:921-979`）。

## 3. 實際 structural objective 與限制式

令 `F` 是上述達 structural exact-pattern 門檻之完整 R/A patterns 對應的 mandatory states、`0` 是 root、`x_v∈{0,1}` 表示節點 `v` 是否被選入 `N`。新 symbolic MILP 實作的是：

\[
\min \sum_{v\in U\setminus(F\cup\{0\})}x_v.
\]

也就是最少的非 root、非 full-observed 額外節點；partial 所選代表與純 connector 都計為 hidden/extra。限制為：

> **名詞邊界**：程式欄名 `minimum_hidden_nodes` 是歷史命名；數學上更精確的說法是 **minimum-extra-state count**。它同時包含 partial group 選中的 representative state 與為了 root connectivity 所需的 intermediate connector，不等於「未觀測的真實 clone 數」，也不是腫瘤細胞數。

1. `x_0=1`、每個 `f∈F` 均 `x_f=1`（`hypercube_exact.py:302-325`）。
2. 每個被選的非根節點至少有一個已選 Hamming-1 predecessor：
   \[
   x_v\leq\sum_{j:v_j=1}x_{v-e_j}.
   \]
   因 predecessor 的 popcount 嚴格減一，反覆回溯必到 root，所以此 local constraint 已保證 root connectivity（`hypercube_exact.py:330-345`）。
3. 每個 partial group 至少被一個已選節點覆蓋：
   \[
   \sum_{v\in G(p)}x_v\geq1\quad\Longleftrightarrow\quad N\cap G(p)\neq\varnothing
   \]
   （`hypercube_exact.py:347-354`）。

在不改變 optimum 的前提下，目前會先做 effective-group duplicate／dominance、mandatory-hit、singleton forcing，再只對 active bits 建 predecessor rows；固定最小 hidden objective 後，使用 sparse exact no-good cuts 列舉不同 vertex sets。`k=1..3` 的 18 個 seeded cases 已與獨立 exhaustive oracle 的完整 optimal-set 集合逐案相同；`k=4` 仍完整列出 24 個 optima、digest 前後一致。這降低 constraint 與 nonzero 數，但「V 個 optimal sets 約需 V+1 次 MILP」的尾端仍存在，故 154-task scaling 仍需真實 pilot。

凍結 solver 的邏輯相同，但實作方式不同：先以 `subcube_members()` 明列各 partial group 的 members，再由 `_covers()` 檢查交集非空、`_is_closed()` 檢查 predecessor closure（`003_tree_enumeration_solver.py:29-42,59-73`）。它由 `e=0,1,…` 逐層枚舉 extra node combinations，第一個可行層即最少 hidden 層（`003_tree_enumeration_solver.py:142-157`）。

這更精確地說是 **rooted directed group-Steiner-like minimum-extra-vertex problem on an oriented Boolean hypercube**。選定 `N` 後，每個非 root 節點只選一條 parent edge，因此 edge 數固定為 `|N|-1`；在 full terminals 固定時，最小 extra vertices 與最小 arborescence size 只差常數。

## 4. 為何不能「每種 completion 分開跑，任一成功就通過」

### 4.1 直接結論

- 若只要某個 completion **可行**就停止：錯，因為先遇到的可行解可能不是全域最少 hidden。
- 若把 partial completion 當成 full terminal：錯，因為它會從 objective 中被免計，改變 hidden 定義。
- 唯一可與 production structural group constraints 等價的 explicit 方法是：先保留達門檻的distinct partial patterns，再枚舉**所有 retained structural groups 的 joint completion assignment**、partial-selected vertex 仍計入 extra cost、每個 assignment 求最小解、最後只保留全域最小 objective 並 deduplicate。這在數學上可等價，但複雜度為所有 retained groups completion 數的乘積，沒有實作價值；低於門檻的informative molecules仍只在likelihood層使用。

若第 `i` 個retained structural partial pattern有 `u_i` 個 `X`，naive joint enumeration 最壞為

\[
\prod_i 2^{u_i}=2^{\sum_i u_i},
\]

而不是單一 pattern 的 `2^u`。

### 4.2 `RRA + AAA + RAX` 是「任一成功」會出錯的直接反例

`RAX` 的一個 `X` 在 full 3-cube 中有兩個 completion：

\[
G(RAX)=\{RAR,RAA\}.
\]

兩者都能形成一棵可行樹，但最小 hidden 不同：

| 強迫使用的 completion | 最少 non-full vertices | 結論 |
|---|---:|---|
| `RAR` | 2 | 可行，但不是全域最佳 |
| `RAA` | 1 | 全域最佳 |

因此若先測到 `RAR`、看到「成功」就停止，會錯過 hidden 少一個的 `RAA` 解。

### 4.3 `AX + XA` 證明多個 partial groups 必須聯合求解

在 `k=2` full cube 中：

```text
G(AX) = {AR, AA}
G(XA) = {RA, AA}
```

4 個 joint completion assignments 會產生 5 筆各 world 的 minimum outputs；跨 worlds 比較同一個全域 hidden objective 並去重後，只剩 3 個合法的全域最小 vertex sets：

```text
[RR, AR, RA]
[RR, AR, AA]
[RR, RA, AA]
```

symbolic solver 與獨立 explicit-joint oracle 均得到 `h*=2`、`V=3`。若把兩條 partial patterns 分開做 first-success，會分別輸出 `[RR,AR]` 與 `[RR,RA]`；前者未命中 `XA`，後者未命中 `AX`，兩者都不是原問題的可行解。這個例子同時證明：必須聯合滿足全部 groups、跨 worlds 比較全域 objective、並對重複 `N` 去重。

## 5. 先列 optimal vertex sets，再列 parent-edge trees

### 5.1 全部 optimal vertex sets

新 solver 先求 objective `h*`，接著固定 objective 等於 `h*`，逐次加入 exact binary no-good cut 排除已找到的 `N`，直到 MILP 回報 `INFEASIBLE`；只有此時 `complete=true` 才能宣稱所有 optimal vertex sets 已列完（`hypercube_exact.py:356-388,466-529`）。若達 `max_sets`、time limit 或非 optimal status，只能標 incomplete/abstain。

凍結 solver 在 `capped=false` 時保留第一個可行 hidden 層的全部 `feasible_N`。但若 `C(pool,e)>per_level_budget`、超過 `extra_cap`，會改用 greedy fallback，不能宣稱全枚舉；此外 `tree_cap=32` 是 stored tree 上限，`trees_complete=false` 時展示清單不等於完整候選集（`003_tree_enumeration_solver.py:138-165,188-209`）。

### 5.2 同一 vertex set 的全部 parent choices

對固定的 `N`，每個非根節點 `x` 可從 `N` 內所有 unit predecessor `x\{j}` 選一個 parent。所有節點 parent choices 的 Cartesian product 就是該 `N` 的全部 arborescences；由於每條邊 popcount 增一，天然無環且 root-connected（`003_tree_enumeration_solver.py:211-220`）。樹數為

\[
T(N)=\prod_{x\in N\setminus\{0\}}|Pred_N(x)|.
\]

新 module 的 `parent_choice_count()` 只計數、**沒有實際輸出 edge combinations**（`hypercube_exact.py:532-546`）；實際 edge enumeration 目前由凍結 solver `_parent_choice_trees()` 提供。故工作流必須分清：

1. `V`：distinct optimal vertex sets；
2. `T`：每個 `V` 下的 parent-edge trees；
3. `T` 多不必然表示 read data 可排序其 edges。

## 6. read-pattern likelihood 如何處理 partial read

對某候選 vertex set `N`，令每個節點的混合比例 `π_v≥0` 且 `Σπ_v=1`。對 pattern `p`，只在 fixed positions 計算 matches/mismatches；`X` 不產生 match 或 mismatch。固定錯誤率 `ε` 下：

\[
q_{pv}=P(p\mid v,\epsilon)
=(1-\epsilon)^{m_{pv}}\epsilon^{d_{pv}}.
\]

這個固定 `ε` 公式是 symbolic pilot 與 M0 engineering baseline 的簡化模型。schema 2.0 M2 primary 會保存每個 fixed R/A call 的 Phred BQ。令 `e=10^{-Q/10}` 是「變成任一錯誤鹼基」的總機率；在明示的 symmetric four-base substitution 假設下，再 condition on observed base ∈ {REF, ALT}：

\[
P(\text{match}\mid R/A)=\frac{1-e}{1-2e/3},\qquad
P(\text{R}\leftrightarrow\text{A flip}\mid R/A)=\frac{e/3}{1-2e/3}.
\]

因此不能直接把 `e` 當成特定 R↔A flip 機率；那會把此錯誤放大約三倍。`O/D/S/L/X` 保留在漏斗，但在目前 primary emission 中作 missing marginalization。

read／molecule group count 為 `n_p` 時，候選 `N` 的 profile log-likelihood 是

\[
\ell(N)=\max_{\pi\in\Delta}
\sum_p n_p\log\left(\sum_{v\in N}\pi_vq_{pv}\right).
\]

這裡 `Σ_v π_v q_pv` 就是對產生該 partial read 的 latent mutation state 做 marginalization；不需列出 X 的所有 completion。M2 primary 的 per-BQ emission 與 certified mixture 串接見 `build_m2_patterns_and_rank.py:595-681`；共用 optimizer 的 warm-start、pairwise refinement與 global certificate 見 `hypercube_exact.py:705-839`。固定錯誤率 sensitivity 另外走同一個 certified objective（`hypercube_exact.py:926-973`）。因此 SLSQP 只作 warm start；`scipy.success` 本身不是權威結論。

此 likelihood 成立需明示至少三個假設：

1. read/molecule 在給定 latent state 後條件獨立；
2. `X` 的缺失機制在條件化 observed mask 後不額外依賴 allele/state；若 allele-specific mapping 或 deletion 造成非隨機缺失，需把 missingness/call class 納入 emission；
3. M0 使用同一個 `ε`；M2 schema 2.0 使用 per-call BQ 的 conditional R/A emission，但仍不是完整的 ONT context／strand／indel／modified-base error calibration。BQ grid 與 fixed-error grid只能作敏感度分析。

### 6.1 為何相同 vertex set 的不同 edges 必須同分

上式只依賴 `N`、pattern 與 mixture weights，不含 parent edge `E`。因此

\[
N_1=N_2\Rightarrow \ell(N_1,E_a)=\ell(N_2,E_b).
\]

這是資料可辨識性限制，不是 optimizer 缺陷。snapshot read 共現與由同批 reads 算出的 VAF 能判斷「哪些 mutation states/clone proportions 更符合資料」，但不能區分同一組 states 間不同的 parent assignment。要分 edge，需要真正 edge-informative 的獨立資訊或額外假設；甲基若只改善 state assignment，仍不自動變成 transition observation。

另一個必須明說的限制是：目前為嚴格兩階段選擇。第一階段只保留 minimum-extra objective 等於 `h*` 的 vertex sets，第二階段才在這個集合內比較 read-pattern likelihood。因此「likelihood winner」的精確意義是 **在 minimum-extra candidate set 內最符合 reads**；現行流程不會把 `h*+1` 或更大但仍可行的結構一起排名，所以不可改寫成「所有可行結構中全域最可能的真實拓撲」。

同理，由相同 reads 得到的 VAF 已由 `π`/pattern counts 使用，不能再當獨立分數加一次，否則 double counting；source-bound method contract明示`same_molecule_vaf_added_as_separate_term=false`。VAF 可作可解釋的 posterior summary、外部校準或與 caller AF 的 sensitivity analysis；若直接加入 joint objective，必須先定義其獨立性與 covariance。這個contract欄位是方法宣告，必須與score construction、producer SHA與negative regression一起驗證，不能把布林值本身當成獨立量測證據。

## 7. 實際執行驗證

### 7.1 單元測試

**輸入**：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`

**執行命令**：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest \
research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py -v
```

**輸出**：stdout；結果摘錄記錄於本節，未建立大型衍生檔。

**實際輸出（exit code 0）**：

```text
test_rra_aaa_rax_unique_minimum ... ok
test_rra_aaa_rax_symbolic_matches_global_completion_world_minimum ... ok
test_ax_xa_requires_joint_coverage_and_deduplicated_global_worlds ... ok
test_symbolic_matches_explicit_for_all_patterns_through_k6 ... ok
test_same_vertex_set_is_edge_invariant ... ok
test_seeded_small_k_matches_independent_exhaustive_oracle ... ok
test_complete_chain_enumeration_keeps_all_24_optima ... ok
test_slsqp_matches_em_log_likelihood ... ok
Ran 23 tests in 1.115s
OK
```

驗證：symbolic membership 與 explicit completion 在 `k=1..6` 全 pattern 一致；`RRA+AAA+RAX` 的 symbolic 解等於跨全部 completion worlds 的全域最小解；`AX+XA` 的joint coverage、global objective與dedup均一致；exact reductions 與 independent exhaustive oracle 的完整 optimal-set 集合一致；同 vertex set 順序不改 likelihood；EM/SLSQP log-likelihood相符。加入兩個固定反例與可重跑exhaustive-reduction audit後，歷史整個research test suite為162/162 PASS（wall 12.96s、peak RSS 90,940KB）。2026-07-17以current test SHA獨立重跑本方法targeted suite為23/23 PASS（1.156 s）；後者不替代真實M2 receipts。

### 7.2 新舊 solver 同例比對

**輸入**：上述兩個 solver，`full=[RRA,AAA]`、`partial=[RAX]`、`k=3`。

**執行命令**：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 - <<'PY'
import importlib.util, json, pathlib, sys
repo = pathlib.Path('/big7_disk/liaoyoyo2001/InterSubMod')
sys.path.insert(0, str(repo / 'research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts'))
from hypercube_exact import SymbolicPattern, enumerate_optimal_vertex_sets, parent_choice_count, state_to_pattern
frozen_path = pathlib.Path('/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/source_bundle/files/imported/003_tree_enumeration_solver.py')
spec = importlib.util.spec_from_file_location('frozen_solver', frozen_path)
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)
p = SymbolicPattern.from_string('RAX')
n = enumerate_optimal_vertex_sets(['RRA', 'AAA'], ['RAX'], 3, max_sets=20)
f = frozen.enumerate_min_trees(['RRA', 'AAA'], ['RAX'], 3, tree_cap=0)
print(json.dumps({
    'completions': [state_to_pattern(v, 3) for v in p.enumerate_completions()],
    'new_hidden': n['objective'],
    'new_complete': n['complete'],
    'new_vertex_sets': [[state_to_pattern(v, 3) for v in vs] for vs in n['vertex_sets']],
    'new_parent_counts': [parent_choice_count(vs) for vs in n['vertex_sets']],
    'frozen_hidden': f['n_hidden'],
    'frozen_n_feasible_N': f['n_feasible_N'],
    'frozen_n_trees': f['n_trees'],
    'frozen_capped': f['capped'],
    'frozen_trees_complete': f['trees_complete'],
    'frozen_trees': f['trees'],
}, indent=2, sort_keys=True))
PY
```

**輸出**：stdout；本文件記錄所有判讀所需的實際值。

**實際輸出（exit code 0）**：

```text
RAX: n_X=1, conceptual completions=[RAR, RAA]
new symbolic MILP: objective_hidden=1, complete=true,
  optimal vertex sets=[[RRR,RRA,RAA,AAA]], parent_choice_count=1
frozen canonical: n_hidden=1, n_feasible_N=1, n_trees=1,
  capped=false, trees_complete=true,
  edges=ROOT→RRA→H_RAA→AAA
```

Step→Verify：

1. 解析 `RAX` → 驗證：一個 `X`、full cube 中恰有 `RAR/RAA` 兩個 members。
2. 求最少 hidden → 驗證：兩個 solver 均為 `1`，且 optimal vertex set 唯一。
3. 列 parent choices → 驗證：只有一種 edge tree，無 recurrence。
4. 強迫 completion 分開求解 → 驗證：`RAR` 最少需 2 個 non-full vertices；`RAA` 只需 1 個，證明 first-success 不保證 global optimum。

第 4 點使用下列完整小型 brute-force 命令（只枚舉 `k=3` 的 8 個 states）：

```bash
python3 - <<'PY'
import itertools
k, root, full = 3, 0, {4, 7}  # RRA, AAA
label = lambda v: ''.join('A' if v & (1 << i) else 'R' for i in range(k))
def closed(N):
    return all(v == 0 or any((v ^ (1 << b)) in N for b in range(k) if v & (1 << b)) for v in N)
def best(forced):
    mandatory = {root, *full, forced}
    optional = set(range(1 << k)) - mandatory
    for r in range(len(optional) + 1):
        feasible = []
        for extra in itertools.combinations(optional, r):
            N = mandatory | set(extra)
            if closed(N):
                feasible.append(N)
        if feasible:
            cost = len(feasible[0] - ({root} | full))
            return cost, [[label(v) for v in sorted(N)] for N in feasible if len(N - ({root} | full)) == cost]
for completion in (2, 6):  # RAR, RAA
    cost, sets = best(completion)
    print(f'forced_completion={label(completion)} minimum_nonfull_vertices={cost} optimal_vertex_sets={sets}')
PY
```

實際輸出（exit code 0）：

```text
forced_completion=RAR minimum_nonfull_vertices=2 optimal_vertex_sets=[...]
forced_completion=RAA minimum_nonfull_vertices=1 optimal_vertex_sets=[['RRR', 'RRA', 'RAA', 'AAA']]
```

### 7.3 既有更廣 pilot 證據（僅作工程驗證）

`results/pilot_v3_identity_bound/pilot_receipt.json` 明示 scope 為 `PILOT_NOT_FINAL_VALIDATION`，並把執行時 runner／solver identity 綁入 receipt；其可重現結果為：

- symbolic exhaustive：9,840 patterns、2,015,538 state checks、0 mismatch（k=1..8）；
- frozen/MILP cross-check：66 vertex-set checks、0 mismatch，另有 9 個 frozen-capped cases 跳過，不能算 pass；
- k9–12：12/12 pilot cases pass，最大 4,096 variables、最大 `0.191109856` s；
- likelihood：同 vertex set/order delta `1.42×10^-14`，partial missing-dimension delta `0`。

這些數字支持實作一致性，但**不是全 7 樣本、全區域的最終生物驗證**。

### 7.4 Presolve／no-good deterministic exhaustive audit

**輸入**：current solver SHA256=`9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`；audit script SHA256=`7d3a1c465b9bec7951670ec2aa8787768eb282b450afcd52b7ff2a4b7123674a`。

**執行命令**：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
/usr/bin/time -v nice -n 10 python3 \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/audit_hypercube_reductions_exhaustive.py \
  --json-output research/20260716_read_linked_hypercube_exact_likelihood_validation/results/hypercube_reduction_exhaustive_v1/receipt.json
```

**輸出**：`results/hypercube_reduction_exhaustive_v1/receipt.json`（SHA256=`a60037e7cf527d82abfc1b676ba7405a7b8a66625df73491c57ad8e7bd7b88b8`）與checksum sidecar；`sha256sum -c` PASS。

| 獨立檢查 | Cases / predicates | Mismatch |
|---|---:|---:|
| k=1..3、全universe masks、≤3 partial-group multiset、≤2 extra mandatory 的presolve cases | 61,340 cases；1,979,356 selected-set predicates | 0 |
| Fixed-cardinality sparse no-good | 23,909 set pairs | 0 |
| General dense no-good | 21,844 set pairs | 0 |

Reduction事件為duplicate 19,918、dominance 1,554、mandatory/forced hit 112,111、singleton forced 4,848；這些是測試事件，不是樣本數。原始group membership由獨立的R/A/X字元比較重建，而非呼叫`SymbolicPattern.compatible`；reducer本身則直接稽核current private `_reduce_partial_groups`。Receipt `all_pass=true`；wall=6.19s、peak RSS=59,368KB。這個audit證明的是上述有限exhaustive scope內的algebraic equivalence；若將來加入capacity、per-read assignment或edge weights，必須重新證明。

## 8. 複雜度與 k>12 的較高效做法

### 8.1 現況複雜度

| 方法 | partial 表示成本 | structural search 主成本 | 限制 |
|---|---:|---:|---|
| naive joint completions | 每 read `2^u`，joint 為 `2^{Σu_i}` | 每個 joint assignment 再解一次 tree | 不應採用 |
| frozen canonical | 每 group materialize `2^{|X∩U|}` members | `Σ_e C(|pool|,e)`；parent trees 為 `Π_x |Pred_N(x)|` | `extra_cap/per_level_budget/tree_cap` 可能截斷 |
| 新 symbolic MILP | 每 pattern O(1) masks | 目前仍有 `2^d` binary vertices，group rows 建構約 O(`G·2^d`) | 避免 completion product，但沒有消除 hypercube 指數性 |

其中 `d=|U|≤k`。所以 symbolic group 是必要改善，但不能宣稱已把 NP-hard/指數性問題消除；列出「全部」optimal sets/edge trees本身也可能有指數級輸出量。

### 8.2 建議的 production 路線（依優先順序）

1. **pattern compression**：相同 `(fixed_mask,alt_mask)` 的 reads 合併 count；structural coverage 每個 unique group 只留一條 constraint，likelihood 保留總 count。
2. **group dominance**：若 `G_a⊆G_b`，覆蓋 `G_a` 已必然覆蓋 `G_b`，structural 層可移除冗餘的 superset constraint；scoring 層仍保留兩者 counts。
3. **implicit/lazy hypercube MILP**：不預先建立全部 `2^d` vertices，先用 full states、partial representatives 與必要 ancestor closure 起始，透過 column generation / branch-and-cut 加入能改善 objective 的 vertices。
4. **BDD/ZDD 或 bitset family 表示**：以 decision diagram 表示 subcubes、ancestor closures 與 no-good families，避免反覆 materialize 大量 members/sets。
5. **read-linkage decomposition**：只在有 molecule bridge 的 locus component 內解；若要再切 component，必須以 zero-bridge 或帶 boundary-state 的 exact decomposition 證明等價。任意 densest-8 截取只能標 local/partial，不能稱完整區域最佳樹。
6. **certified branch-and-bound**：用 closure/coverage lower bounds、dominance memoization、warm-start incumbent；只有 `MIP gap=0` 才稱 objective optimal。
7. **輸出敏感的 tie enumeration**：先求 `h*`，再枚舉 distinct vertex sets；達 time/max-set limit 時回報 lower bound、incumbent、已找到 sets 與 `ENUMERATION_INCOMPLETE`，不可挑第一棵當唯一答案。
8. **vertex-first scoring**：read likelihood 只跑 distinct `N`；不要對相同 `N` 的每個 edge tree重複最佳化。最後才展開 edges，並把同 `N` 下的 edges標成 read-likelihood tie。

對 `k>12`，建議的誠實決策規則是：

- exact solver `OPTIMAL` 且 `gap=0`：可報 objective 已證明；
- optimal-set enumeration 到 infeasible：可報 all vertex sets complete；
- 只取得 incumbent/time limit：可報候選與 bound，不可報穩定最佳；
- 經 exact read-linkage decomposition：可報 component-local exact；
- 任意 cap/取樣：必須標 `PARTIAL/LOCAL_ONLY`，不得併入全區域唯一拓撲比例。

## 9. 最終判讀與方法契約

使用者原理解「`u` 個 X 有 `2^u` 種可能」在 full hypercube 下正確；「合理的都輸出，再由後續資料排序」方向也合理，但需修正為下列可驗證流程：

```text
R/A/X molecules
  → exact R/A/X pattern counts（所有informative molecules保留供likelihood）
  → count≥structural_exact_pattern_minread的distinct symbolic groups
  → 全部retained structural group-cover constraints（不展開joint completions）
  → 求全域 minimum hidden h*
  → 列所有 optimal vertex sets N（需 completeness 證明）
  → 對 distinct N 做 partial/error-aware read likelihood
  → 排序可被 reads 區分的 N
  → 每個 N 才展開 parent-edge trees
  → 同 N 的不同 edges維持 tied/不可辨
  → VAF、甲基與其他資料依獨立性作校準或正式 joint model
```

最關鍵的三條防誤讀規則：

1. **completion 可行 ≠ 全域最小**；不可 first-success。
2. **vertex-set likelihood ≠ edge likelihood**；同 `N` 不可硬排序 parent edges。
3. **同 reads 的 VAF ≠ 獨立證據**；不可把已在 pattern likelihood 使用的資訊再加權一次。
4. **minimum-extra state ≠ hidden clone**；欄名 `minimum_hidden_nodes` 不可直接作為真 clone 數。
5. **winner 是條件式 winner**；它只在 `h=h*` 的全部 minimum-extra vertex sets 內排名，不涵蓋較大可行結構。

---

驗證狀態：本文件的程式與小例子稽核通過；全 7 樣本 quantitative census 由本研究 cycle 其他工作流負責，未在本文件重複宣稱。
