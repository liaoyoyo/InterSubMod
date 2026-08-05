<!--
建立時間: 2026-07-16 08:20 +08:00
目標: 稽核現有 hypercube group-Steiner MILP 在 k>12 的複雜度，設計可保留 global optimality certificate 的高效路徑
處理範圍: 探索性效能 pilot；k=13-16 每 k 兩個合成 seeded cases；不是全樣本、全區域或生物結論
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/benchmark_k_gt12_exact.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot_v3_identity_bound/pilot_receipt.json
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_方法學原始文獻與主張邊界_01.md
-->

# k>12 exact 效率與 certified 策略稽核

> **2026-07-16 10:20 更新／閱讀順序**：本文第3–4節的複雜度盤點是縮減實作前的audit baseline。其後已落地duplicate/dominance、mandatory-hit、singleton forcing、active-bit predecessor、`O(2^u_eff)` sparse compatibility generation、fixed-hidden sparse no-good與`h=0` early-complete；最新驗證見`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_Hypercube_exact_solver縮減與列舉效能稽核_01.md`。目前仍保留`2^m` vertex variables／predecessor rows、fixed-objective row，以及每個distinct optimum約再解一次MILP的主要瓶頸。current solver SHA256=`9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`（相對縮減audit版另有module wording-only修正）。

> **PARTIAL / EXPLORATORY PILOT**：本文件服務 G3/G4/G5，但只回答演算法合約與小型合成效能；不能替代 7 datasets × chr1–22 的 Task-B 全量數據驗證，也不能把數學最優稱為真實 clone ancestry。

用 SCQA：現有方法在固定 active loci 上是可認證的 exact MILP；困難不是名義上的 `k`，而是 `active-k=m` 造成 `2^m` 個 vertex variables，加上 partial-group rows 與逐解 no-good rows 的指數展開。**近期最安全路線是先做 exact active-bit/group dominance reduction，再依 `(m, group數, terminal數, BDD width)` 路由至 MILP、ordinary-terminal DP 或 exact BDD/branch-and-price；任何逾時、非零 gap、width/solution cap 都必須 abstain。**

## 1. 稽核對象與 frozen 證據

| 輸入 | 路徑 | SHA-256 / 狀態 |
|---|---|---|
| exact solver（本文原始複雜度 audit-time snapshot） | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py` | `db85c94f5cd1fda5bb1d4eca9ca64d5df6680829b1b005ac02be9d3c835500bc`；只用來解讀第3節「縮減前 baseline」 |
| exact solver（current） | 同上 | `9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95` |
| current-solver identity-bound pilot | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot_v3_identity_bound/pilot_receipt.json` | `5cab8c37f31908d13b0c921619264598a4e6e3ef29cbcf723d1a6c421b6297c8`；k9–12 12/12 PASS，最大 4,096 variables，最大 `0.191109856` s |
| 新 benchmark | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/benchmark_k_gt12_exact.py` | `bdfcfda2f368f3b8299ecfe72642c3bdf513f749296963d633e26a88c9087440`（正式 receipt 會再次凍結執行時 SHA） |
| 新 unit test | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_benchmark_k_gt12_exact.py` | `b55db6d989e2f71434a0b85e56e775add1d3236857abb0b2394127f69a271958` |
| 原始理論來源 | Mahapatra et al., arXiv:2110.02830 | ordinary-terminal directed-hypercube MSA-DH；**不是**本研究 group-terminal extension 的現成證明 |

版本環境：Python 所見 `SciPy 1.13.1`、`NumPy 1.23.5`；MILP 由 SciPy `milp`／HiGHS 執行，`mip_rel_gap=0`。

## 2. 現有 MILP 到底解什麼

令 active bit mask 為 `M`、`m=|M|`，candidate vertices 為 `V={v:v⊆M}`。root 與 full-read states 是 mandatory set `F`；partial pattern `p_j` 對應相容子立方體 `G_j`。

### 2.1 變數、目標與限制

對每個 `v∈V` 建 binary variable `x_v`：

\[
\min \sum_{v\notin F} x_v
\]

其中 root/full states 由 variable bounds 固定為 1。每個非 root state 有 predecessor constraint：

\[
x_v \le \sum_{u\in Pred(v)}x_u.
\]

每個 partial read group 只有一個「至少一種 completion 被選」的聯合限制：

\[
\sum_{v\in G_j}x_v\ge1.
\]

因此不是逐一把 `2^u` 個 completion 當成 `2^u` 份 observation，也不是「每個 completion 各跑一棵、任何一棵成功就算」。所有 full/partial constraints 在**同一個**最佳化問題中一起滿足。

### 2.2 為何 predecessor constraints 足以表示 rooted tree

每條 predecessor edge 都讓 popcount 減 1；圖是 DAG。只要每個 selected non-root vertex 至少有一個 selected predecessor，反覆往前必到 root。再替每個 non-root vertex 任選一個 selected predecessor，就得到一棵 rooted arborescence。因此：

- MILP 先決定最小 vertex set `V*`；
- 同一個 `V*` 可能有多種 parent-edge assignments；
- read-state likelihood 觀察 state，不直接觀察 mutation edge，所以同 `V*` 不同 edge 應維持並列。

### 2.3 no-good enumeration 合約

第一輪取得 `h*` 後，程式固定 hidden objective=`h*`，每找到一個 binary vertex set 就加一條 exact no-good row。只有在排除已知解後，下一輪 MILP 回傳 `INFEASIBLE`，才能聲稱「所有最優 vertex sets 已完整枚舉」。

- `MAX_SETS_REACHED`：只證明至少找到 `max_sets` 個，不證明總數。
- `LIMIT_REACHED` / hard timeout：可以保存 incumbent 作診斷，但不可稱最優或完整。
- 現行 `time_limit_seconds` 是**每次 MILP solve**，不是整段 enumeration 的 global budget；新 benchmark 額外用 subprocess hard timeout 與 global deadline 封住總時間。

## 3. 靜態複雜度：縮減前 audit baseline 與 current 差異

下表的前五列與「dense no-good baseline」是本文件建立時、exact reductions 尚未落地的複雜度盤點；不可再把 dense no-good 那一列當成 current implementation。current solver 仍有 `2^m` variables、predecessor rows與fixed-`h*` row，但 partial compatibility 改為直接生成相容 indices，且固定 objective 下的 no-good row只列舊解的 extra vertices。

令 `g` 為 partial groups 數、`u_j` 為第 j 個 pattern 在 active bits 中的 X 數、`e` 為已排除的 optimum vertex sets 數：

| 結構 | 數量 / 非零係數 |
|---|---:|
| vertex variables | `n=2^m` |
| predecessor constraints | `2^m−1` |
| predecessor matrix nonzeros | `(2^m−1)+m·2^(m−1)` |
| partial-group rows | `g` rows、約 `Σ_j 2^(u_j)` nonzeros |
| fixed-`h*` row | 最多 `2^m−|F|` nonzeros |
| `e` 條 dense no-good rows（**縮減前 baseline**） | `e·2^m` nonzeros |
| `e` 條 fixed-hidden sparse no-good rows（**current**） | 通常 `e·h*` nonzeros；每條約 `h*`，其中 `h*` 是minimum extra vertices。只有不具 fixed-cardinality equality 的general fallback才回到dense exact exclusion |

原始 audit 的關鍵發現是：`SymbolicPattern` 的 mask 避免在 Python 先產生 completion Cartesian product，但縮減前 `_build_problem` 仍掃描所有 vertices。current solver 已改為每次 MILP build、每個reduced group直接列 `v=alt_mask OR s` 的相容active indices；row本身仍有至多`2^u_eff`個coefficients，且每次all-optimal no-good重建MILP時會重新產生。故改善了row建構，**尚未讓整體 MILP matrix 脫離 `2^m` variables／predecessor rows**。

### 3.1 名義 k 與 active-k 必須分開報

`observed_alt` 模式只保留至少在一個 full/partial pattern 明示 A 的 loci。若 total `k=16`、active `m=8`，variables 仍是 `2^8=256`，不是 65,536。反之只要 active `m=16`，即使 read 很少，variables 就是 65,536。

## 4. 哪些改進仍保證 global optimum

### 4.1 可立即落地的 exact preprocessing

| 方法 | exactness 條件 | 效果 | 結論 |
|---|---|---|---|
| **active-bit compression** | 固定採 `observed-alt` structural contract；沒有任何 full/partial pattern 在該 bit 明示 A，投影後保留 fixed-R／X 相容性 | `2^k→2^m` | 在 **observed-alt model 內完全 exact**；相對 full-universe 只主張保留 optimum objective／至少存在一個投影最優解，不主張保留所有 V-set identities |
| duplicate group compression | structural stage 中相同 `(fixed_mask,alt_mask)` 合併；read count 留在 likelihood | 降 `g` | **完全 exact** |
| mandatory-hit removal | 若 mandatory vertex 已落在 `G_j`，刪除該 group constraint | 降 rows | **完全 exact** |
| group dominance | 若 `G_a⊆G_b`，滿足較嚴格 `G_a` 已保證 `G_b`；刪 `G_b` | 降 rows/NNZ | **完全 exact** |
| sparse compatibility generation | 直接列 `v=alt_mask OR s`、`s⊆(free_mask∩M)`，不再掃描全部 `2^m` vertices | row 建構由 `O(2^m)` 降為 `O(2^u)` | **完全 exact**；row 本身仍有 `2^u` coefficients |
| relevant-ancestor pruning | vertex 必須是某 full terminal 或某可選 group member 的 ancestor；其餘不可能出現在 inclusion-minimal optimum | 降 candidate vertices | **exact**，但需保留所有 group representatives 的 ancestor union |

active-bit 的保守 exactness 邊界：在程式明定的 `observed-alt` universe 中，它是模型定義的一部分，所以 MILP 對該模型完全 exact。若比較更寬的 `all-loci` universe，把「從未明示 A」的 bit 投影回 0 仍可滿足 fixed-R/X groups，且不增加 vertex 數；因此可保留 optimum objective，並保證至少有一個投影最優解。**但 full-universe 中同 objective 的替代 V-set identity 是否全部被保留，未由現有程式或文獻證明，不可主張。**正式方法需把 observed-alt 稱為 structural prior，並在代表性 k 可承受區域做 `all-loci` sensitivity：比較 `h*`、最優 V-set 數／投影 digest 與 downstream likelihood 排名。

group dominance 的 mask 判斷亦可不展開 completions：`G_a⊆G_b` 當且僅當 b 固定的每個 bit 在 a 也固定為同值，即 `(fixed_b & ~fixed_a)==0` 且 `((alt_a xor alt_b)&fixed_b)==0`。這讓 duplicate/dominance preprocessing 對 `2^u` 無依賴。

### 4.2 exact separator / articulation decomposition

可在「pruned directed candidate DAG」找 directed dominator 或小 separator `S`：若所有通往一組 terminals/group members 的 root path 必經 `S`，可在固定 boundary state 後拆子問題。

要保留 global optimum，必須：

1. crossing group 可落在 separator 兩側時，**枚舉所有 boundary/group assignment**；
2. 每個子問題都回傳 gap=0 或 exact DP certificate；
3. 合併時扣除共享 separator vertices，並比較所有 boundary assignments。

只用 undirected articulation、把 crossing group 指派給 read 數最多的一側，或各 block 獨立取一棵再拼接，都是 local heuristic，沒有 global certificate。

### 4.3 ordinary-terminal DP 與 group-terminal 邊界

Mahapatra et al. 對普通 terminal set `R` 給出 LCA（bitwise AND）subset DP：時間約 `O(3^|R|·|R|·k + |R|²k²)`、空間約 `O(2^|R| log(|R|k))`。它對 `k` 是 polynomial，所以當 full terminals 少、沒有 partial groups 時，能避開 explicit `2^k` cube，且是**全域 exact**。

但 partial group 是一組可能 terminals。把每個 group 任選一個 completion 後直接套 DP，若沒有枚舉所有 choices，就不是 exact。可認證的 group extension 有兩條研究路線：

- 用 BDD 表示每個 subcube，對所有 representative choices 做共享狀態的 min-plus DP；
- branch over BDD decisions，使用 ordinary-terminal DP 當 leaf solver與 admissible lower bound，直到所有 branch 被證明或剪除。

這兩條只有在**BDD 不截 width、branch 完整、lower bound 可證**時保留 global optimum；目前尚未實作／證明，不能寫成既有能力。

### 4.4 BDD/ZDD、branch-and-cut / branch-and-price

- 每個 R/A/X subcube 可用 O(k) 決策路徑表示；BDD 可共享大量相同 suffix。
- optimum vertex-set family 可用 ZDD 壓縮，避免每找到一解就新增一條新row並重解MILP。current fixed-hidden sparse no-good已把單列nonzeros由dense `2^m`降到約`h*`，但沒有消除每個distinct optimum約一次solve/build的tail。
- dynamic-column MILP 可只建立目前有用的 vertices／compatibility coefficients；但要稱 exact，必須有 pricing oracle 證明沒有任何 omitted vertex/path 能改善 objective。
- lazy connectivity cuts 也必須把所有 violated cuts 精確分離；只檢查 incumbent 周邊不是證明。

因此，**BDD width cap、beam search、未完成 pricing、只留 top-N columns 都是 local-only**。Exact BDD/ZDD 或 branch-and-price 的優點是「常見 case 很小、最壞仍可 abstain」，不是消除 NP-hard worst case。

### 4.5 `h*` 層的 exact enumeration

建議把工作拆成兩張 certificate：

1. `objective certificate`：找到 `h*` 且 primal=dual、gap=0；
2. `ambiguity certificate`：固定 `h*` 後完整列舉，最後由 infeasibility 或 exact ZDD exhaustion 證明 complete。

可改進項目：persistent solver/warm start、layer-wise DFS + admissible lower bound、instance automorphism 的 orbit enumeration、exact ZDD family。若只跑 `h*` 到 `h*+q`，要明示 parsimony layer；若只找到前 `N` 組，應輸出 `n_optimal_sets ≥ N`，不可寫 `=N`。

## 5. 建議的 certified solver router

1. **Normalize**：壓縮 duplicate patterns、保留 counts 供 likelihood。  
   → 驗證：pattern count mass 在 likelihood 前後守恆。
2. **Exact reduce**：active bits、mandatory-hit、group dominance、relevant ancestors。  
   → 驗證：observed-alt contract 內比對縮減前後 `h*` 與完整 vertex-set digests；另以小 k `all-loci` sensitivity 確認 `h*` 與 projected digests，不能預設 full-universe identities 全同。
3. **Route**：
   - `m≤12`：explicit MILP；
   - no groups 且 `|R|` 小：Mahapatra terminal-subset DP；
   - group BDD width 小：exact BDD/branch-DP；
   - 其餘：branch-and-price/cut，需 dual/pricing certificate；
   - 時間或記憶超限：`ABSTAIN_RESOURCE_LIMIT`。
4. **Enumerate h***：完整則報 exact count；cap 則只報 lower bound。  
   → 驗證：最後 INFEASIBLE/ZDD exhausted 才 `complete=true`。
5. **Likelihood ranking**：每筆 read pattern 一次、partial 位點邊際化；同 V 不同 edge 維持 tie。  
   → 驗證：同 vertex set 的 edge permutation score delta≈0；不可再重複加 derived VAF。

## 6. 明確的 global / local 邊界

| 輸出狀態 | 可以說 | 不可以說 |
|---|---|---|
| MILP/DP gap=0，enumeration complete | 固定模型下 `h*` 與全部最優 V sets 已認證 | 真實生物樹／唯一 parent edge |
| gap=0，enumeration cap | `h*` 已認證；至少有 N 個最優 V sets | 最優 V sets 總數=N |
| time limit 且 gap>0/unknown | 有 incumbent 與 bound；本區域 abstain | incumbent 是最可能／最小樹 |
| BDD width cap、beam/top-N、heuristic partition | local candidate / engineering sensitivity | globally optimal / complete |
| VAF/read likelihood 分出 V sets | 哪些 **state sets** 較符合同一批 reads | 同 V 不同 edges 的祖先方向被 reads 證明 |

## 7. 輕量 benchmark 設計與資源 gate

每個 `k=13..16` 建兩類 deterministic seeded cases，另把 k=13 sparse case 以完全相同 patterns 重跑一次 `all-loci` sensitivity：

- `sparse_active`：total k>12、active-k=8，驗證 exact compression 的變數數固定 256；
- `dense_active`：active-k=k，直接量測 8,192／16,384／32,768／65,536 variables。
- `sparse_all_loci_sensitivity`（只在 k=13）：相同 sparse input 從 256 擴至 8,192 variables；比較 `h*`，並把 all-loci V sets 投影回 observed-alt mask。若任一枚舉有 cap，identity 比較只報 lower bound。

每個 solve/enumeration 放在 disposable subprocess，`OPENBLAS/OMP/MKL/NUMEXPR=1`、child `nice +10`；stage hard timeout=12 s、HiGHS limit=10 s、enumeration cap=8、全 run wall cap=240 s。共 9 cases；逾時輸出 `HARD_TIMEOUT`／`ABSTAIN`，不保留為 certified optimum。

### 7.1 已執行測試

輸入：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/benchmark_k_gt12_exact.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_benchmark_k_gt12_exact.py`

命令：

```bash
python3 -m py_compile \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/benchmark_k_gt12_exact.py \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_benchmark_k_gt12_exact.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_benchmark_k_gt12_exact.py -v
```

實際輸出：`4 tests / 4 PASS / 0 FAIL`，wall `0.001 s`。

subprocess/hard-timeout smoke v2 輸入為 k=3 三個合成 cases（含同 input 的 observed-alt/all-loci sensitivity），輸出 `/tmp/intersubmod_kgt12_smoke_v2/receipt.json`（暫存、不作研究證據）。實際為 3/3 solve `OPTIMAL`、3/3 enumeration 有 `INFEASIBLE` complete certificate、3/3 row contract PASS、wall 2.779 s；sensitivity objective 相同、projected digest 1/1 相交。

### 7.2 k=13–16 正式 pilot 狀態

**截至 2026-07-16 08:37 +08:00：RESOURCE-GATED，正式 pilot 未啟動。** 主機只有 48 logical CPUs，仍有 46 個相關使用者 workers，合計約 3,162% CPU；load average 61.11 / 58.55 / 151.33，1-minute load 已高於可用 CPUs。依任務要求不與大型工作競爭，本次以 `RESOURCE-GATED` 結案；`results/k_gt12_efficiency_pilot/` 經檢查未建立，沒有把 smoke 或預測值冒充正式數據。待 worker 結束且資源 gate 開放後，由主 agent 執行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
nice -n 10 ionice -c2 -n7 \
python3 \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/benchmark_k_gt12_exact.py \
  --output-dir \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/results/k_gt12_efficiency_pilot \
  --k-min 13 --k-max 16 --seed 20260716 \
  --stage-timeout-seconds 12 --solver-time-limit-seconds 10 \
  --enumeration-time-limit-seconds 10 --max-sets 8 \
  --total-wall-limit-seconds 240
```

預期輸出：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/k_gt12_efficiency_pilot/cases.jsonl`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/k_gt12_efficiency_pilot/receipt.json`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/k_gt12_efficiency_pilot/receipt.sha256`

正式 receipt 完成後，本節必補上每 case 的 total-k、active-k、variables、constraints、objective/bound/gap/status、solve runtime、枚舉數／cap／runtime與 timeout/abstain。**目前這些 k=13–16 實測數據不存在；唯一可引用的量化結果仍是既有 k9–12 pilot 與本文件明示為暫存的 k=3 smoke。**

## 8. 結論

1. 現有 MILP 的數學合約是合理的：partial read 是一個 group constraint，全部 reads 聯合求最小 hidden vertex set；不是 completion-by-completion 的「任一成功」。
2. 目前 exactness 的主要限制是 explicit `2^active-k` variables、`Σ2^u_eff` compatibility coefficients，以及每個 optimum set約需再build／solve一次MILP；current fixed-hidden sparse no-good每列約`h*` nonzeros，不再是縮減前的`2^m` dense row。
3. `k>12` 不應一律切掉；應以 active-k、group dominance、terminal/group count 與可認證 solver 狀態判斷。total k 大但 active-k 小可 exact；active-k 大而無 certificate 則 abstain。
4. 最優 objective、完整 V-set enumeration、edge/topology 可辨識性是三個不同層次，必須分開報。
5. 下一個可落地工程優先序：**exact preprocessing → router → global-budget enumeration → BDD/ZDD 或 certified branch-and-price prototype**。在新方法完成小 k 等價驗證前，不應替換現有 canonical solver。

## 參考來源

- Sugyani Mahapatra, Manikandan Narayanan, N. S. Narayanaswamy. *Parameterized Algorithms for the Steiner Arborescence Problem on a Hypercube*. arXiv:2110.02830。ordinary-terminal exact DP 為 `O(3^|R|)` 級、penalty-parameter exact algorithm 為 `O~(36^q)`；本研究的 group-terminal 延伸不在其直接證明範圍。<https://arxiv.org/abs/2110.02830>
