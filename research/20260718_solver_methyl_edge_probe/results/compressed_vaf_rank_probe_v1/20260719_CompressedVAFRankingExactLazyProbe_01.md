<!--
建立時間：2026-07-19
目標：驗證 recurrence-free perfect family 是否能以 DP traceback + certified upper bound，在不先列完整 family 的情況下，重現目前 M2 primary BQ-aware likelihood 的最佳候選與完整 tie class。
處理範圍：Task A exploratory pilot；synthetic m=1..5；hard 09b 因資源互斥未執行。
關聯檔案：receipt.json、receipt.json.sha256、../../scripts/compressed_vaf_rank_probe.py、../../tests/test_compressed_vaf_rank_probe.py
-->

# Compressed VAF/read-AF exact lazy ranking probe

> **PARTIAL / 非 production 證據。** m=1…5 exact oracle 全通過；正式
> hard25 timing 進行期間，`09b` hard case 未執行。因此本報告證明小範圍正確性與
> fail-closed 行為，不宣稱全資料或 production 整體 wall-time 已加速。2026-07-19
> review後，v1 receipt僅保留為當時source-bound小型證據；後續hard execution必須
> 寫入全新`compressed_vaf_rank_probe_v2`，不得覆寫v1。
>
> **P1 supersession：**v1把普通float relaxation剪枝後的結果稱為exact，這個
> machine-certificate claim已撤回。權威小型結果改見
> `compressed_vaf_rank_probe_v2/p1_hardened_r2/`：float-pruned只能是diagnostic，
> authoritative current-endpoint結果必須全候選實評。

## 結論先行

可以在 recurrence-free perfect-family gate 內，以 DP traceback DAG 表示所有
最簡 vertex sets，再用安全 likelihood 上界逐 branch 剪枝。只要
`search_complete=true` 且 `tie_class_complete=true`，即可在**不先 materialize
整個 family**的情況下，回傳與目前 M2 primary BQ-aware profile mixture
likelihood 相同的最佳分數、最佳 vertex-set IDs 與完整 tie class。

小型驗證中，m=1…5 的 score、winner IDs、tie count 對獨立 parent-vector
exhaustive oracle 均為零 mismatch。m=5 的 synthetic case 有120個 logical
candidates，只評分24個、以 certified bound 剪除96個，評分量縮為原本的20%
（5× reduction）。這不是正式 wall-time speedup；hard case尚未量測。

## 同 endpoint 與不同 endpoint

本 probe 保持不變的是：

\[
\max_{\pi\in\Delta}
\sum_i n_i\log\left(\sum_{v\in V}\pi_v Q_{iv}\right)
\]

以及 current ranker 的
`best_log_likelihood - candidate_log_likelihood <= tie_tolerance` tie 規則。

本 probe **不等同**完整 M2 release endpoint：fixed-error sensitivity、
bootstrap、responsibilities、all-candidate TSV、完整 release receipt 都未在此
實作。歷史 `read_af_tree_ordering_multisample.py` 的 ancestor-AF edge heuristic
亦不是上述 likelihood，不能混稱同一分數。

## 演算法與正確性鏈

1. **DP count / traceback DAG**
   - mutation forest以subset mask表示。
   - `forest_count(mask)`與`tree_count(mask)`只儲存子問題數量，不先列出所有解。
   - top-level branches彼此不重疊，且各 branch count總和必須等於完整
     `forest_count(full_mask)`；否則直接丟錯，不允許排名。

2. **Lazy enumeration**
   - 依 branch upper bound由高至低搜尋。
   - 只有需要評分的 branch才逐一 traceback成 vertex set。
   - recurrence-required family在建構期明確拒絕，不能拿此 perfect solver硬算。

3. **安全 branch upper bound**
   - 對 branch內所有可能出現的 states取聯集 \(U\)。
   - 任一 completion \(V_b\) 都有
     \(\operatorname{conv}(V_b)\subseteq\operatorname{conv}(U)\)。
   - 因此在 \(U\) 上擬合相同 mixture後，加上 optimizer KKT gap，得到
     \[
     UB_b \ge \max_{V_b\in b} LL(V_b).
     \]
   - 另有較鬆但仍安全的 row-wise bound：
     \[
     UB_b^{row}
       =\sum_i n_i\log\left(\max_{v\in U}Q_{iv}\right).
     \]

4. **只在嚴格安全條件剪枝**
   - incumbent是已評分可行候選的 likelihood lower bound。
   - 僅當
     \[
     UB_b < incumbent-\text{tie\_tolerance}-\epsilon
     \]
     才能剪掉整個 branch。
   - 所以被剪候選不可能超過 incumbent，也不可能進入其 tie class；之後若找到
     更高 incumbent，這個結論只會更強。

5. **Fail closed**
   - deadline、candidate cap、tie-output cap、非 finite bound、optimizer未通過
     current KKT certificate或count不守恆，均不得發布 authoritative winner。
   - 只有 `search_complete && tie_class_complete` 才回
     `best_log_likelihood`、winner IDs與tie count。
   - 未完成時僅可保留 `diagnostic_incumbent.authoritative=false`。

> **實作界線：**目前prototype只在traceback DAG的**最上層 branches**各算一次
> upper bound，並在該branch內lazy列舉；尚未對每個recursive DP子狀態重算bound，
> 因此是top-level branch-and-bound prototype，不是完整recursive B&B。這不影響
> 已完成case的正確性，但可能讓大型case剪枝不足。

## 小型 exact oracle 結果

| effective m | logical family | 實際評分 | certified prune | tie數 | score差 | winner/tie |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 1 | 1 | 0 | 1 | 0 | PASS |
| 2 | 2 | 1 | 1 | 1 | 0 | PASS |
| 3 | 6 | 2 | 4 | 2 | 0 | PASS |
| 4 | 24 | 6 | 18 | 6 | 0 | PASS |
| 5 | 120 | 24 | 96 | 24 | 0 | PASS |

另有對稱二鏈 control：2個候選皆同分，完整2-way tie正確保留。單元測試亦覆蓋
candidate cap、tie-output cap、recurrence-required rejection與rowwise bound。

### 測試命令

輸入：

- `InterSubMod/research/20260718_solver_methyl_edge_probe/scripts/compressed_vaf_rank_probe.py`
- `InterSubMod/research/20260718_solver_methyl_edge_probe/tests/test_compressed_vaf_rank_probe.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py`

命令：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 -m unittest -v \
research/20260718_solver_methyl_edge_probe/tests/test_compressed_vaf_rank_probe.py
```

實際結果：8/8 PASS，test time 3.862s，wall 4.204s，exit 0。

### Receipt命令與輸出

命令：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/bin/python3 \
research/20260718_solver_methyl_edge_probe/scripts/run_compressed_vaf_rank_probe.py \
--output \
research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v1/receipt.json
```

實際輸出片段：

```json
{
  "elapsed_seconds": 3.331929816864431,
  "hard_status": "NOT_RUN_RESOURCE_MUTEX",
  "small_oracle_pass": true,
  "verdict": "PASS_SMALL_EXACT_ORACLE_HARD_NOT_RUN_RESOURCE_MUTEX"
}
```

輸出：

- `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v1/receipt.json`
- `InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v1/receipt.json.sha256`

兩檔皆為mode `0444`、hardlink count `1`；`sha256sum -c`為PASS。

## 複雜度與無法消失的最壞情況

令 \(m\) 為active mutations、\(F\) 為logical perfect-family size、\(S\) 為實際
評分候選數、\(n\) 為quality-pattern rows。

- DP count：時間 \(O(3^m\,poly(m))\)，空間 \(O(2^m+m^2)\)。
- Lazy traceback：output-sensitive；最壞仍需走訪 \(F\) 個vertex sets。
- Ranking：可寫成
  \[
  O\!\left(
    3^m poly(m)+
    \sum_b Fit(U_b,n)+
    S\cdot Fit(m+1,n)
  \right).
  \]
- 若 bound完全不能剪，\(S=F\)，最壞仍為
  \(\Theta(F\cdot Fit)\)。
- 若有 \(T\) 個候選都落在tie tolerance內，顯式輸出完整tie IDs至少需要
  \(\Omega(T)\)時間與輸出空間。對已知1.22億family，若全部同分，就不可能靠
  一般B&B保證小型完整tie輸出。

目前prototype為了稽核會保留所有已評分 `CandidateScore`，所以記憶體還有
\(O(Sm)\) 項；後續可改成online incumbent + tie buffer，但若完整tie本身很大，
仍需exact compressed tie-class theorem或明確abstain，不能把heuristic冒充exact。

## 現階段可下與不可下的結論

可以下：

- DP traceback + safe B&B在小型 recurrence-free family上可保持current primary
  likelihood的最佳解與完整tie class。
- synthetic m=5可安全少評分96/120候選。
- cap/deadline/tie-output不完整時確實fail closed。

不可下：

- `09b`的122,281,152候選已完成VAF排名。
- hard25或所有樣本已獲得wall-time加速。
- full M2 release endpoint已被替代。
- 任意大family都能在有限時間列出完整tie class。

下一個合法步驟是在formal hard25 timing結束、取得明確資源授權後，才以
`--run-hard --hard-max-candidates 20 --hard-deadline 5`跑單一09b bounded
diagnostic，輸出至全新
`InterSubMod/research/20260718_solver_methyl_edge_probe/results/compressed_vaf_rank_probe_v2/receipt.json`。
v2 hard loader必須通過R3 pointer與sidecar、pointer schema/status、相對路徑
containment、manifest byte/sidecar/semantic digest、authority status、全部frozen
checks與case/unit結構驗證。若回未完成，仍是有效的fail-closed結果，但不是exact
winner。
