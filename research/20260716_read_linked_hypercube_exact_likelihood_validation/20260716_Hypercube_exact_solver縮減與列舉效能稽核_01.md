<!--
建立時間: 2026-07-16 09:20 +08:00
目標: 稽核並安全縮減 Boolean-hypercube rooted group-Steiner exact MILP，保留 partial-group、minimum-hidden 與 all-optimal vertex-set 完整語意
處理範圍: static review + single-core small synthetic（每案 hard limit <=10 秒）；不讀 BAM、不執行 k>=13、不修改 canonical output
關聯檔案: InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py；InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py
任務類型: B Comprehensive validation 的 bounded solver PROBE
研究目標: G3 / G4 / G5
-->

# Hypercube exact solver：語意等價縮減與 all-optimal 列舉效能稽核

> **TL;DR — exact presolve 已安全落地，本稽核當時共享 cycle 108/108 tests PASS；但「每個最佳 vertex set 再解一次 MILP」的漸近 tail 仍存在，因此 code-level GO、154-task scale-level 仍是 PROBE（影響：中，信心：高）。**

> **Post-audit identity note（2026-07-16 10:20）**：其後只修正module docstring，把「完全不materialize completions」改成實際行為——每次MILP build、每個reduced group的sparse row會列出active compatible vertex indices，但不建立completion-wise tree worlds或cross-read Cartesian product；solver邏輯未變。current SHA256=`9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95`，identity-bound pilot receipt=`results/pilot_v3_identity_bound/pilot_receipt.json`（9,840 patterns／2,015,538 checks／0 mismatch）。
>
> 下文的`891c0e...`與108/108皆是**audit-time historical snapshot**，不是current identity／current suite總數。

## 1. 範圍、前提與不做的事

本稽核服務 G3（read-level mutation-state model）、G4（跨樣本可重現）與 G5（可外部驗證）。它不是生物結果驗證，也不會把候選樹升級成真實 clone tree。

啟動研究任務 5 問：

1. Thread D：**間接相關**；本輪只處理 read pattern 的結構 constraint，不處理 methylation。
2. Thread B 撤回範圍：**否**。
3. KDE-corrected：**不適用**。
4. VCF caller AF：**不使用**；結構 objective 不加入 VAF。
5. 長計算／C++／搬檔／NO-GO gate：**只做 bounded Python synthetic**；不讀 BAM、不跑 k≥13、不碰 canonical output。

關鍵假設與 claim boundary：

- ROOT 固定為 `0^k`；邊只允許一次 `R→A`，因此 popcount 嚴格增加，DAG 無 cycle。
- full read state 與 ROOT 是 zero-cost mandatory vertices；其他被選 vertex 都保留原 objective cost。
- partial read 是一個 group constraint，不是先展開多個獨立世界再「任一成功」。
- 只有 `complete=true` 才代表已列完所有 minimum vertex sets；`MAX_SETS_REACHED` 或 time limit 必須 abstain。
- 本稽核只證明 solver 的候選 vertex-set 語意，不證明 parent edge、clone 數、cellular ancestry 或生物真值。

## 2. 實際 partial-read / minimum-hidden 語意

對長度 `k` 的 pattern `p∈{R,A,X}^k`：

- `fixed_mask`：R/A 已觀察的位置。
- `alt_mask`：已觀察為 A 的位置。
- `free_mask`：X 位置。
- 概念上的 full-cube completions 為 `2^u`，其中 `u=#X`。
- 在 observed-ALT structural universe `U` 中，實際有效 completions 是 `2^popcount(free_mask & U)`；X-only、從未觀察到 A 的維度會固定為 R。

一條 partial read 建立一個 subcube group：

`G(p)={v⊆U : ((v xor alt_mask) & fixed_mask)=0}`。

所有 reads 的 groups **一起**進入同一個 MILP，逐一要求：

`sum[v∈G(p)] x_v >= 1`。

另對每個非 ROOT vertex `v`：

`x_v <= sum[x_pred]`，其中 `pred` 是少一個 ALT 的 Hamming-1 predecessor。

因 predecessor 每次嚴格少一個 ALT，這個 local constraint 會遞迴保證 ROOT connectivity。目標仍是：

`min sum[x_v : v 不是 ROOT/full-read mandatory state]`。

求到 minimum hidden/extra 數 `h*` 後，solver 固定 `sum extra x=h*`，保留所有可行的 global minimum vertex sets。後續 read-pattern likelihood 才對 X/latent states 做 marginalization 並排序 vertex sets；相同 vertex set 的不同 parent-edge tree 在 snapshot read/VAF 下不可區分。

## 3. 已落地的 exact reductions

| reduction | 實作 | 為何不改答案 |
|---|---|---|
| duplicate group removal | 相同 effective `(fixed_mask, alt_mask)` 只留一條 constraint | 結構是 Boolean feasibility；相同不等式重複幾次不改 feasible set。read count 仍由 likelihood 使用，未被刪除。 |
| dominated group removal | 若 `G1⊂G2`，只留 `G1` | `G1` 被命中必然使 `G2` 被命中；刪 `G2` 不改 feasible set。 |
| required-hit / singleton forcing | 已被 ROOT/full/forced state 命中的 group 移除；單一 completion 改為 `x_v=1` | 前者是不等式恆真；後者與原 `sum_G x>=1` 在 `|G|=1` 時代數相同。forced partial vertex **仍保留原 objective cost**。 |
| active-bit predecessor loop | predecessor 只巡訪 `U` 中 active bits | vertex 原本就只可能包含 U 中 bits；刪除 inactive-bit 檢查不改 predecessor 集合。變數數仍是 `2^effective-k`。 |
| sparse no-good | 固定 `sum extra x=h*` 後，用 `sum[e∈E_old]x_e<=h*-1` | 每個候選都有恰好 h* 個 extra vertices；少任一舊 extra 才可能是不同集合。與 dense binary Hamming no-good 完全等價。 |
| zero-hidden early certificate | 若 `h*=0`，第一個集合即為唯一集合 | 所有 nonmandatory variables 都被 objective equality 固定為 0；不需再跑一次 infeasible MILP。 |

### 沒有做的「mandatory closure」

除了 singleton group 外，沒有武斷固定某條祖先路徑。完整 Boolean cube 中，popcount >1 的 vertex 通常有多個 predecessor；任選其一會刪掉其他同樣最佳的 topology。只有在未來另有 exact vertex-pruning 後出現唯一 predecessor，才可迭代做 closure。

### 沒有做的 solver replacement

SciPy `milp`/HiGHS 介面沒有在此 workflow 使用 solution-pool/callback/persistent model；因此列舉 V 個最佳 vertex sets 仍需第一個 optimum、V−1 個新解與最後一個 infeasible certificate，最壞仍是 O(V+1) 次 MILP。這次縮減每個模型的 rows/nonzeros，但沒有宣稱消除漸近 tail。

## 4. Step → Verify

1. 鎖定 optimizer-agent 完成後的 stable baseline → 驗證：`hypercube_exact.py` SHA-256=`6576910b52cd8a3463b9bbe3d22fde5d735f5af78ae2c512f0bf4b67031d5fb8`。
2. 以 algebraic proof 篩選 exact reductions → 驗證：每項皆能寫成等價 constraint 或邏輯蕴含；不加入 read count/VAF edge weight。
3. 新增 independent exhaustive oracle → 驗證：k=1..3、18 個 deterministic seeded cases 的 objective 與**完整 vertex-set 集合**逐案一致。
4. 驗證典型多解 → 驗證：`full=AAAA,k=4` 仍 complete 列出 24 個 minimum chains，集合 digest 前後相同。
5. 跑本 cycle 全測試 → 驗證：本稽核當時共享 worktree 108/108 PASS、exit 0、wall 4.10 秒、peak RSS 85,312 KB；後續新增tests不回填成這次歷史執行結果。

## 5. Bounded synthetic 實測

所有命令都在 `/big7_disk/liaoyoyo2001/InterSubMod` 執行，並固定：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
```

### 5.1 重複／支配 groups，k=8

輸入：50 groups = 5 種 pattern 各重複 10 次；`universe_mode=all_loci`；256 variables。

| 指標 | baseline | reduced | 核對 |
|---|---:|---:|---|
| status | OPTIMAL | OPTIMAL | 一致 |
| objective | 5.0 | 5.0 | 一致 |
| variables | 256 | 256 | 一致 |
| constraints | 305 | 257 | −48（−15.74%） |
| groups input → active | 50 → 50 | 50 → 2 | 45 duplicate + 3 dominated |
| 單次 wall | 4.136 s | 3.963 s | 描述性 −4.19%；單次量測，不作穩定 speedup claim |

兩次 solver 可回傳不同的其中一個 optimal incumbent，這不代表候選集合改變；完整性由 algebraic proof 與 exhaustive-oracle tests 驗證。

### 5.2 all-optimal chain，k=4

輸入：full state `AAAA`、無 partial group、`max_sets=32`。

| 指標 | baseline | reduced |
|---|---:|---:|
| objective | 3 | 3 |
| complete | true | true |
| optimal vertex sets | 24 | 24 |
| sorted set digest | `7ec0afac5a6df66f45246df34605a8b546f7c440f867ed52745f6322c277c60e` | 同左 |
| wall | 0.257 s | 0.217 s |

固定 h*=3 時，一條 no-good row 從 16 個 nonzeros 降為 3 個（−81.25%）；這是 exact matrix sparsity 結果，不是以近似法換速度。

### 5.3 hard-limit 記錄

最初的 k=10、200 redundant-group stress 超過本稽核 10 秒上限，已立即終止且未取結果；沒有把 timeout incumbent 當 optimum。之後所有正式比較皆縮至 k≤8 且單案 <10 秒。

## 6. 實際驗證命令與輸出片段

### Targeted solver tests

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /usr/bin/time -f 'ELAPSED=%e RSS_KB=%M EXIT=%x' \
  python3 -m unittest \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py -v
```

實際片段：`Ran 21 tests in 1.509s`、`OK`、`ELAPSED=2.08 RSS_KB=70024 EXIT=0`。

### Full cycle tests

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
/usr/bin/time -v env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python3 -m unittest discover \
  -s research/20260716_read_linked_hypercube_exact_likelihood_validation/tests \
  -p 'test_*.py' -v
```

首次完整執行片段：`Ran 105 tests in 2.739s`、`OK`、elapsed `0:03.43`、maximum RSS `85128 kbytes`、exit status `0`。其他 agent 同步新增 3 tests 後，以當時 source 再跑：`Ran 108 tests in 3.266s`、`OK`、`ELAPSED=4.10 RSS_KB=85312 EXIT=0`。這是audit-time execution receipt，不代表後續共享worktree的current test總數。

### Syntax validation

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
rg --files research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests \
  | rg '\.py$' | xargs -r python3 -m py_compile
```

實際結果：無 stderr、exit 0。

## 7. Verdict 與 154-task 前的 gate

### Code-level：GO

- partial group joint feasibility、minimum-hidden objective、all-optimal vertex-set completeness均保持。
- 新 reductions 有逐項代數證明、seeded exhaustive oracle、典型 24-optima digest 與最終 108-test regression 支持。
- 沒有使用 read count/VAF 當 edge evidence，也沒有刪除 likelihood 中的 read counts。

### Scale-level：PROBE

在 HCC1954 chr22 PS-aware pilot 前，不應宣稱 154 tasks 已解決 solver tail。pilot 必須至少記錄：

- effective-k、partial groups input/active/reduction reasons；
- optimum h*、V enumerated、`complete`/stop reason；
- 每 unit solve/enumeration wall、MILP node/gap；
- `MAX_SETS_REACHED` / time-limit 的 abstain 數；
- candidate-table 規模與 peak RSS。

若 tail 仍由 V 很大主導，下一步應評估可提供 exact solution pool 或 persistent no-good model 的 backend；或以 BDD/ZDD / DP 儲存所有 minimum vertex sets 的壓縮表示。任何 beam search、只留 top-N、或先選一個 partial completion 的作法都只能是 heuristic，必須明示不完整，不能冒充 `complete=true`。

## 8. 交付檔案與 SHA-256

| 檔案 | SHA-256 |
|---|---|
| baseline `hypercube_exact.py`（修改前） | `6576910b52cd8a3463b9bbe3d22fde5d735f5af78ae2c512f0bf4b67031d5fb8` |
| `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`（exact reductions audit-time historical snapshot） | `891c0e469c8b24f37a4e98668a564f46c2cd505ed3da633f9bdd63b6f7294aa7` |
| 同檔案（current；post-audit wording-only correction） | `9dbaf8ec813d709ba55e1c45ca08bcb4220fc15070c764aa2949ab27608bcd95` |
| `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py` | `6872886288389821c82f0b9e5cac1aa46c06c687bfd30c352b189398346f2adf` |

本文件 SHA 應由交付時外部 `sha256sum` 計算，避免 self-referential hash。
