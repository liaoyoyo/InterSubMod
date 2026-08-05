<!--
建立時間: 2026-07-16
目標: 獨立 red-team 稽核 M2 symbolic partial-read Hypercube Steiner candidate generation、read-pattern likelihood ranking、全量彙整與 reuse 契約
處理範圍: 唯讀稽核、合成單元測試、已完成 HCC1954 chr22 pilot 唯讀數據檢查；未啟動 BAM 全量抽取或重運算
關聯檔案: InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py; InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_ranking.py; InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py; InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/extract_lossless_read_linkage.py; InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_m2_patterns_and_rank.py
任務類型: B_COMPREHENSIVE_VALIDATION red-team audit
服務目標: G3 / G4 / G5
狀態: RELEASE_BLOCKED_AT_AUDITED_SNAPSHOT; remediation 修訂後須重跑本文件放行 gate
-->

# M2 ranking 獨立 red-team 稽核

## 0. TL;DR 與結論

**結論：M2 的 symbolic partial-read solver 核心成立，但在本次受稽核 snapshot 中仍有 1 個 P0（PS-aware 分層）與多個 P1；在 PS、effective-k、reuse scope、partial funnel 與可重建 candidate output 通過後，才能將結果宣稱為可靠的 per-HP regional topology comprehensive validation。**

- 影響：高。
- 信心：高。
- 本輪已證實：symbolic group constraint、不建立跨 read Cartesian product、全域 minimum vertex-set enumeration、same vertex/different edge 不可區分、k 超出上限時 fail-closed 的基本契約。
- 本輪發現後已修正：V>1 未初始化 runtime blocker、unique vertex-set aggregate 漏計、primary optimizer nonconvergence 仍執行 bootstrap。
- 重要 snapshot 提醒：25 tests 是對下方「受稽核 SHA」執行；文件建立前 `build_m2_patterns_and_rank.py` 已有後續 remediation 修訂，必須依§ 8 Gate 0 重跑，不可把舊測試結果直接轉移為新版 PASS。

---

## 1. 輸入、版本與稽核範圍

### 1.1 主要輸入

| 類型 | 路徑 |
|---|---|
| M2 ranker | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py` |
| Full runner | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_full_m2_ranking.py` |
| Exact solver | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py` |
| Extractor | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/extract_lossless_read_linkage.py` |
| Lossless contract | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/lossless_read_contract.py` |
| Tests | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_m2_patterns_and_rank.py` 與 `test_hypercube_exact.py` |
| PS pilot calls | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/HCC1954.chr22.molecule_sparse_calls.tsv.gz` |
| PS pilot components | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/HCC1954.chr22.components.tsv.gz` |
| PS pilot membership | `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/HCC1954.chr22.site_component_membership.tsv.gz` |

### 1.2 25 tests 對應的受稽核 SHA

| 檔案 | SHA-256 |
|---|---|
| `build_m2_patterns_and_rank.py` | `ec73d83e21c5af97c1073e0ce05b0d60175b72c9e280e904462289c247a0727a` |
| `run_full_m2_ranking.py` | `f8c6cced22b9738d37500d9a1ea2f6d3c3e4f6da186758f23bc16c32e53633a7` |
| `hypercube_exact.py` | `db85c94f5cd1fda5bb1d4eca9ca64d5df6680829b1b005ac02be9d3c835500bc` |
| `extract_lossless_read_linkage.py` | `3ad014da5a2eadd7220d2bcc015702637d2e7606f068583499fdf782e47e4291` |
| `test_build_m2_patterns_and_rank.py` | `cbfd5d6702e38949346cb16189caee45a86b1e066deea988c548b0a346c3712f` |

### 1.3 建立本文件時的 hash drift

`build_m2_patterns_and_rank.py` 建立本文件時已變更為：

```text
219a21bd780b570515171b31cbc6b4ed761c5a5039f397412fead4d8655f53e6
```

這代表 remediation 已在進行，不代表新版已通過本 red-team audit。放行前應以最終凍結 SHA 重跑全部契約測試與 pilot。

---

## 2. Partial read 實際處理方法

現行數學契約不是「展開所有 completion，只要任一種成功就通過」，而是：

1. 一條有 `u` 個 `X` 的 read，在完整 k-cube 上有 `2^u` 個「概念 completion」。
2. Solver 以一條 symbolic group constraint 表示這條 read：

   \[
   \sum_{v\in G(p)}x_v\ge1
   \]

3. 所有 reads 的 constraints 必須同時成立，不為每條 read 生成情境、也不建立跨 read Cartesian product。
4. 先最小化 root 與 full observed states 之外的額外 vertices。
5. 只列舉全域 minimum vertex sets；到達 `max_vertex_sets=256` 或 timeout 時不強行下結論，而是 abstain。
6. 對完整列舉的 distinct vertex sets 擬合 molecule-pattern mixture likelihood。
7. 同一 vertex set 的不同 parent-edge assignments 不重複評分；snapshot read states 無法區分這些 edges。

但是 solver 使用 `universe_mode=observed_alt`，所以實際 group 大小是：

\[
2^{|free\_mask\cap structural\_ALT\_universe|}
\]

可能小於概念上的 `2^u`。例如 `RAX` 概念上有 2 種；若 X 位點未在 retained structural evidence 出現 ALT，實際 solver 可能只保留 1 種。

### 方法解釋的正確邊界

- 可以說：「對所有 retained symbolic groups 同時求全域 minimum compatible vertex sets」。
- 不可以說：「把每條 partial read 的 `2^u` 種情況都建立成獨立樹，任一成功就輸出」。
- 不可以說：「已找到最高生物機率的所有樹」。現行方法是 lexicographic minimum-parsimony 優先，likelihood 只在 minimum candidates 中排序。

---

## 3. Audit matrix

| 等級 | Finding | 影響 | 建議修正 | 必要測試 |
|---|---|---|---|---|
| **P0** | PS 未進入受稽核 snapshot 的 production component/ranking：extractor 保存 PS，但 bridge 與 unit 只按 HP family | 不同 PS block 的 HP1/HP2 orientation 可相反；混合後 component、pattern、topology 都可能錯 | Primary unit 改為 dataset×chr×PS×HP；missing PS 只進漏斗與 sensitivity；跨 PS 只有獨立 orientation anchor 才可合併 | 合成兩 PS、HP orientation swap，不得合成同一 unit；實際 pilot 需輸出 multi-PS/missing-PS funnel |
| **P1** | `component.k>12` 在 active universe 形成前直接 abstain | k=13 但只有1個 mutation-bearing bit 也會被錯誤丟棄 | 拆成 `component_site_k`、`structural_alt_active_k`、`scoring_alt_k`；exact gate 使用 active k | k=13/effective=1 應 exact；effective=13 才 abstain；aggregate 回報兩種 k |
| **P1** | `observed_alt` 禁止 retained structure 未見的 latent ALT | 可能漏掉低支持 subclone；不能宣稱測試完整 `2^u` completions | 輸出 conceptual/effective completions 與 universe 來源；比較 structural-alt、all-scoring-alt、all-loci sensitivity | `RAX` conceptual=2、effective=1/2 雙模式 golden test |
| **P1** | `minread=3` 是三條完全相同 R/A/X pattern，不是三條相容 reads | partial masks 不同時會碎裂支持 | 改名 `identical_pattern_minread`；執行 1/2/3/5 sensitivity，或設計 compatibility-aware soft constraints | `AXR`、`ARX`、`AXX` 各1條的合成案例必須有預先定義的路由結果 |
| **P1** | 現行目標是 lexicographic「先最少 extra vertices，再在平手中用 likelihood 排序」 | 多一個節點但 likelihood 很高的樹不會被比較；不等於全域最可能生物樹 | 報告明訂 minimum-compatible tree；另做 `h*..h*+q` 或 `LL-λh` sensitivity | 建立 likelihood 明顯支持 nonminimum tree 的合成案例 |
| **P1** | BQ mismatch 把 `e=10^(-Q/10)` 直接當特定 REF↔ALT flip probability | Phred e 是任一錯誤 base 的總機率；現行實際上是 clipped binary-flip surrogate | 降級命名與 claim，或保存 observed base 並使用4-base substitution model | 對稱 `e/3`、條件式 binary、現行 surrogate 三模型 sensitivity |
| **P1** | Full-ranking reuse 未驗 `scope.hp_families` | 只含 HP1 的舊結果可能被 HP1+HP2 run 誤用 | reuse contract 加 expected families 與 basis mode | HP family mismatch 必須拒絕 reuse |
| **P1** | reuse 只驗 extraction receipt，未重新 hash 實際 upstream TSV | calls/components 被修改但 receipt 未改時可能仍 `REUSED_PASS` | reuse 時重驗 ranker 消費檔案的 path confinement、size、SHA | upstream calls/components tamper regression |
| **P1** | 實際 vertex states 與 `MixtureFit.weights` 未輸出，只留不可逆 SHA IDs | 無法回答每區 clone/state 組合、比例，也無法不重跑 solver 獨立重建結果 | 新增 hashed candidate table：states、state role、LL、relative weight、π、edge count、topology class | 只依輸出表不重跑 solver 即可重建 top candidate |
| **P1** | Full aggregate 缺 partial-read 定量漏斗 | 無法量化 full/partial、`u`、`2^u`、effective completions | 按 dataset/basis/threshold 彙整 molecule、RAX-group、BQ-group 三種分母 | partial funnel 各分母與 molecule mass 守恆 |
| **P1** | Topology inclusion 非互斥，且 optimizer abstain 的 class 可能被累計 | 類別總和可超過100%，並混入不可靠結果 | 增加 unique-class counts、ambiguous-set counts、evaluated denominator；abstain 不進正式 topology counts | exclusive sum=unique；unique+ambiguous=evaluated；nonconvergence 不進 topology counts |
| **P1** | Full receipt `all_pass` 主要只看154 scope，未驗 aggregate conservation | 彙整公式錯誤仍可能顯示 PASS | 對每個 dataset×basis×threshold 強制驗 selection、solver、projection、BQ、`T≥V` 守恆 | 任一 aggregate 欄位竄改都必須 fail |
| **P1** | Full receipt 缺每樣本 input molecule/call/HP funnel | 無法清楚呈現各樣本數量與比例，也容易混淆 7 datasets 與 6 biological samples | `by_dataset` 加 input funnel、call-code counts、HP/PS counts；HCC1395 兩個 pipeline 獨立標示 | 各樣本 funnel 加總等於 global，並清楚列出 denominator |
| **P1** | Bootstrap 固定原候選，只重抽 scoring counts | 只是 conditional-on-candidate-set stability，不是整條 tree pipeline 穩定度 | 正確改名；若宣稱全流程 bootstrap，每 replicate 需重跑 threshold 與 candidate generation | 重抽後 candidate set 會變化的合成案例 |
| **P1** | 靜態最壞成本高：每 unit 最多 256 MILP×30秒，另有約 `V×(1+4+20)=25V` 次 SLSQP | 全量可能花數天或無法完成 | active-bit compression、group dominance、先跑 primary/grid、bootstrap 僅對預登記子集 | chr22 pilot 回報 V 分布、wall time、timeout/max-set 比例後才決定全量 |
| **P2** | `no_partial_completions_materialized`、`same_read_vaf_not_added` 等 receipt checks 是 hard-coded `True` | Receipt 本身不能偵測未來 code regression | 移到 method declaration，或改為可由 row/hash/test 驗證的條件 | 故意引入 completion expansion/VAF 欄位的 negative contract test |

---

## 4. PS 問題的實際 pilot 證據

### 4.1 輸入範圍

- Dataset：HCC1954。
- Chromosome：chr22。
- 資料：已完成的 v1 pilot，未讀取當時仍在寫入的 v2 pilot。
- 用途：證明 PS 混合非純理論風險；不可以代替 v2/full 數字。

### 4.2 實際結果

| Metric | 結果 |
|---|---:|
| Sparse molecule rows | 15,671 |
| Distinct PS | 50 |
| HP1 rows | 5,116 |
| HP1 missing PS | 173 / 5,116 = 3.38% |
| HP2 rows | 9,519 |
| HP2 missing PS | 75 / 9,519 = 0.79% |

k>1 pooled-component 映射後，component×HP unit 包含超過一個 PS 的比例：

| Bridge threshold | HP1 | HP2 |
|---:|---:|---:|
| 1 | 2 / 32 = 6.25% | 1 / 43 = 2.33% |
| 2 | 2 / 29 = 6.90% | 1 / 39 = 2.56% |
| 3 | 1 / 28 = 3.57% | 0 / 38 = 0% |
| 5 | 1 / 23 = 4.35% | 0 / 32 = 0% |

Repo 內部方法學 SoT 也明確記載：

- `InterSubMod/docs/method_comparison/20260621_clone_subclone_reconstruction_landscape_and_ism_feasibility_01.md:50,111`：bulk long-read 無跨 phase-set linkage，per-block call 不可無 anchor 串成 lineage。
- `InterSubMod/docs/provenance/ai_sessions/2026/04/20260412_V5_somatic_fallback_getVote_修正與驗證_01.md:547-550`：HP1/HP2 label 在每個 PS block 內任意定向，必須 per-PS orientation correction。
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/lossless_read_contract.py:126-146`：原契約將 pattern 定義為 molecule×HP×PS，受稽核 production ranker 當時與此契約有 drift。

---

## 5. Structural threshold 反例

合成三條 molecules：

```text
AXR × 1
ARX × 1
AXX × 1
```

三條 molecules 都支持第一位 ALT，但因為沒有一個完全相同 pattern 達到 `minread=3`，實際結果是：

```json
{
  "molecule_component_projections": 3,
  "informative_scoring_molecules": 3,
  "structural_retained_molecules": 0,
  "n_structural_pattern_groups": 0,
  "selection_status": "STRUCTURE_ABSTAIN_NO_MINREAD_ALT_PATTERN",
  "candidate_generation_status": "NO_STRUCTURAL_ALT_PATTERN_AT_DECLARED_MINREAD"
}
```

因此報告必須將 `minread` 解釋為「identical R/A/X pattern count」，不可簡寫成「至少3條 reads 支持」。

---

## 6. 執行命令、輸出與實際結果

### 6.1 輕量單元測試

**工作目錄**：`/big7_disk/liaoyoyo2001/InterSubMod`

**執行命令**：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest -v \
research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_hypercube_exact.py \
research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_m2_patterns_and_rank.py
```

**實際 stdout 摘要**：

```text
Ran 25 tests in 14.090s

OK
```

**實際結果**：25 / 25 PASS，exit code 0。

**輸出路徑**：無新分析檔；單元測試僅輸出 stdout/stderr。

### 6.2 PS pilot 唯讀 streaming join

**輸入路徑**：

```text
InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/HCC1954.chr22.molecule_sparse_calls.tsv.gz
InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/HCC1954.chr22.components.tsv.gz
InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m2_pilot_HCC1954_chr22/HCC1954.chr22.site_component_membership.tsv.gz
```

**處理方法**：gzip TSV streaming；以 threshold×component×HP 映射 molecule 的 distinct PS，未寫新檔。

**實際輸出摘要**：

```text
rows=15671 distinct_ps=50
t=1 HP1 k>1_units=32 multiPS=2 fraction=0.062500 maxPS=2
t=1 HP2 k>1_units=43 multiPS=1 fraction=0.023256 maxPS=2
t=2 HP1 k>1_units=29 multiPS=2 fraction=0.068966 maxPS=2
t=2 HP2 k>1_units=39 multiPS=1 fraction=0.025641 maxPS=2
t=3 HP1 k>1_units=28 multiPS=1 fraction=0.035714 maxPS=2
t=3 HP2 k>1_units=38 multiPS=0 fraction=0.000000 maxPS=1
t=5 HP1 k>1_units=23 multiPS=1 fraction=0.043478 maxPS=2
t=5 HP2 k>1_units=32 multiPS=0 fraction=0.000000 maxPS=1
```

**輸出路徑**：無；本文件§ 4 是實際 stdout 結果的持久化紀錄。

### 6.3 唯讀稽核輸出

本次唯一新增檔案：

```text
InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2_ranking獨立red-team稽核_01.md
```

本輪未修改任何 solver、ranker、extractor 或 tests。

---

## 7. 已通過契約與未足之處

### 7.1 已驗證

- `SymbolicPattern.compatible` 在 k≤6 的所有 R/A/X patterns 與 explicit completion 一致。
- Partial read 保持為一條 symbolic group constraint，沒有建立跨 read Cartesian product。
- `RRA + AAA + RAX` 反例得到唯一的全域 minimum vertex set，不是「先找到任一 completion 就成功」。
- Exact enumeration 只在證明無其他 optimum 時標 `complete=true`；max-set/timeout 不會被美化為唯一解。
- 候選先依 vertex set collapse；same vertex/different edge 不會被 read likelihood 錯誤分開。
- 每個 unit 內每條 molecule 有一個 pattern projection；BQ 使用於 fixed R/A calls，X 以 multiplicative-one 邊際化。
- Primary optimizer nonconvergence 後不再執行可宣稱為 complete 的 bootstrap。
- `quality_primary_unique + tied + abstain = units` 已加入守恆測試。

### 7.2 不可由上述 PASS 推論

- 25 tests PASS 不等於 154 tasks 全量數據 PASS。
- BQ surrogate PASS 不等於 ONT error model 已校準。
- Vertex-set likelihood unique 不等於 parent edges unique。
- Coarse topology class unique 不等於 exact topology isomorphism unique。
- Read-level mutation-state mixture 不等於已解決 CN、purity、multiplicity 或跨 HP clone pairing。
- HCC1954 chr22 pilot 不等於 7 datasets×chr1-22 全量結論。

---

## 8. 放行 gate

### Gate 0：版本凍結與重測

- [ ] 凍結 extractor、component builder、ranker、full runner、report builder 與 tests 的 SHA-256。
- [ ] 以凍結 SHA 重跑所有 M2 unit tests，不可沿用舊 hash 的 25/25 PASS。
- [ ] 重跑 PS swap、effective-k、reuse tamper、nonconvergence、aggregate conservation 專用 regression tests。

### Gate 1：PS-aware primary contract（P0）

- [ ] Primary component/unit 明確綁定 dataset×chr×PS×HP family。
- [ ] Missing PS 不可混入 known-PS topology；需獨立漏斗與 sensitivity。
- [ ] 跨 PS 不可僅因 HP 數字相同就合併；必須有外部 orientation anchor 或保持 abstain。
- [ ] 兩 PS orientation-swap synthetic test PASS。
- [ ] 輸出各樣本 PS 漏斗、multi-PS component 數與比例。

### Gate 2：effective mutation-bearing k

- [ ] Exact solver gate 使用 observed-alt active k，不是物理 component site k。
- [ ] k=13/effective=1 可精確求解；k=13/effective=13 維持 fail-closed。
- [ ] 輸出 component k、active k、R-only dimension 數與路由漏斗。

### Gate 3：partial-read 定量與主張邊界

- [ ] 同時輸出 conceptual `2^u` 與 solver-effective completion count。
- [ ] 全量彙整 full/partial/all-X molecules、RAX groups、BQ groups、`u` 分布與 effective-completion 分布。
- [ ] `minread` 文字定義明確寫成 identical-pattern count，並完成 1/2/3/5 sensitivity。
- [ ] 報告明訂只列舉 global minimum candidates，不宣稱已排序所有 nonminimum feasible trees。

### Gate 4：likelihood 與可重建輸出

- [ ] BQ 明確標示為 binary-flip surrogate，或完成4-base/error calibration sensitivity。
- [ ] Same-read VAF 不作為第二個獨立 score。
- [ ] Bootstrap 明確命名 conditional-on-candidate-set，或改為每 replicate 重建 candidates。
- [ ] 持久化 candidate vertex states、mixture weights、likelihood、state roles 與 topology class，不只留 SHA ID。
- [ ] Optimizer abstain 不可進入正式 topology 數量與比例。

### Gate 5：receipt、reuse 與 aggregate

- [ ] Full reuse 驗證 HP families、component bases、thresholds、ranker SHA、upstream receipt SHA 與實際 upstream output hashes。
- [ ] Upstream calls/components 竄改時 reuse 必須 fail。
- [ ] 每個 dataset×basis×threshold 通過 selection、solver、projection、BQ、T/V 與 topology denominator 守恆。
- [ ] 正式互斥比例使用 unique-class counts；inclusion counts 另行標示可重疊。
- [ ] `by_dataset` 含 input molecule/call/HP/PS funnel，並清楚區分 7 datasets / 6 biological samples。

### Gate 6：pilot 效能與資源

- [ ] 先在一個已登記 chromosome pilot 回報 units、V 分布、MILP 次數、SLSQP 次數、wall time、RSS、timeout/max-set 比例。
- [ ] 資源 gate 無使用者衝突作業後才啟動全量。
- [ ] 若 pilot 不可控，先做 active-bit compression、dominance reduction、bootstrap 分離，不盲目啟動 154 tasks。

### Gate 7：全量與報告

- [ ] 7 datasets×chr1-22 = 154 extraction tasks 全數 PASS/REUSED_PASS，無 subset 當作 final。
- [ ] 154 ranking tasks 全數 PASS/REUSED_PASS，每個 child receipt 與 output hash 可追溯。
- [ ] 最終 HTML 只讀取通過上述 gate 的 receipts；任一 M2 full input 缺失時必須 fail final，不可只加註解繼續輸出。
- [ ] 各數字同時列總比例、相對比例與明確 denominator。
- [ ] 獨立 reviewer 復核方法 claim、數據 funnel、拓樸類別互斥性與 7-dataset 範圍後才放行。

---

## 9. 最終 verdict

### 可保留的方法核心

1. 將 partial read 保持為 symbolic subcube constraint。
2. 所有 read groups 同時求解，不建立 completion Cartesian product。
3. 用 exact minimum-compatible vertex-set enumeration 作為可稽核的 parsimony layer。
4. 先依 vertex set collapse，再用 read-pattern likelihood 排序。
5. Same vertex set 的 parent edges 保持 unresolved，不用相同 snapshot reads 虛構 edge 權重。

### 目前不放行的主張

在 Gate 0-7 尚未全數通過前，不放行以下結論：

- 「已完成 7 datasets 的可靠 per-HP clone topology 重建」。
- 「partial read 的所有 `2^u` 情境都已逐一枚舉與排序」。
- 「BQ/VAF 已選出唯一真實 parent-edge topology」。
- 「mixture weights 已等價於可跨 HP 配對的真實 clone fractions」。

**Release decision：BLOCKED，直到最終 remediation SHA 完成獨立重測、PS-aware pilot 與 154-task 全量 receipts 驗證。**
