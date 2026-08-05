<!--
建立時間: 2026-07-23
目標: 獨立驗證 strict two-sSNV region edge、k<=12 block、MINREAD pattern 與 tree solver 實際輸入語意是否一致
處理範圍: HCC1395 chr1-22；strict regions v2；partition v2；authoritative molecule sparse calls；MINREAD=3
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/exact_ps_partition_to_mlhp.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/layered_tree_reconstruction.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/tree_enumeration_solver.py
  - InterSubMod/research/20260723_strict_two_ssnv_phaseblock_validation/agents/audit_tree_input_semantics.py
狀態: VALIDATED READ-ONLY AUDIT；未修改 production
-->

# Strict two-sSNV region 到 tree solver pattern 的語意稽核

**TL;DR：region 建立時的「兩個 sSNV 有 ≥3 molecules 共同 fixed observation」與 solver 真正收到的「某一完整 R/A/X 向量本身 count≥3」不是同一條件。HCC1395 的 76,202 條 strict edges 中，只有 30,344（39.82%）仍以同一 solver-visible pattern pair 出現；11,590 個 adapter tree-input blocks 中，383 個只有 single-fixed-site patterns，1,873 個 supported-pattern graph 無法連通全部 block sites。因此目前可稱「strict read-linked source region＋partial candidate inference」，不可把每個 adapter group 都稱為 block-wide two-sSNV topology。（影響：高；信心：高）**

用 SCQA：

- **Situation**：strict graph 用 pairwise molecule support 建立 source W；k>12 再切成 bounded blocks。
- **Complication**：adapter 不直接把 strict edges 送入 solver，而是重投影 molecule，依「完整 R/A/X vector count≥MINREAD」篩 pattern。
- **Question**：通過 two-sSNV region gate 的 pair evidence，是否完整保留到 solver？
- **Answer**：否。partition boundary 與 exact-vector MINREAD 都會讓 pair relation 在 solver input 消失；solver 仍可對 incomplete patterns 枚舉候選，但那是模型相容解，不等於每個 block 都有 block-wide read linkage。

本任務分類為 **(B) Comprehensive validation**，服務 **G1、G4、G5**。

## 1. 四層定義不能混用

| 層 | 單位 | 實際條件 | 能主張 | 不能直接主張 |
|---|---|---|---|---|
| Strict region edge | endpoint pair | 同 exact PS×HP；`RR+RA+AR+AA≥3` distinct molecules | 兩位點有 pairwise fixed-call linkage | 該 pair 一定會進 solver |
| Bounded block | `B_ij` | source W 經 k≤12 partition | local computational domain | 跨 blocks 的 W topology 已保留 |
| Adapter pattern | 完整 R/A/X vector | **同一 vector** count≥MINREAD=3 | 該 exact／partial pattern 被 solver 看見 | 所有 strict pair evidence 都被保留 |
| Tree solver | pattern presence | full keys＋partial keys；count 不進結構目標 | 最小相容 candidate trees／ambiguity | read-weighted likelihood、cellular clone truth |

程式證據：

- `exact_ps_partition_to_mlhp.py:431-462` 將每個 molecule 投影成 block-wide R/A/X vector。
- `exact_ps_partition_to_mlhp.py:470-480` 只檢查 `k≥2` 與至少一個 vector count≥MINREAD；**沒有要求 pattern 固定 ≥2 sites，也沒有要求 pattern graph connected**。
- `exact_ps_partition_to_mlhp.py:481-517` 將無 X 的 vector 放入 `populations_by_hp`、有 X 的放入 `subread_groups_by_hp`。
- `layered_tree_reconstruction.py:122-138` 把 full dictionary 與 partial pattern keys 送入 `enumerate_min_trees`。
- `tree_enumeration_solver.py:113-122` 對 full dictionary 只使用 keys；partial 也只有 pattern presence。位點只有在某 supported pattern 曾為 `A` 時才進 mutation universe。

所以「region 通過至少兩個 sSNV」只描述上游 raw pairwise linkage；不能自動改寫成「solver 每個 pattern 至少穿過兩個 sSNV」。

## 2. HCC1395 block 到 solver-input 漏斗

本次直接重讀 22／22 chromosome 的 partition、strict edge、extraction receipts 與 molecule sparse calls；四類 receipts／Python-C++ comparison 全部 22／22 PASS。

### 2.1 Block eligibility

| 階段 | blocks | 說明 |
|---|---:|---|
| bounded blocks | 11,712 | partition 全部 child blocks |
| k=1 child blocks | 12 | 排除 tree |
| k≥2 blocks | 11,700 | strict structural candidates |
| 沒有任何 MINREAD-supported vector | 110 | adapter 排除 |
| **依現行 adapter 會成為 tree-input groups** | **11,590** | `k≥2 AND any supported pattern` |

因此「至少一個 MINREAD-supported pattern」在 tree-input denominator 中是 **11,590／11,590（100%）**，但這只是存在一個 vector，不保證該 vector同時固定兩個位點。

### 2.2 Tree-input pattern 語意

| 檢核（分母 11,590） | 通過 | 未通過 | 比例 |
|---|---:|---:|---:|
| 至少一個 supported pattern 固定觀測 ≥2 sSNVs | 11,207 | **383** | 96.70%／3.30% |
| supported patterns 合併後覆蓋 block 全部 sites | 11,282 | **308** | 97.34%／2.66% |
| supported-pattern induced graph 連通全部 sites | 9,717 | **1,873** | 83.84%／16.16% |
| 至少一個 full-coverage supported pattern | 8,274 | 3,316 partial-only | 71.39%／28.61% |
| 至少一個 ALT（mutation-bearing） | 9,624 | 1,966 reference-only | 83.04%／16.96% |
| ALT-supported mutation loci ≥2 | 6,365 | 5,225 僅 0–1 mutation locus | 54.92%／45.08% |
| pattern-connected 且 ALT loci≥2 | **5,009** | 6,581 | **43.22%** |

細節：

- 383 個「only single-fixed-site patterns」中，290 個 patterns 合併後雖覆蓋全部 sites，仍沒有任何 pattern edge；93 個連 coverage 都不完整。
- 1,873 個 pattern-disconnected blocks 中，1,565 個其實覆蓋全部 sites，只是分成兩個以上互不共 fixed-observation 的 pattern components；另 308 個有 sites 完全未被任何 supported pattern fixed-observe。
- 1,490／1,873 個 disconnected blocks 至少有一個 multi-site pattern，但多個 local pattern components 仍沒有連成整個 block。
- 所有 11,590 tree-input blocks 的 upstream strict induced graph 都 connected；因此這 1,873 個 mismatch 是 **MINREAD pattern 表示層失去連通**，不是原 region 沒有 pairwise edge。

現行 11,590 groups 可分成互斥的 downstream claim tiers：

| 建議 tier | blocks | 語意 |
|---|---:|---|
| `TREE_READY_MULTIMUTATION` | **5,009** | pattern graph connected，且至少 2 個 ALT-supported mutation loci |
| `SINGLE_MUTATION_CONTROL` | 2,885 | pattern graph connected，但只有 1 個 ALT locus；只能形成單突變／trivial tree |
| `REFERENCE_ONLY_CONTROL` | 1,823 | pattern graph connected，但 mutation universe=0 |
| `ABSTAIN_PATTERN_DISCONNECTED` | **1,873** | current solver 可枚舉 candidate completions，但不可稱 block-wide read-linked topology |
| `ABSTAIN_PATTERN_UNSUPPORTED` | 110 | k≥2，但沒有 vector 達 MINREAD |
| `ABSTAIN_K1_CHILD_BLOCK` | 12 | 沒有 two-site tree 問題 |
| **合計** | **11,712** | 守恆 |

## 3. Strict pair evidence 在哪裡消失

76,202 條 threshold-3 strict edges 可完整分解：

| disposition | strict edges | 占全部 strict edges | solver 是否仍看見 pair relation |
|---|---:|---:|---|
| 兩端被 partition 分到不同 blocks | **35,008** | **45.94%** | 否；只保留 source component identity |
| 同 block，但該 block pattern unsupported | 3,170 | 4.16% | 否；整個 block 不輸出 |
| 同 tree-input block，但 exact vectors 分散，沒有任何 vector≥3 同時 fixed 兩端 | **7,680** | **10.08%** | 否 |
| 同 tree-input block，至少一個 supported vector fixed 兩端 | **30,344** | **39.82%** | 是 |
| **合計** | **76,202** | **100%** | 守恆 |

45.94% 的 cross-block 是 **edge-weighted** 比例，且集中在 90 個大型 source W；不能誤讀成 45.94% 的 regions 都失敗。判斷分析單位影響時，應優先看下方 source-W 與 block-level 比例。

換另一個分母看：

- 41,194 條 within-block strict edges 中，30,344（73.66%）solver-visible，10,850（26.34%）不再以 pair relation 出現。
- 11,590 個 tree-input blocks 中，2,321 個至少少了一條 internal strict edge；其中 1,873 個因此讓 pattern graph 整體斷開，另 448 個雖少 edge 但仍保持連通。
- 11,462 個 source W 中，2,334（20.36%）至少一條 strict edge 未到 solver-visible pattern；411（3.59%）沒有任何 strict pair relation留在 solver pattern layer。
- 90／90 個 k>12 source W 有 cross-block strict edges；k≤12 W 則為 0。這些大型 W 的 local block trees不能自動拼回 global W tree。

「edge 消失」的精確意思是：molecules 或單點 calls 不一定消失，但 **兩 endpoint 的共同關係**不再出現在任何單一 solver input。

### 3.1 為何 pairwise threshold 與 vector threshold 不等價

假設三個 molecules 在位點 A、B 都 fixed，但對第三位點 C 的狀態不同：

```text
ARX ×1
ARR ×1
ARA ×1
```

在 strict edge 層，A–B pair support=`3`，所以 edge 合格；在 adapter 層，三個完整 vectors 各 count=`1`，全部低於 MINREAD=3，所以沒有任何 supported vector 保留 A–B pair。這不是算術錯誤，而是 aggregation grain 改變。

## 4. Solver 會處理 single-fixed patterns，但「determined」語意較窄

synthetic spot-check：

```text
partial=["AX"], k=2
→ class=determined, n_trees=1, tree ROOT→H_AR, universe=[0]

partial=["AX","XA"], k=2
→ class=ambiguous, n_trees=3, universe=[0,1]

partial=["RX"], k=2
→ class=determined, n_trees=1, empty tree, universe=[]
```

這不是 solver implementation bug；它符合 incomplete-pattern model：

- `AX` 只要求第一個 mutation bit 為 A；
- 第二位點從未有 ALT，故不進 mutation universe；
- solver 的 `determined` 是「在目前 universe 與最小 hidden-node 模型下只有一個 candidate」，不是「一條 read 已共同觀察兩個 sSNVs」。

因此口試不可說：

> 每個送進 solver 的 block 都至少有三條 reads 同時穿過兩個突變位點。

應改為：

> Source region 由 pairwise fixed-call support 建立；partition 後，solver 使用達 MINREAD 的 full／partial R/A/X patterns。部分 block 的 solver patterns 只提供單點或不連通的局部 constraints，因此我們另外檢查 pattern coverage 與 connectivity，並對不足者 abstain。

## 5. Gate 建議

### 5.1 不要只加一個「有 pattern」gate

現行：

```text
k_B >= 2
AND exists exact R/A/X pattern with count >= MINREAD
```

只擋掉 110 blocks，仍放入 383 個 only-single-fixed blocks 與 1,873 個 pattern-disconnected blocks。

### 5.2 論文 primary topology gate

若要主張「block-wide joint mutation-state topology」，至少要求：

```text
k_B >= 2
AND strict_induced_graph_connected
AND exists MINREAD-supported pattern
AND supported_pattern_graph_connected_across_all_block_sites
AND number_of_ALT_supported_loci >= 2
```

HCC1395 在此 gate 下是 **5,009／11,590（43.22%）** tree-input groups。其他 groups 不應刪除資料，而應分流為 single-mutation control、reference control 或 abstain。

### 5.3 更佳的 representation 修正

直接丟棄 1,873 個 disconnected groups 會浪費 upstream raw pair evidence。較好的後續工程方案：

1. 在 adapter 同時輸出 `strict_edge_constraints` 或 pairwise aggregated R/A state counts。
2. 對 pair support 套 MINREAD，而不是先依整個 k-dimensional vector 分桶後再套 MINREAD。
3. solver 明確區分：
   - observed exact／partial genotype groups；
   - upstream pairwise linkage constraints；
   - partition-cut edges（只作 source-W provenance，不進 local topology）。
4. 若暫時不改 solver，至少輸出：
   - `pattern_fixed_site_max`
   - `pattern_site_coverage`
   - `pattern_graph_component_count`
   - `strict_edges_total/internal/cross_block/solver_visible`
   - `claim_role`

### 5.4 k>12 的 claim ceiling

partition 把 35,008 條 strict edges跨 block 分離，並不是 partition 錯誤；它是 k≤12 計算限制的代價。正確敘述是：

> k>12 的 source W 保留完整 read-linked component identity；每個 B 只輸出 local candidate tree。跨 block strict edges 被記錄為 cut evidence，在完成且驗證 stitching 前，不宣稱已重建 W-level global topology。

## 6. 實際執行狀態

`hcc1395_partition_v2` 目前沒有：

- `run_receipt.json`
- `strict_exact_ps_mlhp.json`
- `layered_reconstruction.json`
- topology receipt

而 adapter `exact_ps_partition_to_mlhp.py:564-569` 要求 partition root 必須有 PASS `run_receipt.json`。所以本報告的 11,590 是對 authoritative molecule data **完整重現 current adapter eligibility contract 的 would-be tree input**，不是宣稱 strict-v2 solver 已正式執行完成。

## 7. 可重現證據

輸入：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_partition_v2`
- 其 `strict_regions` symlink 指向 `hcc1395_strict_regions_v2`
- 其 `extraction` symlink 指向 `20260722.../hcc1395_chr1_22_direct_big7_v2`

命令：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  research/20260723_strict_two_ssnv_phaseblock_validation/agents/audit_tree_input_semantics.py \
  --partition-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_partition_v2 \
  --min-read 3
```

實際輸出重點：

```text
blocks: total=11712, k>=2=11700, tree_input=11590
tree input: multisite-pattern=11207, only-single-fixed=383
tree input: coverage-all=11282, pattern-connected=9717, disconnected=1873
strict edges: total=76202, cross-block=35008
within-block strict edges=41194, solver-visible=30344, not-visible=10850
source W: any strict-edge loss=2334, zero solver-visible strict edge=411
all invariants=true
```

Targeted adapter tests：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_partition_to_mlhp.py
```

實際輸出：

```text
4 passed in 0.57s
```

現有 tests 驗證 exact-PS routing、MINREAD 與 k=1 exclusion，但沒有涵蓋：

- only-single-fixed supported patterns；
- pattern graph disconnected；
- strict edge 到 solver-visible pair 的守恆。

## 8. 最終 verdict

**Needs revision before thesis-level topology claim。**

- Region 建立的 two-sSNV pairwise evidence 本身正確。
- Partition 與 adapter 也各自符合目前程式契約。
- 但兩層 threshold 的 aggregation grain 不同，且 current adapter 未對 pattern coverage／connectivity 設 gate。
- 在 representation 修補前，只有 5,009 blocks 適合放入「pattern-connected multi-mutation topology」primary denominator；其餘必須以 control、partial candidate 或 abstain 分層，不可全稱完成的 clone/subclone tree。
