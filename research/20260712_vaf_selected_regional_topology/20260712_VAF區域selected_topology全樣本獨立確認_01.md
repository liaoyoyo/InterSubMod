<!--
建立時間: 2026-07-12 02:20 +0800
目標: 獨立重跑並區分 family-specific read-AF/VAF 可指定的 mutation-labeled exact tree 與 rooted-unlabeled topology shape，以及全 7 datasets 的比例
處理範圍: chr1-22; 7 dataset rows; 6 biological samples; historical layered-v2 engineering snapshot
關聯檔案:
  - research/20260711_read_group_C_tree_T_topology_report/pre-decision-audit.md
  - research/20260711_read_group_C_tree_T_topology_report/20260711_VAF_top050並列樹與拓撲完整重算_01.md
  - research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_census.json
  - research/20260712_vaf_selected_regional_topology/data/20260712_vaf_selected_topology_by_sample.tsv
-->

# VAF 區域 selected tree／topology shape：全樣本獨立確認

> **PARTIAL SCIENTIFIC VALIDATION / PUBLICATION NO-GO**
>
> **TL;DR**：可新增描述性 selection，但必須分兩層。若只問 rooted-unlabeled shape，結構先行再以 VAF 解原 Topo>1 後為 37,039/39,885（92.86%）；若要分辨 mutation-labeled 正／反祖先 exact tree，較嚴格 coverage 是 29,924/39,885（75.03%）。兩者都不是經校準的最可能生物樹或獨立確認。

## 1. Task type、目標與研究問題

- Task type：**B Comprehensive validation**。
- Scope：chr1-22 × 7 dataset rows／6 biological samples；HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 dataset rows。
- 服務目標：G3（read-level evidence）、G4（跨樣本一致性）、G5（可外部稽核）。
- 研究問題：在保留所有 read-compatible candidate trees 的前提下，family-specific read-AF 的 exact ordering-score 第一順位集合，能否為 region 指定單一 mutation-labeled exact tree，或至少指定單一 rooted-unlabeled topology shape？

## 2. 假設、成功條件與失敗條件

### 假設

若 region 內每個 primary HP unit 的最高 exact ordering-score candidate set 只有一個 ordered forest product，則可標為 `unique_first exact tree`；若有多棵 exact trees、但都映射到同一 canonical shape，則只可標為 `selected_one_shape`。若最高分集合跨越多個 shapes，則必須保留為 `tied_different_topology`。

### 成功條件

1. 全 7 dataset rows 都重新枚舉全部 non-capped candidate trees，candidate/shape mismatch=0。
2. Region、unit、summary、checks 與既有 frozen output byte-identical。
3. JSON 排除 `generated_at` 後 semantic-identical；candidate TSV 解壓後 content-identical。
4. `39,885 = 10,832 + 29,053`、`29,053 = 28,183 + 742 + 128` 等分母守恆。
5. 逐樣本五類結果能各自回到該樣本全部 ambiguous complete regions。

### 失敗條件

- 任一 sample candidate/shape mismatch、region key duplicate、分母不守恆或 exact tie 對 `1e-12` 至 `1e-6` tolerance 翻類。
- 把 softmax relative weight 稱為 calibrated posterior，或把 selected tree／shape 稱為 confirmed biological topology。

## 3. 定義

- `unique_first`：最高 exact score 只有一棵 candidate tree。
- `tied_first_same_topology`：最高分有多棵 exact trees，但 canonical topology shape 相同。
- `tied_first_different_topology`：最高分集合跨越多個 topology shapes，仍無法選定區域 Topo。
- `canonical_shape()` 會移除 mutation/node labels，只留下 rooted、directed、unordered branching skeleton；因此同一 shape 仍可能包含不同的 mutation ancestor direction。
- Exact-tree coverage：`structural T=1 + VAF unique_first`。
- Shape coverage：先接受所有原始 `Topo=1`，再只對原始 `Topo>1` 計入 VAF 第一順位集合為同一 shape 者。

因此，`selected shape` 不要求 exact tree 唯一；若要回答「A 是 B 的祖先，還是 B 是 A 的祖先」，必須採 mutation-labeled `unique_first exact tree` 口徑。兩種 selection 都不刪除原候選全集。

### 現行排序公式

對 mutation \(i\)，family-specific read AF 為：

\[
f_i=\frac{ALT_i}{REF_i+ALT_i}
\]

候選樹 \(T\) 的分數為：

\[
S(T)=\sum_{(p\rightarrow c)\in E(T)}\sum_{a\in A(p)}(f_a-f_j)
\]

其中 \(j\) 是 edge 新取得的 mutation，\(A(p)\) 是 parent 已含的 mutations；第一順位集合為
\(\mathcal{T}^*=\{T:S(T)=\max_{T'}S(T')\}\)。程式以整數 read counts 轉成 `Fraction`，只有 exact score 完全相等才算並列。

## 4. 輸入、執行命令與輸出

### 輸入

- Run root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2`
- Manifest：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/input_manifest.json`
- Corrected C/T/Topo report：`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json`
- Legacy VAF summary：`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_region_resolution_summary.json`

### 獨立重算命令

```bash
python3 InterSubMod/research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_top_tie_census.py \
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 \
  --input-manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/input_manifest.json \
  --corrected-report InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json \
  --legacy-vaf-summary InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_region_resolution_summary.json \
  --method-script-dir InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts \
  --output-json /tmp/20260712_vaf_topology_recalc/vaf_top_tie_census.json \
  --output-region-tsv /tmp/20260712_vaf_topology_recalc/vaf_top_tie_regions.tsv \
  --output-unit-tsv /tmp/20260712_vaf_topology_recalc/vaf_top_tie_units.tsv \
  --output-candidate-tsv-gz /tmp/20260712_vaf_topology_recalc/vaf_top_tie_candidates.tsv.gz \
  --output-summary-tsv /tmp/20260712_vaf_topology_recalc/vaf_top_tie_summary.tsv \
  --output-checks-tsv /tmp/20260712_vaf_topology_recalc/vaf_top_tie_checks.tsv \
  --output-report-md /tmp/20260712_vaf_topology_recalc/vaf_top_tie_report.md
```

### 實際輸出片段

```text
OUTPUT JSON -> /tmp/20260712_vaf_topology_recalc/vaf_top_tie_census.json
OUTPUT REGION TSV -> /tmp/20260712_vaf_topology_recalc/vaf_top_tie_regions.tsv
OUTPUT UNIT TSV -> /tmp/20260712_vaf_topology_recalc/vaf_top_tie_units.tsv
OUTPUT CANDIDATE TSV.GZ -> /tmp/20260712_vaf_topology_recalc/vaf_top_tie_candidates.tsv.gz
STATUS -> PASS
```

## 5. 全樣本結果

| Dataset row | Complete | 結構-first＋VAF 單一 shape | Shape coverage | Mutation-labeled exact tree 唯一 | Exact-tree coverage | 原 Topo>1 被 VAF 縮至單 shape |
|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | 6,940 | 6,798 | 97.95% | 6,395 | 92.15% | 3,360/3,502 = 95.95% |
| HCC1395_DORADO | 6,750 | 6,082 | 90.10% | 5,572 | 82.55% | 3,638/4,306 = 84.49% |
| COLO829 | 6,949 | 5,600 | 80.59% | 3,294 | 47.40% | 1,724/3,073 = 56.10% |
| H1437 | 6,984 | 6,622 | 94.82% | 4,643 | 66.48% | 2,038/2,400 = 84.92% |
| H2009 | 5,882 | 5,661 | 96.24% | 4,147 | 70.50% | 2,133/2,354 = 90.61% |
| HCC1937 | 2,557 | 2,542 | 99.41% | 2,343 | 91.63% | 601/616 = 97.56% |
| HCC1954 | 3,823 | 3,734 | 97.67% | 3,530 | 92.34% | 1,569/1,658 = 94.63% |
| **7 dataset rows pooled** | **39,885** | **37,039** | **92.86%** | **29,924** | **75.03%** | **15,063/17,909 = 84.11%** |

完整逐欄結果：`InterSubMod/research/20260712_vaf_selected_regional_topology/data/20260712_vaf_selected_topology_by_sample.tsv`。

### Pooled 守恆

```text
39,885 complete
├─ 21,976 structurally Topo=1
│  ├─ 10,832 T=1 / Topo=1
│  └─ 11,144 T>1 / Topo=1
└─ 17,909 structurally Topo>1
   ├─ 15,063 VAF first set maps to one shape
   ├─  2,205 VAF first set still spans shapes
   └─    641 not evaluable

Final single-shape coverage = 21,976 + 15,063 = 37,039 / 39,885 = 92.86%
Exact-tree coverage = 10,832 + 19,092 VAF unique first = 29,924 / 39,885 = 75.03%
```

另有一個可重現但不建議作主指標的 36,810/39,885（92.29%）：它是 `T=1 + VAF-evaluable T>1 first-set one shape`，會排除 229 個雖無法做 VAF ranking、但結構本來已是 `T>1/Topo=1` 的 regions。正式流程應採上面的結構-first 92.86%。

在 28,183 個可評估 T>1 regions 中，VAF exact unique first 為 19,092/28,183（67.74%）；另有 6,886 個雖 exact trees 並列，但只落在同一 unlabeled shape。

原始 17,909 個 `T>1/Topo>1` regions 中：15,063（84.11%）在 VAF 第一順位集合只剩一種 Topo，2,205 仍跨 Topo，641 不可評估；若只看可評估 17,268 個，為 15,063/17,268（87.23%）。

## 6. 驗證結果

| Artifact | Frozen SHA-256 | Independent rerun | Verdict |
|---|---|---|---|
| region TSV | `a6649bf12dba...` | `a6649bf12dba...` | byte-identical |
| unit TSV | `3629c4c94daa...` | `3629c4c94daa...` | byte-identical |
| summary TSV | `9be8b99aea5b...` | `9be8b99aea5b...` | byte-identical |
| checks TSV | `f48666e17a7b...` | `f48666e17a7b...` | byte-identical |
| JSON excluding `generated_at` | `3bd7fe34c221...` | `3bd7fe34c221...` | semantic-identical |
| candidate TSV decompressed content | `7ebec1d54f4e...` | `7ebec1d54f4e...` | content-identical |

全部既有 aggregate/sample checks PASS；獨立 region group-by 亦重現逐樣本比例。

## 7. 方法判斷與 claim ceiling

### GO

- 方法上可補欄位，但本輪現有 region TSV **尚未 materialize** 實際的 `selected_tree_id` 或 HP-labeled `selected_shape_tuple`；目前只有第一順位樹／shape 數量與分類。
- 建議下一版輸出 `selection_level={exact_tree|unlabeled_shape|unresolved|not_evaluable}`、`n_top_exact`、`n_top_shapes`、`selected_tree_id`、`selected_shape_tuple` 與 `direction_evaluable`。
- 可把 84.11% 報成「VAF 對原 Topo>1 的 single-shape resolving coverage」，92.86% 報成「結構-first＋VAF 的整體 single-shape coverage」。
- 若必須回答正向／反向 mutation ancestry，只能用 75.03% 的 `unique_first exact-tree heuristic coverage`，且必須保留 heuristic 限定詞。

### 不可升格

- 不可寫「92.86% 的正向祖先樹已確認」；unlabeled shape 無法分辨 `A→B` 與 `B→A`。
- 不可把 75.03% 寫成「真實 exact tree 已確認」；它仍只是同一批 read-AF 的 heuristic argmax。
- 不可把 softmax weight 稱為 posterior probability。
- 不可把同一 read-derived evidence 的 ranking 稱為 independent validation。

### 主要限制

1. 21,646/35,457 ranked HP units 的候選樹 ancestry-comparison 數不同，score sum 可能偏向比較項較多的樹。
2. 分數是差值總和；正負 ancestry deltas 可互相抵銷，高總分不保證每一條 edge 都符合方向。
3. 1,738 個 exact-top candidates 沒有 ancestry comparison，direction evidence 不可評估。
4. Raw family read-AF 尚未校正 CN、purity、mutation multiplicity，也沒有 molecule bootstrap。
5. Region 是 HP-labeled ordered forest tuple，不一定是一棵相連的區域樹。
6. COLO829 的 exact-tree coverage 只有 47.40%；在控制 read depth、CN、候選空間前，不解讀成生物差異。
7. 本輪來源仍是 historical layered-v2；clean raw-all layered-v3 完成前維持 publication NO-GO。

## 8. 結論與下一步

**Verdict**：`GO_DESCRIPTIVE_SELECTED_TREE_OR_SHAPE`；`NO_GO_CONFIRMED_TRUE_ANCESTRY`。

下一個可升級步驟不是再調 softmax threshold，而是：

1. 以 normalized score／binomial或 beta-binomial constrained likelihood 重排，處理比較項數與 read-depth uncertainty。
2. 加入 CN、purity、mutation multiplicity 分層或校正。
3. 對 read/molecule bootstrap 報 selected exact tree 與 shape 的穩定率。
4. Materialize 每個 region 的 HP-labeled selected tree／shape tuple，同時保留並列集合與不可評估狀態。
5. 在 clean layered-v3 7/7 aggregate 完成後，以同一 script 全量重生；再與 historical shape 92.86%／exact-tree 75.03% 比較。
