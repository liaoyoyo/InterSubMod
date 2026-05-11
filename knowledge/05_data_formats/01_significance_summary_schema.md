---
id: ism-kb-05-data-formats-significance-summary-schema
name: "significance_summary.csv Schema (59 欄)"
description: "ISM 主要 CSV 輸出的 59 欄完整字典分組摘要：區域識別、統計檢驗、HP/Allele 驗證、品質評估、判定。權威全欄定義在 docs/data_specs/。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "column groupings against docs/data_specs/20260411_significance_summary欄位字典_01.md"
related_ids:
  - ism-kb-05-data-formats-index
  - ism-kb-05-data-formats-master-dataset-schema
  - ism-kb-04-parameters-statistical-methods
  - ism-kb-05-data-formats-per-region-files
  - ism-kb-07-derived-features-hpfinengroups
  - ism-kb-07-derived-features-coverage-multiple
tags: [data-format, csv, schema, significance-summary]
canonical_paths: [05_data_formats/01_significance_summary_schema.md]
alias_paths: []
---

# significance_summary.csv Schema (59 欄)

- 一句結論：ISM C++ 產出的主 CSV，59 欄分 13 群；權威全欄字典在 `docs/data_specs/20260411_significance_summary欄位字典_01.md`
- 適用對象：讀 CSV 做分析、寫後處理 script
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  head -1 /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/*/significance_summary.csv \
    | tr ',' '\n' | nl | head -60
  ```

---

## 生成位置

- **產出函式**：`src/core/RegionProcessor.cpp:1066` (`write_significance_summary()`)
- **輸出路徑**：`<output_dir>/significance_summary.csv`
- **每行代表**：一個 SNV region 的分析結果
- **總欄數**：59

---

## 欄位分群（13 群）

### A. 區域識別（5 欄）
`RegionID`, `Chr`, `Pos`, `Ref`, `Alt`

### B. 基本統計（2 欄）
`NumReads`, `NumCpGs`

### C. 全域檢驗（4 欄）
`GlobalP`, `CramersV`, `HeuristicScore`, `PassedGating`

### D. 距離度量摘要（2 欄）
`PairwiseMeanDist`, `PairwiseMedianDist`

### E. 分群 PERMANOVA（5 欄）
`ClusterPermanovaF`, `ClusterPermanovaP`, `ClusterPermanovaValid`, `ClusterDispersionP`, `ClusterDispersionWarn`

### F. Stage 1: HP Merged（5 欄）
`HPMergedDelta`, `HPMergedP`, `HPMergedSig`, `HP1FamilyN`, `HP2FamilyN`

### G. Stage 2: HP Fine-Grained（4 欄）
`HPFineF`, `HPFineP`, `HPFineSig`, `HPFineNGroups`

### H. Stage 3: Allele（3 欄）
`AlleleDelta`, `AlleleP`, `AlleleSig`

### I. Label × HP/Allele PERMANOVA（10 欄）
`LabelHPPermanovaF/P/Valid`, `LabelHPDispersionP/Warn`, `LabelAllelePermanovaF/P/Valid`, `LabelAlleleDispersionP/Warn`

### J. Stage 4: Unassigned Affinity（5 欄）
`UnassignedAffinity`, `UnassignedAffinityP`, `UnassignedAffinityDir`, `NHP3`, `NHP0`

### K. 判定（3 欄）
`DominantLabel`, `Stability`, `VerificationClass`

### L. 品質評估（7 欄）
`HP_Ratio`, `Potential_LOH`, `Coverage_Multiple`, `Coverage_Category`, `LOH_Subtype`, `Quality_Score`, `Quality_Tier`

### M. 局部檢驗（4 欄）
`LocalBestCluster`, `LocalBestP`, `Significant`, `SuggestFilter`

---

## 重點欄位速查

| 欄位 | 類型 | 用途 | 值域 |
|------|------|------|------|
| `HPFineNGroups` | int | Fine-grained 有效群組數（subclone marker） | 0-4 |
| `HPMergedDelta` | float | HP1-family vs HP2-family 距離差 | 可正可負 |
| `HeuristicScore` | float | 整合顯著性評分 | ≥0 |
| `PairwiseMeanDist` | float | 所有 read pair 的平均距離 | [0, 1] for BERNOULLI |
| `Coverage_Multiple` | float | CN proxy (NumReads / 期望覆蓋度) | >0 |
| `HP_Ratio` | float | HP1-family / effective_hp_reads（Laplace smoothing） | [0, 1] |
| `Potential_LOH` | bool/int | LOH 候選標記 | 0/1 |
| `Significant` | bool | 局部檢驗是否顯著 | 0/1 |

---

## 與統計方法對應

59 欄中，**30 欄**由 [../04_parameters/03_statistical_methods.md](../04_parameters/03_statistical_methods.md) 的方法產生：
- C, E, F, G, H, I, J 群 → PERMANOVA / FFH / 多階段 HP 驗證

---

## 權威全欄字典

本 KB 提供分組摘要；若需每一欄精確定義（dtype、null 處理、計算公式）→ 讀：

👉 **`docs/data_specs/20260411_significance_summary欄位字典_01.md`**

---

## 常用後處理操作

```bash
# 取全欄數
head -1 significance_summary.csv | awk -F, '{print NF}'
# 預期 59

# 篩選 Quality_Tier = "Tier A" 的 regions
awk -F, 'NR==1 || $(NF-2)=="Tier A"' significance_summary.csv

# Python pandas 載入
python3 -c "import pandas as pd; df = pd.read_csv('significance_summary.csv'); print(df.shape)"
# 預期 (N_regions, 59)
```

---

## 相關

- Master dataset（Python 聚合）：[02_master_dataset_schema.md](02_master_dataset_schema.md)
- Per-region 子檔：[03_per_region_files.md](03_per_region_files.md)
- 統計方法：[../04_parameters/03_statistical_methods.md](../04_parameters/03_statistical_methods.md)
- 權威字典：[../../docs/data_specs/20260411_significance_summary欄位字典_01.md](../../docs/data_specs/20260411_significance_summary欄位字典_01.md)
