---
title: 數字口徑定案 — 解釋頁建置前的反捏造前置（§13.0）
doc_type: provenance_resolution
created: 2026-06-12
status: resolved
note: 盤點抓出的同主題不同口徑數字，建頁前以 ground truth 釘死單一口徑。每筆附驗證指令。
---

# 數字口徑定案

> 規則：解釋頁引用這些數字時，用本表定案值 + source_path；不並列舊口徑，不憑記憶。

## ① significance_summary.csv 欄數 = **59**（非 117）

| 證據 | 值 | 來源 |
|---|---|---|
| 實際輸出檔 header 計數（ground truth）| 59 | `head -1 /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv \| tr ',' '\n' \| wc -l` → 59 |
| 權威欄位字典 | 59 欄 | `InterSubMod/docs/data_specs/20260411_significance_summary欄位字典_01.md`（標題「全部 59 欄」）|
| KB schema | 59 欄 | `InterSubMod/knowledge/05_data_formats/01_significance_summary_schema.md` |
| 盤點宣稱的 117 | **查無** | grep "117" 於 `method_comparison/.../01_ism_method_spec_from_source.md` → 0 命中；判定為合成誤差 |

**定案：59 欄。** 解釋頁一律寫 59。

## ② PERMANOVA permutation 次數 = 生產 **99**（函式庫預設 999）

| 層級 | 值 | 來源（file:line）|
|---|---|---|
| 函式庫 config 預設 | 999 | `include/core/StructureTest.hpp:26`、`include/core/LabelTest.hpp:31`（`int n_permutations = 999;`）|
| **生產 pipeline 覆寫**（實跑值）| **99** | `src/core/RegionProcessor.cpp:1573`（`sig_config.structure_config.n_permutations = 99;`）|

**定案：PERMANOVA（StructureTest）實跑 = 99**；p 值解析度下限 ≈ 1/(99+1) = 0.01。
解釋頁敘述：「函式庫預設 999，但生產設定覆寫為 99（RegionProcessor.cpp:1573）」——兩者都標，避免再被誤讀為矛盾。

## ③ 其他（第一批外，但名詞卡 example 可能引）

- **BRCA2 HP-axis Δβ = −0.122**（修正 max-collapse；−0.054 是 buggy 雙列砍半 artifact，禁用）。來源見盤點線 6 key_data。
- **phasing n**：樣本 instance = 7（含 1 個 basecall 變體 HCC1395_DORADO）；**生物 n = 6**。統計用 bio-n=6，標明。
