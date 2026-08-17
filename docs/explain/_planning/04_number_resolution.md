---
title: 數字口徑定案 — 解釋頁建置前的反捏造前置（§13.0）
doc_type: provenance_resolution
created: 2026-06-12
status: resolved_revision_scoped
note: 盤點抓出的同主題不同口徑數字，以 revision-scoped ground truth 定案；historical 與 current 不可混用。每筆附驗證指令。
---

# 數字口徑定案

> 規則：解釋頁引用這些數字時，必須同時標 revision、定案值與 source path；historical 口徑只能作歷史說明，不能當 current contract。

## ① significance_summary.csv 欄數：historical **59**；frozen baseline audited **199**

| 證據 | 值 | 來源 |
|---|---|---|
| historical canonical run（舊 schema）| 59 | `head -1 /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv \| tr ',' '\n' \| wc -l` → 59 |
| historical 欄位字典 | 59 欄 | `InterSubMod/docs/data_specs/20260411_significance_summary欄位字典_01.md`（標題「全部 59 欄」）|
| historical KB snapshot | 59 欄 | `InterSubMod/knowledge/05_data_formats/01_significance_summary_schema.md` |
| **frozen release baseline `ddd8909a`** | **199 欄** | `InterSubMod/include/core/SignificanceCsvColumns.hpp` 與 baseline synthetic E2E schema receipt；`73afaeac-dirty` 僅是 historical byte-equivalent `IN_PROGRESS/PARTIAL` audit，不是 release source |
| 盤點宣稱的 117 | **查無** | grep "117" 於 `method_comparison/.../01_ism_method_spec_from_source.md` → 0 命中；判定為合成誤差 |

**定案：**59 欄只描述 historical/legacy output；本次 frozen `ddd8909a` research contract 是 **199 欄**。公開解釋頁不得把 59 寫成 frozen baseline schema，也不得把本次 snapshot 稱為 production release。

## ② PERMANOVA permutation 次數：historical **99**；frozen baseline research path **999**

| 層級 | 值 | 來源（file:line）|
|---|---|---|
| historical revision `9098f11` 的 StructureTest production override | 99 | `git show 9098f11:src/core/RegionProcessor.cpp`（historical `n_permutations = 99`）|
| historical `9098f11` label-first/default | 999 | `include/core/LabelTest.hpp`／`include/core/StructureTest.hpp` 的 config default |
| **frozen release baseline `ddd8909a` research paths** | **999** | `src/core/RegionProcessor.cpp` 的三個 `n_permutations = 999` 路徑，以及 `include/core/{StructureTest,LabelTest}.hpp`；`73afaeac-dirty` 排除於本次 release source |

**定案：**99 permutations／p-floor 0.01 是 revision-scoped historical StructureTest 行為；本次 frozen `ddd8909a` research contract 是 **999 permutations**，分母 `999+1`，p-floor **0.001**。公開頁若談本次 baseline 只能寫後者，不可將 `73afaeac-dirty` 升格為 release source。

## ③ 其他（第一批外，但名詞卡 example 可能引）

- **BRCA2 HP-axis Δβ = −0.122**（修正 max-collapse；−0.054 是 buggy 雙列砍半 artifact，禁用）。來源見盤點線 6 key_data。
- **phasing n**：樣本 instance = 7（含 1 個 basecall 變體 HCC1395_DORADO）；**生物 n = 6**。統計用 bio-n=6，標明。
