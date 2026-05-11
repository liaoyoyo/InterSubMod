---
id: ism-kb-04-parameters-statistical-methods
name: "Statistical Methods 規格"
description: "ISM 的統計檢驗：Fisher-Freeman-Halton (FFH)、PERMANOVA、PERMDISP2、Cramér's V；含 permutation 次數、閾值、實作位置、多階段 HP 驗證協定。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "statistical methods against include/core/{GlobalTest,StructureTest,LabelTest}.hpp"
related_ids:
  - ism-kb-04-parameters-index
  - ism-kb-04-parameters-distance-metrics
  - ism-kb-05-data-formats-significance-summary-schema
  - ism-kb-04-parameters-config-defaults
  - ism-kb-07-derived-features-hpfinengroups
tags: [parameters, statistics, PERMANOVA, FFH, Cramer-V, significance]
canonical_paths: [04_parameters/03_statistical_methods.md]
alias_paths: []
---

# Statistical Methods 規格

- 一句結論：PERMANOVA 999 permutations；Cramér's V threshold 0.1；多階段 HP 驗證（merged → fine-grained → allele → unassigned）
- 適用對象：理解 significance_summary.csv 欄位者、修改統計流程者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -n "kNumPermutations\|DEFAULT_" /big7_disk/liaoyoyo2001/InterSubMod/include/core/Config.hpp | head -20
  ```

---

## 統計方法速查

| 方法 | 用途 | 實作檔 |
|------|------|--------|
| **FFH (Fisher-Freeman-Halton)** | 列聯表獨立性精確檢定 | `include/core/GlobalTest.hpp` |
| **Cramér's V** | 列聯表效應量 | `include/core/MathUtils.hpp:156-170` |
| **PERMANOVA** | 距離矩陣分群結構檢定 | `include/core/StructureTest.hpp` |
| **PERMDISP2** | 群間離散度同質性檢定 | `src/core/StructureTest.cpp` |
| **Label-First PERMANOVA** | HP/Allele 標籤驗證 | `include/core/LabelTest.hpp` |
| **啟發式評分** | 整合多維度信號用於粗篩 | `src/core/RegionProcessor.cpp` |

---

## 1. Fisher-Freeman-Halton (FFH)

**用途**：2 維以上列聯表的獨立性精確檢定（甲基化狀態 × read 標籤）

**實作**：`include/core/GlobalTest.hpp` + `src/core/GlobalTest.cpp`

**輸出欄位**：
- `GlobalP`：p-value（[0, 1]，越小越顯著）
- `CramersV`：效應量（[0, 1]）
- `PassedGating`：是否通過 gating（預設 p ≤ 0.05 或其他閾值）

**限制**：
- 2×2 表 Cramér's V 有 93% 為零的已知問題（見 `09_conclusions/03_concluded_negative.md` 的 Feature Design R1-R5）

---

## 2. Cramér's V 效應量

**公式**：
```
V = sqrt( chi² / (n × (min(r, c) - 1)) )
```
其中 `r, c` 為列聯表維度，`n` 為總觀察數。

**閾值**：`include/core/GlobalTest.hpp:32` 定義 `min_cramers_v = 0.1`（視為有意義效應）。注意：**Config.hpp 無 `kCramerVThreshold` 常數**，此閾值定義於 `GlobalTest::Config` 結構體。

**取值**：[0, 1]，越大效應越強

---

## 3. PERMANOVA（Permutational MANOVA）

**目的**：檢驗距離矩陣中的分群結構是否源自標籤

**核心統計 (Pseudo-F)**：
```
F = (SS_between / df_between) / (SS_within / df_within)
```
`df_between = k - 1`, `df_within = n - k`（k 群數，n 總樣本數）

**Permutation test**：
- **預設 999 permutations**（`StructureTest.hpp:26` / `LabelTest.hpp:31` 的 struct default）
- ⚠️ **Runtime override**：`src/core/RegionProcessor.cpp:1573` 覆蓋為 **99**，實際 pipeline 執行時只跑 99 次 permutations（p-value 精度上限 0.01）
- 隨機重排 labels，累計 pseudo-F ≥ 觀測值的次數 → p-value

**有效性要求**：
- ≥2 群
- 每群 ≥2 reads
- 否則 `Valid=false`

**實作**：
- SS 計算：`src/core/StructureTest.cpp:140-180`（避免數值消失）
- Permutation loop：同檔案

**輸出欄位**：
- `ClusterPermanovaF`, `ClusterPermanovaP`, `ClusterPermanovaValid`
- `LabelHPPermanovaF/P/Valid`, `LabelAllelePermanovaF/P/Valid`

---

## 4. PERMDISP2（群間離散度檢定）

**目的**：檢驗 PERMANOVA 的前提假設（群間變異度均勻）

**方法**：計算各 reads 到其群中心（centroid）的距離，用 Levene test

**輸出**：
- `ClusterDispersionP`, `ClusterDispersionWarn`
- `LabelHPDispersionP/Warn`, `LabelAlleleDispersionP/Warn`

**警告**：若 p ≤ 0.05 → 離散度不均，PERMANOVA F 可能膨脹（需謹慎解讀）

---

## 5. 多階段 HP 驗證

在 PERMANOVA 之前/之外，做多階段標籤意義驗證：

### Stage 1: HP Merged
合併 `(HP1 + HP1-somatic)` vs `(HP2 + HP2-somatic)`

**檢定**：permutation test Δ = d_between - d_within（999 次）

**條件通過**：`HPMergedP ≤ 0.05 AND HPMergedDelta > 0`

**欄位**：`HPMergedDelta`, `HPMergedP`, `HPMergedSig`

### Stage 2: HP Fine-Grained
四群 PERMANOVA：HP1-G（germline HP1）, HP1-S（somatic HP1）, HP2-G, HP2-S（各需 ≥3 reads）

**關鍵欄位**：`HPFineNGroups`（有效群組數，用於篩選 subclone markers）

**結論**：≥4 群 + NumReads≥80 → TP rate 89.1% / 92.81%（詳見 [../07_derived_features/01_hpfinengroups.md](../07_derived_features/01_hpfinengroups.md)）

### Stage 3: Allele
ALT-supporting vs REF-supporting reads

**注意**：與 caller AF 有 confound，非獨立信號

**欄位**：`AlleleDelta`, `AlleleP`, `AlleleSig`

### Stage 4: Unassigned Affinity
HP3 / HP0（未分配 reads）的歸屬親和度

**欄位**：`UnassignedAffinity`, `UnassignedAffinityP`, `UnassignedAffinityDir`, `NHP3`, `NHP0`

---

## 6. 啟發式評分（Heuristic Score）

**公式**：
```
HeuristicScore = -log10(best_p) + best_v × 2
```

**調整**：
- 若有 dispersion warning → ×0.5
- 若 PERMANOVA 不顯著 → ×0.7

**用途**：整合多維度信號，用於粗篩 regions

---

## 統計輸出欄位對應

59 欄 `significance_summary.csv` 中，統計相關約 30 欄，主要群組：
- **C. 全域檢驗（4 欄）**：GlobalP, CramersV, HeuristicScore, PassedGating
- **E. Cluster PERMANOVA（5 欄）**
- **F. Stage 1 HP Merged（5 欄）**
- **G. Stage 2 Fine-Grained（4 欄）**
- **H. Stage 3 Allele（3 欄）**
- **I. Label × HP/Allele PERMANOVA（10 欄）**
- **J. Stage 4 Unassigned（5 欄）**

完整字典：[../05_data_formats/01_significance_summary_schema.md](../05_data_formats/01_significance_summary_schema.md)

---

## 預設參數速查

| 參數 | 值 | 出處 |
|------|---|------|
| PERMANOVA permutations | 999 (struct default) / **99 (runtime override)** | StructureTest.hpp:26 / RegionProcessor.cpp:1573 |
| Cramér's V threshold | 0.1 | GlobalTest.hpp:32 (`min_cramers_v`) |
| Min reads per group (fine-grained) | 3 | RegionProcessor.cpp |
| Min reads for clustering | 10 | Config.hpp:67 |
| Linkage method | UPGMA | Config.hpp:66 |

---

## 常見陷阱

- ❌ **誤把 Cramér's V=0 當無信號**：2×2 表 93% 為零是結構性限制，不是真無信號
- ❌ **忽略 PERMDISP warning 直接用 PERMANOVA F**：離散度不均會膨脹 F
- ❌ **把 Allele test 當獨立信號**：與 caller AF confound

---

## 相關

- 距離度量（PERMANOVA 輸入）：[02_distance_metrics.md](02_distance_metrics.md)
- 輸出欄位：[../05_data_formats/01_significance_summary_schema.md](../05_data_formats/01_significance_summary_schema.md)
- HPFineNGroups 結論：[../07_derived_features/01_hpfinengroups.md](../07_derived_features/01_hpfinengroups.md)
- 原始碼：`include/core/{GlobalTest,StructureTest,LabelTest}.hpp`
