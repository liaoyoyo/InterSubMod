<!--
建立時間: 2026-04-01 12:00
目標: LOH Evidence Panel 四輪研究審閱文件 — 完整方法、動機、關鍵發現與圖表彙整
處理範圍: Round 1-4 全量 LOH 盤點，7 paired + 7 TO（14 datasets, 748,391 rows）
關聯檔案:
  - docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round2_support_hp0_analysis_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round3_methyl_hp0_filter_01.md
  - docs/plans/2026/03/20260326_LOH盤點執行規格_01.md
  - docs/plans/2026/03/20260327_LOH_round2_execution_spec_01.md
-->

# LOH Evidence Panel 四輪研究審閱

**審閱日期**：2026-04-01
**涵蓋期間**：2026-03-26 -- 2026-03-28（四輪完整執行）
**資料規模**：14 datasets（7 paired + 7 TO），748,391 region rows
**樣本清單**：HCC1395_HKU_5kHz, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954
**核心問題**：LOH（Loss of Heterozygosity）相關特徵是否能作為 ClairS somatic variant calling 的 FP 鑑別或過濾指標？

> **勘誤說明**：Round 1-3 報告中所有 TO 端 HP 相關統計曾受 HP integer tag bug 影響（ReadParser.cpp 中 HP:i:11/21/33 mapping 錯誤），已於 2026-03-30 全面修正。本文引用的 TO 數字均為修正後版本。Paired 端結論不受此 bug 影響。

---

## 研究脈絡與四輪邏輯鏈

四輪研究以逐步深化的方式回答「LOH 是否可用」：

```
Round 1: 跨樣本 LOH 診斷基線（分佈盤點）
    → 發現 HP_Ratio=0.5 假象、paired FP enrichment 1.194x
    → 需要 support quality 分層才能下結論

Round 2: effective_hp Support Tier 分層
    → Tier A(>=30) paired enrichment 1.169x (p=7.2e-7) 是訊號來源
    → HP0 在 paired/TO LOH-like 中方向相反，揭示機制差異

Round 3: HP0 Filter + Methylation 聯合 + Tier 閾值敏感度
    → HP0 filter 假設被否定
    → LOH+HPMergedSig = 7.4x FP（待深度驗證）
    → A>=50 enrichment 最強 (2.02x)
    → 所有 filter F1 delta < 0

Round 4: 基線修正 + 深度驗證
    → F1 基線修正（含 FN 的真實值）
    → LOH+HPMergedSig 7.4x 是 HCC1395 chr8 特異性
    → Tier A(30-49) 與 A+(>=50) 方向相反
```

---

## Round 1：跨樣本 LOH 診斷基線

### 目的

確認 LOH-like 現象在 TP/FP 中的分佈基線。這一輪不是要直接證明 LOH 能提升 F1，而是先盤清楚三件事：

1. LOH-like / HP family imbalance 在 7 個樣本的 paired / TO 中，實際分佈長什麼樣子
2. 這些 LOH-like 現象在 TP / FP 之間，是否有穩定且可泛化的差異
3. Paired 與 TO 在同樣本同位點上，是互相支持還是呈現大量不一致

**分析腳本**：`scripts/analysis/build_loh_round1_cross_sample_audit.py`
**輸出 workspace**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/`

### 方法

#### LOH-like 定義

本輪不直接使用原始 `HP_Ratio`，而是重算更精確的指標：

1. `effective_hp_reads = HP1FamilyN + HP2FamilyN` -- 排除 HP0（unphased）和 HP3（ambiguous）reads
2. `hp_ratio_core = HP1FamilyN / effective_hp_reads` -- 純 phased reads 的比例
3. `core_loh_like = (effective_hp_reads > 0) AND (hp_ratio_core < 0.1 OR hp_ratio_core > 0.9)`

其中 HP1FamilyN = HP1 + HP1-1 reads（HP1 family），HP2FamilyN = HP2 + HP2-1 reads（HP2 family）。這與 C++ LabelTest Stage 1 的 merged HP test 一致。

#### 為什麼用 0.1/0.9 閾值

閾值 0.1/0.9 的意義：一個 haplotype 的 reads 佔 90% 以上，代表幾乎只有一個 allele 在該 region 有 phased reads。這是 LOH 的直接表現 -- 一個等位基因的喪失導致 reads 集中在剩餘的那個等位基因上。沿用 InterSubMod 現有 `Potential_LOH` 定義，但每輪另做 exploratory bins（`<0.05, 0.05-0.1, 0.1-0.2, 0.2-0.8, 0.8-0.9, 0.9-0.95, >0.95`）以驗證閾值是否需要調整。

（定義來源：`docs/plans/2026/03/20260326_LOH盤點執行規格_01.md` 第 4.6 節）

#### HP_Ratio 計算公式

```
hp_ratio_core = HP1FamilyN / (HP1FamilyN + HP2FamilyN)
```

- 當 `HP1FamilyN + HP2FamilyN = 0` 時，`hp_ratio_core` 設為 NaN（不可解讀）
- HP0（unphased）和 HP3（ambiguous）reads **不納入分母**，避免 phase completeness 問題被誤解為生物訊號

#### 資料範圍

Round 1 共納入 14 個 dataset，產出 748,391 個 region rows 與 459,782 個 same-locus paired-vs-TO union rows。

（row 來源：`all_region_rows.tsv.gz`，路徑見 workspace）

### 關鍵發現

#### 1. HP_Ratio=0.5 陷阱：69,807 個 region effective_hp=0

在所有 748,391 rows 中，有 69,807 個 region 的 `effective_hp_reads = 0`（完全沒有 HP1/HP2 family reads）。這些 region 的原始 `HP_Ratio` 全部是 0.5，且 `Potential_LOH` 全部是 False。

**意義**：`HP_Ratio = 0.5` 並不一定代表「haplotype 平衡」，它可能只是「根本沒有有效 HP1/HP2 read」。後續所有 LOH 分析必須先排除 `effective_hp_reads = 0` 的 region（後來定義為 Tier C0）。

（數據來源：`hp_coverage_qc_summary.tsv`）

#### 2. Paired FP enrichment 1.194x

| Mode | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment |
|------|---------|---------|-------------|-------------|-----------------|
| Paired | 325,270 | 3,429 | 29.33% | 35.02% | **1.194x** |
| TO（修正後）| 291,310 | 128,382 | 44.5% | 35.8% | **0.805x** |

Paired 端有輕度 FP enrichment（LOH-like 在 FP 中比 TP 更常見）。TO 端修正後方向相反：LOH-like 在 TP 比 FP 更常見（TP 富集）。

（數據來源：`loh_enrichment_by_sample_mode.tsv`）

#### 3. Paired 的 LOH-like 具有 sample heterogeneity

| Sample | Paired FP/TP Enrichment |
|--------|------------------------|
| HCC1954 | 3.185x |
| H2009 | 2.685x |
| H1437 | 1.795x |
| HCC1937 | 1.505x |
| HCC1395_DORADO | 1.260x |
| COLO829 | 1.155x |
| HCC1395_HKU_5kHz | 1.016x |

有些樣本明顯偏 FP（HCC1954 3.2x），有些幾乎沒有區分力（HCC1395 1.0x）。LOH-like 不能直接寫成單一全域規則。

（數據來源：`loh_enrichment_by_sample_mode.tsv`）

#### 4. TO-only FP 是 paired-vs-TO discordance 的主因

| Concordance | Count | Fraction |
|-------------|-------|----------|
| both_tp | 287,092 | 62.44% |
| to_only_fp | 126,865 | 27.59% |
| paired_only_tp | 38,178 | 8.30% |
| to_only_tp | 4,218 | 0.92% |
| paired_only_fp | 1,912 | 0.42% |
| both_fp | 1,517 | 0.33% |

真正支配 discordance 的是大量 TO-only FP，而非 both_fp。這支持把 TO 視為獨立主線問題。

（數據來源：`same_locus_compare.tsv`）

### 圖表

**Fig01 -- LOH-like Fraction Overview**（各樣本 paired/TO LOH-like 佔比概覽）

![Fig01 LOH-like Fraction Overview](figures/fig01_loh_like_fraction_overview.png)

**Fig02 -- HP Ratio Core Distribution**（hp_ratio_core 分佈，TP vs FP，paired vs TO）

![Fig02 HP Ratio Core Distribution](figures/fig02_hp_ratio_core_distribution.png)

---

## Round 2：Support Quality 分層

### 目的

Round 1 確認了 paired 端的 LOH-like FP enrichment，但 enrichment 包含了大量 `effective_hp_reads` 極低的 region（統計力不足的雜訊）。Round 2 的核心問題是：**在什麼品質條件下，LOH evidence 才是可信的 evidence panel 成員？**

具體假說：排除 `effective_hp_reads < 10` 的雜訊後，Tier A（>=30 reads）的 LOH-like FP/TP enrichment 是否更穩定、更具泛化力？

**分析腳本**：`scripts/analysis/build_loh_round2_support_hp0_analysis.py`
**輸出 workspace**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round2_support_hp0_analysis/`

### 方法

#### Tier 定義

| Tier | effective_hp_reads | 說明 | 設計動機 |
|------|--------------------|------|---------|
| C0 | = 0 | 完全無 HP reads | HP_Ratio=0.5 假象，Round 1 已確認 |
| C | 1-9 | 極弱 support | 1-9 reads 不足以區分 LOH vs 隨機波動 |
| B | 10-29 | 中等 support | 統計力有限，需小心解讀 |
| A | >=30 | 強 support | phasing 統計有足夠 reads 支撐 LOH 判斷 |

**為什麼這樣分 Tier**：Tier 邊界基於 phasing 統計力。在二項分佈下，當 `effective_hp_reads >= 30` 且觀察到 `hp_ratio_core > 0.9` 時，在 null hypothesis（真實比例=0.5）下的 p-value < 1e-4，代表偏離平衡具有統計顯著性。低於 10 reads 時，即使觀察到 `hp_ratio_core = 1.0`，也可能僅是取樣噪聲。

Enrichment 計算：`Enrichment = FP_LOH% / TP_LOH%`。數值 >1 表示 LOH-like 在 FP 中更常見（FP 富集），<1 表示在 TP 中更常見（TP 富集）。

（Tier 定義來源：`docs/plans/2026/03/20260327_LOH_round2_execution_spec_01.md` 第二節）

### 關鍵發現

#### 1. Tier A (>=30) paired enrichment 1.169x (p=7.2e-7) -- 真正的信號來源

| Tier | TP Total | FP Total | TP LOH-like% | FP LOH-like% | FP/TP Enrichment | Fisher p |
|------|---------|---------|-------------|-------------|-----------------|---------|
| C0 | 29,873 | 135 | 0.0% | 0.0% | -- | 1.0 |
| C (1-9) | 2,296 | 50 | 67.5% | 70.0% | 1.037x | 0.76 |
| **B (10-29)** | 25,073 | 1,299 | 38.6% | 34.8% | **0.901x** | 0.006 |
| **A (>=30)** | 268,028 | 1,945 | 31.4% | 36.7% | **1.169x** | **7.2e-7** |

Round 1 整體 paired enrichment = 1.194x，對比顯示 Tier A 幾乎是整體 enrichment 的全部來源。

（數據來源：`core1_tier_enrichment_global.tsv`）

#### 2. Tier B 的「反轉」(0.90x) 是 COLO829 主導的 artifact

Tier B 整體 enrichment = 0.901x（FP 中 LOH-like 比 TP 少），看起來是反向訊號。但 per-sample 分析揭示：

- COLO829 貢獻了 Tier B 中 93% 的 FP（1,206/1,299）
- COLO829 本身 Tier B enrichment = 1.093x（正向）
- 但非 COLO829 樣本的 TP LOH-like fraction 很高（70-77%），而它們幾乎沒有 FP

**結論**：Tier B 整體反轉不是真實的生物訊號，是不同樣本間 TP/FP 比例差異導致的 Simpson's paradox（聚合 artifact）。

（數據來源：`core1_tier_enrichment_by_sample_mode.tsv`）

#### 3. Paired LOH-like HP0 低 vs TO LOH-like HP0 高（方向相反）

| 模式 | LOH-like region HP0 | Non-LOH-like HP0 | 解釋 |
|------|---------------------|-------------------|------|
| Paired | **4.1%** | 9.0% | Real LOH：one-sided phasing 把大多 reads 指向一個 family，HP0 自然少 |
| TO（修正後）| **9.6%** | 4.6% | Partial phasing artifact：TO phasing 產生 LOH-like 訊號時伴隨 phasing 不完全 |

**意義**：Paired LOH-like 是「真實 one-sided phasing」的反映（幾乎所有 reads 都被成功分配到某個 HP family，剩餘 HP0 少）。TO LOH-like 有相當比例是 phasing 不完全產生的（HP0 反而更高），揭示了 paired 與 TO 的 LOH 機制根本不同。

（數據來源：`core2_hp0_by_loh_status.tsv`）

#### 4. HP0 的 Tier 依賴性確認 Tier C 不可信

| Mode | Tier | TP HP0 mean | 說明 |
|------|------|------------|------|
| Paired | C | 40.4% | HP0 極高，Tier C LOH 完全不可信 |
| Paired | A | 3.3% | 可接受的 HP0 水平 |
| TO（修正後）| C0 | 85.3% | 幾乎全是 HP0 |
| TO（修正後）| C | 77.7% | 極高 HP0，完全不可信 |
| TO（修正後）| A | 9.2% | 可接受 |
| TO（修正後）| A+ | 2.7% | 最低，phase 支持品質最好 |

修正後 HP0 單調遞減（tier 越高 -> eff_hp 越高 -> HP0 越低），符合物理預期。

（數據來源：`core2_hp0_by_tier.tsv`）

### 圖表

**Fig01 -- Tier Enrichment Global**（paired vs TO 各 Tier 的 FP/TP enrichment 對比）

![Fig01 Tier Enrichment Global](figures/fig01_tier_enrichment_global.png)

**Fig02 -- Paired Tier Enrichment Heatmap**（各樣本 x Tier 的 paired FP/TP enrichment 熱圖）

![Fig02 Paired Tier Enrichment Heatmap](figures/fig02_paired_tier_enrichment_heatmap.png)

---

## Round 3：HP0 Filter + Methylation 聯合

### 目的

Round 2 發現 TO LOH-like 伴隨高 HP0，直覺上「高 HP0 代表品質差」。Round 3 的核心問題是：

1. **HP0 filter 假設**：hp0_ratio 低的 TO LOH-like region 是否 TP% 更高？（即低 HP0 = 更乾淨的 LOH 訊號？）
2. **LOH x methylation 聯合**：LOH=True + HPMergedSig=True 是否構成更強的 FP 指標？
3. **Tier 閾值穩健性**：A>=30 的 enrichment 在不同閾值下是否穩定？

**分析腳本**：`scripts/analysis/build_loh_round3_methyl_hp0_filter.py`
**輸出 workspace**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round3_methyl_hp0_filter/`

### 方法

#### HP0 filter 假設

在 TO Tier A（>=30 reads）的 LOH-like region 中，按 hp0_ratio 分層（Low <5%, Mid 5-15%, High >=15%），比較各層的 TP%。若假設成立，Low HP0 群體的 TP% 應最高。

#### LOH x HPMergedSig 聯合

在 paired Tier A 中，將 region 按 2x2 矩陣分組（LOH=Y/N x HPMergedSig=Y/N），比較各象限的 FP%。HPMergedSig 代表 HP1-family 與 HP2-family 之間的甲基化差異達統計顯著。

#### Tier 閾值敏感度

以 A>=10, 15, 20, 25, 30, 40, 50 七個閾值重跑 enrichment，檢查結論穩定性。

### 關鍵發現

#### 1. HP0 filter 假設否定：High HP0 TP% > Low HP0 TP%

| HP0 Stratum | n（修正後）| TP% | 方向 |
|------------|---------|-----|------|
| Low (<5%) | 114,229 | 74.8% | 基準 |
| Mid (5-15%) | 15,855 | 75.7% | +0.9pp vs Low |
| High (>=15%) | 17,881 | **76.7%** | **+1.9pp vs Low** |

**方向完全相反**：高 HP0 的 region 反而 TP 比例更高。HP0 filter 假設被數據否定。

原因分析：FP 集中在低 HP0 group（尤其 HCC1954 佔 TO FP 最大來源），代表「phasing 看起來乾淨但仍是 FP」的 region 大量存在。FP 的來源並非 phasing artifact，而是其他因素（如 germline LOH region 中的 artifact）。

（數據來源：Round 3 workspace HP0 threshold sweep 結果）

#### 2. LOH+HPMergedSig 聯合 FP% = 5.61% (vs 純 LOH 0.76%) -- 7.4x 差距

Paired Tier A 中的 2x2 矩陣：

| LOH | HPMergedSig | n_TP | n_FP | FP% |
|-----|-------------|------|------|-----|
| **Y** | **Y** | 1,346 | 80 | **5.61%** |
| Y | N | 82,826 | 634 | 0.76% |
| N | Y | 102,234 | 620 | 0.60% |
| N | N | 81,622 | 611 | 0.74% |

LOH=True + HPMergedSig=True 的 FP 率是純 LOH 的 **7.4 倍**。

**生物學解釋**：這些是 allele-specific methylation (ASM) region -- 兩個等位基因的甲基化本來就不同（正常生物學特性）。LOH 使 reads 集中在一個等位基因，加上甲基化差異顯著，容易把 germline SNP 誤判為 somatic。

但此發現需要 Round 4 深度驗證其普遍性。

（數據來源：Round 3 Core 2 joint feature 分析）

#### 3. 所有 filter 變體的 F1 delta < 0

| 樣本 | Filter | rm_TP% | rm_FP% | F1 delta |
|------|--------|--------|--------|---------|
| H2009 | TierA_LOH | 27.64% | 75.58% | **-0.1601** |
| HCC1937 | TierA_LOH | 53.90% | 82.05% | **-0.3623** |
| HCC1395 | TierA_LOH_HPSig | 1.32% | 11.16% | -0.0055 |
| COLO829 | TierA_LOH_HPSig | 0.02% | 0.04% | -0.0001 |

即使 TierA_LOH_HPSig 是最保守的 filter（移除 TP 最少），F1 delta 仍然是負值。根本限制：paired FP 絕對數量太少（8-2,273 個），移除的 TP 損失始終大於 FP 收益。

**LOH 不能作為 binary filter** -- 這個結論在三輪研究中全部一致。

（數據來源：Round 3 Core 4 filter simulation 結果）

#### 4. Tier 閾值敏感度：A>=50 enrichment 最強 (2.018x)

| A 閾值 | Paired Enrichment | p-value |
|--------|------------------|---------|
| >=10 | 1.122 | 2.45e-06 |
| >=25 | 1.062 | 0.041 |
| >=30 | 1.169 | 7.22e-07 |
| >=40 | 1.580 | 3.57e-38 |
| **>=50** | **2.018** | **2.77e-67** |

Paired enrichment 呈非單調結構：10-25 reads 有局部低谷（enrichment 下降），>=30 後跳升，>=50 達最強 2.02x。A=30 是統計顯著的最低合理門檻，但 A=50 有更強鑑別力。

（數據來源：Round 3 Core 3 tier threshold sensitivity 分析）

### 圖表

**Fig05 -- Filter Simulation F1 Delta**（各樣本各 filter 變體的 F1 變化量）

![Fig05 Filter Simulation F1 Delta](figures/fig05_filter_simulation_f1_delta.png)

---

## Round 4：基線修正與深度驗證

### 目的

Round 3 有兩個需要修正或深度驗證的問題：

1. **F1 基線假象**：Round 3 的 filter simulation F1（如 H2009=0.9997）是在「未計入 FN 的子集上」計算的，不是真實 F1
2. **LOH+HPMergedSig 的 7.4x 是否普遍**：需要 per-sample 與 per-chromosome 分解，確認是全樣本現象還是特定樣本驅動

**分析腳本**：`scripts/analysis/build_loh_round4_final_validation.py`
**輸出 workspace**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260328_loh_round4_final_validation/`

### 方法

#### 基線修正

加入完整 FN（False Negative）數據重算 F1：`F1 = 2 * TP / (2 * TP + FP + FN)`。Round 3 的 F1 只在 `TP + FP` 子集上計算（缺少 FN），導致 F1 被大幅高估。

#### 80 FPs 深度驗證

對 paired Tier A 中 LOH=Y + HPMergedSig=Y 的 80 個 FP 做 per-sample 與 per-chromosome 分解，確認 FP 集中度。

#### 雙層 Tier 拆分

將 Tier A(>=30) 拆分為 A(30-49) 與 A+(>=50)，分別計算 enrichment，確認是否是 uniform enrichment。

### 關鍵發現

#### 1. 真實 F1 基線修正

| 樣本 | Truth Total | TP | FP | FN | Round 3 顯示 F1 | **Round 4 真實 F1** |
|------|-----------|-----|-----|-----|----------------|-------------------|
| H2009 | 142,091 | 132,879 | 85 | 9,212 | ~~0.9997~~ | **0.9662** |
| H1437 | 81,016 | 67,460 | 8 | 13,556 | — | **0.9087** |
| HCC1395_DORADO | 39,447 | 29,877 | 238 | 9,570 | — | **0.8590** |
| HCC1395 | 39,447 | 29,752 | 544 | 9,695 | ~~0.9896~~ | **0.8532** |
| HCC1937 | 16,867 | 12,392 | 195 | 4,475 | — | **0.8414** |
| HCC1954 | 23,030 | 17,864 | 29 | 5,166 | — | **0.8731** |
| COLO829 | 43,192 | 35,184 | 2,273 | 8,008 | ~~0.9689~~ | **0.8725** |

H2009 真實 F1 = 0.9662（非 0.9997），差距 -0.034。但 F1 delta（filter 後的變化量）與 Round 3 幾乎相同，因為 LOH 數據覆蓋率接近 100%。

（數據來源：Round 4 fig01_baseline_f1.png 對應數據）

#### 2. 80 個 LOH+HPSig FP 中 70 個來自 HCC1395（87.5%）、66 個在 chr8（82.5%）

**Per-sample FP 分解**：

| 樣本 | FP 數 | TP 數 | FP% |
|------|-------|-------|-----|
| **HCC1395** | **70** | 393 | **15.1%** |
| H2009 | 3 | 630 | 0.5% |
| HCC1395_DORADO | 3 | 150 | 2.0% |
| HCC1954 | 2 | 22 | 8.3% |
| COLO829 | 1 | 8 | 11.1% |
| HCC1937 | 1 | 25 | 3.8% |
| H1437 | 0 | 118 | 0.0% |

**Per-chromosome FP 分解**：

| 染色體 | FP 數 | 佔比 |
|--------|-------|------|
| **chr8** | **66** | **82.5%** |
| chr3 | 4 | 5.0% |
| chr7 | 4 | 5.0% |
| chr17 | 3 | 3.8% |
| 其他 | 3 | 3.7% |

HCC1395 有 chr8 LOH，且 chr8 同時存在 allele-specific methylation（ASM）。chr8 的 FP 特別集中在 LOH=Y+HPSig=Y 類別中。

（數據來源：Round 4 fig02_fp80_per_sample 與 fig04_chr_distribution_fp 對應數據）

#### 3. 排除 HCC1395/HCC1954 後，7.4x 縮至 1.3x

| 組別 | 包含所有樣本 | 排除 HCC1395+HCC1954 |
|------|------------|---------------------|
| LOH=Y+HPSig=Y FP% | **5.61%** | 0.85% |
| LOH=Y+HPSig=N FP% | 0.76% | 0.64% |
| 相對差距 | **7.4x** | **1.3x** |

**結論**：「LOH + HPMergedSig = 7.4x FP 訊號」是 **HCC1395 chr8 LOH + ASM 的樣本特異現象**，而非普遍規律。

（數據來源：Round 4 報告 4.4 節一致性測試結果）

#### 4. Tier A(30-49) 與 A+(>=50) 方向相反

| Tier | Paired Enrichment | p-value | 方向 |
|------|------------------|---------|------|
| A (30-49) | **0.43x** | 4.1e-79 | **TP 富集（LOH 是 TP 指標）** |
| A+ (>=50) | **2.018x** | 2.8e-67 | **FP 富集（LOH 是 FP 指標）** |
| A combined (>=30) | 1.169x | 7.2e-7 | 兩效應加權平均 |

**A(30-49) 詳細數字**：
- A(30-49) TP 總數 48,088，其中 LOH-like 23,854（49.6%）
- A(30-49) FP 總數 1,065，其中 LOH-like 227（21.3%）
- 在 30-49 reads 的 region 中，LOH-like 是 TP 的可能性更高

**解讀**：
- 30-49 reads（中等 LOH support）：正確偵測到 somatic SNV 位於 LOH region 的正常情況，TP 多 LOH-like
- >=50 reads（超高 LOH support）：古老/克隆性 LOH（如染色體臂 LOH），germline SNP 在此容易被誤判為 somatic，FP 多 LOH-like

Round 2/3 報告的「Tier A(>=30) enrichment 1.169x」需補充：這是 A(30-49) 反向 0.43x 與 A+(>=50) 正向 2.02x 的加權平均，不是 uniform enrichment。

（數據來源：Round 4 報告第 5 節手動計算結果）

#### 5. 全量 benchmark filter simulation 仍全部 F1 delta 負值

| 樣本 | Baseline F1 | TierAplus_LOH_HPSig F1 | delta |
|------|------------|------------------------|-------|
| HCC1395 | 0.8532 | 0.8491 | -0.0041 |
| H2009 | 0.9662 | 0.9641 | -0.0021 |
| COLO829 | 0.8725 | 0.8725 | +/-0.0000 |
| HCC1954 | 0.8731 | 0.8726 | -0.0005 |

最保守的 TierAplus_LOH_HPSig 在 COLO829 達到中性（+/-0.0000），但沒有任何樣本正向受益。

（數據來源：Round 4 fig06_full_benchmark_f1_delta 對應數據）

### 圖表

**Fig01 -- Baseline F1**（含 FN 的真實 benchmark F1）

![Fig01 Baseline F1](figures/fig01_baseline_f1.png)

**Fig05 -- Joint LOH x HPMergedSig Heatmap**（2x2 FP rate heatmap，paired + TO）

![Fig05 Joint Heatmap](figures/fig05_joint_heatmap.png)

**Fig04 -- Chromosomal Distribution of FPs**（80 FPs 的染色體分佈，chr8 82.5%）

![Fig04 Chr Distribution FP](figures/fig04_chr_distribution_fp.png)

---

## 四輪研究總結論

### 結論一覽

| 結論 | 輪次 | 概述 |
|------|------|------|
| HP_Ratio=0.5 不代表平衡 | R1 | 69,807 regions effective_hp=0，Tier C0 必須排除 |
| Paired LOH-like FP enrichment 存在 | R1 | 整體 1.194x，但 sample heterogeneity 大 |
| TO LOH-like 是 TP 富集（非 FP marker） | R1+R2（修正後）| 全樣本一致方向 0.805x |
| Tier A(>=30) 是信號來源 | R2 | 1.169x (p=7.2e-7)，Tier B 反轉是 artifact |
| HP0 在 paired/TO LOH 方向相反 | R2 | Paired LOH: 低 HP0（真實 phasing）；TO LOH: 高 HP0（部分 artifact） |
| HP0 filter 假設否定 | R3 | High HP0 TP%=76.7% > Low HP0 74.8% |
| LOH+HPMergedSig 7.4x 是 HCC1395 chr8 特異 | R3+R4 | 排除 HCC1395/HCC1954 後僅 1.3x |
| 所有 filter F1 delta < 0 | R3+R4 | LOH 不可作 binary filter |
| A(30-49) vs A+(>=50) 方向相反 | R4 | A: 0.43x（TP 富集）；A+: 2.02x（FP 富集） |

### LOH 的最終定位

**LOH 在 paired somatic calling 中的正確用法 -- FP 風險計分因子**：

```
LOH_FP_risk_score (paired) =
    IF effective_hp_reads >= 50 AND core_loh_like:
        weight = 2.02   (Tier A+ = 最強 FP 訊號)
    ELIF effective_hp_reads 30-49 AND core_loh_like:
        weight = 0.43   (Tier A = 反向，是 TP 富集指標)
    ELSE:
        weight = 1.0    (無效應)
```

在 evidence panel 設計中，Tier A（30-49）的 LOH 應賦予**負 FP 風險**（TP 保護），Tier A+（>=50）賦予正 FP 風險。

**TO 端**：LOH 是弱 TP 富集信號（enrichment 0.805x），方向與 paired 相反，需獨立分析框架。

### 方法論注意事項

1. **Sample composition artifact**：多個「全域」enrichment 數字受個別樣本主導（HCC1395 主導 LOH+HPSig；COLO829 主導 Tier B 反轉；HCC1954 主導 TO FP 統計）。任何全域結論都需要 per-sample 驗證。

2. **LOH FP coverage 問題**：HCC1395 的 LOH FP coverage 達 115%（83 個 FP 在 LOH 數據中但不在 benchmark truth scope 內），可能影響該樣本的 enrichment 統計。

3. **Tier enrichment 計算方式**：Round 4 確認 `compute_enrichment_tier` 函式應使用 tier 特定總數作為分母（非全部 paired 資料），Round 4 的手動計算已修正此問題。

---

## 待驗證問題（已驗證 / 已更新）

### 已解決

1. **~~TO LOH TP 富集的機制~~** ✅ **根因：TO phasing 缺乏 normal reference 的系統性 HP 偏移**
   - 方向翻轉完全由 TP LOH rate 驅動：Paired TP 29.3% → TO TP **44.5%**（+51.6%），FP 幾乎不變（35.0% → 35.8%）
   - Shapley 分解：TP rate 變動貢獻 **105.6%**，FP rate 變動僅 -5.4%
   - 機制：287,092 個同位點 TP 中，39,724 個從 Paired non-LOH 翻轉為 TO LOH-like；其中 71.6% 在 TO 的 min(HP1,HP2)=0（單 haplotype 完全無 reads），同位點在 Paired 維持 HP1:HP2 = 8:8
   - TO phasing 缺乏 normal reference → 約 13% TP 位點所有 reads 被分配到單一 haplotype → hp_ratio 極端化 → 自動歸類為 LOH-like
   - 7/7 樣本全部呈相同翻轉方向（Paired FP-enriched, TO TP-enriched），確認為系統性屬性
   - **啟示**：TO LOH-like 信號包含大量 phasing artifact，LOH penalty 在 TO 模式應關閉（已完成 b9eaba7）

2. **~~A(30-49) 為何 TP 富集~~** ✅ **根因：hp_ratio 離散化放大 + Simpson's Paradox**
   - 30-49 reads 時 hp_ratio 更易觸發 LOH-like 閾值（53.2% TP 落在極端值 vs A+ 的 33.1%），放大 TP LOH rate
   - 6/7 樣本 tier-specific 下 A+ 實際優於 A，池化後反轉（Simpson's Paradox）
   - LOH 判別甜蜜點：hp_reads 45-50（enrichment 0.634），>100 時失效（~1.08）

### 尚未解決

3. **HCC1395 chr8 ASM+LOH block 深度分析**：chr8 的 LOH+HPMergedSig 7.4x enrichment 是 HCC1395 特異性，根因尚未確認。定位 P2 待資源。

4. **Tier A+ enrichment 的 genomic 分佈**：A+(>=50) 的 2.02x enrichment 是否集中在特定 chromosome arm 或已知 CNV 區域？定位 P2。

5. **跨 PS block 的 HP1/HP2 一致性**：目前所有 HP-driven 結論僅限於 local phase block 內，Paired 仍缺 variant-level PS export。

---

## 認知門檻補充建議

### 給初次閱讀者的背景知識

1. **LOH（Loss of Heterozygosity）**：腫瘤中一個等位基因的喪失。在 long-read sequencing 中，LOH 表現為 phased reads 集中在一個 haplotype（HP1 或 HP2）。LOH 本身是常見的腫瘤事件，不代表 variant call 的正確或錯誤。

2. **HP_Ratio vs hp_ratio_core**：原始 HP_Ratio 包含 Laplace smoothing 且在 effective_hp=0 時回傳 0.5（假象平衡）。hp_ratio_core 排除 HP0/HP3 reads，只計算有效 phased reads 的比例，是更精確的 LOH 指標。

3. **Enrichment 的解讀**：Enrichment = FP_LOH% / TP_LOH%。數值 1.2x 代表「FP 中有 LOH-like 的比例比 TP 高 20%」。**但 enrichment 高不等於可以做 filter**：當 FP 絕對數量遠小於 TP 時，移除 LOH-like region 會移除大量 TP，F1 反而下降。

4. **Tier 系統的必要性**：讀數量（effective_hp_reads）直接影響 LOH 判斷的統計力。1 vs 0 reads 的 hp_ratio=1.0 不代表 LOH，50 vs 0 reads 的才可能是。Tier 系統把不同統計力的 region 分開分析。

5. **Paired vs TO 的根本差異**：Paired 模式有 normal sample 作為 phasing 參考，LOH-like 較可能反映真實 LOH。TO（tumor-only）缺少 normal 對照，LOH-like 可能是 phasing artifact、tumor purity 效應或真實 LOH 的混合。兩者不可共用同一套結論。

### 圖表閱讀指引

- **Enrichment heatmap**：紅色代表 FP 富集（enrichment > 1），藍色代表 TP 富集（enrichment < 1）。顏色越深代表效應越強。
- **F1 delta bar chart**：所有 bar 都在負值區域，代表所有 filter 都損害整體效能。bar 越短（越接近 0）代表損害越小。
- **Joint heatmap（2x2）**：比較 LOH x HPMergedSig 四個象限的 FP%。右上角（LOH=Y, HPSig=Y）如果特別高，代表聯合特徵有 FP marker 潛力。

---

## 數據與圖表完整索引

### Round 1 數據

| 檔案 | 路徑 |
|------|------|
| all_region_rows.tsv.gz | `observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz` |
| loh_enrichment_by_sample_mode.tsv | 同 workspace |
| hp_coverage_qc_summary.tsv | 同 workspace |
| same_locus_compare.tsv | 同 workspace |
| verificationclass_by_loh_subtype.tsv | 同 workspace |

### Round 2 數據

| 檔案 | 路徑 |
|------|------|
| core1_tier_enrichment_global.tsv | `observation_workspaces/20260327_loh_round2_support_hp0_analysis/` |
| core1_tier_enrichment_by_sample_mode.tsv | 同 workspace |
| core2_hp0_by_loh_status.tsv | 同 workspace |
| core2_hp0_by_tier.tsv | 同 workspace |
| core2_affinity_nhp0gt0_only.tsv | 同 workspace |

### Round 3 數據

| 檔案 | 路徑 |
|------|------|
| HP0 threshold sweep 結果 | `observation_workspaces/20260327_loh_round3_methyl_hp0_filter/` |
| Joint LOH x HPMergedSig 分析 | 同 workspace |
| Tier threshold sensitivity 分析 | 同 workspace |

### Round 4 數據

| 檔案 | 路徑 |
|------|------|
| 全量 benchmark F1 基線 | `observation_workspaces/20260328_loh_round4_final_validation/` |
| 80 FPs per-sample/per-chr 分解 | 同 workspace |
| 雙層 Tier enrichment | 同 workspace |

所有 workspace 根目錄：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/`

### 圖表總索引

| 圖號 | 輪次 | 檔名 | 內容 |
|------|------|------|------|
| R1-01 | R1 | fig01_loh_like_fraction_overview.png | 各樣本 paired/TO LOH-like 佔比概覽 |
| R1-02 | R1 | fig02_hp_ratio_core_distribution.png | hp_ratio_core 分佈 |
| R2-01 | R2 | fig01_tier_enrichment_global.png | 全域 Tier enrichment |
| R2-02 | R2 | fig02_paired_tier_enrichment_heatmap.png | 各樣本 x Tier paired enrichment 熱圖 |
| R3-05 | R3 | fig05_filter_simulation_f1_delta.png | Filter simulation F1 delta |
| R4-01 | R4 | fig01_baseline_f1.png | 含 FN 的真實 benchmark F1 |
| R4-04 | R4 | fig04_chr_distribution_fp.png | 80 FPs 染色體分佈 |
| R4-05 | R4 | fig05_joint_heatmap.png | 2x2 LOH x HPSig FP rate heatmap |
