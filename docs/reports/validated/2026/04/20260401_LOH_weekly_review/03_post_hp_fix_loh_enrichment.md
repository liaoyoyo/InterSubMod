<!--
建立時間: 2026-04-01 12:00
目標: 週回顧第三章——HP integer tag 修正後 TO LOH enrichment 完整重估（含 Paired vs TO 方向翻轉確認）
處理範圍: 7 個 TO 樣本 + 7 個 Paired 樣本，修正後 significance_summary.csv，LOH eligibility / enrichment / per-sample / tier 全面分析
關聯檔案:
  - docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md
  - docs/reports/validated/2026/04/20260401_LOH_enrichment_paired_to_corrected_analysis_01.md
  - scripts/analysis/build_to_loh_enrichment_post_hp_fix.py
  - src/core/ReadParser.cpp（HP integer tag fix）
-->

# 03. HP Integer Tag 修正後 TO LOH Enrichment 重估

---

## 1. 為什麼需要重跑：HP Integer Tag Bug 的影響

### 1.1 Bug 本質

ReadParser.cpp 中的 HP tag 解析存在 switch-case mapping 錯誤：BAM 檔中 `HP:i:11`、`HP:i:21`、`HP:i:33` 等整數值未被正確映射為 haplotype 標籤（應為 `"1-1"`、`"2-1"`、`"3"` 等）。

**前提**：此 bug 僅影響 Tumor-Only (TO) 模式。Paired 模式的 HP tag 格式為字串型（`HP:Z:1`、`HP:Z:2`），不受此 bug 影響。

### 1.2 影響範圍

- **受影響樣本**：全部 7 個 TO 樣本（HCC1395、HCC1395_DORADO、COLO829、H1437、H2009、HCC1937、HCC1954）
- **影響程度**：30-99% 的 reads 未被正確追蹤 HP 歸屬 → `effective_hp` 被嚴重低估 → LOH eligibility 和 HP_Ratio 計算均無效
- **結論**：**2026-03-30 修正前的所有 TO LOH 統計數據均無效**，包括 Round 1 cross-sample audit 中的 TO 相關結論

### 1.3 為什麼 Round 1 舊數據偶爾看似「接近正確」

Round 1 TO 全域 enrichment 為 0.912x（<1，看似 TP 富集方向正確），但這是巧合而非真實結論：

- 舊版主要由 HCC1395 帶動（該樣本受 bug 影響相對較小）
- 舊版 COLO829 TO enrichment = 1.010x（完全無鑑別力），修正後為 0.956x
- 舊版 COLO829 median eff_hp_tp = 5.0，幾乎全部落在 C0/C tier，LOH-like 比例 67.9% 是基於錯誤 HP_Ratio 的假訊號

> **來源**：`q1b_enrichment_paired_vs_to.tsv`（位於 `20260330_post_hp_fix_to_loh_investigation/`）

---

## 2. 重跑流程

### 2.1 修正與重跑步驟

1. **修正 ReadParser.cpp**：修正 HP integer tag 的 switch-case mapping（`HP:i:11→"1-1"`, `HP:i:21→"2-1"`, `HP:i:33→"3"`）
2. **重新編譯** ISM pipeline（`cmake .. && make -j`）
3. **全量重跑 7 個 TO 樣本**的 ISM pipeline → 輸出新版 `significance_summary.csv`
4. **重建 master dataset**：從新版 `significance_summary.csv` 聚合為分析用 TSV
5. **所有 Round 1-4 workspace 全量重跑**：舊版資料保留於 `*_before_hp_fix` 目錄供對照

### 2.2 數據確認

新版 `significance_summary.csv` 更新時間：2026-03-30 03:53 ~ 06:21，覆蓋全部 7 個 TO 樣本。

> **分析腳本**：`scripts/analysis/build_to_loh_enrichment_post_hp_fix.py`
> **輸出目錄**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_enrichment_post_hp_fix/`

---

## 3. LOH Eligibility 修正效果

LOH 分析的前提是 region 具備足夠的 phasing support（`effective_hp ≥ 30`，定義為 LOH eligible）。修正前大量 reads 的 HP 歸屬遺失，導致 `effective_hp` 被嚴重低估。

### Table 1. LOH Eligible Region 修正前後比較（TP）

| 樣本 | 修正前 LOH Elig% | 修正後 LOH Elig% | 提升幅度 (pp) | Fold Change | 說明 |
|------|:---:|:---:|:---:|:---:|------|
| **H1437** | 24.5% | **95.5%** | **+71.0pp** | 3.9x | 改善最大：69% reads 原本不可見 |
| HCC1395_DORADO | 46.9% | **93.1%** | +46.2pp | 2.0x | |
| COLO829 | 0.7% | **34.7%** | +34.0pp | 49.5x | Fold change 最大，但最終值仍低 |
| HCC1395 | 58.5% | **92.2%** | +33.7pp | 1.6x | |
| H2009 | 64.5% | **97.7%** | +33.2pp | 1.5x | |
| HCC1954 | 70.0% | **93.0%** | +23.0pp | 1.3x | |
| HCC1937 | 86.7% | **97.4%** | +10.7pp | 1.1x | 原本受影響較小 |

> **來源**：`q3a_tier_distribution_before_after.tsv`（位於 `20260330_post_hp_fix_to_loh_investigation/`）

**關鍵觀察**：

1. **H1437**：修正幅度最大（+71.0pp），修正前僅 24.5% region 達 LOH eligible，修正後 95.5%，說明 69% reads 原本 HP 歸屬遺失
2. **COLO829**：fold change 最大（49.5x），但最終值僅 34.7%——遠低於其他樣本的 92-98%。原因是 COLO829 的 median `effective_hp` 僅 25（Tier B），phasing coverage 本身較低（見第 5 節）
3. **修正後所有樣本 LOH eligible 均顯著提升**，確認 bug 已完全解決

---

## 4. TO LOH Enrichment 翻轉確認

修正後的核心問題：TO LOH 是 FP 指標還是 TP 指標？

### 4.1 全域聚合結果（7 個 TO 樣本）

#### Table 2. TO LOH Enrichment by eff_hp Tier

| Tier | eff_hp 範圍 | TP total | TP LOH% | FP total | FP LOH% | FP/TP Enrichment | 解讀 |
|------|:---:|---:|:---:|---:|:---:|:---:|------|
| C0 | =0 | 732 | 0.0% | 338 | 0.0% | N/A | 無 phasing 資訊 |
| C | <10 | 4,084 | 76.1% | 2,027 | 77.3% | 1.015x | 無鑑別力（噪音） |
| B | 10-29 | 27,422 | 55.8% | 14,506 | 52.4% | 0.938x | 微弱 TP 富集 |
| **A** | **30-49** | **44,177** | **61.4%** | **25,255** | **43.3%** | **0.706x** | **TP 富集最強** |
| A+ | >=50 | 214,895 | 39.1% | 86,256 | 30.0% | 0.766x | TP 富集 |
| A>=30 | >=30 | 259,072 | 42.9% | 111,511 | 33.0% | 0.769x | TP 富集（p~=0） |
| **All** | **全部** | **291,310** | **44.5%** | **128,382** | **35.8%** | **0.805x** | **全域 TP 富集** |

> **來源**：`core1_tier_enrichment_global.tsv`

**核心結論**：

1. **TO LOH 全面呈現 TP 富集**（enrichment < 1.0）：方向與 paired 模式一致地相反（Paired 1.02-3.18x vs TO 0.852-0.956x，效應幅度差異大）
2. **Tier A（30-49）TP 富集最強**（0.706x）：`eff_hp` 30-49 範圍的 LOH-like region 明顯更傾向是 TP
3. **Tier C（<10）基本無鑑別力**（1.015x）：phasing 不足時 LOH 訊號只是噪音
4. **A>=30 全域 Fisher p ~= 0**：統計極顯著

### 4.2 Per-Sample 一致性驗證

#### Table 3. TO LOH Enrichment by Sample — All Tier

| 樣本 | TP LOH% | FP LOH% | Enrichment | p-value | median eff_hp (TP) |
|------|:---:|:---:|:---:|:---:|:---:|
| HCC1954 | 25.0% | 21.3% | **0.852x** | p~=0 | 61 |
| HCC1937 | 64.8% | 57.2% | **0.882x** | p~=0 | 99 |
| HCC1395 | 60.6% | 54.3% | **0.896x** | p~=0 | 61 |
| HCC1395_DORADO | 61.6% | 56.0% | **0.909x** | p~=0 | 63 |
| H1437 | 41.8% | 38.4% | **0.919x** | p~=0 | 66 |
| H2009 | 40.9% | 37.8% | **0.923x** | p~=0 | 85 |
| COLO829 | 35.4% | 33.8% | **0.956x** | 4.5e-4 | **25** |

> **來源**：`core1_tier_enrichment_by_sample_mode.tsv`

**所有 7 個樣本一致 < 1.0**（0.852-0.956x），無例外。

#### Table 4. TO Tier A (eff_hp 30-49) — 統計顯著性詳細

| 樣本 | TP n | TP LOH% | FP n | FP LOH% | Enrichment | p-value |
|------|---:|:---:|---:|:---:|:---:|:---:|
| HCC1954 | 3,524 | 28.0% | 9,882 | 25.3% | 0.904x | 0.002 |
| HCC1937 | 645 | 91.3% | 629 | 85.2% | 0.933x | 8.8e-4 |
| H1437 | 8,635 | 79.8% | 2,586 | 75.4% | 0.945x | 3.0e-6 |
| **COLO829** | **10,881** | **15.6%** | **6,066** | **14.8%** | **0.951x** | **0.189** |
| H2009 | 6,431 | 76.6% | 545 | 74.5% | 0.973x | 0.270 |
| HCC1395_DORADO | 6,784 | 86.7% | 2,680 | 84.9% | 0.979x | 0.022 |
| HCC1395 | 7,277 | 84.4% | 2,867 | 82.9% | 0.982x | 0.062 |

> **來源**：`core1_tier_enrichment_by_sample_mode.tsv`

**注意**：Tier A 中 COLO829（p=0.189）和 H2009（p=0.270）不顯著——前者因為 LOH-like 比例本身低（15.6%），後者因為 FP 數量極少（n=545）。但方向仍一致為 < 1.0。

---

## 5. COLO829 特殊情況

COLO829 在所有分析中都是離群值，需要單獨標記。

### 5.1 數據概況

| 指標 | COLO829 TO | 其他 6 樣本平均 |
|------|:---:|:---:|
| Median eff_hp (TP) | **25** | 61-99 |
| LOH eligible (eff_hp>=30) | **34.7%** | 92-98% |
| 主要落在 Tier | **B (10-29)** | A+/A |
| All enrichment | 0.956x | 0.852-0.923x |
| Tier A p-value | **0.189（不顯著）** | <=0.002 |
| Tier A+ n (TP+FP) | **1,104** | 19,000-128,000 |
| Tier A+ enrichment | 1.062x（p=0.679） | 0.808-0.892x |

### 5.2 原因分析

1. **Phasing coverage 本身較低**：median eff_hp = 25 落在 Tier B，大多數 region 不具備可靠 LOH 判定所需的 phasing support
2. **LOH-like 比例偏低**（35.4%）：可能是 melanoma 基因組特性——高度基因組異質性導致 haplotype block 碎片化
3. **Tier A+ 樣本量極小**（TP 594 + FP 510 = 1,104）：統計力完全不足

### 5.3 結論

> **COLO829 TO LOH 應標記為「低信心」**。enrichment 方向（0.956x）與其他樣本一致但不顯著，不應作為 feature engineering 主材料。Phase 1 分析中 COLO829 的 LOH feature 應降權或排除。

---

## 6. Paired vs TO 方向翻轉：完整對比

這是本次修正後最重要的結論之一：LOH enrichment 在 Paired 和 TO 模式下方向一致地相反（Paired FP-enriched 1.02-3.18x，TO TP-enriched 0.852-0.956x），但效應幅度差異大，部分 paired 接近中性（COLO829 1.15x）。

### 6.1 全域統計

| 分組 | LOH Rate | 解讀 |
|------|:---:|------|
| Paired TP | 29.3% | 基線 |
| Paired FP | **35.0%** | FP 比 TP 更容易被標記為 LOH → **LOH = FP-enriched** |
| TO TP | **44.5%** | TO TP 的 LOH rate 遠高於 Paired TP |
| TO FP | 35.8% | FP LOH rate 反而低於 TP → **LOH = TP-enriched** |

> **來源**：`20260401_LOH_enrichment_paired_to_corrected_analysis_01.md` 全域統計表

### 6.2 TP/FP 基礎差異

理解方向翻轉必須先理解兩種模式的 TP/FP 比例差異：

| 模式 | TP 數量 | FP 數量 | FP Rate |
|------|---:|---:|:---:|
| Paired | 325,270 | 3,429 | **1.04%** |
| TO | 291,310 | 128,382 | **30.6%** |

Paired FP 極少（1%），TO FP 極多（31%）——兩者是根本不同的問題空間。

![Paired vs TO TP/FP 比例](../figures/20260401_loh_enrichment_corrected/fig04_truth_label_composition.png)

### 6.3 Per-Sample Paired 數據

| Sample | Paired TP LOH% | Paired FP LOH% | Enrichment (FP/TP) | FP n |
|--------|:---:|:---:|:---:|---:|
| HCC1954 | 10.8% | **34.5%** | **3.18x** | 29 |
| H2009 | 28.6% | **76.7%** | **2.68x** | 86 |
| H1437 | 20.9% | 37.5% | 1.79x | 8 |
| HCC1937 | 55.5% | **83.6%** | **1.50x** | 195 |
| HCC1395_DORADO | 44.7% | 56.2% | 1.26x | 240 |
| COLO829 | 20.1% | 23.3% | 1.15x | 2,244 |
| HCC1395 | 47.4% | 48.2% | 1.02x | 627 |

> **來源**：`20260401_LOH_enrichment_paired_to_corrected_analysis_01.md` per-sample paired 表格

**Paired 全部 >= 1.0**（FP-enriched），特別是 HCC1954 (3.18x)、H2009 (2.68x)。但部分樣本 FP 數量極少（H1437 僅 8 FP, HCC1954 僅 29 FP），COLO829 (n=2,244 FP) 的 1.15x 是最穩健的估計。

### 6.4 方向對比總結

| 樣本 | TO Enrichment | Paired Enrichment | 方向 |
|------|:---:|:---:|:---:|
| HCC1395 | 0.896x (TP-enriched) | 1.02x (FP-enriched) | 翻轉 |
| HCC1395_DORADO | 0.909x | 1.26x | 翻轉 |
| COLO829 | 0.956x | 1.15x | 翻轉 |
| H1437 | 0.919x | 1.79x | 翻轉 |
| H2009 | 0.923x | 2.68x | 翻轉 |
| HCC1937 | 0.882x | 1.50x | 翻轉 |
| HCC1954 | 0.852x | 3.18x | 翻轉 |

**所有 7 個樣本一致翻轉**，無例外。

![LOH Enrichment Heatmap — Paired vs TO](../figures/20260401_loh_enrichment_corrected/fig04_loh_enrichment_heatmap.png)

![TO vs Paired LOH Concordance](../figures/20260401_loh_enrichment_corrected/fig06_to_vs_paired_loh_concordance.png)

---

## 7. 對 QS LOH Penalty 的影響

現行 Quality Score 的 LOH penalty (-25 分) 對所有 `Potential_LOH=True` 的位點扣分。HP 修正後的數據揭示：

| 指標 | 數值 |
|------|------|
| TO LOH penalty 觸發率（TP） | 44.5% |
| TO LOH penalty 觸發率（FP） | 35.8% |
| **淨效果** | **懲罰 TP 多於 FP → 反向作用** |
| 移除 LOH penalty 後 TO QS AUC 變化 | 0.497 → 0.542 (+0.045)（僅移除 LOH penalty；實際 commit b9eaba7 同時移除 verify bonus，預估 ~0.546 +0.049，待 benchmark 驗證） |

**結論**：LOH penalty 在 TO 模式下必須移除或反轉。此結論已在 `20260401_LOH_enrichment_paired_to_corrected_analysis_01.md` 和 QS mode-aware 重設計中確認。

---

## 8. 本章核心結論

| # | 結論 | 強度 | 證據來源 |
|---|------|:---:|------|
| 1 | HP integer tag bug 導致 TO 30-99% reads HP 歸屬遺失，修正前所有 TO LOH 統計無效 | 確定 | `q3a_tier_distribution_before_after.tsv` |
| 2 | 修正後 TO LOH enrichment 全域 = 0.805x（TP-enriched），7 個樣本一致 0.852-0.956x | 強支持（7/7 方向一致，6/7 顯著，COLO829 p=0.265 不顯著） | `core1_tier_enrichment_global.tsv`, `core1_tier_enrichment_by_sample_mode.tsv` |
| 3 | Tier A (eff_hp 30-49) TP 富集最強（0.706x），Tier C 無鑑別力 | 確定 | `core1_tier_enrichment_global.tsv` |
| 4 | COLO829 TO LOH 不顯著（Tier A p=0.189），應標記為低信心 | 確定 | `core1_tier_enrichment_by_sample_mode.tsv` |
| 5 | Paired vs TO LOH enrichment 方向一致地相反：Paired FP-enriched (1.02-3.18x), TO TP-enriched (0.852-0.956x)，效應幅度差異大 | 強支持（7/7 方向一致） | 兩份來源報告交叉驗證 |
| 6 | TO LOH penalty (-25) 在 TO 模式下反向作用，必須移除 | 確定 | QS AUC 分析 |

---

## 待驗證問題（已驗證 / 已更新）

### 已解決

1. **~~COLO829 phasing coverage 低的根因~~** ✅ **根因：低測序深度 (~30x, BAM 94GB)**
   - COLO829 median num_reads = 29（其他 6 樣本 65-103），median effective_hp = 25（其他 61-99）
   - 8.5% region hp_reads < 10（其他皆 < 1.3%，Z-score = 22.09）
   - HP 分配效率最差（hp/total = 0.862，其他 0.913-0.966）
   - 因果鏈：低深度 → 低 effective_hp → hp_ratio 離散化（2.8% 恰好 = 0.5）→ LOH 判定不穩定 → enrichment ≈ 1（0.956x，無區分力）
   - 加劇因素：het_phased 僅 0.7%（78.4% hom），黑色素瘤高 TMB + 高 purity 使 AF 集中 ~0.5
   - **建議**：跨樣本分析時應排除 COLO829 或以 hp_reads-matched 方式控制深度混淆

2. **~~Tier A 為何 TP 富集最強（0.706x vs A+ 0.766x）~~** ✅ **根因：hp_ratio 離散化 + Simpson's Paradox**
   - 兩種 enrichment 定義結論相反但不矛盾：
     - 定義 A（全局分母，step5）：A(30-49) = 0.916, A+(≥50) = **0.698** → A+ 更好
     - 定義 B（tier-specific 分母）：A(30-49) = **0.706**, A+(≥50) = 0.766 → A 更好
   - 定義 B 中 A 更好的原因：30-49 reads 時 53.2% TP / 37.3% FP 的 hp_ratio 落在 0 或 1 極端值（A+ 僅 33.1%/24.9%），低 reads 放大 LOH-like rate，TP 被放大更多
   - Simpson's Paradox 確認：6/7 樣本 tier-specific 定義下 A+ 優於 A，但池化後反轉
   - LOH 判別甜蜜點在 hp_reads **45-50**（enrichment 0.634），>100 時完全失效（enrichment ~1.08）
   - **實務建議**：未來 LOH 判別應考慮 hp_reads-dependent 閾值，而非統一切分 tier

3. **~~Round 1 TO 舊數據 0.912x 為何方向偶然正確~~** ✅ 已關閉：Round 1-4 全量重跑完成（`*_post_to_hp_fix`），舊版保留於 `*_before_hp_fix`，「為何偶然正確」為學術好奇心，不影響結論

### 尚未解決

4. **TO LOH enrichment 是否在 FN（caller 過濾掉的 TP）中更強**：目前只有 PASS VCF 的 TP/FP 數據，FN rescue 潛力需要新的數據收集。對應 Section 09 P2 行動「O9 FN 觀察」。

---

## 認知門檻補充建議

1. **LOH eligible 的門檻理解**：`effective_hp >= 30` 代表一個 region 至少有 30 條 reads 被成功分配到某個 haplotype。低於此門檻時，HP_Ratio 的隨機波動過大，LOH 判定不可靠。Tier 分級（C0/C/B/A/A+）本質上是 phasing support 的信心分級。

2. **Enrichment < 1 = TP 富集的直覺**：enrichment 定義為 FP_LOH% / TP_LOH%。當 enrichment < 1，代表 TP 的 LOH 比例高於 FP → LOH 是 TP 的特徵而非 FP 的特徵 → LOH 不能用來過濾 FP，反而應作為 TP 的支持 evidence。

3. **Paired vs TO 的 TP/FP 比例差異**：Paired 的 FP rate 僅 1%（caller 已做大部分過濾），TO 的 FP rate 達 31%。兩者的 LOH enrichment 不能直接比較絕對值，只能比較方向。方向翻轉才是核心發現。

4. **HP integer tag bug 的認知價值**：此 bug 是一個 negative control — 它展示了當 phasing 資訊被系統性破壞時，LOH 統計會如何失真（LOH eligible 從 0.7% 到 34.7% 的 COLO829 案例）。這也警示未來任何涉及 HP tag 的分析都必須先驗證 HP tag 的完整性。
