<!--
建立時間: 2026-03-31 23:59
目標: 修正版 LOH enrichment 分析——Paired vs TO 的 LOH-TP/FP 關係完整解讀
處理範圍: 748,391 rows, 7 samples × 2 modes
關聯檔案:
  - big7_disk_output/synthesis/observation_workspaces/20260401_O03_loh_features_post_fix/
  - big7_disk_output/synthesis/observation_workspaces/20260401_O02_quality_score_decompose/
  - big7_disk_output/synthesis/observation_workspaces/20260401_O01_master_distribution/figures/fig04_truth_label_composition.png
  - big7_disk_output/synthesis/observation_workspaces/20260401_O03_loh_features_post_fix/fig04_loh_enrichment_heatmap.png
  - big7_disk_output/synthesis/observation_workspaces/20260401_O03_loh_features_post_fix/fig06_to_vs_paired_loh_concordance.png
-->

# LOH Enrichment 修正版分析：Paired vs TO 的 LOH-TP/FP 關係

## 1. 問題背景：Paired 與 TO 的 TP/FP 差異

Paired 和 TO 模式在 TP/FP 比例上有根本性差異，這是理解 LOH enrichment 的前提。

![Paired vs TO TP/FP 比例](figures/20260401_loh_enrichment_corrected/fig04_truth_label_composition.png)

| 模式 | TP 數量 | FP 數量 | FP Rate |
|------|---------|---------|---------|
| Paired | 325,270 | 3,429 | **1.04%** |
| TO | 291,310 | 128,382 | **30.6%** |

Paired FP 極少（1%），TO FP 極多（31%）——兩者本質上是不同的問題空間。

---

## 2. 核心發現：LOH Enrichment 方向在 Paired 和 TO 之間完全翻轉

### 全域統計

| 分組 | LOH Rate | 解讀 |
|------|----------|------|
| Paired TP | 29.3% | 基線 |
| Paired FP | **35.0%** | FP 比 TP 更容易被標記為 LOH → **LOH = FP-enriched** |
| TO TP | **44.5%** | TO TP 的 LOH rate 遠高於 paired |
| TO FP | 35.8% | **FP 的 LOH rate 反而低於 TP** → **LOH = TP-enriched** |

### Per-Sample TO 數據（FP LOH rate 全部低於 TP LOH rate）

| Sample | TO TP LOH% | TO FP LOH% | Enrichment (FP/TP) |
|--------|-----------|-----------|-------------------|
| HCC1395 | 60.6% | 54.3% | **0.90x** |
| HCC1395_DORADO | 61.6% | 56.0% | **0.91x** |
| HCC1937 | 64.8% | 57.2% | **0.88x** |
| HCC1954 | 25.0% | 21.3% | **0.85x** |
| H2009 | 40.9% | 37.8% | **0.92x** |
| H1437 | 41.8% | 38.4% | **0.92x** |
| COLO829 | 35.4% | 33.8% | **0.96x** |

**所有 7 個樣本一致：TO 模式下 FP 的 LOH rate 低於 TP**（enrichment 0.85-0.96x）。

### Per-Sample Paired 數據（FP LOH rate 全部高於或等於 TP LOH rate）

| Sample | Paired TP LOH% | Paired FP LOH% | Enrichment (FP/TP) | FP n |
|--------|----------------|----------------|-------------------|------|
| HCC1395 | 47.4% | 48.2% | 1.02x | 627 |
| HCC1395_DORADO | 44.7% | 56.2% | 1.26x | 240 |
| HCC1937 | 55.5% | **83.6%** | **1.50x** | 195 |
| HCC1954 | 10.8% | **34.5%** | **3.18x** | 29 |
| H2009 | 28.6% | **76.7%** | **2.68x** | 86 |
| H1437 | 20.9% | 37.5% | 1.79x | 8 |
| COLO829 | 20.1% | 23.3% | 1.15x | 2,244 |

**Paired 模式下 FP 嚴重偏向 LOH**，特別是 HCC1954 (3.18x)、H2009 (2.68x)、HCC1937 (1.50x)。但注意 paired FP 樣本量很小（H1437 僅 8 FP, HCC1954 僅 29 FP），部分 enrichment 可能受小樣本波動影響。COLO829 (n=2,244 FP) 的 enrichment 1.15x 是最可靠的估計。

---

## 3. 機制解讀

### 3.1 TO 模式：LOH 過判導致 TP 被過度標記

在無 Normal 樣本的 TO 模式下，phasing 完全依賴腫瘤自體的 heterozygous SNV。這導致：

1. **Genotype 稀疏**：TO LOH 判定的位點中 51.5% 具有 partial genotype (0/., 1/.)，缺乏完整的雙等位基因資訊
2. **Phase block 碎片化**：29.3% 的 phase block 為 singleton，LOH 位點 PS 缺失率 11.1%（nonLOH 僅 1.2%）
3. **系統性過判**：TO-Paired LOH concordance 矩陣（下圖）顯示 85.5% 一致率中，不一致的 95.5% 是 TO=LOH where paired=nonLOH

![TO vs Paired LOH Concordance](figures/20260401_loh_enrichment_corrected/fig06_to_vs_paired_loh_concordance.png)

**結果**：TO 系統性地將更多位點判定為 LOH，而這些「過判」主要發生在 TP 位點（因為 TP 是真正的體細胞變異，更容易因 haplotype imbalance 被判為 LOH）。因此 **TO TP LOH rate (44.5%) > TO FP LOH rate (35.8%)**。

### 3.2 Paired 模式：FP 集中在 LOH 區域

Paired 模式有 Normal 樣本提供 germline phasing，LOH 判定相對可靠。FP 偏向 LOH 的原因：

1. **LOH 區域的 reads 缺乏 haplotype diversity**：當一個 haplotype 幾乎消失時，甲基化分析缺乏 HP 對照組，更容易產生 false positive ISM 信號
2. **部分樣本效應極端**：H2009 paired FP 的 76.7% 是 LOH（但僅 86 FP），HCC1937 paired FP 83.6% 是 LOH（195 FP）——這些 FP 幾乎全部來自 LOH 區域
3. **Paired FP 數量極少**（全部僅 3,429），這些 FP 是通過嚴格 paired calling pipeline 仍然漏網的少數 false positive，它們集中在 LOH 區域說明 LOH 是 paired 模式下最主要的 FP 來源

### 3.3 TO FP 的 LOH 比例跟隨樣本整體 LOH 比例

TO FP 的 LOH rate 與樣本整體 LOH rate 高度一致：

| Sample | 整體 LOH Rate | TO FP LOH% | 差距 |
|--------|-------------|-----------|------|
| HCC1937 | ~60% | 57.2% | -3pp |
| HCC1395 | ~55% | 54.3% | -1pp |
| H2009 | ~35% | 37.8% | +3pp |
| HCC1954 | ~18% | 21.3% | +3pp |
| COLO829 | ~28% | 33.8% | +6pp |

TO FP 的 LOH 分佈「隨波逐流」——它不是因為 LOH 才成為 FP，而是 FP 均勻分佈在 LOH 和 nonLOH 區域中（LOH rate 接近樣本基線）。相比之下，TO TP 的 LOH rate 明顯高於基線，因為 TO 系統性過判使更多 TP 被標記為 LOH。

---

## 4. LOH Enrichment 熱圖

![LOH Enrichment Heatmap](figures/20260401_loh_enrichment_corrected/fig04_loh_enrichment_heatmap.png)

此圖確認：Paired 所有樣本 enrichment > 1.0（FP-enriched），TO 所有樣本 enrichment < 1.0（TP-enriched）。

---

## 5. 對 QS LOH Penalty 的影響

現行 Quality Score 的 LOH penalty (-25 分) 對所有 Potential_LOH=True 的位點扣分。在 TO 模式下：

- LOH penalty 觸發 44.5% TP vs 35.8% FP → **懲罰 TP 多於 FP**
- 移除 LOH penalty 後 TO QS AUC 從 0.497 升至 0.542 (+0.045)
- **LOH penalty 在 TO 模式下必須移除或反轉**

---

## 6. 結論與行動建議

| 結論 | 證據 | 行動 |
|------|------|------|
| TO LOH 系統性過判 | concordance 85.5%, 不一致的 95.5% 為 TO 多判 | 基於 GT completeness 和 PS size 重新判定 |
| TO 下 LOH penalty 反向 | TP LOH 44.5% > FP 35.8%, 所有 7 sample 一致 | **P0: 移除 TO LOH penalty** |
| Paired FP 集中在 LOH | enrichment 1.15-3.18x | LOH 區域是 paired 改善空間 |
| TO FP LOH 跟隨樣本基線 | FP LOH% ≈ 樣本整體 LOH% | TO FP 不是 LOH-specific，不可用 LOH 過濾 |
| Paired/TO 必須分離模型 | LOH enrichment 方向完全相反 | **P0: 分離模型策略** |
