<!--
建立時間: 2026-04-06
目標: Step 1.1 HP Imbalance vs LOH.bed 一致性 + Step 1.2 四象限 TP/FP 分佈
處理範圍: 7 samples, TO mode, LongPhase Baseline
關聯檔案:
  - scripts/analysis/build_loh_concordance_analysis.py
-->

# Step 1.1 + 1.2: LOH 雙定義一致性與四象限 TP/FP

## 背景

ISM 使用兩種方式定義 LOH：

| 定義 | 來源 | 粒度 | 方法 |
|------|------|------|------|
| **HP Imbalance** | ISM RegionProcessor | Per-variant | HP_Ratio < 0.1 or > 0.9 |
| **LOH.bed** | LongPhase-TO 輸出 | Region-level BED | Phased genotype ratio |

本分析量化兩者的一致性，並檢驗 LOH 狀態對 TP/FP 的區分貢獻。

---

## Step 1.1: 一致性分析

### Per-sample 一致性

| Sample | HP Imbalance Rate | LOH.bed Rate | Kappa | Jaccard | Sensitivity | Specificity |
|--------|-------------------|-------------|-------|---------|-------------|-------------|
| HCC1395 | 58.8% | 43.8% | 0.688 | 0.731 | 98.9% | 72.5% |
| HCC1395_DORADO | 60.0% | 44.0% | 0.688 | 0.733 | 100.0% | 71.5% |
| COLO829 | 34.8% | 20.9% | 0.662 | 0.601 | 100.0% | 82.4% |
| H1437 | 41.0% | 26.4% | 0.680 | 0.643 | 100.0% | 80.1% |
| H2009 | 40.6% | 24.8% | 0.649 | 0.609 | 99.9% | 78.9% |
| HCC1937 | 61.1% | 50.4% | 0.779 | 0.820 | 99.7% | 78.1% |
| HCC1954 | 22.2% | 6.3% | 0.378 | 0.281 | 99.8% | 83.0% |
| **Overall (TO)** | **41.8%** | **26.7%** | — | — | — | — |

### 四象限分佈

| 象限 | 定義 | 筆數 | 佔比 |
|------|------|------|------|
| Q1 | Both LOH (HP Imbalance + LOH.bed) | 111,932 | 26.7% |
| Q2 | ISM-only (HP Imbalance, 不在 LOH.bed) | 63,610 | 15.2% |
| Q3 | LOH.bed-only (在 LOH.bed, 無 HP Imbalance) | 286 | 0.1% |
| Q4 | Neither (兩者都不是 LOH) | 243,864 | 58.1% |

![Confusion matrix](figures/w1c_01_confusion_matrix.png)

![HP_Ratio distribution per quadrant](figures/w1c_02_hp_ratio_distribution.png)

### 關鍵發現

1. **Sensitivity 近乎完美 (99.7-100%)**：ISM HP Imbalance 幾乎捕捉所有 LOH.bed 區域
2. **Q3 (LOH.bed-only) 僅 286 筆 (0.07%)**：可忽略 → **ISM 是 LOH.bed 的超集 (J1)**
3. **ISM 傾向「過判」**：HP Imbalance rate (41.8%) 顯著高於 LOH.bed rate (26.7%)，差額全部落入 Q2
4. **HCC1954 一致性最低** (kappa=0.378)：HP Imbalance 22.2% vs LOH.bed 6.3%，大量 Q2 過判

![Kappa per sample](figures/w1c_03_kappa_barchart.png)

![HP_Ratio ROC for LOH.bed prediction](figures/w1c_04_hp_ratio_roc.png)

---

## Step 1.2: 四象限 TP/FP 分佈

### 核心數據

| 象限 | TP | FP | FP Rate | FP Enrichment |
|------|-----|-----|---------|---------------|
| Q1 (Both LOH) | 85,143 | 26,789 | **0.239** | 0.86 |
| Q2 (ISM-only) | 44,446 | 19,164 | 0.301 | 1.00 |
| Q3 (LOH.bed-only) | 200 | 86 | 0.301 | 1.09 |
| Q4 (Neither) | 161,521 | 82,343 | **0.338** | 1.06 |

![Quadrant TP/FP stacked bar](figures/w1c_05_quadrant_tp_fp_stacked.png)

![FP rate heatmap per sample × quadrant](figures/w1c_06_fp_rate_heatmap.png)

![FP enrichment](figures/w1c_07_fp_enrichment.png)

### 關鍵發現

1. **Q1 (Both LOH) 的 FP rate 最低 (0.239)**：LOH 區域 TP 比例反而更高（符合 tumor suppressor 失活生物學）
2. **FP enrichment 差異溫和 (0.86-1.09)**：LOH 象限歸屬對 TP/FP 區分的直接貢獻有限 (J3)
3. **Q4 (Neither) FP rate 最高 (0.338)**：非 LOH 區域 FP 比例更高
4. **驗證通過**：四象限 TP+FP 合計 = 每 sample 總數（數學恆等式 PASS）

---

## 操作條件與限制

- **數據版本**: Master dataset 2026-03-27 HP fix rebuild
- **HP Imbalance 定義**: HP_Ratio < 0.1 or > 0.9（ISM RegionProcessor.cpp）
- **LOH.bed**: LongPhase-TO baseline output（region-level BED）
- **限制**: Self-phasing 效應已確認存在（62% TO TP LOH 是 artifact），但本分析不直接量化 self-phasing
