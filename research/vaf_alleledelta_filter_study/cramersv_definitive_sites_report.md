# CramersV 明確位點觀察列表與分析報告

**專案**: InterSubMod - 甲基化輔助體細胞突變檢測  
**分析日期**: 2026-01-23  
**目的**: 分析 CramersV > 0 (明確計算値) 的 TP 和 FP 位點特徵差異

---

## 執行摘要

本報告專注於分析 **CramersV 值明確 (> 0)** 的位點，這些位點代表了 `is_reliable = true` 的統計可靠位點。

| 類別 | 總數 | CramersV > 0 | 比例 |
|------|------|--------------|------|
| TP | 30,476 | 1,874 | 6.15% |
| FP | 4,822 | 111 | 2.30% |

> [!IMPORTANT]
> **關鍵發現**: TP 中有明確 CramersV 的位點比例 (6.15%) 顯著高於 FP (2.30%)，這可能與真實突變位點更可能展現穩定的 ASM 信號有關。

---

## 1. CramersV 區間分布

以下為 CramersV > 0 位點的區間分布：

| 區間 | TP 數量 | TP % | FP 數量 | FP % |
|------|---------|------|---------|------|
| 極低 (0 ~ 0.05) | 4 | 0.21% | 3 | 2.70% |
| 低 (0.05 ~ 0.20) | 8 | 0.43% | 6 | 5.41% |
| 中 (0.20 ~ 0.50) | 174 | 9.28% | 10 | 9.01% |
| 高 (0.50 ~ 0.80) | 254 | 13.55% | 23 | 20.72% |
| 極高 (0.80 ~ 1.00) | 1,434 | 76.52% | 69 | 62.16% |

> [!NOTE]
> TP 位點傾向於有更高的 CramersV 值 (76.52% 處於 0.80-1.00 區間)，而 FP 相對分散。

---

## 2. 統計比較分析

### 2.1 CramersV 統計差異

| 統計量 | TP (n=1,874) | FP (n=111) | 差異 |
|--------|--------------|------------|------|
| Mean | 0.8502 | 0.7550 | TP 更高 |
| Std | 0.2016 | 0.2833 | FP 更分散 |
| Median | 0.9416 | 0.8874 | TP 更高 |
| IQR | [0.82, 1.00] | [0.56, 0.96] | FP 範圍更大 |

**Mann-Whitney U 檢定**: p = 1.44e-04 ***  
**Rank-biserial correlation**: -0.21

### 2.2 多特徵差異比較 (CramersV > 0 位點)

| 特徵 | TP mean | FP mean | TP median | FP median | 顯著性 |
|------|---------|---------|-----------|-----------|--------|
| **AlleleDelta** | 0.0980 | 0.0799 | 0.0940 | 0.0755 | *** |
| **VAF** | 0.4114 | 0.2181 | 0.3975 | 0.1532 | *** |
| **QUAL** | 0.9721 | 0.7190 | 0.9872 | 0.7023 | *** |
| **NumReads** | 75.3 | 109.7 | 70.0 | 94.0 | *** |
| HPMergedDelta | 0.0784 | 0.0648 | 0.0806 | 0.0593 | * |
| HeuristicScore | 20.93 | 19.59 | 21.88 | 21.77 | *** |

> [!WARNING]
> **特殊觀察**: FP 位點具有
> - 明顯較低的 VAF (0.22 vs 0.41)
> - 明顯較低的 QUAL (0.72 vs 0.97)
> - 較高的 NumReads (110 vs 75)

這表明即使在 CramersV 明確的 FP 中，VAF 和 QUAL 仍然是良好的區分特徵。

---

## 3. VerificationClass 分布

### TP (CramersV > 0)

| Class | 數量 | 比例 |
|-------|------|------|
| **Strong** | 1,840 | 98.2% |
| Subclone | 22 | 1.2% |
| Weak | 9 | 0.5% |
| Noise | 3 | 0.2% |

### FP (CramersV > 0)

| Class | 數量 | 比例 |
|-------|------|------|
| **Strong** | 102 | 91.9% |
| Noise | 7 | 6.3% |
| Weak | 2 | 1.8% |

> [!CAUTION]
> 91.9% 的 CramersV > 0 FP 位點被分類為 "Strong"！這意味著 CramersV 本身無法有效區分 TP 和 FP。

---

## 4. 高 CramersV 位點分析 (V ≥ 0.5)

| 指標 | TP (n=1,691) | FP (n=92) |
|------|--------------|-----------|
| AlleleDelta (mean) | 0.1044 | 0.0896 |
| VAF (mean) | 0.4199 | 0.2089 |
| VAF (median) | 0.4125 | 0.1482 |
| QUAL (mean) | 0.9737 | 0.7213 |
| NumReads (mean) | 73.5 | 111.3 |
| PassedGating | 100% | 100% |

---

## 5. 特殊情況分析

### 5.1 強 ASM 信號 (V ≥ 0.5 AND AD > 0.25)

| 類別 | 數量 |
|------|------|
| TP | 28 |
| FP | 3 |

這些位點展現出高統計可靠性和高甲基化差異。

### 5.2 矛盾情況 (V ≥ 0.5 但 AD ≤ 0.1)

| 類別 | 數量 |
|------|------|
| TP | 844 |
| FP | 71 |

> [!NOTE]
> 這些位點顯示高統計關聯性但低絕對差異。可能原因：
> 1. 小樣本統計的高敏感性
> 2. Cluster 與 Label 的完美對應但距離小

### 5.3 低 VAF + 高 V (VAF < 0.24 AND V ≥ 0.5)

| 類別 | 數量 | 說明 |
|------|------|------|
| TP | 196 | 可能為亞克隆突變 |
| FP | 61 | 可能為假 ASM |

> [!IMPORTANT]
> 這 61 個 FP 在目前的過濾條件下 (AD > 0.25 AND V < 0.05 AND VAF < 0.24) 不會被移除，因為它們的 CramersV ≥ 0.5。需要額外策略來處理。

---

## 6. 觀察列表

### 6.1 已生成的 CSV 檔案

| 檔案 | 內容 | 記錄數 |
|------|------|--------|
| [tp_cramersv_definitive.csv](file:///big8_disk/liaoyoyo2001/InterSubMod/analysis/vaf_alleledelta_filter_study/cramersv_sites/tp_cramersv_definitive.csv) | 所有 CramersV > 0 的 TP | 1,874 |
| [fp_cramersv_definitive.csv](file:///big8_disk/liaoyoyo2001/InterSubMod/analysis/vaf_alleledelta_filter_study/cramersv_sites/fp_cramersv_definitive.csv) | 所有 CramersV > 0 的 FP | 111 |
| [tp_high_cramersv.csv](file:///big8_disk/liaoyoyo2001/InterSubMod/analysis/vaf_alleledelta_filter_study/cramersv_sites/tp_high_cramersv.csv) | CramersV ≥ 0.5 的 TP | 1,691 |
| [fp_high_cramersv.csv](file:///big8_disk/liaoyoyo2001/InterSubMod/analysis/vaf_alleledelta_filter_study/cramersv_sites/fp_high_cramersv.csv) | CramersV ≥ 0.5 的 FP | 92 |

### 6.2 Top FP 位點 (高 CramersV)

以下是 CramersV 最高的 FP 位點，值得深入觀察：

| Chr | Pos | VAF | CramersV | AD | QUAL | NumReads | DominantLabel |
|-----|-----|-----|----------|-----|------|----------|---------------|
| chr9 | 41803018 | 0.14 | 1.0 | 0.08 | 0.82 | 156 | hp |
| chr9 | 41800364 | 0.19 | 1.0 | 0.09 | 0.88 | 138 | hp |
| chr9 | 41800163 | 0.12 | 1.0 | 0.08 | 0.79 | 173 | hp |
| chr9 | 41800991 | 0.15 | 1.0 | 0.08 | 0.70 | 162 | hp |
| chr4 | 47643402 | 0.07 | 1.0 | 0.23 | 0.60 | 70 | hp |
| chr8 | 82168704 | 0.25 | 1.0 | 0.26 | 0.98 | 72 | hp |
| chr11 | 55277533 | 0.06 | 1.0 | 0.09 | 0.56 | 97 | hp |

> [!TIP]
> **觀察重點**: chr9:41800163-41803018 區域有多個聚集的 FP，可能是一個 Region-level 的假陽性區域，需特別注意。DominantLabel 多為 "hp"，可能與 Haplotype 差異有關。

---

## 7. 視覺化

### 7.1 CramersV 分布與特徵關係

![CramersV Definitive Analysis](cramersv_sites/cramersv_definitive_analysis.png)

### 7.2 TP vs FP Box Plots

![CramersV Box Plots](cramersv_sites/cramersv_boxplots.png)

---

## 8. 結論與建議

### 8.1 主要發現

1. **TP 中有更多明確 CramersV 的位點** (6.15% vs 2.30%)
2. **TP 的 CramersV 值更高** (mean 0.85 vs 0.76)
3. **VAF 和 QUAL 仍是良好區分特徵** — 即使在 CramersV > 0 的子集中
4. **FP 的 NumReads 較高** — 可能與高覆蓋度區域更易產生假陽性有關
5. **chr9:41.8Mb 區域有聚集 FP** — 需特別關注

### 8.2 對過濾策略的啟示

目前過濾條件 `AD > 0.25 AND V < 0.05 AND VAF < 0.24` 無法觸及 CramersV > 0 的 FP。

可考慮的額外策略：

| 策略 | 潛在 FP 移除 | 潛在 TP 損失 | 優先級 |
|------|-------------|-------------|--------|
| 針對低 VAF + 高 V 添加 QUAL 過濾 | ~40 | ~10 | 高 |
| 區域黑名單 (如 chr9:41.8Mb) | ~15 | 0 | 中 |
| 針對 DominantLabel=hp 的特殊處理 | 待評估 | 待評估 | 低 |

---

**報告生成時間**: 2026-01-23  
**分析腳本**: `analyze_cramersv_definitive.py`  
**輸出目錄**: `/big8_disk/liaoyoyo2001/InterSubMod/analysis/vaf_alleledelta_filter_study/cramersv_sites/`
