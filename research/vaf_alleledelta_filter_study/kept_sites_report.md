# 過濾後剩餘 TP/FP 狀況分析報告

**專案**: InterSubMod - 甲基化輔助體細胞突變檢測  
**分析日期**: 2026-01-23  
**過濾條件**: `AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24` (不含 QUAL)

---

## 執行摘要

本報告分析過濾後**剩餘的 TP 與 FP** 狀況，評估統計顯著性分析與 VerificationClass 分類的變化。

### 過濾效果總結

| 指標 | 過濾前 | 過濾後 | 移除數 | 移除率 |
|------|--------|--------|--------|--------|
| **TP** | 30,476 | 30,470 | 6 | 0.02% |
| **FP** | 4,822 | 4,447 | 375 | 7.78% |
| **FP Strong** | 1,738 | 1,378 | 360 | **20.7%** |

> [!IMPORTANT]
> **關鍵成效**: 過濾條件精準移除了 20.7% 的 FP Strong 位點，同時幾乎不影響 TP (僅移除 0.02%)。

---

## 1. 剩餘位點統計

### 1.1 剩餘 TP 特徵統計

| 特徵 | 平均值 | 中位數 | 標準差 |
|------|--------|--------|--------|
| **AlleleDelta** | 0.0339 | 0.0086 | 0.0641 |
| **CramersV** | 0.0523 | 0.0000 | 0.2103 |
| **VAF** | 0.4940 | 0.4262 | 0.2702 |
| **QUAL** | 0.9672 | 0.9879 | 0.0599 |
| **NumReads** | 71.4 | 66 | 30.6 |

**VerificationClass 分布**:

| Class | 數量 | 比例 |
|-------|------|------|
| Weak | 11,342 | 37.22% |
| Noise | 10,437 | 34.25% |
| Strong | 7,270 | 23.86% |
| Subclone | 1,421 | 4.66% |

**DominantLabel 分布**:

| Label | 數量 | 比例 |
|-------|------|------|
| hp | 12,967 | 42.56% |
| allele | 12,710 | 41.71% |
| none | 4,793 | 15.73% |

### 1.2 剩餘 FP 特徵統計

| 特徵 | 平均值 | 中位數 | 標準差 |
|------|--------|--------|--------|
| **AlleleDelta** | 0.0580 | 0.0326 | 0.0712 |
| **CramersV** | 0.0188 | 0.0000 | 0.1259 |
| **VAF** | 0.1561 | 0.1167 | 0.1343 |
| **QUAL** | 0.6889 | 0.6606 | 0.1407 |
| **NumReads** | 60.2 | 53 | 32.1 |

**VerificationClass 分布**:

| Class | 數量 | 比例 |
|-------|------|------|
| Noise | 1,699 | 38.21% |
| **Strong** | **1,378** | **30.99%** |
| Weak | 1,138 | 25.59% |
| Subclone | 232 | 5.22% |

> [!WARNING]
> **仍有問題**: 剩餘 FP 中仍有 30.99% (1,378 個) 是 Strong 類別，這些是現有方法無法區分的高風險假陽性。

---

## 2. 過濾前後比較

### 2.1 VerificationClass 變化

![Before After Comparison](figures/before_after_comparison.png)

#### TP VerificationClass 變化

| Class | 過濾前 | 過濾後 | 移除 | 移除率 |
|-------|--------|--------|------|--------|
| Strong | 7,275 | 7,270 | 5 | 0.1% |
| Weak | 11,343 | 11,342 | 1 | 0.0% |
| Subclone | 1,421 | 1,421 | 0 | 0.0% |
| Noise | 10,437 | 10,437 | 0 | 0.0% |

#### FP VerificationClass 變化

| Class | 過濾前 | 過濾後 | 移除 | 移除率 |
|-------|--------|--------|------|--------|
| **Strong** | 1,738 | 1,378 | **360** | **20.7%** |
| Weak | 1,153 | 1,138 | 15 | 1.3% |
| Subclone | 232 | 232 | 0 | 0.0% |
| Noise | 1,699 | 1,699 | 0 | 0.0% |

> [!NOTE]
> **觀察**: 過濾條件精準鎖定 FP Strong 類別，因為這些位點符合「高 AlleleDelta + 低 CramersV + 低 VAF」的假 ASM 特徵。

### 2.2 統計顯著性比較

| 指標 | TP | FP |
|------|----|----|
| **Significant** | 1,861 (6.11%) | 101 (2.27%) |
| **PassedGating** | 9,623 (31.58%) | 1,852 (41.65%) |

> [!WARNING]
> FP 的 PassedGating 比例 (41.65%) 高於 TP (31.58%)，這表示現有的 Gating 機制對 FP 的識別能力有限。

---

## 3. 剩餘 FP 深入分析

### 3.1 最需關注的 FP 子群

| 子群 | 數量 | 說明 |
|------|------|------|
| **Strong** | 1,378 | 30.99% of kept FP |
| **高品質 (QUAL > 0.9)** | 508 | 難以被 QUAL 過濾識別 |
| **高 AD + 高 VAF** | 10 | AD>0.25, VAF≥0.24，未被過濾 |

### 3.2 剩餘 FP Strong 特徵

| 特徵 | 平均值 |
|------|--------|
| AlleleDelta | 0.1269 |
| VAF | 0.1489 |
| QUAL | 0.7064 |

**解釋**: 這些 FP Strong 位點的 AlleleDelta (0.13) 比移除的 FP (0.37) 低很多，且 VAF (0.15) 接近門檻，說明它們的特徵與 TP 更相似，難以區分。

---

## 4. 剩餘位點分布圖

### 4.1 剩餘位點概覽

![Kept Sites Overview](figures/kept_sites_overview.png)

### 4.2 VAF vs AlleleDelta 散布圖

![Kept Sites Scatter](figures/kept_sites_scatter.png)

---

## 5. 高風險案例分析

### 5.1 剩餘的高風險 FP (Strong + QUAL > 0.8)

這些是過濾後仍保留的最危險假陽性：

| Chr | Pos | AlleleDelta | CramersV | VAF | QUAL | NumReads |
|-----|-----|-------------|----------|-----|------|----------|
| chr6 | 98149813 | 0.387 | 0.000 | 0.357 | 0.833 | 11 |
| chr12 | 71376990 | 0.305 | 0.000 | 0.948 | 0.966 | 61 |
| chr8 | 82168704 | 0.258 | **1.000** | 0.253 | 0.980 | 72 |
| chr9 | 41253868 | 0.256 | **0.854** | 0.302 | 0.941 | 84 |
| chr4 | 23100385 | 0.250 | 0.000 | 0.238 | 0.962 | 18 |

**為何未被過濾**:
- chr6:98149813: VAF = 0.357 > 0.24 門檻
- chr12:71376990: VAF = 0.948 遠高於門檻
- chr8:82168704: CramersV = 1.0 > 0.05 門檻
- chr9:41253868: CramersV = 0.854 > 0.05 門檻

> [!TIP]
> 這些案例顯示：**高 VAF 或 高 CramersV** 會讓位點躲過過濾，未來可考慮調整條件。

### 5.2 良好 TP 案例 (Strong + QUAL > 0.95 + VAF > 0.4)

對照參考，以下是典型的良好 TP 特徵：

| Chr | Pos | AlleleDelta | CramersV | VAF | QUAL | NumReads |
|-----|-----|-------------|----------|-----|------|----------|
| chr8 | 115093788 | -0.009 | 0.0 | 0.417 | 0.993 | 97 |
| chr8 | 7388709 | 0.000 | 0.0 | 0.980 | 0.990 | 42 |
| chr3 | 68729618 | 0.000 | 0.0 | 0.974 | 0.997 | 69 |
| chr9 | 120019872 | 0.000 | 0.0 | 0.982 | 0.997 | 47 |

**特徵**: 良好 TP 的 AlleleDelta 接近 0，VAF 高，QUAL 高。

---

## 6. 統計驗證方法評估

### 6.1 現有方法的問題

| 問題 | 說明 | 影響 |
|------|------|------|
| **Significant 標籤利用不足** | TP/FP Significant 比例僅 6%/2%，大部分位點無法判定 | 限制過濾能力 |
| **PassedGating 反向** | FP PassedGating 比例 (42%) 高於 TP (32%) | Gating 條件可能需調整 |
| **CramersV 分布偏斜** | 大部分位點 CramersV = 0 | 效果量難以區分 TP/FP |

### 6.2 建議改進方向

1. **調整 PassedGating 條件**: 現有條件讓太多 FP 通過，需重新評估
2. **引入更多特徵**: 考慮 HPMergedDelta、Strand Bias、Read Position Bias
3. **組合 QUAL 過濾**: 加入 QUAL < 0.75 可額外移除 ~3,000 FP

---

## 7. 後續篩選策略建議

### 7.1 激進策略

若要進一步減少 FP：

```python
# 更激進的過濾條件
Remove if: (QUAL < 0.75) OR
           (AlleleDelta > 0.20 AND CramersV < 0.05 AND VAF < 0.30)
```

**預期效果**: 
- 可額外移除 ~600 FP
- 但會多誤刪 ~20 TP

### 7.2 保守策略（目前）

```python
# 目前使用的條件
Remove if: AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24
```

**優點**:
- 誤刪率極低 (0.02%)
- F1 提升 +0.0040

### 7.3 完整策略（推薦用於生產）

```python
# 結合 QUAL 的完整過濾
Remove if: (QUAL < 0.75) OR
           (AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24)
```

**效果**: F1 = 0.8439 (最佳)

---

## 8. 檔案索引

### 8.1 資料檔案

| 檔案 | 說明 |
|------|------|
| [tp_kept_after_filter.csv](data/tp_kept_after_filter.csv) | 過濾後剩餘 TP 完整資料 |
| [fp_kept_after_filter.csv](data/fp_kept_after_filter.csv) | 過濾後剩餘 FP 完整資料 |
| [kept_sites_summary.json](data/kept_sites_summary.json) | 統計摘要 |
| [high_risk_fp_cases.csv](data/high_risk_fp_cases.csv) | 高風險 FP 案例 |
| [good_tp_cases.csv](data/good_tp_cases.csv) | 良好 TP 案例 |

### 8.2 圖表檔案
![kept_sites_overview](figures/kept_sites_overview.png)
![before_after_comparison](figures/before_after_comparison.png)
![kept_sites_scatter](figures/kept_sites_scatter.png)

---

## 9. 結論

### 9.1 過濾效果評估

✅ **過濾有效**: 精準移除 20.7% FP Strong，誤刪率僅 0.02%

### 9.2 仍需關注的問題

⚠️ **剩餘 1,378 個 FP Strong**: 這些位點特徵與 TP 相似，現有方法難以區分

### 9.3 建議後續行動

1. **結合 QUAL 過濾**: 可將 F1 從 0.8195 提升至 0.8439
2. **調整 PassedGating**: FP 通過率過高需處理
3. **機器學習模型**: 訓練分類器學習複雜的 TP/FP 邊界

---

**報告生成時間**: 2026-01-23  
**分析腳本**: [analyze_kept_sites.py](analyze_kept_sites.py)
