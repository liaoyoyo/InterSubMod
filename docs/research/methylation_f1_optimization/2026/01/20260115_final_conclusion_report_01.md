<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 甲基化 F1 研究完整結論報告

**研究時間**: 2026-01-13 ~ 2026-01-14  
**計劃來源**: [20260113_methylation_significance_plan_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/01/20260113_methylation_significance_plan_01.md)

---

## 資料來源

> [!NOTE]
> **VCF 來源**: ClairS 0.4.1 **longphase-s** 模式處理後的結果
> 
> longphase-s 整合了 haplotag (單倍型標記) 資訊，理論上可進一步過濾不合理的呼叫，
> 但本研究發現即使經過 longphase-s 處理，**仍需額外的 QUAL/AF 後過濾**才能達到最佳效果。

| 檔案 | 位置 | 數量 |
|:---|:---|---:|
| TP VCF | `data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz` | 30,490 |
| FP VCF | `data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz` | 4,842 |
| FN VCF | `data/vcf/HCC1395/pileup/filtered_snv_fn.vcf.gz` | 8,957 |

---

## 核心結論

> [!IMPORTANT]
> ### 主要發現
> 1. **VCF 原始特徵 (QUAL、AF) 遠優於甲基化特徵**
> 2. **最佳策略: QUAL < 0.8**，可將 F1 提升 **+2.70%**
> 3. **甲基化 CramersV 區分能力弱** (AUC = 0.52)，但 NumReads 有效 (AUC = 0.63)
> 4. **chr8 染色體** 佔 FP 高聚集區域的 **66.3%**，需特別注意

---

## 驗證結果

### F1-score 計算確認

所有分析均基於 **SEQC2 金標準** 驗證：

| 項目 | 數值 | 說明 |
|:---|---:|:---|
| SEQC2 金標準位點 | 39,447 | HCC1395 高信心體細胞突變 |
| True Positives (TP) | 30,490 | **100% 與 SEQC2 交集** |
| False Positives (FP) | 4,842 | 呼叫但不在 SEQC2 中 |
| False Negatives (FN) | 8,957 | 在 SEQC2 但未被呼叫 |
| **Baseline F1** | **0.8155** | 明恩 SNV 結果 |

---

## 特徵區分能力 (AUC)

### VCF 特徵 vs 甲基化特徵

| 特徵類型 | 特徵名稱 | AUC | 判定 |
|:---|:---|:---:|:---|
| **VCF** | QUAL | **0.9668** | ✅ 極佳 |
| **VCF** | AF | **0.9235** | ✅ 極佳 |
| 甲基化 | NumReads | 0.6303 | ✅ 有效 |
| 甲基化 | GlobalP | 0.5614 | ⚠️ 弱 |
| 甲基化 | CramersV | 0.5194 | ❌ 無效 |
| 甲基化 | HeuristicScore | 0.4437 | ❌ 無效 |

**結論**: VCF 特徵的區分能力是甲基化特徵的 **2-3 倍**。

---

## 過濾策略效果

### 最佳策略比較

| 策略 | 新 F1 | 改善 | 移除 FP | 誤刪 TP |
|:---|---:|---:|---:|---:|
| **QUAL < 0.8** | **0.8424** | **+2.70%** | 75.6% | 3.0% |
| QUAL<0.6 OR AF<0.12 | 0.8412 | +2.57% | 63.1% | 1.8% |
| AF < 0.12 | 0.8361 | +2.06% | 53.6% | 2.3% |
| QUAL<0.5 OR AF<0.10 | 0.8337 | +1.82% | 53.7% | 1.0% |

---

## 甲基化特徵特殊發現

### 1. 顯著性的價值

| 位點類型 | 數量 | TP 比例 |
|:---|---:|---:|
| 甲基化顯著 (Significant=True) | 1,961 | **94.8%** |
| 甲基化非顯著 (Significant=False) | 33,337 | 85.8% |

**解讀**: 若某位點呈現甲基化顯著性，則該位點為 TP 的機率高達 94.8%，比非顯著位點高出約 9 個百分點。

### 2. NumReads 在不同品質區間的效果

| VCF 品質區間 | NumReads AUC | 說明 |
|:---|:---:|:---|
| 低 AF (<0.15) | **0.8403** | ← 價值最高 |
| 中 QUAL (0.5-0.8) | 0.7835 | 次高 |
| 高 QUAL (≥0.8) | 0.5553 | 無效 |

**解讀**: NumReads 在低品質區域的區分能力特別強，可作為邊界案例的輔助判斷依據。

### 3. 可「救援」的風險 TP

- **風險 TP** (QUAL<0.8 OR AF<0.15): 1,503 個
- 其中有 **高甲基化信號** (CramersV>0.3 OR HScore>5): 334 個 (22.2%)
- 這些 TP 可透過甲基化特徵保留，避免被過度過濾

---

## 區域特徵分析

### FP 高聚集區域 (ClusterCount > 20)

| 染色體 | FP 數量 | 佔比 |
|:---|---:|---:|
| **chr8** | 377 | **66.3%** |
| chr9 | 107 | 18.8% |
| chr14 | 44 | 7.7% |
| chr2 | 41 | 7.2% |

**解讀**: chr8 染色體上存在異常高的 FP 聚集，可能與該區域的基因組特性或定序偏差有關，建議在分析時特別留意。

---

## 圖表說明

### Figure 1: QUAL Distribution (TP vs FP)

![QUAL Distribution](../../figures/2026/01/phase1_plots/exp1_1_qual_distribution.png)

**說明**: 此圖顯示 VCF QUAL 分數在 TP 和 FP 之間的分布差異。
- 左圖為密度分布：TP（綠色）集中在高 QUAL 區域 (0.9-1.0)，FP（紅色）分布較廣
- 右圖為箱形圖：清楚顯示 TP 中位數約 0.99，FP 中位數約 0.66
- **結論**: QUAL 分數可有效區分 TP/FP (AUC = 0.9668)

---

### Figure 2: AF Distribution Analysis

![AF Distribution](../../figures/2026/01/phase1_plots/exp1_2_af_distribution.png)

**說明**: 此圖分析 Allele Frequency (AF) 的分布特性。
- 左圖：TP 的 AF 集中在 0.4-0.6（雜合子特徵），FP 集中在低 AF 區域
- 中圖：KDE 平滑曲線更清楚顯示兩者差異
- 右圖：按 AF 範圍計算 TP 比例，低 AF (<0.1) 區域 TP 比例最低
- **結論**: 低 AF 的位點較可能是 FP

---

### Figure 3: Grid Search Results (QUAL and AF Thresholds)

![Grid Search](../../figures/2026/01/phase2_plots/exp2_1_2_grid_search.png)

**說明**: 此圖顯示不同過濾閾值對 F1 的影響。
- 左圖 (QUAL)：閾值 0.8 達到最佳 F1，再高則召回率下降過多
- 右圖 (AF)：閾值 0.12 效果較好
- 紅線為 baseline F1 (0.8155)，綠色條形表示有效策略
- **結論**: QUAL < 0.8 是最佳單一策略

---

### Figure 4: Regional Clustering Analysis

![Regional Clustering](../../figures/2026/01/phase3_plots/regional_clustering_analysis.png)

**說明**: 此圖分析位點的區域聚集特性。
- 左上：FP 的聚集度明顯高於 TP
- 右上：高聚集區域 (>20) 的 TP 比例驟降至約 30%
- 左下：各染色體的 TP/FP 分布
- 右下：chr8 的平均聚集度異常高
- **結論**: 高聚集區域需警惕，chr8 尤其明顯

---

### Figure 5: Methylation Feature Interaction Effects

![Interaction Effects](../../figures/2026/01/phase4_plots/interaction_effects_analysis.png)

**說明**: 此圖探索甲基化特徵與 VCF 特徵的交互效應。
- 左上：QUAL × CramersV 熱圖，顏色表示 TP 比例
- 右上：AF × 顯著性熱圖，顯著位點的 TP 比例較高
- 左下：各特徵 AUC 比較，VCF 特徵明顯優於甲基化
- 右下：NumReads vs HeuristicScore 散點圖
- **結論**: 交互特徵未能提升區分能力

---

### Figure 6: Conditional Methylation Value

![Conditional Value](../../figures/2026/01/phase4_plots/conditional_methylation_value.png)

**說明**: 此圖分析甲基化特徵在不同 VCF 品質區間的價值。
- 左圖熱圖：NumReads 在低 AF 區域 AUC 高達 0.84
- 右圖：各區域最佳甲基化特徵，低 AF/QUAL 區域 NumReads 最有效
- **結論**: 甲基化特徵在邊界案例中有輔助價值

---

### Figure 7: Hybrid Strategy Comparison

![Hybrid Strategy](../../figures/2026/01/phase4_plots/hybrid_strategy_comparison.png)

**說明**: 此圖比較混合過濾策略的效果。
- 各策略相對於 baseline 的 F1 改變百分比
- QUAL<0.8 only 仍是最佳 (+2.70%)
- 加入甲基化「救援」機制的策略略低 (+2.67%)
- **結論**: 純 VCF 過濾已是最佳，甲基化僅供參考

---

## 最終建議

### 主要策略

```
過濾條件: QUAL < 0.8
預期效果: F1 從 0.8155 提升至 0.8424 (+2.70%)
移除 FP: 75.6% (3,659 個)
誤刪 TP: 3.0% (922 個)
```

### 輔助參考

1. **顯著性標記**: 若位點 Significant=True，可信度較高 (94.8% 為 TP)
2. **NumReads**: 在低 AF 邊界案例中可作為額外判斷依據
3. **chr8 警示**: 該區域 FP 聚集程度最高，需額外審視

### 不建議

- ❌ 使用 CramersV 作為過濾條件 (AUC 僅 0.52)
- ❌ 使用 HeuristicScore 作為過濾條件 (AUC 僅 0.44)
- ❌ 僅依賴甲基化顯著性進行過濾 (會大幅降低 F1)

---

## 研究成果總覽

| 項目 | 數量 |
|:---|---:|
| 分析報告 | 4 份 |
| 視覺化圖表 | 13 張 |
| 分析腳本 | 4 個 |
| 測試策略 | 26+ 個 |

### 檔案位置

- 報告: `docs/research/methylation_f1_optimization/2026/01/`
- 圖表: `docs/research/methylation_f1_optimization/figures/2026/01/phase*_plots/`
- 腳本: `scripts/analysis/legacy/methylation_f1_optimization/`
