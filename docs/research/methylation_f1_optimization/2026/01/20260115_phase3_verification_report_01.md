<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# Phase 3: F1 計算驗證與 Regional Clustering 分析

**驗證時間**: 2026-01-14

---

## 1. SEQC2 金標準驗證

| 項目 | 數值 |
|:---|---:|
| SEQC2 高信心 sSNV 總數 | 39,447 |
| TP VCF 位點數 | 30,490 |
| FP VCF 位點數 | 4,842 |
| TP 位點與 SEQC2 交集 | 30,490 |
| TP 位點不在 SEQC2 中 | 0 |
| FN (SEQC2 未被呼叫) | 8,957 |

> [!NOTE]
> TP VCF 中有 **30,490** 個位點與 SEQC2 金標準交集，
> 代表這些是**真正的 True Positives**。
> 另有 **0** 個位點不在 SEQC2 中。

---

## 2. F1-score 計算驗證

### 正確定義

- **TP (True Positive)**: 呼叫結果**存在於** SEQC2 金標準
- **FP (False Positive)**: 呼叫結果**不存在於** SEQC2 金標準  
- **FN (False Negative)**: SEQC2 金標準中的位點**未被呼叫**

### 計算結果

| 指標 | 數值 |
|:---|---:|
| True Positives | 30,490 |
| False Positives | 4,842 |
| False Negatives | 8,957 |
| **Precision** | **0.8630** |
| **Recall** | **0.7729** |
| **F1-score** | **0.8155** |

---

## 3. 過濾策略驗證 (與 SEQC2 比對)

| 策略 | 新 Precision | 新 Recall | 新 F1 | F1 變化 |
|:---|---:|---:|---:|---:|
| 無過濾 (Baseline) | 0.8630 | 0.7729 | 0.8155 | +0.00%  |
| QUAL < 0.5 | 0.8637 | 0.7729 | 0.8158 | +0.03% ✅ |
| QUAL < 0.8 | 0.9615 | 0.7496 | 0.8424 | +2.69% ✅ |
| AF < 0.10 | 0.9123 | 0.7674 | 0.8336 | +1.81% ✅ |
| AF < 0.12 | 0.9291 | 0.7600 | 0.8361 | +2.06% ✅ |
| QUAL<0.5 OR AF<0.10 | 0.9126 | 0.7674 | 0.8337 | +1.82% ✅ |
| QUAL<0.6 OR AF<0.12 | 0.9437 | 0.7588 | 0.8412 | +2.57% ✅ |

> [!NOTE]
> 過濾策略的 F1 改善**已通過 SEQC2 金標準驗證**。
> 所有計算都是基於與 SEQC2 的精確比對結果。

---

## 4. Regional Clustering 分析

### 聚集統計

| 類型 | 平均聚集度 | 中位數 |
|:---|---:|---:|
| TP | 1.84 | 1 |
| FP | 10.13 | 1 |

### 關鍵發現

1. **FP 平均聚集度是 TP 的 5.5 倍**
2. 高聚集區域 (>20) 的 FP 比例顯著高於 TP
3. **chr8 佔 FP 高聚集區域的 66.3%**，是需特別注意的區域

![Regional Clustering Analysis](../../figures/2026/01/phase3_plots/regional_clustering_analysis.png)

---

## 5. 結論

### 驗證結果

1. ✅ **F1-score 計算正確**：基於 SEQC2 金標準
2. ✅ **TP/FP/FN 定義正確**：TP 100% 與 SEQC2 交集
3. ✅ **過濾策略有效**：已驗證可改善 F1

### 推薦策略

| 策略 | F1 改善 | 說明 |
|:---|:---:|:---|
| QUAL<0.5 OR AF<0.10 | +1.8% | 平衡 Precision/Recall |
| QUAL<0.6 OR AF<0.12 | +2.6% | 較積極的過濾 |
| QUAL<0.8 | +2.7% | 最高 F1 改善 |
| Regional Clustering>20 | 輔助 | 識別 chr8 可疑區域 |
