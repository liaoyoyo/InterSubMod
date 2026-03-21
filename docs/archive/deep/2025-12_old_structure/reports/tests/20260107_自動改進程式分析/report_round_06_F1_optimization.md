# 2026-01-07 F1-score 最佳化分析報告

## 結論摘要

| 指標 | 原始值 | 最佳化後 | 變化 |
|:---|:---:|:---:|:---|
| TP | 30,490 | 30,377 | -113 (0.37%) |
| FP | 4,842 | 4,650 | **-192 (3.96%)** |
| FN | 8,960 | 9,073 | +113 |
| **F1-score** | **0.8154** | **0.8157** | **+0.0003 ✅** |

**過濾比例**: FP:TP = 1:0.59 (優於基準 1:0.69)

---

## 最佳策略

### 使用 `SuggestFilter` 欄位

```
過濾條件: LabelDelta > 0.3
```

當 `SuggestFilter = True` 時，建議移除該 SNV 位點。

### 原理

**LabelDelta** 代表標籤驅動的甲基化距離差異：
- 高 LabelDelta (> 0.3) 表示甲基化分群與遺傳標籤 (HP/Allele) 有強烈關聯
- 分析發現 **FP 位點更容易有高 LabelDelta**
- 這可能因為 FP 位於系統性甲基化異質區域 (重複序列、結構變異等)

### 程式碼修改

**檔案**: `src/core/RegionProcessor.cpp`

```cpp
// NEW: Suggest filtering sites with LabelDelta > 0.3 (F1 optimization)
bool suggest_filter = (r.label_delta > 0.3);
```

**CSV 新增欄位**: `SuggestFilter`

---

## 各特徵過濾效果分析

| 特徵 | 閾值 | FP過濾率 | TP損失 | FP:TP 比 | F1 變化 |
|:---|:---|:---:|:---:|:---:|:---:|
| Cramér's V ≥ 0.1 | 0.1 | 97.8% | 93.9% | 0.16 | -0.73 ❌ |
| HeuristicScore ≥ 1 | 1.0 | 54.1% | 68.4% | 0.13 | -0.44 ❌ |
| **LabelDelta > 0.3** | **0.3** | **3.96%** | **0.37%** | **1.69** | **+0.0003 ✅** |
| NumReads ≥ 30 | 30 | 10.1% | 2.9% | 0.56 | -0.01 ❌ |

**只有 LabelDelta > 0.3 的 FP:TP 比 (1.69) 大於 1，適合用於過濾。**

---

## 使用方式

### 在下游工具中過濾

```python
import pandas as pd

df = pd.read_csv("significance_summary.csv")
# 過濾掉 SuggestFilter = True 的位點
filtered_df = df[df['SuggestFilter'] == False]
```

### 預期效果

- F1-score: 0.8154 → 0.8157 (+0.0003)
- 這是目前甲基化分析能達到的最佳改善

---

## 限制與說明

1. **改善幅度有限** (+0.0003)：因為只有 0.37% TP 和 3.96% FP 符合過濾條件
2. **大多數位點無法區分**：~96% 的 FP 和 ~99.6% 的 TP 的 LabelDelta ≤ 0.3
3. **甲基化異質性不是 SNV 真假的主要特徵**：大部分 TP 沒有明顯的 ASM

---

## 輸出檔案位置

| 檔案 | 路徑 |
|:---|:---|
| 分析結果 | `output/bip8_disk_output/20260107_all-with-w5000_6/` |
| significance_summary.csv | 包含新增的 `SuggestFilter` 欄位 |
