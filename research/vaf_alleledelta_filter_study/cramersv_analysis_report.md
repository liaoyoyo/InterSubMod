# CramersV 與 HPMergedSig 計算方法分析報告

**專案**: InterSubMod - 甲基化輔助體細胞突變檢測  
**分析日期**: 2026-01-23  
**目的**: 調查 CramersV 二元分布成因與計算方法合理性

---

## 1. HPMergedSig 意義解釋

### 1.1 定義

**HPMergedSig** (Haplotype Merged Significance) 表示 **HP1 與 HP2 之間甲基化差異是否統計顯著**。

```cpp
// 來源: RegionProcessor.cpp:1108
result.hp_merged_sig = hp_ms.merged.significant;
```

### 1.2 計算流程

1. 將 reads 依據 Haplotype 分為 HP1 和 HP2 兩組
2. 計算兩組之間的甲基化距離 (HPMergedDelta)
3. 進行統計檢定 (如 PERMANOVA)
4. 根據 p-value 判斷是否顯著

### 1.3 生物學意義

| HPMergedSig | 含義 |
|-------------|------|
| **True** | HP1 與 HP2 的甲基化模式存在**統計顯著差異** |
| **False** | 兩個 haplotype 的甲基化模式無顯著差異 |

**解讀**:
- HP 差異可能來自 **Allele-Specific Methylation (ASM)**
- 也可能來自 **印記基因 (Imprinted genes)**

---

## 2. CramersV 計算方法分析

### 2.1 計算公式

```cpp
// 來源: MathUtils.cpp:154-207
double MathUtils::cramers_v(const vector<int>& table, int n_rows, int n_cols, bool& is_reliable) {
    // 1. 計算 chi-square 統計量
    double chi2 = chi_square(table, n_rows, n_cols, min_expected);
    
    // 2. 計算樣本數
    int n = accumulate(table.begin(), table.end(), 0);
    
    // 3. 計算 Cramér's V
    int min_dim = min(n_rows, n_cols);
    double v = sqrt(chi2 / (n * (min_dim - 1)));
    
    return clamp(v, 0.0, 1.0);
}
```

**數學公式**:
```
V = sqrt(χ² / (n × (min(r, c) - 1)))

其中:
- χ² = Σ[(O - E)² / E]  (chi-square 統計量)
- n = 總樣本數
- r, c = contingency table 的行列數
```

### 2.2 Contingency Table 結構

```
          Cluster_0    Cluster_1    Cluster_2    ...
Label_0     a00          a01          a02         ...
Label_1     a10          a11          a12         ...
```

- **Rows**: 甲基化 cluster 標籤 (由 UPGMA clustering 決定)
- **Cols**: Binary label (Ref/Alt 或 HP1/HP2)

### 2.3 可靠性判斷 (is_reliable)

```cpp
// MathUtils.cpp:185-200
int low_expected_count = 0;
for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
        double expected = (row_totals[i] * col_totals[j]) / n_d;
        if (expected < 5.0) {
            ++low_expected_count;
        }
    }
}

// 若超過 20% 的 cells 期望值 < 5，則不可靠
double low_expected_ratio = low_expected_count / total_cells;
is_reliable = (low_expected_ratio <= 0.2);
```

---

## 3. CramersV 二元分布原因分析

### 3.1 觀察到的現象

| 統計 | TP | FP |
|------|----|----|
| V = 0 | 93.85% | 97.70% |
| V > 0 | 6.15% | 2.30% |
| V = 1.0 | 1.60% | 0.39% |
| Mean (V>0) | 0.85 | 0.76 |

**問題**: 幾乎沒有中間值 (0.1 ~ 0.5)，呈現二元分布。

### 3.2 根本原因

#### 原因 1: is_reliable = false 時輸出為 0

```cpp
// RegionProcessor.cpp:1092-1094
double v_alt = sig_result.global_alt.cramers_v_reliable ? sig_result.global_alt.cramers_v : 0.0;
double v_hp = sig_result.global_hp.cramers_v_reliable ? sig_result.global_hp.cramers_v : 0.0;
result.cramers_v = max(v_alt, v_hp);
```

> [!IMPORTANT]
> **關鍵發現**: 當 `cramers_v_reliable = false` 時，CramersV 直接被設為 **0.0**。這是造成大量位點 V=0 的主要原因。

#### 原因 2: 可靠性門檻嚴格

```cpp
// 可靠條件: <= 20% 的 cells 有 expected < 5
is_reliable = (low_expected_ratio <= 0.2);
```

| 情境 | 預期 | 結果 |
|------|------|------|
| 少量 reads | 多個 cells expected < 5 | V = 0 (不可靠) |
| 不平衡的 label | 某些 cells expected 很小 | V = 0 (不可靠) |
| 足夠 reads + 平衡 | 大部分 cells expected >= 5 | V = 實際計算值 |

#### 原因 3: 當 min_dim <= 1 時直接返回 0

```cpp
// MathUtils.cpp:169-173
int min_dim = min(n_rows, n_cols);
if (min_dim <= 1) {
    is_reliable = false;
    return 0.0;
}
```

- 若只有 1 個 cluster，CramersV = 0
- 若只有 1 個 label (全是 Ref 或全是 Alt)，CramersV = 0

### 3.3 為何高 V 值 (>0.5) 也很常見？

當 `is_reliable = true` 時，計算出的 V 值往往很高：

1. **小樣本的統計特性**: 樣本量不大時，一旦有差異，chi-square 會很大
2. **完美分離**: 若 cluster 與 label 完美對應，V ≈ 1.0
3. **過濾效應**: 只有足夠數據的位點才會通過可靠性檢查，這些位點往往有明確的關聯

---

## 4. 方法合理性評估

### 4.1 合理的部分

| 設計 | 評估 |
|------|------|
| 使用 Cramér's V 作為效果量 | ✅ 合理，是標準的列聯表關聯性測度 |
| 設置可靠性門檻 | ✅ 合理，chi-square 在小樣本時不可靠 |
| 結合 p-value 與 V 進行 gating | ✅ 合理，避免僅依賴 p-value |

### 4.2 可能的問題

| 問題 | 說明 | 影響 |
|------|------|------|
| **二元輸出** | V 不是 0 就是高值，失去區分能力 | 無法利用 V 的連續性 |
| **依賴 cluster 數量** | n_rows 由 clustering 決定，不穩定 | 不同位點的 V 不可比 |
| **低覆蓋度直接為 0** | 很多真實差異被忽略 | 可能漏掉有意義的信號 |
| **20% 門檻可能太嚴** | 很多低中等覆蓋的位點被排除 | 大量信息丟失 |

### 4.3 待確認的問題

> [!WARNING]
> 以下問題需要進一步確認：

1. **Cluster 數量如何決定？**
   - UPGMA clustering 的截斷標準是什麼？
   - 不同位點的 cluster 數可能差異很大

2. **expected < 5 的 20% 門檻是否太嚴格？**
   - 標準統計建議是「不超過 20% 的 cells < 5 且沒有 cell < 1」
   - 目前實現可能過於保守

3. **為何不輸出原始 V 值？**
   - 即使不可靠，原始 V 值仍可作為參考
   - 可以設置 flag 標記可靠性，而非歸零

---

## 5. 改進建議

### 5.1 短期改進

| 建議 | 實施方式 |
|------|----------|
| **不要歸零** | 即使 is_reliable=false，仍輸出計算值，另設 flag |
| **輸出原始 V** | 新增欄位 CramersV_raw |
| **調整門檻** | 將 20% 放寬至 30% 或 40% |

### 5.2 中期改進

| 建議 | 實施方式 |
|------|----------|
| **使用 bias-corrected Cramér's V** | 對小樣本進行校正 |
| **考慮樣本量加權** | 報告 (V, n) 對，讓下游可以自行判斷 |
| **多種效果量** | 同時報告 Phi, V, Cohen's w |

### 5.3 長期改進

| 建議 | 實施方式 |
|------|----------|
| **機器學習特徵** | 將原始 contingency table 作為特徵輸入 |
| **不依賴 clustering** | 直接計算 read-level 的關聯性 |

---

## 6. 與過濾條件的關係

### 6.1 目前過濾條件

```python
Remove if: AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24
```

### 6.2 CramersV < 0.05 的含義

由於二元分布，`CramersV < 0.05` 實際上等價於：
- V = 0 (93.85% of TP, 97.70% of FP)
- **或** 0 < V < 0.05 (極少數)

因此這個條件實際上是在篩選：
> 「高 AlleleDelta + 低 VAF + **統計不可靠或無關聯**」的位點

### 6.3 對過濾效果的影響

| 情境 | CramersV | 過濾結果 |
|------|----------|----------|
| 高 AD + 低 VAF + V=0 | < 0.05 | ✅ 被過濾 (假設為 FP) |
| 高 AD + 低 VAF + V>0.5 | >= 0.05 | ❌ 保留 (可能是真 ASM 或 FP) |

> [!TIP]
> **建議**: 若能改進 CramersV 計算，輸出連續值，可以更精確地設置過濾門檻。

---

## 7. 結論

### 7.1 CramersV 二元分布的核心原因

1. **is_reliable = false 時歸零** — 這是主要原因
2. **可靠性門檻嚴格** — 20% cells < 5 的標準排除了大量位點
3. **min_dim <= 1 時直接返回 0** — 單 cluster 或單 label 情況

### 7.2 HPMergedSig 解釋

- 表示 HP1 與 HP2 之間甲基化差異是否統計顯著
- True 表示存在 haplotype-specific methylation

### 7.3 方法是否合理

**整體合理**，但有改進空間：
- 建議不要歸零，而是設置可靠性 flag
- 考慮使用 bias-corrected Cramér's V
- 可同時輸出多種效果量供下游選擇

---

**報告生成時間**: 2026-01-23  
**分析來源**: `src/core/MathUtils.cpp`, `src/core/GlobalTest.cpp`, `src/core/RegionProcessor.cpp`
