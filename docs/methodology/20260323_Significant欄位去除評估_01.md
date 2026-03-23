<!--
建立時間: 2026-03-23
目標: 確認 Significant 欄位是否過擬合、設計替代標準
處理範圍: 7 個樣本 + HCC1395 pileup n=99（significance_summary.csv）
關聯檔案:
  - docs/methodology/20260323_VerificationClass決策樹跨樣本量化_01.md (Section 4 數據來源)
  - src/core/RegionProcessor.cpp (Significant 欄位定義)
  - scripts/pipeline/steps/03_filter_analysis.py
-->

# `Significant` 欄位去除評估

**狀態：分析完成 — 建議移除 `Significant`，改用 `VC != Noise` 作為 downstream 標準**

---

## 背景

`Significant` 是 C++ 輸出至 `significance_summary.csv` 的欄位，由 `RegionProcessor.cpp` 的 Python 層判斷邏輯寫入。
目前在 `03_filter_analysis.py` 與部分下游分析腳本中，`Significant=True` 被用作「有信號」的篩選標準。

---

## 定義比較

| 標準 | 定義 | 特徵 |
|------|------|------|
| `Significant=True` | 更嚴格的複合條件（詳見 RegionProcessor.cpp） | 極保守，捕獲少數高信心位點 |
| `VC != Noise` | Strong OR Weak OR Subclone | 有任何結構性信號的全集 |
| `VC != Noise` 中的 LOH_Subclone 子集 | + `LOH_Subtype = 'LOH_Subclone'` | 可靠 Subclone 分層 |

---

## 跨樣本數據比較

| 樣本 | Sig=True TP% | Sig=True FP% | Sig=True Prec | VC!=Noise TP% | VC!=Noise FP% | VC!=Noise Prec |
|------|-------------|-------------|--------------|--------------|--------------|---------------|
| HCC1395 | **6.2%** | 0.6% | 99.8% | **65.8%** | 79.4% | 97.5% |
| HCC1395_DORADO | **11.9%** | 10.4% | 99.3% | **65.3%** | 41.2% | 99.5% |
| COLO829 | **1.2%** | 1.1% | 94.8% | **33.1%** | 34.8% | 93.7% |
| H1437 | **11.8%** | 12.5% | 100.0% | **66.3%** | 75.0% | 100.0% |
| H2009 | **3.6%** | 0.0% | 100.0% | **41.6%** | 10.5% | 100.0% |
| HCC1937 | **4.9%** | 1.0% | 99.7% | **57.6%** | 10.8% | 99.7% |
| HCC1954 | **6.8%** | 0.0% | 100.0% | **67.5%** | 55.2% | 99.9% |
| HCC1395_pileup_n99 | **6.6%** | 2.2% | 95.1% | **66.6%** | 66.1% | 86.4% |

### 關鍵觀察

1. **`Significant=True` 極度保守**：跨樣本僅捕獲 **1.2~11.9%** 的 TP（平均 ≈6%）
2. **`VC != Noise` 更合理**：捕獲 **33~68%** TP，且 Precision 幾乎相同（93.7~99.9%）
3. **Precision 差距小**：`Significant=True` Precision = 94.8~100%，`VC!=Noise` Precision = 86.4~100%
   → Precision 差距 ≈ 0~9%，但 TP 召回率差距高達 10-60%

---

## Significant 欄位的問題

### 問題 1：捕獲率極低

在 pileup n=99（最完整數據）中：
- `Significant=True`：2017/30475 TP（6.6%），捕獲 5.0% pileup TP
- `VC != Noise`：20287/30475 TP（66.6%）

**`Significant=True` 拋棄了 93% 的 TP 信號。**

### 問題 2：過擬合 pileup 數據格式

`Significant` 欄位的定義包含多重嚴格條件，對 pileup 模式（n=99 reads）效果尤差：
- COLO829 只有 1.2% TP 被標為 Significant，但 33% TP 有結構信號（VC!=Noise）

### 問題 3：FP 辨別力不優

| 樣本 | Sig=True FP% | VC!=Noise FP% | 差異 |
|------|-------------|--------------|------|
| HCC1395 | 0.6% | 79.4% | VC!=Noise 包含更多 FP（但也包含更多 TP） |
| COLO829 | 1.1% | 34.8% | 類似 |

`Significant=True` 雖然 FP 佔比低，但這是因為它捕獲的 TP 太少，整體 Precision 並無顯著提升。

---

## 替代方案

### 推薦：`VC != Noise` 作為主要標準

```python
# 推薦替代
df['has_signal'] = df['VerificationClass'] != 'Noise'

# 更精細分層（可選）
df['signal_strong'] = df['VerificationClass'].isin(['Strong', 'Weak'])
df['signal_loh_subclone'] = (
    (df['VerificationClass'] == 'Subclone') &
    (df['LOH_Subtype'] == 'LOH_Subclone')
)
```

### 降級使用場景（保留 Significant 的場景）

若下游需要極高精確度（如論文展示「high-confidence 位點」）：
```python
# Significant=True 仍可用於「極信心子集」展示
df['high_confidence'] = df['Significant'] == True
# 但不應作為主要 pipeline filter
```

---

## 對 `03_filter_analysis.py` 的影響

目前 `scripts/pipeline/steps/03_filter_analysis.py` 的 `should_filter_variant()` 函式**不直接讀取 `Significant` 欄位**（已確認）。
`Significant` 主要被用於：
- 下游 Python 分析腳本中的手動篩選
- 視覺化輸出的高亮標記

**→ 移除對 `Significant` 的依賴不需修改 `03_filter_analysis.py`。**

---

## 判斷結論

| 問題 | 結論 |
|------|------|
| `Significant=True` 是否過擬合？ | **是** — 捕獲率 1.2~11.9%，遺漏 88~99% TP 信號 |
| 替代標準是否存在？ | **是** — `VC != Noise` 捕獲率 33~68%，Precision 相近 |
| 是否需要修改 C++？ | **否** — 欄位保留，只修改下游 Python 分析邏輯 |
| 是否需要修改 `03_filter_analysis.py`？ | **否** — 該腳本不讀取 `Significant` |
| 建議行動 | 在所有分析腳本中以 `VC != Noise` 取代 `Significant=True` 作為有信號標準 |

**判決：Approved — 移除 `Significant` 作為主要篩選標準，改用 `VC != Noise`。**
