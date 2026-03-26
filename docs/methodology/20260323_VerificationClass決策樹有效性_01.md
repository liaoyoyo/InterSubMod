<!--
建立時間: 2026-03-23
問題類型: C++ 方法 | 統計分類
影響 track: TO | 兩者
狀態: pending_decision
-->

# VerificationClass 決策樹有效性 審查文件

## 問題描述

ISM 的最終輸出 `VerificationClass`（Strong/Weak/Subclone/Noise）是 ISM 對該 region 意義的核心判斷。但數據顯示：

1. **`Significant=True` 的 AUROC 僅 0.511**（TP vs FP 二元分類，3樣本合計）
2. **TP 中 Noise 佔 46.6%**（近半數已知真陽性被認為無意義）
3. **Strong FP 佔 FP 的 16.2%**（最高信心的錯誤判斷）
4. **TP 和 FP 的 class 分佈差異不大**（ISM 目前無法有效區分 TP 和 FP）

這些問題指向 VerificationClass 決策樹邏輯本身存在設計問題。

---

## 量化影響

### VerificationClass 分佈（3樣本合計）

| Class | TP 比例 | FP 比例 | TP/(TP+FP) | 含義 |
|-------|---------|---------|------------|------|
| Strong | 22.1% | 16.2% | ~58% | 最高信心，但仍有 16.2% 是 FP |
| Weak | 28.5% | 24.7% | ~54% | 稍低信心 |
| Subclone | 2.8% | 3.2% | ~47% | 亞克隆，FP 比例稍高 |
| Noise | **46.6%** | **55.8%** | ~46% | 大量 TP 被歸類為無意義 |

**Strong class 精準率約 58%**：意味著 Strong 判斷中有約 42% 是 FP。

### passed_gating 效果

- TP: ClusterPassedGating=True 27.8%（這些能進入後續分析）
- FP: ClusterPassedGating=True 22.0%
- PassedGating 的 AUROC = 0.529（微弱）

**結論**：PassedGating 雖能排除大量 FP，但也同時排除了大量 TP（72.2%）。

---

## VerificationClass 決策樹解析

**程式碼位置**：`src/core/SignificanceAnalyzer.cpp:290-420`

決策邏輯概要（需讀 C++ 確認）：
```
passed_gate=True & sig=True & CramersV≥0.3 → Strong
passed_gate=True & sig=True & CramersV<0.3 → Weak
passed_gate=True & sig=True & label=subclone → Subclone
else → Noise
```

**關鍵問題**：
1. `CramersV` 的閾值 0.3 是如何決定的？是否有數據支持？
2. 當 passed_gate=False 時，直接歸 Noise，即使 AllelePermanovaF 很高
3. 46.6% 的 TP 走到 Noise 路徑：
   - 是因為 passed_gate=False（72.2% TP 不過關）？
   - 還是過關後 CramersV 太低？

---

## 修改選項

### 方案 A：不修改（補文件說明限制）
- **理由**：現況 ISM F1 delta 接近 0，說明分類錯誤不影響最終 calling
- **後續**：文件說明 Noise class 不代表「無信號」，只代表「信號未達閾值」
- **F1 影響**：0

### 方案 B：分析決策樹各分支，找出 TP 被錯誤歸 Noise 的原因
- **目標**：不改邏輯，但量化每個分支的貢獻
- **實作**：
  1. 從 `label_first_metrics.tsv` 分析 Noise TP 的特徵分佈
  2. 確認：是否 passed_gate=False 就能解釋 46.6% Noise TP？
  3. 確認：LabelAllelePermanovaF 高但仍歸 Noise 的 case 有多少？
- **產出**：決策樹流量圖（數字版）
- **成本**：Python 分析，無需改 C++

### 方案 C：放寬 passed_gating 閾值（調整 global_alt.p_value 閾值）
- **程式碼位置**：`src/core/GlobalTest.cpp`（待確認）
- **實作**：將 global_p 閾值從 0.1 調整到 0.15 或 0.2
- **風險**：更多 FP 也會通過 gating，需看 precision/recall 曲線
- **F1 影響**：不確定，需測試

### 方案 D：引入 LabelAllelePermanovaF 作為獨立強化條件
- **邏輯**：即使 passed_gate=False，若 LabelAllelePermanovaF ≥ 10 也視為有信號
- **理由**：LabelAllelePermanovaF AUROC=0.569 是目前最好的特徵
- **實作位置**：`src/core/SignificanceAnalyzer.cpp`（VerificationClass 決策樹）
- **預估 F1 影響**：可能 +0.001~+0.005

---

## 建議下一步（方案 B 優先）

**先執行方案 B**（純數據分析，不改 C++）：

```python
# 分析 Noise TP 的原因
noise_tp = combined[(combined['label']=='TP') & (combined['VerificationClass']=='Noise')]
print(noise_tp['ClusterPassedGating'].value_counts(normalize=True))  # 是因為 passed_gate=False？
print(noise_tp['LabelAllelePermanovaF'].describe())  # 這些 TP 的 Allele F 分佈
```

產出文件：`docs/methodology/20260323_VerificationClass決策樹流量分析_02.md`

---

## 驗收標準

若選 C/D：
- [ ] test-quick 通過
- [ ] test-full F1 delta ≥ +0.001（有實質改進）
- [ ] Strong FP 比例不顯著增加
- [ ] 跨 3 個樣本確認

---

## 用戶決策

**選擇**：[ ] A / [ ] B（先分析）/ [ ] C / [ ] D
**日期**：
**理由**：
