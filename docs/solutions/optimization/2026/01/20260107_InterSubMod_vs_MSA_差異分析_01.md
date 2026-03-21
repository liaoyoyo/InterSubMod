<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# InterSubMod vs MethylSomaticAnalysis (MSA) 結果差異分析

基於 `InterSubMod` (Round 6) 的程式邏輯與一般方法學的差異，以下分析兩者結果不同的原因。

## 1. 程式邏輯與篩選條件差異

最主要的差異來自於 `InterSubMod` 剛剛經過的 6 輪深度優化，包含了極為嚴格的過濾條件，而 `MSA` (假設為基準版本) 可能未包含這些針對性的篩選。

| 特徵 | InterSubMod (Round 6) | MSA (推測基準) | 影響 |
|:---|:---|:---|:---|
| **顯著性定義** | **P ≤ 0.05** AND **V ≥ 0.1** AND **Depth ≥ 30** | 可能僅 P ≤ 0.05 | InterSubMod 會排除大量低效應量 (V < 0.1) 的位點，顯著數大幅減少 |
| **FP 過濾** | **LabelDelta > 0.3** 標記為 `SuggestFilter` | 無此機制 | InterSubMod 能識別並標記 3.96% 的高風險 FP 位點 |
| **距離演算法** | **Bernoulli / Weighted Hamming** | 可能為 Euclidean 或 Hamming | 距離計算不同導致聚類結果不同，影響 Local Test P-value |
| **Gating 機制** | 雙重閾值 (P + V) | 可能僅 P 值 | InterSubMod 提早排除干擾訊號 |

## 2. 結果差異的具體表現

### A. 顯著位點數量 (Significant Sites)
- **InterSubMod**: 由於加入了 `Cramér's V ≥ 0.1` 和 `Depth ≥ 30` 的硬性閾值，**顯著位點數量會顯著少於 MSA**。
- **MSA**: 如果僅依賴 P-value，會包含大量 P 值顯著但效應量極低 (V ≈ 0.0) 的 "噪音" 位點。

### B. 真假陽性比 (TP/FP Ratio)
- **InterSubMod**: 經過 Round 2 優化，TP/FP 比從 0.70 提升至 **2.92**。結果會更乾淨，FP 更少。
- **MSA**: 可能維持在原始水平 (TP/FP ≈ 0.70)，即 FP 數量甚至可能超過 TP。

### C. 邊緣位點的處理
- 對於深度較低 (20-30 reads) 或甲基化差異微弱 (V=0.05-0.1) 的位點，**InterSubMod 會判定為非顯著**，而 MSA 可能判定為顯著。

## 3. 結論與建議

**差異原因總結：**
> `InterSubMod` (Round 6) 是一個**高度特化且保守**的版本，旨在最大化 Precision (減少 FP)。而 `MSA` 可能是一個**廣泛且敏感**的版本，旨在最大化 Recall (減少 FN)。

**如何解讀：**
1. 如果你需要**高信心**的候選位點 → 採用 `InterSubMod` (並過濾 SuggestFilter=True)。
2. 如果你需要**探索性**的分析，不介意高 FP → `MSA` 可能保留了更多微弱訊號。
3. 兩者的差異是**策略性選擇**的結果 (Precision vs Recall 的權衡)。

---

建議檢查 MSA 的輸出設定檔，確認其 P-value 閾值、距離演算法與是否有效應量過濾，即可證實上述分析。
