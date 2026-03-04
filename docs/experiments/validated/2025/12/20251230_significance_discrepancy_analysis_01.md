<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 顯著性判斷差異分析報告
**日期**: 2025-12-30
**目標**: 分析 TP (True Positive) 與 FP (False Positive) 數據集中，算法判斷結果與預期（或視覺判斷）不符的案例，並找出原因。

## 1. 摘要
通過對比 `20251229_w5000_WED_TP` 與 `20251229_w5000_WED_FP` 的結果，我們發現目前的顯著性判斷邏輯（`Significant` column）存在以下主要問題：
1.  **FP 誤判為顯著 (False Positive)**: 算法過度依賴 Fisher's Exact Test 的 P-value，即使在數據稀疏、關聯強度（Cramér's V）不可靠的情況下，仍可能因 P < 0.05 而被標記為顯著。
2.  **TP 誤判為不顯著 (False Negative)**: 算法採用嚴格的 P <= 0.05 閾值，導致部分效應量（Effect Size）極強（V > 0.6）但 P-value 處於邊緣（如 0.052）的案例被剔除。

## 2. 分析方法
- **採樣對象**:
    - **FP 數據集中的顯著案例**: 預期應為不顯著，但算法判讀為 `True`。
    - **TP 數據集中的不顯著案例**: 預期應為顯著，但算法判讀為 `False`。
- **分析工具**: 檢視 `significance_summary.csv` 統計數據、原始 `significance.json` 詳細報告以及 C++ 原始碼 (`RegionProcessor.cpp`, `GlobalTest.cpp`, `SignificanceAnalyzer.cpp`)。

## 3. 案例分析詳細結果

### 3.1 案例 A: FP 被誤判為顯著 (False Positive Claimed Significant)
- **樣本 ID**: `chr4:9649610` (來自 FP 數據集)
- **算法判斷**: `Significant = True`, `GlobalP = 0.001`
- **數據特徵**:
    - **Reads**: 37
    - **Global Alt P-value**: 0.001 (極顯著)
    - **Cramér's V (Effect Size)**: JSON 中為 `1.000` (完美關聯)，但在 CSV 中被輸出為 `0.0000`。
    - **可靠性標記**: `chi_square_reliable = false`
- **原因分析**:
    - 該區域 Reads 數較少 (37)，且可能存在極端分佈（例如某個微小 Cluster 只有 2 條 Reads 且全是 ALT）。
    - 這種情況下，Fisher's Exact Test 會計算出極低的 P-value (0.001)。
    - 然而，由於期望頻數過低，卡方檢驗與 Cramér's V 被標記為「不可靠 (`unreliable`)」。
    - **關鍵問題**: 現行的 `Significant` 判斷邏輯 (`PassedGating && GlobalP <= 0.05`) **完全忽略了可靠性標記**。它僅依據 Fisher P-value，導致這些統計上顯著但數據基礎薄弱（視覺上像噪聲）的案例被保留。

### 3.2 案例 B: 顯著 TP 被誤判為不顯著 (Significant TP Claimed Non-Significant)
- **樣本 ID**: `chr8:84676714` (來自 TP 數據集)
- **算法判斷**: `Significant = False`, `GlobalP = 0.0526`
- **數據特徵**:
    - **Reads**: 63
    - **Global HP P-value**: 0.0526
    - **Cramér's V (Effect Size)**: `0.6975` (強關聯)
    - **Gating**: `PassedGating = True`
- **原因分析**:
    - 該案例的 Cramér's V高達 0.7，顯示單倍型 (HP) 與聚類有極強的關聯，視覺上應能看到明顯分群。
    - 然而，Fisher P-value (0.0526) 略高於硬性閾值 (0.05)。
    - **關鍵問題**: 算法缺乏「彈性」。即使效應量極強，只要 P-value 稍微超標即被剔除。這在樣本數中等（N=63）時容易發生，因樣本數不足以將強效應推至極低 P-value。

## 4. 程式與算法流程檢視
目前的核心判斷邏輯位於 `RegionProcessor.cpp` 的 `write_significance_summary` 函數中：

```cpp
// 當前代碼邏輯
bool is_significant = r.passed_gating && (r.global_p_value <= 0.05);
```

以及 `RegionProcessor.cpp` 在處理不可靠 Cramér's V 時的邏輯：
```cpp
// 若不可靠，則 CSV 輸出 0.0，但不影響 is_significant 判斷
double v_alt = sig_result.global_alt.cramers_v_reliable ? sig_result.global_alt.cramers_v : 0.0;
```

這導致了：
1.  **FP**: 即使 `cramers_v_reliable` 為 `false`，只要 P 小，就判顯著。
2.  **TP**: 即使 `CramersV` 很大，只要 P > 0.05，就判不顯著。

## 5. 結論與建議推測

造成數據與圖片/預期不符的主要原因列表：

| 現象 | 可能原因 | 推測機制 |
| :--- | :--- | :--- |
| **FP 看起來不顯著卻被判顯著** | **數據稀疏導致的統計偽陽性** | 小樣本或極端分佈導致 Fisher P-value 極低，但實際上期望頻數不足，屬於統計學上的「不穩健」結果。視覺上可能表現為僅由 1-2 條 Reads 構成的異常群。 |
| **TP 看起來顯著卻被判不顯著** | **P-value 閾值過於僵化** | 樣本數中等時，強效應 (High V) 可能僅產生邊緣顯著的 P-value (如 0.05-0.10)。硬性切斷導致漏失。 |

### 建議改進方案
1.  **引入可靠性過濾**: 修改 `Significant` 判斷標準，要求 `cramers_v_reliable` 必須為 `true`，或者當 `reliable` 為 `false` 時，要求更嚴格的 P-value (如 0.001) 或要求最小 Cluster Size。
2.  **採納 Heuristic Score**: `HeuristicScore` 已經考慮了 Cramér's V 的加分。可以考慮改用 `HeuristicScore >= Threshold` 作為顯著性標準，或允許 `High CramersV` 救援 `Borderline P-value` 的情況。
3.  **雙重門檻**:
    - `Significant = (P <= 0.05 && Reliable) OR (P <= 0.01 && !Reliable)`
    - 或 `Significant = (P <= 0.05) && (CramersV >= 0.1)` (強制要求最小效應量，過濾掉僅由極大樣本數驅動的微弱效應，或不穩定的極端分佈)。

## 6. 策略優化討論：分群優先 vs 標籤優先

目前的算法邏輯屬於**「分群優先」(Cluster-First)** 策略。
- **流程**: 先計算 Methylation Distance -> 進行 Hierarchical Clustering -> 得到 Cluster Labels -> 檢驗 Cluster 與 HP/Allele Labels 的關聯 (Fisher Test)。
- **優點**: 能夠發現非線性的、未預期的結構。
- **缺點**: 當數據稀疏時，Clustering 容易產生無意義的微小分群 (Small Clusters)，導致 Fisher Test 產生誤導性的低 P-value (FP)。

**思考方向：引入「標籤優先」(Label-First) 策略**
- **流程**: 直接利用已知的 HP 或 Allele Labels 將 Reads 分組 -> 計算組間距離 (Between-Group Distance) 與組內距離 (Within-Group Distance) 的差異 (如 PERMANOVA 或簡單的距離平均值檢定)。
- **預期效果**:
    - **對抗 FP**: 因為不依賴不穩定的 Clustering 結果，若組間距離無顯著差異，即可直接判定無效，避免了「隨機分群剛好對應到標籤」的風險。
    - **救援 TP**: 對於強效應但 P-value 邊緣的案例 (如 effect size 大但 N 小)，Label-First 直接測量距離差異可能比 Fisher Test 的列聯表檢定更敏感且直觀。

| 策略 | 適用場景 | 優勢 | 風險 |
| :--- | :--- | :--- | :--- |
| **Cluster-First (現狀)** | 探索性分析，未知結構，或混合了多種效應的情況 | 能發現潛在的新亞型 | 易受噪聲干擾產生 FP (過擬合) |
| **Label-First (建議)** | 驗證性分析，明確想知道 HP/Allele 是否影響甲基化 | 穩健，直接回答科學問題 | 可能漏掉非標籤驅動的真實結構 |

**建議實作**:
可在現有的 `PassedGating` 邏輯中，加入 Label-First 的快速檢定 (如比較 HP1 vs HP2 的平均距離) 作為輔助條件。

## 7. 最低限制建議 (Minimum Requirements)

為了減少統計偽陽性並提高結果可靠性，建議加入以下最低限制：

### 7.1 樣本數限制
- **Minimum Reads**: 建議設定為 **20** 或 **30** (目前可能有更低的容許值導致不穩定)。
- **Minimum Cluster Size**: 進行 Fisher Test 時，若任一 Cluster 的 Reads 數少於 **5** (或總 Reads 的 10%)，應標記為 `Unreliable` 或直接忽略該測試結果。

### 7.2 效應量限制 (Effect Size Threshold)
- **Cramér's V**:
    - `Unreliable` (卡方檢定條件未滿足) 時，直接視為不顯著，或要求極嚴格的 P-value (例如 1e-5)。
    - 強制要求 `Cramér's V >= 0.2` (或其他經驗值) 才能判為顯著，避免「大樣本下的微弱效應」或「小樣本下的隨機顯著」。

### 7.3 P-value 動態閾值
- 不應使用單一的 0.05。建議根據樣本數 (N) 動態調整，或使用 **Benjamini-Hochberg (FDR)** 校正來控制多重檢定問題。

