# 雙向驗證流程分析與開發策略報告
**日期**: 2025-12-30

## 1. 摘要
為了更準確地判斷甲基化距離矩陣中的結構意義，我們將從單一的「分群優先」(Cluster-First) 策略，轉向 **「分群↔標籤」雙向互證 (Triangulation)** 策略。本報告分析了現有代碼的差距，並提出了整合 Label-First 與 Cluster-First 的完整邏輯流程，以區分強關聯、Subclone、與統計噪聲。

## 2. 現狀代碼分析

### 2.1 現有功能 (Existing Capabilities)
- **StructureTest.cpp**: 已實作 `run_permanova` (PERMANOVA/Adonis) 與 `check_dispersion` (Dispersion Homogeneity)。這為 Label-First 驗證提供了核心算力。
- **SignificanceAnalyzer.cpp**: 目前主要依賴 Fisher's Exact Test (`GlobalTest`)。雖然有提及 PERMANOVA，但在 `analyze_simple` 流程中被 **強制關閉** (`enable_permanova = false`)。更關鍵的是，現有的 PERMANOVA 呼叫邏輯是基於 `cluster_labels` (分群結果)，而非 `HP/Allele` 標籤，這不符合 Label-First 的定義。

### 2.2 缺口 (Gaps)
1.  **Label-First logic 缺失**: 尚未實作「直接使用 HP/Allele 標籤將 Reads 分組並計算距離差異 (Delta)」的功能。
2.  **PERMANOVA 應用錯誤**: 目前僅用於驗證 Cluster 的結構強度，未用於驗證 Label 對距離的解釋力。
3.  **穩定度檢查 (Stability Check) 缺失**: 尚未實作 Bootstrap/Subsampling 來驗證 Clustering 的穩定性，容易受小樣本噪聲影響。
4.  **雙向整合邏輯缺失**: 缺乏將 Label-First 結果與 Cluster-First 結果綜合判斷的決策樹。

## 3. 提出的雙向驗證流程 (Methodology)

### 3.1 核心概念：三角互證
同一個 Region 需同時經過兩條路徑檢驗：
1.  **Path A: Label-First (標籤驗證結構)**
    - 科學問題：已知的單倍型 (HP) 或等位基因 (Allele) 是否導致了甲基化模式的顯著差異？
    - 方法：PERMANOVA (on Labels), Delta Distance (Between - Within), kNN Prediction (Optional).
2.  **Path B: Cluster-First (結構探索標籤)**
    - 科學問題：數據中是否存在自然的次結構 (Substructure)？這些結構是否對應到已知標籤？
    - 方法：Hierarchical Clustering, Stability Check (Bootstrap), Fisher's Test (on Clusters).

### 3.2 判斷邏輯 (Decision Tree)

| 類別 (Class) | Label-First (Path A) | Cluster-First (Path B) | 解釋 |
| :--- | :--- | :--- | :--- |
| **強關聯 (Strong Association)** | **顯著** (High R², Sig P) | **一致且穩定** (High ARI, Stable) | 標籤是驅動甲基化變異的主因。最可信的結果。 |
| **潛在 Subclone (Novel Structure)** | 不顯著 (Low R²) | **顯著且穩定** (Stable Cluster, but Low ARI with Label) | 存在真實的結構，但不完全由 HP/Allele 解釋。可能是細胞亞型或未知的遺傳因子。 |
| **弱關聯/過渡態 (Weak/Transitional)** | 邊緣顯著 (Low R², but Sig P) | 不穩定或部分一致 | 效應存在但不明顯，或是漸層變化。 |
| **統計噪聲/假陽性 (Artifact/Noise)** | 不顯著 | **不穩定** (Unstable Small Clusters) | 隨機分群導致的 Fisher 顯著，應被過濾。 |

## 4. 詳細算法流程設計

### 階段 0: 前置過濾 (Quality Gate)
- Reads >= 20 (or 30)
- Label Group Size >= 5 (or 10%)
- Coverage Check

### 階段 1: Label-First Verification (針對 HP 與 Allele 分別執行)
1.  **分組**: 根據 HP tags (HP1, HP2) 將 Reads 分為兩組 (忽略 Unphased)。
2.  **Delta Calculation**:
    - 計算組內平均距離 (Within-Group Mean Dist)
    - 計算組間平均距離 (Between-Group Mean Dist)
    - $\Delta = \bar{d}_{between} - \bar{d}_{within}$
3.  **Statistic Test**:
    - **PERMANOVA**: 計算 Pseudo-F 與 P-value (Permutation 999次)。
    - **Delta Permutation**: (Optional, 輕量化) 用 Delta 作為統計量進行置換檢定。

### 階段 2: Cluster-First Exploration
1.  **Clustering**: 執行 Hierarchical Clustering (UPGMA/Ward)。
2.  **Stability Check**:
    - Subsampling (e.g., 80% reads) 重複 20-50 次。
    - 計算原始 Cluster 與 Subsample Cluster 的一致性 (e.g., Jaccard Index of co-clustering)。
    - 輸出 `Stability Score` (0-1)。
3.  **Association**: 若 StabilityPass，執行 Fisher's Exact Test 與 Cramér's V。

### 階段 3: 綜合判斷 (Synthesis)
- 結合上述指標，輸出最終分類標籤 (`Strong`, `Subclone`, `Weak`, `Noise`) 與信心分數。

## 5. 輸出資訊擴充
- **CSV/JSON**: 新增 `Label_R2`, `Label_P`, `Delta`, `Cluster_Stability`, `Verification_Class` 等欄位。
- **Visualization**: 
    - Heatmap 標題需標註分類結果。
    - 若 Label-First 顯著，應在 Heatmap 上強制標示 Label 分組線。

## 6. 開發計畫參照
請參閱 `/big8_disk/liaoyoyo2001/InterSubMod/docs/development/plans/2025_12_30/implementation_plan.md`。
