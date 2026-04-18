<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 20251223_Significance_Analysis_Plan

基於先前架構設計與目前程式碼狀態，本文件詳細規劃 **InterSubMod** 專案中「顯著性分析 (Significance Analysis)」模組的實作細節。此模組旨在統計上量化甲基化分群與生物標籤（如 Tumor/Normal、Haplotype、ALT/REF）之間的關聯強度，並驗證 Somatic 變異的真實性。

## 1. 研究目標與背景 (Goal & Context)

### 1.1 問題定義
目前的系統能夠計算 Read 間距離並進行分群，但缺乏一個客觀的統計指標來回答以下問題：
*   **分群是否隨機？** (Is the clustering structure real?)
*   **分群是否與已知生物標籤相關？** (Do clusters correspond to haplotypes or tumor subclones?)
*   **該位點是否為真實的 Somatic 變異？** (Is this SNV supported by a distinct methylation pattern?)

### 1.2 目標功能
實作一套統計檢定流程，針對每個 SNV Region 輸出：
1.  **標籤關聯性 (Label Association)**: 計算分群結果與 HP/Tumor 標籤的 P-value (Fisher's Exact Test / Chi-square)。
2.  **分群品質 (Cluster Quality)**: 計算 Silhouette Score，評估分群的緊密度與分離度。
3.  **綜合評分 (Confidence Score)**: 結合上述指標，為 SNV 給出一個「甲基化支持度」分數。

---

## 2. 統計方法選擇 (Statistical Methods)

### 2.1 關聯性檢定 (Association Tests)

針對每個 Cluster $k$ 與每個 Label $L$ (例如 HP1 vs HP2)，我們建立列聯表 (Contingency Table) 並進行檢定。

#### 2.1.1 Fisher's Exact Test (2x2)
*   **適用場景**: 樣本數較少，或比較兩個特定群體 (e.g., Cluster A vs Rest; Label 1 vs Label 2)。
*   **表格結構**:
    | | Cluster A | Other Clusters |
    |---|---|---|
    | **Label 1 (e.g., HP1)** | $n_{11}$ | $n_{12}$ |
    | **Label 2 (e.g., HP2)** | $n_{21}$ | $n_{22}$ |
*   **輸出**: P-value (雙尾)，衡量 Label 分布是否在 Cluster A 中有顯著差異（富集或排斥）。

#### 2.1.2 Chi-square Test of Independence (RxC)
*   **適用場景**: 總體評估分群結果 (K個 Clusters) 與標籤 (M個 Labels) 是否獨立。
*   **優點**: 提供全域的關聯性評估。
*   **限制**: 當某些格子計數過小 (<5) 時準確度下降，需退回 Fisher 或使用模擬法。本計畫初期以 Fisher 為主，Chi-square 作為全域指標。

#### 2.1.3 多重檢定校正 (Multiple Testing Correction)
*   由於一個 Region 可能有多個 Cluster，且我們同時分析多個 SNV Region，必須控制 False Discovery Rate (FDR)。
*   **方法**:Benjamini-Hochberg (BH) procedure。

### 2.2 分群品質指標 (Clustering Metrics)

#### 2.2.1 Silhouette Score
*   **定義**: 對於每個 Read $i$，計算 $s(i) = \frac{b(i) - a(i)}{\max(a(i), b(i))}$。
    *   $a(i)$: 同群內平均距離。
    *   $b(i)$: 最近鄰群的平均距離。
*   **用途**: 衡量分群的「可信度」。若 P-value 顯著但 Silhouette 很低，可能代表分群邊界模糊，結果需存疑。

---

## 3. 系統架構設計 (System Architecture)

### 3.1 新增資料結構
在 `include/core/DataStructs.hpp` 與 `include/core/Stats.hpp` (新增) 中定義：

```cpp
struct ContingencyTable {
    int n11, n12, n21, n22;
    // ... helper methods ...
};

struct ClusterStats {
    int cluster_id;
    int size;
    
    // Counts
    int count_hp1, count_hp2;
    int count_tumor, count_normal;
    int count_alt, count_ref;
    
    // Statistics
    double p_value_hp;       // Fisher test for HP1 vs HP2
    double p_value_tumor;    // Fisher test for Tumor vs Normal
    double p_value_somatic;  // Fisher test for ALT vs REF
    
    // ...
};

struct AssociationResult {
    int region_id;
    std::vector<ClusterStats> cluster_stats;
    double avg_silhouette_score;
    // ...
};
```

### 3.2 新增統計模組
建立 `src/core/Statistics.cpp` 與 `src/core/SignificanceAnalyzer.cpp`。

*   **Logic**:
    *   `calculate_fisher_exact(n11, n12, n21, n22)`: 實作或引用開源實作 (如 htslib 或 boost，若無依賴則手刻 log-factorial 版本以避免溢位)。
    *   `compute_silhouette(dist_matrix, labels)`: 已在 `HierarchicalClustering` 中有雛形，需提取為通用函式。

### 3.3 整合至 RegionProcessor
在 `RegionProcessor::process_single_region` 流程中加入 S7 步驟：

1.  **Clustering**: 取得 `labels`。
2.  **Analysis**:
    *   呼叫 `SignificanceAnalyzer::analyze(labels, reads_info)`。
    *   計算每個 Cluster 的 HP/Tumor/ALT 計數。
    *   執行 Fisher Test。
3.  **Output**:
    *   將 P-values 寫入 CSV 報告 (`cluster_stats.csv`)。
    *   在 JSONL 中標記主要 Cluster 及其 P-value。

---

## 4. 實作步驟 (Implementation Steps)

### Phase 1: 基礎統計函式庫 (Basic Stats)
- [ ] 實作 `MathUtils::fisher_exact_test`。這需要處理大數階乘的對數運算。
- [ ] 實作 `MathUtils::benjamini_hochberg_correction`。

### Phase 2: 關聯性分析器 (Analyzer)
- [ ] 建立 `SignificanceAnalyzer` class。
- [ ] 實作計數邏輯：給定 reads 分群與 labels，建立 Contingency Table。
- [ ] 整合 Silhouette Score 計算。

### Phase 3: 系統整合與輸出 (Integration)
- [ ] 修改 `RegionProcessor`，在分群後呼叫 Analyzer。
- [ ] 更新 `RegionWriter`，輸出 `cluster_stats.csv`，包含下列欄位：
    `region_id, cluster_id, size, hp1, hp2, p_val_hp, tumor, normal, p_val_tumor, avg_silhouette`
- [ ] CLI 參數：`--compute-significance` (預設開啟)。

---

## 5. 驗證計畫 (Verification Plan)

### 5.1 單元測試
*   **Fisher Test**: 使用標準案例 (如 Wikipedia 例題) 驗證 P-value 計算準確度。
*   **計數邏輯**: 構造假 reads 資料，確認 contingency table 填寫正確。

### 5.2 真實數據驗證 (Bio-validation)
*   **TP/FP Dataset**: 使用 SEQC2 或已知 Somatic 真值集。
*   **預期結果**:
    *   **TP 位點**: 應觀察到顯著的 P-value (e.g., < 0.05) 且 Silhouette 高，代表 ALT reads 形成獨特且一致的甲基化聚類。
    *   **FP 位點**: 應呈現隨機分布 (P-value 不顯著) 或由 Artifact 驅動 (需檢查 MapQ/Strand)。
*   **比較分析**: 繪製 TP 與 FP 的 P-value 分布圖 (QQ plot 或 Density plot)。

## 6. 結論 (Conclusion)
本計畫將補足 InterSubMod 從「探索性分群」到「統計驗證」的關鍵缺口，提供具備生物統計意義的 Somatic 變異與次克隆判讀依據。
