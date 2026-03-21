# 1202_InterSubMod_聚類與演化分析開發計劃

## 1. 研究目標

本研究旨在利用 Read 與 甲基化 (Methylation) 的數值矩陣，計算 Read 間的甲基化距離矩陣，並進行層次式聚類 (Hierarchical Clustering) 分析。隨後，結合 Bootstrap 方法進行驗證，並繪製 Cluster Heatmap。

進一步，我們將結合 Read 的各類標籤 (Tag) 資訊（如 Tumor/Normal、Somatic 位點 ALT/REF、HP Tag 分型、正反股資訊等），分析這些標籤組合在聚類結果中的分佈情形，並透過統計方法驗證其與距離矩陣及演化樹結構的顯著相關性。

最終目標是輸出單個位點的詳細分析結果，以及多個位點的綜合整體分析報告。

---

## 2. 開發階段與詳細規劃

### 階段一：Read-Read 甲基化距離矩陣計算 (Distance Matrix Calculation)

**目標**：從 Read-Methylation 矩陣計算出 Read 兩兩之間的距離。

**輸入**：

* Read-Methylation Matrix (行：Reads, 列：CpG Sites, 值：甲基化狀態 0/1/缺失)

**方法與演算法**：

1. **距離度量 (Distance Metrics)**：
    * **Euclidean Distance**：適用於數值型資料，但在稀疏矩陣中需注意缺失值處理。
    * **Hamming Distance**：適用於二元資料 (0/1)，計算不一致的位點比例。
    * **Pearson Correlation Distance**：$1 - r$，關注變化趨勢而非絕對值。
    * **Jaccard Distance**：關注非零元素的交集與聯集（可能較少用，視情況而定）。
2. **缺失值處理 (Missing Data Handling)**：
    * **Pairwise Deletion**：計算兩 Read 距離時，僅考慮兩者皆有數據的 CpG 位點。需設定最小重疊位點閾值 (Minimum Overlap Threshold)，若重疊位點過少則視為距離無限大或標記為無效。

**程式實作 (C++)**：

* 擴充 `DistanceMatrix` 類別，支援多種距離計算策略。
* 優化計算效能：利用位元運算 (Bitwise operations) 加速二元資料的比較。
* 輸出格式：標準距離矩陣 (CSV/TSV)，包含 Read ID 索引。

### 階段二：層次式聚類與 Bootstrap 驗證 (Hierarchical Clustering & Bootstrap)

**目標**：根據距離矩陣建構演化樹 (Dendrogram) 並評估其穩定性。

**方法與演算法**：

1. **聚類演算法**：
    * **UPGMA (Unweighted Pair Group Method with Arithmetic Mean)**：假設演化速率恆定，適合距離矩陣。
    * **Neighbor Joining (NJ)**：不假設速率恆定，常用於演化樹構建。
    * **Ward's Method**：最小化群內變異數，適合產生緊湊的群集。
    * **Single Linkage / Complete Linkage**：視需求選擇。
2. **Bootstrap 驗證**：
    * **重抽樣 (Resampling)**：對 CpG 位點 (Columns) 進行有放回的重抽樣，產生 $N$ 組新的 Read-Methylation 矩陣。
    * **重複聚類**：對每一組重抽樣資料計算距離矩陣並建構樹狀圖。
    * **一致性計算**：計算原始樹狀圖中每個節點 (Clade) 在 Bootstrap 樹中出現的頻率 (Support Value)。

**程式實作**：

* **C++ 方案**：使用現有庫 (如類似 FastTree 或自行實作基礎 UPGMA/NJ) 進行快速計算。
* **Python 整合方案**：若 C++ 實作複雜度過高，可呼叫 Python `scipy.cluster.hierarchy` 或 `scikit-learn` 進行聚類計算，C++ 負責資料準備與流程控制。

### 階段三：視覺化輸出 (Visualization)

**目標**：產生精緻的 Cluster Heatmap 與演化樹圖片。

**工具選擇**：

* **Python (推薦)**：使用 `seaborn.clustermap` 或 `matplotlib` + `scipy.dendrogram`。Python 在繪圖美觀度與彈性上遠優於 C++ 原生繪圖庫。
* **C++**：若需直接輸出，可考慮 `matplotlib-cpp` (需 Python 環境) 或輸出 SVG 格式，但開發成本較高且美觀度調整不易。

**實作細節**：

* **Heatmap**：顯示 Read (Y軸) x CpG Site (X軸) 的甲基化狀態。
* **Dendrogram**：顯示在 Heatmap 側邊，呈現聚類結構。
* **Annotation Bars**：在 Heatmap 旁添加顏色條，標示 Read 的 Tag 資訊 (Tumor/Normal, HP, Strand 等)。

### 階段四：標籤關聯與統計分析 (Tag Association & Statistical Analysis)

**目標**：驗證 Read 的生物標籤 (Tags) 是否與甲基化聚類結構有顯著關聯。

**分析維度**：

1. **樣本來源**：Tumor vs Normal
2. **Somatic 變異**：ALT (變異) vs REF (參考)
3. **Haplotype (HP)**：HP1, HP2, HP1-1, HP2-1, HP3, Unphased
4. **股向 (Strand)**：Forward (+) vs Reverse (-)

**統計驗證方法**：

1. **PERMANOVA (Permutational Multivariate Analysis of Variance)**：
    * 檢驗不同 Tag 群組 (如 Tumor vs Normal) 在距離矩陣上的中心點是否顯著不同。
    * 使用 `scikit-bio` 或 R 的 `vegan` 套件 (透過 Python 呼叫)。
2. **Mantel Test**：
    * 檢驗「甲基化距離矩陣」與「標籤距離矩陣」之間的相關性。
3. **演化樹拓撲分析**：
    * 計算特定 Tag 在演化樹特定分支 (Clade) 中的富集程度 (Enrichment Analysis, e.g., Fisher's Exact Test)。
    * 計算 **Phylogenetic Signal** (如 Lambda 值)，評估 Tag 分佈是否符合演化樹結構。

### 階段五：綜合分析與報告輸出 (Output & Reporting)

**輸出內容**：

1. **單一變異位點分析 (Single Locus Analysis)**：
    * 針對每個 Somatic SNV 位點視窗。
    * 輸出：原始矩陣、距離矩陣、聚類樹檔 (Newick format)、Heatmap 圖片、統計檢定 P-value 表。
2. **多位點綜合分析 (Multi-Locus Aggregate Analysis)**：
    * 彙整所有位點的統計結果。
    * 分析整體趨勢：例如，是否大部分位點的 Tumor Reads 都聚類在一起？HP1/HP2 是否在甲基化上有顯著分層？

## 3. 程式撰寫與開發注意事項

1. **模組化設計**：
    * 將距離計算、聚類、統計分析拆分為獨立模組或類別，便於替換演算法。
2. **效能考量**：
    * Read 數量可能很大，距離矩陣計算為 $O(N^2)$，需注意記憶體與計算時間。
    * 可考慮使用多執行緒 (OpenMP) 加速距離計算。
3. **Python 與 C++ 的協作**：
    * 建議核心計算 (矩陣生成、距離計算) 由 C++ 完成。
    * 複雜統計與繪圖由 Python 腳本完成 (C++ 透過 `system()` 呼叫或輸出中間檔後由 Python 批次處理)。
4. **資料格式標準化**：
    * 確立 C++ 輸出給 Python 的資料介面 (如 CSV 格式：Header 包含位點資訊，Rows 包含 Read ID, Tags, Methylation String)。

## 4. 待辦事項清單 (To-Do List)

* [ ] **C++**: 實作多種距離演算法 (Euclidean, Hamming, Pearson)。
* [ ] **C++**: 實作 Bootstrap 重抽樣機制。
* [ ] **Python**: 開發繪圖腳本 (接收矩陣與 Tags，繪製帶有 Annotation 的 Clustermap)。
* [ ] **Python**: 開發統計分析模組 (PERMANOVA, Mantel Test)。
* [ ] **Integration**: 整合 Shell/Python 流程，自動化執行「讀取 -> 計算 -> 繪圖 -> 統計 -> 報告」。
