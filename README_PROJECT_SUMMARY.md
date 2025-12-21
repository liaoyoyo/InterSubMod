# InterSubMod 專案全貌與技術總結

## 1. 專案資訊 (Project Info)

*   **專案名稱**: InterSubMod (Inter-Subclonal Methylation Analysis)
*   **核心目標**: 通過整合長讀取 (Long-read) 測序數據中的甲基化 (Methylation) 與體細胞變異 (Somatic SNVs)，解析腫瘤內的亞克隆結構 (Subclonal Structure) 與表觀遺傳異質性 (Epigenetic Heterogeneity)。
*   **技術特點**:
    *   **高效能 C++ 核心**: 採用 C++17 標準，結合 OpenMP 平行運算。
    *   **精確甲基化解析**: 支援 BAM 格式中的 MM/ML 標籤解析，精確定位 CpG 位點的甲基化狀態。
    *   **多樣化距離度量**: 支援 L1 Distance、NHD (Normalized Hamming Distance) 等多種距離算法，用於量化 Read 間的表觀遺傳差異。
    *   **自動化視覺化**: 整合 Python 繪圖工具，自動生成距離熱圖 (Distance Heatmap) 與分群熱圖 (Cluster Heatmap)。

---

## 2. 完整分析流程 (Complete Workflow)

本專案提供一鍵式自動化腳本 (`scripts/run_full_vcf_test.sh`)，串聯從原始 BAM 檔到最終視覺化圖表的完整流程。

### 流程步驟：

1.  **數據讀取與過濾 (Data Ingestion & Filtering)**
    *   輸入：Tumor BAM、Ref Genome、Target SNVs (VCF)。
    *   處理：
        *   根據 VCF 中的 SNV 位點，提取周邊 Reads。
        *   解析 CIGAR 字串，進行座標校正。
        *   過濾低品質 Reads (Mapping Quality, Length)。
        *   識別 Reads 的單倍體型 (Haplotype, HP tag)。

2.  **甲基化特徵提取 (Methylation Extraction)**
    *   解析 MM (Base Modification) 與 ML (Probability) 標籤。
    *   將甲基化機率映射至基因組座標。
    *   建構 Read × CpG 甲基化矩陣。

3.  **距離計算與矩陣建構 (Distance Calculation)**
    *   計算 Read 對 Read (Read-Read) 的甲基化距離。
    *   支援算法：
        *   **L1 Distance**: 曼哈頓距離，適合連續機率值。
        *   **NHD (Normalized Hamming Distance)**: 歸一化漢明距離，適合二值化或離散化數據。

4.  **視覺化與輸出 (Visualization & Output)**
    *   **Distance Heatmap**: 展示 Read 間的相似度結構，包含雙向階層分群 (Hierarchical Clustering)。
    *   **Cluster Heatmap**: 展示 Read 與 CpG 位點的甲基化模式，直觀呈現亞克隆分群。

---

## 3. 核心模組架構 (Core Modules)

### A. C++ 高效運算核心 (`src/`)

| 模組 | 檔案 | 功能描述 |
| :--- | :--- | :--- |
| **BamReader** | `src/core/BamReader.cpp` | 基於 HTSlib，負責高效、線程安全 (Thread-safe) 的 BAM 讀取與索引查詢。 |
| **ReadParser** | `src/core/ReadParser.cpp` | 解析 BAM Record，提取 CIGAR、Mapping Quality、HP Tag 等資訊。 |
| **MethylationParser** | `src/core/MethylationParser.cpp` | **關鍵模組**。解析複雜的 MM/ML 標籤，處理 Delta Encoding，將修飾機率精確映射回 Reference 座標。 |
| **MatrixBuilder** | `src/core/MatrixBuilder.cpp` | 動態建構稀疏甲基化矩陣 (Sparse Methylation Matrix)，優化記憶體使用。 |
| **RegionProcessor** | `src/core/RegionProcessor.cpp` | 管理多執行緒任務調度，平行處理多個 SNV 區域。 |
| **RegionWriter** | `src/io/RegionWriter.cpp` | 輸出標準化的 TSV/CSV 檔案，供後續分析使用。 |

### B. Python 視覺化工具 (`tools/`)

| 工具 | 檔案 | 功能描述 |
| :--- | :--- | :--- |
| **Distance Plotter** | `tools/plot_distance_heatmap.py` | 繪製 **Read-Read 距離矩陣熱圖**。使用 `scipy.cluster.hierarchy` 進行分群，展示 Read 間的親緣關係。 |
| **Cluster Plotter** | `tools/plot_cluster_heatmap.py` | 繪製 **Read-CpG 甲基化模式熱圖**。X 軸為 CpG 位點，Y 軸為分群後的 Reads。直觀顯示甲基化模式的異質性。 |

---

## 4. 輸出檔案結構 (Output Structure)

執行分析後，每個 Region 目錄下包含以下檔案：

*   `metadata.txt`: 區域基本資訊、統計數據 (Read 數、CpG 數、處理時間)。
*   `reads.tsv`: Read 詳細資訊 (ID, Name, Mapping Quality, HP Tag, etc.)。
*   `methylation.csv`: 原始甲基化矩陣 (Rows: Reads, Cols: CpG Sites, Values: 0.0-1.0)。
*   `distance_matrix_[METRIC].csv`: Read-Read 距離矩陣。
*   `distance_heatmap_[METRIC].png`: 距離矩陣熱圖。
*   `cluster_heatmap_[METRIC].png`: 甲基化分群熱圖。

---

## 5. 技術亮點與效能 (Highlights & Performance)

*   **平行化設計**: 使用 OpenMP 實現 Region-level 平行化，在 32 核環境下可達線性加速比，單 Region 平均處理時間 < 300ms。
*   **精準座標映射**: 解決了 ONT 數據中常見的 Insertion/Deletion 對甲基化座標偏移的影響，確保每個 CpG 位點精確對齊。
*   **完整 MM/ML 支援**: 正確處理多種修飾類型共存時的 ML 陣列偏移 (Offset) 問題。

---

## 6. 開發與維護 (Development)

*   **建置**: CMake >= 3.14
*   **依賴**: HTSlib, OpenMP, Eigen3, Python3 (Matplotlib, Seaborn, Scipy, Pandas)
*   **測試**: GoogleTest 框架 (`bin/test_*`) 與 整合測試腳本 (`scripts/verify_output.sh`)。
