# InterSubMod 專案全貌與技術總結

## 1. 專案資訊 (Project Info)

* **專案名稱**: InterSubMod (Inter-Subclonal Methylation Analysis)
* **目前可支持的核心目標**: 以長讀長上的 somatic-SNV 物理分子共現建立**局部、允許 recurrence 的最小突變狀態候選結構**，再以甲基化做 pattern-conditioned association。這不是 confirmed cellular subclone、lineage 或全樣本真實 phylogeny。
* **技術特點**:
  * **高效能 C++ 核心**: 採用 C++17 標準，結合 OpenMP 平行運算。
  * **精確甲基化解析**: 支援 BAM 格式中的 MM/ML 標籤解析，精確定位 CpG 位點的甲基化狀態。
  * **多樣化距離度量**: 支援 L1 Distance、NHD (Normalized Hamming Distance) 等多種距離算法，用於量化 Read 間的表觀遺傳差異。
  * **下游視覺化工具**: Python 工具可由 C++ 輸出的資料生成距離與 read-level methylation cluster 熱圖；C++ binary 本身不會自動渲染 PNG。

---

## 2. 完整分析流程 (Complete Workflow)

`scripts/run_vcf_all_snv.sh` 是**實驗室內部 orchestration**：它能在既定內部資料配置下串接核心運算與 Python 圖表，但依賴未納入 Git 的 `/big7`／`/big8` 絕對路徑。Public clone 目前沒有 runnable fixture，不能把此腳本描述為一鍵 portable workflow。

### 流程步驟

1. **數據讀取與過濾 (Data Ingestion & Filtering)**
    * 輸入：Tumor BAM、Ref Genome、Target SNVs (VCF)。
    * 處理：
        * 根據 VCF 中的 SNV 位點，提取周邊 Reads。
        * 解析 CIGAR 字串，進行座標校正。
        * 過濾低品質 Reads (Mapping Quality, Length)。
        * 識別 Reads 的單倍體型 (Haplotype, HP tag)。

2. **甲基化特徵提取 (Methylation Extraction)**
    * 解析 MM (Base Modification) 與 ML (Probability) 標籤。
    * 將甲基化機率映射至基因組座標。
    * 建構 Read × CpG 甲基化矩陣。

3. **距離計算與矩陣建構 (Distance Calculation)**
    * 計算 Read 對 Read (Read-Read) 的甲基化距離。
    * 支援算法：
        * **L1 Distance**: 曼哈頓距離，適合連續機率值。
        * **NHD (Normalized Hamming Distance)**: 歸一化漢明距離，適合二值化或離散化數據。

4. **視覺化與輸出 (Visualization & Output)**
    * **Distance Heatmap**: 展示 Read 間的相似度結構，包含雙向階層分群 (Hierarchical Clustering)。
    * **Cluster Heatmap**: 展示 Read 與 CpG 位點的甲基化模式與 read-level clusters；這些群不是細胞亞克隆。

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

執行分析後，主要輸出在 `output/YYYYMMDD_vcf_*` 目錄，包含：

* `significance_summary.csv`: 區域層級顯著性摘要（主分析表）。
* `significance_statistics.txt`: 全域統計摘要。
* `full_execution_analysis.log`: 完整執行日誌。
* `filtered_snv_tp/chr*/...`: 依 VCF / 染色體分層的細部輸出與圖表。

---

## 5. 技術亮點與效能 (Highlights & Performance)

* **平行化設計**: 使用 OpenMP 實現 Region-level 平行化；目前沒有鎖定 hardware、commit、
  input distribution、repetitions 與 thread-scaling curve 的公開 benchmark，因此不提供
  speedup 或 per-region latency 數字。
* **精準座標映射**: 解決了 ONT 數據中常見的 Insertion/Deletion 對甲基化座標偏移的影響，確保每個 CpG 位點精確對齊。
* **完整 MM/ML 支援**: 正確處理多種修飾類型共存時的 ML 陣列偏移 (Offset) 問題。

---

## 6. 開發與維護 (Development)

* **建置**: CMake >= 3.14
* **依賴**: HTSlib, OpenMP, Eigen3, Python3 (Matplotlib, Seaborn, Scipy, Pandas)
* **測試**: GoogleTest 框架 (`bin/test_*`) 與 整合測試腳本 (`scripts/verify_output.sh`)。
