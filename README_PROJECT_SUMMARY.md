# InterSubMod 專案全貌與技術總結

## 1. 專案資訊 (Project Info)

* **專案名稱**: InterSubMod (Inter-Subclonal Methylation Analysis)
* **目前可支持的研究目標**: 以長讀長上的 somatic-SNV 物理分子共現建立**局部、允許 recurrence 的最小突變狀態候選結構**，再以甲基化做 pattern-conditioned association。候選結構由獨立的 research exact-PS solver pipeline 產生；`inter_sub_mod` 本體只產生 per-region 甲基化、距離、read clustering 與統計，不輸出該 candidate funnel。兩條線都不能確認 cellular subclone、lineage 或全樣本真實 phylogeny。
* **技術特點**:
  * **C++ 分區平行核心**: 採用 C++17 與 OpenMP；效能需以指定 commit、資料與硬體的 benchmark receipt 評估。
  * **MM/ML 甲基化解析**: 解析 BAM 的 MM/ML 標籤，並依 CIGAR 與 reference 對應 CpG 座標；正確性邊界以測試與具名 fixture 為準。
  * **多樣化距離度量**: 支援 L1 Distance、NHD (Normalized Hamming Distance) 等多種距離算法，用於量化 Read 間的表觀遺傳差異。
  * **下游視覺化工具**: Python 工具可由 C++ 輸出的資料生成距離與 read-level methylation cluster 熱圖；C++ binary 本身不會自動渲染 PNG。

---

## 2. 完整分析流程 (Complete Workflow)

`scripts/run_vcf_all_snv.sh` 是**實驗室內部 orchestration**：它能在既定內部資料配置下串接核心運算與 Python 圖表，但依賴未納入 Git 的 `/big7`／`/big8` 絕對路徑。Public clone 提供 tiny synthetic DEMO fixture 驗證軟體接線，但不是 real-data science；不能因此把內部腳本描述為一鍵 portable workflow。

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

### A. C++ 分析核心 (`src/`)

| 模組 | 檔案 | 功能描述 |
| :--- | :--- | :--- |
| **BamReader** | `src/core/BamReader.cpp` | 基於 HTSlib，負責 BAM 讀取與索引查詢。 |
| **ReadParser** | `src/core/ReadParser.cpp` | 解析 BAM Record，提取 CIGAR、Mapping Quality、HP Tag 等資訊。 |
| **MethylationParser** | `src/core/MethylationParser.cpp` | **關鍵模組**。解析 MM/ML 標籤、處理 Delta Encoding，並依 CIGAR 將修飾機率對應回 Reference 座標。 |
| **MatrixBuilder** | `src/core/MatrixBuilder.cpp` | 動態建構稀疏甲基化矩陣 (Sparse Methylation Matrix)，優化記憶體使用。 |
| **RegionProcessor** | `src/core/RegionProcessor.cpp` | 管理多執行緒任務調度，平行處理多個 SNV 區域。 |
| **RegionWriter** | `src/io/RegionWriter.cpp` | 輸出標準化的 TSV/CSV 檔案，供後續分析使用。 |

### B. Python 視覺化工具 (`tools/`)

| 工具 | 檔案 | 功能描述 |
| :--- | :--- | :--- |
| **Distance Plotter** | `tools/plot_distance_heatmap.py` | 繪製 **Read-Read 距離矩陣熱圖**。使用 `scipy.cluster.hierarchy` 展示 read 間的甲基化相似度與階層分群；這不是細胞關係或演化譜系。 |
| **Cluster Plotter** | `tools/plot_cluster_heatmap.py` | 繪製 **Read-CpG 甲基化模式熱圖**。X 軸為 CpG 位點，Y 軸為分群後的 Reads。直觀顯示甲基化模式的異質性。 |

---

## 4. 輸出檔案結構 (Output Structure)

直接執行 core 時，未指定 `--output-dir` 的 generic root 是 `output/`：

* root summary：`run_params.json`、`run_summary.json`、`significance_summary.csv`、`significance_statistics.txt`、`region_stratification_status.tsv`。
* per-region：`output/<VCF_STEM>/<chr>/<chr_pos>/<chr_start_end>/...`，包含 region 級 JSON/TSV/CSV/Newick 等產物。
* `full_execution_analysis.log`、日期命名目錄與 `filtered_snv_tp` 是 `scripts/run_vcf_all_snv.sh` 及特定 VCF stem 的 internal-wrapper 產物，不是 direct-core generic default。

---

## 5. 技術亮點與效能 (Highlights & Performance)

* **平行化設計**: 使用 OpenMP 實現 Region-level 平行化；目前沒有可重現的 scaling benchmark，因此不宣稱 32 核線性加速或 <300 ms/region。
* **CIGAR-aware 座標對應**: 將 insertion/deletion 納入 read-to-reference 座標轉換；不得在沒有對應測試或 fixture receipt 時宣稱所有 CpG 均完全正確。
* **MM/ML 解析（commit／fixture 範圍）**: 目前實作與測試覆蓋多種修飾類型共存時的 ML 陣列偏移 (Offset)；未經具體 fixture／test receipt 覆蓋的 tag 組合不宣稱已完整驗證。

---

## 6. 開發與維護 (Development)

* **建置**: CMake >= 3.14
* **依賴**: HTSlib, OpenMP, Eigen3, Python3 (Matplotlib, Seaborn, Scipy, Pandas)
* **測試**: GoogleTest 框架 (`bin/test_*`) 與 整合測試腳本 (`scripts/verify_output.sh`)。
