# CLAUDE.md - Claude Code 專案指南

## 專案概述

**InterSubMod (Inter-Subclonal Methylation Analysis)** 是一個生物資訊工具，專門用於透過長讀取 (Long-read) 測序數據偵測腫瘤樣本中的亞克隆結構 (Subclonal Structure)。本專案整合甲基化模式 (Methylation Patterns)、體細胞變異 (Somatic SNVs) 以及單倍體型 (Haplotypes) 來分析表觀遺傳異質性 (Epigenetic Heterogeneity)。

### 核心技術特點

- **高效能 C++17 核心**: 結合 OpenMP 平行運算，單 Region 平均處理時間 < 300ms
- **精確甲基化解析**: 支援 BAM 格式 MM/ML 標籤，精確定位 CpG 位點甲基化狀態
- **多元距離度量**: L1 / L2 / NHD / Bernoulli / Jaccard 等多種距離算法
- **統計顯著性分析**: PERMANOVA、卡方檢驗、Cramér's V 效應量計算
- **自動化視覺化**: Python 工具生成距離熱圖 (Distance Heatmap) 與分群熱圖 (Cluster Heatmap)

---

## 建置命令

```bash
# 配置並建置 (從專案根目錄)
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Debug 模式建置
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# 執行所有測試
cd build && ctest --output-on-failure

# 執行特定測試
./build/tests/test_<name>
```

---

## 專案結構

### C++ 核心模組 (`src/core/`)

| 檔案 | 功能 |
| :--- | :--- |
| `BamReader.cpp` | 基於 HTSlib 的高效 BAM 讀取與索引查詢 |
| `ReadParser.cpp` | 解析 BAM Record (CIGAR、Mapping Quality、HP Tag) |
| `MethylationParser.cpp` | **關鍵模組**：解析 MM/ML 標籤，處理 Delta Encoding |
| `MatrixBuilder.cpp` | 建構稀疏甲基化矩陣 |
| `DistanceMatrix.cpp` | 多種距離度量計算 (NHD/L1/L2/Bernoulli/Jaccard) |
| `HierarchicalClustering.cpp` | 階層分群 (UPGMA/Ward/Single/Complete) |
| `RegionProcessor.cpp` | 多執行緒任務調度，平行處理 SNV 區域 |
| `SignificanceAnalyzer.cpp` | 顯著性分析主控制器 |
| `GlobalTest.cpp` | 卡方關聯測試與 Cramér's V |
| `LabelTest.cpp` | HP/Allele 標籤 PERMANOVA 測試 |
| `StructureTest.cpp` | PERMANOVA 與離散度檢驗 |

### I/O 模組 (`src/io/`)

| 檔案 | 功能 |
| :--- | :--- |
| `RegionWriter.cpp` | 輸出標準化 TSV/CSV 檔案 |
| `TreeWriter.cpp` | 輸出階層樹結構 |

### 工具程式 (`src/utils/`)

| 檔案 | 功能 |
| :--- | :--- |
| `Logger.cpp` | 日誌系統 |
| `FastaReader.cpp` | 參考基因組讀取 |

### Python 視覺化工具 (`tools/`)

| 檔案 | 功能 |
| :--- | :--- |
| `plot_cluster_heatmap.py` | Read-CpG 甲基化模式熱圖 |
| `plot_distance_heatmap.py` | Read-Read 距離矩陣熱圖 |
| `compare_vcf_results.py` | VCF 結果比較分析 |
| `f1_optimization_analysis.py` | F1 分數優化分析 |

### 自動化腳本 (`scripts/`)

| 檔案 | 功能 |
| :--- | :--- |
| `run_vcf_all_snv.sh` | 單一 VCF 分析執行腳本 |
| `run_batch_vcf_analysis.sh` | 批次 VCF 分析工作流程 |
| `verify_output.sh` | 輸出驗證腳本 |

---

## 程式碼規範

- **C++17** 標準
- 遵循 `.clang-format` (Google 風格、120 字元行寬、4 空格縮排)
- 提交前格式化：`clang-format -i <file>`
- **預設回應與註解語言**：繁體中文
- **程式碼註解語言**：英文

---

## 關鍵依賴

| 依賴 | 用途 |
| :--- | :--- |
| HTSlib | BAM 檔案處理 |
| OpenMP | 平行運算 |
| Eigen3 | 線性代數 |
| GoogleTest | 單元測試 |
| jemalloc 5.3.0 | 記憶體分配 |
| Python3 + Matplotlib/Seaborn/Scipy/Pandas | 視覺化 |

---

## 測試

使用 GoogleTest 框架，測試檔案位於 `tests/` 目錄。

```bash
# 執行所有測試
cd build && ctest

# 執行特定測試
./build/tests/test_global_test
./build/tests/test_distance_matrix
./build/tests/test_hierarchical_clustering
```

---

## 常用工作流程

### 1. 完整 VCF 分析 (預設執行命令)

```bash
./scripts/run_vcf_all_snv.sh --mode all-with-w1000
```

### 2. 批次分析 (TP/FP 比較)

```bash
./scripts/run_batch_vcf_analysis.sh
```

### 3. 單點快速驗證

```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification
```

### 4. 新增分析模組

1. 在 `include/core/` 建立標頭檔
2. 在 `src/core/` 建立實作檔
3. 更新 `CMakeLists.txt` (如需要)
4. 在 `tests/` 撰寫單元測試

---

## 輸出檔案結構

每個 Region 目錄包含：

- `metadata.txt`: 區域基本資訊與統計數據
- `reads.tsv`: Read 詳細資訊
- `methylation.csv`: 甲基化矩陣 (Read × CpG)
- `distance_matrix_[METRIC].csv`: Read-Read 距離矩陣
- `significance_summary.csv`: 顯著性分析結果彙總
- `*.png`: 視覺化圖表

---

## 開發重點

1. **甲基化解析精確度**: MM/ML 標籤解碼與 CIGAR 座標校正的正確性
2. **距離計算效能**: OpenMP 平行化與稀疏矩陣最佳化
3. **統計假設檢驗**: PERMANOVA p-value 與 Cramér's V 效應量計算
4. **視覺化品質**: 熱圖標註 (HP tag / Allele) 與分群樹狀圖
5. **批次分析流程**: TP/FP 位點比較與 F1 分數優化

---

## 文件資源

- `README_PROJECT_SUMMARY.md`: 專案完整技術總結
- `QUICKSTART.md`: 快速入門指南
- `docs/`: 完整技術文件 (API、架構、開發、報告)
