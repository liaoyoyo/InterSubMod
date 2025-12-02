# Read-Read 距離矩陣實作報告

**日期**：2025-12-01  
**版本**：v1.0  
**作者**：AI Assistant

---

## 1. 摘要

本次開發完成了 Read-Read 距離矩陣的計算功能，包括：

1. **多種距離度量支援**：NHD（Normalized Hamming Distance）、L1、L2、CORR（Pearson相關性）、JACCARD
2. **正反股分離計算**：支援 strand-aware 的距離矩陣輸出
3. **完整的參數化配置**：所有可能影響計算的細節都可通過參數配置
4. **高效穩定的實作**：使用 OpenMP 平行化，支援大規模資料處理
5. **完整的單元測試**：15 個測試案例全部通過
6. **整合測試驗證**：30,484 個區域成功輸出距離矩陣

---

## 2. 技術實現

### 2.1 核心類別結構

```
include/core/DistanceMatrix.hpp
├── struct DistanceConfig       // 距離計算配置
├── class DistanceMatrix        // 距離矩陣儲存與計算
└── class DistanceCalculator    // Strand-aware 計算包裝器
```

### 2.2 距離度量實作

#### A. Normalized Hamming Distance (NHD) - 預設

適用於二值化甲基化資料（`binary_matrix`）。

$$D_{NHD}(R_i, R_j) = \frac{\sum_{k \in Common} |M_{ik} - M_{jk}|}{|Common|}$$

- **範圍**：[0, 1]
- **特性**：對二值差異敏感，計算效率高

#### B. L1 (Manhattan) Distance

$$D_{L1}(R_i, R_j) = \frac{\sum_{k \in Common} |P_{ik} - P_{jk}|}{|Common|}$$

- **範圍**：[0, 1]（使用原始機率）
- **特性**：對絕對差異敏感

#### C. L2 (Euclidean) Distance

$$D_{L2}(R_i, R_j) = \sqrt{\frac{\sum_{k \in Common} (P_{ik} - P_{jk})^2}{|Common|}}$$

- **範圍**：[0, 1]
- **特性**：對大差異敏感，平滑

#### D. Pearson Correlation Distance (CORR)

$$D_{CORR}(R_i, R_j) = \frac{1 - r(R_i, R_j)}{2}$$

- **範圍**：[0, 1]（正規化後）
- **特性**：測量模式相似性，對偏移不敏感
- **要求**：至少 3 個共同位點

#### E. Jaccard Distance

$$D_{JACCARD}(R_i, R_j) = 1 - \frac{|A \cap B|}{|A \cup B|}$$

- **範圍**：[0, 1]
- **特性**：僅考慮甲基化位點（binary = 1）

### 2.3 缺失值處理

| 參數 | 預設值 | 說明 |
|------|--------|------|
| `min_common_coverage` (C_min) | 3 | 最小共同 CpG 位點數 |
| `nan_strategy` | MAX_DIST | 不足時的處理策略 |
| `max_distance_value` | 1.0 | MAX_DIST 策略使用的值 |

當兩個 reads 的共同有效 CpG 位點數 < C_min 時：
- **MAX_DIST**：設為 `max_distance_value`（預設 1.0）
- **SKIP**：設為 NaN（需聚類演算法支援）

### 2.4 Strand-Aware 計算

將 reads 依 strand 分類後分別計算：

```cpp
std::pair<DistanceMatrix, DistanceMatrix> compute_strand_specific(
    const MethylationMatrix& methyl_mat,
    const std::vector<ReadInfo>& reads
);
```

輸出檔案：
- `distance_matrix.csv` - 所有 reads
- `distance_forward.csv` - Forward strand (+)
- `distance_reverse.csv` - Reverse strand (-)

---

## 3. 參數配置說明

### 3.1 Config 結構新增欄位

```cpp
struct Config {
    // Distance Matrix Configuration
    bool compute_distance_matrix = true;        // 是否計算距離矩陣
    bool output_distance_matrix = true;         // 是否輸出 CSV
    bool output_strand_distance_matrices = true;// 是否輸出 strand-specific
    
    DistanceMetricType distance_metric = NHD;   // 距離度量類型
    int min_common_coverage = 3;                // C_min
    NanDistanceStrategy nan_distance_strategy = MAX_DIST;
    double max_distance_value = 1.0;            // MAX_DIST 的值
    
    bool distance_use_binary = true;            // 使用 binary matrix
    bool distance_pearson_center = true;        // CORR: mean-centered
    bool distance_jaccard_include_unmeth = false; // JACCARD: 是否包含 unmeth
};
```

### 3.2 命令列參數

| 參數 | 類型 | 預設值 | 說明 |
|------|------|--------|------|
| `--compute-distance-matrix` | flag | enabled | 計算距離矩陣 |
| `--no-distance-matrix` | flag | - | 禁用距離矩陣 |
| `--output-distance-matrix` | flag | enabled | 輸出 CSV |
| `--output-strand-distance-matrices` | flag | enabled | 輸出 strand-specific |
| `--distance-metric` | string | NHD | 距離度量：NHD, L1, L2, CORR, JACCARD |
| `--min-common-coverage` | int | 3 | C_min 值 |
| `--nan-distance-strategy` | string | MAX_DIST | MAX_DIST 或 SKIP |
| `--max-distance-value` | double | 1.0 | MAX_DIST 的值 |

### 3.3 使用範例

```bash
# 預設 NHD，C_min=3
./inter_sub_mod --tumor-bam tumor.bam --normal-bam normal.bam \
    --reference ref.fa --vcf snv.vcf --output-dir output

# 使用 L2 距離，C_min=5
./inter_sub_mod ... --distance-metric L2 --min-common-coverage 5

# 使用 Correlation 距離
./inter_sub_mod ... --distance-metric CORR

# 禁用距離矩陣輸出
./inter_sub_mod ... --no-distance-matrix
```

---

## 4. 輸出檔案格式

### 4.1 distance_matrix.csv

對稱 N×N 矩陣，CSV 格式：

```csv
read_id,0,1,2,3,...
0,0.000000,0.333333,1.000000,1.000000,...
1,0.333333,0.000000,0.500000,1.000000,...
2,1.000000,0.500000,0.000000,0.250000,...
...
```

- 第一行：read_id 標題
- 對角線：0（自己與自己的距離）
- 1.0：可能為真實最大距離或因 C_min 不足而設定

### 4.2 distance_stats.txt

統計摘要：

```
Distance Matrix Statistics
==========================

Region ID: 0
Number of reads: 201
Metric: NHD
Min common coverage (C_min): 3

Valid pairs: 1075
Invalid pairs (insufficient overlap): 19025
Valid pair ratio: 5.3%
Average common coverage: 6.62

Distance Statistics:
  Min: 0.0000
  Max: 1.0000
  Mean: 0.9655
  Std Dev: 0.1573
  25th percentile: 1.0000
  Median: 1.0000
  75th percentile: 1.0000
```

---

## 5. 單元測試

### 5.1 測試檔案

`tests/test_distance_matrix.cpp`

### 5.2 測試案例（15 個，全部通過）

| 測試名稱 | 說明 |
|----------|------|
| NHD_BasicComputation | NHD 基本計算正確性 |
| NHD_MinCommonCoverage | C_min 閾值處理 |
| L1_BasicComputation | L1 距離計算 |
| L2_BasicComputation | L2 距離計算 |
| CORR_BasicComputation | 相關性距離計算 |
| JACCARD_BasicComputation | Jaccard 距離計算 |
| StrandSpecific_BasicComputation | Strand 分離計算 |
| DistanceConfig_AllOptions | 配置選項測試 |
| EmptyMatrix | 空矩陣處理 |
| SingleRead | 單一 read 處理 |
| AllMissing | 全缺失值處理 |
| Statistics | 統計計算正確性 |
| WriteCSV | CSV 輸出格式 |
| MetricToString | 度量類型轉換 |
| StringToMetric | 字串解析 |

### 5.3 執行測試

```bash
cd build && ./bin/run_tests --gtest_filter="DistanceMatrix*:DistanceCalculator*"
```

輸出：
```
[==========] Running 15 tests from 2 test suites.
...
[  PASSED  ] 15 tests.
```

---

## 6. 整合測試結果

### 6.1 測試環境

- **模式**：all-with-w1000（±1000bp 窗口，no-filter）
- **執行緒**：32
- **VCF**：HCC1395 filtered_snv_tp.vcf.gz

### 6.2 輸出統計

| 項目 | 數量 |
|------|------|
| Regions 處理 | 30,490 |
| 距離矩陣輸出 | 30,484 |
| Forward strand 矩陣 | 30,484 |
| Reverse strand 矩陣 | 30,484 |

### 6.3 效能觀察

- 每個 region 處理時間：40-80 ms（含距離計算）
- 距離矩陣大小：通常 80-200 reads
- 有效 read pair 比例：5-20%（取決於 CpG 覆蓋）

---

## 7. 架構設計考量

### 7.1 記憶體效率

- 使用 Eigen::MatrixXd 儲存（對稱矩陣只計算上三角）
- N=200 reads: ~320KB 記憶體
- 每個 region 獨立處理，完成後立即輸出

### 7.2 計算效率

- OpenMP 平行化 region 層級
- 單 region 內 pairwise 計算使用 `#pragma omp parallel for`
- O(N² × M) 複雜度，N=reads, M=CpGs

### 7.3 可擴展性

- 新增距離度量只需在 `calculate_distance_impl()` 添加 case
- 支援自定義 NaN 處理策略
- 配置參數化便於實驗比較

---

## 8. 已知限制與未來改進

### 8.1 當前限制

1. **相關性距離穩定性**：當所有值相近時變異數極小，需注意
2. **大矩陣輸出**：N>500 時 CSV 檔案較大，考慮壓縮格式
3. **記憶體使用**：極高深度區域可能需要分批處理

### 8.2 建議改進

1. 支援 binary 矩陣的 SIMD 加速（XOR + popcount）
2. 支援 NPZ 格式輸出
3. 增加 distance histogram 統計
4. 支援動態 C_min（基於 coverage 分佈）

---

## 9. 變更檔案清單

| 檔案 | 變更說明 |
|------|----------|
| `include/core/DistanceMatrix.hpp` | 新增 DistanceConfig, DistanceCalculator |
| `src/core/DistanceMatrix.cpp` | 實作所有距離計算方法 |
| `include/core/Config.hpp` | 新增距離計算配置欄位 |
| `include/utils/ArgParser.hpp` | 新增命令列參數 |
| `include/io/RegionWriter.hpp` | 新增距離矩陣輸出方法 |
| `src/io/RegionWriter.cpp` | 實作距離矩陣輸出 |
| `include/core/RegionProcessor.hpp` | 新增距離計算相關成員 |
| `src/core/RegionProcessor.cpp` | 整合距離計算流程 |
| `tests/test_distance_matrix.cpp` | 新增 15 個單元測試 |
| `CMakeLists.txt` | 添加測試檔案 |

---

## 10. 結論

本次開發成功實作了高效、可配置的 Read-Read 距離矩陣計算功能。通過完整的單元測試和大規模整合測試驗證，確保了計算的正確性和穩定性。所有可能影響計算結果的細節都已參數化，便於後續實驗比較不同設定對聚類分析的影響。

---

**附錄：執行指令範例**

```bash
# 執行 all-with-w1000 模式測試
./scripts/run_full_vcf_test.sh --mode all-with-w1000 --threads 32

# 執行單元測試
cd build && ./bin/run_tests --gtest_filter="DistanceMatrix*"
```

