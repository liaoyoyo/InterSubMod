# Read-Read 甲基化距離矩陣實作說明

## 1. 概述 (Overview)

本模組負責計算同一區域 (Region of Interest) 內，不同 Reads 之間的甲基化模式差異。此距離矩陣是後續階層聚類 (Hierarchical Clustering) 與單倍型分析 (Haplotype Analysis) 的基礎。

**核心目標**：

1. 量化 Read 與 Read 之間的甲基化相似度。
2. 支援多種距離演算法 (Hamming, Euclidean)。
3. **區分正反股 (Strand-aware)**：由於 DNA 甲基化在正反股可能具有不對稱性，且定序過程中的 PCR 偏差可能不同，因此需分別計算正股 (Forward) 與反股 (Reverse) 的距離矩陣。
4. 處理缺失值 (Missing Data)：Reads 覆蓋範圍不同，僅在共同覆蓋的 CpG 位點上計算距離。

## 2. 輸入與輸出 (Input & Output)

### 2.1 輸入資料

- **MethylationMatrix**: 包含 $N$ 個 Reads 與 $M$ 個 CpG sites 的甲基化狀態。
  - `binary_matrix`: $N \times M$ (值為 1, 0, -1)。
  - `raw_matrix`: $N \times M$ (值為 0.0-1.0, NaN)。
- **ReadInfo**: 每個 Read 的元數據，關鍵欄位為 `strand` (FORWARD/REVERSE)。

### 2.2 輸出資料

- **DistanceMatrix**: $N \times N$ 對稱矩陣，儲存 Read 間的距離。
  - 由於需分開正反股，實際操作上會產生兩個矩陣，或在一個大矩陣中將不同股的距離設為無限大 (Infinity)。
  - 建議實作：**產生兩個獨立的 DistanceMatrix 物件**，分別對應 Forward 與 Reverse Reads。

## 3. 演算法詳細設計 (Algorithm Design)

### 3.1 距離度量 (Distance Metrics)

#### A. 正規化漢明距離 (Normalized Hamming Distance, NHD)

適用於二值化資料 (`binary_matrix`)。

$$ D_{Hamming}(R_i, R_j) = \frac{\sum_{k \in Common} |M_{ik} - M_{jk}|}{|Common|} $$

其中 $Common$ 為兩個 Read $i, j$ 都具有有效甲基化數值 (非 -1) 的 CpG 位點集合。

- 若 $M_{ik} = M_{jk}$，距離貢獻為 0。
- 若 $M_{ik} \neq M_{jk}$，距離貢獻為 1。
- 最終除以共同位點數進行正規化，範圍 [0, 1]。

#### B. 歐幾里得距離 (Euclidean Distance)

適用於原始機率資料 (`raw_matrix`)。

$$ D_{Euclidean}(R_i, R_j) = \sqrt{\frac{\sum_{k \in Common} (P_{ik} - P_{jk})^2}{|Common|}} $$

- 同樣需對共同位點數進行正規化 (Root Mean Square Difference)，以避免覆蓋長度影響距離數值。

### 3.2 缺失值處理 (Missing Data Handling)

由於定序深度與覆蓋範圍的變異，任意兩個 Reads 可能只有極少數的共同 CpG 位點。

- **最小覆蓋閾值 (Min Common Coverage)**：設定參數 `min_common_cov` (例如 3 或 5)。
- 若 $|Common| < min\_common\_cov$，則這兩個 Reads 之間的距離視為 **無效 (Invalid/Infinity)**。
- 在聚類時，這些無效距離將導致這兩個 Reads 不會被聚在一起，或使用特殊策略處理。

## 4. 實作流程 (Implementation Steps)

### 步驟 1: 讀取分群 (Strand Splitting)

將 `MethylationMatrix` 中的 Reads 依據 `ReadInfo.strand` 分為兩組索引：

- `forward_indices`: [0, 2, 5, ...]
- `reverse_indices`: [1, 3, 4, ...]

### 步驟 2: 矩陣初始化

建立兩個 $N_f \times N_f$ 與 $N_r \times N_r$ 的距離矩陣，初始化為 NaN 或 MAX_DOUBLE。

### 步驟 3: 平行化計算 (Parallel Computation)

使用 OpenMP 對矩陣的行 (Row) 進行平行化遍歷。

```cpp
// 虛擬碼範例
#pragma omp parallel for schedule(dynamic)
for (int i = 0; i < subset_size; ++i) {
    for (int j = i + 1; j < subset_size; ++j) {
        double dist = compute_pair_distance(read_i, read_j);
        matrix(i, j) = dist;
        matrix(j, i) = dist; // 對稱
    }
}
```

### 步驟 4: 核心距離計算函式

針對每一對 Read，遍歷所有 CpG sites。

**最佳化技巧**：

- **Bitwise Operation**: 若使用二值矩陣，可將每個 Read 的狀態壓縮為 `std::vector<bool>` 或 `boost::dynamic_bitset`，利用 XOR 快速計算差異，AND 計算共同覆蓋。
  - Mask: `valid_mask_i & valid_mask_j` (共同有效位點)
  - Diff: `(value_i ^ value_j) & common_mask`
  - Distance = `popcount(Diff) / popcount(Common)`
- **Early Exit**: 若共同位點數累計尚未達到閾值且已遍歷完，直接回傳無效。

## 5. 程式碼架構 (Code Structure)

建議在 `src/core/DistanceMatrix.cpp` 中實作以下類別與方法：

```cpp
namespace InterSubMod {

class DistanceCalculator {
public:
    // 設定參數
    struct Config {
        DistanceMetricType metric;
        int min_common_coverage;
        double nan_value; // 用於無效距離的數值 (例如 1.0 或 -1.0)
    };

    DistanceCalculator(Config config);

    // 主要入口：計算並回傳兩個矩陣 (Forward, Reverse)
    std::pair<DistanceMatrix, DistanceMatrix> compute_split_matrices(
        const MethylationMatrix& meth_mat,
        const std::vector<ReadInfo>& read_infos
    );

private:
    Config config_;

    // 計算單一子集的距離矩陣
    DistanceMatrix compute_subset(
        const MethylationMatrix& meth_mat,
        const std::vector<int>& row_indices
    );

    // 計算兩列之間的距離
    double calculate_pairwise_distance(
        const Eigen::MatrixXi& binary_data, // 或 MatrixXd
        int row_i,
        int row_j
    );
};

}
```

## 6. 記憶體與效能分析 (Performance & Memory)

- **記憶體**:
  - 距離矩陣為 $N \times N$ double。
  - 若 $N=1000$ (高深度)，大小約 8MB，非常輕量。
  - 主要記憶體消耗仍在 `MethylationMatrix` 本身。
- **計算複雜度**:
  - $O(N^2 \cdot M)$，其中 $N$ 為 Read 數，$M$ 為 CpG 數。
  - 由於 $N$ 通常不大 (<100)，且 $M$ 有限 (<500)，計算非常快速。
  - OpenMP 平行化可線性加速。

## 7. 註解與文件 (Comments & Documentation)

所有函式需包含 Doxygen 風格註解，說明：

- `@param` 輸入參數格式
- `@return` 回傳值意義
- `@note` 特殊處理 (如 NaN 行為)
