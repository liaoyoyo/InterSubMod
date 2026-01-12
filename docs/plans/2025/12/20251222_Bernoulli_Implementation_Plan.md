# 20251222_Bernoulli_Implementation_Plan

基於先前 [可行性分析報告](20251222_Bernoulli.md) 的結論，本文件詳細規劃在 `InterSubMod` 專案中實作 **Bernoulli 距離矩陣計算方法** 的開發步驟、程式碼修改細節與驗證計畫。

## 1. 概述 (Overview)

目標是新增一種距離計算方法，利用 Read 上的原始甲基化機率 (0.0-1.0) 計算期望不一致率，並引入信心加權機制，以降低低信心位點 (p ≈ 0.5) 對分群結果的干擾。

*   **新增方法名稱**: `BERNOULLI`
*   **預期效益**: 提升 Nanopore 甲基化數據分群的魯棒性，減少假陽性分群。

## 2. 方法獨立性與影響分析 (Independence & Impact Analysis)

**Q: 這個實作會影響到底層的數值計算嗎？**
**A: 不會。這是一個完全獨立的距離計算選項。**

*   **底層數據不變**：程式最底層的 `MethylationMatrix` 建構過程（從 BAM 讀取 reads 並解析 MM/ML tags 轉為機率值）**完全保持不變**。系統原本就會儲存所有位點的原始浮點數機率 (`raw_matrix`)。
*   **計算邏輯獨立**：Bernoulli 方法僅是在「計算兩條 Read 之間的距離」這個步驟時，使用不同的公式。它與現有的 NHD (Hamming)、L1、L2 等方法是平行的關係。
*   **參數控制**：完全透過 CLI 參數 `--distance-metric` 控制。
    *   若未指定或指定 `NHD`，程式行為與數值與修改前 **100% 相同**。
    *   只有指定 `--distance-metric BERNOULLI` 時才會執行新的程式碼路徑。

## 3. 系統架構修改 (Architecture Changes)

### 3.1 定義新的距離類型
*   **檔案**: `include/core/Types.hpp`
*   **修改內容**:
    *   在 `enum class DistanceMetricType` 中新增枚舉值 `BERNOULLI`。

```cpp
enum class DistanceMetricType { NHD, L1, L2, CORR, JACCARD, BERNOULLI };
```

### 3.2 更新設定與參數解析
*   **檔案**: `src/core/DistanceMatrix.cpp` (及相關 Header)
*   **修改內容**:
    *   更新 `DistanceCalculator::string_to_metric` 函數，支援從 CLI 參數字串 `"BERNOULLI"` 解析為 `DistanceMetricType::BERNOULLI`。
    *   更新 `DistanceCalculator::metric_to_string` 函數，支援輸出 `"BERNOULLI"` 字串。

## 4. 核心邏輯實作 (Core Logic Implementation)

核心計算邏輯將實作於 `src/core/DistanceMatrix.cpp`。

### 4.1 新增計算函數 `calculate_bernoulli`

在 `DistanceMatrix.cpp` 的匿名命名空間 (Translation unit local) 中新增以下函數：

```cpp
/**
 * @brief Calculate Bernoulli Distance (Expected Disagreement with Confidence Weighting).
 *
 * Formula:
 *   delta(p, q) = p(1-q) + (1-p)q
 *   weight(p) = 2 * |p - 0.5|
 *   w_k = weight(p_i) * weight(p_j)
 *   Dist = sum(w_k * delta_k) / sum(w_k)
 *
 * @param raw_row_i Probability values for read i
 * @param raw_row_j Probability values for read j
 * @param min_cov Minimum number of common valid sites (used for validity check)
 * @param[out] common_count Number of common valid sites (valid means not NaN)
 * @return Distance value [0, 1], or -1.0 if valid pairs < min_cov or total weight is 0
 */
static double calculate_bernoulli(const Eigen::VectorXd& raw_row_i, const Eigen::VectorXd& raw_row_j, int min_cov,
                                  int& common_count) {
    common_count = 0;
    double sum_weighted_diff = 0.0;
    double sum_weights = 0.0;
    const int n = raw_row_i.size();

    auto weight_func = [](double p) {
        return 2.0 * std::abs(p - 0.5);
    };

    for (int k = 0; k < n; ++k) {
        double p_i = raw_row_i(k);
        double p_j = raw_row_j(k);

        // Check for missing data (NaN)
        if (!std::isnan(p_i) && !std::isnan(p_j)) {
            common_count++;
            
            // 1. Calculate Confidence Weights
            double w_i = weight_func(p_i);
            double w_j = weight_func(p_j);
            double w_k = w_i * w_j;

            // Optional: Skip if weight is effectively zero to save compute? 
            // No, keep it simple for now, 0 weight adds 0 to sums.

            // 2. Calculate Expected Disagreement (Bernoulli difference)
            // delta = P(diff) = p_i(1-p_j) + (1-p_i)p_j
            double delta = p_i * (1.0 - p_j) + (1.0 - p_i) * p_j;

            // 3. Accumulate
            sum_weighted_diff += w_k * delta;
            sum_weights += w_k;
        }
    }

    // Standard overlap check
    if (common_count < min_cov) {
        return -1.0;
    }

    // Edge case: If total weight is too small (e.g. all overlapping sites are p=0.5)
    // We treat this as "no information for comparison" -> Invalid distance
    if (sum_weights < 1e-9) {
        return -1.0; 
    }

    // Normalize
    return sum_weighted_diff / sum_weights;
}
```

### 4.2 整合至 Dispatcher

*   **檔案**: `src/core/DistanceMatrix.cpp`
*   **修改函數**: `calculate_distance_impl`
*   **修改內容**:
    *   在 `switch (config.metric)` 中新增 `case DistanceMetricType::BERNOULLI:`。
    *   呼叫 `calculate_bernoulli(mat.raw_matrix.row(row_i), mat.raw_matrix.row(row_j), ...)`。

## 5. 細節與參數注意事項 (Review & Edge Cases)

### 5.1 參數 `min_cov` 的定義
**問題**: `min_cov` 通常指的是「共同覆蓋的位點數」。在 Bernoulli 方法中，若兩個 Read 共同覆蓋了 10 個位點，但這 10 個位點的 $p$ 都是 0.5，則權重和為 0。
**決策**:
*   `min_cov` 依然檢查「非 NaN 的共同位點數」 (Physical Overlap)。
*   新增檢查：若 `sum_weights` 趨近於 0，視為「有效資訊不足」，回傳 -1.0 (Invalid)。這與覆蓋度不足在邏輯上是一致的（無法判斷距離）。

### 5.2 浮點數比較
**細節**: 使用 `1e-9` 作為 epsilon 來判斷 `sum_weights` 是否為 0，避免除以零錯誤。

### 5.3 效能考量
*   雖然涉及浮點數運算，但複雜度仍為 $O(N)$。目前的 OpenMP 平行化架構 (`#pragma omp parallel for`) 對此應能良好支援。
*   `weight_func` 是無狀態的，可考慮 inline。

## 6. 使用者介面 (User Interface)

使用者將透過 CLI 參數啟用此功能：
```bash
# Example usage
./InterSubMod --bam input.bam --vcf somatic.vcf --ref ref.fa --distance-metric bernoulli
```
或是同時計算多種：
```bash
... --distance-metric nhd,bernoulli
```
(現有 `Config::parse_args` 邏輯應已支援逗號分隔解析，只需確認 `string_to_metric` 更新即可)

## 7. 驗證計畫 (Verification Plan)

### 7.1 單元測試 (Unit Logic Check)
透過新增 temporary test code 或在 debug mode 下驗證特定 case：
*   **Case A (完全一致)**: $p_i=[0, 1]$, $p_j=[0, 1]$ → $\delta=[0, 0]$, $D=0$。
*   **Case B (完全相反)**: $p_i=[0, 1]$, $p_j=[1, 0]$ → $\delta=[1, 1]$, $D=1$。
*   **Case C (模糊位點)**: $p_i=[0.5]$, $p_j=[0]$ → $w_i=0$, $w=0$ → 忽略此位點。
*   **Case D (部分模糊)**: $p_i=[0, 0.5]$, $p_j=[1, 0.5]$ → Site 2 權重為 0，只看 Site 1 → $D=1$。

### 7.2 整合測試
執行 `scripts/run_full_vcf_test.sh`，增加 `--mode all-with-w1000` (或自訂 mode)，檢查：
1.  程式是否正常執行無 Crash。
2.  `output/.../distances_BERNOULLI.csv` 是否產生。
3.  比較 `distances_NHD.csv` 與 `distances_BERNOULLI.csv` 的熱圖：
    *   預期：Bernoulli 的熱圖中，模糊區域的雜訊應減少，群聚結構應更清晰（如果數據本身有結構）。
