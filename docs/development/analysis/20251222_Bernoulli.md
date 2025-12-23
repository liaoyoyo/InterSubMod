# Bernoulli 距離矩陣計算方法可行性分析與確認

**日期**: 2025-12-22
**狀態**: 已確認 (Verified)
**目標**: 評估在 InterSubMod 專案中加入「Bernoulli 距離」以改善低信心位點對分群結果影響的合理性與可行性。

---

## 1. 專案目標與背景 (Project Context)

本專案 `InterSubMod` 的核心目標是利用長讀取 (Long Read) 測序數據中的**甲基化圖樣 (Methylation Patterns)** 來識別生物樣本中的**子克隆 (Subclones)**。
*   **關鍵挑戰**：Nanopore 等測序技術產生的甲基化判讀是機率性的 (Probabilistic)，而非絕對的二元值。
*   **現狀問題**：目前的 Hamming Distance (NHD) 或 L1 可能對信心度低 (p ≈ 0.5) 的位點過度敏感，導致這些雜訊位點主導了距離計算，產生虛假的群聚或「假分裂」。

## 2. Bernoulli 距離方法分析 (Method Analysis)

根據先前的筆記與討論，擬採用的「Bernoulli 方法」包含兩個核心組件：**期望不一致率 (Expected Disagreement)** 與 **信心加權 (Confidence Weighting)**。

### 2.1 數學定義

對於兩個 Read $i, j$ 在位點 $k$ 上的甲基化機率 $p_{ik}, p_{jk}$：

1.  **差異量度 (Difference Metric)**:
    採用 Bernoulli 分布的期望不一致機率：
    $$ \delta(p_{ik}, p_{jk}) = p_{ik}(1 - p_{jk}) + (1 - p_{ik})p_{jk} $$
    *   意義：假設兩條 Read 在該位點的真實狀態是根據機率 $p$ 進行 Bernoulli 試驗，$\delta$ 即為兩次試驗結果**不同**的機率。
    *   性質：$\delta(0,1)=1$, $\delta(0.5,0.5)=0.5$, $\delta(0,0)=0$。

2.  **信心權重 (Confidence Weighting)**:
    為了抑制低信心位點的影響，引入權重函數 $c(p)$：
    $$ c(p) = 2 \cdot |p - 0.5| $$
    $$ w_{ijk} = c(p_{ik}) \cdot c(p_{jk}) $$
    *   性質：
        *   $p=0.5 \implies c(p)=0$ (完全忽略)。
        *   $p=0 \text{ or } 1 \implies c(p)=1$ (全權重)。
    *   **結合權重**：只有當**兩條 Read** 在該位點都具有高信心時，該位點的差異才會有顯著貢獻。

3.  **最終距離公式**:
    $$ D_{Bernoulli}(i, j) = \frac{\sum_{k \in K_{ij}} w_{ijk} \cdot \delta(p_{ik}, p_{jk})}{\sum_{k \in K_{ij}} w_{ijk}} $$
    其中 $K_{ij}$ 為兩條 Read 共同覆蓋的 CpG 位點集合。

### 2.2 合理性評估 (Reasonableness)

*   **統計意義 (Statistical Validity)**: **高度合理**。此方法直接利用了原始數據的機率特性，而非將其強制降維為 0/1。$\delta$ 函數在統計上是衡量兩個獨立 Bernoulli 分布差異的自然方式。
*   **雜訊抑制 (Noise Reduction)**: **極佳**。透過 $c(p)$ 權重，有效解決了 Nanopore 數據中常見的模糊區域 (fuzzy regions) 問題。相比於單純的閾值過濾（cutoff），加權法保留了更多資訊且更為平滑。
*   **相容性**: 此方法與現有的 NHD 方法在極端情況下 ($p \in \{0,1\}$) 是等價的，是一種更為泛化的距離定義。

## 3. 系統可行性與實作 (Feasibility & Implementation)

### 3.1 數據結構支援
檢查 `include/core/MethylationMatrix.hpp` 發現：
```cpp
class MethylationMatrix {
public:
    Eigen::MatrixXd raw_matrix;     ///< Raw methylation probabilities (0.0 - 1.0)
    // ...
};
```
*   **結論**: 系統已經儲存了浮點數的原始機率 (`raw_matrix`)，因此**不需要**修改底層數據讀取或儲存邏輯，即可直接實作此算法。

### 3.2 潛在邊界情況 (Edge Cases)
實作時需注意以下情況：
1.  **權重總和為 0**: 若兩條 Read 的共同覆蓋位點全都是 $p=0.5$ (導致 $w=0$)，則分母為 0。
    *   *解決方案*: 此情況應視為「資訊不足以判斷距離」，應回傳 `NaN` 或最大距離 (1.0)，並在後續 Clustering 中將其視為無效邊。
2.  **運算效能**: 涉及浮點數乘法與絕對值運算，比 Hamming 的 XOR 慢，但在現代 CPU 上差異應可忽略不計 ($\approx O(N)$)。

## 4. 總結與建議 (Conclusion)

**確認結果**: 加入 Bernoulli 距離矩陣計算方法是**完全合理且必要的**。它能從根本上提升分群演算法對低品質/模糊數據的魯棒性。

**實作建議步聚**:
1.  在 `DistanceMetricType` (Types.hpp) 中加入 `BERNOULLI`。
2.  在 `DistanceMatrix.cpp` 中實作 `calculate_bernoulli` 函數，依照上述公式。
3.  處理分母為 0 的情況（建議回傳 -1 或 NaN 供上層過濾）。
4.  在 CLI 參數中開放此選項作為進階使用者的選擇。

此文檔確認了該方法的理論基礎與實作可行性，建議立即進行開發。
