# 距離計算方法分析與比較 (Distance Methods Analysis)

本文檔詳細說明 InterSubMod 專案中用於計算 Read-Read 甲基化距離的各種方法，比較其數學定義、特性、限制，並提供使用建議。

## 1. 概述

在甲基化數據分析中，計算兩個 Read 之間的距離（差異度）是進行聚類分析（Clustering）和熱圖繪製（Heatmap）的基礎。由於甲基化測序數據具有稀疏性（Sparsity），即兩個 Read 之間可能只有部分重疊的 CpG 位點，因此所有距離計算方法都必須處理缺失值（Missing Data）。

程式目前實作了以下五種距離計算方法：
1. **NHD** (Normalized Hamming Distance) - 預設方法
2. **L1** (Manhattan / Mean Absolute Difference)
3. **L2** (Euclidean / Root Mean Square Deviation)
4. **CORR** (Pearson Correlation Distance)
5. **JACCARD** (Jaccard Distance)

---

## 2. 詳細方法說明

### 2.1 Normalized Hamming Distance (NHD)
這是目前程式的**預設方法**。

*   **數學定義**：
    $$ D_{NHD} = \frac{\text{差異位點數}}{\text{共同有效位點數}} $$
    只考慮兩個 Read 都有覆蓋到的 CpG 位點（Common Valid Sites）。
    *   輸入：二元化甲基化狀態 (0: 未甲基化, 1: 甲基化)。
*   **特性**：
    *   直觀地表示「不一致率」。
    *   特定於二元數據（Binary Data）。
    *   範圍：[0, 1]。0 表示完全一致，1 表示完全相反。
*   **程式實作細節**：
    *   忽略缺失值 (-1)。
    *   若共同位點數小於設定閾值（`min_cov`），則視為距離無效。

### 2.2 L1 Distance (Manhattan / Mean Absolute Difference)
其實作為「平均絕對誤差」實作。

*   **數學定義**：
    $$ D_{L1} = \frac{\sum |p_i - p_j|}{N} $$
    其中 $N$ 是共同有效位點數。
*   **特性**：
    *   使用浮點數原始甲基化值（若有）。對於二元數據，結果與 NHD 相同。
    *   對異常值（Outliers）的敏感度低於 L2。
    *   範圍：[0, 1]（假設輸入值在 0-1 之間）。
*   **程式實作細節**：
    *   標準化了長度（除以 $N$），因此不會因為重疊區域長短而影響距離，這點非常重要。

### 2.3 L2 Distance (Euclidean / RMSD)
其實作為「均方根誤差」(RMSD) 實作。

*   **數學定義**：
    $$ D_{L2} = \sqrt{\frac{\sum (p_i - p_j)^2}{N}} $$
*   **特性**：
    *   給予較大的差異更大的權重。
    *   幾何意義上的直線距離（經長度歸一化）。
*   **限制**：
    *   對於二元數據，$(0-1)^2 = |0-1|$，因此在純二元輸入下，其單調性與 NHD 相同，但在數學數值上是 NHD 的平方根（$D_{L2} = \sqrt{D_{NHD}}$）。

### 2.4 Pearson Correlation Distance (CORR)
基於皮爾森相關係數。

*   **數學定義**：
    $$ D_{CORR} = \frac{1 - r}{2} $$
    其中 $r$ 是兩個 Read 在重疊位點上的皮爾森相關係數。
*   **特性**：
    *   衡量的是「變化趨勢」的一致性，而非數值絕對大小。
    *   範圍：[0, 1]。0 表示完全正相關 ($r=1$)，0.5 表示無相關 ($r=0$)，1 表示完全負相關 ($r=-1$)。
*   **限制**：
    *   **需要變異量**：如果一個 Read 在重疊區域內的值全為 0 或全為 1（變異數為 0），則無法計算相關係數。
    *   需要較多的重疊位點（程式限制至少 3 個）才有統計意義。
    *   對於甲基化這種常出現全 0 或全 1 片段的數據，容易產生計算失敗或極端值。

### 2.5 Jaccard Distance
集合相似度的度量。

*   **數學定義**：
    $$ D_{Jaccard} = 1 - \frac{|A \cap B|}{|A \cup B|} $$
    其中 $A, B$ 是兩個 Read 中**甲基化位點 (1)** 的集合。
*   **特性**：
    *   **只關注甲基化位點**（預設設定）：認為「共同甲基化」才是相似，而「共同未甲基化」不提供相似度資訊。
    *   這在稀疏矩陣（Sparse Matrix）分析中很常見，但在甲基化分析中可能有爭議（因為未甲基化也是一種明確的生物學狀態）。
*   **程式實作細節**：
    *   程式支援 `include_unmeth` 選項，若開啟則將未甲基化位點也視為特徵，此時行為會接近 NHD。

---

## 3. 比較與優缺點總結

| 方法 | 適用情境 | 優點 | 缺點 |
| :--- | :--- | :--- | :--- |
| **NHD** (預設) | 一般甲基化二元數據 | **最直觀**，物理意義明確（不一致比例），對長度不敏感 | 僅適用於二元數據（0/1） |
| **L1** | 連續型/機率型甲基化值 | 與 NHD 兼容，對異常值魯棒 | 在二元數據下與 NHD 重複 |
| **L2** | 強調大差異的情境 | 數學性質良好 | 在二元數據下是 NHD 的變形，計算開銷稍大 |
| **CORR** | 尋找趨勢相似性 (Pattern) | 忽略絕對數值，關注升降趨勢 | **極不推薦**用於短 Read。因常遇變異數為 0 導致大量無效計算 |
| **JACCARD** | 關注「有甲基化」的特徵 | 忽略「共同缺失」或「共同背景(0)」 | 忽略了「共同未甲基化」的生物學意義 (Unmethylated state is informative) |

---

## 4. 建議 (Recommendations)

### 最佳推薦：NHD (Normalized Hamming Distance)
*   **理由**：
    1.  BS-seq / Nanopore 經過 call methylation 後通常是二元數據 (0/1)。
    2.  NHD 直接反映了兩個 DNA 片段在重疊區域的不一致程度，生物學解釋性最強。
    3.  程式實作已處理了覆蓋度歸一化，不會受 Read 長短影響。

### 次佳選擇：Jaccard (特定情境)
*   **理由**：
    *   如果你只關心「甲基化發生」的事件，而認為未甲基化只是背景，則 Jaccard 較合適。例如在甲基化非常稀疏的區域，尋找共同甲基化的 Read 聚類。

### 不建議：Pearson Correlation (CORR)
*   **理由**：
    *   單條 Read 的甲基化模式通常是階梯狀或片段狀的，且長度短，變異數常為 0（例如：整段都沒甲基化）。這會導致相關係數無法計算或產生誤導性結果（Artifacts）。

### 參數建議
*   **最小覆蓋度 (Min Common Coverage)**：建議至少設為 **5~10**。
    *   若設太低（如 1 或 2），隨機的一致性會導致距離矩陣充滿雜訊。
    *   `tools/plot_distance_heatmap.py` 預設過濾 `min_reads`，底層 C++ 計算也有 `min_cov` 參數。

---

## 5. 程式碼參考
相關實作位於 `src/core/DistanceMatrix.cpp`：
*   `calculate_nhd`: 正規化漢明距離
*   `calculate_l1`: 正規化 L1
*   `calculate_l2`: RMSD
*   `calculate_corr`: 皮爾森相關
*   `calculate_jaccard`: Jaccard 距離
