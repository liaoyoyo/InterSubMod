# UPGMA Clustering 與 CramersV 改進方案評估報告

**專案**: InterSubMod - 甲基化分析  
**分析日期**: 2026-01-23  
**目的**: 評估 UPGMA clustering 截斷標準與 CramersV 計算改進方案

---

## 1. 現行 UPGMA Clustering 截斷機制

### 1.1 截斷標準

```cpp
// RegionProcessor.cpp:1020-1021
auto [best_k, cluster_labels] = TreeCutter::find_optimal_clusters(
    tree, all_dist.dist_matrix, 2, std::min(6, read_list.size() / 2));
```

| 參數 | 值 | 說明 |
|------|-----|------|
| `min_k` | 2 | 最少 2 個 cluster |
| `max_k` | min(6, n_reads/2) | 最多 6 個，但不超過 reads 的一半 |
| **選擇標準** | **Silhouette Score 最高** | 計算每個 k 的平均 Silhouette，選最佳 |

### 1.2 Silhouette Score 計算

```cpp
// HierarchicalClustering.cpp:425-466
for (int k = min_k; k <= max_k; ++k) {
    for (每個 read i) {
        a(i) = 同 cluster 內平均距離
        b(i) = 最近其他 cluster 的平均距離
        s(i) = (b - a) / max(a, b)
    }
    avg_silhouette = mean(s(i))
    選擇最高的 k
}
```

### 1.3 截斷標準評估

| 評估項目 | 結果 | 說明 |
|----------|------|------|
| **方法合理性** | ✅ 合理 | Silhouette 是標準的 cluster 數量選擇方法 |
| **k 範圍限制** | ⚠️ 可調 | max=6 可能太小，某些複雜區域可能需要更多 clusters |
| **穩定性** | ⚠️ 關注 | 小樣本時 Silhouette 可能不穩定 |

---

## 2. UPGMA 替代方案

### 2.1 現有支援的 Linkage 方法

```cpp
switch (config_.method) {
    case LinkageMethod::UPGMA:     // 平均距離 (現行預設)
    case LinkageMethod::WARD:      // 最小化組內方差
    case LinkageMethod::SINGLE:    // 最近點距離
    case LinkageMethod::COMPLETE:  // 最遠點距離
}
```

### 2.2 方法比較

| 方法 | 優點 | 缺點 | 適用情境 |
|------|------|------|----------|
| **UPGMA** | 平衡、穩定 | 假設等速演化 | 甲基化分析 (現行) |
| **Ward** | 最小化組內變異 | 傾向產生相似大小的 cluster | 均衡分群需求 |
| **Single** | 發現細長形狀 | 容易 chaining | 離群點偵測 |
| **Complete** | 產生緊湊 cluster | 傾向等大小 cluster | 避免離群值影響 |

### 2.3 建議測試的替代方案

1. **Ward + Silhouette**: 可能產生更穩定的 cluster 大小
2. **UPGMA + Gap Statistic**: 比 Silhouette 更 robust 的 k 選擇
3. **DBSCAN**: 不需要預設 k，可自動偵測離群點

---

## 3. 固定 k 但改變截斷位置

### 3.1 概念說明

> **問題**: 同樣取 k=3 個 cluster，是否可以有不同的截斷位置？

**答案**: 在現行的 `cut_by_num_clusters` 實現中，截斷位置是**唯一確定的**。

```cpp
// HierarchicalClustering.cpp:360-404
std::vector<int> TreeCutter::cut_by_num_clusters(const Tree& tree, int num_clusters) {
    // 收集所有 merge heights，由高到低排序
    for (const auto& [h, node] : height_nodes) {
        // 找到剛好產生 num_clusters 的切割高度
        auto test_labels = cut_by_distance(tree, h * 2.0 + 0.0001);
        int n_clusters = max_element(test_labels) + 1;
        if (n_clusters >= num_clusters) {
            cut_height = h * 2.0;
            break;
        }
    }
    return cut_by_distance(tree, cut_height);
}
```

### 3.2 替代方案：動態截斷

| 方法 | 說明 | 實現難度 |
|------|------|----------|
| **Gap Statistics** | 比較實際與 null 分布的間隙 | 中 |
| **Elbow Method** | 找 SSE 下降曲線的轉折點 | 低 |
| **Calinski-Harabasz** | 另一種 cluster 質量指標 | 低 |
| **Dynamic Tree Cut** | 依據 tree 結構動態決定 | 高 |

---

## 4. CramersV 原始值輸出

### 4.1 現行問題

```cpp
// RegionProcessor.cpp:1092-1094
double v_alt = sig_result.global_alt.cramers_v_reliable ? sig_result.global_alt.cramers_v : 0.0;
double v_hp = sig_result.global_hp.cramers_v_reliable ? sig_result.global_hp.cramers_v : 0.0;
result.cramers_v = std::max(v_alt, v_hp);
```

> [!WARNING]
> 當 `cramers_v_reliable = false` 時，V 被設為 **0.0**，無法區分「真的無關聯」與「數據不足無法計算」。

### 4.2 建議修改方案

#### 方案 A：新增欄位區分

| 欄位 | 說明 |
|------|------|
| `CramersV_Raw` | 原始計算值，不管可靠性 |
| `CramersV_Reliable` | bool，是否可靠 |
| `CramersV` | 現行邏輯 (可靠時用原始，否則 0) |

#### 方案 B：新增欄位追蹤例外原因

```cpp
enum class CramersVStatus {
    VALID,              // 正常計算
    LOW_EXPECTED,       // >20% cells expected < 5
    SINGLE_CLUSTER,     // 只有 1 個 cluster
    SINGLE_LABEL,       // 只有 1 個 label (全 Ref/全 Alt)
    ZERO_N              // 樣本數為 0
};
```

---

## 5. 建議修改的程式碼區塊

### 5.1 修改 `GlobalTestResult` 結構

```cpp
// GlobalTest.hpp
struct GlobalTestResult {
    // 現有欄位
    double cramers_v = 0.0;
    bool cramers_v_reliable = false;
    
    // 新增欄位
    double cramers_v_raw = 0.0;           // 原始計算值
    std::string cramers_v_status = "";     // 例外狀態字串
};
```

### 5.2 修改 `MathUtils::cramers_v`

```cpp
double MathUtils::cramers_v(const std::vector<int>& table, int n_rows, int n_cols, 
                            bool& is_reliable, std::string& status) {
    if (n == 0) {
        is_reliable = false;
        status = "ZERO_N";
        return 0.0;  // 但仍返回計算值
    }
    
    int min_dim = std::min(n_rows, n_cols);
    if (min_dim <= 1) {
        is_reliable = false;
        status = "SINGLE_DIM";
        return 0.0;
    }
    
    // ... 計算 ...
    
    if (low_expected_ratio > 0.2) {
        is_reliable = false;
        status = "LOW_EXPECTED";
    } else {
        is_reliable = true;
        status = "VALID";
    }
    
    // 不管可靠性都返回計算值
    return std::clamp(v, 0.0, 1.0);
}
```

### 5.3 修改 CSV 輸出

```cpp
// RegionProcessor.cpp - significance_summary.csv
csv_file << "...,CramersV,CramersV_Raw,CramersV_Status,...\n";

csv_file << result.cramers_v << ","
         << result.cramers_v_raw << ","
         << result.cramers_v_status << ",...";
```

---

## 6. 例外狀況分類表

| 狀態碼 | 說明 | CramersV 處理 | 建議後處理 |
|--------|------|---------------|------------|
| `VALID` | 正常計算 | 使用原始值 | 直接使用 |
| `LOW_EXPECTED` | >20% cells expected<5 | 輸出原始值 + flag | 可考慮使用，但需謹慎 |
| `SINGLE_CLUSTER` | 只有 1 個 cluster | 設為 0 | 不適用 V 分析 |
| `SINGLE_LABEL` | 只有 1 個 label | 設為 0 | 不適用 V 分析 |
| `ZERO_N` | 無有效樣本 | 設為 0 | 跳過 |

---

## 7. 總結與建議

### 7.1 UPGMA Clustering

| 項目 | 評估 | 建議 |
|------|------|------|
| **Linkage 方法** | ✅ UPGMA 合理 | 可選擇性測試 Ward |
| **截斷標準** | ⚠️ Silhouette 可能不穩定 | 考慮 Gap Statistic |
| **k 範圍** | ⚠️ max=6 可能太小 | 考慮放寬到 8 或 10 |

### 7.2 CramersV 計算

| 項目 | 評估 | 建議 |
|------|------|------|
| **二元輸出問題** | ❌ 嚴重 | **必須修改**，輸出原始值 |
| **可靠性門檻** | ⚠️ 偏嚴 | 考慮放寬到 30% |
| **例外狀態追蹤** | ❌ 缺失 | **必須新增** status 欄位 |

### 7.3 修改優先順序

1. **高優先**: 新增 `CramersV_Raw` 和 `CramersV_Status` 欄位
2. **中優先**: 調整可靠性門檻 (20% → 30%)
3. **低優先**: 測試替代 Linkage 方法

---

## 8. 實施計劃

### 階段 1: 診斷與資料收集 (1-2 天)

- [ ] 修改程式輸出 `CramersV_Raw` 和 `CramersV_Status`
- [ ] 重新執行分析，收集例外狀態分布
- [ ] 分析 TP vs FP 在各狀態下的差異

### 階段 2: 參數調整 (1 天)

- [ ] 測試放寬 expected < 5 門檻 (20% → 30% → 40%)
- [ ] 評估對 F1-score 的影響

### 階段 3: 方法驗證 (2-3 天)

- [ ] 測試 Ward linkage
- [ ] 測試 Gap Statistic 選擇 k
- [ ] 比較不同組合的效果

---

**報告生成時間**: 2026-01-23  
**相關程式碼**: 
- `src/core/HierarchicalClustering.cpp`
- `src/core/MathUtils.cpp`  
- `src/core/RegionProcessor.cpp`
