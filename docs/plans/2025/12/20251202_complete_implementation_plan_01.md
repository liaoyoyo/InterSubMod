<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# InterSubMod 聚類與演化分析完整實作方案

**版本**：v2.0  
**日期**：2025-12-02  
**基於使用者需求確認版**

---

## 1. 核心決策與實作策略

### 1.1 技術棧選擇

基於使用者需求，採用以下技術方案：

| 階段 | 實作語言 | 工具/套件 | 原因 |
|------|---------|----------|------|
| 距離矩陣計算 | **C++** | 現有 `DistanceMatrix.cpp` | ✅ 已完成，高效能平行化 |
| 層次聚類 | **C++** | **選項見 1.2** | 核心計算，需要高效能 |
| Bootstrap 驗證 | **C++** | 自行實作重抽樣邏輯 | 與聚類整合，避免重複 I/O |
| 演化樹建構 | **C++** | **選項見 1.3** | 研究癌症-甲基化共同演化 |
| Heatmap 視覺化 | **Python** | `seaborn.clustermap` | 精緻美觀，C++ 呼叫 |
| 統計分析 | **Python** | `scikit-bio`, `scipy.stats` | 成熟的統計套件，C++ 呼叫 |

---

### 1.2 C++ 層次聚類套件調研與選擇

#### 選項 A：自行實作 UPGMA / Ward's Method ⭐⭐⭐⭐ (推薦)

**實作複雜度**：中等（200-300 行程式碼）

**優點**：

- 完全掌控演算法，易於整合 Bootstrap
- 僅需距離矩陣作為輸入，已完成計算
- 可直接輸出 Newick 格式演化樹
- 無外部依賴

**缺點**：

- 需要自行驗證正確性

**實作架構**：

```cpp
namespace InterSubMod {

class HierarchicalClustering {
public:
    enum class LinkageMethod {
        UPGMA,       // 適合演化樹
        WARD,        // 適合一般聚類
        SINGLE,
        COMPLETE
    };
    
    struct ClusterNode {
        int left_child;   // -1 for leaf, else node index
        int right_child;
        double height;    // 分支長度
        int node_id;
        std::vector<int> leaf_indices; // 此 clade 包含的 reads
    };
    
    struct Tree {
        std::vector<ClusterNode> nodes;
        std::string to_newick() const;  // 輸出 Newick 格式
        void write_to_file(const std::string& path) const;
    };
    
    HierarchicalClustering(LinkageMethod method);
    
    // 主要入口：從距離矩陣建構樹
    Tree build_tree(const DistanceMatrix& dist_matrix);
    
private:
    LinkageMethod method_;
    double compute_cluster_distance(
        const Eigen::MatrixXd& dist,
        const std::vector<int>& cluster_a,
        const std::vector<int>& cluster_b
    );
};

} // namespace InterSubMod
```

**UPGMA 演算法簡述**：

1. 初始化：每個 read 是一個 cluster
2. 迭代：找到距離最小的兩個 clusters，合併
3. 更新距離矩陣（使用算術平均）
4. 記錄合併高度（分支長度 = distance / 2）
5. 重複直到所有 clusters 合併為一棵樹

**時間複雜度**：$O(N^3)$（未優化），$O(N^2 \log N)$（使用 priority queue）

---

#### 選項 B：使用 `mlpack` (Machine Learning Pack) ⭐⭐⭐

**套件資訊**：

- **網站**：<https://www.mlpack.org/>
- **License**：BSD（可商用）
- **依賴**：Armadillo（線性代數庫），OpenMP

**功能**：

- 提供多種聚類演算法（k-means, DBSCAN 等）
- **但不直接支援距離矩陣輸入的層次聚類**

**評估**：❌ 不適合本專案需求

---

#### 選項 C：使用 `Clustering` 頭文件庫 ⭐⭐

**資訊**：

- GitHub: <https://github.com/davidstutz/clustering> (小型專案)
- 支援：Agglomerative Clustering

**評估**：⚠️ 專案較小，維護度低，不如自行實作

---

#### 選項 D：呼叫 Python `scipy` 後解析結果 ⭐⭐⭐⭐⭐

**流程**：

1. C++ 輸出距離矩陣到檔案
2. 使用 `system()` 或 `pybind11` 呼叫 Python 腳本
3. Python 執行聚類並輸出 Newick
4. C++ 讀取結果繼續後續分析

**優點**：

- `scipy.cluster.hierarchy` 已經過充分測試
- 支援多種 linkage 方法
- 可直接輸出 dendrogram 與 Newick

**缺點**：

- 跨語言呼叫有一定開銷
- 需要 Python 環境

---

#### **🔴 建議方案**

**混合策略（最佳平衡）**：

1. **C++ 自行實作 UPGMA**：用於生產環境的快速聚類
2. **Python `scipy` 作為驗證與視覺化**：確保正確性，並生成 heatmap

**具體工作流程**：

```
C++ 主程式：
1. 計算距離矩陣（已完成）
2. 執行 UPGMA 建構演化樹
3. 執行 Bootstrap（重抽樣 CpG → 重新計算距離 → UPGMA）
4. 輸出 Newick 格式樹檔案 + Bootstrap support values
5. 呼叫 Python 腳本繪製 Clustermap
```

---

### 1.3 演化樹建構方法選擇

基於使用者需求（研究癌症-甲基化共同演化），需要支援演化樹構建。

#### UPGMA vs Neighbor Joining (NJ)

| 方法 | 假設 | 適用情境 | 實作複雜度 |
|------|------|----------|-----------|
| **UPGMA** | 演化速率恆定（molecular clock） | 樣本演化時間相近 | ⭐⭐ (容易) |
| **NJ** | 不假設速率恆定 | 演化速率差異大 | ⭐⭐⭐ (中等) |

**對於甲基化資料**：

- Reads 來自同一時間點的細胞群體
- 不存在真正的「時間演化」
- **UPGMA 較適合**，假設演化速率恆定在這裡可被解釋為「甲基化變化速率相近」

**🔴 建議**：優先實作 UPGMA，若需要可擴展 NJ

---

## 2. Bootstrap 驗證實作方案

### 2.1 Bootstrap 原理

**目的**：評估演化樹拓撲結構的穩定性

**步驟**：

1. **重抽樣**：對原始資料的 CpG 位點（columns）進行有放回抽樣 $N$ 次（例如 100 次）
2. **重複聚類**：每次重抽樣後，重新計算距離矩陣並建構演化樹
3. **計算支持度**：統計原始樹中每個節點（clade）在 Bootstrap 樹中出現的頻率

**範例**：

- 原始樹： `((A,B)95, (C,D)80)100`
- 數字 95, 80, 100 為 Bootstrap support values（百分比）
- 95% 表示在 100 次 Bootstrap 中，有 95 次出現此節點

---

### 2.2 C++ 實作架構

```cpp
namespace InterSubMod {

class BootstrapAnalyzer {
public:
    struct Config {
        int n_iterations = 100;      // Bootstrap 次數
        int random_seed = 42;         // 隨機種子
        int n_threads = 16;           // 平行化執行緒數
    };
    
    struct BootstrapResult {
        Tree original_tree;           // 原始樹
        std::vector<double> support_values; // 每個節點的支持度 [0-100]
        std::vector<Tree> bootstrap_trees; // （可選）儲存所有 Bootstrap 樹
    };
    
    BootstrapAnalyzer(Config config);
    
    // 主要入口
    BootstrapResult run_bootstrap(
        const MethylationMatrix& original_matrix,
        const DistanceConfig& dist_config,
        HierarchicalClustering::LinkageMethod linkage
    );
    
private:
    Config config_;
    std::mt19937 rng_;
    
    // 重抽樣 CpG 位點
    MethylationMatrix resample_columns(
        const MethylationMatrix& original,
        std::vector<int>& sampled_indices
    );
    
    // 比較兩棵樹的拓撲結構
    bool is_clade_present(
        const Tree& tree,
        const std::vector<int>& clade_leaves
    );
};

} // namespace InterSubMod
```

---

### 2.3 平行化策略

**方案 A：OpenMP 平行化 Bootstrap iterations**

```cpp
#pragma omp parallel for schedule(dynamic)
for (int i = 0; i < n_iterations; ++i) {
    // 每個執行緒獨立執行一次 Bootstrap
    MethylationMatrix resampled = resample_columns(original_matrix);
    DistanceMatrix dist = calculate_distance(resampled);
    Tree tree = build_tree(dist);
    
    #pragma omp critical
    {
        bootstrap_trees.push_back(tree);
    }
}
```

**預期加速比**：接近線性（64 執行緒可達到 50-60x 加速）

---

### 2.4 Bootstrap 支持度計算

**演算法**：

1. 對於原始樹的每個內部節點 $n$，記錄其包含的葉節點集合 $L_n$
2. 對於每棵 Bootstrap 樹，檢查是否存在節點包含完全相同的 $L_n$
3. 支持度 = （出現次數 / Bootstrap 總次數）× 100%

**實作注意事項**：

- 葉節點集合需要排序後比較（避免順序差異）
- 僅計算內部節點，葉節點支持度恆為 100%

---

## 3. 多標籤統計分析策略

### 3.1 標籤維度定義

根據使用者需求，分析以下標籤組合：

| 標籤維度 | 類別數 | 優先級 | 分析層級 |
|---------|--------|--------|---------|
| **HP (Haplotype)** | 6-7 種 | ⭐⭐⭐⭐⭐ | **最基本** |
| Tumor/Normal | 2 種 | ⭐⭐⭐⭐ | 次要 |
| ALT/REF | 2 種 | ⭐⭐⭐ | 次要 |
| Strand | 2 種 | ⭐⭐ | 輔助 |

**HP Tag 分類**：

- HP1, HP2：定相的主要單倍型 reads
- HP1-1, HP2-1：原於定相的主要單倍型並有經過ALT 初步認為是突變 reads
- HP3：經過ALT 但沒有定相或是混合兩種單倍型的reads
- Unphased：無法定相的reads
- None：無 HP 標籤的reads

---

### 3.2 分析層次設計

#### Level 1：單維度分析（必須執行）

**分析 HP 標籤與甲基化聚類的關聯**

**方法 1：PERMANOVA (Permutational MANOVA)**

- **問題**：不同 HP 類別的 reads，是否在甲基化距離矩陣上有顯著分離？
- **輸入**：距離矩陣 + HP 分組標籤
- **輸出**：F-statistic, p-value
- **解釋**：p < 0.05 表示 HP 類別顯著影響甲基化模式

**Python 實作**：

```python
from skbio.stats.distance import permanova

# dist_matrix: N x N 距離矩陣
# hp_labels: N 個 reads 的 HP 標籤
result = permanova(dist_matrix, hp_labels, permutations=999)
print(f"F-statistic: {result['test statistic']}")
print(f"p-value: {result['p-value']}")
```

---

**方法 2：演化樹分支富集分析 (Clade Enrichment)**

**目的**：檢驗特定 HP 類別是否在演化樹的某些分支中富集

**步驟**：

1. 識別演化樹的主要分支 (major clades)
2. 對於每個分支，統計各 HP 類別的 read 數量
3. 使用 **Fisher's Exact Test** 檢驗富集顯著性

**範例**：

```
演化樹分支 A 包含 20 個 reads：
- HP1: 15 個
- HP2: 5 個

全部 reads (100 個)：
- HP1: 40 個
- HP2: 60 個

Fisher's Exact Test:
            In Clade A | Not in Clade A
HP1            15      |      25
HP2             5      |      55

p-value < 0.001 → HP1 在此分支中顯著富集
```

---

#### Level 2：雙維度組合分析（條件允許時執行）

**分析 HP × Tumor/Normal 的交互作用**

**方法**：Two-way PERMANOVA

```python
# 建立組合標籤
combined_labels = [f"{hp}_{tumor}" for hp, tumor in zip(hp_labels, tumor_labels)]
result = permanova(dist_matrix, combined_labels, permutations=999)
```

**解釋**：

- 是否在 Tumor 樣本中，HP1 與 HP2 的甲基化差異更大？
- Normal 樣本中是否無差異？

---

#### Level 3：Strand 的整合與分離分析

**策略**：

1. **整合分析**（預設）：
   - 使用完整的 `distance_matrix.csv`（包含所有 reads）
   - 將 Strand 作為統計協變量（covariate），檢驗其是否有顯著效應

2. **分離分析**（當 read 數量足夠時）：
   - 分別對 `distance_forward.csv` 與 `distance_reverse.csv` 執行 HP 分析
   - 比較兩個 strand 的結果是否一致

**判斷標準**：

- 若某 strand 的 read 數 < 20，僅執行整合分析
- 若兩個 strand 的 read 數皆 ≥ 20，執行分離分析

---

### 3.3 統計分析實作架構

**C++ 側**：

- 輸出距離矩陣與 read metadata
- 呼叫 Python 腳本執行統計分析

**Python 側**：

```python
# scripts/statistics/run_statistical_tests.py

import pandas as pd
import numpy as np
from skbio.stats.distance import permanova
from scipy.stats import fisher_exact

def analyze_hp_association(
    distance_matrix_path,
    read_metadata_path,
    output_dir
):
    # 讀取資料
    dist_df = pd.read_csv(distance_matrix_path, index_col=0)
    metadata = pd.read_csv(read_metadata_path, sep='\t')
    
    # PERMANOVA
    permanova_result = permanova(
        dist_df.values,
        metadata['hp_tag'],
        permutations=999
    )
    
    # Clade Enrichment (需先執行聚類取得樹結構)
    tree = read_newick('tree.nwk')
    enrichment_results = []
    for clade in tree.major_clades:
        for hp_class in metadata['hp_tag'].unique():
            contingency_table = build_contingency_table(clade, hp_class, metadata)
            _, p_value = fisher_exact(contingency_table)
            enrichment_results.append({
                'clade_id': clade.id,
                'hp_class': hp_class,
                'p_value': p_value
            })
    
    # 輸出結果
    save_results(permanova_result, enrichment_results, output_dir)
```

---

## 4. 計算量與時間估算

### 4.1 現有效能基準

根據使用者提供的資訊：

- **距離矩陣計算**：40,000 個位點，64 執行緒，約 **40 秒**
- **平均每個位點**：40秒 / 40,000 = **1 ms**

---

### 4.2 聚類計算時間估算

**UPGMA 時間複雜度**：$O(N^2 \log N)$（使用優先佇列優化）

**假設**：

- 平均每個位點有 $N = 100$ reads
- 40,000 個位點

**單一位點聚類時間**：

- 樸素實作 $O(N^3)$：~5 ms
- 優化實作 $O(N^2 \log N)$：~1 ms

**總時間**：40,000 × 1 ms = **40 秒**

**加上距離矩陣計算**：40 + 40 = **80 秒**（約 1.3 分鐘）

---

### 4.3 Bootstrap 時間估算

**Bootstrap 一次迭代的成本**：

- 重抽樣 CpG：< 1 ms（記憶體操作）
- 重新計算距離矩陣：與原始計算相同，~1 ms
- 執行聚類：~1 ms
- **單次 Bootstrap 總計**：~2 ms / 位點

**執行 100 次 Bootstrap**：

- 單一位點：2 ms × 100 = **200 ms**
- 40,000 位點：40,000 × 200 ms = **8,000 秒**（約 2.2 小時）

**平行化加速（64 執行緒）**：

- 理想加速比：50x（考慮同步開銷）
- 實際時間：8,000 / 50 = **160 秒**（約 2.7 分鐘）

---

### 4.4 篩選策略（提升效率）

並非所有位點都需要執行 Bootstrap：

**篩選條件**：

1. Read 數量 ≥ 20（太少的 reads 聚類不穩定）
2. 有效距離對數 ≥ 50%（共同覆蓋的 CpG 位點足夠）
3. HP 標籤多樣性 ≥ 2 種（單一 HP 類別無意義）

**預估篩選後位點數**：~10,000 個（25%）

**優化後 Bootstrap 時間**：160 秒 × 25% = **40 秒**

---

### 4.5 Python 視覺化時間

**Heatmap 繪製**（`seaborn.clustermap`）：

- 單一位點：100 reads × 50 CpGs → ~500 ms
- 若篩選後輸出 10,000 個 heatmaps：10,000 × 0.5s = **5,000 秒**（約 1.4 小時）

**優化**：

- 僅為關鍵位點（如 PERMANOVA p < 0.05）繪製詳細 heatmap
- 預估關鍵位點 ~500 個：500 × 0.5s = **250 秒**（約 4 分鐘）

---

### 4.6 總計算時間匯總

| 步驟 | 全部位點 (40,000) | 篩選後 (10,000) |
|------|------------------|----------------|
| 距離矩陣計算 | 40 秒 | 40 秒 |
| 聚類建樹 | 40 秒 | 10 秒 |
| Bootstrap (100次, 64執行緒) | 160 秒 | 40 秒 |
| 統計分析 (Python) | ~60 秒 | ~15 秒 |
| Heatmap 繪製 (Python) | 5,000 秒 | 250 秒 |
| **總計** | **~1.5 小時** | **~6 分鐘** |

**🔴 建議**：採用篩選策略，僅對關鍵位點執行完整分析

---

## 5. 多重檢驗校正 (Multiple Testing Correction)

### 5.1 問題描述

當執行大量統計檢定時（如 40,000 個位點的 PERMANOVA），即使真實無顯著差異，也會因隨機性出現假陽性（Type I Error）。

**範例**：

- 執行 40,000 次檢定，顯著性水準 $\alpha = 0.05$
- 預期假陽性數量：40,000 × 0.05 = **2,000 個**
- 即使所有位點都無真實效應，仍會有 2,000 個 p < 0.05

---

### 5.2 是否需要校正？

**答案：是，強烈需要**

**原因**：

1. **多重比較問題**：執行上萬次檢定，必然產生大量假陽性
2. **生物學解釋**：未校正的結果會誤導研究結論
3. **發表要求**：期刊審稿人會要求多重檢驗校正

**若不校正的後果**：

- 將大量無意義的位點誤認為顯著
- 後續驗證實驗無法重現結果
- 浪費研究資源

---

### 5.3 校正方法選擇

#### 方法 1：Bonferroni 校正（最保守）

**原理**：將顯著性水準除以檢定次數

$$\alpha_{corrected} = \frac{\alpha}{n} = \frac{0.05}{40000} = 1.25 \times 10^{-6}$$

**優點**：簡單，控制 Family-Wise Error Rate (FWER)
**缺點**：過於保守，可能遺漏真陽性

**適用情境**：探索性研究，需要高信心的候選位點

---

#### 方法 2：Benjamini-Hochberg FDR 校正（推薦）⭐⭐⭐⭐⭐

**原理**：控制 False Discovery Rate（假發現率）

**步驟**：

1. 對所有 p-values 排序：$p_1 \leq p_2 \leq ... \leq p_n$
2. 找到最大的 $i$ 使得 $p_i \leq \frac{i}{n} \times q$（$q$ 為設定的 FDR，如 0.05）
3. 拒絕所有 $p_1, ..., p_i$

**Python 實作**：

```python
from statsmodels.stats.multitest import multipletests

p_values = [...]  # 所有位點的 p-values
reject, p_corrected, _, _ = multipletests(
    p_values, 
    alpha=0.05, 
    method='fdr_bh'
)
```

**優點**：

- 在保證 FDR < 5% 的前提下，最大化檢測力
- 適合 genomics 研究

**缺點**：

- 允許少量假陽性（但比例可控）

**適用情境**：本研究（推薦）

---

#### 方法 3：q-value 方法（進階）

**原理**：估計每個 p-value 對應的 FDR

**套件**：R 的 `qvalue` 或 Python 的 `qvalue`

**適用情境**：需要報告每個位點的 FDR

---

### 5.4 實作建議

**策略**：同時報告原始 p-value 與校正後 p-value

**輸出格式**（統計結果表）：

```
| locus_id       | n_reads | n_hp_classes | permanova_F | p_value | p_adj_fdr | significant |
|----------------|---------|--------------|-------------|---------|-----------|-------------|
| chr19_29283968 | 85      | 3            | 12.5        | 0.001   | 0.023     | True        |
| chr19_29284500 | 120     | 4            | 8.3         | 0.008   | 0.045     | True        |
| chr19_29285000 | 45      | 2            | 2.1         | 0.15    | 0.32      | False       |
```

**判斷標準**：

- `p_adj_fdr < 0.05`：顯著
- `0.05 ≤ p_adj_fdr < 0.1`：邊緣顯著（可報告但需謹慎解釋）
- `p_adj_fdr ≥ 0.1`：不顯著

---

### 5.5 對結果的影響

**預期顯著位點數量變化**：

| 情境 | 原始 p < 0.05 數量 | FDR 校正後顯著數量 | 比例 |
|------|------------------|------------------|------|
| 無真實效應 | 2,000 | ~50 | 2.5% |
| 10% 位點有效應 | 6,000 | ~3,500 | 58% |
| 30% 位點有效應 | 14,000 | ~11,000 | 79% |

**解釋**：

- 若真實存在效應，FDR 校正仍能保留大部分顯著位點
- 僅會過濾掉 p-value 接近 0.05 的邊緣案例

---

## 6. 完整工作流程與輸出

### 6.1 C++ 主程式流程

```cpp
// 偽代碼示意
for (auto& region : all_regions) {
    // Phase 1: 距離矩陣計算（已完成）
    DistanceMatrix dist_all = compute_distance_matrix(region);
    DistanceMatrix dist_fwd = compute_distance_matrix_strand(region, FORWARD);
    DistanceMatrix dist_rev = compute_distance_matrix_strand(region, REVERSE);
    
    // Phase 2: 層次聚類建樹
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree_all = clusterer.build_tree(dist_all);
    Tree tree_fwd = clusterer.build_tree(dist_fwd);
    Tree tree_rev = clusterer.build_tree(dist_rev);
    
    // Phase 3: Bootstrap 驗證（若 reads 數量足夠）
    if (region.reads.size() >= 20) {
        BootstrapAnalyzer bootstrap(config);
        BootstrapResult bs_result = bootstrap.run_bootstrap(
            region.methylation_matrix,
            dist_config,
            LinkageMethod::UPGMA
        );
        tree_all.annotate_support_values(bs_result.support_values);
    }
    
    // Phase 4: 輸出檔案
    write_newick(tree_all, output_dir + "/tree.nwk");
    write_newick(tree_fwd, output_dir + "/tree_forward.nwk");
    write_newick(tree_rev, output_dir + "/tree_reverse.nwk");
    write_distance_matrix(dist_all, output_dir + "/distance_matrix.csv");
    write_methylation_matrix(region, output_dir + "/methylation_matrix.csv");
    write_read_metadata(region, output_dir + "/read_metadata.tsv");
    
    // Phase 5: 呼叫 Python 繪製 Heatmap
    if (region.is_significant) {  // 可先用簡單標準篩選
        call_python_script(
            "scripts/visualization/plot_clustermap.py",
            {output_dir + "/methylation_matrix.csv",
             output_dir + "/tree.nwk",
             output_dir + "/read_metadata.tsv",
             output_dir + "/heatmap.png"}
        );
    }
}

// Phase 6: 呼叫 Python 統計分析（批次處理所有位點）
call_python_script(
    "scripts/statistics/run_permanova.py",
    {all_output_dirs, summary_output_path}
);
```

---

### 6.2 Python 視覺化腳本範例

```python
# scripts/visualization/plot_clustermap.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo

def plot_clustermap(
    methylation_matrix_path,
    tree_path,
    metadata_path,
    output_path
):
    # 讀取資料
    meth_df = pd.read_csv(methylation_matrix_path, index_col=0)
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    
    # 建立 annotation colors
    hp_colors = metadata['hp_tag'].map({
        'HP1': '#F5A3A3', # 淡紅色
        'HP2': '#A3C8F5', # 淡藍色
        'HP1-1': '#8B0000', # 深紅色
        'HP2-1': '#00008B', # 深藍色
        'HP3': '#D6A9E8', # 淡紫色
        'Unphased': '#999999', # 未分型
        'None': '#FFFFFF'
    })
    
    tumor_colors = metadata['tumor_normal'].map({
        'Tumor': '#F7DC6F', # 淡黃色
        'Normal': '#2CA02C' # 淡綠色
    })
    
    row_colors = pd.DataFrame({
        'HP': hp_colors,
        'Tumor': tumor_colors
    })
    
    # 繪製 Clustermap（使用預先計算的樹）
    # 注意：seaborn 會自動執行聚類，若要使用預先計算的樹需自行構建 linkage matrix
    # 此處簡化為讓 seaborn 自動聚類
    g = sns.clustermap(
        meth_df,
        row_colors=row_colors,
        method='average',  # UPGMA
        cmap='RdBu_r',
        cbar_kws={'label': 'Methylation Level'},
        figsize=(12, 10),
        dendrogram_ratio=0.15
    )
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
```

---

### 6.3 輸出檔案結構

```
output/
├── single_locus/
│   ├── chr19_29283968/
│   │   ├── methylation_matrix.csv      # Reads × CpGs 原始矩陣
│   │   ├── distance_matrix.csv          # 完整距離矩陣
│   │   ├── distance_forward.csv         # Forward strand 距離矩陣
│   │   ├── distance_reverse.csv         # Reverse strand 距離矩陣
│   │   ├── tree.nwk                     # Newick 格式演化樹
│   │   ├── tree_forward.nwk
│   │   ├── tree_reverse.nwk
│   │   ├── bootstrap_support.txt        # Bootstrap 支持度詳細報告
│   │   ├── read_metadata.tsv            # Read 標籤資訊
│   │   ├── heatmap.png                  # Cluster Heatmap
│   │   └── stats_summary.txt            # 該位點的統計摘要
│   └── ...
│
└── aggregate/
    ├── all_loci_statistics.tsv          # 所有位點的統計結果
    ├── permanova_results.tsv            # PERMANOVA 結果（含 FDR 校正）
    │   # Columns: locus_id, n_reads, permanova_F, p_value, p_adj_fdr, significant
    ├── enrichment_results.tsv           # Clade enrichment 結果
    ├── significant_loci_list.txt        # 顯著位點清單（FDR < 0.05）
    └── figures/
        ├── p_value_histogram.png        # P-value 分佈直方圖
        ├── volcano_plot.png             # Volcano plot (F-stat vs -log10(p))
        └── hp_tag_distribution.png      # HP 標籤在聚類中的分佈
```

---

## 7. 實作優先級與時程規劃

### Phase 1: 核心聚類功能（1-2 週）

- [x] 距離矩陣計算（已完成）
- [ ] 實作 UPGMA 聚類演算法
- [ ] Newick 格式輸出
- [ ] 單元測試驗證正確性
- [ ] 整合到 `RegionProcessor`

**驗收標準**：能為單一位點輸出正確的演化樹檔案

---

### Phase 2: Bootstrap 驗證（1 週）

- [ ] 實作 CpG 位點重抽樣
- [ ] 平行化 Bootstrap iterations
- [ ] 支持度計算演算法
- [ ] 輸出帶有 Bootstrap values 的 Newick

**驗收標準**：Bootstrap 支持度與 scipy 結果一致

---

### Phase 3: Python 視覺化整合（3-5 天）

- [ ] 開發 `plot_clustermap.py` 腳本
- [ ] 設計 annotation bar 配色
- [ ] C++ 呼叫 Python 機制
- [ ] 批次繪圖能力

**驗收標準**：輸出美觀的 heatmap，標籤清晰可讀

---

### Phase 4: 統計分析（1 週）

- [ ] 實作 PERMANOVA 分析腳本
- [ ] 實作 Clade Enrichment 分析
- [ ] FDR 多重檢驗校正
- [ ] 批次處理所有位點
- [ ] 產生綜合統計報告

**驗收標準**：產生完整的統計結果表，含 FDR 校正

---

### Phase 5: 優化與測試（3-5 天）

- [ ] 效能優化（篩選策略、記憶體管理）
- [ ] 邊界案例測試（少量 reads、單一 HP 類別）
- [ ] 大規模測試（40,000 位點完整執行）
- [ ] 文件撰寫

**驗收標準**：能在合理時間內（< 1 小時）完成 40,000 位點分析

---

## 8. 技術細節與注意事項

### 8.1 UPGMA 實作要點

1. **初始化**：每個 read 是一個 cluster，高度為 0
2. **距離矩陣更新**：合併兩個 clusters 後，新 cluster 到其他 cluster 的距離為算術平均
   $$d(AB, C) = \frac{|A| \cdot d(A, C) + |B| \cdot d(B, C)}{|A| + |B|}$$
3. **分支長度**：新節點高度 = 合併時的距離 / 2

**避免的陷阱**：

- 忘記更新 cluster 大小（影響加權平均）
- 分支長度計算錯誤（應為高度差，非絕對距離）

---

### 8.2 Newick 格式範例

```
((read_0:0.05,read_1:0.05)95:0.03,(read_2:0.08,read_3:0.08)80:0.02)100:0.0;
```

**解析**：

- `read_0`, `read_1` 是葉節點（實際 read IDs）
- `:0.05` 是分支長度
- `)95` 是 Bootstrap 支持度（百分比）
- 最外層 `)100:0.0` 是根節點

---

### 8.3 記憶體管理

**預估記憶體使用**（單一位點，N=100 reads, M=50 CpGs）：

- 甲基化矩陣：100 × 50 × 4 bytes (int32) = 20 KB
- 距離矩陣：100 × 100 × 8 bytes (double) = 80 KB
- 演化樹：< 10 KB
- **總計**：< 200 KB / 位點

**40,000 位點總記憶體**：< 8 GB（完全可接受）

**優化**：每個位點處理完即釋放記憶體，不需全部儲存

---

### 8.4 與現有程式碼整合

**修改檔案**：

1. `include/core/HierarchicalClustering.hpp`（新增）
2. `src/core/HierarchicalClustering.cpp`（新增）
3. `include/core/BootstrapAnalyzer.hpp`（新增）
4. `src/core/BootstrapAnalyzer.cpp`（新增）
5. `include/core/RegionProcessor.hpp`（新增成員變數與方法）
6. `src/core/RegionProcessor.cpp`（整合聚類流程）
7. `include/io/RegionWriter.hpp`（新增 Newick 輸出）
8. `src/io/RegionWriter.cpp`（實作 Newick 輸出）

**Config 新增參數**：

```cpp
struct Config {
    // ... 現有參數
    
    // Clustering Configuration
    bool compute_clustering = true;
    bool compute_bootstrap = true;
    int bootstrap_iterations = 100;
    int bootstrap_min_reads = 20;
    std::string linkage_method = "UPGMA";  // UPGMA, WARD, NJ
    
    // Visualization
    bool generate_heatmaps = true;
    bool heatmap_only_significant = true;
    std::string python_executable = "python3";
};
```

---

## 9. 總結

### 核心技術決策

1. ✅ **C++ 實作 UPGMA 聚類**（自行實作，200-300 行程式碼）
2. ✅ **C++ 實作 Bootstrap**（平行化，64 執行緒可達 50x 加速）
3. ✅ **Python 視覺化與統計**（呼叫成熟套件，精緻美觀）
4. ✅ **FDR 多重檢驗校正**（必須執行，使用 Benjamini-Hochberg 方法）

### 預期效能

- **完整分析時間**（40,000 位點）：
  - 無篩選：~1.5 小時
  - 篩選後（10,000 關鍵位點）：~6 分鐘
- **記憶體使用**：< 8 GB

### 多標籤分析策略

- **Level 1**：HP 標籤（PERMANOVA + Clade Enrichment）
- **Level 2**：HP × Tumor/Normal 組合
- **Level 3**：Strand 整合與分離分析

### 輸出成果

- 每個位點：演化樹（Newick）、距離矩陣、Heatmap、Bootstrap 支持度
- 綜合報告：統計表格（含 FDR 校正）、顯著位點清單、視覺化圖表

---

**下一步行動**：

1. 確認本方案符合研究需求
2. 開始實作 Phase 1（UPGMA 聚類）
3. 設計單元測試案例
