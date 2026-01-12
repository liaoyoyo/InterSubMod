# InterSubMod 聚類與演化分析開發實作指南

**版本**：1.0  
**日期**：2025-12-02  
**文件類型**：詳細開發指南

---

## 目錄

1. [系統架構設計](#1-系統架構設計)
2. [Phase 1: UPGMA 聚類實作](#2-phase-1-upgma-聚類實作)
3. [Phase 2: Bootstrap 驗證實作](#3-phase-2-bootstrap-驗證實作)
4. [Phase 3: Python 視覺化整合](#4-phase-3-python-視覺化整合)
5. [Phase 4: 統計分析實作](#5-phase-4-統計分析實作)
6. [Phase 5: 整合與優化](#6-phase-5-整合與優化)
7. [配置參數說明](#7-配置參數說明)
8. [常見問題與解決方案](#8-常見問題與解決方案)

---

## 1. 系統架構設計

### 1.1 新增的核心類別

```
include/core/
├── HierarchicalClustering.hpp    # 層次聚類與演化樹建構
├── BootstrapAnalyzer.hpp         # Bootstrap 驗證分析
└── TreeStructure.hpp             # 樹狀結構資料型別

src/core/
├── HierarchicalClustering.cpp
├── BootstrapAnalyzer.cpp
└── TreeStructure.cpp

include/io/
└── TreeWriter.hpp                # Newick 格式輸出

src/io/
└── TreeWriter.cpp
```

---

### 1.2 資料流程

```
RegionProcessor::process_region()
  ↓
[已完成] compute_distance_matrix() → DistanceMatrix
  ↓
[新增] HierarchicalClustering::build_tree() → Tree
  ↓
[新增] BootstrapAnalyzer::run_bootstrap() → BootstrapResult
  ↓
[新增] TreeWriter::write_newick() → tree.nwk
  ↓
[新增] call_python_visualization() → heatmap.png
  ↓
[新增] call_python_statistics() → statistics.tsv
```

---

## 2. Phase 1: UPGMA 聚類實作

### 2.1 類別定義

**檔案**：`include/core/HierarchicalClustering.hpp`

```cpp
#pragma once

#include <vector>
#include <string>
#include <memory>
#include "DistanceMatrix.hpp"
#include "TreeStructure.hpp"

namespace InterSubMod {

/**
 * @brief 層次聚類演算法實作 (UPGMA)
 */
class HierarchicalClustering {
public:
    enum class LinkageMethod {
        UPGMA,      ///< Unweighted Pair Group Method with Arithmetic Mean
        WARD,       ///< Ward's minimum variance method
        SINGLE,     ///< Single linkage (最小距離)
        COMPLETE    ///< Complete linkage (最大距離)
    };
    
    /**
     * @brief 建構函式
     * @param method 連結方法
     */
    explicit HierarchicalClustering(LinkageMethod method = LinkageMethod::UPGMA);
    
    /**
     * @brief 從距離矩陣建構演化樹
     * @param dist_matrix 距離矩陣 (N x N)
     * @param read_ids Read ID 列表 (用於葉節點標籤)
     * @return 演化樹結構
     */
    Tree build_tree(
        const DistanceMatrix& dist_matrix,
        const std::vector<std::string>& read_ids
    );
    
private:
    LinkageMethod method_;
    
    /**
     * @brief 計算兩個 cluster 之間的距離
     * @param dist 當前距離矩陣
     * @param cluster_a Cluster A 的元素索引
     * @param cluster_b Cluster B 的元素索引
     * @return 合併後的距離
     */
    double compute_cluster_distance(
        const Eigen::MatrixXd& dist,
        const std::vector<int>& cluster_a,
        const std::vector<int>& cluster_b
    ) const;
    
    /**
     * @brief UPGMA 核心演算法
     */
    Tree build_upgma(
        const Eigen::MatrixXd& dist_matrix,
        const std::vector<std::string>& read_ids
    );
};

} // namespace InterSubMod
```

---

### 2.2 Tree 資料結構

**檔案**：`include/core/TreeStructure.hpp`

```cpp
#pragma once

#include <vector>
#include <string>
#include <memory>

namespace InterSubMod {

/**
 * @brief 樹的節點結構
 */
struct TreeNode {
    int node_id;                    ///< 節點 ID (-1 表示葉節點)
    std::string label;              ///< 節點標籤 (葉節點為 read_id)
    double height;                  ///< 節點高度 (分支長度)
    double bootstrap_support;       ///< Bootstrap 支持度 (0-100)
    
    std::shared_ptr<TreeNode> left;  ///< 左子樹
    std::shared_ptr<TreeNode> right; ///< 右子樹
    
    std::vector<int> leaf_indices;   ///< 此 clade 包含的葉節點索引
    
    /**
     * @brief 判斷是否為葉節點
     */
    bool is_leaf() const {
        return left == nullptr && right == nullptr;
    }
};

/**
 * @brief 演化樹結構
 */
class Tree {
public:
    Tree() = default;
    
    /**
     * @brief 設定樹根
     */
    void set_root(std::shared_ptr<TreeNode> root) { root_ = root; }
    
    /**
     * @brief 取得樹根
     */
    std::shared_ptr<TreeNode> get_root() const { return root_; }
    
    /**
     * @brief 轉換為 Newick 格式字串
     * @param include_bootstrap 是否包含 Bootstrap 支持度
     * @return Newick 格式字串
     */
    std::string to_newick(bool include_bootstrap = true) const;
    
    /**
     * @brief 取得所有內部節點 (用於 Bootstrap 分析)
     */
    std::vector<std::shared_ptr<TreeNode>> get_internal_nodes() const;
    
    /**
     * @brief 標註 Bootstrap 支持度
     */
    void annotate_bootstrap_support(const std::vector<double>& support_values);
    
private:
    std::shared_ptr<TreeNode> root_;
    
    /**
     * @brief 遞迴生成 Newick 字串
     */
    std::string newick_recursive(
        const std::shared_ptr<TreeNode>& node,
        bool include_bootstrap
    ) const;
    
    /**
     * @brief 遞迴收集內部節點
     */
    void collect_internal_nodes(
        const std::shared_ptr<TreeNode>& node,
        std::vector<std::shared_ptr<TreeNode>>& nodes
    ) const;
};

} // namespace InterSubMod
```

---

### 2.3 UPGMA 核心演算法實作

**檔案**：`src/core/HierarchicalClustering.cpp`

```cpp
#include "HierarchicalClustering.hpp"
#include <limits>
#include <algorithm>
#include <queue>

namespace InterSubMod {

HierarchicalClustering::HierarchicalClustering(LinkageMethod method)
    : method_(method) {}

Tree HierarchicalClustering::build_tree(
    const DistanceMatrix& dist_matrix,
    const std::vector<std::string>& read_ids
) {
    switch (method_) {
        case LinkageMethod::UPGMA:
            return build_upgma(dist_matrix.get_matrix(), read_ids);
        // 其他方法的實作可後續擴充
        default:
            throw std::runtime_error("Unsupported linkage method");
    }
}

Tree HierarchicalClustering::build_upgma(
    const Eigen::MatrixXd& dist_matrix,
    const std::vector<std::string>& read_ids
) {
    int n = dist_matrix.rows();
    
    // 初始化：每個 read 是一個 cluster
    std::vector<std::shared_ptr<TreeNode>> clusters;
    std::vector<std::vector<int>> cluster_members;  // 每個 cluster 包含的原始元素索引
    std::vector<int> cluster_sizes;                 // 每個 cluster 的大小
    
    for (int i = 0; i < n; ++i) {
        auto leaf = std::make_shared<TreeNode>();
        leaf->node_id = -1;
        leaf->label = read_ids[i];
        leaf->height = 0.0;
        leaf->bootstrap_support = 100.0;
        leaf->leaf_indices = {i};
        
        clusters.push_back(leaf);
        cluster_members.push_back({i});
        cluster_sizes.push_back(1);
    }
    
    // 複製距離矩陣（會動態更新）
    Eigen::MatrixXd D = dist_matrix;
    
    int next_node_id = n;  // 內部節點 ID 從 n 開始
    
    // 迭代合併 clusters
    while (clusters.size() > 1) {
        // 找到距離最小的兩個 clusters
        double min_dist = std::numeric_limits<double>::max();
        int min_i = -1, min_j = -1;
        
        for (size_t i = 0; i < clusters.size(); ++i) {
            for (size_t j = i + 1; j < clusters.size(); ++j) {
                // 計算這兩個 cluster 之間的距離
                double dist = compute_cluster_distance(
                    D, cluster_members[i], cluster_members[j]
                );
                
                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        
        // 合併 cluster i 和 j
        auto new_node = std::make_shared<TreeNode>();
        new_node->node_id = next_node_id++;
        new_node->label = "";  // 內部節點無標籤
        new_node->height = min_dist / 2.0;  // UPGMA: 高度 = 距離 / 2
        new_node->bootstrap_support = 0.0;  // 初始為 0，後續由 Bootstrap 填入
        new_node->left = clusters[min_i];
        new_node->right = clusters[min_j];
        
        // 合併 leaf_indices
        new_node->leaf_indices = cluster_members[min_i];
        new_node->leaf_indices.insert(
            new_node->leaf_indices.end(),
            cluster_members[min_j].begin(),
            cluster_members[min_j].end()
        );
        
        // 更新 cluster 列表
        std::vector<int> new_members = cluster_members[min_i];
        new_members.insert(
            new_members.end(),
            cluster_members[min_j].begin(),
            cluster_members[min_j].end()
        );
        
        int new_size = cluster_sizes[min_i] + cluster_sizes[min_j];
        
        // 移除舊的 clusters（先移除索引較大的）
        if (min_i > min_j) {
            clusters.erase(clusters.begin() + min_i);
            cluster_members.erase(cluster_members.begin() + min_i);
            cluster_sizes.erase(cluster_sizes.begin() + min_i);
            
            clusters.erase(clusters.begin() + min_j);
            cluster_members.erase(cluster_members.begin() + min_j);
            cluster_sizes.erase(cluster_sizes.begin() + min_j);
        } else {
            clusters.erase(clusters.begin() + min_j);
            cluster_members.erase(cluster_members.begin() + min_j);
            cluster_sizes.erase(cluster_sizes.begin() + min_j);
            
            clusters.erase(clusters.begin() + min_i);
            cluster_members.erase(cluster_members.begin() + min_i);
            cluster_sizes.erase(cluster_sizes.begin() + min_i);
        }
        
        // 加入新 cluster
        clusters.push_back(new_node);
        cluster_members.push_back(new_members);
        cluster_sizes.push_back(new_size);
    }
    
    // 建構 Tree 物件
    Tree tree;
    tree.set_root(clusters[0]);
    return tree;
}

double HierarchicalClustering::compute_cluster_distance(
    const Eigen::MatrixXd& dist,
    const std::vector<int>& cluster_a,
    const std::vector<int>& cluster_b
) const {
    switch (method_) {
        case LinkageMethod::UPGMA: {
            // UPGMA: 算術平均
            double sum = 0.0;
            for (int i : cluster_a) {
                for (int j : cluster_b) {
                    sum += dist(i, j);
                }
            }
            return sum / (cluster_a.size() * cluster_b.size());
        }
        case LinkageMethod::SINGLE: {
            // Single linkage: 最小距離
            double min_dist = std::numeric_limits<double>::max();
            for (int i : cluster_a) {
                for (int j : cluster_b) {
                    min_dist = std::min(min_dist, dist(i, j));
                }
            }
            return min_dist;
        }
        case LinkageMethod::COMPLETE: {
            // Complete linkage: 最大距離
            double max_dist = 0.0;
            for (int i : cluster_a) {
                for (int j : cluster_b) {
                    max_dist = std::max(max_dist, dist(i, j));
                }
            }
            return max_dist;
        }
        default:
            throw std::runtime_error("Unsupported linkage method");
    }
}

} // namespace InterSubMod
```

---

### 2.4 Newick 格式輸出

**檔案**：`src/core/TreeStructure.cpp`

```cpp
#include "TreeStructure.hpp"
#include <sstream>
#include <iomanip>

namespace InterSubMod {

std::string Tree::to_newick(bool include_bootstrap) const {
    if (!root_) {
        return ";";
    }
    return newick_recursive(root_, include_bootstrap) + ";";
}

std::string Tree::newick_recursive(
    const std::shared_ptr<TreeNode>& node,
    bool include_bootstrap
) const {
    if (node->is_leaf()) {
        // 葉節點：返回標籤
        return node->label;
    }
    
    // 內部節點：遞迴處理左右子樹
    std::ostringstream oss;
    oss << "(";
    oss << newick_recursive(node->left, include_bootstrap);
    oss << ",";
    oss << newick_recursive(node->right, include_bootstrap);
    oss << ")";
    
    // 加入 Bootstrap 支持度（若有）
    if (include_bootstrap && node->bootstrap_support > 0.0) {
        oss << std::fixed << std::setprecision(0) << node->bootstrap_support;
    }
    
    // 加入分支長度
    oss << ":" << std::fixed << std::setprecision(6) << node->height;
    
    return oss.str();
}

std::vector<std::shared_ptr<TreeNode>> Tree::get_internal_nodes() const {
    std::vector<std::shared_ptr<TreeNode>> nodes;
    if (root_ && !root_->is_leaf()) {
        collect_internal_nodes(root_, nodes);
    }
    return nodes;
}

void Tree::collect_internal_nodes(
    const std::shared_ptr<TreeNode>& node,
    std::vector<std::shared_ptr<TreeNode>>& nodes
) const {
    if (node->is_leaf()) {
        return;
    }
    
    nodes.push_back(node);
    if (node->left) collect_internal_nodes(node->left, nodes);
    if (node->right) collect_internal_nodes(node->right, nodes);
}

void Tree::annotate_bootstrap_support(const std::vector<double>& support_values) {
    auto internal_nodes = get_internal_nodes();
    if (internal_nodes.size() != support_values.size()) {
        throw std::runtime_error("Bootstrap support values size mismatch");
    }
    
    for (size_t i = 0; i < internal_nodes.size(); ++i) {
        internal_nodes[i]->bootstrap_support = support_values[i];
    }
}

} // namespace InterSubMod
```

---

## 3. Phase 2: Bootstrap 驗證實作

### 3.1 類別定義

**檔案**：`include/core/BootstrapAnalyzer.hpp`

```cpp
#pragma once

#include <vector>
#include <random>
#include "MethylationMatrix.hpp"
#include "DistanceMatrix.hpp"
#include "HierarchicalClustering.hpp"

namespace InterSubMod {

/**
 * @brief Bootstrap 驗證分析器
 */
class BootstrapAnalyzer {
public:
    struct Config {
        int n_iterations = 100;       ///< Bootstrap 迭代次數
        int random_seed = 42;          ///< 隨機種子（可重現）
        int n_threads = 16;            ///< 平行執行緒數
        int min_common_cov = 3;        ///< 最小共同覆蓋 (用於距離計算)
    };
    
    struct BootstrapResult {
        Tree original_tree;                     ///< 原始演化樹
        std::vector<double> support_values;     ///< 每個內部節點的支持度
        int n_successful_iterations;            ///< 成功執行的 Bootstrap 次數
    };
    
    explicit BootstrapAnalyzer(Config config);
    
    /**
     * @brief 執行 Bootstrap 分析
     * @param original_matrix 原始甲基化矩陣
     * @param dist_config 距離計算配置
     * @param linkage 聚類方法
     * @param read_ids Read ID 列表
     * @return Bootstrap 結果
     */
    BootstrapResult run_bootstrap(
        const MethylationMatrix& original_matrix,
        const DistanceConfig& dist_config,
        HierarchicalClustering::LinkageMethod linkage,
        const std::vector<std::string>& read_ids
    );
    
private:
    Config config_;
    std::mt19937 rng_;
    
    /**
     * @brief 對 CpG 位點進行重抽樣
     */
    MethylationMatrix resample_columns(const MethylationMatrix& original);
    
    /**
     * @brief 檢查 clade 是否在樹中存在
     */
    bool is_clade_present(
        const Tree& tree,
        const std::vector<int>& clade_leaves
    ) const;
    
    /**
     * @brief 計算支持度
     */
    std::vector<double> calculate_support_values(
        const Tree& original_tree,
        const std::vector<Tree>& bootstrap_trees
    ) const;
};

} // namespace InterSubMod
```

---

### 3.2 Bootstrap 核心實作

**檔案**：`src/core/BootstrapAnalyzer.cpp`

```cpp
#include "BootstrapAnalyzer.hpp"
#include <omp.h>
#include <algorithm>

namespace InterSubMod {

BootstrapAnalyzer::BootstrapAnalyzer(Config config)
    : config_(config), rng_(config.random_seed) {}

BootstrapAnalyzer::BootstrapResult BootstrapAnalyzer::run_bootstrap(
    const MethylationMatrix& original_matrix,
    const DistanceConfig& dist_config,
    HierarchicalClustering::LinkageMethod linkage,
    const std::vector<std::string>& read_ids
) {
    // Step 1: 建構原始樹
    DistanceCalculator dist_calc(dist_config);
    DistanceMatrix orig_dist = dist_calc.calculate(original_matrix);
    
    HierarchicalClustering clusterer(linkage);
    Tree original_tree = clusterer.build_tree(orig_dist, read_ids);
    
    // Step 2: 平行執行 Bootstrap iterations
    std::vector<Tree> bootstrap_trees;
    bootstrap_trees.reserve(config_.n_iterations);
    
    #pragma omp parallel num_threads(config_.n_threads)
    {
        // 每個執行緒有自己的隨機數生成器
        std::mt19937 thread_rng(config_.random_seed + omp_get_thread_num());
        
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < config_.n_iterations; ++i) {
            try {
                // 重抽樣 CpG 位點
                MethylationMatrix resampled = resample_columns(original_matrix);
                
                // 重新計算距離矩陣
                DistanceMatrix bs_dist = dist_calc.calculate(resampled);
                
                // 建構 Bootstrap 樹
                Tree bs_tree = clusterer.build_tree(bs_dist, read_ids);
                
                #pragma omp critical
                {
                    bootstrap_trees.push_back(bs_tree);
                }
            } catch (const std::exception& e) {
                // 某些 Bootstrap iteration 可能失敗（如距離矩陣無效）
                // 記錄但不中斷
                #pragma omp critical
                {
                    std::cerr << "Bootstrap iteration " << i << " failed: " 
                              << e.what() << std::endl;
                }
            }
        }
    }
    
    // Step 3: 計算支持度
    std::vector<double> support_values = calculate_support_values(
        original_tree,
        bootstrap_trees
    );
    
    // Step 4: 標註支持度到原始樹
    original_tree.annotate_bootstrap_support(support_values);
    
    BootstrapResult result;
    result.original_tree = original_tree;
    result.support_values = support_values;
    result.n_successful_iterations = bootstrap_trees.size();
    
    return result;
}

MethylationMatrix BootstrapAnalyzer::resample_columns(
    const MethylationMatrix& original
) {
    int n_reads = original.binary_matrix.rows();
    int n_cpgs = original.binary_matrix.cols();
    
    // 有放回抽樣 CpG 位點索引
    std::uniform_int_distribution<int> dist(0, n_cpgs - 1);
    std::vector<int> sampled_indices;
    sampled_indices.reserve(n_cpgs);
    
    for (int i = 0; i < n_cpgs; ++i) {
        sampled_indices.push_back(dist(rng_));
    }
    
    // 建構新的矩陣
    MethylationMatrix resampled;
    resampled.binary_matrix.resize(n_reads, n_cpgs);
    resampled.raw_matrix.resize(n_reads, n_cpgs);
    
    for (int i = 0; i < n_reads; ++i) {
        for (int j = 0; j < n_cpgs; ++j) {
            int source_col = sampled_indices[j];
            resampled.binary_matrix(i, j) = original.binary_matrix(i, source_col);
            resampled.raw_matrix(i, j) = original.raw_matrix(i, source_col);
        }
    }
    
    return resampled;
}

std::vector<double> BootstrapAnalyzer::calculate_support_values(
    const Tree& original_tree,
    const std::vector<Tree>& bootstrap_trees
) const {
    auto original_internal = original_tree.get_internal_nodes();
    std::vector<double> support_values(original_internal.size(), 0.0);
    
    int n_bootstrap = bootstrap_trees.size();
    if (n_bootstrap == 0) {
        return support_values;  // 全為 0
    }
    
    // 對每個原始樹的內部節點
    for (size_t i = 0; i < original_internal.size(); ++i) {
        auto& clade_leaves = original_internal[i]->leaf_indices;
        
        // 排序以便比較
        std::vector<int> sorted_leaves = clade_leaves;
        std::sort(sorted_leaves.begin(), sorted_leaves.end());
        
        // 計算在多少棵 Bootstrap 樹中出現
        int count = 0;
        for (const auto& bs_tree : bootstrap_trees) {
            if (is_clade_present(bs_tree, sorted_leaves)) {
                ++count;
            }
        }
        
        support_values[i] = (100.0 * count) / n_bootstrap;
    }
    
    return support_values;
}

bool BootstrapAnalyzer::is_clade_present(
    const Tree& tree,
    const std::vector<int>& clade_leaves
) const {
    auto internal_nodes = tree.get_internal_nodes();
    
    for (const auto& node : internal_nodes) {
        std::vector<int> sorted_node_leaves = node->leaf_indices;
        std::sort(sorted_node_leaves.begin(), sorted_node_leaves.end());
        
        if (sorted_node_leaves == clade_leaves) {
            return true;
        }
    }
    
    return false;
}

} // namespace InterSubMod
```

---

---

## 4. Phase 3: Python 視覺化整合

### 4.1 兩種視覺化類型說明

在 Cluster 分析中，有兩種常見的 Heatmap 視覺化：

**選項 A：Distance-based Cluster Heatmap（推薦）**

- 繪製 **Read × Read 距離矩陣**  
- 直接視覺化聚類結果
- 左側與上方顯示 Dendrogram（演化樹）

**選項 B：Methylation Pattern Heatmap + Dendrogram**

- 繪製 **Read × CpG 甲基化矩陣**
- 展示聚類後的甲基化模式  
- 僅左側顯示 Dendrogram

**本指南提供選項 A 的實作**（更符合聚類分析的直接視覺化需求）

---

### 4.2 Python 腳本：Distance-based Cluster Heatmap

**檔案**：`scripts/visualization/plot_distance_heatmap.py`

```python
#!/usr/bin/env python3
"""
繪製 Read-Read 距離矩陣 Cluster Heatmap
這是真正的聚類視覺化：X軸與Y軸都是Reads，顏色表示距離
"""

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
import argparse

def load_distance_matrix(dist_csv_path):
    """載入 C++ 輸出的距離矩陣"""
    dist_df = pd.read_csv(dist_csv_path, index_col=0)
    return dist_df

def load_metadata(metadata_tsv_path):
    """載入 read metadata"""
    metadata = pd.read_csv(metadata_tsv_path, sep='\t', index_col=0)
    return metadata

def create_annotation_colors(metadata):
    """建立 annotation 顏色映射"""
    row_colors_df = pd.DataFrame(index=metadata.index)
    
    # HP 標籤
    if 'hp' in metadata.columns:
        hp_colors = {
            '0': '#CCCCCC',    # Unphased
            '1': '#E74C3C',    # HP1 - 紅
            '2': '#3498DB',    # HP2 - 藍
            '1-1': '#E74C3C',
            '1-2': '#9B59B6',
            '2-1': '#2ECC71',
            '2-2': '#3498DB',
        }
        row_colors_df['HP'] = metadata['hp'].astype(str).map(hp_colors)
    
    # Tumor/Normal
    if 'is_tumor' in metadata.columns:
        tumor_colors = {
            1: '#E74C3C',   # Tumor - 紅
            0: '#27AE60',   # Normal - 綠
            '1': '#E74C3C',
            '0': '#27AE60'
        }
        row_colors_df['Sample'] = metadata['is_tumor'].map(tumor_colors)
    
    # Strand
    if 'strand' in metadata.columns:
        strand_colors = {
            '+': '#FF6B6B',
            '-': '#4ECDC4',
            'Forward': '#FF6B6B',
            'Reverse': '#4ECDC4'
        }
        row_colors_df['Strand'] = metadata['strand'].map(strand_colors)
    
    # ALT/REF
    if 'alt_support' in metadata.columns:
        alt_colors = {
            'ALT': '#F39C12',
            'REF': '#8E44AD',
            'UNKNOWN': '#BDC3C7'
        }
        row_colors_df['Allele'] = metadata['alt_support'].map(alt_colors)
    
    return row_colors_df

def plot_distance_heatmap(
    dist_matrix_path,
    metadata_path,
    output_path,
    linkage_method='average',
    figsize=(12, 10),
    dpi=150
):
    """
    繪製距離矩陣 Cluster Heatmap
    
    Parameters:
    -----------
    dist_matrix_path : str
        距離矩陣 CSV 檔案路徑 (Reads x Reads)
    metadata_path : str
        Read metadata TSV 檔案路徑
    output_path : str
        輸出圖片路徑
    linkage_method : str
        聚類方法 ('average' = UPGMA, 'ward', 'single', 'complete')
    """
    # 1. 載入資料
    dist_df = load_distance_matrix(dist_matrix_path)
    metadata = load_metadata(metadata_path)
    
    # 2. 確保索引一致
    common_reads = dist_df.index.intersection(metadata.index)
    dist_df = dist_df.loc[common_reads, common_reads]
    metadata = metadata.loc[common_reads]
    
    # 3. 建立 annotation 顏色
    row_colors_df = create_annotation_colors(metadata)
    
    # 4. 計算聚類 linkage
    # 注意：如果 C++ 已輸出 Newick 樹，可從樹檔案讀取
    # 這裡為簡化起見，重新計算 linkage
    dist_condensed = squareform(dist_df.values, checks=False)
    Z = linkage(dist_condensed, method=linkage_method)
    
    # 5. 繪製 Clustermap
    g = sns.clustermap(
        dist_df,
        row_linkage=Z,         # 使用計算好的 linkage
        col_linkage=Z,         # 對稱矩陣，行列使用相同 linkage
        row_colors=row_colors_df,
        col_colors=row_colors_df,
        cmap='viridis',        # 距離用不同配色
        vmin=0, vmax=1,
        figsize=figsize,
        cbar_kws={'label': 'Distance (NHD)'},
        dendrogram_ratio=(0.15, 0.15),  # 左側與上方都顯示樹
        colors_ratio=0.02,
        linewidths=0,
        xticklabels=False,     # Read ID 太多不顯示
        yticklabels=False
    )
    
    # 6. 標題
    n_reads = len(dist_df)
    g.fig.suptitle(
        f'Read-Read Distance Heatmap (N={n_reads}, Linkage={linkage_method.upper()})',
        y=0.98,
        fontsize=14,
        fontweight='bold'
    )
    
    # 7. 軸標籤
    g.ax_heatmap.set_xlabel('Reads (clustered)', fontsize=12)
    g.ax_heatmap.set_ylabel('Reads (clustered)', fontsize=12)
    
    # 8. 儲存
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(g.fig)
    
    print(f"✓ Distance heatmap saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plot Read-Read distance matrix cluster heatmap'
    )
    parser.add_argument('--distance', required=True, 
                       help='Distance matrix CSV (from C++)')
    parser.add_argument('--metadata', required=True, 
                       help='Read metadata TSV')
    parser.add_argument('--output', required=True, 
                       help='Output image path')
    parser.add_argument('--linkage', default='average', 
                       choices=['average', 'ward', 'single', 'complete'],
                       help='Linkage method (default: average=UPGMA)')
    parser.add_argument('--figsize', default='12,10', 
                       help='Figure size as "width,height"')
    parser.add_argument('--dpi', type=int, default=150, 
                       help='Output resolution')
    
    args = parser.parse_args()
    
    figsize = tuple(map(int, args.figsize.split(',')))
    
    plot_distance_heatmap(
        args.distance,
        args.metadata,
        args.output,
        linkage_method=args.linkage,
        figsize=figsize,
        dpi=args.dpi
    )
```

**使用範例**：

```bash
python3 scripts/visualization/plot_distance_heatmap.py \
    --distance output/chr1_877772/chr1_876772_878772/distance/NHD/matrix.csv \
    --metadata output/chr1_877772/chr1_876772_878772/reads/reads.tsv \
    --output output/chr1_877772/chr1_876772_878772/plots/distance_heatmap.png \
    --linkage average \
    --dpi 200
```

---

### 4.3 C++ 整合聚類到 RegionProcessor

**重要**：當前 C++ 實作中，聚類功能已完成但**未整合到 RegionProcessor**。需要修改：

**檔案**：`src/core/RegionProcessor.cpp`

**位置**：在 `process_single_region()` 方法中，距離矩陣輸出後（約 L496）

```cpp
// ==== 現有程式碼：距離矩陣計算與輸出 (L398-496) ====

// 已完成：計算距離矩陣
DistanceMatrix all_dist = dist_calc.compute(meth_mat, read_list);

// 已完成：輸出距離矩陣
if (output_distance_matrix_) {
    std::string region_dir = writer.get_region_dir(...);
    writer.write_distance_matrices(
        region_dir, all_dist, forward_dist, reverse_dist,
        metric, output_strand_distance_matrices_
    );
}

// ==== 新增：執行聚類建樹 ====
if (config_.compute_clustering && result.num_reads >= config_.clustering_min_reads) {
    #include "core/HierarchicalClustering.hpp"
    #include "io/TreeWriter.hpp"
    
    // 準備聚類配置
    HierarchicalClustering::ClusteringConfig cluster_config;
    cluster_config.method = config_.linkage_method;  // UPGMA, WARD, etc.
    cluster_config.min_branch_length = 1e-8;
    
    HierarchicalClustering clusterer(cluster_config);
    
    // 準備 Read IDs
    std::vector<std::string> read_ids;
    for (const auto& r : read_list) {
        read_ids.push_back(r.read_id);
    }
    
    // 建構演化樹（所有 reads）
    Tree tree_all = clusterer.build_tree(all_dist, read_ids);
    
    // 輸出 Newick 格式
    TreeWriter tree_writer;
    std::string tree_path = region_dir + "/tree.nwk";
    tree_writer.write_newick(tree_all, tree_path);
    
    if (log_level_ >= LogLevel::LOG_DEBUG) {
        #pragma omp critical
        {
            std::cout << "  ✓ Tree built: " << tree_all.num_leaves() 
                      << " leaves, saved to " << tree_path << std::endl;
        }
    }
    
    // Strand-specific 樹（若啟用）
    if (output_strand_distance_matrices_) {
        // Forward strand
        if (forward_dist.size() >= 2) {
            std::vector<std::string> fwd_ids;
            for (const auto& r : read_list) {
                if (r.strand == Strand::FORWARD) {
                    fwd_ids.push_back(r.read_id);
                }
            }
            if (fwd_ids.size() >= 2) {
                Tree tree_fwd = clusterer.build_tree(forward_dist, fwd_ids);
                tree_writer.write_newick(tree_fwd, region_dir + "/tree_forward.nwk");
            }
        }
        
        // Reverse strand（同理）
        // ... 
    }
}
```

**必要的 Config 參數修改**：

**檔案**：`include/core/Config.hpp`

```cpp
struct Config {
    // ... 現有參數 ...
    
    // === 新增：Clustering Configuration ===
    bool compute_clustering = true;               ///< 是否執行聚類建樹
    bool output_tree_files = true;                ///< 是否輸出 Newick 樹檔案
    LinkageMethod linkage_method = LinkageMethod::UPGMA;  ///< 連結方法
    int clustering_min_reads = 10;                ///< 最小 Read 數閾值
};
```

---

### 4.4 C++ 呼叫 Python 繪圖（批次模式）

**不建議**：在 `RegionProcessor` 中每個位點立即呼叫 Python（效率低）

**推薦**：使用 Shell 腳本在 C++ 完成後批次呼叫

**檔案**：`scripts/run_full_vcf_test.sh`（已存在，需修改）

```bash
# Step 1: C++ 處理所有位點
./build/bin/inter_sub_mod \
    --vcf data.vcf \
    --output-dir output/ \
    --compute-clustering \    # ← 新增參數
    --threads 64

# Step 2: 批次繪製 Heatmap
python3 tools/plot_cluster_heatmap.py \
    --output-dir output/ \
    --threads 16 \
    --metric NHD \
    --linkage average \
    --plot-type distance      # ← 指定繪製距離矩陣
```

**若必須在 C++ 中呼叫 Python**（不推薦）：

``cpp
void RegionProcessor::call_python_visualization(
    const std::string& region_dir,
    const std::string& region_id
) {
    // 檢查資料檔案是否存在
    std::string dist_csv = region_dir + "/distance/NHD/matrix.csv";
    std::string metadata_tsv = region_dir + "/reads/reads.tsv";

    if (!std::filesystem::exists(dist_csv) || !std::filesystem::exists(metadata_tsv)) {
        logger_.warning("Missing data files for visualization: " + region_id);
        return;
    }
    
    // 建構 Python 命令
    std::string script_path = "scripts/visualization/plot_distance_heatmap.py";
    std::string output_png = region_dir + "/plots/distance_heatmap.png";
    
    // 建立 plots 目錄
    std::filesystem::create_directories(region_dir + "/plots");
    
    std::ostringstream cmd;
    cmd << config_.python_executable << " " << script_path
        << " --distance " << dist_csv
        << " --metadata " << metadata_tsv
        << " --output " << output_png
        << " --linkage average"
        << " --dpi 150";
    
    int ret = std::system(cmd.str().c_str());
    if (ret != 0) {
        logger_.warning("Python visualization failed for: " + region_id);
    }
}

```

---

## 5. Phase 4: 統計分析實作

### 5.1 PERMANOVA 分析腳本

**檔案**：`scripts/statistics/run_permanova.py`

```python
#!/usr/bin/env python3
"""
批次執行 PERMANOVA 統計分析
"""

import pandas as pd
import numpy as np
from skbio.stats.distance import permanova, DistanceMatrix
from statsmodels.stats.multitest import multipletests
import argparse
import glob
import os

def run_permanova_single_locus(
    distance_matrix_path,
    metadata_path,
    group_column='hp_tag',
    permutations=999
):
    """
    對單一位點執行 PERMANOVA
    
    Returns:
    --------
    dict with keys: F_statistic, p_value, n_reads, n_groups
    """
    # 讀取距離矩陣
    dist_df = pd.read_csv(distance_matrix_path, index_col=0)
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    
    # 對齊索引
    common_reads = dist_df.index.intersection(metadata.index)
    if len(common_reads) < 4:
        return None  # 樣本數太少
    
    dist_df = dist_df.loc[common_reads, common_reads]
    metadata = metadata.loc[common_reads]
    
    # 提取分組標籤
    if group_column not in metadata.columns:
        return None
    
    groups = metadata[group_column]
    n_unique_groups = groups.nunique()
    
    if n_unique_groups < 2:
        return None  # 只有一個分組，無法分析
    
    # 建構 DistanceMatrix 物件
    dm = DistanceMatrix(dist_df.values, ids=dist_df.index.tolist())
    
    # 執行 PERMANOVA
    try:
        result = permanova(dm, groups, permutations=permutations)
        
        return {
            'F_statistic': result['test statistic'],
            'p_value': result['p-value'],
            'n_reads': len(common_reads),
            'n_groups': n_unique_groups
        }
    except Exception as e:
        print(f"PERMANOVA failed: {e}")
        return None

def batch_permanova(
    input_dirs,
    output_path,
    group_column='hp_tag'
):
    """
    批次處理多個位點的 PERMANOVA 分析
    """
    results = []
    
    for locus_dir in input_dirs:
        locus_id = os.path.basename(locus_dir)
        
        dist_path = os.path.join(locus_dir, 'distance_matrix.csv')
        metadata_path = os.path.join(locus_dir, 'read_metadata.tsv')
        
        if not os.path.exists(dist_path) or not os.path.exists(metadata_path):
            continue
        
        result = run_permanova_single_locus(
            dist_path,
            metadata_path,
            group_column=group_column
        )
        
        if result:
            results.append({
                'locus_id': locus_id,
                **result
            })
    
    # 彙整結果
    df = pd.DataFrame(results)
    
    # FDR 校正
    if len(df) > 0:
        _, p_adj, _, _ = multipletests(
            df['p_value'],
            alpha=0.05,
            method='fdr_bh'
        )
        df['p_adj_fdr'] = p_adj
        df['significant'] = df['p_adj_fdr'] < 0.05
    
    # 排序並儲存
    df = df.sort_values('p_adj_fdr')
    df.to_csv(output_path, sep='\t', index=False)
    
    print(f"PERMANOVA results saved to: {output_path}")
    print(f"Significant loci (FDR < 0.05): {df['significant'].sum()} / {len(df)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Batch PERMANOVA analysis')
    parser.add_argument('--input-pattern', required=True, 
                        help='Pattern for locus directories (e.g., output/single_locus/*)')
    parser.add_argument('--output', required=True, 
                        help='Output TSV file path')
    parser.add_argument('--group-column', default='hp_tag',
                        help='Column name for grouping (default: hp_tag)')
    
    args = parser.parse_args()
    
    locus_dirs = glob.glob(args.input_pattern)
    
    batch_permanova(
        locus_dirs,
        args.output,
        args.group_column
    )
```

---

## 6. Phase 5: 整合與優化

### 6.1 RegionProcessor 整合

**檔案**：`src/core/RegionProcessor.cpp` (修改 `process_region` 方法)

```cpp
void RegionProcessor::process_region(
    const SomaticSnv& snv,
    int region_start,
    int region_end
) {
    // [現有代碼] 讀取 reads、解析甲基化、建立矩陣...
    
    // === 新增：聚類分析 ===
    if (config_.compute_clustering && methylation_matrix.binary_matrix.rows() >= config_.clustering_min_reads) {
        
        // 1. 建構演化樹
        HierarchicalClustering clusterer(config_.linkage_method);
        Tree tree = clusterer.build_tree(distance_matrix, read_ids);
        
        // 2. Bootstrap 驗證（若啟用）
        if (config_.compute_bootstrap) {
            BootstrapAnalyzer::Config bs_config;
            bs_config.n_iterations = config_.bootstrap_iterations;
            bs_config.n_threads = config_.threads;
            
            BootstrapAnalyzer bs_analyzer(bs_config);
            auto bs_result = bs_analyzer.run_bootstrap(
                methylation_matrix,
                distance_config,
                config_.linkage_method,
                read_ids
            );
            
            tree = bs_result.original_tree;  // 帶有 Bootstrap 支持度的樹
        }
        
        // 3. 輸出演化樹
        TreeWriter tree_writer;
        tree_writer.write_newick(tree, output_dir + "/tree.nwk");
        tree_writer.write_newick(tree_fwd, output_dir + "/tree_forward.nwk");
        tree_writer.write_newick(tree_rev, output_dir + "/tree_reverse.nwk");
        
        // 4. 呼叫 Python 視覺化
        if (config_.generate_heatmaps) {
            bool is_significant = check_if_significant(snv);  // 可事先篩選
            
            if (!config_.heatmap_only_significant || is_significant) {
                call_python_visualization(output_dir, region_id);
            }
        }
    }
}
```

---

### 6.2 篩選策略實作

```cpp
bool RegionProcessor::should_run_full_analysis(
    const MethylationMatrix& matrix,
    const std::vector<ReadInfo>& reads
) const {
    int n_reads = reads.size();
    
    // 條件 1: Read 數量足夠
    if (n_reads < config_.clustering_min_reads) {
        return false;
    }
    
    // 條件 2: HP 標籤多樣性
    std::set<std::string> unique_hp;
    for (const auto& read : reads) {
        unique_hp.insert(read.hp_tag);
    }
    if (unique_hp.size() < 2) {
        return false;
    }
    
    // 條件 3: 有效距離對比例
    // (可從 distance_matrix 的統計資訊計算)
    
    return true;
}
```

---

## 7. 配置參數說明

### 7.1 Config 新增欄位

**檔案**：`include/core/Config.hpp`

```cpp
struct Config {
    // ... 現有參數 ...
    
    // === Clustering 配置 ===
    bool compute_clustering = true;               ///< 是否執行聚類
    int clustering_min_reads = 20;                ///< 最小 Read 數（低於此不執行）
    std::string linkage_method = "UPGMA";         ///< UPGMA, WARD, SINGLE, COMPLETE
    
    // === Bootstrap 配置 ===
    bool compute_bootstrap = true;                ///< 是否執行 Bootstrap
    int bootstrap_iterations = 100;               ///< Bootstrap 迭代次數
    int bootstrap_random_seed = 42;               ///< 隨機種子
    
    // === Visualization 配置 ===
    bool generate_heatmaps = true;                ///< 是否生成 Heatmap
    bool heatmap_only_significant = true;         ///< 僅為顯著位點生成
    std::string python_executable = "python3";    ///< Python 直譯器路徑
    
    // === Statistics 配置 ===
    bool run_statistics = true;                   ///< 是否執行統計分析
    std::string stats_group_column = "hp_tag";    ///< PERMANOVA 分組欄位
    int permanova_permutations = 999;             ///< PERMANOVA 排列次數
    double fdr_alpha = 0.05;                      ///< FDR 顯著性水準
};
```

---

### 7.2 命令列參數

| 參數 | 預設值 | 說明 |
|------|--------|------|
| `--compute-clustering` | true | 執行聚類分析 |
| `--clustering-min-reads` | 20 | 最小 Read 數閾值 |
| `--linkage-method` | UPGMA | 連結方法 |
| `--compute-bootstrap` | true | 執行 Bootstrap |
| `--bootstrap-iterations` | 100 | Bootstrap 次數 |
| `--generate-heatmaps` | true | 生成 Heatmap |
| `--python-exec` | python3 | Python 路徑 |

---

## 8. 常見問題與解決方案

### Q1: UPGMA 結果與 scipy 不一致？

**檢查項目**：

1. 距離矩陣是否完全一致（浮點數精度）
2. 連結方法是否確實為 UPGMA（arithmetic mean）
3. 分支長度計算是否正確（height = distance / 2）

**除錯方法**：

```cpp
// 輸出中間步驟的距離矩陣與合併順序
std::cout << "Merging cluster " << min_i << " and " << min_j 
          << " at distance " << min_dist << std::endl;
```

### Q2: Bootstrap 速度太慢？

**優化方案**：

1. 增加執行緒數：`--threads 64`
2. 降低 Bootstrap 次數（50 次通常已足夠）
3. 啟用篩選策略，僅處理關鍵位點

### Q3: Python 腳本找不到模組？

**解決方法**：

```bash
# 明確指定 conda 環境
--python-exec /path/to/conda/envs/myenv/bin/python3

# 或在腳本開頭加入
import sys
sys.path.insert(0, '/path/to/site-packages')
```

### Q4: 記憶體不足？

**方案**：

1. 處理超高深度位點時，限制平行度：

   ```cpp
   if (n_reads > 500) {
       omp_set_num_threads(1);  // 序列執行
   }
   ```

2. 立即釋放不需要的資料結構

---

**文件狀態**：完成  
**配套文件**：[測試驗證計劃](#)、[研究主軸說明](#)
