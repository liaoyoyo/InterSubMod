#pragma once

#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <Eigen/Dense>
#include "DistanceMatrix.hpp"
#include "TreeStructure.hpp"

namespace InterSubMod {

/**
 * @brief 階層聚類連結方法
 * 
 * 定義不同的聚類合併策略：
 * - UPGMA: 算術平均連結，假設分子鐘（適合甲基化分析）
 * - WARD: 最小變異法，傾向產生大小平衡的群
 * - SINGLE: 單連結法（最近鄰），傾向產生鏈狀結構
 * - COMPLETE: 完全連結法（最遠鄰），傾向產生緊密圓形群
 */
enum class LinkageMethod {
    UPGMA,      ///< Unweighted Pair Group Method with Arithmetic Mean
    WARD,       ///< Ward's minimum variance method
    SINGLE,     ///< Single linkage (最小距離)
    COMPLETE    ///< Complete linkage (最大距離)
};

/**
 * @brief 階層聚類配置
 */
struct ClusteringConfig {
    LinkageMethod method = LinkageMethod::UPGMA;  ///< 聚類方法
    bool optimal_leaf_ordering = false;           ///< 是否優化葉節點排序
    double min_branch_length = 1e-6;              ///< 最小分支長度（避免零長度分支）
};

/**
 * @brief 層次聚類演算法實作
 * 
 * 支援多種連結方法，從距離矩陣建構二元聚類樹。
 * 主要用於甲基化模式聚類分析。
 * 
 * 使用方式：
 * @code
 * HierarchicalClustering clusterer(LinkageMethod::UPGMA);
 * Tree tree = clusterer.build_tree(dist_matrix, read_ids);
 * std::string newick = tree.to_newick();
 * @endcode
 */
class HierarchicalClustering {
public:
    /**
     * @brief 建構函式（簡化版）
     * @param method 連結方法
     */
    explicit HierarchicalClustering(LinkageMethod method = LinkageMethod::UPGMA);
    
    /**
     * @brief 建構函式（完整版）
     * @param config 聚類配置
     */
    explicit HierarchicalClustering(const ClusteringConfig& config);
    
    /**
     * @brief 從 DistanceMatrix 物件建構演化樹
     * 
     * @param dist_matrix 距離矩陣物件
     * @param read_names Read 名稱列表（用於葉節點標籤）
     * @return 演化樹結構
     */
    Tree build_tree(
        const DistanceMatrix& dist_matrix,
        const std::vector<std::string>& read_names
    );
    
    /**
     * @brief 從 Eigen 矩陣建構演化樹
     * 
     * @param dist_matrix N x N 距離矩陣（必須對稱）
     * @param read_names Read 名稱列表
     * @return 演化樹結構
     */
    Tree build_tree(
        const Eigen::MatrixXd& dist_matrix,
        const std::vector<std::string>& read_names
    );
    
    /**
     * @brief 取得目前使用的連結方法
     */
    LinkageMethod get_method() const { return config_.method; }
    
    /**
     * @brief 設定連結方法
     */
    void set_method(LinkageMethod method) { config_.method = method; }
    
    /**
     * @brief 取得配置
     */
    const ClusteringConfig& config() const { return config_; }
    
    /**
     * @brief 將 LinkageMethod 轉為字串
     */
    static std::string method_to_string(LinkageMethod method);
    
    /**
     * @brief 從字串解析 LinkageMethod
     */
    static LinkageMethod string_to_method(const std::string& str);
    
private:
    ClusteringConfig config_;
    
    /**
     * @brief UPGMA 演算法核心實作
     * 
     * 時間複雜度: O(N^3)，空間複雜度: O(N^2)
     */
    Tree build_upgma(
        const Eigen::MatrixXd& dist_matrix,
        const std::vector<std::string>& read_names
    );
    
    /**
     * @brief Ward 演算法核心實作
     */
    Tree build_ward(
        const Eigen::MatrixXd& dist_matrix,
        const std::vector<std::string>& read_names
    );
    
    /**
     * @brief Single Linkage 演算法核心實作
     */
    Tree build_single(
        const Eigen::MatrixXd& dist_matrix,
        const std::vector<std::string>& read_names
    );
    
    /**
     * @brief Complete Linkage 演算法核心實作
     */
    Tree build_complete(
        const Eigen::MatrixXd& dist_matrix,
        const std::vector<std::string>& read_names
    );
    
    /**
     * @brief 通用聚類演算法框架
     * 
     * 使用策略模式，根據不同的距離更新函式執行聚類
     * 
     * @param dist_matrix 距離矩陣
     * @param read_names 節點名稱
     * @param update_distance 距離更新函式
     * @param distance_to_height 距離轉高度函式
     * @return 演化樹
     */
    Tree build_generic(
        const Eigen::MatrixXd& dist_matrix,
        const std::vector<std::string>& read_names,
        std::function<double(
            const Eigen::MatrixXd& D,
            const std::vector<int>& cluster_a,
            const std::vector<int>& cluster_b,
            int size_a, int size_b
        )> compute_cluster_distance,
        std::function<double(double distance)> distance_to_height = [](double d) { return d / 2.0; }
    );
    
    /**
     * @brief 找到距離矩陣中最小的非對角元素
     * 
     * @param D 當前距離矩陣
     * @param active_clusters 活躍的 cluster 索引
     * @param[out] min_i 最小距離的第一個索引
     * @param[out] min_j 最小距離的第二個索引
     * @return 最小距離值
     */
    double find_minimum_distance(
        const Eigen::MatrixXd& D,
        const std::vector<bool>& active,
        int& min_i,
        int& min_j
    );
};

/**
 * @brief 聚類結果切割器
 * 
 * 根據指定的切割條件（距離閾值或群數），從樹狀結構產生 cluster labels
 */
class TreeCutter {
public:
    /**
     * @brief 以距離閾值切割樹
     * 
     * @param tree 演化樹
     * @param distance_threshold 距離閾值
     * @return 每個葉節點的 cluster label (0-based)
     */
    static std::vector<int> cut_by_distance(const Tree& tree, double distance_threshold);
    
    /**
     * @brief 以群數切割樹
     * 
     * @param tree 演化樹
     * @param num_clusters 目標群數
     * @return 每個葉節點的 cluster label (0-based)
     */
    static std::vector<int> cut_by_num_clusters(const Tree& tree, int num_clusters);
    
    /**
     * @brief 計算每個切割下的 silhouette score，選擇最佳群數
     * 
     * @param tree 演化樹
     * @param dist_matrix 原始距離矩陣
     * @param min_k 最小群數
     * @param max_k 最大群數
     * @return pair<最佳群數, 對應的 labels>
     */
    static std::pair<int, std::vector<int>> find_optimal_clusters(
        const Tree& tree,
        const Eigen::MatrixXd& dist_matrix,
        int min_k = 2,
        int max_k = 10
    );
};

} // namespace InterSubMod

