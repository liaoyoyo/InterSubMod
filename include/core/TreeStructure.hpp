#pragma once

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace InterSubMod {

/**
 * @brief 樹的節點結構
 *
 * 用於表示階層聚類產生的二元樹。每個節點可能是：
 * - 葉節點：代表一條 read，有標籤但無子節點
 * - 內部節點：代表一次聚類合併，有左右子節點
 */
struct TreeNode {
    int node_id = -1;                ///< 節點 ID (葉節點 0..N-1, 內部節點 N..2N-2)
    std::string label;               ///< 節點標籤 (葉節點為 read_id 或名稱)
    double height = 0.0;             ///< 節點高度 (從葉節點的距離)
    double branch_length = 0.0;      ///< 到父節點的分支長度
    double bootstrap_support = 0.0;  ///< Bootstrap 支持度 (0-100)

    std::shared_ptr<TreeNode> left;   ///< 左子樹
    std::shared_ptr<TreeNode> right;  ///< 右子樹
    std::weak_ptr<TreeNode> parent;   ///< 父節點 (weak_ptr 避免循環引用)

    std::vector<int> leaf_indices;  ///< 此 clade 包含的葉節點索引 (排序後)

    /**
     * @brief 判斷是否為葉節點
     */
    bool is_leaf() const {
        return left == nullptr && right == nullptr;
    }

    /**
     * @brief 取得此節點下的所有葉節點數量
     */
    int num_leaves() const {
        return static_cast<int>(leaf_indices.size());
    }

    /**
     * @brief 建立葉節點
     */
    static std::shared_ptr<TreeNode> create_leaf(int index, const std::string& name) {
        auto node = std::make_shared<TreeNode>();
        node->node_id = index;
        node->label = name;
        node->height = 0.0;
        node->branch_length = 0.0;
        node->bootstrap_support = 100.0;  // 葉節點支持度預設 100
        node->leaf_indices = {index};
        return node;
    }

    /**
     * @brief 建立內部節點
     */
    static std::shared_ptr<TreeNode> create_internal(int node_id, std::shared_ptr<TreeNode> left_child,
                                                     std::shared_ptr<TreeNode> right_child, double merge_height) {
        auto node = std::make_shared<TreeNode>();
        node->node_id = node_id;
        node->label = "";  // 內部節點無標籤
        node->height = merge_height;
        node->bootstrap_support = 0.0;  // 初始為 0，後續由 Bootstrap 填入

        node->left = left_child;
        node->right = right_child;

        // 計算分支長度
        if (left_child) {
            left_child->branch_length = merge_height - left_child->height;
            left_child->parent = node;
        }
        if (right_child) {
            right_child->branch_length = merge_height - right_child->height;
            right_child->parent = node;
        }

        // 合併 leaf_indices
        node->leaf_indices = left_child ? left_child->leaf_indices : std::vector<int>{};
        if (right_child) {
            node->leaf_indices.insert(node->leaf_indices.end(), right_child->leaf_indices.begin(),
                                      right_child->leaf_indices.end());
        }
        std::sort(node->leaf_indices.begin(), node->leaf_indices.end());

        return node;
    }
};

/**
 * @brief 聚類合併記錄
 *
 * 記錄每次合併的詳細資訊，類似於 scipy 的 linkage matrix 格式
 */
struct MergeRecord {
    int cluster_i;       ///< 第一個被合併的 cluster ID
    int cluster_j;       ///< 第二個被合併的 cluster ID
    double distance;     ///< 合併時的距離
    int new_cluster_id;  ///< 合併後新 cluster 的 ID
    int size;            ///< 新 cluster 包含的葉節點數
};

/**
 * @brief 演化樹結構
 *
 * 封裝完整的樹狀結構，提供：
 * - Newick 格式輸出
 * - 內部節點遍歷
 * - Bootstrap 支持度標註
 * - Clade 比較功能
 */
class Tree {
public:
    Tree() = default;

    /**
     * @brief 設定樹根
     */
    void set_root(std::shared_ptr<TreeNode> root) {
        root_ = root;
    }

    /**
     * @brief 取得樹根
     */
    std::shared_ptr<TreeNode> get_root() const {
        return root_;
    }

    /**
     * @brief 判斷樹是否為空
     */
    bool empty() const {
        return root_ == nullptr;
    }

    /**
     * @brief 取得葉節點數量
     */
    int num_leaves() const {
        return root_ ? root_->num_leaves() : 0;
    }

    /**
     * @brief 取得內部節點數量
     */
    int num_internal_nodes() const {
        return root_ ? static_cast<int>(get_internal_nodes().size()) : 0;
    }

    /**
     * @brief 轉換為 Newick 格式字串
     *
     * @param include_bootstrap 是否包含 Bootstrap 支持度
     * @param include_branch_length 是否包含分支長度
     * @param precision 浮點數精度
     * @return Newick 格式字串
     */
    std::string to_newick(bool include_bootstrap = true, bool include_branch_length = true, int precision = 6) const;

    /**
     * @brief 取得所有內部節點 (用於 Bootstrap 分析)
     *
     * 以 pre-order 順序返回所有內部節點
     */
    std::vector<std::shared_ptr<TreeNode>> get_internal_nodes() const;

    /**
     * @brief 取得所有葉節點
     *
     * 以 in-order 順序返回所有葉節點
     */
    std::vector<std::shared_ptr<TreeNode>> get_leaves() const;

    /**
     * @brief 標註 Bootstrap 支持度
     *
     * @param support_values 每個內部節點的支持度 (與 get_internal_nodes() 順序對應)
     */
    void annotate_bootstrap_support(const std::vector<double>& support_values);

    /**
     * @brief 取得所有 clade 的 leaf_indices 集合
     *
     * 用於 Bootstrap 比較
     */
    std::vector<std::set<int>> get_all_clades() const;

    /**
     * @brief 設定合併記錄 (用於除錯和驗證)
     */
    void set_merge_records(const std::vector<MergeRecord>& records) {
        merge_records_ = records;
    }

    /**
     * @brief 取得合併記錄
     */
    const std::vector<MergeRecord>& get_merge_records() const {
        return merge_records_;
    }

    /**
     * @brief 深拷貝整棵樹
     */
    Tree deep_copy() const;

private:
    std::shared_ptr<TreeNode> root_;
    std::vector<MergeRecord> merge_records_;

    /**
     * @brief 遞迴生成 Newick 字串
     */
    std::string newick_recursive(const std::shared_ptr<TreeNode>& node, bool include_bootstrap,
                                 bool include_branch_length, int precision) const;

    /**
     * @brief 遞迴收集內部節點
     */
    void collect_internal_nodes(const std::shared_ptr<TreeNode>& node,
                                std::vector<std::shared_ptr<TreeNode>>& nodes) const;

    /**
     * @brief 遞迴收集葉節點
     */
    void collect_leaves(const std::shared_ptr<TreeNode>& node, std::vector<std::shared_ptr<TreeNode>>& leaves) const;

    /**
     * @brief 遞迴收集所有 clade
     */
    void collect_clades(const std::shared_ptr<TreeNode>& node, std::vector<std::set<int>>& clades) const;

    /**
     * @brief 遞迴深拷貝節點
     */
    std::shared_ptr<TreeNode> deep_copy_node(const std::shared_ptr<TreeNode>& node) const;
};

}  // namespace InterSubMod
