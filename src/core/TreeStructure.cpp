#include "core/TreeStructure.hpp"
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace InterSubMod {

// ============================================================================
// Tree Methods
// ============================================================================

std::string Tree::to_newick(
    bool include_bootstrap,
    bool include_branch_length,
    int precision
) const {
    if (!root_) {
        return ";";
    }
    
    // 特殊情況：單一葉節點
    if (root_->is_leaf()) {
        std::ostringstream oss;
        oss << root_->label;
        if (include_branch_length && root_->branch_length > 0) {
            oss << ":" << std::fixed << std::setprecision(precision) << root_->branch_length;
        }
        oss << ";";
        return oss.str();
    }
    
    return newick_recursive(root_, include_bootstrap, include_branch_length, precision) + ";";
}

std::string Tree::newick_recursive(
    const std::shared_ptr<TreeNode>& node,
    bool include_bootstrap,
    bool include_branch_length,
    int precision
) const {
    if (!node) {
        return "";
    }
    
    std::ostringstream oss;
    
    if (node->is_leaf()) {
        // 葉節點：返回標籤
        oss << node->label;
    } else {
        // 內部節點：遞迴處理左右子樹
        oss << "(";
        oss << newick_recursive(node->left, include_bootstrap, include_branch_length, precision);
        oss << ",";
        oss << newick_recursive(node->right, include_bootstrap, include_branch_length, precision);
        oss << ")";
        
        // 加入 Bootstrap 支持度（若有）
        if (include_bootstrap && node->bootstrap_support > 0.0) {
            oss << std::fixed << std::setprecision(0) << node->bootstrap_support;
        }
    }
    
    // 加入分支長度
    if (include_branch_length && node->branch_length > 0) {
        oss << ":" << std::fixed << std::setprecision(precision) << node->branch_length;
    }
    
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
    if (!node || node->is_leaf()) {
        return;
    }
    
    // Pre-order: 先加入當前節點
    nodes.push_back(node);
    
    // 遞迴處理子節點
    if (node->left) {
        collect_internal_nodes(node->left, nodes);
    }
    if (node->right) {
        collect_internal_nodes(node->right, nodes);
    }
}

std::vector<std::shared_ptr<TreeNode>> Tree::get_leaves() const {
    std::vector<std::shared_ptr<TreeNode>> leaves;
    if (root_) {
        collect_leaves(root_, leaves);
    }
    return leaves;
}

void Tree::collect_leaves(
    const std::shared_ptr<TreeNode>& node,
    std::vector<std::shared_ptr<TreeNode>>& leaves
) const {
    if (!node) {
        return;
    }
    
    if (node->is_leaf()) {
        leaves.push_back(node);
    } else {
        // In-order: 先左後右
        if (node->left) {
            collect_leaves(node->left, leaves);
        }
        if (node->right) {
            collect_leaves(node->right, leaves);
        }
    }
}

void Tree::annotate_bootstrap_support(const std::vector<double>& support_values) {
    auto internal_nodes = get_internal_nodes();
    
    if (internal_nodes.size() != support_values.size()) {
        throw std::runtime_error(
            "Bootstrap support values size mismatch: expected " + 
            std::to_string(internal_nodes.size()) + 
            ", got " + std::to_string(support_values.size())
        );
    }
    
    for (size_t i = 0; i < internal_nodes.size(); ++i) {
        internal_nodes[i]->bootstrap_support = support_values[i];
    }
}

std::vector<std::set<int>> Tree::get_all_clades() const {
    std::vector<std::set<int>> clades;
    if (root_) {
        collect_clades(root_, clades);
    }
    return clades;
}

void Tree::collect_clades(
    const std::shared_ptr<TreeNode>& node,
    std::vector<std::set<int>>& clades
) const {
    if (!node || node->is_leaf()) {
        return;
    }
    
    // 將此內部節點的 leaf_indices 轉為 set
    std::set<int> clade(node->leaf_indices.begin(), node->leaf_indices.end());
    clades.push_back(clade);
    
    // 遞迴處理子節點
    if (node->left) {
        collect_clades(node->left, clades);
    }
    if (node->right) {
        collect_clades(node->right, clades);
    }
}

Tree Tree::deep_copy() const {
    Tree new_tree;
    if (root_) {
        new_tree.root_ = deep_copy_node(root_);
    }
    new_tree.merge_records_ = merge_records_;
    return new_tree;
}

std::shared_ptr<TreeNode> Tree::deep_copy_node(
    const std::shared_ptr<TreeNode>& node
) const {
    if (!node) {
        return nullptr;
    }
    
    auto new_node = std::make_shared<TreeNode>();
    new_node->node_id = node->node_id;
    new_node->label = node->label;
    new_node->height = node->height;
    new_node->branch_length = node->branch_length;
    new_node->bootstrap_support = node->bootstrap_support;
    new_node->leaf_indices = node->leaf_indices;
    
    if (node->left) {
        new_node->left = deep_copy_node(node->left);
        new_node->left->parent = new_node;
    }
    if (node->right) {
        new_node->right = deep_copy_node(node->right);
        new_node->right->parent = new_node;
    }
    
    return new_node;
}

} // namespace InterSubMod

