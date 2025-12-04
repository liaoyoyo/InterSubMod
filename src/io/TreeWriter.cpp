#include "io/TreeWriter.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <iostream>

namespace InterSubMod {

bool TreeWriter::write_newick(const Tree& tree, const std::string& filepath) const {
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filepath << std::endl;
        return false;
    }
    
    ofs << to_newick_string(tree);
    ofs.close();
    
    return true;
}

std::string TreeWriter::to_newick_string(const Tree& tree) const {
    if (tree.empty()) {
        return ";";
    }
    
    auto root = tree.get_root();
    
    // 使用遞迴生成 Newick 字串
    std::function<std::string(const std::shared_ptr<TreeNode>&)> build_newick;
    build_newick = [&](const std::shared_ptr<TreeNode>& node) -> std::string {
        if (!node) return "";
        
        std::ostringstream oss;
        
        if (node->is_leaf()) {
            // 葉節點：輸出標籤
            oss << process_label(node->label);
        } else {
            // 內部節點：遞迴處理子節點
            oss << "(";
            oss << build_newick(node->left);
            oss << ",";
            oss << build_newick(node->right);
            oss << ")";
            
            // Bootstrap 支持度
            if (options_.include_bootstrap && 
                node->bootstrap_support >= options_.min_bootstrap_to_show) {
                oss << std::fixed << std::setprecision(0) << node->bootstrap_support;
            }
        }
        
        // 分支長度
        if (options_.include_branch_length && node->branch_length > 0) {
            oss << ":" << std::fixed << std::setprecision(options_.precision) 
                << node->branch_length;
        }
        
        return oss.str();
    };
    
    return build_newick(root) + ";";
}

bool TreeWriter::write_linkage_matrix(const Tree& tree, const std::string& filepath) const {
    const auto& records = tree.get_merge_records();
    if (records.empty()) {
        std::cerr << "Warning: No merge records available" << std::endl;
        return false;
    }
    
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filepath << std::endl;
        return false;
    }
    
    // Header
    ofs << "cluster_i\tcluster_j\tdistance\tnew_cluster_id\tsize\n";
    
    // Data
    for (const auto& record : records) {
        ofs << record.cluster_i << "\t"
            << record.cluster_j << "\t"
            << std::fixed << std::setprecision(options_.precision) << record.distance << "\t"
            << record.new_cluster_id << "\t"
            << record.size << "\n";
    }
    
    ofs.close();
    return true;
}

bool TreeWriter::write_tree_stats(const Tree& tree, const std::string& filepath) const {
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filepath << std::endl;
        return false;
    }
    
    ofs << "Tree Statistics\n";
    ofs << "===============\n\n";
    
    if (tree.empty()) {
        ofs << "Tree is empty.\n";
        ofs.close();
        return true;
    }
    
    auto root = tree.get_root();
    int n_leaves = tree.num_leaves();
    int n_internal = tree.num_internal_nodes();
    
    ofs << "Number of leaves (taxa): " << n_leaves << "\n";
    ofs << "Number of internal nodes: " << n_internal << "\n";
    ofs << "Total nodes: " << (n_leaves + n_internal) << "\n";
    ofs << "Tree height (root): " << std::fixed << std::setprecision(6) 
        << root->height << "\n\n";
    
    // Branch length statistics
    auto internal_nodes = tree.get_internal_nodes();
    auto leaves = tree.get_leaves();
    
    std::vector<double> branch_lengths;
    for (const auto& node : leaves) {
        if (node->branch_length > 0) {
            branch_lengths.push_back(node->branch_length);
        }
    }
    for (const auto& node : internal_nodes) {
        if (node->branch_length > 0) {
            branch_lengths.push_back(node->branch_length);
        }
    }
    
    if (!branch_lengths.empty()) {
        std::sort(branch_lengths.begin(), branch_lengths.end());
        double sum = 0.0;
        for (double bl : branch_lengths) sum += bl;
        double mean = sum / branch_lengths.size();
        
        ofs << "Branch Length Statistics:\n";
        ofs << "  Min: " << branch_lengths.front() << "\n";
        ofs << "  Max: " << branch_lengths.back() << "\n";
        ofs << "  Mean: " << mean << "\n";
        ofs << "  Median: " << branch_lengths[branch_lengths.size()/2] << "\n";
        ofs << "  Total tree length: " << sum << "\n\n";
    }
    
    // Bootstrap statistics
    std::vector<double> bootstrap_values;
    for (const auto& node : internal_nodes) {
        if (node->bootstrap_support > 0) {
            bootstrap_values.push_back(node->bootstrap_support);
        }
    }
    
    if (!bootstrap_values.empty()) {
        std::sort(bootstrap_values.begin(), bootstrap_values.end());
        double sum = 0.0;
        for (double bs : bootstrap_values) sum += bs;
        double mean = sum / bootstrap_values.size();
        
        int high_support = std::count_if(bootstrap_values.begin(), bootstrap_values.end(),
                                          [](double bs) { return bs >= 95.0; });
        int medium_support = std::count_if(bootstrap_values.begin(), bootstrap_values.end(),
                                           [](double bs) { return bs >= 75.0 && bs < 95.0; });
        
        ofs << "Bootstrap Support Statistics:\n";
        ofs << "  Min: " << bootstrap_values.front() << "%\n";
        ofs << "  Max: " << bootstrap_values.back() << "%\n";
        ofs << "  Mean: " << std::fixed << std::setprecision(1) << mean << "%\n";
        ofs << "  Median: " << bootstrap_values[bootstrap_values.size()/2] << "%\n";
        ofs << "  High support (>=95%): " << high_support << " nodes\n";
        ofs << "  Medium support (75-95%): " << medium_support << " nodes\n";
    } else {
        ofs << "Bootstrap support: Not available\n";
    }
    
    ofs.close();
    return true;
}

std::string TreeWriter::process_label(const std::string& label) const {
    std::string result = label;
    
    // 替換空格
    if (options_.replace_spaces) {
        std::replace(result.begin(), result.end(), ' ', '_');
    }
    
    // 替換其他特殊字符（Newick 保留字符）
    for (char& c : result) {
        if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '[' || c == ']') {
            c = '_';
        }
    }
    
    // 加引號
    if (options_.quote_labels && !result.empty()) {
        result = "'" + result + "'";
    }
    
    return result;
}

} // namespace InterSubMod

