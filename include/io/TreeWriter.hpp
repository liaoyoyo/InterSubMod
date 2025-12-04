#pragma once

#include <string>
#include <fstream>
#include "core/TreeStructure.hpp"

namespace InterSubMod {

/**
 * @brief 樹狀結構輸出選項
 */
struct TreeOutputOptions {
    bool include_bootstrap = true;       ///< 是否包含 Bootstrap 支持度
    bool include_branch_length = true;   ///< 是否包含分支長度
    int precision = 6;                   ///< 浮點數精度
    double min_bootstrap_to_show = 0.0;  ///< 最小顯示的 Bootstrap 值 (低於此不顯示)
    bool quote_labels = false;           ///< 是否用引號包圍標籤
    bool replace_spaces = true;          ///< 是否替換標籤中的空格為底線
};

/**
 * @brief 樹狀結構檔案輸出器
 * 
 * 支援多種輸出格式：
 * - Newick: 標準的演化樹交換格式
 * - Extended Newick: 包含額外註解的 Newick
 * - 未來可擴充 Nexus, PhyloXML 等格式
 */
class TreeWriter {
public:
    TreeWriter() = default;
    explicit TreeWriter(const TreeOutputOptions& options) : options_(options) {}
    
    /**
     * @brief 將樹寫入 Newick 格式檔案
     * 
     * @param tree 演化樹
     * @param filepath 輸出檔案路徑
     * @return 是否成功
     */
    bool write_newick(const Tree& tree, const std::string& filepath) const;
    
    /**
     * @brief 將樹轉為 Newick 字串
     * 
     * @param tree 演化樹
     * @return Newick 格式字串
     */
    std::string to_newick_string(const Tree& tree) const;
    
    /**
     * @brief 寫入合併記錄 (類似 scipy linkage matrix)
     * 
     * 格式：cluster_i, cluster_j, distance, new_cluster_id, size
     * 
     * @param tree 演化樹（需包含 merge_records）
     * @param filepath 輸出檔案路徑
     * @return 是否成功
     */
    bool write_linkage_matrix(const Tree& tree, const std::string& filepath) const;
    
    /**
     * @brief 寫入樹的統計摘要
     * 
     * @param tree 演化樹
     * @param filepath 輸出檔案路徑
     * @return 是否成功
     */
    bool write_tree_stats(const Tree& tree, const std::string& filepath) const;
    
    /**
     * @brief 取得輸出選項
     */
    const TreeOutputOptions& options() const { return options_; }
    
    /**
     * @brief 設定輸出選項
     */
    void set_options(const TreeOutputOptions& options) { options_ = options; }
    
private:
    TreeOutputOptions options_;
    
    /**
     * @brief 處理標籤（替換特殊字符）
     */
    std::string process_label(const std::string& label) const;
};

} // namespace InterSubMod

