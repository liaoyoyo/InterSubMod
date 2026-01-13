#pragma once

#include <fstream>
#include <string>

#include "core/TreeStructure.hpp"

namespace InterSubMod {

/**
 * @brief Tree structure output options
 */
struct TreeOutputOptions {
    bool include_bootstrap = true;       ///< Whether to include Bootstrap support
    bool include_branch_length = true;   ///< Whether to include branch lengths
    int precision = 6;                   ///< Floating-point precision
    double min_bootstrap_to_show = 0.0;  ///< Minimum Bootstrap value to show (below this not shown)
    bool quote_labels = false;           ///< Whether to quote labels
    bool replace_spaces = true;          ///< Whether to replace spaces in labels with underscores
};

/**
 * @brief Tree structure file writer
 *
 * Supports multiple output formats:
 * - Newick: Standard phylogenetic tree exchange format
 * - Extended Newick: Newick with additional annotations
 * - Future extensible to Nexus, PhyloXML, etc.
 */
class TreeWriter {
public:
    TreeWriter() = default;
    explicit TreeWriter(const TreeOutputOptions& options) : options_(options) {
    }

    /**
     * @brief Write tree to Newick format file
     *
     * @param tree Phylogenetic tree
     * @param filepath Output file path
     * @return Whether successful
     */
    bool write_newick(const Tree& tree, const std::string& filepath) const;

    /**
     * @brief Convert tree to Newick string
     *
     * @param tree Phylogenetic tree
     * @return Newick format string
     */
    std::string to_newick_string(const Tree& tree) const;

    /**
     * @brief Write merge records (similar to scipy linkage matrix)
     *
     * Format: cluster_i, cluster_j, distance, new_cluster_id, size
     *
     * @param tree Phylogenetic tree (must contain merge_records)
     * @param filepath Output file path
     * @return Whether successful
     */
    bool write_linkage_matrix(const Tree& tree, const std::string& filepath) const;

    /**
     * @brief Write tree statistics summary
     *
     * @param tree Phylogenetic tree
     * @param filepath Output file path
     * @return Whether successful
     */
    bool write_tree_stats(const Tree& tree, const std::string& filepath) const;

    /**
     * @brief Get output options
     */
    const TreeOutputOptions& options() const {
        return options_;
    }

    /**
     * @brief Set output options
     */
    void set_options(const TreeOutputOptions& options) {
        options_ = options;
    }

private:
    TreeOutputOptions options_;

    /**
     * @brief Process label (replace special characters)
     */
    std::string process_label(const std::string& label) const;
};

}  // namespace InterSubMod
