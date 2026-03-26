#pragma once

#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "DistanceMatrix.hpp"
#include "TreeStructure.hpp"

namespace InterSubMod {

/**
 * @brief Hierarchical clustering linkage methods
 *
 * Defines different cluster merging strategies:
 * - UPGMA: Arithmetic mean linkage, assumes molecular clock (suitable for methylation analysis)
 * - WARD: Minimum variance method, tends to produce balanced cluster sizes
 * - SINGLE: Single linkage (nearest neighbor), tends to produce chain-like structures
 * - COMPLETE: Complete linkage (furthest neighbor), tends to produce compact spherical clusters
 */
enum class LinkageMethod {
    UPGMA,    ///< Unweighted Pair Group Method with Arithmetic Mean
    WARD,     ///< Ward's minimum variance method
    SINGLE,   ///< Single linkage (minimum distance)
    COMPLETE  ///< Complete linkage (maximum distance)
};

/**
 * @brief Hierarchical clustering configuration
 */
struct ClusteringConfig {
    LinkageMethod method = LinkageMethod::UPGMA;  ///< Clustering method
    bool optimal_leaf_ordering = false;           ///< Whether to optimize leaf node ordering
    double min_branch_length = 1e-6;              ///< Minimum branch length (to avoid zero-length branches)
};

/**
 * @brief Hierarchical clustering algorithm implementation
 *
 * Supports multiple linkage methods, builds binary clustering tree from distance matrix.
 * Primarily used for methylation pattern clustering analysis.
 *
 * Usage:
 * @code
 * HierarchicalClustering clusterer(LinkageMethod::UPGMA);
 * Tree tree = clusterer.build_tree(dist_matrix, read_ids);
 * std::string newick = tree.to_newick();
 * @endcode
 */
class HierarchicalClustering {
public:
    /**
     * @brief Constructor (simplified version)
     * @param method Linkage method
     */
    explicit HierarchicalClustering(LinkageMethod method = LinkageMethod::UPGMA);

    /**
     * @brief Constructor (full version)
     * @param config Clustering configuration
     */
    explicit HierarchicalClustering(const ClusteringConfig& config);

    /**
     * @brief Build phylogenetic tree from DistanceMatrix object
     *
     * @param dist_matrix Distance matrix object
     * @param read_names Read name list (for leaf node labels)
     * @return Phylogenetic tree structure
     */
    Tree build_tree(const DistanceMatrix& dist_matrix, const std::vector<std::string>& read_names);

    /**
     * @brief Build phylogenetic tree from Eigen matrix
     *
     * @param dist_matrix N x N distance matrix (must be symmetric)
     * @param read_names Read name list
     * @return Phylogenetic tree structure
     */
    Tree build_tree(const Eigen::MatrixXd& dist_matrix, const std::vector<std::string>& read_names);

    /**
     * @brief Get current linkage method
     */
    LinkageMethod get_method() const {
        return config_.method;
    }

    /**
     * @brief Set linkage method
     */
    void set_method(LinkageMethod method) {
        config_.method = method;
    }

    /**
     * @brief Get configuration
     */
    const ClusteringConfig& config() const {
        return config_;
    }

    /**
     * @brief Convert LinkageMethod to string
     */
    static std::string method_to_string(LinkageMethod method);

    /**
     * @brief Parse LinkageMethod from string
     */
    static LinkageMethod string_to_method(const std::string& str);

private:
    ClusteringConfig config_;

    /**
     * @brief UPGMA algorithm core implementation
     *
     * Time complexity: O(N^3), Space complexity: O(N^2)
     */
    Tree build_upgma(const Eigen::MatrixXd& dist_matrix, const std::vector<std::string>& read_names);

    /**
     * @brief Ward algorithm core implementation
     */
    Tree build_ward(const Eigen::MatrixXd& dist_matrix, const std::vector<std::string>& read_names);

    /**
     * @brief Single Linkage algorithm core implementation
     */
    Tree build_single(const Eigen::MatrixXd& dist_matrix, const std::vector<std::string>& read_names);

    /**
     * @brief Complete Linkage algorithm core implementation
     */
    Tree build_complete(const Eigen::MatrixXd& dist_matrix, const std::vector<std::string>& read_names);

    /**
     * @brief Generic clustering algorithm framework
     *
     * Uses strategy pattern to perform clustering based on different distance update functions
     *
     * @param dist_matrix Distance matrix
     * @param read_names Node names
     * @param update_distance Distance update function
     * @param distance_to_height Distance to height conversion function
     * @return Phylogenetic tree
     */
    Tree build_generic(
        const Eigen::MatrixXd& dist_matrix, const std::vector<std::string>& read_names,
        std::function<double(const Eigen::MatrixXd& D, const std::vector<int>& cluster_a,
                             const std::vector<int>& cluster_b, int size_a, int size_b)>
            compute_cluster_distance,
        std::function<double(double distance)> distance_to_height = [](double d) { return d / 2.0; });

    /**
     * @brief Find minimum non-diagonal element in distance matrix
     *
     * @param D Current distance matrix
     * @param active_clusters Active cluster indices
     * @param[out] min_i First index of minimum distance
     * @param[out] min_j Second index of minimum distance
     * @return Minimum distance value
     */
    double find_minimum_distance(const Eigen::MatrixXd& D, const std::vector<bool>& active, int& min_i, int& min_j);
};

/**
 * @brief Cluster result cutter
 *
 * Generates cluster labels from tree structure based on specified cutting conditions
 * (distance threshold or number of clusters)
 */
class TreeCutter {
public:
    /**
     * @brief Cut tree by distance threshold
     *
     * @param tree Phylogenetic tree
     * @param distance_threshold Distance threshold
     * @return Cluster label for each leaf node (0-based)
     */
    static std::vector<int> cut_by_distance(const Tree& tree, double distance_threshold);

    /**
     * @brief Cut tree by number of clusters
     *
     * @param tree Phylogenetic tree
     * @param num_clusters Target number of clusters
     * @return Cluster label for each leaf node (0-based)
     */
    static std::vector<int> cut_by_num_clusters(const Tree& tree, int num_clusters);

    /**
     * @brief Compute silhouette score for each cut, select optimal cluster count
     *
     * @param tree Phylogenetic tree
     * @param dist_matrix Original distance matrix
     * @param min_k Minimum number of clusters
     * @param max_k Maximum number of clusters
     * @return pair<optimal cluster count, corresponding labels>
     */
    static std::pair<int, std::vector<int>> find_optimal_clusters(const Tree& tree, const Eigen::MatrixXd& dist_matrix,
                                                                  int min_k = 2, int max_k = 10);
};

}  // namespace InterSubMod
