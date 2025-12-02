#pragma once

#include <vector>
#include <string>
#include <Eigen/Dense>
#include "Types.hpp"
#include "DataStructs.hpp"

namespace InterSubMod {

// Forward declaration
class MethylationMatrix;

/**
 * @brief Configuration for distance matrix calculation.
 * 
 * Allows fine-grained control over distance computation parameters
 * to ensure reproducibility and customization.
 */
struct DistanceConfig {
    DistanceMetricType metric = DistanceMetricType::NHD;  ///< Distance metric
    int min_common_coverage = 5;     ///< Minimum common CpG sites required (C_min)
    NanDistanceStrategy nan_strategy = NanDistanceStrategy::MAX_DIST;  ///< Strategy for invalid pairs
    double max_distance_value = 1.0; ///< Value for MAX_DIST strategy (normalized metrics)
    bool use_binary_matrix = true;   ///< Use binary matrix (true) or raw matrix (false)
    double binary_threshold_high = 0.8;  ///< Threshold for methylated (binary)
    double binary_threshold_low = 0.2;   ///< Threshold for unmethylated (binary)
    int num_threads = 1;             ///< Number of threads for parallel computation
    
    // For correlation distance
    bool pearson_center = true;      ///< Use mean-centered Pearson (true) or raw correlation
    
    // For Jaccard distance
    bool jaccard_include_unmeth = false;  ///< Include unmethylated sites in Jaccard
};

/**
 * @brief Stores pairwise distances between reads within a region.
 * 
 * Used as input for hierarchical clustering. Supports multiple distance
 * metrics and strand-aware computation.
 */
class DistanceMatrix {
public:
    int region_id = -1;
    std::vector<int> read_ids;        ///< Read IDs corresponding to rows/cols
    Eigen::MatrixXd dist_matrix;      ///< Symmetric NxN distance matrix
    int min_common_coverage = 5;      ///< Minimum overlapping sites used
    DistanceMetricType metric_type;   ///< Metric used
    NanDistanceStrategy nan_strategy; ///< NaN strategy used
    
    // Statistics
    int num_valid_pairs = 0;          ///< Number of valid (computed) pairs
    int num_invalid_pairs = 0;        ///< Number of pairs with insufficient overlap
    double avg_common_coverage = 0.0; ///< Average common CpG sites per pair
    
    DistanceMatrix() = default;
    
    /**
     * @brief Check if the matrix is empty (no reads).
     */
    bool empty() const { return read_ids.empty(); }
    
    /**
     * @brief Get the number of reads in the matrix.
     */
    int size() const { return static_cast<int>(read_ids.size()); }
    
    /**
     * @brief Computes pairwise distances from a MethylationMatrix.
     * 
     * @param methyl_mat The input methylation data.
     * @param type Distance metric (e.g., NHD).
     * @param min_cov Minimum common CpG sites required to calculate a valid distance.
     * @param nan_strategy Strategy to handle pairs with insufficient overlap.
     */
    void compute_from_methylation(
        const MethylationMatrix& methyl_mat, 
        DistanceMetricType type, 
        int min_cov, 
        NanDistanceStrategy nan_strategy
    );

    /**
     * @brief Computes pairwise distances with full configuration.
     * 
     * @param methyl_mat The input methylation data.
     * @param config Distance calculation configuration.
     */
    void compute_from_methylation(
        const MethylationMatrix& methyl_mat, 
        const DistanceConfig& config
    );

    /**
     * @brief Computes pairwise distances for a specific subset of reads.
     * 
     * @param methyl_mat The input methylation data.
     * @param row_indices Indices of rows in methyl_mat to include.
     * @param type Distance metric.
     * @param min_cov Minimum common coverage.
     * @param nan_strategy Strategy for missing data.
     */
    void compute_subset(
        const MethylationMatrix& methyl_mat, 
        const std::vector<int>& row_indices, 
        DistanceMetricType type, 
        int min_cov, 
        NanDistanceStrategy nan_strategy
    );
    
    /**
     * @brief Computes pairwise distances for a specific subset with full configuration.
     * 
     * @param methyl_mat The input methylation data.
     * @param row_indices Indices of rows to include.
     * @param config Distance calculation configuration.
     */
    void compute_subset(
        const MethylationMatrix& methyl_mat,
        const std::vector<int>& row_indices,
        const DistanceConfig& config
    );
    
    /**
     * @brief Get distance between two reads by their indices in this matrix.
     */
    double get_distance(int i, int j) const {
        if (i < 0 || i >= size() || j < 0 || j >= size()) return NAN;
        return dist_matrix(i, j);
    }
    
    /**
     * @brief Write the distance matrix to a CSV file.
     * 
     * @param filepath Output file path.
     * @param include_header Include header row with read IDs.
     */
    void write_csv(const std::string& filepath, bool include_header = true) const;
    
    /**
     * @brief Write summary statistics to a text file.
     */
    void write_stats(const std::string& filepath) const;
};

/**
 * @brief Utility class for computing strand-specific distance matrices.
 * 
 * Computes separate distance matrices for forward and reverse strand reads,
 * as methylation patterns may differ between strands.
 */
class DistanceCalculator {
public:
    /**
     * @brief Construct a DistanceCalculator with configuration.
     */
    explicit DistanceCalculator(const DistanceConfig& config) : config_(config) {}
    
    /**
     * @brief Compute distance matrix for all reads.
     * 
     * @param methyl_mat The methylation matrix.
     * @param reads Read information (for strand info).
     * @return The computed distance matrix.
     */
    DistanceMatrix compute(
        const MethylationMatrix& methyl_mat,
        const std::vector<ReadInfo>& reads
    );
    
    /**
     * @brief Compute strand-specific distance matrices.
     * 
     * @param methyl_mat The methylation matrix.
     * @param reads Read information (for strand info).
     * @return Pair of (forward_matrix, reverse_matrix).
     */
    std::pair<DistanceMatrix, DistanceMatrix> compute_strand_specific(
        const MethylationMatrix& methyl_mat,
        const std::vector<ReadInfo>& reads
    );
    
    /**
     * @brief Get the configuration.
     */
    const DistanceConfig& config() const { return config_; }
    
    /**
     * @brief Convert DistanceMetricType to string.
     */
    static std::string metric_to_string(DistanceMetricType type);
    
    /**
     * @brief Convert string to DistanceMetricType.
     */
    static DistanceMetricType string_to_metric(const std::string& str);
    
private:
    DistanceConfig config_;
};

} // namespace InterSubMod
