#pragma once

/**
 * @file LabelTest.hpp
 * @brief Label-based significance testing (Label-First verification)
 *
 * Implements Label-First verification path:
 * - Computes Delta = mean(between-group distance) - mean(within-group distance)
 * - Runs permutation test to assess if labels explain distance structure
 * - Supports multiple label dimensions: HP and Allele
 *
 * This complements the Cluster-First approach by directly testing
 * whether known labels (HP/Allele) correlate with methylation distance.
 */

#include <Eigen/Dense>
#include <random>
#include <string>
#include <vector>

#include "Stats.hpp"
#include "Types.hpp"

namespace InterSubMod {

/**
 * @brief Configuration for label-based tests
 */
struct LabelTestConfig {
    // Permutation settings
    int n_permutations = 999;

    // Minimum requirements
    int min_reads_per_group = 3;  // Minimum reads required in each group
    int min_total_reads = 10;     // Minimum total reads for valid test

    // Significance threshold
    double alpha = 0.05;

    // Random seed
    uint64_t seed = 0;
};

// LabelDeltaResult is defined in Stats.hpp

/**
 * @brief Aggregated result for all label dimensions
 */
struct LabelTestResult {
    // Per-dimension results
    LabelDeltaResult hp_result;     // HP binary (now deprecated, kept for compatibility)
    LabelDeltaResult allele_result; // ALT vs REF

    // NEW: HP multi-group PERMANOVA result
    double hp_permanova_f = 0.0;     // Pseudo-F statistic
    double hp_permanova_p = 1.0;     // P-value from multi-group PERMANOVA
    bool hp_permanova_sig = false;   // Significant at alpha
    int hp_n_groups = 0;             // Number of HP groups with enough reads

    // Best (most significant) dimension
    std::string dominant_dimension;  // "hp", "allele", or "none"
    double best_delta = 0.0;
    double best_p_value = 1.0;
    bool any_significant = false;

    // Overall validity
    bool valid = true;
    std::string invalid_reason;
};

/**
 * @brief Label-based distance structure test
 *
 * Tests whether known biological labels (HP/Allele) explain the
 * methylation distance structure. This is the "Label-First" path
 * in the bidirectional verification framework.
 *
 * Key metric: Delta = d_between - d_within
 *   - Positive Delta indicates labels explain distance structure
 *   - Permutation test assesses statistical significance
 */
class LabelTest {
public:
    explicit LabelTest(const LabelTestConfig& config = LabelTestConfig());

    /**
     * @brief Run all label-based tests
     *
     * Tests HP (HP1 vs HP2) and Allele (ALT vs REF) dimensions.
     *
     * @param dist_matrix N x N distance matrix (may contain NaN)
     * @param full_labels Full label information for each read
     * @return LabelTestResult Aggregated results for all dimensions
     */
    LabelTestResult test_all(const Eigen::MatrixXd& dist_matrix,
                             const std::vector<FullLabel>& full_labels);

    /**
     * @brief Test a specific binary grouping
     *
     * @param dist_matrix N x N distance matrix
     * @param group_labels Binary group assignment (0 or 1, -1 for excluded)
     * @return LabelDeltaResult Delta statistics and permutation p-value
     */
    LabelDeltaResult test_binary_groups(const Eigen::MatrixXd& dist_matrix,
                                        const std::vector<int>& group_labels);

    /**
     * @brief Compute Delta without permutation test (for efficiency)
     *
     * Useful for quick preliminary check before full test.
     *
     * @param dist_matrix N x N distance matrix
     * @param group_labels Binary group assignment
     * @return LabelDeltaResult (without p-value/significance)
     */
    LabelDeltaResult compute_delta(const Eigen::MatrixXd& dist_matrix,
                                   const std::vector<int>& group_labels);

    /**
     * @brief Set random seed for permutation tests
     */
    void set_seed(uint64_t seed);

    const LabelTestConfig& config() const { return config_; }

private:
    LabelTestConfig config_;
    std::mt19937_64 rng_;

    std::vector<int> hp_to_binary_labels(const std::vector<FullLabel>& full_labels);

    /**
     * @brief Convert HP tags to multi-group labels for PERMANOVA
     *
     * Mapping: HP1->0, HP2->1, HP1-1->2, HP2-1->3, HP3->4, HP0/empty->5
     * Returns -1 for reads in groups with < min_reads_per_group members (excluded)
     */
    std::vector<int> hp_to_multigroup_labels(const std::vector<FullLabel>& full_labels);

    /**
     * @brief Convert Allele to binary group labels
     *
     * ALT -> 0, REF -> 1, UNKNOWN -> -1 (excluded)
     */
    std::vector<int> allele_to_binary_labels(const std::vector<FullLabel>& full_labels);

    /**
     * @brief Compute within-group and between-group mean distances
     *
     * Only considers valid (non-NaN) distance pairs where both reads have valid labels.
     */
    void compute_group_distances(const Eigen::MatrixXd& dist_matrix,
                                 const std::vector<int>& group_labels,
                                 double& within_mean, double& between_mean,
                                 int& n_pairs_within, int& n_pairs_between);

    double permutation_test(const Eigen::MatrixXd& dist_matrix,
                            const std::vector<int>& group_labels,
                            double observed_delta);

    /**
     * @brief Run multi-group PERMANOVA test for HP labels
     *
     * Tests if HP labels (HP1, HP2, HP1-1, HP2-1, HP3, HP0) explain distance structure.
     * Uses pseudo-F statistic with permutation test.
     */
    void test_hp_multigroup_permanova(const Eigen::MatrixXd& dist_matrix,
                                      const std::vector<int>& hp_labels,
                                      LabelTestResult& result);
};

}  // namespace InterSubMod
