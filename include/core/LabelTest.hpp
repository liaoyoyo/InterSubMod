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

// LabelDeltaResult, MultiStageHPResult are defined in Stats.hpp

/**
 * @brief Aggregated result for all label dimensions
 *
 * Contains results from multi-stage HP verification (Stages 1, 2, 4)
 * and Allele verification (Stage 3).
 */
struct LabelTestResult {
    // Multi-stage HP verification (NEW)
    MultiStageHPResult hp_multistage;

    // Stage 3: Allele test (ALT vs REF)
    LabelDeltaResult allele_result;

    // Deprecated fields (kept for backward compatibility)
    LabelDeltaResult hp_result;      // Use hp_multistage.merged instead
    double hp_permanova_f = 0.0;     // Use hp_multistage.fine_f instead
    double hp_permanova_p = 1.0;     // Use hp_multistage.fine_p instead
    bool hp_permanova_sig = false;   // Use hp_multistage.fine_sig instead
    int hp_n_groups = 0;             // Use hp_multistage.fine_n_groups instead

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

    // ========================================================================
    // Label Conversion Functions
    // ========================================================================

    /**
     * @brief Convert HP tags to binary labels (deprecated)
     *
     * HP1/HP1-1 -> 0, HP2/HP2-1 -> 1, others -> -1 (excluded)
     */
    std::vector<int> hp_to_binary_labels(const std::vector<FullLabel>& full_labels);

    /**
     * @brief Convert HP tags to merged family labels (Stage 1)
     *
     * HP1-family (HP1 + HP1-1) -> 0
     * HP2-family (HP2 + HP2-1) -> 1
     * HP3, HP0, empty -> -1 (excluded)
     */
    std::vector<int> hp_to_merged_labels(const std::vector<FullLabel>& full_labels);

    /**
     * @brief Convert HP tags to fine-grained labels (Stage 2)
     *
     * HP1 -> 0, HP1-1 -> 1, HP2 -> 2, HP2-1 -> 3
     * HP3, HP0, empty -> -1 (excluded)
     * Returns -1 for groups with < min_reads_per_group members
     */
    std::vector<int> hp_to_fine_labels(const std::vector<FullLabel>& full_labels);

    /**
     * @brief Convert HP tags to multi-group labels for PERMANOVA (deprecated)
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

    // ========================================================================
    // Distance Computation
    // ========================================================================

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

    // ========================================================================
    // Multi-Stage HP Tests (NEW)
    // ========================================================================

    /**
     * @brief Stage 1: HP Family Merged Test
     *
     * Tests binary grouping: (HP1 + HP1-1) vs (HP2 + HP2-1)
     * Uses Delta test with permutation.
     */
    void test_hp_merged(const Eigen::MatrixXd& dist_matrix,
                        const std::vector<FullLabel>& full_labels,
                        MultiStageHPResult& result);

    /**
     * @brief Stage 2: HP Fine-Grained Test
     *
     * Tests 4-group structure: HP1, HP1-1, HP2, HP2-1
     * Uses multi-group PERMANOVA and computes pairwise distances.
     */
    void test_hp_fine_grained(const Eigen::MatrixXd& dist_matrix,
                              const std::vector<FullLabel>& full_labels,
                              MultiStageHPResult& result);

    /**
     * @brief Stage 4: Unassigned Affinity Test
     *
     * Tests if HP3/HP0 reads are closer to HP1-family or HP2-family.
     * Computes affinity score and permutation p-value.
     */
    void test_unassigned_affinity(const Eigen::MatrixXd& dist_matrix,
                                  const std::vector<FullLabel>& full_labels,
                                  const std::vector<int>& merged_labels,
                                  MultiStageHPResult& result);

    /**
     * @brief Run multi-group PERMANOVA test for HP labels (deprecated)
     *
     * Tests if HP labels (HP1, HP2, HP1-1, HP2-1, HP3, HP0) explain distance structure.
     * Uses pseudo-F statistic with permutation test.
     */
    void test_hp_multigroup_permanova(const Eigen::MatrixXd& dist_matrix,
                                      const std::vector<int>& hp_labels,
                                      LabelTestResult& result);
};

}  // namespace InterSubMod
