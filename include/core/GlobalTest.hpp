#pragma once

/**
 * @file GlobalTest.hpp
 * @brief Global association test between cluster labels and biological labels
 *
 * Implements Phase 2 of significance analysis:
 * - Fisher-Freeman-Halton test for RxC contingency tables
 * - Chi-square test (when valid)
 * - Cramér's V effect size
 * - Gating logic to skip costly downstream tests
 */

#include <vector>

#include "FisherExact.hpp"
#include "Stats.hpp"

namespace InterSubMod {

/**
 * @brief Configuration for global association test
 */
struct GlobalTestConfig {
    // Fisher test settings
    FisherConfig fisher_config;

    // Gating: skip local/structure tests if global p > this threshold
    double gating_p_threshold = 0.1;

    // Chi-square validity: require this fraction of cells with expected >= 5
    double chi_square_min_valid_ratio = 0.8;
};

/**
 * @brief Global association test calculator
 *
 * Tests whether cluster membership is associated with biological labels.
 * Priority: ALT/REF > HP > Tumor/Normal
 */
class GlobalTest {
public:
    explicit GlobalTest(const GlobalTestConfig& config = GlobalTestConfig());

    /**
     * @brief Test association between clusters and ALT/REF labels
     *
     * @param cluster_labels Cluster ID for each read (0-indexed)
     * @param full_labels Full label information for each read
     * @return GlobalTestResult
     */
    GlobalTestResult test_allele(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels);

    /**
     * @brief Test association between clusters and HP labels (HP1 vs HP2)
     */
    GlobalTestResult test_hp(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels);

    /**
     * @brief Test association between clusters and Tumor/Normal labels
     */
    GlobalTestResult test_sample(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels);

    /**
     * @brief Run all global tests with priority ordering
     *
     * @param cluster_labels Cluster ID for each read
     * @param full_labels Full label information
     * @return Tuple of (alt_result, hp_result, sample_result)
     */
    std::tuple<GlobalTestResult, GlobalTestResult, GlobalTestResult> test_all(
        const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels);

    /**
     * @brief Get configuration
     */
    const GlobalTestConfig& config() const {
        return config_;
    }

private:
    GlobalTestConfig config_;
    FisherExact fisher_;

    /**
     * @brief Build contingency table from cluster labels and binary labels
     *
     * @param cluster_labels Cluster ID for each read
     * @param binary_labels Binary label (0 or 1) for each read
     * @param valid_mask True if read should be included
     * @return ContingencyTableRxC with clusters as rows, labels as columns
     */
    ContingencyTableRxC build_contingency_table(const std::vector<int>& cluster_labels,
                                                const std::vector<int>& binary_labels,
                                                const std::vector<bool>& valid_mask);

    /**
     * @brief Compute global test result from contingency table
     */
    GlobalTestResult compute_from_table(const ContingencyTableRxC& table);

    /**
     * @brief Apply gating logic
     */
    void apply_gating(GlobalTestResult& result);
};

}  // namespace InterSubMod
