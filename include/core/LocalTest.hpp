#pragma once

/**
 * @file LocalTest.hpp
 * @brief Local (per-cluster) association tests
 *
 * Implements one-vs-rest Fisher's exact tests for each cluster
 * to identify which specific clusters drive the global association.
 */

#include <vector>

#include "FisherExact.hpp"
#include "Stats.hpp"

namespace InterSubMod {

/**
 * @brief Configuration for local tests
 */
struct LocalTestConfig {
    FisherConfig fisher_config;
};

/**
 * @brief Local test calculator
 *
 * Performs one-vs-rest tests for each cluster to identify
 * which clusters are enriched/depleted for specific labels.
 */
class LocalTest {
public:
    explicit LocalTest(const LocalTestConfig& config = LocalTestConfig());

    /**
     * @brief Run local tests for all dimensions
     *
     * @param cluster_labels Cluster ID for each read
     * @param full_labels Full label information
     * @return LocalTestResult with per-cluster statistics
     */
    LocalTestResult run(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels);

    /**
     * @brief Test one-vs-rest for allele dimension
     */
    void test_allele_local(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                           std::vector<ClusterStats>& stats);

    /**
     * @brief Test one-vs-rest for HP dimension
     */
    void test_hp_local(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                       std::vector<ClusterStats>& stats);

    /**
     * @brief Test one-vs-rest for sample dimension
     */
    void test_sample_local(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                           std::vector<ClusterStats>& stats);

private:
    LocalTestConfig config_;
    FisherExact fisher_;

    /**
     * @brief Build 2x2 table for one-vs-rest comparison
     *
     * Table structure:
     *        | This Cluster | Other Clusters |
     * -------|--------------|----------------|
     * Target |      a       |       b        |
     * Other  |      c       |       d        |
     *
     * @param cluster_id Target cluster
     * @param cluster_labels All cluster labels
     * @param target_label Target label value
     * @param binary_labels Binary label (0 or 1)
     * @param valid_mask Validity mask
     */
    ContingencyTable2x2 build_one_vs_rest_table(int cluster_id, const std::vector<int>& cluster_labels,
                                                int target_label, const std::vector<int>& binary_labels,
                                                const std::vector<bool>& valid_mask);

    /**
     * @brief Compute delta proportion: P(target|cluster) - P(target|rest)
     */
    double compute_delta_proportion(const ContingencyTable2x2& table);

    /**
     * @brief Count reads per cluster
     */
    std::vector<ClusterStats> initialize_cluster_stats(const std::vector<int>& cluster_labels,
                                                       const std::vector<FullLabel>& full_labels);
};

}  // namespace InterSubMod
