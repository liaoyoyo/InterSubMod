#pragma once

/**
 * @file SignificanceAnalyzer.hpp
 * @brief Main entry point for significance analysis
 *
 * Integrates all significance tests:
 * - Global association tests (Fisher-Freeman-Halton)
 * - Local tests (One-vs-Rest)
 * - Structure tests (PERMANOVA, Dispersion)
 * - Bootstrap stability
 */

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "Bootstrap.hpp"
#include "GlobalTest.hpp"
#include "LabelTest.hpp"
#include "LocalTest.hpp"
#include "Stats.hpp"
#include "StructureTest.hpp"

namespace InterSubMod {

/**
 * @brief Configuration for SignificanceAnalyzer
 */
struct SignificanceConfig {
    // Enable/disable components
    bool enable_global_test = true;
    bool enable_local_test = true;
    bool enable_permanova = true;
    bool enable_dispersion = true;
    bool enable_bootstrap = false;  // Expensive, disabled by default

    // Component configs
    GlobalTestConfig global_config;
    LocalTestConfig local_config;
    StructureTestConfig structure_config;
    BootstrapConfig bootstrap_config;
    LabelTestConfig label_config;  // NEW: Label-First test config

    // Output options
    std::string run_id;
    std::string vcf_id;
    std::string bam_id;

    // Random seed
    uint64_t seed = 0;
};

/**
 * @brief Main significance analyzer
 *
 * Orchestrates all significance tests for a region.
 */
class SignificanceAnalyzer {
public:
    explicit SignificanceAnalyzer(const SignificanceConfig& config = SignificanceConfig());

    /**
     * @brief Analyze significance for a region
     *
     * @param cluster_labels Cluster assignment for each read
     * @param full_labels Full label information for each read
     * @param dist_matrix Distance matrix
     * @param methylation_matrix Methylation values (reads x CpG sites)
     * @param region_id Region identifier
     * @param anchor_key Anchor key (chrom:pos:ref:alt)
     * @return SignificanceResult
     */
    SignificanceResult analyze(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                               const Eigen::MatrixXd& dist_matrix, const Eigen::MatrixXd& methylation_matrix,
                               int region_id, const std::string& anchor_key);

    /**
     * @brief Simple analysis without structure tests
     *
     * For quick analysis when only contingency tests are needed.
     */
    SignificanceResult analyze_simple(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                                      int region_id, const std::string& anchor_key);

    /**
     * @brief Compute heuristic score
     */
    double compute_heuristic_score(const SignificanceResult& result);

    /**
     * @brief Set seed for all random components
     */
    void set_seed(uint64_t seed);

    const SignificanceConfig& config() const {
        return config_;
    }

private:
    SignificanceConfig config_;

    std::unique_ptr<GlobalTest> global_test_;
    std::unique_ptr<LocalTest> local_test_;
    std::unique_ptr<StructureTest> structure_test_;
    std::unique_ptr<Bootstrap> bootstrap_;
    std::unique_ptr<LabelTest> label_test_;  // NEW: Label-First test

    /**
     * @brief Initialize sub-analyzers with current config
     */
    void init_analyzers();
};

}  // namespace InterSubMod
