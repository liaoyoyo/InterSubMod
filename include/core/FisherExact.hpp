#pragma once

/**
 * @file FisherExact.hpp
 * @brief Fisher's exact test implementation
 *
 * Provides:
 * - 2x2 Fisher's exact test (two-sided, R-compatible)
 * - RxC Fisher-Freeman-Halton test with Monte Carlo
 *
 * Two-sided p-value definition (R compatible):
 * Sum of all probabilities P(table) <= P(observed table)
 */

#include <cstdint>
#include <random>
#include <vector>

#include "Stats.hpp"

namespace InterSubMod {

/**
 * @brief Configuration for Fisher's exact test
 */
struct FisherConfig {
    // Monte Carlo settings
    int min_mc_samples = 2000;
    int max_mc_samples = 10000;

    // Early stopping: stop if 99% CI is completely above or below alpha
    double alpha = 0.05;
    double ci_confidence = 0.99;
    bool use_early_stopping = true;

    // Random seed (0 = use random device)
    uint64_t seed = 0;
};

/**
 * @brief Fisher's exact test calculator
 */
class FisherExact {
public:
    explicit FisherExact(const FisherConfig& config = FisherConfig());

    /**
     * @brief Compute 2x2 Fisher's exact test (two-sided)
     *
     * Implementation strictly follows R's fisher.test():
     * - Two-sided: sum of P(table) for all tables with P <= P(observed)
     * - Uses hypergeometric distribution
     *
     * @param table 2x2 contingency table
     * @return FisherResult with p-value and odds ratio
     */
    FisherResult test_2x2(const ContingencyTable2x2& table);

    /**
     * @brief Compute RxC Fisher-Freeman-Halton test (Monte Carlo)
     *
     * Uses Monte Carlo sampling to estimate p-value for general RxC tables.
     * Implements 99% Binomial CI stopping rule.
     *
     * @param table RxC contingency table
     * @return FisherResult with p-value and CI bounds
     */
    FisherResult test_rxc_mc(const ContingencyTableRxC& table);

    /**
     * @brief Set random seed
     */
    void set_seed(uint64_t seed);

    /**
     * @brief Get current configuration
     */
    const FisherConfig& config() const {
        return config_;
    }

private:
    FisherConfig config_;
    std::mt19937_64 rng_;

    /**
     * @brief Compute log-probability of a 2x2 table under null hypothesis
     */
    static double log_table_prob_2x2(const ContingencyTable2x2& table);

    /**
     * @brief Compute log-probability of an RxC table under null hypothesis
     *
     * P(table | marginals) = prod(row_i!) * prod(col_j!) / (n! * prod(cell_ij!))
     */
    static double log_table_prob_rxc(const ContingencyTableRxC& table);

    /**
     * @brief Generate random table with fixed marginals (Patefield algorithm)
     */
    ContingencyTableRxC sample_table_with_marginals(const std::vector<int>& row_totals,
                                                    const std::vector<int>& col_totals);

    /**
     * @brief Compute 99% Binomial confidence interval for proportion
     *
     * Uses Wilson score interval.
     *
     * @param n_success Number of successes
     * @param n_total Total trials
     * @param confidence Confidence level (default 0.99)
     * @return {lower, upper} bounds
     */
    static std::pair<double, double> binomial_ci(int n_success, int n_total, double confidence = 0.99);
};

}  // namespace InterSubMod
