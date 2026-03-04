#pragma once

/**
 * @file MathUtils.hpp
 * @brief Mathematical utility functions for significance analysis
 *
 * Provides fundamental mathematical operations including:
 * - Log-factorial with Stirling approximation for large n
 * - Log-binomial coefficient
 * - Odds ratio with Haldane-Anscombe correction
 * - Cramér's V effect size
 */

#include <cmath>
#include <cstdint>
#include <vector>

namespace InterSubMod {

/**
 * @brief Mathematical utilities for statistical computations
 *
 * All functions are designed for numerical stability:
 * - Use log-scale operations to avoid overflow
 * - Handle edge cases (zeros, small samples)
 */
class MathUtils {
public:
    /**
     * @brief Compute log(n!) using precomputed table or Stirling approximation
     *
     * For n <= 1024: Uses precomputed lookup table (O(1))
     * For n > 1024: Uses Stirling's approximation with correction terms
     *
     * @param n Non-negative integer
     * @return log(n!) in natural log
     *
     * @note Matches R's lfactorial() within 1e-10 relative error
     */
    static double log_factorial(int n);

    /**
     * @brief Compute log(C(n, k)) = log(n! / (k! * (n-k)!))
     *
     * Uses log-scale arithmetic to avoid overflow for large n.
     *
     * @param n Total count (n >= 0)
     * @param k Selection count (0 <= k <= n)
     * @return log(C(n, k))
     */
    static double log_binomial(int n, int k);

    /**
     * @brief Compute odds ratio with Haldane-Anscombe correction
     *
     * For a 2x2 contingency table:
     *        | Exposed | Not Exposed |
     * -------|---------|-------------|
     * Case   |    a    |      b      |
     * Control|    c    |      d      |
     *
     * OR = (a + 0.5)(d + 0.5) / ((b + 0.5)(c + 0.5))
     *
     * The +0.5 correction (Haldane-Anscombe) prevents division by zero
     * and reduces bias in small samples.
     *
     * @param a Cell (1,1) count
     * @param b Cell (1,2) count
     * @param c Cell (2,1) count
     * @param d Cell (2,2) count
     * @return Odds ratio (always positive)
     */
    static double odds_ratio(int a, int b, int c, int d);

    /**
     * @brief Compute log odds ratio with Haldane-Anscombe correction
     *
     * @return log(OR), can be negative (OR < 1) or positive (OR > 1)
     */
    static double log_odds_ratio(int a, int b, int c, int d);

    /**
     * @brief Compute Cramér's V effect size for contingency table
     *
     * V = sqrt(chi2 / (n * (min(r, c) - 1)))
     *
     * Where:
     * - chi2 is the chi-square statistic
     * - n is total sample size
     * - r, c are number of rows and columns
     *
     * @param table Contingency table (row-major, flattened)
     * @param n_rows Number of rows
     * @param n_cols Number of columns
     * @param[out] is_reliable Set to false if >20% of cells have expected < 5
     * @return Cramér's V in [0, 1]
     */
    static double cramers_v(const std::vector<int>& table, int n_rows, int n_cols, bool& is_reliable);

    /**
     * @brief Compute chi-square statistic for contingency table
     *
     * @param table Contingency table (row-major, flattened)
     * @param n_rows Number of rows
     * @param n_cols Number of columns
     * @param[out] min_expected Minimum expected frequency in any cell
     * @return Chi-square statistic
     */
    static double chi_square(const std::vector<int>& table, int n_rows, int n_cols, double& min_expected);

    /**
     * @brief Compute hypergeometric probability P(X = k)
     *
     * P(X = k) = C(K, k) * C(N-K, n-k) / C(N, n)
     *
     * Where:
     * - N: population size
     * - K: success states in population
     * - n: number of draws
     * - k: observed successes
     *
     * @return Log probability
     */
    static double log_hypergeom_pmf(int k, int N, int K, int n);

    /**
     * @brief Compute hypergeometric probability (non-log)
     */
    static double hypergeom_pmf(int k, int N, int K, int n);

    /**
     * @brief Size of precomputed log-factorial table
     */
    static constexpr int LOG_FACTORIAL_TABLE_SIZE = 1025;

private:
    // Precomputed log-factorial table for n <= 1024
    static const double log_factorial_table_[LOG_FACTORIAL_TABLE_SIZE];

    // Initialize the lookup table
    static const double* init_log_factorial_table();
};

}  // namespace InterSubMod
