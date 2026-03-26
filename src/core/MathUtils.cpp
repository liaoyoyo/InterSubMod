/**
 * @file MathUtils.cpp
 * @brief Implementation of mathematical utility functions
 */

#include "core/MathUtils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace InterSubMod {

// ============================================================================
// Log-Factorial Table Initialization
// ============================================================================

namespace {

// Build log-factorial lookup table at compile time is not feasible,
// so we use a static initialization function
std::vector<double> build_log_factorial_table() {
    std::vector<double> table(MathUtils::LOG_FACTORIAL_TABLE_SIZE);
    table[0] = 0.0;  // log(0!) = log(1) = 0
    for (int i = 1; i < MathUtils::LOG_FACTORIAL_TABLE_SIZE; ++i) {
        table[i] = table[i - 1] + std::log(static_cast<double>(i));
    }
    return table;
}

// Static table initialized once
const std::vector<double>& get_log_factorial_table() {
    static const std::vector<double> table = build_log_factorial_table();
    return table;
}

}  // anonymous namespace

// ============================================================================
// Log-Factorial
// ============================================================================

double MathUtils::log_factorial(int n) {
    if (n < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Use lookup table for small n
    if (n < LOG_FACTORIAL_TABLE_SIZE) {
        return get_log_factorial_table()[n];
    }

    // Stirling's approximation with correction terms for large n
    // log(n!) ≈ n*log(n) - n + 0.5*log(2*pi*n) + 1/(12n) - 1/(360n^3)
    double x = static_cast<double>(n);
    return x * std::log(x) - x + 0.5 * std::log(2.0 * M_PI * x) + 1.0 / (12.0 * x) - 1.0 / (360.0 * x * x * x);
}

// ============================================================================
// Log-Binomial Coefficient
// ============================================================================

double MathUtils::log_binomial(int n, int k) {
    if (k < 0 || k > n || n < 0) {
        return -std::numeric_limits<double>::infinity();  // log(0)
    }
    if (k == 0 || k == n) {
        return 0.0;  // log(1)
    }

    // log(C(n,k)) = log(n!) - log(k!) - log((n-k)!)
    return log_factorial(n) - log_factorial(k) - log_factorial(n - k);
}

// ============================================================================
// Odds Ratio
// ============================================================================

double MathUtils::odds_ratio(int a, int b, int c, int d) {
    // Haldane-Anscombe correction: add 0.5 to all cells
    double a_corr = static_cast<double>(a) + 0.5;
    double b_corr = static_cast<double>(b) + 0.5;
    double c_corr = static_cast<double>(c) + 0.5;
    double d_corr = static_cast<double>(d) + 0.5;

    return (a_corr * d_corr) / (b_corr * c_corr);
}

double MathUtils::log_odds_ratio(int a, int b, int c, int d) {
    // Haldane-Anscombe correction
    double a_corr = static_cast<double>(a) + 0.5;
    double b_corr = static_cast<double>(b) + 0.5;
    double c_corr = static_cast<double>(c) + 0.5;
    double d_corr = static_cast<double>(d) + 0.5;

    return std::log(a_corr) + std::log(d_corr) - std::log(b_corr) - std::log(c_corr);
}

// ============================================================================
// Chi-Square and Cramér's V
// ============================================================================

double MathUtils::chi_square(const std::vector<int>& table, int n_rows, int n_cols, double& min_expected) {
    if (static_cast<int>(table.size()) != n_rows * n_cols) {
        throw std::invalid_argument("Table size does not match dimensions");
    }

    // Compute row totals, column totals, and grand total
    std::vector<int> row_totals(n_rows, 0);
    std::vector<int> col_totals(n_cols, 0);
    int grand_total = 0;

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            int val = table[i * n_cols + j];
            row_totals[i] += val;
            col_totals[j] += val;
            grand_total += val;
        }
    }

    if (grand_total == 0) {
        min_expected = 0.0;
        return 0.0;
    }

    // Compute chi-square statistic
    double chi2 = 0.0;
    min_expected = std::numeric_limits<double>::max();
    double n = static_cast<double>(grand_total);

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            double observed = static_cast<double>(table[i * n_cols + j]);
            double expected = (static_cast<double>(row_totals[i]) * static_cast<double>(col_totals[j])) / n;

            if (expected > 0) {
                min_expected = std::min(min_expected, expected);
                double diff = observed - expected;
                chi2 += (diff * diff) / expected;
            }
        }
    }

    if (min_expected == std::numeric_limits<double>::max()) {
        min_expected = 0.0;
    }

    return chi2;
}

double MathUtils::cramers_v(const std::vector<int>& table, int n_rows, int n_cols, bool& is_reliable) {
    double min_expected;
    double chi2 = chi_square(table, n_rows, n_cols, min_expected);

    // Compute total sample size
    int n = std::accumulate(table.begin(), table.end(), 0);

    if (n == 0) {
        is_reliable = false;
        return 0.0;
    }

    // Check reliability: if any expected frequency < 5
    // For simplicity, we use min_expected < 5 as the criterion
    // A more rigorous check would count the proportion of cells with expected < 5
    int min_dim = std::min(n_rows, n_cols);
    if (min_dim <= 1) {
        is_reliable = false;
        return 0.0;
    }

    // Count cells with expected < 5
    std::vector<int> row_totals(n_rows, 0);
    std::vector<int> col_totals(n_cols, 0);
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            row_totals[i] += table[i * n_cols + j];
            col_totals[j] += table[i * n_cols + j];
        }
    }

    int low_expected_count = 0;
    int total_cells = n_rows * n_cols;
    double n_d = static_cast<double>(n);

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            double expected = (static_cast<double>(row_totals[i]) * static_cast<double>(col_totals[j])) / n_d;
            if (expected < 5.0) {
                ++low_expected_count;
            }
        }
    }

    // Reliable if less than 20% of cells have expected < 5
    double low_expected_ratio = static_cast<double>(low_expected_count) / static_cast<double>(total_cells);
    is_reliable = (low_expected_ratio <= 0.2);

    // Cramér's V
    double v = std::sqrt(chi2 / (n_d * static_cast<double>(min_dim - 1)));

    // Clamp to [0, 1]
    return std::clamp(v, 0.0, 1.0);
}

// ============================================================================
// Hypergeometric Distribution
// ============================================================================

double MathUtils::log_hypergeom_pmf(int k, int N, int K, int n) {
    // P(X = k) = C(K, k) * C(N-K, n-k) / C(N, n)
    // log(P) = log(C(K,k)) + log(C(N-K, n-k)) - log(C(N, n))

    // Bounds check
    int min_k = std::max(0, n - (N - K));
    int max_k = std::min(n, K);

    if (k < min_k || k > max_k) {
        return -std::numeric_limits<double>::infinity();  // log(0)
    }

    return log_binomial(K, k) + log_binomial(N - K, n - k) - log_binomial(N, n);
}

double MathUtils::hypergeom_pmf(int k, int N, int K, int n) {
    double log_p = log_hypergeom_pmf(k, N, K, n);
    if (std::isinf(log_p) && log_p < 0) {
        return 0.0;
    }
    return std::exp(log_p);
}

// Static member definition (not used with the function-based approach, but kept for compatibility)
const double MathUtils::log_factorial_table_[LOG_FACTORIAL_TABLE_SIZE] = {};

double MathUtils::log_sum_exp(const std::vector<double>& log_values) {
    if (log_values.empty()) {
        return -std::numeric_limits<double>::infinity();
    }

    // Find the maximum to use as shift (prevents overflow/underflow)
    double max_val = *std::max_element(log_values.begin(), log_values.end());

    // Guard: if max is -inf, all values are -inf
    if (!std::isfinite(max_val)) {
        return max_val;
    }

    // Stable summation: sum(exp(x_i - max)) then shift back
    double sum = 0.0;
    for (double v : log_values) {
        sum += std::exp(v - max_val);
    }

    return max_val + std::log(sum);
}

}  // namespace InterSubMod
