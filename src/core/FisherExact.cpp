/**
 * @file FisherExact.cpp
 * @brief Implementation of Fisher's exact test
 *
 * Two-sided p-value follows R's fisher.test() definition:
 * Sum of all table probabilities P(table) <= P(observed)
 */

#include "core/FisherExact.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "core/MathUtils.hpp"

namespace InterSubMod {

// ============================================================================
// Constructor
// ============================================================================

FisherExact::FisherExact(const FisherConfig& config) : config_(config) {
    if (config_.seed == 0) {
        std::random_device rd;
        rng_.seed(rd());
    } else {
        rng_.seed(config_.seed);
    }
}

void FisherExact::set_seed(uint64_t seed) {
    config_.seed = seed;
    rng_.seed(seed);
}

// ============================================================================
// 2x2 Fisher's Exact Test
// ============================================================================

double FisherExact::log_table_prob_2x2(const ContingencyTable2x2& table) {
    // P(table | marginals) = C(n1, a) * C(n2, c) / C(N, m1)
    // where n1 = a+b (row1), n2 = c+d (row2), m1 = a+c (col1), N = total

    // Using hypergeometric: P(X = a) = C(m1, a) * C(N-m1, n1-a) / C(N, n1)
    int N = table.grand_total();
    int K = table.col1_total();  // m1 = a + c
    int n = table.row1_total();  // n1 = a + b
    int k = table.a;

    return MathUtils::log_hypergeom_pmf(k, N, K, n);
}

FisherResult FisherExact::test_2x2(const ContingencyTable2x2& table) {
    FisherResult result;

    if (!table.is_valid()) {
        result.p_value = 1.0;
        result.odds_ratio = 1.0;
        result.log_odds_ratio = 0.0;
        return result;
    }

    // Compute odds ratio
    result.odds_ratio = MathUtils::odds_ratio(table.a, table.b, table.c, table.d);
    result.log_odds_ratio = MathUtils::log_odds_ratio(table.a, table.b, table.c, table.d);

    // Fixed marginals
    int N = table.grand_total();
    int K = table.col1_total();  // a + c
    int n = table.row1_total();  // a + b

    // Support of k (value of cell a)
    int k_min = std::max(0, n - (N - K));
    int k_max = std::min(n, K);

    // Observed probability (log scale)
    double log_p_obs = log_table_prob_2x2(table);

    // Two-sided p-value: sum of P(table) for all tables with P <= P(observed)
    // We enumerate all possible values of k and sum probabilities that are <= p_obs
    double p_value = 0.0;

    for (int k = k_min; k <= k_max; ++k) {
        // Construct table with k in cell (0,0)
        ContingencyTable2x2 test_table;
        test_table.a = k;
        test_table.b = n - k;
        test_table.c = K - k;
        test_table.d = (N - K) - (n - k);

        double log_p_k = log_table_prob_2x2(test_table);

        // Two-sided: include if P(k) <= P(observed)
        // Using log scale: log_p_k <= log_p_obs + epsilon (for numerical stability)
        if (log_p_k <= log_p_obs + 1e-10) {
            p_value += std::exp(log_p_k);
        }
    }

    // Clamp to [0, 1]
    result.p_value = std::clamp(p_value, 0.0, 1.0);

    return result;
}

// ============================================================================
// RxC Fisher-Freeman-Halton Test (Monte Carlo)
// ============================================================================

double FisherExact::log_table_prob_rxc(const ContingencyTableRxC& table) {
    // P(table | marginals) ∝ prod(row_i!) * prod(col_j!) / (N! * prod(cell_ij!))

    int n_rows = table.n_rows;
    int n_cols = table.n_cols;

    // Compute row totals
    std::vector<int> row_totals(n_rows, 0);
    std::vector<int> col_totals(n_cols, 0);
    int N = 0;

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            int val = table.get(i, j);
            row_totals[i] += val;
            col_totals[j] += val;
            N += val;
        }
    }

    // log(P) = sum(log(row_i!)) + sum(log(col_j!)) - log(N!) - sum(log(cell_ij!))
    double log_p = 0.0;

    for (int i = 0; i < n_rows; ++i) {
        log_p += MathUtils::log_factorial(row_totals[i]);
    }
    for (int j = 0; j < n_cols; ++j) {
        log_p += MathUtils::log_factorial(col_totals[j]);
    }
    log_p -= MathUtils::log_factorial(N);

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            log_p -= MathUtils::log_factorial(table.get(i, j));
        }
    }

    return log_p;
}

ContingencyTableRxC FisherExact::sample_table_with_marginals(const std::vector<int>& row_totals,
                                                             const std::vector<int>& col_totals) {
    // Patefield's algorithm (1981) for sampling tables with fixed marginals
    // This is a simplified version using sequential filling

    int n_rows = static_cast<int>(row_totals.size());
    int n_cols = static_cast<int>(col_totals.size());

    ContingencyTableRxC result;
    result.n_rows = n_rows;
    result.n_cols = n_cols;
    result.cells.resize(n_rows * n_cols, 0);

    // Working copies of marginals
    std::vector<int> row_remain = row_totals;
    std::vector<int> col_remain = col_totals;

    // Fill cells row by row
    for (int i = 0; i < n_rows - 1; ++i) {
        for (int j = 0; j < n_cols - 1; ++j) {
            // Parameters for hypergeometric sampling
            // We need to fill cell (i,j) given remaining row_remain[i] and col_remain[j]

            int n = row_remain[i];  // items to place in row i
            int K = col_remain[j];  // available in column j
            int N = 0;              // total remaining for columns j to end
            for (int jj = j; jj < n_cols; ++jj) {
                N += col_remain[jj];
            }

            // Sample from hypergeometric
            int k_min = std::max(0, n - (N - K));
            int k_max = std::min(n, K);

            // Build CDF
            std::vector<double> cdf;
            double cumsum = 0.0;
            for (int k = k_min; k <= k_max; ++k) {
                double p = MathUtils::hypergeom_pmf(k, N, K, n);
                cumsum += p;
                cdf.push_back(cumsum);
            }

            // Sample
            std::uniform_real_distribution<double> dist(0.0, cumsum);
            double u = dist(rng_);
            int sampled_k = k_min;
            for (int idx = 0; idx < static_cast<int>(cdf.size()); ++idx) {
                if (u <= cdf[idx]) {
                    sampled_k = k_min + idx;
                    break;
                }
            }

            result.cells[i * n_cols + j] = sampled_k;
            row_remain[i] -= sampled_k;
            col_remain[j] -= sampled_k;
        }
        // Last column of this row is determined
        result.cells[i * n_cols + (n_cols - 1)] = row_remain[i];
        col_remain[n_cols - 1] -= row_remain[i];
        row_remain[i] = 0;
    }

    // Last row is determined by column remainders
    for (int j = 0; j < n_cols; ++j) {
        result.cells[(n_rows - 1) * n_cols + j] = col_remain[j];
    }

    return result;
}

std::pair<double, double> FisherExact::binomial_ci(int n_success, int n_total, double confidence) {
    // Wilson score interval
    if (n_total == 0) {
        return {0.0, 1.0};
    }

    double p_hat = static_cast<double>(n_success) / static_cast<double>(n_total);
    double n = static_cast<double>(n_total);

    // z value for (1 + confidence) / 2 quantile
    // For 99% CI: z ≈ 2.576
    // For 95% CI: z ≈ 1.96
    double z;
    if (confidence >= 0.99) {
        z = 2.576;
    } else if (confidence >= 0.95) {
        z = 1.96;
    } else {
        z = 1.645;
    }

    double z2 = z * z;
    double denom = 1.0 + z2 / n;
    double center = (p_hat + z2 / (2.0 * n)) / denom;
    double margin = (z / denom) * std::sqrt(p_hat * (1.0 - p_hat) / n + z2 / (4.0 * n * n));

    double lower = std::max(0.0, center - margin);
    double upper = std::min(1.0, center + margin);

    return {lower, upper};
}

FisherResult FisherExact::test_rxc_mc(const ContingencyTableRxC& table) {
    FisherResult result;
    result.n_samples = 0;
    result.early_stopped = false;

    if (!table.is_valid()) {
        result.p_value = 1.0;
        return result;
    }

    // Compute marginals
    std::vector<int> row_totals(table.n_rows, 0);
    std::vector<int> col_totals(table.n_cols, 0);

    for (int i = 0; i < table.n_rows; ++i) {
        for (int j = 0; j < table.n_cols; ++j) {
            int val = table.get(i, j);
            row_totals[i] += val;
            col_totals[j] += val;
        }
    }

    // Log probability of observed table
    double log_p_obs = log_table_prob_rxc(table);

    // Monte Carlo sampling
    int n_extreme = 0;  // Count of tables with P <= P(observed)

    for (int iter = 0; iter < config_.max_mc_samples; ++iter) {
        // Sample random table
        ContingencyTableRxC sampled = sample_table_with_marginals(row_totals, col_totals);
        double log_p_sample = log_table_prob_rxc(sampled);

        // Count if as extreme or more extreme
        if (log_p_sample <= log_p_obs + 1e-10) {
            ++n_extreme;
        }

        result.n_samples = iter + 1;

        // Early stopping check after minimum samples
        if (config_.use_early_stopping && result.n_samples >= config_.min_mc_samples) {
            auto [ci_lower, ci_upper] = binomial_ci(n_extreme, result.n_samples, config_.ci_confidence);
            result.ci_lower = ci_lower;
            result.ci_upper = ci_upper;

            // Stop if CI completely above or below alpha
            if (ci_upper < config_.alpha || ci_lower > config_.alpha) {
                result.early_stopped = true;
                break;
            }
        }
    }

    // Final p-value estimate
    result.p_value = static_cast<double>(n_extreme) / static_cast<double>(result.n_samples);
    auto [ci_lower, ci_upper] = binomial_ci(n_extreme, result.n_samples, config_.ci_confidence);
    result.ci_lower = ci_lower;
    result.ci_upper = ci_upper;

    return result;
}

// ============================================================================
// ContingencyTableRxC methods
// ============================================================================

int ContingencyTableRxC::row_total(int row) const {
    int sum = 0;
    for (int j = 0; j < n_cols; ++j) {
        sum += cells[row * n_cols + j];
    }
    return sum;
}

int ContingencyTableRxC::col_total(int col) const {
    int sum = 0;
    for (int i = 0; i < n_rows; ++i) {
        sum += cells[i * n_cols + col];
    }
    return sum;
}

int ContingencyTableRxC::grand_total() const {
    return std::accumulate(cells.begin(), cells.end(), 0);
}

// ============================================================================
// FullLabel methods
// ============================================================================

uint16_t FullLabel::combo_code() const {
    uint16_t code = 0;

    // Bits 0-1: Allele
    switch (allele) {
        case AltSupport::ALT:
            code |= 0;
            break;
        case AltSupport::REF:
            code |= 1;
            break;
        default:
            code |= 2;
            break;
    }

    // Bits 2-4: HP
    uint16_t hp_code = 5;  // Default: unphase
    if (hp_tag == "1") {
        hp_code = 0;
    } else if (hp_tag == "2") {
        hp_code = 1;
    } else if (hp_tag == "1-1") {
        hp_code = 2;
    } else if (hp_tag == "2-1") {
        hp_code = 3;
    } else if (hp_tag == "3") {
        hp_code = 4;
    }
    code |= (hp_code << 2);

    // Bit 5: Strand
    if (strand == Strand::REVERSE) {
        code |= (1 << 5);
    }

    // Bit 6: Sample
    if (!is_tumor) {
        code |= (1 << 6);
    }

    return code;
}

TestLabel FullLabel::get_allele_label() const {
    switch (allele) {
        case AltSupport::ALT:
            return TestLabel::ALLELE_ALT;
        case AltSupport::REF:
            return TestLabel::ALLELE_REF;
        default:
            return TestLabel::ALLELE_UNKNOWN;
    }
}

TestLabel FullLabel::get_hp_label() const {
    if (hp_tag == "1") {
        return TestLabel::HP_1;
    } else if (hp_tag == "2") {
        return TestLabel::HP_2;
    }
    return TestLabel::HP_OTHER;
}

TestLabel FullLabel::get_sample_label() const {
    return is_tumor ? TestLabel::SAMPLE_TUMOR : TestLabel::SAMPLE_NORMAL;
}

}  // namespace InterSubMod
