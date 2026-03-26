/**
 * @file Bootstrap.cpp
 * @brief Implementation of bootstrap stability test
 */

#include "core/Bootstrap.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>

namespace InterSubMod {

Bootstrap::Bootstrap(const BootstrapConfig& config) : config_(config) {
    if (config_.seed == 0) {
        std::random_device rd;
        rng_.seed(rd());
    } else {
        rng_.seed(config_.seed);
    }
}

void Bootstrap::set_seed(uint64_t seed) {
    config_.seed = seed;
    rng_.seed(seed);
}

double Bootstrap::adjusted_rand_index(const std::vector<int>& labels1, const std::vector<int>& labels2) {
    if (labels1.size() != labels2.size() || labels1.empty()) {
        return 0.0;
    }

    int n = static_cast<int>(labels1.size());

    // Build contingency table
    std::map<std::pair<int, int>, int> contingency;
    std::map<int, int> a_sum, b_sum;

    for (int i = 0; i < n; ++i) {
        contingency[{labels1[i], labels2[i]}]++;
        a_sum[labels1[i]]++;
        b_sum[labels2[i]]++;
    }

    // Compute ARI components
    double sum_nij_choose_2 = 0.0;
    for (const auto& [key, nij] : contingency) {
        if (nij >= 2) {
            sum_nij_choose_2 += static_cast<double>(nij) * (nij - 1) / 2.0;
        }
    }

    double sum_ai_choose_2 = 0.0;
    for (const auto& [key, ai] : a_sum) {
        if (ai >= 2) {
            sum_ai_choose_2 += static_cast<double>(ai) * (ai - 1) / 2.0;
        }
    }

    double sum_bj_choose_2 = 0.0;
    for (const auto& [key, bj] : b_sum) {
        if (bj >= 2) {
            sum_bj_choose_2 += static_cast<double>(bj) * (bj - 1) / 2.0;
        }
    }

    double n_choose_2 = static_cast<double>(n) * (n - 1) / 2.0;
    if (n_choose_2 == 0) return 0.0;

    double expected = (sum_ai_choose_2 * sum_bj_choose_2) / n_choose_2;
    double max_index = (sum_ai_choose_2 + sum_bj_choose_2) / 2.0;

    double denominator = max_index - expected;
    if (std::abs(denominator) < 1e-10) {
        return 1.0;  // Perfect agreement
    }

    return (sum_nij_choose_2 - expected) / denominator;
}

Eigen::MatrixXd Bootstrap::resample_columns(const Eigen::MatrixXd& matrix) {
    int n_rows = matrix.rows();
    int n_cols = matrix.cols();

    if (n_cols == 0) return matrix;

    std::uniform_int_distribution<int> dist(0, n_cols - 1);

    Eigen::MatrixXd resampled(n_rows, n_cols);
    for (int j = 0; j < n_cols; ++j) {
        int sampled_col = dist(rng_);
        resampled.col(j) = matrix.col(sampled_col);
    }

    return resampled;
}

Eigen::MatrixXd Bootstrap::compute_l1_distance(const Eigen::MatrixXd& matrix) {
    int n = matrix.rows();
    Eigen::MatrixXd dist_matrix = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int count = 0;
            double sum = 0.0;

            for (int k = 0; k < matrix.cols(); ++k) {
                double vi = matrix(i, k);
                double vj = matrix(j, k);

                if (!std::isnan(vi) && !std::isnan(vj)) {
                    sum += std::abs(vi - vj);
                    ++count;
                }
            }

            double d = (count > 0) ? (sum / count) : std::nan("");
            dist_matrix(i, j) = d;
            dist_matrix(j, i) = d;
        }
    }

    return dist_matrix;
}

BootstrapResult Bootstrap::compute(const Eigen::MatrixXd& methylation_matrix, const std::vector<int>& original_labels,
                                   std::function<std::vector<int>(const Eigen::MatrixXd&)> cluster_func) {
    BootstrapResult result;
    result.n_iterations = 0;
    result.early_stopped = false;

    int n = methylation_matrix.rows();
    if (n < 2 || static_cast<int>(original_labels.size()) != n) {
        result.stable = false;
        return result;
    }

    int consecutive_low = 0;
    int consecutive_high = 0;

    for (int iter = 0; iter < config_.n_iterations; ++iter) {
        // Resample CpG sites
        Eigen::MatrixXd resampled = resample_columns(methylation_matrix);

        // Compute distance matrix
        Eigen::MatrixXd dist = compute_l1_distance(resampled);

        // Re-cluster
        std::vector<int> new_labels = cluster_func(dist);

        // Compute ARI
        double ari = adjusted_rand_index(original_labels, new_labels);
        result.ari_samples.push_back(ari);
        result.n_iterations = iter + 1;

        // Early stopping check
        if (config_.use_early_stopping) {
            if (ari < config_.early_stop_low_threshold) {
                consecutive_low++;
                consecutive_high = 0;
            } else if (ari > config_.early_stop_high_threshold) {
                consecutive_high++;
                consecutive_low = 0;
            } else {
                consecutive_low = 0;
                consecutive_high = 0;
            }

            if (consecutive_low >= config_.early_stop_consecutive ||
                consecutive_high >= config_.early_stop_consecutive) {
                result.early_stopped = true;
                break;
            }
        }
    }

    // Compute statistics
    if (!result.ari_samples.empty()) {
        result.mean_ari =
            std::accumulate(result.ari_samples.begin(), result.ari_samples.end(), 0.0) / result.ari_samples.size();

        double sum_sq = 0.0;
        for (double ari : result.ari_samples) {
            sum_sq += (ari - result.mean_ari) * (ari - result.mean_ari);
        }
        result.std_ari = std::sqrt(sum_sq / result.ari_samples.size());

        result.stable = (result.mean_ari >= 0.2);
    }

    return result;
}

}  // namespace InterSubMod
