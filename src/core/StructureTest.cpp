/**
 * @file StructureTest.cpp
 * @brief Implementation of PERMANOVA and dispersion tests
 */

#include "core/StructureTest.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>

namespace InterSubMod {

StructureTest::StructureTest(const StructureTestConfig& config) : config_(config) {
    if (config_.seed == 0) {
        std::random_device rd;
        rng_.seed(rd());
    } else {
        rng_.seed(config_.seed);
    }
}

void StructureTest::set_seed(uint64_t seed) {
    config_.seed = seed;
    rng_.seed(seed);
}

bool StructureTest::filter_reads_for_complete_matrix(const Eigen::MatrixXd& dist_matrix,
                                                     std::vector<int>& valid_read_indices) {
    int n = dist_matrix.rows();
    if (n == 0) {
        valid_read_indices.clear();
        return false;
    }

    std::vector<bool> is_valid(n, true);
    int n_valid = n;

    while (true) {
        bool has_invalid = false;
        for (int i = 0; i < n && !has_invalid; ++i) {
            if (!is_valid[i]) continue;
            for (int j = i + 1; j < n; ++j) {
                if (!is_valid[j]) continue;
                if (std::isnan(dist_matrix(i, j))) {
                    has_invalid = true;
                    break;
                }
            }
        }

        if (!has_invalid) break;

        if (n_valid <= config_.min_reads_for_permanova) {
            valid_read_indices.clear();
            return false;
        }

        int worst_read = -1;
        int min_degree = n_valid;

        for (int i = 0; i < n; ++i) {
            if (!is_valid[i]) continue;
            int degree = 0;
            for (int j = 0; j < n; ++j) {
                if (i == j || !is_valid[j]) continue;
                if (!std::isnan(dist_matrix(i, j))) ++degree;
            }
            if (degree < min_degree) {
                min_degree = degree;
                worst_read = i;
            }
        }

        if (worst_read >= 0) {
            is_valid[worst_read] = false;
            --n_valid;
        } else {
            break;
        }
    }

    valid_read_indices.clear();
    for (int i = 0; i < n; ++i) {
        if (is_valid[i]) valid_read_indices.push_back(i);
    }

    return static_cast<int>(valid_read_indices.size()) >= config_.min_reads_for_permanova;
}

void StructureTest::compute_ss(const Eigen::MatrixXd& dist_matrix, const std::vector<int>& group_labels,
                               double& ss_total, double& ss_within, double& ss_between) {
    int n = dist_matrix.rows();

    std::map<int, std::vector<int>> groups;
    for (int i = 0; i < n; ++i) {
        groups[group_labels[i]].push_back(i);
    }

    ss_total = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d = dist_matrix(i, j);
            ss_total += d * d;
        }
    }
    ss_total /= static_cast<double>(n);

    ss_within = 0.0;
    for (const auto& [group_id, members] : groups) {
        int n_k = static_cast<int>(members.size());
        if (n_k <= 1) continue;
        double group_ss = 0.0;
        for (size_t i = 0; i < members.size(); ++i) {
            for (size_t j = i + 1; j < members.size(); ++j) {
                double d = dist_matrix(members[i], members[j]);
                group_ss += d * d;
            }
        }
        ss_within += group_ss / static_cast<double>(n_k);
    }

    ss_between = ss_total - ss_within;
}

double StructureTest::compute_pseudo_f(double ss_between, double ss_within, int n, int k) {
    if (k <= 1 || n <= k || ss_within <= 0.0) return 0.0;
    double df_between = static_cast<double>(k - 1);
    double df_within = static_cast<double>(n - k);
    return (ss_between / df_between) / (ss_within / df_within);
}

PermanovaResult StructureTest::run_permanova(const Eigen::MatrixXd& dist_matrix, const std::vector<int>& group_labels) {
    PermanovaResult result;
    result.n_permutations = config_.n_permutations;

    int n = dist_matrix.rows();
    if (n < config_.min_reads_for_permanova || static_cast<int>(group_labels.size()) != n) {
        result.valid = false;
        result.invalid_reason = "insufficient_reads";
        return result;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (std::isnan(dist_matrix(i, j))) {
                result.valid = false;
                result.invalid_reason = "incomplete_distance_matrix";
                return result;
            }
        }
    }

    std::set<int> unique_groups(group_labels.begin(), group_labels.end());
    int k = static_cast<int>(unique_groups.size());

    if (k < 2) {
        result.valid = false;
        result.invalid_reason = "insufficient_groups";
        return result;
    }

    double ss_total, ss_within, ss_between;
    compute_ss(dist_matrix, group_labels, ss_total, ss_within, ss_between);
    result.pseudo_f = compute_pseudo_f(ss_between, ss_within, n, k);

    std::vector<int> permuted_labels = group_labels;
    int n_extreme = 1;

    for (int perm = 0; perm < config_.n_permutations; ++perm) {
        std::shuffle(permuted_labels.begin(), permuted_labels.end(), rng_);
        double perm_ss_total, perm_ss_within, perm_ss_between;
        compute_ss(dist_matrix, permuted_labels, perm_ss_total, perm_ss_within, perm_ss_between);
        double perm_f = compute_pseudo_f(perm_ss_between, perm_ss_within, n, k);
        if (perm_f >= result.pseudo_f) ++n_extreme;
    }

    result.p_value = static_cast<double>(n_extreme) / static_cast<double>(config_.n_permutations + 1);
    result.valid = true;
    return result;
}

std::vector<double> StructureTest::compute_mean_distances_to_centroid(const Eigen::MatrixXd& dist_matrix,
                                                                      const std::vector<int>& group_labels) {
    int n = dist_matrix.rows();
    std::map<int, std::vector<int>> groups;
    for (int i = 0; i < n; ++i) {
        groups[group_labels[i]].push_back(i);
    }

    std::vector<double> mean_distances;
    for (const auto& [group_id, members] : groups) {
        int n_k = static_cast<int>(members.size());
        if (n_k <= 1) {
            mean_distances.push_back(0.0);
            continue;
        }
        double sum_sq_all = 0.0;
        for (size_t i = 0; i < members.size(); ++i) {
            for (size_t j = i + 1; j < members.size(); ++j) {
                double d = dist_matrix(members[i], members[j]);
                sum_sq_all += d * d;
            }
        }
        sum_sq_all *= 2.0;
        double mean_d = std::sqrt(sum_sq_all / (2.0 * n_k * n_k));
        mean_distances.push_back(mean_d);
    }
    return mean_distances;
}

double StructureTest::anova_f_test(const std::vector<double>& values, const std::vector<int>& group_labels) {
    if (values.empty() || values.size() != group_labels.size()) return 0.0;

    int n = static_cast<int>(values.size());
    std::map<int, std::vector<double>> groups;
    for (int i = 0; i < n; ++i) {
        groups[group_labels[i]].push_back(values[i]);
    }

    int k = static_cast<int>(groups.size());
    if (k < 2) return 0.0;

    double grand_mean = std::accumulate(values.begin(), values.end(), 0.0) / n;

    double ss_between = 0.0;
    for (const auto& [gid, vals] : groups) {
        double gm = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
        ss_between += vals.size() * (gm - grand_mean) * (gm - grand_mean);
    }

    double ss_within = 0.0;
    for (const auto& [gid, vals] : groups) {
        double gm = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
        for (double v : vals) ss_within += (v - gm) * (v - gm);
    }

    if (ss_within <= 0.0 || n <= k) return 0.0;
    return (ss_between / (k - 1)) / (ss_within / (n - k));
}

DispersionResult StructureTest::check_dispersion(const Eigen::MatrixXd& dist_matrix,
                                                 const std::vector<int>& group_labels) {
    DispersionResult result;
    int n = dist_matrix.rows();
    if (n < 2 || static_cast<int>(group_labels.size()) != n) return result;

    result.mean_distances_to_centroid = compute_mean_distances_to_centroid(dist_matrix, group_labels);

    std::map<int, std::vector<int>> groups;
    for (int i = 0; i < n; ++i) groups[group_labels[i]].push_back(i);

    std::vector<double> all_distances;
    std::vector<int> all_group_labels;

    for (const auto& [gid, members] : groups) {
        int n_k = static_cast<int>(members.size());
        if (n_k <= 1) continue;
        for (int i : members) {
            double sum_sq = 0.0;
            for (int j : members) {
                if (i == j) continue;
                double d = dist_matrix(i, j);
                sum_sq += d * d;
            }
            all_distances.push_back(std::sqrt(sum_sq / (2.0 * n_k)));
            all_group_labels.push_back(gid);
        }
    }

    result.anova_f = anova_f_test(all_distances, all_group_labels);
    result.anova_p = result.anova_f > 4.0 ? 0.01 : (result.anova_f > 2.5 ? 0.05 : 0.1);
    result.warning = (result.anova_p < config_.dispersion_alpha);

    return result;
}

}  // namespace InterSubMod
