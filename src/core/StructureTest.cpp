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

    // Accumulate per-group within-group and cross-group squared distances in a single pass.
    // Direct formula: ss_between = cross_raw/n + adj_term
    //   where adj_term = Σ_g within_g * (1/n - 1/n_g)
    // This is algebraically equivalent to ss_total - ss_within but avoids
    // catastrophic cancellation when groups have similar spread.
    std::map<int, double> within_raw_per_group;
    for (const auto& [gid, members] : groups) within_raw_per_group[gid] = 0.0;

    double cross_raw = 0.0;
    double total_raw = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d2 = dist_matrix(i, j) * dist_matrix(i, j);
            total_raw += d2;
            if (group_labels[i] == group_labels[j]) {
                within_raw_per_group[group_labels[i]] += d2;
            } else {
                cross_raw += d2;
            }
        }
    }

    ss_total = total_raw / static_cast<double>(n);

    ss_within = 0.0;
    double adj_term = 0.0;
    for (const auto& [gid, members] : groups) {
        int n_g = static_cast<int>(members.size());
        if (n_g <= 0) continue;
        double wg = within_raw_per_group.at(gid);
        ss_within += wg / static_cast<double>(n_g);
        adj_term += wg * (1.0 / static_cast<double>(n) - 1.0 / static_cast<double>(n_g));
    }

    ss_between = cross_raw / static_cast<double>(n) + adj_term;
    if (ss_between < 0.0) ss_between = 0.0;  // clamp floating-point noise to zero
}

double StructureTest::compute_pseudo_f(double ss_between, double ss_within, int n, int k) {
    if (k <= 1 || n <= k) return 0.0;
    if (ss_within <= 0.0) {
        // Perfect separation: no within-group variation.
        // Return large sentinel so permutation test correctly rejects null hypothesis.
        return (ss_between > 0.0) ? 1e9 : 0.0;
    }
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

namespace {
// Regularized incomplete beta function I_x(a, b) (Numerical Recipes betacf/betai).
// Used to evaluate the F-distribution upper-tail p-value analytically (replaces a crude
// 3-bucket F-threshold lookup), so PERMDISP can screen dispersion-driven PERMANOVA hits.
double betacf(double a, double b, double x) {
    const int kMaxIt = 200;
    const double kEps = 3.0e-12;
    const double kFpMin = 1.0e-300;
    const double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < kFpMin) d = kFpMin;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= kMaxIt; ++m) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < kFpMin) d = kFpMin;
        c = 1.0 + aa / c;
        if (std::fabs(c) < kFpMin) c = kFpMin;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < kFpMin) d = kFpMin;
        c = 1.0 + aa / c;
        if (std::fabs(c) < kFpMin) c = kFpMin;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < kEps) break;
    }
    return h;
}

double betai(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double ln_beta = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b);
    const double bt = std::exp(ln_beta + a * std::log(x) + b * std::log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * betacf(a, b, x) / a;
    }
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

// Upper-tail p-value of the F-distribution: P(F_{d1,d2} >= f).
double f_dist_sf(double f, double d1, double d2) {
    if (f <= 0.0 || d1 <= 0.0 || d2 <= 0.0) return 1.0;
    const double x = d2 / (d2 + d1 * f);
    return betai(d2 / 2.0, d1 / 2.0, x);
}
}  // namespace

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
    // Analytic F-distribution upper-tail p: df1 = (#groups with n_k>1) - 1, df2 = N_obs - #groups.
    // all_distances/all_group_labels already exclude singleton groups (n_k<=1), so they define the
    // ANOVA design. Replaces the prior 3-bucket F-threshold lookup with a real p-value.
    const int n_obs = static_cast<int>(all_distances.size());
    const std::set<int> distinct_groups(all_group_labels.begin(), all_group_labels.end());
    const int k_groups = static_cast<int>(distinct_groups.size());
    const int df1 = k_groups - 1;
    const int df2 = n_obs - k_groups;
    result.anova_p = (df1 >= 1 && df2 >= 1) ? f_dist_sf(result.anova_f, df1, df2) : 1.0;
    result.warning = (result.anova_p < config_.dispersion_alpha);

    return result;
}

}  // namespace InterSubMod
