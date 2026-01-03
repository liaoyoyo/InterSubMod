/**
 * @file LabelTest.cpp
 * @brief Implementation of Label-based significance testing
 *
 * Implements the Label-First verification path:
 * 1. Converts biological labels (HP/Allele) to binary group assignments
 * 2. Computes within-group and between-group mean distances
 * 3. Calculates Delta = d_between - d_within
 * 4. Runs permutation test to assess statistical significance
 */

#include "core/LabelTest.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>

namespace InterSubMod {

// ============================================================================
// Constructor
// ============================================================================

LabelTest::LabelTest(const LabelTestConfig& config) : config_(config) {
    if (config_.seed == 0) {
        std::random_device rd;
        rng_.seed(rd());
    } else {
        rng_.seed(config_.seed);
    }
}

void LabelTest::set_seed(uint64_t seed) {
    config_.seed = seed;
    rng_.seed(seed);
}

// ============================================================================
// Main Entry Point
// ============================================================================

LabelTestResult LabelTest::test_all(const Eigen::MatrixXd& dist_matrix,
                                    const std::vector<FullLabel>& full_labels) {
    LabelTestResult result;

    // Basic validation
    int n = dist_matrix.rows();
    if (n < config_.min_total_reads || static_cast<int>(full_labels.size()) != n) {
        result.valid = false;
        result.invalid_reason = "insufficient_reads";
        return result;
    }

    // Test HP dimension using multi-group PERMANOVA (NEW)
    std::vector<int> hp_multi_labels = hp_to_multigroup_labels(full_labels);
    test_hp_multigroup_permanova(dist_matrix, hp_multi_labels, result);

    // Also keep binary HP test for backward compatibility (deprecated)
    std::vector<int> hp_labels = hp_to_binary_labels(full_labels);
    result.hp_result = test_binary_groups(dist_matrix, hp_labels);

    // Test Allele dimension (ALT vs REF) - remains binary
    std::vector<int> allele_labels = allele_to_binary_labels(full_labels);
    result.allele_result = test_binary_groups(dist_matrix, allele_labels);

    // Determine dominant dimension (most significant)
    // Now using HP PERMANOVA p-value instead of binary Delta
    result.any_significant = result.hp_permanova_sig || result.allele_result.significant;

    // Pick dominant dimension based on lowest p-value
    bool hp_valid = (result.hp_n_groups >= 2);
    bool allele_valid = result.allele_result.valid;

    if (hp_valid && allele_valid) {
        // Both valid - pick the one with lower p-value
        if (result.hp_permanova_p <= result.allele_result.p_value) {
            result.dominant_dimension = "hp";
            result.best_delta = result.hp_permanova_f;  // Use pseudo-F as "delta" proxy
            result.best_p_value = result.hp_permanova_p;
        } else {
            result.dominant_dimension = "allele";
            result.best_delta = result.allele_result.delta;
            result.best_p_value = result.allele_result.p_value;
        }
    } else if (hp_valid) {
        result.dominant_dimension = "hp";
        result.best_delta = result.hp_permanova_f;
        result.best_p_value = result.hp_permanova_p;
    } else if (allele_valid) {
        result.dominant_dimension = "allele";
        result.best_delta = result.allele_result.delta;
        result.best_p_value = result.allele_result.p_value;
    } else {
        result.dominant_dimension = "none";
        result.valid = false;
        result.invalid_reason = "no_valid_label_dimension";
    }

    return result;
}

// ============================================================================
// Binary Group Testing
// ============================================================================

LabelDeltaResult LabelTest::test_binary_groups(const Eigen::MatrixXd& dist_matrix,
                                               const std::vector<int>& group_labels) {
    // First compute delta without permutation
    LabelDeltaResult result = compute_delta(dist_matrix, group_labels);

    if (!result.valid) {
        return result;
    }

    // Run permutation test if delta is positive
    if (result.delta > 0) {
        result.p_value = permutation_test(dist_matrix, group_labels, result.delta);
        result.n_permutations = config_.n_permutations;
        result.significant = (result.p_value <= config_.alpha);
    } else {
        // Negative delta means within-group distance > between-group distance
        // This is the opposite of expected biological signal
        result.p_value = 1.0;
        result.n_permutations = 0;
        result.significant = false;
    }

    return result;
}

LabelDeltaResult LabelTest::compute_delta(const Eigen::MatrixXd& dist_matrix,
                                          const std::vector<int>& group_labels) {
    LabelDeltaResult result;
    int n = dist_matrix.rows();

    if (n == 0 || static_cast<int>(group_labels.size()) != n) {
        result.valid = false;
        result.invalid_reason = "matrix_label_size_mismatch";
        return result;
    }

    // Count group sizes
    int count_0 = 0, count_1 = 0;
    for (int label : group_labels) {
        if (label == 0) ++count_0;
        else if (label == 1) ++count_1;
    }

    result.n_group_a = count_0;
    result.n_group_b = count_1;

    // Check minimum requirements
    if (count_0 < config_.min_reads_per_group || count_1 < config_.min_reads_per_group) {
        result.valid = false;
        result.invalid_reason = "insufficient_reads_per_group";
        return result;
    }

    // Compute distances
    compute_group_distances(dist_matrix, group_labels,
                           result.within_mean, result.between_mean,
                           result.n_pairs_within, result.n_pairs_between);

    // Check if we have enough valid pairs
    if (result.n_pairs_within == 0 || result.n_pairs_between == 0) {
        result.valid = false;
        result.invalid_reason = "insufficient_valid_pairs";
        return result;
    }

    // Calculate delta
    result.delta = result.between_mean - result.within_mean;
    result.valid = true;

    return result;
}

// ============================================================================
// Label Conversion
// ============================================================================

std::vector<int> LabelTest::hp_to_binary_labels(const std::vector<FullLabel>& full_labels) {
    std::vector<int> binary_labels(full_labels.size(), -1);  // -1 = excluded

    for (size_t i = 0; i < full_labels.size(); ++i) {
        const std::string& hp = full_labels[i].hp_tag;

        // HP1 variants -> group 0 (includes HP1 from any sample: 1, 1-1)
        if (hp == "1" || hp == "HP1" || hp == "1-1") {
            binary_labels[i] = 0;
        }
        // HP2 variants -> group 1 (includes HP2 from any sample: 2, 2-1)
        else if (hp == "2" || hp == "HP2" || hp == "2-1") {
            binary_labels[i] = 1;
        }
        // All others (0, 3, empty/unphased) -> excluded (-1)
        // HP0 = unphased, HP3 = conflict, these are not used for HP1 vs HP2 comparison
    }

    return binary_labels;
}

std::vector<int> LabelTest::allele_to_binary_labels(const std::vector<FullLabel>& full_labels) {
    std::vector<int> binary_labels(full_labels.size(), -1);  // -1 = excluded

    for (size_t i = 0; i < full_labels.size(); ++i) {
        switch (full_labels[i].allele) {
            case AltSupport::ALT:
                binary_labels[i] = 0;
                break;
            case AltSupport::REF:
                binary_labels[i] = 1;
                break;
            case AltSupport::UNKNOWN:
            default:
                binary_labels[i] = -1;  // Excluded
                break;
        }
    }

    return binary_labels;
}

// ============================================================================
// Distance Computation
// ============================================================================

void LabelTest::compute_group_distances(const Eigen::MatrixXd& dist_matrix,
                                        const std::vector<int>& group_labels,
                                        double& within_mean, double& between_mean,
                                        int& n_pairs_within, int& n_pairs_between) {
    int n = dist_matrix.rows();
    double sum_within = 0.0;
    double sum_between = 0.0;
    n_pairs_within = 0;
    n_pairs_between = 0;

    // Iterate over all unique pairs (i < j)
    for (int i = 0; i < n; ++i) {
        // Skip reads with invalid labels (-1)
        if (group_labels[i] < 0) continue;

        for (int j = i + 1; j < n; ++j) {
            // Skip reads with invalid labels
            if (group_labels[j] < 0) continue;

            double dist = dist_matrix(i, j);

            // Skip NaN distances
            if (std::isnan(dist)) continue;

            if (group_labels[i] == group_labels[j]) {
                // Same group -> within-group distance
                sum_within += dist;
                ++n_pairs_within;
            } else {
                // Different groups -> between-group distance
                sum_between += dist;
                ++n_pairs_between;
            }
        }
    }

    // Compute means
    within_mean = (n_pairs_within > 0) ? (sum_within / n_pairs_within) : 0.0;
    between_mean = (n_pairs_between > 0) ? (sum_between / n_pairs_between) : 0.0;
}

// ============================================================================
// Permutation Test
// ============================================================================

double LabelTest::permutation_test(const Eigen::MatrixXd& dist_matrix,
                                   const std::vector<int>& group_labels,
                                   double observed_delta) {
    int n_extreme = 1;  // Include observed value

    // Create a mutable copy of labels for shuffling
    std::vector<int> permuted_labels = group_labels;

    // Get indices of valid labels (not -1) for shuffling
    std::vector<size_t> valid_indices;
    for (size_t i = 0; i < group_labels.size(); ++i) {
        if (group_labels[i] >= 0) {
            valid_indices.push_back(i);
        }
    }

    // Extract valid labels for shuffling
    std::vector<int> valid_label_values;
    for (size_t idx : valid_indices) {
        valid_label_values.push_back(group_labels[idx]);
    }

    for (int perm = 0; perm < config_.n_permutations; ++perm) {
        // Shuffle only the valid label values
        std::shuffle(valid_label_values.begin(), valid_label_values.end(), rng_);

        // Copy back to permuted_labels
        for (size_t i = 0; i < valid_indices.size(); ++i) {
            permuted_labels[valid_indices[i]] = valid_label_values[i];
        }

        // Compute permuted delta
        double within_mean, between_mean;
        int n_within, n_between;
        compute_group_distances(dist_matrix, permuted_labels,
                               within_mean, between_mean,
                               n_within, n_between);

        double perm_delta = between_mean - within_mean;

        // Count if permuted >= observed
        if (perm_delta >= observed_delta) {
            ++n_extreme;
        }
    }

    // P-value = (n_extreme) / (n_permutations + 1)
    return static_cast<double>(n_extreme) / static_cast<double>(config_.n_permutations + 1);
}

// ============================================================================
// HP Multi-group PERMANOVA
// ============================================================================

std::vector<int> LabelTest::hp_to_multigroup_labels(const std::vector<FullLabel>& full_labels) {
    // First pass: assign raw group labels
    // Mapping: HP1->0, HP2->1, HP1-1->2, HP2-1->3, HP3->4, HP0/empty->5
    std::vector<int> raw_labels(full_labels.size(), -1);
    std::map<int, int> group_counts;

    for (size_t i = 0; i < full_labels.size(); ++i) {
        const std::string& hp = full_labels[i].hp_tag;
        int group = -1;

        if (hp == "1" || hp == "HP1") {
            group = 0;  // HP1
        } else if (hp == "2" || hp == "HP2") {
            group = 1;  // HP2
        } else if (hp == "1-1") {
            group = 2;  // HP1-1
        } else if (hp == "2-1") {
            group = 3;  // HP2-1
        } else if (hp == "3" || hp == "HP3") {
            group = 4;  // HP3
        } else if (hp.empty() || hp == "0" || hp == "HP0" || hp == "unphased") {
            group = 5;  // HP0/unphased
        }

        raw_labels[i] = group;
        if (group >= 0) {
            group_counts[group]++;
        }
    }

    // Second pass: exclude groups with too few reads
    std::vector<int> final_labels(full_labels.size(), -1);
    for (size_t i = 0; i < full_labels.size(); ++i) {
        int group = raw_labels[i];
        if (group >= 0 && group_counts[group] >= config_.min_reads_per_group) {
            final_labels[i] = group;
        }
    }

    return final_labels;
}

void LabelTest::test_hp_multigroup_permanova(const Eigen::MatrixXd& dist_matrix,
                                             const std::vector<int>& hp_labels,
                                             LabelTestResult& result) {
    int n = dist_matrix.rows();

    // Filter out excluded reads (-1) and build valid subset
    std::set<int> unique_groups;
    std::vector<int> valid_indices;
    for (int i = 0; i < n; ++i) {
        if (hp_labels[i] >= 0) {
            unique_groups.insert(hp_labels[i]);
            valid_indices.push_back(i);
        }
    }

    int n_valid = static_cast<int>(valid_indices.size());
    int n_groups = static_cast<int>(unique_groups.size());
    result.hp_n_groups = n_groups;

    // Need at least 2 groups for PERMANOVA
    if (n_groups < 2 || n_valid < config_.min_total_reads) {
        result.hp_permanova_f = 0.0;
        result.hp_permanova_p = 1.0;
        result.hp_permanova_sig = false;
        return;
    }

    // Build subset distance matrix and labels
    Eigen::MatrixXd sub_dist(n_valid, n_valid);
    std::vector<int> sub_labels(n_valid);
    for (int i = 0; i < n_valid; ++i) {
        sub_labels[i] = hp_labels[valid_indices[i]];
        for (int j = 0; j < n_valid; ++j) {
            double d = dist_matrix(valid_indices[i], valid_indices[j]);
            sub_dist(i, j) = std::isnan(d) ? 0.0 : d;  // Replace NaN with 0
        }
    }

    // Compute SS (sum of squares) using PERMANOVA formulas
    auto compute_ss = [&](const std::vector<int>& labels) {
        std::map<int, std::vector<int>> groups;
        for (int i = 0; i < n_valid; ++i) {
            groups[labels[i]].push_back(i);
        }

        double ss_total = 0.0;
        for (int i = 0; i < n_valid; ++i) {
            for (int j = i + 1; j < n_valid; ++j) {
                double d = sub_dist(i, j);
                ss_total += d * d;
            }
        }
        ss_total /= static_cast<double>(n_valid);

        double ss_within = 0.0;
        for (const auto& [group_id, members] : groups) {
            int n_k = static_cast<int>(members.size());
            if (n_k <= 1) continue;
            double group_ss = 0.0;
            for (size_t i = 0; i < members.size(); ++i) {
                for (size_t j = i + 1; j < members.size(); ++j) {
                    double d = sub_dist(members[i], members[j]);
                    group_ss += d * d;
                }
            }
            ss_within += group_ss / static_cast<double>(n_k);
        }

        double ss_between = ss_total - ss_within;
        return std::make_pair(ss_between, ss_within);
    };

    // Compute pseudo-F
    auto compute_pseudo_f = [&](double ss_between, double ss_within, int n, int k) {
        if (k <= 1 || n <= k || ss_within <= 0.0) return 0.0;
        double df_between = static_cast<double>(k - 1);
        double df_within = static_cast<double>(n - k);
        return (ss_between / df_between) / (ss_within / df_within);
    };

    // Observed F
    auto [obs_ss_between, obs_ss_within] = compute_ss(sub_labels);
    double obs_f = compute_pseudo_f(obs_ss_between, obs_ss_within, n_valid, n_groups);
    result.hp_permanova_f = obs_f;

    // Permutation test
    std::vector<int> permuted_labels = sub_labels;
    int n_extreme = 1;  // Include observed

    for (int perm = 0; perm < config_.n_permutations; ++perm) {
        std::shuffle(permuted_labels.begin(), permuted_labels.end(), rng_);
        auto [perm_ss_between, perm_ss_within] = compute_ss(permuted_labels);
        double perm_f = compute_pseudo_f(perm_ss_between, perm_ss_within, n_valid, n_groups);
        if (perm_f >= obs_f) {
            ++n_extreme;
        }
    }

    result.hp_permanova_p = static_cast<double>(n_extreme) / static_cast<double>(config_.n_permutations + 1);
    result.hp_permanova_sig = (result.hp_permanova_p <= config_.alpha);
}

}  // namespace InterSubMod
