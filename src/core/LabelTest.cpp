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
        result.hp_multistage.valid = false;
        result.hp_multistage.invalid_reason = "insufficient_reads";
        return result;
    }

    // ========================================================================
    // Multi-Stage HP Verification (NEW)
    // ========================================================================

    // Stage 1: HP Family Merged Test
    test_hp_merged(dist_matrix, full_labels, result.hp_multistage);

    // Stage 2: HP Fine-Grained Test (4-group PERMANOVA + pairwise distances)
    test_hp_fine_grained(dist_matrix, full_labels, result.hp_multistage);

    // Stage 4: Unassigned Affinity Test (HP3/HP0 → HP1/HP2)
    std::vector<int> merged_labels = hp_to_merged_labels(full_labels);
    test_unassigned_affinity(dist_matrix, full_labels, merged_labels, result.hp_multistage);

    // ========================================================================
    // Backward Compatibility: Populate deprecated fields
    // ========================================================================

    // Copy Stage 1 result to deprecated hp_result
    result.hp_result = result.hp_multistage.merged;

    // Copy Stage 2 result to deprecated PERMANOVA fields
    result.hp_permanova_f = result.hp_multistage.fine_f;
    result.hp_permanova_p = result.hp_multistage.fine_p;
    result.hp_permanova_sig = result.hp_multistage.fine_sig;
    result.hp_n_groups = result.hp_multistage.fine_n_groups;

    // ========================================================================
    // Stage 3: Allele Test (ALT vs REF)
    // ========================================================================

    std::vector<int> allele_labels = allele_to_binary_labels(full_labels);
    result.allele_result = test_binary_groups(dist_matrix, allele_labels);

    // ========================================================================
    // Sample ASM Test (Tumor vs Normal)
    // ========================================================================

    std::vector<int> sample_labels = sample_to_binary_labels(full_labels);
    result.sample_result = test_binary_groups(dist_matrix, sample_labels);

    // ========================================================================
    // Determine Dominant Dimension
    // ========================================================================

    // Use Stage 1 (merged) significance for HP dimension
    bool hp_sig = result.hp_multistage.merged.significant;
    bool allele_sig = result.allele_result.significant;
    result.any_significant = hp_sig || allele_sig;

    // Pick dominant dimension based on lowest p-value
    bool hp_valid = result.hp_multistage.merged.valid;
    bool allele_valid = result.allele_result.valid;

    double hp_p = result.hp_multistage.merged.p_value;
    double allele_p = result.allele_result.p_value;

    if (hp_valid && allele_valid) {
        // Both valid - pick the one with lower p-value
        if (hp_p <= allele_p) {
            result.dominant_dimension = "hp";
            result.best_delta = result.hp_multistage.merged.delta;
            result.best_p_value = hp_p;
        } else {
            result.dominant_dimension = "allele";
            result.best_delta = result.allele_result.delta;
            result.best_p_value = allele_p;
        }
    } else if (hp_valid) {
        result.dominant_dimension = "hp";
        result.best_delta = result.hp_multistage.merged.delta;
        result.best_p_value = hp_p;
    } else if (allele_valid) {
        result.dominant_dimension = "allele";
        result.best_delta = result.allele_result.delta;
        result.best_p_value = allele_p;
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

std::vector<int> LabelTest::hp_to_merged_labels(const std::vector<FullLabel>& full_labels) {
    // Stage 1 grouping:
    // HP1-family (HP1 + HP1-1) -> 0
    // HP2-family (HP2 + HP2-1) -> 1
    // HP3, HP0, empty -> -1 (excluded)
    std::vector<int> labels(full_labels.size(), -1);

    for (size_t i = 0; i < full_labels.size(); ++i) {
        const std::string& hp = full_labels[i].hp_tag;

        if (hp == "1" || hp == "HP1" || hp == "1-1") {
            labels[i] = 0;  // HP1-family
        } else if (hp == "2" || hp == "HP2" || hp == "2-1") {
            labels[i] = 1;  // HP2-family
        }
        // HP3, HP0, empty remain -1 (excluded)
    }

    return labels;
}

std::vector<int> LabelTest::hp_to_fine_labels(const std::vector<FullLabel>& full_labels) {
    // Stage 2 grouping (fine-grained, excludes HP3/HP0):
    // HP1 -> 0, HP1-1 -> 1, HP2 -> 2, HP2-1 -> 3
    // HP3, HP0, empty -> -1 (excluded)

    // First pass: assign raw labels and count
    std::vector<int> raw_labels(full_labels.size(), -1);
    std::map<int, int> group_counts;

    for (size_t i = 0; i < full_labels.size(); ++i) {
        const std::string& hp = full_labels[i].hp_tag;
        int group = -1;

        if (hp == "1" || hp == "HP1") {
            group = 0;  // HP1 (germline)
        } else if (hp == "1-1") {
            group = 1;  // HP1-1 (somatic from HP1)
        } else if (hp == "2" || hp == "HP2") {
            group = 2;  // HP2 (germline)
        } else if (hp == "2-1") {
            group = 3;  // HP2-1 (somatic from HP2)
        }
        // HP3, HP0, empty remain -1 (excluded)

        raw_labels[i] = group;
        if (group >= 0) {
            group_counts[group]++;
        }
    }

    // Second pass: exclude groups with insufficient reads
    std::vector<int> final_labels(full_labels.size(), -1);
    for (size_t i = 0; i < full_labels.size(); ++i) {
        int group = raw_labels[i];
        if (group >= 0 && group_counts[group] >= config_.min_reads_per_group) {
            final_labels[i] = group;
        }
    }

    return final_labels;
}

std::vector<int> LabelTest::lineage_to_labels(const std::vector<FullLabel>& full_labels, const std::string& axis,
                                              const std::string& min_status) {
    std::vector<int> labels(full_labels.size(), -1);  // -1 = excluded

    // Confidence gate. 'P' (partial) carries a subtree-scoped path and 'A' (abstain)
    // carries no topology at all, so admitting them changes what a "group" means.
    const auto admits = [&min_status](char status) {
        if (status == '\0') return false;
        if (min_status == "any") return status != 'A';
        if (min_status == "UM") return status == 'U' || status == 'M';
        return status == 'U';  // default: unique only
    };

    // Assign a dense group index per distinct value, in first-seen order so the
    // numbering is deterministic for a given read ordering.
    std::map<std::string, int> group_of;
    for (size_t i = 0; i < full_labels.size(); ++i) {
        const FullLabel& fl = full_labels[i];
        if (!admits(fl.lineage_status)) continue;

        const std::string* value = nullptr;
        if (axis == "lc") {
            value = &fl.lineage_component;
        } else if (axis == "lu") {
            value = &fl.lineage_block;
        } else if (axis == "lv") {
            value = &fl.lineage_path;
        } else {
            return labels;  // unknown axis: everything excluded
        }
        if (value->empty() || *value == ".") continue;

        auto it = group_of.find(*value);
        if (it == group_of.end()) {
            it = group_of.emplace(*value, static_cast<int>(group_of.size())).first;
        }
        labels[i] = it->second;
    }
    return labels;
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

std::vector<int> LabelTest::sample_to_binary_labels(const std::vector<FullLabel>& full_labels) {
    std::vector<int> labels(full_labels.size());
    for (size_t i = 0; i < full_labels.size(); ++i) {
        labels[i] = full_labels[i].is_tumor ? 0 : 1;  // Tumor=0, Normal=1
    }
    return labels;
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

// ============================================================================
// Multi-Stage HP Tests (NEW)
// ============================================================================

void LabelTest::test_hp_merged(const Eigen::MatrixXd& dist_matrix,
                               const std::vector<FullLabel>& full_labels,
                               MultiStageHPResult& result) {
    // Stage 1: HP Family Merged Test
    // Groups: (HP1 + HP1-1) vs (HP2 + HP2-1)

    std::vector<int> merged_labels = hp_to_merged_labels(full_labels);

    // Count family sizes
    result.n_hp1_family = 0;
    result.n_hp2_family = 0;
    for (int label : merged_labels) {
        if (label == 0) ++result.n_hp1_family;
        else if (label == 1) ++result.n_hp2_family;
    }

    // Run binary delta test
    result.merged = test_binary_groups(dist_matrix, merged_labels);
}

void LabelTest::test_hp_fine_grained(const Eigen::MatrixXd& dist_matrix,
                                     const std::vector<FullLabel>& full_labels,
                                     MultiStageHPResult& result) {
    // Stage 2: HP Fine-Grained Test
    // Groups: HP1(0), HP1-1(1), HP2(2), HP2-1(3)

    std::vector<int> fine_labels = hp_to_fine_labels(full_labels);
    int n = dist_matrix.rows();

    // Count groups
    std::set<int> unique_groups;
    std::vector<int> valid_indices;
    for (int i = 0; i < n; ++i) {
        if (fine_labels[i] >= 0) {
            unique_groups.insert(fine_labels[i]);
            valid_indices.push_back(i);
        }
    }

    // Reset group counts
    result.fine_group_counts = {0, 0, 0, 0};
    for (int label : fine_labels) {
        if (label >= 0 && label < 4) {
            result.fine_group_counts[label]++;
        }
    }
    result.fine_n_groups = static_cast<int>(unique_groups.size());

    int n_valid = static_cast<int>(valid_indices.size());

    // Need at least 2 groups for PERMANOVA
    if (result.fine_n_groups < 2 || n_valid < config_.min_total_reads) {
        result.fine_f = 0.0;
        result.fine_p = 1.0;
        result.fine_sig = false;
        return;
    }

    // Build subset distance matrix
    Eigen::MatrixXd sub_dist(n_valid, n_valid);
    std::vector<int> sub_labels(n_valid);
    for (int i = 0; i < n_valid; ++i) {
        sub_labels[i] = fine_labels[valid_indices[i]];
        for (int j = 0; j < n_valid; ++j) {
            double d = dist_matrix(valid_indices[i], valid_indices[j]);
            sub_dist(i, j) = std::isnan(d) ? 0.0 : d;
        }
    }

    // Compute SS using PERMANOVA formulas
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

        return std::make_pair(ss_total - ss_within, ss_within);
    };

    auto compute_pseudo_f = [](double ss_between, double ss_within, int n, int k) {
        if (k <= 1 || n <= k || ss_within <= 0.0) return 0.0;
        double df_between = static_cast<double>(k - 1);
        double df_within = static_cast<double>(n - k);
        return (ss_between / df_between) / (ss_within / df_within);
    };

    // Observed F
    auto [obs_ss_between, obs_ss_within] = compute_ss(sub_labels);
    result.fine_f = compute_pseudo_f(obs_ss_between, obs_ss_within, n_valid, result.fine_n_groups);

    // Permutation test
    std::vector<int> permuted_labels = sub_labels;
    int n_extreme = 1;

    for (int perm = 0; perm < config_.n_permutations; ++perm) {
        std::shuffle(permuted_labels.begin(), permuted_labels.end(), rng_);
        auto [perm_ss_between, perm_ss_within] = compute_ss(permuted_labels);
        double perm_f = compute_pseudo_f(perm_ss_between, perm_ss_within, n_valid, result.fine_n_groups);
        if (perm_f >= result.fine_f) {
            ++n_extreme;
        }
    }

    result.fine_p = static_cast<double>(n_extreme) / static_cast<double>(config_.n_permutations + 1);
    result.fine_sig = (result.fine_p <= config_.alpha);

    // Compute pairwise mean distances between groups
    // Group indices: HP1(0), HP1-1(1), HP2(2), HP2-1(3)
    std::map<int, std::vector<int>> group_members;
    for (int i = 0; i < n_valid; ++i) {
        group_members[sub_labels[i]].push_back(i);
    }

    auto compute_between_mean = [&](int g1, int g2) -> double {
        if (group_members.find(g1) == group_members.end() ||
            group_members.find(g2) == group_members.end()) {
            return NAN;
        }
        const auto& m1 = group_members[g1];
        const auto& m2 = group_members[g2];
        if (m1.empty() || m2.empty()) return NAN;

        double sum = 0.0;
        int count = 0;
        for (int i : m1) {
            for (int j : m2) {
                sum += sub_dist(i, j);
                ++count;
            }
        }
        return count > 0 ? sum / count : NAN;
    };

    result.fine_pairwise.d_hp1_hp1s = compute_between_mean(0, 1);
    result.fine_pairwise.d_hp1_hp2 = compute_between_mean(0, 2);
    result.fine_pairwise.d_hp1_hp2s = compute_between_mean(0, 3);
    result.fine_pairwise.d_hp1s_hp2 = compute_between_mean(1, 2);
    result.fine_pairwise.d_hp1s_hp2s = compute_between_mean(1, 3);
    result.fine_pairwise.d_hp2_hp2s = compute_between_mean(2, 3);
}

void LabelTest::test_unassigned_affinity(const Eigen::MatrixXd& dist_matrix,
                                         const std::vector<FullLabel>& full_labels,
                                         const std::vector<int>& merged_labels,
                                         MultiStageHPResult& result) {
    // Stage 4: Unassigned Affinity Test
    // Target: HP3, HP0 reads
    // Method: Compute mean distance to HP1-family vs HP2-family

    int n = dist_matrix.rows();

    // Find unassigned (HP3, HP0) and family indices
    std::vector<int> unassigned_indices;
    std::vector<int> hp1_family_indices;
    std::vector<int> hp2_family_indices;

    for (int i = 0; i < n; ++i) {
        const std::string& hp = full_labels[i].hp_tag;

        if (hp == "3" || hp == "HP3") {
            unassigned_indices.push_back(i);
            result.unassigned.n_hp3++;
        } else if (hp.empty() || hp == "0" || hp == "HP0" || hp == "unphased") {
            unassigned_indices.push_back(i);
            result.unassigned.n_hp0++;
        } else if (merged_labels[i] == 0) {
            hp1_family_indices.push_back(i);
        } else if (merged_labels[i] == 1) {
            hp2_family_indices.push_back(i);
        }
    }

    // Check validity
    int n_unassigned = static_cast<int>(unassigned_indices.size());
    if (n_unassigned == 0) {
        result.unassigned.valid = false;
        result.unassigned.invalid_reason = "no_unassigned_reads";
        return;
    }
    if (hp1_family_indices.empty() || hp2_family_indices.empty()) {
        result.unassigned.valid = false;
        result.unassigned.invalid_reason = "insufficient_family_reads";
        return;
    }

    // Compute mean distance from each unassigned read to each family
    auto compute_mean_dist_to_family = [&](int idx, const std::vector<int>& family) {
        double sum = 0.0;
        int count = 0;
        for (int f : family) {
            double d = dist_matrix(idx, f);
            if (!std::isnan(d)) {
                sum += d;
                ++count;
            }
        }
        return count > 0 ? sum / count : NAN;
    };

    double sum_d_hp1 = 0.0;
    double sum_d_hp2 = 0.0;
    int valid_unassigned = 0;

    for (int u : unassigned_indices) {
        double d_hp1 = compute_mean_dist_to_family(u, hp1_family_indices);
        double d_hp2 = compute_mean_dist_to_family(u, hp2_family_indices);

        if (!std::isnan(d_hp1) && !std::isnan(d_hp2)) {
            sum_d_hp1 += d_hp1;
            sum_d_hp2 += d_hp2;
            ++valid_unassigned;
        }
    }

    if (valid_unassigned == 0) {
        result.unassigned.valid = false;
        result.unassigned.invalid_reason = "no_valid_distance_pairs";
        return;
    }

    double mean_d_hp1 = sum_d_hp1 / valid_unassigned;
    double mean_d_hp2 = sum_d_hp2 / valid_unassigned;

    // Compute affinity score: (d_hp2 - d_hp1) / (d_hp1 + d_hp2)
    // Positive = closer to HP1, Negative = closer to HP2
    double denom = mean_d_hp1 + mean_d_hp2;
    result.unassigned.affinity_score = (denom > 0) ? (mean_d_hp2 - mean_d_hp1) / denom : 0.0;

    // Determine direction (threshold = 0.1)
    const double threshold = 0.1;
    if (result.unassigned.affinity_score > threshold) {
        result.unassigned.affinity_direction = "HP1";
    } else if (result.unassigned.affinity_score < -threshold) {
        result.unassigned.affinity_direction = "HP2";
    } else {
        result.unassigned.affinity_direction = "None";
    }

    // Permutation test for affinity score
    double observed_abs_score = std::abs(result.unassigned.affinity_score);
    int n_extreme = 1;

    // Pool HP1 and HP2 family indices for shuffling
    std::vector<int> all_family_indices;
    all_family_indices.insert(all_family_indices.end(), hp1_family_indices.begin(), hp1_family_indices.end());
    all_family_indices.insert(all_family_indices.end(), hp2_family_indices.begin(), hp2_family_indices.end());

    int n_hp1 = static_cast<int>(hp1_family_indices.size());
    int n_total_family = static_cast<int>(all_family_indices.size());

    for (int perm = 0; perm < config_.n_permutations; ++perm) {
        // Shuffle family indices
        std::shuffle(all_family_indices.begin(), all_family_indices.end(), rng_);

        // Split into pseudo HP1 and HP2 families
        std::vector<int> pseudo_hp1(all_family_indices.begin(), all_family_indices.begin() + n_hp1);
        std::vector<int> pseudo_hp2(all_family_indices.begin() + n_hp1, all_family_indices.end());

        // Recompute affinity score
        double perm_sum_d_hp1 = 0.0;
        double perm_sum_d_hp2 = 0.0;
        int perm_valid = 0;

        for (int u : unassigned_indices) {
            double d_hp1 = compute_mean_dist_to_family(u, pseudo_hp1);
            double d_hp2 = compute_mean_dist_to_family(u, pseudo_hp2);

            if (!std::isnan(d_hp1) && !std::isnan(d_hp2)) {
                perm_sum_d_hp1 += d_hp1;
                perm_sum_d_hp2 += d_hp2;
                ++perm_valid;
            }
        }

        if (perm_valid > 0) {
            double perm_mean_d_hp1 = perm_sum_d_hp1 / perm_valid;
            double perm_mean_d_hp2 = perm_sum_d_hp2 / perm_valid;
            double perm_denom = perm_mean_d_hp1 + perm_mean_d_hp2;
            double perm_score = (perm_denom > 0) ? (perm_mean_d_hp2 - perm_mean_d_hp1) / perm_denom : 0.0;

            if (std::abs(perm_score) >= observed_abs_score) {
                ++n_extreme;
            }
        }
    }

    result.unassigned.affinity_p = static_cast<double>(n_extreme) / static_cast<double>(config_.n_permutations + 1);
    result.unassigned.valid = true;
}

}  // namespace InterSubMod
