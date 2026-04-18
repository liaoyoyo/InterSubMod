/**
 * @file LocalTest.cpp
 * @brief Implementation of local (per-cluster) association tests
 */

#include "core/LocalTest.hpp"

#include <algorithm>
#include <set>

#include "core/MathUtils.hpp"

namespace InterSubMod {

// ============================================================================
// Constructor
// ============================================================================

LocalTest::LocalTest(const LocalTestConfig& config) : config_(config), fisher_(config.fisher_config) {
}

// ============================================================================
// Initialize Cluster Stats
// ============================================================================

std::vector<ClusterStats> LocalTest::initialize_cluster_stats(const std::vector<int>& cluster_labels,
                                                              const std::vector<FullLabel>& full_labels) {
    // Find unique clusters
    std::set<int> unique_clusters(cluster_labels.begin(), cluster_labels.end());

    std::vector<ClusterStats> stats;
    for (int cluster_id : unique_clusters) {
        ClusterStats cs;
        cs.cluster_id = cluster_id;
        cs.size = 0;

        // Count reads in this cluster
        for (size_t i = 0; i < cluster_labels.size(); ++i) {
            if (cluster_labels[i] == cluster_id) {
                cs.size++;

                const FullLabel& label = full_labels[i];

                // Count by allele
                switch (label.allele) {
                    case AltSupport::ALT:
                        cs.count_alt++;
                        break;
                    case AltSupport::REF:
                        cs.count_ref++;
                        break;
                    default:
                        cs.count_unknown++;
                        break;
                }

                // Count by HP (fine-grained + backward-compatible collapsed)
                const std::string& hp = label.hp_tag;
                if (hp == "1") {
                    cs.count_hp1++;
                    cs.count_hp_family_1++;
                } else if (hp == "2") {
                    cs.count_hp2++;
                    cs.count_hp_family_2++;
                } else if (hp == "1-1") {
                    cs.count_hp1_1++;
                    cs.count_hp_family_1++;
                    cs.count_hp_other++;
                } else if (hp == "2-1") {
                    cs.count_hp2_1++;
                    cs.count_hp_family_2++;
                    cs.count_hp_other++;
                } else if (hp == "3") {
                    cs.count_hp3++;
                    cs.count_hp_other++;
                } else {
                    cs.count_hp0++;
                    cs.count_hp_other++;
                }

                // Count by sample
                if (label.is_tumor) {
                    cs.count_tumor++;
                } else {
                    cs.count_normal++;
                }
            }
        }

        stats.push_back(cs);
    }

    return stats;
}

// ============================================================================
// Build One-vs-Rest Table
// ============================================================================

ContingencyTable2x2 LocalTest::build_one_vs_rest_table(int cluster_id, const std::vector<int>& cluster_labels,
                                                       int target_label, const std::vector<int>& binary_labels,
                                                       const std::vector<bool>& valid_mask) {
    ContingencyTable2x2 table;

    for (size_t i = 0; i < cluster_labels.size(); ++i) {
        if (!valid_mask[i]) continue;

        bool in_cluster = (cluster_labels[i] == cluster_id);
        bool has_target = (binary_labels[i] == target_label);

        if (in_cluster && has_target) {
            table.a++;
        } else if (in_cluster && !has_target) {
            table.c++;
        } else if (!in_cluster && has_target) {
            table.b++;
        } else {
            table.d++;
        }
    }

    return table;
}

// ============================================================================
// Compute Delta Proportion
// ============================================================================

double LocalTest::compute_delta_proportion(const ContingencyTable2x2& table) {
    // P(target | this cluster) - P(target | other clusters)
    // = a / (a + c) - b / (b + d)

    int in_cluster_total = table.a + table.c;
    int out_cluster_total = table.b + table.d;

    if (in_cluster_total == 0 || out_cluster_total == 0) {
        return 0.0;
    }

    double p_in = static_cast<double>(table.a) / static_cast<double>(in_cluster_total);
    double p_out = static_cast<double>(table.b) / static_cast<double>(out_cluster_total);

    return p_in - p_out;
}

// ============================================================================
// Test Allele Local
// ============================================================================

void LocalTest::test_allele_local(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                                  std::vector<ClusterStats>& stats) {
    size_t n = cluster_labels.size();
    std::vector<int> binary_labels(n);
    std::vector<bool> valid_mask(n);

    for (size_t i = 0; i < n; ++i) {
        TestLabel label = full_labels[i].get_allele_label();
        if (label == TestLabel::ALLELE_ALT) {
            binary_labels[i] = 1;  // Target = ALT
            valid_mask[i] = true;
        } else if (label == TestLabel::ALLELE_REF) {
            binary_labels[i] = 0;
            valid_mask[i] = true;
        } else {
            binary_labels[i] = -1;
            valid_mask[i] = false;
        }
    }

    for (auto& cs : stats) {
        ContingencyTable2x2 table =
            build_one_vs_rest_table(cs.cluster_id, cluster_labels, 1, binary_labels, valid_mask);

        cs.fisher_allele = fisher_.test_2x2(table);
        cs.delta_proportion_alt = compute_delta_proportion(table);
    }
}

// ============================================================================
// Test HP Local
// ============================================================================

void LocalTest::test_hp_local(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                              std::vector<ClusterStats>& stats) {
    size_t n = cluster_labels.size();
    std::vector<int> binary_labels(n);
    std::vector<bool> valid_mask(n);

    for (size_t i = 0; i < n; ++i) {
        TestLabel label = full_labels[i].get_hp_label();
        if (label == TestLabel::HP_1) {
            binary_labels[i] = 1;  // Target = HP1
            valid_mask[i] = true;
        } else if (label == TestLabel::HP_2) {
            binary_labels[i] = 0;
            valid_mask[i] = true;
        } else {
            binary_labels[i] = -1;
            valid_mask[i] = false;
        }
    }

    for (auto& cs : stats) {
        ContingencyTable2x2 table =
            build_one_vs_rest_table(cs.cluster_id, cluster_labels, 1, binary_labels, valid_mask);

        cs.fisher_hp = fisher_.test_2x2(table);
    }
}

// ============================================================================
// Test HP Family Local
// ============================================================================

void LocalTest::test_hp_family_local(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                                     std::vector<ClusterStats>& stats) {
    size_t n = cluster_labels.size();
    std::vector<int> binary_labels(n);
    std::vector<bool> valid_mask(n);

    for (size_t i = 0; i < n; ++i) {
        TestLabel label = full_labels[i].get_hp_family_label();
        if (label == TestLabel::HP_FAMILY_1) {
            binary_labels[i] = 1;  // Target = HP-family-1
            valid_mask[i] = true;
        } else if (label == TestLabel::HP_FAMILY_2) {
            binary_labels[i] = 0;
            valid_mask[i] = true;
        } else {
            binary_labels[i] = -1;
            valid_mask[i] = false;
        }
    }

    for (auto& cs : stats) {
        ContingencyTable2x2 table =
            build_one_vs_rest_table(cs.cluster_id, cluster_labels, 1, binary_labels, valid_mask);

        cs.fisher_hp_family = fisher_.test_2x2(table);
    }
}

// ============================================================================
// Test Sample Local
// ============================================================================

void LocalTest::test_sample_local(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels,
                                  std::vector<ClusterStats>& stats) {
    size_t n = cluster_labels.size();
    std::vector<int> binary_labels(n);
    std::vector<bool> valid_mask(n, true);

    for (size_t i = 0; i < n; ++i) {
        binary_labels[i] = full_labels[i].is_tumor ? 1 : 0;  // Target = Tumor
    }

    for (auto& cs : stats) {
        ContingencyTable2x2 table =
            build_one_vs_rest_table(cs.cluster_id, cluster_labels, 1, binary_labels, valid_mask);

        cs.fisher_sample = fisher_.test_2x2(table);
    }
}

// ============================================================================
// Run All Local Tests
// ============================================================================

LocalTestResult LocalTest::test_all(const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels) {
    LocalTestResult result;

    if (cluster_labels.empty()) {
        return result;
    }

    // Initialize cluster statistics
    result.cluster_stats = initialize_cluster_stats(cluster_labels, full_labels);

    // Run tests for each dimension
    test_allele_local(cluster_labels, full_labels, result.cluster_stats);
    test_hp_local(cluster_labels, full_labels, result.cluster_stats);
    test_hp_family_local(cluster_labels, full_labels, result.cluster_stats);
    test_sample_local(cluster_labels, full_labels, result.cluster_stats);

    // Find best cluster (lowest p-value across all dimensions)
    result.best_p_value = 1.0;
    result.best_cluster_id = -1;

    for (const auto& cs : result.cluster_stats) {
        // Check allele Fisher
        if (cs.fisher_allele.p_value < result.best_p_value) {
            result.best_p_value = cs.fisher_allele.p_value;
            result.best_log_odds_ratio = cs.fisher_allele.log_odds_ratio;
            result.best_cluster_id = cs.cluster_id;
            result.best_dimension = "allele";
        }

        // Check HP Fisher (pure germline)
        if (cs.fisher_hp.p_value < result.best_p_value) {
            result.best_p_value = cs.fisher_hp.p_value;
            result.best_log_odds_ratio = cs.fisher_hp.log_odds_ratio;
            result.best_cluster_id = cs.cluster_id;
            result.best_dimension = "hp";
        }

        // Check HP Family Fisher
        if (cs.fisher_hp_family.p_value < result.best_p_value) {
            result.best_p_value = cs.fisher_hp_family.p_value;
            result.best_log_odds_ratio = cs.fisher_hp_family.log_odds_ratio;
            result.best_cluster_id = cs.cluster_id;
            result.best_dimension = "hp_family";
        }

        // Check Sample Fisher
        if (cs.fisher_sample.p_value < result.best_p_value) {
            result.best_p_value = cs.fisher_sample.p_value;
            result.best_log_odds_ratio = cs.fisher_sample.log_odds_ratio;
            result.best_cluster_id = cs.cluster_id;
            result.best_dimension = "sample";
        }
    }

    return result;
}

}  // namespace InterSubMod
