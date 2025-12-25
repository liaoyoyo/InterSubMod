/**
 * @file GlobalTest.cpp
 * @brief Implementation of global association tests
 */

#include "core/GlobalTest.hpp"

#include <algorithm>
#include <map>
#include <set>

#include "core/MathUtils.hpp"

namespace InterSubMod {

// ============================================================================
// Constructor
// ============================================================================

GlobalTest::GlobalTest(const GlobalTestConfig& config) : config_(config), fisher_(config.fisher_config) {
}

// ============================================================================
// Build Contingency Table
// ============================================================================

ContingencyTableRxC GlobalTest::build_contingency_table(const std::vector<int>& cluster_labels,
                                                        const std::vector<int>& binary_labels,
                                                        const std::vector<bool>& valid_mask) {
    // Determine number of clusters
    std::set<int> unique_clusters;
    for (size_t i = 0; i < cluster_labels.size(); ++i) {
        if (valid_mask[i]) {
            unique_clusters.insert(cluster_labels[i]);
        }
    }

    // Map cluster IDs to consecutive indices
    std::map<int, int> cluster_to_idx;
    int idx = 0;
    for (int c : unique_clusters) {
        cluster_to_idx[c] = idx++;
    }

    int n_clusters = static_cast<int>(unique_clusters.size());
    int n_labels = 2;  // Binary labels (0 and 1)

    ContingencyTableRxC table;
    table.n_rows = n_clusters;
    table.n_cols = n_labels;
    table.cells.resize(n_clusters * n_labels, 0);

    // Fill table
    for (size_t i = 0; i < cluster_labels.size(); ++i) {
        if (!valid_mask[i]) continue;

        int row = cluster_to_idx[cluster_labels[i]];
        int col = binary_labels[i];  // 0 or 1

        if (col >= 0 && col < n_labels) {
            table.cells[row * n_labels + col]++;
        }
    }

    return table;
}

// ============================================================================
// Compute Result from Table
// ============================================================================

GlobalTestResult GlobalTest::compute_from_table(const ContingencyTableRxC& table) {
    GlobalTestResult result;

    if (!table.is_valid() || table.n_rows < 1) {
        result.valid = false;
        result.invalid_reason = "empty_or_invalid_table";
        return result;
    }

    // Fisher-Freeman-Halton test (Monte Carlo)
    result.fisher_ffh = fisher_.test_rxc_mc(table);

    // Chi-square test
    std::vector<int> cells = table.cells;
    double min_expected;
    result.chi_square = MathUtils::chi_square(cells, table.n_rows, table.n_cols, min_expected);

    // Chi-square reliability check
    // Count proportion of cells with expected >= 5
    int n = table.grand_total();
    if (n > 0) {
        std::vector<int> row_totals(table.n_rows, 0);
        std::vector<int> col_totals(table.n_cols, 0);

        for (int i = 0; i < table.n_rows; ++i) {
            for (int j = 0; j < table.n_cols; ++j) {
                int val = table.get(i, j);
                row_totals[i] += val;
                col_totals[j] += val;
            }
        }

        int valid_cells = 0;
        int total_cells = table.n_rows * table.n_cols;
        double n_d = static_cast<double>(n);

        for (int i = 0; i < table.n_rows; ++i) {
            for (int j = 0; j < table.n_cols; ++j) {
                double expected = (static_cast<double>(row_totals[i]) * static_cast<double>(col_totals[j])) / n_d;
                if (expected >= 5.0) {
                    ++valid_cells;
                }
            }
        }

        double valid_ratio = static_cast<double>(valid_cells) / static_cast<double>(total_cells);
        result.chi_square_reliable = (valid_ratio >= config_.chi_square_min_valid_ratio);

        // Chi-square p-value (approximation using chi-square distribution)
        // For now, we skip p-value computation and rely on Fisher test
        // In production, you'd use a chi-square CDF
        result.chi_square_p = result.fisher_ffh.p_value;  // Placeholder
    }

    // Cramér's V
    result.cramers_v = MathUtils::cramers_v(cells, table.n_rows, table.n_cols, result.cramers_v_reliable);

    result.valid = true;
    return result;
}

// ============================================================================
// Apply Gating
// ============================================================================

void GlobalTest::apply_gating(GlobalTestResult& result) {
    // If p-value > threshold, mark as not passing gate
    result.passed_gate = (result.fisher_ffh.p_value <= config_.gating_p_threshold);
}

// ============================================================================
// Test Allele Association
// ============================================================================

GlobalTestResult GlobalTest::test_allele(const std::vector<int>& cluster_labels,
                                         const std::vector<FullLabel>& full_labels) {
    size_t n = cluster_labels.size();
    std::vector<int> binary_labels(n);
    std::vector<bool> valid_mask(n);

    for (size_t i = 0; i < n; ++i) {
        TestLabel label = full_labels[i].get_allele_label();
        if (label == TestLabel::ALLELE_ALT) {
            binary_labels[i] = 0;
            valid_mask[i] = true;
        } else if (label == TestLabel::ALLELE_REF) {
            binary_labels[i] = 1;
            valid_mask[i] = true;
        } else {
            binary_labels[i] = -1;
            valid_mask[i] = false;  // Exclude UNKNOWN
        }
    }

    ContingencyTableRxC table = build_contingency_table(cluster_labels, binary_labels, valid_mask);
    GlobalTestResult result = compute_from_table(table);
    apply_gating(result);

    return result;
}

// ============================================================================
// Test HP Association
// ============================================================================

GlobalTestResult GlobalTest::test_hp(const std::vector<int>& cluster_labels,
                                     const std::vector<FullLabel>& full_labels) {
    size_t n = cluster_labels.size();
    std::vector<int> binary_labels(n);
    std::vector<bool> valid_mask(n);

    for (size_t i = 0; i < n; ++i) {
        TestLabel label = full_labels[i].get_hp_label();
        if (label == TestLabel::HP_1) {
            binary_labels[i] = 0;
            valid_mask[i] = true;
        } else if (label == TestLabel::HP_2) {
            binary_labels[i] = 1;
            valid_mask[i] = true;
        } else {
            // HP_OTHER: exclude from HP1 vs HP2 comparison
            binary_labels[i] = -1;
            valid_mask[i] = false;
        }
    }

    ContingencyTableRxC table = build_contingency_table(cluster_labels, binary_labels, valid_mask);
    GlobalTestResult result = compute_from_table(table);
    apply_gating(result);

    return result;
}

// ============================================================================
// Test Sample Association
// ============================================================================

GlobalTestResult GlobalTest::test_sample(const std::vector<int>& cluster_labels,
                                         const std::vector<FullLabel>& full_labels) {
    size_t n = cluster_labels.size();
    std::vector<int> binary_labels(n);
    std::vector<bool> valid_mask(n, true);

    for (size_t i = 0; i < n; ++i) {
        TestLabel label = full_labels[i].get_sample_label();
        binary_labels[i] = (label == TestLabel::SAMPLE_TUMOR) ? 0 : 1;
    }

    ContingencyTableRxC table = build_contingency_table(cluster_labels, binary_labels, valid_mask);
    GlobalTestResult result = compute_from_table(table);
    apply_gating(result);

    return result;
}

// ============================================================================
// Test All
// ============================================================================

std::tuple<GlobalTestResult, GlobalTestResult, GlobalTestResult> GlobalTest::test_all(
    const std::vector<int>& cluster_labels, const std::vector<FullLabel>& full_labels) {
    return {test_allele(cluster_labels, full_labels), test_hp(cluster_labels, full_labels),
            test_sample(cluster_labels, full_labels)};
}

}  // namespace InterSubMod
