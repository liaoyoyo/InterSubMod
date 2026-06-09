/**
 * @file PerCpgAsm.cpp
 * @brief Per-CpG ASM detection and epiallele heterogeneity metrics
 *
 * Literature references:
 * - Per-CpG Fisher: DAMEfinder (Orjuela 2020), pycoMeth (Snajder 2023)
 * - NME: CPEL (Jenkinson 2020, Nature Communications)
 * - Epipolymorphism: methclone (Li 2014), Metheor (Lee 2023)
 */

#include "core/PerCpgAsm.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace InterSubMod {

// ============================================================================
// BH-FDR Correction
// ============================================================================

std::vector<double> bh_fdr_correction(const std::vector<double>& p_values) {
    int n = static_cast<int>(p_values.size());
    std::vector<double> fdr(n, std::numeric_limits<double>::quiet_NaN());

    // Collect valid (non-NaN) indices and p-values
    std::vector<int> valid_idx;
    valid_idx.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (!std::isnan(p_values[i])) {
            valid_idx.push_back(i);
        }
    }
    int m = static_cast<int>(valid_idx.size());
    if (m == 0) return fdr;

    // Sort valid indices by p-value (ascending)
    std::sort(valid_idx.begin(), valid_idx.end(),
              [&p_values](int a, int b) { return p_values[a] < p_values[b]; });

    // BH correction: q_i = p_i * m / rank_i, then enforce monotonicity from right
    std::vector<double> q(m);
    for (int i = 0; i < m; ++i) {
        q[i] = p_values[valid_idx[i]] * static_cast<double>(m) / static_cast<double>(i + 1);
    }

    // Enforce monotonicity: q[i] = min(q[i], q[i+1]) scanning right to left
    double running_min = q[m - 1];
    q[m - 1] = std::min(q[m - 1], 1.0);
    running_min = q[m - 1];
    for (int i = m - 2; i >= 0; --i) {
        q[i] = std::min(q[i], running_min);
        q[i] = std::min(q[i], 1.0);
        running_min = q[i];
    }

    // Write back
    for (int i = 0; i < m; ++i) {
        fdr[valid_idx[i]] = q[i];
    }
    return fdr;
}

// ============================================================================
// NME (Normalized Methylation Entropy)
// ============================================================================

double compute_nme(const Eigen::MatrixXi& binary_rows, int min_cpgs) {
    int n_reads = binary_rows.rows();
    int n_cpgs = binary_rows.cols();

    if (n_reads < 1 || n_cpgs < min_cpgs) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // For each read, build a binary pattern string from valid CpGs
    // then count unique patterns
    std::map<std::string, int> pattern_counts;
    int valid_reads = 0;

    for (int i = 0; i < n_reads; ++i) {
        std::string pattern;
        pattern.reserve(n_cpgs);
        bool has_valid = false;
        for (int j = 0; j < n_cpgs; ++j) {
            int val = binary_rows(i, j);
            if (val == -1) {
                pattern += 'N';
            } else {
                pattern += (val == 1 ? '1' : '0');
                has_valid = true;
            }
        }
        if (has_valid) {
            pattern_counts[pattern]++;
            valid_reads++;
        }
    }

    if (valid_reads < 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Shannon entropy H = -sum(p_i * log2(p_i))
    double H = 0.0;
    for (const auto& kv : pattern_counts) {
        double p = static_cast<double>(kv.second) / static_cast<double>(valid_reads);
        if (p > 0.0) {
            H -= p * std::log2(p);
        }
    }

    // H_max = log2(min(n_reads, 2^n_cpgs))
    // For large n_cpgs, 2^n_cpgs >> n_reads, so H_max = log2(n_reads)
    int valid_cpg_count = 0;
    for (int j = 0; j < n_cpgs; ++j) {
        bool any_valid = false;
        for (int i = 0; i < n_reads; ++i) {
            if (binary_rows(i, j) != -1) {
                any_valid = true;
                break;
            }
        }
        if (any_valid) valid_cpg_count++;
    }

    if (valid_cpg_count < min_cpgs) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // 2^n_cpgs can overflow for large n, but log2(2^n) = n
    double log2_max_patterns = static_cast<double>(valid_cpg_count);
    double log2_n_reads = std::log2(static_cast<double>(valid_reads));
    double H_max = std::min(log2_n_reads, log2_max_patterns);

    if (H_max <= 0.0) {
        return 0.0;  // Only 1 read or 0 valid CpGs → no entropy possible
    }

    return H / H_max;
}

// ============================================================================
// Epipolymorphism
// ============================================================================

double compute_epipolymorphism(const Eigen::MatrixXi& binary_rows, int window_size) {
    int n_reads = binary_rows.rows();
    int n_cpgs = binary_rows.cols();

    if (n_reads < 1 || n_cpgs < window_size) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Identify columns with at least one valid value
    std::vector<int> valid_cols;
    valid_cols.reserve(n_cpgs);
    for (int j = 0; j < n_cpgs; ++j) {
        for (int i = 0; i < n_reads; ++i) {
            if (binary_rows(i, j) != -1) {
                valid_cols.push_back(j);
                break;
            }
        }
    }

    int n_valid = static_cast<int>(valid_cols.size());
    if (n_valid < window_size) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Sliding window across valid columns
    int n_windows = n_valid - window_size + 1;
    double epipoly_sum = 0.0;
    int valid_windows = 0;

    for (int w = 0; w < n_windows; ++w) {
        // Build pattern for each read in this window
        std::map<int, int> pattern_counts;
        int valid_reads = 0;

        for (int i = 0; i < n_reads; ++i) {
            int pattern = 0;
            bool all_valid = true;
            for (int k = 0; k < window_size; ++k) {
                int col = valid_cols[w + k];
                int val = binary_rows(i, col);
                if (val == -1) {
                    all_valid = false;
                    break;
                }
                pattern |= (val << k);
            }
            if (all_valid) {
                pattern_counts[pattern]++;
                valid_reads++;
            }
        }

        if (valid_reads < 2) continue;

        // Epipolymorphism = 1 - sum(p_i^2)
        double sum_p2 = 0.0;
        for (const auto& kv : pattern_counts) {
            double p = static_cast<double>(kv.second) / static_cast<double>(valid_reads);
            sum_p2 += p * p;
        }

        epipoly_sum += (1.0 - sum_p2);
        valid_windows++;
    }

    if (valid_windows == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return epipoly_sum / static_cast<double>(valid_windows);
}

// ============================================================================
// Main: compute_per_cpg_asm
// ============================================================================

PerCpgAsmResult compute_per_cpg_asm(const MethylationMatrix& meth_mat,
                                     const std::vector<ReadInfo>& read_list,
                                     int min_reads_per_group) {
    PerCpgAsmResult result;

    int n_reads = meth_mat.num_reads();
    int n_cpgs = meth_mat.num_sites();

    if (n_reads < 2 || n_cpgs < 1) return result;

    // Build HP-family labels: 0 = HP1-family, 1 = HP2-family, -1 = other
    std::vector<int> hp_labels(n_reads, -1);
    std::vector<int> hp1_indices, hp2_indices;

    for (int i = 0; i < n_reads; ++i) {
        const std::string& hp = read_list[i].hp_tag;
        if (hp == "1" || hp == "HP1" || hp == "1-1") {
            hp_labels[i] = 0;
            hp1_indices.push_back(i);
        } else if (hp == "2" || hp == "HP2" || hp == "2-1") {
            hp_labels[i] = 1;
            hp2_indices.push_back(i);
        }
    }

    int n_hp1 = static_cast<int>(hp1_indices.size());
    int n_hp2 = static_cast<int>(hp2_indices.size());

    if (n_hp1 < min_reads_per_group || n_hp2 < min_reads_per_group) {
        return result;  // Insufficient HP group sizes
    }

    result.valid = true;

    // ========================================================================
    // 1. Per-CpG Fisher's exact test
    // ========================================================================

    FisherConfig fc;
    fc.seed = 42;
    FisherExact fisher(fc);

    std::vector<double> raw_pvalues(n_cpgs, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> deltas(n_cpgs, std::numeric_limits<double>::quiet_NaN());
    int n_tested = 0;

    for (int j = 0; j < n_cpgs; ++j) {
        // Count meth/unmeth for each HP group at this CpG
        int hp1_meth = 0, hp1_unmeth = 0;
        int hp2_meth = 0, hp2_unmeth = 0;

        for (int idx : hp1_indices) {
            double val = meth_mat.raw_matrix(idx, j);
            if (std::isnan(val)) continue;
            if (val > 0.5) hp1_meth++;
            else hp1_unmeth++;
        }
        for (int idx : hp2_indices) {
            double val = meth_mat.raw_matrix(idx, j);
            if (std::isnan(val)) continue;
            if (val > 0.5) hp2_meth++;
            else hp2_unmeth++;
        }

        int hp1_total = hp1_meth + hp1_unmeth;
        int hp2_total = hp2_meth + hp2_unmeth;

        if (hp1_total < 1 || hp2_total < 1) continue;

        // 2x2 table: rows = HP groups, cols = meth/unmeth
        ContingencyTable2x2 table;
        table.a = hp1_meth;
        table.b = hp1_unmeth;
        table.c = hp2_meth;
        table.d = hp2_unmeth;

        FisherResult fr = fisher.test_2x2(table);
        raw_pvalues[j] = fr.p_value;

        // Delta = HP1 meth fraction - HP2 meth fraction
        double hp1_frac = static_cast<double>(hp1_meth) / static_cast<double>(hp1_total);
        double hp2_frac = static_cast<double>(hp2_meth) / static_cast<double>(hp2_total);
        deltas[j] = hp1_frac - hp2_frac;

        n_tested++;
    }

    result.fisher_n_tested = n_tested;

    if (n_tested > 0) {
        // BH-FDR correction
        std::vector<double> fdr = bh_fdr_correction(raw_pvalues);

        int n_sig = 0;
        double max_neg_log_fdr = 0.0;

        for (int j = 0; j < n_cpgs; ++j) {
            if (std::isnan(fdr[j])) continue;
            if (fdr[j] < 0.05) {
                n_sig++;
            }
            double neg_log_fdr = -std::log10(std::max(fdr[j], 1e-300));
            if (neg_log_fdr > max_neg_log_fdr) {
                max_neg_log_fdr = neg_log_fdr;
            }
        }

        result.fisher_n_sig = n_sig;
        result.fisher_frac_sig = static_cast<double>(n_sig) / static_cast<double>(n_tested);
        result.fisher_max_neg_log_fdr = max_neg_log_fdr;
    }

    // ========================================================================
    // 2. NME per HP group
    // ========================================================================

    // Extract binary submatrices for each HP group
    auto extract_rows = [&](const std::vector<int>& indices) -> Eigen::MatrixXi {
        Eigen::MatrixXi sub(static_cast<int>(indices.size()), n_cpgs);
        for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
            sub.row(i) = meth_mat.binary_matrix.row(indices[i]);
        }
        return sub;
    };

    Eigen::MatrixXi hp1_binary = extract_rows(hp1_indices);
    Eigen::MatrixXi hp2_binary = extract_rows(hp2_indices);

    result.nme_hp1 = compute_nme(hp1_binary, 2);
    result.nme_hp2 = compute_nme(hp2_binary, 2);

    if (!std::isnan(result.nme_hp1) && !std::isnan(result.nme_hp2)) {
        result.entropy_imbalance = std::abs(result.nme_hp1 - result.nme_hp2);
    }

    // ========================================================================
    // 3. Epipolymorphism per HP group
    // ========================================================================

    result.epipoly_hp1 = compute_epipolymorphism(hp1_binary, 4);
    result.epipoly_hp2 = compute_epipolymorphism(hp2_binary, 4);

    if (!std::isnan(result.epipoly_hp1) && !std::isnan(result.epipoly_hp2)) {
        result.epipoly_delta = std::abs(result.epipoly_hp1 - result.epipoly_hp2);
    }

    return result;
}

}  // namespace InterSubMod
