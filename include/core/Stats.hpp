#pragma once

/**
 * @file Stats.hpp
 * @brief Statistical data structures for significance analysis
 *
 * Defines core data structures for:
 * - Contingency tables
 * - Cluster statistics
 * - Label encoding (dual-layer design)
 * - Test results
 */

#include <cstdint>
#include <string>
#include <vector>

#include "Types.hpp"

namespace InterSubMod {

// ============================================================================
// Label Encoding (Dual-Layer Design)
// ============================================================================

/**
 * @brief Low-dimensional labels for statistical tests
 *
 * These labels are used in Fisher/Chi-square tests to avoid
 * sparse RxC tables that would hurt Monte Carlo performance.
 */
enum class TestLabel : uint8_t {
    // Allele dimension
    ALLELE_ALT = 0,
    ALLELE_REF = 1,
    ALLELE_UNKNOWN = 2,

    // HP dimension (collapsed)
    HP_1 = 10,
    HP_2 = 11,
    HP_OTHER = 12,  // HP1-1, HP2-1, HP3, unphase

    // Sample dimension
    SAMPLE_TUMOR = 20,
    SAMPLE_NORMAL = 21
};

/**
 * @brief Full label for debugging and traceability
 *
 * High-dimensional combo code: HP × Allele × Strand × Sample
 */
struct FullLabel {
    AltSupport allele = AltSupport::UNKNOWN;
    std::string hp_tag;  // Original HP tag string
    Strand strand = Strand::UNKNOWN;
    bool is_tumor = true;

    /**
     * @brief Generate combo code
     *
     * Encoding:
     * - Bits 0-1: Allele (0=ALT, 1=REF, 2=UNKNOWN)
     * - Bits 2-4: HP (0=HP1, 1=HP2, 2=HP1-1, 3=HP2-1, 4=HP3, 5=unphase)
     * - Bit 5: Strand (0=FWD, 1=REV)
     * - Bit 6: Sample (0=Tumor, 1=Normal)
     */
    uint16_t combo_code() const;

    /**
     * @brief Get collapsed TestLabel for a specific dimension
     */
    TestLabel get_allele_label() const;
    TestLabel get_hp_label() const;
    TestLabel get_sample_label() const;
};

// ============================================================================
// Contingency Tables
// ============================================================================

/**
 * @brief 2x2 contingency table
 *
 *           | Category B | Not B |
 * ----------|------------|-------|
 * Category A|     a      |   b   | row1
 * Not A     |     c      |   d   | row2
 *           |   col1     | col2  |
 */
struct ContingencyTable2x2 {
    int a = 0;  // (row1, col1)
    int b = 0;  // (row1, col2)
    int c = 0;  // (row2, col1)
    int d = 0;  // (row2, col2)

    int row1_total() const {
        return a + b;
    }
    int row2_total() const {
        return c + d;
    }
    int col1_total() const {
        return a + c;
    }
    int col2_total() const {
        return b + d;
    }
    int grand_total() const {
        return a + b + c + d;
    }

    bool is_valid() const {
        return grand_total() > 0;
    }
};

/**
 * @brief Generic RxC contingency table (row-major storage)
 */
struct ContingencyTableRxC {
    std::vector<int> cells;
    int n_rows = 0;
    int n_cols = 0;

    int get(int row, int col) const {
        return cells[row * n_cols + col];
    }
    void set(int row, int col, int val) {
        cells[row * n_cols + col] = val;
    }

    int row_total(int row) const;
    int col_total(int col) const;
    int grand_total() const;

    bool is_valid() const {
        return grand_total() > 0 && n_rows > 0 && n_cols > 0;
    }
};

// ============================================================================
// Test Results
// ============================================================================

/**
 * @brief Result of a Fisher's exact test
 */
struct FisherResult {
    double p_value = 1.0;
    double odds_ratio = 1.0;
    double log_odds_ratio = 0.0;

    // For Monte Carlo methods
    int n_samples = 0;      // Number of MC samples used
    double ci_lower = 0.0;  // 95% CI lower bound for p-value
    double ci_upper = 1.0;  // 95% CI upper bound for p-value
    bool early_stopped = false;

    bool is_significant(double alpha = 0.05) const {
        return p_value <= alpha;
    }
};

/**
 * @brief Result of global association test
 */
struct GlobalTestResult {
    // Fisher-Freeman-Halton test (Monte Carlo)
    FisherResult fisher_ffh;

    // Chi-square test (if valid)
    double chi_square = 0.0;
    double chi_square_p = 1.0;
    bool chi_square_reliable = false;  // False if >20% cells have expected < 5

    // Effect size
    double cramers_v = 0.0;
    bool cramers_v_reliable = false;

    // Gating
    bool passed_gate = true;

    // Validity
    bool valid = true;
    std::string invalid_reason;
};

/**
 * @brief Statistics for a single cluster
 */
struct ClusterStats {
    int cluster_id = -1;
    int size = 0;

    // Counts by label
    int count_alt = 0;
    int count_ref = 0;
    int count_unknown = 0;
    int count_hp1 = 0;
    int count_hp2 = 0;
    int count_hp_other = 0;
    int count_tumor = 0;
    int count_normal = 0;

    // One-vs-Rest test results
    FisherResult fisher_allele;  // ALT vs REF
    FisherResult fisher_hp;      // HP1 vs HP2
    FisherResult fisher_sample;  // Tumor vs Normal

    double delta_proportion_alt = 0.0;  // P(ALT|this) - P(ALT|rest)
};

/**
 * @brief Result of local (per-cluster) tests
 */
struct LocalTestResult {
    std::vector<ClusterStats> cluster_stats;

    // Best cluster (lowest p-value)
    int best_cluster_id = -1;
    double best_p_value = 1.0;
    double best_log_odds_ratio = 0.0;
    std::string best_dimension;  // "allele", "hp", or "sample"
};

/**
 * @brief Result of PERMANOVA test
 */
struct PermanovaResult {
    double pseudo_f = 0.0;
    double p_value = 1.0;
    int n_permutations = 0;

    bool valid = true;
    std::string invalid_reason;
};

/**
 * @brief Result of dispersion check
 */
struct DispersionResult {
    std::vector<double> mean_distances_to_centroid;  // Per cluster
    double anova_f = 0.0;
    double anova_p = 1.0;
    bool warning = false;  // True if significant heterogeneity
};

/**
 * @brief Result of bootstrap stability test
 */
struct BootstrapResult {
    double mean_ari = 0.0;
    double std_ari = 0.0;
    std::vector<double> ari_samples;
    int n_iterations = 0;
    bool early_stopped = false;
    bool stable = true;  // True if mean_ari >= 0.2
};

/**
 * @brief Aggregated significance result for a region
 */
struct SignificanceResult {
    // Join keys
    std::string run_id;
    std::string vcf_id;
    std::string bam_id;
    int region_id = -1;
    std::string anchor_key;  // chrom:pos:ref:alt

    // Validity
    bool valid_flag = true;
    std::string invalid_reason;

    // Quality metrics
    int n_reads = 0;
    int n_reads_valid = 0;
    int n_cpg = 0;
    double effective_overlap_median = 0.0;
    double mapq_mean = 0.0;

    // Global test results
    GlobalTestResult global_alt;
    GlobalTestResult global_hp;
    GlobalTestResult global_sample;

    // Local test results
    LocalTestResult local_result;

    // Structure tests
    PermanovaResult permanova;
    DispersionResult dispersion;
    BootstrapResult bootstrap;

    // Final judgment
    bool passed_gate = true;
    bool dispersion_warning = false;
    double heuristic_score = 0.0;
    bool score_ambiguous = false;

    // Random seed used
    uint64_t seed = 0;
};

// ============================================================================
// Label Codebook
// ============================================================================

/**
 * @brief Entry in label codebook
 */
struct LabelCodeEntry {
    uint16_t code;
    std::string allele_str;
    std::string hp_str;
    std::string strand_str;
    std::string sample_str;
};

/**
 * @brief Schema version for output
 */
struct SignificanceSchema {
    std::string version = "1.0.0";
    std::vector<std::string> column_names;
    std::vector<std::string> column_types;
};

}  // namespace InterSubMod
