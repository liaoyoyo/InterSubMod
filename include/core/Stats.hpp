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

#include <array>
#include <cmath>
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

    // HP dimension (collapsed — pure germline only)
    HP_1 = 10,
    HP_2 = 11,
    HP_OTHER = 12,  // HP1-1, HP2-1, HP3, unphase

    // HP dimension (family-level — germline + somatic merged)
    HP_FAMILY_1 = 13,  // HP1 + HP1-1
    HP_FAMILY_2 = 14,  // HP2 + HP2-1

    // HP dimension (fine-grained — 4-group)
    HP_FINE_1 = 15,    // Pure germline HP1
    HP_FINE_1_1 = 16,  // Somatic on HP1
    HP_FINE_2 = 17,    // Pure germline HP2
    HP_FINE_2_1 = 18,  // Somatic on HP2

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

    // Lineage axes, populated from the BAM aux tags written by LongLineage's
    // tag_bam. Empty when the input never went through that step, in which case
    // the lineage axes are simply skipped rather than treated as one big group.
    std::string lineage_component;  // lc:Z
    std::string lineage_block;      // lu:Z
    std::string lineage_path;       // lv:Z, e.g. "HP2-1-1" or "HP2-1+" (subtree)
    char lineage_status = '\0';     // ls:A  U | M | P | A ; '\0' when absent

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
    TestLabel get_hp_family_label() const;
    TestLabel get_hp_fine_label() const;
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

    // Multi-group info (for HP fine-grained test)
    int n_valid_groups = 0;  ///< Number of valid groups (after sparse exclusion)

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
    int count_hp_other = 0;  // Backward-compatible: HP1-1 + HP2-1 + HP3 + HP0
    int count_hp1_1 = 0;     // Somatic on haplotype 1
    int count_hp2_1 = 0;     // Somatic on haplotype 2
    int count_hp3 = 0;       // Unresolvable ALT haplotype
    int count_hp0 = 0;       // Unphased
    int count_hp_family_1 = 0;  // HP1 + HP1-1
    int count_hp_family_2 = 0;  // HP2 + HP2-1
    int count_tumor = 0;
    int count_normal = 0;

    // One-vs-Rest test results
    FisherResult fisher_allele;      // ALT vs REF
    FisherResult fisher_hp;          // HP1 vs HP2 (pure germline)
    FisherResult fisher_hp_family;   // HP1-family vs HP2-family
    FisherResult fisher_sample;      // Tumor vs Normal

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
    double r_squared = 0.0;
    int n_permutations_realized = 0;
    std::string permutation_mode = "unrestricted";
    std::string evaluation_status = "EVALUABLE";

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
    int n_permutations = 0;
    int n_permutations_realized = 0;
    std::string permutation_mode = "analytic_f";
    std::string evaluation_status = "EVALUABLE";
    bool valid = true;
    std::string invalid_reason;
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
 * @brief Result of Label-First Delta test for one label dimension
 *
 * Measures the difference between between-group and within-group distances
 * to assess if a label (HP/Allele) explains the distance structure.
 */
struct LabelDeltaResult {
    // Distance statistics
    double within_mean = 0.0;   // Mean pairwise distance within same label
    double between_mean = 0.0;  // Mean pairwise distance between different labels
    double delta = 0.0;         // between_mean - within_mean

    // Permutation test results
    double p_value = 1.0;       // Permutation p-value
    int n_permutations = 0;     // Actual number of permutations completed
    bool significant = false;   // delta > 0 and p_value <= alpha

    // Group statistics
    int n_group_a = 0;          // Size of group A (e.g., HP1 or ALT)
    int n_group_b = 0;          // Size of group B (e.g., HP2 or REF)
    int n_pairs_within = 0;     // Number of within-group pairs used
    int n_pairs_between = 0;    // Number of between-group pairs used

    // Validity
    bool valid = true;
    std::string invalid_reason;
};

// ============================================================================
// Multi-Stage HP Verification Structures (NEW)
// ============================================================================

/**
 * @brief Pairwise distances between HP fine-grained groups
 *
 * Used in Stage 2 (HP Fine-Grained Test) to understand
 * the distance structure between germline and somatic subgroups.
 */
struct HPFinePairwise {
    double d_hp1_hp1s = NAN;   // d(HP1, HP1-1): Same haplotype, germline vs somatic
    double d_hp1_hp2 = NAN;    // d(HP1, HP2): Different haplotype, both germline
    double d_hp1_hp2s = NAN;   // d(HP1, HP2-1): Different haplotype, germline vs somatic
    double d_hp1s_hp2 = NAN;   // d(HP1-1, HP2): Different haplotype, somatic vs germline
    double d_hp1s_hp2s = NAN;  // d(HP1-1, HP2-1): Different haplotype, both somatic
    double d_hp2_hp2s = NAN;   // d(HP2, HP2-1): Same haplotype, germline vs somatic
};

/**
 * @brief Result of unassigned reads affinity test (Stage 4)
 *
 * Determines if HP3/HP0 reads are closer to HP1-family or HP2-family
 * based on mean distance comparison.
 */
struct UnassignedAffinityResult {
    // Affinity score: (d_to_hp2 - d_to_hp1) / (d_to_hp1 + d_to_hp2)
    // Positive = closer to HP1-family, Negative = closer to HP2-family
    double affinity_score = 0.0;
    double affinity_p = 1.0;                   // Permutation p-value
    std::string affinity_direction = "None";   // "HP1", "HP2", or "None"

    // Group counts
    int n_hp3 = 0;   // Number of HP3 reads (somatic unassignable)
    int n_hp0 = 0;   // Number of HP0/unphased reads

    // Validity
    bool valid = false;
    std::string invalid_reason;
};

/**
 * @brief Multi-stage HP verification result
 *
 * Contains results from all HP-related verification stages:
 * - Stage 1: HP Family Merged Test (HP1+HP1-1 vs HP2+HP2-1)
 * - Stage 2: HP Fine-Grained Test (4-group PERMANOVA)
 * - Stage 4: Unassigned Affinity Test (HP3/HP0 → HP1/HP2)
 */
struct MultiStageHPResult {
    // Stage 1: HP Family Merged Test
    LabelDeltaResult merged;    // Binary delta test: HP1-family vs HP2-family
    int n_hp1_family = 0;       // HP1 + HP1-1 count
    int n_hp2_family = 0;       // HP2 + HP2-1 count

    // Stage 2: HP Fine-Grained Test (4-group PERMANOVA)
    double fine_f = 0.0;                        // Pseudo-F statistic
    double fine_p = 1.0;                        // Permutation p-value
    bool fine_sig = false;                      // Significant at alpha
    int fine_n_groups = 0;                      // Number of valid groups (max 4)
    std::array<int, 4> fine_group_counts = {};  // [HP1, HP1-1, HP2, HP2-1]
    HPFinePairwise fine_pairwise;               // Pairwise mean distances

    // Stage 4: Unassigned Affinity Test
    UnassignedAffinityResult unassigned;

    // Overall validity
    bool valid = true;
    std::string invalid_reason;
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

    // Global test results (Cluster-First)
    GlobalTestResult global_alt;
    GlobalTestResult global_hp;          // Layer 2: Pure germline HP1 vs HP2
    GlobalTestResult global_hp_family;   // Layer 1: HP1-family vs HP2-family
    GlobalTestResult global_hp_fine;     // Layer 3: Fine-grained 4-group
    GlobalTestResult global_sample;

    // Local test results
    LocalTestResult local_result;

    // Structure tests
    PermanovaResult permanova;
    DispersionResult dispersion;
    BootstrapResult bootstrap;

    // Multi-Stage HP Verification (NEW)
    MultiStageHPResult hp_multistage;

    // Label-First verification results
    // Note: label_hp is deprecated, use hp_multistage.merged instead
    LabelDeltaResult label_hp;      // HP1 vs HP2 delta test (deprecated)
    LabelDeltaResult label_allele;  // ALT vs REF delta test (Stage 3)
    LabelDeltaResult label_sample;  // Tumor vs Normal delta test (Sample ASM)
    PermanovaResult label_hp_permanova;       // HP merged labels on distance matrix
    DispersionResult label_hp_dispersion;
    bool label_hp_dispersion_warning = false;
    PermanovaResult label_allele_permanova;   // ALT/REF labels on distance matrix
    DispersionResult label_allele_dispersion;
    bool label_allele_dispersion_warning = false;
    std::string dominant_label_dimension;  // "hp", "allele", or "none"
    bool label_significant = false;

    // Bidirectional verification classification
    std::string verification_class = "Noise";  // "Strong", "Subclone", "Weak", or "Noise"

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

// ============================================================================
// Quality Assessment Structures (NEW - Multi-Layer Validation)
// ============================================================================

/**
 * @brief Coverage category for CNV detection
 */
enum class CoverageCategory : uint8_t {
    CNV_LOSS_CANDIDATE = 0,  // < 0.5x expected coverage
    LOW_COVERAGE = 1,        // 0.5-0.8x expected coverage
    NORMAL = 2,              // 0.8-1.2x expected coverage
    ELEVATED = 3,            // 1.2-1.5x expected coverage
    CNV_GAIN_CANDIDATE = 4,  // 1.5-2.0x expected coverage
    HIGH_COPY = 5            // > 2.0x expected coverage
};

/**
 * @brief LOH subtype classification
 *
 * Classifies LOH events based on their verification class to
 * understand the nature of allelic imbalance.
 */
enum class LOHSubtype : uint8_t {
    NONE = 0,            // Not LOH (balanced HP)
    LOH_NOISE = 1,       // LOH + Noise verification class
    LOH_WEAK = 2,        // LOH + Weak verification class
    LOH_STRONG = 3,      // LOH + Strong verification class (may be real signal)
    LOH_SUBCLONE = 4     // LOH + Subclone verification class (high priority)
};

/**
 * @brief Quality tier for confidence assessment
 */
enum class QualityTier : uint8_t {
    HIGH = 0,     // Quality score >= 70
    MEDIUM = 1,   // Quality score 40-69
    LOW = 2       // Quality score < 40
};

/**
 * @brief Quality assessment flags for multi-layer validation
 *
 * Contains derived metrics for LOH detection, coverage validation,
 * and overall quality scoring to enable more accurate classification.
 */
struct QualityFlags {
    // HP balance metrics
    double hp_ratio = 0.5;         // HP1/(HP1+HP2), range [0,1]
    bool potential_loh = false;    // True if hp_ratio < 0.1 or > 0.9

    // Coverage metrics (normalized to expected 75x)
    double coverage_multiple = 1.0;     // NumReads/75.0
    CoverageCategory coverage_category = CoverageCategory::NORMAL;

    // LOH classification
    LOHSubtype loh_subtype = LOHSubtype::NONE;

    // Overall quality
    float quality_score = 50.0f;        // 0-100 composite score
    QualityTier quality_tier = QualityTier::MEDIUM;

    // Confidence adjustments applied
    float coverage_confidence_adj = 0.0f;   // Adjustment from coverage
    float loh_confidence_adj = 0.0f;        // Adjustment from LOH status

    /**
     * @brief Get coverage category as string
     */
    static std::string coverage_category_to_string(CoverageCategory cat) {
        switch (cat) {
            case CoverageCategory::CNV_LOSS_CANDIDATE: return "CNV_Loss";
            case CoverageCategory::LOW_COVERAGE: return "Low";
            case CoverageCategory::NORMAL: return "Normal";
            case CoverageCategory::ELEVATED: return "Elevated";
            case CoverageCategory::CNV_GAIN_CANDIDATE: return "CNV_Gain";
            case CoverageCategory::HIGH_COPY: return "High_Copy";
            default: return "Unknown";
        }
    }

    /**
     * @brief Get LOH subtype as string
     */
    static std::string loh_subtype_to_string(LOHSubtype subtype) {
        switch (subtype) {
            case LOHSubtype::NONE: return "None";
            case LOHSubtype::LOH_NOISE: return "LOH_Noise";
            case LOHSubtype::LOH_WEAK: return "LOH_Weak";
            case LOHSubtype::LOH_STRONG: return "LOH_Strong";
            case LOHSubtype::LOH_SUBCLONE: return "LOH_Subclone";
            default: return "Unknown";
        }
    }

    /**
     * @brief Get quality tier as string
     */
    static std::string quality_tier_to_string(QualityTier tier) {
        switch (tier) {
            case QualityTier::HIGH: return "High";
            case QualityTier::MEDIUM: return "Medium";
            case QualityTier::LOW: return "Low";
            default: return "Unknown";
        }
    }
};

}  // namespace InterSubMod
