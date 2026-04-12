#pragma once

#include <string>
#include <vector>

namespace InterSubMod {

// Forward declaration
struct RegionResult;

/**
 * @brief Per-region ASM feature profile for cross-region subclone analysis.
 *
 * Aggregates key metrics from RegionResult into a compact feature vector
 * suitable for clustering regions into subclone groups.
 */
struct RegionAsmProfile {
    int region_id;
    std::string chrom;
    uint32_t pos;

    // Core ASM features
    double hp_asm_delta;        ///< HP merged delta (haplotype ASM strength)
    double sample_asm_delta;    ///< Sample ASM delta (tumor vs normal)
    double hp_fine_f;           ///< Fine-grained heterogeneity F-statistic
    double allele_delta;        ///< Allele (ALT vs REF) delta

    // HP and coverage features
    double hp_ratio;            ///< HP1/(HP1+HP2)
    double coverage_multiple;   ///< Coverage relative to diploid expectation
    bool potential_loh;         ///< ISM hp_ratio-based LOH flag

    // LOH BED annotation
    bool loh_bed_overlap;       ///< External LOH BED overlap
    std::string loh_source;     ///< "none", "bed_only", "ratio_only", "both"

    // Verification class
    std::string verification_class; ///< "Strong", "Subclone", "Weak", "Noise"

    // Quality
    float quality_score;
    int num_reads;
    int num_cpgs;
};

/**
 * @brief Summary statistics for a single subclone group.
 */
struct SubcloneSummary {
    int subclone_id;
    int n_regions;

    // Mean feature values across regions in this subclone
    double mean_hp_asm;
    double mean_sample_asm;
    double mean_hp_ratio;
    double mean_hp_fine_f;
    double mean_coverage_multiple;

    // LOH statistics
    int n_loh_bed;          ///< Regions overlapping LOH BED
    int n_loh_ratio;        ///< Regions with ISM hp_ratio LOH
    int n_loh_both;         ///< Regions with both LOH sources

    // Verification class distribution
    int n_strong;
    int n_subclone;
    int n_weak;
    int n_noise;
    std::string dominant_verification_class;
};

/**
 * @brief Result of cross-region subclone analysis.
 */
struct SubcloneResult {
    int n_subclones;
    std::vector<int> region_assignments;    ///< Per-region subclone assignment (-1 = unassigned)
    std::vector<SubcloneSummary> summaries; ///< Per-subclone summary statistics
    bool valid;                             ///< Whether analysis was successful

    // Overall sample-level statistics
    int total_regions_analyzed;
    int total_regions_assigned;
    double mean_hp_asm_all;
    double mean_sample_asm_all;
    double loh_fraction;                    ///< Fraction of regions with any LOH evidence

    SubcloneResult() : n_subclones(0), valid(false), total_regions_analyzed(0),
                       total_regions_assigned(0), mean_hp_asm_all(0.0),
                       mean_sample_asm_all(0.0), loh_fraction(0.0) {}
};

/**
 * @brief Cross-region subclone structure analyzer.
 *
 * Aggregates per-region ASM profiles from all processed regions and identifies
 * subclone groups based on epigenetic similarity patterns. Uses a combination
 * of LOH status, HP ASM, and sample ASM to classify regions into subclones.
 *
 * Strategy:
 * 1. Extract valid RegionAsmProfiles from RegionResults
 * 2. Stratify by LOH status (LOH vs non-LOH have fundamentally different ASM patterns)
 * 3. Within each stratum, cluster by (hp_asm_delta, sample_asm_delta, hp_fine_f)
 * 4. Assign subclone IDs and compute per-subclone summary statistics
 */
class SubcloneAnalyzer {
public:
    /**
     * @brief Analyze cross-region subclone structure.
     *
     * @param results Vector of RegionResults from process_all_regions()
     * @param min_regions Minimum regions required for analysis (default: 50)
     * @return SubcloneResult with assignments and summaries
     */
    SubcloneResult analyze(const std::vector<RegionResult>& results, int min_regions = 50) const;

    /**
     * @brief Write subclone analysis report to file.
     *
     * @param result SubcloneResult to report
     * @param output_path Path for output report file
     */
    static void write_report(const SubcloneResult& result, const std::string& output_path);

    /**
     * @brief Write per-region subclone assignments to TSV.
     *
     * @param profiles Region profiles with assignments
     * @param assignments Per-region subclone IDs
     * @param output_path Path for output TSV file
     */
    static void write_assignments_tsv(const std::vector<RegionAsmProfile>& profiles,
                                       const std::vector<int>& assignments,
                                       const std::string& output_path);

private:
    /**
     * @brief Extract valid RegionAsmProfiles from RegionResults.
     *
     * Filters out failed regions, regions with too few reads/CpGs.
     */
    std::vector<RegionAsmProfile> extract_profiles(const std::vector<RegionResult>& results) const;

    /**
     * @brief Simple stratification-based subclone assignment.
     *
     * Groups regions by LOH status and ASM pattern:
     * - Group 0: Non-LOH, low ASM (normal diploid)
     * - Group 1: Non-LOH, high HP ASM (epigenetic heterogeneity)
     * - Group 2: LOH regions (LOH-driven pattern)
     * - Group 3: High sample ASM (tumor-specific changes)
     *
     * @param profiles Valid region profiles
     * @return Per-profile subclone assignments
     */
    std::vector<int> assign_subclones(const std::vector<RegionAsmProfile>& profiles) const;

    /**
     * @brief Compute summary statistics for each subclone group.
     */
    std::vector<SubcloneSummary> compute_summaries(const std::vector<RegionAsmProfile>& profiles,
                                                    const std::vector<int>& assignments,
                                                    int n_subclones) const;
};

}  // namespace InterSubMod
