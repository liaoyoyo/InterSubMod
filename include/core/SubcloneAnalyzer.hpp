#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace InterSubMod {

// Forward declaration
struct RegionResult;

constexpr int REGION_STRATIFICATION_SCHEMA_VERSION = 1;
constexpr std::size_t REGION_STRATUM_SLOT_COUNT = 4;
constexpr int DEFAULT_MIN_REGIONS_FOR_STRATIFICATION = 50;

/**
 * @brief Fixed, non-biological region strata used by the cross-region analysis.
 *
 * These values describe region-level ASM patterns. They must never be interpreted
 * as cellular subclones.
 */
enum class RegionStratum : int {
    BaselineLowAsm = 0,
    HighHpAsm = 1,
    LohFlagged = 2,
    HighSampleAsm = 3,
};

enum class RegionStratificationStatus {
    Valid,
    InsufficientRegions,
    NotApplicableTumorOnly,
    Failed,
};

const char* region_stratum_label(RegionStratum stratum);
const char* region_stratum_reason(RegionStratum stratum);
const char* region_stratification_status_string(RegionStratificationStatus status);

/**
 * @brief Per-region ASM feature profile for region stratification.
 */
struct RegionAsmProfile {
    std::size_t result_index = 0;
    int region_id = -1;
    std::string chrom;
    uint32_t pos = 0;

    double hp_asm_delta = 0.0;
    double sample_asm_delta = 0.0;
    double hp_fine_f = 0.0;
    double allele_delta = 0.0;
    double hp_ratio = 0.0;
    double coverage_multiple = 0.0;
    bool potential_loh = false;
    bool loh_bed_overlap = false;
    std::string loh_source;

    // verification_class is retained as the exact legacy-v1 compatibility alias.
    std::string verification_class;
    std::string verification_class_current;
    std::string verification_class_legacy;
    int verification_schema_version = 0;

    float quality_score = 0.0F;
    int num_reads = 0;
    int num_cpgs = 0;
};

struct RegionStratumAssignment {
    std::size_t result_index = 0;
    int region_id = -1;
    RegionStratum stratum = RegionStratum::BaselineLowAsm;
};

struct AssignmentValidationResult {
    bool valid = false;
    std::string reason;
};

/**
 * @brief Summary statistics for one of the four fixed region-stratum slots.
 */
struct SubcloneSummary {
    // Deprecated name retained for source compatibility; value is the region-stratum ID.
    int subclone_id = -1;
    int n_regions = 0;

    double mean_hp_asm = 0.0;
    double mean_sample_asm = 0.0;
    double mean_hp_ratio = 0.0;
    double mean_hp_fine_f = 0.0;
    double mean_coverage_multiple = 0.0;

    int n_loh_bed = 0;
    int n_loh_ratio = 0;
    int n_loh_both = 0;

    // Deprecated legacy-v1 mirrors retained for source compatibility.
    int n_strong = 0;
    int n_subclone = 0;
    int n_weak = 0;
    int n_noise = 0;
    std::string dominant_verification_class;

    std::array<int, 11> verification_class_v2_counts{};
    int verification_class_v2_unknown = 0;
    std::array<int, 4> verification_class_legacy_counts{};
    int verification_class_legacy_unknown = 0;
    std::string dominant_verification_classes_v2;
    std::string dominant_verification_classes_legacy;
};

/**
 * @brief Result of fixed-slot cross-region ASM region stratification.
 */
struct SubcloneResult {
    // Deprecated name retained for source compatibility; this is occupied-stratum count.
    int n_subclones = 0;
    int n_occupied_region_strata = 0;
    std::vector<int> region_assignments;
    std::vector<RegionStratumAssignment> assignments;
    std::vector<SubcloneSummary> summaries;
    bool valid = false;
    RegionStratificationStatus status = RegionStratificationStatus::Failed;
    std::string reason = "RUN_IN_PROGRESS";
    int warning_count = 0;

    int total_regions_analyzed = 0;
    int total_regions_assigned = 0;
    double mean_hp_asm_all = 0.0;
    double mean_sample_asm_all = 0.0;
    double loh_fraction = 0.0;
};

bool is_region_stratification_eligible(const RegionResult& result);
RegionStratum assign_region_stratum(const RegionAsmProfile& profile);

/**
 * @brief Validate the complete assignment set without mutating RegionResult.
 */
AssignmentValidationResult validate_region_stratification_assignments(
    const std::vector<RegionResult>& results,
    const std::vector<RegionStratumAssignment>& assignments,
    std::size_t expected_assignment_count);

/**
 * @brief Validate all assignments and commit them only if the complete set is valid.
 */
AssignmentValidationResult commit_region_stratification_assignments(
    std::vector<RegionResult>& results,
    const std::vector<RegionStratumAssignment>& assignments,
    std::size_t expected_assignment_count);

/**
 * @brief Compatibility facade for the historical SubcloneAnalyzer API.
 *
 * The implementation now performs deterministic region stratification. The class
 * name is retained for one compatibility release.
 */
class SubcloneAnalyzer {
public:
    SubcloneResult analyze(const std::vector<RegionResult>& results,
                           int min_regions = DEFAULT_MIN_REGIONS_FOR_STRATIFICATION) const;

    static void write_report(const SubcloneResult& result, const std::string& output_path);

    static void write_assignments_tsv(const std::vector<RegionAsmProfile>& profiles,
                                      const std::vector<int>& assignments,
                                      const std::string& output_path);

private:
    std::vector<RegionAsmProfile> extract_profiles(const std::vector<RegionResult>& results) const;
    std::vector<RegionStratumAssignment> assign_strata(const std::vector<RegionAsmProfile>& profiles) const;
    std::vector<SubcloneSummary> compute_summaries(
        const std::vector<RegionAsmProfile>& profiles,
        const std::vector<RegionStratumAssignment>& assignments) const;
};

}  // namespace InterSubMod
