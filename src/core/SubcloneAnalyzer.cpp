#include "core/SubcloneAnalyzer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_set>

#include "core/RegionProcessor.hpp"

namespace InterSubMod {
namespace {

constexpr double HIGH_HP_ASM_THRESHOLD = 0.05;
constexpr double HIGH_SAMPLE_ASM_THRESHOLD = 0.10;
constexpr int MIN_READS_FOR_PROFILE = 10;
constexpr int MIN_CPGS_FOR_PROFILE = 3;

constexpr std::array<const char*, REGION_STRATUM_SLOT_COUNT> REGION_STRATUM_LABELS = {
    "BaselineLowAsm", "HighHpAsm", "LohFlagged", "HighSampleAsm"};
constexpr std::array<const char*, REGION_STRATUM_SLOT_COUNT> REGION_STRATUM_REASONS = {
    "BASELINE_LOW_ASM", "HIGH_HP_ASM", "LOH_FLAGGED", "HIGH_SAMPLE_ASM"};

constexpr std::array<const char*, 11> VERIFICATION_V2_CLASSES = {
    "Strong_Bidirectional", "ClusterFirstOnly", "LOH-Structure", "MultiGroupNoLabel",
    "LabelShift", "PermanovaLocation", "StructureNoLabel", "DispersionStructure",
    "Noise_Uniform", "Noise_Chaotic", "Noise_Uncorrelated"};
constexpr std::array<const char*, 4> VERIFICATION_LEGACY_CLASSES = {
    "Strong", "Subclone", "Weak", "Noise"};

template <std::size_t N>
int class_index(const std::string& value, const std::array<const char*, N>& values) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (value == values[i]) return static_cast<int>(i);
    }
    return -1;
}

template <std::size_t N>
std::string dominant_classes(const std::array<int, N>& counts,
                             const std::array<const char*, N>& names,
                             int unknown_count,
                             const char* unknown_name) {
    const int max_count = std::max(*std::max_element(counts.begin(), counts.end()), unknown_count);
    if (max_count <= 0) return "";

    std::ostringstream joined;
    bool first = true;
    for (std::size_t i = 0; i < counts.size(); ++i) {
        if (counts[i] != max_count) continue;
        if (!first) joined << ";";
        joined << names[i];
        first = false;
    }
    if (unknown_count == max_count) {
        if (!first) joined << ";";
        joined << unknown_name;
    }
    return joined.str();
}

bool valid_stratum(RegionStratum stratum) {
    const int value = static_cast<int>(stratum);
    return value >= 0 && value < static_cast<int>(REGION_STRATUM_SLOT_COUNT);
}

}  // namespace

const char* region_stratum_label(RegionStratum stratum) {
    if (!valid_stratum(stratum)) return "Unknown";
    return REGION_STRATUM_LABELS[static_cast<std::size_t>(static_cast<int>(stratum))];
}

const char* region_stratum_reason(RegionStratum stratum) {
    if (!valid_stratum(stratum)) return "INVALID_STRATUM";
    return REGION_STRATUM_REASONS[static_cast<std::size_t>(static_cast<int>(stratum))];
}

const char* region_stratification_status_string(RegionStratificationStatus status) {
    switch (status) {
        case RegionStratificationStatus::Valid:
            return "VALID";
        case RegionStratificationStatus::InsufficientRegions:
            return "INSUFFICIENT_REGIONS";
        case RegionStratificationStatus::NotApplicableTumorOnly:
            return "NOT_APPLICABLE_TUMOR_ONLY";
        case RegionStratificationStatus::Failed:
            return "FAILED";
    }
    return "FAILED";
}

bool is_region_stratification_eligible(const RegionResult& result) {
    return result.success && result.num_reads >= MIN_READS_FOR_PROFILE &&
           result.num_cpgs >= MIN_CPGS_FOR_PROFILE && result.significance_computed;
}

RegionStratum assign_region_stratum(const RegionAsmProfile& profile) {
    // Precedence is part of the schema contract: LOH > sample ASM > HP ASM > baseline.
    if (profile.potential_loh || profile.loh_bed_overlap) return RegionStratum::LohFlagged;
    if (std::abs(profile.sample_asm_delta) > HIGH_SAMPLE_ASM_THRESHOLD) {
        return RegionStratum::HighSampleAsm;
    }
    if (std::abs(profile.hp_asm_delta) > HIGH_HP_ASM_THRESHOLD) return RegionStratum::HighHpAsm;
    return RegionStratum::BaselineLowAsm;
}

AssignmentValidationResult validate_region_stratification_assignments(
    const std::vector<RegionResult>& results,
    const std::vector<RegionStratumAssignment>& assignments,
    std::size_t expected_assignment_count) {
    if (assignments.size() != expected_assignment_count) {
        return {false, "ASSIGNMENT_COUNT_MISMATCH"};
    }

    std::size_t eligible_count = 0;
    for (const auto& result : results) {
        if (is_region_stratification_eligible(result)) ++eligible_count;
    }
    if (eligible_count != expected_assignment_count) {
        return {false, "ELIGIBLE_REGION_COUNT_MISMATCH"};
    }

    std::unordered_set<std::size_t> seen_indices;
    seen_indices.reserve(assignments.size());
    for (const auto& assignment : assignments) {
        if (assignment.result_index >= results.size()) return {false, "RESULT_INDEX_OUT_OF_RANGE"};
        if (!seen_indices.insert(assignment.result_index).second) {
            return {false, "DUPLICATE_RESULT_INDEX"};
        }
        const auto& result = results[assignment.result_index];
        if (!is_region_stratification_eligible(result)) return {false, "INELIGIBLE_ASSIGNMENT"};
        if (result.region_id != assignment.region_id) return {false, "REGION_ID_MISMATCH"};
        if (!valid_stratum(assignment.stratum)) return {false, "INVALID_STRATUM"};
    }

    for (std::size_t i = 0; i < results.size(); ++i) {
        if (is_region_stratification_eligible(results[i]) && seen_indices.count(i) == 0) {
            return {false, "ELIGIBLE_REGION_UNASSIGNED"};
        }
    }
    return {true, "OK"};
}

AssignmentValidationResult commit_region_stratification_assignments(
    std::vector<RegionResult>& results,
    const std::vector<RegionStratumAssignment>& assignments,
    std::size_t expected_assignment_count) {
    const auto validation =
        validate_region_stratification_assignments(results, assignments, expected_assignment_count);
    if (!validation.valid) return validation;

    // Mutation begins only after the complete assignment set has passed validation.
    for (const auto& assignment : assignments) {
        auto& result = results[assignment.result_index];
        const int stratum_id = static_cast<int>(assignment.stratum);
        result.region_stratification_schema_version = REGION_STRATIFICATION_SCHEMA_VERSION;
        result.region_stratum_id = stratum_id;
        result.region_stratum_label = region_stratum_label(assignment.stratum);
        result.region_stratum_reason = region_stratum_reason(assignment.stratum);
        result.subclone_id = stratum_id;  // Exact one-release compatibility alias.
    }
    return {true, "OK"};
}

SubcloneResult SubcloneAnalyzer::analyze(const std::vector<RegionResult>& results, int min_regions) const {
    SubcloneResult result;
    result.summaries = compute_summaries({}, {});  // Always publish four fixed slots.

    const auto profiles = extract_profiles(results);
    result.total_regions_analyzed = static_cast<int>(profiles.size());

    // Empty is never valid, even when a caller explicitly sets min_regions to zero.
    if (profiles.empty()) {
        result.status = RegionStratificationStatus::InsufficientRegions;
        result.reason = "NO_ELIGIBLE_REGIONS";
        return result;
    }
    if (static_cast<int>(profiles.size()) < std::max(0, min_regions)) {
        result.status = RegionStratificationStatus::InsufficientRegions;
        result.reason = "BELOW_MIN_REGIONS";
        return result;
    }

    result.assignments = assign_strata(profiles);
    result.region_assignments.reserve(result.assignments.size());
    std::array<bool, REGION_STRATUM_SLOT_COUNT> occupied{};
    for (const auto& assignment : result.assignments) {
        const int id = static_cast<int>(assignment.stratum);
        result.region_assignments.push_back(id);
        occupied[static_cast<std::size_t>(id)] = true;
    }
    result.n_occupied_region_strata =
        static_cast<int>(std::count(occupied.begin(), occupied.end(), true));
    result.n_subclones = result.n_occupied_region_strata;
    result.total_regions_assigned = static_cast<int>(result.assignments.size());
    result.summaries = compute_summaries(profiles, result.assignments);

    double sum_hp = 0.0;
    double sum_sample = 0.0;
    int n_loh = 0;
    for (const auto& profile : profiles) {
        sum_hp += std::abs(profile.hp_asm_delta);
        sum_sample += std::abs(profile.sample_asm_delta);
        if (profile.potential_loh || profile.loh_bed_overlap) ++n_loh;
    }
    const double denominator = static_cast<double>(profiles.size());
    result.mean_hp_asm_all = sum_hp / denominator;
    result.mean_sample_asm_all = sum_sample / denominator;
    result.loh_fraction = static_cast<double>(n_loh) / denominator;

    for (const auto& summary : result.summaries) {
        result.warning_count += summary.verification_class_v2_unknown;
        result.warning_count += summary.verification_class_legacy_unknown;
    }
    result.status = RegionStratificationStatus::Valid;
    result.reason = result.warning_count == 0 ? "OK" : "VALID_WITH_WARNINGS";
    result.valid = true;
    return result;
}

std::vector<RegionAsmProfile> SubcloneAnalyzer::extract_profiles(
    const std::vector<RegionResult>& results) const {
    std::vector<RegionAsmProfile> profiles;
    profiles.reserve(results.size());

    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];
        if (!is_region_stratification_eligible(result)) continue;

        RegionAsmProfile profile;
        profile.result_index = i;
        profile.region_id = result.region_id;
        profile.hp_asm_delta = result.hp_merged_delta;
        profile.sample_asm_delta = result.sample_asm_delta;
        profile.hp_fine_f = result.hp_fine_f;
        profile.allele_delta = result.allele_delta;
        profile.hp_ratio = result.hp_ratio;
        profile.coverage_multiple = result.coverage_multiple;
        profile.potential_loh = result.potential_loh;
        profile.loh_bed_overlap = result.loh_bed_overlap;
        profile.loh_source = result.loh_source;
        profile.verification_class = result.verification_class_legacy;
        profile.verification_class_current = result.verification_class;
        profile.verification_class_legacy = result.verification_class_legacy;
        profile.verification_schema_version = result.verification_schema_version;
        profile.quality_score = result.quality_score;
        profile.num_reads = result.num_reads;
        profile.num_cpgs = result.num_cpgs;
        profiles.push_back(std::move(profile));
    }
    return profiles;
}

std::vector<RegionStratumAssignment> SubcloneAnalyzer::assign_strata(
    const std::vector<RegionAsmProfile>& profiles) const {
    std::vector<RegionStratumAssignment> assignments;
    assignments.reserve(profiles.size());
    for (const auto& profile : profiles) {
        assignments.push_back(
            {profile.result_index, profile.region_id, assign_region_stratum(profile)});
    }
    return assignments;
}

std::vector<SubcloneSummary> SubcloneAnalyzer::compute_summaries(
    const std::vector<RegionAsmProfile>& profiles,
    const std::vector<RegionStratumAssignment>& assignments) const {
    std::vector<SubcloneSummary> summaries(REGION_STRATUM_SLOT_COUNT);
    for (std::size_t i = 0; i < summaries.size(); ++i) summaries[i].subclone_id = static_cast<int>(i);

    for (std::size_t i = 0; i < profiles.size(); ++i) {
        if (i >= assignments.size() || !valid_stratum(assignments[i].stratum)) continue;
        auto& summary = summaries[static_cast<std::size_t>(static_cast<int>(assignments[i].stratum))];
        const auto& profile = profiles[i];
        ++summary.n_regions;
        summary.mean_hp_asm += std::abs(profile.hp_asm_delta);
        summary.mean_sample_asm += std::abs(profile.sample_asm_delta);
        summary.mean_hp_ratio += profile.hp_ratio;
        summary.mean_hp_fine_f += profile.hp_fine_f;
        summary.mean_coverage_multiple += profile.coverage_multiple;

        if (profile.loh_source == "bed_only") {
            ++summary.n_loh_bed;
        } else if (profile.loh_source == "ratio_only") {
            ++summary.n_loh_ratio;
        } else if (profile.loh_source == "both") {
            ++summary.n_loh_bed;
            ++summary.n_loh_ratio;
            ++summary.n_loh_both;
        }

        const int current_index =
            profile.verification_schema_version == VERIFICATION_SCHEMA_VERSION
                ? class_index(profile.verification_class_current, VERIFICATION_V2_CLASSES)
                : -1;
        if (current_index >= 0) {
            ++summary.verification_class_v2_counts[static_cast<std::size_t>(current_index)];
        } else {
            ++summary.verification_class_v2_unknown;
        }

        const int legacy_index = class_index(profile.verification_class_legacy, VERIFICATION_LEGACY_CLASSES);
        if (legacy_index >= 0) {
            ++summary.verification_class_legacy_counts[static_cast<std::size_t>(legacy_index)];
        } else {
            ++summary.verification_class_legacy_unknown;
        }
    }

    for (auto& summary : summaries) {
        if (summary.n_regions > 0) {
            const double denominator = static_cast<double>(summary.n_regions);
            summary.mean_hp_asm /= denominator;
            summary.mean_sample_asm /= denominator;
            summary.mean_hp_ratio /= denominator;
            summary.mean_hp_fine_f /= denominator;
            summary.mean_coverage_multiple /= denominator;
        }

        summary.n_strong = summary.verification_class_legacy_counts[0];
        summary.n_subclone = summary.verification_class_legacy_counts[1];
        summary.n_weak = summary.verification_class_legacy_counts[2];
        summary.n_noise = summary.verification_class_legacy_counts[3];
        summary.dominant_verification_classes_v2 =
            dominant_classes(summary.verification_class_v2_counts, VERIFICATION_V2_CLASSES,
                             summary.verification_class_v2_unknown, "UnknownCurrentClass");
        summary.dominant_verification_classes_legacy =
            dominant_classes(summary.verification_class_legacy_counts, VERIFICATION_LEGACY_CLASSES,
                             summary.verification_class_legacy_unknown, "UnknownLegacyClass");
        summary.dominant_verification_class = summary.dominant_verification_classes_legacy;
    }
    return summaries;
}

void SubcloneAnalyzer::write_report(const SubcloneResult& result, const std::string& output_path) {
    std::ofstream output(output_path);
    if (!output) {
        std::cerr << "[SubcloneAnalyzer] Failed to open compatibility report: " << output_path << "\n";
        return;
    }
    output << "Cross-region ASM region stratification\n";
    output << "Status: " << region_stratification_status_string(result.status) << "\n";
    output << "Reason: " << result.reason << "\n";
    if (!result.valid) return;
    output << "Occupied region strata: " << result.n_occupied_region_strata << "\n";
    for (const auto& summary : result.summaries) {
        if (summary.n_regions == 0) continue;
        output << "Region stratum " << summary.subclone_id << " ("
               << region_stratum_label(static_cast<RegionStratum>(summary.subclone_id)) << ")\n";
        output << "  Regions: " << summary.n_regions << "\n";
    }
}

void SubcloneAnalyzer::write_assignments_tsv(const std::vector<RegionAsmProfile>& profiles,
                                              const std::vector<int>& assignments,
                                              const std::string& output_path) {
    std::ofstream output(output_path);
    if (!output) {
        std::cerr << "[SubcloneAnalyzer] Failed to open compatibility assignments: " << output_path << "\n";
        return;
    }
    output << "RegionID\tRegionStratum_ID\tRegionStratum_Label\tRegionStratificationSchemaVersion\n";
    for (std::size_t i = 0; i < profiles.size(); ++i) {
        const int id = i < assignments.size() ? assignments[i] : -1;
        const char* label = id >= 0 && id < static_cast<int>(REGION_STRATUM_SLOT_COUNT)
                                ? region_stratum_label(static_cast<RegionStratum>(id))
                                : "Unassigned";
        output << profiles[i].region_id << "\t" << id << "\t" << label << "\t"
               << REGION_STRATIFICATION_SCHEMA_VERSION << "\n";
    }
}

}  // namespace InterSubMod
