#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <string>
#include <vector>

#include "core/RegionProcessor.hpp"
#include "core/SubcloneAnalyzer.hpp"

using namespace InterSubMod;

namespace {

RegionResult eligible_region(int id, double hp_delta = 0.0, double sample_delta = 0.0) {
    RegionResult result;
    result.region_id = id;
    result.success = true;
    result.num_reads = 10;
    result.num_cpgs = 3;
    result.significance_computed = true;
    result.hp_merged_delta = hp_delta;
    result.sample_asm_delta = sample_delta;
    result.verification_class = "Noise_Uncorrelated";
    result.verification_class_legacy = "Noise";
    result.loh_source = "none";
    return result;
}

std::vector<RegionResult> make_regions(int count) {
    std::vector<RegionResult> results;
    results.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) results.push_back(eligible_region(i));
    return results;
}

int sum_current(const SubcloneSummary& summary) {
    return std::accumulate(summary.verification_class_v2_counts.begin(),
                           summary.verification_class_v2_counts.end(), 0) +
           summary.verification_class_v2_unknown;
}

int sum_legacy(const SubcloneSummary& summary) {
    return std::accumulate(summary.verification_class_legacy_counts.begin(),
                           summary.verification_class_legacy_counts.end(), 0) +
           summary.verification_class_legacy_unknown;
}

}  // namespace

TEST(RegionStratifierTest, EmptyGuardWinsEvenWhenMinimumIsZero) {
    SubcloneAnalyzer analyzer;
    const auto result = analyzer.analyze({}, 0);
    EXPECT_FALSE(result.valid);
    EXPECT_EQ(result.status, RegionStratificationStatus::InsufficientRegions);
    EXPECT_EQ(result.reason, "NO_ELIGIBLE_REGIONS");
    EXPECT_EQ(result.total_regions_assigned, 0);
    EXPECT_EQ(result.n_occupied_region_strata, 0);
    ASSERT_EQ(result.summaries.size(), REGION_STRATUM_SLOT_COUNT);
}

TEST(RegionStratifierTest, EligibilityBoundaryAtFortyNineAndFiftyIsExact) {
    SubcloneAnalyzer analyzer;
    const auto below = analyzer.analyze(make_regions(49));
    EXPECT_FALSE(below.valid);
    EXPECT_EQ(below.reason, "BELOW_MIN_REGIONS");
    EXPECT_EQ(below.total_regions_assigned, 0);
    EXPECT_EQ(below.n_occupied_region_strata, 0);

    const auto exact = analyzer.analyze(make_regions(50));
    EXPECT_TRUE(exact.valid);
    EXPECT_EQ(exact.total_regions_assigned, 50);
    EXPECT_EQ(exact.assignments.size(), 50U);
    EXPECT_EQ(exact.n_occupied_region_strata, 1);
    EXPECT_EQ(exact.n_subclones, exact.n_occupied_region_strata);
}

TEST(RegionStratifierTest, OccupiedCountIsDistinctAndStorageAlwaysHasFourSlots) {
    SubcloneAnalyzer analyzer;
    auto only_three = std::vector<RegionResult>{eligible_region(0, 0.0, 0.11)};
    auto result = analyzer.analyze(only_three, 1);
    ASSERT_TRUE(result.valid);
    EXPECT_EQ(result.n_occupied_region_strata, 1);
    EXPECT_EQ(result.n_subclones, 1);
    ASSERT_EQ(result.summaries.size(), 4U);
    EXPECT_EQ(result.summaries[3].n_regions, 1);

    auto zero_and_three = std::vector<RegionResult>{eligible_region(0), eligible_region(1, 0.0, 0.11)};
    result = analyzer.analyze(zero_and_three, 1);
    ASSERT_TRUE(result.valid);
    EXPECT_EQ(result.n_occupied_region_strata, 2);
    ASSERT_EQ(result.summaries.size(), 4U);
    EXPECT_EQ(result.summaries[0].n_regions, 1);
    EXPECT_EQ(result.summaries[3].n_regions, 1);

    auto only_loh = std::vector<RegionResult>{eligible_region(2, 0.20, 0.20)};
    only_loh[0].potential_loh = true;
    only_loh[0].loh_source = "ratio_only";
    result = analyzer.analyze(only_loh, 1);
    ASSERT_TRUE(result.valid);
    EXPECT_EQ(result.n_occupied_region_strata, 1);
    EXPECT_EQ(result.summaries[2].n_regions, 1);
}

TEST(RegionStratifierTest, PrecedenceAndStrictThresholdsRemainUnchanged) {
    RegionAsmProfile profile;
    profile.hp_asm_delta = 0.05;
    profile.sample_asm_delta = 0.10;
    EXPECT_EQ(assign_region_stratum(profile), RegionStratum::BaselineLowAsm);

    profile.hp_asm_delta = 0.050001;
    EXPECT_EQ(assign_region_stratum(profile), RegionStratum::HighHpAsm);
    profile.sample_asm_delta = 0.100001;
    EXPECT_EQ(assign_region_stratum(profile), RegionStratum::HighSampleAsm);
    profile.potential_loh = true;
    EXPECT_EQ(assign_region_stratum(profile), RegionStratum::LohFlagged);
}

TEST(RegionStratifierTest, SharedEligibilityPreservesOriginalResultIndices) {
    std::vector<RegionResult> results;
    results.push_back(eligible_region(10));
    results.back().num_reads = 9;
    results.push_back(eligible_region(11));
    results.push_back(eligible_region(12));
    results.back().num_cpgs = 2;
    results.push_back(eligible_region(13));
    results.push_back(eligible_region(14));
    results.back().significance_computed = false;

    EXPECT_FALSE(is_region_stratification_eligible(results[0]));
    EXPECT_TRUE(is_region_stratification_eligible(results[1]));
    EXPECT_FALSE(is_region_stratification_eligible(results[2]));
    EXPECT_TRUE(is_region_stratification_eligible(results[3]));
    EXPECT_FALSE(is_region_stratification_eligible(results[4]));

    SubcloneAnalyzer analyzer;
    const auto analysis = analyzer.analyze(results, 2);
    ASSERT_TRUE(analysis.valid);
    ASSERT_EQ(analysis.assignments.size(), 2U);
    EXPECT_EQ(analysis.assignments[0].result_index, 1U);
    EXPECT_EQ(analysis.assignments[0].region_id, 11);
    EXPECT_EQ(analysis.assignments[1].result_index, 3U);
    EXPECT_EQ(analysis.assignments[1].region_id, 13);
}

TEST(RegionStratifierTest, ValidateAllThenCommitAllPreventsPartialWriteback) {
    auto results = make_regions(2);
    const std::vector<RegionStratumAssignment> duplicate = {
        {0, 0, RegionStratum::BaselineLowAsm}, {0, 0, RegionStratum::HighHpAsm}};
    auto validation = commit_region_stratification_assignments(results, duplicate, 2);
    EXPECT_FALSE(validation.valid);
    EXPECT_EQ(validation.reason, "DUPLICATE_RESULT_INDEX");
    EXPECT_EQ(results[0].region_stratum_id, -1);
    EXPECT_EQ(results[1].region_stratum_id, -1);

    const std::vector<RegionStratumAssignment> mismatch = {
        {0, 999, RegionStratum::BaselineLowAsm}, {1, 1, RegionStratum::HighHpAsm}};
    validation = commit_region_stratification_assignments(results, mismatch, 2);
    EXPECT_FALSE(validation.valid);
    EXPECT_EQ(validation.reason, "REGION_ID_MISMATCH");
    EXPECT_EQ(results[0].region_stratum_id, -1);
    EXPECT_EQ(results[1].region_stratum_id, -1);

    const std::vector<RegionStratumAssignment> out_of_range = {
        {0, 0, RegionStratum::BaselineLowAsm}, {2, 1, RegionStratum::HighHpAsm}};
    validation = commit_region_stratification_assignments(results, out_of_range, 2);
    EXPECT_FALSE(validation.valid);
    EXPECT_EQ(validation.reason, "RESULT_INDEX_OUT_OF_RANGE");
    EXPECT_EQ(results[0].region_stratum_id, -1);

    const std::vector<RegionStratumAssignment> valid = {
        {0, 0, RegionStratum::BaselineLowAsm}, {1, 1, RegionStratum::HighHpAsm}};
    validation = commit_region_stratification_assignments(results, valid, 2);
    ASSERT_TRUE(validation.valid);
    EXPECT_EQ(results[0].region_stratum_id, 0);
    EXPECT_EQ(results[1].region_stratum_id, 1);
    EXPECT_EQ(results[0].subclone_id, results[0].region_stratum_id);
    EXPECT_EQ(results[1].subclone_id, results[1].region_stratum_id);
}

TEST(RegionStratifierTest, CurrentAndLegacyLeafCountersConserveUnknownsAndTies) {
    const std::array<const char*, 11> current_classes = {
        "Strong_Bidirectional", "ClusterFirstOnly", "LOH-Structure", "MultiGroupNoLabel",
        "LabelShift", "PermanovaLocation", "StructureNoLabel", "DispersionStructure",
        "Noise_Uniform", "Noise_Chaotic", "Noise_Uncorrelated"};
    const std::array<const char*, 4> legacy_classes = {"Strong", "Subclone", "Weak", "Noise"};

    std::vector<RegionResult> results;
    for (std::size_t i = 0; i < current_classes.size(); ++i) {
        auto result = eligible_region(static_cast<int>(i));
        result.verification_class = current_classes[i];
        result.verification_class_legacy = legacy_classes[i % legacy_classes.size()];
        results.push_back(std::move(result));
    }
    auto unknown = eligible_region(11);
    unknown.verification_class = "FutureCurrentClass";
    unknown.verification_class_legacy = "FutureLegacyClass";
    results.push_back(std::move(unknown));

    SubcloneAnalyzer analyzer;
    const auto analysis = analyzer.analyze(results, 1);
    ASSERT_TRUE(analysis.valid);
    ASSERT_EQ(analysis.summaries.size(), 4U);
    const auto& summary = analysis.summaries[0];
    EXPECT_EQ(summary.n_regions, 12);
    for (int count : summary.verification_class_v2_counts) EXPECT_EQ(count, 1);
    EXPECT_EQ(summary.verification_class_v2_unknown, 1);
    EXPECT_EQ(sum_current(summary), summary.n_regions);
    EXPECT_EQ(sum_legacy(summary), summary.n_regions);
    EXPECT_EQ(summary.verification_class_legacy_counts[0], 3);
    EXPECT_EQ(summary.verification_class_legacy_counts[1], 3);
    EXPECT_EQ(summary.verification_class_legacy_counts[2], 3);
    EXPECT_EQ(summary.verification_class_legacy_counts[3], 2);
    EXPECT_EQ(summary.verification_class_legacy_unknown, 1);
    EXPECT_EQ(analysis.warning_count, 2);
    EXPECT_EQ(analysis.reason, "VALID_WITH_WARNINGS");
    EXPECT_EQ(summary.dominant_verification_classes_v2,
              "Strong_Bidirectional;ClusterFirstOnly;LOH-Structure;MultiGroupNoLabel;LabelShift;"
              "PermanovaLocation;StructureNoLabel;DispersionStructure;Noise_Uniform;Noise_Chaotic;"
              "Noise_Uncorrelated;UnknownCurrentClass");
    EXPECT_EQ(summary.dominant_verification_classes_legacy, "Strong;Subclone;Weak");
}

TEST(RegionStratifierTest, CurrentClassRequiresExplicitSchemaVersionTwo) {
    auto result = eligible_region(0);
    result.verification_class = "Strong_Bidirectional";
    result.verification_schema_version = 1;
    result.verification_class_legacy = "Strong";

    SubcloneAnalyzer analyzer;
    const auto analysis = analyzer.analyze({result}, 1);
    ASSERT_TRUE(analysis.valid);
    const auto& summary = analysis.summaries[0];
    EXPECT_EQ(summary.verification_class_v2_counts[0], 0);
    EXPECT_EQ(summary.verification_class_v2_unknown, 1);
    EXPECT_EQ(summary.verification_class_legacy_counts[0], 1);
    EXPECT_EQ(summary.dominant_verification_classes_v2, "UnknownCurrentClass");
    EXPECT_EQ(analysis.warning_count, 1);
}
