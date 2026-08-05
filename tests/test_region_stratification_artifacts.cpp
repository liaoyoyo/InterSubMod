#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

#include "core/RegionProcessor.hpp"

using namespace InterSubMod;

namespace {

std::string unique_output_dir() {
    const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
    return "/tmp/intersubmod_region_artifacts_" + std::to_string(static_cast<long long>(::getpid())) +
           "_" + std::to_string(tick);
}

std::string read_all(const std::string& path) {
    std::ifstream input(path);
    std::ostringstream content;
    content << input.rdbuf();
    return content.str();
}

int newline_count(const std::string& text) {
    return static_cast<int>(std::count(text.begin(), text.end(), '\n'));
}

}  // namespace

TEST(RegionStratificationArtifactsTest, StatusHasExactHeaderAndRfc3339UtcTimestamp) {
    const std::string output_dir = unique_output_dir();
    RegionStratificationStatusRecord status;
    status.status = RegionStratificationStatus::Failed;
    status.reason = "RUN_IN_PROGRESS";
    std::string error;
    ASSERT_TRUE(write_region_stratification_status_atomic(output_dir, status, error)) << error;

    const std::string content = read_all(output_dir + "/region_stratification_status.tsv");
    EXPECT_EQ(content.substr(0, content.find('\n')),
              "status\treason\tschema_version\teligible_region_count\tmin_regions_required\t"
              "assignment_count\tn_occupied_region_strata\twarning_count\tgenerated_at");
    EXPECT_NE(content.find("FAILED\tRUN_IN_PROGRESS\t1\t0\t50\t0\t0\t0\t"), std::string::npos);
    const std::regex timestamp("[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z");
    EXPECT_TRUE(std::regex_search(content, timestamp));
}

TEST(RegionStratificationArtifactsTest, ValidArtifactsUseStableKeyAndOnlyOccupiedSections) {
    const std::string output_dir = unique_output_dir();
    RegionStratificationArtifactRow row;
    row.result_index = 7;
    row.region_id = 42;
    row.chromosome = "chr17";
    row.position = 7579472;
    row.reference = 'G';
    row.alternate = 'A';
    row.stratum_id = 3;
    row.stratum_label = "HighSampleAsm";
    row.stratum_reason = "HIGH_SAMPLE_ASM";

    std::vector<SubcloneSummary> summaries(REGION_STRATUM_SLOT_COUNT);
    for (std::size_t i = 0; i < summaries.size(); ++i) summaries[i].subclone_id = static_cast<int>(i);
    summaries[3].n_regions = 1;
    summaries[3].verification_class_v2_counts[0] = 1;
    summaries[3].verification_class_legacy_counts[0] = 1;
    summaries[3].dominant_verification_classes_v2 = "Strong_Bidirectional";
    summaries[3].dominant_verification_classes_legacy = "Strong";

    RegionStratificationStatusRecord status;
    status.status = RegionStratificationStatus::Valid;
    status.reason = "OK";
    status.eligible_region_count = 1;
    status.min_regions_required = 1;
    status.assignment_count = 1;
    status.n_occupied_region_strata = 1;

    std::string error;
    ASSERT_TRUE(write_region_stratification_artifacts_atomic(
        output_dir, {row}, summaries, status, error)) << error;

    const std::string assignments = read_all(output_dir + "/region_stratification.tsv");
    EXPECT_EQ(newline_count(assignments), 2);
    EXPECT_NE(assignments.find("42\tchr17\t7579472\tG\tA\t3\tHighSampleAsm\tHIGH_SAMPLE_ASM\t1"),
              std::string::npos);

    const std::string summary = read_all(output_dir + "/region_stratification_summary.txt");
    EXPECT_EQ(summary.find("Cross-region ASM region stratification"), 0U);
    EXPECT_NE(summary.find("Status: VALID"), std::string::npos);
    EXPECT_NE(summary.find("RegionStratum 3: HighSampleAsm"), std::string::npos);
    EXPECT_EQ(summary.find("RegionStratum 0:"), std::string::npos);

    EXPECT_EQ(read_all(output_dir + "/subclone_structure.txt"),
              "DEPRECATED: This file does not contain cellular subclones.\n"
              "Canonical report: region_stratification_summary.txt\n"
              "Status: VALID\n"
              "Reason: OK\n");

    for (const auto& entry : std::filesystem::directory_iterator(output_dir)) {
        EXPECT_EQ(entry.path().filename().string().find(".tmp."), std::string::npos);
    }
}

TEST(RegionStratificationArtifactsTest, NonValidRunsOverwritePriorValidArtifactsWithoutStaleRows) {
    const std::string output_dir = unique_output_dir();
    RegionStratificationArtifactRow stale_row;
    stale_row.region_id = 1;
    stale_row.chromosome = "chr1";
    stale_row.position = 100;
    stale_row.reference = 'C';
    stale_row.alternate = 'T';
    stale_row.stratum_id = 0;
    stale_row.stratum_label = "BaselineLowAsm";
    stale_row.stratum_reason = "BASELINE_LOW_ASM";
    std::vector<SubcloneSummary> summaries(REGION_STRATUM_SLOT_COUNT);
    for (std::size_t i = 0; i < summaries.size(); ++i) summaries[i].subclone_id = static_cast<int>(i);
    summaries[0].n_regions = 1;

    RegionStratificationStatusRecord valid;
    valid.status = RegionStratificationStatus::Valid;
    valid.reason = "OK";
    valid.eligible_region_count = 1;
    valid.min_regions_required = 1;
    valid.assignment_count = 1;
    valid.n_occupied_region_strata = 1;
    std::string error;
    ASSERT_TRUE(write_region_stratification_artifacts_atomic(
        output_dir, {stale_row}, summaries, valid, error)) << error;

    const std::vector<std::pair<RegionStratificationStatus, std::string>> states = {
        {RegionStratificationStatus::InsufficientRegions, "BELOW_MIN_REGIONS"},
        {RegionStratificationStatus::NotApplicableTumorOnly, "NORMAL_BAM_REQUIRED"},
        {RegionStratificationStatus::Failed, "ASSIGNMENT_VALIDATION_FAILED"}};
    for (const auto& [state, reason] : states) {
        RegionStratificationStatusRecord status;
        status.status = state;
        status.reason = reason;
        ASSERT_TRUE(write_region_stratification_artifacts_atomic(
            output_dir, {}, summaries, status, error)) << error;

        const std::string assignments = read_all(output_dir + "/region_stratification.tsv");
        EXPECT_EQ(newline_count(assignments), 1);
        EXPECT_EQ(assignments.find("chr1"), std::string::npos);
        const std::string summary = read_all(output_dir + "/region_stratification_summary.txt");
        EXPECT_EQ(summary.find("RegionStratum 0:"), std::string::npos);
        EXPECT_NE(summary.find("Reason: " + reason), std::string::npos);
        const std::string stub = read_all(output_dir + "/subclone_structure.txt");
        EXPECT_EQ(stub.find("DEPRECATED: This file does not contain cellular subclones."), 0U);
        EXPECT_NE(stub.find("Reason: " + reason), std::string::npos);
    }
}

TEST(RegionStratificationArtifactsTest, ContradictoryStatusAndRowsFailBeforePublication) {
    const std::string output_dir = unique_output_dir();
    RegionStratificationStatusRecord status;
    status.status = RegionStratificationStatus::Valid;
    status.reason = "OK";
    status.eligible_region_count = 1;
    status.min_regions_required = 1;
    status.assignment_count = 1;
    status.n_occupied_region_strata = 1;

    std::vector<SubcloneSummary> summaries(REGION_STRATUM_SLOT_COUNT);
    for (std::size_t i = 0; i < summaries.size(); ++i) summaries[i].subclone_id = static_cast<int>(i);
    std::string error;
    EXPECT_FALSE(write_region_stratification_artifacts_atomic(
        output_dir, {}, summaries, status, error));
    EXPECT_EQ(error, "REGION_ARTIFACT_ROW_COUNT_MISMATCH");
    EXPECT_FALSE(std::filesystem::exists(output_dir));

    status.status = RegionStratificationStatus::Failed;
    EXPECT_FALSE(write_region_stratification_status_atomic(output_dir, status, error));
    EXPECT_EQ(error, "NONVALID_REGION_STATUS_HAS_ASSIGNMENTS");
    EXPECT_FALSE(std::filesystem::exists(output_dir));
}

TEST(RegionStratificationArtifactsTest, AtomicRenameFailureRemovesOwnedTempFile) {
    const std::string output_dir = unique_output_dir();
    const std::string final_path = output_dir + "/region_stratification_status.tsv";
    ASSERT_TRUE(std::filesystem::create_directories(final_path));

    RegionStratificationStatusRecord status;
    status.status = RegionStratificationStatus::Failed;
    status.reason = "RUN_IN_PROGRESS";
    std::string error;
    EXPECT_FALSE(write_region_stratification_status_atomic(output_dir, status, error));
    EXPECT_EQ(error.find("ATOMIC_RENAME_FAILED:"), 0U);

    for (const auto& entry : std::filesystem::directory_iterator(output_dir)) {
        EXPECT_EQ(entry.path().filename().string().find(".tmp."), std::string::npos);
    }
    EXPECT_TRUE(std::filesystem::is_directory(final_path));
}
