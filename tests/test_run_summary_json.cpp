// Tests for RegionProcessor::write_run_summary_json (added 2026-08-06).
//
// print_summary() writes the same statistics to the log, where they are lost
// once the terminal closes. run_summary.json persists them so a finished run
// can be audited or compared without re-parsing logs. These tests pin down the
// accumulation rules — in particular that failed regions contribute nothing.

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "core/Config.hpp"
#include "core/RegionProcessor.hpp"

using namespace InterSubMod;

namespace {

/// Minimal config pointing at a scratch output directory. The RegionProcessor
/// constructor only copies fields and parses the VCF stem, so no real BAM or
/// reference is required to exercise the summary writer.
Config make_config(const std::string& out_dir) {
    Config c;
    c.output_dir = out_dir;
    c.somatic_vcf_path = "dummy_cohort.vcf";
    c.threads = 4;
    c.compute_distance_matrix = true;
    return c;
}

RegionResult make_result(bool success, int reads, int cpgs, int fwd, int rev, int filtered, double ms,
                         int valid_pairs = 0, int invalid_pairs = 0, double avg_cov = 0.0) {
    RegionResult r;
    r.success = success;
    r.num_reads = reads;
    r.num_cpgs = cpgs;
    r.num_forward_reads = fwd;
    r.num_reverse_reads = rev;
    r.num_filtered_reads = filtered;
    r.elapsed_ms = ms;
    r.num_valid_pairs = valid_pairs;
    r.num_invalid_pairs = invalid_pairs;
    r.avg_common_coverage = avg_cov;
    return r;
}

std::string read_file(const std::string& path) {
    std::ifstream in(path);
    if (!in.good()) {
        return "";
    }
    std::stringstream buf;
    buf << in.rdbuf();
    return buf.str();
}

/// Removes the scratch directory and everything under it.
void cleanup(const std::string& dir) {
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

}  // namespace

TEST(RunSummaryJsonTest, CountsSucceededAndFailedRegionsSeparately) {
    const std::string dir = "run_summary_test_counts";
    cleanup(dir);

    Config config = make_config(dir);
    RegionProcessor processor(config);

    std::vector<RegionResult> results = {
        make_result(true, 100, 10, 60, 40, 5, 12.5),
        make_result(true, 200, 20, 120, 80, 7, 20.0),
        make_result(false, 999, 999, 999, 999, 999, 999.0),  // must contribute nothing
    };

    std::string err;
    ASSERT_TRUE(processor.write_run_summary_json(results, 500.0, err)) << err;

    const std::string json = read_file(dir + "/run_summary.json");
    ASSERT_FALSE(json.empty());

    EXPECT_NE(json.find("\"total\": 3"), std::string::npos);
    EXPECT_NE(json.find("\"succeeded\": 2"), std::string::npos);
    EXPECT_NE(json.find("\"failed\": 1"), std::string::npos);

    // 100 + 200, with the failed region's 999 excluded.
    EXPECT_NE(json.find("\"total\": 300"), std::string::npos);
    EXPECT_NE(json.find("\"forward_strand\": 180"), std::string::npos);
    EXPECT_NE(json.find("\"reverse_strand\": 120"), std::string::npos);
    EXPECT_NE(json.find("\"filtered_out\": 12"), std::string::npos);

    cleanup(dir);
}

TEST(RunSummaryJsonTest, RecordsWallClockAndThreadCount) {
    const std::string dir = "run_summary_test_timing";
    cleanup(dir);

    Config config = make_config(dir);
    config.threads = 8;
    RegionProcessor processor(config);

    std::vector<RegionResult> results = {make_result(true, 10, 1, 5, 5, 0, 40.0)};

    std::string err;
    ASSERT_TRUE(processor.write_run_summary_json(results, 1234.5, err)) << err;

    const std::string json = read_file(dir + "/run_summary.json");
    EXPECT_NE(json.find("\"wall_clock\": 1234.5"), std::string::npos);
    EXPECT_NE(json.find("\"summed_region_time\": 40"), std::string::npos);
    EXPECT_NE(json.find("\"threads\": 8"), std::string::npos);

    cleanup(dir);
}

TEST(RunSummaryJsonTest, AggregatesDistanceMatrixPairStatistics) {
    const std::string dir = "run_summary_test_distance";
    cleanup(dir);

    Config config = make_config(dir);
    RegionProcessor processor(config);

    std::vector<RegionResult> results = {
        make_result(true, 50, 5, 25, 25, 0, 10.0, /*valid=*/30, /*invalid=*/10, /*avg_cov=*/6.0),
        make_result(true, 50, 5, 25, 25, 0, 10.0, /*valid=*/50, /*invalid=*/10, /*avg_cov=*/8.0),
    };

    std::string err;
    ASSERT_TRUE(processor.write_run_summary_json(results, 100.0, err)) << err;

    const std::string json = read_file(dir + "/run_summary.json");
    EXPECT_NE(json.find("\"regions_with_matrix\": 2"), std::string::npos);
    EXPECT_NE(json.find("\"valid_pairs\": 80"), std::string::npos);
    EXPECT_NE(json.find("\"invalid_pairs_insufficient_overlap\": 20"), std::string::npos);
    // 80 / (80 + 20) = 0.8
    EXPECT_NE(json.find("\"valid_pair_ratio\": 0.8"), std::string::npos);
    // (6.0 + 8.0) / 2
    EXPECT_NE(json.find("\"mean_common_cpg_coverage\": 7"), std::string::npos);

    cleanup(dir);
}

TEST(RunSummaryJsonTest, HandlesEmptyResultSetWithoutDivisionByZero) {
    const std::string dir = "run_summary_test_empty";
    cleanup(dir);

    Config config = make_config(dir);
    RegionProcessor processor(config);

    std::vector<RegionResult> results;

    std::string err;
    ASSERT_TRUE(processor.write_run_summary_json(results, 0.0, err)) << err;

    const std::string json = read_file(dir + "/run_summary.json");
    EXPECT_NE(json.find("\"total\": 0"), std::string::npos);
    EXPECT_NE(json.find("\"succeeded\": 0"), std::string::npos);
    // Means must be 0, not nan/inf.
    EXPECT_EQ(json.find("nan"), std::string::npos);
    EXPECT_EQ(json.find("inf"), std::string::npos);

    cleanup(dir);
}

TEST(RunSummaryJsonTest, ReportsErrorWhenOutputDirCannotBeCreated) {
    // A regular file occupying the output directory path: directory creation
    // must fail and the failure must surface rather than be swallowed.
    const std::string blocker = "run_summary_blocking_file";
    {
        std::ofstream ofs(blocker);
        ofs << "not a directory";
    }

    Config config = make_config(blocker);
    RegionProcessor processor(config);

    std::vector<RegionResult> results = {make_result(true, 1, 1, 1, 0, 0, 1.0)};

    std::string err;
    EXPECT_FALSE(processor.write_run_summary_json(results, 1.0, err));
    EXPECT_FALSE(err.empty());

    std::remove(blocker.c_str());
}
