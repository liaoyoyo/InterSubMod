#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "core/RegionProcessor.hpp"

using namespace InterSubMod;

namespace {

RegionResult make_region(int num_reads, bool success = true) {
    RegionResult r;
    r.num_reads = num_reads;
    r.success = success;
    return r;
}

}  // namespace

TEST(CoverageMultipleTest, BasicDivision) {
    EXPECT_DOUBLE_EQ(compute_coverage_multiple(75, 75.0), 1.0);
    EXPECT_DOUBLE_EQ(compute_coverage_multiple(150, 75.0), 2.0);
    EXPECT_DOUBLE_EQ(compute_coverage_multiple(30, 60.0), 0.5);
}

TEST(CoverageMultipleTest, DefaultExpectedIs75) {
    EXPECT_DOUBLE_EQ(compute_coverage_multiple(75), 1.0);
    EXPECT_DOUBLE_EQ(compute_coverage_multiple(150), 2.0);
}

TEST(EstimateDiploidCoverageTest, FallbackOnInsufficientData) {
    std::vector<RegionResult> results;
    for (int i = 0; i < 50; ++i) {
        results.push_back(make_region(60));
    }
    EXPECT_DOUBLE_EQ(estimate_diploid_coverage(results), 75.0);
}

TEST(EstimateDiploidCoverageTest, FallbackOnEmptyInput) {
    std::vector<RegionResult> empty_results;
    EXPECT_DOUBLE_EQ(estimate_diploid_coverage(empty_results), 75.0);
}

TEST(EstimateDiploidCoverageTest, FallbackOnLowReadCounts) {
    // All regions have num_reads <= 5, so none enter the estimation
    std::vector<RegionResult> results;
    for (int i = 0; i < 200; ++i) {
        results.push_back(make_region(3));
    }
    EXPECT_DOUBLE_EQ(estimate_diploid_coverage(results), 75.0);
}

TEST(EstimateDiploidCoverageTest, ModeNear60) {
    // Construct a distribution where mode sits around NumReads=60
    std::vector<RegionResult> results;
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(60.0, 8.0);
    for (int i = 0; i < 1000; ++i) {
        int nr = static_cast<int>(std::max(10.0, dist(rng)));
        results.push_back(make_region(nr));
    }
    double est = estimate_diploid_coverage(results);
    // KDE should land close to mode of 60; allow wide tolerance due to p80 truncation + histogram binning
    EXPECT_GT(est, 45.0);
    EXPECT_LT(est, 80.0);
}

TEST(EstimateDiploidCoverageTest, ModeNear90) {
    // Higher sequencing depth sample
    std::vector<RegionResult> results;
    std::mt19937 rng(123);
    std::normal_distribution<double> dist(90.0, 10.0);
    for (int i = 0; i < 1000; ++i) {
        int nr = static_cast<int>(std::max(15.0, dist(rng)));
        results.push_back(make_region(nr));
    }
    double est = estimate_diploid_coverage(results);
    EXPECT_GT(est, 70.0);
    EXPECT_LT(est, 110.0);
}

TEST(EstimateDiploidCoverageTest, FailedRegionsSkipped) {
    // Mix of successful regions at mode 60 + many failed regions with extreme values
    std::vector<RegionResult> results;
    std::mt19937 rng(7);
    std::normal_distribution<double> dist(60.0, 8.0);
    for (int i = 0; i < 500; ++i) {
        int nr = static_cast<int>(std::max(10.0, dist(rng)));
        results.push_back(make_region(nr, true));
    }
    // Add failed regions with extreme values that should be ignored
    for (int i = 0; i < 500; ++i) {
        results.push_back(make_region(200, false));
    }
    double est = estimate_diploid_coverage(results);
    EXPECT_GT(est, 45.0);
    EXPECT_LT(est, 80.0);
}

TEST(RegionResultTest, DiploidCoverageUsedDefault) {
    RegionResult r;
    EXPECT_DOUBLE_EQ(r.diploid_coverage_used, 75.0);
}
