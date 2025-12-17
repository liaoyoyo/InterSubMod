/**
 * @file test_distance_matrix.cpp
 * @brief Unit tests for DistanceMatrix and DistanceCalculator
 *
 * Tests cover:
 * 1. Various distance metrics (NHD, L1, L2, CORR, JACCARD)
 * 2. Missing value handling
 * 3. Minimum common coverage threshold
 * 4. Strand-specific distance calculation
 * 5. Edge cases (empty matrices, single read, etc.)
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

#include "core/DataStructs.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/MethylationMatrix.hpp"

using namespace InterSubMod;

// ============================================================================
// Test Fixtures
// ============================================================================

class DistanceMatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple 4x5 methylation matrix
        // Rows: 4 reads, Columns: 5 CpG sites
        //
        // Read 0 (FORWARD): [1,    1,    0,    0,    NaN]
        // Read 1 (FORWARD): [1,    0,    0,    NaN,  NaN]  -> Similar to Read 0
        // Read 2 (REVERSE): [0,    0,    1,    1,    1]    -> Different from Read 0
        // Read 3 (REVERSE): [NaN,  NaN,  1,    1,    1]    -> Similar to Read 2

        meth_mat.region_id = 0;
        meth_mat.read_ids = {0, 1, 2, 3};
        meth_mat.cpg_ids = {0, 1, 2, 3, 4};

        int n_reads = 4;
        int n_cpgs = 5;

        meth_mat.raw_matrix = Eigen::MatrixXd(n_reads, n_cpgs);
        meth_mat.binary_matrix = Eigen::MatrixXi(n_reads, n_cpgs);

        // Read 0: [1, 1, 0, 0, NaN]
        meth_mat.raw_matrix(0, 0) = 0.95;
        meth_mat.raw_matrix(0, 1) = 0.90;
        meth_mat.raw_matrix(0, 2) = 0.10;
        meth_mat.raw_matrix(0, 3) = 0.05;
        meth_mat.raw_matrix(0, 4) = NAN;

        meth_mat.binary_matrix(0, 0) = 1;
        meth_mat.binary_matrix(0, 1) = 1;
        meth_mat.binary_matrix(0, 2) = 0;
        meth_mat.binary_matrix(0, 3) = 0;
        meth_mat.binary_matrix(0, 4) = -1;

        // Read 1: [1, 0, 0, NaN, NaN]
        meth_mat.raw_matrix(1, 0) = 0.85;
        meth_mat.raw_matrix(1, 1) = 0.15;
        meth_mat.raw_matrix(1, 2) = 0.10;
        meth_mat.raw_matrix(1, 3) = NAN;
        meth_mat.raw_matrix(1, 4) = NAN;

        meth_mat.binary_matrix(1, 0) = 1;
        meth_mat.binary_matrix(1, 1) = 0;
        meth_mat.binary_matrix(1, 2) = 0;
        meth_mat.binary_matrix(1, 3) = -1;
        meth_mat.binary_matrix(1, 4) = -1;

        // Read 2: [0, 0, 1, 1, 1]
        meth_mat.raw_matrix(2, 0) = 0.10;
        meth_mat.raw_matrix(2, 1) = 0.05;
        meth_mat.raw_matrix(2, 2) = 0.90;
        meth_mat.raw_matrix(2, 3) = 0.95;
        meth_mat.raw_matrix(2, 4) = 0.85;

        meth_mat.binary_matrix(2, 0) = 0;
        meth_mat.binary_matrix(2, 1) = 0;
        meth_mat.binary_matrix(2, 2) = 1;
        meth_mat.binary_matrix(2, 3) = 1;
        meth_mat.binary_matrix(2, 4) = 1;

        // Read 3: [NaN, NaN, 1, 1, 1]
        meth_mat.raw_matrix(3, 0) = NAN;
        meth_mat.raw_matrix(3, 1) = NAN;
        meth_mat.raw_matrix(3, 2) = 0.92;
        meth_mat.raw_matrix(3, 3) = 0.88;
        meth_mat.raw_matrix(3, 4) = 0.95;

        meth_mat.binary_matrix(3, 0) = -1;
        meth_mat.binary_matrix(3, 1) = -1;
        meth_mat.binary_matrix(3, 2) = 1;
        meth_mat.binary_matrix(3, 3) = 1;
        meth_mat.binary_matrix(3, 4) = 1;

        // Set up ReadInfo for strand testing
        read_infos.resize(4);
        read_infos[0].read_id = 0;
        read_infos[0].strand = Strand::FORWARD;
        read_infos[1].read_id = 1;
        read_infos[1].strand = Strand::FORWARD;
        read_infos[2].read_id = 2;
        read_infos[2].strand = Strand::REVERSE;
        read_infos[3].read_id = 3;
        read_infos[3].strand = Strand::REVERSE;
    }

    MethylationMatrix meth_mat;
    std::vector<ReadInfo> read_infos;
};

// ============================================================================
// Basic NHD Tests
// ============================================================================

TEST_F(DistanceMatrixTest, NHD_BasicComputation) {
    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(meth_mat, DistanceMetricType::NHD,
                                      2,  // min_common_coverage
                                      NanDistanceStrategy::MAX_DIST);

    ASSERT_EQ(dist_mat.size(), 4);

    // Distance to self should be 0
    EXPECT_DOUBLE_EQ(dist_mat.get_distance(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(dist_mat.get_distance(1, 1), 0.0);

    // Read 0 vs Read 1: Common sites [0, 1, 2], differences at site 1
    // NHD = 1/3 ≈ 0.333
    double d01 = dist_mat.get_distance(0, 1);
    EXPECT_NEAR(d01, 1.0 / 3.0, 0.001);

    // Read 0 vs Read 2: Common sites [0, 1, 2, 3], all different at 0,1,2,3
    // NHD = 4/4 = 1.0
    double d02 = dist_mat.get_distance(0, 2);
    EXPECT_NEAR(d02, 1.0, 0.001);

    // Read 2 vs Read 3: Common sites [2, 3, 4], no differences
    // NHD = 0/3 = 0.0
    double d23 = dist_mat.get_distance(2, 3);
    EXPECT_NEAR(d23, 0.0, 0.001);

    // Matrix should be symmetric
    EXPECT_DOUBLE_EQ(dist_mat.get_distance(0, 1), dist_mat.get_distance(1, 0));
    EXPECT_DOUBLE_EQ(dist_mat.get_distance(0, 2), dist_mat.get_distance(2, 0));
}

TEST_F(DistanceMatrixTest, NHD_MinCommonCoverage) {
    DistanceMatrix dist_mat;

    // With min_common_coverage = 4, Read 0 vs Read 1 should be invalid
    // (only 3 common sites)
    dist_mat.compute_from_methylation(meth_mat, DistanceMetricType::NHD,
                                      4,  // min_common_coverage
                                      NanDistanceStrategy::MAX_DIST);

    // Read 0 vs Read 1: Only 3 common sites < 4, should be MAX_DIST (1.0)
    double d01 = dist_mat.get_distance(0, 1);
    EXPECT_NEAR(d01, 1.0, 0.001);

    // Read 0 vs Read 2: 4 common sites >= 4, should be computed
    double d02 = dist_mat.get_distance(0, 2);
    EXPECT_NEAR(d02, 1.0, 0.001);  // All different
}

// ============================================================================
// L1 (Manhattan) Distance Tests
// ============================================================================

TEST_F(DistanceMatrixTest, L1_BasicComputation) {
    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(meth_mat, DistanceMetricType::L1, 2, NanDistanceStrategy::MAX_DIST);

    ASSERT_EQ(dist_mat.size(), 4);

    // Distance to self should be 0
    EXPECT_NEAR(dist_mat.get_distance(0, 0), 0.0, 0.001);

    // Read 0 vs Read 1: Common sites [0, 1, 2]
    // L1 = mean(|0.95-0.85|, |0.90-0.15|, |0.10-0.10|) = mean(0.10, 0.75, 0.0) = 0.283
    double d01 = dist_mat.get_distance(0, 1);
    EXPECT_NEAR(d01, 0.283, 0.01);
}

// ============================================================================
// L2 (Euclidean) Distance Tests
// ============================================================================

TEST_F(DistanceMatrixTest, L2_BasicComputation) {
    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(meth_mat, DistanceMetricType::L2, 2, NanDistanceStrategy::MAX_DIST);

    ASSERT_EQ(dist_mat.size(), 4);

    // Distance to self should be 0
    EXPECT_NEAR(dist_mat.get_distance(0, 0), 0.0, 0.001);

    // Read 2 vs Read 3: Common sites [2, 3, 4]
    // L2 = sqrt(mean((0.90-0.92)^2, (0.95-0.88)^2, (0.85-0.95)^2))
    //    = sqrt(mean(0.0004, 0.0049, 0.01)) = sqrt(0.00513) ≈ 0.0716
    double d23 = dist_mat.get_distance(2, 3);
    EXPECT_NEAR(d23, 0.0716, 0.01);
}

// ============================================================================
// CORR (Correlation) Distance Tests
// ============================================================================

TEST_F(DistanceMatrixTest, CORR_BasicComputation) {
    // Create a special matrix for correlation testing
    MethylationMatrix corr_mat;
    corr_mat.region_id = 0;
    corr_mat.read_ids = {0, 1};
    corr_mat.cpg_ids = {0, 1, 2, 3};

    int n_reads = 2;
    int n_cpgs = 4;

    corr_mat.raw_matrix = Eigen::MatrixXd(n_reads, n_cpgs);
    corr_mat.binary_matrix = Eigen::MatrixXi(n_reads, n_cpgs);

    // Read 0: [0.1, 0.3, 0.7, 0.9] - increasing pattern
    corr_mat.raw_matrix(0, 0) = 0.1;
    corr_mat.raw_matrix(0, 1) = 0.3;
    corr_mat.raw_matrix(0, 2) = 0.7;
    corr_mat.raw_matrix(0, 3) = 0.9;

    // Read 1: [0.2, 0.4, 0.6, 0.8] - similar increasing pattern
    corr_mat.raw_matrix(1, 0) = 0.2;
    corr_mat.raw_matrix(1, 1) = 0.4;
    corr_mat.raw_matrix(1, 2) = 0.6;
    corr_mat.raw_matrix(1, 3) = 0.8;

    for (int i = 0; i < n_reads; ++i) {
        for (int j = 0; j < n_cpgs; ++j) {
            corr_mat.binary_matrix(i, j) = corr_mat.raw_matrix(i, j) >= 0.5 ? 1 : 0;
        }
    }

    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(corr_mat, DistanceMetricType::CORR,
                                      3,  // Need at least 3 for correlation
                                      NanDistanceStrategy::MAX_DIST);

    ASSERT_EQ(dist_mat.size(), 2);

    // Read 0 vs Read 1: Both have increasing patterns, should be positively correlated
    // CORR distance = (1 - r) / 2
    // If r is high (close to 1), distance should be close to 0
    double d01 = dist_mat.get_distance(0, 1);
    EXPECT_LT(d01, 0.2);  // Should be strongly positively correlated
}

// ============================================================================
// JACCARD Distance Tests
// ============================================================================

TEST_F(DistanceMatrixTest, JACCARD_BasicComputation) {
    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(meth_mat, DistanceMetricType::JACCARD, 2, NanDistanceStrategy::MAX_DIST);

    ASSERT_EQ(dist_mat.size(), 4);

    // Read 2 vs Read 3: Methylated sites in common
    // Read 2: sites 2,3,4 are methylated (binary=1)
    // Read 3: sites 2,3,4 are methylated, but 0,1 are missing
    // Common valid sites: 2,3,4, all methylated in both
    // Intersection = 3, Union = 3, Jaccard = 3/3 = 1
    // Jaccard distance = 1 - 1 = 0
    double d23 = dist_mat.get_distance(2, 3);
    EXPECT_NEAR(d23, 0.0, 0.001);

    // Read 0 vs Read 2: No overlap in methylated sites
    // Read 0 methylated at 0,1; Read 2 methylated at 2,3,4
    // Intersection = 0, Union = 4 (sites 0,1,2,3 have valid data in both)
    // Wait, for Jaccard we only count methylated sites
    // Jaccard distance = 1 - 0/4 = 1.0
    double d02 = dist_mat.get_distance(0, 2);
    EXPECT_NEAR(d02, 1.0, 0.001);
}

// ============================================================================
// Strand-Specific Distance Tests
// ============================================================================

TEST_F(DistanceMatrixTest, StrandSpecific_BasicComputation) {
    DistanceConfig config;
    config.metric = DistanceMetricType::NHD;
    config.min_common_coverage = 2;
    config.nan_strategy = NanDistanceStrategy::MAX_DIST;

    DistanceCalculator calc(config);

    auto [forward_mat, reverse_mat] = calc.compute_strand_specific(meth_mat, read_infos);

    // Forward matrix should have reads 0 and 1
    ASSERT_EQ(forward_mat.size(), 2);

    // Reverse matrix should have reads 2 and 3
    ASSERT_EQ(reverse_mat.size(), 2);

    // Forward reads 0 vs 1: NHD = 1/3 ≈ 0.333
    double d_fwd = forward_mat.get_distance(0, 1);
    EXPECT_NEAR(d_fwd, 1.0 / 3.0, 0.001);

    // Reverse reads 2 vs 3: NHD = 0/3 = 0
    double d_rev = reverse_mat.get_distance(0, 1);
    EXPECT_NEAR(d_rev, 0.0, 0.001);
}

// ============================================================================
// DistanceConfig Tests
// ============================================================================

TEST_F(DistanceMatrixTest, DistanceConfig_AllOptions) {
    DistanceConfig config;
    config.metric = DistanceMetricType::NHD;
    config.min_common_coverage = 3;
    config.nan_strategy = NanDistanceStrategy::MAX_DIST;
    config.max_distance_value = 1.0;
    config.use_binary_matrix = true;
    config.num_threads = 2;

    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(meth_mat, config);

    ASSERT_EQ(dist_mat.size(), 4);
    EXPECT_EQ(dist_mat.metric_type, DistanceMetricType::NHD);
    EXPECT_EQ(dist_mat.min_common_coverage, 3);
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_F(DistanceMatrixTest, EmptyMatrix) {
    MethylationMatrix empty_mat;
    empty_mat.region_id = 0;
    empty_mat.raw_matrix = Eigen::MatrixXd(0, 0);
    empty_mat.binary_matrix = Eigen::MatrixXi(0, 0);

    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(empty_mat, DistanceMetricType::NHD, 2, NanDistanceStrategy::MAX_DIST);

    EXPECT_TRUE(dist_mat.empty());
    EXPECT_EQ(dist_mat.size(), 0);
}

TEST_F(DistanceMatrixTest, SingleRead) {
    MethylationMatrix single_mat;
    single_mat.region_id = 0;
    single_mat.read_ids = {0};
    single_mat.cpg_ids = {0, 1, 2};
    single_mat.raw_matrix = Eigen::MatrixXd(1, 3);
    single_mat.binary_matrix = Eigen::MatrixXi(1, 3);
    single_mat.raw_matrix(0, 0) = 0.9;
    single_mat.raw_matrix(0, 1) = 0.1;
    single_mat.raw_matrix(0, 2) = 0.5;
    single_mat.binary_matrix(0, 0) = 1;
    single_mat.binary_matrix(0, 1) = 0;
    single_mat.binary_matrix(0, 2) = -1;

    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(single_mat, DistanceMetricType::NHD, 2, NanDistanceStrategy::MAX_DIST);

    EXPECT_EQ(dist_mat.size(), 1);
    EXPECT_DOUBLE_EQ(dist_mat.get_distance(0, 0), 0.0);
}

TEST_F(DistanceMatrixTest, AllMissing) {
    MethylationMatrix missing_mat;
    missing_mat.region_id = 0;
    missing_mat.read_ids = {0, 1};
    missing_mat.cpg_ids = {0, 1, 2};
    missing_mat.raw_matrix = Eigen::MatrixXd(2, 3);
    missing_mat.binary_matrix = Eigen::MatrixXi(2, 3);

    // All values are missing
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            missing_mat.raw_matrix(i, j) = NAN;
            missing_mat.binary_matrix(i, j) = -1;
        }
    }

    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(missing_mat, DistanceMetricType::NHD, 1, NanDistanceStrategy::MAX_DIST);

    // Should return MAX_DIST since no common valid sites
    double d01 = dist_mat.get_distance(0, 1);
    EXPECT_NEAR(d01, 1.0, 0.001);
}

// ============================================================================
// Statistics Tests
// ============================================================================

TEST_F(DistanceMatrixTest, Statistics) {
    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(meth_mat, DistanceMetricType::NHD, 2, NanDistanceStrategy::MAX_DIST);

    // 4 reads -> 6 pairs
    int total_pairs = (4 * 3) / 2;  // C(4,2) = 6
    EXPECT_EQ(dist_mat.num_valid_pairs + dist_mat.num_invalid_pairs, total_pairs);

    // Average common coverage should be positive
    EXPECT_GT(dist_mat.avg_common_coverage, 0.0);
}

// ============================================================================
// Metric Conversion Tests
// ============================================================================

TEST(DistanceCalculatorTest, MetricToString) {
    EXPECT_EQ(DistanceCalculator::metric_to_string(DistanceMetricType::NHD), "NHD");
    EXPECT_EQ(DistanceCalculator::metric_to_string(DistanceMetricType::L1), "L1");
    EXPECT_EQ(DistanceCalculator::metric_to_string(DistanceMetricType::L2), "L2");
    EXPECT_EQ(DistanceCalculator::metric_to_string(DistanceMetricType::CORR), "CORR");
    EXPECT_EQ(DistanceCalculator::metric_to_string(DistanceMetricType::JACCARD), "JACCARD");
}

TEST(DistanceCalculatorTest, StringToMetric) {
    EXPECT_EQ(DistanceCalculator::string_to_metric("NHD"), DistanceMetricType::NHD);
    EXPECT_EQ(DistanceCalculator::string_to_metric("nhd"), DistanceMetricType::NHD);
    EXPECT_EQ(DistanceCalculator::string_to_metric("HAMMING"), DistanceMetricType::NHD);
    EXPECT_EQ(DistanceCalculator::string_to_metric("L1"), DistanceMetricType::L1);
    EXPECT_EQ(DistanceCalculator::string_to_metric("MANHATTAN"), DistanceMetricType::L1);
    EXPECT_EQ(DistanceCalculator::string_to_metric("L2"), DistanceMetricType::L2);
    EXPECT_EQ(DistanceCalculator::string_to_metric("EUCLIDEAN"), DistanceMetricType::L2);
    EXPECT_EQ(DistanceCalculator::string_to_metric("CORR"), DistanceMetricType::CORR);
    EXPECT_EQ(DistanceCalculator::string_to_metric("PEARSON"), DistanceMetricType::CORR);
    EXPECT_EQ(DistanceCalculator::string_to_metric("JACCARD"), DistanceMetricType::JACCARD);
}

// ============================================================================
// CSV Output Test
// ============================================================================

TEST_F(DistanceMatrixTest, WriteCSV) {
    DistanceMatrix dist_mat;
    dist_mat.compute_from_methylation(meth_mat, DistanceMetricType::NHD, 2, NanDistanceStrategy::MAX_DIST);

    // Write to temp file
    std::string temp_path = "/tmp/test_distance_matrix.csv";
    dist_mat.write_csv(temp_path);

    // Read back and verify
    std::ifstream ifs(temp_path);
    ASSERT_TRUE(ifs.is_open());

    std::string header;
    std::getline(ifs, header);
    EXPECT_EQ(header.substr(0, 7), "read_id");  // Header should start with "read_id"

    ifs.close();
    std::remove(temp_path.c_str());
}
