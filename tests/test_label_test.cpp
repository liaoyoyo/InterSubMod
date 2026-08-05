/**
 * @file test_label_test.cpp
 * @brief Unit tests for LabelTest — previously 0% coverage (Task 2.2)
 *
 * Tests cover:
 * 1. compute_delta() — well-separated, null model, edge cases
 * 2. test_binary_groups() — significance detection, insufficient data
 * 3. test_all() — HP and allele dimensions, combined scenarios
 * 4. Label conversion edge cases (all-same-group, empty, single read)
 */

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <vector>

#include "core/LabelTest.hpp"
#include "core/Stats.hpp"
#include "core/Types.hpp"

using namespace InterSubMod;

// ============================================================================
// Helpers
// ============================================================================

/**
 * Build a symmetric distance matrix with:
 *   within-group pairs (i,j both in same half) = within_val
 *   between-group pairs = between_val
 *   diagonal = 0
 * n reads: first half = group 0, second half = group 1
 */
static Eigen::MatrixXd make_block_dist(int n, double within_val, double between_val) {
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(n, n);
    int half = n / 2;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            bool same_group = (i < half) == (j < half);
            double val = same_group ? within_val : between_val;
            m(i, j) = val;
            m(j, i) = val;
        }
    }
    return m;
}

/** group_labels: first half = 0, second half = 1 */
static std::vector<int> make_binary_labels(int n) {
    int half = n / 2;
    std::vector<int> labels(n);
    for (int i = 0; i < n; ++i) labels[i] = (i < half) ? 0 : 1;
    return labels;
}

/** Build FullLabel with hp_tag and allele */
static FullLabel make_label(const std::string& hp, AltSupport allele = AltSupport::UNKNOWN) {
    FullLabel lbl;
    lbl.hp_tag = hp;
    lbl.allele = allele;
    lbl.strand = Strand::FORWARD;
    lbl.is_tumor = true;
    return lbl;
}

// ============================================================================
// Test Fixture
// ============================================================================

class LabelTestFixture : public ::testing::Test {
protected:
    void SetUp() override {
        LabelTestConfig cfg;
        cfg.n_permutations = 199;   // Fast for tests
        cfg.min_reads_per_group = 3;
        cfg.min_total_reads = 6;
        cfg.alpha = 0.05;
        cfg.seed = 42;
        label_test = std::make_unique<LabelTest>(cfg);
    }

    std::unique_ptr<LabelTest> label_test;
};

// ============================================================================
// compute_delta() Tests
// ============================================================================

TEST_F(LabelTestFixture, ComputeDelta_WellSeparated) {
    // 6 reads: 3 group-0, 3 group-1
    // within = 0.1, between = 0.9 → delta = 0.8
    auto dist = make_block_dist(6, 0.1, 0.9);
    auto labels = make_binary_labels(6);

    LabelDeltaResult result = label_test->compute_delta(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_NEAR(result.within_mean, 0.1, 1e-6);
    EXPECT_NEAR(result.between_mean, 0.9, 1e-6);
    EXPECT_NEAR(result.delta, 0.8, 1e-6);
    EXPECT_EQ(result.n_group_a, 3);
    EXPECT_EQ(result.n_group_b, 3);
}

TEST_F(LabelTestFixture, ComputeDelta_NullModel) {
    // All distances equal → delta should be 0
    Eigen::MatrixXd dist = Eigen::MatrixXd::Constant(6, 6, 0.5);
    dist.diagonal().setZero();
    auto labels = make_binary_labels(6);

    LabelDeltaResult result = label_test->compute_delta(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_NEAR(result.delta, 0.0, 1e-6);
}

TEST_F(LabelTestFixture, ComputeDelta_AllSameGroup_Invalid) {
    // All reads have same group label → no between-group pairs
    auto dist = make_block_dist(6, 0.1, 0.9);
    std::vector<int> labels(6, 0);  // All group 0

    LabelDeltaResult result = label_test->compute_delta(dist, labels);

    EXPECT_FALSE(result.valid);
}

TEST_F(LabelTestFixture, ComputeDelta_InsufficientGroupSize_Invalid) {
    // Group A has only 1 read, below min_reads_per_group=3
    auto dist = make_block_dist(4, 0.1, 0.9);
    std::vector<int> labels = {0, 1, 1, 1};  // 1 in group 0, 3 in group 1

    LabelDeltaResult result = label_test->compute_delta(dist, labels);

    EXPECT_FALSE(result.valid);
}

TEST_F(LabelTestFixture, ComputeDelta_WithExcludedReads) {
    // Label -1 means excluded from test
    auto dist = make_block_dist(8, 0.1, 0.9);
    // reads 0,1,2 = group 0; reads 4,5,6 = group 1; reads 3,7 excluded
    std::vector<int> labels = {0, 0, 0, -1, 1, 1, 1, -1};

    LabelDeltaResult result = label_test->compute_delta(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_NEAR(result.delta, 0.8, 0.01);
}

// ============================================================================
// test_binary_groups() Tests
// ============================================================================

TEST_F(LabelTestFixture, BinaryGroups_WellSeparated_Significant) {
    // Clear separation should yield significant p-value
    auto dist = make_block_dist(8, 0.05, 0.95);
    auto labels = make_binary_labels(8);

    LabelDeltaResult result = label_test->test_binary_groups(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_GT(result.delta, 0.0);
    EXPECT_LT(result.p_value, 0.05);
    EXPECT_TRUE(result.significant);
    EXPECT_GT(result.n_permutations, 0);
}

TEST_F(LabelTestFixture, BinaryGroups_NullModel_NotSignificant) {
    // Uniform distances → permutation p-value should not be small
    Eigen::MatrixXd dist = Eigen::MatrixXd::Constant(8, 8, 0.5);
    dist.diagonal().setZero();
    auto labels = make_binary_labels(8);

    LabelDeltaResult result = label_test->test_binary_groups(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_NEAR(result.delta, 0.0, 1e-6);
    EXPECT_FALSE(result.significant);
}

TEST_F(LabelTestFixture, BinaryGroups_EmptyInput_Invalid) {
    Eigen::MatrixXd dist(0, 0);
    std::vector<int> labels;

    LabelDeltaResult result = label_test->test_binary_groups(dist, labels);

    EXPECT_FALSE(result.valid);
}

TEST_F(LabelTestFixture, BinaryGroups_SingleRead_Invalid) {
    Eigen::MatrixXd dist = Eigen::MatrixXd::Zero(1, 1);
    std::vector<int> labels = {0};

    LabelDeltaResult result = label_test->test_binary_groups(dist, labels);

    EXPECT_FALSE(result.valid);
}

// ============================================================================
// test_all() — HP dimension tests
// ============================================================================

TEST_F(LabelTestFixture, TestAll_HP_WellSeparated) {
    // 8 reads: 4 HP1, 4 HP2, well separated
    int n = 8;
    auto dist = make_block_dist(n, 0.05, 0.95);

    std::vector<FullLabel> labels;
    for (int i = 0; i < 4; ++i) labels.push_back(make_label("1"));   // HP1
    for (int i = 0; i < 4; ++i) labels.push_back(make_label("2"));   // HP2

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_TRUE(result.hp_multistage.valid);
    // Stage 1 merged test should detect separation
    EXPECT_GT(result.hp_multistage.merged.delta, 0.0);
    EXPECT_TRUE(result.hp_multistage.merged.significant);
}

TEST_F(LabelTestFixture, TestAll_HP_AllSameTag_Invalid) {
    // All reads have HP1 → no HP2 group → HP test invalid
    int n = 6;
    Eigen::MatrixXd dist = make_block_dist(n, 0.1, 0.9);

    std::vector<FullLabel> labels(n, make_label("1"));

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_FALSE(result.hp_multistage.merged.valid);
}

TEST_F(LabelTestFixture, TestAll_HP_UnphasedOnly_Invalid) {
    // All reads have unphased HP → HP test invalid
    int n = 6;
    Eigen::MatrixXd dist = make_block_dist(n, 0.1, 0.9);

    std::vector<FullLabel> labels(n, make_label("0"));  // HP0 = unphased

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_FALSE(result.hp_multistage.merged.valid);
}

// ============================================================================
// test_all() — Allele dimension tests
// ============================================================================

TEST_F(LabelTestFixture, TestAll_Allele_WellSeparated) {
    // 10 reads: 5 ALT, 5 REF, well separated (n=6 only has C(6,3)=20 permutations;
    // 10 reads gives C(10,5)=252 so alpha=0.05 is reachable)
    int n = 10;
    auto dist = make_block_dist(n, 0.05, 0.95);

    std::vector<FullLabel> labels;
    for (int i = 0; i < 5; ++i) labels.push_back(make_label("0", AltSupport::ALT));
    for (int i = 0; i < 5; ++i) labels.push_back(make_label("0", AltSupport::REF));

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_TRUE(result.allele_result.valid);
    EXPECT_GT(result.allele_result.delta, 0.0);
    EXPECT_TRUE(result.allele_result.significant);
}

TEST_F(LabelTestFixture, TestAll_Allele_AllUnknown_Invalid) {
    // All alleles UNKNOWN → allele test invalid
    int n = 6;
    Eigen::MatrixXd dist = make_block_dist(n, 0.1, 0.9);

    std::vector<FullLabel> labels(n, make_label("0", AltSupport::UNKNOWN));

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_FALSE(result.allele_result.valid);
}

// ============================================================================
// test_all() — Mixed HP + Allele
// ============================================================================

TEST_F(LabelTestFixture, TestAll_BothDimensions) {
    // 8 reads: both HP and allele present
    int n = 8;
    auto dist = make_block_dist(n, 0.05, 0.95);

    std::vector<FullLabel> labels;
    for (int i = 0; i < 4; ++i) labels.push_back(make_label("1", AltSupport::ALT));
    for (int i = 0; i < 4; ++i) labels.push_back(make_label("2", AltSupport::REF));

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_TRUE(result.valid);
    // At least one significant dimension
    EXPECT_TRUE(result.any_significant);
    EXPECT_FALSE(result.dominant_dimension.empty());
    EXPECT_NE(result.dominant_dimension, "none");
}

TEST_F(LabelTestFixture, TestAll_EmptyInput_Invalid) {
    Eigen::MatrixXd dist(0, 0);
    std::vector<FullLabel> labels;

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_FALSE(result.valid);
}

// ============================================================================
// HP sub-group label mapping edge cases
// ============================================================================

TEST_F(LabelTestFixture, TestAll_HP_FamilyMerge) {
    // HP1 + HP1-1 should merge into same family group
    // HP2 + HP2-1 should merge into same family group
    int n = 8;
    auto dist = make_block_dist(n, 0.05, 0.95);

    std::vector<FullLabel> labels;
    labels.push_back(make_label("1"));    // HP1 → family 0
    labels.push_back(make_label("1"));
    labels.push_back(make_label("1-1")); // HP1-1 → family 0 (same as HP1)
    labels.push_back(make_label("1-1"));
    labels.push_back(make_label("2"));    // HP2 → family 1
    labels.push_back(make_label("2"));
    labels.push_back(make_label("2-1")); // HP2-1 → family 1 (same as HP2)
    labels.push_back(make_label("2-1"));

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_TRUE(result.valid);
    // HP1-family = 4 (HP1 + HP1-1), HP2-family = 4
    EXPECT_EQ(result.hp_multistage.n_hp1_family, 4);
    EXPECT_EQ(result.hp_multistage.n_hp2_family, 4);
    EXPECT_TRUE(result.hp_multistage.merged.valid);
}

TEST_F(LabelTestFixture, TestAll_HP_InsufficientFamilySize_Invalid) {
    // HP1-family has only 1 read (below min_reads_per_group=3)
    int n = 6;
    Eigen::MatrixXd dist = make_block_dist(n, 0.1, 0.9);

    std::vector<FullLabel> labels;
    labels.push_back(make_label("1"));   // Only 1 in HP1-family
    for (int i = 0; i < 5; ++i) labels.push_back(make_label("2"));

    LabelTestResult result = label_test->test_all(dist, labels);

    EXPECT_FALSE(result.hp_multistage.merged.valid);
}
