/**
 * @file test_phylo_labeler.cpp
 * @brief Unit tests for PhyloLabeler (phylo-v4.1 hierarchical cluster labeling, Phase B1 coarse).
 *
 * Covers: BERNOULLI distance recompute correctness, two-clear-group split (sensitivity),
 * pure-noise no-spurious-split (specificity ~0% FP), and small-n / degenerate edge cases.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

#include "core/PhyloLabeler.hpp"

using namespace InterSubMod;

namespace {

// Build a reads x CpG methylation matrix; group A fully methylated, group B fully unmethylated.
Eigen::MatrixXd two_group_meth(int per_group, int C) {
    int n = 2 * per_group;
    Eigen::MatrixXd M(n, C);
    for (int i = 0; i < n; ++i) {
        double v = (i < per_group) ? 0.95 : 0.05;  // confident, opposite states
        for (int c = 0; c < C; ++c) M(i, c) = v;
    }
    return M;
}

}  // namespace

TEST(PhyloLabeler, BernoulliIdenticalDeterministicIsZero) {
    // Identical reads at p=1.0 (deterministic methylated) -> delta=0 -> distance 0.
    PhyloLabeler lab;
    Eigen::MatrixXd M(2, 10);
    M.setConstant(1.0);
    Eigen::MatrixXd D = lab.bernoulli_dist(M);
    EXPECT_NEAR(D(0, 1), 0.0, 1e-9);
    EXPECT_NEAR(D(0, 0), 0.0, 1e-9);
}

TEST(PhyloLabeler, BernoulliIdenticalProbabilisticIsExpectedDisagreement) {
    // Identical reads at p=0.95 -> expected disagreement = 2*0.95*0.05 = 0.095 (NOT zero;
    // this is the BERNOULLI "expected disagreement" semantics, matching Python bernoulli_dist).
    PhyloLabeler lab;
    Eigen::MatrixXd M(2, 10);
    M.setConstant(0.95);
    Eigen::MatrixXd D = lab.bernoulli_dist(M);
    EXPECT_NEAR(D(0, 1), 0.095, 1e-3);
}

TEST(PhyloLabeler, BernoulliOppositeIsLarge) {
    PhyloLabeler lab;
    Eigen::MatrixXd M(2, 10);
    for (int c = 0; c < 10; ++c) {
        M(0, c) = 0.99;
        M(1, c) = 0.01;
    }
    Eigen::MatrixXd D = lab.bernoulli_dist(M);
    EXPECT_GT(D(0, 1), 0.9);  // near-maximal disagreement
}

TEST(PhyloLabeler, BernoulliInvalidWhenTooFewCommon) {
    PhyloConfig cfg;
    cfg.min_common = 3;
    PhyloLabeler lab(cfg);
    Eigen::MatrixXd M(2, 2);  // only 2 CpGs < min_common=3
    M.setConstant(0.9);
    Eigen::MatrixXd D = lab.bernoulli_dist(M);
    EXPECT_LT(D(0, 1), 0.0);  // -1 invalid
}

TEST(PhyloLabeler, TwoClearGroupsSplitIntoTwo) {
    PhyloLabeler lab;
    Eigen::MatrixXd M = two_group_meth(10, 30);  // 20 reads, 30 CpG, two opposite groups
    Eigen::MatrixXd D = lab.bernoulli_dist(M);
    PhyloResult r = lab.label(D, M);
    ASSERT_TRUE(r.valid);
    EXPECT_EQ(r.n_groups, 2);
    EXPECT_EQ(r.n_outlier, 0);
    // each read assigned a non-empty label
    EXPECT_EQ(static_cast<int>(r.labels.size()), 20);
}

TEST(PhyloLabeler, PureNoiseNoSpuriousSplit) {
    // Independent per-read Bernoulli with per-CpG rates: marginal structure, zero read structure.
    std::mt19937_64 rng(12345);
    int n = 40, C = 76;
    std::uniform_real_distribution<double> ur(0.1, 0.9);
    std::vector<double> rates(C);
    for (int c = 0; c < C; ++c) rates[c] = ur(rng);
    Eigen::MatrixXd M(n, C);
    for (int i = 0; i < n; ++i)
        for (int c = 0; c < C; ++c) {
            std::bernoulli_distribution b(rates[c]);
            M(i, c) = b(rng) ? 1.0 : 0.0;
        }
    PhyloLabeler lab;
    Eigen::MatrixXd D = lab.bernoulli_dist(M);
    PhyloResult r = lab.label(D, M);
    ASSERT_TRUE(r.valid);
    EXPECT_LE(r.n_groups, 1);  // null95 gate: noise must not split (specificity ~0% FP)
}

TEST(PhyloLabeler, TooFewReadsInvalid) {
    PhyloLabeler lab;
    Eigen::MatrixXd M = two_group_meth(2, 30);  // 4 reads < 2*min_sz=6
    Eigen::MatrixXd D = lab.bernoulli_dist(M);
    PhyloResult r = lab.label(D, M);
    EXPECT_FALSE(r.valid);
}
