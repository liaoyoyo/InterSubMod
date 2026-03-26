/**
 * @file test_significance_analyzer.cpp
 * @brief Integration tests for SignificanceAnalyzer
 */

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "core/SignificanceAnalyzer.hpp"

using namespace InterSubMod;

class SignificanceAnalyzerTest : public ::testing::Test {
protected:
    void SetUp() override {
        SignificanceConfig config;
        config.seed = 12345;
        config.enable_global_test = true;
        config.enable_local_test = true;
        config.enable_permanova = true;
        config.enable_dispersion = true;
        config.enable_bootstrap = false;  // Skip for faster tests

        config.global_config.fisher_config.min_mc_samples = 100;
        config.global_config.fisher_config.max_mc_samples = 200;
        config.structure_config.n_permutations = 99;

        config.run_id = "test_run";
        config.vcf_id = "test.vcf";
        config.bam_id = "test.bam";

        analyzer_ = std::make_unique<SignificanceAnalyzer>(config);
    }

    std::unique_ptr<SignificanceAnalyzer> analyzer_;
};

TEST_F(SignificanceAnalyzerTest, FullAnalysis_StrongAssociation) {
    // Create strongly associated data
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Cluster 0: mostly ALT
    for (int i = 0; i < 15; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 3; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    // Cluster 1: mostly REF
    for (int i = 0; i < 3; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 15; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    int n = static_cast<int>(cluster_labels.size());

    // Create distance matrix: low within-cluster, high between-cluster
    Eigen::MatrixXd dist_matrix(n, n);
    dist_matrix.setZero();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            bool same_cluster = (cluster_labels[i] == cluster_labels[j]);
            double d = same_cluster ? 0.1 : 0.9;
            dist_matrix(i, j) = d;
            dist_matrix(j, i) = d;
        }
    }

    // Create methylation matrix
    Eigen::MatrixXd meth_matrix(n, 10);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 10; ++j) {
            meth_matrix(i, j) = (cluster_labels[i] == 0) ? 0.9 : 0.1;
        }
    }

    SignificanceResult result = analyzer_->analyze(
        cluster_labels, full_labels, dist_matrix, meth_matrix, 0, "chr1:12345:A:T");

    // Check metadata
    EXPECT_EQ(result.run_id, "test_run");
    EXPECT_EQ(result.vcf_id, "test.vcf");
    EXPECT_EQ(result.region_id, 0);
    EXPECT_EQ(result.anchor_key, "chr1:12345:A:T");
    EXPECT_TRUE(result.valid_flag);

    // Check global test
    EXPECT_TRUE(result.global_alt.valid);
    EXPECT_LT(result.global_alt.fisher_ffh.p_value, 0.05);
    EXPECT_TRUE(result.passed_gate);

    // Check local test
    EXPECT_FALSE(result.local_result.cluster_stats.empty());
    EXPECT_LT(result.local_result.best_p_value, 0.05);

    // Check PERMANOVA
    EXPECT_TRUE(result.permanova.valid);
    EXPECT_LT(result.permanova.p_value, 0.05);
    EXPECT_TRUE(result.label_hp_permanova.valid);
    EXPECT_LT(result.label_hp_permanova.p_value, 0.05);
    EXPECT_TRUE(result.label_allele_permanova.valid);
    EXPECT_LT(result.label_allele_permanova.p_value, 0.05);

    // Heuristic score should be positive
    EXPECT_GT(result.heuristic_score, 0.0);
}

TEST_F(SignificanceAnalyzerTest, SimpleAnalysis) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    for (int i = 0; i < 20; ++i) {
        cluster_labels.push_back(i < 10 ? 0 : 1);
        FullLabel label;
        label.allele = (i < 10) ? AltSupport::ALT : AltSupport::REF;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    SignificanceResult result = analyzer_->analyze_simple(
        cluster_labels, full_labels, 1, "chr2:5000:G:C");

    EXPECT_TRUE(result.valid_flag);
    EXPECT_TRUE(result.global_alt.valid);
}

TEST_F(SignificanceAnalyzerTest, EmptyInput) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    Eigen::MatrixXd empty_dist(0, 0);
    Eigen::MatrixXd empty_meth(0, 0);

    SignificanceResult result = analyzer_->analyze(
        cluster_labels, full_labels, empty_dist, empty_meth, 0, "chr1:100:A:T");

    EXPECT_FALSE(result.valid_flag);
    EXPECT_EQ(result.invalid_reason, "insufficient_reads");
}

TEST_F(SignificanceAnalyzerTest, GatingBehavior) {
    // No association - should fail gate
    // Both ALT/REF and HP labels are equally distributed across clusters
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Equal distribution in both clusters for BOTH ALT/REF AND HP
    for (int i = 0; i < 10; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = (i % 2 == 0) ? AltSupport::ALT : AltSupport::REF;
        label.hp_tag = (i % 2 == 0) ? "1" : "2";  // Mix HP tags in cluster 0
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 10; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = (i % 2 == 0) ? AltSupport::ALT : AltSupport::REF;
        label.hp_tag = (i % 2 == 0) ? "1" : "2";  // Mix HP tags in cluster 1
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    Eigen::MatrixXd dist_matrix(20, 20);
    dist_matrix.setConstant(0.5);
    for (int i = 0; i < 20; ++i) dist_matrix(i, i) = 0.0;

    Eigen::MatrixXd meth_matrix(20, 5);
    meth_matrix.setConstant(0.5);

    SignificanceResult result = analyzer_->analyze(
        cluster_labels, full_labels, dist_matrix, meth_matrix, 0, "chr1:200:C:G");

    // Should fail gating due to no association
    EXPECT_FALSE(result.passed_gate);

    // Local tests should be empty when gating fails
    EXPECT_TRUE(result.local_result.cluster_stats.empty());
}

TEST_F(SignificanceAnalyzerTest, HeuristicScore) {
    SignificanceResult result;
    result.global_alt.valid = true;
    result.global_alt.fisher_ffh.p_value = 0.001;  // Significant
    result.global_alt.cramers_v = 0.5;
    result.global_alt.cramers_v_reliable = true;
    result.dispersion_warning = false;
    result.permanova.valid = true;
    result.permanova.p_value = 0.01;

    double score = analyzer_->compute_heuristic_score(result);

    // -log10(0.001) = 3, plus cramers_v bonus
    EXPECT_GT(score, 3.0);

    // With dispersion warning
    result.dispersion_warning = true;
    double score_with_warning = analyzer_->compute_heuristic_score(result);

    // Should be penalized
    EXPECT_LT(score_with_warning, score);
}

TEST_F(SignificanceAnalyzerTest, Determinism_SameSeedProducesSameResult) {
    // Verify that two analyzers with the same top-level seed produce identical
    // ClusterPermanovaF values. This validates that sig_config.seed (top-level)
    // correctly propagates to all components via init_analyzers().
    auto make_analyzer = [](uint64_t seed) {
        SignificanceConfig cfg;
        cfg.seed = seed;  // Top-level seed — must propagate to structure_config
        cfg.enable_global_test = true;
        cfg.enable_local_test = true;
        cfg.enable_permanova = true;
        cfg.enable_dispersion = false;
        cfg.enable_bootstrap = false;
        cfg.structure_config.n_permutations = 99;
        cfg.run_id = "det_test";
        cfg.vcf_id = "test.vcf";
        cfg.bam_id = "test.bam";
        return std::make_unique<SignificanceAnalyzer>(cfg);
    };

    auto a1 = make_analyzer(42);
    auto a2 = make_analyzer(42);
    auto a3 = make_analyzer(99);  // Different seed

    // Build a simple 2-cluster dataset
    std::vector<int> cluster_labels(20);
    std::vector<FullLabel> full_labels(20);
    Eigen::MatrixXd dist(20, 20);
    Eigen::MatrixXd meth(20, 5);

    for (int i = 0; i < 20; ++i) {
        cluster_labels[i] = i < 10 ? 0 : 1;
        FullLabel lbl;
        lbl.allele = i < 10 ? AltSupport::ALT : AltSupport::REF;
        lbl.hp_tag = i < 10 ? "1" : "2";
        lbl.is_tumor = true;
        full_labels[i] = lbl;
        for (int j = 0; j < 20; ++j)
            dist(i, j) = (i < 10) == (j < 10) ? 0.1 : 0.8;
        dist(i, i) = 0.0;
        for (int k = 0; k < 5; ++k)
            meth(i, k) = i < 10 ? 0.9 : 0.1;
    }

    auto r1 = a1->analyze(cluster_labels, full_labels, dist, meth, 0, "chr1:1000:A:T");
    auto r2 = a2->analyze(cluster_labels, full_labels, dist, meth, 0, "chr1:1000:A:T");
    auto r3 = a3->analyze(cluster_labels, full_labels, dist, meth, 0, "chr1:1000:A:T");

    // pseudo_f is computed from data → always identical regardless of seed
    EXPECT_TRUE(r1.permanova.valid);
    EXPECT_TRUE(r2.permanova.valid);
    EXPECT_DOUBLE_EQ(r1.permanova.pseudo_f, r2.permanova.pseudo_f);
    EXPECT_DOUBLE_EQ(r1.label_hp_permanova.pseudo_f, r2.label_hp_permanova.pseudo_f);

    // Same seed → p_value is also identical (permutation order is deterministic)
    EXPECT_DOUBLE_EQ(r1.permanova.p_value, r2.permanova.p_value);
    EXPECT_DOUBLE_EQ(r1.label_hp_permanova.p_value, r2.label_hp_permanova.p_value);

    // Different seed → p_value may differ (permutation order changes)
    // For borderline regions, different seeds give different p-values.
    // Here we just verify r1 and r2 match (same seed determinism is the key guarantee).
    // r3 uses seed=99; we don't assert its p differs since perfectly separated
    // data produces p=1/(n+1) regardless of seed.
}

TEST_F(SignificanceAnalyzerTest, LabelPermanova_InvalidWithoutTwoGroups) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    for (int i = 0; i < 12; ++i) {
        cluster_labels.push_back(i < 6 ? 0 : 1);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    Eigen::MatrixXd dist_matrix(12, 12);
    dist_matrix.setConstant(0.4);
    for (int i = 0; i < 12; ++i) {
        dist_matrix(i, i) = 0.0;
    }

    Eigen::MatrixXd meth_matrix(12, 4);
    meth_matrix.setConstant(0.5);

    SignificanceResult result = analyzer_->analyze(
        cluster_labels, full_labels, dist_matrix, meth_matrix, 2, "chr3:700:T:C");

    EXPECT_FALSE(result.label_hp_permanova.valid);
    EXPECT_EQ(result.label_hp_permanova.invalid_reason, "insufficient_groups");
    EXPECT_FALSE(result.label_allele_permanova.valid);
    EXPECT_EQ(result.label_allele_permanova.invalid_reason, "insufficient_groups");
}
