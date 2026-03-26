/**
 * @file test_structure_bootstrap.cpp
 * @brief Unit tests for StructureTest and Bootstrap
 */

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "core/Bootstrap.hpp"
#include "core/StructureTest.hpp"

using namespace InterSubMod;

// ============================================================================
// StructureTest Tests
// ============================================================================

class StructureTestTest : public ::testing::Test {
protected:
    void SetUp() override {
        StructureTestConfig config;
        config.seed = 12345;
        config.n_permutations = 99;  // Fewer for faster tests
        config.min_reads_for_permanova = 5;

        structure_test_ = std::make_unique<StructureTest>(config);
    }

    std::unique_ptr<StructureTest> structure_test_;
};

TEST_F(StructureTestTest, FilterReads_CompleteMatrix) {
    // Complete matrix with 6 reads - no filtering needed
    Eigen::MatrixXd dist(6, 6);
    dist.setZero();
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            dist(i, j) = 0.1 * (i + j);
            dist(j, i) = dist(i, j);
        }
    }

    std::vector<int> valid_indices;
    bool success = structure_test_->filter_reads_for_complete_matrix(dist, valid_indices);

    EXPECT_TRUE(success);
    EXPECT_EQ(valid_indices.size(), 6);
}

TEST_F(StructureTestTest, FilterReads_WithNans) {
    // Matrix with 7 reads, some NaN values
    Eigen::MatrixXd dist(7, 7);
    dist.setZero();

    // Fill with valid distances first
    for (int i = 0; i < 7; ++i) {
        for (int j = i + 1; j < 7; ++j) {
            dist(i, j) = 0.1 * std::abs(i - j);
            dist(j, i) = dist(i, j);
        }
    }

    // Add NaN for one problematic read (read 3)
    dist(3, 0) = NAN; dist(0, 3) = NAN;
    dist(3, 1) = NAN; dist(1, 3) = NAN;

    std::vector<int> valid_indices;
    bool success = structure_test_->filter_reads_for_complete_matrix(dist, valid_indices);

    // Should remove read 3, leaving 6 valid reads
    EXPECT_TRUE(success);
    EXPECT_GE(valid_indices.size(), 5);

    // Check no NaN in filtered submatrix
    for (size_t i = 0; i < valid_indices.size(); ++i) {
        for (size_t j = i + 1; j < valid_indices.size(); ++j) {
            EXPECT_FALSE(std::isnan(dist(valid_indices[i], valid_indices[j])));
        }
    }
}

TEST_F(StructureTestTest, FilterReads_TooManyNans) {
    // Matrix with too many NaNs - can't get enough valid reads
    Eigen::MatrixXd dist(5, 5);
    dist << 0, NAN, NAN, NAN, NAN,
            NAN, 0, NAN, NAN, NAN,
            NAN, NAN, 0, NAN, NAN,
            NAN, NAN, NAN, 0, NAN,
            NAN, NAN, NAN, NAN, 0;

    std::vector<int> valid_indices;
    bool success = structure_test_->filter_reads_for_complete_matrix(dist, valid_indices);

    EXPECT_FALSE(success);
}

TEST_F(StructureTestTest, Permanova_SignificantDifference) {
    // Two clearly separated groups
    Eigen::MatrixXd dist(10, 10);
    dist.setZero();

    // Group 0: reads 0-4 are close to each other
    // Group 1: reads 5-9 are close to each other
    // Between groups: large distance

    for (int i = 0; i < 10; ++i) {
        for (int j = i + 1; j < 10; ++j) {
            bool same_group = (i < 5 && j < 5) || (i >= 5 && j >= 5);
            double d = same_group ? 0.1 : 0.9;
            dist(i, j) = d;
            dist(j, i) = d;
        }
    }

    std::vector<int> labels = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    PermanovaResult result = structure_test_->run_permanova(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_GT(result.pseudo_f, 1.0);  // Should have some separation
    EXPECT_LT(result.p_value, 0.05);  // Should be significant
}

TEST_F(StructureTestTest, Permanova_NoDifference) {
    // Random labels - no real structure
    Eigen::MatrixXd dist(10, 10);
    dist.setZero();

    // Uniform distances
    for (int i = 0; i < 10; ++i) {
        for (int j = i + 1; j < 10; ++j) {
            dist(i, j) = 0.5;
            dist(j, i) = 0.5;
        }
    }

    std::vector<int> labels = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    PermanovaResult result = structure_test_->run_permanova(dist, labels);

    EXPECT_TRUE(result.valid);
    EXPECT_GT(result.p_value, 0.05);  // Should not be significant
}

TEST_F(StructureTestTest, Permanova_InvalidInput) {
    Eigen::MatrixXd dist(3, 3);  // Too few reads
    dist.setZero();
    std::vector<int> labels = {0, 0, 1};

    PermanovaResult result = structure_test_->run_permanova(dist, labels);

    EXPECT_FALSE(result.valid);
}

TEST_F(StructureTestTest, Dispersion_Homogeneous) {
    // Both groups have similar dispersion
    Eigen::MatrixXd dist(10, 10);
    dist.setZero();

    for (int i = 0; i < 10; ++i) {
        for (int j = i + 1; j < 10; ++j) {
            dist(i, j) = 0.3;
            dist(j, i) = 0.3;
        }
    }

    std::vector<int> labels = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    DispersionResult result = structure_test_->check_dispersion(dist, labels);

    EXPECT_FALSE(result.warning);  // No warning for homogeneous dispersion
    EXPECT_EQ(result.mean_distances_to_centroid.size(), 2);
}

TEST_F(StructureTestTest, Dispersion_Heterogeneous) {
    // Group 0: tight cluster, Group 1: loose cluster
    Eigen::MatrixXd dist(10, 10);
    dist.setZero();

    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            dist(i, j) = 0.05;  // Very tight
            dist(j, i) = 0.05;
        }
    }

    for (int i = 5; i < 10; ++i) {
        for (int j = i + 1; j < 10; ++j) {
            dist(i, j) = 0.8;  // Very loose
            dist(j, i) = 0.8;
        }
    }

    // Between groups
    for (int i = 0; i < 5; ++i) {
        for (int j = 5; j < 10; ++j) {
            dist(i, j) = 0.5;
            dist(j, i) = 0.5;
        }
    }

    std::vector<int> labels = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    DispersionResult result = structure_test_->check_dispersion(dist, labels);

    // Check that mean distances are different between groups
    EXPECT_EQ(result.mean_distances_to_centroid.size(), 2);
    // The first group (tight) should have smaller mean distance than second (loose)
    EXPECT_LT(result.mean_distances_to_centroid[0], result.mean_distances_to_centroid[1]);
}

// ============================================================================
// Bootstrap Tests
// ============================================================================

class BootstrapTest : public ::testing::Test {
protected:
    void SetUp() override {
        BootstrapConfig config;
        config.seed = 12345;
        config.n_iterations = 20;  // Fewer for faster tests
        config.use_early_stopping = false;

        bootstrap_ = std::make_unique<Bootstrap>(config);
    }

    std::unique_ptr<Bootstrap> bootstrap_;
};

TEST_F(BootstrapTest, ARI_PerfectAgreement) {
    std::vector<int> labels1 = {0, 0, 0, 1, 1, 1};
    std::vector<int> labels2 = {0, 0, 0, 1, 1, 1};

    double ari = Bootstrap::adjusted_rand_index(labels1, labels2);
    EXPECT_NEAR(ari, 1.0, 0.01);
}

TEST_F(BootstrapTest, ARI_RelabeledClusters) {
    // Same structure, different labels
    std::vector<int> labels1 = {0, 0, 0, 1, 1, 1};
    std::vector<int> labels2 = {1, 1, 1, 0, 0, 0};  // Swapped labels

    double ari = Bootstrap::adjusted_rand_index(labels1, labels2);
    EXPECT_NEAR(ari, 1.0, 0.01);  // ARI is label-agnostic
}

TEST_F(BootstrapTest, ARI_DifferentClustering) {
    std::vector<int> labels1 = {0, 0, 0, 1, 1, 1};
    std::vector<int> labels2 = {0, 1, 0, 1, 0, 1};  // Different structure

    double ari = Bootstrap::adjusted_rand_index(labels1, labels2);
    EXPECT_LT(ari, 0.5);  // Low similarity
}

TEST_F(BootstrapTest, ARI_AllSameCluster) {
    std::vector<int> labels1 = {0, 0, 0, 0, 0, 0};
    std::vector<int> labels2 = {0, 0, 0, 0, 0, 0};

    double ari = Bootstrap::adjusted_rand_index(labels1, labels2);
    EXPECT_NEAR(ari, 1.0, 0.01);
}

TEST_F(BootstrapTest, Bootstrap_StableCluster) {
    // Create methylation matrix with clear structure
    Eigen::MatrixXd meth(10, 20);

    // Group 0: high methylation
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 20; ++j) {
            meth(i, j) = 0.9;
        }
    }

    // Group 1: low methylation
    for (int i = 5; i < 10; ++i) {
        for (int j = 0; j < 20; ++j) {
            meth(i, j) = 0.1;
        }
    }

    std::vector<int> original_labels = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

    // Simple clustering function: split by mean methylation
    auto cluster_func = [](const Eigen::MatrixXd& dist) -> std::vector<int> {
        int n = dist.rows();
        std::vector<int> labels(n);
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            int count = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j && !std::isnan(dist(i, j))) {
                    sum += dist(i, j);
                    ++count;
                }
            }
            labels[i] = (count > 0 && sum / count < 0.5) ? 0 : 1;
        }
        return labels;
    };

    BootstrapResult result = bootstrap_->compute(meth, original_labels, cluster_func);

    // Check that bootstrap ran
    EXPECT_GT(result.n_iterations, 0);
    EXPECT_FALSE(result.ari_samples.empty());
}

TEST_F(BootstrapTest, Bootstrap_EarlyStopping) {
    BootstrapConfig config;
    config.seed = 42;
    config.n_iterations = 100;
    config.use_early_stopping = true;
    config.early_stop_consecutive = 5;
    config.early_stop_high_threshold = 0.9;

    Bootstrap bootstrap(config);

    // Create strongly structured data
    Eigen::MatrixXd meth(8, 10);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 10; ++j) meth(i, j) = 0.95;
    }
    for (int i = 4; i < 8; ++i) {
        for (int j = 0; j < 10; ++j) meth(i, j) = 0.05;
    }

    std::vector<int> labels = {0, 0, 0, 0, 1, 1, 1, 1};

    auto cluster_func = [](const Eigen::MatrixXd& dist) -> std::vector<int> {
        int n = dist.rows();
        std::vector<int> labels(n);
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                if (!std::isnan(dist(i, j))) sum += dist(i, j);
            }
            labels[i] = (sum / n < 0.4) ? 0 : 1;
        }
        return labels;
    };

    BootstrapResult result = bootstrap.compute(meth, labels, cluster_func);

    // May or may not early stop depending on random sampling
    EXPECT_GT(result.n_iterations, 0);
}

// ============================================================================
// PERMANOVA SS Accuracy Tests (Task 1.4 — catastrophic cancellation fix)
// ============================================================================

// Helper: build a symmetric distance matrix from a lambda
static Eigen::MatrixXd make_dist(int n, std::function<double(int, int)> fn) {
    Eigen::MatrixXd D(n, n);
    D.setZero();
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            D(i, j) = D(j, i) = fn(i, j);
    return D;
}

TEST(PermanovaSSTest, WellSeparated_SmallWithin_PseudoFIsLarge) {
    // Small within-group distances (0.1), large between-group (10.0).
    // Clear group structure → pseudo-F >> 1.
    const int ng = 10, n = 2 * ng;
    auto D = make_dist(n, [&](int i, int j) {
        return ((i < ng) == (j < ng)) ? 0.1 : 10.0;
    });
    std::vector<int> labels(n);
    for (int i = 0; i < n; ++i) labels[i] = (i < ng) ? 0 : 1;

    StructureTestConfig cfg;
    cfg.seed = 1;
    cfg.n_permutations = 99;
    cfg.min_reads_for_permanova = 5;
    StructureTest st(cfg);

    auto res = st.run_permanova(D, labels);
    EXPECT_TRUE(res.valid);
    EXPECT_GT(res.pseudo_f, 10.0) << "Well-separated groups should yield large pseudo-F";
    EXPECT_GE(res.pseudo_f, 0.0);
}

TEST(PermanovaSSTest, WellSeparated_ZeroWithin_PseudoFIsLarge) {
    // Perfectly separated: all within-group distances = 0, between = 1.
    // ss_within = 0 → old code returned 0 (division-by-zero guard).
    // Fixed code must return a large sentinel so permutation test is meaningful.
    const int ng = 10, n = 2 * ng;
    auto D = make_dist(n, [&](int i, int j) {
        return ((i < ng) == (j < ng)) ? 0.0 : 1.0;
    });
    std::vector<int> labels(n);
    for (int i = 0; i < n; ++i) labels[i] = (i < ng) ? 0 : 1;

    StructureTestConfig cfg;
    cfg.seed = 1;
    cfg.n_permutations = 99;
    cfg.min_reads_for_permanova = 5;
    StructureTest st(cfg);

    auto res = st.run_permanova(D, labels);
    EXPECT_TRUE(res.valid);
    // With fix: pseudo_f should be very large (perfect separation)
    EXPECT_GT(res.pseudo_f, 1e6) << "Zero within-group variance should yield very large pseudo-F";
}

TEST(PermanovaSSTest, NullModel_PseudoFIsNearOne) {
    // All pairwise distances identical (= 1.0) across 2 equal groups.
    // No group signal → pseudo-F ≈ 1.
    const int ng = 10, n = 2 * ng;
    auto D = make_dist(n, [](int, int) { return 1.0; });
    std::vector<int> labels(n);
    for (int i = 0; i < n; ++i) labels[i] = (i < ng) ? 0 : 1;

    StructureTestConfig cfg;
    cfg.seed = 1;
    cfg.n_permutations = 99;
    cfg.min_reads_for_permanova = 5;
    StructureTest st(cfg);

    auto res = st.run_permanova(D, labels);
    EXPECT_TRUE(res.valid);
    EXPECT_NEAR(res.pseudo_f, 1.0, 0.5) << "Null model pseudo-F should be near 1";
    EXPECT_GE(res.pseudo_f, 0.0);
}

TEST(PermanovaSSTest, CatastrophicCancellation_PseudoFNonNegative) {
    // n=40, k=2, ng=20. Large within distances (1e4) and between distances
    // chosen so that true ss_between ≈ 0 (between_dist = sqrt(0.95)*within_dist).
    //
    // Analysis: for balanced K=2 PERMANOVA,
    //   ss_between = 10*D_c² - 9.5*D_w²
    // Setting D_c = sqrt(0.95)*D_w gives ss_between_true = 0.
    //
    // In double precision, sqrt(0.95)*1e4 is represented as a slightly smaller
    // number, making D_c² < 0.95*D_w².  Then:
    //   ss_total = (within_raw + cross_raw)/n  < ss_within (rounding pushes cross down)
    // so the old formula produces ss_between < 0 → pseudo_f < 0.
    //
    // The fix: direct formula + clamp ss_between ≥ 0.
    const int ng = 20, n = 2 * ng;
    const double D_w = 1e4;
    const double D_c = std::sqrt(0.95) * D_w;  // nearest double → D_c² slightly < 0.95*D_w²

    auto D = make_dist(n, [&](int i, int j) {
        return ((i < ng) == (j < ng)) ? D_w : D_c;
    });
    std::vector<int> labels(n);
    for (int i = 0; i < n; ++i) labels[i] = (i < ng) ? 0 : 1;

    StructureTestConfig cfg;
    cfg.seed = 1;
    cfg.n_permutations = 99;
    cfg.min_reads_for_permanova = 5;
    StructureTest st(cfg);

    auto res = st.run_permanova(D, labels);
    EXPECT_TRUE(res.valid);
    EXPECT_GE(res.pseudo_f, 0.0)
        << "Direct ss_between formula must not produce negative pseudo-F";
}
