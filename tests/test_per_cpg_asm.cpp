/**
 * @file test_per_cpg_asm.cpp
 * @brief Unit tests for PerCpgAsm — per-CpG Fisher, NME, Epipolymorphism
 *
 * Tests cover:
 * 1. bh_fdr_correction() — basic, all NaN, monotonicity
 * 2. compute_nme() — uniform, single pattern, edge cases
 * 3. compute_epipolymorphism() — uniform, single pattern, window size
 * 4. compute_per_cpg_asm() — full pipeline, insufficient data, no ASM, clear ASM
 */

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <vector>

#include "core/PerCpgAsm.hpp"

using namespace InterSubMod;

// ============================================================================
// Helpers
// ============================================================================

static MethylationMatrix make_meth_mat(const Eigen::MatrixXd& raw, const Eigen::MatrixXi& binary) {
    MethylationMatrix mm;
    mm.region_id = 0;
    mm.raw_matrix = raw;
    mm.binary_matrix = binary;
    mm.read_ids.resize(raw.rows());
    mm.cpg_ids.resize(raw.cols());
    for (int i = 0; i < raw.rows(); ++i) mm.read_ids[i] = i;
    for (int j = 0; j < raw.cols(); ++j) mm.cpg_ids[j] = j;
    return mm;
}

static std::vector<ReadInfo> make_reads(int n, const std::vector<std::string>& hp_tags) {
    std::vector<ReadInfo> reads(n);
    for (int i = 0; i < n; ++i) {
        reads[i].read_id = i;
        reads[i].hp_tag = hp_tags[i];
        reads[i].is_tumor = true;
        reads[i].alt_support = AltSupport::UNKNOWN;
        reads[i].strand = Strand::FORWARD;
    }
    return reads;
}

// ============================================================================
// BH-FDR Correction Tests
// ============================================================================

TEST(PerCpgAsmTest, BhFdr_Basic) {
    std::vector<double> pvals = {0.01, 0.04, 0.03, 0.20};
    auto fdr = bh_fdr_correction(pvals);

    ASSERT_EQ(fdr.size(), 4u);
    // All should be valid (non-NaN)
    for (double q : fdr) {
        EXPECT_FALSE(std::isnan(q));
        EXPECT_GE(q, 0.0);
        EXPECT_LE(q, 1.0);
    }
    // Smallest p-value → smallest FDR
    EXPECT_LE(fdr[0], fdr[1]);
}

TEST(PerCpgAsmTest, BhFdr_AllNaN) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> pvals = {nan, nan, nan};
    auto fdr = bh_fdr_correction(pvals);
    for (double q : fdr) {
        EXPECT_TRUE(std::isnan(q));
    }
}

TEST(PerCpgAsmTest, BhFdr_SingleValue) {
    std::vector<double> pvals = {0.03};
    auto fdr = bh_fdr_correction(pvals);
    ASSERT_EQ(fdr.size(), 1u);
    EXPECT_NEAR(fdr[0], 0.03, 1e-10);  // m=1, rank=1, q = p*1/1 = p
}

TEST(PerCpgAsmTest, BhFdr_MixedNaN) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> pvals = {0.01, nan, 0.04};
    auto fdr = bh_fdr_correction(pvals);
    EXPECT_FALSE(std::isnan(fdr[0]));
    EXPECT_TRUE(std::isnan(fdr[1]));
    EXPECT_FALSE(std::isnan(fdr[2]));
}

// ============================================================================
// NME Tests
// ============================================================================

TEST(PerCpgAsmTest, NME_AllSamePattern) {
    // All reads have identical binary pattern → H=0 → NME=0
    Eigen::MatrixXi binary(6, 4);
    binary.setOnes();  // All 1s
    double nme = compute_nme(binary, 2);
    EXPECT_NEAR(nme, 0.0, 1e-10);
}

TEST(PerCpgAsmTest, NME_TwoEqualPatterns) {
    // 2 reads with opposite patterns → H = log2(2) = 1
    Eigen::MatrixXi binary(2, 3);
    binary.row(0) << 1, 1, 1;
    binary.row(1) << 0, 0, 0;
    double nme = compute_nme(binary, 2);
    // H = 1.0, H_max = log2(min(2, 2^3)) = log2(2) = 1.0
    EXPECT_NEAR(nme, 1.0, 1e-10);
}

TEST(PerCpgAsmTest, NME_MaxEntropy) {
    // 4 reads, 2 CpGs, all 4 patterns → maximum entropy
    Eigen::MatrixXi binary(4, 2);
    binary.row(0) << 0, 0;
    binary.row(1) << 0, 1;
    binary.row(2) << 1, 0;
    binary.row(3) << 1, 1;
    double nme = compute_nme(binary, 2);
    // H = log2(4) = 2.0, H_max = log2(min(4, 4)) = 2.0, NME = 1.0
    EXPECT_NEAR(nme, 1.0, 1e-10);
}

TEST(PerCpgAsmTest, NME_InsufficientCpgs) {
    Eigen::MatrixXi binary(5, 1);
    binary.setOnes();
    double nme = compute_nme(binary, 2);  // min_cpgs=2, but only 1
    EXPECT_TRUE(std::isnan(nme));
}

TEST(PerCpgAsmTest, NME_EmptyMatrix) {
    Eigen::MatrixXi binary(0, 4);
    double nme = compute_nme(binary, 2);
    EXPECT_TRUE(std::isnan(nme));
}

TEST(PerCpgAsmTest, NME_WithMissing) {
    // Some missing values (-1) should still work
    Eigen::MatrixXi binary(4, 3);
    binary.row(0) << 1, 1, -1;
    binary.row(1) << 1, 1, -1;
    binary.row(2) << 0, 0, -1;
    binary.row(3) << 0, 0, -1;
    double nme = compute_nme(binary, 2);
    EXPECT_FALSE(std::isnan(nme));
    EXPECT_GE(nme, 0.0);
    EXPECT_LE(nme, 1.0);
}

// ============================================================================
// Epipolymorphism Tests
// ============================================================================

TEST(PerCpgAsmTest, Epipoly_AllSamePattern) {
    // All reads same pattern → all p_i^2 = 1 → epipoly = 0
    Eigen::MatrixXi binary(10, 5);
    binary.setOnes();
    double ep = compute_epipolymorphism(binary, 4);
    EXPECT_NEAR(ep, 0.0, 1e-10);
}

TEST(PerCpgAsmTest, Epipoly_TwoEqualPatterns) {
    // 10 reads, 5 CpGs, half 0000x and half 1111x
    Eigen::MatrixXi binary(10, 5);
    for (int i = 0; i < 5; ++i) {
        binary.row(i) << 0, 0, 0, 0, 0;
    }
    for (int i = 5; i < 10; ++i) {
        binary.row(i) << 1, 1, 1, 1, 1;
    }
    double ep = compute_epipolymorphism(binary, 4);
    // In each 4-CpG window: 5 reads "0000" + 5 reads "1111"
    // p1=0.5, p2=0.5 → epipoly = 1 - 0.25 - 0.25 = 0.5
    EXPECT_NEAR(ep, 0.5, 1e-10);
}

TEST(PerCpgAsmTest, Epipoly_InsufficientCpgs) {
    Eigen::MatrixXi binary(10, 3);  // Only 3 CpGs, need 4
    binary.setOnes();
    double ep = compute_epipolymorphism(binary, 4);
    EXPECT_TRUE(std::isnan(ep));
}

TEST(PerCpgAsmTest, Epipoly_SingleRead) {
    Eigen::MatrixXi binary(1, 5);
    binary.setOnes();
    double ep = compute_epipolymorphism(binary, 4);
    // Only 1 read → valid_reads < 2 in each window → NaN
    EXPECT_TRUE(std::isnan(ep));
}

// ============================================================================
// Full Pipeline Tests: compute_per_cpg_asm
// ============================================================================

TEST(PerCpgAsmTest, FullPipeline_InsufficientReads) {
    // Only 3 reads per HP group (below min_reads_per_group=5)
    Eigen::MatrixXd raw(6, 4);
    raw.setConstant(0.8);
    Eigen::MatrixXi binary(6, 4);
    binary.setOnes();
    auto mm = make_meth_mat(raw, binary);

    std::vector<std::string> tags = {"1", "1", "1", "2", "2", "2"};
    auto reads = make_reads(6, tags);

    auto result = compute_per_cpg_asm(mm, reads, 5);
    EXPECT_FALSE(result.valid);
}

TEST(PerCpgAsmTest, FullPipeline_NoASM) {
    // 12 reads, 5 CpGs, all methylated → no ASM
    int n = 12, c = 5;
    Eigen::MatrixXd raw(n, c);
    raw.setConstant(0.9);
    Eigen::MatrixXi binary(n, c);
    binary.setOnes();
    auto mm = make_meth_mat(raw, binary);

    std::vector<std::string> tags = {"1", "1", "1", "1", "1", "1",
                                      "2", "2", "2", "2", "2", "2"};
    auto reads = make_reads(n, tags);

    auto result = compute_per_cpg_asm(mm, reads, 5);
    EXPECT_TRUE(result.valid);
    // No ASM: all reads same methylation → Fisher p=1 for all CpGs
    EXPECT_EQ(result.fisher_n_sig, 0);
    EXPECT_NEAR(result.fisher_frac_sig, 0.0, 1e-10);
    EXPECT_EQ(result.fisher_n_tested, c);
    // NME: all same pattern → 0
    EXPECT_NEAR(result.nme_hp1, 0.0, 1e-10);
    EXPECT_NEAR(result.nme_hp2, 0.0, 1e-10);
    EXPECT_NEAR(result.entropy_imbalance, 0.0, 1e-10);
    // Epipoly: all same pattern → 0
    EXPECT_NEAR(result.epipoly_hp1, 0.0, 1e-10);
    EXPECT_NEAR(result.epipoly_hp2, 0.0, 1e-10);
    EXPECT_NEAR(result.epipoly_delta, 0.0, 1e-10);
}

TEST(PerCpgAsmTest, FullPipeline_ClearASM) {
    // HP1 all methylated, HP2 all unmethylated → strong per-CpG ASM
    int n = 12, c = 5;
    Eigen::MatrixXd raw(n, c);
    Eigen::MatrixXi binary(n, c);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < c; ++j) {
            raw(i, j) = 0.9;    // HP1 methylated
            binary(i, j) = 1;
        }
    }
    for (int i = 6; i < 12; ++i) {
        for (int j = 0; j < c; ++j) {
            raw(i, j) = 0.1;    // HP2 unmethylated
            binary(i, j) = 0;
        }
    }
    auto mm = make_meth_mat(raw, binary);

    std::vector<std::string> tags = {"1", "1", "1", "1", "1", "1",
                                      "2", "2", "2", "2", "2", "2"};
    auto reads = make_reads(n, tags);

    auto result = compute_per_cpg_asm(mm, reads, 5);
    EXPECT_TRUE(result.valid);
    // Strong ASM: all CpGs should be significant
    EXPECT_EQ(result.fisher_n_sig, c);
    EXPECT_NEAR(result.fisher_frac_sig, 1.0, 1e-10);
    EXPECT_GT(result.fisher_max_neg_log_fdr, 0.0);
    // NME: each HP has a single pattern → NME=0 for both
    EXPECT_NEAR(result.nme_hp1, 0.0, 1e-10);
    EXPECT_NEAR(result.nme_hp2, 0.0, 1e-10);
    EXPECT_NEAR(result.entropy_imbalance, 0.0, 1e-10);
}

TEST(PerCpgAsmTest, FullPipeline_MixedHP_Tags) {
    // Test HP1-1 and HP2-1 tags (somatic) are correctly mapped to families
    int n = 12, c = 4;
    Eigen::MatrixXd raw(n, c);
    Eigen::MatrixXi binary(n, c);
    raw.setConstant(0.5);
    binary.setZero();

    auto mm = make_meth_mat(raw, binary);

    // Mix of germline and somatic HP tags
    std::vector<std::string> tags = {"1", "HP1", "1-1", "1", "1", "1",
                                      "2", "HP2", "2-1", "2", "2", "2"};
    auto reads = make_reads(n, tags);

    auto result = compute_per_cpg_asm(mm, reads, 5);
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.fisher_n_tested, c);
}

TEST(PerCpgAsmTest, FullPipeline_NaN_RawValues) {
    // Some NaN in raw_matrix → should skip those CpGs for those reads
    int n = 10, c = 4;
    Eigen::MatrixXd raw(n, c);
    Eigen::MatrixXi binary(n, c);

    raw.setConstant(0.9);
    binary.setOnes();

    // Set some values to NaN
    raw(0, 0) = std::numeric_limits<double>::quiet_NaN();
    raw(1, 1) = std::numeric_limits<double>::quiet_NaN();
    binary(0, 0) = -1;
    binary(1, 1) = -1;

    auto mm = make_meth_mat(raw, binary);

    std::vector<std::string> tags = {"1", "1", "1", "1", "1",
                                      "2", "2", "2", "2", "2"};
    auto reads = make_reads(n, tags);

    auto result = compute_per_cpg_asm(mm, reads, 5);
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.fisher_n_tested, c);  // All CpGs should still be testable
}

TEST(PerCpgAsmTest, FullPipeline_EntropyImbalance) {
    // HP1: diverse patterns, HP2: uniform → entropy imbalance > 0
    int c = 4;
    int n_hp1 = 8, n_hp2 = 8;
    int n = n_hp1 + n_hp2;

    Eigen::MatrixXd raw(n, c);
    Eigen::MatrixXi binary(n, c);

    // HP1: 4 different patterns (diverse)
    int patterns[4][4] = {{0,0,0,0}, {0,0,1,1}, {1,1,0,0}, {1,1,1,1}};
    for (int i = 0; i < n_hp1; ++i) {
        int p = i % 4;
        for (int j = 0; j < c; ++j) {
            binary(i, j) = patterns[p][j];
            raw(i, j) = patterns[p][j] == 1 ? 0.9 : 0.1;
        }
    }
    // HP2: all same pattern (uniform)
    for (int i = n_hp1; i < n; ++i) {
        for (int j = 0; j < c; ++j) {
            binary(i, j) = 1;
            raw(i, j) = 0.9;
        }
    }

    auto mm = make_meth_mat(raw, binary);
    std::vector<std::string> tags = {"1", "1", "1", "1", "1", "1", "1", "1",
                                      "2", "2", "2", "2", "2", "2", "2", "2"};
    auto reads = make_reads(n, tags);

    auto result = compute_per_cpg_asm(mm, reads, 5);
    EXPECT_TRUE(result.valid);

    // HP1 should have higher NME than HP2 (diverse vs uniform)
    EXPECT_GT(result.nme_hp1, result.nme_hp2);
    // Entropy imbalance should be positive
    EXPECT_GT(result.entropy_imbalance, 0.0);

    // HP1 epipolymorphism should be higher (more diverse patterns)
    EXPECT_GT(result.epipoly_hp1, result.epipoly_hp2);
    EXPECT_GT(result.epipoly_delta, 0.0);
}
