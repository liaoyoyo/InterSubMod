/**
 * @file test_significance_stats.cpp
 * @brief Unit tests for MathUtils and FisherExact
 *
 * Tests:
 * 1. log_factorial accuracy vs R lfactorial()
 * 2. log_binomial correctness
 * 3. odds_ratio with Haldane correction
 * 4. chi_square and cramers_v
 * 5. 2x2 Fisher's exact test vs R fisher.test()
 * 6. RxC Fisher-Freeman-Halton Monte Carlo
 */

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "core/FisherExact.hpp"
#include "core/MathUtils.hpp"
#include "core/Stats.hpp"

using namespace InterSubMod;

// ============================================================================
// MathUtils Tests
// ============================================================================

class MathUtilsTest : public ::testing::Test {
protected:
    // R reference values computed with:
    // lfactorial(0:20)
    const std::vector<double> r_log_factorial = {
        0.000000000,   // 0
        0.000000000,   // 1
        0.693147181,   // 2
        1.791759469,   // 3
        3.178053830,   // 4
        4.787491743,   // 5
        6.579251212,   // 6
        8.525161361,   // 7
        10.604602903,  // 8
        12.801827480,  // 9
        15.104412573,  // 10
        17.502307846,  // 11
        19.987214496,  // 12
        22.552163853,  // 13
        25.191221182,  // 14
        27.899271383,  // 15
        30.671860106,  // 16
        33.505073450,  // 17
        36.395445208,  // 18
        39.339884188,  // 19
        42.335616461   // 20
    };
};

TEST_F(MathUtilsTest, LogFactorial_SmallValues) {
    for (int n = 0; n <= 20; ++n) {
        double result = MathUtils::log_factorial(n);
        EXPECT_NEAR(result, r_log_factorial[n], 1e-6)
            << "Mismatch at n=" << n;
    }
}

TEST_F(MathUtilsTest, LogFactorial_LargeValues) {
    // R: lfactorial(100) = 363.7393756
    EXPECT_NEAR(MathUtils::log_factorial(100), 363.7393756, 1e-4);

    // R: lfactorial(1000) = 5912.128178
    EXPECT_NEAR(MathUtils::log_factorial(1000), 5912.128178, 1e-2);

    // R: lfactorial(10000) = 82108.92784
    EXPECT_NEAR(MathUtils::log_factorial(10000), 82108.92784, 1.0);
}

TEST_F(MathUtilsTest, LogFactorial_EdgeCases) {
    EXPECT_DOUBLE_EQ(MathUtils::log_factorial(0), 0.0);
    EXPECT_DOUBLE_EQ(MathUtils::log_factorial(1), 0.0);
    EXPECT_TRUE(std::isnan(MathUtils::log_factorial(-1)));
}

TEST_F(MathUtilsTest, LogBinomial_BasicCases) {
    // C(5, 2) = 10, log(10) = 2.302585
    EXPECT_NEAR(MathUtils::log_binomial(5, 2), std::log(10.0), 1e-10);

    // C(10, 5) = 252
    EXPECT_NEAR(MathUtils::log_binomial(10, 5), std::log(252.0), 1e-10);

    // C(n, 0) = 1
    EXPECT_DOUBLE_EQ(MathUtils::log_binomial(100, 0), 0.0);

    // C(n, n) = 1
    EXPECT_DOUBLE_EQ(MathUtils::log_binomial(100, 100), 0.0);
}

TEST_F(MathUtilsTest, LogBinomial_EdgeCases) {
    // k > n should return -inf
    EXPECT_TRUE(std::isinf(MathUtils::log_binomial(5, 10)));
    EXPECT_LT(MathUtils::log_binomial(5, 10), 0);

    // Negative values
    EXPECT_TRUE(std::isinf(MathUtils::log_binomial(-1, 0)));
}

TEST_F(MathUtilsTest, OddsRatio_BasicCases) {
    // Standard 2x2 table: a=10, b=5, c=3, d=12
    // OR = (10.5 * 12.5) / (5.5 * 3.5) = 131.25 / 19.25 ≈ 6.818
    double or_val = MathUtils::odds_ratio(10, 5, 3, 12);
    EXPECT_NEAR(or_val, 6.818, 0.01);
}

TEST_F(MathUtilsTest, OddsRatio_WithZeros) {
    // Zero cells should not cause division by zero
    double or_val = MathUtils::odds_ratio(10, 0, 0, 5);
    EXPECT_TRUE(std::isfinite(or_val));
    EXPECT_GT(or_val, 0);
}

TEST_F(MathUtilsTest, LogOddsRatio_Symmetry) {
    // log(OR(a,b,c,d)) = -log(OR(b,a,d,c))
    double log_or1 = MathUtils::log_odds_ratio(10, 5, 3, 12);
    double log_or2 = MathUtils::log_odds_ratio(5, 10, 12, 3);
    EXPECT_NEAR(log_or1, -log_or2, 1e-10);
}

TEST_F(MathUtilsTest, ChiSquare_IndependentData) {
    // Perfectly independent data
    // Expected counts equal observed counts -> chi2 = 0
    std::vector<int> table = {25, 25, 25, 25};  // 2x2, all equal
    double min_exp;
    double chi2 = MathUtils::chi_square(table, 2, 2, min_exp);
    EXPECT_DOUBLE_EQ(chi2, 0.0);
    EXPECT_DOUBLE_EQ(min_exp, 25.0);
}

TEST_F(MathUtilsTest, ChiSquare_DependentData) {
    // Strong association: a=20, b=0, c=0, d=20
    std::vector<int> table = {20, 0, 0, 20};
    double min_exp;
    double chi2 = MathUtils::chi_square(table, 2, 2, min_exp);
    // chi2 = 40 for this case
    EXPECT_NEAR(chi2, 40.0, 0.01);
}

TEST_F(MathUtilsTest, CramersV_Range) {
    // Cramér's V should be in [0, 1]
    std::vector<int> table = {10, 5, 3, 12};
    bool reliable;
    double v = MathUtils::cramers_v(table, 2, 2, reliable);
    EXPECT_GE(v, 0.0);
    EXPECT_LE(v, 1.0);
}

TEST_F(MathUtilsTest, CramersV_Reliability) {
    // Small counts should be marked unreliable
    std::vector<int> table = {1, 1, 1, 1};  // All cells have expected = 1 < 5
    bool reliable;
    MathUtils::cramers_v(table, 2, 2, reliable);
    EXPECT_FALSE(reliable);

    // Large counts should be reliable
    std::vector<int> large_table = {100, 100, 100, 100};
    MathUtils::cramers_v(large_table, 2, 2, reliable);
    EXPECT_TRUE(reliable);
}

TEST_F(MathUtilsTest, HypergeomPMF_BasicCases) {
    // Drawing k=2 successes from n=5 draws, with K=10 success states in N=20
    // P(X=2) should be a valid probability
    double p = MathUtils::hypergeom_pmf(2, 20, 10, 5);
    EXPECT_GT(p, 0.0);
    EXPECT_LT(p, 1.0);

    // Sum of all probabilities should equal 1
    double sum = 0.0;
    for (int k = 0; k <= 5; ++k) {
        sum += MathUtils::hypergeom_pmf(k, 20, 10, 5);
    }
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

// ============================================================================
// Fisher Exact Test
// ============================================================================

class FisherExactTest : public ::testing::Test {
protected:
    void SetUp() override {
        FisherConfig config;
        config.seed = 12345;  // Fixed seed for reproducibility
        fisher_ = std::make_unique<FisherExact>(config);
    }

    std::unique_ptr<FisherExact> fisher_;
};

TEST_F(FisherExactTest, Fisher2x2_LadyTastingTea) {
    // Classic Lady Tasting Tea example
    // R: fisher.test(matrix(c(3,1,1,3), nrow=2))$p.value = 0.4857143

    ContingencyTable2x2 table;
    table.a = 3;
    table.b = 1;
    table.c = 1;
    table.d = 3;

    FisherResult result = fisher_->test_2x2(table);

    // R reference: 0.4857143
    EXPECT_NEAR(result.p_value, 0.4857143, 1e-5);
}

TEST_F(FisherExactTest, Fisher2x2_SignificantAssociation) {
    // Strong association
    // R: fisher.test(matrix(c(10,0,0,10), nrow=2))$p.value ≈ 0.0000497

    ContingencyTable2x2 table;
    table.a = 10;
    table.b = 0;
    table.c = 0;
    table.d = 10;

    FisherResult result = fisher_->test_2x2(table);

    // Should be highly significant
    EXPECT_LT(result.p_value, 0.001);
}

TEST_F(FisherExactTest, Fisher2x2_NoAssociation) {
    // No association (proportional rows)
    // R: fisher.test(matrix(c(10,10,10,10), nrow=2))$p.value = 1

    ContingencyTable2x2 table;
    table.a = 10;
    table.b = 10;
    table.c = 10;
    table.d = 10;

    FisherResult result = fisher_->test_2x2(table);

    EXPECT_NEAR(result.p_value, 1.0, 1e-10);
}

TEST_F(FisherExactTest, Fisher2x2_OddsRatio) {
    ContingencyTable2x2 table;
    table.a = 10;
    table.b = 5;
    table.c = 3;
    table.d = 12;

    FisherResult result = fisher_->test_2x2(table);

    // OR should be > 1 (positive association)
    EXPECT_GT(result.odds_ratio, 1.0);
    EXPECT_GT(result.log_odds_ratio, 0.0);
}

TEST_F(FisherExactTest, Fisher2x2_EmptyTable) {
    ContingencyTable2x2 table;  // All zeros

    FisherResult result = fisher_->test_2x2(table);

    EXPECT_DOUBLE_EQ(result.p_value, 1.0);
}

TEST_F(FisherExactTest, FisherRxC_BasicTest) {
    // 3x2 table
    ContingencyTableRxC table;
    table.n_rows = 3;
    table.n_cols = 2;
    table.cells = {10, 5, 3, 12, 8, 7};

    FisherConfig config;
    config.seed = 42;
    config.min_mc_samples = 1000;
    config.max_mc_samples = 5000;

    FisherExact fisher(config);
    FisherResult result = fisher.test_rxc_mc(table);

    // P-value should be valid
    EXPECT_GE(result.p_value, 0.0);
    EXPECT_LE(result.p_value, 1.0);
    EXPECT_GT(result.n_samples, 0);
}

TEST_F(FisherExactTest, FisherRxC_StrongAssociation) {
    // Strong diagonal pattern
    ContingencyTableRxC table;
    table.n_rows = 2;
    table.n_cols = 2;
    table.cells = {20, 1, 1, 20};

    FisherConfig config;
    config.seed = 42;
    config.min_mc_samples = 1000;

    FisherExact fisher(config);
    FisherResult result = fisher.test_rxc_mc(table);

    // Should be significant
    EXPECT_LT(result.p_value, 0.01);
}

TEST_F(FisherExactTest, FisherRxC_EarlyStopping) {
    // Very clear case should trigger early stopping
    ContingencyTableRxC table;
    table.n_rows = 2;
    table.n_cols = 2;
    table.cells = {50, 0, 0, 50};

    FisherConfig config;
    config.seed = 42;
    config.min_mc_samples = 500;
    config.max_mc_samples = 10000;
    config.use_early_stopping = true;

    FisherExact fisher(config);
    FisherResult result = fisher.test_rxc_mc(table);

    // Should stop early because p-value CI is clearly < 0.05
    EXPECT_TRUE(result.early_stopped);
    EXPECT_LT(result.n_samples, config.max_mc_samples);
}

// ============================================================================
// Stats Structures Tests
// ============================================================================

TEST(StatsTest, ContingencyTable2x2_Totals) {
    ContingencyTable2x2 table;
    table.a = 10;
    table.b = 5;
    table.c = 3;
    table.d = 12;

    EXPECT_EQ(table.row1_total(), 15);
    EXPECT_EQ(table.row2_total(), 15);
    EXPECT_EQ(table.col1_total(), 13);
    EXPECT_EQ(table.col2_total(), 17);
    EXPECT_EQ(table.grand_total(), 30);
    EXPECT_TRUE(table.is_valid());
}

TEST(StatsTest, ContingencyTableRxC_Operations) {
    ContingencyTableRxC table;
    table.n_rows = 2;
    table.n_cols = 3;
    table.cells = {1, 2, 3, 4, 5, 6};

    EXPECT_EQ(table.get(0, 0), 1);
    EXPECT_EQ(table.get(1, 2), 6);
    EXPECT_EQ(table.row_total(0), 6);   // 1+2+3
    EXPECT_EQ(table.row_total(1), 15);  // 4+5+6
    EXPECT_EQ(table.col_total(0), 5);   // 1+4
    EXPECT_EQ(table.col_total(1), 7);   // 2+5
    EXPECT_EQ(table.col_total(2), 9);   // 3+6
    EXPECT_EQ(table.grand_total(), 21);
}

TEST(StatsTest, FullLabel_ComboCode) {
    FullLabel label;
    label.allele = AltSupport::ALT;
    label.hp_tag = "1";
    label.strand = Strand::FORWARD;
    label.is_tumor = true;

    uint16_t code = label.combo_code();

    // allele=0, hp=0, strand=0, sample=0 -> code = 0
    EXPECT_EQ(code, 0);

    label.allele = AltSupport::REF;
    label.hp_tag = "2";
    label.strand = Strand::REVERSE;
    label.is_tumor = false;

    code = label.combo_code();
    // allele=1, hp=1<<2=4, strand=1<<5=32, sample=1<<6=64
    // code = 1 + 4 + 32 + 64 = 101
    EXPECT_EQ(code, 101);
}

TEST(StatsTest, FullLabel_TestLabels) {
    FullLabel label;
    label.allele = AltSupport::ALT;
    label.hp_tag = "1";
    label.is_tumor = true;

    EXPECT_EQ(label.get_allele_label(), TestLabel::ALLELE_ALT);
    EXPECT_EQ(label.get_hp_label(), TestLabel::HP_1);
    EXPECT_EQ(label.get_sample_label(), TestLabel::SAMPLE_TUMOR);

    label.hp_tag = "1-1";
    EXPECT_EQ(label.get_hp_label(), TestLabel::HP_OTHER);

    label.hp_tag = "unphase";
    EXPECT_EQ(label.get_hp_label(), TestLabel::HP_OTHER);
}
