/**
 * @file test_global_local.cpp
 * @brief Unit tests for GlobalTest and LocalTest
 *
 * Tests:
 * 1. Global association tests (allele, HP, sample)
 * 2. Gating logic
 * 3. Local one-vs-rest tests
 * 4. Best cluster identification
 */

#include <gtest/gtest.h>

#include <vector>

#include "core/GlobalTest.hpp"
#include "core/LocalTest.hpp"

using namespace InterSubMod;

// ============================================================================
// Test Fixtures
// ============================================================================

class GlobalTestTest : public ::testing::Test {
protected:
    void SetUp() override {
        GlobalTestConfig config;
        config.fisher_config.seed = 12345;
        config.fisher_config.min_mc_samples = 500;
        config.fisher_config.max_mc_samples = 1000;
        config.gating_p_threshold = 0.1;

        global_test_ = std::make_unique<GlobalTest>(config);
    }

    // Helper to create test data with strong allele association
    void create_strong_allele_association(std::vector<int>& cluster_labels,
                                           std::vector<FullLabel>& full_labels) {
        // 2 clusters: Cluster 0 has mostly ALT, Cluster 1 has mostly REF
        cluster_labels.clear();
        full_labels.clear();

        // Cluster 0: 15 ALT, 2 REF
        for (int i = 0; i < 15; ++i) {
            cluster_labels.push_back(0);
            FullLabel label;
            label.allele = AltSupport::ALT;
            label.hp_tag = "1";
            label.is_tumor = true;
            full_labels.push_back(label);
        }
        for (int i = 0; i < 2; ++i) {
            cluster_labels.push_back(0);
            FullLabel label;
            label.allele = AltSupport::REF;
            label.hp_tag = "1";
            label.is_tumor = true;
            full_labels.push_back(label);
        }

        // Cluster 1: 2 ALT, 15 REF
        for (int i = 0; i < 2; ++i) {
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
    }

    // Helper to create data with no association
    void create_no_association(std::vector<int>& cluster_labels,
                               std::vector<FullLabel>& full_labels) {
        cluster_labels.clear();
        full_labels.clear();

        // Cluster 0: 10 ALT, 10 REF
        for (int i = 0; i < 10; ++i) {
            cluster_labels.push_back(0);
            FullLabel label;
            label.allele = AltSupport::ALT;
            label.hp_tag = "1";
            label.is_tumor = true;
            full_labels.push_back(label);
        }
        for (int i = 0; i < 10; ++i) {
            cluster_labels.push_back(0);
            FullLabel label;
            label.allele = AltSupport::REF;
            label.hp_tag = "2";
            label.is_tumor = true;
            full_labels.push_back(label);
        }

        // Cluster 1: 10 ALT, 10 REF (same proportions)
        for (int i = 0; i < 10; ++i) {
            cluster_labels.push_back(1);
            FullLabel label;
            label.allele = AltSupport::ALT;
            label.hp_tag = "1";
            label.is_tumor = true;
            full_labels.push_back(label);
        }
        for (int i = 0; i < 10; ++i) {
            cluster_labels.push_back(1);
            FullLabel label;
            label.allele = AltSupport::REF;
            label.hp_tag = "2";
            label.is_tumor = true;
            full_labels.push_back(label);
        }
    }

    std::unique_ptr<GlobalTest> global_test_;
};

// ============================================================================
// GlobalTest Tests
// ============================================================================

TEST_F(GlobalTestTest, AlleleTest_StrongAssociation) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;
    create_strong_allele_association(cluster_labels, full_labels);

    GlobalTestResult result = global_test_->test_allele(cluster_labels, full_labels);

    EXPECT_TRUE(result.valid);
    EXPECT_LT(result.fisher_ffh.p_value, 0.01);  // Should be significant
    EXPECT_TRUE(result.passed_gate);
    EXPECT_GT(result.cramers_v, 0.3);  // Moderate to strong effect
}

TEST_F(GlobalTestTest, AlleleTest_NoAssociation) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;
    create_no_association(cluster_labels, full_labels);

    GlobalTestResult result = global_test_->test_allele(cluster_labels, full_labels);

    EXPECT_TRUE(result.valid);
    EXPECT_GT(result.fisher_ffh.p_value, 0.1);  // Should not be significant
    EXPECT_FALSE(result.passed_gate);  // Should fail gating
}

TEST_F(GlobalTestTest, HPTest_Works) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;
    create_strong_allele_association(cluster_labels, full_labels);

    // This test data also has HP1 in cluster 0 and HP2 in cluster 1
    GlobalTestResult result = global_test_->test_hp(cluster_labels, full_labels);

    EXPECT_TRUE(result.valid);
    // Should also find HP association since HP is confounded with allele in this data
    EXPECT_LT(result.fisher_ffh.p_value, 0.05);
}

TEST_F(GlobalTestTest, GatingLogic) {
    GlobalTestConfig config;
    config.gating_p_threshold = 0.05;  // Strict threshold
    config.fisher_config.seed = 42;
    config.fisher_config.min_mc_samples = 500;

    GlobalTest strict_test(config);

    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;
    create_no_association(cluster_labels, full_labels);

    GlobalTestResult result = strict_test.test_allele(cluster_labels, full_labels);

    // P-value should be high for no association
    EXPECT_GT(result.fisher_ffh.p_value, 0.05);
    // Should fail strict gating
    EXPECT_FALSE(result.passed_gate);
}

TEST_F(GlobalTestTest, EmptyInput) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    GlobalTestResult result = global_test_->test_allele(cluster_labels, full_labels);

    EXPECT_FALSE(result.valid);
}

TEST_F(GlobalTestTest, TestAll) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;
    create_strong_allele_association(cluster_labels, full_labels);

    auto [alt_result, hp_result, sample_result] = global_test_->test_all(cluster_labels, full_labels);

    EXPECT_TRUE(alt_result.valid);
    EXPECT_TRUE(hp_result.valid);
    EXPECT_TRUE(sample_result.valid);
}

// ============================================================================
// LocalTest Tests
// ============================================================================

class LocalTestTest : public ::testing::Test {
protected:
    void SetUp() override {
        LocalTestConfig config;
        config.fisher_config.seed = 12345;
        local_test_ = std::make_unique<LocalTest>(config);
    }

    std::unique_ptr<LocalTest> local_test_;
};

TEST_F(LocalTestTest, BasicLocalTest) {
    // Create simple data: Cluster 0 enriched for ALT
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Cluster 0: 18 ALT, 2 REF
    for (int i = 0; i < 18; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 2; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    // Cluster 1: 2 ALT, 18 REF
    for (int i = 0; i < 2; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 18; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    LocalTestResult result = local_test_->test_all(cluster_labels, full_labels);

    // Should have 2 clusters
    EXPECT_EQ(result.cluster_stats.size(), 2);

    // Best cluster should be identified
    EXPECT_NE(result.best_cluster_id, -1);
    EXPECT_LT(result.best_p_value, 0.001);  // Should be highly significant
}

TEST_F(LocalTestTest, ClusterCounts) {
    std::vector<int> cluster_labels = {0, 0, 0, 1, 1};
    std::vector<FullLabel> full_labels(5);

    for (int i = 0; i < 3; ++i) {
        full_labels[i].allele = AltSupport::ALT;
        full_labels[i].hp_tag = "1";
        full_labels[i].is_tumor = true;
    }
    for (int i = 3; i < 5; ++i) {
        full_labels[i].allele = AltSupport::REF;
        full_labels[i].hp_tag = "2";
        full_labels[i].is_tumor = false;
    }

    LocalTestResult result = local_test_->test_all(cluster_labels, full_labels);

    // Cluster 0 should have 3 reads
    bool found_cluster_0 = false;
    for (const auto& cs : result.cluster_stats) {
        if (cs.cluster_id == 0) {
            EXPECT_EQ(cs.size, 3);
            EXPECT_EQ(cs.count_alt, 3);
            EXPECT_EQ(cs.count_ref, 0);
            EXPECT_EQ(cs.count_hp1, 3);
            EXPECT_EQ(cs.count_tumor, 3);
            found_cluster_0 = true;
        }
    }
    EXPECT_TRUE(found_cluster_0);
}

TEST_F(LocalTestTest, DeltaProportion) {
    // Cluster 0: all ALT, Cluster 1: all REF
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    for (int i = 0; i < 10; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 10; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    LocalTestResult result = local_test_->test_all(cluster_labels, full_labels);

    // Cluster 0 should have delta_proportion_alt = 1.0 - 0.0 = 1.0
    for (const auto& cs : result.cluster_stats) {
        if (cs.cluster_id == 0) {
            EXPECT_NEAR(cs.delta_proportion_alt, 1.0, 0.01);
        }
    }
}

TEST_F(LocalTestTest, BestDimensionIdentified) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Strong HP association, weak allele association
    for (int i = 0; i < 20; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = (i < 10) ? AltSupport::ALT : AltSupport::REF;  // Mixed
        label.hp_tag = "1";  // All HP1
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 20; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = (i < 10) ? AltSupport::ALT : AltSupport::REF;  // Mixed
        label.hp_tag = "2";  // All HP2
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    LocalTestResult result = local_test_->test_all(cluster_labels, full_labels);

    // HP should be the best dimension since allele is balanced
    EXPECT_EQ(result.best_dimension, "hp");
}

// ============================================================================
// FullLabel Method Tests (HP Family / HP Fine)
// ============================================================================

TEST(FullLabelTest, HPFamilyLabel_Mapping) {
    FullLabel label;

    // HP1 → HP_FAMILY_1
    label.hp_tag = "1";
    EXPECT_EQ(label.get_hp_family_label(), TestLabel::HP_FAMILY_1);

    // HP1-1 → HP_FAMILY_1
    label.hp_tag = "1-1";
    EXPECT_EQ(label.get_hp_family_label(), TestLabel::HP_FAMILY_1);

    // HP2 → HP_FAMILY_2
    label.hp_tag = "2";
    EXPECT_EQ(label.get_hp_family_label(), TestLabel::HP_FAMILY_2);

    // HP2-1 → HP_FAMILY_2
    label.hp_tag = "2-1";
    EXPECT_EQ(label.get_hp_family_label(), TestLabel::HP_FAMILY_2);

    // HP3 → HP_OTHER
    label.hp_tag = "3";
    EXPECT_EQ(label.get_hp_family_label(), TestLabel::HP_OTHER);

    // HP0 → HP_OTHER
    label.hp_tag = "0";
    EXPECT_EQ(label.get_hp_family_label(), TestLabel::HP_OTHER);

    // Empty → HP_OTHER
    label.hp_tag = "";
    EXPECT_EQ(label.get_hp_family_label(), TestLabel::HP_OTHER);
}

TEST(FullLabelTest, HPFineLabel_Mapping) {
    FullLabel label;

    label.hp_tag = "1";
    EXPECT_EQ(label.get_hp_fine_label(), TestLabel::HP_FINE_1);

    label.hp_tag = "1-1";
    EXPECT_EQ(label.get_hp_fine_label(), TestLabel::HP_FINE_1_1);

    label.hp_tag = "2";
    EXPECT_EQ(label.get_hp_fine_label(), TestLabel::HP_FINE_2);

    label.hp_tag = "2-1";
    EXPECT_EQ(label.get_hp_fine_label(), TestLabel::HP_FINE_2_1);

    label.hp_tag = "3";
    EXPECT_EQ(label.get_hp_fine_label(), TestLabel::HP_OTHER);

    label.hp_tag = "0";
    EXPECT_EQ(label.get_hp_fine_label(), TestLabel::HP_OTHER);
}

// ============================================================================
// GlobalTest HP Family / HP Fine Tests
// ============================================================================

TEST_F(GlobalTestTest, HPFamilyTest_StrongAssociation) {
    // Create data where HP family (HP1+HP1-1 vs HP2+HP2-1) strongly associates with clusters,
    // but pure HP (HP1 vs HP2 only) has fewer reads and weaker power
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Cluster 0: mix of HP1 and HP1-1 (all HP-family-1)
    for (int i = 0; i < 8; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 7; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1-1";  // Somatic on HP1
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    // Cluster 1: mix of HP2 and HP2-1 (all HP-family-2)
    for (int i = 0; i < 8; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 7; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2-1";  // Somatic on HP2
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    GlobalTestResult result = global_test_->test_hp_family(cluster_labels, full_labels);

    EXPECT_TRUE(result.valid);
    EXPECT_LT(result.fisher_ffh.p_value, 0.01);
    EXPECT_TRUE(result.passed_gate);
    EXPECT_GT(result.cramers_v, 0.5);
}

TEST_F(GlobalTestTest, HPFineTest_FourGroups) {
    // Create data with 4 HP fine groups separated across clusters
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Cluster 0: HP1 and HP1-1
    for (int i = 0; i < 6; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 6; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1-1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    // Cluster 1: HP2 and HP2-1
    for (int i = 0; i < 6; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    for (int i = 0; i < 6; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2-1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    GlobalTestResult result = global_test_->test_hp_fine(cluster_labels, full_labels);

    EXPECT_TRUE(result.valid);
    EXPECT_LT(result.fisher_ffh.p_value, 0.05);
    EXPECT_EQ(result.n_valid_groups, 4);
}

TEST_F(GlobalTestTest, HPFineTest_SparseGroup) {
    // Create data where HP1-1 and HP2-1 have too few reads (< min_reads_per_group=3)
    // Only HP1 and HP2 remain as valid groups
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Cluster 0: 10 HP1, 1 HP1-1 (sparse)
    for (int i = 0; i < 10; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = "1-1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    // Cluster 1: 10 HP2, 1 HP2-1 (sparse)
    for (int i = 0; i < 10; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2";
        label.is_tumor = true;
        full_labels.push_back(label);
    }
    {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = AltSupport::REF;
        label.hp_tag = "2-1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    GlobalTestResult result = global_test_->test_hp_fine(cluster_labels, full_labels, 3);

    // Should still be valid with 2 groups (HP1 and HP2)
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.n_valid_groups, 2);
    EXPECT_LT(result.fisher_ffh.p_value, 0.05);
}

TEST_F(GlobalTestTest, HPFineTest_AllSparse_Invalid) {
    // All 4 groups have < min_reads → should be invalid
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Only 1 read per HP group
    for (const auto& hp : {"1", "1-1", "2", "2-1"}) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = AltSupport::ALT;
        label.hp_tag = hp;
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    GlobalTestResult result = global_test_->test_hp_fine(cluster_labels, full_labels, 3);

    // Less than 2 valid groups → invalid
    EXPECT_FALSE(result.valid);
}

TEST_F(GlobalTestTest, TestAllExpanded) {
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;
    create_strong_allele_association(cluster_labels, full_labels);

    auto expanded = global_test_->test_all_expanded(cluster_labels, full_labels);

    EXPECT_TRUE(expanded.alt.valid);
    EXPECT_TRUE(expanded.hp.valid);
    EXPECT_TRUE(expanded.hp_family.valid);
    EXPECT_TRUE(expanded.sample.valid);
    // hp_fine may or may not be valid depending on sparse groups
}

// ============================================================================
// LocalTest HP Family Tests
// ============================================================================

TEST_F(LocalTestTest, HPFamilyLocal_FineGrainedCounts) {
    // Verify fine-grained HP counting in initialize_cluster_stats
    std::vector<int> cluster_labels = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
    std::vector<FullLabel> full_labels(12);

    // Cluster 0: 2 HP1, 2 HP1-1, 1 HP3, 1 HP0
    full_labels[0].hp_tag = "1";   full_labels[0].allele = AltSupport::ALT; full_labels[0].is_tumor = true;
    full_labels[1].hp_tag = "1";   full_labels[1].allele = AltSupport::ALT; full_labels[1].is_tumor = true;
    full_labels[2].hp_tag = "1-1"; full_labels[2].allele = AltSupport::ALT; full_labels[2].is_tumor = true;
    full_labels[3].hp_tag = "1-1"; full_labels[3].allele = AltSupport::ALT; full_labels[3].is_tumor = true;
    full_labels[4].hp_tag = "3";   full_labels[4].allele = AltSupport::ALT; full_labels[4].is_tumor = true;
    full_labels[5].hp_tag = "0";   full_labels[5].allele = AltSupport::ALT; full_labels[5].is_tumor = true;

    // Cluster 1: 2 HP2, 2 HP2-1, 1 HP3, 1 HP0
    full_labels[6].hp_tag = "2";   full_labels[6].allele = AltSupport::REF; full_labels[6].is_tumor = true;
    full_labels[7].hp_tag = "2";   full_labels[7].allele = AltSupport::REF; full_labels[7].is_tumor = true;
    full_labels[8].hp_tag = "2-1"; full_labels[8].allele = AltSupport::REF; full_labels[8].is_tumor = true;
    full_labels[9].hp_tag = "2-1"; full_labels[9].allele = AltSupport::REF; full_labels[9].is_tumor = true;
    full_labels[10].hp_tag = "3";  full_labels[10].allele = AltSupport::REF; full_labels[10].is_tumor = true;
    full_labels[11].hp_tag = "0";  full_labels[11].allele = AltSupport::REF; full_labels[11].is_tumor = true;

    LocalTestResult result = local_test_->test_all(cluster_labels, full_labels);

    EXPECT_EQ(result.cluster_stats.size(), 2);

    for (const auto& cs : result.cluster_stats) {
        if (cs.cluster_id == 0) {
            EXPECT_EQ(cs.count_hp1, 2);
            EXPECT_EQ(cs.count_hp1_1, 2);
            EXPECT_EQ(cs.count_hp_family_1, 4);  // HP1 + HP1-1
            EXPECT_EQ(cs.count_hp_family_2, 0);
            EXPECT_EQ(cs.count_hp3, 1);
            EXPECT_EQ(cs.count_hp0, 1);
        } else if (cs.cluster_id == 1) {
            EXPECT_EQ(cs.count_hp2, 2);
            EXPECT_EQ(cs.count_hp2_1, 2);
            EXPECT_EQ(cs.count_hp_family_1, 0);
            EXPECT_EQ(cs.count_hp_family_2, 4);  // HP2 + HP2-1
            EXPECT_EQ(cs.count_hp3, 1);
            EXPECT_EQ(cs.count_hp0, 1);
        }
    }
}

TEST_F(LocalTestTest, HPFamilyLocal_BestDimension) {
    // All reads use "1-1" / "2-1" (somatic HP tags), so pure HP test has no valid reads,
    // but HP family test should still find strong association
    std::vector<int> cluster_labels;
    std::vector<FullLabel> full_labels;

    // Cluster 0: all HP1-1
    for (int i = 0; i < 15; ++i) {
        cluster_labels.push_back(0);
        FullLabel label;
        label.allele = (i % 2 == 0) ? AltSupport::ALT : AltSupport::REF;  // Mixed allele
        label.hp_tag = "1-1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    // Cluster 1: all HP2-1
    for (int i = 0; i < 15; ++i) {
        cluster_labels.push_back(1);
        FullLabel label;
        label.allele = (i % 2 == 0) ? AltSupport::ALT : AltSupport::REF;  // Mixed allele
        label.hp_tag = "2-1";
        label.is_tumor = true;
        full_labels.push_back(label);
    }

    LocalTestResult result = local_test_->test_all(cluster_labels, full_labels);

    // HP family should be the best dimension (not "hp" since pure has no valid reads)
    EXPECT_EQ(result.best_dimension, "hp_family");
    EXPECT_LT(result.best_p_value, 0.01);
}
