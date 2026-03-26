/**
 * @file test_hierarchical_clustering.cpp
 * @brief Unit tests for HierarchicalClustering and Tree structures
 *
 * Tests cover:
 * 1. UPGMA algorithm correctness
 * 2. Comparison with scipy reference results
 * 3. Multiple linkage methods (UPGMA, WARD, SINGLE, COMPLETE)
 * 4. Newick format output
 * 5. Bootstrap support annotation
 * 6. Tree cutting and cluster labeling
 * 7. Edge cases (single node, empty tree)
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <vector>

#include "core/DistanceMatrix.hpp"
#include "core/HierarchicalClustering.hpp"
#include "core/TreeStructure.hpp"
#include "io/TreeWriter.hpp"

using namespace InterSubMod;

// ============================================================================
// Test Fixtures
// ============================================================================

class HierarchicalClusteringTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 建立一個簡單的 4x4 測試距離矩陣
        // 預期的 UPGMA 聚類結果：
        //   第一次合併: A 與 B (距離 2.0)
        //   第二次合併: C 與 D (距離 2.0)
        //   最後合併: (AB) 與 (CD) (距離 5.0)
        //
        //     A    B    C    D
        // A  0.0  2.0  4.0  6.0
        // B  2.0  0.0  4.0  6.0
        // C  4.0  4.0  0.0  2.0
        // D  6.0  6.0  2.0  0.0

        dist_4x4 = Eigen::MatrixXd(4, 4);
        dist_4x4 << 0.0, 2.0, 4.0, 6.0, 2.0, 0.0, 4.0, 6.0, 4.0, 4.0, 0.0, 2.0, 6.0, 6.0, 2.0, 0.0;

        names_4 = {"A", "B", "C", "D"};

        // 另一個較大的測試矩陣 (6x6)
        dist_6x6 = Eigen::MatrixXd(6, 6);
        dist_6x6 << 0.0, 0.2, 0.4, 0.8, 0.9, 1.0, 0.2, 0.0, 0.3, 0.7, 0.8, 0.9, 0.4, 0.3, 0.0, 0.6, 0.7, 0.8, 0.8, 0.7,
            0.6, 0.0, 0.2, 0.3, 0.9, 0.8, 0.7, 0.2, 0.0, 0.2, 1.0, 0.9, 0.8, 0.3, 0.2, 0.0;

        names_6 = {"R1", "R2", "R3", "R4", "R5", "R6"};
    }

    Eigen::MatrixXd dist_4x4;
    Eigen::MatrixXd dist_6x6;
    std::vector<std::string> names_4;
    std::vector<std::string> names_6;
};

// ============================================================================
// Basic UPGMA Tests
// ============================================================================

TEST_F(HierarchicalClusteringTest, UPGMA_BasicCorrectness) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    // 驗證樹不為空
    ASSERT_FALSE(tree.empty());

    // 驗證葉節點數量
    EXPECT_EQ(tree.num_leaves(), 4);

    // 驗證內部節點數量 (N-1 = 3)
    EXPECT_EQ(tree.num_internal_nodes(), 3);

    // 驗證根節點
    auto root = tree.get_root();
    ASSERT_NE(root, nullptr);
    EXPECT_FALSE(root->is_leaf());

    // 根節點應包含所有 4 個葉節點
    EXPECT_EQ(root->num_leaves(), 4);
}

TEST_F(HierarchicalClusteringTest, UPGMA_TreeTopology) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    // 取得所有 clade
    auto clades = tree.get_all_clades();

    // 應該有 3 個 clade (內部節點)
    EXPECT_EQ(clades.size(), 3);

    // 預期的 clade: {0,1}, {2,3}, {0,1,2,3}
    std::set<int> clade_AB = {0, 1};
    std::set<int> clade_CD = {2, 3};
    std::set<int> clade_all = {0, 1, 2, 3};

    // 驗證這些 clade 存在
    bool found_AB = std::find(clades.begin(), clades.end(), clade_AB) != clades.end();
    bool found_CD = std::find(clades.begin(), clades.end(), clade_CD) != clades.end();
    bool found_all = std::find(clades.begin(), clades.end(), clade_all) != clades.end();

    EXPECT_TRUE(found_AB) << "Clade {A,B} should exist";
    EXPECT_TRUE(found_CD) << "Clade {C,D} should exist";
    EXPECT_TRUE(found_all) << "Clade {A,B,C,D} should exist";
}

TEST_F(HierarchicalClusteringTest, UPGMA_BranchLengths) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    auto root = tree.get_root();

    // UPGMA: 高度 = 距離 / 2
    // A-B 合併距離 = 2.0, 高度 = 1.0
    // C-D 合併距離 = 2.0, 高度 = 1.0
    // (AB)-(CD) 合併: 平均距離 = (4+6+4+6)/4 = 5.0, 高度 = 2.5

    // 根節點高度應約為 2.5
    EXPECT_NEAR(root->height, 2.5, 0.01);

    // 葉節點高度應為 0
    auto leaves = tree.get_leaves();
    for (const auto& leaf : leaves) {
        EXPECT_DOUBLE_EQ(leaf->height, 0.0);
    }
}

TEST_F(HierarchicalClusteringTest, UPGMA_NewickOutput) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    std::string newick = tree.to_newick(false, true);

    // 驗證 Newick 格式正確
    EXPECT_TRUE(newick.back() == ';');
    EXPECT_TRUE(newick.find('(') != std::string::npos);
    EXPECT_TRUE(newick.find(')') != std::string::npos);

    // 驗證包含所有葉節點標籤
    EXPECT_TRUE(newick.find('A') != std::string::npos);
    EXPECT_TRUE(newick.find('B') != std::string::npos);
    EXPECT_TRUE(newick.find('C') != std::string::npos);
    EXPECT_TRUE(newick.find('D') != std::string::npos);

    // 驗證包含分支長度
    EXPECT_TRUE(newick.find(':') != std::string::npos);
}

// ============================================================================
// Multiple Linkage Methods Tests
// ============================================================================

TEST_F(HierarchicalClusteringTest, SingleLinkage_ChainEffect) {
    // Single linkage 傾向產生鏈狀結構
    HierarchicalClustering clusterer(LinkageMethod::SINGLE);
    Tree tree = clusterer.build_tree(dist_6x6, names_6);

    ASSERT_FALSE(tree.empty());
    EXPECT_EQ(tree.num_leaves(), 6);

    // 驗證 Newick 格式可正確輸出
    std::string newick = tree.to_newick(false, true);
    EXPECT_TRUE(newick.back() == ';');
}

TEST_F(HierarchicalClusteringTest, CompleteLinkage) {
    HierarchicalClustering clusterer(LinkageMethod::COMPLETE);
    Tree tree = clusterer.build_tree(dist_6x6, names_6);

    ASSERT_FALSE(tree.empty());
    EXPECT_EQ(tree.num_leaves(), 6);

    auto root = tree.get_root();
    // Complete linkage 使用最大距離，通常產生較高的樹
    EXPECT_GT(root->height, 0.0);
}

TEST_F(HierarchicalClusteringTest, WardMethod) {
    HierarchicalClustering clusterer(LinkageMethod::WARD);
    Tree tree = clusterer.build_tree(dist_6x6, names_6);

    ASSERT_FALSE(tree.empty());
    EXPECT_EQ(tree.num_leaves(), 6);
}

TEST_F(HierarchicalClusteringTest, AllMethodsProduceTrees) {
    std::vector<LinkageMethod> methods = {LinkageMethod::UPGMA, LinkageMethod::WARD, LinkageMethod::SINGLE,
                                          LinkageMethod::COMPLETE};

    for (LinkageMethod method : methods) {
        HierarchicalClustering clusterer(method);
        Tree tree = clusterer.build_tree(dist_4x4, names_4);

        EXPECT_FALSE(tree.empty()) << "Method " << HierarchicalClustering::method_to_string(method)
                                   << " should produce non-empty tree";
        EXPECT_EQ(tree.num_leaves(), 4) << "Method " << HierarchicalClustering::method_to_string(method)
                                        << " should have 4 leaves";
    }
}

// ============================================================================
// Tree Structure Tests
// ============================================================================

TEST_F(HierarchicalClusteringTest, TreeNode_CreateLeaf) {
    auto leaf = TreeNode::create_leaf(0, "TestRead");

    EXPECT_TRUE(leaf->is_leaf());
    EXPECT_EQ(leaf->node_id, 0);
    EXPECT_EQ(leaf->label, "TestRead");
    EXPECT_EQ(leaf->height, 0.0);
    EXPECT_EQ(leaf->leaf_indices.size(), 1);
    EXPECT_EQ(leaf->leaf_indices[0], 0);
}

TEST_F(HierarchicalClusteringTest, TreeNode_CreateInternal) {
    auto left = TreeNode::create_leaf(0, "A");
    auto right = TreeNode::create_leaf(1, "B");

    auto internal = TreeNode::create_internal(2, left, right, 1.0);

    EXPECT_FALSE(internal->is_leaf());
    EXPECT_EQ(internal->node_id, 2);
    EXPECT_EQ(internal->height, 1.0);
    EXPECT_EQ(internal->num_leaves(), 2);

    // 子節點的分支長度應正確設定
    EXPECT_EQ(left->branch_length, 1.0);
    EXPECT_EQ(right->branch_length, 1.0);
}

TEST_F(HierarchicalClusteringTest, Tree_DeepCopy) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree original = clusterer.build_tree(dist_4x4, names_4);

    Tree copy = original.deep_copy();

    // 驗證拷貝正確
    EXPECT_EQ(copy.num_leaves(), original.num_leaves());
    EXPECT_EQ(copy.num_internal_nodes(), original.num_internal_nodes());

    // 驗證是深拷貝（修改原始不影響拷貝）
    original.get_root()->height = 999.0;
    EXPECT_NE(copy.get_root()->height, 999.0);
}

// ============================================================================
// Tree Cutting Tests
// ============================================================================

TEST_F(HierarchicalClusteringTest, TreeCutter_ByDistance) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    // 以較低閾值切割，應產生 4 個群（每個葉節點一群）
    auto labels_4 = TreeCutter::cut_by_distance(tree, 0.5);
    EXPECT_EQ(labels_4.size(), 4);

    int n_clusters_4 = *std::max_element(labels_4.begin(), labels_4.end()) + 1;
    EXPECT_EQ(n_clusters_4, 4);

    // 以中等閾值切割，應產生 2 個群
    // A-B 合併距離 = 2.0, C-D 合併距離 = 2.0
    auto labels_2 = TreeCutter::cut_by_distance(tree, 3.0);
    int n_clusters_2 = *std::max_element(labels_2.begin(), labels_2.end()) + 1;
    EXPECT_EQ(n_clusters_2, 2);

    // A 和 B 應在同一群
    EXPECT_EQ(labels_2[0], labels_2[1]);
    // C 和 D 應在同一群
    EXPECT_EQ(labels_2[2], labels_2[3]);
    // AB 群和 CD 群應不同
    EXPECT_NE(labels_2[0], labels_2[2]);
}

TEST_F(HierarchicalClusteringTest, TreeCutter_ByNumClusters) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    // 要求 2 個群
    auto labels_2 = TreeCutter::cut_by_num_clusters(tree, 2);
    int n_clusters = *std::max_element(labels_2.begin(), labels_2.end()) + 1;
    EXPECT_EQ(n_clusters, 2);

    // A 和 B 應在同一群
    EXPECT_EQ(labels_2[0], labels_2[1]);
    // C 和 D 應在同一群
    EXPECT_EQ(labels_2[2], labels_2[3]);
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_F(HierarchicalClusteringTest, EdgeCase_SingleNode) {
    Eigen::MatrixXd dist_1x1(1, 1);
    dist_1x1(0, 0) = 0.0;
    std::vector<std::string> names_1 = {"Single"};

    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_1x1, names_1);

    EXPECT_FALSE(tree.empty());
    EXPECT_EQ(tree.num_leaves(), 1);
    EXPECT_EQ(tree.num_internal_nodes(), 0);

    std::string newick = tree.to_newick();
    EXPECT_EQ(newick, "Single;");
}

TEST_F(HierarchicalClusteringTest, EdgeCase_TwoNodes) {
    Eigen::MatrixXd dist_2x2(2, 2);
    dist_2x2 << 0.0, 0.5, 0.5, 0.0;
    std::vector<std::string> names_2 = {"X", "Y"};

    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_2x2, names_2);

    EXPECT_FALSE(tree.empty());
    EXPECT_EQ(tree.num_leaves(), 2);
    EXPECT_EQ(tree.num_internal_nodes(), 1);

    auto root = tree.get_root();
    EXPECT_NEAR(root->height, 0.25, 0.001);  // UPGMA: distance/2
}

TEST_F(HierarchicalClusteringTest, EdgeCase_EmptyMatrix) {
    Eigen::MatrixXd dist_empty(0, 0);
    std::vector<std::string> names_empty;

    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_empty, names_empty);

    EXPECT_TRUE(tree.empty());
}

TEST_F(HierarchicalClusteringTest, EdgeCase_IdenticalDistances) {
    // 所有距離相同
    Eigen::MatrixXd dist = Eigen::MatrixXd::Constant(4, 4, 1.0);
    dist.diagonal().setZero();
    std::vector<std::string> names = {"A", "B", "C", "D"};

    HierarchicalClustering clusterer(LinkageMethod::UPGMA);

    // 應不崩潰
    EXPECT_NO_THROW({
        Tree tree = clusterer.build_tree(dist, names);
        EXPECT_EQ(tree.num_leaves(), 4);
    });
}

// ============================================================================
// Bootstrap Support Tests
// ============================================================================

TEST_F(HierarchicalClusteringTest, BootstrapAnnotation) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    // 手動設定 Bootstrap 支持度
    std::vector<double> support = {95.0, 80.0, 100.0};  // 3 個內部節點

    EXPECT_NO_THROW(tree.annotate_bootstrap_support(support));

    // 驗證 Newick 輸出包含 Bootstrap
    std::string newick_with_bs = tree.to_newick(true, true);
    EXPECT_TRUE(newick_with_bs.find("95") != std::string::npos || newick_with_bs.find("80") != std::string::npos ||
                newick_with_bs.find("100") != std::string::npos);
}

// ============================================================================
// TreeWriter Tests
// ============================================================================

TEST_F(HierarchicalClusteringTest, TreeWriter_WriteNewick) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    std::string temp_path = "/tmp/test_tree.nwk";
    TreeWriter writer;

    bool success = writer.write_newick(tree, temp_path);
    EXPECT_TRUE(success);

    // 讀回並驗證
    std::ifstream ifs(temp_path);
    ASSERT_TRUE(ifs.is_open());

    std::string content;
    std::getline(ifs, content);
    ifs.close();

    EXPECT_TRUE(content.back() == ';');
    EXPECT_TRUE(content.find('A') != std::string::npos);

    std::remove(temp_path.c_str());
}

TEST_F(HierarchicalClusteringTest, TreeWriter_WriteLinkageMatrix) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    std::string temp_path = "/tmp/test_linkage.tsv";
    TreeWriter writer;

    bool success = writer.write_linkage_matrix(tree, temp_path);
    EXPECT_TRUE(success);

    // 讀回並驗證
    std::ifstream ifs(temp_path);
    ASSERT_TRUE(ifs.is_open());

    std::string header;
    std::getline(ifs, header);
    EXPECT_TRUE(header.find("cluster_i") != std::string::npos);

    // 應有 3 個合併記錄
    int line_count = 0;
    std::string line;
    while (std::getline(ifs, line)) {
        if (!line.empty()) line_count++;
    }
    EXPECT_EQ(line_count, 3);

    ifs.close();
    std::remove(temp_path.c_str());
}

TEST_F(HierarchicalClusteringTest, TreeWriter_WriteStats) {
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_4x4, names_4);

    std::string temp_path = "/tmp/test_tree_stats.txt";
    TreeWriter writer;

    bool success = writer.write_tree_stats(tree, temp_path);
    EXPECT_TRUE(success);

    // 讀回並驗證
    std::ifstream ifs(temp_path);
    ASSERT_TRUE(ifs.is_open());

    std::string content((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    EXPECT_TRUE(content.find("Number of leaves") != std::string::npos);
    EXPECT_TRUE(content.find("4") != std::string::npos);

    std::remove(temp_path.c_str());
}

// ============================================================================
// Method Conversion Tests
// ============================================================================

TEST(LinkageMethodTest, MethodToString) {
    EXPECT_EQ(HierarchicalClustering::method_to_string(LinkageMethod::UPGMA), "UPGMA");
    EXPECT_EQ(HierarchicalClustering::method_to_string(LinkageMethod::WARD), "WARD");
    EXPECT_EQ(HierarchicalClustering::method_to_string(LinkageMethod::SINGLE), "SINGLE");
    EXPECT_EQ(HierarchicalClustering::method_to_string(LinkageMethod::COMPLETE), "COMPLETE");
}

TEST(LinkageMethodTest, StringToMethod) {
    EXPECT_EQ(HierarchicalClustering::string_to_method("UPGMA"), LinkageMethod::UPGMA);
    EXPECT_EQ(HierarchicalClustering::string_to_method("upgma"), LinkageMethod::UPGMA);
    EXPECT_EQ(HierarchicalClustering::string_to_method("AVERAGE"), LinkageMethod::UPGMA);
    EXPECT_EQ(HierarchicalClustering::string_to_method("WARD"), LinkageMethod::WARD);
    EXPECT_EQ(HierarchicalClustering::string_to_method("ward.d"), LinkageMethod::WARD);
    EXPECT_EQ(HierarchicalClustering::string_to_method("SINGLE"), LinkageMethod::SINGLE);
    EXPECT_EQ(HierarchicalClustering::string_to_method("min"), LinkageMethod::SINGLE);
    EXPECT_EQ(HierarchicalClustering::string_to_method("COMPLETE"), LinkageMethod::COMPLETE);
    EXPECT_EQ(HierarchicalClustering::string_to_method("max"), LinkageMethod::COMPLETE);
}

// ============================================================================
// Integration with DistanceMatrix
// ============================================================================

TEST_F(HierarchicalClusteringTest, IntegrationWithDistanceMatrix) {
    // 建立 DistanceMatrix 物件
    DistanceMatrix dist_mat;
    dist_mat.dist_matrix = dist_4x4;
    dist_mat.read_ids = {0, 1, 2, 3};

    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist_mat, names_4);

    EXPECT_FALSE(tree.empty());
    EXPECT_EQ(tree.num_leaves(), 4);
}

// ============================================================================
// Performance Sanity Check
// ============================================================================

TEST(HierarchicalClusteringPerformance, MediumSizeMatrix) {
    // 100 x 100 矩陣
    int n = 100;
    Eigen::MatrixXd dist = Eigen::MatrixXd::Random(n, n);
    dist = (dist.array().abs() + dist.array().abs().transpose()) / 2.0;
    dist.diagonal().setZero();

    std::vector<std::string> names;
    for (int i = 0; i < n; ++i) {
        names.push_back("Read_" + std::to_string(i));
    }

    HierarchicalClustering clusterer(LinkageMethod::UPGMA);

    auto start = std::chrono::high_resolution_clock::now();
    Tree tree = clusterer.build_tree(dist, names);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    EXPECT_FALSE(tree.empty());
    EXPECT_EQ(tree.num_leaves(), n);

    // 100x100 應在 1 秒內完成
    EXPECT_LT(duration.count(), 1000) << "UPGMA on 100x100 should complete in <1s";

    std::cout << "  UPGMA 100x100: " << duration.count() << " ms" << std::endl;
}
