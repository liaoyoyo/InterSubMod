#include "core/HierarchicalClustering.hpp"
#include <limits>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <queue>
#include <numeric>

namespace InterSubMod {

// ============================================================================
// HierarchicalClustering Implementation
// ============================================================================

HierarchicalClustering::HierarchicalClustering(LinkageMethod method) {
    config_.method = method;
}

HierarchicalClustering::HierarchicalClustering(const ClusteringConfig& config)
    : config_(config) {}

Tree HierarchicalClustering::build_tree(
    const DistanceMatrix& dist_matrix,
    const std::vector<std::string>& read_names
) {
    return build_tree(dist_matrix.dist_matrix, read_names);
}

Tree HierarchicalClustering::build_tree(
    const Eigen::MatrixXd& dist_matrix,
    const std::vector<std::string>& read_names
) {
    // 驗證輸入
    int n = dist_matrix.rows();
    if (n != dist_matrix.cols()) {
        throw std::runtime_error("Distance matrix must be square");
    }
    if (n != static_cast<int>(read_names.size())) {
        throw std::runtime_error("Number of names must match matrix dimension");
    }
    
    // 特殊情況處理
    if (n == 0) {
        return Tree();
    }
    if (n == 1) {
        Tree tree;
        tree.set_root(TreeNode::create_leaf(0, read_names[0]));
        return tree;
    }
    
    // 根據方法選擇演算法
    switch (config_.method) {
        case LinkageMethod::UPGMA:
            return build_upgma(dist_matrix, read_names);
        case LinkageMethod::WARD:
            return build_ward(dist_matrix, read_names);
        case LinkageMethod::SINGLE:
            return build_single(dist_matrix, read_names);
        case LinkageMethod::COMPLETE:
            return build_complete(dist_matrix, read_names);
        default:
            throw std::runtime_error("Unsupported linkage method");
    }
}

Tree HierarchicalClustering::build_upgma(
    const Eigen::MatrixXd& dist_matrix,
    const std::vector<std::string>& read_names
) {
    // UPGMA: 距離 = 兩群所有成員對的平均距離
    auto compute_distance = [](
        const Eigen::MatrixXd& D,
        const std::vector<int>& cluster_a,
        const std::vector<int>& cluster_b,
        int /*size_a*/, int /*size_b*/
    ) -> double {
        double sum = 0.0;
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                sum += D(i, j);
            }
        }
        return sum / (cluster_a.size() * cluster_b.size());
    };
    
    // UPGMA 高度 = 距離 / 2
    auto distance_to_height = [](double d) { return d / 2.0; };
    
    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_single(
    const Eigen::MatrixXd& dist_matrix,
    const std::vector<std::string>& read_names
) {
    // Single linkage: 距離 = 最近點對距離
    auto compute_distance = [](
        const Eigen::MatrixXd& D,
        const std::vector<int>& cluster_a,
        const std::vector<int>& cluster_b,
        int /*size_a*/, int /*size_b*/
    ) -> double {
        double min_dist = std::numeric_limits<double>::max();
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                min_dist = std::min(min_dist, D(i, j));
            }
        }
        return min_dist;
    };
    
    // Single linkage 高度 = 距離 / 2
    auto distance_to_height = [](double d) { return d / 2.0; };
    
    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_complete(
    const Eigen::MatrixXd& dist_matrix,
    const std::vector<std::string>& read_names
) {
    // Complete linkage: 距離 = 最遠點對距離
    auto compute_distance = [](
        const Eigen::MatrixXd& D,
        const std::vector<int>& cluster_a,
        const std::vector<int>& cluster_b,
        int /*size_a*/, int /*size_b*/
    ) -> double {
        double max_dist = 0.0;
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                max_dist = std::max(max_dist, D(i, j));
            }
        }
        return max_dist;
    };
    
    // Complete linkage 高度 = 距離 / 2
    auto distance_to_height = [](double d) { return d / 2.0; };
    
    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_ward(
    const Eigen::MatrixXd& dist_matrix,
    const std::vector<std::string>& read_names
) {
    // Ward's method: 最小化合併後的群內變異增量
    // 使用 Lance-Williams 公式的特殊形式
    // 注意：Ward 需要平方距離
    
    auto compute_distance = [](
        const Eigen::MatrixXd& D,
        const std::vector<int>& cluster_a,
        const std::vector<int>& cluster_b,
        int size_a, int size_b
    ) -> double {
        // 計算兩群中心之間的距離
        // Ward 距離增量 = (n_a * n_b / (n_a + n_b)) * d(centroid_a, centroid_b)^2
        // 這裡我們用群間平均距離來近似
        double sum = 0.0;
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                double d = D(i, j);
                sum += d * d;  // Ward 使用平方距離
            }
        }
        double avg_sq_dist = sum / (cluster_a.size() * cluster_b.size());
        
        // Ward 增量公式
        double n = static_cast<double>(size_a + size_b);
        double weight = (static_cast<double>(size_a) * static_cast<double>(size_b)) / n;
        return weight * avg_sq_dist;
    };
    
    // Ward 高度 = sqrt(距離增量) / 2，但通常直接用增量作為高度
    auto distance_to_height = [](double d) { return std::sqrt(d) / 2.0; };
    
    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_generic(
    const Eigen::MatrixXd& dist_matrix,
    const std::vector<std::string>& read_names,
    std::function<double(
        const Eigen::MatrixXd& D,
        const std::vector<int>& cluster_a,
        const std::vector<int>& cluster_b,
        int size_a, int size_b
    )> compute_cluster_distance,
    std::function<double(double distance)> distance_to_height
) {
    int n = dist_matrix.rows();
    
    // 初始化：每個 read 是一個獨立的 cluster
    std::vector<std::shared_ptr<TreeNode>> nodes;
    std::vector<std::vector<int>> cluster_members;  // 每個 cluster 包含的原始元素索引
    std::vector<int> cluster_sizes;
    std::vector<bool> active(n, true);  // 活躍的 cluster
    
    for (int i = 0; i < n; ++i) {
        auto leaf = TreeNode::create_leaf(i, read_names[i]);
        nodes.push_back(leaf);
        cluster_members.push_back({i});
        cluster_sizes.push_back(1);
    }
    
    // 合併記錄
    std::vector<MergeRecord> merge_records;
    
    int next_node_id = n;  // 內部節點 ID 從 n 開始
    int active_count = n;
    
    // 迭代合併，直到只剩一個 cluster
    while (active_count > 1) {
        // 找到距離最小的兩個活躍 cluster
        double min_dist = std::numeric_limits<double>::max();
        int min_i = -1, min_j = -1;
        
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            if (!active[i]) continue;
            for (int j = i + 1; j < static_cast<int>(nodes.size()); ++j) {
                if (!active[j]) continue;
                
                double dist = compute_cluster_distance(
                    dist_matrix,
                    cluster_members[i],
                    cluster_members[j],
                    cluster_sizes[i],
                    cluster_sizes[j]
                );
                
                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        
        if (min_i < 0 || min_j < 0) {
            // 所有剩餘的 cluster 之間距離為無窮大
            break;
        }
        
        // 計算合併高度
        double merge_height = distance_to_height(min_dist);
        
        // 確保高度單調遞增
        double max_child_height = std::max(nodes[min_i]->height, nodes[min_j]->height);
        if (merge_height < max_child_height + config_.min_branch_length) {
            merge_height = max_child_height + config_.min_branch_length;
        }
        
        // 建立新的內部節點
        auto new_node = TreeNode::create_internal(
            next_node_id,
            nodes[min_i],
            nodes[min_j],
            merge_height
        );
        
        // 記錄合併
        MergeRecord record;
        record.cluster_i = min_i;
        record.cluster_j = min_j;
        record.distance = min_dist;
        record.new_cluster_id = next_node_id;
        record.size = cluster_sizes[min_i] + cluster_sizes[min_j];
        merge_records.push_back(record);
        
        // 更新 cluster 資訊
        std::vector<int> new_members = cluster_members[min_i];
        new_members.insert(
            new_members.end(),
            cluster_members[min_j].begin(),
            cluster_members[min_j].end()
        );
        int new_size = cluster_sizes[min_i] + cluster_sizes[min_j];
        
        // 標記舊的 cluster 為非活躍
        active[min_i] = false;
        active[min_j] = false;
        
        // 加入新 cluster
        nodes.push_back(new_node);
        cluster_members.push_back(new_members);
        cluster_sizes.push_back(new_size);
        active.push_back(true);
        
        next_node_id++;
        active_count--;
    }
    
    // 找到最後一個活躍的 cluster 作為根節點
    std::shared_ptr<TreeNode> root = nullptr;
    for (int i = static_cast<int>(nodes.size()) - 1; i >= 0; --i) {
        if (active[i]) {
            root = nodes[i];
            break;
        }
    }
    
    // 建構 Tree 物件
    Tree tree;
    tree.set_root(root);
    tree.set_merge_records(merge_records);
    
    return tree;
}

double HierarchicalClustering::find_minimum_distance(
    const Eigen::MatrixXd& D,
    const std::vector<bool>& active,
    int& min_i,
    int& min_j
) {
    double min_dist = std::numeric_limits<double>::max();
    min_i = -1;
    min_j = -1;
    
    int n = D.rows();
    for (int i = 0; i < n; ++i) {
        if (!active[i]) continue;
        for (int j = i + 1; j < n; ++j) {
            if (!active[j]) continue;
            if (D(i, j) < min_dist) {
                min_dist = D(i, j);
                min_i = i;
                min_j = j;
            }
        }
    }
    
    return min_dist;
}

std::string HierarchicalClustering::method_to_string(LinkageMethod method) {
    switch (method) {
        case LinkageMethod::UPGMA: return "UPGMA";
        case LinkageMethod::WARD: return "WARD";
        case LinkageMethod::SINGLE: return "SINGLE";
        case LinkageMethod::COMPLETE: return "COMPLETE";
        default: return "UNKNOWN";
    }
}

LinkageMethod HierarchicalClustering::string_to_method(const std::string& str) {
    std::string upper = str;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    
    if (upper == "UPGMA" || upper == "AVERAGE") return LinkageMethod::UPGMA;
    if (upper == "WARD" || upper == "WARD.D" || upper == "WARD.D2") return LinkageMethod::WARD;
    if (upper == "SINGLE" || upper == "MIN") return LinkageMethod::SINGLE;
    if (upper == "COMPLETE" || upper == "MAX") return LinkageMethod::COMPLETE;
    
    // 預設為 UPGMA
    return LinkageMethod::UPGMA;
}

// ============================================================================
// TreeCutter Implementation
// ============================================================================

std::vector<int> TreeCutter::cut_by_distance(const Tree& tree, double distance_threshold) {
    if (tree.empty()) {
        return {};
    }
    
    auto root = tree.get_root();
    int n_leaves = root->num_leaves();
    std::vector<int> labels(n_leaves, 0);
    
    // BFS 遍歷樹，找到高度低於閾值的節點作為 cluster 根
    int current_label = 0;
    std::queue<std::shared_ptr<TreeNode>> queue;
    queue.push(root);
    
    while (!queue.empty()) {
        auto node = queue.front();
        queue.pop();
        
        // 高度轉換為距離（height 是到葉節點的距離，合併距離是 height * 2）
        double merge_distance = node->height * 2.0;
        
        if (node->is_leaf() || merge_distance <= distance_threshold) {
            // 這個節點及其所有葉節點屬於同一 cluster
            for (int leaf_idx : node->leaf_indices) {
                labels[leaf_idx] = current_label;
            }
            current_label++;
        } else {
            // 繼續向下探索
            if (node->left) queue.push(node->left);
            if (node->right) queue.push(node->right);
        }
    }
    
    return labels;
}

std::vector<int> TreeCutter::cut_by_num_clusters(const Tree& tree, int num_clusters) {
    if (tree.empty()) {
        return {};
    }
    
    auto root = tree.get_root();
    int n_leaves = root->num_leaves();
    
    if (num_clusters <= 0) {
        num_clusters = 1;
    }
    if (num_clusters >= n_leaves) {
        // 每個葉節點是一個 cluster
        std::vector<int> labels(n_leaves);
        std::iota(labels.begin(), labels.end(), 0);
        return labels;
    }
    
    // 收集所有合併高度，從高到低排序
    std::vector<std::pair<double, std::shared_ptr<TreeNode>>> height_nodes;
    auto internal_nodes = tree.get_internal_nodes();
    for (auto& node : internal_nodes) {
        height_nodes.emplace_back(node->height, node);
    }
    std::sort(height_nodes.begin(), height_nodes.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // 找到切割高度：使得產生恰好 num_clusters 個群
    // 從最高處開始，逐漸降低切割高度
    double cut_height = height_nodes.empty() ? 0.0 : height_nodes[0].first + 1.0;
    
    // 二分搜尋找到正確的切割高度
    for (const auto& [h, node] : height_nodes) {
        // 在高度 h 處切割會把這個節點分開
        // 計算當前切割會產生多少個 cluster
        auto test_labels = cut_by_distance(tree, h * 2.0 + 0.0001);  // 略高於 h
        int n_clusters = *std::max_element(test_labels.begin(), test_labels.end()) + 1;
        
        if (n_clusters >= num_clusters) {
            cut_height = h * 2.0;
            break;
        }
    }
    
    return cut_by_distance(tree, cut_height);
}

std::pair<int, std::vector<int>> TreeCutter::find_optimal_clusters(
    const Tree& tree,
    const Eigen::MatrixXd& dist_matrix,
    int min_k,
    int max_k
) {
    if (tree.empty()) {
        return {0, {}};
    }
    
    int n = tree.num_leaves();
    if (max_k > n) max_k = n;
    if (min_k < 2) min_k = 2;
    if (min_k > max_k) min_k = max_k;
    
    double best_score = -2.0;  // Silhouette 範圍是 [-1, 1]
    int best_k = min_k;
    std::vector<int> best_labels;
    
    for (int k = min_k; k <= max_k; ++k) {
        auto labels = cut_by_num_clusters(tree, k);
        
        // 計算 Silhouette Score
        double total_silhouette = 0.0;
        int valid_count = 0;
        
        for (int i = 0; i < n; ++i) {
            int cluster_i = labels[i];
            
            // 計算 a(i): 同群內的平均距離
            double a_i = 0.0;
            int count_a = 0;
            for (int j = 0; j < n; ++j) {
                if (j != i && labels[j] == cluster_i) {
                    a_i += dist_matrix(i, j);
                    count_a++;
                }
            }
            if (count_a > 0) a_i /= count_a;
            
            // 計算 b(i): 最近其他群的平均距離
            double b_i = std::numeric_limits<double>::max();
            for (int c = 0; c < k; ++c) {
                if (c == cluster_i) continue;
                
                double avg_dist = 0.0;
                int count_b = 0;
                for (int j = 0; j < n; ++j) {
                    if (labels[j] == c) {
                        avg_dist += dist_matrix(i, j);
                        count_b++;
                    }
                }
                if (count_b > 0) {
                    avg_dist /= count_b;
                    b_i = std::min(b_i, avg_dist);
                }
            }
            
            // 計算 silhouette
            if (count_a > 0 && b_i < std::numeric_limits<double>::max()) {
                double s_i = (b_i - a_i) / std::max(a_i, b_i);
                total_silhouette += s_i;
                valid_count++;
            }
        }
        
        double avg_silhouette = valid_count > 0 ? total_silhouette / valid_count : -1.0;
        
        if (avg_silhouette > best_score) {
            best_score = avg_silhouette;
            best_k = k;
            best_labels = labels;
        }
    }
    
    return {best_k, best_labels};
}

} // namespace InterSubMod

