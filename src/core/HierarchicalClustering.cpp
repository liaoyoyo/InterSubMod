#include "core/HierarchicalClustering.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>

namespace InterSubMod {

// ============================================================================
// HierarchicalClustering Implementation
// ============================================================================

HierarchicalClustering::HierarchicalClustering(LinkageMethod method) {
    config_.method = method;
}

HierarchicalClustering::HierarchicalClustering(const ClusteringConfig& config) : config_(config) {
}

Tree HierarchicalClustering::build_tree(const DistanceMatrix& dist_matrix, const std::vector<std::string>& read_names) {
    return build_tree(dist_matrix.dist_matrix, read_names);
}

Tree HierarchicalClustering::build_tree(const Eigen::MatrixXd& dist_matrix,
                                        const std::vector<std::string>& read_names) {
    // Validate input
    int n = dist_matrix.rows();
    if (n != dist_matrix.cols()) {
        throw std::runtime_error("Distance matrix must be square");
    }
    if (n != static_cast<int>(read_names.size())) {
        throw std::runtime_error("Number of names must match matrix dimension");
    }

    // Handle special cases
    if (n == 0) {
        return Tree();
    }
    if (n == 1) {
        Tree tree;
        tree.set_root(TreeNode::create_leaf(0, read_names[0]));
        return tree;
    }

    // Select algorithm based on method
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

Tree HierarchicalClustering::build_upgma(const Eigen::MatrixXd& dist_matrix,
                                         const std::vector<std::string>& read_names) {
    // UPGMA: distance = average distance of all pairs between two clusters
    auto compute_distance = [](const Eigen::MatrixXd& D, const std::vector<int>& cluster_a,
                               const std::vector<int>& cluster_b, int /*size_a*/, int /*size_b*/
                               ) -> double {
        double sum = 0.0;
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                sum += D(i, j);
            }
        }
        return sum / (cluster_a.size() * cluster_b.size());
    };

    // UPGMA height = distance / 2
    auto distance_to_height = [](double d) { return d / 2.0; };

    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_single(const Eigen::MatrixXd& dist_matrix,
                                          const std::vector<std::string>& read_names) {
    // Single linkage: distance = nearest pair distance
    auto compute_distance = [](const Eigen::MatrixXd& D, const std::vector<int>& cluster_a,
                               const std::vector<int>& cluster_b, int /*size_a*/, int /*size_b*/
                               ) -> double {
        double min_dist = std::numeric_limits<double>::max();
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                min_dist = std::min(min_dist, D(i, j));
            }
        }
        return min_dist;
    };

    // Single linkage height = distance / 2
    auto distance_to_height = [](double d) { return d / 2.0; };

    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_complete(const Eigen::MatrixXd& dist_matrix,
                                            const std::vector<std::string>& read_names) {
    // Complete linkage: distance = farthest pair distance
    auto compute_distance = [](const Eigen::MatrixXd& D, const std::vector<int>& cluster_a,
                               const std::vector<int>& cluster_b, int /*size_a*/, int /*size_b*/
                               ) -> double {
        double max_dist = 0.0;
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                max_dist = std::max(max_dist, D(i, j));
            }
        }
        return max_dist;
    };

    // Complete linkage height = distance / 2
    auto distance_to_height = [](double d) { return d / 2.0; };

    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_ward(const Eigen::MatrixXd& dist_matrix,
                                        const std::vector<std::string>& read_names) {
    // Ward's method: minimize within-cluster variance increase after merge
    // Using special form of Lance-Williams formula
    // Note: Ward requires squared distances

    auto compute_distance = [](const Eigen::MatrixXd& D, const std::vector<int>& cluster_a,
                               const std::vector<int>& cluster_b, int size_a, int size_b) -> double {
        // Calculate distance between two cluster centroids
        // Ward distance increment = (n_a * n_b / (n_a + n_b)) * d(centroid_a, centroid_b)^2
        // Here we approximate using average inter-cluster distance
        double sum = 0.0;
        for (int i : cluster_a) {
            for (int j : cluster_b) {
                double d = D(i, j);
                sum += d * d;  // Ward uses squared distances
            }
        }
        double avg_sq_dist = sum / (cluster_a.size() * cluster_b.size());

        // Ward increment formula
        double n = static_cast<double>(size_a + size_b);
        double weight = (static_cast<double>(size_a) * static_cast<double>(size_b)) / n;
        return weight * avg_sq_dist;
    };

    // Ward height = sqrt(distance_increment) / 2, but often increment is used directly as height
    auto distance_to_height = [](double d) { return std::sqrt(d) / 2.0; };

    return build_generic(dist_matrix, read_names, compute_distance, distance_to_height);
}

Tree HierarchicalClustering::build_generic(
    const Eigen::MatrixXd& dist_matrix, const std::vector<std::string>& read_names,
    std::function<double(const Eigen::MatrixXd& D, const std::vector<int>& cluster_a, const std::vector<int>& cluster_b,
                         int size_a, int size_b)>
        compute_cluster_distance,
    std::function<double(double distance)> distance_to_height) {
    int n = dist_matrix.rows();

    // Initialize: each read is a separate cluster
    std::vector<std::shared_ptr<TreeNode>> nodes;
    std::vector<std::vector<int>> cluster_members;  // Original element indices in each cluster
    std::vector<int> cluster_sizes;
    std::vector<bool> active(n, true);  // Active clusters

    for (int i = 0; i < n; ++i) {
        auto leaf = TreeNode::create_leaf(i, read_names[i]);
        nodes.push_back(leaf);
        cluster_members.push_back({i});
        cluster_sizes.push_back(1);
    }

    // Merge records
    std::vector<MergeRecord> merge_records;

    int next_node_id = n;  // Internal node IDs start from n
    int active_count = n;

    // Iteratively merge until only one cluster remains
    while (active_count > 1) {
        // Find two active clusters with minimum distance
        double min_dist = std::numeric_limits<double>::max();
        int min_i = -1, min_j = -1;

        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            if (!active[i]) continue;
            for (int j = i + 1; j < static_cast<int>(nodes.size()); ++j) {
                if (!active[j]) continue;

                double dist = compute_cluster_distance(dist_matrix, cluster_members[i], cluster_members[j],
                                                       cluster_sizes[i], cluster_sizes[j]);

                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        if (min_i < 0 || min_j < 0) {
            // All remaining clusters have infinite distance between them
            break;
        }

        // Calculate merge height
        double merge_height = distance_to_height(min_dist);

        // Ensure height is monotonically increasing
        double max_child_height = std::max(nodes[min_i]->height, nodes[min_j]->height);
        if (merge_height < max_child_height + config_.min_branch_length) {
            merge_height = max_child_height + config_.min_branch_length;
        }

        // Create new internal node
        auto new_node = TreeNode::create_internal(next_node_id, nodes[min_i], nodes[min_j], merge_height);

        // Record merge
        MergeRecord record;
        record.cluster_i = min_i;
        record.cluster_j = min_j;
        record.distance = min_dist;
        record.new_cluster_id = next_node_id;
        record.size = cluster_sizes[min_i] + cluster_sizes[min_j];
        merge_records.push_back(record);

        // Update cluster information
        std::vector<int> new_members = cluster_members[min_i];
        new_members.insert(new_members.end(), cluster_members[min_j].begin(), cluster_members[min_j].end());
        int new_size = cluster_sizes[min_i] + cluster_sizes[min_j];

        // Mark old clusters as inactive
        active[min_i] = false;
        active[min_j] = false;

        // Add new cluster
        nodes.push_back(new_node);
        cluster_members.push_back(new_members);
        cluster_sizes.push_back(new_size);
        active.push_back(true);

        next_node_id++;
        active_count--;
    }

    // Find last active cluster as root node
    std::shared_ptr<TreeNode> root = nullptr;
    for (int i = static_cast<int>(nodes.size()) - 1; i >= 0; --i) {
        if (active[i]) {
            root = nodes[i];
            break;
        }
    }

    // Construct Tree object
    Tree tree;
    tree.set_root(root);
    tree.set_merge_records(merge_records);

    return tree;
}

double HierarchicalClustering::find_minimum_distance(const Eigen::MatrixXd& D, const std::vector<bool>& active,
                                                     int& min_i, int& min_j) {
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
        case LinkageMethod::UPGMA:
            return "UPGMA";
        case LinkageMethod::WARD:
            return "WARD";
        case LinkageMethod::SINGLE:
            return "SINGLE";
        case LinkageMethod::COMPLETE:
            return "COMPLETE";
        default:
            return "UNKNOWN";
    }
}

LinkageMethod HierarchicalClustering::string_to_method(const std::string& str) {
    std::string upper = str;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    if (upper == "UPGMA" || upper == "AVERAGE") return LinkageMethod::UPGMA;
    if (upper == "WARD" || upper == "WARD.D" || upper == "WARD.D2") return LinkageMethod::WARD;
    if (upper == "SINGLE" || upper == "MIN") return LinkageMethod::SINGLE;
    if (upper == "COMPLETE" || upper == "MAX") return LinkageMethod::COMPLETE;

    // Default to UPGMA
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

    // BFS traversal to find nodes below threshold as cluster roots
    int current_label = 0;
    std::queue<std::shared_ptr<TreeNode>> queue;
    queue.push(root);

    while (!queue.empty()) {
        auto node = queue.front();
        queue.pop();

        // Convert height to distance (height is distance to leaves, merge distance is height * 2)
        double merge_distance = node->height * 2.0;

        if (node->is_leaf() || merge_distance <= distance_threshold) {
            // This node and all its leaves belong to the same cluster
            for (int leaf_idx : node->leaf_indices) {
                labels[leaf_idx] = current_label;
            }
            current_label++;
        } else {
            // Continue exploring downward
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
        // Each leaf is a separate cluster
        std::vector<int> labels(n_leaves);
        std::iota(labels.begin(), labels.end(), 0);
        return labels;
    }

    // Collect all merge heights, sort from high to low
    std::vector<std::pair<double, std::shared_ptr<TreeNode>>> height_nodes;
    auto internal_nodes = tree.get_internal_nodes();
    for (auto& node : internal_nodes) {
        height_nodes.emplace_back(node->height, node);
    }
    std::sort(height_nodes.begin(), height_nodes.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

    // Find cut height to produce exactly num_clusters groups
    // Start from the highest point, gradually lower the cut height
    double cut_height = height_nodes.empty() ? 0.0 : height_nodes[0].first + 1.0;

    // Binary search to find the correct cut height
    for (const auto& [h, node] : height_nodes) {
        // Cutting at height h will split this node
        // Calculate how many clusters the current cut will produce
        auto test_labels = cut_by_distance(tree, h * 2.0 + 0.0001);  // slightly above h
        int n_clusters = *std::max_element(test_labels.begin(), test_labels.end()) + 1;

        if (n_clusters >= num_clusters) {
            cut_height = h * 2.0;
            break;
        }
    }

    return cut_by_distance(tree, cut_height);
}

std::pair<int, std::vector<int>> TreeCutter::find_optimal_clusters(const Tree& tree, const Eigen::MatrixXd& dist_matrix,
                                                                   int min_k, int max_k) {
    if (tree.empty()) {
        return {0, {}};
    }

    int n = tree.num_leaves();
    if (max_k > n) max_k = n;
    if (min_k < 2) min_k = 2;
    if (min_k > max_k) min_k = max_k;

    double best_score = -2.0;  // Silhouette range is [-1, 1]
    int best_k = min_k;
    std::vector<int> best_labels;

    for (int k = min_k; k <= max_k; ++k) {
        auto labels = cut_by_num_clusters(tree, k);

        // Calculate Silhouette Score
        double total_silhouette = 0.0;
        int valid_count = 0;

        for (int i = 0; i < n; ++i) {
            int cluster_i = labels[i];

            // Calculate a(i): average distance within same cluster
            double a_i = 0.0;
            int count_a = 0;
            for (int j = 0; j < n; ++j) {
                if (j != i && labels[j] == cluster_i) {
                    a_i += dist_matrix(i, j);
                    count_a++;
                }
            }
            if (count_a > 0) a_i /= count_a;

            // Calculate b(i): average distance to nearest other cluster
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

            // Calculate silhouette
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

}  // namespace InterSubMod
