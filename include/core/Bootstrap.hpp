#pragma once

/**
 * @file Bootstrap.hpp
 * @brief Bootstrap stability test for clustering
 *
 * Tests cluster stability by resampling CpG sites and
 * measuring Adjusted Rand Index (ARI) between original
 * and bootstrap cluster assignments.
 */

#include <Eigen/Dense>
#include <random>
#include <vector>

#include "Stats.hpp"

namespace InterSubMod {

/**
 * @brief Configuration for bootstrap test
 */
struct BootstrapConfig {
    int n_iterations = 200;

    // Early stopping: stop if ARI is consistently low or high
    int early_stop_consecutive = 10;
    double early_stop_low_threshold = 0.2;
    double early_stop_high_threshold = 0.9;
    bool use_early_stopping = true;

    // Random seed
    uint64_t seed = 0;
};

/**
 * @brief Bootstrap stability test
 */
class Bootstrap {
public:
    explicit Bootstrap(const BootstrapConfig& config = BootstrapConfig());

    /**
     * @brief Compute bootstrap stability test
     *
     * @param methylation_matrix N x M matrix (reads x CpG sites), NaN for missing
     * @param original_labels Original cluster assignments
     * @param cluster_func Function to re-cluster given a distance matrix
     * @return BootstrapResult
     */
    BootstrapResult compute(const Eigen::MatrixXd& methylation_matrix, const std::vector<int>& original_labels,
                            std::function<std::vector<int>(const Eigen::MatrixXd&)> cluster_func);

    [[deprecated("Use compute() instead")]]
    BootstrapResult run(const Eigen::MatrixXd& methylation_matrix, const std::vector<int>& original_labels,
                        std::function<std::vector<int>(const Eigen::MatrixXd&)> cluster_func) {
        return compute(methylation_matrix, original_labels, cluster_func);
    }

    /**
     * @brief Compute Adjusted Rand Index between two clusterings
     *
     * ARI measures clustering similarity, adjusted for chance.
     * ARI = 1: perfect agreement
     * ARI = 0: random agreement
     * ARI < 0: worse than random
     *
     * @param labels1 First clustering
     * @param labels2 Second clustering
     * @return ARI in [-1, 1]
     */
    static double adjusted_rand_index(const std::vector<int>& labels1, const std::vector<int>& labels2);

    void set_seed(uint64_t seed);

    const BootstrapConfig& config() const {
        return config_;
    }

private:
    BootstrapConfig config_;
    std::mt19937_64 rng_;

    /**
     * @brief Resample CpG sites (columns) with replacement
     */
    Eigen::MatrixXd resample_columns(const Eigen::MatrixXd& matrix);

    /**
     * @brief Compute pairwise L1 distance matrix from methylation values
     */
    Eigen::MatrixXd compute_l1_distance(const Eigen::MatrixXd& matrix);
};

}  // namespace InterSubMod
