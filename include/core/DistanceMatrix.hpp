#pragma once

#include <vector>
#include <Eigen/Dense>
#include "Types.hpp"

namespace InterSubMod {

/**
 * @brief Stores pairwise distances between reads within a region.
 * Used as input for hierarchical clustering.
 */
class DistanceMatrix {
public:
    int region_id;
    std::vector<int> read_ids;   ///< Read IDs corresponding to rows/cols
    Eigen::MatrixXd dist_matrix; ///< Symmetric NxN distance matrix
    int min_common_coverage;     ///< Minimum overlapping sites required
    DistanceMetricType metric_type; ///< Metric used (e.g., NHD)

    /**
     * @brief Computes pairwise distances from a MethylationMatrix.
     * 
     * @param methyl_mat The input methylation data.
     * @param type Distance metric (e.g., NHD).
     * @param min_cov Minimum common CpG sites required to calculate a valid distance.
     * @param nan_strategy Strategy to handle pairs with insufficient overlap (e.g., MAX_DIST).
     */
    void compute_from_methylation(const class MethylationMatrix& methyl_mat, DistanceMetricType type, int min_cov, NanDistanceStrategy nan_strategy);
};

} // namespace InterSubMod
