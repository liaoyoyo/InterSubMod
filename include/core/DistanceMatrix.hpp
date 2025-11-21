#pragma once

#include <vector>
#include <Eigen/Dense>
#include "Types.hpp"

namespace InterSubMod {

class DistanceMatrix {
public:
    int region_id;
    std::vector<int> read_ids;
    Eigen::MatrixXd dist_matrix; // Symmetric N x N
    int min_common_coverage;
    DistanceMetricType metric_type;

    // Methods
    void compute_from_methylation(const class MethylationMatrix& methyl_mat, DistanceMetricType type, int min_cov, NanDistanceStrategy nan_strategy);
};

} // namespace InterSubMod

