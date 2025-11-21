#pragma once

#include <vector>
#include <Eigen/Dense>
#include "DataStructs.hpp"

namespace InterSubMod {

class MethylationMatrix {
public:
    int region_id;
    std::vector<int> read_ids; // Row indices mapping to global Read ID
    std::vector<int> cpg_ids;  // Col indices mapping to global CpG ID
    
    Eigen::MatrixXd raw_matrix;     // 0.0 - 1.0 or NaN
    Eigen::MatrixXi binary_matrix;  // 1, 0, or -1 (missing)
    
    int num_reads() const { return read_ids.size(); }
    int num_sites() const { return cpg_ids.size(); }
    
    // Methods to be implemented
    void build(const std::vector<ReadInfo>& reads, const std::vector<CpGSite>& sites);
};

} // namespace InterSubMod

