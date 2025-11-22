#pragma once

#include <vector>
#include <Eigen/Dense>
#include "DataStructs.hpp"

namespace InterSubMod {

/**
 * @brief Represents the methylation status of Reads x CpG Sites within a Region.
 * 
 * Stores both raw probability values (0.0-1.0) and binary calls (0/1/-1).
 * Rows correspond to Reads, Columns correspond to CpG Sites.
 */
class MethylationMatrix {
public:
    int region_id;
    std::vector<int> read_ids; ///< Maps row index to global Read ID
    std::vector<int> cpg_ids;  ///< Maps column index to global CpG ID
    
    Eigen::MatrixXd raw_matrix;     ///< Raw methylation probabilities (0.0 - 1.0), NaN for missing
    Eigen::MatrixXi binary_matrix;  ///< Binary methylation status: 1 (meth), 0 (unmeth), -1 (missing)
    
    int num_reads() const { return read_ids.size(); }
    int num_sites() const { return cpg_ids.size(); }
    
    /**
     * @brief Builds the matrix from a set of reads and CpG sites.
     * 
     * Iterates through reads, maps their MM/ML tags to the provided CpG sites,
     * and fills the raw and binary matrices.
     */
    void build(const std::vector<ReadInfo>& reads, const std::vector<CpGSite>& sites);
};

} // namespace InterSubMod
