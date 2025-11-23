#include "core/MatrixBuilder.hpp"
#include <algorithm>
#include <set>
#include <stdexcept>

namespace InterSubMod {

int MatrixBuilder::add_read(const ReadInfo& read_info, const std::vector<MethylCall>& methyl_calls) {
    if (finalized_) {
        throw std::runtime_error("MatrixBuilder::add_read: Cannot add reads after finalize()");
    }
    
    int read_id = reads_.size();
    reads_.push_back(read_info);
    
    // Store methylation calls in temporary map
    auto& methyl_map = read_methyl_map_[read_id];
    for (const auto& call : methyl_calls) {
        methyl_map[call.ref_pos] = call.probability;
    }
    
    return read_id;
}

void MatrixBuilder::finalize() {
    if (finalized_) {
        return;  // Already finalized
    }
    
    // 1. Collect all unique CpG positions
    std::set<int32_t> unique_cpgs;
    for (const auto& [read_id, methyl_map] : read_methyl_map_) {
        for (const auto& [cpg_pos, prob] : methyl_map) {
            unique_cpgs.insert(cpg_pos);
        }
    }
    
    // 2. Convert to sorted vector (set is already sorted)
    cpg_positions_.assign(unique_cpgs.begin(), unique_cpgs.end());
    
    // 3. Create position -> column index map
    std::map<int32_t, int> pos_to_col;
    for (size_t col = 0; col < cpg_positions_.size(); col++) {
        pos_to_col[cpg_positions_[col]] = col;
    }
    
    // 4. Allocate matrix
    int num_rows = reads_.size();
    int num_cols = cpg_positions_.size();
    matrix_.resize(num_rows);
    for (int r = 0; r < num_rows; r++) {
        matrix_[r].resize(num_cols, -1.0);  // -1.0 indicates no coverage
    }
    
    // 5. Fill in methylation values
    for (const auto& [read_id, methyl_map] : read_methyl_map_) {
        for (const auto& [cpg_pos, prob] : methyl_map) {
            int col = pos_to_col[cpg_pos];
            matrix_[read_id][col] = prob;
        }
    }
    
    // 6. Clear temporary storage
    read_methyl_map_.clear();
    
    finalized_ = true;
}

void MatrixBuilder::clear() {
    reads_.clear();
    read_methyl_map_.clear();
    cpg_positions_.clear();
    matrix_.clear();
    finalized_ = false;
}

} // namespace InterSubMod

