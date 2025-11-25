#include "core/MatrixBuilder.hpp"
#include <algorithm>
#include <set>
#include <stdexcept>

namespace InterSubMod {

// Assuming MethylCall struct and ReadInfo struct are defined elsewhere and ref_pos in MethylCall is uint32_t
// Assuming MatrixBuilder class definition has:
// std::vector<ReadInfo> reads_;
// std::vector<std::map<uint32_t, double>> read_methyl_map_; // Changed key type
// std::vector<uint32_t> cpg_positions_; // Changed element type
// std::vector<std::vector<double>> matrix_;
// bool finalized_ = false;

int MatrixBuilder::add_read(const ReadInfo& read_info, const std::vector<MethylCall>& methyl_calls) {
    if (finalized_) {
        throw std::runtime_error("MatrixBuilder::add_read: Cannot add reads after finalize()");
    }
    
    int read_id = static_cast<int>(reads_.size());
    reads_.push_back(read_info);
    
    // Ensure vector is large enough
    if (read_methyl_data_.size() <= static_cast<size_t>(read_id)) {
        read_methyl_data_.resize(read_id + 1);
    }

    // Store methylation calls in temporary vector
    std::vector<std::pair<int32_t, float>>& calls = read_methyl_data_[read_id];
    calls.reserve(methyl_calls.size());
    
    for (const auto& call : methyl_calls) {
        calls.emplace_back(call.ref_pos, call.probability);
    }
    
    return read_id;
}

void MatrixBuilder::finalize() {
    if (finalized_) {
        return;  // Already finalized
    }
    
    // 1. Collect all unique CpG positions
    // We iterate through all reads and collect positions
    // 1. Collect all unique CpG positions
    // We iterate through all reads and collect positions into a single vector
    std::vector<int32_t> all_positions;
    size_t total_calls = 0;
    for (const auto& calls : read_methyl_data_) {
        total_calls += calls.size();
    }
    all_positions.reserve(total_calls);

    for (const auto& calls : read_methyl_data_) {
        for (const auto& p : calls) {
            all_positions.push_back(p.first);
        }
    }
    
    // 2. Sort and remove duplicates to get unique positions
    std::sort(all_positions.begin(), all_positions.end());
    auto last = std::unique(all_positions.begin(), all_positions.end());
    all_positions.erase(last, all_positions.end());
    
    cpg_positions_ = std::move(all_positions);
    
    // 3. Map position to column index
    std::unordered_map<int32_t, int> pos_to_col;
    for (size_t i = 0; i < cpg_positions_.size(); i++) {
        pos_to_col[cpg_positions_[i]] = i;
    }
    
    // 4. Allocate matrix
    int num_rows = reads_.size();
    int num_cols = cpg_positions_.size();
    matrix_.resize(num_rows);
    
    // Initialize with -1.0 (NaN/No Coverage)
    // Using parallel loop here could help if matrix is huge
    for (int r = 0; r < num_rows; r++) {
        matrix_[r].resize(num_cols, -1.0);
    }
    
    // 5. Fill in methylation values
    for (size_t read_id = 0; read_id < read_methyl_data_.size(); ++read_id) {
        const auto& calls = read_methyl_data_[read_id];
        std::vector<double>& row = matrix_[read_id];
        for (const auto& p : calls) {
            int32_t pos = p.first;
            float prob = p.second;
            
            auto it = pos_to_col.find(pos);
            if (it != pos_to_col.end()) {
                row[it->second] = prob;
            }
        }
    }
    
    // 6. Clear temporary storage
    read_methyl_data_.clear();
    
    finalized_ = true;
}

void MatrixBuilder::clear() {
    reads_.clear();
    read_methyl_data_.clear();

    cpg_positions_.clear();
    matrix_.clear();
    finalized_ = false;
}

} // namespace InterSubMod

