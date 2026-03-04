#pragma once

#include <cstdint>
#include <map>
#include <vector>

#include "core/DataStructs.hpp"
#include "core/MethylationParser.hpp"

namespace InterSubMod {

/**
 * @brief Build and manage Read x CpG methylation matrix
 *
 * This class is responsible for:
 * 1. Collecting methylation information from all reads
 * 2. Building unique CpG site list (sorted)
 * 3. Constructing sparse matrix (using NaN for missing data)
 * 4. Providing matrix access and query interface
 *
 * Matrix dimensions:
 * - Rows: reads (in order of reading)
 * - Cols: unique CpG sites (sorted by genomic coordinates)
 *
 * Value meanings:
 * - [0.0, 1.0]: Methylation probability
 * - -1.0: This read does not cover this CpG site (using -1 to represent NaN)
 *
 * Note: Uses vector<vector<double>> instead of Eigen to avoid dependencies
 */
class MatrixBuilder {
public:
    MatrixBuilder() = default;

    /**
     * @brief Add methylation information for a read
     *
     * @param read_info Basic read information (for recording read metadata)
     * @param methyl_calls All methylation calls for this read
     * @return read_id (row index in matrix)
     */
    int add_read(const ReadInfo& read_info, const std::vector<MethylCall>& methyl_calls);

    /**
     * @brief Complete data collection and build final matrix
     *
     * This method will:
     * 1. Sort all unique CpG sites
     * 2. Allocate matrix space (rows x cols)
     * 3. Fill in methylation values
     * 4. Set uncovered positions to NaN
     */
    void finalize();

    /**
     * @brief Get final methylation matrix (read-only)
     * @return Matrix reference (rows=reads, cols=CpGs)
     */
    const std::vector<std::vector<double>>& get_matrix() const {
        return matrix_;
    }

    /**
     * @brief Get all Read information (by row order)
     */
    const std::vector<ReadInfo>& get_reads() const {
        return reads_;
    }

    /**
     * @brief Get all CpG site coordinates (by column order, sorted)
     */
    const std::vector<int32_t>& get_cpg_positions() const {
        return cpg_positions_;
    }

    /**
     * @brief Get matrix dimension information
     */
    int num_reads() const {
        return reads_.size();
    }
    int num_cpgs() const {
        return cpg_positions_.size();
    }

    /**
     * @brief Clear all data (for processing next region)
     */
    void clear();

private:
    std::vector<ReadInfo> reads_;  ///< Read metadata (row indices)

    // Temporary storage: read_id -> vector of (pos, prob)
    // Changed from map<int, map> to map<int, vector> for performance
    // Changed from map<int, map> to vector<vector> for performance (O(1) access)
    std::vector<std::vector<std::pair<int32_t, float>>> read_methyl_data_;

    // Final data
    std::vector<int32_t> cpg_positions_;       ///< Sorted unique CpG positions (column indices)
    std::vector<std::vector<double>> matrix_;  ///< Final matrix (rows x cols), -1.0 = no coverage

    bool finalized_ = false;
};

}  // namespace InterSubMod
