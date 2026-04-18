#pragma once

/**
 * @file NormalBaseline.hpp
 * @brief Normal methylation baseline for per-CpG reference values
 *
 * Provides two capabilities for Phase 2A (Dual-BAM analysis):
 * 1. Per-CpG normal methylation statistics (mean, std, coverage)
 * 2. Residual matrix computation (raw - normal_mean) for pure somatic HP ASM
 *
 * Biological motivation:
 * - FP sites have germline ASM (mean=-0.0205, p=1.16e-14)
 * - Subtracting normal baseline removes this germline signal
 * - Remaining HP ASM on residual matrix captures purely somatic differences
 */

#include <Eigen/Dense>
#include <vector>

namespace InterSubMod {

/**
 * @brief Per-CpG methylation statistics from normal reads
 *
 * Built from normal reads in the merged Tumor+Normal matrix.
 * Used to quantify tumor-vs-normal differences and compute residuals.
 */
struct NormalBaseline {
    std::vector<double> per_cpg_mean;      ///< Mean methylation per CpG from normal reads
    std::vector<double> per_cpg_std;       ///< Std deviation per CpG from normal reads
    std::vector<int> per_cpg_coverage;     ///< Number of normal reads covering each CpG
    int total_normal_reads = 0;            ///< Total normal reads contributing
    double overall_mean = 0.0;             ///< Mean of per_cpg_mean (overall normal methylation)
    double mean_coverage = 0.0;            ///< Mean per-CpG coverage from normal reads
    bool valid = false;                    ///< True if >= min_normal_reads normal reads
};

/**
 * @brief Build normal baseline from merged methylation matrix
 *
 * Extracts normal reads (identified by is_tumor_mask) and computes
 * per-CpG methylation statistics.
 *
 * @param raw_matrix       Merged methylation matrix (reads x CpGs), NaN = no coverage
 * @param is_tumor_mask    True for tumor reads, false for normal reads (aligned with rows)
 * @param min_normal_reads Minimum normal reads for valid baseline (default: 5)
 * @param min_cpg_coverage Minimum normal reads per CpG for valid mean (default: 3)
 * @return NormalBaseline with per-CpG statistics
 */
NormalBaseline build_normal_baseline(const Eigen::MatrixXd& raw_matrix,
                                     const std::vector<bool>& is_tumor_mask,
                                     int min_normal_reads = 5,
                                     int min_cpg_coverage = 3);

/**
 * @brief Compute residual matrix by subtracting normal baseline
 *
 * For each cell: residual = raw - normal_mean (if both valid).
 * Cells where normal baseline has insufficient coverage remain unchanged.
 * NaN cells remain NaN.
 *
 * @param raw_matrix Original methylation matrix (reads x CpGs)
 * @param baseline   Normal baseline (must be valid)
 * @return Residual matrix (same dimensions as raw_matrix)
 */
Eigen::MatrixXd compute_residual_matrix(const Eigen::MatrixXd& raw_matrix,
                                        const NormalBaseline& baseline);

}  // namespace InterSubMod
