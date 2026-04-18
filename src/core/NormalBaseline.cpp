#include "core/NormalBaseline.hpp"

#include <cmath>

namespace InterSubMod {

NormalBaseline build_normal_baseline(const Eigen::MatrixXd& raw_matrix,
                                     const std::vector<bool>& is_tumor_mask,
                                     int min_normal_reads,
                                     int min_cpg_coverage) {
    NormalBaseline baseline;

    int n_reads = static_cast<int>(raw_matrix.rows());
    int n_cpgs = static_cast<int>(raw_matrix.cols());

    if (n_reads == 0 || n_cpgs == 0) return baseline;

    // Collect normal read indices
    std::vector<int> normal_indices;
    normal_indices.reserve(n_reads);
    for (int i = 0; i < n_reads; ++i) {
        if (i < static_cast<int>(is_tumor_mask.size()) && !is_tumor_mask[i]) {
            normal_indices.push_back(i);
        }
    }

    baseline.total_normal_reads = static_cast<int>(normal_indices.size());
    if (baseline.total_normal_reads < min_normal_reads) {
        return baseline;  // valid remains false
    }

    // Compute per-CpG statistics from normal reads
    baseline.per_cpg_mean.resize(n_cpgs, NAN);
    baseline.per_cpg_std.resize(n_cpgs, NAN);
    baseline.per_cpg_coverage.resize(n_cpgs, 0);

    double total_mean_sum = 0.0;
    int valid_cpg_count = 0;
    double total_coverage_sum = 0.0;

    for (int j = 0; j < n_cpgs; ++j) {
        double sum = 0.0;
        double sum_sq = 0.0;
        int count = 0;

        for (int idx : normal_indices) {
            double val = raw_matrix(idx, j);
            if (!std::isnan(val)) {
                sum += val;
                sum_sq += val * val;
                count++;
            }
        }

        baseline.per_cpg_coverage[j] = count;
        total_coverage_sum += count;

        if (count >= min_cpg_coverage) {
            double mean = sum / count;
            baseline.per_cpg_mean[j] = mean;
            total_mean_sum += mean;
            valid_cpg_count++;

            if (count > 1) {
                double variance = (sum_sq - sum * sum / count) / (count - 1);
                baseline.per_cpg_std[j] = std::sqrt(std::max(0.0, variance));
            } else {
                baseline.per_cpg_std[j] = 0.0;
            }
        }
        // Otherwise per_cpg_mean[j] and per_cpg_std[j] remain NAN
    }

    // Compute overall statistics
    if (valid_cpg_count > 0) {
        baseline.overall_mean = total_mean_sum / valid_cpg_count;
    }
    baseline.mean_coverage = (n_cpgs > 0) ? total_coverage_sum / n_cpgs : 0.0;

    baseline.valid = true;
    return baseline;
}

Eigen::MatrixXd compute_residual_matrix(const Eigen::MatrixXd& raw_matrix,
                                        const NormalBaseline& baseline) {
    Eigen::MatrixXd residual = raw_matrix;

    int n_reads = static_cast<int>(raw_matrix.rows());
    int n_cpgs = static_cast<int>(raw_matrix.cols());
    int baseline_size = static_cast<int>(baseline.per_cpg_mean.size());

    for (int j = 0; j < n_cpgs && j < baseline_size; ++j) {
        double normal_mean = baseline.per_cpg_mean[j];
        if (std::isnan(normal_mean)) continue;

        for (int i = 0; i < n_reads; ++i) {
            double val = raw_matrix(i, j);
            if (!std::isnan(val)) {
                residual(i, j) = val - normal_mean;
            }
        }
    }

    return residual;
}

}  // namespace InterSubMod
