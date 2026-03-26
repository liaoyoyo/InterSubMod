#include "core/DistanceMatrix.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "core/MethylationMatrix.hpp"

namespace InterSubMod {

// ============================================================================
// Distance Calculation Helper Functions
// ============================================================================

/**
 * @brief Calculate Normalized Hamming Distance between two rows.
 *
 * NHD = (number of differences) / (number of common valid sites)
 *
 * @param binary_row_i Binary values for read i (-1 = missing)
 * @param binary_row_j Binary values for read j (-1 = missing)
 * @param min_cov Minimum common coverage required
 * @param[out] common_count Number of common valid sites
 * @return Distance value, or -1.0 if insufficient overlap
 */
static double calculate_nhd(const Eigen::VectorXi& binary_row_i, const Eigen::VectorXi& binary_row_j, int min_cov,
                            int& common_count) {
    common_count = 0;
    int diff_count = 0;
    const int n = binary_row_i.size();

    for (int k = 0; k < n; ++k) {
        int val_i = binary_row_i(k);
        int val_j = binary_row_j(k);

        // Both must be valid (not -1)
        if (val_i != -1 && val_j != -1) {
            common_count++;
            if (val_i != val_j) {
                diff_count++;
            }
        }
    }

    if (common_count < min_cov) {
        return -1.0;
    }

    return static_cast<double>(diff_count) / static_cast<double>(common_count);
}

/**
 * @brief Calculate L1 (Manhattan) Distance between two rows.
 *
 * L1 = mean(|p_i - p_j|) for common valid sites
 */
static double calculate_l1(const Eigen::VectorXd& raw_row_i, const Eigen::VectorXd& raw_row_j, int min_cov,
                           int& common_count) {
    common_count = 0;
    double sum_diff = 0.0;
    const int n = raw_row_i.size();

    for (int k = 0; k < n; ++k) {
        double val_i = raw_row_i(k);
        double val_j = raw_row_j(k);

        if (!std::isnan(val_i) && !std::isnan(val_j)) {
            common_count++;
            sum_diff += std::abs(val_i - val_j);
        }
    }

    if (common_count < min_cov) {
        return -1.0;
    }

    return sum_diff / static_cast<double>(common_count);
}

/**
 * @brief Calculate L2 (Euclidean) Distance between two rows.
 *
 * L2 = sqrt(mean((p_i - p_j)^2)) for common valid sites
 */
static double calculate_l2(const Eigen::VectorXd& raw_row_i, const Eigen::VectorXd& raw_row_j, int min_cov,
                           int& common_count) {
    common_count = 0;
    double sum_sq_diff = 0.0;
    const int n = raw_row_i.size();

    for (int k = 0; k < n; ++k) {
        double val_i = raw_row_i(k);
        double val_j = raw_row_j(k);

        if (!std::isnan(val_i) && !std::isnan(val_j)) {
            common_count++;
            double diff = val_i - val_j;
            sum_sq_diff += diff * diff;
        }
    }

    if (common_count < min_cov) {
        return -1.0;
    }

    return std::sqrt(sum_sq_diff / static_cast<double>(common_count));
}

/**
 * @brief Calculate Pearson Correlation Distance between two rows.
 *
 * CORR = 1 - pearson_correlation(p_i, p_j)
 *
 * Note: Requires at least 3 common sites for meaningful correlation.
 */
static double calculate_corr(const Eigen::VectorXd& raw_row_i, const Eigen::VectorXd& raw_row_j, int min_cov,
                             int& common_count, bool center = true) {
    // Collect common valid values
    std::vector<double> vals_i, vals_j;
    const int n = raw_row_i.size();

    for (int k = 0; k < n; ++k) {
        double val_i = raw_row_i(k);
        double val_j = raw_row_j(k);

        if (!std::isnan(val_i) && !std::isnan(val_j)) {
            vals_i.push_back(val_i);
            vals_j.push_back(val_j);
        }
    }

    common_count = static_cast<int>(vals_i.size());

    // Need at least min_cov and at least 3 for correlation
    if (common_count < min_cov || common_count < 3) {
        return -1.0;
    }

    // Calculate means
    double mean_i = std::accumulate(vals_i.begin(), vals_i.end(), 0.0) / common_count;
    double mean_j = std::accumulate(vals_j.begin(), vals_j.end(), 0.0) / common_count;

    // Calculate correlation
    double sum_prod = 0.0;
    double sum_sq_i = 0.0;
    double sum_sq_j = 0.0;

    for (int k = 0; k < common_count; ++k) {
        double di = center ? (vals_i[k] - mean_i) : vals_i[k];
        double dj = center ? (vals_j[k] - mean_j) : vals_j[k];
        sum_prod += di * dj;
        sum_sq_i += di * di;
        sum_sq_j += dj * dj;
    }

    // Handle edge cases
    if (sum_sq_i < 1e-10 || sum_sq_j < 1e-10) {
        // One or both have zero variance - return max distance
        return 1.0;
    }

    double corr = sum_prod / (std::sqrt(sum_sq_i) * std::sqrt(sum_sq_j));

    // Clamp to [-1, 1] due to floating point errors
    corr = std::max(-1.0, std::min(1.0, corr));

    // Distance = 1 - correlation (range [0, 2])
    // Normalize to [0, 1] for consistency
    return (1.0 - corr) / 2.0;
}

/**
 * @brief Calculate Jaccard Distance between two rows.
 *
 * Jaccard = 1 - |A ∩ B| / |A ∪ B|
 * where A, B are sets of methylated sites (binary = 1).
 *
 * @param include_unmeth If true, also consider unmethylated sites (binary = 0)
 */
static double calculate_jaccard(const Eigen::VectorXi& binary_row_i, const Eigen::VectorXi& binary_row_j, int min_cov,
                                int& common_count, bool include_unmeth = false) {
    int intersection = 0;
    int union_size = 0;
    common_count = 0;
    const int n = binary_row_i.size();

    for (int k = 0; k < n; ++k) {
        int val_i = binary_row_i(k);
        int val_j = binary_row_j(k);

        // Skip if either is missing
        if (val_i == -1 || val_j == -1) {
            continue;
        }

        common_count++;

        if (include_unmeth) {
            // Consider both methylated (1) and unmethylated (0) as features
            if (val_i == val_j) {
                intersection++;
            }
            union_size++;  // Always increment for valid pairs
        } else {
            // Only consider methylated sites
            bool in_i = (val_i == 1);
            bool in_j = (val_j == 1);

            if (in_i || in_j) {
                union_size++;
                if (in_i && in_j) {
                    intersection++;
                }
            }
        }
    }

    if (common_count < min_cov) {
        return -1.0;
    }

    if (union_size == 0) {
        // No methylated sites in either read
        return 0.0;  // Consider them identical
    }

    double jaccard_index = static_cast<double>(intersection) / static_cast<double>(union_size);
    return 1.0 - jaccard_index;
}

/**
 * @brief Calculate Bernoulli Distance (Expected Disagreement with Confidence Weighting).
 *
 * This method calculates distance using raw probability values and applies
 * confidence weighting to reduce the impact of low-confidence sites (p ≈ 0.5).
 *
 * Formula:
 *   delta(p, q) = p(1-q) + (1-p)q  (Expected disagreement rate)
 *   weight(p) = 2 * |p - 0.5|      (Confidence weight)
 *   w_k = weight(p_i) * weight(p_j)
 *   Dist = sum(w_k * delta_k) / sum(w_k)
 *
 * @param raw_row_i Probability values for read i (NaN = missing)
 * @param raw_row_j Probability values for read j (NaN = missing)
 * @param min_cov Minimum number of common valid sites
 * @param[out] common_count Number of common valid sites (non-NaN)
 * @return Distance value [0, 1], or -1.0 if insufficient overlap or zero total weight
 */
static double calculate_bernoulli(const Eigen::VectorXd& raw_row_i, const Eigen::VectorXd& raw_row_j, int min_cov,
                                  int& common_count) {
    common_count = 0;
    double sum_weighted_diff = 0.0;
    double sum_weights = 0.0;
    const int n = raw_row_i.size();

    // Confidence weight function: high confidence at p=0 or p=1, zero at p=0.5
    auto weight_func = [](double p) {
        return 2.0 * std::abs(p - 0.5);
    };

    for (int k = 0; k < n; ++k) {
        double p_i = raw_row_i(k);
        double p_j = raw_row_j(k);

        // Check for missing data (NaN)
        if (!std::isnan(p_i) && !std::isnan(p_j)) {
            common_count++;

            // 1. Calculate Confidence Weights
            double w_i = weight_func(p_i);
            double w_j = weight_func(p_j);
            double w_k = w_i * w_j;

            // 2. Calculate Expected Disagreement (Bernoulli difference)
            // delta = P(states differ) = p_i(1-p_j) + (1-p_i)p_j
            double delta = p_i * (1.0 - p_j) + (1.0 - p_i) * p_j;

            // 3. Accumulate weighted difference
            sum_weighted_diff += w_k * delta;
            sum_weights += w_k;
        }
    }

    // Check minimum coverage requirement
    if (common_count < min_cov) {
        return -1.0;
    }

    // Edge case: If total weight is too small (all overlapping sites are p≈0.5)
    // Treat as "insufficient information" -> Invalid distance
    if (sum_weights < 1e-9) {
        return -1.0;
    }

    // Normalize and return
    return sum_weighted_diff / sum_weights;
}

/**
 * @brief Generic distance calculation dispatcher.
 */
static double calculate_distance_impl(const MethylationMatrix& mat, int row_i, int row_j, const DistanceConfig& config,
                                      int& common_count) {
    switch (config.metric) {
        case DistanceMetricType::NHD: {
            Eigen::VectorXi vec_i = mat.binary_matrix.row(row_i);
            Eigen::VectorXi vec_j = mat.binary_matrix.row(row_j);
            return calculate_nhd(vec_i, vec_j, config.min_common_coverage, common_count);
        }

        case DistanceMetricType::L1: {
            Eigen::VectorXd vec_i = mat.raw_matrix.row(row_i);
            Eigen::VectorXd vec_j = mat.raw_matrix.row(row_j);
            return calculate_l1(vec_i, vec_j, config.min_common_coverage, common_count);
        }

        case DistanceMetricType::L2: {
            Eigen::VectorXd vec_i = mat.raw_matrix.row(row_i);
            Eigen::VectorXd vec_j = mat.raw_matrix.row(row_j);
            return calculate_l2(vec_i, vec_j, config.min_common_coverage, common_count);
        }

        case DistanceMetricType::CORR: {
            Eigen::VectorXd vec_i = mat.raw_matrix.row(row_i);
            Eigen::VectorXd vec_j = mat.raw_matrix.row(row_j);
            return calculate_corr(vec_i, vec_j, config.min_common_coverage, common_count, config.pearson_center);
        }

        case DistanceMetricType::JACCARD: {
            Eigen::VectorXi vec_i = mat.binary_matrix.row(row_i);
            Eigen::VectorXi vec_j = mat.binary_matrix.row(row_j);
            return calculate_jaccard(vec_i, vec_j, config.min_common_coverage, common_count,
                                     config.jaccard_include_unmeth);
        }

        case DistanceMetricType::BERNOULLI: {
            Eigen::VectorXd vec_i = mat.raw_matrix.row(row_i);
            Eigen::VectorXd vec_j = mat.raw_matrix.row(row_j);
            return calculate_bernoulli(vec_i, vec_j, config.min_common_coverage, common_count);
        }

        default:
            common_count = 0;
            return -1.0;
    }
}

// Legacy helper for backward compatibility
static double calculate_distance(const MethylationMatrix& mat, int row_i, int row_j, DistanceMetricType type,
                                 int min_cov) {
    DistanceConfig config;
    config.metric = type;
    config.min_common_coverage = min_cov;
    int common_count;
    return calculate_distance_impl(mat, row_i, row_j, config, common_count);
}

// ============================================================================
// DistanceMatrix Methods
// ============================================================================

void DistanceMatrix::compute_from_methylation(const MethylationMatrix& methyl_mat, DistanceMetricType type, int min_cov,
                                              NanDistanceStrategy nan_strat) {
    DistanceConfig config;
    config.metric = type;
    config.min_common_coverage = min_cov;
    config.nan_strategy = nan_strat;
    compute_from_methylation(methyl_mat, config);
}

void DistanceMatrix::compute_from_methylation(const MethylationMatrix& methyl_mat, const DistanceConfig& config) {
    std::vector<int> all_indices(methyl_mat.num_reads());
    std::iota(all_indices.begin(), all_indices.end(), 0);
    compute_subset(methyl_mat, all_indices, config);
}

void DistanceMatrix::compute_subset(const MethylationMatrix& methyl_mat, const std::vector<int>& row_indices,
                                    DistanceMetricType type, int min_cov, NanDistanceStrategy nan_strat) {
    DistanceConfig config;
    config.metric = type;
    config.min_common_coverage = min_cov;
    config.nan_strategy = nan_strat;
    compute_subset(methyl_mat, row_indices, config);
}

void DistanceMatrix::compute_subset(const MethylationMatrix& methyl_mat, const std::vector<int>& row_indices,
                                    const DistanceConfig& config) {
    // Initialize metadata
    this->region_id = methyl_mat.region_id;
    this->metric_type = config.metric;
    this->min_common_coverage = config.min_common_coverage;
    this->nan_strategy = config.nan_strategy;

    const int n = static_cast<int>(row_indices.size());

    // Copy read IDs
    this->read_ids.resize(n);
    for (int i = 0; i < n; ++i) {
        this->read_ids[i] = methyl_mat.read_ids[row_indices[i]];
    }

    // Initialize distance matrix
    this->dist_matrix = Eigen::MatrixXd::Zero(n, n);

    // Determine NaN replacement value
    double nan_val = config.max_distance_value;
    if (config.nan_strategy == NanDistanceStrategy::SKIP) {
        nan_val = NAN;
    }

    // Statistics accumulators (thread-safe)
    std::atomic<int> valid_pairs{0};
    std::atomic<int> invalid_pairs{0};
    std::atomic<long long> total_common_coverage{0};

    // Set number of threads
    int num_threads = config.num_threads > 0 ? config.num_threads : omp_get_max_threads();

// Parallel computation of upper triangle
#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int common_count = 0;
            double dist = calculate_distance_impl(methyl_mat, row_indices[i], row_indices[j], config, common_count);

            if (dist < 0) {
                // Invalid pair
                dist = nan_val;
                invalid_pairs++;
            } else {
                valid_pairs++;
                total_common_coverage += common_count;
            }

            this->dist_matrix(i, j) = dist;
            this->dist_matrix(j, i) = dist;
        }
    }

    // Store statistics
    this->num_valid_pairs = valid_pairs.load();
    this->num_invalid_pairs = invalid_pairs.load();

    if (this->num_valid_pairs > 0) {
        this->avg_common_coverage = static_cast<double>(total_common_coverage.load()) / this->num_valid_pairs;
    }
}

void DistanceMatrix::write_csv(const std::string& filepath, bool include_header) const {
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filepath << std::endl;
        return;
    }

    // Header
    if (include_header) {
        ofs << "read_id";
        for (int id : read_ids) {
            ofs << "," << id;
        }
        ofs << "\n";
    }

    // Data
    const int n = size();
    for (int i = 0; i < n; ++i) {
        ofs << read_ids[i];
        for (int j = 0; j < n; ++j) {
            ofs << ",";
            double val = dist_matrix(i, j);
            if (std::isnan(val)) {
                ofs << "NA";
            } else {
                ofs << std::fixed << std::setprecision(6) << val;
            }
        }
        ofs << "\n";
    }

    ofs.close();
}

void DistanceMatrix::write_stats(const std::string& filepath) const {
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        return;
    }

    ofs << "Distance Matrix Statistics\n";
    ofs << "==========================\n\n";
    ofs << "Region ID: " << region_id << "\n";
    ofs << "Number of reads: " << size() << "\n";
    ofs << "Metric: " << DistanceCalculator::metric_to_string(metric_type) << "\n";
    ofs << "Min common coverage (C_min): " << min_common_coverage << "\n";
    ofs << "\n";
    ofs << "Valid pairs: " << num_valid_pairs << "\n";
    ofs << "Invalid pairs (insufficient overlap): " << num_invalid_pairs << "\n";

    int total_pairs = (size() * (size() - 1)) / 2;
    if (total_pairs > 0) {
        double valid_ratio = 100.0 * num_valid_pairs / total_pairs;
        ofs << "Valid pair ratio: " << std::fixed << std::setprecision(1) << valid_ratio << "%\n";
    }

    ofs << "Average common coverage: " << std::fixed << std::setprecision(2) << avg_common_coverage << "\n";

    // Distance statistics (excluding diagonal and invalid)
    if (num_valid_pairs > 0) {
        std::vector<double> valid_distances;
        const int n = size();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double d = dist_matrix(i, j);
                if (!std::isnan(d)) {
                    valid_distances.push_back(d);
                }
            }
        }

        if (!valid_distances.empty()) {
            std::sort(valid_distances.begin(), valid_distances.end());

            double sum = std::accumulate(valid_distances.begin(), valid_distances.end(), 0.0);
            double mean = sum / valid_distances.size();

            double sq_sum = 0.0;
            for (double d : valid_distances) {
                sq_sum += (d - mean) * (d - mean);
            }
            double std_dev = std::sqrt(sq_sum / valid_distances.size());

            ofs << "\nDistance Statistics:\n";
            ofs << "  Min: " << std::fixed << std::setprecision(4) << valid_distances.front() << "\n";
            ofs << "  Max: " << std::fixed << std::setprecision(4) << valid_distances.back() << "\n";
            ofs << "  Mean: " << std::fixed << std::setprecision(4) << mean << "\n";
            ofs << "  Std Dev: " << std::fixed << std::setprecision(4) << std_dev << "\n";

            // Percentiles
            size_t n_dist = valid_distances.size();
            ofs << "  25th percentile: " << valid_distances[n_dist / 4] << "\n";
            ofs << "  Median: " << valid_distances[n_dist / 2] << "\n";
            ofs << "  75th percentile: " << valid_distances[3 * n_dist / 4] << "\n";
        }
    }

    ofs.close();
}

// ============================================================================
// DistanceCalculator Methods
// ============================================================================

DistanceMatrix DistanceCalculator::compute(const MethylationMatrix& methyl_mat, const std::vector<ReadInfo>& reads) {
    DistanceMatrix result;
    result.compute_from_methylation(methyl_mat, config_);
    return result;
}

std::pair<DistanceMatrix, DistanceMatrix> DistanceCalculator::compute_strand_specific(
    const MethylationMatrix& methyl_mat, const std::vector<ReadInfo>& reads) {
    // Collect indices by strand
    std::vector<int> forward_indices;
    std::vector<int> reverse_indices;

    for (size_t i = 0; i < reads.size(); ++i) {
        if (reads[i].strand == Strand::FORWARD) {
            forward_indices.push_back(static_cast<int>(i));
        } else if (reads[i].strand == Strand::REVERSE) {
            reverse_indices.push_back(static_cast<int>(i));
        }
        // UNKNOWN strand reads are excluded from strand-specific analysis
    }

    DistanceMatrix forward_matrix;
    DistanceMatrix reverse_matrix;

    if (!forward_indices.empty()) {
        forward_matrix.compute_subset(methyl_mat, forward_indices, config_);
    }

    if (!reverse_indices.empty()) {
        reverse_matrix.compute_subset(methyl_mat, reverse_indices, config_);
    }

    return {forward_matrix, reverse_matrix};
}

std::string DistanceCalculator::metric_to_string(DistanceMetricType type) {
    switch (type) {
        case DistanceMetricType::NHD:
            return "NHD";
        case DistanceMetricType::L1:
            return "L1";
        case DistanceMetricType::L2:
            return "L2";
        case DistanceMetricType::CORR:
            return "CORR";
        case DistanceMetricType::JACCARD:
            return "JACCARD";
        case DistanceMetricType::BERNOULLI:
            return "BERNOULLI";
        default:
            return "UNKNOWN";
    }
}

DistanceMetricType DistanceCalculator::string_to_metric(const std::string& str) {
    std::string upper = str;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    if (upper == "NHD" || upper == "HAMMING") return DistanceMetricType::NHD;
    if (upper == "L1" || upper == "MANHATTAN") return DistanceMetricType::L1;
    if (upper == "L2" || upper == "EUCLIDEAN") return DistanceMetricType::L2;
    if (upper == "CORR" || upper == "CORRELATION" || upper == "PEARSON") return DistanceMetricType::CORR;
    if (upper == "JACCARD") return DistanceMetricType::JACCARD;
    if (upper == "BERNOULLI") return DistanceMetricType::BERNOULLI;

    // Default
    return DistanceMetricType::NHD;
}

}  // namespace InterSubMod
