#pragma once

/**
 * @file StructureTest.hpp
 * @brief Structural validation tests (PERMANOVA, Dispersion)
 *
 * Implements Phase 3 of significance analysis:
 * - PERMANOVA (Permutational MANOVA) for testing group differences
 * - Dispersion check to detect heterogeneous within-group variance
 * - Read filtering for complete distance matrix
 */

#include <Eigen/Dense>
#include <random>
#include <vector>

#include "Stats.hpp"

namespace InterSubMod {

enum class PermutationMode {
    kUnrestricted,
    kWithinStrata,
};

/**
 * @brief Exchangeability contract for permutation-based structure tests
 *
 * For kWithinStrata, labels are shuffled only among samples sharing the same
 * integer stratum. If no stratum contains at least two different labels, the
 * test is returned as NOT_EVALUABLE instead of silently falling back to an
 * unrestricted permutation.
 */
struct PermutationOptions {
    PermutationMode mode = PermutationMode::kUnrestricted;
    std::vector<int> strata;
};

/**
 * @brief Configuration for structure tests
 */
struct StructureTestConfig {
    // PERMANOVA settings
    int n_permutations = 999;

    // Dispersion check
    double dispersion_alpha = 0.05;  // Threshold for dispersion warning

    // Read filtering
    int min_reads_for_permanova = 5;  // Nmin

    // Random seed
    uint64_t seed = 0;
};

/**
 * @brief Structure test calculator
 *
 * Performs PERMANOVA and dispersion check on distance matrices.
 * Requires valid (complete) distance matrix - invalid pairs must be filtered first.
 */
class StructureTest {
public:
    explicit StructureTest(const StructureTestConfig& config = StructureTestConfig());

    /**
     * @brief Filter reads to obtain complete distance matrix
     *
     * Uses greedy pruning: repeatedly removes read with lowest valid_degree
     * until no invalid pairs remain or n_reads < min_reads.
     *
     * @param dist_matrix Distance matrix (NaN for invalid pairs)
     * @param[out] valid_read_indices Indices of reads to keep
     * @return True if filtering succeeded (enough reads remain)
     */
    bool filter_reads_for_complete_matrix(const Eigen::MatrixXd& dist_matrix, std::vector<int>& valid_read_indices);

    /**
     * @brief Run PERMANOVA test
     *
     * Tests whether group centroids differ significantly.
     * Requires complete distance matrix (no NaN values).
     *
     * @param dist_matrix N x N distance matrix
     * @param group_labels Group assignment for each sample
     * @return PermanovaResult
     */
    PermanovaResult run_permanova(const Eigen::MatrixXd& dist_matrix, const std::vector<int>& group_labels);

    /**
     * @brief Run PERMANOVA with an explicit exchangeability contract
     *
     * @param dist_matrix N x N complete distance matrix
     * @param group_labels Group assignment for each sample
     * @param permutation_options Unrestricted or within-stratum permutation
     * @return PermanovaResult including R-squared and realized permutation audit
     */
    PermanovaResult run_permanova(const Eigen::MatrixXd& dist_matrix, const std::vector<int>& group_labels,
                                  const PermutationOptions& permutation_options);

    /**
     * @brief Check dispersion homogeneity
     *
     * Computes mean distance to centroid for each group.
     * Warns if groups have significantly different dispersions.
     *
     * @param dist_matrix N x N distance matrix
     * @param group_labels Group assignment
     * @return DispersionResult
     */
    DispersionResult check_dispersion(const Eigen::MatrixXd& dist_matrix, const std::vector<int>& group_labels);

    /**
     * @brief Run permutation-based PERMDISP with an exchangeability contract
     *
     * Unlike the legacy analytic overload, distances to the permuted group
     * centroids are recomputed for every permutation.
     */
    DispersionResult check_dispersion(const Eigen::MatrixXd& dist_matrix, const std::vector<int>& group_labels,
                                      const PermutationOptions& permutation_options);

    /**
     * @brief Set random seed
     */
    void set_seed(uint64_t seed);

    const StructureTestConfig& config() const {
        return config_;
    }

private:
    StructureTestConfig config_;
    std::mt19937_64 rng_;

    /**
     * @brief Compute sum of squared distances within and between groups
     *
     * SS_Total = sum of all pairwise squared distances / N
     * SS_Within = sum of squared distances within each group
     * SS_Between = SS_Total - SS_Within
     */
    void compute_ss(const Eigen::MatrixXd& dist_matrix, const std::vector<int>& group_labels, double& ss_total,
                    double& ss_within, double& ss_between);

    /**
     * @brief Compute pseudo-F statistic
     *
     * F = (SS_Between / (k-1)) / (SS_Within / (N-k))
     */
    double compute_pseudo_f(double ss_between, double ss_within, int n, int k);

    /**
     * @brief Compute mean distance to centroid for each group
     *
     * Uses the distance-based formula without explicit centroid coordinates.
     */
    std::vector<double> compute_mean_distances_to_centroid(const Eigen::MatrixXd& dist_matrix,
                                                           const std::vector<int>& group_labels);

    /**
     * @brief Compute one distance-to-centroid value for every sample
     */
    std::vector<double> compute_distances_to_centroid_per_sample(const Eigen::MatrixXd& dist_matrix,
                                                                 const std::vector<int>& group_labels);

    /**
     * @brief Perform one-way ANOVA F-test on dispersions
     */
    double anova_f_test(const std::vector<double>& values, const std::vector<int>& group_labels);

    /**
     * @brief Validate and apply the requested permutation restriction
     */
    bool validate_permutation_options(int n, const std::vector<int>& group_labels,
                                      const PermutationOptions& permutation_options, std::string& invalid_reason);
    void permute_labels(std::vector<int>& labels, const PermutationOptions& permutation_options);
    static const char* permutation_mode_name(PermutationMode mode);
};

}  // namespace InterSubMod
