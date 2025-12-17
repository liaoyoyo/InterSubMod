#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "Types.hpp"

namespace InterSubMod {

/**
 * @brief Configuration structure holding all runtime parameters.
 *
 * Stores paths to input/output files and various analysis thresholds.
 * Validated by both CLI11 (basic checks) and internal validate() method (complex logic).
 */
struct Config {
    // Input/Output
    std::string tumor_bam_path;         ///< Path to Tumor BAM file (Required)
    std::string normal_bam_path;        ///< Path to Normal BAM file (Optional)
    std::string reference_fasta_path;   ///< Path to Reference FASTA (Required)
    std::string somatic_vcf_path;       ///< Path to Somatic VCF (Required)
    std::string output_dir = "output";  ///< Output directory for results
    std::string pmd_bed_path;           ///< Path to PMD annotation BED (Optional)

    // Global Parameters
    int window_size_bp = 1000;   ///< Analysis window size around somatic SNV (±bp)
    int min_mapq = 20;           ///< Minimum Mapping Quality to consider a read
    int min_read_length = 1000;  ///< Minimum Read Length (bp)
    int min_base_quality = 20;   ///< Minimum Base Quality for SNV/CpG sites

    double binary_methyl_high = 0.8;  ///< Threshold for methylated (1) call
    double binary_methyl_low = 0.2;   ///< Threshold for unmethylated (0) call

    int min_site_coverage = 5;    ///< Minimum reads covering a CpG site to keep it
    int min_common_coverage = 3;  ///< Minimum common CpG sites to calculate distance (C_min)

    NanDistanceStrategy nan_distance_strategy = NanDistanceStrategy::MAX_DIST;     ///< Strategy for missing distances
    std::vector<DistanceMetricType> distance_metrics = {DistanceMetricType::NHD};  ///< Distance metrics to use

    bool pmd_gating = true;  ///< Whether to exclude CpG sites in PMDs
    int threads = 16;        ///< Number of threads for parallel processing

    // Distance Matrix Configuration
    bool compute_distance_matrix = true;           ///< Whether to compute read-read distance matrix
    bool output_distance_matrix = true;            ///< Whether to output distance matrix to CSV
    bool output_strand_distance_matrices = true;   ///< Whether to output strand-specific distance matrices
    double max_distance_value = 1.0;               ///< Value for MAX_DIST strategy (normalized metrics)
    bool distance_use_binary = true;               ///< Use binary matrix (true) or raw matrix (false)
    bool distance_pearson_center = true;           ///< Use mean-centered Pearson for CORR metric
    bool distance_jaccard_include_unmeth = false;  ///< Include unmethylated sites in Jaccard

    // Hierarchical Clustering Configuration
    bool compute_clustering = true;        ///< Whether to perform hierarchical clustering
    bool output_tree_files = true;         ///< Whether to output Newick tree files
    std::string linkage_method = "UPGMA";  ///< Linkage method: UPGMA, WARD, SINGLE, COMPLETE
    int clustering_min_reads = 10;         ///< Minimum reads required for clustering
    bool output_linkage_matrix = true;     ///< Whether to output linkage matrix CSV

    // Logging and Debug
    LogLevel log_level = LogLevel::LOG_INFO;  ///< Logging verbosity level
    std::string debug_output_dir = "";        ///< Directory for debug outputs (filtered reads, etc.)
                                              ///< If empty, defaults to output_dir/debug
    bool output_filtered_reads = false;       ///< Output filtered reads in debug mode
    bool no_filter_output = false;            ///< If true, output all reads without filtering (for verification)

    /**
     * @brief Validates configuration logic and file formats.
     *
     * Performs checks that CLI11 cannot handle, such as:
     * - Logical relationships (e.g., high threshold > low threshold)
     * - Deep file format verification using htslib (BAM headers, VCF format, FASTA index)
     *
     * @return true if configuration is valid, false otherwise.
     */
    bool validate() const;

    /**
     * @brief Prints the current configuration to stdout.
     */
    void print() const;

    /**
     * @brief Returns the effective debug output directory.
     */
    std::string get_debug_output_dir() const {
        if (!debug_output_dir.empty()) {
            return debug_output_dir;
        }
        return output_dir + "/debug";
    }

    /**
     * @brief Check if debug mode is enabled.
     */
    bool is_debug() const {
        return log_level >= LogLevel::LOG_DEBUG;
    }
};

}  // namespace InterSubMod
