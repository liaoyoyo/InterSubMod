#pragma once

#include <string>
#include <vector>
#include <iostream>
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
    std::string tumor_bam_path;       ///< Path to Tumor BAM file (Required)
    std::string normal_bam_path;      ///< Path to Normal BAM file (Optional)
    std::string reference_fasta_path; ///< Path to Reference FASTA (Required)
    std::string somatic_vcf_path;     ///< Path to Somatic VCF (Required)
    std::string output_dir = "output";///< Output directory for results
    std::string pmd_bed_path;         ///< Path to PMD annotation BED (Optional)

    // Global Parameters
    int window_size_bp = 1000;        ///< Analysis window size around somatic SNV (±bp)
    int min_mapq = 20;                ///< Minimum Mapping Quality to consider a read
    int min_read_length = 1000;       ///< Minimum Read Length (bp)
    int min_base_quality = 20;        ///< Minimum Base Quality for SNV/CpG sites
    
    double binary_methyl_high = 0.8;  ///< Threshold for methylated (1) call
    double binary_methyl_low = 0.2;   ///< Threshold for unmethylated (0) call
    
    int min_site_coverage = 5;       ///< Minimum reads covering a CpG site to keep it
    int min_common_coverage = 3;      ///< Minimum common CpG sites to calculate distance
    
    NanDistanceStrategy nan_distance_strategy = NanDistanceStrategy::MAX_DIST; ///< Strategy for missing distances
    DistanceMetricType distance_metric = DistanceMetricType::NHD;              ///< Distance metric to use
    
    bool pmd_gating = true;           ///< Whether to exclude CpG sites in PMDs
    int threads = 16;                  ///< Number of threads for parallel processing

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
};

} // namespace InterSubMod
