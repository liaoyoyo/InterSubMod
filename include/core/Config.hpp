#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "Types.hpp"

namespace InterSubMod {

struct Config {
    // Input/Output
    std::string tumor_bam_path;
    std::string normal_bam_path;
    std::string reference_fasta_path;
    std::string somatic_vcf_path;
    std::string output_dir = "output";
    std::string pmd_bed_path;

    // Global Parameters
    int window_size_bp = 1000;
    int min_mapq = 20;
    int min_read_length = 1000;
    int min_base_quality = 20;
    
    double binary_methyl_high = 0.8;
    double binary_methyl_low = 0.2;
    
    int min_site_coverage = 10;
    int min_common_coverage = 5; // C_min
    
    NanDistanceStrategy nan_distance_strategy = NanDistanceStrategy::MAX_DIST;
    DistanceMetricType distance_metric = DistanceMetricType::NHD;
    
    bool pmd_gating = true;
    int threads = 1; // Default to 1, auto can be handled in parser

    // Validation method
    bool validate() const;
    
    // Print current config
    void print() const;
};

} // namespace InterSubMod

