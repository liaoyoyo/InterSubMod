#include "core/Config.hpp"
#include "core/DistanceMatrix.hpp"
#include <iostream>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

namespace InterSubMod {

bool Config::validate() const {
    bool valid = true;

    if (tumor_bam_path.empty()) {
        std::cerr << "Error: Tumor BAM path is required." << std::endl;
        valid = false;
    } else {
        // Verify if it's a valid BAM/CRAM/SAM
        samFile* fp = sam_open(tumor_bam_path.c_str(), "r");
        if (fp == NULL) {
             std::cerr << "Error: Cannot open Tumor BAM file: " << tumor_bam_path << std::endl;
             valid = false;
        } else {
            hts_idx_t* idx = sam_index_load(fp, tumor_bam_path.c_str());
            if (idx == NULL) {
                 std::cerr << "Warning: Tumor BAM index not found. Random access may fail." << std::endl;
                 // Not strictly invalid for basic opening, but good to warn
            } else {
                hts_idx_destroy(idx);
            }
            
            sam_hdr_t* hdr = sam_hdr_read(fp);
            if (hdr == NULL) {
                std::cerr << "Error: Cannot read header from Tumor BAM file." << std::endl;
                valid = false;
            } else {
                sam_hdr_destroy(hdr);
            }
            sam_close(fp);
        }
    }
    
    if (!normal_bam_path.empty()) {
        samFile* fp = sam_open(normal_bam_path.c_str(), "r");
        if (fp == NULL) {
             std::cerr << "Error: Cannot open Normal BAM file: " << normal_bam_path << std::endl;
             valid = false;
        } else {
             sam_hdr_t* hdr = sam_hdr_read(fp);
             if (hdr == NULL) {
                 std::cerr << "Error: Cannot read header from Normal BAM file." << std::endl;
                 valid = false;
             } else {
                 sam_hdr_destroy(hdr);
             }
             sam_close(fp);
        }
    }

    if (reference_fasta_path.empty()) {
        std::cerr << "Error: Reference FASTA path is required." << std::endl;
        valid = false;
    } else {
        // Verify FASTA index (.fai)
        faidx_t* fai = fai_load(reference_fasta_path.c_str());
        if (fai == NULL) {
            std::cerr << "Error: Cannot load Reference FASTA (or .fai index missing): " << reference_fasta_path << std::endl;
            valid = false;
        } else {
            fai_destroy(fai);
        }
    }

    if (somatic_vcf_path.empty()) {
        std::cerr << "Error: Somatic VCF path is required." << std::endl;
        valid = false;
    } else {
        // Verify VCF/BCF
        vcfFile* fp = vcf_open(somatic_vcf_path.c_str(), "r");
        if (fp == NULL) {
            std::cerr << "Error: Cannot open Somatic VCF file: " << somatic_vcf_path << std::endl;
            valid = false;
        } else {
            bcf_hdr_t* hdr = bcf_hdr_read(fp);
            if (hdr == NULL) {
                std::cerr << "Error: Cannot read header from Somatic VCF file." << std::endl;
                valid = false;
            } else {
                bcf_hdr_destroy(hdr);
            }
            vcf_close(fp);
        }
    }

    if (window_size_bp <= 0) {
        std::cerr << "Error: window_size_bp must be positive." << std::endl;
        valid = false;
    }

    if (binary_methyl_high <= binary_methyl_low) {
        std::cerr << "Error: binary_methyl_high must be greater than binary_methyl_low." << std::endl;
        valid = false;
    }

    if (binary_methyl_high > 1.0 || binary_methyl_high < 0.0 || 
        binary_methyl_low > 1.0 || binary_methyl_low < 0.0) {
        std::cerr << "Error: Methylation thresholds must be between 0.0 and 1.0." << std::endl;
        valid = false;
    }

    return valid;
}

void Config::print() const {
    std::cout << "--- Configuration ---" << std::endl;
    std::cout << "Tumor BAM: " << tumor_bam_path << std::endl;
    std::cout << "Normal BAM: " << (normal_bam_path.empty() ? "None" : normal_bam_path) << std::endl;
    std::cout << "Reference: " << reference_fasta_path << std::endl;
    std::cout << "Somatic VCF: " << somatic_vcf_path << std::endl;
    std::cout << "Output Dir: " << output_dir << std::endl;
    std::cout << "Window Size: " << window_size_bp << " bp" << std::endl;
    std::cout << "Min MapQ: " << min_mapq << std::endl;
    std::cout << "Min Read Length: " << min_read_length << std::endl;
    std::cout << "Methylation Thresholds: Low=" << binary_methyl_low << ", High=" << binary_methyl_high << std::endl;
    std::cout << "Threads: " << threads << std::endl;
    std::cout << "Distance Metrics: ";
    for (size_t i = 0; i < distance_metrics.size(); ++i) {
        std::cout << (i > 0 ? ", " : "") << DistanceCalculator::metric_to_string(distance_metrics[i]);
    }
    std::cout << std::endl;
    std::cout << "---------------------" << std::endl;
}

} // namespace InterSubMod
