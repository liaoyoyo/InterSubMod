#include "core/Config.hpp"

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <iostream>
#include <memory>

#include "core/DistanceMatrix.hpp"

namespace InterSubMod {

bool Config::validate() const {
    bool valid = true;

    // HTSlib RAII helpers (lambdas required because HTSlib uses macros/function pointers)
    auto close_sam = [](samFile* f) { if (f) sam_close(f); };
    auto close_idx = [](hts_idx_t* i) { if (i) hts_idx_destroy(i); };
    auto close_hdr = [](sam_hdr_t* h) { if (h) sam_hdr_destroy(h); };

    if (tumor_bam_path.empty()) {
        std::cerr << "Error: Tumor BAM path is required." << std::endl;
        valid = false;
    } else {
        // RAII guard: all HTSlib resources released automatically on scope exit
        std::unique_ptr<samFile, decltype(close_sam)> fp(
            sam_open(tumor_bam_path.c_str(), "r"), close_sam);
        if (!fp) {
            std::cerr << "Error: Cannot open Tumor BAM file: " << tumor_bam_path << std::endl;
            valid = false;
        } else {
            std::unique_ptr<hts_idx_t, decltype(close_idx)> idx(
                sam_index_load(fp.get(), tumor_bam_path.c_str()), close_idx);
            if (!idx) {
                std::cerr << "Warning: Tumor BAM index not found. Random access may fail." << std::endl;
            }

            std::unique_ptr<sam_hdr_t, decltype(close_hdr)> hdr(
                sam_hdr_read(fp.get()), close_hdr);
            if (!hdr) {
                std::cerr << "Error: Cannot read header from Tumor BAM file." << std::endl;
                valid = false;
            }
        }  // fp, idx, hdr automatically released here
    }

    if (!normal_bam_path.empty()) {
        std::unique_ptr<samFile, decltype(close_sam)> fp(
            sam_open(normal_bam_path.c_str(), "r"), close_sam);
        if (!fp) {
            std::cerr << "Error: Cannot open Normal BAM file: " << normal_bam_path << std::endl;
            valid = false;
        } else {
            std::unique_ptr<sam_hdr_t, decltype(close_hdr)> hdr(
                sam_hdr_read(fp.get()), close_hdr);
            if (!hdr) {
                std::cerr << "Error: Cannot read header from Normal BAM file." << std::endl;
                valid = false;
            }
        }  // fp, hdr automatically released here
    }

    if (reference_fasta_path.empty()) {
        std::cerr << "Error: Reference FASTA path is required." << std::endl;
        valid = false;
    } else {
        // Verify FASTA index (.fai)
        faidx_t* fai = fai_load(reference_fasta_path.c_str());
        if (fai == NULL) {
            std::cerr << "Error: Cannot load Reference FASTA (or .fai index missing): " << reference_fasta_path
                      << std::endl;
            valid = false;
        } else {
            fai_destroy(fai);
        }
    }

    if (somatic_vcf_path.empty()) {
        std::cerr << "Error: Somatic VCF path is required." << std::endl;
        valid = false;
    } else {
        auto close_vcf = [](vcfFile* f) { if (f) vcf_close(f); };
        auto close_bcf_hdr = [](bcf_hdr_t* h) { if (h) bcf_hdr_destroy(h); };

        // Verify VCF/BCF - RAII guards ensure cleanup on all exit paths
        std::unique_ptr<vcfFile, decltype(close_vcf)> fp(
            vcf_open(somatic_vcf_path.c_str(), "r"), close_vcf);
        if (!fp) {
            std::cerr << "Error: Cannot open Somatic VCF file: " << somatic_vcf_path << std::endl;
            valid = false;
        } else {
            std::unique_ptr<bcf_hdr_t, decltype(close_bcf_hdr)> hdr(
                bcf_hdr_read(fp.get()), close_bcf_hdr);
            if (!hdr) {
                std::cerr << "Error: Cannot read header from Somatic VCF file." << std::endl;
                valid = false;
            }
        }  // fp, hdr automatically released here
    }

    if (window_size_bp <= 0) {
        std::cerr << "Error: window_size_bp must be positive." << std::endl;
        valid = false;
    }

    if (binary_methyl_high <= binary_methyl_low) {
        std::cerr << "Error: binary_methyl_high must be greater than binary_methyl_low." << std::endl;
        valid = false;
    }

    if (binary_methyl_high > 1.0 || binary_methyl_high < 0.0 || binary_methyl_low > 1.0 || binary_methyl_low < 0.0) {
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
    std::cout << "Full Read Span: " << (use_full_read_span ? "Enabled" : "Disabled") << std::endl;
    std::cout << "---------------------" << std::endl;
}

}  // namespace InterSubMod
