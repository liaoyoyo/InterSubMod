#include "core/Config.hpp"

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <system_error>

#include "core/DistanceMatrix.hpp"

// Build identity injected by CMake (see CMakeLists.txt). Falls back to
// "unknown" so the binary still builds outside the project's CMake setup.
#ifdef ISM_GIT_COMMIT
#define ISM_GIT_COMMIT_STR ISM_GIT_COMMIT
#else
#define ISM_GIT_COMMIT_STR "unknown"
#endif

#ifdef ISM_BUILD_TYPE
#define ISM_BUILD_TYPE_STR ISM_BUILD_TYPE
#else
#define ISM_BUILD_TYPE_STR "unknown"
#endif

namespace InterSubMod {

namespace {

/// Escapes a string for embedding in a JSON document.
std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '\n':
                out += "\\n";
                break;
            case '\r':
                out += "\\r";
                break;
            case '\t':
                out += "\\t";
                break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned char>(c));
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

const char* nan_strategy_to_string(NanDistanceStrategy s) {
    switch (s) {
        case NanDistanceStrategy::MAX_DIST:
            return "MAX_DIST";
        case NanDistanceStrategy::SKIP:
            return "SKIP";
    }
    return "UNKNOWN";
}

const char* log_level_to_string(LogLevel l) {
    switch (l) {
        case LogLevel::LOG_FATAL:
            return "fatal";
        case LogLevel::LOG_ERROR:
            return "error";
        case LogLevel::LOG_WARN:
            return "warn";
        case LogLevel::LOG_INFO:
            return "info";
        case LogLevel::LOG_DEBUG:
            return "debug";
        case LogLevel::LOG_TRACE:
            return "trace";
    }
    return "unknown";
}

/// Current local time as an ISO-8601 string with UTC offset.
std::string iso8601_now() {
    const std::time_t now = std::time(nullptr);
    std::tm tm_buf{};
    if (localtime_r(&now, &tm_buf) == nullptr) {
        return "unknown";
    }
    char buf[40];
    if (std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S%z", &tm_buf) == 0) {
        return "unknown";
    }
    return std::string(buf);
}

}  // namespace

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
    // Fields below were previously computed but never shown; a finished run gave no
    // way to tell whether e.g. --germline-hp-only or --permute-hp-seed had been set.
    std::cout << "Min Base Quality: " << min_base_quality << std::endl;
    std::cout << "Min Common Coverage (C_min): " << min_common_coverage << std::endl;
    std::cout << "NaN Distance Strategy: " << nan_strategy_to_string(nan_distance_strategy) << std::endl;
    std::cout << "Expected Coverage: " << expected_coverage << (expected_coverage == 0.0 ? " (auto-estimate)" : "")
              << std::endl;
    std::cout << "Germline HP Only: " << (germline_hp_only ? "Enabled" : "Disabled") << std::endl;
    std::cout << "Permute HP Seed: " << permute_hp_seed << (permute_hp_seed == 0 ? " (off)" : "") << std::endl;
    std::cout << "Log Level: " << log_level_to_string(log_level) << std::endl;
    std::cout << "(full parameter set is written to " << output_dir << "/run_params.json)" << std::endl;
    std::cout << "---------------------" << std::endl;
}

std::string Config::to_json() const {
    std::ostringstream o;
    o << std::boolalpha;

    o << "{\n";
    o << "  \"schema_name\": \"intersubmod.run_params\",\n";
    o << "  \"schema_version\": \"1.0.0\",\n";

    // --- Build identity: which binary produced this run ---
    o << "  \"build\": {\n";
    o << "    \"git_commit\": \"" << json_escape(ISM_GIT_COMMIT_STR) << "\",\n";
    o << "    \"build_type\": \"" << json_escape(ISM_BUILD_TYPE_STR) << "\",\n";
    o << "    \"compiled_at\": \"" << __DATE__ << " " << __TIME__ << "\"\n";
    o << "  },\n";

    // --- Run identity ---
    o << "  \"run\": {\n";
    o << "    \"started_at\": \"" << iso8601_now() << "\"\n";
    o << "  },\n";

    // --- Input / output ---
    o << "  \"io\": {\n";
    o << "    \"tumor_bam_path\": \"" << json_escape(tumor_bam_path) << "\",\n";
    o << "    \"normal_bam_path\": \"" << json_escape(normal_bam_path) << "\",\n";
    o << "    \"reference_fasta_path\": \"" << json_escape(reference_fasta_path) << "\",\n";
    o << "    \"somatic_vcf_path\": \"" << json_escape(somatic_vcf_path) << "\",\n";
    o << "    \"output_dir\": \"" << json_escape(output_dir) << "\",\n";
    o << "    \"pmd_bed_path\": \"" << json_escape(pmd_bed_path) << "\",\n";
    o << "    \"loh_bed_path\": \"" << json_escape(loh_bed_path) << "\"\n";
    o << "  },\n";

    // --- Windowing, filtering, parallelism ---
    o << "  \"region\": {\n";
    o << "    \"window_size_bp\": " << window_size_bp << ",\n";
    o << "    \"use_full_read_span\": " << use_full_read_span << ",\n";
    o << "    \"threads\": " << threads << "\n";
    o << "  },\n";

    o << "  \"read_filter\": {\n";
    o << "    \"min_mapq\": " << min_mapq << ",\n";
    o << "    \"min_read_length\": " << min_read_length << ",\n";
    o << "    \"min_base_quality\": " << min_base_quality << ",\n";
    o << "    \"no_filter_output\": " << no_filter_output << "\n";
    o << "  },\n";

    // --- Methylation binarization and CpG coverage ---
    o << "  \"methylation\": {\n";
    o << "    \"binary_methyl_high\": " << binary_methyl_high << ",\n";
    o << "    \"binary_methyl_low\": " << binary_methyl_low << ",\n";
    o << "    \"min_site_coverage\": " << min_site_coverage << ",\n";
    o << "    \"pmd_gating\": " << pmd_gating << "\n";
    o << "  },\n";

    // --- Distance matrix ---
    o << "  \"distance\": {\n";
    o << "    \"metrics\": [";
    for (size_t i = 0; i < distance_metrics.size(); ++i) {
        o << (i > 0 ? ", " : "") << "\"" << DistanceCalculator::metric_to_string(distance_metrics[i]) << "\"";
    }
    o << "],\n";
    o << "    \"min_common_coverage\": " << min_common_coverage << ",\n";
    o << "    \"nan_distance_strategy\": \"" << nan_strategy_to_string(nan_distance_strategy) << "\",\n";
    o << "    \"max_distance_value\": " << max_distance_value << ",\n";
    o << "    \"compute_distance_matrix\": " << compute_distance_matrix << ",\n";
    o << "    \"output_distance_matrix\": " << output_distance_matrix << ",\n";
    o << "    \"output_strand_distance_matrices\": " << output_strand_distance_matrices << ",\n";
    o << "    \"distance_use_binary\": " << distance_use_binary << ",\n";
    o << "    \"distance_pearson_center\": " << distance_pearson_center << ",\n";
    o << "    \"distance_jaccard_include_unmeth\": " << distance_jaccard_include_unmeth << "\n";
    o << "  },\n";

    // --- Clustering ---
    o << "  \"clustering\": {\n";
    o << "    \"compute_clustering\": " << compute_clustering << ",\n";
    o << "    \"linkage_method\": \"" << json_escape(linkage_method) << "\",\n";
    o << "    \"clustering_min_reads\": " << clustering_min_reads << ",\n";
    o << "    \"output_tree_files\": " << output_tree_files << ",\n";
    o << "    \"output_linkage_matrix\": " << output_linkage_matrix << "\n";
    o << "  },\n";

    // --- Coverage normalization and haplotype handling ---
    // These three change results substantially yet were invisible in the old
    // stdout summary; they are the main reason this file exists.
    o << "  \"coverage_and_haplotype\": {\n";
    o << "    \"expected_coverage\": " << expected_coverage << ",\n";
    o << "    \"expected_coverage_mode\": \"" << (expected_coverage == 0.0 ? "auto_kde" : "user_specified")
      << "\",\n";
    o << "    \"germline_hp_only\": " << germline_hp_only << ",\n";
    o << "    \"permute_hp_seed\": " << permute_hp_seed << ",\n";
    o << "    \"permute_hp_enabled\": " << (permute_hp_seed != 0) << "\n";
    o << "  },\n";

    // --- Logging ---
    o << "  \"logging\": {\n";
    o << "    \"log_level\": \"" << log_level_to_string(log_level) << "\",\n";
    o << "    \"debug_output_dir\": \"" << json_escape(get_debug_output_dir()) << "\",\n";
    o << "    \"output_filtered_reads\": " << output_filtered_reads << "\n";
    o << "  }\n";

    o << "}\n";
    return o.str();
}

bool Config::write_run_params_json(std::string& error_message) const {
    error_message.clear();

    std::error_code ec;
    std::filesystem::create_directories(output_dir, ec);
    if (ec) {
        error_message = "cannot create output directory '" + output_dir + "': " + ec.message();
        return false;
    }

    const std::string path = output_dir + "/run_params.json";
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) {
        error_message = "cannot open '" + path + "' for writing";
        return false;
    }

    out << to_json();
    out.flush();
    if (!out) {
        error_message = "write failed for '" + path + "'";
        return false;
    }

    return true;
}

}  // namespace InterSubMod
