#pragma once

#include "core/Config.hpp"
#include "vendor/CLI11.hpp"
#include <iostream>
#include <map>

namespace InterSubMod {
namespace Utils {

/**
 * @brief Command-line argument parser wrapper around CLI11.
 */
class ArgParser {
public:
    /**
     * @brief Parses command line arguments and populates the Config object.
     * 
     * Uses CLI11 to handle argument parsing, type conversion, and basic validation
     * (e.g., file existence, numeric ranges).
     * 
     * @param argc Argument count.
     * @param argv Argument values.
     * @param config Reference to Config object to populate.
     * @return true if parsing was successful and execution should continue.
     * @return false if parsing failed or help was requested (execution should stop).
     */
    static bool parse(int argc, char** argv, Config& config) {
        CLI::App app{"InterSubMod - Read-level methylation and somatic variant analysis"};

        // Input/Output
        app.add_option("-t,--tumor-bam", config.tumor_bam_path, "Path to Tumor BAM (Required)")
            ->required()
            ->check(CLI::ExistingFile);
            
        app.add_option("-n,--normal-bam", config.normal_bam_path, "Path to Normal BAM (Optional)")
            ->check(CLI::ExistingFile);
            
        app.add_option("-r,--reference", config.reference_fasta_path, "Path to Reference FASTA (Required)")
            ->required()
            ->check(CLI::ExistingFile);
            
        app.add_option("-v,--vcf", config.somatic_vcf_path, "Path to Somatic VCF (Required)")
            ->required()
            ->check(CLI::ExistingFile);
            
        app.add_option("-o,--output-dir", config.output_dir, "Output Directory (Default: output)");

        // Parameters
        app.add_option("-w,--window-size", config.window_size_bp, "Window size in bp (Default: 1000)")
            ->check(CLI::PositiveNumber);
            
        app.add_option("-j,--threads", config.threads, "Number of threads (Default: 1)")
            ->check(CLI::PositiveNumber);

        // Methylation Thresholds (custom check)
        app.add_option("--methyl-high", config.binary_methyl_high, "Binary methylation high threshold")
            ->check(CLI::Range(0.0, 1.0));
            
        app.add_option("--methyl-low", config.binary_methyl_low, "Binary methylation low threshold")
            ->check(CLI::Range(0.0, 1.0));

        // Filter parameters
        app.add_option("--min-mapq", config.min_mapq, "Minimum mapping quality (Default: 20)")
            ->check(CLI::Range(0, 60));
            
        app.add_option("--min-read-length", config.min_read_length, "Minimum read length in bp (Default: 1000)")
            ->check(CLI::PositiveNumber);
            
        app.add_option("--min-base-quality", config.min_base_quality, "Minimum base quality at SNV (Default: 20)")
            ->check(CLI::Range(0, 93));

        // Logging and Debug
        std::string log_level_str = "info";
        app.add_option("--log-level", log_level_str, 
            "Logging level: error, warn, info, debug (Default: info)")
            ->check(CLI::IsMember({"error", "warn", "info", "debug"}, CLI::ignore_case));
            
        app.add_option("--debug-output-dir", config.debug_output_dir,
            "Directory for debug outputs (Default: <output-dir>/debug)");
            
        app.add_flag("--output-filtered-reads", config.output_filtered_reads,
            "Output filtered reads with reasons in debug mode");
            
        app.add_flag("--no-filter", config.no_filter_output,
            "Output all reads without filtering (for verification purposes)");

        try {
            app.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            // If help is requested (ret=0) or error occurs (ret>0), we print message and return false.
            app.exit(e);
            return false;
        }

        // Convert log level string to enum
        static const std::map<std::string, LogLevel> log_level_map = {
            {"error", LogLevel::LOG_ERROR},
            {"warn", LogLevel::LOG_WARN},
            {"info", LogLevel::LOG_INFO},
            {"debug", LogLevel::LOG_DEBUG}
        };
        
        std::string log_lower = log_level_str;
        std::transform(log_lower.begin(), log_lower.end(), log_lower.begin(), ::tolower);
        auto it = log_level_map.find(log_lower);
        if (it != log_level_map.end()) {
            config.log_level = it->second;
        }
        
        // Auto-enable output_filtered_reads in debug mode
        if (config.log_level >= LogLevel::LOG_DEBUG) {
            config.output_filtered_reads = true;
        }

        return true;
    }
};

} // namespace Utils
} // namespace InterSubMod
