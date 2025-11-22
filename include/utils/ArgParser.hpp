#pragma once

#include "core/Config.hpp"
#include "vendor/CLI11.hpp"
#include <iostream>

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

        try {
            app.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            // If help is requested (ret=0) or error occurs (ret>0), we print message and return false.
            return app.exit(e) == 0 && false; // This logic is slightly flawed in catch block
                                              // app.exit() prints the message.
                                              // We want main() to exit.
                                              // Let's just return false.
            return false;
        }

        return true;
    }
};

} // namespace Utils
} // namespace InterSubMod
