#include <chrono>
#include <iostream>

#include "core/Config.hpp"
#include "core/RegionProcessor.hpp"
#include "utils/ArgParser.hpp"
#include "utils/ResourceMonitor.hpp"
// #include "core/SomaticSnv.hpp"

int main(int argc, char** argv) {
    InterSubMod::Utils::ResourceMonitor monitor;

    InterSubMod::Config config;

    if (!InterSubMod::Utils::ArgParser::parse(argc, argv, config)) {
        return 1;  // Parse failed or help printed
    }

    std::cout << "Validating configuration..." << std::endl;
    if (!config.validate()) {
        std::cerr << "Configuration validation failed." << std::endl;
        return 1;
    }

    config.print();

    // Print debug mode status
    if (config.is_debug()) {
        std::cout << "\n=== DEBUG MODE ENABLED ===" << std::endl;
        std::cout << "Filtered reads will be logged to: " << config.get_debug_output_dir() << std::endl;
        if (config.no_filter_output) {
            std::cout << "No-filter mode: All reads will be output without filtering" << std::endl;
        }
        std::cout << "==========================\n" << std::endl;
    }

    std::cout << "Configuration valid. Starting analysis..." << std::endl;

    try {
        // Use the new Config-based constructor
        InterSubMod::RegionProcessor processor(config);

        std::cout << "[1] Loading SNVs from VCF..." << std::endl;
        int num_snvs = processor.load_snvs_from_vcf(config.somatic_vcf_path);

        if (num_snvs == 0) {
            std::cerr << "No SNVs loaded. Exiting." << std::endl;
            return 1;
        }

        std::cout << "[2] Processing " << num_snvs << " regions..." << std::endl;
        auto t_start = std::chrono::high_resolution_clock::now();

        auto results = processor.process_all_regions(0);  // Process all

        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        std::cout << "[3] Analysis Complete." << std::endl;
        processor.print_summary(results);

        std::cout << "Total Wall-clock time: " << total_time << " ms" << std::endl;
        std::cout << "Output directory: " << config.output_dir << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }

    monitor.print_stats("Total Execution");

    return 0;
}
