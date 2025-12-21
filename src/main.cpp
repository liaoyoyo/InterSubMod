#include <chrono>
#include <iostream>

#include "core/Config.hpp"
#include "core/RegionProcessor.hpp"
#include "utils/ArgParser.hpp"
#include "utils/Logger.hpp"
#include "utils/ResourceMonitor.hpp"
// #include "core/SomaticSnv.hpp"

int main(int argc, char** argv) {
    InterSubMod::Utils::ResourceMonitor monitor;

    InterSubMod::Config config;

    if (!InterSubMod::Utils::ArgParser::parse(argc, argv, config)) {
        return 1;  // Parse failed or help printed
    }

    // Configure Logger
    auto& logger = InterSubMod::Utils::Logger::instance();
    logger.set_log_level(config.log_level);
    if (!config.debug_output_dir.empty()) {
        // Optional: log to a file in debug dir
        // logger.set_log_file(config.debug_output_dir + "/inter_sub_mod.log");
    }

    if (!config.validate()) {
        LOG_ERROR("Configuration validation failed.");
        return 1;
    }

    config.print();

    // Print debug mode status
    if (config.is_debug()) {
        LOG_INFO("\n=== DEBUG MODE ENABLED ===");
        LOG_INFO("Filtered reads will be logged to: " + config.get_debug_output_dir());
        if (config.no_filter_output) {
            LOG_INFO("No-filter mode: All reads will be output without filtering");
        }
        LOG_INFO("==========================\n");
    }

    LOG_INFO("Configuration valid. Starting analysis...");

    try {
        // Use the new Config-based constructor
        InterSubMod::RegionProcessor processor(config);

        InterSubMod::Utils::ScopedLogger main_scope("Main Execution");

        LOG_INFO("[1] Loading SNVs from VCF...");
        int num_snvs = processor.load_snvs_from_vcf(config.somatic_vcf_path);

        if (num_snvs == 0) {
            LOG_ERROR("No SNVs loaded. Exiting.");
            return 1;
        }

        LOG_INFO("[2] Processing " + std::to_string(num_snvs) + " regions...");
        auto t_start = std::chrono::high_resolution_clock::now();

        auto results = processor.process_all_regions(0);  // Process all

        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        std::cout << "[3] Analysis Complete." << std::endl;
        processor.print_summary(results);

        LOG_INFO("Total Wall-clock time: " + std::to_string(total_time) + " ms");
        LOG_INFO("Output directory: " + config.output_dir);

    } catch (const std::exception& e) {
        LOG_ERROR("Fatal error: " + std::string(e.what()));
        return 1;
    }

    monitor.print_stats("Total Execution");

    return 0;
}
