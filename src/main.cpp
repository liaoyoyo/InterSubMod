#include <iostream>
#include <chrono>
#include "core/Config.hpp"
#include "utils/ArgParser.hpp"
#include "utils/ResourceMonitor.hpp"
#include "core/RegionProcessor.hpp"
// #include "core/SomaticSnv.hpp"

int main(int argc, char** argv) {
    InterSubMod::Utils::ResourceMonitor monitor;

    InterSubMod::Config config;
    
    if (!InterSubMod::Utils::ArgParser::parse(argc, argv, config)) {
        return 1; // Parse failed or help printed
    }

    std::cout << "Validating configuration..." << std::endl;
    if (!config.validate()) {
        std::cerr << "Configuration validation failed." << std::endl;
        return 1;
    }

    config.print();

    std::cout << "Configuration valid. Starting analysis..." << std::endl;
    
    try {
        InterSubMod::RegionProcessor processor(
            config.tumor_bam_path,
            config.normal_bam_path,
            config.reference_fasta_path,
            config.output_dir,
            config.threads,
            config.window_size_bp
        );

        std::cout << "[1] Loading SNVs from VCF..." << std::endl;
        int num_snvs = processor.load_snvs_from_vcf(config.somatic_vcf_path);
        
        if (num_snvs == 0) {
            std::cerr << "No SNVs loaded. Exiting." << std::endl;
            return 1;
        }

        std::cout << "[2] Processing " << num_snvs << " regions..." << std::endl;
        auto t_start = std::chrono::high_resolution_clock::now();
        
        auto results = processor.process_all_regions(0); // Process all

        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        std::cout << "[3] Analysis Complete." << std::endl;
        processor.print_summary(results);
        
        std::cout << "Total Wall-clock time: " << total_time << " ms" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
    
    monitor.print_stats("Total Execution");

    return 0;
}
