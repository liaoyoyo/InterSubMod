#include <chrono>
#include <iostream>

#include "core/RegionProcessor.hpp"

using namespace InterSubMod;

int main(int argc, char** argv) {
    std::cout << "=== Phase 4 & 5 Test: Parallel Processing of SNVs ===" << std::endl << std::endl;

    try {
        // Configuration with defaults
        std::string tumor_bam = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam";
        std::string normal_bam = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam";
        std::string ref_fasta = "/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa";
        std::string output_dir = "/big8_disk/liaoyoyo2001/InterSubMod/output";
        std::string snv_table = "/big8_disk/liaoyoyo2001/InterSubMod/data/test_snvs_32.tsv";

        int num_threads = 4;
        int32_t window_size = 2000;  // ±2000 bp

        // Override from command line args
        if (argc > 1) {
            snv_table = argv[1];
        }
        if (argc > 2) {
            output_dir = argv[2];
        }
        if (argc > 3) {
            num_threads = std::stoi(argv[3]);
        }

        std::cout << "[1] Initializing RegionProcessor..." << std::endl;
        std::cout << "  - SNV Table: " << snv_table << std::endl;
        std::cout << "  - Threads: " << num_threads << std::endl;
        std::cout << "  - Window size: ±" << window_size << " bp" << std::endl;
        std::cout << "  - Output: " << output_dir << std::endl << std::endl;

        RegionProcessor processor(tumor_bam, normal_bam, ref_fasta, output_dir, num_threads, window_size);

        std::cout << "[2] Loading SNV table..." << std::endl;
        int num_snvs = processor.load_snvs(snv_table);
        std::cout << "✓ Loaded " << num_snvs << " SNVs" << std::endl << std::endl;

        if (num_snvs == 0) {
            std::cerr << "✗ No SNVs loaded, exiting" << std::endl;
            return 1;
        }

        // Show first few SNVs
        std::cout << "First 5 SNVs:" << std::endl;
        const auto& snvs = processor.get_snvs();
        for (size_t i = 0; i < std::min(snvs.size(), size_t(5)); i++) {
            std::cout << "  " << i << ". chr" << snvs[i].chr_id << ":" << snvs[i].pos << " " << snvs[i].ref_base << ">"
                      << snvs[i].alt_base << std::endl;
        }
        std::cout << std::endl;

        std::cout << "[3] Processing all regions with " << num_threads << " threads..." << std::endl;
        auto t_start = std::chrono::high_resolution_clock::now();

        auto results = processor.process_all_regions(0);  // Process all loaded SNVs

        auto t_end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        std::cout << "\n[4] Results Summary" << std::endl;
        processor.print_summary(results);

        std::cout << "\n[5] Performance Metrics" << std::endl;
        std::cout << "Wall-clock time: " << total_time << " ms" << std::endl;
        std::cout << "Speedup estimate: " << (total_time / (total_time / num_threads)) << "x" << std::endl;
        std::cout << "Parallel efficiency: " << (100.0 / num_threads) << "%" << std::endl;

        // Show failed regions if any
        int failed = 0;
        for (const auto& r : results) {
            if (!r.success) {
                if (failed == 0) {
                    std::cout << "\n[6] Failed Regions:" << std::endl;
                }
                std::cout << "  Region " << r.region_id << ": " << r.error_message << std::endl;
                failed++;
            }
        }

        if (failed == 0) {
            std::cout << "\n✓ All regions processed successfully!" << std::endl;
        } else {
            std::cout << "\n⚠ " << failed << " regions failed" << std::endl;
        }

        std::cout << "\n✓ Phase 4 & 5 tests completed!" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "✗ Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
