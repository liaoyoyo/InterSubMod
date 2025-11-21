#pragma once

#include "core/Config.hpp"
#include <getopt.h>
#include <iostream>
#include <cstdlib>
#include <vector>

namespace InterSubMod {
namespace Utils {

class ArgParser {
public:
    static bool parse(int argc, char** argv, Config& config) {
        static struct option long_options[] = {
            {"tumor-bam", required_argument, 0, 't'},
            {"normal-bam", required_argument, 0, 'n'},
            {"reference", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},
            {"output-dir", required_argument, 0, 'o'},
            {"window-size", required_argument, 0, 'w'},
            {"threads", required_argument, 0, 'j'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int opt;
        int option_index = 0;

        while ((opt = getopt_long(argc, argv, "t:n:r:v:o:w:j:h", long_options, &option_index)) != -1) {
            switch (opt) {
                case 't': config.tumor_bam_path = optarg; break;
                case 'n': config.normal_bam_path = optarg; break;
                case 'r': config.reference_fasta_path = optarg; break;
                case 'v': config.somatic_vcf_path = optarg; break;
                case 'o': config.output_dir = optarg; break;
                case 'w': config.window_size_bp = std::atoi(optarg); break;
                case 'j': config.threads = std::atoi(optarg); break;
                case 'h':
                    print_usage(argv[0]);
                    return false;
                case '?':
                    return false;
                default:
                    return false;
            }
        }
        return true;
    }

    static void print_usage(const char* program_name) {
        std::cout << "Usage: " << program_name << " [options]\n"
                  << "Options:\n"
                  << "  -t, --tumor-bam <path>      Path to Tumor BAM (Required)\n"
                  << "  -r, --reference <path>      Path to Reference FASTA (Required)\n"
                  << "  -v, --vcf <path>            Path to Somatic VCF (Required)\n"
                  << "  -n, --normal-bam <path>     Path to Normal BAM (Optional)\n"
                  << "  -o, --output-dir <path>     Output Directory (Default: output)\n"
                  << "  -w, --window-size <int>     Window size in bp (Default: 1000)\n"
                  << "  -j, --threads <int>         Number of threads (Default: 1)\n"
                  << "  -h, --help                  Show this help message\n";
    }
};

} // namespace Utils
} // namespace InterSubMod

