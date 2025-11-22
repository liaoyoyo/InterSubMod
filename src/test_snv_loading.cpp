#include "core/SomaticSnv.hpp"
#include "utils/Logger.hpp"
#include <iostream>
#include <string>

using namespace InterSubMod;

int main(int argc, char* argv[]) {
    // Setup Logger
    Utils::Logger::instance().set_log_level(Utils::LogLevel::L_DEBUG);
    Utils::Logger::info("Starting SNV Loading Test...");

    if (argc < 2) {
        Utils::Logger::error("Usage: test_snv_loading <vcf_path>");
        return 1;
    }

    std::string vcf_path = argv[1];
    
    ChromIndex chrom_index;
    SomaticSnvTable snv_table;

    Utils::Logger::info("Loading VCF from: " + vcf_path);
    
    if (!snv_table.load_from_vcf(vcf_path, chrom_index)) {
        Utils::Logger::error("Failed to load VCF file.");
        return 1;
    }

    Utils::Logger::info("Successfully loaded SNV Table.");
    Utils::Logger::info("Total SNVs: " + std::to_string(snv_table.size()));

    // Print first 10 SNVs as verification
    const auto& snvs = snv_table.all();
    int limit = 10;
    int count = 0;
    
    std::cout << "\n--- First " << limit << " SNVs ---\n";
    std::cout << "ID\tChr\tPos\tRef\tAlt\tVAF\n";
    
    for (const auto& snv : snvs) {
        if (count >= limit) break;
        std::cout << snv.snv_id << "\t"
                  << chrom_index.get_name(snv.chr_id) << "\t"
                  << snv.pos << "\t"
                  << snv.ref_base << "\t"
                  << snv.alt_base << "\t"
                  << snv.somatic_conf << "\n";
        count++;
    }
    std::cout << "--------------------------\n";

    // Optional: Save to TSV
    std::string output_tsv = "test_snv_output.tsv";
    if (snv_table.save_to_tsv(output_tsv, chrom_index)) {
         Utils::Logger::info("Saved SNV table to " + output_tsv);
    }

    return 0;
}

