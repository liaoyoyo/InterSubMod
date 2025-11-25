#include "io/RegionWriter.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

namespace InterSubMod {

RegionWriter::RegionWriter(const std::string& output_dir)
    : output_dir_(output_dir) {
    // Create output directory if it doesn't exist
    mkdir(output_dir_.c_str(), 0755);
}

std::string RegionWriter::create_region_dir(const std::string& chr_name, int32_t snv_pos, int32_t region_start, int32_t region_end) {
    // Level 1: output_dir/chrName_snvPos
    std::ostringstream level1;
    level1 << output_dir_ << "/" << chr_name << "_" << snv_pos;
    std::string level1_dir = level1.str();
    mkdir(level1_dir.c_str(), 0755);

    // Level 2: output_dir/chrName_snvPos/chrName_start_end
    std::ostringstream level2;
    level2 << level1_dir << "/" << chr_name << "_" << region_start << "_" << region_end;
    std::string level2_dir = level2.str();
    mkdir(level2_dir.c_str(), 0755);

    return level2_dir;
}

void RegionWriter::write_region(
    const SomaticSnv& snv,
    const std::string& chr_name,
    int region_id,
    int32_t region_start,
    int32_t region_end,
    const std::vector<ReadInfo>& reads,
    const std::vector<int32_t>& cpg_positions,
    const std::vector<std::vector<double>>& matrix,
    double elapsed_ms,
    double peak_memory_mb
) {
    std::string region_dir = create_region_dir(chr_name, snv.pos, region_start, region_end);
    
    write_metadata(region_dir, snv, chr_name, region_id, region_start, region_end,
                   reads.size(), cpg_positions.size(), elapsed_ms, peak_memory_mb);
    
    write_reads(region_dir, reads, chr_name);
    write_cpg_sites(region_dir, chr_name, cpg_positions);
    write_matrix_csv(region_dir, matrix, cpg_positions);
}

void RegionWriter::write_metadata(
    const std::string& region_dir,
    const SomaticSnv& snv,
    const std::string& chr_name,
    int region_id,
    int32_t region_start,
    int32_t region_end,
    int num_reads,
    int num_cpgs,
    double elapsed_ms,
    double peak_memory_mb
) {
    std::ofstream ofs(region_dir + "/metadata.txt");
    
    ofs << "Region ID: " << region_id << "\n";
    ofs << "Region: " << chr_name << ":" << region_start << "-" << region_end << "\n";
    ofs << "Region Size: " << (region_end - region_start + 1) << " bp\n";
    ofs << "\n";
    ofs << "SNV ID: " << snv.snv_id << "\n";
    ofs << "SNV Position: " << chr_name << ":" << snv.pos << "\n";
    ofs << "SNV: " << snv.ref_base << " -> " << snv.alt_base << "\n";
    ofs << "SNV Quality: " << snv.qual << "\n";
    ofs << "\n";
    ofs << "Num Reads: " << num_reads << "\n";
    ofs << "Num CpG Sites: " << num_cpgs << "\n";
    ofs << "Matrix Dimensions: " << num_reads << " × " << num_cpgs << "\n";
    ofs << "\n";
    ofs << "Processing Time: " << std::fixed << std::setprecision(2) << elapsed_ms << " ms\n";
    ofs << "Peak Memory: " << std::fixed << std::setprecision(2) << peak_memory_mb << " MB\n";
    
    ofs.close();
}

void RegionWriter::write_reads(
    const std::string& region_dir,
    const std::vector<ReadInfo>& reads,
    const std::string& chr_name
) {
    std::ofstream ofs(region_dir + "/reads.tsv");
    
    // Header
    ofs << "read_id\tread_name\tchr\tstart\tend\tmapq\thp\talt_support\tis_tumor\n";
    
    // Data
    for (const auto& read : reads) {
        ofs << read.read_id << "\t"
            << read.read_name << "\t"
            << chr_name << "\t"
            << read.align_start << "\t"
            << read.align_end << "\t"
            << read.mapq << "\t"
            << read.hp_tag << "\t";
        
        // AltSupport enum
        switch (read.alt_support) {
            case AltSupport::REF:     ofs << "REF"; break;
            case AltSupport::ALT:     ofs << "ALT"; break;
            case AltSupport::UNKNOWN: ofs << "UNKNOWN"; break;
        }
        
        ofs << "\t"
            << (read.is_tumor ? "1" : "0") << "\n";
    }
    
    ofs.close();
}

void RegionWriter::write_cpg_sites(
    const std::string& region_dir,
    const std::string& chr_name,
    const std::vector<int32_t>& cpg_positions
) {
    std::ofstream ofs(region_dir + "/cpg_sites.tsv");
    
    // Header
    ofs << "cpg_id\tchr\tposition\n";
    
    // Data
    for (size_t i = 0; i < cpg_positions.size(); i++) {
        ofs << i << "\t"
            << chr_name << "\t"
            << cpg_positions[i] << "\n";
    }
    
    ofs.close();
}

void RegionWriter::write_matrix_csv(
    const std::string& region_dir,
    const std::vector<std::vector<double>>& matrix,
    const std::vector<int32_t>& cpg_positions
) {
    std::ofstream ofs(region_dir + "/methylation.csv");
    
    // Header: CpG positions
    ofs << "read_id";
    for (auto pos : cpg_positions) {
        ofs << "," << pos;
    }
    ofs << "\n";
    
    // Data: each row is a read
    for (size_t r = 0; r < matrix.size(); r++) {
        ofs << r;
        for (size_t c = 0; c < matrix[r].size(); c++) {
            ofs << ",";
            if (matrix[r][c] < 0.0) {  // -1.0 indicates no coverage
                ofs << "NA";
            } else {
                ofs << std::fixed << std::setprecision(4) << matrix[r][c];
            }
        }
        ofs << "\n";
    }
    
    ofs.close();
}

} // namespace InterSubMod

