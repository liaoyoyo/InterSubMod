#include <iostream>
#include <chrono>
#include <unordered_set>
#include "core/BamReader.hpp"
#include "core/ReadParser.hpp"
#include "core/MethylationParser.hpp"
#include "core/SomaticSnv.hpp"
#include "core/MatrixBuilder.hpp"
#include "utils/FastaReader.hpp"
#include "io/RegionWriter.hpp"

using namespace InterSubMod;

int main() {
    std::cout << "=== Phase 3 Test: Matrix Building & Output ===" << std::endl << std::endl;
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    try {
        // Configuration
        std::string bam_path = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam";
        std::string ref_path = "/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa";
        std::string output_dir = "/big8_disk/liaoyoyo2001/InterSubMod/output";
        
        std::string test_chr = "chr17";
        uint32_t snv_pos = 7578000;  // SNV position
        int32_t window = 2000;       // ±2000 bp
        int32_t region_start = snv_pos - window;
        int32_t region_end = snv_pos + window;
        
        std::cout << "[1] Opening files..." << std::endl;
        BamReader bam_reader(bam_path);
        FastaReader fasta_reader(ref_path);
        std::cout << "✓ Files opened" << std::endl << std::endl;
        
        std::cout << "[2] Fetching reads from " << test_chr << ":" << region_start << "-" << region_end << std::endl;
        auto reads = bam_reader.fetch_reads(test_chr, region_start, region_end);
        std::cout << "✓ Fetched " << reads.size() << " reads" << std::endl << std::endl;
        
        std::cout << "[3] Fetching reference sequence..." << std::endl;
        std::string ref_seq = fasta_reader.fetch_sequence(test_chr, region_start, region_end);
        std::cout << "✓ Fetched " << ref_seq.length() << " bp" << std::endl << std::endl;
        
        std::cout << "[4] Parsing reads and methylation..." << std::endl;
        ReadFilterConfig config;
        config.require_mm_ml = false;  // Don't require for filtering
        config.min_base_quality = 5;   // Lower threshold for test data
        
        ReadParser read_parser(config);
        MethylationParser methyl_parser;
        MatrixBuilder matrix_builder;
        
        // Create a test SNV
        SomaticSnv test_snv;
        test_snv.snv_id = 0;
        test_snv.chr_id = 17;
        test_snv.pos = snv_pos;
        test_snv.ref_base = 'A'; // Adjusted based on BAM inspection
        test_snv.alt_base = 'T';
        test_snv.qual = 100.0f;
        
        int passed = 0;
        int with_methyl = 0;
        std::unordered_set<std::string> processed_read_names;
        
        for (auto* b : reads) {
            if (read_parser.should_keep(b)) {
                ReadInfo info = read_parser.parse(b, passed, true, test_snv, ref_seq, region_start);
                
                if (info.alt_support == AltSupport::UNKNOWN) {
                    continue;
                }
                
                // Skip duplicates
                if (processed_read_names.find(info.read_name) != processed_read_names.end()) {
                    continue;
                }
                processed_read_names.insert(info.read_name);
                
                auto methyl_calls = methyl_parser.parse_read(b, ref_seq, region_start);
                
                matrix_builder.add_read(info, methyl_calls);
                
                passed++;
                if (!methyl_calls.empty()) {
                    with_methyl++;
                }
            }
        }
        
        std::cout << "✓ Parsed " << passed << " reads" << std::endl;
        std::cout << "  - " << with_methyl << " reads with methylation" << std::endl << std::endl;
        
        std::cout << "[5] Building matrix..." << std::endl;
        matrix_builder.finalize();
        
        const auto& matrix = matrix_builder.get_matrix();
        const auto& cpg_pos = matrix_builder.get_cpg_positions();
        
        std::cout << "✓ Matrix dimensions: " << matrix.size() << " reads × " << (matrix.empty() ? 0 : matrix[0].size()) << " CpGs" << std::endl;
        if (!cpg_pos.empty()) {
            std::cout << "  - CpG range: " << cpg_pos.front() << " - " << cpg_pos.back() << std::endl << std::endl;
        } else {
             std::cout << "  - CpG range: None" << std::endl << std::endl;
        }
        
        // Calculate stats
        int na_count = 0;    // -1.0 values (no coverage)
        int zero_count = 0;
        int nonzero_count = 0;
        
        for (size_t r = 0; r < matrix.size(); r++) {
            for (size_t c = 0; c < matrix[r].size(); c++) {
                double val = matrix[r][c];
                if (val < 0.0) {  // -1.0 indicates no coverage
                    na_count++;
                } else if (val == 0.0) {
                    zero_count++;
                } else {
                    nonzero_count++;
                }
            }
        }
        
        int total_cells = matrix.size() * (matrix.empty() ? 0 : matrix[0].size());
        std::cout << "  Matrix statistics:" << std::endl;
        std::cout << "    - NA (no coverage): " << na_count << std::endl;
        std::cout << "    - Zero methylation: " << zero_count << std::endl;
        std::cout << "    - Non-zero methylation: " << nonzero_count << std::endl;
        std::cout << "    - Sparsity: " << (100.0 * na_count / total_cells) << "%" << std::endl << std::endl;
        
        std::cout << "[6] Writing output..." << std::endl;
        
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        
        RegionWriter writer(output_dir);
        writer.write_region(
            test_snv,
            test_chr,
            0,  // region_id
            region_start,
            region_end,
            matrix_builder.get_reads(),
            cpg_pos,
            matrix,
            elapsed_ms,
            0.0  // memory tracking not implemented yet
        );
        
        std::cout << "✓ Output written to " << output_dir << "/region_0000/" << std::endl << std::endl;
        
        std::cout << "[7] Summary" << std::endl;
        std::cout << "  Processing time: " << elapsed_ms << " ms" << std::endl;
        std::cout << "  Reads processed: " << passed << std::endl;
        std::cout << "  CpG sites found: " << cpg_pos.size() << std::endl;
        
        size_t matrix_mem = 0;
        for (const auto& row : matrix) {
            matrix_mem += row.size() * sizeof(double);
        }
        std::cout << "  Matrix size: " << (matrix_mem / 1024.0 / 1024.0) << " MB" << std::endl;
        
        // Cleanup
        for (auto* r : reads) {
            bam_destroy1(r);
        }
        
        std::cout << "\n✓ All Phase 3 tests completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

