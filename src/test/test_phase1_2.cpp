#include <iostream>
#include <iomanip>
#include "core/BamReader.hpp"
#include "core/ReadParser.hpp"
#include "core/MethylationParser.hpp"
#include "core/SomaticSnv.hpp"
#include "utils/FastaReader.hpp"

using namespace InterSubMod;

int main(int argc, char** argv) {
    std::cout << "=== Phase 1 & 2 Functionality Test ===" << std::endl;
    
    // Configuration
    // tumor.bam has MM/ML tags (confirmed by samtools inspection)
    std::string bam_path = "/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam";
    std::string ref_path = "/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa";
    std::string test_chr = "chr17";
    int32_t test_start = 7577000;
    int32_t test_end = 7579000;
    
    try {
        // Test 1: BamReader
        std::cout << "\n[Test 1] BamReader" << std::endl;
        std::cout << "Opening BAM: " << bam_path << std::endl;
        BamReader bam_reader(bam_path);
        std::cout << "✓ BAM opened successfully" << std::endl;
        
        // Test 2: FastaReader
        std::cout << "\n[Test 2] FastaReader" << std::endl;
        std::cout << "Opening FASTA: " << ref_path << std::endl;
        FastaReader fasta_reader(ref_path);
        std::cout << "✓ FASTA opened successfully" << std::endl;
        
        // Test 3: Fetch reads
        std::cout << "\n[Test 3] Fetch Reads from " << test_chr << ":" << test_start << "-" << test_end << std::endl;
        auto reads = bam_reader.fetch_reads(test_chr, test_start, test_end);
        std::cout << "✓ Fetched " << reads.size() << " reads" << std::endl;
        
        if (reads.empty()) {
            std::cout << "⚠ No reads in test region, exiting early" << std::endl;
            return 0;
        }
        
        // Test 4: Fetch reference sequence
        std::cout << "\n[Test 4] Fetch Reference Sequence" << std::endl;
        std::string ref_seq = fasta_reader.fetch_sequence(test_chr, test_start, test_end);
        std::cout << "✓ Fetched " << ref_seq.length() << " bp sequence" << std::endl;
        std::cout << "  First 50 bp: " << ref_seq.substr(0, 50) << std::endl;
        
        // Count CpG sites in reference
        int cpg_count = 0;
        for (size_t i = 0; i + 1 < ref_seq.size(); i++) {
            if (ref_seq[i] == 'C' && ref_seq[i+1] == 'G') {
                cpg_count++;
            }
        }
        std::cout << "  CpG sites in reference: " << cpg_count << std::endl;
        
        // Test 5: ReadParser & MethylationParser
        std::cout << "\n[Test 5] ReadParser & MethylationParser (full pipeline)" << std::endl;
        
        // Configure parser to NOT require MM/ML tags for filtering
        // (we'll check for them separately to see which reads have them)
        ReadFilterConfig config;
        config.require_mm_ml = false;
        
        ReadParser read_parser(config);
        MethylationParser methyl_parser;
        
        // Create a dummy SNV for testing
        SomaticSnv test_snv;
        test_snv.snv_id = 0;
        test_snv.chr_id = 1;  // chr17
        test_snv.pos = (test_start + test_end) / 2;  // Middle of region
        test_snv.ref_base = 'C';
        test_snv.alt_base = 'T';
        
        int passed_filter = 0;
        int has_mm_ml = 0;
        int total_methyl_calls = 0;
        
        std::cout << "Processing first 10 reads..." << std::endl;
        
        for (size_t i = 0; i < std::min(reads.size(), size_t(10)); i++) {
            bam1_t* b = reads[i];
            
            if (read_parser.should_keep(b)) {
                passed_filter++;
                
                // Parse read info
                ReadInfo info = read_parser.parse(b, i, true, test_snv, ref_seq, test_start);
                
                std::cout << "\nRead " << i << ": " << info.read_name << std::endl;
                std::cout << "  Position: " << test_chr << ":" << info.align_start << "-" << info.align_end << std::endl;
                std::cout << "  MAPQ: " << info.mapq << std::endl;
                std::cout << "  HP tag: " << info.hp_tag << std::endl;
                std::cout << "  AltSupport: ";
                switch (info.alt_support) {
                    case AltSupport::ALT: std::cout << "ALT"; break;
                    case AltSupport::REF: std::cout << "REF"; break;
                    case AltSupport::UNKNOWN: std::cout << "UNKNOWN"; break;
                }
                std::cout << std::endl;
                
                // Test 6: MethylationParser
                auto methyl_calls = methyl_parser.parse_read(b, ref_seq, test_start);
                if (!methyl_calls.empty()) {
                    has_mm_ml++;
                    total_methyl_calls += methyl_calls.size();
                    
                    std::cout << "  Methylation: " << methyl_calls.size() << " CpG sites" << std::endl;
                    
                    // Show first few methylation calls
                    for (size_t j = 0; j < std::min(methyl_calls.size(), size_t(5)); j++) {
                        std::cout << "    CpG@" << methyl_calls[j].ref_pos 
                                  << ": " << std::fixed << std::setprecision(3) 
                                  << methyl_calls[j].probability << std::endl;
                    }
                    if (methyl_calls.size() > 5) {
                        std::cout << "    ... (" << (methyl_calls.size() - 5) << " more)" << std::endl;
                    }
                }
            }
        }
        
        // Summary
        std::cout << "\n=== Summary ===" << std::endl;
        std::cout << "Total reads fetched: " << reads.size() << std::endl;
        std::cout << "Passed filter: " << passed_filter << " / 10 processed" << std::endl;
        std::cout << "Reads with MM/ML tags: " << has_mm_ml << std::endl;
        std::cout << "Total methylation calls: " << total_methyl_calls << std::endl;
        if (has_mm_ml > 0) {
            std::cout << "Average CpGs per read: " << (total_methyl_calls / (float)has_mm_ml) << std::endl;
        }
        
        // Cleanup
        for (auto* r : reads) {
            bam_destroy1(r);
        }
        
        std::cout << "\n✓ All tests completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

