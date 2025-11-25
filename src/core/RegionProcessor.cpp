#include "core/RegionProcessor.hpp"
#include <fstream>
#include <sstream>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <unordered_set>

namespace InterSubMod {

RegionProcessor::RegionProcessor(
    const std::string& tumor_bam_path,
    const std::string& normal_bam_path,
    const std::string& ref_fasta_path,
    const std::string& output_dir,
    int num_threads,
    int32_t window_size
) : tumor_bam_path_(tumor_bam_path),
    normal_bam_path_(normal_bam_path),
    ref_fasta_path_(ref_fasta_path),
    output_dir_(output_dir),
    num_threads_(num_threads),
    window_size_(window_size) {
    
    // Set OpenMP threads
    omp_set_num_threads(num_threads_);
    
    std::cout << "RegionProcessor initialized with " << num_threads_ 
              << " threads, window_size=±" << window_size_ << "bp" << std::endl;
}

int RegionProcessor::load_snvs(const std::string& snv_table_path) {
    std::ifstream ifs(snv_table_path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open SNV table: " << snv_table_path << std::endl;
        return 0;
    }
    
    snvs_.clear();
    std::string line;
    int line_num = 0;
    
    // Skip header if present
    std::getline(ifs, line);
    line_num++;
    
    // Check if first line is header
    bool has_header = (line.find("chr") != std::string::npos || 
                      line.find("pos") != std::string::npos);
    
    if (!has_header) {
        // First line is data, parse it
        std::istringstream iss(line);
        std::string chr_str;
        uint32_t pos;
        char ref, alt;
        float qual = 0.0f;
        
        if (iss >> chr_str >> pos >> ref >> alt) {
            iss >> qual;  // Optional
            
            SomaticSnv snv;
            snv.snv_id = 0;
            snv.chr_id = chrom_index_.get_or_create_id(chr_str);
            snv.pos = pos;
            snv.ref_base = ref;
            snv.alt_base = alt;
            snv.qual = qual;
            
            snvs_.push_back(snv);
        }
    }
    
    // Parse remaining lines
    while (std::getline(ifs, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string chr_str;
        uint32_t pos;
        char ref, alt;
        float qual = 0.0f;
        
        if (iss >> chr_str >> pos >> ref >> alt) {
            iss >> qual;  // Optional
            
            SomaticSnv snv;
            snv.snv_id = snvs_.size();
            snv.chr_id = chrom_index_.get_or_create_id(chr_str);
            snv.pos = pos;
            snv.ref_base = ref;
            snv.alt_base = alt;
            snv.qual = qual;
            
            snvs_.push_back(snv);
        } else {
            std::cerr << "Failed to parse SNV at line " << line_num << ": " << line << std::endl;
        }
    }
    
    ifs.close();
    std::cout << "Loaded " << snvs_.size() << " SNVs from " << snv_table_path << std::endl;
    return snvs_.size();
}

int RegionProcessor::load_snvs_from_vcf(const std::string& vcf_path) {
    SomaticSnvTable table;
    if (table.load_from_vcf(vcf_path, chrom_index_)) {
        snvs_ = table.all();
        std::cout << "Loaded " << snvs_.size() << " SNVs from VCF: " << vcf_path << std::endl;
        return snvs_.size();
    }
    
    std::cerr << "Failed to load SNVs from VCF: " << vcf_path << std::endl;
    return 0;
}

std::vector<RegionResult> RegionProcessor::process_all_regions(int max_snvs) {
    int num_to_process = (max_snvs > 0 && max_snvs < static_cast<int>(snvs_.size())) 
                         ? max_snvs : snvs_.size();
    
    std::cout << "Processing " << num_to_process << " regions with " 
              << num_threads_ << " threads..." << std::endl;
    
    std::vector<RegionResult> results(num_to_process);
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    // OpenMP parallel loop
    // OpenMP parallel region to manage thread-local resources
    #pragma omp parallel
    {
        // Initialize thread-local readers once per thread
        // This avoids opening/closing the BAM file and loading the index for every region
        BamReader tumor_reader(tumor_bam_path_);
        FastaReader ref_reader(ref_fasta_path_);
        
        // Optional: Normal BAM reader if path is provided
        std::unique_ptr<BamReader> normal_reader;
        if (!normal_bam_path_.empty()) {
            normal_reader = std::make_unique<BamReader>(normal_bam_path_);
        }

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < num_to_process; i++) {
            const auto& snv = snvs_[i];
            std::string chr_name = chrom_index_.get_name(snv.chr_id);
            
            // Logging (protected by critical section)
            #pragma omp critical
            {
                std::cout << "[Thread " << omp_get_thread_num() << "] Processing region " 
                          << i << " (SNV " << chr_name << ":" << snv.pos << ")" << std::endl;
            }
            
            // Process the region using the thread-local readers
            // Note: We pass the readers by reference
            results[i] = process_single_region(snv, i, tumor_reader, ref_reader);
            
            #pragma omp critical
            {
                if (results[i].success) {
                    std::cout << "[Thread " << omp_get_thread_num() << "] ✓ Region " << i 
                              << " completed: " << results[i].num_reads << " reads, " 
                              << results[i].num_cpgs << " CpGs, " 
                              << results[i].elapsed_ms << " ms" << std::endl;
                } else {
                    std::cerr << "[Thread " << omp_get_thread_num() << "] ✗ Region " << i 
                              << " failed: " << results[i].error_message << std::endl;
                }
            }
        }
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    double total_elapsed = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    
    std::cout << "All regions processed in " << total_elapsed << " ms (" 
              << (total_elapsed / num_to_process) << " ms/region)" << std::endl;
    
    return results;
}

RegionResult RegionProcessor::process_single_region(
    const SomaticSnv& snv, 
    int region_id,
    BamReader& bam_reader,
    FastaReader& fasta_reader
) {
    RegionResult result;
    result.region_id = region_id;
    result.snv_id = snv.snv_id;
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    try {
        // Initialize parsers (lightweight, can be recreated or reused)
        ReadParser read_parser;
        MethylationParser methyl_parser;
        MatrixBuilder matrix_builder;
        
        // Define region window
        // Ensure we don't go below 1 or beyond chromosome end (checked later)
        int32_t region_start = static_cast<int32_t>(snv.pos) - static_cast<int32_t>(window_size_);
        int32_t region_end = static_cast<int32_t>(snv.pos) + static_cast<int32_t>(window_size_);
        
        // Get chromosome name and length
        std::string chr_name = chrom_index_.get_name(snv.chr_id);
        int32_t chr_length = fasta_reader.get_chr_length(chr_name);

        // Clamp coordinates
        if (region_start < 1) region_start = 1;
        if (chr_length > 0 && region_end > chr_length) region_end = chr_length;
        
        // Fetch reads from BAM
        // This uses the thread-local reader, which is much faster than re-opening
        auto reads = bam_reader.fetch_reads(chr_name, region_start, region_end);
        
        // Fetch reference sequence for this window
        // Needed for CIGAR parsing and CpG verification
        std::string ref_seq = fasta_reader.fetch_sequence(chr_name, region_start, region_end);
        
        if (ref_seq.empty()) {
            throw std::runtime_error("Failed to fetch reference sequence");
        }
        
        // Process each read
        int read_count = 0;
        std::unordered_set<std::string> processed_read_names;

        for (auto* b : reads) {
            // 1. Filter and Parse Read Info
            if (read_parser.should_keep(b)) {
                ReadInfo info = read_parser.parse(b, read_count, true, snv, ref_seq, region_start);
                
                // Skip reads that don't cover the SNV or have low quality at the SNV site
                if (info.alt_support == AltSupport::UNKNOWN) {
                    continue;
                }
                
                // Optional: Duplicate read name check
                // if (processed_read_names.find(info.read_name) != processed_read_names.end()) {
                //     continue;
                // }
                processed_read_names.insert(info.read_name);

                // 2. Parse Methylation (MM/ML tags)
                auto methyl_calls = methyl_parser.parse_read(b, ref_seq, region_start);
                
                // 3. Add to Matrix Builder
                matrix_builder.add_read(info, methyl_calls);
                read_count++;
            }
        }
        
        // Finalize matrix construction (sort CpGs, fill NaNs)
        matrix_builder.finalize();
        
        result.num_reads = matrix_builder.num_reads();
        result.num_cpgs = matrix_builder.num_cpgs();
        
        // Write output to disk
        RegionWriter writer(output_dir_);
        writer.write_region(
            snv,
            chr_name,
            region_id,
            region_start,
            region_end,
            matrix_builder.get_reads(),
            matrix_builder.get_cpg_positions(),
            matrix_builder.get_matrix(),
            0.0,  // elapsed_ms will be set below
            0.0   // peak_memory_mb not tracked yet
        );
        
        // Cleanup BAM records
        for (auto* r : reads) {
            bam_destroy1(r);
        }
        
        result.success = true;
        
    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = e.what();
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    
    return result;
}

void RegionProcessor::print_summary(const std::vector<RegionResult>& results) const {
    int success_count = 0;
    int total_reads = 0;
    int total_cpgs = 0;
    double total_time = 0.0;
    
    for (const auto& r : results) {
        if (r.success) {
            success_count++;
            total_reads += r.num_reads;
            total_cpgs += r.num_cpgs;
            total_time += r.elapsed_ms;
        }
    }
    
    std::cout << "\n=== Processing Summary ===" << std::endl;
    std::cout << "Total regions: " << results.size() << std::endl;
    std::cout << "Successful: " << success_count << std::endl;
    std::cout << "Failed: " << (results.size() - success_count) << std::endl;
    std::cout << "Total reads processed: " << total_reads << std::endl;
    std::cout << "Total CpG sites found: " << total_cpgs << std::endl;
    std::cout << "Total processing time: " << total_time << " ms" << std::endl;
    std::cout << "Average time per region: " << (total_time / results.size()) << " ms" << std::endl;
    std::cout << "Average reads per region: " << (total_reads / static_cast<double>(success_count)) << std::endl;
    std::cout << "Average CpGs per region: " << (total_cpgs / static_cast<double>(success_count)) << std::endl;
}

} // namespace InterSubMod

