#pragma once

#include <string>
#include <vector>
#include <memory>
#include <htslib/sam.h>

namespace InterSubMod {

/**
 * @brief RAII wrapper for BAM file reading with HTSlib.
 * 
 * This class provides thread-safe access to BAM files for querying reads
 * in specific genomic regions. Each thread should maintain its own instance
 * to avoid file pointer contention.
 * 
 * Usage:
 *   BamReader reader("sample.bam");
 *   auto reads = reader.fetch_reads("chr17", 7577000, 7578000);
 *   // ... process reads ...
 *   for (auto* r : reads) bam_destroy1(r);
 */
class BamReader {
public:
    /**
     * @brief Constructs a BAM reader for the specified file.
     * @param bam_path Path to the BAM file (must have a .bai index).
     * @param n_threads Number of decompression threads (default 1).
     * @throws std::runtime_error if file cannot be opened or indexed.
     */
    explicit BamReader(const std::string& bam_path, int n_threads = 1);
    
    /**
     * @brief Destructor - releases all HTSlib resources.
     */
    ~BamReader();
    
    // Disable copy, allow move
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;
    BamReader(BamReader&&) noexcept;
    BamReader& operator=(BamReader&&) noexcept;
    
    /**
     * @brief Fetches all reads overlapping the specified region.
     * 
     * @param chr Chromosome name (must match BAM header).
     * @param start 0-based inclusive start position.
     * @param end 0-based exclusive end position.
     * @return Vector of bam1_t pointers. Caller must free with bam_destroy1().
     * 
     * @note Returns empty vector if chromosome not found or region invalid.
     * @note Reads are allocated with bam_dup1(), caller owns the memory.
     */
    std::vector<bam1_t*> fetch_reads(const std::string& chr, int32_t start, int32_t end);
    
    /**
     * @brief Gets the BAM header for chromosome name lookups.
     * @return Pointer to the BAM header (valid until BamReader is destroyed).
     */
    const bam_hdr_t* get_header() const { return hdr_; }
    
    /**
     * @brief Checks if the BAM file was successfully opened.
     * @return true if file is open and ready for queries.
     */
    bool is_open() const { return fp_ != nullptr; }
    
    /**
     * @brief Gets the path of the opened BAM file.
     */
    const std::string& get_path() const { return bam_path_; }

private:
    std::string bam_path_;
    samFile* fp_;
    bam_hdr_t* hdr_;
    hts_idx_t* idx_;
};

} // namespace InterSubMod

