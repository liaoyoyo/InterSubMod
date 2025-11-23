#pragma once

#include <string>
#include <htslib/faidx.h>

namespace InterSubMod {

/**
 * @brief RAII wrapper for FASTA file reading with HTSlib.
 * 
 * This class provides access to reference genome sequences using faidx.
 * Sequences are fetched on-demand to minimize memory usage.
 * 
 * Thread-safety: faidx operations are generally thread-safe for read operations,
 * but each thread should ideally have its own FastaReader instance.
 * 
 * Usage:
 *   FastaReader fasta("hg38.fa");
 *   std::string seq = fasta.fetch_sequence("chr17", 7577000, 7579000);
 */
class FastaReader {
public:
    /**
     * @brief Constructs a FASTA reader for the specified file.
     * @param fasta_path Path to the FASTA file (must have a .fai index).
     * @throws std::runtime_error if file cannot be opened or indexed.
     */
    explicit FastaReader(const std::string& fasta_path);
    
    /**
     * @brief Destructor - releases faidx resources.
     */
    ~FastaReader();
    
    // Disable copy, allow move
    FastaReader(const FastaReader&) = delete;
    FastaReader& operator=(const FastaReader&) = delete;
    FastaReader(FastaReader&&) noexcept;
    FastaReader& operator=(FastaReader&&) noexcept;
    
    /**
     * @brief Fetches a subsequence from the reference genome.
     * 
     * @param chr Chromosome name (e.g., "chr17", "17").
     * @param start 0-based inclusive start position.
     * @param end 0-based exclusive end position.
     * @return Uppercase DNA sequence, or empty string if region invalid.
     * 
     * @note The returned sequence is always uppercase (A/C/G/T/N).
     * @note If start >= end or chromosome not found, returns empty string.
     */
    std::string fetch_sequence(const std::string& chr, int32_t start, int32_t end);
    
    /**
     * @brief Gets the length of a chromosome.
     * @param chr Chromosome name.
     * @return Length in bp, or -1 if chromosome not found.
     */
    int64_t get_chr_length(const std::string& chr) const;
    
    /**
     * @brief Checks if the FASTA file was successfully loaded.
     */
    bool is_loaded() const { return fai_ != nullptr; }
    
    /**
     * @brief Gets the path of the loaded FASTA file.
     */
    const std::string& get_path() const { return fasta_path_; }

private:
    std::string fasta_path_;
    faidx_t* fai_;
};

} // namespace InterSubMod

