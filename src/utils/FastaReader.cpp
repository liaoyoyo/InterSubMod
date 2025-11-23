#include "utils/FastaReader.hpp"
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace InterSubMod {

FastaReader::FastaReader(const std::string& fasta_path)
    : fasta_path_(fasta_path), fai_(nullptr) {
    
    // Load FASTA index
    fai_ = fai_load(fasta_path.c_str());
    if (!fai_) {
        throw std::runtime_error("Failed to load FASTA index: " + fasta_path + ".fai");
    }
}

FastaReader::~FastaReader() {
    if (fai_) {
        fai_destroy(fai_);
    }
}

FastaReader::FastaReader(FastaReader&& other) noexcept
    : fasta_path_(std::move(other.fasta_path_)),
      fai_(other.fai_) {
    other.fai_ = nullptr;
}

FastaReader& FastaReader::operator=(FastaReader&& other) noexcept {
    if (this != &other) {
        if (fai_) {
            fai_destroy(fai_);
        }
        fasta_path_ = std::move(other.fasta_path_);
        fai_ = other.fai_;
        other.fai_ = nullptr;
    }
    return *this;
}

std::string FastaReader::fetch_sequence(
    const std::string& chr, 
    int32_t start, 
    int32_t end
) {
    if (!fai_) {
        return "";
    }
    
    // Validate range
    if (start < 0 || end <= start) {
        return "";
    }
    
    // Fetch sequence from FASTA
    // faidx_fetch_seq expects 0-based start and (end-1) for inclusive end
    int len = 0;
    char* seq = faidx_fetch_seq(fai_, chr.c_str(), start, end - 1, &len);
    
    if (!seq || len <= 0) {
        if (seq) free(seq);
        return "";  // Chromosome not found or invalid region
    }
    
    // Convert to std::string and uppercase
    std::string result(seq, len);
    free(seq);
    
    // Convert to uppercase
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    
    return result;
}

int64_t FastaReader::get_chr_length(const std::string& chr) const {
    if (!fai_) {
        return -1;
    }
    
    return faidx_seq_len(fai_, chr.c_str());
}

} // namespace InterSubMod

