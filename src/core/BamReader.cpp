#include "core/BamReader.hpp"

#include <sstream>
#include <stdexcept>

namespace InterSubMod {

BamReader::BamReader(const std::string &bam_path, int n_threads)
    : bam_path_(bam_path), fp_(nullptr), hdr_(nullptr), idx_(nullptr) {
    // Open BAM file
    fp_ = sam_open(bam_path.c_str(), "r");
    if (!fp_) {
        throw std::runtime_error("Failed to open BAM file: " + bam_path);
    }

    // Set decompression threads if requested
    if (n_threads > 1) {
        if (hts_set_threads(fp_, n_threads) != 0) {
            sam_close(fp_);
            throw std::runtime_error("Failed to set threads for BAM: " + bam_path);
        }
    }

    // Read header
    hdr_ = sam_hdr_read(fp_);
    if (!hdr_) {
        sam_close(fp_);
        throw std::runtime_error("Failed to read BAM header: " + bam_path);
    }

    // Load index
    idx_ = sam_index_load(fp_, bam_path.c_str());
    if (!idx_) {
        bam_hdr_destroy(hdr_);
        sam_close(fp_);
        throw std::runtime_error("Failed to load BAM index (.bai): " + bam_path);
    }
}

BamReader::~BamReader() {
    if (idx_) hts_idx_destroy(idx_);
    if (hdr_) bam_hdr_destroy(hdr_);
    if (fp_) sam_close(fp_);
}

BamReader::BamReader(BamReader &&other) noexcept
    : bam_path_(std::move(other.bam_path_)), fp_(other.fp_), hdr_(other.hdr_), idx_(other.idx_) {
    other.fp_ = nullptr;
    other.hdr_ = nullptr;
    other.idx_ = nullptr;
}

BamReader &BamReader::operator=(BamReader &&other) noexcept {
    if (this != &other) {
        // Clean up current resources
        if (idx_) hts_idx_destroy(idx_);
        if (hdr_) bam_hdr_destroy(hdr_);
        if (fp_) sam_close(fp_);

        // Move from other
        bam_path_ = std::move(other.bam_path_);
        fp_ = other.fp_;
        hdr_ = other.hdr_;
        idx_ = other.idx_;

        other.fp_ = nullptr;
        other.hdr_ = nullptr;
        other.idx_ = nullptr;
    }
    return *this;
}

std::vector<bam1_t *> BamReader::fetch_reads(const std::string &chr, int32_t start, int32_t end) {
    std::vector<bam1_t *> reads;

    if (!fp_ || !hdr_ || !idx_) {
        return reads;  // Not initialized
    }

    // Build region string "chr:start-end"
    std::ostringstream region_ss;
    region_ss << chr << ":" << start << "-" << end;
    std::string region_str = region_ss.str();

    // Create iterator for region query
    hts_itr_t *iter = sam_itr_querys(idx_, hdr_, region_str.c_str());
    if (!iter) {
        // Region not found or invalid - return empty vector
        return reads;
    }

    // Fetch all reads in region
    bam1_t *b = bam_init1();
    int ret;
    while ((ret = sam_itr_next(fp_, iter, b)) >= 0) {
        // Duplicate the read (caller will own the memory)
        bam1_t *copy = bam_dup1(b);
        reads.push_back(copy);
    }

    // Clean up
    bam_destroy1(b);
    hts_itr_destroy(iter);

    // Check for read errors (ret < -1 indicates error)
    if (ret < -1) {
        // Error occurred during iteration
        for (auto *r : reads) {
            bam_destroy1(r);
        }
        reads.clear();
    }

    return reads;
}

}  // namespace InterSubMod
