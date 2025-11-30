#pragma once

#include <htslib/sam.h>
#include "core/DataStructs.hpp"
#include "core/SomaticSnv.hpp"
#include <utility>

namespace InterSubMod {

/**
 * @brief Configuration for read filtering criteria.
 */
struct ReadFilterConfig {
    int min_mapq = 20;           ///< Minimum mapping quality
    int min_read_length = 1000;  ///< Minimum read length in bp
    int min_base_quality = 20;   ///< Minimum base quality for SNV support
    bool require_mm_ml = true;   ///< Require MM and ML tags
};

/**
 * @brief Result of alt support determination with filter reason.
 */
struct AltSupportResult {
    AltSupport support;
    FilterReason filter_reason;
    
    AltSupportResult(AltSupport s = AltSupport::UNKNOWN, FilterReason r = FilterReason::NONE)
        : support(s), filter_reason(r) {}
};

/**
 * @brief Parser for extracting ReadInfo from BAM records.
 * 
 * This class handles:
 * - Filtering reads based on quality criteria
 * - Extracting basic information (QNAME, FLAG, MAPQ, positions)
 * - Parsing phasing tags (HP, PS)
 * - Determining strand orientation (forward/reverse)
 * - Determining ALT/REF support at SNV positions (requires CIGAR parsing)
 * 
 * Thread-safe: This class is stateless and can be used from multiple threads.
 */
class ReadParser {
public:
    /**
     * @brief Constructs a parser with the given filter configuration.
     */
    explicit ReadParser(const ReadFilterConfig& config = {});
    
    /**
     * @brief Checks if a read passes all filtering criteria.
     * 
     * Filters based on:
     * - FLAG: removes secondary, supplementary, duplicate, unmapped
     * - MAPQ: requires >= min_mapq
     * - Length: requires >= min_read_length
     * - Tags: requires MM and ML if configured
     * 
     * @param b BAM record to check.
     * @return true if read should be kept, false if filtered out.
     */
    bool should_keep(const bam1_t* b) const;
    
    /**
     * @brief Checks if a read passes filtering and returns the reason if not.
     * 
     * @param b BAM record to check.
     * @return pair of (should_keep, filter_reason). If should_keep is true, filter_reason is NONE.
     */
    std::pair<bool, FilterReason> should_keep_with_reason(const bam1_t* b) const;
    
    /**
     * @brief Parses a BAM record into a ReadInfo structure.
     * 
     * @param b BAM record to parse.
     * @param read_id Internal read ID to assign.
     * @param is_tumor true if from tumor BAM, false if from normal.
     * @param anchor_snv The SNV position for determining ALT support.
     * @param ref_seq Reference sequence covering the read's alignment region.
     * @param ref_start_pos 0-based start position of ref_seq.
     * @return Populated ReadInfo structure.
     * 
     * @note AltSupport determination requires valid ref_seq covering the SNV.
     */
    ReadInfo parse(
        const bam1_t* b,
        int read_id,
        bool is_tumor,
        const SomaticSnv& anchor_snv,
        const std::string& ref_seq,
        int32_t ref_start_pos
    ) const;
    
    /**
     * @brief Creates a FilteredReadInfo from a BAM record for debug logging.
     * 
     * @param b BAM record.
     * @param is_tumor true if from tumor BAM.
     * @param reasons Filter reasons that caused this read to be filtered.
     * @return Populated FilteredReadInfo structure.
     */
    FilteredReadInfo create_filtered_info(
        const bam1_t* b,
        bool is_tumor,
        FilterReason reasons
    ) const;
    
    /**
     * @brief Determines strand orientation from BAM FLAG.
     * 
     * @param b BAM record.
     * @return Strand::FORWARD if on positive strand, Strand::REVERSE if on negative strand.
     */
    static Strand determine_strand(const bam1_t* b);
    
    /**
     * @brief Gets the filter configuration.
     */
    const ReadFilterConfig& get_config() const { return config_; }

private:
    ReadFilterConfig config_;
    
    /**
     * @brief Determines if a read supports ALT, REF, or is UNKNOWN at SNV position.
     * 
     * This requires:
     * 1. Read must cover the SNV position
     * 2. Base quality at SNV >= min_base_quality
     * 3. CIGAR traversal to find the read offset corresponding to SNV
     * 
     * @return AltSupport::ALT, REF, or UNKNOWN.
     */
    AltSupport determine_alt_support(
        const bam1_t* b,
        const SomaticSnv& snv,
        const std::string& ref_seq,
        int32_t ref_start_pos
    ) const;
    
    /**
     * @brief Determines alt support with detailed filter reason for debug mode.
     * 
     * @return AltSupportResult containing both support status and filter reason.
     */
    AltSupportResult determine_alt_support_with_reason(
        const bam1_t* b,
        const SomaticSnv& snv,
        const std::string& ref_seq,
        int32_t ref_start_pos
    ) const;
};

} // namespace InterSubMod

