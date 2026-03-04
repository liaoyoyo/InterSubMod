#pragma once

#include <cstdint>
#include <string>

#include "Types.hpp"

namespace InterSubMod {

/**
 * @brief Represents a single CpG site in the genome.
 */
struct CpGSite {
    int cpg_id;                ///< Unique internal ID
    int chr_id;                ///< Chromosome ID (mapped via ChromIndex)
    uint32_t pos;              ///< 1-based genomic position
    bool in_pmd;               ///< True if site is within a Partially Methylated Domain
    bool in_repressive_state;  ///< True if site is in a repressive chromatin state
    bool accessible;           ///< True if site is accessible (ATAC-seq peak)
};

/**
 * @brief Simplified read information extracted from BAM.
 */
struct ReadInfo {
    int read_id;             ///< Unique internal read ID
    std::string read_name;   ///< Original read name (QNAME)
    int chr_id;              ///< Chromosome ID
    int32_t align_start;     ///< Alignment start position (0-based)
    int32_t align_end;       ///< Alignment end position (0-based)
    int mapq;                ///< Mapping Quality
    std::string hp_tag;      ///< Haplotype tag (HP): "1", "2", "1-1", "2-1", "unphase", etc.
    bool is_tumor;           ///< True if from Tumor BAM, False if from Normal BAM
    AltSupport alt_support;  ///< Support for somatic variant (ALT, REF, or UNKNOWN)
    Strand strand;           ///< Strand orientation (FORWARD/+ or REVERSE/-)
};

/**
 * @brief Information about a read that was filtered out.
 *
 * Used in debug mode to record why reads were excluded from analysis.
 */
struct FilteredReadInfo {
    std::string read_name;  ///< Original read name (QNAME)
    int chr_id;             ///< Chromosome ID
    int32_t align_start;    ///< Alignment start position (0-based)
    int32_t align_end;      ///< Alignment end position (0-based)
    int mapq;               ///< Mapping Quality
    Strand strand;          ///< Strand orientation
    FilterReason reasons;   ///< Bitwise OR of all filter reasons
    bool is_tumor;          ///< True if from Tumor BAM

    FilteredReadInfo()
        : chr_id(-1),
          align_start(0),
          align_end(0),
          mapq(0),
          strand(Strand::UNKNOWN),
          reasons(FilterReason::NONE),
          is_tumor(true) {
    }
};

}  // namespace InterSubMod
