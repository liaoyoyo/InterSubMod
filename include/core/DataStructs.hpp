#pragma once

#include <string>
#include <cstdint>
#include "Types.hpp"

namespace InterSubMod {

/**
 * @brief Represents a single CpG site in the genome.
 */
struct CpGSite {
    int cpg_id;             ///< Unique internal ID
    int chr_id;             ///< Chromosome ID (mapped via ChromIndex)
    int32_t pos;            ///< 1-based genomic position
    bool in_pmd;            ///< True if site is within a Partially Methylated Domain
    bool in_repressive_state; ///< True if site is in a repressive chromatin state
    bool accessible;        ///< True if site is accessible (ATAC-seq peak)
};

/**
 * @brief Simplified read information extracted from BAM.
 */
struct ReadInfo {
    int read_id;            ///< Unique internal read ID
    std::string read_name;  ///< Original read name (QNAME)
    int chr_id;             ///< Chromosome ID
    int32_t align_start;    ///< Alignment start position (0-based)
    int32_t align_end;      ///< Alignment end position (0-based)
    int mapq;               ///< Mapping Quality
    int hp_tag;             ///< Haplotype tag (HP): 0=Unknown, 1=H1, 2=H2
    bool is_tumor;          ///< True if from Tumor BAM, False if from Normal BAM
    AltSupport alt_support; ///< Support for somatic variant (ALT, REF, or UNKNOWN)
};

} // namespace InterSubMod
