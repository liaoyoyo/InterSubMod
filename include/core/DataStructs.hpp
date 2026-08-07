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
    std::string hp_tag;      ///< Haplotype tag (HP): "1","2"=germline; "1-1","2-1"=germline-anchored somatic; "3"=HP3 (touches a somatic ALT, germline phasing unconfirmed); "0"=unphased/germline-conflict. NOTE: ReadParser emits "0" for the unphased case, never the literal string "unphase" (verified genome-wide 2026-07-04, V1).
                             ///< May be demoted to "0" when Config::germline_hp_only is set.
    std::string hp_tag_raw;  ///< Raw HP tag before any germline-only demotion (audit). Mirrors hp_tag when flag is off.
    bool is_tumor;           ///< True if from Tumor BAM, False if from Normal BAM
    AltSupport alt_support;  ///< Support for somatic variant (ALT, REF, or UNKNOWN)
    Strand strand;           ///< Strand orientation (FORWARD/+ or REVERSE/-)

    // ── Lineage tags (written by pipeline/lineage/bin/tag_bam) ────────────────
    // SAM spec reserves tags containing lowercase letters for local use.
    // Empty string means the tag was absent on the alignment.
    int64_t phase_set;       ///< PS:i phase set; -1 when absent
    std::string lineage_component;  ///< lc:Z unit_id of the read-linked component
    std::string lineage_block;      ///< lu:Z block_id
    std::string lineage_path;       ///< lv:Z hierarchical path, e.g. "HP2-1-1"
    std::string lineage_pattern;    ///< lp:Z observed R/A/X pattern for the block
    std::string mutation_order;     ///< lo:Z acquired positions along the path, ">"-separated
    char lineage_status;     ///< ls:A  'U' unique | 'M' tie | 'P' partial | 'A' abstain | '\0' absent

    /// True when this read carries a topologically unique lineage assignment.
    /// Anything else must not be treated as a single-vertex conclusion.
    bool has_unique_lineage() const { return lineage_status == 'U' && !lineage_path.empty(); }
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
