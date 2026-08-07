#include "core/ReadParser.hpp"

#include <htslib/sam.h>

#include <limits>

#include "utils/Logger.hpp"

namespace InterSubMod {

ReadParser::ReadParser(const ReadFilterConfig& config) : config_(config) {
}

Strand ReadParser::determine_strand(const bam1_t* b) {
    // Check BAM_FREVERSE flag (0x10) to determine strand
    // If set, read is mapped to the reverse strand
    if (b->core.flag & BAM_FREVERSE) {
        return Strand::REVERSE;
    }
    return Strand::FORWARD;
}

bool ReadParser::should_keep(const bam1_t* b) const {
    auto [keep, reason] = should_keep_with_reason(b);
    return keep;
}

std::pair<bool, FilterReason> ReadParser::should_keep_with_reason(const bam1_t* b) const {
    FilterReason reasons = FilterReason::NONE;

    // Check FLAG - filter out unwanted reads
    uint16_t flag = b->core.flag;
    if (flag & BAM_FSECONDARY) {
        reasons |= FilterReason::FLAG_SECONDARY;
    }
    if (flag & BAM_FSUPPLEMENTARY) {
        reasons |= FilterReason::FLAG_SUPPLEMENTARY;
    }
    if (flag & BAM_FDUP) {
        reasons |= FilterReason::FLAG_DUPLICATE;
    }
    if (flag & BAM_FUNMAP) {
        reasons |= FilterReason::FLAG_UNMAPPED;
    }

    // Check MAPQ
    if (b->core.qual < config_.min_mapq) {
        reasons |= FilterReason::LOW_MAPQ;
    }

    // Check read length
    int read_len = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
    if (read_len < config_.min_read_length) {
        reasons |= FilterReason::SHORT_READ;
    }

    // Check for MM/ML tags if required
    if (config_.require_mm_ml) {
        uint8_t* mm_aux = bam_aux_get(b, "MM");
        uint8_t* ml_aux = bam_aux_get(b, "ML");
        if (!mm_aux) {
            reasons |= FilterReason::MISSING_MM_TAG;
        }
        if (!ml_aux) {
            reasons |= FilterReason::MISSING_ML_TAG;
        }
    }

    bool keep = (reasons == FilterReason::NONE);
    return {keep, reasons};
}

FilteredReadInfo ReadParser::create_filtered_info(const bam1_t* b, bool is_tumor, FilterReason reasons) const {
    FilteredReadInfo info;
    info.read_name = bam_get_qname(b);

    // Handle paired-end suffixes
    if (b->core.flag & BAM_FPAIRED) {
        if (b->core.flag & BAM_FREAD1) {
            info.read_name += "/1";
        } else if (b->core.flag & BAM_FREAD2) {
            info.read_name += "/2";
        }
    }

    info.chr_id = b->core.tid;
    info.align_start = b->core.pos;
    info.align_end = bam_endpos(b);
    info.mapq = b->core.qual;
    info.strand = determine_strand(b);
    info.is_tumor = is_tumor;
    info.reasons = reasons;

    return info;
}

ReadInfo ReadParser::parse(const bam1_t* b, int read_id, bool is_tumor, const SomaticSnv& anchor_snv,
                           const std::string& ref_seq, int32_t ref_start_pos) const {
    ReadInfo info;

    // Basic information
    info.read_id = read_id;
    info.read_name = bam_get_qname(b);

    // Handle paired-end suffixes
    if (b->core.flag & BAM_FPAIRED) {
        if (b->core.flag & BAM_FREAD1) {
            info.read_name += "/1";
        } else if (b->core.flag & BAM_FREAD2) {
            info.read_name += "/2";
        }
    }

    info.chr_id = anchor_snv.chr_id;
    info.align_start = b->core.pos;  // 0-based
    info.align_end = bam_endpos(b);  // 0-based, exclusive
    info.mapq = b->core.qual;
    info.is_tumor = is_tumor;

    // Determine strand orientation from BAM FLAG
    info.strand = determine_strand(b);

    // Extract HP tag (haplotype). Parse the raw value first; apply any germline-only
    // demotion downstream, but always keep the raw value in hp_tag_raw for audit.
    std::string hp_raw = "0";  // Default: unknown/unphased
    uint8_t* hp_aux = bam_aux_get(b, "HP");
    if (hp_aux) {
        char type = hp_aux[0];
        if (type == 'Z' || type == 'H') {
            // String format (longphase-s): HP:Z:1, HP:Z:2, HP:Z:1-1, HP:Z:2-1, HP:Z:3
            hp_raw = bam_aux2Z(hp_aux);
        } else if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I') {
            // Integer format (longphase-to): 1=HP1, 2=HP2, 11=HP1-1, 21=HP2-1, 33=HP3
            // Map to canonical string format used throughout LabelTest
            int hp_int = bam_aux2i(hp_aux);
            switch (hp_int) {
                case 1:  hp_raw = "1";   break;
                case 2:  hp_raw = "2";   break;
                case 11: hp_raw = "1-1"; break;
                case 21: hp_raw = "2-1"; break;
                case 33: hp_raw = "3";   break;
                default:
                    // Unrecognized HP integer: kept verbatim but it maps to HP_OTHER downstream and is
                    // silently dropped from all HP grouping. Surface it so unexpected phasing tags are visible.
                    hp_raw = std::to_string(hp_int);
                    LOG_DEBUG("ReadParser: unrecognized HP integer tag " + std::to_string(hp_int) +
                              " (expected 1/2/11/21/33); read excluded from HP grouping");
                    break;
            }
        }
    }
    info.hp_tag_raw = hp_raw;
    // Self-phasing fallback: treat somatic HP tags as unphased.
    // LongPhase-TO marks somatic-only phase blocks with HP:i:11/21/33; these reflect
    // the somatic variant phasing itself (circular dependency) rather than germline
    // haplotypes. When this flag is on we demote them to "0" so downstream HP-based
    // features only trust the germline phasing (HP:i:1/2).
    if (config_.germline_hp_only && (hp_raw == "1-1" || hp_raw == "2-1" || hp_raw == "3")) {
        info.hp_tag = "0";
    } else {
        info.hp_tag = hp_raw;
    }

    // Extract PS (phase set) and lineage tags written by pipeline/lineage/bin/tag_bam.
    // These are optional: BAMs that never went through tag_bam simply leave them empty,
    // so every downstream consumer must tolerate absence rather than assume presence.
    extract_lineage_tags(b, info);

    // Determine ALT support
    info.alt_support = determine_alt_support(b, anchor_snv, ref_seq, ref_start_pos);

    return info;
}

void ReadParser::extract_lineage_tags(const bam1_t* b, ReadInfo& info) {
    info.phase_set = -1;
    info.lineage_status = '\0';

    if (uint8_t* ps_aux = bam_aux_get(b, "PS"); ps_aux != nullptr) {
        const char type = ps_aux[0];
        if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I') {
            info.phase_set = bam_aux2i(ps_aux);
        }
    }

    // lc / lu / lv / lp / lo are all Z-type strings
    const std::pair<const char*, std::string*> string_tags[] = {
        {"lc", &info.lineage_component},
        {"lu", &info.lineage_block},
        {"lv", &info.lineage_path},
        {"lp", &info.lineage_pattern},
        {"lo", &info.mutation_order},
    };
    for (const auto& [tag, dest] : string_tags) {
        uint8_t* aux = bam_aux_get(b, tag);
        if (aux != nullptr && (aux[0] == 'Z' || aux[0] == 'H')) {
            const char* v = bam_aux2Z(aux);
            if (v != nullptr) *dest = v;
        }
    }

    // ls is an A-type single character
    if (uint8_t* ls_aux = bam_aux_get(b, "ls"); ls_aux != nullptr && ls_aux[0] == 'A') {
        info.lineage_status = bam_aux2A(ls_aux);
    }

    // Invariant from HIERARCHICAL_TAG_SPEC: a lineage path is only meaningful when its
    // confidence status travelled with it. Drop an orphaned path rather than let a
    // downstream axis silently treat it as a confirmed assignment.
    if (!info.lineage_path.empty() && info.lineage_status == '\0') {
        LOG_DEBUG("read " + info.read_name + " carries lv without ls; dropping lv");
        info.lineage_path.clear();
        info.mutation_order.clear();
    }
}

AltSupport ReadParser::determine_alt_support(const bam1_t* b, const SomaticSnv& snv, const std::string& ref_seq,
                                             int32_t ref_start_pos) const {
    return determine_alt_support_with_reason(b, snv, ref_seq, ref_start_pos).support;
}

AltSupportResult ReadParser::determine_alt_support_with_reason(const bam1_t* b, const SomaticSnv& snv,
                                                               const std::string& ref_seq [[maybe_unused]],
                                                               int32_t ref_start_pos [[maybe_unused]]) const {
    // SNV position (convert 1-based to 0-based).
    // Guard against pos==0 (invalid 1-based) and values that would overflow int32_t.
    if (snv.pos == 0 || snv.pos > static_cast<uint32_t>(std::numeric_limits<int32_t>::max())) {
        return AltSupportResult(AltSupport::UNKNOWN, FilterReason::SNV_NOT_COVERED);
    }
    int32_t snv_pos_0based = static_cast<int32_t>(snv.pos) - 1;

    // Check if read covers the SNV position
    int32_t read_start = b->core.pos;
    int32_t read_end = bam_endpos(b);

    if (snv_pos_0based < read_start || snv_pos_0based >= read_end) {
        return AltSupportResult(AltSupport::UNKNOWN, FilterReason::SNV_NOT_COVERED);
    }

    // Traverse CIGAR to find read offset at SNV position
    int32_t ref_pos = read_start;  // Current reference position
    int seq_pos = 0;               // Current sequence position
    int read_offset = -1;          // Sequence offset at SNV position

    uint32_t* cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        switch (op) {
            case BAM_CMATCH:  // M
            case BAM_CEQUAL:  // =
            case BAM_CDIFF:   // X
                // Match/mismatch: consumes both ref and seq
                // Check if the SNV position falls within this CIGAR operation
                if (ref_pos <= snv_pos_0based && snv_pos_0based < ref_pos + len) {
                    // Found it! Calculate the offset in the read sequence
                    read_offset = seq_pos + (snv_pos_0based - ref_pos);
                    goto found_offset;
                }
                ref_pos += len;
                seq_pos += len;
                break;

            case BAM_CINS:        // I
            case BAM_CSOFT_CLIP:  // S
                // Insertion/soft clip: consumes only seq
                // These bases are in the read but not in the reference
                seq_pos += len;
                break;

            case BAM_CDEL:       // D
            case BAM_CREF_SKIP:  // N
                // Deletion/skip: consumes only ref
                // These positions are in the reference but skipped in the read
                if (ref_pos <= snv_pos_0based && snv_pos_0based < ref_pos + len) {
                    // SNV falls within a deletion - so the read does not support ALT or REF
                    return AltSupportResult(AltSupport::UNKNOWN, FilterReason::SNV_IN_DELETION);
                }
                ref_pos += len;
                break;

            case BAM_CHARD_CLIP:  // H
                // Hard clip: consumes nothing (removed from read sequence)
                break;

            default:
                // Unknown CIGAR operation
                break;
        }
    }

found_offset:
    if (read_offset == -1) {
        return AltSupportResult(AltSupport::UNKNOWN, FilterReason::SNV_NOT_COVERED);
    }

    // Check base quality at SNV position
    uint8_t* qual = bam_get_qual(b);
    if (qual[read_offset] < config_.min_base_quality) {
        return AltSupportResult(AltSupport::UNKNOWN, FilterReason::LOW_BASE_QUALITY);
    }

    // Get the base at SNV position
    uint8_t* seq = bam_get_seq(b);
    char base = seq_nt16_str[bam_seqi(seq, read_offset)];

    // Compare with REF and ALT
    if (base == snv.alt_base) {
        return AltSupportResult(AltSupport::ALT, FilterReason::NONE);
    } else if (base == snv.ref_base) {
        return AltSupportResult(AltSupport::REF, FilterReason::NONE);
    } else {
        return AltSupportResult(AltSupport::UNKNOWN, FilterReason::NOT_REF_OR_ALT);
    }
}

}  // namespace InterSubMod
