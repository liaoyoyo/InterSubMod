#include "core/ReadParser.hpp"

#include <htslib/sam.h>

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

    // Extract HP tag (haplotype)
    info.hp_tag = "0";  // Default: unknown
    uint8_t* hp_aux = bam_aux_get(b, "HP");
    if (hp_aux) {
        char type = hp_aux[0];
        if (type == 'Z' || type == 'H') {
            info.hp_tag = bam_aux2Z(hp_aux);
        } else if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I') {
            info.hp_tag = std::to_string(bam_aux2i(hp_aux));
        }
    }

    // Determine ALT support
    info.alt_support = determine_alt_support(b, anchor_snv, ref_seq, ref_start_pos);

    return info;
}

AltSupport ReadParser::determine_alt_support(const bam1_t* b, const SomaticSnv& snv, const std::string& ref_seq,
                                             int32_t ref_start_pos) const {
    return determine_alt_support_with_reason(b, snv, ref_seq, ref_start_pos).support;
}

AltSupportResult ReadParser::determine_alt_support_with_reason(const bam1_t* b, const SomaticSnv& snv,
                                                               const std::string& ref_seq [[maybe_unused]],
                                                               int32_t ref_start_pos [[maybe_unused]]) const {
    // SNV position (convert 1-based to 0-based)
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
