#include "core/MethylationParser.hpp"

#include <algorithm>
#include <charconv>
#include <cstring>

namespace InterSubMod {

// Helper function to convert little-endian bytes to uint32
static inline uint32_t le_to_u32(const uint8_t* bytes) {
    return static_cast<uint32_t>(bytes[0]) | (static_cast<uint32_t>(bytes[1]) << 8) |
           (static_cast<uint32_t>(bytes[2]) << 16) | (static_cast<uint32_t>(bytes[3]) << 24);
}

std::vector<MethylCall> MethylationParser::parse_read(const bam1_t* b, const std::string& ref_seq,
                                                      int32_t ref_start_pos) {
    std::vector<MethylCall> calls;

    // Retrieve MM and ML tags
    uint8_t* mm_aux = bam_aux_get(b, "MM");
    uint8_t* ml_aux = bam_aux_get(b, "ML");

    if (!mm_aux || !ml_aux) {
        return calls;  // No methylation tags
    }

    const char* mm_str = bam_aux2Z(mm_aux);

    // Parse ML array
    // Format: 'B' (array type), 'C' (uint8), [4 bytes length], [data...]
    if (ml_aux[0] != 'B' || ml_aux[1] != 'C') {
        return calls;  // Invalid ML format
    }

    uint32_t ml_len = le_to_u32(ml_aux + 2);
    const uint8_t* ml_data = ml_aux + 6;

    // Parse MM tag to get delta-encoded positions and ML offset
    // MM format: "C+h?,<deltas>;C+m?,<deltas>;"
    // Each modification type has its own delta list
    // ML array contains all probabilities in order: [C+h? probs..., C+m? probs...]

    int ml_offset = 0;
    std::vector<int> deltas;

    // Pointer-based block scan: avoids std::string copy and istringstream.
    // Each ';'-delimited block is processed as a [block_start, block_end) range.
    bool found_target = false;
    const char* p     = mm_str;

    while (*p) {
        const char* block_start = p;
        while (*p && *p != ';') ++p;
        const char* block_end = p;
        if (*p == ';') ++p;

        if (block_start == block_end) continue;

        // Check if this block is "C+m?" (target modification)
        const char* target = "C+m?";
        const size_t tlen  = 4;
        bool is_target = (static_cast<size_t>(block_end - block_start) >= tlen) &&
                         std::memcmp(block_start, target, tlen) == 0;

        if (is_target) {
            deltas      = parse_mm_tag(mm_str, "C+m?");
            found_target = true;
            break;
        }

        // Not our target: count deltas (commas after the first ',') to advance ml_offset.
        // Block format: "X+y?,d1,d2,...", so delta_count = comma_count_in_delta_part + 1.
        const char* comma =
            static_cast<const char*>(std::memchr(block_start, ',', static_cast<size_t>(block_end - block_start)));
        if (comma) {
            int delta_count = 0;
            for (const char* q = comma + 1; q < block_end; ++q) {
                if (*q == ',') ++delta_count;
            }
            delta_count++;  // +1 for last value (no trailing comma)
            ml_offset += delta_count;
        }
    }

    if (!found_target) {
        return calls;  // C+m? not found
    }

    if (deltas.empty()) {
        return calls;  // No 5mC modifications found
    }

    // Verify we have enough ML data
    if (static_cast<size_t>(ml_offset + deltas.size()) > ml_len) {
        return calls;  // Not enough ML data
    }

    // Build sequence-to-reference mapping
    // This maps each base in the read sequence to its genomic coordinate
    std::vector<int32_t> seq_to_ref = build_seq_to_ref_map(b);

    // Iterate through read sequence to find modified bases
    // Note: MM tag refers to the SEQ in the BAM file.
    // - If read is mapped to Forward strand: SEQ matches Ref. We look for 'C' and "C+m?".
    // - If read is mapped to Reverse strand: SEQ is Reverse Complemented.
    //   Original 'C's become 'G's in SEQ. We must look for 'G'.
    //   CRITICAL: MM tags list deltas in 5'->3' order of the ORIGINAL read.
    //   BAM SEQ for reverse strand is 3'->5' of original read (RevComp).
    //   Therefore, for reverse strand, we must iterate BAM SEQ BACKWARDS (Len-1 -> 0)
    //   to match the MM tag order.

    uint8_t* seq = bam_get_seq(b);
    int seq_len = b->core.l_qseq;
    bool is_reverse = bam_is_rev(b);

    int base_count = 0;    // Count of target bases found
    size_t delta_idx = 0;  // Current index in deltas array
    int next_target = deltas.empty() ? -1 : deltas[0];

    if (is_reverse) {
        // REVERSE STRAND LOGIC
        // Iterate BACKWARDS (seq_len-1 -> 0) to match 5'->3' of original read
        // Target 'G' (which is 'C' in original read)
        char target_base = 'G';

        for (int seq_idx = seq_len - 1; seq_idx >= 0; seq_idx--) {
            char base = seq_nt16_str[bam_seqi(seq, seq_idx)];

            if (base == target_base) {
                if (base_count == next_target) {
                    // Match!
                    int32_t ref_pos_0based = seq_to_ref[seq_idx];

                    if (ref_pos_0based >= 0) {
                        int ref_offset = ref_pos_0based - ref_start_pos;

                        // Validate: Ref should be G, Pred should be C (CpG)
                        if (ref_offset > 0 && static_cast<size_t>(ref_offset) < ref_seq.size() &&
                            ref_seq[ref_offset] == 'G' && ref_seq[ref_offset - 1] == 'C') {
                            float prob = ml_data[ml_offset + delta_idx] / 255.0f;
                            // Report at C position (ref_pos - 1)
                            calls.emplace_back(ref_pos_0based,
                                               prob);  // ref_pos_0based - 1 + 1 (1-based) = ref_pos_0based
                        }
                    }

                    // Next delta
                    delta_idx++;
                    if (delta_idx < deltas.size()) {
                        next_target += deltas[delta_idx] + 1;
                    } else {
                        next_target = -1;
                    }
                }
                base_count++;
            }
        }
    } else {
        // FORWARD STRAND LOGIC
        // Iterate FORWARDS (0 -> seq_len-1)
        // Target 'C'
        char target_base = 'C';

        for (int seq_idx = 0; seq_idx < seq_len; seq_idx++) {
            char base = seq_nt16_str[bam_seqi(seq, seq_idx)];

            if (base == target_base) {
                if (base_count == next_target) {
                    // Match!
                    int32_t ref_pos_0based = seq_to_ref[seq_idx];

                    if (ref_pos_0based >= 0) {
                        int ref_offset = ref_pos_0based - ref_start_pos;

                        // Validate: Ref should be C, Next should be G (CpG)
                        if (ref_offset >= 0 && static_cast<size_t>(ref_offset + 1) < ref_seq.size() &&
                            ref_seq[ref_offset] == 'C' &&
                            ref_seq[ref_offset + 1] == 'G') {  // Simple is_cpg_site inline

                            float prob = ml_data[ml_offset + delta_idx] / 255.0f;
                            calls.emplace_back(ref_pos_0based + 1, prob);
                        }
                    }

                    // Next delta
                    delta_idx++;
                    if (delta_idx < deltas.size()) {
                        next_target += deltas[delta_idx] + 1;
                    } else {
                        next_target = -1;
                    }
                }
                base_count++;
            }
        }
    }

    return calls;
}

std::vector<int> MethylationParser::parse_mm_tag(const char* mm_str, const std::string& mod_code) {
    std::vector<int> deltas;

    if (!mm_str) return deltas;

    // Locate modification code without creating std::string — avoids copy allocation
    const char* found = std::strstr(mm_str, mod_code.c_str());
    if (!found) return deltas;

    // Skip past modification code and mandatory leading comma
    const char* p = found + mod_code.length();
    if (*p != ',') return deltas;
    ++p;

    // Parse comma-separated integers until ';' or end of string.
    // std::from_chars (C++17): no allocation, no exceptions, stops at first non-digit.
    deltas.reserve(64);
    const char* str_end = mm_str + std::strlen(mm_str);

    while (p < str_end && *p != ';') {
        int val;
        auto [next, ec] = std::from_chars(p, str_end, val);
        if (ec != std::errc{}) break;
        deltas.push_back(val);
        p = next;
        if (p < str_end && *p == ',') ++p;
    }

    return deltas;
}

std::vector<int32_t> MethylationParser::build_seq_to_ref_map(const bam1_t* b) {
    std::vector<int32_t> seq_to_ref(b->core.l_qseq, -1);

    int32_t ref_pos = b->core.pos;  // Current reference position (0-based)
    int seq_pos = 0;                // Current sequence position

    uint32_t* cigar = bam_get_cigar(b);

    for (uint32_t i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        switch (op) {
            case BAM_CMATCH:  // M
            case BAM_CEQUAL:  // =
            case BAM_CDIFF:   // X
                // Match/mismatch: both ref and seq advance
                for (int j = 0; j < len; j++) {
                    if (seq_pos < b->core.l_qseq) {
                        seq_to_ref[seq_pos] = ref_pos;
                    }
                    seq_pos++;
                    ref_pos++;
                }
                break;

            case BAM_CINS:        // I
            case BAM_CSOFT_CLIP:  // S
                // Insertion/soft clip: only seq advances, ref positions remain -1
                seq_pos += len;
                break;

            case BAM_CDEL:       // D
            case BAM_CREF_SKIP:  // N
                // Deletion/skip: only ref advances
                ref_pos += len;
                break;

            case BAM_CHARD_CLIP:  // H
                // Hard clip: does not appear in SEQ, no action needed
                break;

            default:
                // Unknown CIGAR operation, skip
                break;
        }
    }

    return seq_to_ref;
}

bool MethylationParser::is_cpg_site(const std::string& ref_seq, size_t offset) {
    if (offset + 1 >= ref_seq.size()) {
        return false;
    }
    return (ref_seq[offset] == 'C' && ref_seq[offset + 1] == 'G');
}

}  // namespace InterSubMod
