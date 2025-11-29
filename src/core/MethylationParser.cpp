#include "core/MethylationParser.hpp"
#include <sstream>
#include <algorithm>

namespace InterSubMod {

// Helper function to convert little-endian bytes to uint32
static inline uint32_t le_to_u32(const uint8_t* bytes) {
    return static_cast<uint32_t>(bytes[0]) |
           (static_cast<uint32_t>(bytes[1]) << 8) |
           (static_cast<uint32_t>(bytes[2]) << 16) |
           (static_cast<uint32_t>(bytes[3]) << 24);
}

std::vector<MethylCall> MethylationParser::parse_read(
    const bam1_t* b,
    const std::string& ref_seq,
    int32_t ref_start_pos
) {
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
    
    // First, count deltas for all modification types before "C+m?"
    std::string mm_string(mm_str);
    size_t cm_pos = mm_string.find("C+m?");
    
    if (cm_pos != std::string::npos) {
        // Count commas before "C+m?" to find ML offset
        for (size_t i = 0; i < cm_pos; i++) {
            if (mm_string[i] == ',') {
                ml_offset++;
            }
        }
        // Subtract semicolons (they don't count towards ML offset)
        for (size_t i = 0; i < cm_pos; i++) {
            if (mm_string[i] == ';') {
                ml_offset--;
            }
        }
        
        // Now parse the deltas for "C+m?"
        deltas = parse_mm_tag(mm_str, "C+m?");
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
    //   Original 'C's become 'G's in SEQ. We must look for 'G' and "G-m?" (or "G+m?").
    
    uint8_t* seq = bam_get_seq(b);
    int seq_len = b->core.l_qseq;
    bool is_reverse = bam_is_rev(b);
    
    // Determine target base and modification code based on strand
    char target_base = is_reverse ? 'G' : 'C';
    
    // For reverse strand, we might need to parse a different delta list if the tool produces "G-m?"
    // However, we already parsed "C+m?" above. 
    // TODO: Verify if the upstream tool (Dorado) produces "C+m?" for all reads (referring to original strand)
    // or if it produces "G-m?" for reverse reads.
    // Standard htslib/modkit usually normalizes. 
    // If deltas were found for "C+m?", we assume they apply to the canonical base 'C' on the original strand.
    // BUT, if the read is reverse, the 'C' on original strand is 'G' in BAM SEQ.
    // So we iterate 'G's in SEQ and apply the deltas.
    
    // If no deltas found for C+m?, check G-m? (for completeness)
    if (deltas.empty() && is_reverse) {
        deltas = parse_mm_tag(mm_str, "G-m?");
        if (deltas.empty()) deltas = parse_mm_tag(mm_str, "G+m?");
        
        // Recalculate ML offset if needed (complex, skipping for now assuming C+m? is standard)
    }

    int base_count = 0;    // Count of target bases in read sequence
    size_t delta_idx = 0;  // Current index in deltas array
    int next_target = deltas.empty() ? -1 : deltas[0];
    
    for (int seq_idx = 0; seq_idx < seq_len; seq_idx++) {
        char base = seq_nt16_str[bam_seqi(seq, seq_idx)];
        
        if (base == target_base) {
            if (base_count == next_target) {
                // This base has methylation information
                int32_t ref_pos_0based = seq_to_ref[seq_idx];
                
                if (ref_pos_0based >= 0) {  // Not an insertion
                    // Calculate offset in ref_seq
                    int ref_offset = ref_pos_0based - ref_start_pos;
                    
                    // Validate bounds and CpG context
                    if (ref_offset >= 0 && static_cast<size_t>(ref_offset) < ref_seq.size()) {
                        bool valid_cpg = false;
                        
                        if (is_reverse) {
                            // Reverse strand: We are looking at a 'G' in BAM SEQ.
                            // This corresponds to a 'C' on the reverse strand.
                            // A CpG on reverse strand is 3'-GC-5' (which is 5'-CG-3').
                            // So on the Forward Reference, we see 5'-CG-3'.
                            // The 'G' we are looking at is the G of the CpG pair.
                            // So we check if the PREVIOUS base in Ref is 'C'.
                            if (ref_seq[ref_offset] == 'G') {
                                if (ref_offset > 0 && ref_seq[ref_offset - 1] == 'C') {
                                    valid_cpg = true;
                                }
                            }
                        } else {
                            // Forward strand: We are looking at 'C'.
                            // Check if NEXT base in Ref is 'G'.
                            if (ref_seq[ref_offset] == 'C' && is_cpg_site(ref_seq, ref_offset)) {
                                valid_cpg = true;
                            }
                        }

                        if (valid_cpg) {
                            // Valid CpG site - use ml_offset to get correct probability
                            float prob = ml_data[ml_offset + delta_idx] / 255.0f;
                            
                            // For reverse strand, the CpG site coordinate is usually reported 
                            // as the 'C' position on the forward strand (standard VCF/bedMethyl convention).
                            // If we are at 'G' (pos i), the 'C' is at pos i-1.
                            // Let's standardize to reporting the 'C' position of the CpG dinucleotide.
                            int32_t report_pos = is_reverse ? (ref_pos_0based - 1) : ref_pos_0based;
                            
                            calls.emplace_back(report_pos + 1, prob);  // Convert to 1-based
                        }
                    }
                }
                
                // Move to next delta
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
    
    return calls;
}

std::vector<int> MethylationParser::parse_mm_tag(const char* mm_str, const std::string& mod_code) {
    std::vector<int> deltas;
    
    if (!mm_str) {
        return deltas;
    }
    
    std::string mm(mm_str);
    
    // Find the modification code (e.g., "C+m?")
    size_t pos = mm.find(mod_code);
    if (pos == std::string::npos) {
        return deltas;  // Modification type not found
    }
    
    // Skip past the modification code and comma
    pos += mod_code.length();
    if (pos >= mm.length() || mm[pos] != ',') {
        return deltas;  // No deltas following (shouldn't happen for valid MM)
    }
    pos++;  // Skip comma
    
    // Parse comma-separated integers until ';' or end of string
    std::string remaining = mm.substr(pos);
    std::istringstream ss(remaining);
    std::string token;
    
    while (std::getline(ss, token, ',')) {
        // Stop if we hit a semicolon (start of next modification type)
        size_t semi_pos = token.find(';');
        if (semi_pos != std::string::npos) {
            token = token.substr(0, semi_pos);
            if (!token.empty()) {
                try {
                    deltas.push_back(std::stoi(token));
                } catch (...) {
                    // Invalid number, skip
                }
            }
            break;
        }
        
        if (!token.empty()) {
            try {
                deltas.push_back(std::stoi(token));
            } catch (...) {
                // Invalid number, stop parsing
                break;
            }
        }
    }
    
    return deltas;
}

std::vector<int32_t> MethylationParser::build_seq_to_ref_map(const bam1_t* b) {
    std::vector<int32_t> seq_to_ref(b->core.l_qseq, -1);
    
    int32_t ref_pos = b->core.pos;  // Current reference position (0-based)
    int seq_pos = 0;                 // Current sequence position
    
    uint32_t* cigar = bam_get_cigar(b);
    
    for (uint32_t i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        
        switch (op) {
            case BAM_CMATCH:    // M
            case BAM_CEQUAL:    // =
            case BAM_CDIFF:     // X
                // Match/mismatch: both ref and seq advance
                for (int j = 0; j < len; j++) {
                    if (seq_pos < b->core.l_qseq) {
                        seq_to_ref[seq_pos] = ref_pos;
                    }
                    seq_pos++;
                    ref_pos++;
                }
                break;
                
            case BAM_CINS:      // I
            case BAM_CSOFT_CLIP: // S
                // Insertion/soft clip: only seq advances, ref positions remain -1
                seq_pos += len;
                break;
                
            case BAM_CDEL:      // D
            case BAM_CREF_SKIP:  // N
                // Deletion/skip: only ref advances
                ref_pos += len;
                break;
                
            case BAM_CHARD_CLIP: // H
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

} // namespace InterSubMod

