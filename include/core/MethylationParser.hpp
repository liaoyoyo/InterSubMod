#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <htslib/sam.h>

namespace InterSubMod {

/**
 * @brief Represents a single methylation call at a genomic position.
 */
struct MethylCall {
    int32_t ref_pos;      ///< 1-based genomic position (hg38 coordinates)
    float probability;    ///< Methylation probability [0.0, 1.0]
    
    MethylCall() : ref_pos(0), probability(0.0f) {}
    MethylCall(int32_t pos, float prob) : ref_pos(pos), probability(prob) {}
};

/**
 * @brief Parser for extracting methylation information from BAM MM/ML tags.
 * 
 * This class handles the complex task of:
 * 1. Parsing MM tag (delta-encoded modification positions)
 * 2. Extracting ML tag (modification probabilities)
 * 3. Mapping read sequence positions to reference genome coordinates
 * 4. Validating CpG context in the reference sequence
 * 
 * MM Tag Format: "C+m?,skip1,skip2,...;C+h?,..." 
 *   - "C+m?" indicates 5-methylcytosine modifications
 *   - Numbers are delta-encoded skip counts between modified bases
 * 
 * ML Tag Format: Array of uint8 values [0-255]
 *   - Each value corresponds to a modification in MM
 *   - Probability = ML[i] / 255.0
 * 
 * Thread-safe: This class is stateless and can be used from multiple threads.
 */
class MethylationParser {
public:
    /**
     * @brief Parses methylation information from a BAM record.
     * 
     * @param b BAM record containing MM and ML tags.
     * @param ref_seq Reference sequence covering the read's alignment.
     * @param ref_start_pos 0-based start position of ref_seq in genome.
     * @return Vector of MethylCall structures (may be empty if no valid calls).
     * 
     * @note Only returns calls that:
     *       1. Are in CpG context (verified against ref_seq)
     *       2. Map to valid reference positions (not in insertions)
     *       3. Have matching MM/ML array lengths
     */
    std::vector<MethylCall> parse_read(
        const bam1_t* b,
        const std::string& ref_seq,
        int32_t ref_start_pos
    );

private:
    /**
     * @brief Parses MM tag and extracts delta-encoded skip counts.
     * 
     * Searches for the specified modification code (e.g., "C+m?") and
     * extracts the comma-separated skip counts that follow it.
     * 
     * @param mm_str MM tag string (e.g., "C+m?,3,5,0,2;").
     * @param mod_code Modification code to search for (default "C+m?" for 5mC).
     * @return Vector of delta skip counts, empty if mod_code not found.
     * 
     * @note Delta encoding means: skip_count[i] is the number of unmodified
     *       bases of the target type (e.g., 'C') between modifications.
     */
    std::vector<int> parse_mm_tag(const char* mm_str, const std::string& mod_code = "C+m?");
    
    /**
     * @brief Builds a mapping from read sequence index to reference position.
     * 
     * Traverses the CIGAR string to determine which read bases align to
     * which reference positions. Insertions map to -1 (not in reference).
     * 
     * @param b BAM record.
     * @return Vector where seq_to_ref[i] = reference position for seq[i],
     *         or -1 if seq[i] is an insertion.
     */
    std::vector<int32_t> build_seq_to_ref_map(const bam1_t* b);
    
    /**
     * @brief Checks if a position in the reference sequence is a CpG site.
     * 
     * A valid CpG site has 'C' at offset and 'G' at offset+1.
     * 
     * @param ref_seq Reference sequence (uppercase).
     * @param offset 0-based offset within ref_seq.
     * @return true if position is a CpG dinucleotide.
     */
    bool is_cpg_site(const std::string& ref_seq, size_t offset);
};

} // namespace InterSubMod

