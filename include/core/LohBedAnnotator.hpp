#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace InterSubMod {

/// Classification of LOH evidence source
enum class LohSource : uint8_t {
    NONE,        ///< No LOH evidence
    BED_ONLY,    ///< LOH from external BED file only
    RATIO_ONLY,  ///< LOH from ISM hp_ratio only
    BOTH         ///< LOH from both BED and hp_ratio
};

/// Convert LohSource to string
inline std::string loh_source_to_string(LohSource src) {
    switch (src) {
        case LohSource::NONE: return "none";
        case LohSource::BED_ONLY: return "bed_only";
        case LohSource::RATIO_ONLY: return "ratio_only";
        case LohSource::BOTH: return "both";
    }
    return "none";
}

/**
 * @brief Loads and queries LOH regions from a BED file.
 *
 * BED format: chrom start end [annotation]
 * Coordinates are 0-based, half-open [start, end).
 * ISM uses 1-based positions; query methods handle conversion.
 *
 * Regions are sorted per-chromosome for binary search lookup.
 */
class LohBedAnnotator {
public:
    LohBedAnnotator() = default;

    /**
     * @brief Load LOH regions from a BED file.
     * @param bed_path Path to BED file
     * @return Number of regions loaded, or -1 on error
     */
    int load(const std::string& bed_path);

    /**
     * @brief Check if a 1-based position overlaps any LOH region.
     * @param chrom Chromosome name (e.g., "chr19")
     * @param pos 1-based genomic position
     * @return true if position falls within a LOH region
     */
    bool overlaps(const std::string& chrom, int32_t pos) const;

    /**
     * @brief Get the annotation string for the overlapping LOH region.
     * @return Annotation or empty string if no overlap
     */
    std::string get_annotation(const std::string& chrom, int32_t pos) const;

    /**
     * @brief Classify LOH source combining BED overlap and ISM hp_ratio.
     * @param chrom Chromosome name
     * @param pos 1-based position
     * @param hp_ratio_loh Whether ISM hp_ratio indicates LOH
     * @return LohSource classification
     */
    LohSource classify(const std::string& chrom, int32_t pos, bool hp_ratio_loh) const;

    /// Number of loaded regions
    int num_regions() const { return total_regions_; }

    /// Whether any regions are loaded
    bool loaded() const { return total_regions_ > 0; }

private:
    struct BedRegion {
        int32_t start;       ///< 0-based start
        int32_t end;         ///< 0-based end (exclusive)
        std::string annotation;
    };

    std::unordered_map<std::string, std::vector<BedRegion>> regions_;
    int total_regions_ = 0;

    /// Binary search for overlap (assumes sorted regions)
    const BedRegion* find_overlap(const std::string& chrom, int32_t pos_0based) const;
};

}  // namespace InterSubMod
