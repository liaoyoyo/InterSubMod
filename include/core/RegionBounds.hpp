#pragma once

/**
 * @file RegionBounds.hpp
 * @brief Pure data struct for genomic region coordinates.
 *
 * Header-only. No HTSlib or heavy dependencies.
 */

#include <algorithm>
#include <cstdint>
#include <string>

namespace InterSubMod {

/**
 * @brief Genomic coordinate window for a single analysis region.
 *
 * Coordinates are 1-based, inclusive on both ends, matching the SNV
 * position convention used throughout this codebase.
 */
struct RegionBounds {
    std::string chr;
    int32_t start = 0;  ///< 1-based, inclusive
    int32_t end   = 0;  ///< 1-based, inclusive

    /** True if both endpoints are positive and end >= start (allows single-base regions). */
    bool is_valid() const noexcept { return start >= 1 && end >= start; }
    /** Number of bases covered (inclusive both ends: [start, end]). */
    int32_t width() const noexcept { return end - start + 1; }
};

/**
 * @brief Compute the initial analysis window centered on a SNV position.
 *
 * @param chr_name    Chromosome name (e.g. "chr17")
 * @param snv_pos     SNV 1-based position (must be > 0)
 * @param window_size Half-width in bp; final window = [pos-w, pos+w]
 * @param chr_length  Chromosome length in bp; 0 means unknown (no end-clamping)
 */
inline RegionBounds compute_region_bounds(const std::string& chr_name,
                                           int32_t            snv_pos,
                                           int32_t            window_size,
                                           int32_t            chr_length) noexcept {
    RegionBounds b;
    b.chr   = chr_name;
    b.start = std::max(1, snv_pos - window_size);
    b.end   = (chr_length > 0) ? std::min(chr_length, snv_pos + window_size)
                                : (snv_pos + window_size);
    return b;
}

}  // namespace InterSubMod
