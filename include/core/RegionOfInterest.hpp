#pragma once

#include <cstdint>
#include <vector>

namespace InterSubMod {

/**
 * @brief Defines an analysis window centered around a Somatic SNV.
 * This is the unit of parallelization for the analysis.
 */
struct RegionOfInterest {
    int region_id;          ///< Unique Region ID
    int snv_id;             ///< ID of the anchor Somatic SNV
    int chr_id;             ///< Chromosome ID
    int32_t win_start_pos;  ///< Window Start (1-based, inclusive)
    int32_t win_end_pos;    ///< Window End (1-based, inclusive)
};

} // namespace InterSubMod
