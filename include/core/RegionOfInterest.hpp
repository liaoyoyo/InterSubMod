#pragma once

#include <cstdint>
#include <vector>

namespace InterSubMod {

struct RegionOfInterest {
    int region_id;
    int snv_id; // anchor SNV
    int chr_id;
    int32_t win_start_pos;
    int32_t win_end_pos;
};

} // namespace InterSubMod

