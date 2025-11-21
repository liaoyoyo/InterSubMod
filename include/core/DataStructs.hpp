#pragma once

#include <string>
#include <cstdint>
#include "Types.hpp"

namespace InterSubMod {

struct CpGSite {
    int cpg_id;
    int chr_id;
    int32_t pos;
    bool in_pmd;
    bool in_repressive_state;
    bool accessible;
};

struct ReadInfo {
    int read_id;
    std::string read_name;
    int chr_id;
    int32_t align_start;
    int32_t align_end;
    int mapq;
    int hp_tag; // 0=Unknown, 1=H1, 2=H2
    bool is_tumor;
    AltSupport alt_support;
};

} // namespace InterSubMod

