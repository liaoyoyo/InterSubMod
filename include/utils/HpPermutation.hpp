#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "core/DataStructs.hpp"

namespace InterSubMod {

/**
 * @brief R-SELFREF permutation null: shuffle HP tags among reads within a region.
 *
 * Shuffles the (hp_tag, hp_tag_raw) pair as a unit across reads, preserving the per-region HP
 * marginal (same multiset of tags) but breaking the read<->tag self-phasing link. Used to build a
 * null distribution for HP-based discrimination metrics (does the observed same-HP structure exceed
 * what random HP assignment would produce?).
 *
 * Deterministic: the same (seed, chrom, pos) always yields the same permutation, so a run is
 * reproducible. Uses an explicit integer region hash (not std::hash, which is implementation-defined
 * across builds). No-op when seed <= 0 or the region has fewer than 2 reads.
 *
 * @param reads  Per-region reads; hp_tag/hp_tag_raw are permuted in place.
 * @param seed   Permutation seed (>0 to permute; <=0 is a no-op).
 * @param chrom  Region chromosome name (part of the per-region seed).
 * @param pos    Region anchor position (part of the per-region seed).
 */
inline void permute_hp_tags_per_region(std::vector<ReadInfo>& reads, int seed, const std::string& chrom, int pos) {
    if (reads.size() < 2 || seed <= 0) {
        return;
    }
    uint32_t region_seed = static_cast<uint32_t>(seed) * 2654435761u + static_cast<uint32_t>(pos);
    for (char c : chrom) {
        region_seed = region_seed * 31u + static_cast<unsigned char>(c);
    }
    std::mt19937 rng(region_seed);
    std::vector<std::pair<std::string, std::string>> tags;
    tags.reserve(reads.size());
    for (const auto& r : reads) {
        tags.emplace_back(r.hp_tag, r.hp_tag_raw);
    }
    std::shuffle(tags.begin(), tags.end(), rng);
    for (std::size_t i = 0; i < reads.size(); ++i) {
        reads[i].hp_tag = tags[i].first;
        reads[i].hp_tag_raw = tags[i].second;
    }
}

}  // namespace InterSubMod
