#include "core/LohBedAnnotator.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace InterSubMod {

int LohBedAnnotator::load(const std::string& bed_path) {
    std::ifstream ifs(bed_path);
    if (!ifs.is_open()) {
        std::cerr << "[LohBedAnnotator] Failed to open BED file: " << bed_path << "\n";
        return -1;
    }

    regions_.clear();
    total_regions_ = 0;

    std::string line;
    int line_num = 0;
    while (std::getline(ifs, line)) {
        ++line_num;
        // Skip empty lines and comment lines
        if (line.empty() || line[0] == '#') continue;
        // Skip browser/track lines
        if (line.compare(0, 7, "browser") == 0 || line.compare(0, 5, "track") == 0) continue;

        std::istringstream ss(line);
        std::string chrom;
        int32_t start, end;

        if (!(ss >> chrom >> start >> end)) {
            std::cerr << "[LohBedAnnotator] Skipping malformed line " << line_num << ": " << line << "\n";
            continue;
        }

        // Validate coordinates
        if (start < 0 || end < 0 || start >= end) {
            std::cerr << "[LohBedAnnotator] Skipping invalid coordinates at line " << line_num
                      << ": start=" << start << " end=" << end << "\n";
            continue;
        }

        BedRegion region;
        region.start = start;
        region.end = end;

        // Read optional annotation (rest of line)
        std::string annotation;
        if (std::getline(ss, annotation)) {
            // Trim leading whitespace/tab
            size_t pos = annotation.find_first_not_of(" \t");
            if (pos != std::string::npos) {
                region.annotation = annotation.substr(pos);
            }
        }

        regions_[chrom].push_back(region);
        ++total_regions_;
    }

    // Sort regions per chromosome by start position for binary search
    for (auto& [chrom, vec] : regions_) {
        std::sort(vec.begin(), vec.end(), [](const BedRegion& a, const BedRegion& b) {
            return a.start < b.start || (a.start == b.start && a.end < b.end);
        });
    }

    return total_regions_;
}

bool LohBedAnnotator::overlaps(const std::string& chrom, int32_t pos) const {
    // Convert 1-based ISM position to 0-based BED coordinate
    int32_t pos_0based = pos - 1;
    return find_overlap(chrom, pos_0based) != nullptr;
}

std::string LohBedAnnotator::get_annotation(const std::string& chrom, int32_t pos) const {
    int32_t pos_0based = pos - 1;
    const BedRegion* region = find_overlap(chrom, pos_0based);
    if (region) {
        return region->annotation;
    }
    return "";
}

LohSource LohBedAnnotator::classify(const std::string& chrom, int32_t pos, bool hp_ratio_loh) const {
    bool bed_overlap = overlaps(chrom, pos);

    if (bed_overlap && hp_ratio_loh) return LohSource::BOTH;
    if (bed_overlap) return LohSource::BED_ONLY;
    if (hp_ratio_loh) return LohSource::RATIO_ONLY;
    return LohSource::NONE;
}

const LohBedAnnotator::BedRegion* LohBedAnnotator::find_overlap(const std::string& chrom,
                                                                  int32_t pos_0based) const {
    auto it = regions_.find(chrom);
    if (it == regions_.end()) return nullptr;

    const auto& vec = it->second;
    if (vec.empty()) return nullptr;

    // Binary search: find last region where start <= pos_0based
    // Then check if pos_0based < end (i.e., pos is within [start, end))
    auto upper = std::upper_bound(vec.begin(), vec.end(), pos_0based,
                                  [](int32_t p, const BedRegion& r) { return p < r.start; });

    // upper points to first region with start > pos_0based
    // We need to check the region before it (if any)
    if (upper == vec.begin()) return nullptr;

    --upper;
    // Check if pos_0based falls within [start, end)
    if (pos_0based >= upper->start && pos_0based < upper->end) {
        return &(*upper);
    }

    return nullptr;
}

}  // namespace InterSubMod
