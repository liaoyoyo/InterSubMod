#pragma once

#include <cstdint>
#include <string>

namespace InterSubMod {

enum class DistanceMetricType {
    NHD,
    L1,
    L2,
    CORR,
    JACCARD
};

enum class NanDistanceStrategy {
    MAX_DIST,
    SKIP
};

enum class AltSupport {
    ALT,
    REF,
    UNKNOWN
};

/**
 * @brief Strand orientation of a read.
 * 
 * Determined by BAM FLAG bit 0x10 (BAM_FREVERSE):
 * - FORWARD: Read maps to the forward/positive strand
 * - REVERSE: Read maps to the reverse/negative strand
 * - UNKNOWN: Strand cannot be determined (should be rare)
 */
enum class Strand : uint8_t {
    FORWARD = 0,  ///< Forward strand (+)
    REVERSE = 1,  ///< Reverse strand (-)
    UNKNOWN = 2   ///< Unknown strand
};

/**
 * @brief Log level for controlling output verbosity.
 */
enum class LogLevel {
    LOG_ERROR = 0,    ///< Only errors
    LOG_WARN = 1,     ///< Errors and warnings
    LOG_INFO = 2,     ///< Normal operational messages
    LOG_DEBUG = 3     ///< Detailed debug output including filtered reads
};

/**
 * @brief Reasons why a read may be filtered out.
 * 
 * Used in debug mode to track and report why reads were excluded.
 */
enum class FilterReason : uint16_t {
    NONE = 0,
    FLAG_SECONDARY = 1 << 0,      ///< Secondary alignment
    FLAG_SUPPLEMENTARY = 1 << 1,  ///< Supplementary alignment
    FLAG_DUPLICATE = 1 << 2,      ///< PCR/optical duplicate
    FLAG_UNMAPPED = 1 << 3,       ///< Unmapped read
    LOW_MAPQ = 1 << 4,            ///< MAPQ below threshold
    SHORT_READ = 1 << 5,          ///< Read length below threshold
    MISSING_MM_TAG = 1 << 6,      ///< Missing MM tag
    MISSING_ML_TAG = 1 << 7,      ///< Missing ML tag
    SNV_NOT_COVERED = 1 << 8,     ///< Read does not cover SNV position
    SNV_IN_DELETION = 1 << 9,     ///< SNV falls in deletion region
    LOW_BASE_QUALITY = 1 << 10,   ///< Base quality at SNV below threshold
    NOT_REF_OR_ALT = 1 << 11      ///< Base at SNV is neither REF nor ALT
};

/**
 * @brief Enable bitwise operations for FilterReason.
 */
inline FilterReason operator|(FilterReason a, FilterReason b) {
    return static_cast<FilterReason>(static_cast<uint16_t>(a) | static_cast<uint16_t>(b));
}

inline FilterReason operator&(FilterReason a, FilterReason b) {
    return static_cast<FilterReason>(static_cast<uint16_t>(a) & static_cast<uint16_t>(b));
}

inline FilterReason& operator|=(FilterReason& a, FilterReason b) {
    a = a | b;
    return a;
}

inline bool has_flag(FilterReason flags, FilterReason check) {
    return (static_cast<uint16_t>(flags) & static_cast<uint16_t>(check)) != 0;
}

/**
 * @brief Convert FilterReason to human-readable string.
 */
inline std::string filter_reason_to_string(FilterReason reason) {
    if (reason == FilterReason::NONE) return "NONE";
    
    std::string result;
    if (has_flag(reason, FilterReason::FLAG_SECONDARY)) result += "SECONDARY,";
    if (has_flag(reason, FilterReason::FLAG_SUPPLEMENTARY)) result += "SUPPLEMENTARY,";
    if (has_flag(reason, FilterReason::FLAG_DUPLICATE)) result += "DUPLICATE,";
    if (has_flag(reason, FilterReason::FLAG_UNMAPPED)) result += "UNMAPPED,";
    if (has_flag(reason, FilterReason::LOW_MAPQ)) result += "LOW_MAPQ,";
    if (has_flag(reason, FilterReason::SHORT_READ)) result += "SHORT_READ,";
    if (has_flag(reason, FilterReason::MISSING_MM_TAG)) result += "MISSING_MM,";
    if (has_flag(reason, FilterReason::MISSING_ML_TAG)) result += "MISSING_ML,";
    if (has_flag(reason, FilterReason::SNV_NOT_COVERED)) result += "SNV_NOT_COVERED,";
    if (has_flag(reason, FilterReason::SNV_IN_DELETION)) result += "SNV_IN_DELETION,";
    if (has_flag(reason, FilterReason::LOW_BASE_QUALITY)) result += "LOW_BASE_QUALITY,";
    if (has_flag(reason, FilterReason::NOT_REF_OR_ALT)) result += "NOT_REF_OR_ALT,";
    
    // Remove trailing comma
    if (!result.empty() && result.back() == ',') {
        result.pop_back();
    }
    return result;
}

} // namespace InterSubMod
