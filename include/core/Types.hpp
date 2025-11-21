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

} // namespace InterSubMod
