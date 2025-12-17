# Investigation of Distance Heatmap "No Intersection" Behavior

**Date**: 2025-12-17
**Status**: Verified / Expected Behavior

## Issue Description

User observed that in the distance heatmap (specifically `distance_heatmap.png`), many read pairs display a distance of approximately 1.0 (indicating maximum difference), even when those reads have "completely no intersection methylation" (i.e., no overlapping CpG sites).

## Investigation Findings

We analyzed the source code responsible for distance matrix calculation, specifically:

- `src/core/DistanceMatrix.cpp`
- `include/core/DistanceMatrix.hpp`
- `include/core/Config.hpp`

### 1. Distance Calculation Logic

The distance calculation functions (e.g., `calculate_nhd`, `calculate_l1`, etc.) enforce a minimum overlap requirement.

```cpp
// src/core/DistanceMatrix.cpp
if (common_count < min_cov) {
    return -1.0;
}
```

If the number of shared valid CpG sites (`common_count`) is less than `min_common_coverage` (default: 3), the function returns `-1.0`, indicating an "invalid" or "insufficiently overlapping" pair. This naturally includes pairs with **zero** intersection.

### 2. Handling Invalid Pairs

In the main matrix computation loop, these "invalid" return values are caught and replaced:

```cpp
// src/core/DistanceMatrix.cpp
if (dist < 0) {
    // Invalid pair
    dist = nan_val;
    invalid_pairs++;
}
```

### 3. Default Configuration

The value of `nan_val` determines the final distance assigned to non-overlapping pairs. This is controlled by the configuration:

- `NanDistanceStrategy` defaults to `MAX_DIST` (in `DistanceMatrix.hpp`).
- `max_distance_value` defaults to `1.0` (in `Config.hpp`).

Therefore, any read pair with insufficient overlap (including no intersection) is automatically assigned a distance of **1.0**.

## Conclusion

The observed behavior is **reasonable** and **working as designed**.

### Why is this reasonable?

1. **Clustering Stability**: If non-overlapping reads were assigned a distance of 0 (identical) or an intermediate value, they might erroneously cluster with reads they have no relationship with.
2. **Conservative Approach**: Assigning the maximum distance (1.0) ensures that we only cluster reads together when there is positive evidence of similarity. "No evidence" is treated conservatively as "potentially different".
3. **Visual Separation**: In the heatmap, this forces non-overlapping reads to be pushed apart or pushed to the edges of clusters, rather than forming false tight clusters.

### Interpretation Guide

When interpreting the distance heatmaps:

- **Distance ≈ 0.0 (Dark Blue)**: Strong evidence of identical methylation patterns.
- **Distance ≈ 1.0 (Bright Yellow/Red)**: Indicates **EITHER**:
    1. Strong evidence of different methylation patterns (high NHD), **OR**
    2. **Insufficient information/overlap** to judge similarity.
