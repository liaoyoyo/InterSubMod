# Implementation Plan - Bidirectional Verification Workflow
**Date**: 2025-12-30
**Goal**: Implement a "Label-First" vs "Cluster-First" bidirectional verification system for methylation distance matrices.

## 1. Overview
The current system prioritizes clustering (Cluster-First) which is prone to False Positives in sparse data. We will implement a Label-First path (verifying if labels explain distance structure) and integrate it with an improved Cluster-First path (adding stability checks).

## 2. Methodology
- **Label-First**: Compute `Delta` (Between-Within distance) and run PERMANOVA using HP/Allele labels directly.
- **Cluster-First**: Add `Bootstrap/Subsampling` to `HierarchicalClustering` to assess stability before Significance Analysis.
- **Integration**: `SignificanceAnalyzer` will run both paths and classify regions into `Strong`, `Subclone`, `Weak`, or `Artifact`.

## 3. Component Changes

### 3.1 `core/StructureTest.cpp` & `.hpp`
- **[NEW] `compute_delta_and_test`**:
    - Input: `dist_matrix`, `labels`.
    - Logic: Calculate actual `Delta = mean(between) - mean(within)`. Run permutation test (shuffle labels) to get P-value for Delta. This is often more robust and intuitive than PERMANOVA for simple 2-group comparisons.
    - Output: `Delta`, `P_value`, `WithinMean`, `BetweenMean`.

### 3.2 `core/SignificanceAnalyzer.cpp` & `.hpp`
- **[MODIFY] `analyze` / `analyze_simple`**:
    - **Enable Label-First**:
        - Extract HP labels and Allele labels from `full_labels`.
        - Check overlapping/balance (Minimum Requirements).
        - Call `structure_test_->run_permanova` AND `compute_delta` using these labels.
    - **Integrate Results**:
        - Store Label-First results (R2, Delta, P) in `SignificanceResult`.
        - Implement the classification logic (Decision Tree from Analysis Report).

### 3.3 `core/HierarchicalClustering.cpp` (or new `ClusteringStability.cpp`)
- **[NEW] `assess_stability`**:
    - Input: `dist_matrix`, `original_labels`.
    - Logic:
        - Subsample reads (e.g., 80%).
        - Re-run UPGMA.
        - Cut tree at `k` (or same height).
        - Compute Jaccard/ARI between original clusters (restricted to subsampled items) and new clusters.
        - Repeat ~20-50 times.
    - Output: `StabilityScore` (Mean Jaccard/ARI).

### 3.4 `core/RegionProcessor.cpp`
- **[MODIFY] Output Generation**:
    - Update `write_significance_summary` (CSV) to include new columns: `Label_P`, `Label_Delta`, `Stability`, `Class`.
    - Update `significance.json` writer to include detailed structure test results.

## 4. Execution Steps

### Step 1: Label-First Core (C++)
1.  Modify `StructureTest` to add `compute_delta`.
2.  Modify `SignificanceAnalyzer` to prepare HP/Allele specific label vectors and call `StructureTest` methods on them.
3.  Target: `core/StructureTest.cpp`, `core/SignificanceAnalyzer.cpp`.

### Step 2: Stability Check (C++)
1.  Implement lightweight stability check in `Clustering` module.
2.  Integrate into `SignificanceAnalyzer` (before or parallel to GlobalTest).
3.  Target: `core/HierarchicalClustering.cpp`.

### Step 3: Synthesis & Output
1.  Implement the `Class` determination logic in `SignificanceAnalyzer`.
2.  Update `RegionProcessor` to write new fields.
3.  Target: `core/RegionProcessor.cpp`.

## 5. Verification Plan
- **Unit Test**: Create synthetic distance matrices:
    - Case A: Perfect separation by label (Expect: High Delta, Low P, Stable).
    - Case B: Random noise (Expect: Low Delta, High P, Unstable).
    - Case C: Strong cluster orthogonal to label (Expect: Low Delta/P for Label, High Stability for Cluster).
- **Integration Test**: Run on the identified discrepancy samples (`chr4:9649610` FP, `chr8:84676714` TP).
    - Verify `chr4:9649610` (FP) gets classified as `Noise` or `Weak` (due to low stability or low Delta significance).
    - Verify `chr8:84676714` (TP) gets classified as `Strong` (if Delta is significant despite marginal Fisher P).
