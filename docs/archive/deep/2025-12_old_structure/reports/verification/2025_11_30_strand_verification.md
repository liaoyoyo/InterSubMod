# Verification Report: Strand-Aware Processing and Coordinate Logic

**Date:** 2025-11-30
**Scope:** Verify strand handling, 5mC extraction, and coordinate calculation in `InterSubMod`.

## 1. Executive Summary

The verification confirms that the codebase correctly handles strand detection and 5mC extraction. **The "1-bp offset" between forward and reverse strands is explicitly handled by the code to align reverse strand calls to the forward strand CpG start position.** This means the output files (`methylation.csv`, etc.) present aggregated CpG sites with unified coordinates, eliminating the offset in the final output.

## 2. Detailed Findings

### 2.1 Strand Detection

* **File:** `src/core/ReadParser.cpp`
* **Logic:** Uses `b->core.flag & BAM_FREVERSE` (0x10) to determine strand.
* **Status:** **Correct.** This is the standard method for determining read orientation in BAM files.

### 2.2 5mC Extraction and Target Base

* **File:** `src/core/MethylationParser.cpp`
* **Logic:**
  * **Forward Strand:** The code looks for Cytosine ('C') in the BAM sequence.
  * **Reverse Strand:** The code looks for Guanine ('G') in the BAM sequence.
* **Analysis:**
  * Reads mapped to the reverse strand have their sequence reverse-complemented in the BAM file (per SAM specification).
  * A methylated Cytosine on the reverse strand pairs with a Guanine on the forward strand.
  * Therefore, in the BAM sequence (which matches the forward reference), this position appears as a 'G'.
  * **Status:** **Correct.** The logic correctly targets the base corresponding to the methylated Cytosine on the original strand.

### 2.3 Coordinate Calculation and 1-bp Offset

* **File:** `src/core/MethylationParser.cpp`
* **Logic:**
  * **Forward Strand:** Reports `ref_pos` (0-based) converted to 1-based.
  * **Reverse Strand:** Reports `ref_pos - 1` (0-based) converted to 1-based.
* **Analysis:**
  * A CpG dinucleotide on the forward strand is `5'-C(i)G(i+1)-3'`.
  * The corresponding CpG on the reverse strand is `3'-G(i)C(i+1)-5'`, which reads as `5'-C(i+1)G(i)-3'` on the reverse strand.
  * The methylated Cytosine on the reverse strand is at position `i+1`.
  * The code identifies the 'G' at position `i+1` (in forward coordinates) and subtracts 1, resulting in position `i`.
* **Result:**
  * Both forward and reverse strand methylation calls for the same CpG site are reported at position `i` (the position of the Forward Cytosine).
  * **Status:** **Intentionally Aligned.** The code effectively "merges" the strands by shifting the reverse strand coordinate. This is standard behavior for CpG-level analysis but technically reports the "Forward C" position for reverse strand events.

### 2.4 Output Verification

* **Files:** `methylation_forward.csv` vs `methylation_reverse.csv`
* **Observation:** The column headers (genomic positions) are identical in both files.
* **Conclusion:** The 1-bp offset handling is working as implemented. There is no discrepancy in the output files.

## 3. Potential Changes

If the intention is to report the **exact genomic coordinate of the methylated base** (rather than the CpG site start), the following change is required:

* **Location:** `src/core/MethylationParser.cpp`, line 162.
* **Current Code:**

    ```cpp
    int32_t report_pos = is_reverse ? (ref_pos_0based - 1) : ref_pos_0based;
    ```

* **Proposed Change (for exact C coordinate):**

    ```cpp
    // Report the actual coordinate of the methylated base (which is at ref_pos_0based for Reverse strand 'G' which corresponds to 'C')
    // Wait, if we are at 'G' (pos i+1), the C on reverse strand is at i+1.
    // So we should just report ref_pos_0based.
    int32_t report_pos = ref_pos_0based;
    ```

**Recommendation:**

* **Keep current behavior** if the goal is to analyze CpG sites as units (most common for methylation analysis).
* **Apply change** only if strand-specific resolution at the exact base level is required (e.g., for non-CpG methylation or specific strand asymmetry studies).

## 4. CIGAR Parsing

* **File:** `src/core/ReadParser.cpp`
* **Status:** **Correct.** The CIGAR parsing logic correctly handles matches, insertions, deletions, and soft-clips to map read sequence indices to reference coordinates.
