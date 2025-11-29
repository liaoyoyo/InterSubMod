# Strand-Aware Processing Design

## 1. Overview

This document details the design for implementing strand-aware processing in the InterSubMod system. As identified in the verification analysis (`docs/reports/verification/2025_11_30_verification_analysis.md`), distinguishing between forward and reverse strands is critical for filtering sequencing artifacts and verifying biological signals.

## 2. Motivation

* **Artifact Removal**: Sequencing errors often exhibit strand bias (appearing only on one strand).
* **Biological Validity**: Hemimethylation events can be detected by analyzing strands separately.
* **Quality Control**: A robust filter requiring support from both strands (e.g., >3 reads on forward, >4 reads on reverse) significantly reduces false positives.

## 3. Data Structure Changes

### 3.1 `ReadInfo` Struct

Add a `strand` field to the `ReadInfo` structure in `include/core/DataStructs.hpp`.

```cpp
enum class Strand : uint8_t {
    FORWARD = 0,
    REVERSE = 1,
    UNKNOWN = 2
};

struct ReadInfo {
    // ... existing fields ...
    Strand strand;
};
```

### 3.2 Output Format (`reads.tsv`)

Update the output format to include a `strand` column.

```tsv
read_id read_name chr start end mapq hp alt_support is_tumor strand
0 read_1 chr19 100 200 60 1 ALT 1 +
1 read_2 chr19 110 210 60 1 REF 1 -
```

## 4. Implementation Logic

### 4.1 Parsing (`ReadParser.cpp`)

* Check `bam1_t::core.flag`.
* If `flag & BAM_FREVERSE`, set `strand = Strand::REVERSE`.
* Otherwise, set `strand = Strand::FORWARD`.

### 4.2 Filtering Strategy

* **Strict Mode**: Require minimum read support on *both* strands for a variant/methylation pattern to be considered valid.
* **Loose Mode**: Allow single-strand support but flag it.

## 5. Verification Plan

1. Run `run_full_vcf_test.sh`.
2. Check `reads.tsv` for the new `strand` column.
3. Verify that the distribution of strands matches `samtools view` statistics.
