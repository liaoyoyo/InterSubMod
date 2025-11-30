# Strand-Aware Processing & Debug Mode Implementation Report

**Date**: 2025-11-30  
**Version**: 1.0  
**Status**: Completed

## 1. Executive Summary

This report documents the implementation of strand-aware processing and debug logging features for the InterSubMod project, as specified in the verification documents.

### Key Accomplishments
- ✅ Added `strand` field to `ReadInfo` structure and all outputs
- ✅ Implemented strand-separated methylation matrices (`methylation_forward.csv`, `methylation_reverse.csv`)
- ✅ Added `--log-level` parameter supporting ERROR/WARN/INFO/DEBUG levels
- ✅ Implemented filtered read logging with detailed filter reasons
- ✅ Added `--no-filter` mode for verification purposes
- ✅ Created two test modes: baseline and all-with-window
- ✅ Successfully executed both test modes with complete output

## 2. Implementation Details

### 2.1 Strand Information

#### Data Structure Changes

**`include/core/Types.hpp`**
```cpp
enum class Strand : uint8_t {
    FORWARD = 0,  ///< Forward strand (+)
    REVERSE = 1,  ///< Reverse strand (-)
    UNKNOWN = 2   ///< Unknown strand
};
```

**`include/core/DataStructs.hpp`**
```cpp
struct ReadInfo {
    // ... existing fields ...
    Strand strand;  ///< Strand orientation (FORWARD/+ or REVERSE/-)
};
```

#### Strand Determination Logic

In `ReadParser.cpp`, strand is determined from BAM FLAG bit 0x10:
```cpp
Strand ReadParser::determine_strand(const bam1_t* b) {
    if (b->core.flag & BAM_FREVERSE) {
        return Strand::REVERSE;
    }
    return Strand::FORWARD;
}
```

#### Output Format Changes

**`reads.tsv`** - New `strand` column added:
```tsv
read_id  read_name  chr  start  end  mapq  hp  alt_support  is_tumor  strand
0        read_1     chr19 100    200  60    1   ALT          1         +
1        read_2     chr19 110    210  60    1   REF          1         -
```

**`metadata.txt`** - Now includes strand statistics:
```
Num Reads: 48
  Forward Strand (+): 21
  Reverse Strand (-): 27
```

### 2.2 Strand-Separated Methylation Matrices

Two new output files are generated:

1. **`methylation_forward.csv`** - Contains only forward strand reads
2. **`methylation_reverse.csv`** - Contains only reverse strand reads

Both files include:
- `read_id`: New sequential ID within the strand-specific matrix
- `original_read_id`: Reference to the full matrix row index
- CpG columns with methylation probabilities

### 2.3 Debug Mode & Logging

#### New Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--log-level` | string | info | Logging level: error, warn, info, debug |
| `--output-filtered-reads` | flag | false | Output filtered reads with reasons |
| `--no-filter` | flag | false | Output all reads without filtering |
| `--debug-output-dir` | string | `<output>/debug` | Directory for debug outputs |

#### Filter Reasons

The system tracks the following filter reasons:

| Reason | Description |
|--------|-------------|
| `SECONDARY` | Secondary alignment |
| `SUPPLEMENTARY` | Supplementary alignment |
| `DUPLICATE` | PCR/optical duplicate |
| `UNMAPPED` | Unmapped read |
| `LOW_MAPQ` | MAPQ below threshold |
| `SHORT_READ` | Read length below threshold |
| `MISSING_MM` | Missing MM tag |
| `MISSING_ML` | Missing ML tag |
| `SNV_NOT_COVERED` | Read does not cover SNV position |
| `SNV_IN_DELETION` | SNV falls in deletion region |
| `LOW_BASE_QUALITY` | Base quality at SNV below threshold |
| `NOT_REF_OR_ALT` | Base at SNV is neither REF nor ALT |

#### Filtered Reads Output

**`filtered_reads.tsv`**:
```tsv
read_name  chr    start     end       mapq  strand  is_tumor  filter_reasons
read_1     chr19  10048296  10100752  60    -       1         SNV_NOT_COVERED
```

## 3. Test Execution Results

### 3.1 Baseline Mode

**Command:**
```bash
./scripts/run_full_vcf_test.sh --mode baseline --threads 32
```

**Output Directory:** `output/vcf_baseline_1130`

**Results:**
| Metric | Value |
|--------|-------|
| Total Regions | 30,490 |
| Total Reads | 2,176,225 |
| Forward Strand (+) | 1,086,952 |
| Reverse Strand (-) | 1,089,273 |
| Filtered Reads | 1,568,851 |
| Total CpG Sites | 602,912 |
| Wall-clock Time | 98.05 s |
| Average Reads/Region | 71.4 |

### 3.2 All-with-Window Mode

**Command:**
```bash
./scripts/run_full_vcf_test.sh --mode all-with-window --threads 32
```

**Output Directory:** `output/vcf_all_w2000_1130`

**Settings:**
- Window size: ±2000 bp
- No filtering (--no-filter)

**Results:**
| Metric | Value |
|--------|-------|
| Total Regions | 30,490 |
| Total Reads | 4,133,855 |
| Forward Strand (+) | 2,071,818 |
| Reverse Strand (-) | 2,062,037 |
| Filtered Reads | 0 |
| Total CpG Sites | 1,196,722 |
| Wall-clock Time | 81.27 s |
| Average Reads/Region | 135.6 |

### 3.3 Comparison Analysis

| Metric | Baseline | All-with-Window | Ratio |
|--------|----------|-----------------|-------|
| Reads | 2,176,225 | 4,133,855 | 1.90x |
| CpG Sites | 602,912 | 1,196,722 | 1.98x |
| Avg Reads/Region | 71.4 | 135.6 | 1.90x |

The all-with-window mode captured approximately **1.9x more reads** due to:
1. Larger window size (±2000 bp vs ±1000 bp)
2. No filtering applied (includes UNKNOWN alt_support reads)

## 4. Output Directory Structure

```
output/vcf_baseline_1130/
├── chr19_10084871/
│   └── chr19_10083871_10085871/
│       ├── metadata.txt              # Region & SNV info with strand stats
│       ├── reads.tsv                 # Read list with strand column
│       ├── cpg_sites.tsv             # CpG positions
│       ├── methylation.csv           # Full methylation matrix
│       ├── methylation_forward.csv   # Forward strand matrix
│       ├── methylation_reverse.csv   # Reverse strand matrix
│       └── filtered_reads.tsv        # Filtered reads with reasons
├── chr19_10085272/
│   └── ...
└── full_execution_analysis.log       # Execution log
```

## 5. Code Changes Summary

### Modified Files

| File | Changes |
|------|---------|
| `include/core/Types.hpp` | Added `Strand`, `LogLevel`, `FilterReason` enums |
| `include/core/DataStructs.hpp` | Added `strand` to `ReadInfo`, new `FilteredReadInfo` struct |
| `include/core/Config.hpp` | Added log level and debug parameters |
| `include/core/ReadParser.hpp` | Added `determine_strand()`, filter reason tracking |
| `src/core/ReadParser.cpp` | Implemented strand parsing and filter reasons |
| `include/io/RegionWriter.hpp` | Added strand matrix and filtered reads output |
| `src/io/RegionWriter.cpp` | Implemented new output methods |
| `include/core/RegionProcessor.hpp` | Added debug mode support |
| `src/core/RegionProcessor.cpp` | Integrated strand and debug features |
| `include/utils/ArgParser.hpp` | Added new CLI parameters |
| `src/main.cpp` | Updated to use Config-based constructor |
| `scripts/run_full_vcf_test.sh` | Added baseline and all-with-window modes |

## 6. Usage Examples

### Baseline Mode (Normal Filtering + Debug)
```bash
./build/bin/inter_sub_mod \
    --tumor-bam data/bam/tumor.bam \
    --normal-bam data/bam/normal.bam \
    --reference data/ref/hg38.fa \
    --vcf data/vcf/snvs.vcf.gz \
    --output-dir output/baseline \
    --threads 32 \
    --log-level debug \
    --output-filtered-reads
```

### All-with-Window Mode (No Filtering)
```bash
./build/bin/inter_sub_mod \
    --tumor-bam data/bam/tumor.bam \
    --normal-bam data/bam/normal.bam \
    --reference data/ref/hg38.fa \
    --vcf data/vcf/snvs.vcf.gz \
    --output-dir output/all_reads \
    --threads 32 \
    --window-size 2000 \
    --log-level debug \
    --no-filter
```

## 7. Verification Checklist

- [x] Read verification documentation and implemented according to specs
- [x] All outputs include consistent strand information (+/-)
- [x] Default output separates reads by strand (forward/reverse matrices)
- [x] Debug mode outputs filtered reads with detailed reasons
- [x] `run_full_vcf_test.sh` successfully executes:
  - [x] Baseline mode → `output/vcf_baseline_1130`
  - [x] All-with-window mode → `output/vcf_all_w2000_1130`
- [x] Log content sufficient for debugging and verification

## 8. Conclusion

The implementation successfully addresses all requirements specified in the verification documents:

1. **Strand Information**: Complete implementation with forward/reverse strand tracking
2. **Debug Mode**: Comprehensive logging with filter reason tracking
3. **Test Modes**: Both baseline and all-with-window modes functional

The system is ready for production use with enhanced debugging capabilities.

---

**References:**
- `docs/reports/verification/2025_11_30_verification_analysis.md`
- `docs/reports/verification/2025_11_30_output_integrity_check.md`
- `docs/architecture/strand_aware_processing_design.md`

