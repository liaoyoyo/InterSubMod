# Full VCF Test Execution Report

## 1. Introduction

This report documents the execution and results of the full VCF test for the InterSubMod project. The goal is to verify the system's ability to process a real-world somatic VCF file, generating methylation matrices and clustering results for all identified somatic SNVs.

## 2. Test Environment

* **Date**: 2025-11-24
* **Script**: `scripts/run_full_vcf_test.sh`
* **Input VCF**: `/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz`
* **Threads**: 64
* **Output Directory**: `/big8_disk/liaoyoyo2001/InterSubMod/output/full_vcf_test`

## 3. Process

The test follows a three-step process defined in the `run_full_vcf_test.sh` script:

1. **SNV Table Generation**:
    * Extracts chromosome, position, reference, alternate allele, and quality from the input VCF.
    * Generates a tab-separated value (TSV) file `full_snvs.tsv`.

2. **Execution**:
    * Runs the `test_phase4_5` executable.
    * Parameters: Input TSV, Output Directory, Thread Count.
    * The executable performs the following for each SNV:
        * Defines a Region of Interest (ROI) around the SNV (±2000bp).
        * Fetches reads from BAM files (Tumor/Normal) overlapping the ROI.
        * Parses methylation (MM/ML tags) and builds a Read x CpG matrix.
        * Writes output files for each region.

3. **Verification**:
    * Checks if the number of output region directories matches the number of input SNVs.

## 4. Technologies Used

The core system is built using high-performance C++ technologies:

* **Language**: C++17
* **Parallelism**: OpenMP (for processing multiple regions concurrently)
* **Bioinformatics Libraries**: HTSlib (for efficient VCF and BAM file handling)
* **Math/Matrix**: Eigen3 (for high-performance matrix operations and distance calculations)
* **Memory Management**: Jemalloc (linked statically for optimized memory allocation in multi-threaded environments)
* **Build System**: CMake

## 5. Results

* **Total SNVs Processed**: 30,490
* **Total Execution Time**: 3 minutes 49 seconds
* **Success Rate**: 100% (30,490/30,490 regions generated)

### Output Verification

* **Expected Regions**: 30,490
* **Generated Regions**: 30,490
* **Output Structure**:
    For each SNV, a directory structure is created: `output_dir/chr_pos/chr_start_end/`.
    Each directory contains:
  * `metadata.txt`: Summary of the region, SNV, and processing stats.
  * `reads.tsv`: Detailed information for each read (ID, mapping quality, HP tag, Alt support).
  * `cpg_sites.tsv`: List of CpG sites within the region.
  * `methylation.csv`: The Read x CpG methylation matrix (values 0.0-1.0, or NA).

* **Note**: Clustering and association analysis results (e.g., trees, cluster labels) are not currently output by the `test_phase4_5` executable in this configuration.
