#!/bin/bash
# Script to run ALL SNV tests from a real VCF file

# Defaults
VCF_PATH="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output/full_vcf_test"
THREADS=64
TUMOR_BAM="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL_BAM="/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF_FASTA="/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
LOG_FILE="${OUTPUT_DIR}/full_execution_analysis.log"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v|--vcf)
            VCF_PATH="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -o|--out)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Update log file path if output dir changed
LOG_FILE="${OUTPUT_DIR}/full_execution_analysis.log"

echo "=== InterSubMod FULL VCF Test ==="
echo "VCF: ${VCF_PATH}"
echo "Tumor BAM: ${TUMOR_BAM}"
echo "Normal BAM: ${NORMAL_BAM}"
echo "Reference: ${REF_FASTA}"
echo "Output Dir: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Log File: ${LOG_FILE}"
echo "---------------------------------"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Step 1: Run Execution
echo ""
echo "[1] Running InterSubMod..."
EXECUTABLE="/big8_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod"

if [ ! -f "${EXECUTABLE}" ]; then
    echo "Error: Executable not found at ${EXECUTABLE}"
    echo "Please build the project first."
    exit 1
fi

# Construct Command
CMD="${EXECUTABLE} \
    --tumor-bam ${TUMOR_BAM} \
    --normal-bam ${NORMAL_BAM} \
    --reference ${REF_FASTA} \
    --vcf ${VCF_PATH} \
    --output-dir ${OUTPUT_DIR} \
    --threads ${THREADS}"

echo "Command to run (for IDE Debugging):"
echo "${CMD}"
echo ""

# Use /usr/bin/time for resource monitoring
echo "    Starting execution with ${THREADS} threads..."
echo "    (This may take a few minutes...)"
{
    /usr/bin/time -v ${CMD}
} 2>&1 | tee "${LOG_FILE}"

echo ""
echo "=== Analysis Log Saved to: ${LOG_FILE} ==="
echo "Done."

