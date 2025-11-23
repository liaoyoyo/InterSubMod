#!/bin/bash
# Script to run ALL SNV tests from a real VCF file

# Defaults
VCF_PATH="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output_full"
THREADS=64
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
echo "Output Dir: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Log File: ${LOG_FILE}"
echo "---------------------------------"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Step 1: Generate Full SNV Table
echo "[1] Generating full SNV table from VCF..."
TEMP_TSV="${OUTPUT_DIR}/full_snvs.tsv"

# Write Header
echo -e "chr\tpos\tref\talt\tqual" > "${TEMP_TSV}"

# Extract ALL lines from VCF (excluding header)
if [ -f "${VCF_PATH}" ]; then
    if [[ "${VCF_PATH}" == *.gz ]]; then
        CAT_CMD="zgrep"
    else
        CAT_CMD="grep"
    fi
    
    # Extract all variants
    ${CAT_CMD} -v "^#" "${VCF_PATH}" | awk '{print $1, $2, $4, $5, $6}' >> "${TEMP_TSV}"
    
    NUM_SNVS=$(($(wc -l < "${TEMP_TSV}") - 1))
    echo "    Generated: ${TEMP_TSV}"
    echo "    Total SNVs to process: ${NUM_SNVS}"
else
    echo "Error: VCF file not found at ${VCF_PATH}"
    exit 1
fi

# Step 2: Run Execution
echo ""
echo "[2] Running InterSubMod..."
EXECUTABLE="/big8_disk/liaoyoyo2001/InterSubMod/build/bin/test_phase4_5"

if [ ! -f "${EXECUTABLE}" ]; then
    echo "Error: Executable not found at ${EXECUTABLE}"
    echo "Please build the project first."
    exit 1
fi

# Use /usr/bin/time for resource monitoring
echo "    Starting execution with ${THREADS} threads..."
echo "    (This may take a few minutes...)"
{
    /usr/bin/time -v "${EXECUTABLE}" "${TEMP_TSV}" "${OUTPUT_DIR}" "${THREADS}"
} 2>&1 | tee "${LOG_FILE}"

# Step 3: Basic Verification
echo ""
echo "[3] Verifying results..."
REGION_COUNT=$(ls -d "${OUTPUT_DIR}/region_"* 2>/dev/null | wc -l)

echo "    Expected regions: ${NUM_SNVS}"
echo "    Found regions:    ${REGION_COUNT}"

if [ "${REGION_COUNT}" -eq "${NUM_SNVS}" ]; then
    echo "    ✓ Region count matches."
else
    echo "    ⚠ Region count mismatch!"
fi

echo ""
echo "=== Analysis Log Saved to: ${LOG_FILE} ==="
echo "Done."

