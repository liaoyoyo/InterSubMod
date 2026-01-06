#!/bin/bash
# Script to run Batch VCF Analysis
#
# This script executes run_vcf_all_snv.sh for multiple VCF inputs in parallel if possible,
# and then runs a Python tool to analyze and compare the results.
#
# Usage:
#   ./scripts/run_batch_vcf_analysis.sh
#

set -e

# ============================================================================
# Configuration
# ============================================================================

# Default settings
THREADS=120
MODE="all-with-w5000"
METRICS="BERNOULLI"
PLOT_TYPE="all"

# Output Base Directory
# Note: The script will create a subdirectory {YYYYMMDD}_MODE_{SEQ} inside this base
BASE_OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output"

# Input VCFs (Label VCF_Path)
# defined as associative array or just parallel arrays. 
# Using parallel arrays for simplicity in bash.
LABELS=("filtered_snv_tp" "filtered_snv_fp")
VCF_PATHS=(
    "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
    "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz"
)

# Scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINGLE_RUN_SCRIPT="${SCRIPT_DIR}/run_vcf_all_snv.sh"
PYTHON_TOOL="${SCRIPT_DIR}/../tools/compare_vcf_results.py"

# Date string
DATE_STR=$(date +%Y%m%d)

# ============================================================================
# Output Directory Setup
# ============================================================================

# Ensure base exists
mkdir -p "${BASE_OUTPUT_DIR}"

# Generate unique output directory name: {YYYYMMDD}_MODE_{SEQ}
# e.g., 20260106_all-with-w5000_1
SEQ=1
while true; do
    OUTPUT_DIR="${BASE_OUTPUT_DIR}/${DATE_STR}_${MODE}_${SEQ}"
    if [ ! -d "${OUTPUT_DIR}" ]; then
        break
    fi
    SEQ=$((SEQ + 1))
done

mkdir -p "${OUTPUT_DIR}"
echo "=== Batch VCF Analysis ==="
echo "Output Directory: ${OUTPUT_DIR}"
echo "Mode: ${MODE}"
echo "Threads: ${THREADS}"
echo "Metrics: ${METRICS}"
echo "--------------------------"

# ============================================================================
# Execute VCF Runs
# ============================================================================

PIDS=()
SUB_DIRS=()

for i in "${!LABELS[@]}"; do
    LABEL="${LABELS[$i]}"
    VCF="${VCF_PATHS[$i]}"
    SUB_OUTPUT_DIR="${OUTPUT_DIR}/${LABEL}"
    
    echo ""
    echo "[Batch] Starting analysis for: ${LABEL}"
    echo "        VCF: ${VCF}"
    echo "        Dir: ${SUB_OUTPUT_DIR}"
    
    # Check if VCF exists
    if [ ! -f "${VCF}" ]; then
        echo "Error: VCF file not found: ${VCF}"
        continue
    fi
    
    SUB_DIRS+=("${SUB_OUTPUT_DIR}")
    
    # Run the single analysis script in background
    # We pass the calculated sub-directory as the output dir for this run
    "${SINGLE_RUN_SCRIPT}" \
        --vcf "${VCF}" \
        --out "${SUB_OUTPUT_DIR}" \
        --threads "${THREADS}" \
        --mode "${MODE}" \
        --metrics "${METRICS}" \
        --plot-type "${PLOT_TYPE}" \
        > "${OUTPUT_DIR}/${LABEL}.log" 2>&1 &
        
    PIDS+=($!)
    echo "        Job launched with PID: $!"
done

echo ""
echo "[Batch] Waiting for ${#PIDS[@]} jobs to complete..."
echo "        (Logs are being saved to ${OUTPUT_DIR}/<label>.log)"

# Wait for all background jobs
wait

echo ""
echo "[Batch] All VCF analysis jobs completed."

# ============================================================================
# Run Comparison Analysis
# ============================================================================

ANALYSIS_DIR="${OUTPUT_DIR}/analysis"
mkdir -p "${ANALYSIS_DIR}"

echo ""
echo "[Batch] Running Python Comparison Analysis..."
echo "        Tool: ${PYTHON_TOOL}"

# Check python dependencies
if python3 -c "import pandas, seaborn" 2>/dev/null; then
    # Construct arguments for python script
    # --labels L1 L2 --paths P1 P2
    
    # Need to verify if the output directories were actually created (script might have failed)
    VALID_LABELS=()
    VALID_PATHS=()
    
    for i in "${!LABELS[@]}"; do
        LABEL="${LABELS[$i]}"
        PATH="${SUB_DIRS[$i]}"
        
        if [ -d "${PATH}" ]; then
            VALID_LABELS+=("${LABEL}")
            VALID_PATHS+=("${PATH}")
        else
            echo "Warning: Output directory for ${LABEL} missing, skipping analysis for this one."
        fi
    done
    
    if [ ${#VALID_LABELS[@]} -gt 0 ]; then
        python3 "${PYTHON_TOOL}" \
            --output-dir "${ANALYSIS_DIR}" \
            --labels "${VALID_LABELS[@]}" \
            --paths "${VALID_PATHS[@]}"
    else
        echo "Error: No valid output directories to analyze."
    fi
else
    echo "Warning: Python dependencies (pandas, seaborn) missing. Skipping comparison plotting."
fi

echo ""
echo "=== Batch Analysis Workflow Complete ==="
echo "Results located at: ${OUTPUT_DIR}"

