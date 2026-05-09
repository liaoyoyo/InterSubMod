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

# v1.8 T1-2: enforce per-process TMPDIR + free-space pre-flight (prevent /tmp disasters).
# See InterSubMod/scripts/lib/tmpdir_guard.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib/tmpdir_guard.sh"

# ============================================================================
# Configuration
# ============================================================================

# Default settings
THREADS=120
MODE="all-with-w5000"
METRICS="BERNOULLI"
PLOT_TYPE="no"
RUN_COMPARISON=true
LOG_LEVEL="info"

# Caller mode: "to" (ClairS-TO, tumor-only) or "paired" (ClairS paired, tumor+normal)
# IMPORTANT: TO mode 使用 ClairS-TO VCF + V3-Fixed haplotag BAM
#            Paired mode 使用 ClairS paired pileup VCF + paired haplotag BAM
#            兩者 VCF 不可混用（TO 無 pileup 模型，TP/FP 數量與性質不同）
CALLER_MODE="to"

# Output Base Directory
# Note: The script will create a subdirectory {YYYYMMDD}_MODE_{SEQ} inside this base
BASE_OUTPUT_DIR="/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output"

# --- ClairS-TO (tumor-only) VCFs ---
TO_LABELS=("clairsto_tp" "clairsto_fp")
TO_VCF_PATHS=(
    "/big7_disk/liaoyoyo2001/longphase-to-mod/output/clairsto_v3fixed_work/clairsto_tp.vcf.gz"
    "/big7_disk/liaoyoyo2001/longphase-to-mod/output/clairsto_v3fixed_work/clairsto_fp.vcf.gz"
)

# --- ClairS paired VCFs ---
PAIRED_LABELS=("filtered_snv_tp" "filtered_snv_fp")
PAIRED_VCF_PATHS=(
    "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
    "/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz"
)

# Scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINGLE_RUN_SCRIPT="${SCRIPT_DIR}/run_vcf_all_snv.sh"
PYTHON_TOOL="${SCRIPT_DIR}/../tools/compare_vcf_results.py"
PYTHON_EXEC="/usr/bin/python3"

# Argument parsing
while [[ $# -gt 0 ]]; do
    case $1 in
        --caller-mode)
            CALLER_MODE="$2"
            shift 2
            ;;
        --skip-comparison)
            RUN_COMPARISON=false
            shift
            ;;
        --plot-type)
            PLOT_TYPE="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --log-level)
            LOG_LEVEL="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            shift
            ;;
    esac
done

# Set VCF paths based on caller mode
if [[ "$CALLER_MODE" == "to" ]]; then
    LABELS=("${TO_LABELS[@]}")
    VCF_PATHS=("${TO_VCF_PATHS[@]}")
elif [[ "$CALLER_MODE" == "paired" ]]; then
    LABELS=("${PAIRED_LABELS[@]}")
    VCF_PATHS=("${PAIRED_VCF_PATHS[@]}")
else
    echo "Error: Unknown caller mode '$CALLER_MODE'. Use 'to' or 'paired'." >&2
    exit 1
fi

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
echo "Caller Mode: ${CALLER_MODE}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Mode: ${MODE}"
echo "Threads: ${THREADS}"
echo "Metrics: ${METRICS}"
echo "Log Level: ${LOG_LEVEL}"
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
        --caller-mode "${CALLER_MODE}" \
        --vcf "${VCF}" \
        --out "${SUB_OUTPUT_DIR}" \
        --threads "${THREADS}" \
        --mode "${MODE}" \
        --metrics "${METRICS}" \
        --log-level "${LOG_LEVEL}" \
        --plot-type "${PLOT_TYPE}" \
        > "${OUTPUT_DIR}/${LABEL}.log" 2>&1 &
        
    PIDS+=($!)
    echo "        Job launched with PID: $!"
done

echo ""
echo "[Batch] Waiting for ${#PIDS[@]} jobs to complete..."
echo "        (Logs are being saved to ${OUTPUT_DIR}/<label>.log)"

# Wait for all background jobs and track failures
FAILED=0
for i in "${!PIDS[@]}"; do
    PID="${PIDS[$i]}"
    LABEL="${LABELS[$i]}"
    if ! wait "$PID"; then
        FAILED=$((FAILED + 1))
        echo "[Batch] WARNING: Job '${LABEL}' (PID ${PID}) failed"
    fi
done

echo ""
if [ $FAILED -gt 0 ]; then
    echo "[Batch] WARNING: ${FAILED}/${#PIDS[@]} job(s) failed. Check individual logs for details."
else
    echo "[Batch] All ${#PIDS[@]} VCF analysis jobs completed successfully."
fi

# ============================================================================
# Run Comparison Analysis
# ============================================================================

ANALYSIS_DIR="${OUTPUT_DIR}/analysis"
mkdir -p "${ANALYSIS_DIR}"

if [ "$RUN_COMPARISON" = true ]; then
    echo ""
    echo "[Batch] Running Python Comparison Analysis..."
    echo "        Tool: ${PYTHON_TOOL}"

    # Check python dependencies
    if "${PYTHON_EXEC}" -c "import pandas, seaborn" 2>/dev/null; then
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
            "${PYTHON_EXEC}" "${PYTHON_TOOL}" \
                --output-dir "${ANALYSIS_DIR}" \
                --labels "${VALID_LABELS[@]}" \
                --paths "${VALID_PATHS[@]}"
        else
            echo "Error: No valid output directories to analyze."
        fi
    else
        echo "Warning: Python dependencies (pandas, seaborn) missing. Skipping comparison plotting."
    fi
else
    echo ""
    echo "[Batch] Skipping Comparison Analysis (--skip-comparison used)."
fi

echo ""
echo "=== Batch Analysis Workflow Complete ==="
echo "Results located at: ${OUTPUT_DIR}"

