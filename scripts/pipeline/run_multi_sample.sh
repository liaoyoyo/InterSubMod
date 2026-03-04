#!/bin/bash
# run_multi_sample.sh - Run benchmark pipeline for multiple samples sequentially
#
# Usage:
#   ./scripts/pipeline/run_multi_sample.sh [SAMPLE1] [SAMPLE2] ...
#   ./scripts/pipeline/run_multi_sample.sh --dry-run [SAMPLE1] ...
#
# Example:
#   ./scripts/pipeline/run_multi_sample.sh HCC1395 HCC1395_DORADO

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SCRIPT="${SCRIPT_DIR}/run_benchmark.sh"

# Default samples if none provided
DEFAULT_SAMPLES=("HCC1395" "HCC1395_DORADO")

# Parse arguments
SAMPLES=()
DRY_RUN=false
MODE="s-pure"

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN=true; shift ;;
        --mode)    MODE="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--dry-run] [--mode MODE] [SAMPLE1] [SAMPLE2] ..."
            exit 0
            ;;
        *)
            if [[ "$1" == -* ]]; then
                echo "Unknown argument: $1" >&2
                exit 1
            fi
            SAMPLES+=("$1")
            shift
            ;;
    esac
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    SAMPLES=("${DEFAULT_SAMPLES[@]}")
fi

echo "================================================================="
echo "  Multi-Sample Benchmark Pipeline"
echo "================================================================="
echo "Samples: ${SAMPLES[*]}"
echo "Mode:    ${MODE}"
echo "Dry-run: ${DRY_RUN}"
echo "-----------------------------------------------------------------"

for sample in "${SAMPLES[@]}"; do
    echo "[Multi-Sample] Starting pipeline for: ${sample}"
    
    CMD=("${PIPELINE_SCRIPT}" --mode "${MODE}" --sample "${sample}")
    
    if [[ "${DRY_RUN}" == true ]]; then
        CMD+=(--dry-run)
    fi

    # Execute
    echo "[Multi-Sample] Command: ${CMD[*]}"
    
    if [[ "${DRY_RUN}" == true ]]; then
        "${CMD[@]}"
    else
        # Run with error handling
        if "${CMD[@]}"; then
            echo "[Multi-Sample] SUCCESS: ${sample}"
        else
            echo "[Multi-Sample] FAILED: ${sample}"
            # Decide whether to stop or continue. For now, we stop on failure to avoid cascading issues.
            exit 1
        fi
    fi
    
    echo "-----------------------------------------------------------------"
done

echo "[Multi-Sample] All samples completed successfully."
