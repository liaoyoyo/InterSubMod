#!/bin/bash
# quick_validate.sh - Single-sample quick validation (~5 min)
#
# Runs ISM on one sample (default: HCC1395) using existing haplotagged BAMs.
# Skips LongPhase-S to save time. Produces metrics.json for comparison.
#
# Usage:
#   ./scripts/validation/runners/quick_validate.sh --experiment-id EXP_xxx
#   ./scripts/validation/runners/quick_validate.sh --sample COLO829 --experiment-id EXP_xxx
#   ./scripts/validation/runners/quick_validate.sh --dry-run --experiment-id EXP_xxx

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../validate_config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

SAMPLE="${QUICK_SAMPLE}"
MODE="${QUICK_MODE}"
EXPERIMENT_ID=""
DRY_RUN=false
LOCAL_THREADS="${THREADS}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)           SAMPLE="$2"; shift 2 ;;
        --mode)             MODE="$2"; shift 2 ;;
        --experiment-id)    EXPERIMENT_ID="$2"; shift 2 ;;
        --threads)          LOCAL_THREADS="$2"; shift 2 ;;
        --dry-run)          DRY_RUN=true; shift ;;
        -h|--help)
            cat <<'HELP'
Usage: quick_validate.sh [OPTIONS]

Options:
  --sample SAMPLE         Sample to validate (default: HCC1395)
  --mode MODE             Pipeline mode: s-pure-pileup, s-pure, to-pure, to-full
  --experiment-id ID      Experiment identifier (required)
  --threads N             Override thread count
  --dry-run               Print commands without executing
  -h, --help              Show this help
HELP
            exit 0
            ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${EXPERIMENT_ID}" ]]; then
    EXPERIMENT_ID="$(generate_experiment_id)"
    log_validate "Auto-generated experiment ID: ${EXPERIMENT_ID}"
fi

# ============================================================================
# Pre-flight Checks
# ============================================================================

log_validate "Quick validation: sample=${SAMPLE}, mode=${MODE}"

# Check disk space
if ! require_min_free_gb "/big7_disk" "${VALIDATE_MIN_FREE_GB}" "quick_validate"; then
    exit 1
fi

# Determine output directory
EXPERIMENT_OUTPUT="${EXPERIMENT_OUTPUT_ROOT}/${EXPERIMENT_ID}"
mkdir -p "${EXPERIMENT_OUTPUT}"

# ============================================================================
# Run Benchmark Pipeline (skip LongPhase)
# ============================================================================

BENCHMARK_CMD=(
    "${PIPELINE_DIR}/run_benchmark.sh"
    --mode "${MODE}"
    --sample "${SAMPLE}"
    --threads "${LOCAL_THREADS}"
    --skip-longphase
    --skip-cleanup
    --output-root "${EXPERIMENT_OUTPUT}"
    --run-tag "${EXPERIMENT_ID}"
)

if [[ "${DRY_RUN}" == "true" ]]; then
    BENCHMARK_CMD+=(--dry-run)
fi

log_validate "Running: ${BENCHMARK_CMD[*]}"

RUN_START=$(date +%s)

if "${BENCHMARK_CMD[@]}"; then
    RUN_STATUS="success"
else
    RUN_STATUS="failed"
    log_validate "Benchmark pipeline FAILED for ${SAMPLE}"
fi

RUN_END=$(date +%s)
RUN_TIME=$(( RUN_END - RUN_START ))

log_validate "Quick validation completed: status=${RUN_STATUS}, time=${RUN_TIME}s"

# ============================================================================
# Output Summary
# ============================================================================

echo "EXPERIMENT_ID=${EXPERIMENT_ID}"
echo "EXPERIMENT_OUTPUT=${EXPERIMENT_OUTPUT}"
echo "SAMPLE=${SAMPLE}"
echo "MODE=${MODE}"
echo "RUN_STATUS=${RUN_STATUS}"
echo "RUN_TIME=${RUN_TIME}"

if [[ "${RUN_STATUS}" != "success" ]]; then
    exit 1
fi
