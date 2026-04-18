#!/bin/bash
# full_validate.sh - Multi-sample parallel validation
#
# Runs ISM on all 7 samples in parallel using existing haplotagged BAMs.
# Each sample runs as a background job; waits for all to complete.
#
# Usage:
#   ./scripts/validation/runners/full_validate.sh --experiment-id EXP_xxx
#   ./scripts/validation/runners/full_validate.sh --experiment-id EXP_xxx --samples "HCC1395,COLO829"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../validate_config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

EXPERIMENT_ID=""
MODE="${FULL_MODE}"
CUSTOM_SAMPLES=""
DRY_RUN=false
LOCAL_THREADS="${THREADS}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --experiment-id)    EXPERIMENT_ID="$2"; shift 2 ;;
        --mode)             MODE="$2"; shift 2 ;;
        --samples)          CUSTOM_SAMPLES="$2"; shift 2 ;;
        --threads)          LOCAL_THREADS="$2"; shift 2 ;;
        --dry-run)          DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: full_validate.sh --experiment-id ID [--samples S1,S2,...] [--mode MODE]"
            exit 0
            ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${EXPERIMENT_ID}" ]]; then
    EXPERIMENT_ID="$(generate_experiment_id)"
fi

# Determine samples
if [[ -n "${CUSTOM_SAMPLES}" ]]; then
    IFS=',' read -ra SAMPLES <<< "${CUSTOM_SAMPLES}"
else
    SAMPLES=("${ALL_SAMPLES[@]}")
fi

# ============================================================================
# Pre-flight
# ============================================================================

EXPERIMENT_OUTPUT="${EXPERIMENT_OUTPUT_ROOT}/${EXPERIMENT_ID}"
mkdir -p "${EXPERIMENT_OUTPUT}"

TOTAL_NEEDED=$(( ${#SAMPLES[@]} * ESTIMATED_GB_PER_SAMPLE ))
if ! require_min_free_gb "/big7_disk" "${VALIDATE_MIN_FREE_GB}" "full_validate (need ~${TOTAL_NEEDED}GB)"; then
    exit 1
fi

# Calculate per-job threads
NJOBS=${#SAMPLES[@]}
MAX_PARALLEL=${MAX_PARALLEL_VALIDATE}
if (( NJOBS < MAX_PARALLEL )); then
    MAX_PARALLEL=${NJOBS}
fi
THREADS_PER_JOB=$(( LOCAL_THREADS / MAX_PARALLEL ))
if (( THREADS_PER_JOB < 4 )); then
    THREADS_PER_JOB=4
fi

log_validate "Full validation: ${NJOBS} samples, ${MAX_PARALLEL} parallel, ${THREADS_PER_JOB} threads/job"

# ============================================================================
# Run All Samples
# ============================================================================

PIDS=()
SAMPLES_RUNNING=()
FAILED_SAMPLES=()
RUNNING_COUNT=0

for sample in "${SAMPLES[@]}"; do
    # Wait if we've hit max parallel
    while (( RUNNING_COUNT >= MAX_PARALLEL )); do
        # Wait for any child to finish
        wait -n 2>/dev/null || true
        # Recount running
        RUNNING_COUNT=0
        for pid in "${PIDS[@]}"; do
            if kill -0 "${pid}" 2>/dev/null; then
                (( RUNNING_COUNT++ )) || true
            fi
        done
    done

    LOG_FILE="${EXPERIMENT_OUTPUT}/${sample}_validate.log"

    BENCHMARK_CMD=(
        "${PIPELINE_DIR}/run_benchmark.sh"
        --mode "${MODE}"
        --sample "${sample}"
        --threads "${THREADS_PER_JOB}"
        --skip-longphase
        --skip-cleanup
        --output-root "${EXPERIMENT_OUTPUT}"
        --run-tag "${EXPERIMENT_ID}"
    )

    if [[ "${DRY_RUN}" == "true" ]]; then
        BENCHMARK_CMD+=(--dry-run)
    fi

    log_validate "Starting ${sample}... (log: ${LOG_FILE})"

    "${BENCHMARK_CMD[@]}" >"${LOG_FILE}" 2>&1 &
    PIDS+=($!)
    SAMPLES_RUNNING+=("${sample}")
    (( RUNNING_COUNT++ )) || true
done

# ============================================================================
# Wait for All Jobs
# ============================================================================

log_validate "Waiting for ${#PIDS[@]} jobs to complete..."

for i in "${!PIDS[@]}"; do
    pid="${PIDS[$i]}"
    sample="${SAMPLES_RUNNING[$i]}"

    if wait "${pid}"; then
        log_validate "  [OK] ${sample} (pid=${pid})"
    else
        log_validate "  [FAIL] ${sample} (pid=${pid})"
        FAILED_SAMPLES+=("${sample}")
    fi
done

# ============================================================================
# Summary
# ============================================================================

TOTAL=${#SAMPLES[@]}
FAILED=${#FAILED_SAMPLES[@]}
PASSED=$(( TOTAL - FAILED ))

log_validate "Full validation complete: ${PASSED}/${TOTAL} passed"

if (( FAILED > 0 )); then
    log_validate "Failed samples: ${FAILED_SAMPLES[*]}"
    # Don't exit with error — partial results may still be useful for comparison
    log_validate "WARNING: ${FAILED} sample(s) failed. Partial results available."
fi

echo "EXPERIMENT_ID=${EXPERIMENT_ID}"
echo "EXPERIMENT_OUTPUT=${EXPERIMENT_OUTPUT}"
echo "SAMPLES_TOTAL=${TOTAL}"
echo "SAMPLES_PASSED=${PASSED}"
echo "SAMPLES_FAILED=${FAILED}"
