#!/bin/bash
# validate.sh - Master entrypoint for single-change validation
#
# Orchestrates: build → run ISM → extract metrics → compare → report → log
#
# Usage:
#   ./scripts/validation/validate.sh --hypothesis "Description" --mode quick
#   ./scripts/validation/validate.sh --hypothesis "Description" --mode full
#   ./scripts/validation/validate.sh --hypothesis "Description" --mode full --samples "HCC1395,COLO829"
#   ./scripts/validation/validate.sh --skip-build --mode quick  # reuse existing build

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/validate_config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

HYPOTHESIS=""
HYPOTHESIS_ID=""
VALIDATION_MODE="quick"
CUSTOM_SAMPLES=""
TRACK=""
PIPELINE_MODE=""
SKIP_BUILD=false
SKIP_COMPARE=false
DRY_RUN=false
LOCAL_THREADS="${THREADS}"

show_help() {
    cat <<'HELP'
Usage: validate.sh [OPTIONS]

Single-change validation framework for InterSubMod.
Builds, runs ISM on canonical samples, compares metrics to baseline.

Options:
  --hypothesis TEXT       Description of the change being tested (required)
  --hypothesis-id ID     Hypothesis identifier (optional, e.g., H_NEW_001)
  --mode MODE            Validation mode: quick (1 sample) or full (7 samples)
  --track TRACK          Analysis track: paired (default) or to (tumor-only)
  --samples LIST         Comma-separated sample list (overrides mode default)
  --pipeline-mode MODE   ISM pipeline mode (overrides --track default)
  --threads N            Override thread count
  --skip-build           Skip build step (reuse existing binary)
  --skip-compare         Skip comparison step (just run ISM)
  --dry-run              Print commands without executing
  -h, --help             Show this help

Examples:
  # Quick feedback during development (paired mode)
  ./scripts/validation/validate.sh \
    --hypothesis "Increase binary_methyl_high from 0.8 to 0.9" \
    --mode quick

  # Full cross-sample validation (paired mode)
  ./scripts/validation/validate.sh \
    --hypothesis "Add normal methylation reference" \
    --hypothesis-id 2A-1 \
    --mode full

  # Tumor-only mode validation
  ./scripts/validation/validate.sh \
    --hypothesis "Improve TO quality score" \
    --track to \
    --mode full
HELP
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --hypothesis)       HYPOTHESIS="$2"; shift 2 ;;
        --hypothesis-id)    HYPOTHESIS_ID="$2"; shift 2 ;;
        --mode)             VALIDATION_MODE="$2"; shift 2 ;;
        --track)            TRACK="$2"; shift 2 ;;
        --samples)          CUSTOM_SAMPLES="$2"; shift 2 ;;
        --pipeline-mode)    PIPELINE_MODE="$2"; shift 2 ;;
        --threads)          LOCAL_THREADS="$2"; shift 2 ;;
        --skip-build)       SKIP_BUILD=true; shift ;;
        --skip-compare)     SKIP_COMPARE=true; shift ;;
        --dry-run)          DRY_RUN=true; shift ;;
        -h|--help)          show_help ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${HYPOTHESIS}" ]]; then
    echo "[ERROR] --hypothesis is required" >&2
    show_help
fi

# Resolve track → pipeline mode (--pipeline-mode overrides --track)
if [[ -z "${PIPELINE_MODE}" ]]; then
    if [[ -n "${TRACK}" ]]; then
        PIPELINE_MODE="${TRACK_PIPELINE_MODE[${TRACK}]:-}"
        if [[ -z "${PIPELINE_MODE}" ]]; then
            echo "[ERROR] Unknown track: ${TRACK} (use: paired, to)" >&2
            exit 1
        fi
    else
        PIPELINE_MODE="${FULL_MODE}"
        TRACK="paired"
    fi
elif [[ -z "${TRACK}" ]]; then
    # Infer track from pipeline-mode
    case "${PIPELINE_MODE}" in
        to-pure|to-full|to_pure|to_full) TRACK="to" ;;
        *) TRACK="paired" ;;
    esac
fi

# ============================================================================
# Initialization
# ============================================================================

EXPERIMENT_ID="$(generate_experiment_id)"
EXPERIMENT_OUTPUT="${EXPERIMENT_OUTPUT_ROOT}/${EXPERIMENT_ID}"
WALL_START=$(date +%s)

# Capture git state
eval "$(get_git_state)"
# COMMIT, BRANCH, DIRTY, DIFF_STAT are now set

log_validate "============================================================"
log_validate "EXPERIMENT: ${EXPERIMENT_ID}"
log_validate "Hypothesis: ${HYPOTHESIS}"
log_validate "Mode: ${VALIDATION_MODE} | Track: ${TRACK} | Pipeline: ${PIPELINE_MODE}"
log_validate "Git: ${BRANCH}@${COMMIT} (dirty=${DIRTY})"
log_validate "============================================================"

mkdir -p "${EXPERIMENT_OUTPUT}"

# Determine samples to run
if [[ -n "${CUSTOM_SAMPLES}" ]]; then
    IFS=',' read -ra SAMPLES_TO_RUN <<< "${CUSTOM_SAMPLES}"
elif [[ "${VALIDATION_MODE}" == "quick" ]]; then
    SAMPLES_TO_RUN=("${QUICK_SAMPLE}")
elif [[ "${VALIDATION_MODE}" == "full" ]]; then
    SAMPLES_TO_RUN=("${ALL_SAMPLES[@]}")
else
    echo "[ERROR] Unknown validation mode: ${VALIDATION_MODE}" >&2
    exit 1
fi

log_validate "Samples: ${SAMPLES_TO_RUN[*]}"

# ============================================================================
# Step 1: Build & Test
# ============================================================================

BUILD_STATUS="skipped"
BUILD_TIME=0
TEST_STATUS="skipped"

if [[ "${SKIP_BUILD}" != "true" ]]; then
    log_validate "Step 1: Build & unit test..."

    BUILD_OUTPUT=$("${SCRIPT_DIR}/runners/build_and_test.sh" 2>&1) || {
        log_validate "BUILD OR TEST FAILED — aborting"
        echo "${BUILD_OUTPUT}"
        exit 1
    }

    # Parse output
    eval "$(echo "${BUILD_OUTPUT}" | grep -E '^(BUILD_|TEST_)')"
    log_validate "Build: ${BUILD_STATUS} (${BUILD_TIME}s) | Tests: ${TEST_STATUS}"
else
    log_validate "Step 1: Skipped (--skip-build)"
fi

# ============================================================================
# Step 2: Run ISM on samples
# ============================================================================

log_validate "Step 2: Running ISM on ${#SAMPLES_TO_RUN[@]} sample(s)..."

if [[ "${VALIDATION_MODE}" == "full" && ${#SAMPLES_TO_RUN[@]} -gt 1 ]]; then
    # Full mode: run in parallel
    "${SCRIPT_DIR}/runners/full_validate.sh" \
        --experiment-id "${EXPERIMENT_ID}" \
        --mode "${PIPELINE_MODE}" \
        --samples "$(IFS=,; echo "${SAMPLES_TO_RUN[*]}")" \
        --threads "${LOCAL_THREADS}" \
        ${DRY_RUN:+--dry-run} || {
        log_validate "FULL VALIDATION FAILED"
        exit 1
    }
else
    # Quick mode or single sample
    for sample in "${SAMPLES_TO_RUN[@]}"; do
        "${SCRIPT_DIR}/runners/quick_validate.sh" \
            --sample "${sample}" \
            --mode "${PIPELINE_MODE}" \
            --experiment-id "${EXPERIMENT_ID}" \
            --threads "${LOCAL_THREADS}" \
            ${DRY_RUN:+--dry-run} || {
            log_validate "VALIDATION FAILED for ${sample}"
            exit 1
        }
    done
fi

# ============================================================================
# Step 3: Extract Metrics
# ============================================================================

log_validate "Step 3: Extracting metrics..."
METRICS_BUNDLE="${EXPERIMENT_OUTPUT}/metrics_bundle.json"

python3 "${SCRIPT_DIR}/metrics/compute_metrics.py" \
    --experiment-dir "${EXPERIMENT_OUTPUT}" \
    --output "${METRICS_BUNDLE}"

# ============================================================================
# Step 4: Compare to Baseline
# ============================================================================

COMPARISON_FILE="${EXPERIMENT_OUTPUT}/comparison.json"
REGRESSION_FILE="${EXPERIMENT_OUTPUT}/regression_check.json"
VERDICT="unknown"
VERDICT_REASON=""

if [[ "${SKIP_COMPARE}" == "true" ]]; then
    log_validate "Step 4: Skipped (--skip-compare)"
elif ! has_baseline; then
    log_validate "Step 4: No baseline found — skipping comparison"
    log_validate "  Run: python3 scripts/validation/baseline/snapshot_baseline.py"
    VERDICT="no_baseline"
    VERDICT_REASON="No baseline to compare against"
else
    log_validate "Step 4: Comparing to baseline..."

    # Determine canonical mode name for comparison
    CANONICAL_MODE="$(canonical_mode_name "${PIPELINE_MODE}")"

    python3 "${SCRIPT_DIR}/metrics/compare_metrics.py" \
        --baseline "${BASELINE_STORE}" \
        --experiment "${METRICS_BUNDLE}" \
        --mode "${CANONICAL_MODE}" \
        --validation-mode "${VALIDATION_MODE}" \
        --output "${COMPARISON_FILE}"

    # Extract verdict from comparison
    VERDICT=$(python3 -c "import json; print(json.load(open('${COMPARISON_FILE}'))['verdict'])" 2>/dev/null || echo "error")
    VERDICT_REASON=$(python3 -c "import json; print(json.load(open('${COMPARISON_FILE}')).get('verdict_reason', ''))" 2>/dev/null || echo "")

    # Step 4b: Regression check
    log_validate "Step 4b: Regression check..."
    python3 "${SCRIPT_DIR}/metrics/regression_check.py" \
        --comparison "${COMPARISON_FILE}" \
        --threshold "${F1_REGRESSION_THRESHOLD}" \
        --output "${REGRESSION_FILE}" || true
fi

# ============================================================================
# Step 5: Generate Report
# ============================================================================

log_validate "Step 5: Generating report..."
REPORT_FILE="${EXPERIMENT_OUTPUT}/experiment_report.md"

python3 "${SCRIPT_DIR}/reporting/generate_report.py" \
    --experiment-id "${EXPERIMENT_ID}" \
    --hypothesis "${HYPOTHESIS}" \
    --hypothesis-id "${HYPOTHESIS_ID}" \
    --git-commit "${COMMIT}" \
    --git-branch "${BRANCH}" \
    --git-dirty "${DIRTY}" \
    --validation-mode "${VALIDATION_MODE}" \
    --verdict "${VERDICT}" \
    --verdict-reason "${VERDICT_REASON}" \
    --metrics-bundle "${METRICS_BUNDLE}" \
    --comparison "${COMPARISON_FILE}" \
    --output "${REPORT_FILE}" 2>/dev/null || log_validate "Report generation failed (non-critical)"

# ============================================================================
# Step 6: Update Evidence Ledger
# ============================================================================

log_validate "Step 6: Updating evidence ledger..."

WALL_END=$(date +%s)
WALL_TIME=$(( WALL_END - WALL_START ))

python3 "${SCRIPT_DIR}/ledger/experiment_registry.py" \
    --experiment-id "${EXPERIMENT_ID}" \
    --hypothesis "${HYPOTHESIS}" \
    --hypothesis-id "${HYPOTHESIS_ID}" \
    --git-commit "${COMMIT}" \
    --git-branch "${BRANCH}" \
    --validation-mode "${VALIDATION_MODE}" \
    --verdict "${VERDICT}" \
    --wall-time "${WALL_TIME}" \
    --metrics-bundle "${METRICS_BUNDLE}" \
    --comparison "${COMPARISON_FILE}" \
    --ledger "${EVIDENCE_LEDGER}" 2>/dev/null || log_validate "Ledger update failed (non-critical)"

# ============================================================================
# Final Summary
# ============================================================================

echo ""
echo "================================================================"
echo "  EXPERIMENT COMPLETE: ${EXPERIMENT_ID}"
echo "  Hypothesis: ${HYPOTHESIS}"
echo "  Mode: ${VALIDATION_MODE} | Track: ${TRACK} | Pipeline: ${PIPELINE_MODE}"
echo "  Samples: ${SAMPLES_TO_RUN[*]}"
echo "  Git: ${BRANCH}@${COMMIT}"
echo "  Wall time: ${WALL_TIME}s"
echo ""
echo "  VERDICT: ${VERDICT}"
echo "  Reason:  ${VERDICT_REASON}"
echo ""
echo "  Output: ${EXPERIMENT_OUTPUT}"
echo "  Report: ${REPORT_FILE}"
echo "================================================================"
