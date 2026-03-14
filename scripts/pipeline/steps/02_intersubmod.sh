#!/bin/bash
# 02_intersubmod.sh - Run InterSubMod methylation analysis on TP and FP VCFs
#
# Runs InterSubMod sequentially: first on TP VCF, then on FP VCF.
# Produces significance_summary.csv for each.
#
# Usage:
#   ./scripts/pipeline/steps/02_intersubmod.sh \
#     --tagged-bam BAM --normal-bam BAM --tp-vcf VCF --fp-vcf VCF \
#     --output-dir DIR [--threads N] [--dry-run]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

TAGGED_BAM=""
NORMAL_BAM_ARG=""
TP_VCF=""
FP_VCF=""
ISM_OUTPUT_DIR=""
LOCAL_THREADS="${THREADS}"
LOCAL_WINDOW_SIZE="${WINDOW_SIZE}"
LOCAL_DISTANCE_METRIC="${DISTANCE_METRIC}"
LOCAL_LOG_LEVEL="${LOG_LEVEL}"
LOCAL_PLOT_TYPE="${PLOT_TYPE}"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --tagged-bam)       TAGGED_BAM="$2"; shift 2 ;;
        --normal-bam)       NORMAL_BAM_ARG="$2"; shift 2 ;;
        --tp-vcf)           TP_VCF="$2"; shift 2 ;;
        --fp-vcf)           FP_VCF="$2"; shift 2 ;;
        --output-dir)       ISM_OUTPUT_DIR="$2"; shift 2 ;;
        --threads)          LOCAL_THREADS="$2"; shift 2 ;;
        --window-size)      LOCAL_WINDOW_SIZE="$2"; shift 2 ;;
        --distance-metric)  LOCAL_DISTANCE_METRIC="$2"; shift 2 ;;
        --log-level)        LOCAL_LOG_LEVEL="$2"; shift 2 ;;
        --plot-type)        LOCAL_PLOT_TYPE="$2"; shift 2 ;;
        --dry-run)          DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: $0 --tagged-bam BAM --normal-bam BAM --tp-vcf VCF --fp-vcf VCF --output-dir DIR [OPTIONS]" >&2
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${TAGGED_BAM}" ]] || [[ -z "${TP_VCF}" ]] || [[ -z "${FP_VCF}" ]] || [[ -z "${ISM_OUTPUT_DIR}" ]]; then
    log_error "--tagged-bam, --tp-vcf, --fp-vcf, and --output-dir are all required."
    exit 1
fi

# ============================================================================
# Resolve InterSubMod Binary
# ============================================================================

ISM_BIN=$(resolve_intersubmod_bin) || exit 1

# ============================================================================
# Input Validation
# ============================================================================

log_info "[Step 02] InterSubMod methylation analysis"
log_info "  Binary:          ${ISM_BIN}"
log_info "  Tagged BAM:      ${TAGGED_BAM}"
log_info "  Normal BAM:      ${NORMAL_BAM_ARG:-<not set>}"
log_info "  TP VCF:          ${TP_VCF}"
log_info "  FP VCF:          ${FP_VCF}"
log_info "  Window size:     ${LOCAL_WINDOW_SIZE}"
log_info "  Distance metric: ${LOCAL_DISTANCE_METRIC}"
log_info "  Threads:         ${LOCAL_THREADS}"
log_info "  Output dir:      ${ISM_OUTPUT_DIR}"

if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${TAGGED_BAM}" "Tagged BAM"
    validate_file "${TP_VCF}" "TP VCF"
    validate_file "${FP_VCF}" "FP VCF"
    validate_file "${REFERENCE}" "Reference FASTA"
fi

# ============================================================================
# Build Common Arguments (as array)
# ============================================================================
# Use ALL threads for each task (sequential execution maximizes resource usage)
log_info "  Running TP and FP sequentially (${LOCAL_THREADS} threads each, full CPU utilization)"

# Common args for InterSubMod (shared between TP and FP)
ISM_COMMON_ARGS=(
    --tumor-bam "${TAGGED_BAM}"
    --reference "${REFERENCE}"
    --threads "${LOCAL_THREADS}"
    --window-size "${LOCAL_WINDOW_SIZE}"
    --distance-metric "${LOCAL_DISTANCE_METRIC}"
    --log-level "${LOCAL_LOG_LEVEL}"
    --output-filtered-reads
)

# Add normal BAM if provided
if [[ -n "${NORMAL_BAM_ARG}" ]]; then
    ISM_COMMON_ARGS+=(--normal-bam "${NORMAL_BAM_ARG}")
fi

# ============================================================================
# Build TP/FP commands
# ============================================================================

TP_OUTPUT="${ISM_OUTPUT_DIR}/intersubmod_tp"
FP_OUTPUT="${ISM_OUTPUT_DIR}/intersubmod_fp"
mkdir -p "${TP_OUTPUT}" "${FP_OUTPUT}"

TP_CMD=(
    "${ISM_BIN}"
    "${ISM_COMMON_ARGS[@]}"
    --vcf "${TP_VCF}"
    --output-dir "${TP_OUTPUT}"
)

FP_CMD=(
    "${ISM_BIN}"
    "${ISM_COMMON_ARGS[@]}"
    --vcf "${FP_VCF}"
    --output-dir "${FP_OUTPUT}"
)

# ============================================================================
# Execute TP then FP sequentially (full CPU utilization per task)
# ============================================================================

if [[ "${DRY_RUN}" == true ]]; then
    log_info "[DRY-RUN] InterSubMod (TP):"
    echo "  ${TP_CMD[*]}" >&2
    log_info "[DRY-RUN] InterSubMod (FP):"
    echo "  ${FP_CMD[*]}" >&2
else
    TP_LOG="${ISM_OUTPUT_DIR}/intersubmod_tp.log"
    FP_LOG="${ISM_OUTPUT_DIR}/intersubmod_fp.log"

    # --- Run TP first (all threads) ---
    log_info "  Running InterSubMod on TP VCF (${LOCAL_THREADS} threads)..."
    TP_START=$(date +%s)
    /usr/bin/time -v "${TP_CMD[@]}" > "${TP_LOG}" 2>&1
    TP_EXIT=$?
    TP_END=$(date +%s)
    TP_ELAPSED=$(( TP_END - TP_START ))

    if [[ ${TP_EXIT} -ne 0 ]]; then
        log_error "  InterSubMod TP failed with exit code ${TP_EXIT}. See: ${TP_LOG}"
        exit 1
    fi
    if [[ -f "${TP_OUTPUT}/significance_summary.csv" ]]; then
        TP_REGIONS=$(tail -n +2 "${TP_OUTPUT}/significance_summary.csv" | wc -l)
        python3 "${SCRIPT_DIR}/../../analysis/extract_label_first_metrics.py" \
            --summary-csv "${TP_OUTPUT}/significance_summary.csv" \
            --output-tsv "${TP_OUTPUT}/label_first_metrics.tsv"
        log_info "  TP analysis complete: ${TP_REGIONS} regions (${TP_ELAPSED}s)"
    else
        log_error "  TP significance_summary.csv not found!"
        exit 1
    fi

    # --- Run FP next (all threads) ---
    log_info "  Running InterSubMod on FP VCF (${LOCAL_THREADS} threads)..."
    FP_START=$(date +%s)
    /usr/bin/time -v "${FP_CMD[@]}" > "${FP_LOG}" 2>&1
    FP_EXIT=$?
    FP_END=$(date +%s)
    FP_ELAPSED=$(( FP_END - FP_START ))

    if [[ ${FP_EXIT} -ne 0 ]]; then
        log_error "  InterSubMod FP failed with exit code ${FP_EXIT}. See: ${FP_LOG}"
        exit 1
    fi
    if [[ -f "${FP_OUTPUT}/significance_summary.csv" ]]; then
        FP_REGIONS=$(tail -n +2 "${FP_OUTPUT}/significance_summary.csv" | wc -l)
        python3 "${SCRIPT_DIR}/../../analysis/extract_label_first_metrics.py" \
            --summary-csv "${FP_OUTPUT}/significance_summary.csv" \
            --output-tsv "${FP_OUTPUT}/label_first_metrics.tsv"
        log_info "  FP analysis complete: ${FP_REGIONS} regions (${FP_ELAPSED}s)"
    else
        log_error "  FP significance_summary.csv not found!"
        exit 1
    fi
fi

log_info "[Step 02] InterSubMod analysis complete."
