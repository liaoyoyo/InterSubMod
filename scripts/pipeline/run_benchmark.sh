#!/bin/bash
# run_benchmark.sh - Main benchmark pipeline orchestrator
#
# Runs the complete methylation filter validation pipeline:
#   00_prepare_germline  → Resolve germline phased VCF
#   01_longphase_s       → Haplotagging + TP/FP VCF generation
#   02_intersubmod       → Methylation analysis (TP + FP)
#   03_filter_analysis   → Filter condition analysis + F1 calculation
#   04_cleanup           → Remove large intermediate files
#
# Usage:
#   ./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395
#   ./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395 --dry-run
#   ./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395 --skip-longphase
#   ./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395 --skip-cleanup

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

MODE="s-pure"
SAMPLE="HCC1395"
DRY_RUN=false
SKIP_LONGPHASE=false
SKIP_INTERSUBMOD=false
SKIP_FILTER_ANALYSIS=false
SKIP_CLEANUP=false
COMPRESS_REGIONS=false
LOCAL_THREADS="${THREADS}"
LOCAL_INDEX_THREADS=""
MIN_FREE_GB="${MIN_FREE_GB_DEFAULT}"
VCF_SOURCE="output"
RUN_TAG=""
FORCE_RUN_METHYLATION=false
AUTO_SKIPPED_METHYLATION=false
METH_MM_TAGS="NA"
METH_ML_TAGS="NA"

# Override paths (for --skip-longphase with existing outputs)
EXISTING_TAGGED_BAM=""
EXISTING_TP_VCF=""
EXISTING_FP_VCF=""

show_help() {
    cat <<'HELP'
Usage: run_benchmark.sh [OPTIONS]

Options:
  --mode MODE           Verification mode (default: s-pure)
  --sample SAMPLE       Sample name (default: HCC1395)
  --threads N           Override thread count
  --index-threads N     Threads used for BAM indexing / sort-heavy I/O
  --min-free-gb N       Minimum free space required on output volume before run
  --output-root DIR     Override canonical output root
  --run-tag TEXT        Optional suffix appended to canonical run id
  --dry-run             Print commands without executing
  --skip-longphase      Skip LongPhase-S (use existing tagged BAM + TP/FP VCFs)
  --skip-intersubmod    Skip InterSubMod (use existing significance_summary.csv)
  --force-run-methylation  Force InterSubMod even when MM/ML tags are missing
  --skip-cleanup        Keep tagged BAM after completion
  --compress-regions    Compress per-region output directories
  --tagged-bam PATH     Path to existing tagged BAM (with --skip-longphase)
  --tp-vcf PATH         Path to existing TP VCF (with --skip-longphase)
  --fp-vcf PATH         Path to existing FP VCF (with --skip-longphase)
  -h, --help            Show this help
HELP
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)             MODE="$2"; shift 2 ;;
        --sample)           SAMPLE="$2"; shift 2 ;;
        --threads)          LOCAL_THREADS="$2"; shift 2 ;;
        --index-threads)    LOCAL_INDEX_THREADS="$2"; shift 2 ;;
        --min-free-gb)      MIN_FREE_GB="$2"; shift 2 ;;
        --output-root)      OUTPUT_ROOT="$2"; shift 2 ;;
        --run-tag)          RUN_TAG="$2"; shift 2 ;;
        --dry-run)          DRY_RUN=true; shift ;;
        --skip-longphase)   SKIP_LONGPHASE=true; shift ;;
        --skip-intersubmod) SKIP_INTERSUBMOD=true; shift ;;
        --force-run-methylation) FORCE_RUN_METHYLATION=true; shift ;;
        --skip-cleanup)     SKIP_CLEANUP=true; shift ;;
        --compress-regions) COMPRESS_REGIONS=true; shift ;;
        --tagged-bam)       EXISTING_TAGGED_BAM="$2"; shift 2 ;;
        --tp-vcf)           EXISTING_TP_VCF="$2"; shift 2 ;;
        --fp-vcf)           EXISTING_FP_VCF="$2"; shift 2 ;;
        --vcf-source)       VCF_SOURCE="$2"; shift 2 ;;
        -h|--help)          show_help ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ============================================================================
# Setup
# ============================================================================

# Load sample configuration
eval "$(get_sample_config "${SAMPLE}")"

if [[ -z "${LOCAL_INDEX_THREADS}" ]]; then
    LOCAL_INDEX_THREADS="$(recommended_index_threads "${LOCAL_THREADS}")"
fi

# Create output directory
DATE_STR=$(date +%Y%m%d)
CANONICAL_MODE="$(canonical_mode_name "${MODE}")"
CALLER_MODEL="$(caller_model_name "${VCF_SOURCE}")"
SAMPLE_OUTPUT_DIR="$(build_canonical_run_base "${SAMPLE}" "${CANONICAL_MODE}" "${CALLER_MODEL}" "${DATE_STR}" "${RUN_TAG}")"

# Avoid overwriting existing output
SEQ=1
FINAL_OUTPUT_DIR="${SAMPLE_OUTPUT_DIR}"
while [[ -d "${FINAL_OUTPUT_DIR}" ]]; do
    FINAL_OUTPUT_DIR="${SAMPLE_OUTPUT_DIR}_${SEQ}"
    SEQ=$((SEQ + 1))
done
SAMPLE_OUTPUT_DIR="${FINAL_OUTPUT_DIR}"

mkdir -p "${SAMPLE_OUTPUT_DIR}"

# Setup logging (capture both stdout and stderr to log file)
RUN_LOG="${SAMPLE_OUTPUT_DIR}/run.log"
exec > >(tee -a "${RUN_LOG}") 2> >(tee -a "${RUN_LOG}" >&2)

# ============================================================================
# Banner
# ============================================================================

echo "=================================================================" >&2
echo "  Benchmark Pipeline - Methylation Filter Validation" >&2
echo "=================================================================" >&2
log_info "Sample:      ${SAMPLE}"
log_info "Mode:        ${MODE}"
log_info "Canonical:   ${CANONICAL_MODE}"
log_info "Caller mdl:  ${CALLER_MODEL}"
log_info "Threads:     ${LOCAL_THREADS}"
log_info "Index thr:   ${LOCAL_INDEX_THREADS}"
log_info "Min free GB: ${MIN_FREE_GB}"
log_info "Output dir:  ${SAMPLE_OUTPUT_DIR}"
log_info "Dry-run:     ${DRY_RUN}"
log_info "Skip LongPhase-S:  ${SKIP_LONGPHASE}"
log_info "Skip InterSubMod:  ${SKIP_INTERSUBMOD}"
log_info "Force methylation: ${FORCE_RUN_METHYLATION}"
log_info "Skip cleanup:      ${SKIP_CLEANUP}"
echo "-----------------------------------------------------------------" >&2
print_disk_space
if [[ "${DRY_RUN}" != true ]]; then
    require_min_free_gb "${SAMPLE_OUTPUT_DIR}" "${MIN_FREE_GB}" "run_benchmark:${SAMPLE}:${CANONICAL_MODE}"
fi

PIPELINE_START=$(date +%s)

# ============================================================================
# Step 1: LongPhase-S (or use existing outputs)
# ============================================================================

LONGPHASE_DIR="${SAMPLE_OUTPUT_DIR}/longphase_s"

if [[ "${SKIP_LONGPHASE}" == true ]]; then
    log_info "[Pipeline] Skipping LongPhase-S (using existing outputs)"

    # Resolve paths
    if [[ -n "${EXISTING_TAGGED_BAM}" ]]; then
        TAGGED_BAM="${EXISTING_TAGGED_BAM}"
    else
        TAGGED_BAM="$(find_canonical_artifact "${SAMPLE}" "${CANONICAL_MODE}" "longphase_s/${SAMPLE}_tagged.bam" || true)"
        if [[ -z "${TAGGED_BAM}" ]] && [[ "${SAMPLE}" == "HCC1395" ]]; then
            TAGGED_BAM="/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam"
        fi
    fi

    if [[ -n "${EXISTING_TP_VCF}" ]]; then
        TP_VCF="${EXISTING_TP_VCF}"
    else
        TP_VCF="$(find_canonical_artifact "${SAMPLE}" "${CANONICAL_MODE}" "longphase_s/filtered_snv_tp.vcf.gz" || true)"
        if [[ -z "${TP_VCF}" ]] && [[ "${SAMPLE}" == "HCC1395" ]]; then
            TP_VCF="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
        fi
    fi

    if [[ -n "${EXISTING_FP_VCF}" ]]; then
        FP_VCF="${EXISTING_FP_VCF}"
    else
        FP_VCF="$(find_canonical_artifact "${SAMPLE}" "${CANONICAL_MODE}" "longphase_s/filtered_snv_fp.vcf.gz" || true)"
        if [[ -z "${FP_VCF}" ]] && [[ "${SAMPLE}" == "HCC1395" ]]; then
            FP_VCF="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz"
        fi
    fi

    validate_file "${TAGGED_BAM}" "Existing tagged BAM"
    validate_file "${TP_VCF}" "Existing TP VCF"
    validate_file "${FP_VCF}" "Existing FP VCF"

    TP_COUNT=$(zgrep -c -v "^#" "${TP_VCF}" || echo 0)
    FP_COUNT=$(zgrep -c -v "^#" "${FP_VCF}" || echo 0)
    log_info "  Using existing TP VCF: ${TP_VCF} (${TP_COUNT} variants)"
    log_info "  Using existing FP VCF: ${FP_VCF} (${FP_COUNT} variants)"
else
    LONGPHASE_ARGS=(
        --sample "${SAMPLE}"
        --output-dir "${SAMPLE_OUTPUT_DIR}"
        --threads "${LOCAL_THREADS}"
        --index-threads "${LOCAL_INDEX_THREADS}"
    )
    if [[ "${DRY_RUN}" == true ]]; then
        LONGPHASE_ARGS+=(--dry-run)
    fi

    # Resolve VCF source
    TARGET_VCF=""
    if [[ "${VCF_SOURCE}" == "pileup" ]]; then
        TARGET_VCF="${SOMATIC_VCF_PILEUP}"
        log_info "Using VCF Source: PILEUP (${TARGET_VCF})"
    elif [[ "${VCF_SOURCE}" == "indel" ]]; then
        TARGET_VCF="${SOMATIC_VCF_INDEL}"
        log_info "Using VCF Source: INDEL (${TARGET_VCF})"
    else
        TARGET_VCF="${SOMATIC_VCF}"
        log_info "Using VCF Source: OUTPUT (${TARGET_VCF})"
    fi

    # Pass specific VCF to LongPhase-S
    LONGPHASE_ARGS+=(--somatic-vcf "${TARGET_VCF}")

    "${SCRIPT_DIR}/steps/01_longphase_s.sh" "${LONGPHASE_ARGS[@]}"

    TAGGED_BAM="${LONGPHASE_DIR}/${SAMPLE}_tagged.bam"
    TP_VCF="${LONGPHASE_DIR}/filtered_snv_tp.vcf.gz"
    FP_VCF="${LONGPHASE_DIR}/filtered_snv_fp.vcf.gz"

    if [[ "${DRY_RUN}" != true ]]; then
        # Read counts from the variant_counts.txt generated by step 01
        if [[ -f "${LONGPHASE_DIR}/variant_counts.txt" ]]; then
            source "${LONGPHASE_DIR}/variant_counts.txt"
        else
            TP_COUNT=$(zgrep -c -v "^#" "${TP_VCF}" || echo 0)
            FP_COUNT=$(zgrep -c -v "^#" "${FP_VCF}" || echo 0)
        fi
    fi
fi

# ============================================================================
# Methylation Guard: Check MM/ML tags before InterSubMod
# ============================================================================

if [[ "${DRY_RUN}" != true ]] && [[ "${FORCE_RUN_METHYLATION}" != true ]] && [[ "${SKIP_INTERSUBMOD}" != true ]]; then
    validate_file "${TAGGED_BAM}" "Tagged BAM"
    METH_TAG_COUNTS="$(count_mm_ml_tags "${TAGGED_BAM}" 1000)"
    METH_MM_TAGS="${METH_TAG_COUNTS%,*}"
    METH_ML_TAGS="${METH_TAG_COUNTS#*,}"

    if [[ "${METH_MM_TAGS}" -eq 0 ]] || [[ "${METH_ML_TAGS}" -eq 0 ]]; then
        AUTO_SKIPPED_METHYLATION=true
        SKIP_INTERSUBMOD=true
        SKIP_FILTER_ANALYSIS=true
        log_warn "[Guard] Tagged BAM lacks MM/ML tags (MM=${METH_MM_TAGS}, ML=${METH_ML_TAGS})."
        log_warn "[Guard] Auto-skip InterSubMod + filter analysis to avoid invalid methylation conclusions."

        cat > "${SAMPLE_OUTPUT_DIR}/methylation_guard_status.tsv" <<EOF
sample	tagged_bam	mm_tags	ml_tags	auto_skip_intersubmod	auto_skip_filter	reason
${SAMPLE}	${TAGGED_BAM}	${METH_MM_TAGS}	${METH_ML_TAGS}	true	true	missing_mm_ml_tags
EOF
    else
        log_info "[Guard] MM/ML tags detected (MM=${METH_MM_TAGS}, ML=${METH_ML_TAGS}); methylation steps enabled."
    fi
fi

# ============================================================================
# Step 2: InterSubMod Methylation Analysis
# ============================================================================

if [[ "${SKIP_INTERSUBMOD}" == true ]]; then
    if [[ "${AUTO_SKIPPED_METHYLATION}" == true ]]; then
        log_info "[Pipeline] Skipping InterSubMod (auto-guard: missing MM/ML tags)"
    else
        log_info "[Pipeline] Skipping InterSubMod (using existing outputs)"
    fi
else
    ISM_ARGS=(
        --tagged-bam "${TAGGED_BAM}"
        --normal-bam "${NORMAL_BAM}"
        --tp-vcf "${TP_VCF}"
        --fp-vcf "${FP_VCF}"
        --output-dir "${SAMPLE_OUTPUT_DIR}"
        --threads "${LOCAL_THREADS}"
    )
    if [[ "${DRY_RUN}" == true ]]; then
        ISM_ARGS+=(--dry-run)
    fi

    "${SCRIPT_DIR}/steps/02_intersubmod.sh" "${ISM_ARGS[@]}"
fi

# ============================================================================
# Step 3: Filter Analysis + F1 Calculation
# ============================================================================

TP_SUMMARY="${SAMPLE_OUTPUT_DIR}/intersubmod_tp/significance_summary.csv"
FP_SUMMARY="${SAMPLE_OUTPUT_DIR}/intersubmod_fp/significance_summary.csv"

if [[ "${SKIP_FILTER_ANALYSIS}" != true ]] && [[ "${SKIP_INTERSUBMOD}" == true ]]; then
    if [[ ! -f "${TP_SUMMARY}" ]] || [[ ! -f "${FP_SUMMARY}" ]]; then
        SKIP_FILTER_ANALYSIS=true
        log_warn "[Pipeline] InterSubMod summaries not found; skipping filter analysis."
    fi
fi

if [[ "${SKIP_FILTER_ANALYSIS}" == true ]]; then
    log_info "[Pipeline] Skipping filter analysis."
elif [[ "${DRY_RUN}" == true ]]; then
    log_info "[DRY-RUN] Filter analysis:"
    echo "  python3 ${SCRIPT_DIR}/steps/03_filter_analysis.py \\" >&2
    echo "    --tp-summary ${TP_SUMMARY} \\" >&2
    echo "    --fp-summary ${FP_SUMMARY} \\" >&2
    echo "    --tp-count <TP_COUNT> --fp-count <FP_COUNT> --truth-total ${TRUTH_TOTAL} \\" >&2
    echo "    --output-dir ${SAMPLE_OUTPUT_DIR} \\" >&2
    echo "    --sample ${SAMPLE} --mode ${MODE}" >&2
else
    log_info "[Pipeline] Running filter analysis..."
    python3 "${SCRIPT_DIR}/steps/03_filter_analysis.py" \
        --tp-summary "${TP_SUMMARY}" \
        --fp-summary "${FP_SUMMARY}" \
        --tp-count "${TP_COUNT}" \
        --fp-count "${FP_COUNT}" \
        --truth-total "${TRUTH_TOTAL}" \
        --tp-vcf-file "${TP_VCF}" \
        --fp-vcf-file "${FP_VCF}" \
        --output-dir "${SAMPLE_OUTPUT_DIR}" \
        --sample "${SAMPLE}" \
        --mode "${MODE}"
fi

# ============================================================================
# Step 4: Cleanup
# ============================================================================

if [[ "${SKIP_CLEANUP}" == true ]] || [[ "${SKIP_LONGPHASE}" == true ]]; then
    log_info "[Pipeline] Skipping cleanup (--skip-cleanup or --skip-longphase)"
else
    CLEANUP_ARGS=(--output-dir "${SAMPLE_OUTPUT_DIR}")
    if [[ "${DRY_RUN}" == true ]]; then
        CLEANUP_ARGS+=(--dry-run)
    fi
    if [[ "${COMPRESS_REGIONS}" == true ]]; then
        CLEANUP_ARGS+=(--compress-regions)
    fi

    "${SCRIPT_DIR}/steps/04_cleanup.sh" "${CLEANUP_ARGS[@]}"
fi

# ============================================================================
# Summary
# ============================================================================

PIPELINE_END=$(date +%s)
PIPELINE_ELAPSED=$((PIPELINE_END - PIPELINE_START))
PIPELINE_MIN=$((PIPELINE_ELAPSED / 60))
PIPELINE_SEC=$((PIPELINE_ELAPSED % 60))

echo "" >&2
echo "=================================================================" >&2
echo "  Pipeline Complete" >&2
echo "=================================================================" >&2
log_info "Total time: ${PIPELINE_MIN}m ${PIPELINE_SEC}s"
log_info "Output:     ${SAMPLE_OUTPUT_DIR}"
log_info "Log:        ${RUN_LOG}"

if [[ "${DRY_RUN}" != true ]] && [[ -f "${SAMPLE_OUTPUT_DIR}/metrics.json" ]]; then
    log_info "Metrics:"
    cat "${SAMPLE_OUTPUT_DIR}/metrics.json" >&2
fi

echo "=================================================================" >&2
