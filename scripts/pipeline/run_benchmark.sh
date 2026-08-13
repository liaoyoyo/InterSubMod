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
#   ./scripts/pipeline/run_benchmark.sh --site-profile config/site-profile.local.json --plan
#   ./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395 --dry-run
#   ./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395 --skip-longphase
#   ./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395 --skip-cleanup

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
ORIGINAL_ARGV=("$@")

# ============================================================================
# Argument Parsing
# ============================================================================

MODE="s-pure"
SAMPLE="HCC1395"
DRY_RUN=false
PLAN_ONLY=false
SITE_PROFILE=""
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
CLI_OUTPUT_ROOT=""
FORCE_RUN_METHYLATION=false
AUTO_SKIPPED_METHYLATION=false
METH_MM_TAGS="NA"
METH_ML_TAGS="NA"

# Override paths (for --skip-longphase with existing outputs)
EXISTING_TAGGED_BAM=""
EXISTING_TP_VCF=""
EXISTING_FP_VCF=""
EXISTING_TAGGED_BAM_SHA256=""
EXISTING_TP_VCF_SHA256=""
EXISTING_FP_VCF_SHA256=""
THREADS_PROVENANCE="PROFILE_OR_LEGACY_DEFAULT"
INDEX_THREADS_PROVENANCE="DERIVED_FROM_THREADS"
MIN_FREE_GB_PROVENANCE="DEFAULT"
OUTPUT_ROOT_PROVENANCE="PROFILE_OR_LEGACY_DEFAULT"
RUN_TAG_PROVENANCE="DEFAULT_EMPTY"

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
  --site-profile FILE   Load machine paths from a validated site profile
  --run-tag TEXT        Optional suffix appended to canonical run id
  --plan                Print the effective read-only execution plan and exit
  --dry-run             Print science commands; still writes run/provenance/log scaffolding
  --skip-longphase      Skip LongPhase-S (use existing tagged BAM + TP/FP VCFs)
  --skip-intersubmod    Skip InterSubMod (use existing significance_summary.csv)
  --force-run-methylation  Force InterSubMod even when MM/ML tags are missing
  --skip-cleanup        Keep tagged BAM after completion
  --compress-regions    Compress per-region output directories
  --tagged-bam PATH     Path to existing tagged BAM (with --skip-longphase)
  --tp-vcf PATH         Path to existing TP VCF (with --skip-longphase)
  --fp-vcf PATH         Path to existing FP VCF (with --skip-longphase)
  --tagged-bam-sha256 SHA256  Expected SHA-256 for --tagged-bam (profile mode)
  --tp-vcf-sha256 SHA256      Expected SHA-256 for --tp-vcf (profile mode)
  --fp-vcf-sha256 SHA256      Expected SHA-256 for --fp-vcf (profile mode)
  -h, --help            Show this help
HELP
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)             MODE="$2"; shift 2 ;;
        --sample)           SAMPLE="$2"; shift 2 ;;
        --threads)          LOCAL_THREADS="$2"; THREADS_PROVENANCE="CLI_OVERRIDE"; shift 2 ;;
        --index-threads)    LOCAL_INDEX_THREADS="$2"; INDEX_THREADS_PROVENANCE="CLI_OVERRIDE"; shift 2 ;;
        --min-free-gb)      MIN_FREE_GB="$2"; MIN_FREE_GB_PROVENANCE="CLI_OVERRIDE"; shift 2 ;;
        --output-root)      CLI_OUTPUT_ROOT="$2"; OUTPUT_ROOT_PROVENANCE="CLI_OVERRIDE"; shift 2 ;;
        --run-tag)          RUN_TAG="$2"; RUN_TAG_PROVENANCE="CLI_OVERRIDE"; shift 2 ;;
        --site-profile)     SITE_PROFILE="$2"; shift 2 ;;
        --plan)             PLAN_ONLY=true; shift ;;
        --dry-run)          DRY_RUN=true; shift ;;
        --skip-longphase)   SKIP_LONGPHASE=true; shift ;;
        --skip-intersubmod) SKIP_INTERSUBMOD=true; shift ;;
        --force-run-methylation) FORCE_RUN_METHYLATION=true; shift ;;
        --skip-cleanup)     SKIP_CLEANUP=true; shift ;;
        --compress-regions) COMPRESS_REGIONS=true; shift ;;
        --tagged-bam)       EXISTING_TAGGED_BAM="$2"; shift 2 ;;
        --tp-vcf)           EXISTING_TP_VCF="$2"; shift 2 ;;
        --fp-vcf)           EXISTING_FP_VCF="$2"; shift 2 ;;
        --tagged-bam-sha256) EXISTING_TAGGED_BAM_SHA256="$2"; shift 2 ;;
        --tp-vcf-sha256)     EXISTING_TP_VCF_SHA256="$2"; shift 2 ;;
        --fp-vcf-sha256)     EXISTING_FP_VCF_SHA256="$2"; shift 2 ;;
        --vcf-source)       VCF_SOURCE="$2"; shift 2 ;;
        -h|--help)          show_help ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

validate_benchmark_mode "${MODE}" || exit $?
validate_vcf_source "${VCF_SOURCE}" || exit $?
if ! is_safe_run_component "${SAMPLE}"; then
    echo "[ERROR] Unsafe --sample component: ${SAMPLE}" >&2
    exit 2
fi
if [[ -n "${RUN_TAG}" ]] && ! is_safe_run_component "${RUN_TAG}"; then
    echo "[ERROR] Unsafe --run-tag; use 1-64 ASCII letters, digits, dot, underscore, or hyphen." >&2
    exit 2
fi

# ============================================================================
# Setup
# ============================================================================

# A profile is a standalone authority. It must never inherit a missing role
# from get_sample_config or another machine's /big* defaults.
if [[ -n "${SITE_PROFILE}" ]]; then
    SITE_PROFILE="$(realpath -e -- "${SITE_PROFILE}")" || {
        echo "[ERROR] Site profile not found: ${SITE_PROFILE}" >&2
        exit 3
    }
    load_site_profile_config "${SITE_PROFILE}" "${SAMPLE}" || exit $?
else
    eval "$(get_sample_config "${SAMPLE}")"
fi
if [[ -n "${CLI_OUTPUT_ROOT}" ]]; then
    OUTPUT_ROOT="${CLI_OUTPUT_ROOT}"
fi
# The historical cleanup stage directly removes intermediates. Portable profile
# mode is archive-first: cleanup stays disabled until a separately reviewed
# archive target and receipt contract are supplied.
if [[ "${SITE_PROFILE_ACTIVE}" == true ]] && [[ "${SKIP_CLEANUP}" != true ]]; then
    SKIP_CLEANUP=true
    echo "[WARN] Portable profile mode forces --skip-cleanup; stage 04 is REPRODUCIBLE_LEGACY, not SUPPORTED." >&2
fi
if path_has_symlink_component "${OUTPUT_ROOT}"; then
    echo "[ERROR] Output root contains a symlink component: ${OUTPUT_ROOT}" >&2
    exit 3
fi
if [[ "${SITE_PROFILE_ACTIVE}" == true ]] && ! path_is_within "${OUTPUT_ROOT}" "${BIG7_OUTPUT_ROOT}"; then
    echo "[ERROR] Profile-mode output root escapes data_roots.output: ${OUTPUT_ROOT}" >&2
    exit 3
fi

TARGET_VCF="$(select_somatic_vcf "${MODE}" "${VCF_SOURCE}")" || exit $?
if [[ -z "${TARGET_VCF}" ]]; then
    echo "[ERROR] Selected VCF role resolved to an empty path." >&2
    exit 3
fi

# Portable skip mode never guesses a "latest" directory. The operator must
# bind all three existing artifacts by explicit path and expected SHA-256.
if [[ "${SITE_PROFILE_ACTIVE}" == true ]] && [[ "${SKIP_LONGPHASE}" == true ]]; then
    TAGGED_BAM="$(resolve_governed_existing_artifact "${EXISTING_TAGGED_BAM}" "${EXISTING_TAGGED_BAM_SHA256}" "Existing tagged BAM")" || exit $?
    TP_VCF="$(resolve_governed_existing_artifact "${EXISTING_TP_VCF}" "${EXISTING_TP_VCF_SHA256}" "Existing TP VCF")" || exit $?
    FP_VCF="$(resolve_governed_existing_artifact "${EXISTING_FP_VCF}" "${EXISTING_FP_VCF_SHA256}" "Existing FP VCF")" || exit $?
fi

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
if ! path_is_within "${SAMPLE_OUTPUT_DIR}" "${OUTPUT_ROOT}"; then
    echo "[ERROR] Computed run directory escapes output root: ${SAMPLE_OUTPUT_DIR}" >&2
    exit 3
fi

if [[ "${PLAN_ONLY}" == true ]]; then
    ENABLED_STAGES=()
    [[ "${SKIP_LONGPHASE}" != true ]] && ENABLED_STAGES+=(LongPhase-S)
    [[ "${SKIP_INTERSUBMOD}" != true ]] && ENABLED_STAGES+=(InterSubMod)
    if [[ "${SKIP_FILTER_ANALYSIS}" != true ]] && [[ "${SKIP_INTERSUBMOD}" != true ]]; then
        ENABLED_STAGES+=(filter-analysis)
    fi
    if [[ "${SKIP_CLEANUP}" != true ]] && [[ "${SKIP_LONGPHASE}" != true ]]; then
        ENABLED_STAGES+=(cleanup)
    fi
    STAGE_LIST="none"
    if (( ${#ENABLED_STAGES[@]} > 0 )); then
        STAGE_LIST="$(IFS=,; echo "${ENABLED_STAGES[*]}")"
    fi
    echo "[INPUT] site_profile=${SITE_PROFILE:-BUILT_IN_FALLBACK}"
    if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
        echo "[INPUT] site_profile_sha256=${SITE_PROFILE_SHA256}"
        echo "[PROVENANCE] profile_resolution=single_parent_parse"
    fi
    echo "[INPUT] sample=${SAMPLE} mode=${MODE} vcf_source=${VCF_SOURCE}"
    echo "[PROVENANCE] threads=${THREADS_PROVENANCE}:${LOCAL_THREADS}"
    echo "[PROVENANCE] index_threads=${INDEX_THREADS_PROVENANCE}:${LOCAL_INDEX_THREADS}"
    echo "[PROVENANCE] min_free_gb=${MIN_FREE_GB_PROVENANCE}:${MIN_FREE_GB}"
    echo "[PROVENANCE] output_root=${OUTPUT_ROOT_PROVENANCE}:${OUTPUT_ROOT}"
    echo "[PROVENANCE] run_tag=${RUN_TAG_PROVENANCE}:${RUN_TAG:-<empty>}"
    printf '[PROVENANCE] command='
    printf '%q ' "${0}" "${ORIGINAL_ARGV[@]}"
    printf '\n'
    echo "[PLAN] project_root=${PROJECT_ROOT_DEFAULT}"
    echo "[PLAN] reference=${REFERENCE}"
    echo "[PLAN] tumor_bam=${TUMOR_BAM}"
    echo "[PLAN] normal_bam=${NORMAL_BAM}"
    echo "[PLAN] selected_vcf=${TARGET_VCF}"
    if [[ "${SKIP_LONGPHASE}" == true ]]; then
        echo "[PLAN] existing_tagged_bam=${TAGGED_BAM:-AUTO_DISCOVERY_LEGACY}"
        echo "[PLAN] existing_tp_vcf=${TP_VCF:-AUTO_DISCOVERY_LEGACY}"
        echo "[PLAN] existing_fp_vcf=${FP_VCF:-AUTO_DISCOVERY_LEGACY}"
        if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
            echo "[PROVENANCE] existing_tagged_bam=CLI_PATH+EXPECTED_SHA256:${EXISTING_TAGGED_BAM_SHA256,,}"
            echo "[PROVENANCE] existing_tp_vcf=CLI_PATH+EXPECTED_SHA256:${EXISTING_TP_VCF_SHA256,,}"
            echo "[PROVENANCE] existing_fp_vcf=CLI_PATH+EXPECTED_SHA256:${EXISTING_FP_VCF_SHA256,,}"
        fi
    fi
    echo "[PLAN] longphase_s=${LONGPHASE_S_BIN}"
    echo "[PLAN] intersubmod=${INTERSUBMOD_BIN}"
    echo "[PLAN] output_dir=${SAMPLE_OUTPUT_DIR}"
    echo "[PLAN] stages=${STAGE_LIST}"
    echo "[RESULT] plan_only=true side_effects=none"
    exit 0
fi

# Gate real portable execution before creating any run output. The preflight
# command is read-only, emits one JSON document, and scopes real-data checks to
# the selected sample/dataset.
SITE_PREFLIGHT_JSON=""
if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
    CURRENT_PROFILE_SHA256="$(sha256sum -- "${SITE_PROFILE}" | awk '{print $1}')"
    if [[ "${CURRENT_PROFILE_SHA256}" != "${SITE_PROFILE_SHA256}" ]]; then
        echo "[ERROR] Site profile changed after it was parsed; no run output was created." >&2
        exit 3
    fi
    SITE_PREFLIGHT_JSON="${SITE_PROFILE_PREFLIGHT_JSON}"
    SITE_PREFLIGHT_PASS="$(python3 -c 'import json,sys; print(str(json.load(sys.stdin)["pass"]).lower())' <<< "${SITE_PREFLIGHT_JSON}")" || exit 3
    if [[ "${SITE_PREFLIGHT_PASS}" != true ]]; then
        echo "[ERROR] Real-data site preflight failed (exit 5); no run output was created." >&2
        echo "${SITE_PREFLIGHT_JSON}" >&2
        exit 5
    fi
fi

mkdir -p "${SAMPLE_OUTPUT_DIR}"
RUN_SENTINEL="${SAMPLE_OUTPUT_DIR}/.intersubmod-run-root"
printf 'schema_version=1\nrun_root=%s\nallowed_root=%s\n' \
    "$(realpath -e -- "${SAMPLE_OUTPUT_DIR}")" "$(realpath -m -- "${OUTPUT_ROOT}")" > "${RUN_SENTINEL}"

# Freeze the exact profile used for this run and fail closed on the complete
# real-data doctor before any scientific stage is launched.
if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
    PROFILE_PROVENANCE_DIR="${SAMPLE_OUTPUT_DIR}/provenance"
    mkdir -p "${PROFILE_PROVENANCE_DIR}"
    {
        printf 'command='
        printf '%q ' "${0}" "${ORIGINAL_ARGV[@]}"
        printf '\nthreads=%s\nindex_threads=%s\nmin_free_gb=%s\noutput_root=%s\nrun_tag=%s\n' \
            "${THREADS_PROVENANCE}:${LOCAL_THREADS}" \
            "${INDEX_THREADS_PROVENANCE}:${LOCAL_INDEX_THREADS}" \
            "${MIN_FREE_GB_PROVENANCE}:${MIN_FREE_GB}" \
            "${OUTPUT_ROOT_PROVENANCE}:${OUTPUT_ROOT}" \
            "${RUN_TAG_PROVENANCE}:${RUN_TAG:-<empty>}"
    } > "${PROFILE_PROVENANCE_DIR}/command_receipt.txt"
    SITE_PROFILE_LOCK_PATH="${PROFILE_PROVENANCE_DIR}/site-profile.locked.json"
    install -m 0444 -- "${SITE_PROFILE}" "${SITE_PROFILE_LOCK_PATH}"
    LOCKED_PROFILE_SHA256="$(sha256sum -- "${SITE_PROFILE_LOCK_PATH}" | awk '{print $1}')"
    if [[ "${LOCKED_PROFILE_SHA256}" != "${SITE_PROFILE_SHA256}" ]]; then
        echo "[ERROR] Run-local profile copy failed SHA-256 verification." >&2
        exit 3
    fi
    printf 'schema_version\tprofile_source\tprofile_sha256\tlocked_profile\tsample\n1\t%s\t%s\t%s\t%s\n' \
        "${SITE_PROFILE}" "${SITE_PROFILE_SHA256}" "${SITE_PROFILE_LOCK_PATH}" "${SAMPLE}" \
        > "${PROFILE_PROVENANCE_DIR}/profile_lock_receipt.tsv"
    printf '%s\n' "${SITE_PREFLIGHT_JSON}" > "${PROFILE_PROVENANCE_DIR}/real_data_preflight_receipt.json"
    SITE_PROFILE_PARENT_LOCKED=true
    export SITE_PROFILE_PARENT_LOCKED SITE_PROFILE_LOCK_PATH SITE_PROFILE_SHA256

    if [[ "${SKIP_LONGPHASE}" == true ]]; then
        verify_artifact_sha256 "${TAGGED_BAM}" "${EXISTING_TAGGED_BAM_SHA256}" "Existing tagged BAM" || exit $?
        verify_artifact_sha256 "${TP_VCF}" "${EXISTING_TP_VCF_SHA256}" "Existing TP VCF" || exit $?
        verify_artifact_sha256 "${FP_VCF}" "${EXISTING_FP_VCF_SHA256}" "Existing FP VCF" || exit $?
    fi
fi

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
if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
    log_info "Profile SHA:       ${SITE_PROFILE_SHA256}"
    log_info "Profile lock:      ${SITE_PROFILE_LOCK_PATH}"
    log_info "Profile loading:   single parent parse; child steps verify locked SHA"
fi
echo "-----------------------------------------------------------------" >&2
print_disk_space "${SAMPLE_OUTPUT_DIR}"
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

    # Legacy mode retains historical discovery for compatibility. Profile mode
    # was already resolved and SHA-bound before planning or output creation.
    if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
        :
    elif [[ -n "${EXISTING_TAGGED_BAM}" ]]; then
        TAGGED_BAM="${EXISTING_TAGGED_BAM}"
    else
        TAGGED_BAM="$(find_canonical_artifact "${SAMPLE}" "${CANONICAL_MODE}" "longphase_s/${SAMPLE}_tagged.bam" || true)"
        if [[ "${SITE_PROFILE_ACTIVE}" != true ]] && [[ -z "${TAGGED_BAM}" ]] && [[ "${SAMPLE}" == "HCC1395" ]]; then
            TAGGED_BAM="/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam"
        fi
    fi

    if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
        :
    elif [[ -n "${EXISTING_TP_VCF}" ]]; then
        TP_VCF="${EXISTING_TP_VCF}"
    else
        TP_VCF="$(find_canonical_artifact "${SAMPLE}" "${CANONICAL_MODE}" "longphase_s/filtered_snv_tp.vcf.gz" || true)"
        if [[ "${SITE_PROFILE_ACTIVE}" != true ]] && [[ -z "${TP_VCF}" ]] && [[ "${SAMPLE}" == "HCC1395" ]]; then
            TP_VCF="/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz"
        fi
    fi

    if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
        :
    elif [[ -n "${EXISTING_FP_VCF}" ]]; then
        FP_VCF="${EXISTING_FP_VCF}"
    else
        FP_VCF="$(find_canonical_artifact "${SAMPLE}" "${CANONICAL_MODE}" "longphase_s/filtered_snv_fp.vcf.gz" || true)"
        if [[ "${SITE_PROFILE_ACTIVE}" != true ]] && [[ -z "${FP_VCF}" ]] && [[ "${SAMPLE}" == "HCC1395" ]]; then
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
        --somatic-vcf "${TARGET_VCF}"
    )
    if [[ "${DRY_RUN}" == true ]]; then
        LONGPHASE_ARGS+=(--dry-run)
    fi

    log_info "Using VCF source ${VCF_SOURCE}: ${TARGET_VCF}"

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
    CLEANUP_ARGS=(--output-dir "${SAMPLE_OUTPUT_DIR}" --allowed-root "${OUTPUT_ROOT}")
    if [[ -n "${SITE_PROFILE}" ]]; then
        CLEANUP_ARGS+=(--site-profile "${SITE_PROFILE}" --sample "${SAMPLE}")
    fi
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
