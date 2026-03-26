#!/bin/bash
# Resume a paired canonical run after LongPhase-S using existing tagged BAM + _sc.vcf.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

SAMPLE=""
CANONICAL_MODE=""
RUN_DIR=""
LOCAL_THREADS="${THREADS}"
MIN_FREE_GB="${MIN_FREE_GB_DEFAULT}"
DRY_RUN=false
SKIP_METHOD_VALIDATION=false
SKIP_DASHBOARD=false

show_help() {
    cat <<'EOF'
Usage: resume_paired_run_from_tagged_outputs.sh [OPTIONS]

Options:
  --sample NAME            Sample name
  --canonical-mode MODE    paired_full or paired_pileup
  --run-dir DIR            Existing canonical run directory
  --threads N              Threads for InterSubMod
  --min-free-gb N          Minimum free disk on /big7_disk before resuming
  --skip-method-validation Skip validate_method_design.py
  --skip-dashboard         Skip build_round_dashboard.py
  --dry-run                Print commands without executing
  -h, --help               Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --canonical-mode) CANONICAL_MODE="$2"; shift 2 ;;
        --run-dir) RUN_DIR="$2"; shift 2 ;;
        --threads) LOCAL_THREADS="$2"; shift 2 ;;
        --min-free-gb) MIN_FREE_GB="$2"; shift 2 ;;
        --skip-method-validation) SKIP_METHOD_VALIDATION=true; shift ;;
        --skip-dashboard) SKIP_DASHBOARD=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) show_help; exit 0 ;;
        *)
            log_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "${SAMPLE}" || -z "${CANONICAL_MODE}" || -z "${RUN_DIR}" ]]; then
    log_error "--sample, --canonical-mode, and --run-dir are required."
    exit 1
fi

case "${CANONICAL_MODE}" in
    paired_full)
        PIPELINE_MODE="s-pure"
        ;;
    paired_pileup)
        PIPELINE_MODE="s-pure-pileup"
        ;;
    *)
        log_error "Unsupported canonical mode: ${CANONICAL_MODE}"
        exit 1
        ;;
esac

case "${SAMPLE}" in
    HCC1395)
        PLATFORM="ONT_5kHz"
        ;;
    HCC1395_DORADO)
        PLATFORM="ONT_Dorado"
        ;;
    *)
        PLATFORM="ONT"
        ;;
esac

case "${SAMPLE}" in
    HCC1395|HCC1395_DORADO)
        TRUTH_SET="SEQC2_v1.2.1_HC_SNV"
        ;;
    COLO829)
        TRUTH_SET="NYGC"
        ;;
    *)
        TRUTH_SET="orthogonal-tools"
        ;;
esac

run_cmd() {
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] $*" >&2
    else
        "$@"
    fi
}

eval "$(get_sample_config "${SAMPLE}")"

case "${CANONICAL_MODE}" in
    paired_full)
        CALLER_VCF_PATH="${SOMATIC_VCF}"
        ;;
    paired_pileup)
        CALLER_VCF_PATH="${SOMATIC_VCF_PILEUP}"
        ;;
esac

RUN_DIR="$(realpath "${RUN_DIR}")"
LONGPHASE_DIR="${RUN_DIR}/longphase_s"
TAGGED_BAM="${LONGPHASE_DIR}/${SAMPLE}_tagged.bam"
TAGGED_BAI="${LONGPHASE_DIR}/${SAMPLE}_tagged.bam.bai"
SC_VCF="${LONGPHASE_DIR}/${SAMPLE}_tagged_sc.vcf"
TP_VCF="${LONGPHASE_DIR}/filtered_snv_tp.vcf.gz"
FP_VCF="${LONGPHASE_DIR}/filtered_snv_fp.vcf.gz"
COUNTS_FILE="${LONGPHASE_DIR}/variant_counts.txt"
HAPLOTAG_QC_LONGPHASE="${LONGPHASE_DIR}/haplotag_qc.tsv"
HAPLOTAG_QC_ROOT="${RUN_DIR}/haplotag_qc.tsv"
RUN_CONTEXT="${RUN_DIR}/run_context.json"
ROUND_CONTEXT="${RUN_DIR}/round_context.json"
RESUME_LOG="${RUN_DIR}/resume_from_tagged.log"
TP_SUMMARY="${RUN_DIR}/intersubmod_tp/significance_summary.csv"
FP_SUMMARY="${RUN_DIR}/intersubmod_fp/significance_summary.csv"
METRICS_JSON="${RUN_DIR}/metrics.json"
BENCHMARK_TSV="${RUN_DIR}/benchmark_comparison.tsv"
BENCHMARK_MD="${RUN_DIR}/benchmark_comparison.md"
METHOD_VALIDATION_TSV="${RUN_DIR}/method_design_validation.tsv"
AGREEMENT_TSV="${RUN_DIR}/label_cluster_agreement.tsv"
ROUND_SUMMARY_MD="${RUN_DIR}/round_summary.md"

exec > >(tee -a "${RESUME_LOG}") 2> >(tee -a "${RESUME_LOG}" >&2)

log_info "[Resume paired] Sample=${SAMPLE} CanonicalMode=${CANONICAL_MODE}"
log_info "[Resume paired] Run dir=${RUN_DIR}"
log_info "[Resume paired] Threads=${LOCAL_THREADS} MinFreeGB=${MIN_FREE_GB}"
print_disk_space "/big7_disk"
if [[ "${DRY_RUN}" != true ]]; then
    require_min_free_gb "${RUN_DIR}" "${MIN_FREE_GB}" "resume_paired:${SAMPLE}:${CANONICAL_MODE}"
fi

validate_dir "${RUN_DIR}" "Canonical run dir"
validate_dir "${LONGPHASE_DIR}" "LongPhase output dir"
validate_file "${TAGGED_BAM}" "Tagged BAM"
validate_file "${TAGGED_BAI}" "Tagged BAM index"
validate_file "${SC_VCF}" "Somatic calling VCF"
validate_file "${TRUTH_VCF}" "Truth VCF"
if [[ -n "${TRUTH_BED}" ]]; then
    validate_file "${TRUTH_BED}" "Truth BED"
fi

if [[ -f "${HAPLOTAG_QC_LONGPHASE}" ]]; then
    run_cmd cp -f "${HAPLOTAG_QC_LONGPHASE}" "${HAPLOTAG_QC_ROOT}"
else
    log_warn "[Resume paired] haplotag_qc.tsv missing; continuing without per-sample haplotag QC."
fi

if [[ -f "${TP_VCF}" && -f "${FP_VCF}" && -f "${COUNTS_FILE}" ]]; then
    log_info "[Resume paired] Reusing existing TP/FP split outputs"
else
    log_info "[Resume paired] Splitting _sc.vcf into TP/FP"
    SPLIT_CMD=(
        bash
        "${PROJECT_ROOT}/scripts/pipeline/utils/benchmark_split_snv_vcf.sh"
        --input-vcf "${SC_VCF}"
        --truth-vcf "${TRUTH_VCF}"
        --output-dir "${LONGPHASE_DIR}"
    )
    if [[ -n "${TRUTH_BED}" ]]; then
        SPLIT_CMD+=(--truth-bed "${TRUTH_BED}")
    fi
    run_cmd "${SPLIT_CMD[@]}"
fi

if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${TP_VCF}" "TP VCF"
    validate_file "${FP_VCF}" "FP VCF"
    validate_file "${COUNTS_FILE}" "Variant counts"
    # shellcheck disable=SC1090
    source "${COUNTS_FILE}"
else
    TP_COUNT=0
    FP_COUNT=0
fi

cat > "${RUN_CONTEXT}" <<EOF
{
  "created_at": "$(date '+%Y-%m-%d %H:%M:%S %Z')",
  "sample": "${SAMPLE}",
  "sample_set": "complete_matrix_resume",
  "platform": "${PLATFORM}",
  "analysis_mode": "${PIPELINE_MODE}",
  "truth_set": "${TRUTH_SET}",
  "truth_vcf": "${TRUTH_VCF}",
  "truth_bed": "${TRUTH_BED}",
  "truth_total": ${TRUTH_TOTAL},
  "caller_name": "ClairS",
  "somatic_vcf": "${CALLER_VCF_PATH}",
  "source_run_dir": "${RUN_DIR}",
  "resume_from_tagged": true
}
EOF
run_cmd cp -f "${RUN_CONTEXT}" "${ROUND_CONTEXT}"

if [[ -f "${TP_SUMMARY}" && -f "${FP_SUMMARY}" ]]; then
    log_info "[Resume paired] Reusing existing InterSubMod summaries"
else
    log_info "[Resume paired] Running InterSubMod"
    run_cmd bash "${PROJECT_ROOT}/scripts/pipeline/steps/02_intersubmod.sh" \
        --tagged-bam "${TAGGED_BAM}" \
        --normal-bam "${NORMAL_BAM}" \
        --tp-vcf "${TP_VCF}" \
        --fp-vcf "${FP_VCF}" \
        --output-dir "${RUN_DIR}" \
        --threads "${LOCAL_THREADS}"
fi

if [[ -f "${METRICS_JSON}" ]]; then
    log_info "[Resume paired] Reusing existing metrics.json"
else
    log_info "[Resume paired] Running filter analysis"
    run_cmd python3 "${PROJECT_ROOT}/scripts/pipeline/steps/03_filter_analysis.py" \
        --tp-summary "${TP_SUMMARY}" \
        --fp-summary "${FP_SUMMARY}" \
        --tp-count "${TP_COUNT}" \
        --fp-count "${FP_COUNT}" \
        --truth-total "${TRUTH_TOTAL}" \
        --tp-vcf-file "${TP_VCF}" \
        --fp-vcf-file "${FP_VCF}" \
        --output-dir "${RUN_DIR}" \
        --sample "${SAMPLE}" \
        --mode "${PIPELINE_MODE}"
fi

log_info "[Resume paired] Building benchmark comparison"
run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/compare_benchmark_f1.py" \
    --run-dir "${RUN_DIR}" \
    --context-json "${RUN_CONTEXT}" \
    --output-tsv "${BENCHMARK_TSV}" \
    --output-md "${BENCHMARK_MD}"

if [[ "${SKIP_METHOD_VALIDATION}" != true ]]; then
    if [[ -f "${METHOD_VALIDATION_TSV}" && -f "${AGREEMENT_TSV}" ]]; then
        log_info "[Resume paired] Reusing existing method validation outputs"
    else
        log_info "[Resume paired] Running method design validation"
        run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/validate_method_design.py" \
            --summary-csv "${TP_SUMMARY}" \
            --region-root "${RUN_DIR}/intersubmod_tp" \
            --haplotag-qc "${HAPLOTAG_QC_ROOT}" \
            --sample "${SAMPLE}" \
            --output-dir "${RUN_DIR}"
    fi
fi

if [[ "${SKIP_DASHBOARD}" != true ]]; then
    log_info "[Resume paired] Building round dashboard"
    run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/build_round_dashboard.py" \
        --sample-dir "${RUN_DIR}"
fi

validate_file "${BENCHMARK_TSV}" "Benchmark comparison TSV"
validate_file "${ROUND_SUMMARY_MD}" "Round summary"
log_info "[Resume paired] Complete: ${RUN_DIR}"
