#!/bin/bash
# Resume an existing tumor-only pilot run in place using tagged BAM + step04 benchmark outputs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

RUN_DIR=""
LOCAL_THREADS="${THREADS}"
MIN_FREE_GB="${MIN_FREE_GB_DEFAULT}"
SKIP_METHOD_VALIDATION=false
SKIP_DASHBOARD=false
DRY_RUN=false

show_help() {
    cat <<'EOF'
Usage: resume_to_pilot_from_tagged_outputs.sh --run-dir DIR [OPTIONS]

Options:
  --run-dir DIR               Existing tumor-only pilot round directory
  --threads N                 Threads for InterSubMod
  --min-free-gb N             Minimum free disk on /big7_disk before resuming
  --skip-method-validation    Skip validate_method_design.py
  --skip-dashboard            Skip build_round_dashboard.py
  --dry-run                   Print commands without executing
  -h, --help                  Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
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

if [[ -z "${RUN_DIR}" ]]; then
    log_error "--run-dir is required."
    exit 1
fi

run_cmd() {
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] $*" >&2
    else
        "$@"
    fi
}

RUN_DIR="$(realpath "${RUN_DIR}")"
ROUND_CONTEXT="${RUN_DIR}/round_context.json"
STEP02_DIR="${RUN_DIR}/step02_benchmark_clairs_to"
STEP03_DIR="${RUN_DIR}/step03_longphase_to"
STEP04_DIR="${RUN_DIR}/step04_benchmark_longphase_to"
STEP05_DIR="${RUN_DIR}/step05_intersubmod"

TAGGED_BAM="${STEP03_DIR}/tumor_tagged.bam"
TAGGED_BAI="${TAGGED_BAM}.bai"
HAPLOTAG_QC="${RUN_DIR}/haplotag_qc.tsv"
TP_VCF="${STEP04_DIR}/filtered_snv_tp.vcf.gz"
FP_VCF="${STEP04_DIR}/filtered_snv_fp.vcf.gz"
COUNTS_FILE="${STEP04_DIR}/variant_counts.txt"
TP_SUMMARY="${STEP05_DIR}/intersubmod_tp/significance_summary.csv"
FP_SUMMARY="${STEP05_DIR}/intersubmod_fp/significance_summary.csv"
METRICS_JSON="${RUN_DIR}/metrics.json"
BENCHMARK_TSV="${RUN_DIR}/benchmark_comparison.tsv"
BENCHMARK_MD="${RUN_DIR}/benchmark_comparison.md"
METHOD_VALIDATION_TSV="${RUN_DIR}/method_design_validation.tsv"
AGREEMENT_TSV="${RUN_DIR}/label_cluster_agreement.tsv"
ROUND_SUMMARY_MD="${RUN_DIR}/round_summary.md"
TO_STAGE_METRICS="${RUN_DIR}/to_stage_metrics.json"
RESUME_LOG="${RUN_DIR}/resume_to_from_tagged.log"

exec > >(tee -a "${RESUME_LOG}") 2> >(tee -a "${RESUME_LOG}" >&2)

validate_dir "${RUN_DIR}" "TO round dir"
validate_file "${ROUND_CONTEXT}" "round_context.json"
validate_file "${TAGGED_BAM}" "Tagged BAM"
validate_file "${TAGGED_BAI}" "Tagged BAM index"
validate_file "${TP_VCF}" "LongPhase-TO TP VCF"
validate_file "${FP_VCF}" "LongPhase-TO FP VCF"
validate_file "${COUNTS_FILE}" "LongPhase-TO variant counts"

if [[ "${DRY_RUN}" != true ]]; then
    require_min_free_gb "${RUN_DIR}" "${MIN_FREE_GB}" "resume_to_pilot:${RUN_DIR}"
fi

read_context_field() {
    local key="$1"
    python3 - "$ROUND_CONTEXT" "$key" <<'PY'
import json
import sys
path, key = sys.argv[1], sys.argv[2]
with open(path, "r", encoding="utf-8") as fh:
    payload = json.load(fh)
value = payload.get(key, "")
if value is None:
    value = ""
print(value)
PY
}

SAMPLE="$(read_context_field sample)"
TRUTH_TOTAL="$(read_context_field truth_total)"
ANALYSIS_MODE="$(read_context_field analysis_mode)"
if [[ -z "${ANALYSIS_MODE}" ]]; then
    ANALYSIS_MODE="to-pure"
fi

source "${COUNTS_FILE}"
LP_PASS_TOTAL="${PASS_TOTAL}"
LP_TP_COUNT="${TP_COUNT}"
LP_FP_COUNT="${FP_COUNT}"
LP_FN_COUNT="${FN_COUNT:-0}"

BASELINE_COUNTS="${STEP02_DIR}/variant_counts.txt"
if [[ -f "${BASELINE_COUNTS}" ]]; then
    source "${BASELINE_COUNTS}"
    BASELINE_PASS_TOTAL="${PASS_TOTAL}"
    BASELINE_TP_COUNT="${TP_COUNT}"
    BASELINE_FP_COUNT="${FP_COUNT}"
    BASELINE_FN_COUNT="${FN_COUNT:-0}"
else
    BASELINE_PASS_TOTAL=0
    BASELINE_TP_COUNT=0
    BASELINE_FP_COUNT=0
    BASELINE_FN_COUNT=0
fi

if [[ ! -f "${HAPLOTAG_QC}" ]]; then
    log_warn "[Resume TO] haplotag_qc.tsv missing; continuing without per-sample haplotag QC."
fi

if [[ -f "${TP_SUMMARY}" && -f "${FP_SUMMARY}" ]]; then
    log_info "[Resume TO] Reusing existing InterSubMod summaries"
else
    log_info "[Resume TO] Running InterSubMod"
    run_cmd bash "${PROJECT_ROOT}/scripts/pipeline/steps/02_intersubmod.sh" \
        --tagged-bam "${TAGGED_BAM}" \
        --tp-vcf "${TP_VCF}" \
        --fp-vcf "${FP_VCF}" \
        --output-dir "${STEP05_DIR}" \
        --threads "${LOCAL_THREADS}"
fi

if [[ -f "${METRICS_JSON}" ]]; then
    log_info "[Resume TO] Reusing existing metrics.json"
else
    log_info "[Resume TO] Running filter analysis"
    run_cmd python3 "${PROJECT_ROOT}/scripts/pipeline/steps/03_filter_analysis.py" \
        --tp-summary "${TP_SUMMARY}" \
        --fp-summary "${FP_SUMMARY}" \
        --tp-count "${LP_TP_COUNT}" \
        --fp-count "${LP_FP_COUNT}" \
        --truth-total "${TRUTH_TOTAL}" \
        --tp-vcf-file "${TP_VCF}" \
        --fp-vcf-file "${FP_VCF}" \
        --sample "${SAMPLE}" \
        --mode "${ANALYSIS_MODE}" \
        --purity "100" \
        --output-dir "${RUN_DIR}"
fi

if [[ "${SKIP_METHOD_VALIDATION}" != true ]]; then
    if [[ -f "${METHOD_VALIDATION_TSV}" && -f "${AGREEMENT_TSV}" ]]; then
        log_info "[Resume TO] Reusing existing method validation outputs"
    else
        log_info "[Resume TO] Running method design validation"
        VALIDATE_CMD=(
            python3
            "${PROJECT_ROOT}/scripts/analysis/validate_method_design.py"
            --summary-csv "${TP_SUMMARY}"
            --summary-csv "${FP_SUMMARY}"
            --region-root "${STEP05_DIR}/intersubmod_tp"
            --region-root "${STEP05_DIR}/intersubmod_fp"
            --sample "${SAMPLE}"
            --output-dir "${RUN_DIR}"
        )
        if [[ -f "${HAPLOTAG_QC}" ]]; then
            VALIDATE_CMD+=(--haplotag-qc "${HAPLOTAG_QC}")
        fi
        run_cmd "${VALIDATE_CMD[@]}"
    fi
fi

log_info "[Resume TO] Building benchmark comparison"
run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/compare_benchmark_f1.py" \
    --run-dir "${RUN_DIR}" \
    --context-json "${ROUND_CONTEXT}" \
    --baseline-counts-file "${STEP04_DIR}/variant_counts.txt" \
    --baseline-method "LongPhase-TO" \
    --output-tsv "${BENCHMARK_TSV}" \
    --output-md "${BENCHMARK_MD}"

log_info "[Resume TO] Writing to_stage_metrics.json"
python3 - "${RUN_DIR}" "${BASELINE_PASS_TOTAL}" "${BASELINE_TP_COUNT}" "${BASELINE_FP_COUNT}" "${BASELINE_FN_COUNT}" "${LP_PASS_TOTAL}" "${LP_TP_COUNT}" "${LP_FP_COUNT}" "${LP_FN_COUNT}" <<'PY'
import json
import sys
from pathlib import Path

round_root = Path(sys.argv[1])
baseline_calls, baseline_tp, baseline_fp, baseline_fn = map(int, sys.argv[2:6])
lp_calls, lp_tp, lp_fp, lp_fn = map(int, sys.argv[6:10])
metrics_path = round_root / "metrics.json"
summary_path = round_root / "to_stage_metrics.json"

def calc(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "calls_total": tp + fp if tp + fp else 0,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
    }

payload = {
    "baseline_clairs_to": calc(baseline_tp, baseline_fp, baseline_fn),
    "longphase_to": calc(lp_tp, lp_fp, lp_fn),
}
if metrics_path.exists():
    with metrics_path.open("r", encoding="utf-8") as fh:
        payload["intersubmod"] = json.load(fh).get("filtered", {})
with summary_path.open("w", encoding="utf-8") as fh:
    json.dump(payload, fh, indent=2, ensure_ascii=False)
    fh.write("\n")
PY

if [[ "${SKIP_DASHBOARD}" != true ]]; then
    log_info "[Resume TO] Building round dashboard"
    run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/build_round_dashboard.py" \
        --sample-dir "${RUN_DIR}"
fi

validate_file "${METRICS_JSON}" "metrics.json"
validate_file "${BENCHMARK_TSV}" "benchmark_comparison.tsv"
validate_file "${ROUND_SUMMARY_MD}" "round_summary.md"
validate_file "${TO_STAGE_METRICS}" "to_stage_metrics.json"
log_info "[Resume TO] Complete: ${RUN_DIR}"
