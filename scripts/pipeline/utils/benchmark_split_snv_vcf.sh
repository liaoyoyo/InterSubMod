#!/bin/bash
# benchmark_split_snv_vcf.sh - Scope PASS SNV calls to truth regions and split TP/FP/FN with bcftools isec
#
# Produces:
#   - filtered_snv_tp.vcf.gz / filtered_snv_fp.vcf.gz / filtered_snv_fn.vcf.gz
#   - optional plain tp.vcf / fp.vcf / fn.vcf
#   - variant_counts.txt
#
# Usage:
#   ./scripts/pipeline/utils/benchmark_split_snv_vcf.sh \
#     --input-vcf calls.vcf.gz \
#     --truth-vcf truth.vcf.gz \
#     --truth-bed truth.bed \
#     --output-dir out_dir

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

INPUT_VCF=""
TRUTH_VCF=""
TRUTH_BED=""
OUTPUT_DIR=""
FILTER_TAG="PASS"
WRITE_PLAIN=true
KEEP_TMP=true
LOCAL_BCFTOOLS="${BCFTOOLS}"
LOCAL_BGZIP="${BGZIP:-$(command -v bgzip || true)}"
LOCAL_TABIX="${TABIX:-$(command -v tabix || true)}"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-vcf)   INPUT_VCF="$2"; shift 2 ;;
        --truth-vcf)   TRUTH_VCF="$2"; shift 2 ;;
        --truth-bed)   TRUTH_BED="$2"; shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
        --filter-tag)  FILTER_TAG="$2"; shift 2 ;;
        --bcftools)    LOCAL_BCFTOOLS="$2"; shift 2 ;;
        --no-plain)    WRITE_PLAIN=false; shift ;;
        --dry-run)     DRY_RUN=true; shift ;;
        -h|--help)
            cat <<'HELP' >&2
Usage: benchmark_split_snv_vcf.sh --input-vcf FILE --truth-vcf FILE --truth-bed BED --output-dir DIR [OPTIONS]

Options:
  --filter-tag TAG   FILTER tag to keep from input VCF (default: PASS)
  --bcftools PATH    Override bcftools binary
  --no-plain         Do not emit uncompressed tp.vcf / fp.vcf / fn.vcf
  --dry-run          Print commands without executing
HELP
            exit 0
            ;;
        *)
            log_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "${INPUT_VCF}" || -z "${TRUTH_VCF}" || -z "${OUTPUT_DIR}" ]]; then
    log_error "--input-vcf, --truth-vcf and --output-dir are required."
    exit 1
fi

if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${INPUT_VCF}" "Input VCF"
    validate_file "${TRUTH_VCF}" "Truth VCF"
    if [[ -n "${TRUTH_BED}" ]]; then
        validate_file "${TRUTH_BED}" "Truth BED"
    fi
fi

mkdir -p "${OUTPUT_DIR}"

INPUT_SCOPED="${OUTPUT_DIR}/input_scoped_pass.vcf.gz"
TRUTH_SCOPED="${OUTPUT_DIR}/truth_scoped.vcf.gz"
ISEC_DIR="${OUTPUT_DIR}/isec_tmp"
TP_VCF_GZ="${OUTPUT_DIR}/filtered_snv_tp.vcf.gz"
FP_VCF_GZ="${OUTPUT_DIR}/filtered_snv_fp.vcf.gz"
FN_VCF_GZ="${OUTPUT_DIR}/filtered_snv_fn.vcf.gz"

run_cmd() {
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] $*" >&2
    else
        "$@"
    fi
}

prepare_bgzip_vcf() {
    local source_vcf="$1"
    local output_vcf_gz="$2"

    if [[ -z "${LOCAL_BGZIP}" || -z "${LOCAL_TABIX}" ]]; then
        log_error "bgzip/tabix not found, cannot prepare indexed VCF for: ${source_vcf}"
        exit 1
    fi

    log_info "[benchmark_split] Compressing plain VCF with bgzip: ${source_vcf}"
    run_cmd "${LOCAL_BGZIP}" -f -c "${source_vcf}" > "${output_vcf_gz}"
    run_cmd "${LOCAL_TABIX}" -f -p vcf "${output_vcf_gz}"
}

prepare_input_for_region_query() {
    local source_vcf="$1"
    local prepared_vcf="$2"

    if [[ "${source_vcf}" == *.gz ]]; then
        echo "${source_vcf}"
        return 0
    fi

    prepare_bgzip_vcf "${source_vcf}" "${prepared_vcf}"
    echo "${prepared_vcf}"
}

QUERY_INPUT_VCF="$(prepare_input_for_region_query "${INPUT_VCF}" "${OUTPUT_DIR}/input_region_query.vcf.gz")"

INPUT_VIEW_CMD=("${LOCAL_BCFTOOLS}" "view")
if [[ -n "${FILTER_TAG}" ]]; then
    INPUT_VIEW_CMD+=("-f" "${FILTER_TAG}")
fi
if [[ -n "${TRUTH_BED}" ]]; then
    INPUT_VIEW_CMD+=("-R" "${TRUTH_BED}")
fi
INPUT_VIEW_CMD+=("${QUERY_INPUT_VCF}" "-Oz" "-o" "${INPUT_SCOPED}")

TRUTH_VIEW_CMD=("${LOCAL_BCFTOOLS}" "view")
if [[ -n "${TRUTH_BED}" ]]; then
    TRUTH_VIEW_CMD+=("-R" "${TRUTH_BED}")
fi
TRUTH_VIEW_CMD+=("${TRUTH_VCF}" "-Oz" "-o" "${TRUTH_SCOPED}")

ISEC_CMD=("${LOCAL_BCFTOOLS}" "isec" "-p" "${ISEC_DIR}" "${INPUT_SCOPED}" "${TRUTH_SCOPED}")

log_info "[benchmark_split] Input VCF: ${INPUT_VCF}"
log_info "[benchmark_split] Truth VCF: ${TRUTH_VCF}"
log_info "[benchmark_split] Truth BED: ${TRUTH_BED:-<none>}"
log_info "[benchmark_split] Output dir: ${OUTPUT_DIR}"

run_cmd "${INPUT_VIEW_CMD[@]}"
run_cmd "${LOCAL_BCFTOOLS}" index -f "${INPUT_SCOPED}"
run_cmd "${TRUTH_VIEW_CMD[@]}"
run_cmd "${LOCAL_BCFTOOLS}" index -f "${TRUTH_SCOPED}"

run_cmd mkdir -p "${ISEC_DIR}"
run_cmd "${ISEC_CMD[@]}"

# bcftools isec output:
#   0000 = input only (FP)
#   0001 = truth only (FN)
#   0002 = common from input (TP)
#   0003 = common from truth
run_cmd "${LOCAL_BCFTOOLS}" view -Oz -o "${TP_VCF_GZ}" "${ISEC_DIR}/0002.vcf"
run_cmd "${LOCAL_BCFTOOLS}" index -f "${TP_VCF_GZ}"
run_cmd "${LOCAL_BCFTOOLS}" view -Oz -o "${FP_VCF_GZ}" "${ISEC_DIR}/0000.vcf"
run_cmd "${LOCAL_BCFTOOLS}" index -f "${FP_VCF_GZ}"
run_cmd "${LOCAL_BCFTOOLS}" view -Oz -o "${FN_VCF_GZ}" "${ISEC_DIR}/0001.vcf"
run_cmd "${LOCAL_BCFTOOLS}" index -f "${FN_VCF_GZ}"

if [[ "${WRITE_PLAIN}" == true ]]; then
    run_cmd "${LOCAL_BCFTOOLS}" view -Ov -o "${OUTPUT_DIR}/tp.vcf" "${TP_VCF_GZ}"
    run_cmd "${LOCAL_BCFTOOLS}" view -Ov -o "${OUTPUT_DIR}/fp.vcf" "${FP_VCF_GZ}"
    run_cmd "${LOCAL_BCFTOOLS}" view -Ov -o "${OUTPUT_DIR}/fn.vcf" "${FN_VCF_GZ}"
fi

if [[ "${DRY_RUN}" != true ]]; then
    PASS_TOTAL=$("${LOCAL_BCFTOOLS}" view -H "${INPUT_SCOPED}" | wc -l)
    TP_COUNT=$("${LOCAL_BCFTOOLS}" view -H "${TP_VCF_GZ}" | wc -l)
    FP_COUNT=$("${LOCAL_BCFTOOLS}" view -H "${FP_VCF_GZ}" | wc -l)
    FN_COUNT=$("${LOCAL_BCFTOOLS}" view -H "${FN_VCF_GZ}" | wc -l)
    TRUTH_TOTAL=$("${LOCAL_BCFTOOLS}" view -H "${TRUTH_SCOPED}" | wc -l)

    cat > "${OUTPUT_DIR}/variant_counts.txt" <<EOF
PASS_TOTAL=${PASS_TOTAL}
TP_COUNT=${TP_COUNT}
FP_COUNT=${FP_COUNT}
FN_COUNT=${FN_COUNT}
TRUTH_TOTAL=${TRUTH_TOTAL}
EOF

    log_info "[benchmark_split] PASS_TOTAL=${PASS_TOTAL}"
    log_info "[benchmark_split] TP=${TP_COUNT} FP=${FP_COUNT} FN=${FN_COUNT}"

fi

log_info "[benchmark_split] Complete."
