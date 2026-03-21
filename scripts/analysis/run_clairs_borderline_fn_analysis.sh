#!/bin/bash
# run_clairs_borderline_fn_analysis.sh
#
# Pool B FN analysis for ClairS borderline rescue experiment.
#
# Steps:
#   1. Extract non-PASS variants from ClairS VCF (exclude NonSomatic)
#   2. Scope to truth BED region
#   3. Intersect with truth VCF via bcftools isec → Pool B FN / Pool B non-FN
#   4. Run Python analysis + rescue rule evaluation
#   5. (Optional) Check methylation coverage for combined rule evaluation
#
# Usage:
#   ./scripts/analysis/run_clairs_borderline_fn_analysis.sh \
#     --snv-vcf /path/to/step01_clairs_to/snv.vcf.gz \
#     --truth-vcf /path/to/truth.vcf.gz \
#     --truth-bed /path/to/truth.bed \
#     --significance-csv /path/to/step05_intersubmod/significance_summary.csv \
#     --baseline-tp 28396 \
#     --baseline-fp 11843 \
#     --truth-total 39447 \
#     --sample HCC1395 \
#     --output-dir /path/to/output/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

# ── Defaults ──────────────────────────────────────────────────────────────────
SNV_VCF=""
TRUTH_VCF=""
TRUTH_BED=""
SIGNIFICANCE_CSV=""
BASELINE_TP=""
BASELINE_FP=""
TRUTH_TOTAL=""
SAMPLE="HCC1395"
OUTPUT_DIR=""
PYTHON="${PYTHON:-python3}"
BCFTOOLS_BIN="${BCFTOOLS:-$(command -v bcftools 2>/dev/null || true)}"
BGZIP_BIN="${BGZIP:-$(command -v bgzip 2>/dev/null || true)}"
TABIX_BIN="${TABIX:-$(command -v tabix 2>/dev/null || true)}"

# ── Argument Parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --snv-vcf)         SNV_VCF="$2";          shift 2 ;;
        --truth-vcf)       TRUTH_VCF="$2";         shift 2 ;;
        --truth-bed)       TRUTH_BED="$2";         shift 2 ;;
        --significance-csv) SIGNIFICANCE_CSV="$2"; shift 2 ;;
        --baseline-tp)     BASELINE_TP="$2";       shift 2 ;;
        --baseline-fp)     BASELINE_FP="$2";       shift 2 ;;
        --truth-total)     TRUTH_TOTAL="$2";       shift 2 ;;
        --sample)          SAMPLE="$2";            shift 2 ;;
        --output-dir)      OUTPUT_DIR="$2";        shift 2 ;;
        --python)          PYTHON="$2";            shift 2 ;;
        --bcftools)        BCFTOOLS_BIN="$2";      shift 2 ;;
        -h|--help)
            cat <<'HELP'
Usage: run_clairs_borderline_fn_analysis.sh [OPTIONS]

Required:
  --snv-vcf PATH          ClairS-TO SNV VCF (contains PASS + non-PASS)
  --truth-vcf PATH        SEQC2 truth set VCF (bgzipped + tabixed)
  --truth-bed PATH        High-confidence BED regions
  --baseline-tp INT       ClairS-TO PASS baseline TP count
  --baseline-fp INT       ClairS-TO PASS baseline FP count
  --truth-total INT       Total truth SNV count
  --output-dir PATH       Output directory

Optional:
  --significance-csv PATH InterSubMod significance_summary.csv for combined rules
  --sample NAME           Sample name (default: HCC1395)
  --python PATH           Python executable (default: python3)
  --bcftools PATH         bcftools binary override
HELP
            exit 0
            ;;
        *)
            log_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# ── Validation ─────────────────────────────────────────────────────────────────
for var in SNV_VCF TRUTH_VCF TRUTH_BED BASELINE_TP BASELINE_FP TRUTH_TOTAL OUTPUT_DIR; do
    if [[ -z "${!var}" ]]; then
        log_error "--${var//_/-} is required"
        exit 1
    fi
done
validate_file "${SNV_VCF}" "SNV VCF"
validate_file "${TRUTH_VCF}" "Truth VCF"
validate_file "${TRUTH_BED}" "Truth BED"
if [[ -n "${SIGNIFICANCE_CSV}" ]]; then
    validate_file "${SIGNIFICANCE_CSV}" "Significance CSV"
fi
if [[ -z "${BCFTOOLS_BIN}" ]]; then
    log_error "bcftools not found. Set BCFTOOLS env or use --bcftools."
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
WORK_DIR="${OUTPUT_DIR}/work"
mkdir -p "${WORK_DIR}"

log_info "=== ClairS Pool B FN Analysis ==="
log_info "  Sample:       ${SAMPLE}"
log_info "  SNV VCF:      ${SNV_VCF}"
log_info "  Truth VCF:    ${TRUTH_VCF}"
log_info "  Truth BED:    ${TRUTH_BED}"
log_info "  Baseline:     TP=${BASELINE_TP}, FP=${BASELINE_FP}, truth_total=${TRUTH_TOTAL}"
log_info "  Output dir:   ${OUTPUT_DIR}"

# ─── Step 1: Extract non-PASS, non-NonSomatic variants ────────────────────────
log_info "[Step 1] Extracting Pool B candidates (non-PASS excluding NonSomatic)..."
POOL_B_CANDIDATES="${WORK_DIR}/pool_b_candidates.vcf.gz"

# Include all non-PASS that do NOT have NonSomatic in their FILTER field
"${BCFTOOLS_BIN}" view \
    -e 'FILTER=="PASS" || FILTER~"NonSomatic"' \
    "${SNV_VCF}" \
    -O z -o "${POOL_B_CANDIDATES}"
"${TABIX_BIN}" -f "${POOL_B_CANDIDATES}"

POOL_B_COUNT=$("${BCFTOOLS_BIN}" view -H "${POOL_B_CANDIDATES}" 2>/dev/null | wc -l)
log_info "  Pool B candidates (non-PASS, non-NonSomatic): ${POOL_B_COUNT}"

# ─── Step 2: Scope to truth BED ───────────────────────────────────────────────
log_info "[Step 2] Scoping Pool B to truth high-confidence BED regions..."
POOL_B_SCOPED="${WORK_DIR}/pool_b_scoped.vcf.gz"

"${BCFTOOLS_BIN}" view \
    -R "${TRUTH_BED}" \
    "${POOL_B_CANDIDATES}" \
    -O z -o "${POOL_B_SCOPED}"
"${TABIX_BIN}" -f "${POOL_B_SCOPED}"

SCOPED_COUNT=$("${BCFTOOLS_BIN}" view -H "${POOL_B_SCOPED}" 2>/dev/null | wc -l)
log_info "  Pool B after BED scoping: ${SCOPED_COUNT}"

# Ensure truth VCF is tabixed
TRUTH_SCOPED="${WORK_DIR}/truth_scoped.vcf.gz"
"${BCFTOOLS_BIN}" view \
    -R "${TRUTH_BED}" \
    "${TRUTH_VCF}" \
    -O z -o "${TRUTH_SCOPED}"
"${TABIX_BIN}" -f "${TRUTH_SCOPED}"

# ─── Step 3: bcftools isec → Pool B FN vs Pool B non-FN ───────────────────────
log_info "[Step 3] Running bcftools isec to identify Pool B FN vs non-FN..."
ISEC_DIR="${WORK_DIR}/isec"
rm -rf "${ISEC_DIR}"

"${BCFTOOLS_BIN}" isec \
    -c none \
    -p "${ISEC_DIR}" \
    "${POOL_B_SCOPED}" \
    "${TRUTH_SCOPED}"

# bcftools isec output:
#   0000.vcf = in pool_b_scoped only       → Pool B non-FN (called but not truth: false rescues if restored)
#   0001.vcf = in truth_scoped only        → (irrelevant: truth not in pool_b)
#   0002.vcf = in both pool_b AND truth    → Pool B FN (non-PASS but IS truth TP)
#   0003.vcf = in both (from truth side)   → same intersection from truth side

POOL_B_FN_VCF="${OUTPUT_DIR}/pool_b_fn.vcf.gz"
POOL_B_NON_FN_VCF="${OUTPUT_DIR}/pool_b_non_fn.vcf.gz"

# Pool B FN = in pool_b_scoped AND in truth (file 0002 = left side of intersection)
"${BCFTOOLS_BIN}" view "${ISEC_DIR}/0002.vcf" -O z -o "${POOL_B_FN_VCF}"
"${TABIX_BIN}" -f "${POOL_B_FN_VCF}"

# Pool B non-FN = in pool_b_scoped only (file 0000 = left only)
"${BCFTOOLS_BIN}" view "${ISEC_DIR}/0000.vcf" -O z -o "${POOL_B_NON_FN_VCF}"
"${TABIX_BIN}" -f "${POOL_B_NON_FN_VCF}"

FN_COUNT=$("${BCFTOOLS_BIN}" view -H "${POOL_B_FN_VCF}" 2>/dev/null | wc -l)
NON_FN_COUNT=$("${BCFTOOLS_BIN}" view -H "${POOL_B_NON_FN_VCF}" 2>/dev/null | wc -l)
log_info "  Pool B FN (in truth): ${FN_COUNT}"
log_info "  Pool B non-FN (not in truth): ${NON_FN_COUNT}"

# ─── Go/No-Go Check ───────────────────────────────────────────────────────────
if [[ "${FN_COUNT}" -lt 30 ]]; then
    log_info ""
    log_info "=== Go/No-Go: NO-GO ==="
    log_info "  Pool B FN (${FN_COUNT}) < 30 threshold."
    log_info "  Recommendation: Pool B is not an effective FN source."
    log_info "  Consider exploring Pool A-light (phasing intermediates)."
    log_info ""
fi

# ─── Step 4: Python Analysis ───────────────────────────────────────────────────
log_info "[Step 4] Running Python analysis..."

PYTHON_ARGS=(
    "${SCRIPT_DIR}/analyze_clairs_borderline_fn.py"
    --pool-b-fn-vcf "${POOL_B_FN_VCF}"
    --pool-b-non-fn-vcf "${POOL_B_NON_FN_VCF}"
    --baseline-tp "${BASELINE_TP}"
    --baseline-fp "${BASELINE_FP}"
    --truth-total "${TRUTH_TOTAL}"
    --sample "${SAMPLE}"
    --output-dir "${OUTPUT_DIR}"
)

if [[ -n "${SIGNIFICANCE_CSV}" ]]; then
    PYTHON_ARGS+=(--significance-csv "${SIGNIFICANCE_CSV}")
fi

cd "${SCRIPT_DIR}" && "${PYTHON}" "${PYTHON_ARGS[@]}"

# ─── Step 5: FILTER distribution quick stats ─────────────────────────────────
log_info "[Step 5] FILTER distribution of Pool B FN..."
"${BCFTOOLS_BIN}" view -H "${POOL_B_FN_VCF}" 2>/dev/null | \
    awk '{print $7}' | sort | uniq -c | sort -rn | head -20 > "${OUTPUT_DIR}/pool_b_fn_filter_dist.txt"
log_info "  Written to: ${OUTPUT_DIR}/pool_b_fn_filter_dist.txt"

# ─── Summary ──────────────────────────────────────────────────────────────────
log_info ""
log_info "=== Pool B Analysis Complete ==="
log_info "  Pool B FN:     ${FN_COUNT}"
log_info "  Pool B non-FN: ${NON_FN_COUNT}"
log_info "  Output:        ${OUTPUT_DIR}/"
log_info "  Key files:"
log_info "    pool_b_fn.vcf.gz             - Pool B FN variants (truth TP filtered by ClairS)"
log_info "    pool_b_non_fn.vcf.gz         - Pool B non-FN variants (not in truth)"
log_info "    pool_b_fn_stats.tsv          - Overall stats + Go/No-Go"
log_info "    pool_b_fn_distribution.tsv   - GQ/AF/FILTER distribution"
log_info "    pool_b_caller_rescue_rules.tsv  - Caller-only rescue rules"
log_info "    pool_b_combined_rescue_rules.tsv - Caller+Methyl combined rules"
log_info "    pool_b_joined_features.tsv   - Full feature table"
log_info "    pool_b_summary.md            - Human-readable summary"
