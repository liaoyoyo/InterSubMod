#!/bin/bash
# 01_longphase_s.sh - Run LongPhase-S somatic_haplotag
#
# Produces: tagged BAM, somatic calling VCF (_sc.vcf), benchmark metrics
# Then splits somatic VCF into TP/FP by intersecting with truth set.
#
# Usage:
#   ./scripts/pipeline/steps/01_longphase_s.sh \
#     --sample HCC1395 --output-dir DIR [--threads N] [--dry-run]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

SAMPLE=""
SAMPLE_OUTPUT_DIR=""
LOCAL_THREADS="${THREADS}"
LOCAL_INDEX_THREADS=""
OVERRIDE_SOMATIC_VCF=""
OVERRIDE_TUMOR_BAM=""
DRY_RUN=false

bcftools_sort_supports_threads() {
    "${BCFTOOLS}" sort --help 2>&1 | grep -q -- "--threads"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)     SAMPLE="$2"; shift 2 ;;
        --output-dir) SAMPLE_OUTPUT_DIR="$2"; shift 2 ;;
        --threads)    LOCAL_THREADS="$2"; shift 2 ;;
        --index-threads) LOCAL_INDEX_THREADS="$2"; shift 2 ;;
        --somatic-vcf) OVERRIDE_SOMATIC_VCF="$2"; shift 2 ;;
        --tumor-bam)   OVERRIDE_TUMOR_BAM="$2"; shift 2 ;;
        --dry-run)    DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: $0 --sample NAME --output-dir DIR [--threads N] [--somatic-vcf FILE] [--tumor-bam FILE] [--dry-run]" >&2
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${SAMPLE}" ]] || [[ -z "${SAMPLE_OUTPUT_DIR}" ]]; then
    log_error "Both --sample and --output-dir are required."
    exit 1
fi

if [[ -z "${LOCAL_INDEX_THREADS}" ]]; then
    LOCAL_INDEX_THREADS="$(recommended_index_threads "${LOCAL_THREADS}")"
fi

# ============================================================================
# Load Sample Configuration
# ============================================================================

eval "$(get_sample_config "${SAMPLE}")"

if [[ -n "${OVERRIDE_SOMATIC_VCF:-}" ]]; then
    log_info "Overriding Somatic VCF: ${OVERRIDE_SOMATIC_VCF}"
    SOMATIC_VCF="${OVERRIDE_SOMATIC_VCF}"
fi

if [[ -n "${OVERRIDE_TUMOR_BAM:-}" ]]; then
    log_info "Overriding Tumor BAM: ${OVERRIDE_TUMOR_BAM}"
    TUMOR_BAM="${OVERRIDE_TUMOR_BAM}"
fi

LONGPHASE_OUTPUT_DIR="${SAMPLE_OUTPUT_DIR}/longphase_s"
mkdir -p "${LONGPHASE_OUTPUT_DIR}"

# ============================================================================
# Input Validation
# ============================================================================

log_info "[Step 01] LongPhase-S somatic_haplotag for ${SAMPLE}"
log_info "  Tumor BAM:    ${TUMOR_BAM}"
log_info "  Normal BAM:   ${NORMAL_BAM}"
log_info "  Somatic VCF:  ${SOMATIC_VCF}"
log_info "  Truth VCF:    ${TRUTH_VCF}"
log_info "  Truth BED:    ${TRUTH_BED:-"(none - whole genome)"}"
log_info "  Threads:      ${LOCAL_THREADS}"
log_info "  Index thr:    ${LOCAL_INDEX_THREADS}"
log_info "  Output dir:   ${LONGPHASE_OUTPUT_DIR}"

validate_file "${TUMOR_BAM}" "Tumor BAM"
validate_file "${NORMAL_BAM}" "Normal BAM"
validate_file "${SOMATIC_VCF}" "Somatic VCF"
validate_file "${TRUTH_VCF}" "Truth VCF"
if [[ -n "${TRUTH_BED}" ]]; then
    validate_file "${TRUTH_BED}" "Truth BED"
fi
validate_file "${REFERENCE}" "Reference FASTA"
validate_file "${LONGPHASE_S_BIN}" "LongPhase-S binary"

# ============================================================================
# Step 1: Prepare Germline Phased VCF
# ============================================================================

log_info "  Resolving germline phased VCF..."
GERMLINE_ARGS=(--germline-dir "${GERMLINE_PHASED_DIR}" --output-dir "${LONGPHASE_OUTPUT_DIR}")
if [[ "${DRY_RUN}" == true ]]; then
    GERMLINE_ARGS+=(--dry-run)
fi

GERMLINE_PHASED_VCF=$("${SCRIPT_DIR}/00_prepare_germline.sh" "${GERMLINE_ARGS[@]}")
log_info "  Germline VCF: ${GERMLINE_PHASED_VCF}"

if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${GERMLINE_PHASED_VCF}" "Germline phased VCF"
fi

# ============================================================================
# Step 2: Run LongPhase-S
# ============================================================================

OUTPUT_PREFIX="${LONGPHASE_OUTPUT_DIR}/${SAMPLE}_tagged"

# Pre-filter somatic VCF to PASS-only to avoid LongPhase-S crash
# on LowQual variants (e.g., chrX:72880028 with QUAL=0.000)
PASS_SOMATIC_VCF="${LONGPHASE_OUTPUT_DIR}/somatic_pass.vcf.gz"

if [[ "${DRY_RUN}" == true ]]; then
    log_info "[DRY-RUN] bcftools view -f PASS ${SOMATIC_VCF} -Oz -o ${PASS_SOMATIC_VCF}"
else
    log_info "  Pre-filtering somatic VCF to PASS variants only..."
    
    # Handle VCF input (gzip or plain) and fix GQ header type (Int -> Float) for pileup VCFs
    # Some ClairS pileup VCFs contain decimal GQ values, which bcftools rejects when header is Integer.
    if [[ "${SOMATIC_VCF}" == *.gz ]]; then
        CAT_CMD="zcat"
    else
        CAT_CMD="cat"
    fi

    HEADER_SNIPPET="$("${CAT_CMD}" "${SOMATIC_VCF}" | sed -n '1,200p')"
    if grep -q 'ID=GQ,Number=1,Type=Integer' <<< "${HEADER_SNIPPET}"; then
        log_info "  Detected GQ Integer header; patching to Float for bcftools compatibility."
    fi

    "${CAT_CMD}" "${SOMATIC_VCF}" | \
    sed 's/ID=GQ,Number=1,Type=Integer/ID=GQ,Number=1,Type=Float/' | \
    "${BCFTOOLS}" view -f PASS -Oz -o "${PASS_SOMATIC_VCF}"
    
    "${BCFTOOLS}" index "${PASS_SOMATIC_VCF}"
    PASS_COUNT=$("${BCFTOOLS}" view -H "${PASS_SOMATIC_VCF}" | wc -l)
    log_info "  PASS somatic variants: ${PASS_COUNT}"
fi

LONGPHASE_CMD=(
    "${LONGPHASE_S_BIN}" somatic_haplotag
    -s "${GERMLINE_PHASED_VCF}"
    -b "${NORMAL_BAM}"
    --tumor-snv-file "${PASS_SOMATIC_VCF}"
    --tumor-bam-file "${TUMOR_BAM}"
    -r "${REFERENCE}"
    -t "${LOCAL_THREADS}"
    --tagSupplementary
    -q 20
    --output-somatic-vcf
    --truth-vcf "${TRUTH_VCF}"
    -o "${OUTPUT_PREFIX}"
)
# Add --truth-bed only if BED file is specified
if [[ -n "${TRUTH_BED}" ]]; then
    LONGPHASE_CMD+=(--truth-bed "${TRUTH_BED}")
fi

if [[ "${DRY_RUN}" == true ]]; then
    log_info "[DRY-RUN] LongPhase-S command:"
    echo "  ${LONGPHASE_CMD[*]}" >&2
else
    log_info "  Running LongPhase-S..."
    print_disk_space

    LONGPHASE_LOG="${LONGPHASE_OUTPUT_DIR}/longphase_s.log"

    /usr/bin/time -v "${LONGPHASE_CMD[@]}" 2>&1 | tee "${LONGPHASE_LOG}"

    log_info "  LongPhase-S completed. Log: ${LONGPHASE_LOG}"
fi

# ============================================================================
# Step 3: Index Tagged BAM
# ============================================================================

TAGGED_BAM="${OUTPUT_PREFIX}.bam"

if [[ "${DRY_RUN}" == true ]]; then
    log_info "[DRY-RUN] samtools index ${TAGGED_BAM}"
else
    if [[ -f "${TAGGED_BAM}" ]]; then
        log_info "  Indexing tagged BAM..."
        "${SAMTOOLS}" index -@ "${LOCAL_INDEX_THREADS}" "${TAGGED_BAM}"
        log_info "  Tagged BAM indexed: ${TAGGED_BAM}.bai"

        HAPLOTAG_QC_TSV="${LONGPHASE_OUTPUT_DIR}/haplotag_qc.tsv"
        log_info "  Summarizing HP tags..."
        python3 "${SCRIPT_DIR}/../../analysis/haplotag_qc.py" \
            --bam "${TAGGED_BAM}" \
            --sample "${SAMPLE}" \
            --output-tsv "${HAPLOTAG_QC_TSV}" \
            --samtools "${SAMTOOLS}"
        log_info "  Haplotag QC: ${HAPLOTAG_QC_TSV}"
    else
        log_error "Tagged BAM not found: ${TAGGED_BAM}"
        exit 1
    fi
fi

# ============================================================================
# Step 4: Generate TP/FP VCFs by intersecting with truth set
# ============================================================================

SOMATIC_SC_VCF="${OUTPUT_PREFIX}_sc.vcf"

if [[ "${DRY_RUN}" == true ]]; then
    log_info "[DRY-RUN] Generating TP/FP VCFs from somatic calling output"
    log_info "[DRY-RUN] bcftools view -f PASS ${SOMATIC_SC_VCF} | bcftools sort | bcftools view -Oz > somatic_pass.vcf.gz"
    log_info "[DRY-RUN] bcftools isec -p isec_output somatic_pass.vcf.gz truth.vcf.gz -R truth.bed"
else
    if [[ ! -f "${SOMATIC_SC_VCF}" ]]; then
        log_error "Somatic calling VCF not found: ${SOMATIC_SC_VCF}"
        exit 1
    fi

    log_info "  Splitting somatic VCF into TP/FP..."

    ISEC_DIR="${LONGPHASE_OUTPUT_DIR}/isec_tmp"
    PASS_VCF="${LONGPHASE_OUTPUT_DIR}/somatic_pass.vcf.gz"
    SORT_ARGS=()
    if bcftools_sort_supports_threads; then
        SORT_ARGS+=(--threads "${LOCAL_INDEX_THREADS}")
    else
        log_warn "  bcftools sort does not support --threads; falling back to single-process sort."
    fi

    # Extract PASS variants from somatic calling VCF
    "${BCFTOOLS}" view -f PASS "${SOMATIC_SC_VCF}" | \
        "${BCFTOOLS}" sort "${SORT_ARGS[@]}" -Oz -o "${PASS_VCF}"
    "${BCFTOOLS}" index "${PASS_VCF}"

    # Intersect with truth VCF (within HC regions if BED provided)
    mkdir -p "${ISEC_DIR}"
    ISEC_CMD=("${BCFTOOLS}" isec)
    if [[ -n "${TRUTH_BED}" ]]; then
        ISEC_CMD+=(-R "${TRUTH_BED}")
    fi
    ISEC_CMD+=(-p "${ISEC_DIR}" "${PASS_VCF}" "${TRUTH_VCF}")
    "${ISEC_CMD[@]}"

    # isec output:
    #   0000.vcf = in somatic PASS only (FP)
    #   0001.vcf = in truth only (FN, not used here)
    #   0002.vcf = in both from somatic (TP)
    #   0003.vcf = in both from truth (not used)

    TP_VCF="${LONGPHASE_OUTPUT_DIR}/filtered_snv_tp.vcf.gz"
    FP_VCF="${LONGPHASE_OUTPUT_DIR}/filtered_snv_fp.vcf.gz"

    "${BCFTOOLS}" view -Oz -o "${TP_VCF}" "${ISEC_DIR}/0002.vcf"
    "${BCFTOOLS}" index "${TP_VCF}"

    "${BCFTOOLS}" view -Oz -o "${FP_VCF}" "${ISEC_DIR}/0000.vcf"
    "${BCFTOOLS}" index "${FP_VCF}"

    # Count and report
    TP_COUNT=$(zgrep -c -v "^#" "${TP_VCF}" || echo 0)
    FP_COUNT=$(zgrep -c -v "^#" "${FP_VCF}" || echo 0)
    PASS_TOTAL=$(zgrep -c -v "^#" "${PASS_VCF}" || echo 0)

    log_info "  PASS variants: ${PASS_TOTAL}"
    log_info "  TP: ${TP_COUNT}, FP: ${FP_COUNT}"
    log_info "  TP VCF: ${TP_VCF}"
    log_info "  FP VCF: ${FP_VCF}"

    # Save counts for later use
    cat > "${LONGPHASE_OUTPUT_DIR}/variant_counts.txt" <<EOF
PASS_TOTAL=${PASS_TOTAL}
TP_COUNT=${TP_COUNT}
FP_COUNT=${FP_COUNT}
EOF

    # Cleanup isec temp
    rm -rf "${ISEC_DIR}"
    rm -f "${PASS_VCF}" "${PASS_VCF}.csi" "${PASS_VCF}.tbi"
fi

log_info "[Step 01] LongPhase-S complete."
