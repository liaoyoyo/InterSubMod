#!/bin/bash
# 00_prepare_germline.sh - Prepare germline phased VCF for LongPhase-S
#
# Searches for merged germline VCF in Clair3 output directory.
# If only per-chromosome VCFs exist, merges them with bcftools concat.
#
# Usage:
#   source scripts/pipeline/config.sh
#   ./scripts/pipeline/steps/00_prepare_germline.sh --germline-dir DIR --output-dir DIR
#
# Output:
#   Prints the resolved germline phased VCF path to stdout.
#   All other messages go to stderr.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

GERMLINE_DIR=""
OUTPUT_DIR=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --germline-dir) GERMLINE_DIR="$2"; shift 2 ;;
        --output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
        --dry-run)      DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: $0 --germline-dir DIR --output-dir DIR [--dry-run]" >&2
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${GERMLINE_DIR}" ]] || [[ -z "${OUTPUT_DIR}" ]]; then
    log_error "Both --germline-dir and --output-dir are required."
    exit 1
fi

validate_dir "${GERMLINE_DIR}" "Germline phased VCF directory"

# ============================================================================
# Main Logic
# ============================================================================

log_info "Preparing germline phased VCF..."
log_info "  Germline dir: ${GERMLINE_DIR}"

# Strategy 1: Look for phased per-chromosome VCFs from Clair3 phase_output
# These contain actual phased genotypes (0|1 with PS tag)
PHASED_DIR="${GERMLINE_DIR}/tmp/phase_output/phase_vcf"
if [[ -d "${PHASED_DIR}" ]]; then
    PER_CHR_VCFS=()
    for chr_num in $(seq 1 22) X Y; do
        chr_vcf="${PHASED_DIR}/phased_chr${chr_num}.vcf.gz"
        if [[ -f "${chr_vcf}" ]]; then
            PER_CHR_VCFS+=("${chr_vcf}")
        fi
    done

    if [[ ${#PER_CHR_VCFS[@]} -gt 0 ]]; then
        log_info "  Found ${#PER_CHR_VCFS[@]} phased per-chromosome VCFs in ${PHASED_DIR}"

        mkdir -p "${OUTPUT_DIR}"
        OUTPUT_VCF="${OUTPUT_DIR}/germline_phased_merged.vcf.gz"

        # Check if we already merged
        if [[ -f "${OUTPUT_VCF}" ]]; then
            log_info "  Using existing merged phased VCF: ${OUTPUT_VCF}"
            echo "${OUTPUT_VCF}"
            exit 0
        fi

        if [[ "${DRY_RUN}" == true ]]; then
            echo "[DRY-RUN] ${BCFTOOLS} concat -a ${PER_CHR_VCFS[*]} -Oz -o ${OUTPUT_VCF}" >&2
            echo "[DRY-RUN] ${BCFTOOLS} index ${OUTPUT_VCF}" >&2
            echo "${OUTPUT_VCF}"
        else
            log_info "  Merging phased per-chromosome VCFs..."
            "${BCFTOOLS}" concat --naive "${PER_CHR_VCFS[@]}" -Oz -o "${OUTPUT_VCF}"
            "${BCFTOOLS}" index "${OUTPUT_VCF}"
            log_info "  Merged phased VCF created: ${OUTPUT_VCF}"
            echo "${OUTPUT_VCF}"
        fi
        exit 0
    fi
fi

# Strategy 2: Check for existing merge_output.vcf.gz (may or may not be phased)
MERGED_VCF="${GERMLINE_DIR}/merge_output.vcf.gz"
if [[ -f "${MERGED_VCF}" ]]; then
    log_info "  Found existing merged VCF: ${MERGED_VCF}"
    log_warn "  WARNING: merge_output.vcf.gz may lack phasing. Prefer phased VCFs if available."
    # Verify index exists
    if [[ ! -f "${MERGED_VCF}.tbi" ]] && [[ ! -f "${MERGED_VCF}.csi" ]]; then
        log_warn "  Index not found, creating..."
        if [[ "${DRY_RUN}" == true ]]; then
            echo "[DRY-RUN] ${BCFTOOLS} index ${MERGED_VCF}" >&2
        else
            "${BCFTOOLS}" index "${MERGED_VCF}"
        fi
    fi
    echo "${MERGED_VCF}"
    exit 0
fi

# Strategy 3: Check for pileup or full_alignment VCF
PHASED_VCF="${GERMLINE_DIR}/pileup.vcf.gz"
if [[ -f "${PHASED_VCF}" ]]; then
    log_info "  Found pileup VCF: ${PHASED_VCF}"
    echo "${PHASED_VCF}"
    exit 0
fi

PHASED_VCF="${GERMLINE_DIR}/full_alignment.vcf.gz"
if [[ -f "${PHASED_VCF}" ]]; then
    log_info "  Found full_alignment VCF: ${PHASED_VCF}"
    echo "${PHASED_VCF}"
    exit 0
fi

# Strategy 4: Check for per-chromosome VCFs directly in germline dir
PER_CHR_VCFS=()
for chr_num in $(seq 1 22) X Y; do
    chr_vcf="${GERMLINE_DIR}/chr${chr_num}.vcf.gz"
    if [[ -f "${chr_vcf}" ]]; then
        PER_CHR_VCFS+=("${chr_vcf}")
    fi
done

if [[ ${#PER_CHR_VCFS[@]} -gt 0 ]]; then
    log_info "  Found ${#PER_CHR_VCFS[@]} per-chromosome VCFs. Merging..."

    mkdir -p "${OUTPUT_DIR}"
    OUTPUT_VCF="${OUTPUT_DIR}/germline_merged.vcf.gz"

    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] ${BCFTOOLS} concat -a ${PER_CHR_VCFS[*]} -Oz -o ${OUTPUT_VCF}" >&2
        echo "[DRY-RUN] ${BCFTOOLS} index ${OUTPUT_VCF}" >&2
        echo "${OUTPUT_VCF}"
    else
        "${BCFTOOLS}" concat -a "${PER_CHR_VCFS[@]}" -Oz -o "${OUTPUT_VCF}"
        "${BCFTOOLS}" index "${OUTPUT_VCF}"
        log_info "  Merged VCF created: ${OUTPUT_VCF}"
        echo "${OUTPUT_VCF}"
    fi
    exit 0
fi

log_error "No VCF files found in ${GERMLINE_DIR}"
log_error "Searched: phase_output/phase_vcf/phased_chr*.vcf.gz, merge_output.vcf.gz, pileup.vcf.gz, chr*.vcf.gz"
exit 1

