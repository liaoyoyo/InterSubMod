#!/bin/bash
# 04_cleanup.sh - Clean up large intermediate files to save disk space
#
# Removes tagged BAM and index files (typically ~264GB).
# Optionally compresses per-region InterSubMod output.
#
# Usage:
#   ./scripts/pipeline/steps/04_cleanup.sh --output-dir DIR [--compress-regions] [--dry-run]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

CLEANUP_OUTPUT_DIR=""
COMPRESS_REGIONS=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir)       CLEANUP_OUTPUT_DIR="$2"; shift 2 ;;
        --compress-regions) COMPRESS_REGIONS=true; shift ;;
        --dry-run)          DRY_RUN=true; shift ;;
        -h|--help)
            echo "Usage: $0 --output-dir DIR [--compress-regions] [--dry-run]" >&2
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "${CLEANUP_OUTPUT_DIR}" ]]; then
    log_error "--output-dir is required."
    exit 1
fi

# ============================================================================
# Main Cleanup
# ============================================================================

log_info "[Step 04] Cleanup for: ${CLEANUP_OUTPUT_DIR}"
print_disk_space

TOTAL_FREED=0

# Remove tagged BAM and index
LONGPHASE_DIR="${CLEANUP_OUTPUT_DIR}/longphase_s"
if [[ -d "${LONGPHASE_DIR}" ]]; then
    for bam_file in "${LONGPHASE_DIR}"/*_tagged.bam; do
        if [[ -f "${bam_file}" ]]; then
            BAM_SIZE=$(stat -c%s "${bam_file}" 2>/dev/null || echo 0)
            BAM_SIZE_GB=$(echo "scale=1; ${BAM_SIZE} / 1073741824" | bc)

            if [[ "${DRY_RUN}" == true ]]; then
                log_info "[DRY-RUN] Would remove: ${bam_file} (${BAM_SIZE_GB} GB)"
            else
                log_info "  Removing tagged BAM: ${bam_file} (${BAM_SIZE_GB} GB)"
                rm -f "${bam_file}"
                rm -f "${bam_file}.bai"
                if [[ ! -f "${bam_file}" ]]; then
                     log_info "    Removal successful."
                     TOTAL_FREED=$((TOTAL_FREED + BAM_SIZE))
                else
                     log_error "    Failed to remove: ${bam_file}"
                fi
            fi
        fi
    done
fi

# Remove somatic calling VCF (intermediate)
for sc_vcf in "${LONGPHASE_DIR}"/*_sc.vcf; do
    if [[ -f "${sc_vcf}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[DRY-RUN] Would remove: ${sc_vcf}"
        else
            SC_SIZE=$(stat -c%s "${sc_vcf}" 2>/dev/null || echo 0)
            log_info "  Removing somatic calling VCF: ${sc_vcf}"
            rm -f "${sc_vcf}"
            TOTAL_FREED=$((TOTAL_FREED + SC_SIZE))
        fi
    fi
done

# Optionally compress per-region output directories
if [[ "${COMPRESS_REGIONS}" == true ]]; then
    for ism_dir in "${CLEANUP_OUTPUT_DIR}"/intersubmod_tp "${CLEANUP_OUTPUT_DIR}"/intersubmod_fp; do
        if [[ -d "${ism_dir}" ]]; then
            # Keep significance_summary.csv but compress per-region dirs
            REGION_DIRS=$(find "${ism_dir}" -maxdepth 1 -type d -name "chr*" 2>/dev/null | wc -l)
            if [[ ${REGION_DIRS} -gt 0 ]]; then
                ARCHIVE="${ism_dir}/regions.tar.gz"
                if [[ "${DRY_RUN}" == true ]]; then
                    log_info "[DRY-RUN] Would compress ${REGION_DIRS} region dirs in ${ism_dir}"
                else
                    log_info "  Compressing ${REGION_DIRS} region dirs in ${ism_dir}..."
                    (cd "${ism_dir}" && tar -czf regions.tar.gz chr*/ && rm -rf chr*/)
                    log_info "  Compressed to: ${ARCHIVE}"
                fi
            fi
        fi
    done
fi

# Report
if [[ "${DRY_RUN}" != true ]]; then
    FREED_GB=$(echo "scale=1; ${TOTAL_FREED} / 1073741824" | bc)
    log_info "  Total space freed: ${FREED_GB} GB"
    print_disk_space
fi

log_info "[Step 04] Cleanup complete."
