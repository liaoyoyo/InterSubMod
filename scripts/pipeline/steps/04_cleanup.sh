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
ALLOWED_ROOT=""
SITE_PROFILE=""
SAMPLE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir)       CLEANUP_OUTPUT_DIR="$2"; shift 2 ;;
        --compress-regions) COMPRESS_REGIONS=true; shift ;;
        --dry-run)          DRY_RUN=true; shift ;;
        --allowed-root)     ALLOWED_ROOT="$2"; shift 2 ;;
        --site-profile)     SITE_PROFILE="$2"; shift 2 ;;
        --sample)           SAMPLE="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 --output-dir DIR [--compress-regions] [--dry-run]" >&2
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -n "${SITE_PROFILE}" ]]; then
    if [[ -z "${SAMPLE}" ]]; then
        log_error "--sample is required with --site-profile."
        exit 2
    fi
    load_site_profile_config "${SITE_PROFILE}" "${SAMPLE}" || exit $?
fi
ALLOWED_ROOT="${ALLOWED_ROOT:-${OUTPUT_ROOT}}"

if [[ -z "${CLEANUP_OUTPUT_DIR}" ]]; then
    log_error "--output-dir is required."
    exit 1
fi

# Destructive cleanup is allowed only for a pipeline-created, non-symlink run
# directory whose sentinel binds the exact canonical run and allowed root.
if [[ ! -d "${CLEANUP_OUTPUT_DIR}" ]] || path_has_symlink_component "${CLEANUP_OUTPUT_DIR}"; then
    log_error "Cleanup output must be an existing non-symlink directory: ${CLEANUP_OUTPUT_DIR}"
    exit 3
fi
CLEANUP_REAL="$(realpath -e -- "${CLEANUP_OUTPUT_DIR}")"
ALLOWED_REAL="$(realpath -m -- "${ALLOWED_ROOT}")"
if ! path_is_within "${CLEANUP_REAL}" "${ALLOWED_REAL}" || [[ "${CLEANUP_REAL}" == "${ALLOWED_REAL}" ]]; then
    log_error "Cleanup output is not a child of allowed root: ${CLEANUP_REAL}"
    exit 3
fi
SENTINEL="${CLEANUP_OUTPUT_DIR}/.intersubmod-run-root"
if [[ ! -f "${SENTINEL}" ]] || [[ -L "${SENTINEL}" ]]; then
    log_error "Cleanup sentinel is missing or unsafe: ${SENTINEL}"
    exit 3
fi
if ! grep -Fxq "schema_version=1" "${SENTINEL}" \
    || ! grep -Fxq "run_root=${CLEANUP_REAL}" "${SENTINEL}" \
    || ! grep -Fxq "allowed_root=${ALLOWED_REAL}" "${SENTINEL}"; then
    log_error "Cleanup sentinel does not bind the canonical run/allowed root."
    exit 3
fi

# ============================================================================
# Main Cleanup
# ============================================================================

log_info "[Step 04] Cleanup for: ${CLEANUP_OUTPUT_DIR}"
print_disk_space "${CLEANUP_REAL}"

TOTAL_FREED=0

# Remove tagged BAM and index
LONGPHASE_DIR="${CLEANUP_OUTPUT_DIR}/longphase_s"
if [[ -d "${LONGPHASE_DIR}" ]]; then
    if [[ -L "${LONGPHASE_DIR}" ]] || ! path_is_within "${LONGPHASE_DIR}" "${CLEANUP_REAL}"; then
        log_error "Unsafe LongPhase directory: ${LONGPHASE_DIR}"
        exit 3
    fi
    for bam_file in "${LONGPHASE_DIR}"/*_tagged.bam; do
        if [[ -f "${bam_file}" ]]; then
            if [[ -L "${bam_file}" ]] || ! path_is_within "${bam_file}" "${CLEANUP_REAL}"; then
                log_error "Refusing symlink or escaped BAM cleanup target: ${bam_file}"
                exit 3
            fi
            if [[ -L "${bam_file}.bai" ]]; then
                log_error "Refusing symlink BAM index cleanup target: ${bam_file}.bai"
                exit 3
            fi
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
        if [[ -L "${sc_vcf}" ]] || ! path_is_within "${sc_vcf}" "${CLEANUP_REAL}"; then
            log_error "Refusing symlink or escaped VCF cleanup target: ${sc_vcf}"
            exit 3
        fi
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
            if [[ -L "${ism_dir}" ]] || ! path_is_within "${ism_dir}" "${CLEANUP_REAL}"; then
                log_error "Unsafe InterSubMod output directory: ${ism_dir}"
                exit 3
            fi
            # Keep significance_summary.csv but compress per-region dirs
            REGION_DIRS=$(find "${ism_dir}" -maxdepth 1 -type d -name "chr*" 2>/dev/null | wc -l)
            if [[ ${REGION_DIRS} -gt 0 ]]; then
                ARCHIVE="${ism_dir}/regions.tar.gz"
                if [[ "${DRY_RUN}" == true ]]; then
                    log_info "[DRY-RUN] Would compress ${REGION_DIRS} region dirs in ${ism_dir}"
                else
                    log_info "  Compressing ${REGION_DIRS} region dirs in ${ism_dir}..."
                    if [[ -e "${ARCHIVE}" ]]; then
                        log_error "Refusing to overwrite existing region archive: ${ARCHIVE}"
                        exit 3
                    fi
                    mapfile -d '' REGION_PATHS < <(find "${ism_dir}" -mindepth 1 -maxdepth 1 -type d -name 'chr*' -print0)
                    REGION_NAMES=()
                    for region_path in "${REGION_PATHS[@]}"; do
                        if [[ -L "${region_path}" ]] || ! path_is_within "${region_path}" "${ism_dir}"; then
                            log_error "Unsafe region directory: ${region_path}"
                            exit 3
                        fi
                        REGION_NAMES+=("$(basename -- "${region_path}")")
                    done
                    (cd "${ism_dir}" && tar -czf regions.tar.gz -- "${REGION_NAMES[@]}")
                    for region_path in "${REGION_PATHS[@]}"; do
                        rm -rf -- "${region_path}"
                    done
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
    print_disk_space "${CLEANUP_REAL}"
fi

log_info "[Step 04] Cleanup complete."
