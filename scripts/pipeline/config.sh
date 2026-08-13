#!/bin/bash
# config.sh - Global configuration for benchmark pipeline
#
# All paths, parameters, and sample definitions are centralized here.
# Source this file from other pipeline scripts.

# ============================================================================
# Tool Paths
# ============================================================================

PROJECT_ROOT_DEFAULT="${INTERSUBMOD_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
SITE_PROFILE_ACTIVE="${SITE_PROFILE_ACTIVE:-false}"
LONGPHASE_S_BIN="${LONGPHASE_S_BIN:-/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s}"
INTERSUBMOD_BIN="${INTERSUBMOD_BIN:-${PROJECT_ROOT_DEFAULT}/build/bin/inter_sub_mod}"
INTERSUBMOD_BIN_FALLBACK="${INTERSUBMOD_BIN_FALLBACK:-/big8_disk/liaoyoyo2001/Knowledge/codebase/InterSubMod/build/bin/inter_sub_mod}"
REFERENCE="${REFERENCE:-/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta}"
SAMTOOLS="${SAMTOOLS:-/usr/local/bin/samtools}"
BCFTOOLS="${BCFTOOLS:-/usr/bin/bcftools}"

# ============================================================================
# Execution Parameters
# ============================================================================

if command -v nproc >/dev/null; then
    THREADS=$(nproc)
else
    THREADS=64
fi
WINDOW_SIZE=5000
DISTANCE_METRIC="BERNOULLI"
LOG_LEVEL="info"
PLOT_TYPE="no"
MIN_FREE_GB_DEFAULT="${MIN_FREE_GB_DEFAULT:-800}"
MAX_PARALLEL_DEFAULT="${MAX_PARALLEL_DEFAULT:-4}"
INDEX_THREADS_CAP_DEFAULT="${INDEX_THREADS_CAP_DEFAULT:-8}"

# ============================================================================
# Output Root Directory
# ============================================================================

BIG7_OUTPUT_ROOT="${BIG7_OUTPUT_ROOT:-/big7_disk/liaoyoyo2001/big7_disk_output}"
CANONICAL_OUTPUT_ROOT="${CANONICAL_OUTPUT_ROOT:-${BIG7_OUTPUT_ROOT}/canonical}"
SYNTHESIS_OUTPUT_ROOT="${SYNTHESIS_OUTPUT_ROOT:-${BIG7_OUTPUT_ROOT}/synthesis}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${CANONICAL_OUTPUT_ROOT}}"

# ============================================================================
# Sample Configurations
# ============================================================================

# HCC1395 sample (uses ONT_5khz BAM which has MM/ML methylation tags)
HCC1395_TUMOR_BAM="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam"
HCC1395_NORMAL_BAM="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam"
HCC1395_SOMATIC_VCF="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz"
HCC1395_GERMLINE_PHASED_DIR="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/tmp/clair3_output/clair3_normal_output/"
HCC1395_TRUTH_VCF="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
HCC1395_TRUTH_BED="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"

# ============================================================================
# Helper Functions
# ============================================================================

# Resolve InterSubMod binary (prefer project build, fallback to Knowledge)
resolve_intersubmod_bin() {
    if [[ -x "${INTERSUBMOD_BIN}" ]]; then
        echo "${INTERSUBMOD_BIN}"
    elif [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then
        echo "[ERROR] Profile-declared InterSubMod binary is unavailable: ${INTERSUBMOD_BIN}" >&2
        return 1
    elif [[ -x "${INTERSUBMOD_BIN_FALLBACK}" ]]; then
        echo "[CONFIG] Using fallback InterSubMod binary: ${INTERSUBMOD_BIN_FALLBACK}" >&2
        echo "${INTERSUBMOD_BIN_FALLBACK}"
    else
        echo "[ERROR] No InterSubMod binary found." >&2
        echo "[ERROR]   Primary: ${INTERSUBMOD_BIN}" >&2
        echo "[ERROR]   Fallback: ${INTERSUBMOD_BIN_FALLBACK}" >&2
        return 1
    fi
}

# Load the complete profile as a standalone authority. In profile mode no
# built-in sample or /big* fallback is consulted for missing values.
load_site_profile_config() {
    local profile="$1"
    local sample="$2"
    local loader="${PROJECT_ROOT_DEFAULT}/scripts/site/site_profile.py"
    local runtime_project_root
    runtime_project_root="$(realpath -e -- "${PROJECT_ROOT_DEFAULT}")" || return 3
    if [[ ! -f "${profile}" ]]; then
        log_error "Site profile not found: ${profile}"
        return 3
    fi
    if [[ ! -f "${loader}" ]]; then
        log_error "Site profile loader not found: ${loader}"
        return 3
    fi
    local profile_sha_before profile_sha_after
    profile_sha_before="$(sha256sum -- "${profile}" | awk '{print $1}')" || return 3
    local profile_payload assignments
    profile_payload="$(python3 - "${loader}" "${profile}" "${sample}" <<'PY'
import importlib.util
import json
import sys
from pathlib import Path

loader_path = Path(sys.argv[1])
profile_path = Path(sys.argv[2])
sample = sys.argv[3]
spec = importlib.util.spec_from_file_location("site_profile_runtime", loader_path)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)
profile = module.load_profile(profile_path)
dataset_id = module.dataset_id_for_sample(profile, sample)
receipt = module.inspect_profile(
    profile,
    include_real_data=True,
    profile_path=profile_path,
    dataset_ids={dataset_id},
)
print(json.dumps({
    "shell": module.shell_assignments(profile, sample),
    "preflight": receipt,
}, ensure_ascii=False, sort_keys=True))
PY
    )" || return $?
    profile_sha_after="$(sha256sum -- "${profile}" | awk '{print $1}')" || return 3
    if [[ "${profile_sha_before}" != "${profile_sha_after}" ]]; then
        log_error "Site profile changed while it was being parsed; refusing mixed-authority configuration."
        return 3
    fi
    assignments="$(python3 -c 'import json,sys; print(json.load(sys.stdin)["shell"])' <<< "${profile_payload}")" || return 3
    SITE_PROFILE_PREFLIGHT_JSON="$(python3 -c 'import json,sys; print(json.dumps(json.load(sys.stdin)["preflight"], ensure_ascii=False, sort_keys=True))' <<< "${profile_payload}")" || return 3
    eval "${assignments}"
    local declared_project_root
    declared_project_root="$(realpath -e -- "${PROJECT_ROOT_DEFAULT}")" || {
        log_error "Profile project_root does not exist: ${PROJECT_ROOT_DEFAULT}"
        return 3
    }
    if [[ "${declared_project_root}" != "${runtime_project_root}" ]]; then
        log_error "Profile project_root does not match the executing clone: declared=${declared_project_root} runtime=${runtime_project_root}"
        return 3
    fi
    SITE_PROFILE_SOURCE="$(realpath -e -- "${profile}")" || return 3
    SITE_PROFILE_SHA256="${profile_sha_after}"
    SITE_PROFILE_ACTIVE=true
    export SITE_PROFILE_ACTIVE SITE_PROFILE_SOURCE SITE_PROFILE_SHA256 SITE_PROFILE_PREFLIGHT_JSON
    export PROJECT_ROOT_DEFAULT REFERENCE REFERENCE_GENOME_BUILD CONTIG_NAMING CONTIG_SCOPE
    export SAMTOOLS BCFTOOLS LONGPHASE_S_BIN INTERSUBMOD_BIN
    export BIG7_OUTPUT_ROOT CANONICAL_OUTPUT_ROOT OUTPUT_ROOT
    export PLATFORM_LABEL TRUTH_SET_LABEL TRUTH_TOTAL TUMOR_BAM TUMOR_BAM_INDEX
    export NORMAL_BAM NORMAL_BAM_INDEX SOMATIC_VCF SOMATIC_VCF_INDEX
    export SOMATIC_VCF_PILEUP SOMATIC_VCF_PILEUP_INDEX SOMATIC_VCF_INDEL SOMATIC_VCF_INDEL_INDEX
    export TO_SOMATIC_VCF TO_SOMATIC_VCF_INDEX TO_SOMATIC_VCF_PILEUP TO_SOMATIC_VCF_PILEUP_INDEX
    export TO_SOMATIC_VCF_INDEL TO_SOMATIC_VCF_INDEL_INDEX
    export GERMLINE_PHASED_DIR TRUTH_VCF TRUTH_VCF_INDEX TRUTH_BED
}

# Child steps in a portable run consume the already-resolved environment. They
# verify the run-local immutable profile instead of re-reading the operator's
# mutable source file.
verify_parent_profile_lock() {
    if [[ "${SITE_PROFILE_PARENT_LOCKED:-false}" != true ]]; then
        return 0
    fi
    if [[ -z "${SITE_PROFILE_LOCK_PATH:-}" ]] || [[ -z "${SITE_PROFILE_SHA256:-}" ]]; then
        log_error "Profile parent lock metadata is incomplete."
        return 3
    fi
    if [[ ! -f "${SITE_PROFILE_LOCK_PATH}" ]]; then
        log_error "Locked site profile is missing: ${SITE_PROFILE_LOCK_PATH}"
        return 3
    fi
    local locked_sha
    locked_sha="$(sha256sum -- "${SITE_PROFILE_LOCK_PATH}" | awk '{print $1}')" || return 3
    if [[ "${locked_sha}" != "${SITE_PROFILE_SHA256}" ]]; then
        log_error "Locked site profile SHA-256 mismatch: expected ${SITE_PROFILE_SHA256}, got ${locked_sha}"
        return 3
    fi
}

validate_expected_sha256() {
    local value="$1"
    local label="$2"
    if [[ ! "${value}" =~ ^[0-9a-fA-F]{64}$ ]]; then
        log_error "${label} requires a 64-character SHA-256 value."
        return 3
    fi
}

resolve_governed_existing_artifact() {
    local requested_path="$1"
    local expected_sha="$2"
    local label="$3"
    validate_expected_sha256 "${expected_sha}" "${label}" || return $?
    if [[ -z "${requested_path}" ]]; then
        log_error "${label} path is required in profile mode with --skip-longphase."
        return 3
    fi
    local resolved
    resolved="$(realpath -e -- "${requested_path}")" || {
        log_error "${label} not found: ${requested_path}"
        return 3
    }
    if [[ ! -f "${resolved}" ]] || [[ ! -s "${resolved}" ]]; then
        log_error "${label} must be a non-empty regular file: ${resolved}"
        return 3
    fi
    echo "${resolved}"
}

verify_artifact_sha256() {
    local path="$1"
    local expected_sha="$2"
    local label="$3"
    local actual_sha
    actual_sha="$(sha256sum -- "${path}" | awk '{print $1}')" || return 3
    if [[ "${actual_sha}" != "${expected_sha,,}" ]]; then
        log_error "${label} SHA-256 mismatch: expected ${expected_sha,,}, got ${actual_sha}"
        return 3
    fi
}

is_safe_run_component() {
    [[ "$1" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$ ]]
}

path_is_within() {
    local candidate root
    candidate="$(realpath -m -- "$1")" || return 1
    root="$(realpath -m -- "$2")" || return 1
    [[ "${candidate}" == "${root}" || "${candidate}" == "${root}/"* ]]
}

path_has_symlink_component() {
    local path="$1"
    if [[ "${path}" != /* ]]; then
        path="${PWD}/${path}"
    fi
    local current="/"
    local component
    local -a components=()
    IFS='/' read -r -a components <<< "${path}"
    for component in "${components[@]}"; do
        [[ -z "${component}" ]] && continue
        current="${current%/}/${component}"
        [[ -L "${current}" ]] && return 0
    done
    return 1
}

# Get sample configuration by name
# Usage: eval "$(get_sample_config HCC1395)"
get_sample_config() {
    local sample="$1"
    case "${sample}" in
        HCC1395)
            cat <<'SAMPLE_EOF'
PLATFORM_LABEL="ONT_5kHz"
TRUTH_SET_LABEL="SEQC2_v1.2.1_HC_SNV"
TUMOR_BAM="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam"
NORMAL_BAM="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/indel.vcf.gz"
TO_SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
TRUTH_TOTAL=39447
SAMPLE_EOF
            ;;
        HCC1395_DORADO)
            cat <<'SAMPLE_EOF'
PLATFORM_LABEL="ONT_Dorado"
TRUTH_SET_LABEL="SEQC2_v1.2.1_HC_SNV"
TUMOR_BAM="/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam"
NORMAL_BAM="/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/indel.vcf.gz"
TO_SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/snv.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
TRUTH_TOTAL=39447
SAMPLE_EOF
            ;;
        COLO829)
            cat <<'SAMPLE_EOF'
PLATFORM_LABEL="ONT_PAO"
TRUTH_SET_LABEL="NYGC"
TUMOR_BAM="/big8_disk/data/COLO829/ONT_PAO/PAO29420.bam"
NORMAL_BAM="/big8_disk/data/COLO829/ONT_PAO/PAO33946.bam"
SOMATIC_VCF="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/indel.vcf.gz"
TO_SOMATIC_VCF_PILEUP="/big8_disk/data/COLO829/ONT_PAO/ClairS_TO_v0_3_0/snv.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz"
TRUTH_BED=""
TRUTH_TOTAL=41427
SAMPLE_EOF
            ;;
        H1437)
            cat <<'SAMPLE_EOF'
PLATFORM_LABEL="ONT"
TRUTH_SET_LABEL="orthogonal-tools"
TUMOR_BAM="/big8_disk/data/H1437/ONT/H1437.bam"
NORMAL_BAM="/big8_disk/data/H1437/ONT/H1437BL.bam"
SOMATIC_VCF="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/indel.vcf.gz"
TO_SOMATIC_VCF_PILEUP="/big8_disk/data/H1437/ONT/ClairS_TO_v0_3_0/snv.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark_somatic-only.vcf.gz"
TRUTH_BED="/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark.bed"
TRUTH_TOTAL=90129
SAMPLE_EOF
            ;;
        H2009)
            cat <<'SAMPLE_EOF'
PLATFORM_LABEL="ONT"
TRUTH_SET_LABEL="orthogonal-tools"
TUMOR_BAM="/big8_disk/data/H2009/ONT/H2009.bam"
NORMAL_BAM="/big8_disk/data/H2009/ONT/H2009BL.bam"
SOMATIC_VCF="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/indel.vcf.gz"
TO_SOMATIC_VCF_PILEUP="/big8_disk/data/H2009/ONT/ClairS_TO_v0_3_0/snv.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark_somatic-only.vcf.gz"
TRUTH_BED="/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark.bed"
TRUTH_TOTAL=168529
SAMPLE_EOF
            ;;
        HCC1937)
            cat <<'SAMPLE_EOF'
PLATFORM_LABEL="ONT"
TRUTH_SET_LABEL="orthogonal-tools"
TUMOR_BAM="/big8_disk/data/HCC1937/ONT/HCC1937.bam"
NORMAL_BAM="/big8_disk/data/HCC1937/ONT/HCC1937BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/indel.vcf.gz"
TO_SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1937/ONT/ClairS_TO_v0_3_0/snv.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark_somatic-only.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark.bed"
TRUTH_TOTAL=60691
SAMPLE_EOF
            ;;
        HCC1954)
            cat <<'SAMPLE_EOF'
PLATFORM_LABEL="ONT"
TRUTH_SET_LABEL="orthogonal-tools"
TUMOR_BAM="/big8_disk/data/HCC1954/ONT/HCC1954.bam"
NORMAL_BAM="/big8_disk/data/HCC1954/ONT/HCC1954BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/indel.vcf.gz"
TO_SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1954/ONT/ClairS_TO_v0_3_0/snv.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark_somatic-only.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark.bed"
TRUTH_TOTAL=26567
SAMPLE_EOF
            ;;
        *)
            echo "[ERROR] Unknown sample: ${sample}" >&2
            return 1
            ;;
    esac
}

canonical_mode_name() {
    local mode="${1:-}"
    case "${mode}" in
        s-pure) echo "paired_full" ;;
        s-pure-pileup) echo "paired_pileup" ;;
        to-pure|to_pure) echo "to_pileup" ;;
        to-full|to_full) echo "to_full" ;;
        *) echo "${mode}" ;;
    esac
}

caller_model_name() {
    local source="${1:-output}"
    case "${source}" in
        pileup) echo "pileup" ;;
        indel) echo "indel" ;;
        output|full|"") echo "full" ;;
        *) echo "${source}" ;;
    esac
}

validate_benchmark_mode() {
    case "${1:-}" in
        s-pure|s-pure-pileup|to-pure|to_pure|to-full|to_full) return 0 ;;
        *) log_error "Unsupported --mode: ${1:-<empty>} (expected s-pure, s-pure-pileup, to-pure, or to-full)"; return 2 ;;
    esac
}

validate_vcf_source() {
    case "${1:-}" in
        output|pileup|indel) return 0 ;;
        *) log_error "Unsupported --vcf-source: ${1:-<empty>} (expected output, pileup, or indel)"; return 2 ;;
    esac
}

select_somatic_vcf() {
    local mode="$1"
    local source="$2"
    validate_benchmark_mode "${mode}" || return $?
    validate_vcf_source "${source}" || return $?
    case "${mode}:${source}" in
        s-pure:output|s-pure-pileup:output) echo "${SOMATIC_VCF}" ;;
        s-pure:pileup|s-pure-pileup:pileup) echo "${SOMATIC_VCF_PILEUP}" ;;
        s-pure:indel|s-pure-pileup:indel) echo "${SOMATIC_VCF_INDEL}" ;;
        to-pure:output|to_pure:output|to-full:output|to_full:output)
            if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then echo "${TO_SOMATIC_VCF}"; else echo "${TO_SOMATIC_VCF:-${SOMATIC_VCF}}"; fi ;;
        to-pure:pileup|to_pure:pileup|to-full:pileup|to_full:pileup) echo "${TO_SOMATIC_VCF_PILEUP}" ;;
        to-pure:indel|to_pure:indel|to-full:indel|to_full:indel)
            if [[ "${SITE_PROFILE_ACTIVE}" == true ]]; then echo "${TO_SOMATIC_VCF_INDEL}"; else echo "${TO_SOMATIC_VCF_INDEL:-${SOMATIC_VCF_INDEL}}"; fi ;;
        *) return 2 ;;
    esac
}

build_canonical_run_base() {
    local sample="$1"
    local canonical_mode="$2"
    local caller_model="$3"
    local date_str="${4:-$(date +%Y%m%d)}"
    local run_tag="${5:-}"
    local run_id="${date_str}_${sample}_${canonical_mode}_${caller_model}"
    if [[ -n "${run_tag}" ]]; then
        run_id="${run_id}_${run_tag}"
    fi
    echo "${OUTPUT_ROOT}/${sample}/${canonical_mode}/${run_id}"
}

latest_canonical_run_dir() {
    local sample="$1"
    local canonical_mode="$2"
    local sample_root="${OUTPUT_ROOT}/${sample}/${canonical_mode}"
    if [[ ! -d "${sample_root}" ]]; then
        return 0
    fi
    find "${sample_root}" -mindepth 1 -maxdepth 1 -type d | sort | tail -n 1
}

find_canonical_artifact() {
    local sample="$1"
    local canonical_mode="$2"
    local relative_path="$3"
    local latest_run
    latest_run="$(latest_canonical_run_dir "${sample}" "${canonical_mode}")"
    if [[ -n "${latest_run}" ]] && [[ -f "${latest_run}/${relative_path}" ]]; then
        echo "${latest_run}/${relative_path}"
    fi
}

# Validate that a file exists
validate_file() {
    local filepath="$1"
    local description="$2"
    if [[ ! -f "${filepath}" ]]; then
        echo "[ERROR] ${description} not found: ${filepath}" >&2
        return 1
    fi
}

# Validate that a directory exists
validate_dir() {
    local dirpath="$1"
    local description="$2"
    if [[ ! -d "${dirpath}" ]]; then
        echo "[ERROR] ${description} not found: ${dirpath}" >&2
        return 1
    fi
}

# Log with timestamp
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*" >&2
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*" >&2
}

recommended_index_threads() {
    local worker_threads="${1:-${THREADS}}"
    local cap="${2:-${INDEX_THREADS_CAP_DEFAULT}}"
    local floor=4

    if (( worker_threads < 1 )); then
        worker_threads=1
    fi

    local recommended=$(( worker_threads / 2 ))
    if (( recommended < floor )); then
        recommended="${floor}"
    fi
    if (( recommended > cap )); then
        recommended="${cap}"
    fi
    if (( recommended > worker_threads )); then
        recommended="${worker_threads}"
    fi

    echo "${recommended}"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

# Print disk space info
print_disk_space() {
    local path="${1:-/big7_disk}"
    local available
    available=$(df -h "${path}" | awk 'NR==2 {print $4}')
    log_info "Available disk space on ${path}: ${available}"
}

free_disk_gb() {
    local path="${1:-/big7_disk}"
    df -BG "${path}" | awk 'NR==2 {gsub(/G/, "", $4); print $4}'
}

require_min_free_gb() {
    local path="${1:-/big7_disk}"
    local min_free_gb="${2:-${MIN_FREE_GB_DEFAULT}}"
    local label="${3:-job}"
    local available_gb

    available_gb="$(free_disk_gb "${path}")"
    if [[ -z "${available_gb}" ]]; then
        log_warn "[Guard] Unable to determine free disk space for ${path}; skip threshold check."
        return 0
    fi

    if (( available_gb < min_free_gb )); then
        log_error "[Guard] ${label} requires >= ${min_free_gb}G free on ${path}, only ${available_gb}G available."
        return 1
    fi

    log_info "[Guard] Disk threshold passed for ${label}: ${available_gb}G free on ${path} (threshold=${min_free_gb}G)."
}

# Sample first N alignments and count MM/ML methylation tags.
# Output format: "MM_COUNT,ML_COUNT"
count_mm_ml_tags() {
    local bam_path="$1"
    local sample_n="${2:-1000}"
    local samtools_bin="${SAMTOOLS:-samtools}"
    local sampled_reads mm_count ml_count

    # Suppress expected SIGPIPE noise from early-stop sampling.
    sampled_reads="$("${samtools_bin}" view "${bam_path}" 2>/dev/null | head -n "${sample_n}" || true)"
    mm_count="$(printf '%s\n' "${sampled_reads}" | grep -c "MM:Z:" || true)"
    ml_count="$(printf '%s\n' "${sampled_reads}" | grep -c "ML:B:C" || true)"

    echo "${mm_count},${ml_count}"
}

# Return success if both MM and ML tags are present in sampled reads.
has_mm_ml_tags() {
    local bam_path="$1"
    local sample_n="${2:-1000}"
    local counts mm_count ml_count

    counts="$(count_mm_ml_tags "${bam_path}" "${sample_n}")"
    mm_count="${counts%,*}"
    ml_count="${counts#*,}"

    [[ "${mm_count}" -gt 0 && "${ml_count}" -gt 0 ]]
}
