#!/bin/bash
# config.sh - Global configuration for benchmark pipeline
#
# All paths, parameters, and sample definitions are centralized here.
# Source this file from other pipeline scripts.

# ============================================================================
# Tool Paths
# ============================================================================

PROJECT_ROOT_DEFAULT="/big7_disk/liaoyoyo2001/InterSubMod"
LONGPHASE_S_BIN="${LONGPHASE_S_BIN:-/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s}"
INTERSUBMOD_BIN="${INTERSUBMOD_BIN:-${PROJECT_ROOT_DEFAULT}/build/bin/inter_sub_mod}"
INTERSUBMOD_BIN_FALLBACK="/big8_disk/liaoyoyo2001/Knowledge/codebase/InterSubMod/build/bin/inter_sub_mod"
REFERENCE="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
SAMTOOLS="/usr/local/bin/samtools"
BCFTOOLS="/usr/bin/bcftools"

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
