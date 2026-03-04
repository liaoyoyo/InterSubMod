#!/bin/bash
# config.sh - Global configuration for benchmark pipeline
#
# All paths, parameters, and sample definitions are centralized here.
# Source this file from other pipeline scripts.

# ============================================================================
# Tool Paths
# ============================================================================

LONGPHASE_S_BIN="/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s"
INTERSUBMOD_BIN="/big8_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod"
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

# ============================================================================
# Output Root Directory
# ============================================================================

OUTPUT_ROOT="/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output"

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
TUMOR_BAM="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam"
NORMAL_BAM="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/indel.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
TRUTH_TOTAL=39447
SAMPLE_EOF
            ;;
        HCC1395_DORADO)
            cat <<'SAMPLE_EOF'
TUMOR_BAM="/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam"
NORMAL_BAM="/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/indel.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
TRUTH_TOTAL=39447
SAMPLE_EOF
            ;;
        COLO829)
            cat <<'SAMPLE_EOF'
TUMOR_BAM="/big8_disk/data/COLO829/ONT_PAO/PAO29420.bam"
NORMAL_BAM="/big8_disk/data/COLO829/ONT_PAO/PAO33946.bam"
SOMATIC_VCF="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/indel.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz"
TRUTH_BED=""
TRUTH_TOTAL=41427
SAMPLE_EOF
            ;;
        H1437)
            cat <<'SAMPLE_EOF'
TUMOR_BAM="/big8_disk/data/H1437/ONT/H1437.bam"
NORMAL_BAM="/big8_disk/data/H1437/ONT/H1437BL.bam"
SOMATIC_VCF="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/indel.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark_somatic-only.vcf.gz"
TRUTH_BED="/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark.bed"
TRUTH_TOTAL=90129
SAMPLE_EOF
            ;;
        H2009)
            cat <<'SAMPLE_EOF'
TUMOR_BAM="/big8_disk/data/H2009/ONT/H2009.bam"
NORMAL_BAM="/big8_disk/data/H2009/ONT/H2009BL.bam"
SOMATIC_VCF="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/indel.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark_somatic-only.vcf.gz"
TRUTH_BED="/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark.bed"
TRUTH_TOTAL=168529
SAMPLE_EOF
            ;;
        HCC1937)
            cat <<'SAMPLE_EOF'
TUMOR_BAM="/big8_disk/data/HCC1937/ONT/HCC1937.bam"
NORMAL_BAM="/big8_disk/data/HCC1937/ONT/HCC1937BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/indel.vcf.gz"
GERMLINE_PHASED_DIR="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/tmp/clair3_output/clair3_normal_output/"
TRUTH_VCF="/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark_somatic-only.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark.bed"
TRUTH_TOTAL=60691
SAMPLE_EOF
            ;;
        HCC1954)
            cat <<'SAMPLE_EOF'
TUMOR_BAM="/big8_disk/data/HCC1954/ONT/HCC1954.bam"
NORMAL_BAM="/big8_disk/data/HCC1954/ONT/HCC1954BL.bam"
SOMATIC_VCF="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/output.vcf.gz"
SOMATIC_VCF_PILEUP="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SOMATIC_VCF_INDEL="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/indel.vcf.gz"
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

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

# Print disk space info
print_disk_space() {
    local path="${1:-/big8_disk}"
    local available
    available=$(df -h "${path}" | awk 'NR==2 {print $4}')
    log_info "Available disk space on ${path}: ${available}"
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
