#!/bin/bash
# run_purity_and_standard_verification.sh
#
# 目的：
# 1) 自動掃描指定 sample 的 subsample 所有 t*_n* purity 版本
# 2) 逐一執行 longphase-s -> intersubmod -> filter analysis
# 3) 產出不覆蓋的 run 目錄與總結狀態表
# 4) 額外檢查 MM/ML 標籤可用性，避免誤判「流程成功但無甲基化訊號」

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../pipeline"
PURITY_OUTPUT_ROOT="/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup"
SAMPLE="HCC1395"
BASE_DATA_DIR=""
THREADS=112
RUN_TAG="$(date +%Y%m%d_%H%M%S)"
DRY_RUN=false
ONLY_SUBDIR=""
CLEANUP_TAGGED_BAM=true
INCLUDE_ZERO_TUMOR=false
SOMATIC_VCF_REL_PATH="ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"
SAMPLE_BAM_PREFIX=""
REQUIRE_MM_ML_FOR_METHYLATION=true
MM_ML_SAMPLE_READS=1000

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

usage() {
    cat <<EOF
Usage:
  $0 [--sample NAME] [--base-data-dir DIR] [--bam-prefix PREFIX] [--somatic-vcf-rel RELPATH] [--threads N] [--run-tag TAG] [--output-root DIR] [--only-subdir tXX_nYY] [--dry-run] [--keep-tagged-bam] [--include-zero-tumor]

Options:
  --sample NAME        樣本名稱（對應 pipeline config，預設: ${SAMPLE}）
  --base-data-dir DIR  subsample 根目錄（預設依 sample 自動對應）
  --bam-prefix PREFIX  subsample BAM 檔名前綴（預設: sample 去掉後綴，例如 HCC1395_DORADO -> HCC1395）
  --somatic-vcf-rel P  各 subsample 下 somatic VCF 相對路徑（預設: ${SOMATIC_VCF_REL_PATH}）
  --threads N           執行緒數 (預設: ${THREADS})
  --run-tag TAG         自訂 run tag，避免覆蓋輸出
  --output-root DIR     輸出根目錄 (預設: ${PURITY_OUTPUT_ROOT})
  --only-subdir NAME    僅跑指定 subsample (例如 t40_n10)
  --dry-run             只印出將執行內容，不實際執行
  --keep-tagged-bam     不刪除 longphase 產生的 tagged BAM
  --include-zero-tumor  包含 t00_n*（預設會跳過）
  --allow-no-mm-ml      即使 BAM 缺 MM/ML 仍執行 InterSubMod/Filter（預設會跳過）
  --mm-ml-sample-reads N  MM/ML 抽樣讀數（預設: ${MM_ML_SAMPLE_READS}）
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --base-data-dir) BASE_DATA_DIR="$2"; shift 2 ;;
        --bam-prefix) SAMPLE_BAM_PREFIX="$2"; shift 2 ;;
        --somatic-vcf-rel) SOMATIC_VCF_REL_PATH="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --run-tag) RUN_TAG="$2"; shift 2 ;;
        --output-root) PURITY_OUTPUT_ROOT="$2"; shift 2 ;;
        --only-subdir) ONLY_SUBDIR="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --keep-tagged-bam) CLEANUP_TAGGED_BAM=false; shift ;;
        --include-zero-tumor) INCLUDE_ZERO_TUMOR=true; shift ;;
        --allow-no-mm-ml) REQUIRE_MM_ML_FOR_METHYLATION=false; shift ;;
        --mm-ml-sample-reads) MM_ML_SAMPLE_READS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

source "${PIPELINE_DIR}/config.sh"
eval "$(get_sample_config "${SAMPLE}")"

if [[ -z "${BASE_DATA_DIR}" ]]; then
    case "${SAMPLE}" in
        HCC1395)
            BASE_DATA_DIR="/big8_disk/data/HCC1395/ONT/subsample"
            ;;
        HCC1395_DORADO)
            BASE_DATA_DIR="/big8_disk/data/HCC1395/ONT_Dorado/subsample"
            ;;
        COLO829)
            BASE_DATA_DIR="/big8_disk/data/COLO829/ONT_PAO/subsample"
            ;;
        H1437)
            BASE_DATA_DIR="/big8_disk/data/H1437/ONT/subsample"
            ;;
        H2009)
            BASE_DATA_DIR="/big8_disk/data/H2009/ONT/subsample"
            ;;
        HCC1937)
            BASE_DATA_DIR="/big8_disk/data/HCC1937/ONT/subsample"
            ;;
        HCC1954)
            BASE_DATA_DIR="/big8_disk/data/HCC1954/ONT/subsample"
            ;;
        *)
            log "[ERROR] Unknown sample for auto base-data-dir: ${SAMPLE}. Please pass --base-data-dir."
            exit 1
            ;;
    esac
fi

if [[ -z "${SAMPLE_BAM_PREFIX}" ]]; then
    SAMPLE_BAM_PREFIX="${SAMPLE%%_*}"
fi

STATUS_DIR="${PURITY_OUTPUT_ROOT}/${SAMPLE}/purity_runs/${RUN_TAG}"
mkdir -p "${STATUS_DIR}"
STATUS_TSV="${STATUS_DIR}/purity_status.tsv"

if [[ ! -f "${STATUS_TSV}" ]]; then
cat > "${STATUS_TSV}" <<EOF
run_tag	subsample	estimated_purity	tumor_bam	somatic_vcf	source_mm_tags	source_ml_tags	step01	step02	step03	cleanup	tp_regions_analyzed	fp_regions_analyzed	baseline_f1	filtered_f1	f1_delta	notes	output_dir
EOF
fi

count_mm_ml() {
    local bam_file="$1"
    local sampled="${MM_ML_SAMPLE_READS}"
    local sampled_reads mm_count ml_count

    sampled_reads="$(samtools view "${bam_file}" 2>/dev/null | head -n "${sampled}" || true)"
    mm_count="$(printf '%s\n' "${sampled_reads}" | grep -c "MM:Z:" || true)"
    ml_count="$(printf '%s\n' "${sampled_reads}" | grep -c "ML:B:C" || true)"
    echo "${mm_count},${ml_count}"
}

estimate_purity_pct() {
    local subdir="$1"
    if [[ "${subdir}" =~ ^t([0-9]+)_n([0-9]+)$ ]]; then
        local t="${BASH_REMATCH[1]}"
        local n="${BASH_REMATCH[2]}"
        local total=$((t + n))
        if [[ "${total}" -eq 0 ]]; then
            echo "NA"
        else
            python3 - <<PY
t=${t}
n=${n}
print(f"{(t/(t+n))*100:.1f}")
PY
        fi
    else
        echo "NA"
    fi
}

run_filter_analysis_wrapper() {
    local run_dir="$1"
    local purity_label="$2"
    local lp_dir="${run_dir}/longphase_s"
    local tp_dir="${run_dir}/intersubmod_tp"
    local fp_dir="${run_dir}/intersubmod_fp"
    local tp_summary="${tp_dir}/significance_summary.csv"
    local fp_summary="${fp_dir}/significance_summary.csv"
    local tp_vcf="${lp_dir}/filtered_snv_tp.vcf.gz"
    local fp_vcf="${lp_dir}/filtered_snv_fp.vcf.gz"

    if [[ ! -f "${lp_dir}/variant_counts.txt" ]]; then
        log "[ERROR] variant_counts.txt not found: ${lp_dir}/variant_counts.txt"
        return 1
    fi
    # shellcheck source=/dev/null
    source "${lp_dir}/variant_counts.txt"

    if [[ ! -f "${tp_summary}" ]] || [[ ! -f "${fp_summary}" ]]; then
        log "[ERROR] significance_summary.csv missing under ${run_dir}"
        return 1
    fi

    python3 "${PIPELINE_DIR}/steps/03_filter_analysis.py" \
        --tp-summary "${tp_summary}" \
        --fp-summary "${fp_summary}" \
        --tp-count "${TP_COUNT}" \
        --fp-count "${FP_COUNT}" \
        --truth-total "${TRUTH_TOTAL}" \
        --output-dir "${run_dir}" \
        --tp-vcf-file "${tp_vcf}" \
        --fp-vcf-file "${fp_vcf}" \
        --sample "${SAMPLE}" \
        --mode "s-pure-pileup" \
        --purity "${purity_label}"
}

collect_metrics_fields() {
    local metrics_file="$1"
    python3 - "${metrics_file}" <<'PY'
import json, sys
p = sys.argv[1]
with open(p) as f:
    m = json.load(f)
tp = m.get("tp_regions_analyzed", "NA")
fp = m.get("fp_regions_analyzed", "NA")
b = m.get("baseline", {})
fl = m.get("filtered", {})
imp = m.get("improvement", {})
print("\t".join([
    str(tp),
    str(fp),
    str(b.get("f1", "NA")),
    str(fl.get("f1", "NA")),
    str(imp.get("f1_delta", "NA"))
]))
PY
}

log "=== Start purity validation run ==="
log "SAMPLE: ${SAMPLE}"
log "RUN_TAG: ${RUN_TAG}"
log "OUTPUT_ROOT: ${PURITY_OUTPUT_ROOT}"
log "THREADS: ${THREADS}"
log "BASE_DATA_DIR: ${BASE_DATA_DIR}"
log "SAMPLE_BAM_PREFIX: ${SAMPLE_BAM_PREFIX}"
log "SOMATIC_VCF_REL_PATH: ${SOMATIC_VCF_REL_PATH}"
log "STATUS_TSV: ${STATUS_TSV}"
log "CLEANUP_TAGGED_BAM: ${CLEANUP_TAGGED_BAM}"
log "INCLUDE_ZERO_TUMOR: ${INCLUDE_ZERO_TUMOR}"
log "REQUIRE_MM_ML_FOR_METHYLATION: ${REQUIRE_MM_ML_FOR_METHYLATION}"
log "MM_ML_SAMPLE_READS: ${MM_ML_SAMPLE_READS}"

mapfile -t SUBDIRS < <(find "${BASE_DATA_DIR}" -mindepth 1 -maxdepth 1 -type d -name 't*_n*' -printf '%f\n' | sort)

if [[ ${#SUBDIRS[@]} -eq 0 ]]; then
    log "[ERROR] No subsample directories found under ${BASE_DATA_DIR}"
    exit 1
fi

for subdir in "${SUBDIRS[@]}"; do
    if [[ -n "${ONLY_SUBDIR}" ]] && [[ "${subdir}" != "${ONLY_SUBDIR}" ]]; then
        continue
    fi

    if [[ "${subdir}" =~ ^t00_n[0-9]+$ ]] && [[ "${INCLUDE_ZERO_TUMOR}" == false ]]; then
        log "=== Skipping ${subdir} (zero-tumor purity) ==="
        continue
    fi

    purity_pct="$(estimate_purity_pct "${subdir}")"
    tumor_bam="${BASE_DATA_DIR}/${subdir}/${SAMPLE_BAM_PREFIX}_${subdir}.bam"
    somatic_vcf="${BASE_DATA_DIR}/${subdir}/${SOMATIC_VCF_REL_PATH}"
    output_dir="${PURITY_OUTPUT_ROOT}/${SAMPLE}/purity_${subdir}_${RUN_TAG}"

    step01="SKIP"
    step02="SKIP"
    step03="SKIP"
    cleanup_status="SKIP"
    notes=""
    source_mm="NA"
    source_ml="NA"
    tp_regions="NA"
    fp_regions="NA"
    baseline_f1="NA"
    filtered_f1="NA"
    f1_delta="NA"
    missing_mm_ml=false

    log "=== Processing ${subdir} (estimated purity=${purity_pct}%) ==="

    if [[ ! -f "${tumor_bam}" ]]; then
        notes="missing_tumor_bam"
        log "[WARN] ${notes}: ${tumor_bam}"
        echo -e "${RUN_TAG}\t${subdir}\t${purity_pct}\t${tumor_bam}\t${somatic_vcf}\t${source_mm}\t${source_ml}\t${step01}\t${step02}\t${step03}\t${cleanup_status}\t${tp_regions}\t${fp_regions}\t${baseline_f1}\t${filtered_f1}\t${f1_delta}\t${notes}\t${output_dir}" >> "${STATUS_TSV}"
        continue
    fi

    if [[ ! -f "${somatic_vcf}" ]]; then
        notes="missing_somatic_vcf"
        log "[WARN] ${notes}: ${somatic_vcf}"
        echo -e "${RUN_TAG}\t${subdir}\t${purity_pct}\t${tumor_bam}\t${somatic_vcf}\t${source_mm}\t${source_ml}\t${step01}\t${step02}\t${step03}\t${cleanup_status}\t${tp_regions}\t${fp_regions}\t${baseline_f1}\t${filtered_f1}\t${f1_delta}\t${notes}\t${output_dir}" >> "${STATUS_TSV}"
        continue
    fi

    mm_ml="$(count_mm_ml "${tumor_bam}")"
    source_mm="${mm_ml%,*}"
    source_ml="${mm_ml#*,}"
    if [[ "${source_mm}" -eq 0 ]] || [[ "${source_ml}" -eq 0 ]]; then
        missing_mm_ml=true
        notes="source_bam_missing_mm_ml"
        log "[WARN] ${notes}: ${tumor_bam} (MM=${source_mm}, ML=${source_ml})"
    fi

    mkdir -p "${output_dir}"

    if [[ "${DRY_RUN}" == true ]]; then
        log "[DRY-RUN] would run steps for ${subdir} into ${output_dir}"
        step01="DRY_RUN"
        step02="DRY_RUN"
        step03="DRY_RUN"
    else
        if "${PIPELINE_DIR}/steps/01_longphase_s.sh" \
            --sample "${SAMPLE}" \
            --output-dir "${output_dir}" \
            --threads "${THREADS}" \
            --somatic-vcf "${somatic_vcf}" \
            --tumor-bam "${tumor_bam}"; then
            step01="OK"
        else
            step01="FAIL"
            notes="${notes};step01_failed"
        fi

        if [[ "${step01}" == "OK" ]]; then
            if [[ "${missing_mm_ml}" == true ]] && [[ "${REQUIRE_MM_ML_FOR_METHYLATION}" == true ]]; then
                step02="SKIP_MMML"
                step03="SKIP_MMML"
                notes="${notes};skip_methylation_missing_mm_ml"
                log "[WARN] Skip InterSubMod/Filter for ${subdir}: missing MM/ML tags in source BAM."
            fi
        fi

        if [[ "${step01}" == "OK" ]] && [[ "${step02}" != "SKIP_MMML" ]]; then
            tagged_bam="${output_dir}/longphase_s/${SAMPLE}_tagged.bam"
            tp_vcf="${output_dir}/longphase_s/filtered_snv_tp.vcf.gz"
            fp_vcf="${output_dir}/longphase_s/filtered_snv_fp.vcf.gz"
            if "${PIPELINE_DIR}/steps/02_intersubmod.sh" \
                --tagged-bam "${tagged_bam}" \
                --tp-vcf "${tp_vcf}" \
                --fp-vcf "${fp_vcf}" \
                --output-dir "${output_dir}" \
                --threads "${THREADS}"; then
                step02="OK"
            else
                step02="FAIL"
                notes="${notes};step02_failed"
            fi
        fi

        if [[ "${step02}" == "OK" ]]; then
            if run_filter_analysis_wrapper "${output_dir}" "${purity_pct}"; then
                step03="OK"
            else
                step03="FAIL"
                notes="${notes};step03_failed"
            fi
        fi
    fi

    metrics_file="${output_dir}/metrics.json"
    if [[ -f "${metrics_file}" ]]; then
        metrics_fields="$(collect_metrics_fields "${metrics_file}")"
        tp_regions="$(echo "${metrics_fields}" | cut -f1)"
        fp_regions="$(echo "${metrics_fields}" | cut -f2)"
        baseline_f1="$(echo "${metrics_fields}" | cut -f3)"
        filtered_f1="$(echo "${metrics_fields}" | cut -f4)"
        f1_delta="$(echo "${metrics_fields}" | cut -f5)"
        if [[ "${tp_regions}" == "0" ]] && [[ "${fp_regions}" == "0" ]]; then
            notes="${notes};no_intersubmod_regions_analyzed"
        fi
    elif [[ "${step03}" == "SKIP_MMML" ]] || [[ "${step03}" == "DRY_RUN" ]] || [[ "${step03}" == "SKIP" ]]; then
        true
    else
        notes="${notes};metrics_missing"
    fi

    tagged_bam="${output_dir}/longphase_s/${SAMPLE}_tagged.bam"
    tagged_bam_index="${tagged_bam}.bai"
    if [[ "${DRY_RUN}" == true ]]; then
        cleanup_status="DRY_RUN"
    elif [[ "${CLEANUP_TAGGED_BAM}" == true ]] && [[ "${step02}" != "FAIL" ]] && [[ -f "${tagged_bam}" ]]; then
        rm -f "${tagged_bam}" "${tagged_bam_index}"
        cleanup_status="REMOVED"
    elif [[ "${CLEANUP_TAGGED_BAM}" == false ]]; then
        cleanup_status="KEPT"
    else
        cleanup_status="SKIP"
    fi

    echo -e "${RUN_TAG}\t${subdir}\t${purity_pct}\t${tumor_bam}\t${somatic_vcf}\t${source_mm}\t${source_ml}\t${step01}\t${step02}\t${step03}\t${cleanup_status}\t${tp_regions}\t${fp_regions}\t${baseline_f1}\t${filtered_f1}\t${f1_delta}\t${notes}\t${output_dir}" >> "${STATUS_TSV}"
done

python3 "${SCRIPT_DIR}/build_purity_aware_tables.py" \
    --status-tsv "${STATUS_TSV}" \
    --sample "${SAMPLE}" \
    --output-dir "${STATUS_DIR}"

log "All tasks complete. Status table: ${STATUS_TSV}"
