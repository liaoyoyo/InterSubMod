#!/bin/bash
# run_longphase_to_intersubmod_pilot.sh - Tumor-only pilot:
#   ClairS-TO -> benchmark baseline -> LongPhase-TO phase -> benchmark ->
#   LongPhase-TO haplotag -> InterSubMod -> filter analysis -> round summary

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

LONGPHASE_TO_BIN="/big8_disk/liaoyoyo2001/longphase-to/longphase-to"
ARCHIVE_PENDING_DELETE="/big7_disk/liaoyoyo2001/big7_disk_output/archive_pending_delete"

SAMPLE="HCC1395"
PLATFORM="ONT_5kHz"
ANALYSIS_MODE="to-pure"
TUMOR_BAM="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam"
REFERENCE_FASTA="${REFERENCE}"
TRUTH_VCF="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_BED="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
TRUTH_TOTAL="39447"
TRUTH_SET="SEQC2_v1.2.1_HC_SNV"
CALLER_NAME="ClairS-TO"
CALLER_MODE="clairs_to_ssrs"
CALLER_MODEL="ont_r10_dorado_sup_5khz_ssrs"
THREADS_TO_USE="${THREADS}"
INDEX_THREADS_TO_USE=""
MIN_FREE_GB="${MIN_FREE_GB_DEFAULT}"
RUN_CALLER=true
SKIP_PHASE=false
SKIP_HAPLOTAG=false
SKIP_INTERSUBMOD=false
SKIP_FILTER_ANALYSIS=false
ARCHIVE_TAGGED_BAM=false
TAG_PREFIX_OVERRIDE=""
DRY_RUN=false
TEST_FOCUS="HCC1395 5kHz tumor-only LongPhase-TO + InterSubMod pilot"
CHANGES="建立 tumor-only pilot，驗證 5kHz TO 是否可重現或強化低 VAF + 高 AlleleDelta / label-first 觀察"
CONCLUSION="待正式執行完成後補寫"
NEXT_STEP="若 5kHz TO 接線成功，完成比較後將 tagged BAM 移至 Archive_pending_delete，再切 HCC1395_DORADO"
OUTPUT_ROOT="${SYNTHESIS_OUTPUT_ROOT}/research_rounds"
CALLER_EXISTING_VCF=""

USER_SET_PLATFORM=false
USER_SET_TUMOR_BAM=false
USER_SET_TRUTH_VCF=false
USER_SET_TRUTH_BED=false
USER_SET_TRUTH_TOTAL=false
USER_SET_TRUTH_SET=false
USER_SET_CALLER_EXISTING_VCF=false

PON_FILE_DEFAULT="/big8_disk/data/PON/clairs-to_databases/1000g-pon.sites.vcf.gz,/big8_disk/data/PON/clairs-to_databases/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz"
STRICT_PON_FILE_DEFAULT="/big8_disk/data/PON/clairs-to_databases/dbsnp.b138.non-somatic.sites.vcf.gz,/big8_disk/data/PON/clairs-to_databases/gnomad.r2.1.af-ge-0.001.sites.vcf.gz"
PON_FILE="${PON_FILE_DEFAULT}"
STRICT_PON_FILE="${STRICT_PON_FILE_DEFAULT}"

show_help() {
    cat <<'HELP'
Usage: run_longphase_to_intersubmod_pilot.sh [OPTIONS]

Options:
  --sample NAME              Sample label (default: HCC1395)
  --platform NAME            Platform label (default: ONT_5kHz)
  --tumor-bam PATH           Tumor BAM with MM/ML
  --truth-vcf PATH           Override truth VCF
  --truth-bed PATH           Override truth BED (can be empty)
  --truth-total N            Override truth total count for round context
  --truth-set TEXT           Override truth set label for round context
  --caller-existing-vcf PATH Reuse existing ClairS-TO VCF instead of running caller
  --skip-caller              Skip ClairS-TO and require --caller-existing-vcf
  --skip-phase               Skip LongPhase-TO phase
  --skip-haplotag            Skip LongPhase-TO haplotag
  --skip-intersubmod         Skip InterSubMod step
  --skip-filter-analysis     Skip 03_filter_analysis.py
  --archive-tagged-bam       Move tagged BAM and index to Archive_pending_delete after run
  --tag-prefix PATH          Override LongPhase-TO haplotag output prefix
  --threads N                Threads for caller / longphase / intersubmod
  --index-threads N          Threads used for BAM indexing / sort-heavy I/O
  --min-free-gb N            Minimum free space required on output volume before run
  --output-root DIR          Override output root
  --test-focus TEXT          Override round_context test_focus
  --changes TEXT             Override round_context changes
  --conclusion TEXT          Override round_context conclusion
  --next-step TEXT           Override round_context next_step
  --dry-run                  Print commands only
  -h, --help                 Show help
HELP
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --platform) PLATFORM="$2"; USER_SET_PLATFORM=true; shift 2 ;;
        --tumor-bam) TUMOR_BAM="$2"; USER_SET_TUMOR_BAM=true; shift 2 ;;
        --truth-vcf) TRUTH_VCF="$2"; USER_SET_TRUTH_VCF=true; shift 2 ;;
        --truth-bed) TRUTH_BED="$2"; USER_SET_TRUTH_BED=true; shift 2 ;;
        --truth-total) TRUTH_TOTAL="$2"; USER_SET_TRUTH_TOTAL=true; shift 2 ;;
        --truth-set) TRUTH_SET="$2"; USER_SET_TRUTH_SET=true; shift 2 ;;
        --caller-existing-vcf) CALLER_EXISTING_VCF="$2"; USER_SET_CALLER_EXISTING_VCF=true; RUN_CALLER=false; shift 2 ;;
        --skip-caller) RUN_CALLER=false; shift ;;
        --skip-phase) SKIP_PHASE=true; shift ;;
        --skip-haplotag) SKIP_HAPLOTAG=true; shift ;;
        --skip-intersubmod) SKIP_INTERSUBMOD=true; shift ;;
        --skip-filter-analysis) SKIP_FILTER_ANALYSIS=true; shift ;;
        --archive-tagged-bam) ARCHIVE_TAGGED_BAM=true; shift ;;
        --tag-prefix) TAG_PREFIX_OVERRIDE="$2"; shift 2 ;;
        --threads) THREADS_TO_USE="$2"; shift 2 ;;
        --index-threads) INDEX_THREADS_TO_USE="$2"; shift 2 ;;
        --min-free-gb) MIN_FREE_GB="$2"; shift 2 ;;
        --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
        --test-focus) TEST_FOCUS="$2"; shift 2 ;;
        --changes) CHANGES="$2"; shift 2 ;;
        --conclusion) CONCLUSION="$2"; shift 2 ;;
        --next-step) NEXT_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) show_help ;;
        *)
            log_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

PLATFORM_OVERRIDE="${PLATFORM}"
TUMOR_BAM_OVERRIDE="${TUMOR_BAM}"
TRUTH_VCF_OVERRIDE="${TRUTH_VCF}"
TRUTH_BED_OVERRIDE="${TRUTH_BED}"
TRUTH_TOTAL_OVERRIDE="${TRUTH_TOTAL}"
TRUTH_SET_OVERRIDE="${TRUTH_SET}"
CALLER_EXISTING_VCF_OVERRIDE="${CALLER_EXISTING_VCF}"

if sample_cfg="$(get_sample_config "${SAMPLE}" 2>/dev/null)"; then
    eval "${sample_cfg}"
    if [[ "${USER_SET_PLATFORM}" == true ]]; then
        PLATFORM="${PLATFORM_OVERRIDE}"
    elif [[ -n "${PLATFORM_LABEL:-}" ]]; then
        PLATFORM="${PLATFORM_LABEL}"
    fi
    if [[ "${USER_SET_TUMOR_BAM}" == true ]]; then
        TUMOR_BAM="${TUMOR_BAM_OVERRIDE}"
    fi
    if [[ "${USER_SET_TRUTH_VCF}" == true ]]; then
        TRUTH_VCF="${TRUTH_VCF_OVERRIDE}"
    fi
    if [[ "${USER_SET_TRUTH_BED}" == true ]]; then
        TRUTH_BED="${TRUTH_BED_OVERRIDE}"
    fi
    if [[ "${USER_SET_TRUTH_TOTAL}" == true ]]; then
        TRUTH_TOTAL="${TRUTH_TOTAL_OVERRIDE}"
    fi
    if [[ "${USER_SET_TRUTH_SET}" == true ]]; then
        TRUTH_SET="${TRUTH_SET_OVERRIDE}"
    elif [[ -n "${TRUTH_SET_LABEL:-}" ]]; then
        TRUTH_SET="${TRUTH_SET_LABEL}"
    fi
    if [[ "${RUN_CALLER}" == false && "${USER_SET_CALLER_EXISTING_VCF}" == true ]]; then
        CALLER_EXISTING_VCF="${CALLER_EXISTING_VCF_OVERRIDE}"
    elif [[ "${RUN_CALLER}" == false && -n "${TO_SOMATIC_VCF_PILEUP:-}" ]]; then
        CALLER_EXISTING_VCF="${TO_SOMATIC_VCF_PILEUP}"
    fi
fi

if [[ "${RUN_CALLER}" == false && -z "${CALLER_EXISTING_VCF}" ]]; then
    log_error "--skip-caller / --caller-existing-vcf requires an existing caller VCF."
    exit 1
fi

if [[ -z "${INDEX_THREADS_TO_USE}" ]]; then
    INDEX_THREADS_TO_USE="$(recommended_index_threads "${THREADS_TO_USE}")"
fi

validate_file "${TUMOR_BAM}" "Tumor BAM"
validate_file "${REFERENCE_FASTA}" "Reference FASTA"
validate_file "${TRUTH_VCF}" "Truth VCF"
if [[ -n "${TRUTH_BED}" ]]; then
    validate_file "${TRUTH_BED}" "Truth BED"
fi
validate_file "${LONGPHASE_TO_BIN}" "LongPhase-TO binary"

DATE_STR=$(date +%Y%m%d)
ROUND_ROOT_BASE="${OUTPUT_ROOT}/${DATE_STR}_${SAMPLE,,}_to_pilot"
ROUND_ROOT="${ROUND_ROOT_BASE}"
SEQ=1
while [[ -e "${ROUND_ROOT}" ]]; do
    ROUND_ROOT="${ROUND_ROOT_BASE}_${SEQ}"
    SEQ=$((SEQ + 1))
done
mkdir -p "${ROUND_ROOT}"

RUN_LOG="${ROUND_ROOT}/run.log"
exec > >(tee -a "${RUN_LOG}") 2> >(tee -a "${RUN_LOG}" >&2)

STEP01_DIR="${ROUND_ROOT}/step01_clairs_to"
STEP02_DIR="${ROUND_ROOT}/step02_benchmark_clairs_to"
STEP03_DIR="${ROUND_ROOT}/step03_longphase_to"
STEP04_DIR="${ROUND_ROOT}/step04_benchmark_longphase_to"
STEP05_DIR="${ROUND_ROOT}/step05_intersubmod"
mkdir -p "${STEP01_DIR}" "${STEP02_DIR}" "${STEP03_DIR}" "${STEP04_DIR}" "${STEP05_DIR}"

ROUND_CONTEXT="${ROUND_ROOT}/round_context.json"
cat > "${ROUND_CONTEXT}" <<EOF
{
  "created_at": "$(date '+%Y-%m-%d %H:%M:%S %Z')",
  "sample": "${SAMPLE}",
  "sample_set": "tumor_only_pilot",
  "platform": "${PLATFORM}",
  "analysis_mode": "${ANALYSIS_MODE}",
  "truth_set": "${TRUTH_SET}",
  "truth_vcf": "${TRUTH_VCF}",
  "truth_bed": "${TRUTH_BED}",
  "truth_total": ${TRUTH_TOTAL},
  "caller_name": "${CALLER_NAME}",
  "somatic_vcf": "",
  "test_focus": "${TEST_FOCUS}",
  "changes": "${CHANGES}",
  "source_run_dir": "${ROUND_ROOT}",
  "conclusion": "${CONCLUSION}",
  "next_step": "${NEXT_STEP}"
}
EOF

log_info "[TO pilot] Sample=${SAMPLE} Platform=${PLATFORM}"
log_info "[TO pilot] Output root=${ROUND_ROOT}"
log_info "[TO pilot] Threads=${THREADS_TO_USE} IndexThreads=${INDEX_THREADS_TO_USE} MinFreeGB=${MIN_FREE_GB}"
print_disk_space "/big7_disk"
if [[ "${DRY_RUN}" != true ]]; then
    require_min_free_gb "${ROUND_ROOT}" "${MIN_FREE_GB}" "run_longphase_to_intersubmod_pilot:${SAMPLE}"
fi

run_cmd() {
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] $*" >&2
    else
        "$@"
    fi
}

update_context_somatic_vcf() {
    local vcf_path="$1"
    python3 - "$ROUND_CONTEXT" "$vcf_path" <<'PY'
import json
import sys
path, vcf = sys.argv[1], sys.argv[2]
with open(path, "r", encoding="utf-8") as fh:
    payload = json.load(fh)
payload["somatic_vcf"] = vcf
with open(path, "w", encoding="utf-8") as fh:
    json.dump(payload, fh, indent=2, ensure_ascii=False)
    fh.write("\n")
PY
}

# Step 1: ClairS-TO
if [[ "${RUN_CALLER}" == true ]]; then
    log_info "[Step 01] Running ClairS-TO caller"
    CALLER_LOG="${STEP01_DIR}/run_clairs_to.log"
    CALLER_CMD=(
        docker run --rm
        -v /big7_disk:/big7_disk
        -v /big8_disk:/big8_disk
        -v /bip8_disk:/bip8_disk
        -u "$(id -u):$(id -g)"
        hkubal/clairs-to:v0.3.0
        /opt/bin/run_clairs_to
        --tumor_bam_fn "${TUMOR_BAM}"
        --ref_fn "${REFERENCE_FASTA}"
        --threads "${THREADS_TO_USE}"
        --platform "${CALLER_MODEL}"
        --output_dir "${STEP01_DIR}"
    )
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] ${CALLER_CMD[*]}" >&2
    else
        /usr/bin/time -v "${CALLER_CMD[@]}" 2>&1 | tee "${CALLER_LOG}"
    fi
    CALLER_VCF="${STEP01_DIR}/snv.vcf.gz"
else
    log_info "[Step 01] Reusing existing ClairS-TO VCF"
    CALLER_VCF="${CALLER_EXISTING_VCF}"
fi

if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${CALLER_VCF}" "ClairS-TO VCF"
fi
update_context_somatic_vcf "${CALLER_VCF}"

# Step 2: Benchmark caller baseline
log_info "[Step 02] Benchmark ClairS-TO baseline"
run_cmd bash "${PROJECT_ROOT}/scripts/pipeline/utils/benchmark_split_snv_vcf.sh" \
    --input-vcf "${CALLER_VCF}" \
    --truth-vcf "${TRUTH_VCF}" \
    --truth-bed "${TRUTH_BED}" \
    --output-dir "${STEP02_DIR}"

BASELINE_COUNTS="${STEP02_DIR}/variant_counts.txt"
if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${BASELINE_COUNTS}" "Baseline variant counts"
    source "${BASELINE_COUNTS}"
    BASELINE_TP_COUNT="${TP_COUNT}"
    BASELINE_FP_COUNT="${FP_COUNT}"
    BASELINE_FN_COUNT="${FN_COUNT}"
    BASELINE_PASS_TOTAL="${PASS_TOTAL}"
else
    BASELINE_TP_COUNT=0
    BASELINE_FP_COUNT=0
    BASELINE_FN_COUNT=0
    BASELINE_PASS_TOTAL=0
fi

PHASE_PREFIX="${STEP03_DIR}/tumor_phased"
PHASED_VCF="${PHASE_PREFIX}.vcf"
LOH_BED="${PHASE_PREFIX}_LOH.bed"
TAG_PREFIX="${STEP03_DIR}/tumor_tagged"
if [[ -n "${TAG_PREFIX_OVERRIDE}" ]]; then
    TAG_PREFIX="${TAG_PREFIX_OVERRIDE}"
fi
mkdir -p "$(dirname "${TAG_PREFIX}")"
TAGGED_BAM="${TAG_PREFIX}.bam"

# Step 3: LongPhase-TO phase
if [[ "${SKIP_PHASE}" != true ]]; then
    log_info "[Step 03] Running LongPhase-TO phase"
    PHASE_CMD=(
        "${LONGPHASE_TO_BIN}" phase
        -s "${CALLER_VCF}"
        -b "${TUMOR_BAM}"
        -r "${REFERENCE_FASTA}"
        --ont
        --caller "${CALLER_MODE}"
        --pon-file "${PON_FILE}"
        --strict-pon-file "${STRICT_PON_FILE}"
        --loh
        -t "${THREADS_TO_USE}"
        -o "${PHASE_PREFIX}"
    )
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] ${PHASE_CMD[*]}" >&2
    else
        /usr/bin/time -v "${PHASE_CMD[@]}" 2>&1 | tee "${STEP03_DIR}/longphase_to_phase.log"
    fi
fi
if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${PHASED_VCF}" "LongPhase-TO phased VCF"
fi

# Step 4: Benchmark LongPhase-TO phased VCF
log_info "[Step 04] Benchmark LongPhase-TO phased VCF"
run_cmd bash "${PROJECT_ROOT}/scripts/pipeline/utils/benchmark_split_snv_vcf.sh" \
    --input-vcf "${PHASED_VCF}" \
    --truth-vcf "${TRUTH_VCF}" \
    --truth-bed "${TRUTH_BED}" \
    --output-dir "${STEP04_DIR}"

LONGPHASE_COUNTS="${STEP04_DIR}/variant_counts.txt"
if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${LONGPHASE_COUNTS}" "LongPhase-TO variant counts"
    source "${LONGPHASE_COUNTS}"
    LP_PASS_TOTAL="${PASS_TOTAL}"
    LP_TP_COUNT="${TP_COUNT}"
    LP_FP_COUNT="${FP_COUNT}"
    LP_FN_COUNT="${FN_COUNT}"
else
    LP_PASS_TOTAL=0
    LP_TP_COUNT=0
    LP_FP_COUNT=0
    LP_FN_COUNT=0
fi

# Step 5: LongPhase-TO haplotag + QC
if [[ "${SKIP_HAPLOTAG}" != true ]]; then
    log_info "[Step 05] Running LongPhase-TO haplotag"
    HAPLOTAG_CMD=(
        "${LONGPHASE_TO_BIN}" haplotag
        -s "${PHASED_VCF}"
        -b "${TUMOR_BAM}"
        -r "${REFERENCE_FASTA}"
        --tagSupplementary
        -t "${THREADS_TO_USE}"
        -o "${TAG_PREFIX}"
    )
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY-RUN] ${HAPLOTAG_CMD[*]}" >&2
    else
        /usr/bin/time -v "${HAPLOTAG_CMD[@]}" 2>&1 | tee "${STEP03_DIR}/longphase_to_haplotag.log"
    fi
fi

if [[ "${DRY_RUN}" != true ]]; then
    validate_file "${TAGGED_BAM}" "LongPhase-TO tagged BAM"
    run_cmd "${SAMTOOLS}" index -@ "${INDEX_THREADS_TO_USE}" "${TAGGED_BAM}"
    run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/haplotag_qc.py" \
        --bam "${TAGGED_BAM}" \
        --sample "${SAMPLE}" \
        --output-tsv "${ROUND_ROOT}/haplotag_qc.tsv" \
        --samtools "${SAMTOOLS}"
fi

# Step 6: InterSubMod
if [[ "${SKIP_INTERSUBMOD}" != true ]]; then
    log_info "[Step 06] Running InterSubMod"
    run_cmd "${PROJECT_ROOT}/scripts/pipeline/steps/02_intersubmod.sh" \
        --tagged-bam "${TAGGED_BAM}" \
        --tp-vcf "${STEP04_DIR}/filtered_snv_tp.vcf.gz" \
        --fp-vcf "${STEP04_DIR}/filtered_snv_fp.vcf.gz" \
        --output-dir "${STEP05_DIR}" \
        --threads "${THREADS_TO_USE}"
fi

# Step 7: Filter analysis (InterSubMod F1)
if [[ "${SKIP_FILTER_ANALYSIS}" != true ]]; then
    log_info "[Step 07] Running methylation filter analysis"
    run_cmd python3 "${PROJECT_ROOT}/scripts/pipeline/steps/03_filter_analysis.py" \
        --tp-summary "${STEP05_DIR}/intersubmod_tp/significance_summary.csv" \
        --fp-summary "${STEP05_DIR}/intersubmod_fp/significance_summary.csv" \
        --tp-count "${LP_TP_COUNT}" \
        --fp-count "${LP_FP_COUNT}" \
        --truth-total "${TRUTH_TOTAL}" \
        --tp-vcf-file "${STEP04_DIR}/filtered_snv_tp.vcf.gz" \
        --fp-vcf-file "${STEP04_DIR}/filtered_snv_fp.vcf.gz" \
        --sample "${SAMPLE}" \
        --mode "${ANALYSIS_MODE}" \
        --purity "100" \
        --output-dir "${ROUND_ROOT}"
fi

# Step 8: Design validation and benchmark summary
if [[ -f "${STEP05_DIR}/intersubmod_tp/significance_summary.csv" && -f "${STEP05_DIR}/intersubmod_fp/significance_summary.csv" ]]; then
    run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/validate_method_design.py" \
        --summary-csv "${STEP05_DIR}/intersubmod_tp/significance_summary.csv" \
        --summary-csv "${STEP05_DIR}/intersubmod_fp/significance_summary.csv" \
        --region-root "${STEP05_DIR}/intersubmod_tp" \
        --region-root "${STEP05_DIR}/intersubmod_fp" \
        --haplotag-qc "${ROUND_ROOT}/haplotag_qc.tsv" \
        --sample "${SAMPLE}" \
        --output-dir "${ROUND_ROOT}"
fi

if [[ -f "${ROUND_ROOT}/metrics.json" ]]; then
    run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/compare_benchmark_f1.py" \
        --run-dir "${ROUND_ROOT}" \
        --context-json "${ROUND_CONTEXT}" \
        --baseline-counts-file "${STEP04_DIR}/variant_counts.txt" \
        --baseline-method "LongPhase-TO" \
        --output-tsv "${ROUND_ROOT}/benchmark_comparison.tsv" \
        --output-md "${ROUND_ROOT}/benchmark_comparison.md"
fi

if [[ -f "${ROUND_ROOT}/metrics.json" ]]; then
    run_cmd python3 "${PROJECT_ROOT}/scripts/analysis/build_round_dashboard.py" \
        --sample-dir "${ROUND_ROOT}"
fi

# Step 9: Append explicit TO benchmark comparison note
python3 - "${ROUND_ROOT}" "${BASELINE_PASS_TOTAL}" "${BASELINE_TP_COUNT}" "${BASELINE_FP_COUNT}" "${BASELINE_FN_COUNT}" "${LP_PASS_TOTAL}" "${LP_TP_COUNT}" "${LP_FP_COUNT}" "${LP_FN_COUNT}" <<'PY'
import json
import math
import sys
from pathlib import Path

round_root = Path(sys.argv[1])
baseline_calls, baseline_tp, baseline_fp, baseline_fn = map(int, sys.argv[2:6])
lp_calls, lp_tp, lp_fp, lp_fn = map(int, sys.argv[6:10])
metrics_path = round_root / "metrics.json"
summary_path = round_root / "to_stage_metrics.json"

def calc(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "calls_total": tp + fp if tp + fp else 0,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
    }

payload = {
    "baseline_clairs_to": calc(baseline_tp, baseline_fp, baseline_fn),
    "longphase_to": calc(lp_tp, lp_fp, lp_fn),
}
if metrics_path.exists():
    with metrics_path.open("r", encoding="utf-8") as fh:
        payload["intersubmod"] = json.load(fh).get("filtered", {})
with summary_path.open("w", encoding="utf-8") as fh:
    json.dump(payload, fh, indent=2, ensure_ascii=False)
    fh.write("\n")
PY

# Step 10: Optional archive move for tagged BAM
if [[ "${ARCHIVE_TAGGED_BAM}" == true && "${DRY_RUN}" != true ]]; then
    mkdir -p "${ARCHIVE_PENDING_DELETE}"
    ARCHIVE_TARGET_PREFIX="${ARCHIVE_PENDING_DELETE}/$(basename "${TAGGED_BAM}")"
    if [[ -e "${ARCHIVE_TARGET_PREFIX}" ]]; then
        ARCHIVE_TARGET_PREFIX="${ARCHIVE_PENDING_DELETE}/$(basename "${TAG_PREFIX}")_$(date +%Y%m%d_%H%M%S).bam"
    fi
    log_info "[Step 10] Moving tagged BAM to archive pending delete"
    mv "${TAGGED_BAM}" "${ARCHIVE_TARGET_PREFIX}"
    if [[ -f "${TAGGED_BAM}.bai" ]]; then
        mv "${TAGGED_BAM}.bai" "${ARCHIVE_TARGET_PREFIX}.bai"
    fi
    echo -e "original_path\tarchived_path" > "${ROUND_ROOT}/archive_pending_delete.tsv"
    echo -e "${TAGGED_BAM}\t${ARCHIVE_TARGET_PREFIX}" >> "${ROUND_ROOT}/archive_pending_delete.tsv"
fi

log_info "[TO pilot] Completed: ${ROUND_ROOT}"
