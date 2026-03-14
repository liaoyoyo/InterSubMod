#!/bin/bash
# run_pure_research_round.sh - Orchestrate pure-sample research round bundles.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

SAMPLE_SET="main_5khz"
MODE="s-pure"
LOCAL_THREADS="${THREADS}"
VCF_SOURCE="output"
RUN_PIPELINE=false
DRY_RUN=false
SKIP_CALLER_INPUT=false
DISTANCE_METHOD="BERNOULLI"
WINDOW_BP="5000"
READ_FILTER_PROFILE="stable_default"
METHYLATION_MODE="analog_full"
CLUSTER_METHOD="UPGMA"
LABEL_MODE="hp+allele"
MIN_READS="20"
TEST_FOCUS=""
CHANGES=""
CONCLUSION=""
NEXT_STEP=""
OUTPUT_DIR=""
CUSTOM_SAMPLE=""
declare -a RUN_DIRS=()

show_help() {
    cat <<'EOF'
Usage: run_pure_research_round.sh [OPTIONS]

Options:
  --sample-set NAME        main_5khz | dorado_validation | other_pure | all_pure | custom
  --sample NAME            Required when --sample-set custom
  --mode MODE              Pipeline mode (default: s-pure)
  --run-pipeline           Execute pipeline/run_benchmark.sh before analysis
  --run-dir DIR            Existing source run directory (repeatable)
  --skip-caller-input      Skip caller-input benchmark row for faster round aggregation
  --output-dir DIR         Research round output dir
  --threads N              Threads for pipeline execution
  --vcf-source SRC         output | pileup | indel
  --distance-method NAME   Context value for dashboard/sensitivity output
  --window-bp N            Context value for dashboard/sensitivity output
  --read-filter-profile X  Context value for dashboard/sensitivity output
  --methylation-mode X     Context value for dashboard/sensitivity output
  --cluster-method X       Context value for dashboard/sensitivity output
  --label-mode X           Context value for dashboard/sensitivity output
  --min-reads N            Context value for dashboard/sensitivity output
  --test-focus TEXT        Short experiment focus
  --changes TEXT           Short change summary
  --conclusion TEXT        Initial conclusion text
  --next-step TEXT         Suggested next step
  --dry-run                Print commands without executing
  -h, --help               Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-set) SAMPLE_SET="$2"; shift 2 ;;
        --sample) CUSTOM_SAMPLE="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        --run-pipeline) RUN_PIPELINE=true; shift ;;
        --run-dir) RUN_DIRS+=("$2"); shift 2 ;;
        --skip-caller-input) SKIP_CALLER_INPUT=true; shift ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) LOCAL_THREADS="$2"; shift 2 ;;
        --vcf-source) VCF_SOURCE="$2"; shift 2 ;;
        --distance-method) DISTANCE_METHOD="$2"; shift 2 ;;
        --window-bp) WINDOW_BP="$2"; shift 2 ;;
        --read-filter-profile) READ_FILTER_PROFILE="$2"; shift 2 ;;
        --methylation-mode) METHYLATION_MODE="$2"; shift 2 ;;
        --cluster-method) CLUSTER_METHOD="$2"; shift 2 ;;
        --label-mode) LABEL_MODE="$2"; shift 2 ;;
        --min-reads) MIN_READS="$2"; shift 2 ;;
        --test-focus) TEST_FOCUS="$2"; shift 2 ;;
        --changes) CHANGES="$2"; shift 2 ;;
        --conclusion) CONCLUSION="$2"; shift 2 ;;
        --next-step) NEXT_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

resolve_samples() {
    case "${SAMPLE_SET}" in
        main_5khz) echo "HCC1395" ;;
        dorado_validation) echo "HCC1395_DORADO" ;;
        other_pure) echo "COLO829 H1437 H2009 HCC1937 HCC1954" ;;
        all_pure) echo "HCC1395 HCC1395_DORADO COLO829 H1437 H2009 HCC1937 HCC1954" ;;
        custom)
            if [[ -z "${CUSTOM_SAMPLE}" ]]; then
                echo "[ERROR] --sample is required when --sample-set custom" >&2
                exit 1
            fi
            echo "${CUSTOM_SAMPLE}"
            ;;
        *)
            echo "[ERROR] Unknown sample set: ${SAMPLE_SET}" >&2
            exit 1
            ;;
    esac
}

infer_platform() {
    local sample="$1"
    if [[ "${sample}" == "HCC1395_DORADO" ]]; then
        echo "ONT_Dorado"
    elif [[ "${sample}" == "HCC1395" ]]; then
        echo "ONT_5kHz"
    else
        echo "ONT"
    fi
}

latest_run_dir() {
    local sample="$1"
    local sample_root="${OUTPUT_ROOT}/${MODE}/${sample}"
    if [[ ! -d "${sample_root}" ]]; then
        return 0
    fi
    find "${sample_root}" -mindepth 1 -maxdepth 1 -type d | sort | tail -n 1
}

copy_if_exists() {
    local src="$1"
    local dst="$2"
    if [[ -f "${src}" ]]; then
        cp "${src}" "${dst}"
    fi
}

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -z "${OUTPUT_DIR}" ]]; then
    OUTPUT_DIR="${OUTPUT_ROOT}/research_rounds/${TIMESTAMP}_${SAMPLE_SET}"
fi
mkdir -p "${OUTPUT_DIR}"

SAMPLES=($(resolve_samples))

if [[ "${RUN_PIPELINE}" == true ]]; then
    for sample in "${SAMPLES[@]}"; do
        cmd=(
            "${PROJECT_ROOT}/scripts/pipeline/run_benchmark.sh"
            --mode "${MODE}"
            --sample "${sample}"
            --threads "${LOCAL_THREADS}"
            --vcf-source "${VCF_SOURCE}"
        )
        if [[ "${DRY_RUN}" == true ]]; then
            cmd+=(--dry-run)
            echo "[DRY-RUN] ${cmd[*]}" >&2
        else
            "${cmd[@]}"
            discovered="$(latest_run_dir "${sample}")"
            if [[ -n "${discovered}" ]]; then
                RUN_DIRS+=("${discovered}")
            fi
        fi
    done
fi

if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
    for sample in "${SAMPLES[@]}"; do
        discovered="$(latest_run_dir "${sample}")"
        if [[ -n "${discovered}" ]]; then
            RUN_DIRS+=("${discovered}")
        fi
    done
fi

if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
    echo "[ERROR] No source run directories available." >&2
    exit 1
fi

MANIFEST="${OUTPUT_DIR}/round_manifest.tsv"
printf "sample\tplatform\tsource_run_dir\tsample_bundle_dir\n" > "${MANIFEST}"

declare -a SAMPLE_BUNDLES=()

for run_dir in "${RUN_DIRS[@]}"; do
    if [[ ! -d "${run_dir}" ]]; then
        echo "[WARN] Skip missing run dir: ${run_dir}" >&2
        continue
    fi

    sample="$(basename "$(dirname "${run_dir}")")"
    platform="$(infer_platform "${sample}")"
    sample_bundle="${OUTPUT_DIR}/${sample}"
    mkdir -p "${sample_bundle}"
    SAMPLE_BUNDLES+=("${sample_bundle}")

    eval "$(get_sample_config "${sample}")"

    python3 - <<PY
import json
from pathlib import Path

payload = {
    "created_at": "${TIMESTAMP}",
    "sample_set": "${SAMPLE_SET}",
    "sample": "${sample}",
    "platform": "${platform}",
    "analysis_mode": "${MODE}",
    "region_scope": "${MODE}",
    "metric_family": "paired_pure_research",
    "distance_metric": "${DISTANCE_METHOD}",
    "window_bp": "${WINDOW_BP}",
    "read_filter_profile": "${READ_FILTER_PROFILE}",
    "methylation_mode": "${METHYLATION_MODE}",
    "cluster_method": "${CLUSTER_METHOD}",
    "label_mode": "${LABEL_MODE}",
    "min_reads": "${MIN_READS}",
    "test_focus": "${TEST_FOCUS}",
    "changes": "${CHANGES}",
    "conclusion": "${CONCLUSION}",
    "next_step": "${NEXT_STEP}",
    "source_run_dir": "${run_dir}",
    "somatic_vcf": "${SOMATIC_VCF}",
    "truth_vcf": "${TRUTH_VCF}",
    "truth_bed": "${TRUTH_BED}",
    "caller_name": "",
}
Path("${sample_bundle}/round_context.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\\n", encoding="utf-8")
PY

    copy_if_exists "${run_dir}/metrics.json" "${sample_bundle}/metrics.json"
    copy_if_exists "${run_dir}/longphase_s/haplotag_qc.tsv" "${sample_bundle}/haplotag_qc.tsv"

    COMPARE_CMD=(
        python3 "${SCRIPT_DIR}/compare_benchmark_f1.py"
        --run-dir "${run_dir}"
        --context-json "${sample_bundle}/round_context.json"
        --output-tsv "${sample_bundle}/benchmark_comparison.tsv"
        --output-md "${sample_bundle}/benchmark_comparison.md"
    )
    if [[ "${SKIP_CALLER_INPUT}" == true ]]; then
        COMPARE_CMD+=(--no-caller-input)
    fi
    "${COMPARE_CMD[@]}"

    python3 "${SCRIPT_DIR}/validate_method_design.py" \
        --summary-csv "${run_dir}/intersubmod_tp/significance_summary.csv" \
        --summary-csv "${run_dir}/intersubmod_fp/significance_summary.csv" \
        --region-root "${run_dir}/intersubmod_tp" \
        --region-root "${run_dir}/intersubmod_fp" \
        --haplotag-qc "${run_dir}/longphase_s/haplotag_qc.tsv" \
        --sample "${sample}" \
        --min-reads "${MIN_READS}" \
        --output-dir "${sample_bundle}"

    python3 "${SCRIPT_DIR}/build_round_dashboard.py" \
        --sample-dir "${sample_bundle}"

    printf "%s\t%s\t%s\t%s\n" "${sample}" "${platform}" "${run_dir}" "${sample_bundle}" >> "${MANIFEST}"
done

if [[ ${#SAMPLE_BUNDLES[@]} -gt 0 ]]; then
    MATRIX_CMD=(python3 "${SCRIPT_DIR}/methodology_sensitivity_matrix.py")
    for bundle in "${SAMPLE_BUNDLES[@]}"; do
        MATRIX_CMD+=(--run-dir "${bundle}")
    done
    MATRIX_CMD+=(--output-tsv "${OUTPUT_DIR}/methodology_sensitivity.tsv")
    "${MATRIX_CMD[@]}"
fi

python3 - "${OUTPUT_DIR}" <<'PY'
from pathlib import Path
import csv
import sys

root = Path(sys.argv[1])
manifest = root / "round_manifest.tsv"
lines = ["# Research Round", "", f"- Round dir: `{root}`", ""]
if manifest.exists():
    lines.append("## Samples")
    lines.append("")
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            bundle = Path(row["sample_bundle_dir"])
            rel = bundle.relative_to(root)
            lines.append(f"- [{sample}]({rel}/round_summary.md)")
    lines.append("")
if (root / "methodology_sensitivity.tsv").exists():
    lines.append("- [methodology_sensitivity.tsv](methodology_sensitivity.tsv)")
root.joinpath("README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
PY

echo "[run_pure_research_round] Output dir: ${OUTPUT_DIR}" >&2
