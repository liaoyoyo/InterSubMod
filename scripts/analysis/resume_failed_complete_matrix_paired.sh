#!/bin/bash
# Scan complete_matrix paired runs that failed after LongPhase-S and resume them in place.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

CANONICAL_ROOT="${CANONICAL_OUTPUT_ROOT}"
RUN_TAG_FILTER="complete_matrix"
MAX_PARALLEL="${MAX_PARALLEL_DEFAULT}"
THREADS_PER_JOB=""
MIN_FREE_GB="${MIN_FREE_GB_DEFAULT}"
BATCH_ROOT=""
DRY_RUN=false

show_help() {
    cat <<'EOF'
Usage: resume_failed_complete_matrix_paired.sh [OPTIONS]

Options:
  --canonical-root DIR   Canonical root to scan
  --run-tag TEXT         Only scan run directories containing this text (default: complete_matrix)
  --max-parallel N       Concurrent resume jobs (default: config MAX_PARALLEL_DEFAULT)
  --threads-per-job N    InterSubMod threads per resume job
  --min-free-gb N        Minimum free disk on /big7_disk before each job
  --batch-root DIR       Override resume batch metadata root
  --dry-run              Print queue only
  -h, --help             Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --canonical-root) CANONICAL_ROOT="$2"; shift 2 ;;
        --run-tag) RUN_TAG_FILTER="$2"; shift 2 ;;
        --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
        --threads-per-job) THREADS_PER_JOB="$2"; shift 2 ;;
        --min-free-gb) MIN_FREE_GB="$2"; shift 2 ;;
        --batch-root) BATCH_ROOT="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) show_help; exit 0 ;;
        *)
            log_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "${THREADS_PER_JOB}" ]]; then
    THREADS_PER_JOB=$(( THREADS / MAX_PARALLEL ))
    if (( THREADS_PER_JOB < 1 )); then
        THREADS_PER_JOB=1
    fi
fi

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -z "${BATCH_ROOT}" ]]; then
    BATCH_ROOT="${SYNTHESIS_OUTPUT_ROOT}/batch_runs/${TIMESTAMP}_resume_failed_paired"
fi
mkdir -p "${BATCH_ROOT}/logs"

BATCH_LOG="${BATCH_ROOT}/batch.log"
QUEUE_TSV="${BATCH_ROOT}/job_queue.tsv"
STATUS_TSV="${BATCH_ROOT}/job_status.tsv"

declare -a JOB_IDS=()
declare -a JOB_SAMPLES=()
declare -a JOB_MODES=()
declare -a JOB_RUN_DIRS=()
declare -a JOB_COMMANDS=()
declare -a JOB_LOGS=()
declare -a JOB_STATUSES=()
declare -a JOB_STARTS=()
declare -a JOB_ENDS=()
declare -a JOB_EXITS=()
declare -A PID_TO_INDEX=()

log_batch() {
    local msg="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${msg}" | tee -a "${BATCH_LOG}" >&2
}

append_job() {
    local sample="$1"
    local mode="$2"
    local run_dir="$3"
    local job_id="resume_${sample}_${mode}"
    local log_path="${BATCH_ROOT}/logs/${job_id}.log"
    local cmd="cd ${PROJECT_ROOT} && ./scripts/analysis/resume_paired_run_from_tagged_outputs.sh --sample ${sample} --canonical-mode ${mode} --run-dir ${run_dir} --threads ${THREADS_PER_JOB} --min-free-gb ${MIN_FREE_GB}"

    JOB_IDS+=("${job_id}")
    JOB_SAMPLES+=("${sample}")
    JOB_MODES+=("${mode}")
    JOB_RUN_DIRS+=("${run_dir}")
    JOB_COMMANDS+=("${cmd}")
    JOB_LOGS+=("${log_path}")
    JOB_STATUSES+=("queued")
    JOB_STARTS+=("")
    JOB_ENDS+=("")
    JOB_EXITS+=("")
}

write_queue_manifest() {
    printf "job_id\tsample\tcanonical_mode\trun_dir\tthreads\tlog_path\tcommand\n" > "${QUEUE_TSV}"
    local idx
    for idx in "${!JOB_IDS[@]}"; do
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${JOB_IDS[$idx]}" \
            "${JOB_SAMPLES[$idx]}" \
            "${JOB_MODES[$idx]}" \
            "${JOB_RUN_DIRS[$idx]}" \
            "${THREADS_PER_JOB}" \
            "${JOB_LOGS[$idx]}" \
            "${JOB_COMMANDS[$idx]}" >> "${QUEUE_TSV}"
    done
}

write_status_manifest() {
    printf "job_id\tsample\tcanonical_mode\tstatus\texit_code\tstarted_at\tended_at\trun_dir\tlog_path\n" > "${STATUS_TSV}"
    local idx
    for idx in "${!JOB_IDS[@]}"; do
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${JOB_IDS[$idx]}" \
            "${JOB_SAMPLES[$idx]}" \
            "${JOB_MODES[$idx]}" \
            "${JOB_STATUSES[$idx]}" \
            "${JOB_EXITS[$idx]}" \
            "${JOB_STARTS[$idx]}" \
            "${JOB_ENDS[$idx]}" \
            "${JOB_RUN_DIRS[$idx]}" \
            "${JOB_LOGS[$idx]}" >> "${STATUS_TSV}"
    done
}

launch_job() {
    local idx="$1"
    JOB_STATUSES[$idx]="running"
    JOB_STARTS[$idx]="$(date '+%Y-%m-%d %H:%M:%S')"
    write_status_manifest

    log_batch "Launch ${JOB_IDS[$idx]}: ${JOB_COMMANDS[$idx]}"
    (
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] START ${JOB_IDS[$idx]}"
        echo "[COMMAND] ${JOB_COMMANDS[$idx]}"
        bash -lc "${JOB_COMMANDS[$idx]}"
    ) > "${JOB_LOGS[$idx]}" 2>&1 &

    PID_TO_INDEX[$!]="${idx}"
}

reap_one_job() {
    while true; do
        local pid
        for pid in "${!PID_TO_INDEX[@]}"; do
            if ! kill -0 "${pid}" 2>/dev/null; then
                local idx="${PID_TO_INDEX[$pid]}"
                local exit_code=0
                if wait "${pid}"; then
                    exit_code=0
                    JOB_STATUSES[$idx]="completed"
                else
                    exit_code=$?
                    JOB_STATUSES[$idx]="failed"
                fi
                JOB_EXITS[$idx]="${exit_code}"
                JOB_ENDS[$idx]="$(date '+%Y-%m-%d %H:%M:%S')"
                unset PID_TO_INDEX["${pid}"]
                write_status_manifest
                log_batch "Finish ${JOB_IDS[$idx]}: status=${JOB_STATUSES[$idx]} exit=${exit_code}"
                return 0
            fi
        done
        sleep 10
    done
}

while IFS=$'\t' read -r sample mode run_dir; do
    [[ -z "${sample}" ]] && continue
    append_job "${sample}" "${mode}" "${run_dir}"
done < <(
    python3 - "${CANONICAL_ROOT}" "${RUN_TAG_FILTER}" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
run_tag = sys.argv[2]
for run_dir in sorted(root.glob(f"*/paired_*/*{run_tag}*")):
    if not run_dir.is_dir():
        continue
    sample = run_dir.parts[-3]
    mode = run_dir.parts[-2]
    lp = run_dir / "longphase_s"
    tagged = lp / f"{sample}_tagged.bam"
    sc_vcf = lp / f"{sample}_tagged_sc.vcf"
    tp_vcf = lp / "filtered_snv_tp.vcf.gz"
    fp_vcf = lp / "filtered_snv_fp.vcf.gz"
    metrics = run_dir / "metrics.json"
    benchmark_tsv = run_dir / "benchmark_comparison.tsv"
    round_summary = run_dir / "round_summary.md"
    method_validation = run_dir / "method_design_validation.tsv"
    agreement = run_dir / "label_cluster_agreement.tsv"
    if not tagged.exists() or not sc_vcf.exists():
        continue
    if (
        tp_vcf.exists()
        and fp_vcf.exists()
        and metrics.exists()
        and benchmark_tsv.exists()
        and round_summary.exists()
        and method_validation.exists()
        and agreement.exists()
    ):
        continue
    print(f"{sample}\t{mode}\t{run_dir}")
PY
)

write_queue_manifest
write_status_manifest

log_batch "Resume batch root: ${BATCH_ROOT}"
log_batch "Queue manifest: ${QUEUE_TSV}"
log_batch "Status manifest: ${STATUS_TSV}"
log_batch "Threads per job=${THREADS_PER_JOB}, max_parallel=${MAX_PARALLEL}, min_free_gb=${MIN_FREE_GB}"

if [[ ${#JOB_IDS[@]} -eq 0 ]]; then
    log_batch "No resumable paired complete_matrix runs found."
    exit 0
fi

if [[ "${DRY_RUN}" == true ]]; then
    local_idx=0
    for local_idx in "${!JOB_IDS[@]}"; do
        log_batch "DRY-RUN ${JOB_IDS[$local_idx]} => ${JOB_COMMANDS[$local_idx]}"
    done
    exit 0
fi

require_min_free_gb "/big7_disk" "${MIN_FREE_GB}" "resume_failed_paired_batch"

for idx in "${!JOB_IDS[@]}"; do
    while (( ${#PID_TO_INDEX[@]} >= MAX_PARALLEL )); do
        reap_one_job
    done
    launch_job "${idx}"
done

while (( ${#PID_TO_INDEX[@]} > 0 )); do
    reap_one_job
done

log_batch "Resume batch complete."
