#!/bin/bash
# Run the full runnable big7 experiment matrix and keep tagged BAM outputs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

MAX_PARALLEL="${MAX_PARALLEL_DEFAULT}"
THREADS_PER_JOB=""
INDEX_THREADS_PER_JOB=""
MIN_FREE_GB="${MIN_FREE_GB_DEFAULT}"
RUN_TAG="complete_matrix"
RUN_POSTPROCESS=true
DRY_RUN=false
BATCH_ROOT=""
INCLUDE_SAMPLES_CSV=""
EXCLUDE_SAMPLES_CSV=""
EXCLUDE_JOB_IDS_CSV=""

PAIRED_SAMPLES=(HCC1395 HCC1395_DORADO COLO829 H1437 H2009 HCC1937 HCC1954)
TO_SAMPLES=(HCC1395 HCC1395_DORADO COLO829 H1437 H2009 HCC1937 HCC1954)

show_help() {
    cat <<'EOF'
Usage: run_complete_big7_experiments.sh [OPTIONS]

Options:
  --max-parallel N      Number of concurrent jobs (default: config MAX_PARALLEL_DEFAULT)
  --threads-per-job N   Threads per job; default = nproc / max-parallel
  --index-threads N     Threads used for BAM indexing / sort-heavy I/O
  --min-free-gb N       Minimum free disk on /big7_disk required per job
  --run-tag TEXT        Tag suffix for paired canonical runs
  --batch-root DIR      Override batch metadata root
  --include-samples CSV Run only these samples (comma-separated)
  --exclude-samples CSV Skip these samples (comma-separated)
  --exclude-job-ids CSV Skip these specific job ids (comma-separated)
  --no-postprocess      Skip paired round aggregation after jobs finish
  --dry-run             Print queue and commands only
  -h, --help            Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
        --threads-per-job) THREADS_PER_JOB="$2"; shift 2 ;;
        --index-threads) INDEX_THREADS_PER_JOB="$2"; shift 2 ;;
        --min-free-gb) MIN_FREE_GB="$2"; shift 2 ;;
        --run-tag) RUN_TAG="$2"; shift 2 ;;
        --batch-root) BATCH_ROOT="$2"; shift 2 ;;
        --include-samples) INCLUDE_SAMPLES_CSV="$2"; shift 2 ;;
        --exclude-samples) EXCLUDE_SAMPLES_CSV="$2"; shift 2 ;;
        --exclude-job-ids) EXCLUDE_JOB_IDS_CSV="$2"; shift 2 ;;
        --no-postprocess) RUN_POSTPROCESS=false; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) show_help; exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
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

if [[ -z "${INDEX_THREADS_PER_JOB}" ]]; then
    INDEX_THREADS_PER_JOB="$(recommended_index_threads "${THREADS_PER_JOB}")"
fi

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -z "${BATCH_ROOT}" ]]; then
    BATCH_ROOT="${SYNTHESIS_OUTPUT_ROOT}/batch_runs/${TIMESTAMP}_complete_matrix"
fi

mkdir -p "${BATCH_ROOT}/logs"

BATCH_LOG="${BATCH_ROOT}/batch.log"
QUEUE_TSV="${BATCH_ROOT}/job_queue.tsv"
STATUS_TSV="${BATCH_ROOT}/job_status.tsv"
SUMMARY_MD="${BATCH_ROOT}/batch_summary.md"
POSTPROCESS_TSV="${BATCH_ROOT}/postprocess.tsv"

declare -a JOB_IDS=()
declare -a JOB_TYPES=()
declare -a JOB_SAMPLES=()
declare -a JOB_MODES=()
declare -a JOB_COMMANDS=()
declare -a JOB_LOGS=()
declare -a JOB_STATUSES=()
declare -a JOB_STARTS=()
declare -a JOB_ENDS=()
declare -a JOB_EXITS=()

declare -a POSTPROCESS_NAMES=()
declare -a POSTPROCESS_COMMANDS=()
declare -a POSTPROCESS_LOGS=()
declare -a POSTPROCESS_OUTPUTS=()
declare -a POSTPROCESS_STATUSES=()

declare -A PID_TO_INDEX=()

csv_has_token() {
    local csv="${1:-}"
    local value="${2:-}"
    local token

    [[ -z "${csv}" ]] && return 1
    IFS=',' read -r -a tokens <<< "${csv}"
    for token in "${tokens[@]}"; do
        if [[ "${token}" == "${value}" ]]; then
            return 0
        fi
    done
    return 1
}

log_batch() {
    local msg="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${msg}" | tee -a "${BATCH_LOG}" >&2
}

append_job() {
    local job_id="$1"
    local job_type="$2"
    local sample="$3"
    local mode="$4"
    local cmd="$5"
    local log_path="$6"

    if [[ -n "${INCLUDE_SAMPLES_CSV}" ]] && ! csv_has_token "${INCLUDE_SAMPLES_CSV}" "${sample}"; then
        return 0
    fi
    if csv_has_token "${EXCLUDE_SAMPLES_CSV}" "${sample}"; then
        return 0
    fi
    if csv_has_token "${EXCLUDE_JOB_IDS_CSV}" "${job_id}"; then
        return 0
    fi

    JOB_IDS+=("${job_id}")
    JOB_TYPES+=("${job_type}")
    JOB_SAMPLES+=("${sample}")
    JOB_MODES+=("${mode}")
    JOB_COMMANDS+=("${cmd}")
    JOB_LOGS+=("${log_path}")
    JOB_STATUSES+=("queued")
    JOB_STARTS+=("")
    JOB_ENDS+=("")
    JOB_EXITS+=("")
}

write_queue_manifest() {
    printf "job_id\tjob_type\tsample\tmode\tthreads_per_job\tmin_free_gb\tlog_path\tcommand\n" > "${QUEUE_TSV}"
    local idx
    for idx in "${!JOB_IDS[@]}"; do
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${JOB_IDS[$idx]}" \
            "${JOB_TYPES[$idx]}" \
            "${JOB_SAMPLES[$idx]}" \
            "${JOB_MODES[$idx]}" \
            "${THREADS_PER_JOB}" \
            "${MIN_FREE_GB}" \
            "${JOB_LOGS[$idx]}" \
            "${JOB_COMMANDS[$idx]}" >> "${QUEUE_TSV}"
    done
}

write_status_manifest() {
    printf "job_id\tjob_type\tsample\tmode\tstatus\texit_code\tstarted_at\tended_at\tlog_path\n" > "${STATUS_TSV}"
    local idx
    for idx in "${!JOB_IDS[@]}"; do
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${JOB_IDS[$idx]}" \
            "${JOB_TYPES[$idx]}" \
            "${JOB_SAMPLES[$idx]}" \
            "${JOB_MODES[$idx]}" \
            "${JOB_STATUSES[$idx]}" \
            "${JOB_EXITS[$idx]}" \
            "${JOB_STARTS[$idx]}" \
            "${JOB_ENDS[$idx]}" \
            "${JOB_LOGS[$idx]}" >> "${STATUS_TSV}"
    done
}

launch_job() {
    local idx="$1"
    local job_id="${JOB_IDS[$idx]}"
    local cmd="${JOB_COMMANDS[$idx]}"
    local log_path="${JOB_LOGS[$idx]}"

    JOB_STATUSES[$idx]="running"
    JOB_STARTS[$idx]="$(date '+%Y-%m-%d %H:%M:%S')"
    write_status_manifest

    log_batch "Launch ${job_id}: ${cmd}"
    (
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] START ${job_id}"
        echo "[COMMAND] ${cmd}"
        bash -lc "${cmd}"
    ) > "${log_path}" 2>&1 &

    PID_TO_INDEX[$!]="${idx}"
}

reap_one_job() {
    while true; do
        local pid
        for pid in "${!PID_TO_INDEX[@]}"; do
            if ! kill -0 "${pid}" 2>/dev/null; then
                local idx="${PID_TO_INDEX[$pid]}"
                local job_id="${JOB_IDS[$idx]}"
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
                log_batch "Finish ${job_id}: status=${JOB_STATUSES[$idx]} exit=${exit_code}"
                return 0
            fi
        done
        sleep 10
    done
}

run_postprocess_step() {
    local name="$1"
    local output_dir="$2"
    local cmd="$3"
    local log_path="$4"
    local status="completed"

    log_batch "Postprocess ${name}: ${cmd}"
    if [[ "${DRY_RUN}" == true ]]; then
        status="dry_run"
    else
        if ! bash -lc "${cmd}" > "${log_path}" 2>&1; then
            status="failed"
        fi
    fi

    POSTPROCESS_NAMES+=("${name}")
    POSTPROCESS_COMMANDS+=("${cmd}")
    POSTPROCESS_LOGS+=("${log_path}")
    POSTPROCESS_OUTPUTS+=("${output_dir}")
    POSTPROCESS_STATUSES+=("${status}")
}

write_postprocess_manifest() {
    printf "name\tstatus\toutput_dir\tlog_path\tcommand\n" > "${POSTPROCESS_TSV}"
    local idx
    for idx in "${!POSTPROCESS_NAMES[@]}"; do
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "${POSTPROCESS_NAMES[$idx]}" \
            "${POSTPROCESS_STATUSES[$idx]}" \
            "${POSTPROCESS_OUTPUTS[$idx]}" \
            "${POSTPROCESS_LOGS[$idx]}" \
            "${POSTPROCESS_COMMANDS[$idx]}" >> "${POSTPROCESS_TSV}"
    done
}

write_summary() {
    {
        echo "# Complete Big7 Experiment Batch"
        echo
        echo "- Batch root: \`${BATCH_ROOT}\`"
        echo "- Queue manifest: \`${QUEUE_TSV}\`"
        echo "- Status manifest: \`${STATUS_TSV}\`"
        echo "- Threads per job: \`${THREADS_PER_JOB}\`"
        echo "- Max parallel jobs: \`${MAX_PARALLEL}\`"
        echo "- Index threads per job: \`${INDEX_THREADS_PER_JOB}\`"
        echo "- Min free disk guard: \`${MIN_FREE_GB}G\`"
        echo "- Closeout command: \`cd ${PROJECT_ROOT} && ./scripts/analysis/run_complete_matrix_closeout.sh\`"
        echo
        echo "## Job Status"
        echo
        echo "| Job | Type | Sample | Mode | Status | Exit | Log |"
        echo "| --- | --- | --- | --- | --- | --- | --- |"
        local idx
        for idx in "${!JOB_IDS[@]}"; do
            echo "| ${JOB_IDS[$idx]} | ${JOB_TYPES[$idx]} | ${JOB_SAMPLES[$idx]} | ${JOB_MODES[$idx]} | ${JOB_STATUSES[$idx]} | ${JOB_EXITS[$idx]} | ${JOB_LOGS[$idx]} |"
        done
        if [[ ${#POSTPROCESS_NAMES[@]} -gt 0 ]]; then
            echo
            echo "## Postprocess"
            echo
            echo "| Name | Status | Output | Log |"
            echo "| --- | --- | --- | --- |"
            for idx in "${!POSTPROCESS_NAMES[@]}"; do
                echo "| ${POSTPROCESS_NAMES[$idx]} | ${POSTPROCESS_STATUSES[$idx]} | ${POSTPROCESS_OUTPUTS[$idx]} | ${POSTPROCESS_LOGS[$idx]} |"
            done
        fi
    } > "${SUMMARY_MD}"
}

for sample in "${PAIRED_SAMPLES[@]}"; do
    append_job \
        "paired_full_${sample}" \
        "paired" \
        "${sample}" \
        "paired_full" \
        "cd ${PROJECT_ROOT} && ./scripts/pipeline/run_benchmark.sh --sample ${sample} --mode s-pure --threads ${THREADS_PER_JOB} --index-threads ${INDEX_THREADS_PER_JOB} --min-free-gb ${MIN_FREE_GB} --run-tag ${RUN_TAG} --skip-cleanup" \
        "${BATCH_ROOT}/logs/paired_full_${sample}.log"

    append_job \
        "paired_pileup_${sample}" \
        "paired" \
        "${sample}" \
        "paired_pileup" \
        "cd ${PROJECT_ROOT} && ./scripts/pipeline/run_benchmark.sh --sample ${sample} --mode s-pure-pileup --vcf-source pileup --threads ${THREADS_PER_JOB} --index-threads ${INDEX_THREADS_PER_JOB} --min-free-gb ${MIN_FREE_GB} --run-tag ${RUN_TAG} --skip-cleanup" \
        "${BATCH_ROOT}/logs/paired_pileup_${sample}.log"
done

for sample in "${TO_SAMPLES[@]}"; do
    eval "$(get_sample_config "${sample}")"
    if [[ -z "${TO_SOMATIC_VCF_PILEUP:-}" ]]; then
        continue
    fi
    append_job \
        "to_pileup_${sample}" \
        "tumor_only" \
        "${sample}" \
        "to_pileup" \
        "cd ${PROJECT_ROOT} && ./scripts/analysis/run_longphase_to_intersubmod_pilot.sh --sample ${sample} --platform ${PLATFORM_LABEL} --tumor-bam ${TUMOR_BAM} --truth-vcf ${TRUTH_VCF} --truth-bed '${TRUTH_BED}' --truth-total ${TRUTH_TOTAL} --truth-set ${TRUTH_SET_LABEL} --caller-existing-vcf ${TO_SOMATIC_VCF_PILEUP} --threads ${THREADS_PER_JOB} --index-threads ${INDEX_THREADS_PER_JOB} --min-free-gb ${MIN_FREE_GB}" \
        "${BATCH_ROOT}/logs/to_pileup_${sample}.log"
done

write_queue_manifest
write_status_manifest

log_batch "Batch root: ${BATCH_ROOT}"
log_batch "Queue manifest: ${QUEUE_TSV}"
log_batch "Status manifest: ${STATUS_TSV}"
log_batch "Threads per job=${THREADS_PER_JOB}, index_threads=${INDEX_THREADS_PER_JOB}, max_parallel=${MAX_PARALLEL}, min_free_gb=${MIN_FREE_GB}"

if [[ "${DRY_RUN}" == true ]]; then
    local_idx=0
    for local_idx in "${!JOB_IDS[@]}"; do
        log_batch "DRY-RUN ${JOB_IDS[$local_idx]} => ${JOB_COMMANDS[$local_idx]}"
    done
    write_summary
    exit 0
fi

require_min_free_gb "/big7_disk" "${MIN_FREE_GB}" "complete_big7_batch"

for idx in "${!JOB_IDS[@]}"; do
    while (( ${#PID_TO_INDEX[@]} >= MAX_PARALLEL )); do
        reap_one_job
    done
    launch_job "${idx}"
done

while (( ${#PID_TO_INDEX[@]} > 0 )); do
    reap_one_job
done

paired_failed=false
for idx in "${!JOB_IDS[@]}"; do
    if [[ "${JOB_TYPES[$idx]}" == "paired" && "${JOB_STATUSES[$idx]}" != "completed" ]]; then
        paired_failed=true
        break
    fi
done

if [[ "${RUN_POSTPROCESS}" == true && "${paired_failed}" == false ]]; then
    FULL_ROUND_DIR="${SYNTHESIS_OUTPUT_ROOT}/research_rounds/${TIMESTAMP}_all_pure_paired_full_complete_matrix"
    PILEUP_ROUND_DIR="${SYNTHESIS_OUTPUT_ROOT}/research_rounds/${TIMESTAMP}_all_pure_paired_pileup_complete_matrix"

    run_postprocess_step \
        "paired_full_round" \
        "${FULL_ROUND_DIR}" \
        "cd ${PROJECT_ROOT} && ./scripts/analysis/run_pure_research_round.sh --sample-set all_pure --mode s-pure --output-dir ${FULL_ROUND_DIR} --threads ${THREADS}" \
        "${BATCH_ROOT}/logs/post_paired_full_round.log"

    run_postprocess_step \
        "paired_pileup_round" \
        "${PILEUP_ROUND_DIR}" \
        "cd ${PROJECT_ROOT} && ./scripts/analysis/run_pure_research_round.sh --sample-set all_pure --mode s-pure-pileup --output-dir ${PILEUP_ROUND_DIR} --threads ${THREADS}" \
        "${BATCH_ROOT}/logs/post_paired_pileup_round.log"
fi

write_postprocess_manifest
write_summary
log_batch "Batch summary: ${SUMMARY_MD}"
