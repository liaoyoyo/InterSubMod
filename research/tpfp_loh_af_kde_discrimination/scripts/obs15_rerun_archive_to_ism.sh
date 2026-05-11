#!/bin/bash
# obs15_rerun_archive_to_ism.sh
# Re-run ISM pipeline for 5 archive TO pilots to regenerate significance_summary.csv
# with KDE (Diploid_Coverage_Used) + Phase 2 features (NHP_Somatic, Fisher_Frac_Sig, etc.)
#
# Strategy:
#   1. Backup each pilot's existing intersubmod_*, metrics.json to *_archive_v1_20260422
#   2. Run resume_to_pilot_from_tagged_outputs.sh in parallel (5 pilots, 9 threads each)
#   3. --skip-method-validation --skip-dashboard to keep rerun focused on ISM + filter analysis

set -euo pipefail

ARCHIVE_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots"
PROJECT_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
RESUME_SCRIPT="${PROJECT_ROOT}/scripts/analysis/resume_to_pilot_from_tagged_outputs.sh"

PILOTS=(
    "20260315_hcc1395_dorado_to_pilot"
    "20260318_h1437_to_pilot_fastresume"
    "20260318_h2009_to_pilot_fastresume"
    "20260318_hcc1937_to_pilot_fastresume"
    "20260318_hcc1954_to_pilot_fastresume"
)

THREADS_PER_JOB=9
LOG_DIR="${PROJECT_ROOT}/research/tpfp_loh_af_kde_discrimination/logs/obs15_ism_rerun_20260422"
BACKUP_SUFFIX="_archive_v1_20260422"
mkdir -p "${LOG_DIR}"

echo "==> obs15: ISM rerun for 5 archive TO pilots"
echo "    Threads per job: ${THREADS_PER_JOB} (total ~$((THREADS_PER_JOB * 5)))"
echo "    Log dir: ${LOG_DIR}"
echo

# ============================================================================
# Phase 1: Backup existing artifacts (force rerun)
# ============================================================================
echo "==> Phase 1: Backup existing artifacts"
for pilot in "${PILOTS[@]}"; do
    run_dir="${ARCHIVE_ROOT}/${pilot}"
    step05="${run_dir}/step05_intersubmod"
    echo "  [${pilot}]"
    # Backup intersubmod_tp and intersubmod_fp
    for sub in intersubmod_tp intersubmod_fp; do
        if [[ -d "${step05}/${sub}" ]]; then
            target="${step05}/${sub}${BACKUP_SUFFIX}"
            if [[ -d "${target}" ]]; then
                echo "    Already backed up: ${sub}${BACKUP_SUFFIX}"
            else
                mv "${step05}/${sub}" "${target}"
                echo "    ✓ backed up ${sub} → ${sub}${BACKUP_SUFFIX}"
            fi
        fi
    done
    # Backup metrics.json, benchmark_comparison, round_summary, method_design_validation etc.
    for f in metrics.json benchmark_comparison.tsv benchmark_comparison.md round_summary.md \
             to_stage_metrics.json method_design_validation.tsv label_cluster_agreement.tsv; do
        src="${run_dir}/${f}"
        if [[ -f "${src}" ]]; then
            # Construct backup name: file.ext.backup or file_backup.ext
            base="${f%.*}"
            ext="${f##*.}"
            target="${run_dir}/${base}${BACKUP_SUFFIX}.${ext}"
            if [[ -f "${target}" ]]; then
                echo "    Already backed up: ${f}"
            else
                mv "${src}" "${target}"
                echo "    ✓ backed up ${f}"
            fi
        fi
    done
done
echo

# ============================================================================
# Phase 2: Parallel rerun
# ============================================================================
echo "==> Phase 2: Parallel ISM rerun (${THREADS_PER_JOB} threads each)"
declare -A PIDS
for pilot in "${PILOTS[@]}"; do
    run_dir="${ARCHIVE_ROOT}/${pilot}"
    logfile="${LOG_DIR}/${pilot}.log"
    echo "  launching ${pilot} (log: ${logfile})"
    bash "${RESUME_SCRIPT}" \
        --run-dir "${run_dir}" \
        --threads "${THREADS_PER_JOB}" \
        --skip-method-validation \
        --skip-dashboard \
        > "${logfile}" 2>&1 &
    PIDS["${pilot}"]=$!
done

echo "  PIDs: ${PIDS[@]}"
echo

# ============================================================================
# Phase 3: Wait and collect exit codes
# ============================================================================
echo "==> Phase 3: Waiting for all jobs to complete"
FAILED=()
for pilot in "${PILOTS[@]}"; do
    pid="${PIDS[${pilot}]}"
    if wait "${pid}"; then
        echo "  ✓ ${pilot} complete"
    else
        echo "  ✗ ${pilot} FAILED (exit $?)"
        FAILED+=("${pilot}")
    fi
done

echo
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo "==> ALL 5 PILOTS SUCCESSFUL"
    echo "    New significance_summary.csv paths:"
    for pilot in "${PILOTS[@]}"; do
        echo "      ${ARCHIVE_ROOT}/${pilot}/step05_intersubmod/intersubmod_tp/significance_summary.csv"
    done
    exit 0
else
    echo "==> FAILURES: ${FAILED[*]}"
    exit 1
fi
