#!/bin/bash
# obs15b_resume_failed_pilots.sh
# Recovery driver after obs15 was killed by session compaction at ~12:06.
#
# Status at recovery time (2026-04-22 12:14):
#   DORADO                : TP ✓ + FP ✓ (fully done, skip)
#   HCC1937               : TP ✓, FP partial  → FP-only rerun + downstream
#   HCC1954               : TP ✓, FP partial  → FP-only rerun + downstream
#   H1437                 : TP partial        → full resume_to_pilot
#   H2009                 : TP partial        → full resume_to_pilot
#
# Detachment: caller must invoke with `setsid nohup ... & disown` so this
# script survives parent shell exit (session compaction).

set -euo pipefail

PROJECT_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
ARCHIVE_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots"
ISM_BIN="${PROJECT_ROOT}/build/bin/inter_sub_mod"
REFERENCE="/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"
RESUME_SCRIPT="${PROJECT_ROOT}/scripts/analysis/resume_to_pilot_from_tagged_outputs.sh"
EXTRACT_METRICS="${PROJECT_ROOT}/scripts/analysis/extract_label_first_metrics.py"

LOG_DIR="${PROJECT_ROOT}/research/tpfp_loh_af_kde_discrimination/logs/obs15b_resume_20260422"
mkdir -p "${LOG_DIR}"

THREADS=9
WIN=5000
METRIC=BERNOULLI

# ----------------------------------------------------------------------
# Helper: FP-only rerun for pilots that have TP but not FP
# ----------------------------------------------------------------------
run_fp_only() {
    local pilot="$1"
    local run_dir="${ARCHIVE_ROOT}/${pilot}"
    local step05="${run_dir}/step05_intersubmod"
    local tagged_bam="${run_dir}/step03_longphase_to/tumor_tagged.bam"
    local fp_vcf="${run_dir}/step04_benchmark_longphase_to/filtered_snv_fp.vcf.gz"
    local fp_out="${step05}/intersubmod_fp"
    local fp_log="${step05}/intersubmod_fp.log"
    local logfile="${LOG_DIR}/${pilot}.fp_only.log"

    {
        echo "[$(date '+%F %T')] FP-only rerun: ${pilot}"
        # Wipe partial FP output
        if [[ -d "${fp_out}" ]]; then
            echo "  removing partial ${fp_out}"
            rm -rf "${fp_out}"
        fi
        mkdir -p "${fp_out}"

        echo "  launching inter_sub_mod FP..."
        /usr/bin/time -v "${ISM_BIN}" \
            --tumor-bam "${tagged_bam}" \
            --reference "${REFERENCE}" \
            --threads "${THREADS}" \
            --window-size "${WIN}" \
            --distance-metric "${METRIC}" \
            --log-level info \
            --output-filtered-reads \
            --vcf "${fp_vcf}" \
            --output-dir "${fp_out}" \
            > "${fp_log}" 2>&1

        if [[ ! -f "${fp_out}/significance_summary.csv" ]]; then
            echo "  ✗ FAILED: no significance_summary.csv"
            return 1
        fi

        local n_fp
        n_fp=$(tail -n +2 "${fp_out}/significance_summary.csv" | wc -l)
        echo "  ✓ FP complete: ${n_fp} regions"

        echo "  extracting label_first_metrics..."
        python3 "${EXTRACT_METRICS}" \
            --summary-csv "${fp_out}/significance_summary.csv" \
            --output-tsv "${fp_out}/label_first_metrics.tsv"

        echo "  FP-only rerun done for ${pilot}"
    } > "${logfile}" 2>&1
    return $?
}

# ----------------------------------------------------------------------
# Helper: full resume_to_pilot for pilots that need TP+FP
# ----------------------------------------------------------------------
run_full_resume() {
    local pilot="$1"
    local run_dir="${ARCHIVE_ROOT}/${pilot}"
    local step05="${run_dir}/step05_intersubmod"
    local logfile="${LOG_DIR}/${pilot}.full_resume.log"

    {
        echo "[$(date '+%F %T')] Full resume: ${pilot}"
        # Wipe any partial outputs so resume restarts cleanly
        for sub in intersubmod_tp intersubmod_fp; do
            if [[ -d "${step05}/${sub}" && ! -f "${step05}/${sub}/significance_summary.csv" ]]; then
                echo "  removing partial ${step05}/${sub}"
                rm -rf "${step05}/${sub}"
            fi
        done

        bash "${RESUME_SCRIPT}" \
            --run-dir "${run_dir}" \
            --threads "${THREADS}" \
            --skip-method-validation \
            --skip-dashboard

        echo "  ✓ full resume done for ${pilot}"
    } > "${logfile}" 2>&1
    return $?
}

# ======================================================================
# Launch all 4 in parallel
# ======================================================================
echo "==> obs15b: parallel resume at $(date '+%F %T')"
echo "    Log dir: ${LOG_DIR}"

declare -A PIDS

# FP-only pilots
for pilot in 20260318_hcc1937_to_pilot_fastresume 20260318_hcc1954_to_pilot_fastresume; do
    run_fp_only "${pilot}" &
    PIDS["${pilot}"]=$!
    echo "  launched FP-only ${pilot} (PID ${PIDS[${pilot}]})"
done

# Full-resume pilots
for pilot in 20260318_h1437_to_pilot_fastresume 20260318_h2009_to_pilot_fastresume; do
    run_full_resume "${pilot}" &
    PIDS["${pilot}"]=$!
    echo "  launched full ${pilot} (PID ${PIDS[${pilot}]})"
done

echo
echo "==> Waiting for all 4 pilots"
FAILED=()
for pilot in "${!PIDS[@]}"; do
    pid="${PIDS[${pilot}]}"
    if wait "${pid}"; then
        echo "  ✓ ${pilot} done"
    else
        echo "  ✗ ${pilot} FAILED (pid ${pid})"
        FAILED+=("${pilot}")
    fi
done

echo
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo "==> ALL 4 PILOTS RECOVERED at $(date '+%F %T')"
    exit 0
else
    echo "==> FAILURES: ${FAILED[*]}"
    exit 1
fi
