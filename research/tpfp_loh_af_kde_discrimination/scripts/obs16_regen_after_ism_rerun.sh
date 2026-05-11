#!/bin/bash
# obs16_regen_after_ism_rerun.sh
# Post-ISM-rerun driver: rebuild master dataset + all downstream figures.
#
# Runs in order (each step dep on previous):
#   1. step0_build_master.py      (ng_kde_rescaling)  -> updates master TSV
#   2. step4_tpfp_discrimination.py (ng_kde_rescaling) -> tpfp_cube_coarse.tsv
#   3. step4b_tpfp_tomode.py      (ng_kde_rescaling)  -> tpfp_cube_*_TO.tsv
#   4. step6_tpfp_detailed.py     (ng_kde_rescaling)  -> tpfp_per_sample_scheme_full.tsv
#   5. step7_visualize_detailed.py (ng_kde_rescaling) -> fig_v2_7, fig_v2_8...
#   6. obs09-12                    (tpfp_loh_af_kde_discrimination) -> L1/L2/L3/consistency

set -euo pipefail

PROJECT_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
KDE_DIR="${PROJECT_ROOT}/research/ng_kde_rescaling/scripts"
OBS_DIR="${PROJECT_ROOT}/research/tpfp_loh_af_kde_discrimination/scripts"

LOG_DIR="${PROJECT_ROOT}/research/tpfp_loh_af_kde_discrimination/logs/obs16_regen_20260422"
mkdir -p "${LOG_DIR}"

echo "==> obs16: downstream regen after ISM rerun"

cd "${PROJECT_ROOT}"

for step in \
    "${KDE_DIR}/step0_build_master.py" \
    "${KDE_DIR}/step4_tpfp_discrimination.py" \
    "${KDE_DIR}/step4b_tpfp_tomode.py" \
    "${KDE_DIR}/step6_tpfp_detailed.py" \
    "${KDE_DIR}/step7_visualize_detailed.py" \
    ; do
    name=$(basename "${step}")
    log="${LOG_DIR}/${name%.py}.log"
    echo "==> Running ${name}"
    if python3 "${step}" 2>&1 | tee "${log}"; then
        echo "    ✓ ${name} OK"
    else
        echo "    ✗ ${name} FAILED (see ${log})"
        exit 1
    fi
done

# obs09-12 will be invoked separately after user reviews fig_v2_7
echo
echo "==> Core regen complete. fig_v2_7 updated at:"
echo "    ${PROJECT_ROOT}/docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled/fig_v2_7_baseline_tp_fp_fold.png"
echo "    (symlinked in ${PROJECT_ROOT}/research/tpfp_loh_af_kde_discrimination/figures/existing/)"
