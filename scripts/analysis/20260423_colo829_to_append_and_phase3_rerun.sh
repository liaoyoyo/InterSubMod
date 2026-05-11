#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# COLO829 TO — append to merged master + rerun Phase 3 S1-S7 synthesis
#
# Purpose:
#   When the background-running COLO829 TO full pipeline finishes, this script
#   appends the new ISM TP+FP rows to the existing merged master TSV and re-runs
#   the Phase 3 synthesis script to produce a 7/7 TO cross-sample view.
#
# Idempotent:
#   - Detects existence + non-emptiness of step05 significance_summary.csv TP/FP
#   - Output TSV is date-stamped; rerunning overwrites it safely
#   - Temporary backup of original merged master is restored on success/failure
#
# Usage:
#   bash scripts/analysis/20260423_colo829_to_append_and_phase3_rerun.sh
#   bash scripts/analysis/20260423_colo829_to_append_and_phase3_rerun.sh --check-only
#
# Date: 2026-04-23
# ----------------------------------------------------------------------------
set -euo pipefail

# ----- Paths (absolute, following InterSubMod conventions) ------------------
readonly REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
readonly PILOT_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot"
readonly STEP05_DIR="${PILOT_ROOT}/step05_intersubmod"
readonly TP_CSV="${STEP05_DIR}/intersubmod_tp/significance_summary.csv"
readonly FP_CSV="${STEP05_DIR}/intersubmod_fp/significance_summary.csv"

readonly DATA_DIR="${REPO_ROOT}/research/ng_kde_rescaling/data"
readonly BASE_MASTER="${DATA_DIR}/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
readonly NEW_MASTER="${DATA_DIR}/merged_7samples_paired_full_plus_7to_full_20260423.tsv.gz"

readonly PHASE3_SCRIPT="${REPO_ROOT}/scripts/analysis/20260423_phase3_synthesis.py"
readonly PHASE3_OUT="${REPO_ROOT}/docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis"

readonly LOG_PREFIX="[colo829_to_append]"

# ----- Argument parsing -----------------------------------------------------
CHECK_ONLY="false"
if [[ "${1:-}" == "--check-only" ]]; then
    CHECK_ONLY="true"
fi

# ----- Step 1: Wait/check for pipeline output -------------------------------
echo "${LOG_PREFIX} Step 1: checking COLO829 TO pipeline outputs"
echo "  TP_CSV=${TP_CSV}"
echo "  FP_CSV=${FP_CSV}"

if [[ ! -f "${TP_CSV}" ]]; then
    echo "${LOG_PREFIX} ERROR: ${TP_CSV} not found" >&2
    echo "${LOG_PREFIX} Pipeline still running or failed. Check ${PILOT_ROOT}/run.log" >&2
    exit 1
fi
if [[ ! -f "${FP_CSV}" ]]; then
    echo "${LOG_PREFIX} ERROR: ${FP_CSV} not found" >&2
    exit 1
fi

# Non-empty check (at least 2 rows = header + 1 data)
tp_lines=$(wc -l < "${TP_CSV}")
fp_lines=$(wc -l < "${FP_CSV}")
if [[ "${tp_lines}" -lt 2 ]] || [[ "${fp_lines}" -lt 2 ]]; then
    echo "${LOG_PREFIX} ERROR: TP has ${tp_lines} lines, FP has ${fp_lines} lines (both must be >= 2)" >&2
    exit 1
fi
echo "${LOG_PREFIX} TP rows=$((tp_lines - 1)), FP rows=$((fp_lines - 1))"

if [[ "${CHECK_ONLY}" == "true" ]]; then
    echo "${LOG_PREFIX} --check-only: pipeline outputs present. Exiting."
    exit 0
fi

# ----- Step 2: Verify Diploid_Coverage_Used non-null (new KDE binary) -------
echo "${LOG_PREFIX} Step 2: verifying KDE-fixed binary output (Diploid_Coverage_Used)"
dcu_check=$(python3 -c "
import pandas as pd
tp = pd.read_csv('${TP_CSV}', nrows=100)
if 'Diploid_Coverage_Used' not in tp.columns:
    print('MISSING')
    exit(1)
vals = tp['Diploid_Coverage_Used'].dropna().unique()
if len(vals) == 0:
    print('EMPTY')
    exit(1)
print(f'OK baseline={vals[0]}')
")
echo "  ${dcu_check}"

# ----- Step 3: Append + produce new master ---------------------------------
echo "${LOG_PREFIX} Step 3: appending COLO829 TO rows to new master"

python3 <<'PYEOF'
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Paths (mirror bash vars)
REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
DATA_DIR = REPO_ROOT / "research/ng_kde_rescaling/data"
TP_CSV = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod/intersubmod_tp/significance_summary.csv")
FP_CSV = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod/intersubmod_fp/significance_summary.csv")
BASE_MASTER = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"
NEW_MASTER = DATA_DIR / "merged_7samples_paired_full_plus_7to_full_20260423.tsv.gz"

# Import utils_io for load_sample_run + SAMPLE_ORDER + loh_side
sys.path.insert(0, str(REPO_ROOT / "research/ng_kde_rescaling/scripts"))
from utils_io import load_sample_run, SAMPLE_ORDER, loh_side

print("[py] loading base master:", BASE_MASTER)
base = pd.read_csv(BASE_MASTER, sep="\t", low_memory=False)
print(f"[py]   base rows={len(base)}, cols={len(base.columns)}")
print(f"[py]   existing (sample,mode) pairs: {sorted(base[['sample','mode']].drop_duplicates().apply(tuple, axis=1).tolist())}")

# Safety: if COLO829/to_pileup already in base, skip to avoid duplication
already = ((base["sample"] == "COLO829") & (base["mode"] == "to_pileup")).sum()
if already > 0:
    print(f"[py] WARNING: base master already has {already} COLO829 to_pileup rows; overwriting them.")
    base = base[~((base["sample"] == "COLO829") & (base["mode"] == "to_pileup"))]

print("[py] loading COLO829 TO TP+FP via load_sample_run")
df_co = load_sample_run(TP_CSV, FP_CSV, "COLO829", "to_pileup", "new_direct_20260423")

# Derive analysis columns (mirror step0_build_master.py § "Derive analysis columns")
df_co["CovM_used"] = df_co["Coverage_Multiple"].astype(float)
if "Diploid_Coverage_Used" in df_co.columns:
    df_co["baseline_x"] = float(df_co["Diploid_Coverage_Used"].median())
else:
    df_co["baseline_x"] = np.nan
df_co["sample_order"] = SAMPLE_ORDER.index("COLO829") if "COLO829" in SAMPLE_ORDER else 6
df_co["loh_side"] = df_co.apply(loh_side, axis=1)
# AlleleDelta == AF (per step0_build_master)
df_co["AF"] = df_co["AlleleDelta"].clip(lower=0.0, upper=1.0)
df_co["AF_bin10"] = pd.cut(df_co["AF"], bins=np.arange(0.0, 1.01, 0.1), include_lowest=True)
df_co["AF_class"] = df_co["AF"].apply(
    lambda x: "Extreme" if (pd.notna(x) and (x < 0.1 or x > 0.9))
    else ("Near-half" if (pd.notna(x) and 0.4 <= x <= 0.6)
          else ("Intermediate" if pd.notna(x) else "Unknown"))
)

print(f"[py] COLO829 TO rows={len(df_co)}, baseline_x={df_co['baseline_x'].iloc[0]}, kde_status={df_co['kde_status'].iloc[0]}")
print(f"[py] AF_class distribution: {df_co['AF_class'].value_counts().to_dict()}")
print(f"[py] HPFineNGroups distribution: {df_co['HPFineNGroups'].value_counts().sort_index().head(10).to_dict()}")

# Align columns: union of base + new, fill missing with NaN
base_cols = list(base.columns)
new_cols = list(df_co.columns)
union_cols = base_cols + [c for c in new_cols if c not in base_cols]

for c in union_cols:
    if c not in base.columns:
        base[c] = np.nan
    if c not in df_co.columns:
        df_co[c] = np.nan

# Reorder both to the union
base = base[union_cols]
df_co = df_co[union_cols]

# Row-append
combined = pd.concat([base, df_co], ignore_index=True)
print(f"[py] combined rows={len(combined)} (base {len(base)} + COLO829_TO {len(df_co)})")

# Sanity: (sample, mode) coverage
coverage = combined[["sample", "mode"]].drop_duplicates().sort_values(["sample", "mode"])
print("[py] final (sample, mode) coverage:")
print(coverage.to_string(index=False))

# Write new master (gzip-compressed TSV)
combined.to_csv(NEW_MASTER, sep="\t", index=False, compression="gzip")
print(f"[py] wrote {NEW_MASTER}  ({NEW_MASTER.stat().st_size / 1e6:.1f} MB)")
PYEOF

echo "${LOG_PREFIX} Step 3 complete: ${NEW_MASTER}"

# ----- Step 4: Run Phase 3 synthesis ----------------------------------------
# Strategy: phase3_synthesis.py hard-codes SRC to the base master path. We
# temporarily swap in the new master via symlink (backup+restore), then
# run phase3. Script remains idempotent because both files are immutable inputs.
echo "${LOG_PREFIX} Step 4: running Phase 3 synthesis against new master"

BASE_BAK="${BASE_MASTER}.bak.colo829_append_$(date +%Y%m%d_%H%M%S)"
echo "${LOG_PREFIX}   backing up base master -> ${BASE_BAK}"
cp -p "${BASE_MASTER}" "${BASE_BAK}"

# Cleanup restore on exit (success or failure)
restore_master() {
    if [[ -f "${BASE_BAK}" ]]; then
        echo "${LOG_PREFIX}   restoring base master from ${BASE_BAK}"
        mv -f "${BASE_BAK}" "${BASE_MASTER}"
    fi
}
trap restore_master EXIT

# Swap in new master
cp -f "${NEW_MASTER}" "${BASE_MASTER}"
echo "${LOG_PREFIX}   swapped in ${NEW_MASTER} at ${BASE_MASTER}"

# Run phase3 synthesis
echo "${LOG_PREFIX}   invoking ${PHASE3_SCRIPT}"
python3 "${PHASE3_SCRIPT}"

# Restore handled by trap
echo "${LOG_PREFIX} Step 4 complete: phase3 outputs under ${PHASE3_OUT}"

# ----- Step 5: Report --------------------------------------------------------
echo ""
echo "${LOG_PREFIX} ============================================================"
echo "${LOG_PREFIX} DONE — summary"
echo "${LOG_PREFIX} ------------------------------------------------------------"
echo "${LOG_PREFIX}   new master  : ${NEW_MASTER}"
echo "${LOG_PREFIX}   phase3 out  : ${PHASE3_OUT}/summary.json"
echo "${LOG_PREFIX}                 ${PHASE3_OUT}/s1s7_per_sample.tsv"
echo "${LOG_PREFIX}                 ${PHASE3_OUT}/s1s7_heatmap_tp_rate.png"
echo "${LOG_PREFIX}                 ${PHASE3_OUT}/s1s7_heatmap_fold.png"
echo "${LOG_PREFIX}                 ${PHASE3_OUT}/s3_cross_sample_wilcoxon.json"
echo "${LOG_PREFIX}                 ${PHASE3_OUT}/ng_gap_aggregate.json"
echo "${LOG_PREFIX} ============================================================"
echo "${LOG_PREFIX} COLO829 TO S1-S7 contribution:"
python3 - <<PYEOF
import json
from pathlib import Path
p = Path("${PHASE3_OUT}/s1s7_per_sample.tsv")
if p.exists():
    import pandas as pd
    df = pd.read_csv(p, sep="\t")
    co = df[(df["sample"] == "COLO829") & (df["mode"] == "to_pileup")]
    if len(co) == 0:
        print("  (no COLO829/to_pileup rows in phase3 output)")
    else:
        print(co[["scheme","n","k_tp","tp_rate","fold_vs_baseline"]].to_string(index=False))
    # vs cohort TO
    print("")
    print("Other TO samples per-scheme median for reference:")
    to_other = df[(df["mode"] == "to_pileup") & (df["sample"] != "COLO829")]
    if len(to_other) > 0:
        pivot = to_other.groupby("scheme")["tp_rate"].median()
        print(pivot.to_string())
PYEOF

echo "${LOG_PREFIX} DONE"
