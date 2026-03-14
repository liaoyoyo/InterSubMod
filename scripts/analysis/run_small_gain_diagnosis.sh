#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-/big8_disk/liaoyoyo2001/InterSubMod_runs/output/pure_tumor_evaluation/pure_tumor_eval_20260305_143418}"
DATA_DIR="${BASE_DIR}/data"

if [[ ! -f "${DATA_DIR}/per_variant_decision.tsv.gz" ]]; then
    echo "[ERROR] Missing ${DATA_DIR}/per_variant_decision.tsv.gz" >&2
    exit 1
fi
if [[ ! -f "${DATA_DIR}/run_manifest.tsv" ]]; then
    echo "[ERROR] Missing ${DATA_DIR}/run_manifest.tsv" >&2
    exit 1
fi

python scripts/analysis/small_gain_diagnosis.py \
    --input "${DATA_DIR}/per_variant_decision.tsv.gz" \
    --manifest "${DATA_DIR}/run_manifest.tsv" \
    --out-dir "${DATA_DIR}"

echo "[DONE] small_gain_diagnosis completed under ${DATA_DIR}"

