#!/bin/bash
# run_pure_tumor_evaluation.sh - One-command wrapper for pure tumor methylation evaluation.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/pure_tumor_methylation_evaluation.py"

if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "[ERROR] Python script not found: ${PY_SCRIPT}" >&2
    exit 1
fi

python3 "${PY_SCRIPT}" "$@"
