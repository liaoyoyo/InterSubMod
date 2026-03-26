#!/bin/bash
# Build the final big7 complete-matrix closeout package.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/scripts/pipeline/config.sh"

OUTPUT_DIR="${SYNTHESIS_OUTPUT_ROOT}/final_closeout"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

python3 "${SCRIPT_DIR}/build_complete_matrix_closeout.py" \
    --canonical-root "${CANONICAL_OUTPUT_ROOT}" \
    --synthesis-root "${SYNTHESIS_OUTPUT_ROOT}" \
    --output-dir "${OUTPUT_DIR}" \
    "${EXTRA_ARGS[@]}"
