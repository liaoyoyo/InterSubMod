#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP_DIR="$(cd "${TEST_DIR}/.." && pwd)"
BUILD_DIR="${TEST_DIR}/build"
mkdir -p "${BUILD_DIR}"

CXX="${CXX:-g++}"
CXXFLAGS=(-std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror)

"${CXX}" "${CXXFLAGS[@]}" \
    "${TEST_DIR}/test_strict_endpoint_graph.cpp" \
    -o "${BUILD_DIR}/test_strict_endpoint_graph"

"${CXX}" "${CXXFLAGS[@]}" \
    "${CPP_DIR}/strict_endpoint_graph_verify.cpp" \
    -o "${BUILD_DIR}/strict_endpoint_graph_verify"

"${BUILD_DIR}/test_strict_endpoint_graph"

"${BUILD_DIR}/strict_endpoint_graph_verify" \
    --input "${TEST_DIR}/strict_endpoint_graph.fixture" \
    --threshold 2 \
    --edges-output "${BUILD_DIR}/edges_threshold2.tsv" \
    --components-output "${BUILD_DIR}/components_threshold2.tsv" \
    --receipt-output "${BUILD_DIR}/receipt_threshold2.json"

diff -u "${TEST_DIR}/edges_threshold2.expected" "${BUILD_DIR}/edges_threshold2.tsv"
diff -u "${TEST_DIR}/components_threshold2.expected" "${BUILD_DIR}/components_threshold2.tsv"

"${BUILD_DIR}/strict_endpoint_graph_verify" \
    --input "${TEST_DIR}/strict_endpoint_graph.fixture" \
    --threshold 3 \
    --edges-output "${BUILD_DIR}/edges_threshold3.tsv" \
    --components-output "${BUILD_DIR}/components_threshold3.tsv" \
    --receipt-output "${BUILD_DIR}/receipt_threshold3.json"

awk -F '\t' 'NR > 1 && $15 != 0 { exit 1 } END { if (NR != 7) exit 1 }' "${BUILD_DIR}/edges_threshold3.tsv"
awk -F '\t' 'NR > 1 && ($7 != 1 || $10 != 0 || $12 != "false" || $13 != "ABSTAIN_SINGLETON_UNLINKED") { exit 1 } END { if (NR != 12) exit 1 }' \
    "${BUILD_DIR}/components_threshold3.tsv"

python3 - "${BUILD_DIR}/receipt_threshold2.json" <<'PY'
import json
import pathlib
import sys

receipt = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert receipt["threshold"] == 2
assert receipt["counts"]["containers"] == 4
assert receipt["counts"]["qualifying_edges"] == 6
assert receipt["counts"]["components"] == 5
assert receipt["invariants"]["only_fixed_RA_endpoints_add_support"] is True
PY

if "${BUILD_DIR}/strict_endpoint_graph_verify" \
    --input "${TEST_DIR}/duplicate_molecule.fixture" \
    --threshold 2 \
    --edges-output "${BUILD_DIR}/duplicate_edges.tsv" \
    --components-output "${BUILD_DIR}/duplicate_components.tsv" \
    --receipt-output "${BUILD_DIR}/duplicate_receipt.json" \
    >"${BUILD_DIR}/duplicate_stdout.log" 2>"${BUILD_DIR}/duplicate_stderr.log"; then
    echo "duplicate molecule fixture unexpectedly passed" >&2
    exit 1
fi
grep -F "duplicate molecule_id: m01" "${BUILD_DIR}/duplicate_stderr.log" >/dev/null

echo "run_tests.sh: PASS"
echo "binary=${BUILD_DIR}/strict_endpoint_graph_verify"
echo "edges=${BUILD_DIR}/edges_threshold2.tsv"
echo "components=${BUILD_DIR}/components_threshold2.tsv"
echo "receipt=${BUILD_DIR}/receipt_threshold2.json"
