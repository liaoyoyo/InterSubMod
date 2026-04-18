#!/bin/bash
# build_and_test.sh - Build InterSubMod and run unit tests
#
# Wraps cmake + make + ctest into a single validation step.
# Returns exit code 0 only if both build and tests pass.
#
# Usage:
#   ./scripts/validation/runners/build_and_test.sh
#   ./scripts/validation/runners/build_and_test.sh --skip-tests
#   ./scripts/validation/runners/build_and_test.sh --build-type Debug

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../validate_config.sh"

# ============================================================================
# Argument Parsing
# ============================================================================

SKIP_TESTS=false
LOCAL_BUILD_TYPE="${BUILD_TYPE}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-tests)       SKIP_TESTS=true; shift ;;
        --build-type)       LOCAL_BUILD_TYPE="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: build_and_test.sh [--skip-tests] [--build-type Release|Debug]"
            exit 0
            ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ============================================================================
# Build
# ============================================================================

log_validate "Starting build (${LOCAL_BUILD_TYPE})..."
BUILD_START=$(date +%s)

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake "${PROJECT_ROOT}" \
    -DCMAKE_BUILD_TYPE="${LOCAL_BUILD_TYPE}" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    2>&1 | tail -5

NPROC=$(nproc 2>/dev/null || echo 8)
make -j"${NPROC}" 2>&1 | tail -20

BUILD_END=$(date +%s)
BUILD_TIME=$(( BUILD_END - BUILD_START ))

# Verify binary exists
ISM_BIN="$(resolve_intersubmod_bin)"
if [[ -z "${ISM_BIN}" ]]; then
    log_validate "BUILD FAILED: No InterSubMod binary found"
    exit 1
fi

log_validate "Build succeeded in ${BUILD_TIME}s: ${ISM_BIN}"

# ============================================================================
# Unit Tests
# ============================================================================

if [[ "${SKIP_TESTS}" == "true" ]]; then
    log_validate "Skipping unit tests (--skip-tests)"
    echo "BUILD_STATUS=success"
    echo "BUILD_TIME=${BUILD_TIME}"
    echo "TEST_STATUS=skipped"
    echo "TEST_COUNT=0"
    exit 0
fi

log_validate "Running unit tests..."
TEST_START=$(date +%s)

cd "${BUILD_DIR}"
if ctest --output-on-failure --timeout 120 2>&1; then
    TEST_STATUS="pass"
else
    TEST_STATUS="fail"
fi

TEST_END=$(date +%s)
TEST_TIME=$(( TEST_END - TEST_START ))

# Count tests
TEST_COUNT=$(ctest -N 2>/dev/null | grep -c "Test #" || echo "0")

log_validate "Tests: ${TEST_STATUS} (${TEST_COUNT} tests in ${TEST_TIME}s)"

# Output summary
echo "BUILD_STATUS=success"
echo "BUILD_TIME=${BUILD_TIME}"
echo "TEST_STATUS=${TEST_STATUS}"
echo "TEST_TIME=${TEST_TIME}"
echo "TEST_COUNT=${TEST_COUNT}"

if [[ "${TEST_STATUS}" != "pass" ]]; then
    log_validate "UNIT TESTS FAILED — aborting validation"
    exit 1
fi
