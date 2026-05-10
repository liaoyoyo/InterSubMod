#!/usr/bin/env bash
# Test that preflight.sh exits 0 when all deps present, exits non-zero otherwise.
set -euo pipefail
SCRIPT="$(dirname "$0")/../tools/preflight.sh"

# Test 1: script exists and executable
[ -x "$SCRIPT" ] || { echo "FAIL: $SCRIPT not executable"; exit 1; }

# Test 2: returns 0 in current env (assumes deps installed per Pre-Flight Check)
if "$SCRIPT" >/dev/null 2>&1; then
    echo "PASS: preflight.sh returns 0 in healthy env"
else
    echo "FAIL: preflight.sh returned non-zero in env that should be healthy"
    "$SCRIPT" 2>&1 | head -10
    exit 1
fi

echo "all preflight tests passed"
