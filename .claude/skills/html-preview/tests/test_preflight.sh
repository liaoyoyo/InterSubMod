#!/usr/bin/env bash
# Test that preflight.sh exits 0 when all deps present.
set -euo pipefail
SCRIPT="$(dirname "$0")/../tools/preflight.sh"

[ -x "$SCRIPT" ] || { echo "FAIL: $SCRIPT not executable"; exit 1; }

if "$SCRIPT" >/dev/null 2>&1; then
    echo "PASS: preflight.sh returns 0 in healthy env"
else
    echo "FAIL: preflight.sh returned non-zero"
    "$SCRIPT" 2>&1 | head -10
    exit 1
fi

echo "all preflight tests passed"
