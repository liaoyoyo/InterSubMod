#!/usr/bin/env bash
# Test suite for context_size_audit.sh (TDD)
# Verifies default mode, snapshot mode, and graceful handling of missing files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUDIT_SCRIPT="${SCRIPT_DIR}/../context_size_audit.sh"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

PASS=0
FAIL=0

pass() { PASS=$((PASS + 1)); echo "  PASS: $1"; }
fail() { FAIL=$((FAIL + 1)); echo "  FAIL: $1" >&2; }

check() {
    local name="$1"
    local cond="$2"
    if eval "$cond"; then
        pass "$name"
    else
        fail "$name (cond: $cond)"
    fi
}

echo "== Test 1: default mode runs, exits 0, output contains TOTAL =="
OUT=$(bash "$AUDIT_SCRIPT" 2>&1) && RC=0 || RC=$?
check "exit code 0" "[ $RC -eq 0 ]"
check "output contains TOTAL" "echo \"\$OUT\" | grep -q 'TOTAL'"

echo "== Test 2: output contains all 6 mandatory files (CURRENT_FOCUS) =="
check "mentions CURRENT_FOCUS" "echo \"\$OUT\" | grep -q 'CURRENT_FOCUS'"
check "mentions experiments/INDEX" "echo \"\$OUT\" | grep -q 'experiments/INDEX'"
check "mentions docs/README" "echo \"\$OUT\" | grep -q 'docs/README'"

echo "== Test 3: output contains landscape files =="
check "mentions research_landscape" "echo \"\$OUT\" | grep -q 'research_landscape'"

echo "== Test 4: --snapshot mode creates JSON file in .claude/state/ =="
LABEL="testlabel_$$"
STATE_DIR="${REPO_ROOT}/.claude/state"
BEFORE_COUNT=$(find "$STATE_DIR" -maxdepth 1 -name "context_audit_${LABEL}_*.json" 2>/dev/null | wc -l)
SNAP_OUT=$(bash "$AUDIT_SCRIPT" --snapshot "$LABEL" 2>&1) && SNAP_RC=0 || SNAP_RC=$?
check "snapshot exit code 0" "[ $SNAP_RC -eq 0 ]"
AFTER_COUNT=$(find "$STATE_DIR" -maxdepth 1 -name "context_audit_${LABEL}_*.json" 2>/dev/null | wc -l)
check "JSON file created" "[ $AFTER_COUNT -gt $BEFORE_COUNT ]"

JSON_FILE=$(find "$STATE_DIR" -maxdepth 1 -name "context_audit_${LABEL}_*.json" 2>/dev/null | sort | tail -1)

echo "== Test 5: JSON file is valid =="
if [ -n "$JSON_FILE" ] && [ -f "$JSON_FILE" ]; then
    if python3 -m json.tool "$JSON_FILE" > /dev/null 2>&1; then
        pass "JSON valid"
    else
        fail "JSON invalid: $JSON_FILE"
    fi
    # Cleanup test artifact
    rm -f "$JSON_FILE"
else
    fail "no JSON file located for cleanup/validation"
fi

echo "== Test 6: missing CLAUDE.md path → reports MISSING but exits 0 =="
FAKE_PATH="/tmp/nonexistent_claude_md_$$.md"
MISS_OUT=$(CLAUDE_MD_OVERRIDE="$FAKE_PATH" bash "$AUDIT_SCRIPT" 2>&1) && MISS_RC=0 || MISS_RC=$?
check "missing-file exit code 0" "[ $MISS_RC -eq 0 ]"
check "output contains MISSING" "echo \"\$MISS_OUT\" | grep -q 'MISSING'"

echo
echo "== SUMMARY =="
echo "PASS: $PASS  FAIL: $FAIL"
[ $FAIL -eq 0 ]
