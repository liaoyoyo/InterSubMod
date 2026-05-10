#!/usr/bin/env bash
# test_skill_usage_logger.sh - TDD test for skill_usage_logger.sh
# Tests PreToolUse hook that logs Skill tool invocations.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGGER="${SCRIPT_DIR}/../skill_usage_logger.sh"

PASS=0
FAIL=0

pass() {
    echo "PASS: $1"
    PASS=$((PASS + 1))
}

fail() {
    echo "FAIL: $1"
    FAIL=$((FAIL + 1))
}

run_test() {
    local name="$1"
    local payload="$2"
    local tmp_dir
    tmp_dir="$(mktemp -d)"
    local log_file="${tmp_dir}/skill_usage.log"
    local exit_code

    set +e
    echo "$payload" | SKILL_USAGE_LOG_DIR="$tmp_dir" bash "$LOGGER" 2>/dev/null
    exit_code=$?
    set -e

    echo "$tmp_dir|$log_file|$exit_code"
}

# ----------------------------------------------------------------------------
# Test 1: valid Skill payload writes correct log line
# ----------------------------------------------------------------------------
echo "--- Test 1: valid Skill payload ---"
T1_DIR="$(mktemp -d)"
T1_LOG="${T1_DIR}/skill_usage.log"
PAYLOAD_1='{"tool_input":{"skill":"confirmation-protocol"},"session_id":"test-sess-1"}'
echo "$PAYLOAD_1" | SKILL_USAGE_LOG_DIR="$T1_DIR" bash "$LOGGER"
T1_EXIT=$?

if [ "$T1_EXIT" -eq 0 ]; then
    pass "Test 1a: exit code 0"
else
    fail "Test 1a: exit code $T1_EXIT (expected 0)"
fi

if [ -f "$T1_LOG" ]; then
    pass "Test 1b: log file exists"
else
    fail "Test 1b: log file not created at $T1_LOG"
fi

if [ -f "$T1_LOG" ] && grep -q "confirmation-protocol" "$T1_LOG"; then
    pass "Test 1c: log contains skill name"
else
    fail "Test 1c: log missing skill name 'confirmation-protocol'"
fi

if [ -f "$T1_LOG" ] && grep -q "test-sess-1" "$T1_LOG"; then
    pass "Test 1d: log contains session_id"
else
    fail "Test 1d: log missing session_id 'test-sess-1'"
fi

if [ -f "$T1_LOG" ] && grep -qE '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z' "$T1_LOG"; then
    pass "Test 1e: log line has ISO8601 UTC timestamp"
else
    fail "Test 1e: log line missing ISO8601 UTC timestamp"
fi

rm -rf "$T1_DIR"

# ----------------------------------------------------------------------------
# Test 2: payload without skill field exits 0 and writes nothing
# ----------------------------------------------------------------------------
echo "--- Test 2: payload without skill field ---"
T2_DIR="$(mktemp -d)"
T2_LOG="${T2_DIR}/skill_usage.log"
PAYLOAD_2='{"tool_input":{"command":"ls"}}'
set +e
echo "$PAYLOAD_2" | SKILL_USAGE_LOG_DIR="$T2_DIR" bash "$LOGGER"
T2_EXIT=$?
set -e

if [ "$T2_EXIT" -eq 0 ]; then
    pass "Test 2a: exit code 0 on missing skill"
else
    fail "Test 2a: exit code $T2_EXIT (expected 0)"
fi

if [ ! -f "$T2_LOG" ] || [ ! -s "$T2_LOG" ]; then
    pass "Test 2b: no log line written for non-Skill payload"
else
    fail "Test 2b: log file unexpectedly written: $(cat "$T2_LOG")"
fi

rm -rf "$T2_DIR"

# ----------------------------------------------------------------------------
# Test 3: invalid JSON exits 0 (don't break Claude Code's hook chain)
# ----------------------------------------------------------------------------
echo "--- Test 3: invalid JSON ---"
T3_DIR="$(mktemp -d)"
set +e
echo "not json at all" | SKILL_USAGE_LOG_DIR="$T3_DIR" bash "$LOGGER" 2>/dev/null
T3_EXIT=$?
set -e

if [ "$T3_EXIT" -eq 0 ]; then
    pass "Test 3a: exit code 0 on invalid JSON"
else
    fail "Test 3a: exit code $T3_EXIT (expected 0)"
fi

rm -rf "$T3_DIR"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
echo ""
echo "================================"
echo "  PASS: $PASS"
echo "  FAIL: $FAIL"
echo "================================"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
