#!/usr/bin/env bash
# test_skill_usage_audit.sh - TDD test for skill_usage_audit.sh
# Tests weekly audit report that summarizes Skill tool usage.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUDIT="${SCRIPT_DIR}/../skill_usage_audit.sh"

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

# Helper: build a fake skills directory with N fake skill subdirs
mk_fake_skills_dir() {
    local dir="$1"
    shift
    mkdir -p "$dir"
    for name in "$@"; do
        mkdir -p "${dir}/${name}"
        : > "${dir}/${name}/SKILL.md"
    done
}

# Compute "now - days" in ISO8601 UTC (GNU date)
iso_days_ago() {
    local n="$1"
    date -u -d "${n} days ago" +"%Y-%m-%dT%H:%M:%SZ"
}

# ----------------------------------------------------------------------------
# Test 1: --days 7 filters out entries older than 7 days
# ----------------------------------------------------------------------------
echo "--- Test 1: --days 7 filters older entries ---"
T1_DIR="$(mktemp -d)"
T1_LOG="${T1_DIR}/skill_usage.log"
T1_SKILLS="${T1_DIR}/skills"
mk_fake_skills_dir "$T1_SKILLS" \
    confirmation-protocol methodology-audit results-report \
    auc-confound-guard check-staleness data-audit doc-standards \
    init-research

# Build 5 entries spanning 10 days; only the 3 most recent should be in the
# 7-day window.
TS_RECENT=$(iso_days_ago 0)
TS_2D=$(iso_days_ago 2)
TS_5D=$(iso_days_ago 5)
TS_8D=$(iso_days_ago 8)
TS_10D=$(iso_days_ago 10)

{
    echo "${TS_RECENT}|confirmation-protocol|sess-A|/tmp/x"
    echo "${TS_2D}|methodology-audit|sess-A|/tmp/x"
    echo "${TS_5D}|results-report|sess-B|/tmp/x"
    echo "${TS_8D}|auc-confound-guard|sess-C|/tmp/x"
    echo "${TS_10D}|doc-standards|sess-D|/tmp/x"
} > "$T1_LOG"

set +e
T1_OUT=$(SKILLS_DIR="$T1_SKILLS" bash "$AUDIT" --days 7 --log-file "$T1_LOG" 2>&1)
T1_EXIT=$?
set -e

if [ "$T1_EXIT" -eq 0 ]; then
    pass "Test 1a: exit code 0"
else
    fail "Test 1a: exit code $T1_EXIT (expected 0)"
fi

if echo "$T1_OUT" | grep -q "Total invocations: 3"; then
    pass "Test 1b: counts only 3 in-window entries"
else
    fail "Test 1b: expected 'Total invocations: 3', got: $(echo "$T1_OUT" | grep -E 'Total invocations' || echo NONE)"
fi

if echo "$T1_OUT" | grep -q "confirmation-protocol"; then
    pass "Test 1c: top-table includes recent skill"
else
    fail "Test 1c: top-table missing confirmation-protocol"
fi

if echo "$T1_OUT" | grep -qE "auc-confound-guard\s+[0-9]"; then
    fail "Test 1d: out-of-window skill listed in top-triggered"
else
    pass "Test 1d: 8-day-old entry excluded from top-triggered"
fi

if echo "$T1_OUT" | grep -qE "Unique skills triggered: 3 / 8"; then
    pass "Test 1e: unique-skills count 3 / 8 (fake skills dir size)"
else
    fail "Test 1e: expected 'Unique skills triggered: 3 / 8', got: $(echo "$T1_OUT" | grep -E 'Unique skills triggered' || echo NONE)"
fi

rm -rf "$T1_DIR"

# ----------------------------------------------------------------------------
# Test 2: empty log file prints empty message and exits 0
# ----------------------------------------------------------------------------
echo "--- Test 2: empty log file ---"
T2_DIR="$(mktemp -d)"
T2_LOG="${T2_DIR}/skill_usage.log"
T2_SKILLS="${T2_DIR}/skills"
mk_fake_skills_dir "$T2_SKILLS" foo bar baz
: > "$T2_LOG"  # touch empty file

set +e
T2_OUT=$(SKILLS_DIR="$T2_SKILLS" bash "$AUDIT" --log-file "$T2_LOG" 2>&1)
T2_EXIT=$?
set -e

if [ "$T2_EXIT" -eq 0 ]; then
    pass "Test 2a: exit code 0 on empty log"
else
    fail "Test 2a: exit code $T2_EXIT (expected 0)"
fi

if echo "$T2_OUT" | grep -qi "empty"; then
    pass "Test 2b: empty-log message printed"
else
    fail "Test 2b: expected 'empty' in output, got: $T2_OUT"
fi

rm -rf "$T2_DIR"

# ----------------------------------------------------------------------------
# Test 3: missing log file prints not-found message and exits 0
# ----------------------------------------------------------------------------
echo "--- Test 3: missing log file ---"
T3_DIR="$(mktemp -d)"
T3_LOG="${T3_DIR}/does_not_exist.log"
T3_SKILLS="${T3_DIR}/skills"
mk_fake_skills_dir "$T3_SKILLS" foo bar

set +e
T3_OUT=$(SKILLS_DIR="$T3_SKILLS" bash "$AUDIT" --log-file "$T3_LOG" 2>&1)
T3_EXIT=$?
set -e

if [ "$T3_EXIT" -eq 0 ]; then
    pass "Test 3a: exit code 0 on missing log"
else
    fail "Test 3a: exit code $T3_EXIT (expected 0)"
fi

if echo "$T3_OUT" | grep -qi "not found"; then
    pass "Test 3b: not-found message printed"
else
    fail "Test 3b: expected 'not found' in output, got: $T3_OUT"
fi

rm -rf "$T3_DIR"

# ----------------------------------------------------------------------------
# Test 4: 0-trigger skills section lists untriggered skills
# ----------------------------------------------------------------------------
echo "--- Test 4: 0-trigger skills section ---"
T4_DIR="$(mktemp -d)"
T4_LOG="${T4_DIR}/skill_usage.log"
T4_SKILLS="${T4_DIR}/skills"
mk_fake_skills_dir "$T4_SKILLS" \
    confirmation-protocol methodology-audit results-report \
    auc-confound-guard check-staleness data-audit doc-standards \
    init-research conclude-research review-evidence

TS_NOW=$(iso_days_ago 0)
{
    echo "${TS_NOW}|confirmation-protocol|sess-A|/tmp/x"
    echo "${TS_NOW}|methodology-audit|sess-A|/tmp/x"
} > "$T4_LOG"

set +e
T4_OUT=$(SKILLS_DIR="$T4_SKILLS" bash "$AUDIT" --days 7 --log-file "$T4_LOG" 2>&1)
T4_EXIT=$?
set -e

if [ "$T4_EXIT" -eq 0 ]; then
    pass "Test 4a: exit code 0"
else
    fail "Test 4a: exit code $T4_EXIT (expected 0)"
fi

if echo "$T4_OUT" | grep -qi "0-trigger skills"; then
    pass "Test 4b: 0-trigger section header present"
else
    fail "Test 4b: '0-trigger skills' header missing"
fi

# Expect the 8 untriggered skills to appear (10 total - 2 triggered = 8)
if echo "$T4_OUT" | grep -q "auc-confound-guard" \
    && echo "$T4_OUT" | grep -q "doc-standards" \
    && echo "$T4_OUT" | grep -q "review-evidence"; then
    pass "Test 4c: 0-trigger list includes untriggered skills"
else
    fail "Test 4c: 0-trigger list missing expected entries"
fi

# Check for "(8 total" in the section header
if echo "$T4_OUT" | grep -qE "0-trigger skills \(8 total"; then
    pass "Test 4d: 0-trigger header reports correct count (8)"
else
    fail "Test 4d: expected '(8 total' in 0-trigger header, got: $(echo "$T4_OUT" | grep -E '0-trigger' || echo NONE)"
fi

rm -rf "$T4_DIR"

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
