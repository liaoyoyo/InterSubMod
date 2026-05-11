#!/usr/bin/env bash
# skill_usage_audit.sh - Weekly Skill tool usage audit
#
# Reads .claude/state/skill_usage.log (written by skill_usage_logger.sh)
# and prints:
#   - Top-triggered skills (count + last-triggered timestamp)
#   - 0-trigger skills (candidates for description audit per Anthropic K9)
#
# Log line format (one per Skill invocation):
#   <ISO8601-UTC>|<skill_name>|<session_id>|<cwd>
#
# Usage:
#   bash skill_usage_audit.sh [--days N] [--log-file PATH]
#
# Defaults:
#   --days     7
#   --log-file ${SKILL_USAGE_LOG_DIR:-<repo_root>/.claude/state}/skill_usage.log
#
# Test/override env:
#   SKILLS_DIR  override the skills directory used for 0-trigger comparison
#               (default: <repo_root>/.claude/skills)

set -euo pipefail

# ---------------------------------------------------------------------------
# Repo root resolution (works whether called from repo or anywhere)
# ---------------------------------------------------------------------------
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if REPO_ROOT="$(git -C "$SCRIPT_PATH" rev-parse --show-toplevel 2>/dev/null)"; then
    :
else
    # Fallback: scripts/analysis is two levels under repo root
    REPO_ROOT="$(cd "${SCRIPT_PATH}/../.." && pwd)"
fi

# ---------------------------------------------------------------------------
# Defaults + arg parsing
# ---------------------------------------------------------------------------
DAYS=7
DEFAULT_LOG_DIR="${SKILL_USAGE_LOG_DIR:-${REPO_ROOT}/.claude/state}"
LOG_FILE="${DEFAULT_LOG_DIR}/skill_usage.log"
SKILLS_DIR="${SKILLS_DIR:-${REPO_ROOT}/.claude/skills}"

while [ $# -gt 0 ]; do
    case "$1" in
        --days)
            DAYS="$2"
            shift 2
            ;;
        --log-file)
            LOG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            sed -n '2,20p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Graceful handling: missing or empty log
# ---------------------------------------------------------------------------
if [ ! -f "$LOG_FILE" ]; then
    echo "Log file not found at ${LOG_FILE} -- has any Skill been invoked since the logger was registered?"
    exit 0
fi

if [ ! -s "$LOG_FILE" ]; then
    echo "Log file empty -- no Skill invocations recorded yet"
    exit 0
fi

# ---------------------------------------------------------------------------
# Compute time window
# ---------------------------------------------------------------------------
CUTOFF=$(date -u -d "${DAYS} days ago" +"%Y-%m-%dT%H:%M:%SZ")
PERIOD_START=$(date -u -d "${DAYS} days ago" +"%Y-%m-%dT00:00:00Z")
PERIOD_END=$(date -u +"%Y-%m-%dT23:59:59Z")

# ---------------------------------------------------------------------------
# Enumerate skills directory (each subdir = one skill)
# ---------------------------------------------------------------------------
ALL_SKILLS_FILE="$(mktemp)"
trap 'rm -f "$ALL_SKILLS_FILE" "$TRIGGERED_FILE" "$AGG_FILE" 2>/dev/null || true' EXIT

if [ -d "$SKILLS_DIR" ]; then
    # Only include directories (skip stray files)
    find "$SKILLS_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
        | sort -u > "$ALL_SKILLS_FILE"
fi
TOTAL_SKILLS=$(wc -l < "$ALL_SKILLS_FILE" | tr -d ' ')

# ---------------------------------------------------------------------------
# Aggregate log: filter by cutoff, count per skill, track last-seen + sessions
# ---------------------------------------------------------------------------
AGG_FILE="$(mktemp)"
TRIGGERED_FILE="$(mktemp)"

# awk emits two outputs to separate files via shell redirection:
#  - per-skill aggregate to AGG_FILE: <count>\t<last_ts>\t<skill>
#  - just the unique-triggered skill names to TRIGGERED_FILE
# Plus stdout summary line: "<total>\t<unique>\t<sessions>"
SUMMARY=$(awk -F'|' \
    -v cutoff="$CUTOFF" \
    -v agg_file="$AGG_FILE" \
    -v trig_file="$TRIGGERED_FILE" \
    '
    NF < 2 { next }                       # skip malformed lines
    $1 < cutoff { next }                  # outside time window (ISO8601 lex order)
    {
        skill = $2
        ts    = $1
        sess  = ($3 == "") ? "unknown" : $3
        count[skill]++
        if (last[skill] == "" || ts > last[skill]) last[skill] = ts
        sessions[sess] = 1
        total++
    }
    END {
        unique = 0
        for (s in count) {
            unique++
            printf "%d\t%s\t%s\n", count[s], last[s], s >> agg_file
            print s >> trig_file
        }
        sess_n = 0
        for (s in sessions) sess_n++
        printf "%d\t%d\t%d\n", total+0, unique, sess_n
    }
    ' "$LOG_FILE")

TOTAL=$(echo "$SUMMARY" | awk '{print $1}')
UNIQUE=$(echo "$SUMMARY" | awk '{print $2}')
SESSIONS=$(echo "$SUMMARY" | awk '{print $3}')

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
echo "=== Skill Usage Audit (last ${DAYS} days) ==="
echo "Period: ${PERIOD_START} to ${PERIOD_END}"
echo "Log file: ${LOG_FILE}"
echo "Total invocations: ${TOTAL}"
echo "Unique skills triggered: ${UNIQUE} / ${TOTAL_SKILLS}"
echo "Sessions: ${SESSIONS}"
echo ""

# ---------------------------------------------------------------------------
# Top-triggered table (sorted by count desc, then last-triggered desc)
# ---------------------------------------------------------------------------
echo "--- Top triggered skills ---"
if [ "$TOTAL" -eq 0 ]; then
    echo "(none in the selected window)"
else
    printf "%-32s %-7s %s\n" "Skill" "Count" "Last triggered"
    # sort: -k1nr (count desc), -k2r (timestamp desc)
    sort -t$'\t' -k1,1nr -k2,2r "$AGG_FILE" \
        | awk -F'\t' '{ printf "%-32s %-7s %s\n", $3, $1, $2 }'
fi
echo ""

# ---------------------------------------------------------------------------
# 0-trigger skills (registered but not invoked in window)
# ---------------------------------------------------------------------------
if [ "$TOTAL_SKILLS" -gt 0 ]; then
    # comm requires sorted inputs; both are sort -u'd already
    sort -u "$TRIGGERED_FILE" -o "$TRIGGERED_FILE" 2>/dev/null || true
    ZERO_TMP="$(mktemp)"
    comm -23 "$ALL_SKILLS_FILE" "$TRIGGERED_FILE" > "$ZERO_TMP"
    ZERO_COUNT=$(wc -l < "$ZERO_TMP" | tr -d ' ')
    echo "--- 0-trigger skills (${ZERO_COUNT} total, candidates for description audit) ---"
    if [ "$ZERO_COUNT" -eq 0 ]; then
        echo "(all registered skills triggered at least once)"
    else
        cat "$ZERO_TMP"
    fi
    rm -f "$ZERO_TMP"
else
    echo "--- 0-trigger skills (0 total, candidates for description audit) ---"
    echo "(SKILLS_DIR not found or empty: ${SKILLS_DIR})"
fi

exit 0
