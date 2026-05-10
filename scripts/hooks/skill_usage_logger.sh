#!/usr/bin/env bash
# skill_usage_logger.sh - PreToolUse hook (matcher: Skill)
# Logs every Skill tool invocation to .claude/state/skill_usage.log for
# weekly audit (Anthropic best-practice: identify under-triggered skills).
#
# Log line format:
#   <ISO8601-UTC>|<skill_name>|<session_id>|<cwd>
#
# Behavior:
#   - Reads PreToolUse JSON payload from stdin
#   - Extracts .tool_input.skill (skill name) and .session_id
#   - If skill field missing/null, exit 0 silently (other matchers may fire)
#   - On any unexpected error, exit 0 (never break Claude Code's hook chain)
#
# Override log location for tests via SKILL_USAGE_LOG_DIR env var.

set -uo pipefail

# Resolve log directory: env override > project default
DEFAULT_LOG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/.claude/state"
LOG_DIR="${SKILL_USAGE_LOG_DIR:-$DEFAULT_LOG_DIR}"
LOG_FILE="${LOG_DIR}/skill_usage.log"

# Read stdin payload
INPUT=$(cat)

# Extract skill name; if jq fails or field missing, exit 0 silently
SKILL=$(echo "$INPUT" | jq -r '.tool_input.skill // empty' 2>/dev/null || true)
[ -z "$SKILL" ] && exit 0

# Extract session_id (default "unknown" if missing)
SESSION_ID=$(echo "$INPUT" | jq -r '.session_id // "unknown"' 2>/dev/null || echo "unknown")
[ -z "$SESSION_ID" ] && SESSION_ID="unknown"

# Build log line: timestamp | skill | session | cwd
TS=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
CWD=$(pwd)
LINE="${TS}|${SKILL}|${SESSION_ID}|${CWD}"

# Ensure log dir exists; on failure, exit 0 silently (don't break hook chain)
mkdir -p "$LOG_DIR" 2>/dev/null || exit 0
echo "$LINE" >> "$LOG_FILE" 2>/dev/null || exit 0

exit 0
