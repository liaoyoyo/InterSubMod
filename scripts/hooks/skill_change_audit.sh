#!/bin/bash
# skill_change_audit.sh — PostToolUse hook (matcher: Edit|Write)
#
# Purpose: When .claude/skills/<name>/SKILL.md is modified, append an audit entry
#          to docs/postmortems/skill_audit_<YYYYMM>.log. Prevents 6-month drift
#          by accumulating change history for monthly review.
#
# Source:  plan §F.4 P3 Living-update mechanism (scientific-rigor §11 line 550)
# Trigger: PostToolUse on Edit|Write, file path matches .claude/skills/*/SKILL.md
#          or .claude/skills/<name>/(references|scripts|evals.json)

INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

[ -z "$FILE_PATH" ] && exit 0

# Match skill-related paths
if ! echo "$FILE_PATH" | grep -qE '\.claude/skills/[^/]+/(SKILL\.md|evals\.json|references/.*|scripts/.*)$'; then
    exit 0
fi

# Extract skill name (between .claude/skills/ and next /)
SKILL_NAME=$(echo "$FILE_PATH" | sed -E 's|.*/\.claude/skills/([^/]+)/.*|\1|')
FILE_NAME=$(basename "$FILE_PATH")

# Audit log path (monthly rotation)
YYYYMM=$(date +%Y%m)
LOG_DIR="/big7_disk/liaoyoyo2001/InterSubMod/docs/postmortems"
LOG_FILE="${LOG_DIR}/skill_audit_${YYYYMM}.log"

mkdir -p "$LOG_DIR" 2>/dev/null

# Compute file size + line count (best effort, file may not exist if Write-deleted)
SIZE_BYTES=0
LINE_COUNT=0
if [ -f "$FILE_PATH" ]; then
    SIZE_BYTES=$(stat -c %s "$FILE_PATH" 2>/dev/null || echo 0)
    LINE_COUNT=$(wc -l < "$FILE_PATH" 2>/dev/null || echo 0)
fi

# Append audit entry (tab-separated)
TIMESTAMP=$(date -Iseconds)
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // "unknown"' 2>/dev/null)
printf '%s\tskill=%s\tfile=%s\ttool=%s\tbytes=%s\tlines=%s\tpath=%s\n' \
    "$TIMESTAMP" "$SKILL_NAME" "$FILE_NAME" "$TOOL_NAME" "$SIZE_BYTES" "$LINE_COUNT" "$FILE_PATH" \
    >> "$LOG_FILE"

# Notify user (concise)
echo "[Skill Audit] /${SKILL_NAME} modified: ${FILE_NAME} (${LINE_COUNT} lines)"

# Monthly recommendation: if log > 30 entries, suggest review
if [ -f "$LOG_FILE" ]; then
    ENTRY_COUNT=$(wc -l < "$LOG_FILE" 2>/dev/null || echo 0)
    if [ "$ENTRY_COUNT" -gt 30 ]; then
        echo "  ⚠ Audit log has ${ENTRY_COUNT} entries this month — consider monthly skill review"
        echo "  Log: InterSubMod/docs/postmortems/skill_audit_${YYYYMM}.log"
    fi
fi

exit 0
