#!/bin/bash
# skill_registry_sync.sh — PostToolUse hook (matcher: Edit|Write)
#
# Purpose: When an index doc (CLAUDE.md / skills/README.md) is edited, recompute
#          the ACTUAL on-disk skill + agent counts and compare against the counts
#          the doc claims. Warns on mismatch. This is the "delete-direction" guard
#          that creation_guard.sh (Write-only, additions) does not cover — together
#          they close the bidirectional drift gap that caused the 43/44/45/46 +
#          phantom-grill-me incidents.
#
# Source:  Community gap G5 (registry-sync) + 2026-05-30 audit recommendation.
# Trigger: PostToolUse on Edit|Write, file path is CLAUDE.md or skills/README.md.
# Exit:    Always 0 (advisory).

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
SKILLS_DIR="${ROOT}/.claude/skills"
AGENTS_DIR="${ROOT}/.claude/agents"

INPUT=$(cat 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)
[ -z "$FILE_PATH" ] && exit 0

# Only act on the two index docs
IS_CLAUDE=0
IS_README=0
echo "$FILE_PATH" | grep -qE '/\.claude/CLAUDE\.md$' && IS_CLAUDE=1
echo "$FILE_PATH" | grep -qE '/\.claude/skills/README\.md$' && IS_README=1
[ "$IS_CLAUDE" -eq 0 ] && [ "$IS_README" -eq 0 ] && exit 0

# Ground truth from disk
ACTUAL_SKILLS=$(find "$SKILLS_DIR" -name SKILL.md 2>/dev/null | wc -l | tr -d ' ')
ACTUAL_AGENTS=$(ls "$AGENTS_DIR"/*.md 2>/dev/null | wc -l | tr -d ' ')

WARN=0

# Compare the GRAND TOTAL (max claimed number) against disk truth.
# Rationale: category subtotals (e.g. "Meta — 8 skills") are always <= the grand
# total, so taking the MAX claimed number avoids false-positives on subtotals.
if [ -f "$FILE_PATH" ]; then
    MAX_SKILL=$(grep -oE '[0-9]+ 個 (SKILL\.md|skill)|[0-9]+ skills' "$FILE_PATH" 2>/dev/null | grep -oE '[0-9]+' | sort -n | tail -1)
    if [ -n "$MAX_SKILL" ] && [ "$MAX_SKILL" != "$ACTUAL_SKILLS" ]; then
        echo "[registry-sync] ⚠ $(basename "$FILE_PATH") 宣稱最大 skill 總數 '${MAX_SKILL}'，但磁碟實際 ${ACTUAL_SKILLS} 個 SKILL.md — 請校正總數（分類小計可 <總數）"
        WARN=1
    fi
    MAX_AGENT=$(grep -oE '[0-9]+ 個 (project )?(sub-?)?agent' "$FILE_PATH" 2>/dev/null | grep -oE '[0-9]+' | sort -n | tail -1)
    if [ -n "$MAX_AGENT" ] && [ "$MAX_AGENT" != "$ACTUAL_AGENTS" ]; then
        echo "[registry-sync] ⚠ $(basename "$FILE_PATH") 宣稱最大 agent 數 '${MAX_AGENT}'，但磁碟實際 ${ACTUAL_AGENTS} 個 .claude/agents/*.md — 請校正"
        WARN=1
    fi
fi

if [ "$WARN" -eq 0 ]; then
    echo "[registry-sync] ✓ $(basename "$FILE_PATH") 計數與磁碟一致（${ACTUAL_SKILLS} skills / ${ACTUAL_AGENTS} agents）"
fi

exit 0
