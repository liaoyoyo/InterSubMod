#!/bin/bash
# creation_guard.sh — PreToolUse hook (matcher: Write)
#
# Purpose: Pre-flight duplicate / drift check before creating a new SKILL.md or
#          agent definition. Catches: (a) overwriting an existing skill, (b) near-
#          duplicate skill names, (c) reminds to sync the CLAUDE.md §3 skill count.
#          Prevents the index-drift class that caused the 44-vs-45 / phantom grill-me
#          incident (fixed 2026-05-28).
#
# Source:  Community gap G5 — claude-research creation-guard / skill-preflight.
# Trigger: PreToolUse on Write, target path matches a skill or agent definition.
# Exit:    Always 0 (advisory — surfaces warnings, never blocks creation).

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
SKILLS_DIR="${ROOT}/.claude/skills"

INPUT=$(cat 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

[ -z "$FILE_PATH" ] && exit 0

# Only act on new skill SKILL.md or agent .md creations
IS_SKILL=0
IS_AGENT=0
echo "$FILE_PATH" | grep -qE '\.claude/skills/[^/]+/SKILL\.md$' && IS_SKILL=1
echo "$FILE_PATH" | grep -qE '\.claude/agents/[^/]+\.md$' && IS_AGENT=1
[ "$IS_SKILL" -eq 0 ] && [ "$IS_AGENT" -eq 0 ] && exit 0

if [ "$IS_SKILL" -eq 1 ]; then
    NEW_NAME=$(echo "$FILE_PATH" | sed -E 's|.*/\.claude/skills/([^/]+)/SKILL\.md$|\1|')

    # (a) Overwrite check
    if [ -f "$FILE_PATH" ]; then
        echo "[creation-guard] ⚠ /${NEW_NAME} 已存在 — 這是 OVERWRITE 既有 SKILL.md，非新建。確認是否有意覆寫。"
    else
        echo "[creation-guard] 新建 skill: /${NEW_NAME}"
    fi

    # (b) Near-duplicate name check (shared token with existing skill)
    if [ -d "$SKILLS_DIR" ]; then
        for existing in "$SKILLS_DIR"/*/; do
            ex=$(basename "$existing")
            [ "$ex" = "$NEW_NAME" ] && continue
            # split on '-' and look for a shared meaningful token (len>=4)
            for tok in $(echo "$NEW_NAME" | tr '-' ' '); do
                [ ${#tok} -lt 4 ] && continue
                if echo "$ex" | grep -q "$tok"; then
                    echo "[creation-guard] 🔍 名稱近似既有 /${ex}（共享 token '${tok}'）— 確認非重複功能。"
                    break
                fi
            done
        done
    fi

    # (c) Count-sync reminder
    CUR_COUNT=$(find "$SKILLS_DIR" -name SKILL.md 2>/dev/null | wc -l | tr -d ' ')
    AFTER=$((CUR_COUNT + (1 - $([ -f "$FILE_PATH" ] && echo 1 || echo 0))))
    echo "[creation-guard] 目前 ${CUR_COUNT} 個 SKILL.md → 新建後 ${AFTER}。記得同步 CLAUDE.md §3 計數與分類。"
fi

if [ "$IS_AGENT" -eq 1 ]; then
    NEW_NAME=$(basename "$FILE_PATH" .md)
    if [ -f "$FILE_PATH" ]; then
        echo "[creation-guard] ⚠ agent ${NEW_NAME} 已存在 — OVERWRITE 確認。"
    else
        echo "[creation-guard] 新建 agent: ${NEW_NAME} — 確認與既有 agent 無職責重疊。"
    fi
fi

exit 0
