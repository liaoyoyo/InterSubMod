#!/bin/bash
# evidence_ledger_sync.sh
# Triggered by PostToolUse when evidence_ledger.jsonl is edited.
# Reminds the user to sync MEMORY.md, experiments/INDEX.md, CURRENT_FOCUS.md.
#
# Detection: reads Claude Code hook stdin JSON for tool_input.file_path.
# Action: reminder-only (non-blocking).

set -u

PROJECT_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
LEDGER="research/autoresearch/evidence_ledger.jsonl"
MEMORY_INDEX="/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md"
EXP_INDEX="$PROJECT_ROOT/docs/experiments/INDEX.md"
FOCUS="$PROJECT_ROOT/docs/CURRENT_FOCUS.md"

# Read hook JSON from stdin; silently bail if jq not available or parse fails.
FILE_PATH="$(jq -r '.tool_input.file_path // empty' 2>/dev/null || true)"

# Only act when evidence_ledger.jsonl was edited
case "$FILE_PATH" in
    */evidence_ledger.jsonl|"$LEDGER")
        ;;
    *)
        exit 0
        ;;
esac

cd "$PROJECT_ROOT" 2>/dev/null || exit 0

# Extract last entry's verdict and hypothesis_id for context
LAST_LINE="$(tail -1 "$LEDGER" 2>/dev/null || true)"
[ -z "$LAST_LINE" ] && exit 0

VERDICT="$(echo "$LAST_LINE" | jq -r '.verdict // "unknown"' 2>/dev/null)"
HID="$(echo "$LAST_LINE" | jq -r '.hypothesis_id // "unknown"' 2>/dev/null)"
TIMESTAMP="$(echo "$LAST_LINE" | jq -r '.timestamp // ""' 2>/dev/null)"

# Determine if sync-worthy (any conclusive verdict)
case "$VERDICT" in
    positive*|negative|NO-GO|no-go|concluded*|positive_pilot|positive_validated)
        SYNC_NEEDED=1
        ;;
    *)
        SYNC_NEEDED=0
        ;;
esac

if [ "$SYNC_NEEDED" -eq 1 ]; then
    echo ""
    echo "[EvidenceLedger Sync] 新結論寫入：${HID} (verdict=${VERDICT}, ts=${TIMESTAMP})"
    echo "  建議同步更新以下索引（依實際適用範圍）："
    echo "    1. MEMORY.md       → 新增 Concluded 條目或更新現有"
    echo "       路徑：$MEMORY_INDEX"
    echo "    2. experiments/INDEX.md → 更新方向狀態"
    echo "       路徑：$EXP_INDEX"
    echo "    3. CURRENT_FOCUS.md → 移除已完成項目、更新阻塞"
    echo "       路徑：$FOCUS"
    echo ""

    # Detect staleness: files older than evidence_ledger by >3 days
    LEDGER_MTIME=$(stat -c %Y "$LEDGER" 2>/dev/null || echo 0)
    for f in "$MEMORY_INDEX" "$EXP_INDEX" "$FOCUS"; do
        if [ -f "$f" ]; then
            F_MTIME=$(stat -c %Y "$f" 2>/dev/null || echo 0)
            DIFF=$(( LEDGER_MTIME - F_MTIME ))
            if [ "$DIFF" -gt 259200 ]; then  # 3 days = 259200s
                DAYS=$(( DIFF / 86400 ))
                echo "  ⚠ $(basename "$f") 距離 ledger 已超過 ${DAYS} 天未更新"
            fi
        fi
    done
fi

exit 0
