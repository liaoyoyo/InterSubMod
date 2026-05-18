#!/bin/bash
# evidence_read_tracker.sh — Records when evidence_ledger.jsonl is read.
#
# Pairs with verify_gate.sh — tracks Read events on evidence_ledger so
# verify_gate can confirm "recent evidence consulted".
# Trigger: PostToolUse on Read
# Tracks:  .claude/.evidence-reads (timestamp per line)

set -o pipefail

INPUT=$(cat 2>/dev/null || echo "{}")
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // "unknown"' 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

[ "$TOOL_NAME" != "Read" ] && exit 0
[ -z "$FILE_PATH" ] && exit 0

# Match evidence_ledger.jsonl or hypothesis_queue.json
case "$FILE_PATH" in
    */evidence_ledger.jsonl|*/hypothesis_queue.json|*/MEMORY.md)
        TRACKER="/big7_disk/liaoyoyo2001/InterSubMod/.claude/.evidence-reads"
        echo "$(date -Iseconds)" >> "$TRACKER"
        # Keep only last 100 entries
        tail -100 "$TRACKER" > "$TRACKER.tmp" 2>/dev/null && mv "$TRACKER.tmp" "$TRACKER" 2>/dev/null
        ;;
esac

exit 0
