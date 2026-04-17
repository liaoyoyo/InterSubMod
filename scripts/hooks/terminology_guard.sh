#!/bin/bash
# terminology_guard.sh
# Triggered by PostToolUse when editing research/autoresearch/*.json(l) or hypothesis_queue.json.
# Warns (non-blocking) if non-canonical verdict/track/phase values are used.

set -u

PROJECT_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
TERM_FILE="$PROJECT_ROOT/docs/standards/research_terminology.json"

FILE_PATH="$(jq -r '.tool_input.file_path // empty' 2>/dev/null || true)"

# Scope: only evidence_ledger.jsonl, hypothesis_queue.json, or cycles/*.json
case "$FILE_PATH" in
    */evidence_ledger.jsonl|*/hypothesis_queue.json|*/cycles/*/cycle.json)
        ;;
    *)
        exit 0
        ;;
esac

[ ! -f "$TERM_FILE" ] && exit 0
[ ! -f "$FILE_PATH" ] && exit 0

# Build canonical value lists from terminology JSON
CANONICAL_VERDICTS=$(jq -r '.verdict.canonical_values[]' "$TERM_FILE" 2>/dev/null | tr '\n' '|' | sed 's/|$//')
CANONICAL_TRACKS=$(jq -r '.pipeline_track.canonical_values[]' "$TERM_FILE" 2>/dev/null | tr '\n' '|' | sed 's/|$//')

# Extract actual values from the edited file (best-effort; supports JSONL and JSON)
if [[ "$FILE_PATH" == *.jsonl ]]; then
    FOUND_VERDICTS=$(jq -r '.verdict // empty' "$FILE_PATH" 2>/dev/null | sort -u)
    FOUND_TRACKS=$(jq -r '.pipeline_track // empty' "$FILE_PATH" 2>/dev/null | sort -u)
else
    FOUND_VERDICTS=$(jq -r '.. | .verdict? // empty' "$FILE_PATH" 2>/dev/null | sort -u)
    FOUND_TRACKS=$(jq -r '.. | .pipeline_track? // empty' "$FILE_PATH" 2>/dev/null | sort -u)
fi

VIOLATIONS=0
WARNINGS=""

for v in $FOUND_VERDICTS; do
    [ -z "$v" ] && continue
    if ! echo "$v" | grep -qE "^(${CANONICAL_VERDICTS})\$"; then
        WARNINGS+="  verdict: '$v' 非 canonical（允許值：$CANONICAL_VERDICTS）\n"
        VIOLATIONS=$((VIOLATIONS + 1))
    fi
done

for t in $FOUND_TRACKS; do
    [ -z "$t" ] && continue
    if ! echo "$t" | grep -qE "^(${CANONICAL_TRACKS})\$"; then
        WARNINGS+="  pipeline_track: '$t' 非 canonical（允許值：$CANONICAL_TRACKS）\n"
        VIOLATIONS=$((VIOLATIONS + 1))
    fi
done

if [ "$VIOLATIONS" -gt 0 ]; then
    echo ""
    echo "[Terminology Guard] $FILE_PATH 偵測到 $VIOLATIONS 處非 canonical 命名："
    echo -e "$WARNINGS"
    echo "  規範：docs/standards/research_terminology.json"
    echo "  （警告不阻擋；視需要更正）"
    echo ""
fi

exit 0
