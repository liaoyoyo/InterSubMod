#!/bin/bash
# verify_gate.sh — Default-FAIL evidence gate for validated reports.
#
# Pattern: cwc-long-running-agents Default-FAIL contract + evidence-gated write
# Source:  P4 Top 3 high-ROI fix (T2)
# Trigger: PreToolUse on Edit|Write to docs/reports/validated/** or docs/experiments/concluded/**
# Behavior:
#   - If target file is in validated/concluded path AND user has NOT read evidence_ledger.jsonl in last N seconds
#     → emit WARNING + exit 0 (advisory, NOT blocking) per "advisory by default" principle
#   - If user explicitly bypasses, log to .claude/.evidence-reads-bypass
#   - Otherwise transparent (exit 0 silent)
#
# This is a SOFT gate (advisory only) — not Hard Gate (exit 2 blocking) to avoid
# disrupting normal workflow. Future strengthen to exit 2 if pattern reveals abuse.
#
# Tracking file: .claude/.evidence-reads (one timestamp per line)
# Bypass log:    docs/postmortems/evidence_gate_bypass.log

set -o pipefail

INPUT=$(cat 2>/dev/null || echo "{}")
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // "unknown"' 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

# Only act on Edit/Write
[ "$TOOL_NAME" != "Edit" ] && [ "$TOOL_NAME" != "Write" ] && exit 0
[ -z "$FILE_PATH" ] && exit 0

# Only gate writes to validated reports or concluded experiments
case "$FILE_PATH" in
    */docs/reports/validated/*|*/docs/experiments/concluded/*|*/research/*/concluded/*)
        ;;
    *)
        exit 0  # not a gated path, transparent
        ;;
esac

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
READS_TRACKER="$ROOT/.claude/.evidence-reads"
BYPASS_LOG="$ROOT/docs/postmortems/evidence_gate_bypass.log"
LEDGER="$ROOT/research/autoresearch/evidence_ledger.jsonl"

# How recently must evidence_ledger have been read? (seconds)
RECENCY_WINDOW=3600  # 1 hour

# Check if reads tracker exists and recent
recent_read=0
if [ -f "$READS_TRACKER" ]; then
    last_read=$(tail -1 "$READS_TRACKER" 2>/dev/null)
    if [ -n "$last_read" ]; then
        last_epoch=$(date -d "$last_read" +%s 2>/dev/null || echo 0)
        now_epoch=$(date +%s)
        age=$((now_epoch - last_epoch))
        if [ "$age" -lt "$RECENCY_WINDOW" ]; then
            recent_read=1
        fi
    fi
fi

if [ "$recent_read" -eq 0 ]; then
    # Emit advisory warning
    cat <<EOF
[verify_gate] ⚠ Writing to validated/concluded path: $(basename "$FILE_PATH")
  → No recent evidence_ledger.jsonl read (within ${RECENCY_WINDOW}s).
  → Suggestion: Read InterSubMod/research/autoresearch/evidence_ledger.jsonl first
    to verify your claim is consistent with prior evidence (cwc-long-running-agents pattern).
  → To suppress: Read evidence_ledger.jsonl or touch InterSubMod/.claude/.evidence-reads
  → This is ADVISORY (exit 0). Future may upgrade to Hard Gate if pattern shows abuse.
EOF
    # Log bypass
    mkdir -p "$(dirname "$BYPASS_LOG")" 2>/dev/null
    {
        printf '%s\tfile=%s\ttool=%s\tlast_read=%s\n' \
            "$(date -Iseconds)" "$FILE_PATH" "$TOOL_NAME" "${last_read:-never}"
    } >> "$BYPASS_LOG"
fi

exit 0
