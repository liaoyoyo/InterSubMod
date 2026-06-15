#!/bin/bash
# subagent_completion_logger.sh — SubagentStop hook 記錄 cost / cache / status.
#
# Purpose: 量化 sub-agent ROI（P3 audit Subagent-Logger）
#          - 哪些 agent 真有用（呼叫頻率 / 成功率）
#          - cache hit rate vs 主 agent
#          - cost-per-task 分析
# Trigger: SubagentStop event
# Output: docs/postmortems/subagent_completion_<YYYYMM>.log (monthly rotation)

set -o pipefail

LOG_DIR="/big7_disk/liaoyoyo2001/InterSubMod/docs/postmortems"
YYYYMM=$(date +%Y%m)
LOG_FILE="${LOG_DIR}/subagent_completion_${YYYYMM}.log"
mkdir -p "$LOG_DIR" 2>/dev/null

# SubagentStop hook stdin format (Anthropic spec):
#   { "session_id": "...", "transcript_path": "...", "stop_hook_active": ..., ... }
# May also include subagent_type, total_cost_usd, etc. depending on CLI version
INPUT=$(cat)

TIMESTAMP=$(date -Iseconds)
SUBAGENT_TYPE=$(echo "$INPUT" | jq -r '.subagent_type // .agent_type // "unknown"' 2>/dev/null)
SESSION_ID=$(echo "$INPUT" | jq -r '.session_id // "—"' 2>/dev/null)
TRANSCRIPT_PATH=$(echo "$INPUT" | jq -r '.transcript_path // "—"' 2>/dev/null)
COST=$(echo "$INPUT" | jq -r '.total_cost_usd // .cost // 0' 2>/dev/null)

# Try to extract cache stats from transcript if path provided
CACHE_READ=0
CACHE_CREATE=0
INPUT_TOKENS=0
OUTPUT_TOKENS=0
TURN_COUNT=0
if [ -f "$TRANSCRIPT_PATH" ]; then
    STATS=$(python3 -c "
import json, sys
cr = cc = it = ot = tc = 0
try:
    with open('$TRANSCRIPT_PATH') as f:
        for line in f:
            try:
                d = json.loads(line)
                u = d.get('message', {}).get('usage', {})
                if u:
                    tc += 1
                    cr += u.get('cache_read_input_tokens', 0)
                    cc += u.get('cache_creation_input_tokens', 0)
                    it += u.get('input_tokens', 0)
                    ot += u.get('output_tokens', 0)
            except: pass
except: pass
print(f'{cr} {cc} {it} {ot} {tc}')
" 2>/dev/null)
    if [ -n "$STATS" ]; then
        read -r CACHE_READ CACHE_CREATE INPUT_TOKENS OUTPUT_TOKENS TURN_COUNT <<< "$STATS"
    fi
fi

# Cache hit rate
HIT_RATE="—"
TOTAL_IN=$((CACHE_READ + CACHE_CREATE + INPUT_TOKENS))
if [ "$TOTAL_IN" -gt 0 ]; then
    HIT_RATE=$(python3 -c "print(f'{$CACHE_READ * 100 / $TOTAL_IN:.1f}')")
fi

# Append entry (tab-separated)
printf '%s\tagent=%s\tsession=%s\tcost=%s\tturns=%s\tcache_hit=%s%%\tin=%s\tcr=%s\tcc=%s\tout=%s\n' \
    "$TIMESTAMP" "$SUBAGENT_TYPE" "${SESSION_ID:0:12}" "$COST" "$TURN_COUNT" "$HIT_RATE" \
    "$INPUT_TOKENS" "$CACHE_READ" "$CACHE_CREATE" "$OUTPUT_TOKENS" \
    >> "$LOG_FILE"

# Notify
echo "[SubAgent Logger] $SUBAGENT_TYPE done — cost=\$$COST turns=$TURN_COUNT cache_hit=${HIT_RATE}%"

# Output-token advisory (2026-06-15 audit D3-3 — wire the >3K flag CLAUDE.md §8 promised but never echoed).
# §8 return-contract soft target ~1-2K; >3K = ensure full detail was LANDED to an OUTPUT_DIR file
# (concise-emit), not just returned inline where it risks the output-token cap (friction #1).
if [ "${OUTPUT_TOKENS:-0}" -gt 3000 ]; then
    echo "  ⚠ out=${OUTPUT_TOKENS} tokens >3K — confirm this agent's detail was written to a file & only a {status,metrics,path} summary returned (§8 concise-emit)"
fi

# Monthly recommendation
ENTRY_COUNT=$(wc -l < "$LOG_FILE" 2>/dev/null || echo 0)
if [ "$ENTRY_COUNT" -gt 50 ]; then
    echo "  ℹ ${ENTRY_COUNT} subagent invocations this month — review effectiveness in:"
    echo "    InterSubMod/docs/postmortems/subagent_completion_${YYYYMM}.log"
fi

exit 0
