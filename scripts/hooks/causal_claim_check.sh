#!/bin/bash
# causal_claim_check.sh — Advisory lint for unsubstantiated causal claims in validated reports
#
# Pattern: /scientific-rigor §4 DAG enforcement (Pearl Causal Inference)
# Source:  Plan F.3 落地 (2026-05-18)
# Trigger: PostToolUse on Edit|Write to docs/reports/validated/** or docs/experiments/validated/**
# Behavior:
#   - Scan for causal language (導致 / 由於 / 因此 / cause / due to / leads to)
#   - If found AND no DAG/confound/control mention → emit WARNING
#   - exit 0 (advisory, NOT blocking)
#
# Goal: Prevent causal claims without DAG / confound control (L2 Collider Bias 防線)
# Reference: InterSubMod/.claude/skills/scientific-rigor/SKILL.md §4 (DAG + Pearl)
#            InterSubMod/.claude/skills/known-pitfalls/SKILL.md P-01 (L2 Collider Bias)

set -o pipefail

INPUT=$(cat 2>/dev/null || echo "{}")
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // "unknown"' 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

[ "$TOOL_NAME" != "Edit" ] && [ "$TOOL_NAME" != "Write" ] && exit 0
[ -z "$FILE_PATH" ] && exit 0
case "$FILE_PATH" in
  *.md) ;;
  *) exit 0 ;;
esac
case "$FILE_PATH" in
  */docs/reports/validated/*|*/docs/experiments/validated/*|*/docs/experiments/concluded/*) ;;
  *) exit 0 ;;
esac
[ ! -f "$FILE_PATH" ] && exit 0

causal_count=$(grep -ciE "導致|由於|因此|因為.*所以|cause |caused |due to|leads? to|results? in" "$FILE_PATH" 2>/dev/null | tr -d ' ')
dag_count=$(grep -ciE "DAG|confound|collider|mediator|對照組|control group|residualize|residual.*OLS|within.group" "$FILE_PATH" 2>/dev/null | tr -d ' ')
causal_count=${causal_count:-0}
dag_count=${dag_count:-0}

# Threshold: >3 causal claims and <1 DAG/confound mention → warn
if [ "$causal_count" -gt 3 ] && [ "$dag_count" -lt 1 ]; then
  echo "[causal-claim-check] WARN: $FILE_PATH has $causal_count causal-language claims but $dag_count DAG/confound mentions. Consider InterSubMod/.claude/skills/scientific-rigor/SKILL.md §4 (DAG) + InterSubMod/.claude/skills/known-pitfalls/SKILL.md P-01 (collider bias)." >&2
fi

exit 0
