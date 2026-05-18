#!/bin/bash
# evidence_level_lint.sh — Advisory lint for evidence tier ribbon in validated reports
#
# Pattern: /scientific-rigor §2 evidence tier enforcement
# Source:  Plan F.3 落地 (2026-05-18)
# Trigger: PostToolUse on Edit|Write to docs/reports/validated/** or docs/experiments/validated/**
# Behavior:
#   - Scan written .md for AUC/F1 number patterns (e.g., F1=0.X, AUC=0.X, +0.0XX)
#   - If number found AND no L1-L5 / ⭐ / Cohen ribbon nearby → emit WARNING
#   - exit 0 (advisory, NOT blocking) — non-disruptive
#
# Goal: Prevent reports from publishing claims without evidence tier (Replication Crisis 教訓)
# Reference: InterSubMod/.claude/skills/scientific-rigor/SKILL.md §2 (L1-L5 ladder)

set -o pipefail

INPUT=$(cat 2>/dev/null || echo "{}")
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // "unknown"' 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

# Only act on Edit/Write to .md
[ "$TOOL_NAME" != "Edit" ] && [ "$TOOL_NAME" != "Write" ] && exit 0
[ -z "$FILE_PATH" ] && exit 0
case "$FILE_PATH" in
  *.md) ;;
  *) exit 0 ;;
esac

# Only validated reports + concluded experiments + postmortems
case "$FILE_PATH" in
  */docs/reports/validated/*|*/docs/experiments/validated/*|*/docs/experiments/concluded/*|*/docs/postmortems/*) ;;
  *) exit 0 ;;
esac

[ ! -f "$FILE_PATH" ] && exit 0

# Count claim patterns vs evidence tier ribbon
claim_count=$(grep -ciE "(AUC|F1|delta|recall|precision)[ =:]+[0-9]\.[0-9]+|[+-]0\.0[0-9]+" "$FILE_PATH" 2>/dev/null | tr -d ' ')
tier_count=$(grep -ciE "L[1-5][ \t]|⭐|Cohen|95% CI|95% 信賴" "$FILE_PATH" 2>/dev/null | tr -d ' ')
claim_count=${claim_count:-0}
tier_count=${tier_count:-0}

# Threshold: >5 claims and <2 tier ribbons → warn
if [ "$claim_count" -gt 5 ] && [ "$tier_count" -lt 2 ]; then
  echo "[evidence-level-lint] WARN: $FILE_PATH has $claim_count metric claims but only $tier_count evidence tier ribbon (L1-L5 / Cohen / CI). Consider InterSubMod/.claude/skills/scientific-rigor/SKILL.md §2 + §3." >&2
fi

exit 0
