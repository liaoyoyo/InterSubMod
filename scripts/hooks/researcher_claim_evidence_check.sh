#!/bin/bash
# researcher_claim_evidence_check.sh — Detect L3 推測 / unverified claim language and prompt §2 Evidence Tier.
#
# Pattern: 對齊 feedback_researcher_claim_needs_empirical_verification.md +
#          /scientific-rigor §2 Evidence Tier L1-L5
# Source:  P4 Architect E7 / 業界對應 OpenAI structured docs + Walking Labs L09
# Trigger: PostToolUse on Edit|Write to .md files
# Behavior: Scan written content for hedge / speculation patterns; if found, emit
#           advisory reminder to apply L1-L5 evidence ladder annotation.

set -o pipefail

INPUT=$(cat 2>/dev/null || echo "{}")
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // ""' 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // ""' 2>/dev/null)

# Only act on Edit/Write
[ "$TOOL_NAME" != "Edit" ] && [ "$TOOL_NAME" != "Write" ] && exit 0
[ -z "$FILE_PATH" ] && exit 0

# Only check .md files in research / docs/reports / docs/experiments
case "$FILE_PATH" in
    */docs/reports/*|*/docs/experiments/*|*/research/*/00_INDEX.md|*/research/*/REFLECTION.md)
        ;;
    *)
        exit 0  # Not a target path
        ;;
esac

# Skip if file doesn't exist (might be in middle of write)
[ ! -f "$FILE_PATH" ] && exit 0

# Speculation / hedge patterns (case-insensitive)
HEDGE_PATTERNS='(我推測|我認為|probably|likely|maybe|可能是|should be|大概|似乎|seems|appears to|might be|可能會|或許|presumably|conceivably|tentatively)'

# Count matches using grep -o + wc -l (more robust than grep -c)
MATCH_COUNT=$(grep -oiE "$HEDGE_PATTERNS" "$FILE_PATH" 2>/dev/null | wc -l | tr -d ' ')
CLAIM_NO_REF=$(grep -oiE '(顯著|significant|breakthrough|first to|prove[sd]?|定論|確認)' "$FILE_PATH" 2>/dev/null | wc -l | tr -d ' ')
HAS_TIER=$(grep -oE '(L[1-5] |⭐+|Evidence Tier|證據分級|Cohen)' "$FILE_PATH" 2>/dev/null | wc -l | tr -d ' ')

# Fallback to 0 if empty
MATCH_COUNT=${MATCH_COUNT:-0}
CLAIM_NO_REF=${CLAIM_NO_REF:-0}
HAS_TIER=${HAS_TIER:-0}

# Decision: warn if hedges > 2 AND tier annotation missing
if [ "$MATCH_COUNT" -gt 2 ] && [ "$HAS_TIER" -lt 1 ]; then
    cat <<EOF
[researcher_claim_check] ⚠ $(basename "$FILE_PATH"):
  Found ${MATCH_COUNT} hedge/speculation pattern(s) (我推測 / probably / maybe 等)
  but no Evidence Tier annotation (L1-L5 / ⭐⭐⭐ / Cohen ribbon).

  Per /scientific-rigor §2 + feedback_researcher_claim_needs_empirical_verification:
    - L1 完全佐證 (多 dataset + 機制清楚 + 反例排除)
    - L2 部分佐證 (主數據支持)
    - L3 推論 (邏輯合理 + 間接證據) ← 「我推測 / probably」屬此級
    - L4 假設 (待驗證)
    - L5 矛盾 (有反例)

  Suggestion: 對 ${MATCH_COUNT} 個 hedge 句子標註對應 tier，或補實證使升 L1。
EOF
fi

# Also warn on strong claim without evidence ref pattern
if [ "$CLAIM_NO_REF" -gt 3 ] && [ "$HAS_TIER" -lt 1 ]; then
    echo "[researcher_claim_check] ⚠ ${CLAIM_NO_REF} strong claim verb(s) (顯著/breakthrough/prove) without Evidence Tier — apply /scientific-rigor §3 Effect Size ribbon (Cohen's d / CI / sample size)."
fi

exit 0
