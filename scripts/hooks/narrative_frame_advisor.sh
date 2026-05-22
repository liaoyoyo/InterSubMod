#!/usr/bin/env bash
# narrative_frame_advisor.sh — UserPromptSubmit hook advisor
#
# Detect keyword in user prompt → emit hookSpecificOutput JSON advising
# the model to apply /narrative-frame. Not enforcing; only nudging.
# User can opt out via "不用框架" / "skip framework".
#
# Trigger event: UserPromptSubmit
# Behavior: emit JSON-wrapped advisory if keyword match; silent otherwise. Always exit 0.
#
# V1 review fix (2026-05-22):
#   F-I2-01: parse Claude Code stdin JSON envelope + emit hookSpecificOutput JSON
#   F-I2-02: tighten keyword 教 / 寫 to multi-char tokens (avoid false positive on 教授 / 寫程式)

# Read stdin JSON envelope; fallback to env var for manual testing
PAYLOAD=$(cat 2>/dev/null || true)
if command -v jq >/dev/null 2>&1 && [[ -n "$PAYLOAD" ]]; then
  PROMPT=$(echo "$PAYLOAD" | jq -r '.prompt // empty' 2>/dev/null)
else
  PROMPT=""
fi
# Fallback chain: jq-parsed prompt → env var → raw payload (for manual test)
[[ -z "$PROMPT" ]] && PROMPT="${CLAUDE_USER_PROMPT:-$PAYLOAD}"

# Skip empty
[[ -z "$PROMPT" ]] && exit 0

# Skip if user explicitly opted out
if echo "$PROMPT" | grep -qE "不用框架|不要框架|skip framework|no framework|raw answer"; then
  exit 0
fi

# Skip if already in framework context (avoid recursion)
if echo "$PROMPT" | grep -qE "/narrative-frame|narrative-frame N[1-6]|N1 場景識別"; then
  exit 0
fi

# Keyword detection — multi-char tokens only (F-I2-02 fix)
# Chinese: avoid single-char regex like 教/寫 that match 教授/寫程式 (false positive)
ZH_KEYWORDS="整理|報告|説明|說明|彙報|總結|介紹|講解|解釋|簡報|整合|分析|對比|比較|pitch|答辯|呈現|教我|教學|教練|教過|教你|寫報告|寫週報|寫文章|寫摘要|寫結論|寫論文"
EN_KEYWORDS="explain|summarize|summarise|report|pitch|present|teach|walk through|integrate|outline|tell me about|breakdown"

MATCH=""
if echo "$PROMPT" | grep -qE "$ZH_KEYWORDS"; then
  MATCH="zh"
elif echo "$PROMPT" | grep -qiE "$EN_KEYWORDS"; then
  MATCH="en"
fi

[[ -z "$MATCH" ]] && exit 0

# Length / complexity heuristic
# Skip very short prompts (likely simple Q&A or single command)
# Threshold: 12 chars (Chinese-friendly; 1 char ≈ 1 word for CJK)
PROMPT_LEN=${#PROMPT}
if [[ $PROMPT_LEN -lt 12 ]]; then
  exit 0
fi

# Emit hookSpecificOutput JSON envelope (Claude Code UserPromptSubmit standard)
cat <<'EOF'
{"hookSpecificOutput":{"hookEventName":"UserPromptSubmit","additionalContext":"[narrative-frame advisor]\n本回覆涉及整理 / 報告 / 説明類內容；建議套用敘述框架以減少理解負擔。\n\n推薦做法（依複雜度 Tier）:\n- Tier 1 簡單問答: skip framework — 直接答\n- Tier 2 中度整理 (200-500 字 / 跨 2-3 概念): 回覆首行聲明 framework（如「用 PREP：」/「用 SCQA：」）\n- Tier 3 複雜回覆 (≥500 字 / 多文件): 跑完整 /narrative-frame N1-N6 (場景識別 → 框架推薦 → 萃取 → 套用 → 自審)\n- Tier 4 重大 paper-scope / NO-GO: Tier 3 + 對齊 /scientific-rigor §2-§7 (每 claim 標 L1-L5)\n\n快速查框架對照: InterSubMod/.claude/skills/narrative-frame/references/scenario_to_framework.md\n完整 catalog 50+: InterSubMod/.claude/skills/narrative-frame/references/framework_catalog.md\n\n如不需要框架，可直接告知「不用框架」 — 本 advisor 將 skip。"}}
EOF

exit 0
