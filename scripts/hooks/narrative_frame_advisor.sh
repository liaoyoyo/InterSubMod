#!/usr/bin/env bash
# narrative_frame_advisor.sh — UserPromptSubmit hook advisor
#
# Detect keyword in user prompt → inject advisory recommending /narrative-frame.
# Not enforcing; only nudging. User can ignore or say "不用框架".
#
# Trigger event: UserPromptSubmit
# Behavior: stdout advisory if keyword match; silent otherwise. Always exit 0.

# Read user prompt from stdin (Claude Code hook standard)
PROMPT="${CLAUDE_USER_PROMPT:-$(cat)}"

# Skip empty
[[ -z "$PROMPT" ]] && exit 0

# Skip if prompt too short (< 10 chars likely command not narrative request)
[[ ${#PROMPT} -lt 10 ]] && exit 0

# Skip if user explicitly opted out
if echo "$PROMPT" | grep -qE "不用框架|不要框架|skip framework|no framework|raw answer"; then
  exit 0
fi

# Skip if already in framework context (avoid recursion)
if echo "$PROMPT" | grep -qE "/narrative-frame|narrative-frame N[1-6]|N1 場景識別"; then
  exit 0
fi

# Keyword detection (Chinese + English)
ZH_KEYWORDS="整理|報告|説明|說明|彙報|總結|介紹|講解|解釋|簡報|教|寫|整合|分析|對比|比較|pitch|答辯|交給|呈現"
EN_KEYWORDS="explain|summarize|summarise|report|pitch|present|teach|walk through|integrate|outline|tell me about|breakdown"

MATCH=""
if echo "$PROMPT" | grep -qE "$ZH_KEYWORDS"; then
  MATCH="zh"
elif echo "$PROMPT" | grep -qiE "$EN_KEYWORDS"; then
  MATCH="en"
fi

[[ -z "$MATCH" ]] && exit 0

# Length / complexity heuristic
# Only trigger if prompt suggests output likely ≥200 chars or ≥2 concepts
PROMPT_LEN=${#PROMPT}

# Skip very short prompts (likely simple Q&A)
# Threshold: 12 chars (Chinese-friendly; 1 char ≈ 1 word for CJK)
if [[ $PROMPT_LEN -lt 12 ]]; then
  exit 0
fi

# Emit advisory
cat <<'EOF'

[narrative-frame advisor]
本回覆涉及整理 / 報告 / 説明類內容；建議套用敘述框架以減少理解負擔。

推薦做法（依複雜度）:
- Tier 2（200-500 字）: 回覆首行聲明 framework（如「用 PREP：」）+ 結構化內容
- Tier 3（≥500 字 / 多文件）: 跑完整 /narrative-frame N1-N6（場景識別 → 框架推薦 → 萃取 → 套用 → 自審）

快速查框架對照: InterSubMod/.claude/skills/narrative-frame/references/scenario_to_framework.md
完整 catalog 50+: InterSubMod/.claude/skills/narrative-frame/references/framework_catalog.md

如不需要框架，可直接告知「不用框架」 — 本 advisor 將 skip。

EOF

exit 0
