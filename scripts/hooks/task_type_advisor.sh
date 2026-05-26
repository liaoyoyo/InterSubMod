#!/usr/bin/env bash
# task_type_advisor.sh — UserPromptSubmit hook for task type classification
#
# Purpose (A6, 2026-05-26 落地 — 5/24 incident postmortem):
#   Detect task type keyword in user prompt → emit hookSpecificOutput JSON
#   advising the model to classify task type (A-F) before any further action.
#   Triggers active recall of feedback_observation_scope_default_comprehensive
#   + feedback_task_first_then_doc_then_plan memories.
#
# Trigger event: UserPromptSubmit
# Behavior: emit JSON-wrapped advisory if keyword match; silent otherwise. Always exit 0.
#
# Mirror of narrative_frame_advisor.sh JSON envelope pattern (V1 fix verified).

# Read stdin JSON envelope; fallback to env var for manual testing
PAYLOAD=$(cat 2>/dev/null || true)
if command -v jq >/dev/null 2>&1 && [[ -n "$PAYLOAD" ]]; then
  PROMPT=$(echo "$PAYLOAD" | jq -r '.prompt // empty' 2>/dev/null)
else
  PROMPT=""
fi
[[ -z "$PROMPT" ]] && PROMPT="${CLAUDE_USER_PROMPT:-$PAYLOAD}"

# Skip empty / short
[[ -z "$PROMPT" ]] && exit 0
PROMPT_LEN=${#PROMPT}
if [[ $PROMPT_LEN -lt 15 ]]; then
  exit 0
fi

# Skip if user explicitly opted out
if echo "$PROMPT" | grep -qE "不分類|skip task type|raw answer|不用分類"; then
  exit 0
fi

# Skip if already in task type context (avoid recursion)
if echo "$PROMPT" | grep -qE "task type [A-F]|task_type=[A-F]|§15\.3"; then
  exit 0
fi

# Task type keyword detection — high-recall keywords for B/C/D (full scope mandatory)
# These map to AGENTS.md §15.3 6 categories.
B_KEYWORDS="完整驗證|完整觀察|完整 scope|完整跑|全部驗證|全 chr|全樣本|對應已有觀察|reviewer response|rebuttal|final validation|comprehensive"
C_KEYWORDS="release|deploy|merge to main|production|對外發布|上線"
D_KEYWORDS="handoff|交付|HKU|collaborator|external|對外|給 PI|給教授|給 reviewer|外部團隊"
A_KEYWORDS="pilot|探索|試跑|小規模|快速驗證|sanity check|先看"
E_KEYWORDS="bugfix|hotfix|regress|urgent|排除|急修"
F_KEYWORDS="demo|示意|教學|show|範例|illustration"

# Detect highest-priority match (B/C/D >> A/E/F because full-scope mandates are more critical)
MATCHED_TYPE=""
MATCHED_KW=""
if echo "$PROMPT" | grep -qiE "$D_KEYWORDS"; then
  MATCHED_TYPE="D (Handoff to external)"
  MATCHED_KW=$(echo "$PROMPT" | grep -oiE "$D_KEYWORDS" | head -1)
elif echo "$PROMPT" | grep -qiE "$C_KEYWORDS"; then
  MATCHED_TYPE="C (Production deployment)"
  MATCHED_KW=$(echo "$PROMPT" | grep -oiE "$C_KEYWORDS" | head -1)
elif echo "$PROMPT" | grep -qiE "$B_KEYWORDS"; then
  MATCHED_TYPE="B (Comprehensive validation)"
  MATCHED_KW=$(echo "$PROMPT" | grep -oiE "$B_KEYWORDS" | head -1)
elif echo "$PROMPT" | grep -qiE "$E_KEYWORDS"; then
  MATCHED_TYPE="E (Hotfix / Bugfix)"
  MATCHED_KW=$(echo "$PROMPT" | grep -oiE "$E_KEYWORDS" | head -1)
elif echo "$PROMPT" | grep -qiE "$A_KEYWORDS"; then
  MATCHED_TYPE="A (Exploratory pilot)"
  MATCHED_KW=$(echo "$PROMPT" | grep -oiE "$A_KEYWORDS" | head -1)
elif echo "$PROMPT" | grep -qiE "$F_KEYWORDS"; then
  MATCHED_TYPE="F (Demo / Illustration)"
  MATCHED_KW=$(echo "$PROMPT" | grep -oiE "$F_KEYWORDS" | head -1)
fi

[[ -z "$MATCHED_TYPE" ]] && exit 0

# Log to telemetry
YYYYMM=$(date +%Y%m)
LOG="/big7_disk/liaoyoyo2001/InterSubMod/docs/postmortems/task_type_recall_${YYYYMM}.log"
mkdir -p "$(dirname "$LOG")" 2>/dev/null
printf '%s\t%s\t%s\n' "$(date -Iseconds)" "$MATCHED_TYPE" "$MATCHED_KW" >> "$LOG"

# Emit hookSpecificOutput JSON envelope
cat <<EOF
{"hookSpecificOutput":{"hookEventName":"UserPromptSubmit","additionalContext":"[task-type advisor — 2026-05-26 5/24 incident A6 落地]\n\n偵測到 keyword 「${MATCHED_KW}」→ 預估 task type = **${MATCHED_TYPE}**\n\n**啟動前必做**（per AGENTS.md §15.3 + CLAUDE.md §0）:\n1. 確認 task type — 若預估不對請改判\n2. 對齊 default scope:\n   - A pilot: subset OK\n   - **B / C / D: 全基因組 + 全樣本 mandatory**\n   - E hotfix: 最小可重現\n   - F demo: subset OK + DEMO 標\n3. AskUserQuestion 必含 scope 維度（即使 deliverable 已 ack）\n4. subset 必標 partial flag（/doc-standards 3-locations）\n5. 「見樹也見林」四層（aggregate / canonical / outlier / well-explained）\n\n**主動 recall 必讀**:\n- feedback_observation_scope_default_comprehensive (Step 1-4 + 6 類對照表)\n- feedback_task_first_then_doc_then_plan (Step 1 必含 task type classification)\n\n如預估錯誤或不需 gate，可説「不分類」/「skip task type」 — 本 advisor 將 skip。"}}
EOF

exit 0
