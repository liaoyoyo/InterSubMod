# 風格 Hooks 遷移清單（14 個 baseline）

> 來源：`/big7_disk/liaoyoyo2001/InterSubMod/scripts/hooks/`
> 目標：`<new-project>/scripts/hooks/`
> **重要**：所有 hook 含 hardcoded 絕對路徑 `/big7_disk/liaoyoyo2001/InterSubMod`，拷貝後**必跑** `portable/hook_path_rewrite.sh` 替換為 `${CLAUDE_PROJECT_DIR}` 或 `$(pwd)`。

---

## 拷貝指令

```bash
SRC=/big7_disk/liaoyoyo2001/InterSubMod/scripts/hooks
DST=<new-project>/scripts/hooks

PORTABLE_HOOKS=(
  # SessionStart
  session_start_inject_focus.sh

  # UserPromptSubmit
  md_path_format_rule.sh
  narrative_frame_advisor.sh
  task_type_advisor.sh

  # PreToolUse
  verify_gate.sh                  # Edit|Write Default-FAIL evidence gate
  pre_commit_compile_check.sh     # Bash C++ commit 編譯（非 C++ 專案可不啟用）
  no_binary_commit.sh             # Bash binary 阻擋
  skill_usage_logger.sh           # Skill 使用統計

  # PostToolUse
  external_input_sanitizer.sh     # WebFetch injection 偵測
  evidence_read_tracker.sh        # Read 追蹤（搭配 verify_gate）
  memory_recall_logger.sh         # Read 記憶引用量化
  skill_change_audit.sh           # Edit|Write skill 變動 log

  # SubagentStop
  subagent_completion_logger.sh

  # Manual / 觸發後執行
  cache_telemetry.sh              # prompt cache hit 統計
  allow_list_audit.sh             # permissions audit
)

mkdir -p "$DST"
for h in "${PORTABLE_HOOKS[@]}"; do
  if [ -f "$SRC/$h" ]; then
    cp "$SRC/$h" "$DST/"
    chmod +x "$DST/$h"
    echo "  cp $h"
  else
    echo "  MISSING: $h"
  fi
done

# 路徑改寫（必跑）
bash <new-project>/docs/references/migration/portable/hook_path_rewrite.sh "$DST"
```

---

## settings.local.json hooks 配置（拷貝至新專案）

新專案 `.claude/settings.local.json` 的 `hooks` 區段：

```json
{
  "hooks": {
    "SessionStart": [
      {"hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/session_start_inject_focus.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]}
    ],
    "UserPromptSubmit": [
      {"hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/md_path_format_rule.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]},
      {"hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/narrative_frame_advisor.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]},
      {"hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/task_type_advisor.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]}
    ],
    "PreToolUse": [
      {"matcher": "Bash", "hooks": [{"type": "command", "command": "jq -r '.tool_input.command // empty' | grep -qE '(rm -rf|git reset --hard|git push --force)' && echo '[PreHook] 危險操作警告！' || true"}]},
      {"matcher": "Bash", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/pre_commit_compile_check.sh 2>/dev/null || exit 0"}]},
      {"matcher": "Bash", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/no_binary_commit.sh"}]},
      {"matcher": "Skill", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/skill_usage_logger.sh"}]},
      {"matcher": "Edit|Write", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/verify_gate.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]}
    ],
    "PostToolUse": [
      {"matcher": "Edit|Write", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/skill_change_audit.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]},
      {"matcher": "Read", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/evidence_read_tracker.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]},
      {"matcher": "Read", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/memory_recall_logger.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]},
      {"matcher": "WebFetch", "hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/external_input_sanitizer.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]}
    ],
    "SubagentStop": [
      {"hooks": [{"type": "command", "command": "echo '[SubagentStop] 子代理任務完成，請檢查輸出結果並確認是否符合預期' || true"}]},
      {"hooks": [{"type": "command", "command": "bash ${CLAUDE_PROJECT_DIR}/scripts/hooks/subagent_completion_logger.sh 2>>${CLAUDE_PROJECT_DIR}/docs/postmortems/hook_failures.log || true"}]}
    ],
    "Stop": [
      {"hooks": [{"type": "command", "command": "printf '[Hook] 會話即將結束，請確認：\\n1. 是否已撰寫執行報告\\n2. 是否有未完成任務需要記錄\\n' || true"}]}
    ]
  }
}
```

---

## Hook 用途速查

| Hook | Event | 一行用途 | 必要性 |
|------|-------|---------|--------|
| `session_start_inject_focus.sh` | SessionStart | 注入 CURRENT_FOCUS.md 等價檔到 context | ⭐⭐⭐ |
| `md_path_format_rule.sh` | UserPromptSubmit | 注入 .md 路徑前綴規則提醒 | ⭐⭐ |
| `narrative_frame_advisor.sh` | UserPromptSubmit | 偵測「整理 / 報告 / 説明」keyword → 推薦 framework | ⭐⭐⭐ |
| `task_type_advisor.sh` | UserPromptSubmit | 偵測 6 類 task type keyword + 注入 advisory | ⭐⭐⭐ |
| `verify_gate.sh` | PreToolUse Edit\|Write | Default-FAIL evidence gate（未 Read 過不准 Edit） | ⭐⭐⭐ |
| `pre_commit_compile_check.sh` | PreToolUse Bash | C++ commit 前強制編譯（非 C++ 專案可不啟用） | ⭐⭐ (條件) |
| `no_binary_commit.sh` | PreToolUse Bash | binary commit 阻擋 | ⭐⭐⭐ |
| `external_input_sanitizer.sh` | PostToolUse WebFetch | WebFetch 結果 injection 偵測 | ⭐⭐⭐ |
| `evidence_read_tracker.sh` | PostToolUse Read | Read 追蹤（搭配 verify_gate） | ⭐⭐ |
| `memory_recall_logger.sh` | PostToolUse Read | 記憶引用率量化 | ⭐⭐ |
| `skill_change_audit.sh` | PostToolUse Edit\|Write | skill 檔案變動月度 log | ⭐⭐ |
| `skill_usage_logger.sh` | PreToolUse Skill | skill 使用統計 | ⭐⭐ |
| `subagent_completion_logger.sh` | SubagentStop | subagent cost / cache / artifact 紀錄 | ⭐⭐ |
| `cache_telemetry.sh` | manual | prompt cache hit 統計 | ⭐ |
| `allow_list_audit.sh` | manual | permissions allow-list audit | ⭐ |

---

## 不遷移的 hooks（16 個 ISM 綁定）

`knowledge_check.sh` `kb_freshness_warn.sh` `kb_schema_check.sh` `kb_sot_guard.sh` `research_direction_guard.sh` `pipeline_block_check.sh` `evidence_ledger_sync.sh` `pre_tier_upgrade_check.sh` `cpp_edit_guard.sh` `compile_success_clear.sh` `post_cpp_commit_invalidate.sh` `terminology_guard.sh` `standalone_trigger.sh` `trigger_routing.sh` `researcher_claim_evidence_check.sh` `evidence_level_lint.sh` `causal_claim_check.sh` `compact_test.sh` `md_link_check.sh`

→ 綁定 ISM Knowledge MCP / evidence_ledger / 7-Phase cycle / C++ 結構；新專案不適用。

但其中以下兩個**概念可重用**（新專案若有對應 artifact 可仿製）：
- `evidence_ledger_sync.sh` — 若新任務有「事實 ledger」追蹤需求
- `researcher_claim_evidence_check.sh` `evidence_level_lint.sh` `causal_claim_check.sh` — 若新任務有「claim 證據分級」需求
