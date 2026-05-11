---
id: ism-kb-00-governance-hooks-and-automation
name: "Hooks 與自動化設置"
description: "`.claude/settings.local.json` 內 10+ hooks 配置說明：觸發時機、動作、阻擋級別；包含 C++ commit 編譯檢查、危險命令警告、待編譯標記、evidence ledger sync 等。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "hooks against .claude/settings.local.json HEAD"
related_ids:
  - ism-kb-00-governance-index
  - ism-kb-06-workflows-cpp-change-pdd
  - ism-kb-06-workflows-build-and-test
tags: [governance, hooks, automation, commit-block, pre-tool-use]
canonical_paths: [00_governance/08_hooks_and_automation.md]
alias_paths: []
---

# Hooks 與自動化設置

- 一句結論：10+ hooks 自動執行於 UserPromptSubmit/PreToolUse/PostToolUse/SubagentStop/Stop；C++ 未編譯 commit 會被 Hard Gate 阻擋
- 適用對象：理解「為何 git commit 被擋」、新增 hook 前
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  jq '.hooks | keys' /big7_disk/liaoyoyo2001/InterSubMod/.claude/settings.local.json
  ```

---

## 全 Hooks 對照表

| Hook 類型 | Matcher | 腳本 / 動作 | 阻擋級別 |
|-----------|---------|-------------|----------|
| `UserPromptSubmit` | — | `scripts/hooks/knowledge_check.sh`（知識庫相關主題提醒）| 提醒 |
| `UserPromptSubmit` | — | `scripts/hooks/research_direction_guard.sh`（**已關閉**研究方向警告）| 警告 |
| `PreToolUse` | Bash | **危險命令警告**（`rm -rf` / `git reset --hard` / `git push --force`）| 提醒 |
| `PreToolUse` | Bash | `scripts/hooks/pre_commit_compile_check.sh`（**C++ 未編譯阻擋 commit**）| 🔴 **阻擋 `exit 2`** |
| `PostToolUse` | Edit/Write | `scripts/hooks/cpp_edit_guard.sh`（編輯 `.cpp`/`.hpp`/`.h` 建立待編譯標記）| 標記 |
| `PostToolUse` | Edit/Write | `scripts/hooks/trigger_routing.sh`（依檔案類型建議對應流程）| 提醒 |
| `PostToolUse` | Edit/Write | `scripts/hooks/evidence_ledger_sync.sh`（evidence_ledger.jsonl 同步檢查）| 提醒 |
| `PostToolUse` | Edit/Write | `scripts/hooks/terminology_guard.sh`（術語一致性檢查）| 提醒 |
| `PostToolUse` | Edit/Write | `scripts/hooks/md_link_check.sh`（Markdown link 有效性）| 提醒 |
| `PostToolUse` | Bash | `scripts/hooks/compile_success_clear.sh`（make 成功 → 清除待編譯標記）| 清除 |
| `PostToolUse` | Bash | `git commit` 後提醒確認測試與文檔 | 提醒 |
| `SubagentStop` | — | 子代理完成提醒檢查結果 | 提醒 |
| `Stop` | — | 會話結束提醒撰寫執行報告 | 提醒 |
| `UserPromptSubmit` | — | **`scripts/hooks/kb_freshness_warn.sh`**（10_research_status/ 快照 >14 天 → 提醒） | 提醒 |
| `PreToolUse` | Bash | **`scripts/hooks/kb_schema_check.sh`**（git commit 前 knowledge/ 變動跑 3 驗證腳本，失敗阻擋） | 🔴 **阻擋 `exit 2`** |
| `PostToolUse` | Edit/Write | **`scripts/hooks/kb_sot_guard.sh`**（偵測 ΔF1 SoT 受保護數字寫入非 SoT 文件 → 提醒連結） | 提醒 |

---

## 🔴 Hard Gate：KB schema 違規阻擋 commit（v0.5+）

**腳本**：`scripts/hooks/kb_schema_check.sh`
**觸發**：`git commit` 前（僅當 staged 檔案含 `knowledge/**/*.md` 時）
**動作**：依序跑 3 驗證腳本：
1. `validate_frontmatter.py` — frontmatter schema 合規
2. `check_related_ids_symmetry.py` — related_ids 雙向對稱
3. `check_canonical_paths.py` — canonical_paths 存在

**阻擋條件**：任一失敗 → `exit 2`
**繞過**：`git commit --no-verify`（不建議；僅在 v0.5 容錯期允許遷移舊 `source_type`）

---

## ⏰ KB Freshness 警示（v0.5+）

**腳本**：`scripts/hooks/kb_freshness_warn.sh`
**觸發**：`UserPromptSubmit`
**監控檔案**：
- `10_research_status/01_current_focus_snapshot.md`
- `10_research_status/02_active_hypotheses.md`
- `10_research_status/04_blockers_and_risks.md`
- `10_research_status/05_next_milestones.md`
- `09_conclusions/05_hypothesis_queue_snapshot.md`

**閾值**：`last_verified` > 14 天 → 印提醒
**阻擋**：不阻擋

**解除方式**：更新對應文件 `last_verified` 為今日 + `verified_scope` 描述當次驗證

---

## 🎯 SoT 破窗偵測（v0.5+）

**腳本**：`scripts/hooks/kb_sot_guard.sh`
**觸發**：`PostToolUse` / `Edit|Write`
**偵測模式**：若 `tool_input.new_string` 或 `tool_input.content` 含以下 ΔF1 SoT 受保護數字：
- `+0.0112` / `0.9762`（CL-002 paired_full）
- `-0.0206` / `0.9650`（CL-003 TO）

**豁免**：
- `03_pipelines/05_f1_baseline_canonical.md`（SoT 本尊）
- `CHANGELOG.md`（版本紀錄）

**動作**：印提醒「請連結 SoT 文件而非複製數字」；不阻擋
**目的**：避免 v0.1 時代「ΔF1 22 處重複」的 SoT 破窗問題重演

---

## 🔴 Hard Gate：C++ 未編譯阻擋 commit

**腳本**：`scripts/hooks/pre_commit_compile_check.sh`
**觸發**：`git commit` 前
**動作**：若有 `.cpp`/`.hpp` 待編譯標記存在 → `exit 2` 阻擋

**解除方法**：
```bash
cd build && make -j$(nproc)
# 成功後 compile_success_clear.sh 會自動清標記
```

**繞過警告**：`--no-verify` 會**不**執行 hook（僅在 Hard Gate 明示除外時用；見 [06_workflows/07_cpp_change_pdd.md](../06_workflows/07_cpp_change_pdd.md)）

---

## 📋 C++ 編輯標記流程

```
1. Edit src/core/XX.cpp
   └─ PostToolUse: cpp_edit_guard.sh → 建立 .claude/pending_compile 標記

2. cd build && make -j
   └─ PostToolUse(Bash): compile_success_clear.sh → 清除標記

3. git commit
   └─ PreToolUse: pre_commit_compile_check.sh
      ├─ 標記存在 → exit 2 阻擋
      └─ 標記不存在 → 通過
```

---

## 🎯 Evidence Ledger 同步檢查

**觸發**：Edit/Write `research/autoresearch/evidence_ledger.jsonl` 時
**檢查項**：JSON 欄位完整性、cycle_id 唯一、artifacts_path 存在

**違反時**：提醒訊息（非阻擋），由作者決定是否修正

---

## 術語守衛（terminology_guard.sh）

**觸發**：Edit/Write 時
**檢查**：專案術語一致性（例：避免 ISM vs InterSubMod 混用、HP tag 命名）

---

## Markdown 連結檢查（md_link_check.sh）

**觸發**：Edit/Write `.md` 時
**檢查**：`[text](path)` 是否指向存在檔案
**KB 用途**：補充 `scripts/check_canonical_paths.py`（後者僅檢 frontmatter）

---

## 新增 Hook 的流程

1. 撰寫 `scripts/hooks/<name>.sh`（標準 stdin 讀 JSON）
2. 編輯 `.claude/settings.local.json` 的 `hooks` 區段
3. 測試：執行對應 tool，檢視 hook 輸出
4. **謹慎**：阻擋型 hook（`exit 2`）須多輪測試，避免誤擋

---

## 關閉既有 Hook

編輯 `.claude/settings.local.json`，移除對應 `hooks[]` 條目，或改為 `|| true` 結尾（讓 hook 永不失敗）。

**已關閉**：`research_direction_guard.sh`（研究方向警告，目前 noop）

---

## 調試 Hook

```bash
# 模擬 tool 輸入
echo '{"tool_input": {"command": "git commit -m test"}}' | \
  bash /big7_disk/liaoyoyo2001/InterSubMod/scripts/hooks/pre_commit_compile_check.sh
echo "Exit code: $?"
```

---

## 相關

- PDD 6 steps（與 C++ commit Hard Gate 緊密關聯）：[../06_workflows/07_cpp_change_pdd.md](../06_workflows/07_cpp_change_pdd.md)
- Build/test：[../06_workflows/01_build_and_test.md](../06_workflows/01_build_and_test.md)
- Evidence ledger：[../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)
- 原始配置：`.claude/settings.local.json` hooks 區段
- Hook 腳本：`scripts/hooks/*.sh`
