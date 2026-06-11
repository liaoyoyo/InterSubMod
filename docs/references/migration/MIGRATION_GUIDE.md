# Claude Code 工作流與風格遷移指南

> **目的**：把當前環境的「**做事方式 + 風格偏好**」遷移到另一台新機 Claude Code，**不帶 InterSubMod 專案研究內容**。
> **適用**：新任務可能不同（不一定是癌症生信），但 confirmation 矩陣、暫停判斷、報告格式、HTML 整理、敘述框架、task type gate 等做事規則要完整保留。
> **撰寫日期**：2026-05-26
> **Framework**：Inventory → Classify (3-layer) → Migration How-to → Verification Checklist

---

## §1 環境盤點（當前實況）

### 1.1 目錄與檔案分佈

| 層級 | 路徑 | 內容 | 跨機遷移性 |
|------|------|------|-----------|
| **用戶級（user）** | `~/.claude/` (實際 `/bip7_disk/liaoyoyo2001/.claude/`) | settings.json (全域)、plugins/、daemon、history | ✅ 部分可移植（去專案路徑後） |
| **專案級（project）** | `<project>/.claude/` | CLAUDE.md、settings.local.json、skills/、agents/、commands/、rules/ | ⚠ 混合：規則可移、研究 skill 不可移 |
| **專案腳本** | `<project>/scripts/hooks/` | 35 個 .sh hook scripts（被 settings.local.json 引用） | ⚠ 混合：通用可移、ISM 綁定不可移 |
| **專案根目錄** | `<project>/AGENTS.md` `<project>/docs/CURRENT_FOCUS.md` | 跨 agent governance、live working state | ⚠ AGENTS.md 部分風格可抽 |
| **記憶（per-project）** | `~/.claude/projects/<encoded-path>/memory/` | 98 個 .md 含 MEMORY.md 索引；其中 41 feedback (風格) | ✅ 風格 41 個全可移；project/reference 不移 |

### 1.2 當前數量速覽

- **Skills**: 44 個（專案內 `.claude/skills/`）+ 全域 plugin skills（superpowers / hookify / plugin-dev 等）
- **Agents**: 17 個（專案內 `.claude/agents/`）
- **Commands**: 8 個（build / git-* / test-* / validate）— 全 ISM 流程
- **Hooks**: 30 個 hook scripts 跨 22 matchers（settings.local.json 定義；scripts 在 `scripts/hooks/`）
- **Rules**: 4 個 path-scoped rules
- **Memory**: 98 個（41 feedback 風格 + 38 project 研究 + 19 reference）

---

## §2 三層分類：可遷移 vs 部分遷移 vs 不遷移

### 🟢 第一層：**完全可遷移**（核心做事/風格）

#### 2.1.1 用戶級設定（`~/.claude/settings.json`）

```jsonc
{
  "language": "繁體中文",
  "effortLevel": "xhigh",
  "skipDangerousModePermissionPrompt": true,
  "remoteControlAtStartup": true,
  "agentPushNotifEnabled": true,
  "enabledPlugins": { ... }  // 21 plugins (見 §2.1.2)
}
```

**遷移動作**：直接 cp。permissions allow-list 需清洗（見 §3）。

#### 2.1.2 全域 plugins（21 個，已啟用）

| Plugin | 功能 | 風格相關性 |
|--------|------|------------|
| `superpowers` | brainstorming / TDD / debugging / verification-before-completion | ⭐⭐⭐ 核心做事 |
| `hookify` | conversation 行為偵測 → hook 生成 | ⭐⭐ |
| `plugin-dev` | 自製 plugin 工作流 | ⭐⭐ |
| `skill-creator` | skill 創建/優化 | ⭐⭐⭐ |
| `claude-md-management` | CLAUDE.md audit 工具 | ⭐⭐ |
| `claude-code-setup` | 推薦 Claude Code automation | ⭐⭐ |
| `pr-review-toolkit` | code review / silent-failure-hunter 等 | ⭐⭐ |
| `code-review` `code-simplifier` `feature-dev` `commit-commands` `frontend-design` | 程式開發風格 | ⭐⭐ |
| `clangd-lsp` `agent-sdk-dev` `context7` `playground` `ralph-loop` `security-guidance` | 工具類 | ⭐ |
| `github` `sourcegraph` `anthropic-skills` `stripe` `example-skills` | (已停用或專案特定) | ⭐ |

**遷移動作**：新機 `/plugin install <name>` 重灌；marketplace 為 `claude-plugins-official` + `anthropic-agent-skills`。

#### 2.1.3 通用做事 skills（`.claude/skills/`）— 17 個

可直接拷貝到新專案 `.claude/skills/`，無 InterSubMod 內容：

| 類別 | Skill | 功能 |
|------|-------|------|
| 元方法論 | `confirmation-protocol` | Hard Gate / Gate / Review / FYI 4-tier 確認協議 |
| 元方法論 | `scientific-rigor` | 證據分級 L1-L5 / Effect Size / Pre-registration |
| 元方法論 | `fast-learning-coach` | 快速深度學習導師（費曼+主動回想+80/20） |
| 元方法論 | `problem-framing-ideation` | 5W1H + gap analysis 收斂假說 |
| 元方法論 | `pre-decision-audit` | 決定前審查（≤30min 7 outputs + Cynefin gate） |
| 元方法論 | `known-pitfalls` | 已知陷阱清單（內容需清洗，留通用陷阱） |
| 元方法論 | `infra-ops` | 磁碟/基礎設施 preflight |
| 敘述框架 | `narrative-frame` ⭐ | 主入口 + 50+ framework catalog（取代固定範本） |
| 報告生成 | `weekly-report` | 週進度多主題報告（thin wrapper → Multi-Thread-Narrative） |
| 報告生成 | `structured-tech-report` | 13 段技術報告（thin wrapper → A3+ADR+Postmortem） |
| 報告生成 | `results-report` | 實驗結果報告（thin wrapper → Data-Showcase） |
| 報告生成 | `report` | AI 會話報告（thin wrapper → AI-Session-Companion） |
| 報告生成 | `myPPT` | PPT 場景路由器 |
| 報告生成 | `pptx-build` | PPTX 製作（含 §20 主軸聚焦 + Tier 1/2/3 分流） |
| 報告生成 | `html-report-build` | LLM-direct HTML 3 模式（report/slide/standalone PI） |
| 報告生成 | `image-gen` `image-vision-check` | 圖片生成 + 品質審核（codex AI + cairo 雙軌） |
| 文件管理 | `doc-standards` | 文件命名規範 / partial flag ribbon |
| 文件管理 | `memory-consolidation` | 記憶生命週期管理 |
| 文件管理 | `data-audit` | 檔案組織完整性檢核 |
| 文件管理 | `citation-verification` | 學術引用 .bib 驗證 |
| 視覺化 | `html-preview` | (已 deprecated → html-report-build 替代) |

#### 2.1.4 通用 rule（`.claude/rules/`）

| Rule | 內容 | 動作 |
|------|------|------|
| `opus47-behavior.md` | Opus 4.7 模型執行特性（literal 指令、少 tool calls 等） | ✅ 直接拷貝 |

#### 2.1.5 通用 hooks（`scripts/hooks/`）— 14 個

可拷貝後**改路徑**重用（hook 內 hardcoded `/big7_disk/liaoyoyo2001/InterSubMod/...`）：

| Hook | Event | 功能 | 通用性 |
|------|-------|------|--------|
| `session_start_inject_focus.sh` | SessionStart | 從 CURRENT_FOCUS.md 注入 live state | ⭐⭐⭐ 概念通用 |
| `md_path_format_rule.sh` | UserPromptSubmit | .md 路徑前綴規則 | ⭐⭐ |
| `narrative_frame_advisor.sh` | UserPromptSubmit | keyword 偵測 → 推薦 framework | ⭐⭐⭐ |
| `task_type_advisor.sh` | UserPromptSubmit | 6 類 task type keyword 偵測 + 注入 advisory | ⭐⭐⭐ |
| `verify_gate.sh` | PreToolUse Edit\|Write | Default-FAIL evidence gate | ⭐⭐⭐ |
| `pre_commit_compile_check.sh` | PreToolUse Bash | C++ commit 前強制編譯（C++ 專案通用） | ⭐⭐ |
| `no_binary_commit.sh` | PreToolUse Bash | binary commit 阻擋 | ⭐⭐⭐ |
| `external_input_sanitizer.sh` | PostToolUse WebFetch | injection 偵測 | ⭐⭐⭐ |
| `subagent_completion_logger.sh` | SubagentStop | cost / cache / artifact 紀錄 | ⭐⭐ |
| `evidence_read_tracker.sh` | PostToolUse Read | Read 追蹤（搭配 verify_gate） | ⭐⭐ |
| `memory_recall_logger.sh` | PostToolUse Read | 記憶引用率量化 | ⭐⭐ |
| `skill_change_audit.sh` | PostToolUse Edit\|Write | skill 變動月度 log | ⭐⭐ |
| `skill_usage_logger.sh` | PreToolUse Skill | skill 使用統計 | ⭐⭐ |
| `cache_telemetry.sh` | manual | prompt cache hit 統計 | ⭐⭐ |
| `allow_list_audit.sh` | manual | settings allow-list 158 entries audit | ⭐ |

#### 2.1.6 風格 feedback memory — 41 個（核心做事規範）

完整清單見 `portable/feedback_memory_to_migrate.md`。重點分群：

- **暫停判斷 / 確認協議**（5 個）：`execution_mode_hierarchy`、`strategy_then_per_item_confirmation`、`full_auto_parallel_execution`、`evidence_driven_iteration_workflow`、`task_first_then_doc_then_plan`
- **task type / scope**（3 個）：`observation_scope_default_comprehensive`、`risk_structured_iterative`、`small_scale_validation_first`
- **方法論嚴謹**（6 個）：`cynefin_domain_classification`、`productive_failure_reopen_threshold`、`existing_artifacts_must_verify`、`outside_claim_must_query_kb`、`researcher_claim_needs_empirical_verification`、`L2_collider_bias`、`pooled_ols_residualization_trap`、`spatial_autocorrelation_confound`
- **PPT / 簡報設計**（11 個）：`design_principles_canonical`、`ppt_minimal_visual_first`、`pi_report_html_stack`、`pptx_*` 系列 8 個
- **報告 / 寫作工作流**（5 個）：`weekly_report_workflow`、`md_path_prefix_rule`、`paper_positioning_de_prioritized`、`skill_md_must_state_dependencies_and_diagnostics`、`project_subfolder_structure`
- **視覺化**（3 個）：`matplotlib_cjk_font_rule`、`visual_inspection_requirement`、`figure_layout_standard`
- **基礎設施**（1 個）：`tmp_disk_full_pipeline_pitfall`
- **PI 報告 HTML stack**（1 個）：`pi_report_html_stack`

#### 2.1.7 通用 agents（`.claude/agents/`）

可移植 agents（不綁定 ISM）：`architect.md`、`developer.md`、`optimizer.md`、`reviewer.md`、`tester.md`、`researcher.md`、`parallel-analysis.md`、`parallel-benchmark.md`、`headless-research.md`、`research-orchestrator.md`、`narrative-organizer.md`、`literature-reviewer.md`、`paper-miner.md`、`evaluator.md`、`release.md`

ISM 綁定 agent：`methodology-reviewer.md`（ISM C++ 方法）

---

### 🟡 第二層：**部分遷移**（清洗後可用）

#### 2.2.1 `CLAUDE.md`（風格部分 vs 專案部分混合）

| 段落 | 性質 | 動作 |
|------|------|------|
| §0 任務類型 gate（6 類 task type） | ✅ 完全可移 | 拷貝後改 keyword（如 HKU 改成新場景） |
| §1 確認協議與暫停判斷矩陣 | ✅ 完全可移 | 直接拷貝 |
| §2 Opus 4.7 模型特性備註 | ✅ 完全可移 | 直接拷貝 |
| §3 Skills 分類索引 | ⚠ 部分可移 | 留通用類別、刪研究專用類別 |
| §4 Hooks 概覽 | ⚠ 部分可移 | 留 14 個通用 hook 描述 |
| §5 Rules 載入狀態 | ⚠ 部分可移 | 留 opus47-behavior，刪 cpp-build/workflow-commands/output-structure |
| §6 Working State Pointer | ✅ 概念可移 | 改名為新專案的 live state 檔案 |
| §7 Context 壓縮保留指令 | ✅ 完全可移 | 直接拷貝 |
| §8 Subagent 觸發 | ✅ 完全可移 | 直接拷貝 |
| §9 Agent 上下文 3 入口分工 | ✅ 概念可移 | 改檔案路徑 |
| §10 回應分級機制 | ✅ 完全可移 | 拷貝；引用 AGENTS.md §15 |
| §11 敘述框架預設啟用 | ✅ 完全可移 | 直接拷貝 |

→ 已產出清洗版模板：`portable/CLAUDE_md_template.md`

#### 2.2.2 `AGENTS.md`（風格部分可抽）

| 段落 | 性質 | 動作 |
|------|------|------|
| 跨 agent 語言規範 | ✅ 可移 | 拷貝 |
| Step→Verify 格式 | ✅ 可移 | 拷貝 |
| 假設陳述規則 | ✅ 可移 | 拷貝 |
| 證據敘述標準 | ✅ 可移 | 拷貝 |
| §15 回應分級機制 | ✅ 可移 | 拷貝 |
| §15.3 task type 6 類完整定義 | ✅ 可移 | 拷貝 |
| 專案研究目標 / KB 義務 / pipeline 規範 | ❌ 不移 | 重寫 |

#### 2.2.3 `~/.claude/settings.json` permissions allow-list

158 條 allow entries 含大量 ISM 絕對路徑（`/big7_disk/liaoyoyo2001/InterSubMod/...`、`/big8_disk/...`、`longphase-to-mod/...`）。

**清洗動作**：
- 留下：`Bash(ls:*) Bash(grep:*) Bash(find:*) Bash(awk:*) Bash(python3:*) Bash(git add:*) Bash(git commit:*) Bash(git checkout:*) Bash(git merge:*) Bash(clang-format:*) WebSearch Bash(gh repo:*) Bash(gh api:*) WebFetch(domain:github.com) WebFetch(domain:raw.githubusercontent.com)` 等通用條目
- **刪除**：所有 `/big7_disk/` `/big8_disk/` `/bip7_disk/` 絕對路徑、ISM 研究腳本路徑、HCC1395/COLO829/etc 樣本路徑
- 模板：`portable/settings_template.json`

#### 2.2.4 通用但需改路徑的 hooks（14 個 → 含 ISM 路徑）

所有 hook scripts 內含 `cd /big7_disk/liaoyoyo2001/InterSubMod && ...` 或 `LOG_FILE="docs/postmortems/..."` 類 hardcoded path。

**清洗動作**：
- 用 `sed -i 's|/big7_disk/liaoyoyo2001/InterSubMod|$PROJECT_ROOT|g' *.sh`
- 改用 `${CLAUDE_PROJECT_DIR:-$(pwd)}` 取代絕對路徑
- 模板 / 清洗腳本：`portable/hooks_to_migrate.md` + `portable/hook_path_rewrite.sh`

---

### 🔴 第三層：**不遷移**（純 InterSubMod 研究內容）

#### 2.3.1 研究專用 skills（17 個 — 純 ISM 工作流）

`auc-confound-guard` `feature-layered-observation` `multi-sample-consistency` `pivot-direction` `inject-hypothesis` `init-research` `review-evidence` `observation-analysis` `cycle-init` `research-loop` `check-staleness` `run-evaluator` `conclude-research` `validation-protocol` `cycle-state` `research-dashboard` `research-context-loader` `provenance-tier-audit` `methodology-audit` `verification-loop` `cpp-change` `results-analysis` `implementation-notes`

> 注：`results-analysis` `cpp-change` `verification-loop` `methodology-audit` 是 C++ 研究通用但內含 ISM 流程，新專案若是 C++ 可參考改寫，**不直接拷貝**。

#### 2.3.2 研究專用 rules

`cpp-build.md` `workflow-commands.md` `output-structure.md` — 含 ISM 路徑與腳本名稱

#### 2.3.3 研究專用 commands（8 個）

`build / test-quick / test-full / test-data / git-start / git-finish / git-commit / validate` — 全 ISM pipeline shortcut

#### 2.3.4 研究專用 hooks（16 個 — 綁定 ISM artifact）

`knowledge_check.sh` `kb_freshness_warn.sh` `kb_schema_check.sh` `kb_sot_guard.sh` `research_direction_guard.sh` `pipeline_block_check.sh` `evidence_ledger_sync.sh` `pre_tier_upgrade_check.sh` `cpp_edit_guard.sh` `compile_success_clear.sh` `post_cpp_commit_invalidate.sh` `terminology_guard.sh` `standalone_trigger.sh` `trigger_routing.sh` `researcher_claim_evidence_check.sh` `evidence_level_lint.sh` `causal_claim_check.sh` `compact_test.sh` `md_link_check.sh`

#### 2.3.5 研究記憶（57 個 — 純 ISM artifact）

- `project_*` 系列 38 個（active research / pending / concluded）
- `reference_*` 系列 19 個（ISM source code 行號等）

#### 2.3.6 ISM 專案根目錄

- `AGENTS.md` 內 InterSubMod 研究目標 / KB 義務
- `docs/CURRENT_FOCUS.md` live 主軸
- `state/` `research/` `docs/experiments/` `docs/reports/` 全部

---

## §3 遷移 How-to（新機 step-by-step）

### Step 1：準備 Claude Code 環境

```bash
# 1.1 新機安裝 Claude Code（依官方 docs）
# 1.2 確認版本
claude --version  # 需 Opus 4.7 支援

# 1.3 建立 home 目錄
mkdir -p ~/.claude
```

### Step 2：拷貝用戶級 settings + 安裝 plugins

```bash
# 2.1 從舊機拷貝 settings template（已清洗）
scp <old-machine>:/big7_disk/liaoyoyo2001/InterSubMod/docs/references/migration/portable/settings_template.json ~/.claude/settings.json

# 2.2 在新機啟動 claude 後，逐一安裝 plugins
claude
> /plugin install superpowers@claude-plugins-official
> /plugin install hookify@claude-plugins-official
> /plugin install plugin-dev@claude-plugins-official
> /plugin install skill-creator@claude-plugins-official
> /plugin install claude-md-management@claude-plugins-official
> /plugin install claude-code-setup@claude-plugins-official
> /plugin install pr-review-toolkit@claude-plugins-official
> /plugin install frontend-design@claude-plugins-official
> /plugin install code-review@claude-plugins-official
> /plugin install code-simplifier@claude-plugins-official
> /plugin install commit-commands@claude-plugins-official
> /plugin install feature-dev@claude-plugins-official
> /plugin install agent-sdk-dev@claude-plugins-official
> /plugin install context7@claude-plugins-official
> /plugin install playground@claude-plugins-official
> /plugin install ralph-loop@claude-plugins-official
> /plugin install security-guidance@claude-plugins-official
> /plugin install clangd-lsp@claude-plugins-official   # 若新任務含 C++
```

### Step 3：在新專案內建立 `.claude/` 框架

```bash
cd <new-project>
mkdir -p .claude/{skills,rules,hooks,agents,commands}
mkdir -p scripts/hooks
mkdir -p docs/{references,postmortems}
```

### Step 4：拷貝風格 skills（17 個）

```bash
# 在舊機把 portable skills 打包
cd /big7_disk/liaoyoyo2001/InterSubMod/.claude/skills
SKILLS_PORTABLE="confirmation-protocol scientific-rigor fast-learning-coach problem-framing-ideation pre-decision-audit known-pitfalls infra-ops narrative-frame weekly-report structured-tech-report results-report report myPPT pptx-build html-report-build image-gen image-vision-check doc-standards memory-consolidation data-audit citation-verification"
tar czf ~/portable_skills.tar.gz $SKILLS_PORTABLE

# 拷貝到新機
scp ~/portable_skills.tar.gz <new-machine>:<new-project>/.claude/skills/
ssh <new-machine> "cd <new-project>/.claude/skills && tar xzf portable_skills.tar.gz && rm portable_skills.tar.gz"
```

**清洗**：拷貝後跑 `grep -rn "InterSubMod\|HCC1395\|COLO829\|big7_disk\|big8_disk" .claude/skills/` 找出殘留專案 reference，手動審視（多在 known-pitfalls 案例中）。

### Step 5：拷貝通用 hooks（14 個）+ 改路徑

```bash
HOOKS_PORTABLE="session_start_inject_focus md_path_format_rule narrative_frame_advisor task_type_advisor verify_gate pre_commit_compile_check no_binary_commit external_input_sanitizer subagent_completion_logger evidence_read_tracker memory_recall_logger skill_change_audit skill_usage_logger cache_telemetry allow_list_audit"

cd /big7_disk/liaoyoyo2001/InterSubMod/scripts/hooks
for h in $HOOKS_PORTABLE; do cp ${h}.sh /tmp/portable_hooks/; done
tar czf ~/portable_hooks.tar.gz -C /tmp portable_hooks

# 新機拷貝後執行路徑改寫
ssh <new-machine> "cd <new-project>/scripts/hooks && tar xzf ~/portable_hooks.tar.gz --strip-components=1"
ssh <new-machine> "cd <new-project>/scripts/hooks && sed -i 's|/big7_disk/liaoyoyo2001/InterSubMod|$(pwd | head -c -13)|g' *.sh"
# 或統一改用 \${CLAUDE_PROJECT_DIR}
```

### Step 6：拷貝 CLAUDE.md（已清洗 template）

```bash
cp portable/CLAUDE_md_template.md <new-project>/.claude/CLAUDE.md
# 編輯：
#   - §0 task type 6 類保留，keyword 改成新場景的（如 HKU → 新對外 stakeholder）
#   - §6 Working State Pointer 改成新專案的 live state 檔案
#   - §3 Skills 索引按新專案實際 skill 重列
#   - §4 Hooks 概覽按新專案實際 hook 重列
```

### Step 7：建立新專案的 settings.local.json

```bash
# 從 portable/settings_template_local.json 起手
cp portable/settings_template_local.json <new-project>/.claude/settings.local.json
# 編輯：
#   - hooks 區段內路徑替換為新專案路徑
#   - enabledMcpjsonServers 按需移除（knowledge / biorxiv / ensembl 為 ISM 專用）
#   - permissions allow-list 按新專案實際工具加入
```

### Step 8：拷貝風格 memory（41 個 feedback_*）

```bash
# 新機 Claude Code 用戶記憶會自動建立 projects/<encoded-path>/memory/ 目錄
# 路徑編碼：用戶 home → /bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-<project_name>/memory/

# 1. 找出新機的 memory 目錄
NEW_MEM_DIR=$(ls -d ~/.claude/projects/*<new-project>*/memory/ | head -1)
mkdir -p $NEW_MEM_DIR

# 2. 拷貝 41 個 feedback memory（清單見 portable/feedback_memory_to_migrate.md）
OLD_MEM=/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory
FEEDBACK_LIST=$(grep -lE "^type: feedback" $OLD_MEM/feedback_*.md)
for f in $FEEDBACK_LIST; do cp $f $NEW_MEM_DIR/; done

# 3. 手動建立 MEMORY.md 索引（注意 <200 行限制）— 從 portable/MEMORY_template.md 起手
cp portable/MEMORY_template.md $NEW_MEM_DIR/MEMORY.md
```

**清洗**：feedback memory 內若引用 ISM 事件（如「5/24 HKU handoff」「Cycle 4 LR 失敗」），用戶可自行重新表述為通用 lesson。

### Step 9：可選 — 拷貝通用 rules + agents

```bash
# rules
cp portable/opus47-behavior.md <new-project>/.claude/rules/

# agents（按需）
cd /big7_disk/liaoyoyo2001/InterSubMod/.claude/agents
AGENTS_PORTABLE="architect developer optimizer reviewer tester researcher parallel-analysis parallel-benchmark headless-research research-orchestrator narrative-organizer literature-reviewer paper-miner evaluator release"
for a in $AGENTS_PORTABLE; do cp ${a}.md <new-project>/.claude/agents/; done
```

### Step 10：驗證新機行為

```bash
cd <new-project>
claude
> 隨意輸入測試 prompt，例如「整理今天做的事 → 應該觸發 narrative-frame advisor」
> /confirmation-protocol  # 確認 skill 載入
> 試一個 task type keyword 如「驗證新功能」→ 應該觸發 task_type_advisor 偵測 B validation
```

---

## §4 遷移後驗證 Checklist

- [ ] `claude /help` 顯示語言為繁體中文（用戶 settings.json `language` 生效）
- [ ] `/confirmation-protocol` skill 載入成功 → 4-tier 矩陣可查
- [ ] `/scientific-rigor` skill 載入 → §2-§7 完整
- [ ] `/narrative-frame` skill 載入 → catalog 50+ framework
- [ ] UserPromptSubmit hook 觸發：輸入「整理 X」應看到 `[narrative-frame advisor]` 提示
- [ ] UserPromptSubmit hook 觸發：輸入「完整驗證 X」應看到 `[task-type advisor]` 偵測 B validation
- [ ] SessionStart hook 觸發：開啟 session 應自動注入 working state（若 CURRENT_FOCUS.md 存在）
- [ ] PreToolUse Bash hook 阻擋：嘗試 `rm -rf /tmp/foo` → 應看到 [PreHook] 警告
- [ ] PostToolUse Edit|Write hook：編輯 .md 後應觸發 evidence_level_lint / md_link_check（若有開）
- [ ] PostToolUse SubagentStop hook：spawn 子代理結束後應看到 [SubagentStop] 提醒
- [ ] Plugin 載入：`/superpowers:brainstorming` / `/superpowers:using-superpowers` 可用
- [ ] Memory 載入：開新 session，問「我之前對 PPT 的偏好」應從 feedback_ppt_* 系列召回

---

## §5 三層 portable 附件

詳細模板存於 `portable/` 子目錄：

| 檔案 | 內容 |
|------|------|
| `portable/CLAUDE_md_template.md` | 已清洗的 CLAUDE.md 風格模板（§0-§11 留風格，刪 ISM） |
| `portable/settings_template.json` | 用戶級 `~/.claude/settings.json` 清洗版（無 ISM 路徑） |
| `portable/settings_template_local.json` | 專案級 `.claude/settings.local.json` 清洗版（hooks 配置 + 路徑替換點） |
| `portable/skills_to_migrate.md` | 17 個可遷移 skills 詳細清單 + tar 指令 |
| `portable/hooks_to_migrate.md` | 14 個可遷移 hooks 詳細清單 + 路徑改寫腳本 |
| `portable/feedback_memory_to_migrate.md` | 41 個風格 feedback memory 清單 + 拷貝指令 |
| `portable/MEMORY_template.md` | 新機 MEMORY.md 索引起手模板 |
| `portable/hook_path_rewrite.sh` | 自動把 hook 內絕對路徑改為 `${CLAUDE_PROJECT_DIR}` |

---

## §6 維護建議

### 6.1 哪些檔案應該定期跟舊機 sync？

| 檔案 | sync 頻率 | 機制 |
|------|----------|------|
| `~/.claude/settings.json`（plugin 啟用清單） | 安裝新 plugin 時 | 手動 git diff |
| `.claude/skills/`（17 個風格 skill） | 每 2 週 | git pull from 共用 repo |
| `.claude/hooks/` / `scripts/hooks/`（14 個風格 hook） | 修改後 | git pull |
| `feedback_*` memory（41 個） | 修改後 | 手動 cp（per-project memory 無 git） |

### 6.2 建議建立 dotfiles repo

把以下打包成獨立 git repo（如 `claude-code-portable-config`）：

```
claude-code-portable-config/
├── settings.json                  # ~/.claude/settings.json 模板
├── CLAUDE.md                      # 風格模板
├── skills/                        # 17 個風格 skill
├── hooks/                         # 14 個風格 hook
├── rules/opus47-behavior.md
├── agents/                        # 15 個風格 agent
├── memory/feedback_*.md           # 41 個風格 memory
├── MEMORY_template.md
├── settings_local_template.json   # 專案級 settings 模板
└── install.sh                     # 一鍵 install 腳本
```

新機只需 `git clone` + `bash install.sh <new-project-dir>` → 自動拷貝 + 路徑改寫 + plugin 安裝。

### 6.3 風格 drift 偵測

建議每月跑一次 `diff` 比對：

```bash
diff -r /bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-<project>/memory/ \
        ~/claude-code-portable-config/memory/ | grep feedback
```

→ 偵測新出現的 feedback memory 是否該回 push 進 dotfiles repo。

---

## §7 業界對照

| 工具 | 對應概念 |
|------|---------|
| Cline | Memory Bank `activeContext.md` ≈ 本專案 `CURRENT_FOCUS.md` |
| Cursor | `.cursorrules` ≈ `.claude/rules/` (path-scoped globs) |
| Aider | `CONVENTIONS.md` ≈ `AGENTS.md` (跨 agent governance) |
| Continue.dev | `config.json` ≈ `~/.claude/settings.json` |
| GitHub Copilot | `~/.config/github-copilot/` ≈ `~/.claude/` |

**SoT 原則**（Anthropic + OpenAI 跨平台共識）：
- One source of truth per scope（不要 5 入口，3 入口 = governance / agent-specific / live state）
- Path-scoped rules 優於 always-loaded
- Memory 為跨 session 持久層，不要塞當前 task 細節

---

**結尾 next-step**：先看 §2-§3 確認分類正確，再依 Step 1-10 執行；過程中若有 InterSubMod-specific 內容未清洗乾淨，回報修正模板。
