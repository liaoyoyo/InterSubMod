# CLAUDE.md — Claude Code 行為規範

> **此檔職責**：僅約束 Claude Code 特定行為（確認矩陣、Skills/Hooks 機制、Opus 4.8 特性）。
> **跨 agent governance（專案目標、研究主軸、Step→Verify、KB 義務、Output 規範等）→ `InterSubMod/AGENTS.md`**。
> **🗺️ 任何任務的主工作流地圖（agentic-loop 四層 + git 多工/分支/commit/合併四決策表 + 鉤子機制）→ `InterSubMod/docs/references/20260611_master_workflow_architecture_01.md`**（2026-06-11 收斂為唯一 SoT；動 git 前查其 §3 四表）。
> Claude Code 啟動時 concatenate 載入兩檔。

---

## §0 任務類型 gate（2026-05-26 5/24 incident postmortem 落地 — prerequisite）

**規則**: 任何 conversation 啟動時 / 接到新任務時，**先做 task type classification** 再進 §1 暫停判斷矩陣。task type 決定 default scope 與 deliverable 義務。

**6 類速查**（完整定義見 `InterSubMod/AGENTS.md §15.3`）:

| Type | Keyword | Default Scope |
|------|---------|--------------|
| A pilot | pilot / 探索 / 試 / 先看 | Subset OK |
| **B validation** | 完整 / final / validation / 驗證 | **全基因組 + 全樣本** |
| **C production** | release / deploy / merge | **全 scope + cross-platform** |
| **D handoff** | handoff / 交付 / HKU / external | **全 scope + 説明文件** |
| E hotfix | bugfix / hotfix / urgent | 最小可重現 |
| F demo | demo / 示意 / 教學 | Subset OK + DEMO 標 |

**強制行為**:
1. **啟動先分類** — 不可省略；模糊 → AskUserQuestion 必問 task type
2. **AskUserQuestion 強制含 scope 維度** — 觀察/報告類任務問卷必含「scope 全 vs subset？」
3. **subset 必標 partial flag** — `/doc-standards` ribbon
4. **「見樹也見林」四層** — aggregate / canonical / extreme outlier / well-explained

**主動 recall**: 啟動 prompt 含 keyword（完整 / 報告 / 整理 / handoff / 驗證 / final / 對外）→ 必 recall `feedback_observation_scope_default_comprehensive` 與 `feedback_task_first_then_doc_then_plan` 兩條 memory。

**業界對照**: Cynefin domain + DACI/DECIDE + Bland Assumption Map。完整 governance → `InterSubMod/AGENTS.md §15.3`。

---

## §1 確認協議與暫停判斷（Claude Code 特定）

**模式觸發詞**: `互動模式`（預設） / `全自動` / `auto`
**級別**: Hard Gate（必停）> Gate（互動停/自動過）> Review（互動停/自動過）> FYI

**Hard Gate 必停場景**（不可逆操作永遠 🔴）：
1. **刪除檔案**
2. **C++ commit**（pre-commit hook 強制編譯檢查）
3. **研究方向 NO-GO 判定**

**暫停判斷矩陣**（影響度 × 信心度）：
| 影響 ＼ 信心 | 高 | 中 | 低 |
|---|---|---|---|
| **低**（可逆/<10min） | 🟢 一行告知 | 🟡 列假設 | 🟠 節點暫停 |
| **中**（10min-1h） | 🟡 列假設 | 🟠 節點暫停 | 🔴 立即暫停 |
| **高**（>1h/影響結論） | 🟠 節點暫停 | 🔴 立即暫停 | 🔴 立即暫停 |

**可逆性 override**: 刪檔／C++ commit／NO-GO／覆寫 evidence_ledger 永遠 🔴。

**長計算例外**: 用戶**當輪明示**啟動長計算（「跑全量」「平行 7 樣本」）→ 不二次確認，直接執行 + 一行告知。模型自判需長計算 → 仍 🔴。

**一行告知格式**: `[決策]（影響: 低/中/高, 信心: 低/中/高, 理由: 一句）`

完整維度定義 / 場景對照 → `/confirmation-protocol` skill。

---

## §2 Opus 4.8 模型特性備註（2026-05-30 由 4.7 對齊；model = claude-opus-4-8）

> **Opus 4.8 literal 特性（4.8 更甚）**：模型不會推斷未明講需求、不會悄悄泛化指令。模糊輸入會被按字面執行 — **假設陳述是唯一屏障**。

- **預設不 spawn subagent**（單回合優先；大 fan-out 用 Dynamic Workflow，見 §8）
- **預設少 tool calls**（reasoning 解決優先於反覆讀檔）
- **回應長度動態化**（不再固定冗長）
- **首 turn 規格完整度檢核**：規格缺項 ≥2 且高影響 → 必須回問
- **Effort 比前代更關鍵**：project 已 pin `xhigh`；遇淺推理 raise effort 而非 prompt around it；深推理用 `ultrathink` 關鍵字
- **Adaptive thinking 預設 OFF**（取代 4.7 的 budget_tokens / display:summarized）

完整 Step→Verify 格式 / 假設陳述規則 / 證據敘述 → **`InterSubMod/AGENTS.md` §6**。
完整模型執行特性 → **`.claude/rules/opus48-behavior.md`**（條件式載入；取代已停用的 opus47-behavior.md — 詳見 §5）

---

## §3 Skills 分類索引（49 個 SKILL.md，Claude Code 特定 — 2026-06-19 +task-graph）

> **drift 修正紀錄**：
> - 2026-05-20：原寫 45 → 一度改 44（誤判 `grill-me/` 為 phantom +1）。
> - **2026-05-30 實測校正**：`find .claude/skills -name SKILL.md | wc -l` = **45**；46（含 deprecated html-preview）→ git rm `/html-preview` 後 = 45 全 active。校正機制：`creation_guard.sh`（Write 新 SKILL.md 提醒）+ `skill_registry_sync.sh`（編 README/CLAUDE.md 時比對磁碟實際計數）。
> - **2026-05-30 清理**：`grill-me` 曾是 dangling symlink（→ `../../.agents/skills/grill-me`，target 不存在），**已 `git rm` 移除**。現 `find -maxdepth 1 -type d` 與 `find -name SKILL.md` 一致 = 45（無 orphan symlink）。
> - 新增 12 個未分類 skills 進對應類別（2026-05-20）。
> - **2026-05-31**：新增 `/harness-health`（元方法論 9→10）→ 45→**46**。磁碟現有 **2** 個 Dynamic Workflow（`cross_sample_benchmark.js` + 新增 `harness_audit_2026.js`；§8 未硬編碼計數，故無數字需改）。`harness_health.py` 燈 #1 持續監看此計數 drift。
> - **2026-06-05**：新增 `/methods-example`（視覺化 3→4）→ 46→**47**。方法解釋圖 generate+verify 整體 skill（物件庫 primitive + 案例模板 case + data_ref 注入 renderer〔缺 verified 真值 refuse〕+ 圖例細節 verify loop）；C1 BRCA2 Δβ pilot 已驗證可跑。
> - **2026-06-19**：新增 `/task-graph`（視覺化 5→6）→ 48→**49**（49 = 47 + verify-workstation〔2026-06-15〕+ task-graph）。研究任務有向圖層：`state/tasks/graph.json` 機械真值 → 自動產生 `TASKS.md`（留底）+ `tasks_board.html`（確認+顯示驗證主介面）。補 §盤點出的唯一大洞 D2（跨任務有向依賴）+ 程式流程圖（CWL-style 必需/可選 I/O）+ 主任務 `<details>` 階層；HTML 複用 `build_workstation.py`（§13-A 反捏造）、流程圖為純手刻 SVG 零依賴。`scripts/tasks/task_graph.py` = validate/check/ready/render/render-html。
> **重複交叉位置**：`/feature-layered-observation`（P3 + 研究專用）/ `/multi-sample-consistency`（P4 + 研究專用）/ `/pre-decision-audit`（元方法論 + 三層樓 pre）/ `/run-evaluator`（P5 + 三層樓 post）— 在多分類列出表示同 skill 多角色。

- **元方法論（10）**: `/confirmation-protocol` `/known-pitfalls` `/cycle-state` `/research-context-loader` `/fast-learning-coach` `/scientific-rigor` `/pre-decision-audit` `/problem-framing-ideation` `/provenance-tier-audit` `/harness-health` ⭐ 新（2026-05-31；harness 自我稽核 8 燈儀表板，read-only `scripts/harness_health.py`；2026-06-03 +燈#7 memory-drift +燈#8 doc-path-currency）
- **7-Phase Waterfall（7）**: P0 `/cycle-init` → P1 `/research-loop` → P2 `/check-staleness` → P3 `/feature-layered-observation` → P4 `/multi-sample-consistency` → P5 `/run-evaluator` → P6 `/conclude-research`
- **程式開發（4）**: `/cpp-change` `/methodology-audit` `/infra-ops` `/verification-loop`
- **文件管理（5）**: `/doc-standards` `/data-audit` `/memory-consolidation` `/citation-verification` `/pipeline-manifest` ⭐ 新（reproducibility provenance DAG；與 data-audit 分工：data-audit 查組織、pipeline-manifest 查 script→figure 因果鏈）
- **報告生成 retrospective（7）**: `/weekly-report` → `/pptx-build` / `/html-report-build` / `/results-report` / `/structured-tech-report` / `/report` / `/myPPT`
- **視覺化（6）**: `/image-gen` `/image-vision-check` `/research-dashboard` `/methods-example` `/verify-workstation` `/task-graph` ⭐ 新（2026-06-19；研究任務有向圖層 — `state/tasks/graph.json` 機械真值 → 自動產生 TASKS.md 留底 + tasks_board.html 確認介面；主任務 `<details>` 階層 + 程式流程圖〔CWL-style 必需/可選 I/O，純手刻 SVG〕+ 完整依賴 DAG + 細節缺資訊清單；複用 build_workstation.py；driver `scripts/tasks/task_graph.py`；補盤點出的唯一大洞 D2 跨任務有向依賴；是人工分派地圖非自動執行器）（verify-workstation 2026-06-15；觀察驗證介面 — 一批結果生成互動逐項判讀工作站：人工判讀 localStorage + JSON/CSV 匯出 + 嵌觀察圖 + 修正過程紀錄 changelog + §13-A 注入/refuse-on-missing；通用 generator `tools/build_workstation.py`，worked 實例 chr2:18M〔inline SVG〕+ HCC1395 ASM〔30,350 位點外部 PNG，script 102〕）（methods-example 2026-06-05；html-preview 2026-05-30 移除→/html-report-build）
- **研究專用（8）**: `/auc-confound-guard` `/feature-layered-observation` `/multi-sample-consistency` `/pivot-direction` `/inject-hypothesis` `/init-research` `/review-evidence` `/observation-analysis`
- **資料分析 / 驗證（2）**: `/results-analysis` `/validation-protocol`
- **假說驗證三層樓（pre → process → post）⭐ 2026-05-19**:
  - **Pre** (entry-point ≤30min): `/pre-decision-audit` (7 outputs + Cynefin gate + 5-dim credibility + GO/PROBE/NO-GO verdict)
  - **Process** (spec 實作中 live): `/implementation-notes` (4 sections + Lore — 設計決定 / 偏離 / 折衷 / 未決)
  - **Post** (P5 cycle 結束後): `/run-evaluator` (tier ⭐1-5 + 6 risk components)
- **敘述框架庫 ⭐ 2026-05-22**: `/narrative-frame`（主入口 + 50+ catalog；取代既有 7 報告類 skill 固定範本）
  - 7 thin wrapper skill: `/structured-tech-report`（→ A3+ADR+Postmortem-hybrid）/ `/weekly-report`（→ Multi-Thread-Narrative）/ `/pptx-build`（→ Audience-Scenario-Pitch）/ `/results-report`（→ Data-Showcase）/ `/conclude-research`（→ Verdict-Pyramid）/ `/report`（→ AI-Session-Companion）/ `/myPPT`（場景路由）
  - 1 sub-agent: `narrative-organizer`（≥3 文件並行萃取 + cross-file 主題聚類）
  - 1 hook: `narrative_frame_advisor.sh`（UserPromptSubmit 偵測 keyword 推薦套 framework）
  - 與 `/pre-decision-audit` decision 層正交（decision 決定要不要做 + narrative-frame 決定怎麼講）

---

## §4 Hooks 概覽（Claude Code 特定 — 2026-05-18 P4 完整收尾）

依 `InterSubMod/.claude/settings.local.json` 完整定義（含 SessionStart / UserPromptSubmit / PreToolUse / PostToolUse / SubagentStop / Stop / **PreCompact** 7 個 events；**45 hook scripts**，2026-06-23 實測 `ls scripts/hooks/*.sh | wc -l`（2026-06-23 +`script_lint_advisor.sh`〔.py/.sh post-edit advisory lint，治理稽核補〕 +`palette_drift_advisor.sh`〔ISM renderer 色盤漂移 advisory〕；2026-06-15 +`provenance_stamp_advisor.sh`）；2026-06-11 +`git_branch_commit_guard.sh`〔主線保護 exit-2 + 並行警告〕；2026-06-09 +`health_drift_advisor.sh`〔C7〕 +`concurrent_session_advisor.sh`〔§G〕；2026-06-01 +`number_provenance_check.sh`）。

**2026-06-09 新增 2 個 hook**（loop-engineering ADR C7 + git governance §G — 修 invocation-dependence 漂移 + 並行 session 互撞）:
- `health_drift_advisor.sh`（SessionStart — **變動觸發的唯讀漂移 advisor**：偵測 `.claude/{skills,agents,rules,hooks,workflows}` + settings/CLAUDE.md/AGENTS.md 自上次 marker 後是否變動 → nudge 跑 `/harness-health`；順帶提醒 active cycle >7d 未推進。**advisory-only 永遠 exit 0、不跑 harness_health.py 本體（零延遲/不 hang/不 spam snapshot）、marker debounce ~once/change-batch**。非 `/loop`/`/goal`/cron。設計依據：`InterSubMod/docs/plans/20260609_loop_engineering_research_cycle_architecture_review_01.md` §5）
- `concurrent_session_advisor.sh`（SessionStart — **並行 session 互撞 advisor**：偵測主 repo transcript dir（`/bip7_disk/.../projects/-big7-...-InterSubMod`）>1 活躍 .jsonl（worktree session 的 transcript dir 名不同故自動排除）→ 提醒並行 session 各開 git worktree，避免共用 HEAD 讓 commit 落錯 branch。讀 stdin transcript_path 精準排除自己；**advisory-only exit 0、fail-OPEN**。落地 2026-06-09 跨 session commit 污染事故；見 git governance §G）

**2026-05-28~30 新增 4 個 hook**（社群 gap + 搜尋紀律）:
- `search_scope_guard.sh`（PreToolUse Bash — **exit 2 阻擋** `grep -r .` / `find .` 無 maxdepth / `du */` 等不尊重 .gitignore 的 repo-root 遞迴搜尋；§12 搜尋紀律強制層）
- `precompact_autosave.sh`（PreCompact 事件 — 壓縮前自動 dump active cycle + CURRENT_FOCUS 快照到 `state/compact_snapshots/`；對齊 §7 保留指令；gap G6）
- `creation_guard.sh`（PreToolUse Write — 新建 SKILL.md/agent 前 dedup + 計數同步提醒；防 §3 計數 drift；gap G5）
- `skill_registry_sync.sh`（PostToolUse Edit|Write — 編輯 README/CLAUDE.md 時比對「磁碟實際 skill/agent 計數 vs 文件宣稱最大數」，雙向 drift 守衛；2026-05-30 gap G5 補完）

**2026-06-01 新增 1 個 hook**（數據捏造防呆 — 落地 `20260601_fabricated_metric_in_html_preview_postmortem` A4/A5/A7）:
- `number_provenance_check.sh`（PreToolUse Edit|Write `*.md`/`*.html` — **分級 anti-fabrication gate**：抽報告內「metric 形數字」(`AUC=` `p=` `%` `Δ` `≥2 位小數`)，逐一去 bounded 來源（frontmatter `data_sources:` / 同層 `_assets/` / 同目錄 `.json|.tsv|.txt|.csv`）grep；**validated/pi_reports 路徑找不到來源 → exit 2 阻擋**；其他路徑 → advisory exit 0 提醒。fail-OPEN（python 缺/解析錯 → exit 0，絕不擋所有 Write — 與 neutering bug 本質不同）。override：內文含 `<!-- provenance-verified: 理由 -->`。三層防線見 §13。)

**2026-06-15 新增 1 個 hook**（工作流稽核 D6-2 落地 — provenance stamp 自動 nudge）:
- `provenance_stamp_advisor.sh`（PostToolUse Edit|Write — **advisory exit 0**：audit/status/盤點/manifest/稽核/inventory 類 `.md`/`.html` 缺 `build_branch:` stamp → 印 `scripts/provenance_stamp.sh` 可貼 stamp。窄 scope（只這幾類）+ fail-OPEN 永不擋 Write。把 provenance_stamp.sh 從 manual-only → 自動提醒，防 P-17 跨-worktree 幻覺。⚠ settings wiring 已加但因 settings.local.json 有 65 行未提交 churn（4 個已宣稱 wired hook 的 wiring）+ 並行 session，暫**不 commit settings**，live 生效；待 settings 乾淨時連同提交。)

**Hard Gate hooks**（不可繞過 — 真 `exit 2` 阻擋；2026-06-01 audit = **4 核心** + search/tier/number-provenance 三個條件式 exit-2）:
- `pre_commit_compile_check.sh`（C++ commit 必編譯）
- `kb_schema_check.sh`（KB 寫入前 schema 檢核）
- `pipeline_block_check.sh`（長 pipeline 磁碟檢核）
- `no_binary_commit.sh`（commit binary 阻擋）
- `search_scope_guard.sh`（exit 2 阻擋 repo-root 遞迴搜尋；§12，亦 exit-2 但屬搜尋紀律類）
- `pre_tier_upgrade_check.sh`（**2026-05-31 wired** — Edit/Write `state/cycles/*/state.json` 含 ⭐4/5 tier 但無 `evaluation.json`(verdict=approve_tier) 或 override 註解 → exit 2；INDEX/CURRENT_FOCUS 散文路徑降 advisory 不擋。研究誠信 gate）
- `number_provenance_check.sh`（**2026-06-01 wired** — validated/pi_reports 報告含無來源 metric → exit 2；in_progress → advisory。數據誠信 gate，見上 + §13）
- `git_branch_commit_guard.sh`（**2026-06-11 wired** — PreToolUse Bash，commit 前：直接 commit 到 `main`/`master`/`develop` → exit 2 阻擋（本 repo 一律 feature branch）；偵測並行 session → advisory warn。fail-OPEN〔非 commit/解析錯 → exit 0〕；只在 `commit` 為 git 子命令時觸發。落地主工作流 §3-C/§4③）

> ⚠ **2026-05-31 校正 + 修復**：
> - `kb_sot_guard.sh`（F1 SoT）與 `verify_gate.sh`（evidence gate）**本就是 advisory**（全 `exit 0` / 自宣告 SOFT），非 exit-2；文件曾誤列為 Hard Gate。
> - 🔴 **修復 neutering bug**：`pre_commit_compile_check` + `kb_schema_check` 的 settings wiring 先前是 `2>/dev/null || exit 0` → `|| exit 0` 把 script 的 `exit 2` 吃成 `exit 0` → **兩個 Hard Gate 實際從不阻擋**（C++ compile gate de-facto 失效一段時間）。已移除 mask（改 `1>&2` 讓阻擋原因 surface）。
> - 🔴 **修復 lifecycle bug**：`compile_success_clear.sh` 讀錯欄位 `.tool_result.exit_code`（應 `.tool_response`）→ jq 永遠 fallback `"1"` → marker 從不被清 → 自 2026-05-10 累積 32-檔 stale marker。已改 `.tool_response` + 多重 success 訊號 cascade + 清掉 stale marker。
> - 磁碟實況（修復後）= **6 個真 exit-2**：pre_commit_compile_check / kb_schema_check / pipeline_block_check / no_binary_commit（4 核心）+ search_scope_guard + pre_tier_upgrade_check。`/harness-health` 燈 #2 永久監看 neutering 復發。

**已落地 hooks**（2026-05-18 P0-P4 新增 11 個）:
- `session_start_inject_focus.sh` ✅（SessionStart CURRENT_FOCUS 注入，commit ee648fb）
- `md_path_format_rule.sh` ✅（UserPromptSubmit 路徑前綴）
- `skill_change_audit.sh` ✅（PostToolUse skill 變動月度 log）
- `verify_gate.sh` ✅（PreToolUse Edit|Write Default-FAIL evidence gate）
- `evidence_read_tracker.sh` ✅（PostToolUse Read 追蹤）
- `cache_telemetry.sh` ✅（manual telemetry, 證實 96.8% hit）
- `subagent_completion_logger.sh` ✅（SubagentStop cost/cache 紀錄）
- `compact_test.sh` ✅（manual /compact 後 preservation 驗證）
- `allow_list_audit.sh` ✅（manual 158 entries audit）
- `researcher_claim_evidence_check.sh` ✅（PostToolUse 偵測 hedge 語言）
- `memory_recall_logger.sh` ✅（PostToolUse 引用率量化）
- `external_input_sanitizer.sh` ✅（PostToolUse WebFetch injection 偵測）
- `task_type_advisor.sh` ✅ **2026-05-26 新增**（UserPromptSubmit task type 6 類 keyword 偵測 → 注入 advisory + telemetry log；A6 落地 5/24 incident postmortem）

---

## §5 Rules 載入狀態（path-scoped 已落地 2026-05-18）

| Rule 檔案 | 狀態 | globs 條件 |
|----------|------|-----------|
| `.claude/rules/cpp-build.md` | **條件式載入** ✅ | `src/**/*.cpp`, `src/**/*.hpp`, `include/**/*.hpp`, `include/**/*.h`, `tests/**/*.cpp`, `CMakeLists.txt` |
| `.claude/rules/opus48-behavior.md` | **條件式載入** ✅（2026-05-30 取代 opus47）| `.claude/skills/**/SKILL.md`, `.claude/skills/**/*.json`, `.claude/rules/**/*.md`, `.claude/skills/**/playbook.md`, `.claude/skills/**/prompts/*.md` |
| ~~`.claude/rules/opus47-behavior.md`~~ | **已 git rm 2026-05-30**（由 opus48-behavior.md 取代）| — |
| `.claude/rules/workflow-commands.md` | **條件式載入** ✅ | `scripts/**/*.sh`, `scripts/**/*.py` |
| `.claude/rules/output-structure.md` | **條件式載入** ✅ | `output/**/*`, `results/**/*` |

✅ **2026-05-18 P3 audit M2 落地**：4 個 rules 全部 path-scoped（採 `globs:` Cursor spec，與 `paths:` Anthropic skill spec 等價）。
✅ **省 ~5KB always-loaded**：rules 改為條件式後，僅在編輯對應 path 時才載入。
✅ **C++ Hard Gate 保留**：即使 `cpp-build.md` 變條件式（編輯 src/ 才載入），C++ commit Hard Gate（§1 + pre-commit hook）**仍適用** — Hard Gate 不依賴 rule 載入，而是不可逆操作的絕對規則

---

## §6 Working State Pointer（動態狀態 — 2026-05-18 SessionStart hook 已落地）

**當前主軸**: 由 SessionStart hook 從 `InterSubMod/docs/CURRENT_FOCUS.md` 自動注入（截 ~3000 chars / ~500 tokens）✅
**Pending Items**: 見 MEMORY.md 索引
**研究方向變更**: 觸發 `/pivot-direction`
**詳細上下文載入**: 觸發 `/research-context-loader`（Tier 1 / Tier 2 / Tier 3 深度）

✅ **SessionStart hook 狀態**: **已落地** (commit ee648fb, 2026-05-18)
✅ **mtime > 7 天**: hook 自動加 STALE WARNING（CURRENT_FOCUS.md weekly cadence 強制）
⚠ **何時不需詳細上下文**: 純 code edits（make/test/commit）、單檔 doc 寫作、簡單問答

⚠ **變動頻率上限：週級** — CURRENT_FOCUS.md > 7 天未更新 → hook 主動提醒。

> **業界對照**: `docs/CURRENT_FOCUS.md` 概念等同 [Cline Memory Bank `activeContext.md`](https://docs.cline.bot/features/memory-bank)（live working state，最高變動層）。命名保留 `CURRENT_FOCUS.md` 以維持既有 20+ 引用，不重新命名。

> **🔑 雙層 SoT 分工（2026-06-09 明文化 — loop-engineering ADR §8 resolution）**：避免「機械 cycle 狀態機」與「敘述層」漂移（曾凍在 06-02 vs 工作跑到 06-08）：
> - **多週論文主軸（如 G6/G1）的權威 SoT = `CURRENT_FOCUS.md` + `evidence_ledger.jsonl` + memory**（敘述層；變動頻率高、粒度大）。
> - **`state/cycles/` + `active.json` 狀態機 = 僅用於真正的短驗證 hypothesis cycle**（P0→P6 被實際驅動者）；**不再把多週主軸 backfill 註冊成 live cycle 假裝在跑**。
> - active.json 內現存的 G6/G1 = 06-02 backfill placeholder，其真實進度看 CURRENT_FOCUS（C7 `health_drift_advisor` 會在 >7d 未推進時提醒檢視/archive）。
> 設計依據：`InterSubMod/docs/plans/20260609_loop_engineering_research_cycle_architecture_review_01.md`。

---

## §7 Context 壓縮保留指令（Claude Code `/compact`）

`/compact` 必保留:
- 架構決策 + 理由 / 未解決 blockers / 涉及檔案路徑
- 用戶限制條件 / 假說 ID + 驗證層級 / 未完成步驟
- 本檔 §1-§7 規則 + `AGENTS.md` §1-§N governance

**可安全壓縮**: 中間計算 / 已通過完整測試輸出 / 冗長工具呼叫結果。

**跨 session 交接骨架**: `InterSubMod/docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md`

---

## §8 Subagent 觸發 + Dynamic Workflow 路由（2026-05-30 更新）

> **🔴 並行 session = git worktree 標準作業（2026-06-25 架構稽核落地；原「建議」升「標準」）**：多個 Claude Code session 同時開時，**每個 session 各開自己的 git worktree**，否則共用 HEAD → commit 交錯落錯 branch + settings.local.json live 編輯互相 clobber（2026-06-25 in-the-wild：另一 session 移除本 session 的 hook wiring + palette_drift 重工）。
> - **官方原生機制**（`https://code.claude.com/docs/en/worktrees`）：CLI `claude --worktree <name>`（worktree 落 `.claude/worktrees/worktree-<name>/`）；對話中用 `EnterWorktree` 工具；`worktree.baseRef:"head"` 帶未推 commit；無變更自動清理；desktop app **每 session 自動開 worktree**。gitignored config 用 `.worktreeinclude` 自動複製（tracked 檔如 settings.local.json 本就在每 worktree 各一份、不互撞）。
> - `concurrent_session_advisor.sh`（§4，SessionStart）偵測 >1 活躍 .jsonl 時提醒；**advisory 不強制**（fail-OPEN 故意不擋單人多視窗）。
>
> **並行度上限（外部實證；Augment Code 2026 multi-agent fail 研究）**：並行 subagent/session **建議 ≤2-4**（再多 overhead>benefit，~37% 失敗為協調類）；**dependent task 串行 merge 非同時 merge**（pipeline() 的 barrier 紀律）。

**預設不 spawn**（單回合優先）。**明示觸發語**才啟動：
- 跨樣本平行 benchmark → `parallel-benchmark`
- 大範圍探索（>3 query）→ `Explore`
- PR 審查 → `pr-review-toolkit`
- 夜間/離線研究 → `headless-research`

> **`.claude/agents/` 實有 18 個 project agent**（2026-05-30 實測 `ls .claude/agents/*.md`）：architect / developer / optimizer / tester / researcher / reviewer / evaluator / methodology-reviewer / parallel-benchmark / parallel-analysis / headless-research / narrative-organizer / research-orchestrator / literature-reviewer / paper-miner / release / **reproducibility-audit** ⭐ / **security-reviewer** ⭐（後 2 為 2026-05-30 新增 fresh-context read-only agent）。新增/刪除 agent 必同步本清單（`creation_guard.sh` 會在 Write agent 時提醒）。
>
> **新增的 read-only 機械型 agent 用 `model: haiku` 降 token**：`research-orchestrator`（純路由）+ `narrative-organizer`（純萃取）；驗證/coding agent（evaluator/reviewer/methodology-reviewer/architect/developer 等）維持 `inherit` 不降智力。

**任務切割原則**：可切割且需要大量 context 處理 → 啟動 agent + **必須清楚回報主 agent** + **科學工程紀錄成文件**供檢核驗證。

**Sub-agent return-contract（2026-05-30）**：subagent 回主 agent 時 — (a) full detail 落地成 OUTPUT_DIR 下 .md/.json；(b) 回 parent **只**回 `{status, key metrics/verdict, anomalies, path-to-landed-doc}`；(c) ~1-2K token soft target（**非 hard cap** — 勿截斷 parallel-benchmark 的多樣本 canonical 表）。`subagent_completion_logger.sh` 已捕 OUTPUT_TOKENS 可作 >3K advisory flag。

**Concise-emit（2026-06-11 — 對抗 output-token 上限 session 丟失）**：任何單次回覆/工具輸出**勿超 output-token 上限**（曾整 session 只剩 error 訊息）；大量內容**寫檔回 path** 而非塞進回覆；長 analysis 切「每單元落檔→下一批」。此規則把上述「detail 落檔」明確綁到 output-token 錯誤類（主工作流 §5）。
> **4-狀態 status enum（2026-06-02 借鑑 superpowers subagent-driven-development）**：`status` 用 `DONE` / `DONE_WITH_CONCERNS`（非阻擋疑慮，附說明）/ `NEEDS_CONTEXT`（缺資訊，parent 補後 re-dispatch）/ `BLOCKED`（需升級或重構任務）四值，取代模糊的 done/fail。多階段 review **spec-compliance 先過，才進 code-quality**（禁顛倒 — generator/evaluator 分離的明確排序）。

### Dynamic Workflow 路由規則（2026-05-30 落地 — Opus 4.8 Workflow 工具）

| 條件 | 用什麼 | 理由 |
|------|--------|------|
| **大規模 fan-out + 無資料依賴 + 終態明確**（跨 7 樣本 benchmark / 文獻 cross-check / NO-GO 前多角度 stress-test）| **Dynamic Workflow**（prompt 含 `workflow` keyword 觸發 / 存 `.claude/workflows/`）| plan 從 context 移到 script、resumable、adversarial 內建 |
| **含 Hard Gate**（C++ commit 必編譯 / 刪檔 / NO-GO 判定 / evidence_ledger 覆寫）| **維持主 agent 編排，絕不包進 workflow** | ⚠ workflow subagent 一律 acceptEdits、**繞過 exit-2 hook**（pre_commit_compile_check / kb_schema_check / pre_tier_upgrade_check）、無 mid-run 暫停 |
| **互動探索 / 需中間結果決定下一步** | 主回合 or 既有 sub-agent | workflow 無 mid-run user input |
| **單樣本 pilot / 5 分鐘小活** | 主回合直接做 | workflow 10× token，殺雞用牛刀 |
| **產數字的長 compute**（ISM C++ 全跑 / BAM 處理 / `scripts/run_*.sh` / 任何結果無法 mid-run Read-back 的重跑）| **絕不放進 workflow `agent()` step**；主回合 `Bash(run_in_background)` 或 `scripts/run_*.sh` 跑 → 落檔 → Read 驗 → **才**用 workflow fan-out 匯總「已落檔結果」 | workflow step 短壽 + acceptEdits 黑箱 + 無 mid-run Read-back → 中斷則半完成無法驗 = 撞 §13.0 捏造根因；§1 長計算須主回合可見啟動 |

> **長 job × workflow 鐵則（2026-06-03）**：workflow `agent()` step 是「**輕量驗證 / 讀檔匯總**」單元，**非長 compute 容器**。判準＝該 step 結果能否 **mid-run Read-back 驗證**（用可驗證性，非分鐘數——「10min」只是粗估）。**標準 hybrid**：① 主回合背景跑 compute 落 `.json/.tsv` → ② Read 回真值確認非 error/未完成（§13.0 step 2）→ ③ **才** workflow fan-out 匯總。⚠ ≥3 樣本「實跑腳本」用 `parallel-benchmark` agent，但須先 resource preflight（N 個 30-45min BAM 平行恐撞 CPU/mem → 預設改循序背景跑）；workflow 只接「讀已算結果」階段。實證：cross-sample ASM workflow 4 樣本曾被殺改背景補跑（memory `project_cross_sample_asm_reproducibility`）。**判斷錯誤屬路由類（非 §13 數字捏造類）→ 純文字規則即足，不加機械守衛**（除非日後出現實際 in-the-wild 事故）。

> `/effort ultracode` 不設長期預設（每任務自動編排 workflow + 略過 launch 提問 → 繞過確認矩陣 + 放大 token）；改逐案 `workflow` keyword 觸發，符合 §1「長計算需當輪明示」。

---

## §9 Agent 上下文 3 入口分工（2026-05-18 P4 E5 — 從 5 降到 3）

| 入口 | 權威範圍 | Claude Code 載入方式 |
|------|--------|--------|
| `InterSubMod/AGENTS.md` | 跨 agent governance（語言、結構、build、目標、KB 義務、研究輸出組織、回應分級 §15）| 自動 concatenate |
| `.claude/CLAUDE.md`（本檔）| Claude Code 模式特定（確認矩陣、Skills、Hooks、Rules、subagent）| 自動載入 |
| `InterSubMod/docs/CURRENT_FOCUS.md` | live 主軸、阻塞、active cycle、剛 audit / fix 紀錄 | **SessionStart hook 自動注入** ✅（commit ee648fb 已落地）|

**降級後備用入口**（不在主流程）：
- `InterSubMod/docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` — **降為 reference manual**；新內容直接寫 CURRENT_FOCUS §11 Agent Harness 狀態
- `InterSubMod/research/autoresearch/research_direction.md` — **降為 hypothesis backlog 索引**；不再作為對話啟動入口

⚠ **2026-05-18 變更理由**: Researcher agent 與 Architect agent 雙重 cross-validate 結論「Anthropic single source of truth + Walking Labs L04 'one giant instruction file fails' + OpenAI 'map not 1000-page manual'」三條業界共識指向 — 5 入口稀釋 SoT，降到 3 入口（governance / agent-specific / live state）對齊。

---

## §10 回應分級機制

詳見 `InterSubMod/AGENTS.md §15`（跨 agent governance；Claude Code 完全遵循）。
本檔不重複內容以避免 drift。

---

## §11 敘述框架預設啟用（2026-05-22 落地）

> **目的**：減少用戶理解負擔 — 所有「整理 / 報告 / 説明」AI 回覆預設套敘述框架。

### Tier 對應啟用條件（4-tier，與 SKILL.md §1 + AGENTS.md §15.2 同步）

| Tier | 條件 | 行為 |
|------|------|------|
| **Tier 1** | factual lookup / single-line answer / `<12 字 prompt` / yes-no | **skip** framework — 直接答 |
| **Tier 2** | 200-500 字 / 跨 2-3 概念 | 回覆**首行**聲明 framework（如「用 PREP：」）+ 結構化內容 |
| **Tier 3** | ≥500 字 / 跨 ≥3 概念 / 多文件統整 / 結論性報告 | 完整跑 `/narrative-frame` N1-N6 + structured output + source mapping + 業界源 footer |
| **Tier 4** | 重大 paper-scope / NO-GO 判定 / tier 升級至 ⭐4-5 | Tier 3 + 對齊 `/scientific-rigor` §2-§7（每 claim 標 L1-L5 + DAG + Pre-reg 對照表）|

### 觸發 keyword（中英）

UserPromptSubmit hook `narrative_frame_advisor.sh` 偵測：
- **中**：整理 / 報告 / 説明 / 彙報 / 總結 / 介紹 / 講解 / 解釋 / 簡報 / 教 / 寫 / 整合 / 分析 / 對比 / 比較 / pitch / 答辯
- **英**：explain / summarize / report / pitch / present / teach / walk through / integrate / outline / breakdown

### 用戶 override

- 説「不用框架」/「skip framework」 → advisor skip + AI 不套用
- 説「換 X framework」 → 直接跳 N5 重套
- 説「用 Tier 2 inline」 → 強制走 Tier 2（即使 ≥500 字）

### 與 thin wrapper skill 關係

- 既有 7 報告類 skill（`/structured-tech-report` / `/weekly-report` / `/pptx-build` 等）已 thin-wrapper 化
- 用戶説「整理 commit」 → 路由 `/structured-tech-report` thin wrapper → 預設 framework = A3+ADR+Postmortem-hybrid（即原 13 段）
- 用戶説「換 framework」 → 跳 `/narrative-frame` N1-N6 動態挑

### 50+ framework catalog

`InterSubMod/.claude/skills/narrative-frame/references/framework_catalog.md`（10 大類 SoT）
快速查表：`InterSubMod/.claude/skills/narrative-frame/references/scenario_to_framework.md`

⚠ **變動頻率上限**：新框架進 catalog → 必同步 scenario_to_framework + business_sources；每月 audit 一次 drift。

---

## §12 搜尋紀律（Search Discipline — 2026-05-30 落地）

> **背景**：本 repo 極大（build/output/data/_deps/research/figures/assets + 3 個 worktree 副本，**均已 .gitignore**）。大規模遞迴搜尋是「卡住 / 高消耗」主因。

**核心原則：判斷先 → 精準搜尋。**

| 工具 | 尊重 .gitignore？ | 用法 |
|------|:---:|------|
| **Grep / Glob 工具**（ripgrep）| ✅ **自動跳過 9 重目錄** | **廣搜首選** — 直接用，無需手動排除 |
| Bash `grep -r` / `find` / `du` | ❌ 掃全 repo 含 GB 級 data/output/worktree | **僅限 scope 到子目錄**，或加 `-maxdepth 1` |

**規則**：
1. **廣搜一律用 Grep/Glob 工具**（非 Bash grep -r）— 它們已自動排除重目錄。
2. Bash `grep -r`/`find` **必 scope 到輕目錄**（`.claude/ src include docs scripts state tests tools templates`），禁止對 repo root（`.`）遞迴。
3. `find` 淺掃加 `-maxdepth 1/2`；深掃必加 `-path ... -prune` 排除重目錄。
4. **預先排除目錄**（重，勿掃）：`build output data _deps research figures assets slides .claude/worktrees` + binary（`*.bam *.vcf.gz *.png`）。
5. 真要全掃 → 指令加 `ALLOW_FULL_SCAN` 顯式 opt-in。

**強制層**：`search_scope_guard.sh`（PreToolUse Bash，**exit 2 阻擋** `grep -r .` / `find .` 無 maxdepth / `du */`）。

---

## §13 數據誠信 — 三層防捏造（2026-06-01 落地）

> **背景**：2026-06-01 事件 — AI 在報告/HTML 把「預期數字」當真值寫入，分析其實未完成/失敗，方向還相反。postmortem §9/§9b 證實**同 session 內純文字規則失效三次** → 只有機械防線 + 結構紀律有效。完整：`InterSubMod/docs/postmortems/20260601_fabricated_metric_in_html_preview_postmortem.md`、memory `feedback_no_fabricated_numbers_in_reports`。

### 🔴 §13.0 最高優先鐵則 — 「先有驗證過的數字，才可撰寫文件」（precondition，不可協商）

> **三次捏造的精確根因（postmortem §9b 定案）**：捏造 = **「在同一個 tool-call batch 裡，同時發出『產生數字的指令』(Bash/script) 與『寫入該數字的指令』(Write/Edit 報告)」**。平行 batch 讓 Write 拿不到當批還沒回傳的數字 → AI 用記憶/預期補。**不是忘記規則，是 batch 結構問題。**

**絕對撰寫前置（任何含數字的報告 / HTML / slide / 整理）**：
1. **分析必先完整跑完** → 數字**全部落檔**（.json/.tsv/.txt）。
2. **Read 讀回真值** + **確認非 error/INCONCLUSIVE/未完成** + **確認限制與適用範圍**（樣本數、單樣本？confound？）。
3. **數字到齊且驗證過** → **才**開始撰寫 / 整理文件。
4. **物理隔離（鐵律）**：撰寫報告的 `Write`/`Edit` 與產生數字的分析（`Bash`/script）**永不放同一個 tool-call batch**。先送分析 batch → 等回傳 → Read → **下一個** batch 才 Write 報告。
5. 數字未齊/未驗證 → 該處寫 `{{待填}}` 或整段不寫，**絕不填「預期值/合理範圍」**。

**Pre-write 自我檢查（撰寫前逐項確認，全 ✓ 才動筆）**：
- [ ] 這份報告要用的每個數字，**現在**都能在某個檔案 grep 到？
- [ ] 每個數字都是**這一輪剛 Read 回來**的真值（非記憶/預期）？
- [ ] 分析腳本回傳是 success（非 error / INCONCLUSIVE / 檔案不存在）？
- [ ] 限制 / 適用範圍 / confound 已確認並會在報告誠實標註？
- [ ] 撰寫的 Write 與分析的 Bash **不在同一 batch**？

> **機械後盾**：違反此鐵則的產物會被 §13 三層攔（A 由構造使手打不可能 / B gate 抓無來源數字 / C 收尾溯源表）。但**結構紀律（§13.0）是第一線**——機械層是 backstop，不是免死金牌。

**核心問句**：報告裡每個數字，問「**這個數字現在能在哪個檔案 grep 到？**」grep 不到 = 捏造。

| 層 | 機制 | 何時 | 強弱 |
|----|------|------|------|
| **A 由構造防止**（最優先）| `scripts/fill_report.py <template> <data.json> -o <out>`：數字從 data.json 注入，缺 key 直接 refuse 不 render。含 metric 報告**一律 template+data 注入，不手打**（延續 `harness_health.py` 從不捏造的模式）| 所有資料密集 HTML/報告 | ★★★ 物理上無法捏造 |
| **B 寫入時 gate**（backstop）| `number_provenance_check.sh`（§4）：抽 metric 形數字 → grep bounded 來源。validated/pi → exit 2 擋；其他 → advisory | 手寫不可免時 | ★★ 抓手打捏造 |
| **C 任務結束溯源表**（紀錄依據）| `python3 scripts/number_provenance.py audit <report> [--sources ...]`：產「metric→檔案:行+狀態」表 | validated / PI / handoff 收尾前 | ★★ 紀錄 + 邊界攔截 |

**鐵則（與模型無關，純文字規則不夠所以有機械層）**：
1. **序列依賴**：`分析 → 寫檔(.json/.tsv) → Read 讀回 → 才 Write 報告`。**「跑分析」與「寫報告」絕不同批平行**。
2. **未完成留白**：寫 `{{待填}}` 或整段不寫，**絕不填預期值**。
3. **聲明來源**：報告 frontmatter 加 `data_sources: <path>,<path>` 讓 B 層 gate 找得到；或數字放同層 `_assets/`。
4. **override**：數字確實來自他處（如引用另一 validated 報告）→ 內文加 `<!-- provenance-verified: 來源説明 -->`。
5. **平衡（真實×消耗×錯誤率）**：三層全確定性 grep（flat-rate≈0）；B 只抓 metric 形數字（跳年份/章節號）降假陽性；衍生數字假陰性靠 A 層補。**不裝 LLM-judge 偵測**（EMNLP 2025：對數字/量詞最不準，三軸全輸）。

### 🔴 §13.7 完成宣稱通用 gate（2026-06-02 借鑑 superpowers verification-before-completion）

> **§13.0 數字鐵則是此通用鐵則的特例**：不只「數字」，**任何完成/成功宣稱**都需 fresh 驗證證據。

**鐵律：NO COMPLETION CLAIMS WITHOUT FRESH VERIFICATION EVIDENCE** — 宣稱「done / 修好了 / tests pass / build ok / 已落地 / 已驗證」前，跑 5-step gate：
1. **IDENTIFY**：什麼指令/檔案能證明此宣稱？
2. **RUN**：跑**完整**指令（fresh，非憑記憶）。
3. **READ**：讀完整輸出 + 查 exit code + 數失敗數。
4. **VERIFY**：輸出真的支持宣稱嗎？
5. **CLAIM**：確認後才宣稱。**跳任一步 = 説謊不是驗證**。

**紅旗語言（宣稱前出現 = 停下驗證）**：`should / probably / seems to / 應該 / 大概`、`Done! / Perfect! / 搞定`（未驗證的滿足感）、**信任 subagent/postmortem 回報而不獨立驗證**（本輪實例：postmortem 宣稱「P-15 陷阱卡已落地」但 grep 不到 → 驗證才發現根本沒做）。半自動對應 hook `researcher_claim_evidence_check.sh`（hedge 偵測）。
