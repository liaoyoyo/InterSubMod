# CLAUDE.md — Claude Code 行為規範

> **此檔職責**：僅約束 Claude Code 特定行為（確認矩陣、Skills/Hooks 機制、Opus 4.8 特性）。
> **跨 agent governance（專案目標、研究主軸、Step→Verify、KB 義務、Output 規範等）→ `InterSubMod/AGENTS.md`**。
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

## §3 Skills 分類索引（46 個 SKILL.md，Claude Code 特定 — 2026-05-28 drift 再校正）

> **drift 修正紀錄**：
> - 2026-05-20：原寫 45 → 一度改 44（誤判 `grill-me/` 為 phantom +1）。
> - **2026-05-28 實測校正**：`find .claude/skills -name SKILL.md | wc -l` = **46**（先前一度誤記 44）；46 = 原 45（含 deprecated 的 `/html-preview` 仍計數）+ 新增 `/pipeline-manifest`（社群 gap G1）。校正機制：新增 `creation_guard.sh` hook 在 Write 新 SKILL.md 時提醒同步本計數。
> - **2026-05-30 補正**：`grill-me` 並非「完全不存在」— 它是一個 **dangling symlink**（`.claude/skills/grill-me -> ../../.agents/skills/grill-me`，target 已不存在）。因無 `SKILL.md`，**不計入 46**（故 `find -maxdepth 1 -type d` = 47 但 `find -name SKILL.md` = 46，差 1 即此 symlink）。Claude Code 按 SKILL.md 列舉故不受影響。物理清除屬 Hard Gate 刪檔（待 ack `git rm .claude/skills/grill-me`）。
> - 新增 12 個未分類 skills 進對應類別（2026-05-20）。
> **重複交叉位置**：`/feature-layered-observation`（P3 + 研究專用）/ `/multi-sample-consistency`（P4 + 研究專用）/ `/pre-decision-audit`（元方法論 + 三層樓 pre）/ `/run-evaluator`（P5 + 三層樓 post）— 在多分類列出表示同 skill 多角色。

- **元方法論（9）**: `/confirmation-protocol` `/known-pitfalls` `/cycle-state` `/research-context-loader` `/fast-learning-coach` `/scientific-rigor` `/pre-decision-audit` ⭐ 新 `/problem-framing-ideation` `/provenance-tier-audit`
- **7-Phase Waterfall（7）**: P0 `/cycle-init` → P1 `/research-loop` → P2 `/check-staleness` → P3 `/feature-layered-observation` → P4 `/multi-sample-consistency` → P5 `/run-evaluator` → P6 `/conclude-research`
- **程式開發（4）**: `/cpp-change` `/methodology-audit` `/infra-ops` `/verification-loop`
- **文件管理（5）**: `/doc-standards` `/data-audit` `/memory-consolidation` `/citation-verification` `/pipeline-manifest` ⭐ 新（reproducibility provenance DAG；與 data-audit 分工：data-audit 查組織、pipeline-manifest 查 script→figure 因果鏈）
- **報告生成 retrospective（7）**: `/weekly-report` → `/pptx-build` / `/html-report-build` / `/results-report` / `/structured-tech-report` / `/report` / `/myPPT`
- **視覺化（4）**: `/html-preview` `/image-gen` `/image-vision-check` `/research-dashboard`
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

依 `InterSubMod/.claude/settings.local.json` 完整定義（含 SessionStart / UserPromptSubmit / PreToolUse / PostToolUse / SubagentStop / Stop / **PreCompact** 7 個 events；**38 hook scripts**，2026-05-30 實測 `ls scripts/hooks/*.sh | wc -l`）。

**2026-05-28~30 新增 4 個 hook**（社群 gap + 搜尋紀律）:
- `search_scope_guard.sh`（PreToolUse Bash — **exit 2 阻擋** `grep -r .` / `find .` 無 maxdepth / `du */` 等不尊重 .gitignore 的 repo-root 遞迴搜尋；§12 搜尋紀律強制層）
- `precompact_autosave.sh`（PreCompact 事件 — 壓縮前自動 dump active cycle + CURRENT_FOCUS 快照到 `state/compact_snapshots/`；對齊 §7 保留指令；gap G6）
- `creation_guard.sh`（PreToolUse Write — 新建 SKILL.md/agent 前 dedup + 計數同步提醒；防 §3 計數 drift；gap G5）
- `skill_registry_sync.sh`（PostToolUse Edit|Write — 編輯 README/CLAUDE.md 時比對「磁碟實際 skill/agent 計數 vs 文件宣稱最大數」，雙向 drift 守衛；2026-05-30 gap G5 補完）

**Hard Gate hooks**（不可繞過 — `exit 2` 阻擋）:
- `pre_commit_compile_check.sh`（C++ commit 必編譯）
- `kb_schema_check.sh`（KB 寫入前 schema 檢核）
- `pipeline_block_check.sh`（長 pipeline 磁碟檢核）
- `no_binary_commit.sh`（commit binary 阻擋）
- `kb_sot_guard.sh`（F1 SoT 數字保護）

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
| ~~`.claude/rules/opus47-behavior.md`~~ | **DEPRECATED**（已移除 globs，不再 auto-load；待 Hard Gate `git rm`）| — |
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

### Dynamic Workflow 路由規則（2026-05-30 落地 — Opus 4.8 Workflow 工具）

| 條件 | 用什麼 | 理由 |
|------|--------|------|
| **大規模 fan-out + 無資料依賴 + 終態明確**（跨 7 樣本 benchmark / 文獻 cross-check / NO-GO 前多角度 stress-test）| **Dynamic Workflow**（prompt 含 `workflow` keyword 觸發 / 存 `.claude/workflows/`）| plan 從 context 移到 script、resumable、adversarial 內建 |
| **含 Hard Gate**（C++ commit 必編譯 / 刪檔 / NO-GO 判定 / evidence_ledger 覆寫）| **維持主 agent 編排，絕不包進 workflow** | ⚠ workflow subagent 一律 acceptEdits、**繞過 exit-2 hook**（pre_commit_compile_check / kb_sot_guard）、無 mid-run 暫停 |
| **互動探索 / 需中間結果決定下一步** | 主回合 or 既有 sub-agent | workflow 無 mid-run user input |
| **單樣本 pilot / 5 分鐘小活** | 主回合直接做 | workflow 10× token，殺雞用牛刀 |

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
