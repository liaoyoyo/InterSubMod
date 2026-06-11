# CLAUDE.md — Claude Code 行為規範（風格 portable 模板）

> **此檔職責**：約束 Claude Code 特定行為（確認矩陣、Skills/Hooks 機制、Opus 4.7 特性）。
> **跨 agent governance（專案目標、output 規範等）→ `AGENTS.md`**。
> Claude Code 啟動時 concatenate 載入兩檔。
>
> **此為遷移模板**：保留所有「做事 / 風格」規則，刪除原 InterSubMod 專案特定內容。
> **新機 setup 步驟**：拷貝後依 §3 §4 §5 §6 替換為新專案實際 skill/hook/rule 清單；§0 §1 §2 §7 §8 §11 通常直接保留即可。

---

## §0 任務類型 gate

**規則**: 任何 conversation 啟動時 / 接到新任務時，**先做 task type classification** 再進 §1 暫停判斷矩陣。task type 決定 default scope 與 deliverable 義務。

**6 類速查**:

| Type | Keyword | Default Scope |
|------|---------|--------------|
| A pilot | pilot / 探索 / 試 / 先看 | Subset OK |
| **B validation** | 完整 / final / validation / 驗證 | **全範圍 + 全樣本** |
| **C production** | release / deploy / merge | **全 scope + cross-platform** |
| **D handoff** | handoff / 交付 / external | **全 scope + 説明文件** |
| E hotfix | bugfix / hotfix / urgent | 最小可重現 |
| F demo | demo / 示意 / 教學 | Subset OK + DEMO 標 |

**強制行為**:
1. **啟動先分類** — 不可省略；模糊 → AskUserQuestion 必問 task type
2. **AskUserQuestion 強制含 scope 維度** — 觀察/報告類任務問卷必含「scope 全 vs subset？」
3. **subset 必標 partial flag** — 文件 ribbon 標記
4. **「見樹也見林」四層** — aggregate / canonical / extreme outlier / well-explained

**主動 recall**: 啟動 prompt 含 keyword（完整 / 報告 / 整理 / handoff / 驗證 / final / 對外）→ 必 recall `feedback_observation_scope_default_comprehensive` 與 `feedback_task_first_then_doc_then_plan` 兩條 memory。

**業界對照**: Cynefin domain + DACI/DECIDE + Bland Assumption Map。

---

## §1 確認協議與暫停判斷

**模式觸發詞**: `互動模式`（預設） / `全自動` / `auto`
**級別**: Hard Gate（必停）> Gate（互動停/自動過）> Review（互動停/自動過）> FYI

**Hard Gate 必停場景**（不可逆操作永遠 🔴）：
1. **刪除檔案**
2. **編譯型語言 commit**（pre-commit hook 強制編譯檢查；若新專案是 C++/Rust/Go 等保留此項）
3. **重大方向 NO-GO 判定**

**暫停判斷矩陣**（影響度 × 信心度）：
| 影響 ＼ 信心 | 高 | 中 | 低 |
|---|---|---|---|
| **低**（可逆/<10min） | 🟢 一行告知 | 🟡 列假設 | 🟠 節點暫停 |
| **中**（10min-1h） | 🟡 列假設 | 🟠 節點暫停 | 🔴 立即暫停 |
| **高**（>1h/影響結論） | 🟠 節點暫停 | 🔴 立即暫停 | 🔴 立即暫停 |

**可逆性 override**: 刪檔／commit／NO-GO／覆寫持久 ledger 永遠 🔴。

**長計算例外**: 用戶**當輪明示**啟動長計算（「跑全量」「平行 N 樣本」）→ 不二次確認，直接執行 + 一行告知。模型自判需長計算 → 仍 🔴。

**一行告知格式**: `[決策]（影響: 低/中/高, 信心: 低/中/高, 理由: 一句）`

完整維度定義 → `/confirmation-protocol` skill。

---

## §2 Opus 4.7 模型特性備註

> **Opus 4.7 literal 特性**：模型不會推斷未明講需求、不會悄悄泛化指令。模糊輸入會被按字面執行 — **假設陳述是唯一屏障**。

- **預設不 spawn subagent**（單回合優先）
- **預設少 tool calls**（reasoning 解決優先於反覆讀檔）
- **回應長度動態化**（不再固定冗長）
- **首 turn 規格完整度檢核**：規格缺項 ≥2 且高影響 → 必須回問

完整模型執行特性 → `.claude/rules/opus47-behavior.md`（path-scoped 條件式載入）

---

## §3 Skills 分類索引（風格 portable 17 個 + 新增空槽）

> ⚠ **新機調整**：以下為「風格 portable」17 個 baseline；新專案內加入 domain-specific skill 後請更新此索引。

- **元方法論（7）**: `/confirmation-protocol` `/scientific-rigor` `/fast-learning-coach` `/problem-framing-ideation` `/pre-decision-audit` `/known-pitfalls` `/infra-ops`
- **敘述框架（1）**: `/narrative-frame`（主入口 + 50+ catalog；取代固定範本）
- **報告生成 retrospective（7）**: `/weekly-report` `/structured-tech-report` `/results-report` `/report` `/myPPT` `/pptx-build` `/html-report-build`
- **視覺化（2）**: `/image-gen` `/image-vision-check`
- **文件管理（4）**: `/doc-standards` `/memory-consolidation` `/data-audit` `/citation-verification`
- **域 specific（待新專案填入）**: `<TODO: 依新任務性質填入 domain skills>`

---

## §4 Hooks 概覽（風格 portable 14 個）

依 `.claude/settings.local.json` 完整定義（6 events: SessionStart / UserPromptSubmit / PreToolUse / PostToolUse / SubagentStop / Stop）

**Hard Gate hooks**（不可繞過 — `exit 2` 阻擋）:
- `pre_commit_compile_check.sh`（編譯型語言 commit 必編譯；新專案非 C++ 可註解）
- `no_binary_commit.sh`（binary 阻擋）

**通用 advisory hooks**:
- `session_start_inject_focus.sh` — SessionStart 注入 live working state（讀新專案的 CURRENT_FOCUS.md 等價檔）
- `md_path_format_rule.sh` — UserPromptSubmit 路徑前綴規則
- `narrative_frame_advisor.sh` — UserPromptSubmit keyword 偵測 → 推薦 framework
- `task_type_advisor.sh` — UserPromptSubmit 6 類 task type keyword 偵測 + 注入 advisory
- `verify_gate.sh` — PreToolUse Edit|Write Default-FAIL evidence gate
- `external_input_sanitizer.sh` — PostToolUse WebFetch injection 偵測
- `subagent_completion_logger.sh` — SubagentStop cost/cache 紀錄
- `evidence_read_tracker.sh` `memory_recall_logger.sh` `skill_change_audit.sh` `skill_usage_logger.sh` `cache_telemetry.sh` `allow_list_audit.sh`（telemetry/audit 類）

---

## §5 Rules 載入狀態（path-scoped）

| Rule 檔案 | 狀態 | globs 條件 |
|----------|------|-----------|
| `.claude/rules/opus47-behavior.md` | **條件式載入** | `.claude/skills/**/SKILL.md`, `.claude/skills/**/*.json`, `.claude/rules/**/*.md` |
| `<TODO 新專案 build rule>` | TBD | TBD |
| `<TODO 新專案 workflow rule>` | TBD | TBD |

**Hard Gate 保留**: 即使 rule 變條件式（編輯對應 path 才載入），編譯型語言 commit Hard Gate（§1 + pre-commit hook）**仍適用** — Hard Gate 不依賴 rule 載入，而是不可逆操作的絕對規則。

---

## §6 Working State Pointer（動態狀態）

**當前主軸**: 由 SessionStart hook 從 `<new-project>/docs/CURRENT_FOCUS.md`（或新專案命名的等價檔）自動注入 ~3000 chars / ~500 tokens
**Pending Items**: 見 MEMORY.md 索引
**詳細上下文載入**: 觸發 `/research-context-loader` 或新專案的 context loader skill

⚠ **何時不需詳細上下文**: 純 code edits、單檔 doc 寫作、簡單問答
⚠ **變動頻率上限：週級** — CURRENT_FOCUS.md > 7 天未更新 → hook 主動提醒

> **業界對照**: `docs/CURRENT_FOCUS.md` 概念等同 [Cline Memory Bank `activeContext.md`](https://docs.cline.bot/features/memory-bank)（live working state，最高變動層）。

---

## §7 Context 壓縮保留指令（Claude Code `/compact`）

`/compact` 必保留:
- 架構決策 + 理由 / 未解決 blockers / 涉及檔案路徑
- 用戶限制條件 / 假說 ID + 驗證層級 / 未完成步驟
- 本檔 §0-§11 規則 + `AGENTS.md` 跨 agent governance

**可安全壓縮**: 中間計算 / 已通過完整測試輸出 / 冗長工具呼叫結果。

---

## §8 Opus 4.7 Subagent 觸發

**預設不 spawn**（單回合優先）。**明示觸發語**才啟動：
- 跨樣本平行 benchmark → `parallel-benchmark`
- 大範圍探索（>3 query）→ `Explore` / `general-purpose`
- PR 審查 → `pr-review-toolkit:code-reviewer`
- 夜間/離線研究 → `headless-research`

**任務切割原則**：可切割且需要大量 context 處理 → 啟動 agent + **必須清楚回報主 agent** + **科學工程紀錄成文件**供檢核驗證。

---

## §9 Agent 上下文 3 入口分工

| 入口 | 權威範圍 | Claude Code 載入方式 |
|------|--------|--------|
| `<new-project>/AGENTS.md` | 跨 agent governance（語言、結構、build、目標、output 規範、回應分級）| 自動 concatenate |
| `.claude/CLAUDE.md`（本檔）| Claude Code 模式特定（確認矩陣、Skills、Hooks、Rules、subagent）| 自動載入 |
| `<new-project>/docs/CURRENT_FOCUS.md` | live 主軸、阻塞、active task | **SessionStart hook 自動注入** |

⚠ **設計理由**: 三入口 = governance / agent-specific / live state，對齊 Anthropic single source of truth + OpenAI 'map not 1000-page manual' 共識。

---

## §10 回應分級機制

詳見 `AGENTS.md §15`（跨 agent governance；Claude Code 完全遵循）。
本檔不重複內容以避免 drift。

---

## §11 敘述框架預設啟用

> **目的**：減少用戶理解負擔 — 所有「整理 / 報告 / 説明」AI 回覆預設套敘述框架。

### Tier 對應啟用條件（4-tier）

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

- 報告類 skill（`/structured-tech-report` / `/weekly-report` / `/pptx-build` 等）已 thin-wrapper 化
- 用戶説「整理 commit」 → 路由 `/structured-tech-report` thin wrapper → 預設 framework = A3+ADR+Postmortem-hybrid（即原 13 段）
- 用戶説「換 framework」 → 跳 `/narrative-frame` N1-N6 動態挑

### 50+ framework catalog

`.claude/skills/narrative-frame/references/framework_catalog.md`（10 大類 SoT）
快速查表：`.claude/skills/narrative-frame/references/scenario_to_framework.md`

⚠ **變動頻率上限**：新框架進 catalog → 必同步 scenario_to_framework + business_sources；每月 audit 一次 drift。
