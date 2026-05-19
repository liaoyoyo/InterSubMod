# CLAUDE.md — Claude Code 行為規範

> **此檔職責**：僅約束 Claude Code 特定行為（確認矩陣、Skills/Hooks 機制、Opus 4.7 特性）。
> **跨 agent governance（專案目標、研究主軸、Step→Verify、KB 義務、Output 規範等）→ `InterSubMod/AGENTS.md`**。
> Claude Code 啟動時 concatenate 載入兩檔。

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

## §2 Opus 4.7 模型特性備註

> **Opus 4.7 literal 特性**：模型不會推斷未明講需求、不會悄悄泛化指令。模糊輸入會被按字面執行 — **假設陳述是唯一屏障**。

- **預設不 spawn subagent**（單回合優先）
- **預設少 tool calls**（reasoning 解決優先於反覆讀檔）
- **回應長度動態化**（不再固定冗長）
- **首 turn 規格完整度檢核**：規格缺項 ≥2 且高影響 → 必須回問

完整 Step→Verify 格式 / 假設陳述規則 / 證據敘述 → **`InterSubMod/AGENTS.md` §6**。
完整模型執行特性 → `.claude/rules/opus47-behavior.md`（**目前永遠載入**；frontmatter 加入後變條件式 — 詳見 §5）

---

## §3 Skills 分類索引（45 個 SKILL.md，Claude Code 特定 — 2026-05-19 新增 /pre-decision-audit）

- **元方法論**: `/confirmation-protocol` `/known-pitfalls` `/cycle-state` `/grill-me` `/research-context-loader` `/fast-learning-coach` `/scientific-rigor` `/pre-decision-audit` ⭐ 新
- **7-Phase Waterfall**: P0 `/cycle-init` → P1 `/research-loop` → P2 `/check-staleness` → P3 `/feature-layered-observation` → P4 `/multi-sample-consistency` → P5 `/run-evaluator` → P6 `/conclude-research`
- **程式開發**: `/cpp-change` `/methodology-audit` `/infra-ops` `/verification-loop`
- **文件管理**: `/doc-standards` `/data-audit` `/memory-consolidation` `/citation-verification`
- **報告生成（retrospective）**: `/weekly-report` → `/pptx-build` / `/html-report-build` / `/results-report` / `/structured-tech-report` / `/report`
- **假說驗證三層樓（pre → process → post）⭐ 新 (2026-05-19)**:
  - **Pre** (entry-point ≤30min): `/pre-decision-audit` (7 outputs + Cynefin gate + 5-dim credibility + GO/PROBE/NO-GO verdict)
  - **Process** (spec 實作中 live): `/implementation-notes` (4 sections + Lore — 設計決定 / 偏離 / 折衷 / 未決)
  - **Post** (P5 cycle 結束後): `/run-evaluator` (tier ⭐1-5 + 6 risk components)
- **研究專用**: `/auc-confound-guard` `/feature-layered-observation` `/multi-sample-consistency` `/pivot-direction` `/inject-hypothesis`

---

## §4 Hooks 概覽（Claude Code 特定 — 2026-05-18 P4 完整收尾）

依 `InterSubMod/.claude/settings.local.json` 完整定義（含 SessionStart / UserPromptSubmit / PreToolUse / PostToolUse / SubagentStop / Stop 6 個 events；**30 hook scripts** 跨 22 matchers）。

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

---

## §5 Rules 載入狀態（path-scoped 已落地 2026-05-18）

| Rule 檔案 | 狀態 | globs 條件 |
|----------|------|-----------|
| `.claude/rules/cpp-build.md` | **條件式載入** ✅ | `src/**/*.cpp`, `src/**/*.hpp`, `include/**/*.hpp`, `include/**/*.h`, `tests/**/*.cpp`, `CMakeLists.txt` |
| `.claude/rules/opus47-behavior.md` | **條件式載入** ✅ | `.claude/skills/**/SKILL.md`, `.claude/skills/**/*.json`, `.claude/rules/**/*.md`, `.claude/skills/**/playbook.md`, `.claude/skills/**/prompts/*.md` |
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

## §8 Opus 4.7 Subagent 觸發

**預設不 spawn**（單回合優先）。**明示觸發語**才啟動：
- 跨樣本平行 benchmark → `parallel-benchmark`
- 大範圍探索（>3 query）→ `Explore`
- PR 審查 → `pr-review-toolkit`
- 夜間/離線研究 → `headless-research`

**任務切割原則**：可切割且需要大量 context 處理 → 啟動 agent + **必須清楚回報主 agent** + **科學工程紀錄成文件**供檢核驗證。

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
