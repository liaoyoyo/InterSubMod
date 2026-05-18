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

## §3 Skills 分類索引（44 個，Claude Code 特定）

- **元方法論**: `/confirmation-protocol` `/known-pitfalls` `/cycle-state` `/grill-me` `/research-context-loader` `/fast-learning-coach` `/scientific-rigor`
- **7-Phase Waterfall**: P0 `/cycle-init` → P1 `/research-loop` → P2 `/check-staleness` → P3 `/feature-layered-observation` → P4 `/multi-sample-consistency` → P5 `/run-evaluator` → P6 `/conclude-research`
- **程式開發**: `/cpp-change` `/methodology-audit` `/infra-ops` `/verification-loop`
- **文件管理**: `/doc-standards` `/data-audit` `/memory-consolidation` `/citation-verification`
- **報告生成**: `/weekly-report` → `/pptx-build` / `/html-report-build` / `/results-report` / `/structured-tech-report`
- **研究專用**: `/auc-confound-guard` `/feature-layered-observation` `/multi-sample-consistency` `/pivot-direction` `/inject-hypothesis`

---

## §4 Hooks 概覽（Claude Code 特定）

依 `InterSubMod/.claude/settings.local.json` 完整定義（含 UserPromptSubmit / PreToolUse / PostToolUse / SubagentStop / Stop 多個 matcher 與 command；實際 hook script 數量 20+，每個 matcher 可對應多個 command）。

**Hard Gate hooks**（不可繞過）:
- `pre_commit_compile_check.sh`（C++ commit 必編譯，`exit 2` 阻擋）

**規劃新增**（Phase 5 部署）:
- `session_start_inject_focus.sh`（SessionStart 注入 CURRENT_FOCUS）
- `md_path_format_rule.sh`（強化 `InterSubMod/...` 前綴檢核）

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

## §6 Working State Pointer（動態狀態）

**當前主軸**: 由 SessionStart hook 從 `InterSubMod/docs/CURRENT_FOCUS.md` 注入（≤ 500 tokens）
**Pending Items**: 見 MEMORY.md 索引
**研究方向變更**: 觸發 `/pivot-direction`
**詳細上下文載入**: 觸發 `/research-context-loader`（Tier 1 / Tier 2 / Tier 3 深度）

⚠ **SessionStart hook 狀態**: **待建立**（Phase 5）
⚠ **hook 建立前的暫行做法**: AI 必須主動執行 `Read InterSubMod/docs/CURRENT_FOCUS.md`
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

## §9 Agent 上下文 5 入口分工（Claude Code 角度）

| 入口 | 權威範圍 | Claude Code 載入方式 |
|------|--------|--------|
| `InterSubMod/AGENTS.md` | 跨 agent governance（語言、結構、build、目標、KB 義務、研究輸出組織）| 自動 concatenate |
| `.claude/CLAUDE.md`（本檔）| Claude Code 模式特定（確認矩陣、Skills、Hooks、Rules、subagent）| 自動載入 |
| `InterSubMod/docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` | 研究壓縮上下文、重要數據、任務順序、待決策矩陣 | 由 AI 主動 Read（觸發：對話開始 + 跨 session 交接）|
| `InterSubMod/docs/CURRENT_FOCUS.md` | live 主軸、阻塞 | SessionStart hook 注入（待建立）+ AI 主動 Read |
| `InterSubMod/research/autoresearch/research_direction.md` | AutoResearch 候選 queue | **僅候選**，不作自動執行觸發 |

---

## §10 回應分級機制

詳見 `InterSubMod/AGENTS.md §15`（跨 agent governance；Claude Code 完全遵循）。
本檔不重複內容以避免 drift。
