---
title: claude-mem / claude-md-management 深挖借鑑 + 互動 dashboard 準則 + 參考 HTML genre 裁決
date: 2026-06-03
type: governance / borrow-analysis
status: active
evidence: workflow wf_b676a2c8-b48（5 agents：claude-mem web + claude-md web + master_draft 分析 + 綜合 + critique）；critique 親驗磁碟事實（39 hook/110 memory/18 agent/46 skill）；主 agent 親讀 concept board + design_principles:255 drift
data_sources:
  - .claude/skills/html-report-build/references/design_principles.md
  - state/focus_board.html
  - scripts/hooks/skill_registry_sync.sh
---

# 借鑑深挖 — 他們有何我們缺的更好細節？

> **一句話**：claude-mem 與 claude-md-management **整套都不裝**（會造成競爭 memory authority / 把治理散文誤判 bloat / async LLM 壓縮是新捏造面）。但**確實有真缺口值得學**——以「擴充既有 hook/skill」方式借概念。critique verdict = `needs_trim`，已套用：真做 2-3 個零成本項，砍掉 component library（N=1 YAGNI）。

## 1. claude-mem（thedotmack）— 它做得比我們好的
架構：自動 hook 捕捉 → SQLite+Chroma 雙存 → async worker AI 壓縮 → 2-stage 檢索（search ID → hydrate）。
| 它更好的 | 我們現況 | 借鑑 verdict |
|---------|---------|------------|
| **自動 async 捕捉 tool-trail**（PostToolUse 立即寫 queue）| **零 capture path**，全靠 AI 記得手寫 | 🟡 pilot（session_trail hook，見下；**但 critique 警告：這不修 fabricated-numbers/P-15，那是 claim-time 已由 §13.7 蓋**）|
| 2-stage 漸進檢索（ID→hydrate）| SessionStart 固定注入 slab，無 query 過濾 | ⚪ skip（我們已有 file 兩層 MEMORY.md↔memory/；110 檔 Grep 夠）|
| watermark idempotent 一致性 | MEMORY.md↔memory/ drift 只在手跑修 | ✅ **借**（升 harness_health 燈#7）|
| 向量語意檢索（Chroma）| Grep/keyword | ⚪ skip（110 檔不值 daemon+第二 store）|
- **整套裝？❌** worker+SQLite+Chroma = 競爭 memory authority + 非 git-tracked 破 provenance + async LLM 壓縮繞過 §13。

## 2. claude-md-management（官方 Anthropic plugin）— 它做得比我們好的
架構：skill `claude-md-improver`（6 維 rubric 評分 + report-first）+ `/revise-claude-md`（session-end 反射）。純 prompt-as-code，propose→approve。
| 它更好的 | 我們現況 | 借鑑 verdict |
|---------|---------|------------|
| **Currency/Red-Flag：內文宣稱 vs 磁碟**（deleted-file refs / 失效指令）| §4 手打「39 hook/18 agent」，skill_registry_sync **只查 skill+agent 不查 hook** | ✅ **借**（hook-count + 路徑存在性偵測，~5 行）|
| 6 維內容品質 rubric + A-F grade | harness_health 6 燈全是機械一致性，零個評散文品質 | 🟡 later（改寫成治理型變體當燈#8）|
| anti-bloat（What-NOT-to-Add + Conciseness 扣分）| CLAUDE.md 365 行持續加 §，零 bloat 偵測 | 🟡 later（編輯前自查段，人工）|
| report-first gate（改前必先出 report）| 直接 Edit CLAUDE.md | 🟡 later（沿用 confirmation-protocol）|
- **整套裝？❌** 其 rubric 把 CLAUDE.md 當 codebase-context → 把 §0-§13 治理散文判 verbose/bloat 提刪、Commands=0 給誤導 D/F、Phase5 自動 Apply 繞過 confirmation-protocol。

## 3. 參考 HTML genre 裁決（master_draft + concept board）
| 檔 | 作 dashboard 參考 | 作報告/概念參考 |
|----|:---:|:---:|
| `master_draft.standalone.html`（週報）| ❌ **不該學**（print-first 長敘述、sticky TOC、多 table → 對快掃 dashboard 是噪音）| ✅ **報告 genre 金標範本**（Assertion-Evidence/5秒測試/CRAP/誠實邊界/provenance footer 全達標）|
| `20260603_目標概念整理…雙軸.standalone.html`（概念討論板）| 🟡 部分（Anthropic clay 配色 + provenance footer 可借）| ✅ **概念/討論板 exemplar**（token 系統 + SVG 概念圖 + checkbox 討論 + §13 provenance）|
| `state/focus_board.html`（焦點板）| ✅ **dashboard 的 canonical reference impl** | ❌（live 狀態非敘述）|
> **關鍵 genre 紀律（critique）**：報告 vs dashboard 是互補集，**各自登記各自 genre 的 reference，勿跨借錯方向**。focus_board 不該學週報的 TOC/多 table（它 L0-L3 嵌套已是導航）；週報不該學 since-diff（async 報告無 temporal）。

## 4. 互動 dashboard 設計準則「如何建立」— 關鍵更正
**Rule 14 已存在**（design_principles.md，2026-06-03 加，8 子條 + Shneiderman mantra + dashboard 5 秒 pre-publish + focus_board mapping）。問題不是「要不要建」，是「讓它從散文變可套」。
- ✅ **已修**：Rule 14.7 行 255 stale drift（「(規劃中)since-diff」→ 實際 focus_board sincebar 已實作）。
- ❌ **CUT（critique YAGNI）**：抽 4 個 dashboard component（collapsible/stepper/since-diff/density）進 components/ —— **N=1 dashboard 無第二消費者，不值抽象**。既有 11 component 是靠 report/slide 重複使用掙得位置；dashboard 只有 focus_board 一個。**第 2 個 dashboard 真的出現才抽**。
- ✅ focus_board.html = Rule 14 canonical reference impl（SKILL.md 指向它，不 copy 避免 drift）。

## 5. Critique-收斂行動清單（restraint-confirmed）
**真做（零~低成本）**：
1. ✅ **已做**：修 Rule 14.7 doc drift（1 行真相）。
2. 🟢 **建議做**：harness_health 燈#7 = memory index↔files drift（升 /memory-consolidation Step-4 手動 grep 為持續；read-only 零 LLM）。
3. 🟢 **建議做**：skill_registry_sync / harness_health 燈#1 擴 hook-count ground-truth（`ls scripts/hooks/*.sh`）+ CLAUDE.md 明確 .sh/.py 路徑存在性 grep（~5 行；preventive 防 phantom-count）。
4. 📋 登記：master_draft + concept board → `/report`、`/weekly-report` examples 當 genre exemplar。

**砍/緩**：
- ❌ CUT dashboard component library（N=1 YAGNI）。
- ⚪ memory_query.py = flat skip（110 檔 Grep 夠；0-hit 指標未量測，defer trigger 是 theater）。
- 🟡 session_trail capture hook = **降為審慎 pilot**（不是 fabrication 修法——那是 §13.7 claim-time gate 的事；它是 audit infra，價值真但 speculative）。若做：capture-only 最小版 + 觀察一週 + distill 必經人批准。

> 對齊 [[feedback_harness_restraint_over_adoption]]：整體 0 外部工具裝；真 ROI = 修 1 處 drift + 2 個零成本偵測燈，非裝 memory/CLAUDE.md 維護器。
