---
title: 外部工具/skill/plugin 採用評估（restraint-first）— 29 項清單裁決
date: 2026-06-03
type: governance / tool-adoption-decision
status: active
evidence: L1（grep 親查 domain + 既裝清單）+ 對齊 [[feedback_harness_restraint_over_adoption]]
data_sources:
  - .claude/skills/
  - docs/references/migration/COMMUNITY_COMPARISON.md
---

# 外部工具採用評估 — 29 項清單裁決

> **一句話結論**：**現在採用 0 項**。清單約 70% 是 web/UI/mobile/video/browser 工具 = 對 C++/Python 癌症基因體 CLI 專案硬 domain 錯位；其餘多數**已裝**或**會與既有系統衝突**。唯一值得「之後論文期再看」= `pdf`（文獻 PDF）。此裁決與 [[feedback_harness_restraint_over_adoption]] 的 21-agent 研究（12 high-fit → 4 存活，且皆是流程紀律借鑑非裝工具）完全一致。

> 🔴 **2026-06-03 修正（用戶正確反駁 + §13.7）**：本表原把 `frontend-design`/`web-design-guidelines`/`ui-skills` 一律判 **SKIP-domain**，**過於武斷**。儀表板 + PI report + weekly + pptx **全是 HTML 產出層** = 真實 web/UI domain，design/UI/JS 對「**這一層**」有用（純 web-app/mobile/browser/video 工具才是真 mismatch）。**grep 親查**：既有 `design_principles.md`(13 條) 為 print/figure 導向，**互動 web-UI patterns 0 覆蓋**（progressive disclosure / collapsible / 互動 / 資訊密度 / responsive / scannability / affordance 全缺）= 真實缺口。**restraint 解（已落地）**：擴充既有設計系統 **Rule 14（Interactive Dashboard / Web-UI Patterns）** + 套用到 focus_board（collapsible / 密度切換 / since-last-view 進退 diff），**而非裝外部 skill**（保持單一 curated SoT + 零依賴）。→ 修正後仍 **0 外部 skill 採用**，但「能力」已補。`frontend-design`(anthropics) 改判 🟡 PROBE-low（若日後要更廣 web-component 系統可當參考；目前 Rule 14 已覆蓋 dashboard 需求）。

## 決定性事實（grep 親查 2026-06-03）
- 專案語言：**68 cpp/hpp 檔**；**ZERO** `.ts` / `package.json` / `.tsx` / `.jsx` / 前端 / mobile / video。→ 任何 web/UI/mobile/TS/video/browser 工具 = domain 錯位。
- 既裝（system context + grep 確認）：**superpowers**（含 brainstorming）、**context7 MCP**、**knowledge MCP**、**pr-review-toolkit**（code-reviewer agent）、**code-simplifier**、**feature-dev** agents、**security-reviewer** agent（2026-05-30）、46 skills、18 agents、39 hooks。

## 裁決表

| # | 工具 | 裁決 | 理由 |
|---|------|------|------|
| 1 | find-skills | 🟡 PROBE-low | skill 發現；harness 已用 CLAUDE.md §3 索引 + research-orchestrator 路由。skill 數爆增才需 |
| 2 | frontend-design (anthropics) | ❌ SKIP-domain | 無前端 |
| 3 | web-design-guidelines | ❌ SKIP-domain | 無 web |
| 4 | brainstorming (superpowers) | ✅ ALREADY | superpowers 已裝 + `/problem-framing-ideation` |
| 5 | remotion-best-practices | ❌ SKIP-domain | 影片/React 渲染 |
| 6 | agent-browser | ❌ SKIP-domain | 瀏覽器自動化；無 web 目標（IGV 是桌面非瀏覽器）|
| 7 | browser-use | ❌ SKIP-domain | 同上 |
| 8 | sleek-design-mobile-apps | ❌ SKIP-domain | mobile |
| 9 | ibelick/ui-skills | ❌ SKIP-domain | UI |
| 10 | **pdf (anthropics)** | 🟡 **PROBE-med（唯一較相關）** | 文獻 PDF 抽取（論文期）；但已有 researcher + knowledge MCP + Zotero；G6 在 P4 尚非寫作期 → 列 **P2 paper-phase 採購單，不提前** |
| 11 | skill-creator (anthropics) | 🟡 PROBE-low | 已有 `SKILL_FRONTMATTER_SPEC.md` + `creation_guard.sh`；skill 撰寫頻繁才需 |
| 12 | code-review-expert | ✅ ALREADY | pr-review-toolkit + reviewer agent + `/cpp-change` |
| 13 | Superpowers plugin | ✅ ALREADY | 已裝（SessionStart 確認）|
| 14 | claude-mem | 🔴 SKIP-CONFLICT | **與既有檔案式記憶系統衝突**（MEMORY.md + memory/ + memory_recall_logger + `/memory-consolidation`）；雙記憶系統 = drift |
| 15 | security-guidance | ✅ ALREADY | security-reviewer agent（2026-05-30 新增）|
| 16 | code-review | ✅ ALREADY | pr-review-toolkit |
| 17 | pr-review-toolkit | ✅ ALREADY | 已裝 |
| 18 | commit-commands | ⚪ SKIP-low | release agent + git 治理 + commit 紀律已覆蓋；加值低 |
| 19 | code-simplifier | ✅ ALREADY | 已裝 |
| 20 | claude-md-management | 🔴 SKIP-CONFLICT | CLAUDE.md 為**重度手工治理**（§0-§13 + skill_registry_sync 漂移守衛）；auto-maintainer 會與精心結構打架 |
| 21 | explanatory-output-style | ⚪ SKIP-personal | 全域改輸出；用戶為專家研究者非學習者；要的話 `/output-style` 即可，非 harness 改進 |
| 22 | learning-output-style | ⚪ SKIP-personal | 同上 |
| 23 | frontend-design (重複) | ❌ SKIP-domain | 同 #2 |
| 24 | claude-code-setup | ⚪ SKIP | 專案已建立，非新 scaffolding 需求；研究專案用 `/init-research` |
| 25 | feature-dev | ✅ ALREADY | feature-dev:* agents 已裝 |
| 26 | typescript-lsp | ❌ SKIP-domain | **ZERO TypeScript**；C++ 要 LSP 用 clangd（另議），非 TS-LSP |
| 27 | security-guidance (重複) | ✅ ALREADY | 同 #15 |
| 28 | context7 | ✅ ALREADY | MCP 已裝（`mcp__plugin_context7_context7`）|
| 29 | playwright | ❌ SKIP-domain | 瀏覽器測試；無前端可測 |
| — | claude-plugins-official / skillsmp.com | 📋 REFERENCE | 瀏覽可，**勿 bulk install**；熱門度 ≠ 適配（restraint memory：多數熱門是 overlap/炒作/正中已知失敗類）|

## 統計
- ✅ ALREADY HAVE：**9**（brainstorming / code-review-expert / Superpowers / security-guidance / code-review / pr-review-toolkit / code-simplifier / feature-dev / context7）
- ❌ SKIP domain-mismatch：**10**（frontend×2 / web-design / remotion / agent-browser / browser-use / mobile / ui-skills / typescript-lsp / playwright）
- 🔴 SKIP conflict：**2**（claude-mem / claude-md-management）
- ⚪ SKIP low/personal：**4**（commit-commands / 2×output-style / claude-code-setup）
- 🟡 PROBE-later：**3**（pdf=P2 論文期 / skill-creator=low / find-skills=low）
- **現在採用：0**

## Bottom line
1. **現在不裝任何一個。** 與 restraint 原則 + 21-agent 前例一致：harness 已超前 2026 基線，真 ROI 在修 latent bug + drift（如本 session 抓出的 active.json 空落差、queue 計數誤報），不在裝工具。
2. **唯一書籤 `pdf`** — 待 G6 進論文寫作期（P2）再評估文獻 PDF 抽取需求；屆時先 grep researcher/knowledge/Zotero 是否已夠。
3. **明確避免** claude-mem（記憶系統衝突）、claude-md-management（與手工治理打架）、所有 web/UI/mobile/TS/video 工具（domain 錯位）。
4. 未來遇「熱門 skill 清單」一律先跑此流程：grep domain-fit → grep 既裝覆蓋 → 對抗驗證 headline → 預設 SKIP。

> 銜接 `docs/references/migration/COMMUNITY_COMPARISON.md`（Round 1/2 Claude Code 研究 plugin + Round 3 自主科研框架）。本次為「具名工具清單」裁決。
