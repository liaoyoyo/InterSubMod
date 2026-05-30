<!--
建立時間: 2026-05-18
目標: .claude/skills/ 視覺化分層索引 — 46 個 skill 的場景對應 + 6 大類 + 依賴 hub 圖
處理範圍: M2 (41 skills 視覺分層) 的方案 2 實作
關聯檔案:
  - InterSubMod/.claude/CLAUDE.md §3 (Skills 分類索引，本檔擴展視覺化)
  - InterSubMod/AGENTS.md §17 (4 層導航)
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md §11 (12 步協作圖)
-->

# InterSubMod Claude Code Skills — 視覺化分層索引

> `.claude/skills/` 下 **46 個 skill** 的分層導覽（45 active + 1 deprecated `/html-preview`；2026-05-28 drift 校正）。
> CLAUDE.md §3 已含簡短分類；本檔提供**場景對應 + 依賴 hub 圖 + 觸發頻率**等視覺化補充。

---

## §1 使用導覽（場景 → skill）

| 你想做什麼 | 用哪個 skill | 階段 |
|---------|------------|-----|
| 開始新研究專案 | `/init-research`（多週級）/ `/cycle-init`（單一假說 cycle）| P0 |
| 學新主題 | `/fast-learning-coach`（費曼+主動回想 5 步）| — |
| 改 C++ 程式碼 | `/methodology-audit` → `/cpp-change` → `/verification-loop` | dev |
| 驗證假說 | `/validation-protocol` L1-L4 + `/auc-confound-guard` 3-gate | P3-P4 |
| 寫週報 | `/weekly-report` → `/pptx-build` 或 `/html-report-build` | output |
| 結論宣告（任何「更好」claim）| **`/scientific-rigor`** ⭐（Hub）| 元 |
| 確認暫停時機 | `/confirmation-protocol`（4 層 Hard Gate / Gate / Review / FYI）| 元 |
| 查歷史教訓（防重複犯錯）| `/known-pitfalls`（P-01~P-14 反例庫）| 元 |
| 新增 .md 檔案 | `/doc-standards`（命名 + 元數據 + 圖片規則）| docs |
| 整理 MEMORY.md | `/memory-consolidation`（掃描 + 合併 + 降級）| docs |
| Spaced recall（30d / 90d 後）| `/scientific-rigor §10.1`（觸發 `/memory-consolidation`）| 元 |
| NEGATIVE postmortem | `/scientific-rigor §9.2` + `InterSubMod/templates/postmortem.md` | 元 |
| Pre-registration（事先註冊假設）| `/scientific-rigor §7.1` + `InterSubMod/templates/research_index.md` | P1 |
| 跨樣本一致性 | `/multi-sample-consistency`（7 樣本驗證）| P4 |
| Tier 升級評估 | `/run-evaluator`（retraction risk score）| P5 |
| 收尾研究方向 | `/conclude-research`（INDEX + MEMORY + ledger 同步）| P6 |
| 找 cycle 狀態 | `/cycle-state`（read-only dashboard）| meta |
| Pipeline 異常 | `/infra-ops`（preflight + disk-full + OOM）| ops |

---

## §2 6 大類分層

### §2.1 元方法論（Meta — 8 skills）

所有任務的**通用基礎方法**。

| Skill | 一句話 | 頻率 |
|-------|------|----|
| `/scientific-rigor` ⭐ | 元 skill — 證據分級 + DAG + Pre-reg + Postmortem + 啟發式工作流映射 | high |
| `/confirmation-protocol` | 4 層確認協議（Hard Gate / Gate / Review / FYI）| high |
| `/known-pitfalls` | AI 已知陷阱清單（P-01~P-14 反例庫）| high |
| `/pre-decision-audit` | entry-point ≤30min 決策審計（7 outputs + Cynefin + 5-dim credibility + GO/PROBE/NO-GO）| mid |
| `/cycle-state` | 跨 cycle dashboard（read-only）| mid |
| `/research-context-loader` | 3-tier landscape 上下文按需載入 | mid |
| `/implementation-notes` | spec 實作中 live 紀錄（設計決定 / 偏離 / 折衷 / 未決）| mid |
| `/fast-learning-coach` | 費曼 + 帕雷托 + 主動回想 + 間隔重複 5 步學新主題 | low |

### §2.2 7-Phase Resilient Waterfall（Framework — 7 skills）

研究 cycle 階段控制。

| Phase | Skill | 一句話 |
|-------|-------|------|
| P0 REGISTER | `/cycle-init` | 短週期 hypothesis cycle 初始化 |
| P1 PLAN | `/research-loop` | 觀察→定向→假設→驗證設計 |
| P2 PRECHECK | `/check-staleness` | binary / dataset / upstream 新鮮度 |
| P3 PILOT | `/feature-layered-observation` | 特徵 6 層分層觀察 |
| P4 GENERALIZE | `/multi-sample-consistency` | 7 樣本 cross-sample 一致性 |
| P5 EVALUATE | `/run-evaluator` | retraction risk 評分 + tier 升級 |
| P6 CONCLUDE | `/conclude-research` | 收尾 + INDEX/MEMORY/ledger 同步 |

### §2.3 程式開發（Development — 4 skills）

C++ 修改 + pipeline 維護。

| Skill | 一句話 |
|-------|------|
| `/cpp-change` | 6 步 PDD C++ 修改協議（baseline → impl → test → feature → review → verify）|
| `/methodology-audit` | 改動前 3-4 方案評估（Step 1-5 IDENTIFY/QUANTIFY/OPTIONS/WRITE/DECISION）|
| `/verification-loop` | 程式碼 build/test 6 phase 驗證（build/type/lint/test/security/diff）|
| `/infra-ops` | pipeline preflight + disk-full / OOM 診斷 |

### §2.4 文件管理（Docs — 5 skills）

.md 規範 + 記憶系統 + 可重現性。

| Skill | 一句話 |
|-------|------|
| `/doc-standards` | 檔名 + 元數據 + 圖片 + 目錄結構規範 |
| `/data-audit` | 研究輸出 6 項組織檢核（圖片連結 / INDEX 覆蓋 / 命名 / 元數據 / gitignore / 散落）|
| `/memory-consolidation` | 記憶生命週期（掃描 / 合併 / 降級 / 索引 < 200 行）|
| `/citation-verification` | 學術引用 WebSearch + Google Scholar 必驗 |
| `/pipeline-manifest` | script→inputs→outputs→figures/tables provenance DAG + orphan 偵測（與 data-audit 分工：查因果鏈）|

### §2.5 報告生成（Reporting — 8 active + 1 deprecated）

從週報母稿到 PPT / HTML 輸出。敘述框架由 `/narrative-frame` 統一挑選。

| Skill | 一句話 | 場景 |
|-------|------|----|
| `/narrative-frame` | 全域敘述框架挑選 + 套用 + 自審（50+ catalog；7 報告 skill 之上層動態挑選）| 框架入口 |
| `/myPPT` | PPT 場景識別**總入口**（路由到 weekly-report / pptx-build）| 路由 |
| `/weekly-report` | W1-W7 週報母稿生成（17 段分流）| 週報 |
| `/pptx-build` | P1-P5 PPTX 製作 + 6 報告模板識別 + Vision check | 給教授 |
| `/html-report-build` | 自包含 HTML 報告 / standalone PI-view / slide deck | 給 PI |
| `/results-report` | 實驗結果決策導向報告 | 單實驗 |
| `/structured-tech-report` | 13 段技術報告（Toyota A3 + ADR + SRE postmortem）| 工程修改 |
| `/report` | AI 對話執行報告（session-end）| 對話結尾 |
| ~~`/html-preview`~~ | **DEPRECATED 2026-05-13** — 取代為 /html-report-build（待移除）| — |

### §2.6 研究專用（Research-specific — 13 skills）

ISM 研究域特定方法。

| Skill | 一句話 |
|-------|------|
| `/auc-confound-guard` | 三關 AUC confound 驗證（within-group OLS / AF-bin / permutation）|
| `/validation-protocol` | L1-L4 假說驗證梯度 + 4-track coverage |
| `/pivot-direction` | 快速切換研究方向 |
| `/inject-hypothesis` | 注入新假說到 queue |
| `/research-dashboard` | 專案級看板（方向 / 假說佇列 / 阻塞）|
| `/provenance-tier-audit` | 跨 cycle 證據鏈一致性審計（system-level）|
| `/problem-framing-ideation` | P0 前置 5W1H + gap analysis（< 2h）|
| `/review-evidence` | 查閱歷史結論（hypothesis_queue + ledger）|
| `/observation-analysis` | 標準化 O-系列觀察分析腳本 |
| `/init-research` | 多週級長期專案初始化（research/<topic>/ scaffolding）|
| `/image-gen` | 雙軌生圖（AI 軌 + cairo 軌）|
| `/image-vision-check` | 6 維度圖片品檢（4/6 過）|
| `/results-analysis` | 嚴格統計分析 + 圖表 |

---

## §3 依賴關係圖（Hub 結構）

```
                      ┌───────────────────────┐
                      │  /scientific-rigor ⭐  │ Hub (元方法論層)
                      │  §1-§11.6              │
                      └───────────┬───────────┘
                                  │
        ┌─────────────┬───────────┼───────────┬─────────────┐
        │             │           │           │             │
        ▼             ▼           ▼           ▼             ▼
  /known-pitfalls  /method   /validation /auc-confound  /verification
  P-01~P-14        -audit    -protocol   -guard         -loop
  (反例庫)         (消融)    (對照組)    (3-gate)       (程式碼級)
        │             │           │           │             │
        └─────────────┴───────────┼───────────┴─────────────┘
                                  │
                      ┌───────────▼───────────┐
                      │ /fast-learning-coach  │ §10.3 互補
                      │ (學新主題)            │ (學 vs 用)
                      └───────────────────────┘

  其他重要 hub:
  - /confirmation-protocol  (4 層 Hard Gate 決策樹)
  - /research-context-loader (3-tier landscape 載入)
  - /myPPT                   (報告路由總入口)
```

**Cross-Reference 全雙向** — 8 個 skill 各自含「與 /scientific-rigor 元方法論的關係」段（commits `c1dde00` + `dce837f`）。

---

## §4 觸發頻率分布

| Tier | 觸發 | Skills | 行為 |
|------|------|-------|------|
| **High freq**（每 session 多次）| AI 自動觸發 | `/scientific-rigor` / `/confirmation-protocol` / `/known-pitfalls` / `/cpp-change` / `/doc-standards` | description USE WHEN 廣 |
| **Mid freq**（週級）| 條件觸發 | `/memory-consolidation` / `/weekly-report` / `/research-context-loader` / `/validation-protocol` / `/cycle-state` | 特定階段或週期 |
| **Low freq**（月級或顯式）| 顯式呼叫 | `/init-research` / `/provenance-tier-audit` / `/pipeline-manifest` / `/fast-learning-coach` / `/citation-verification` | 用戶明示觸發 |

---

## §5 配套 templates

| Template | 用途 | 對應 skill |
|---------|------|----------|
| `InterSubMod/templates/postmortem.md` | SRE-style Blameless NEGATIVE postmortem | `/scientific-rigor §9.2` + `/conclude-research` |
| `InterSubMod/templates/research_index.md` | 新研究 `00_INDEX.md` 含 Pre-reg 3 欄 + G1-G5 + reproducibility 7 項 + 子目錄結構 | `/scientific-rigor §7.1` + `/init-research` |

---

## §6 版本 / 演進

- **目前**: 46 skills（45 active + 1 deprecated `/html-preview`）— 2026-05-28 drift 校正（補列 implementation-notes / narrative-frame / pre-decision-audit / pipeline-manifest，索引移除 grill-me 條目）
- **grill-me 註記**（2026-05-30 補正）: `.claude/skills/grill-me` 為 **dangling symlink**（target `../../.agents/skills/grill-me` 已不存在），無 SKILL.md 故不計入 46。索引已不列；物理清除待 Hard Gate ack。
- **2026-05-17 建立**: `/scientific-rigor` 元方法論層（commit `42217cf`）
- **2026-05-13 deprecated**: `/html-preview`（Python middleware）→ 取代為 `/html-report-build`（LLM-direct）
- **2026-05-18 落地**: F.4 P2 Cynefin → `/confirmation-protocol §Cynefin 域對照`；F.4 P2 Productive Failure → `/scientific-rigor §8.3.1`
- **未來規劃**: F.4 P3 Living-update hook（`skill_change_audit.sh` PostToolUse on `.claude/skills/`）— 詳見 plan `~/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md`

---

## §7 與外部 plugin skills 的關係

本檔僅列 InterSubMod 自製 skill。Claude Code 啟動時還會載入 plugin-provided skills（透過 `system-reminder` 列出），常用如：
- `superpowers:*`（14 個工程技能：brainstorming / TDD / debugging 等）
- `feature-dev:*`（code-architect / code-explorer / code-reviewer）
- `pr-review-toolkit:*`（PR 審查 7 agent）
- `skill-creator` / `plugin-dev:*`（meta：建 skill 與 plugin）

**Plugin skills 用 `plugin:skill` 命名格式**；本檔自製 skill 用單一 slash `/<name>`。

---

## §8 設計原則（呼應業界共識）

依 2026-05-17 業界對照（researcher）整理：

| 原則 | 來源 | 本專案實現 |
|------|----|---------|
| Progressive disclosure | Anthropic Skills 2025-10 | SKILL.md only 載入 metadata；觸發後讀 body |
| ≤200 行 / 檔（SKILL.md）| Anthropic Claude Code docs | 各 SKILL.md 平均 100-300 行 |
| Flat skill library | Voyager / OpenHands / AutoGen | `.claude/skills/<name>/` 平面結構 |
| Meta-skill 層 | **本專案創新**（非業界共識）| `/scientific-rigor` 整合 7 個底層 skill |
| Cognitive Load 7±2 限制 | Miller 1956 / Paas 2020 CLT | `/scientific-rigor §0.5` 最小可用子集 |
| Blameless postmortem | Google SRE | `templates/postmortem.md` |
| Double-loop learning | Argyris 1977 | `/scientific-rigor §11.6` 3 問 |

---

> **本檔不可自動載入** — 是 reference 文件，需用戶或 AI 明示 Read 才會進入 context。
> **更新時機**: 新 skill 建立、deprecated、結構大改時同步本檔。
