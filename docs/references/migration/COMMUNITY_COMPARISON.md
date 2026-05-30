# 社群 Claude Code 研究 plugin vs InterSubMod 內部 harness — 比較與改進 Roadmap

> **日期**：2026-05-28
> **目的**：記錄社群著名研究/學術 plugin 與本專案 44+ skill harness 的逐項差異、已採納的改進、以及論文期採購清單。
> **證據級別**：L1 = 讀實際 SKILL.md 全文 / L2 = 讀該 repo 官方 docs / L3 = README 推論。
> **方法**：對齊 `/scientific-rigor §2` 證據分級；每條 claim 標 tier。

---

## 1. 比較對象

| Repo | 規模 | 領域貼近度 | 主力證據級別 |
|------|------|-----------|-------------|
| [flonat/claude-research](https://github.com/flonat/claude-research) | 48 skills / 6 agents / 9 hooks / 8 rules | PhD 學術 workflow（通用） | **L2**（讀了 docs/skills.md + docs/rules.md） |
| [matsengrp/plugins](https://github.com/matsengrp/plugins) | 11 agents / 1 command | 生信實驗室（最貼近） | L3 |
| [K-Dense-AI/claude-scientific-writer](https://github.com/K-Dense-AI/claude-scientific-writer) | 19+ skills | 科學寫作 | L3 |
| [imbad0202/academic-research-skills](https://github.com/imbad0202/academic-research-skills) | 4-stage 多 agent 套件 | 學術全流程（orchestrator） | L3 |

**一句話結論**：本 harness **研究過程嚴謹度**為核心（強）、論文產出為空白；社群 plugin **論文產出全鏈**為核心、過程嚴謹度多僅 1 條 rule。兩者互補 > 競爭。

---

## 2. 維度對照（強 / 相等 / 缺）

### 2.1 本 harness 明顯強或相等

| 維度 | 本 harness（L1） | 社群對應（L2） | 判定 |
|------|------|------|------|
| 科學嚴謹方法論 | `scientific-rigor` 603 行 11 章（L1-L5 tier / DAG / Effect Size / Pre-reg / Reflexion / Productive Failure / SRE postmortem） | claude-research **僅 Rule 1「Design Before Results」單條** + `devils-advocate` skill | **本 harness 強（decisive）** |
| claim 證據追蹤 | evidence_ledger + `provenance-tier-audit` + 5 hook gate | `[LEARN]` tag + promise-checker hook | **本 harness 強** |
| 研究流程狀態機 | 7-Phase Waterfall + `state/cycles/` + `run-evaluator` 6-component | `init-project-research` + `handoff`，無 phase gate | **本 harness 強** |
| 任務分級/確認 | 6 類 task-type gate + 影響×信心矩陣 + Hard Gate | Rule 8 Scope + Rule 6 Plan | **本 harness 強** |
| memory 生命週期 | `memory-consolidation` | `consolidate-memory` + `memory-cleanup` | 相等 |
| session 交接 | `session_start_inject_focus` + CURRENT_FOCUS | `handoff` + `update-focus` | 相等 |
| code review | pr-review-toolkit 5 agents | `code-review` 平行 specialist | 相等（領域差：C++ vs R/Julia） |

### 2.2 社群有、本 harness 實際缺（已 L1/L2 驗證）

| # | 空缺能力 | 社群來源 | 本 harness 現狀(L1) | 急迫度 | 處置 |
|---|---|---|---|---|---|
| **G1** | pipeline-manifest（script→figure provenance DAG） | claude-research `pipeline-manifest` | grep 零命中；data-audit 不做因果鏈 | 🔴 高 | ✅ **已建 `/pipeline-manifest`**（2026-05-28） |
| **G2** | LaTeX 工具鏈（autofix/proofread/health-check/beamer） | claude-research `latex*` `proofread` `beamer-deck`；matsengrp `scientific-tex-editor` | 完全沒有 | 🔴 高（論文期） | 🟡 P2 論文期裝 |
| **G3** | .bib key 交叉驗證（missing/unused/typo'd \cite） | claude-research `bib-validate`；matsengrp `journal-submission-checker` | citation-verification 驗存在性+claim 真實性，不做 key↔\cite 對照 | 🟡 中（互補） | 🟡 P2 論文期裝 |
| **G4** | referee 回覆工作流（PDF 意見→DAG 修訂） | claude-research `parse-reviews` `strategic-revision` `referee2-reviewer` | 完全沒有 | 🔴 高（revision 期） | 🟡 P2 投稿後裝 |
| **G5** | skill 建立前 dedup 防 drift | claude-research `creation-guard` `skill-preflight` | 無自動檢查 → 已發生 44/45 drift | 🟡 中 | ✅ **已建 `creation_guard.sh` hook**（2026-05-28） |
| **G6** | /compact 前後狀態保全 | claude-research `precompact-autosave.py` | 僅 compact_test.sh（manual） | 🟡 中 | ✅ **已建 `precompact_autosave.sh` hook**（2026-05-28） |
| **G7** | 獨立 devils-advocate / multi-perspective | claude-research `devils-advocate` `multi-perspective` | adversarial 散落 6 skill，無獨立紅隊 | 🟡 中 | ⏸ 暫不自建（已散落覆蓋，ROI 低） |
| **G8** | Snakemake reproducible pipeline | matsengrp `snakemake-pipeline-expert` | bash + C++ 自建 | 🟡 中 | ⏸ 可選評估（browse 即可） |

---

## 3. 已採納改進 log（2026-05-28）

| 改進 | 對應空缺 | Artifact | 性質 |
|------|---------|----------|------|
| 修 CLAUDE.md §3 計數 drift（44→45，grill-me 確認不存在） | G5 衍生 | `.claude/CLAUDE.md §3` | 衛生 |
| 新建 `/pipeline-manifest` skill | G1 | `InterSubMod/.claude/skills/pipeline-manifest/SKILL.md` | 做事 |
| 新建 `precompact_autosave.sh` hook（PreCompact 事件） | G6 | `InterSubMod/scripts/hooks/precompact_autosave.sh` | 韌性 |
| 新建 `creation_guard.sh` hook（PreToolUse Write） | G5 | `InterSubMod/scripts/hooks/creation_guard.sh` | 做事 |

---

## 4. P2 論文期採購清單（動筆寫論文那天再裝）

> 觸發條件：研究進入論文撰寫 / 投稿 / revision 階段（目前 `feedback_paper_positioning_de_prioritized` 降優先中）。

```bash
# LaTeX + bib + referee 全鏈（最對症 G2/G3/G4）
/plugin marketplace add flonat/claude-research
# → 挑裝：latex-autofix / proofread / bib-validate / parse-reviews / strategic-revision
#         + agents: referee2-reviewer / paper-critic

# 生信實驗室 LaTeX 投稿鏈（matsengrp，領域最近）
/plugin marketplace add matsengrp/plugins
# → 挑裝 agents: scientific-tex-editor / tex-grammar-checker /
#                pdf-proof-reader / journal-submission-checker
#   可選現在評估: snakemake-pipeline-expert (G8)

# poster + grant（會議/申請才需要）
/plugin marketplace add K-Dense-AI/claude-scientific-writer
# → 挑裝: latex-posters / research-grants
```

**裝法原則**：`/plugin browse` 後**挑單 skill 安裝，勿全包**（避免與本 harness 既有 7-Phase / scientific-rigor / narrative-frame 重疊打架）。

**不裝**：academic-research-skills 整套（13+12+7+10 agent orchestrator，與本 harness 高度重疊，借鏡其 Devil's Advocate 概念即可）。

---

## 5. 驗證限制聲明

- 社群 plugin 判定主力為 **L2（官方 docs）**，少數 **L3（README）**；未 git clone 讀每個 skill 原始碼。對「缺口識別」目的 L2/L3 足夠（本 harness 明確無 latex/bib-key/referee/manifest）。
- 本 harness 判定全 **L1**（讀實際 SKILL.md）。
- 未做：實際安裝任一社群 plugin 跑行為測試（屬 L1 行為驗證，超出本次「比較差異」範圍）。

---

## 5b. Round 2 — 設定 + 生態系審計落地（2026-05-30，Opus 4.8 Dynamic Workflow 6-路研究）

> 6-路 workflow 審計（permission docs / model+effort / industry plugins / dynamic workflows / agents-subagents / internal skills）後落地。完整 raw findings: workflow run `wf_03d29c90-316`。

### 設定變更（已落地 + 官方文件 L1 佐證）

| 設定 | 變更 | 官方依據 |
|------|------|---------|
| `defaultMode` | project 層 `dontAsk` **移除** → 繼承 user `auto` | v2.1.142 起 project/local `auto` 被忽略；移除是唯一正解（code.claude.com/docs/en/permission-modes）|
| `model` | pin `claude-opus-4-8[1m]` | alias 會隨 Opus 4.9 漂移污染跨 cycle benchmark；完整 ID 鎖版（model-config）|
| `effortLevel` | project 層釘 `xhigh` | issue #49076：/model picker 默默改 effortLevel；project 層每次啟動重套防護 |

> 強制限制排查（L1 實證）：無 managed-settings.json、無 permission-forcing env var → 限制純來自 project config，可改。

### Dynamic Workflows 定位（互補非取代，已寫進 CLAUDE.md §8）

- **3 場景改用**：跨樣本 benchmark（`.claude/workflows/cross_sample_benchmark.js` 已建模板）/ 文獻 cross-check（`/deep-research`）/ NO-GO 前 stress-test
- **絕不用**：含 Hard Gate 環節（workflow subagent 一律 acceptEdits、繞過 exit-2 hook）、互動探索、5 分鐘小活
- `/effort ultracode` 不設長期預設（繞過確認矩陣 + 放大 token）

### agents 變更（已落地）

- **新增 2 fresh-context agent**：`reproducibility-audit`（包 /pipeline-manifest DAG）+ `security-reviewer`（C++ 記憶體 + pipeline 路徑注入）→ 16→18 agent
- **read-only 機械 agent 加 `model: haiku`**：research-orchestrator（路由）+ narrative-organizer（萃取）；驗證/coding agent 維持 inherit
- **evaluator 三合一合併（deferred）**：需刪 3 檔（Hard Gate）+ 改 20 檔 cross-ref（風險高）→ 待明確確認

### skills drift 修正（已落地）

- README.md 43→46 + 補列 implementation-notes/narrative-frame/pre-decision-audit/pipeline-manifest + 移除 phantom grill-me
- 5 個 SKILL.md 的 grill-me forward-link → 改指 /scientific-rigor §11.6 或 /pre-decision-audit
- 新增 `skill_registry_sync.sh` hook（雙向 drift 守衛：編 README/CLAUDE.md 時比對磁碟實際計數）

### 業界 plugin/skill（top 3）

| # | 項目 | 採用 | source |
|---|------|------|--------|
| 1 | 原生 Dynamic Workflows + /deep-research | 不裝 plugin，已內建 | code.claude.com/docs/en/workflows |
| 2 | **K-Dense scientific-agent-skills VCF 子集**（26.5k★ MIT, v2.45.0 5/28）| `npx skills add` 後**只留** pysam/Ensembl VEP/ClinVar/COSMIC/TileDB-VCF（不整包 140）| github.com/K-Dense-AI/scientific-agent-skills |
| 3 | 借鏡 passport.yaml provenance schema（CC-BY-NC，只取設計）| 概念加進 evidence_ledger | imbad0202/academic-research-skills |

**不裝**：wshobson/agents、barkain orchestration、2389 plugins（與既有 18 agent 重疊；原生 workflow 出後價值更低）。

### Deferred / Resolved（2026-05-30 對抗式 impact-audit workflow `wf_4fccc199-341` 裁決）

- **evaluator 三合一 → RESOLVED = KEEP-SEPARATE**（不合併）。決定性理由（L1，已對證磁碟）：3 agent 的 `tools` 權限不可參數化 —— evaluator `Read,Grep,Glob,NotebookRead`（機械保證無 Write/Bash = Fresh-Context Evaluator 核心）vs reviewer `+Write,Bash(python3)` vs methodology-reviewer `+Bash(python3,ls)`。合併取並集 → evaluator 喪失 no-write 保證；取交集 → 另二者功能殘廢。HYBRID 共用前言去重屬低價值 churn，不採。
- **D1 移除 deprecated `/html-preview` → 待 Hard Gate ack**：`git rm -r .claude/skills/html-preview/` + 8 檔 lockstep（CLAUDE.md §3 46→45 + 視覺化(4)→(3) / README ×4 / image-gen:25 / image-vision-check:15,25 + playbook:33 / myPPT/evals.json:55,68）。單獨 commit。
- **grill-me dangling symlink → 待 Hard Gate ack**：`.claude/skills/grill-me -> ../../.agents/skills/grill-me`（target 不存在）；`git rm .claude/skills/grill-me`。可與 D1 同批。
- **P1 hook 計數 36→37（已修）** + **P2 grill-me「完全不存在」假陳述（已改為 dangling symlink 正確描述）** — 2026-05-30 非破壞性修正完成。

---

## 6. Lineage

- 比較執行：Claude Code Opus 4.8 session（2026-05-28 round 1 + 2026-05-30 round 2）
- 來源 WebFetch：上述 4 個 GitHub repo + code.claude.com/docs 官方文件（round 2）
- 內部 skill：`InterSubMod/.claude/skills/{scientific-rigor,citation-verification,run-evaluator}/SKILL.md`（L1 讀取）
- Round 2 workflow run: `wf_03d29c90-316`（7 agents, 538k tokens）
