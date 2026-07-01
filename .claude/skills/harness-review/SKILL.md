---
name: harness-review
description: 週度流程/架構 meta-review + restraint gap 分析。把「這週進度與問題 → 對照既有覆蓋只留真 gap → harness 健康 → 網路新實務 delta → 修正清單」固化成可一鍵重跑的流程。USE WHEN 用戶說「確認這週進度與問題」「是否需要修正或改進」「網路最新架構新知」「持續疊代學習改進」「週度盤點」「harness review」、跑完 /insights 想把摩擦點落地、月度自我稽核架構層。SKIP WHEN 研究內容週報（用 /weekly-report）、純 config drift 11 燈（用 /harness-health）、單一 bug fix、研究假說驗證、純 build/commit。
globs:
  - ".claude/skills/**/SKILL.md"
  - ".claude/rules/**/*.md"
  - "scripts/hooks/**/*.sh"
---

# harness-review — 週度流程/架構 meta-review + restraint gap

> **一句話**：每週把「進度 + 問題 + 新想法/insights 建議」**逐項對照既有 skill/hook/rule 是否已實作**，只留**真 gap** 修，避免重造輪子。restraint 是本 skill 的靈魂。

## 定位 / Phase

- **Phase**: cross-cutting 元方法論層（非 P0-P6 專屬）；retrospective 週度/月度觸發。
- **與既有分工（正交，勿混）**：
  - `/harness-health` = **config drift 11 燈**（機械、read-only、grep 磁碟 count/hook/tier/staleness）。
  - `/weekly-report` = **研究內容週報**（實驗結果、證據鏈；用 narrative-frame）。
  - `/harness-review`（本）= **流程/架構週度 meta-review + restraint gap**（把新想法對照既有覆蓋、產修正清單）。三者互補。

## 何時用 / SKIP

- **USE**：「確認這週進度與問題」「是否需修正/改進」「網路最新架構」「持續疊代」「週度盤點」；跑完 `/insights` 想落地摩擦點。
- **SKIP**：研究內容週報（`/weekly-report`）；純 11 燈（`/harness-health`）；單 bug fix；純 build/commit/docs。

## 5 步流程

### ① 本週進度摘要
```bash
git log --since='<7d前>' --until='<今>' --oneline research/subclonal-reconstruction-202606 | head -40
```
+ 讀 memory `MEMORY.md`（最新條目）+ `docs/CURRENT_FOCUS.md`。歸納：完成的 milestone + 收斂結論。
> 大範圍/多面向 → 可派 `Explore` agent（medium）平行掃 git+memory+postmortem，回結構化摘要（避免主 context 塞爆）。

### ② 本週問題
讀 `docs/postmortems/hook_failures.log` tail + `research/autoresearch/evidence_ledger.jsonl` tail + `git status --short | wc -l`（並行髒度/worktree 數）。歸納：blocker、被抓的捏造/overclaim、多-session 干擾、stale doc。

### ③ 🔴 Restraint gap 檢查（核心）
把「新想法 / insights 建議 / 網路新實務」**逐項**對照既有：
```bash
# 對每個候選改進，grep 既有覆蓋（skills/hooks/rules/scripts/memory）
grep -rl "<關鍵字>" .claude/skills .claude/rules scripts/hooks scripts 2>/dev/null | head
```
- 每項標 **ALREADY-COVERED（引 path）** 或 **GENUINE-GAP**。
- **只有 GENUINE-GAP 才進修正清單**。ALREADY-COVERED 明標 SKIP + 引用既有實作（示範 restraint）。
- 依據 memory `feedback_harness_restraint_over_adoption`（裝前查覆蓋 + 對抗驗證 gap）。
> 大量候選 → 派 `Explore` agent 做逐項 grep 對照（見本 session 06-25 週 review：8 insights 建議 → 7 已實作、1 真 gap）。

### ④ 網路新實務 delta（選配）
`WebSearch` 查該領域最新（agentic harness / verification / memory / rate-limit…）。
- 多數新知常**驗證既有設計**（tiered memory=MEMORY.md、context compaction=concise-emit）→ 標「已對應」。
- 🔴 網路結果的 arXiv ID / 論文數字**引用進任何 SoT 前必 `/citation-verification`**（future-dated / 幻覺 ID 風險）。

### ⑤ 產「真 gap + 修正建議」清單
輸出表：`| gap | 根因 | 修正 | 影響/工時 | 是否撞 DEAD tombstone |`。
- 依 CLAUDE.md §1 暫停矩陣 + task-type 判 scope。
- 大改動走對應協議（C++ → `/cpp-change`；新 skill → 同步 §3 計數）。

## Restraint 紀律（本 skill 靈魂）

1. **裝前查覆蓋**：任何「加 skill/hook/rule/工具」前先 §③ grep 既有；命中即 SKIP。
2. **7/8 已覆蓋是常態**：成熟 harness 下，多數「建議」已實作；找出**真 gap** 才是價值。
3. **不重造既有**：report skill / 反捏造 / PERMDISP / anti-overclaim / test-lint / fan-out cap 均已存在（勿再造）。
4. **純文字規則優先於機械守衛**：路由/操作類判斷用文字規則即足；機械 hook 只在**實際 in-the-wild 事故**後加（§13 三層是唯一因反覆捏造才加機械層的例外）。

## 🔴 DEAD tombstone（不可提議）

- 全自動研究迴圈 / 排程 heavy compute / external board 當 canonical memory（memory `project_loop_engineering_harness_review`）。
- 本 skill 本身 = **user-invoked 手動**（非 auto-loop）→ 不撞 tombstone。

## 輸出格式

Tier 3 結構化（narrative-frame 可套 A3/PDCA）：**① 進度 → ② 問題 → ③ gap 對照表（多數 SKIP、少數真 gap）→ ④ 網路 delta → ⑤ 修正清單 + 建議**。落地修正走各自協議 + `/harness-health` 複驗。

## 參考
- restraint: memory `feedback_harness_restraint_over_adoption`
- 主工作流 SoT: `docs/references/20260611_master_workflow_architecture_01.md`
- tombstone: memory `project_loop_engineering_harness_review`
- 落地實例: 本 skill 由 2026-07-01 週度 review 產出（8 insights 建議 → 7 ALREADY-COVERED / 1 真 gap〔CJK 字型〕）。
