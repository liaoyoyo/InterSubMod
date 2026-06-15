---
name: results-report
description: 實驗結果報告撰寫。將完成的實驗 artifact 轉為結構化決策導向研究報告。USE WHEN：「寫實驗報告」「summarize results」「實驗複盤」「results report」。前置：先用 /results-analysis 完成統計分析。輸出到 docs/experiments/*.md。decision section 必繼承 InterSubMod/.claude/skills/scientific-rigor/SKILL.md §2-§7（升 tier 對齊 evidence ladder、結論不超 pre-reg 範圍避 HARKing、NEGATIVE 走 §9.2 Postmortem）。SKIP WHEN 統計分析未完成（先用 results-analysis）、多主題週進度（用 weekly-report）、工程修改報告（用 structured-tech-report）、研究收尾（用 conclude-research）、AI session log（用 report）。
version: 0.1.0
tags: [Research, Reporting, Experiments, Obsidian]
user-invocable: true
---

> **⚠ 2026-05-22 thin wrapper migration**: 本 skill 為 `/narrative-frame` skill 的 **thin wrapper**。
> 預設 framework = **Data-Showcase**（Hypothesis → Data → Stats → Caveats）。
> 等同 `/narrative-frame apply Data-Showcase`。
> 用戶可隨時換 framework：對話中説「換 SCQA」或直接走 `/narrative-frame N1-N6` 動態挑。
> Catalog: `InterSubMod/.claude/skills/narrative-frame/references/framework_catalog.md` §11 hybrid。

# Results Report

Write the **complete post-experiment summary report** after analysis artifacts are ready.

This skill is for the stage **after** `results-analysis`.

## Role boundary

### `results-analysis` does
- strict statistics,
- real figures,
- figure interpretation scaffolding,
- stats appendix.

### `results-report` does
- complete experiment wrap-up report,
- decision-oriented narrative,
- figure-by-figure interpretation inside a coherent structure,
- limitations, failure cases, and next actions,
- Obsidian write-back into `Results/Reports/`.

Do not replace strict analysis with confident prose. If the analysis bundle is missing, first identify the blocker and request or produce the missing bundle.

## Default output

The default report is an **internal research report**, not manuscript prose.

It should be named as:

```text
YYYY-MM-DD--{experiment-line}--r{round}--{purpose}.md
```

Example:
- `2026-03-18--freezing--r03--transfer-summary.md`
- `2026-03-18--contrastive-adversarial--r02--ablation-report.md`

The note title should be:

```text
{Experiment Line} / Round {N} / {Purpose} / {YYYY-MM-DD}
```

Read `references/report-naming.md` before finalizing the filename or note title.

## Required frontmatter

```yaml
---
type: results-report
date: 2026-03-18
experiment_line: freezing
round: 3
purpose: transfer-summary
status: active
source_artifacts:
  - analysis-output/analysis-report.md
  - analysis-output/stats-appendix.md
linked_experiments:
  - Experiments/Freezing-Study.md
linked_results:
  - Results/Freezing-vs-Adapter.md
---
```

## Default report structure

The report must include all sections below.

1. **Executive Summary**
2. **Experiment Identity and Decision Context**
3. **Setup and Evaluation Protocol**
4. **Main Findings**
5. **Statistical Validation**
6. **Figure-by-Figure Interpretation**
7. **Failure Cases / Negative Results / Limitations**
8. **What Changed Our Belief**
9. **Next Actions**
10. **Artifact and Reproducibility Index**

Read `references/report-structure.md` before writing.

## Workflow

### 1. Confirm the report object

Lock these fields first:
- date,
- experiment line,
- round,
- purpose,
- linked experiment note,
- linked durable result note if one already exists.

If round is unknown, do not silently invent a semantic round. Use `r00` only as a temporary placeholder and state that it should be normalized later.

### 2. Read the strict analysis bundle

Minimum required inputs:
- `analysis-report.md`
- `stats-appendix.md`
- `figure-catalog.md`
- actual figures, if available

If these are missing, either generate them first with `results-analysis` or explicitly state which claims cannot be supported.

### 3. Write the report as a decision object

This report is not a transcript of outputs.

Each section must answer a real question:
- What did we test?
- What changed numerically?
- What is actually supported?
- What failed or remains uncertain?
- What should we do next?

Read `references/decision-oriented-analysis.md` for the expected reasoning depth.

### 4. Interpret figures inside the report

Do not only attach figures.

For each main figure:
- introduce why it is included,
- state the key observation,
- explain the supported interpretation,
- explain the decision implication.

Read `references/figure-interpretation.md` and `references/statistical-completeness.md` as needed.

### 5. Choose the write target explicitly

If the current repo is bound to an Obsidian project knowledge base:
- create or update `Results/Reports/{report-name}.md`,
- link back to the relevant `Experiments/` note,
- update the matching canonical `Results/` note when a durable conclusion is now supported,
- append a short trace to today's `Daily/` note,
- update `.claude/project-memory/<project_id>.md`.

If the repo is **not** bound:
- write the report as a local markdown artifact in the requested output location or next to the analysis bundle,
- keep the same filename contract,
- explicitly say that no Obsidian write-back was attempted.

Use `obsidian-project-memory` conventions only for bound repos. Internal experiment reports belong in `Results/Reports/`, not `Writing/`.

### 6. End with explicit next actions

The report must end with operational decisions, for example:
- stop a weak branch,
- schedule one missing ablation,
- promote a stable finding into manuscript-facing writing,
- update the active plan.

## Required quality bar

- The report must be dateable, searchable, and attributable to one experiment line and one round.
- The report must cite actual evidence from the analysis bundle.
- The report must include negative results when they matter.
- The report must separate stable conclusion from tentative interpretation.
- The report must say what changed in project belief and what should happen next.
- **Provenance（§13-C，2026-06-15 audit D6-3）**：報告含 metric → 寫完跑 `python3 scripts/number_provenance.py audit <report.md>`，把「metric→來源檔:行」表納入；任一 metric 無源 → 回填真值或刪。資料密集表優先 `scripts/fill_report.py` template+data 注入（§13-A，缺值 refuse），不手打。

## 嚴謹度繼承（/scientific-rigor）

results-report 是「decision object」— 比 results-analysis 更需要嚴謹度（影響升 tier 決策）。必繼承 `InterSubMod/.claude/skills/scientific-rigor/SKILL.md`：

- **§2 證據分級**: 報告每個 finding 必標 L1-L5 ribbon + 升 tier 建議必對齊 tier
- **§3 Effect Size**: F1/AUC delta 必含 Cohen ribbon + 95% CI（results-analysis 已抽 stats，此處只需引用）
- **§4 DAG**: 因果 claim 必 reference DAG 路徑（`InterSubMod/docs/concepts/DAG/<topic>.md`）
- **§7 Pre-reg**: decision section 必對照 pre-reg `H_預測 / decision_threshold`，禁止 post-hoc HARKing
- **§9.2 Blameless Postmortem**: NEGATIVE / NO-GO 報告必走 SRE 5-段格式

**最小可用子集**:
- Upgrade tier ⭐4-5 報告: §2 + §3 + §4 + §7 全跑（強制觸發 evaluator agent 7-check 通過）
- ⭐3 中段報告: §2 + §3
- 草稿 / 內部 dry-run: §2

**Decision 紅旗**（NEEDS_WORK 訊號）:
- 「mostly positive」「partial」「marginal」未量化為 effect size
- 「跨樣本 X 一致」未列 n / CI
- 結論超出 pre-reg 範圍（HARKing）

## Reference files

Load only what is needed:
- `references/report-structure.md`
- `references/report-naming.md`
- `references/figure-interpretation.md`
- `references/statistical-completeness.md`
- `references/decision-oriented-analysis.md`
- `references/EVIDENCE-PROPAGATION.md`
- `examples/example-results-report.md`
