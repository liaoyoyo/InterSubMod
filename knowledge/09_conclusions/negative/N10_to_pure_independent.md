---
id: ism-kb-09-conclusions-negative-10-to-pure-independent
name: "N10: TO-pure Independent Modeling NEGATIVE"
description: "2026-04-08 LOSO 7-fold 4 models × 419,692 TO regions；HP-free ISM AUC ~0.53（隨機）；caller_af alone AUC=0.654 超越全 ISM LR；ISM+Caller RF 增量僅 +0.030；Per-sample FP rate 8.6× spread。"
status: archived
last_verified: 2026-04-08
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_to_pure_independent_modeling_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, to-mode, loso, caller-af, hp-free]
canonical_paths: [09_conclusions/negative/N10_to_pure_independent.md]
alias_paths: []
---

# N10: TO-pure Independent Modeling NEGATIVE

## Status
`concluded-negative`（TO 獨立 ISM-based filter 方向關閉）。

## Context
TO 模式 self-phasing 汙染 HP-dependent 特徵。假說：是否可以**只用 HP-free ISM 特徵**訓練獨立 model（繞過 self-phasing），或 ISM+Caller 是否顯著超越 caller 單獨。

## Decision-tested
4 種 model 配置（A: HP-free only、B: All-ISM、C: Caller-only、D: ISM+Caller）× LR/RF 是否提供 F1 增益。

## Method
- LOSO 7-fold cross-validation
- 7 TO samples，419,692 regions
- Per-sample FP rate 異質性分析

## Result
| Model | Method | Mean AUC |
|-------|--------|----------|
| D: ISM+Caller | RF | **0.662** |
| D: ISM+Caller | LR | 0.636 |
| B: All-ISM | RF | 0.635 |
| C: Caller-only | LR | **0.632** |
| B: All-ISM | LR | 0.601 |
| A: HP-free | RF | 0.535 |
| A: HP-free | LR | 0.529 |

- **HP-free ISM = random** (AUC ~0.53)
- **Caller-only = 0.632**（超越 All-ISM LR 0.601）
- **ISM 增量**：+0.004 (LR), +0.030 (RF)
- **caller_af alone AUC = 0.654** — single best feature
- Per-sample FP rate 0.087 (H2009) 到 0.746 (HCC1954) — 8.6× spread

## Consequences
- **根因**：TO self-phasing 汙染 HP-dependent 特徵；HP-free 特徵有 identifiability 問題（germline ASM 無法與 somatic heterogeneity 區分，見 N15）。
- **替代方向**：TO 立即使用 caller_af 作 best discriminator。ISM 作為 TO filter 需等 PON-only phasing 上游修正後重評估。
- **決策**：2026-04-08 關閉 TO-pure 獨立 ISM modeling 方向。

## References
- MEMORY：`project_to_pure_independent_modeling_negative.md`
- Analysis：`.../20260408_to_pure_independent_modeling/analysis_summary.txt`
- Results：`.../20260408_to_pure_independent_modeling/loso_results.tsv`
