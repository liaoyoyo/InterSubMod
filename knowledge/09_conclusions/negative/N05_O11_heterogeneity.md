---
id: ism-kb-09-conclusions-negative-05-o11-heterogeneity
name: "N05: O11 Methylation Heterogeneity NEGATIVE (n_reads confound)"
description: "2026-03-31 within-group methylation heterogeneity 假說否決；epipolymorphism 初始 AUC=0.845 完全為 n_reads confound（殘差化後 0.530）；6 特徵 AUC 0.509-0.578。"
status: archived
last_verified: 2026-03-31
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_O11_heterogeneity_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, methylation, heterogeneity, confound, n-reads]
canonical_paths: [09_conclusions/negative/N05_O11_heterogeneity.md]
alias_paths: []
---

# N05: O11 Methylation Heterogeneity 假說 NEGATIVE

## Status
`concluded-negative` (within-group heterogeneity 方向正式關閉)。

## Context
假說：germline ASM（mQTL 驅動）在同 haplotype 組內應比 cancer ASM（stochastic）更一致（更低 heterogeneity），因此 within-group heterogeneity 可區分 germline/somatic variants。

## Decision-tested
Epipolymorphism 等 within-group heterogeneity 特徵是否有獨立於 n_reads 的 TP/FP 區分力。

## Method
- 6 個 heterogeneity 特徵（epipolymorphism 等）
- 殘差化 n_reads + num_cpgs
- Read-count matched bin [81-120] 對照驗證

## Result
- Raw epipolymorphism AUC = 0.845（疑似突破）
- **n_reads confound 暴露**：TP regions 有 1.87× 更多 reads (median 157 vs 92)、n_reads AUC(TP)=0.926、Epipolymorphism ↔ n_reads Spearman r=0.79
- 殘差化後：所有 6 特徵 AUC = 0.509-0.578（近隨機）
- Read-count matched bin [81-120] 中 epipolymorphism AUC=0.560

## Consequences
- **根因**：初始 AUC 0.845 完全是 mechanical artifact（reads 多→heterogeneity 定義自動高）。TP/FP read count 不平衡（paired manifest 中 n_reads 與 truth 高度相關）導致虛假信號。
- **替代方向**：未來所有 region-level TP/FP 分析必須控制 n_reads；germline/somatic methylation 區分的唯一剩餘路徑是 Phase 2 Normal Methylation Reference（後已被 N02 R1-Global 否決）。
- **決策**：2026-03-31 within-group heterogeneity 方向正式關閉。

## References
- MEMORY：`project_O11_heterogeneity_negative.md`
- 文獻調查：`docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md`（ROCIT、MethylBERT、mQTL 等 25 篇）
- O11 workspace：`observation_workspaces/20260331_O11_methylation_heterogeneity/`
