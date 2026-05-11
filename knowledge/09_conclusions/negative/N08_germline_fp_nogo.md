---
id: ism-kb-09-conclusions-negative-08-germline-fp-nogo
name: "N08: TO Germline FP Identification NO-GO（G1-G7 60+ 特徵）"
description: "2026-04-05 G1-G7 七模組 48 圖；60+ VCF 特徵全 AUC<0.64；FP removal at ≤2% TP loss = 0%；根因 PoN 已移除 99.48% germline，殘餘為罕見 germline 本質似 somatic。"
status: archived
last_verified: 2026-04-05
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_germline_fp_identification_nogo.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, germline, vcf-features, pon, to-mode]
canonical_paths: [09_conclusions/negative/N08_germline_fp_nogo.md]
alias_paths: []
---

# N08: TO Germline FP Identification NO-GO

## Status
`concluded-negative`（TO single-sample post-hoc 特徵探索方向正式關閉）。

## Context
TO 模式殘餘 germline FP 是主要 F1 損失來源。假說：雖然 PoN 已移除大部分 germline，剩餘 germline FP 仍可透過 VCF 特徵（GT、AF、strand bias、CpG context、trinucleotide 等）組合識別。

## Decision-tested
VCF 層級 60+ 特徵是否存在組合可達 AUC ≥0.70 或 FP removal ≥10% at ≤2% TP loss。

## Method
- 7 模組 G1-G7 × 48 圖表 × 516K variants × 40+ 新提取 VCF 特徵 × 7 samples
- G2 Germline Signatures (H+GT+AF)
- G3 AF Architecture
- G4 Strand Bias
- G5 CpG Context + transitions
- G6 LR LOSO（7-fold）
- G7 Bootstrap

## Result
- G2 Combined AUC = 0.507（隨機）
- G3 AF raw AUC = 0.566（AF reversal confirmed，高 AF = 更多 FP）
- G4 Strand Bias AUC = 0.502（預期 negative — germline FP 是真 variant）
- G5 is_transition AUC = 0.604（最佳單一特徵）；CpG+Ti combined = 0.616
- G6 LR LOSO Mean AUC = 0.638（0.491-0.725）
- G7 Bootstrap AUC = 0.503 [0.501, 0.504]
- **FP removal at ≤2% TP loss = 0%**

## Consequences
- **根因**：PoN 已移除 99.48% germline FP。殘餘的是 population database 覆蓋不到的罕見 germline variants — 它們與 somatic TP 本質相似（true biallelic variants, real strand support, real AF distribution）。
- **替代方向**：TO FP 過濾已窮盡。正確方向 = **FN rescue + PoN 改進**；結合 O1-O13 的 60+ 特徵 landscape 已形成 CL-008 Beyond-AUC ≤0.58 ceiling 骨幹。
- **決策**：2026-04-05 關閉 TO single-sample post-hoc 特徵探索方向。

## References
- MEMORY：`project_germline_fp_identification_nogo.md`
