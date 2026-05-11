---
id: ism-kb-09-conclusions-negative-18-ism-af-discriminability
name: "N18: ISM Methylation vs AF Gradient — 候選池偏差"
description: "2026-04-01 HCC1395_5kHz 候選池 vs 全量對比；PassedGating AF≥0.9 富集 0.45×（候選池反效果）→ 1.07×（全量中性）；AF≥0.9 VerificationClass Noise ≥87% 全樣本普遍；ISM 在 LOH (AF≥0.9) 完全失效。"
status: archived
last_verified: 2026-04-01
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_ism_af_discriminability.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, af-gradient, passedgating, candidate-pool-bias, loh]
canonical_paths: [09_conclusions/negative/N18_ism_af_discriminability.md]
alias_paths: []
---

# N18: ISM Methylation vs AF Discriminability Gradient — 候選池偏差修正

## Status
`concluded-negative`（ISM 對 TO FP 過濾效果在全量統計中非常有限 1.0-1.5×）。

## Context
HCC1395 5kHz TO 分析（2026-03-22）初始觀察：ISM 甲基化判別力有 AF 梯度依賴（AF 0.6-0.9 PassedGating 富集 2.00×，AF ≥0.9 反效果 0.45×）。需在跨 7 樣本全量（非 rescue 候選池）驗證是否普遍。

## Decision-tested
PassedGating 的 AF 梯度效果（特別是 AF≥0.9 的 0.45× 反效果）在全量統計是否成立。

## Method
- HCC1395_5kHz（rescue 候選池）vs HCC1395_DORADO（全量）對比
- 7 樣本跨 AF bin 富集統計（全量 ISM 輸出 significance_summary.csv）
- VerificationClass 分佈跨樣本驗證
- LOH-germline-escape 癌症類型特異性分析

## Result
- **候選池偏差暴露**：
  - AF≥0.9 PassedGating 富集：候選池 **0.45× (反效果)** → 全量 **1.07× (中性)**
  - AF 0.6-0.9 富集：候選池 2.00× → 全量 1.13×
  - 0.6-0.9 TP_Strong/FP_Strong：候選池 40%/7% → 全量 40%/35%
- **普遍確認的發現（所有 7 樣本）**：
  - AF≥0.9 VerificationClass Noise ≥ 87%（TP 和 FP 皆然，87-100%）
  - AF≥0.9 PassedGating 通過率 <15% 且 TP/FP 無差距
  - ISM 在 LOH 區域完全失效
- **樣本特異性**：COLO829 QS 中位數=35（其他=75），ISM 幾乎無效；HCC1937 低 AF ISM 反效果；H1437/H2009 AF≥0.9 最高正效果 1.5-1.6×
- **LOH-germline-escape**：HCC1937 (BRCA1 mut) FP≥0.9% 45.9% vs TP 13.1%（**3.50× 富集**）— BRCA/HR-deficient 乳腺癌特有

## Consequences
- **根因**：救援候選池選 borderline FN 候選，FP 代表低品質 FP，導致 ISM 特徵差異被放大。全量統計梯度平坦（0.9-1.6×）。
- **替代方向**：
  - 引用 ISM AF 梯度數據時必須區分「候選池 vs 全量」（全量更可信）
  - 不建議單純依靠 ISM PassedGating 過濾 TO FP
  - ISM 在 LOH（AF≥0.9）完全無效普遍認定
  - COLO829（黑色素瘤）ISM 最無效
- **決策**：2026-04-01 PassedGating AF≥0.9 反效果假說否決。

## References
- MEMORY：`project_ism_af_discriminability.md`
