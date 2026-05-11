---
id: ism-kb-09-conclusions-negative-07-o13-cross-region
name: "N07: O13/O13v2 Cross-Region Methylation Correlation NEGATIVE"
description: "2026-04-01 TO 跨區域甲基化 correlation NEGATIVE；表面 FP-FP r=0.739 > TP-TP 0.623 是 shared read count confound 驅動（FP-FP shared=36 vs TP-TP=21）；雙重分層 + residualize + matching 後方向不一致；ML/CpG-SNP 補充驗證全 < 0.565。"
status: archived
last_verified: 2026-04-01
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_O13_cross_region_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, cross-region, methylation, shared-reads, confound]
canonical_paths: [09_conclusions/negative/N07_O13_cross_region.md]
alias_paths: []
---

# N07: O13/O13v2 Cross-Region Methylation Correlation NEGATIVE

## Status
`concluded-negative`（含所有補充驗證 Gap 1-3）；置信度 9/10。

## Context
假說：同 subclone 的 TP-TP pairs 跨區域甲基化 correlation > FP-FP pairs（因 FP 多為 germline variants，mQTL 驅動的 ASM pattern 較隨機）。

## Decision-tested
跨區域甲基化 correlation 是否有獨立於 shared read count confound 的 TP/FP 區分力。

## Method
- 數據：8,400 pairs，5,756 valid (shared≥5)，7 samples
- O13v1：雙重分層（6 strata）
- O13v2：OLS residualized + propensity matched (n=500) + per-sample
- Gap 1-3 補充：20+ 未測試 ISM 欄位 / ML 組合 / CpG-SNP proxy

## Result
1. **O13v1 表面結果**：FP-FP median r=0.739 > TP-TP 0.623（p=8.38e-13），**方向反轉**
2. **Root cause confound**：FP-FP median shared reads=36 vs TP-TP=21
3. **O13v2 嚴格驗證**：
   - 雙重分層：3 TP>FP + 3 FP>TP（方向不一致）
   - OLS residualized：delta=0.028, p=0.464, d=-0.071
   - Propensity matched：delta=-0.037, p=0.403
   - Per-sample：resid 3/7 TP>FP（無一致方向）
4. **補充驗證**：20+ 未測試 ISM 欄位 pooled AUC < 0.565；ML 組合 additive < +0.013 over caller；CpG-SNP proxy AUC=0.449（反轉）
5. **生物學**：HP concordance=1.0（shared reads 同 haplotype）；ALT concordance TP-TP=0.700 vs FP-FP=0.994（variant calling 特性，非甲基化）

## Consequences
- **根因**：shared read count 是 pair-level correlation 的主要驅動；ISM C++ 架構為單區域設計，無原生跨區域能力。
- **替代方向**：與 N05（O11 within-group）、N06（O12 LOH scenario）構成完整證據鏈 — TO 甲基化在所有維度（within-region, LOH, cross-region）+ ML 組合 + 未測試欄位 + CpG-SNP 均 NEGATIVE。多點甲基化關聯性方向正式關閉。
- **決策**：2026-04-01 關閉跨區域甲基化方向。

## References
- MEMORY：`project_O13_cross_region_negative.md`
- O13v1 workspace：`observation_workspaces/20260331_O13_cross_region_correlation/`
- O13v2 workspace：`observation_workspaces/20260401_O13v2_cross_region_confound_control/`
