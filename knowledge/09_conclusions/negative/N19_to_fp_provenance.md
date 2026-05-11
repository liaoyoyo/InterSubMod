---
id: ism-kb-09-conclusions-negative-19-to-fp-provenance
name: "N19: TO FP Provenance — ISM 無法進一步過濾"
description: "2026-03-22 HCC1395 5kHz TO 全量 FP 分析；11,598 殘餘 PASS FP 中 98.6% (11,430) 是 raw_absent（ClairS-TO 在無 normal 下 call 出、ClairS-paired 從未輸出）；GQ/QS 分佈 TP/FP 完全相同。"
status: archived
last_verified: 2026-03-28
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_to_fp_provenance.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, to-fp, provenance, raw-absent, clairs-to, fn-rescue]
canonical_paths: [09_conclusions/negative/N19_to_fp_provenance.md]
alias_paths: []
---

# N19: TO FP Provenance — FP 過濾方向關閉

## Status
`concluded-negative`（TO FP 過濾方向關閉；正確方向為 FN rescue）。

## Context
TO 模式 FP 是主要 F1 損失。需先了解 FP 來源結構，才能判斷 ISM 過濾是否可行。將 TO PASS FP 分為 5 類 provenance（raw_absent、persistent、raw_filter、longphase_s、longphase_w）。

## Decision-tested
TO PASS FP 的 provenance 分佈是否允許 ISM 過濾策略（即 FP 集中在某類可識別 provenance）。

## Method
- HCC1395 5kHz TO 全量 FP 分析
- FP provenance 5 類分類
- GQ、Quality_Score 分佈 TP vs FP 比對
- ISM 覆蓋率量化

## Result
- TO 殘餘 PASS FP：11,598 個
  - **raw_absent = 11,430 (98.6%)** — ClairS-TO 無 normal 下 call 出、ClairS-paired 從未輸出
  - persistent = 87
  - raw_filter = 63
  - longphase_s = 18
- **GQ 和 Quality_Score 分佈在 TP/FP 之間完全相同** — 任何 filter 都同等傷害 TP
- 當前最佳 rescue：H012 GQ>=3 → delta_F1 = +0.009365（F1=0.722062）
- ISM 覆蓋率：773/11,051 FN = 7%

## Consequences
- **根因**：TO FP 本質是 ClairS-TO calling model 在無 normal 情況下的系統性假陽性；這些 variants 在 VCF-level 與真 somatic 不可區分。
- **替代方向**：TO F1 提升的唯一有效路徑是 **FN rescue（擴大 ISM 覆蓋率，當前 7% → 更大）**，而非 FP 過濾。
- **決策**：2026-03-28 不再嘗試 TO FP 過濾方向。

## References
- MEMORY：`project_to_fp_provenance.md`
- Report：`docs/reports/validated/2026/03/20260322_TO_FP_provenance_methylation_analysis_01.md`
- Data：`research/fp_provenance/20260322_hcc1395_5khz_to/`
