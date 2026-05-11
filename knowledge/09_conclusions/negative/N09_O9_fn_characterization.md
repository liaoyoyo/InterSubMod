---
id: ism-kb-09-conclusions-negative-09-o9-fn-characterization
name: "N09: O9 FN Characterization NO-GO"
description: "2026-04-08 122,790 FN regions × 7 samples；HP-free methylation 全 AUC<0.53；TO Quality_Score INVERTED AUC=0.338；FN VerificationClass ≈ TP（55%/23%/18% vs TP 相同比例）；FN 在 methylation space 與 TP 無法區分。"
status: archived
last_verified: 2026-04-08
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_O9_fn_characterization_nogo.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, fn-rescue, methylation, caller, quality-score]
canonical_paths: [09_conclusions/negative/N09_O9_fn_characterization.md]
alias_paths: []
---

# N09: O9 FN Characterization NO-GO

## Status
`concluded-negative`（ISM-based FN rescue 方向關閉）。

## Context
TO F1 提升的可能路徑之一：透過 ISM 甲基化特徵 rescue caller 漏 call 的 FN（將 FN 重新劃入 TP）。假說：若 FN 是真的 somatic variants，其甲基化 pattern 應與 TP 相似而與非 variant 不同。

## Decision-tested
ISM 甲基化特徵能否在 read-level 識別 FN 為類似 TP 的 somatic variant（即 FN ≠ TP 時提供 rescue 訊號）。

## Method
- 數據：122,790 FN regions (44,415 paired + 78,375 TO) across 7 samples
- HP-free methylation features（NumReads, CramersV, PairwiseMeanDist）AUC on FN vs TP
- Per-sample QS、AlleleDelta、LabelAllelePermanovaF
- VerificationClass 比例對比

## Result
1. **HP-free methylation AUC all < 0.53** — NumReads=0.507, CramersV=0.507, PairwiseMeanDist=0.492
2. **Strongest signal is AF proxy** — LabelAllelePermanovaF=0.664 (paired), AlleleDelta=0.642
3. **TO Quality_Score INVERTED** — AUC=0.338, FN>TP, Cohen's d=-0.671（LOH penalty 反而懲罰 TP）
4. **FN VerificationClass ≈ TP** — 55% Noise, 23% Weak, 18% Strong（同 TP 比例）

## Consequences
- **根因**：FN 在 methylation space 與 TP 完全相同（identical pattern）。FN 被 caller 漏 call 是因為 low AF / borderline evidence，非甲基化差異。
- **替代方向**：不追求 ISM-based FN rescue。若 FN rescue 需要，應用 caller features（AF, DP, GQ）或 read-level evidence 如 haplotype consistency。
- **決策**：2026-04-08 關閉 ISM FN rescue 方向。

## References
- MEMORY：`project_O9_fn_characterization_nogo.md`
- Report：`.../20260408_O09_fn_characterization/O9_FN_characterization_report.md`
- Manifest：`.../20260408_O09_fn_characterization/fn_manifest.tsv.gz` (122,790 rows)
