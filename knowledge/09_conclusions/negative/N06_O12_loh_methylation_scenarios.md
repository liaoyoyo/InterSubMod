---
id: ism-kb-09-conclusions-negative-06-o12-loh-methylation-scenarios
name: "N06: O12 LOH Methylation Scenario Discrimination NEGATIVE"
description: "2026-03-31 TO LOH 三場景甲基化區分 NEGATIVE；AlleleDelta raw AUC 0.66 完全被 caller_af confound 驅動（AF-bin 內 <0.59）；發現 L2 residualization collider bias 機制；175,542 rows × 7 samples 驗證。"
status: archived
last_verified: 2026-04-01
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_O12_loh_methylation_scenarios.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, loh, methylation, collider-bias, residualization]
canonical_paths: [09_conclusions/negative/N06_O12_loh_methylation_scenarios.md]
alias_paths: []
---

# N06: O12 LOH Methylation Scenario Discrimination NEGATIVE

## Status
`concluded-negative`；附帶重要方法論發現 L2 collider bias。

## Context
假說：TO LOH 區域的三種生物學場景（LOH→sSNV、sSNV→LOH、germline LOH）因不同的 allelic structure 而有可區分的甲基化 pattern。

## Decision-tested
AlleleDelta、CramersV、HPMergedSig 等 HP/甲基化特徵是否能在 TO LOH 區域區分 TP vs FP。

## Method
- 數據：175,542 TO LOH rows，9 existing + 13 novel features，7 samples × 4-level confound
- L1 → L2 (residualize on AF) → L3 (AF-bin stratification) 三層驗證
- Novel features：region_methyl_mean/low_fraction、spatial features (transition_count, autocorrelation, block_length_mean)

## Result
1. **AlleleDelta**：Raw per-sample AUC 最高 0.66 (HCC1937)，AlleleDelta↔AF Spearman r = -0.23 to -0.57，AF-bin 內全部 <0.59
2. **CramersV/HeuristicScore/HPMergedSig/PassedGating**：L2 AUC 高達 0.80（**collider bias**）— 這些特徵在 TO 中近乎常數，residualizing on 強 AF confound 產生虛假 residual-outcome 相關；AF-bin (L3) AUC = 0.50
3. **Novel region features**：2-4/7 samples L2 > 0.55，但 < 5/7 → 不穩定
4. **Spatial features**：全部 < 0.55

## Consequences
- **根因**：AlleleDelta ≡ AF confound；常數特徵 residualized on 強 confound 產生 collider bias 虛假信號。
- **替代方向**：TO LOH 最佳策略 = caller features only (AF, GQ, DP)；甲基化區分唯一路徑 = Phase 2 Normal Methylation Reference（後已被 N02 R1-Global 否決）。
- **方法論發現（影響所有後續研究）**：
  - L2 residualization collider bias 偵測方式：L2 >> max(L3) + 0.10
  - 解法：以 AF-bin stratification (L3) 為可靠標準
  - Pooled OLS residualization 是 data snooping（N12 進一步確認）
- **決策**：2026-03-31 TO LOH 甲基化區分方向正式關閉。

## References
- MEMORY：`project_O12_loh_methylation_scenarios.md`
- O12 workspace：`observation_workspaces/20260331_O12_loh_methylation_scenarios/`
