---
id: ism-kb-09-conclusions-negative-14-wave3-loh-stratification
name: "N14: Wave 3 LOH 分層分析確定性關閉"
description: "2026-04-06 6 模組 J11-J16 判定；Non-LOH max AUC=0.643；LOH Voting AUC=0.577<0.58；cnLOH PairwiseMeanDist 0.587 為 Simpson's Paradox（per-sample mean 0.50）；CramersV 被 NumReads confound；AlleleDelta=0.556 為 LOH 唯一 confound-free 但不足。"
status: archived
last_verified: 2026-04-06
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_wave3_loh_stratification_closure.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, loh, wave3, simpsons-paradox, cnloh, cramersv]
canonical_paths: [09_conclusions/negative/N14_wave3_loh_stratification.md]
alias_paths: []
---

# N14: Wave 3 LOH 分層分析關閉

## Status
`concluded-negative`（LOH stratification 作為 TP/FP 區分路線正式關閉）。

## Context
Wave 1-2 已測試 LOH 單特徵區分力。Wave 3 探索多特徵組合、cnLOH 子區域、特徵方向修正是否提升 LOH 區分力。

## Decision-tested
1. LOH 多特徵組合是否達 AUC ≥ 0.58
2. cnLOH 子區域是否提供獨立信號
3. 特徵反向 (FP>TP) 修正後是否可用

## Method
- Wave 3 六模組 × 7 samples
- Voting / stacking 多特徵組合
- cnLOH per-sample 分析
- Residualize NumReads 對 CramersV 的影響

## Result
| # | 判定 | 證據 |
|---|------|------|
| J11 | Non-LOH 區分力同樣有限 | max AUC=0.643 (HPFineNGroups, read count proxy) |
| J13 | LOH 多特徵組合不可行 | best Voting AUC=0.577 < 0.58 |
| J14 | cnLOH PairwiseMeanDist 0.587 是 Simpson's Paradox | per-sample mean AUC=0.50, 5/7 一致性 |
| J15 | CramersV 被 NumReads confound | AUC 0.511→0.464 after residualize |
| J16 | AlleleDelta 是 LOH 唯一 confound-free 信號 | AUC=0.556, 7/7, 但不足作 filter/combo |

## Consequences
- **根因**：LOH 區域 TP/FP 均受 self-phasing / CN / AF 結構性驅動，甲基化無獨立訊號；Simpson's Paradox — cnLOH 整體 AUC 0.587 被 H2009 (N=31,970, 48.6%) 主導。
- **新發現**：CramersV NumReads confound（首次發現 LOH 區域 CramersV AUC 大部分來自 NumReads，殘差化後低於隨機）。
- **ISM Code 影響**：
  - Code-1 (LOH-Aware dimension switching) 降為 P2（AlleleDelta AUC=0.556 不足）
  - Code-2 (QS TO 重設計) 維持 P1，方向改為「消除反向特徵」而非「加入新正向特徵」
- **替代方向**：不再投入資源在 LOH-based TP/FP 區分；未來轉向 read-level 或 FN rescue（後者已被 N09 否決）。
- **決策**：2026-04-06 LOH stratification 路線關閉。

## References
- MEMORY：`project_wave3_loh_stratification_closure.md`
- Report：`docs/reports/validated/2026/04/20260406_LOH雙定義交叉分析報告/07-09`
- 視覺檢視：`20260406_visual_inspection/`
