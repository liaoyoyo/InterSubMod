---
id: ism-kb-09-conclusions-negative-12-beyond-auc-exhaustion
name: "N12: Beyond-AUC 7 方法綜合驗證 — 特徵空間耗盡確認"
description: "2026-04-09 M1-M7 七種互補統計方法 × 多層分層 × 25 特徵 × 748,391 regions；pure methylation 全 AUC ≤0.58；Pooled OLS residualization data snooping 陷阱發現；HPFineNGroups somatic heterogeneity marker 確認。"
status: archived
last_verified: 2026-04-09
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_beyond_auc_exhaustion_confirmed.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, beyond-auc, residualization, pooled-ols, ceiling]
canonical_paths: [09_conclusions/negative/N12_beyond_auc_exhaustion.md]
alias_paths: []
---

# N12: Beyond-AUC 7 方法特徵空間耗盡確認（CL-008 ceiling）

## Status
`concluded-negative`（ISM 甲基化特徵探索正式關閉，CL-008 ≤0.58 ceiling 確立為 ⭐5）。

## Context
歷經 O1-O13、G1-G7、LOH Wave 1-3、Fine-Pairwise 等 15+ 研究均以 AUC-ROC 為唯一判定。在正式關閉甲基化特徵空間前，需排除 AUC 遺漏信號（calibration、PR-AUC、nonlinear combination、bootstrap、tail enrichment）的可能性。

## Decision-tested
7 種互補統計方法（AUC、PR-AUC、calibration/BSS、residualization、bootstrap CI、nonlinear ML、tail enrichment）是否揭示 AUC 未捕捉的信號。

## Method
- M1-M7 × 多層分層（LOH×AF×Sample）× 25 特徵 × 748,391 regions
- 45+ 篇文獻交叉查核（germline mQTL vs somatic passenger SNV ASM）

## Result
1. **Pure methylation 特徵全部 AUC ≤ 0.58**（pooled），within-group residualized 後更低
2. **PR-AUC 未揭示 AUC 遺漏信號**（Richardson 2024 確認 AUC 穩健）
3. **Pooled OLS residualization 是 data snooping** — confounder-only AUC 0.66-0.85 完全解釋 residualized 提升
4. **Calibration (BSS) 全為負值** — 無校準價值
5. **M5 Bootstrap 200 iterations** — 22 features CI_lo > 0.50 但均 < 0.58
6. **M6 TO ISM+Caller > Caller-only +0.081**（但 sample-dependent；HCC1937/HCC1954 主導，COLO829 反轉；增量來自 HP-dependent 特徵）
7. **M6 Paired ISM 增量僅 +0.006**（不顯著）
8. **M7 零特徵 tail 5% enrichment > 3×**（max 1.25× HPFineNGroups），AUC 未遺漏尾部信號
9. **HPFineNGroups TO non-LOH low AF AUC=0.72, 6/7 samples**（COLO829 at chance）— somatic heterogeneity marker

## Consequences
- **根因**：germline mQTL ASM >> somatic passenger SNV ASM（3-6×）為文獻普遍規律；無直接先例用甲基化區分 variant origin。
- **方法論發現**：Pooled OLS `feature ~ confounders` 在 TP/FP confounders 分布不同時會 leak；**Within-group OLS 是正確方法**，影響所有先前 residualization 結論。
- **替代方向**：ISM 甲基化探索關閉；後續方向 Phase 2A Normal Methylation Reference（後 N02 否決 F1 面向）、CpG selection、somatic heterogeneity characterization。TO HP-dependent 非線性增量 (+0.081) 是 HP-tag 信號非甲基化。
- **決策**：2026-04-09 CL-008 ≤0.58 ceiling 確立。

## References
- MEMORY：`project_beyond_auc_exhaustion_confirmed.md`
- Report：`docs/experiments/finalized/2026/04/20260409_beyond_auc_comprehensive_validation_01.md`
- Literature：`docs/references/20260409_beyond_auc_biology_literature_review.md`
- Data：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260409_beyond_auc_validation/`
