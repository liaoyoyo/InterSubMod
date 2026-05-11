---
id: ism-kb-09-conclusions-negative-02-r1-global-sample-asm
name: "N02: R1-Global Sample ASM 全域驗證 NEGATIVE"
description: "2026-04-21 HCC1395 TO 全 40,237 variants 驗證；SampleASM_Delta residualized AUC 0.527 [0.520, 0.533]，F pilot subset 0.643 為 pre-selection overfit。CL-025a ⭐3→⭐2；CL-008 ≤0.58 ceiling 升至 ⭐5。"
status: archived
last_verified: 2026-04-21
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_r1_global_sample_asm_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, sample-asm, normal-bam, pre-selection-overfit]
canonical_paths: [09_conclusions/negative/N02_r1_global_sample_asm.md]
alias_paths: []
---

# N02: R1-Global Sample ASM 全域驗證 NEGATIVE（Phase 2 TO arm）

## Status
`concluded-negative` (F1-filter)；characterization-only 保留但無 POSITIVE claim。

## Context
Phase 2A Normal Methylation Reference 假說：Sample ASM (Normal vs Tumor 甲基化差異) 可作為 TP/FP filter。F pilot subset (NG=4+AF<0.4+NR≥80+NonLOH, n≈800) 顯示 SampleASM_Delta AUC 0.643，疑似突破 CL-008 ≤0.58 ceiling。需在 global scale 驗證。

## Decision-tested
Phase 2 Normal BAM 整合的 Sample ASM / HPFineF / Normal_HP_Delta 等特徵在全域（非 F pilot subset）是否保留 AUC ≥ 0.58 的區分力。

## Method
- HCC1395 TO arm 全 40,239 variants（TP=28,396 / FP=11,841），Phase 2 pipeline 48.3 min wall time
- Residualized AUC 對 NumReads + Coverage_Multiple（within-group OLS 正確方法）
- Pre-filter subset vs global scale 對比
- Global ∩ F pilot (n=6,965) 交叉驗證

## Result
- **Pre-filter subset AUC**：SampleASM_Delta 0.643、NormalBaseline_Coverage 0.687、HPFineF 0.582
- **Global residualized AUC**：SampleASM_Delta **0.527 [0.520, 0.533]**、Normal_HP_Delta 0.513、HPFineF 0.527
- 全 15 特徵 upper CI 無一超過 0.54
- Global ∩ F pilot subset AUC 0.51-0.53，與全域一致

## Consequences
- **根因**：Pre-selection overfit。F pilot 4-way filter 把 800 個高 TP rate (~90%) region 選出後，特徵的區分力來自「這 800 行」的幸運組合而非生物訊號。
- **替代方向**：Part C.2 Phase 2A Normal Methylation Reference 作 F1-filter 方向 NEGATIVE；Paper positioning D.1.B subclonal heterogeneity biomarker 路線須改以 characterization 為主軸；Pivot 回五大研究目標（read-level epigenetic characterization 主軸）。
- **方法論教訓**：Pre-selection 永遠要 global scale 驗證。未來任何基於 F pilot/NG subset/Zone subset 的 characterization 聲稱必須在 global 40K+ rows 重算 AUC + CI。
- **決策**：2026-04-21 CL-025a ⭐3→⭐2；CL-008 ⭐4→⭐5。

## References
- MEMORY：`project_r1_global_sample_asm_negative.md`
- Finding：`research/F_hpfinengroups_deepening/observations/step8_r1_global_to_arm.md`
- Registry：`docs/reports/research_landscape/10_Research_Chain_Registry.md` CL-025a / CL-008
- Upstream：`output/hcc1395_normal_pilot/TO/`
- Downstream：`output/hcc1395_normal_pilot_global/TO/significance_summary.csv`
