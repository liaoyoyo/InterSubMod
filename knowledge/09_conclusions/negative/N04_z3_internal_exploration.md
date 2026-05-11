---
id: ism-kb-09-conclusions-negative-04-z3-internal-exploration
name: "N04: Z3 (TO LOH+extreme AF+NG≤1) 內部二階特徵探索 NEGATIVE"
description: "2026-04-18 Step 1/2.5/3 完成；12 特徵 × 7 樣本無任何 |AUC|≥0.60 於 ≥3 樣本；AF∈[0.4,0.6] × CN × NGroups germline 分層僅 1/7 樣本；甲基化特徵 NG≤1 下結構性退化 ~0.50。"
status: archived
last_verified: 2026-04-18
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_z3_internal_exploration_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, zone-aware, z3, amplicon, hcc1954]
canonical_paths: [09_conclusions/negative/N04_z3_internal_exploration.md]
alias_paths: []
---

# N04: Z3 Internal Feature Exploration NEGATIVE

## Status
`concluded-negative` (Zone-Aware Framework 定位維持 characterization only)。

## Context
Zone-Aware Framework 中 Z3（TO LOH + extreme AF + NGroups≤1）被認為是最 FP-enriched zone。假說 H-Z3d：在 Z3 內用 AF × CN × NGroups 二階分層可分離 germline residue。

## Decision-tested
1. Z3 內部是否存在跨樣本穩定的二階特徵（|AUC|≥0.60 於 ≥3 樣本）
2. AF∈[0.4,0.6] band 是否可用 CN+NG 分層識別 germline pattern

## Method
- Step 1：12 特徵 × 7 樣本 AUC 掃描
- Step 2.5：AF×CN×NGroups 三維分層（H-Z3d 假說測試）
- Step 3：HCC1954 amplicon FP 結構分析

## Result
- **Step 1**：無特徵在 ≥3 樣本達 |AUC|≥0.60
- **Step 2.5**：僅 1/7 樣本符合 germline pattern
- 甲基化特徵（AlleleDelta/HPFineF/HPMergedDelta）在 NG≤1 結構性退化至 ~0.50
- **HCC1954 特殊機制**（已知）：Z3 FP 集中 chr5 (43%) + chr8 (29%) + chr17 (13%)（HER2/MYC amplicon），FP NumReads=55 vs TP=37（p=4.7e-9）、Coverage_Multiple FP=0.73 vs TP=0.49 → CNV amplicon artifact 驅動，已被 F pilot `NG=4+AF<0.4+NR≥80+NonLOH` canonical filter 覆蓋

## Consequences
- **根因**：Z3 本身已是最濃縮的 FP zone，內部特徵已無二階區分能力；甲基化特徵在 NG≤1 下因 HP 分化不足而結構性失效。
- **替代方向**：Zone-Aware Framework 定位不變（characterization only）；HCC1954 amplicon-FP 機制是唯一可深入的方向，但需 SEQC2 F1 panel 驗證，預期增益 <+0.005。
- **決策**：2026-04-18 關閉 Z3 二階分層 / AF germline separation 方向。

## References
- MEMORY：`project_z3_internal_exploration_negative.md`
- Report：`docs/experiments/in_progress/2026/04/20260418_Z3_Internal_Feature_Exploration_01.md`
- Scripts：`research/z3_internal_feature_exploration/scripts/step{1,2_5,3}_*.py`
