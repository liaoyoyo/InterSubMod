---
id: ism-kb-09-conclusions-negative-15-feature-design-r1r5
name: "N15: ISM Feature Design R1-R5 — Identifiability 非設計問題"
description: "2026-04-07 R1-R5 系統研究；CramersV 93% 為零（2×2 限制，N≥3 失效但 HPFineNGroups 已克服 +0.125-0.139）；HPMergedDelta N=4 反向（均值壓縮）；結構清楚子集 AUC 全下降；Bootstrap 32-site 僅 29% p<0.05。"
status: archived
last_verified: 2026-04-07
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_feature_design_limitations_r1r5.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, feature-design, identifiability, cramersv, hpmergeddelta]
canonical_paths: [09_conclusions/negative/N15_feature_design_r1r5.md]
alias_paths: []
---

# N15: ISM Feature Design Limitations R1-R5

## Status
`concluded-negative`（設計修正不能解決，根因是 identifiability）。

## Context
用戶質疑 CramersV/HPMergedDelta 等特徵設計是否不足以捕捉 TP/FP 差異。R1-R5 系統研究需定量確認「是特徵設計問題還是根本不可分」。

## Decision-tested
修改特徵定義（擴展到 N≥3、多群設計、結構清楚子集、更嚴格統計）是否突破 CL-008 ≤0.58 ceiling。

## Method
- R1-R5 五個子研究
- CramersV 2×2 vs HPFineNGroups（擴展 contingency）
- HPMergedDelta 跨 N 群分析
- 結構清楚子集（CramersV>0, N≥3, 佔 8.7%）per-subset AUC
- Bootstrap 32-site p-value 分析
- Excess groups = HPFineN - expected(LOH) 跨子集比較

## Result
1. **CramersV 93% 為零**：2×2 contingency 框架在 N≥3 位點系統性失效。但 HPFineNGroups 已克服（AUC +0.125-0.139）→ 不需修改 C++ 核心
2. **HPMergedDelta 多群時反向**：N=4 時 TP<FP by 0.006 → 均值壓縮設計缺陷。HPFineNGroups 不受此影響
3. **結構清楚子集 AUC 全部下降** → FP 的 germline ASM 在結構清楚時更強
4. **Bootstrap 32-site**：HPFineNGroups 只有 29% 的抽樣達 p<0.05 → 小樣本肉眼檢視容易產生過度樂觀印象
5. **Excess groups**：跨子集統一 +0.059，子集內 = 常數偏移，無新信號

## Consequences
- **根因**：不是特徵設計問題，是 **identifiability 問題**。germline ASM 與 somatic heterogeneity 的信號來源混淆，更好的統計方法也無法解決。
- **替代方向**：不再投入特徵 re-design 投資；接受 CL-008 ≤0.58 ceiling；研究主軸改為 characterization（HPFineNGroups 可作 somatic heterogeneity marker 但非 filter）。
- **決策**：2026-04-07 特徵設計改進路線關閉。

## References
- MEMORY：`project_feature_design_limitations_r1r5.md`
