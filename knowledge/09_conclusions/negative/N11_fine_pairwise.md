---
id: ism-kb-09-conclusions-negative-11-fine-pairwise
name: "N11: Fine-Pairwise Distance NEGATIVE — ISM 特徵空間耗盡"
description: "2026-04-08 748,391 regions × 7 samples × 2 modes；Paired 6 pairwise 距離全 AUC<0.50（反轉，germline ASM > somatic ASM）；TO 最高 0.579；LOH 層 AUC=0.132（極端反轉）；self-phasing 指紋確認。"
status: archived
last_verified: 2026-04-08
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_fine_pairwise_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, fine-pairwise, methylation, self-phasing, feature-exhaustion]
canonical_paths: [09_conclusions/negative/N11_fine_pairwise.md]
alias_paths: []
---

# N11: Fine-Pairwise Distance NEGATIVE

## Status
`concluded-negative`（ISM 甲基化特徵空間正式耗盡）。

## Context
HP 四群組（HP1、HP2、HP1-1、HP2-2）6 種 pairwise 距離（d(HP1,HP2)、d(HP1-1,HP2) 等）為 ISM 最後一組尚未暴露的內部特徵。假說：within/cross-haplotype methylation distance 是否有 TP/FP 信號。

## Decision-tested
6 pairwise 距離中是否任一在 paired 或 TO 模式達 AUC ≥ 0.58。

## Method
- 748,391 regions (Paired 328,699 + TO 419,692)
- 7 samples × 2 modes
- POOLED + LOH 層 + per-sample 多層驗證
- 多 Agent 驗證（數據抽驗 + 一致性檢查）

## Result
- **Paired POOLED**：全部 6 距離 AUC < 0.50（反轉），最低 d(HP1-1,HP2) = 0.347
- **TO POOLED**：最高 d(HP1,HP2) = 0.579（<0.58 門檻）
- **LOH 層 Paired**：d(HP1-1,HP2) AUC = 0.132（極端反轉）
- NaN 率：Paired 61-87%, TO 64-76%（需四群組各≥3 reads）
- **Self-phasing 指紋**：Paired HP1-1 (somatic) 65% >> HP1 (germline) 39%（正確分離）；TO HP1 64% ≈ HP1-1 61%（群組邊界模糊）
- 多 Agent 驗證 6/6 一致性通過

## Consequences
- **根因（Paired 反轉機制）**：FP (germline variants) 標記 cis-regulatory mQTL → inter-haplotype 甲基化差異天然大於 TP 的 stochastic somatic ASM。與 SNV-ASM 驗證的 FP>>TP 結論（stringent ASM 3.5×）完全一致。
- **替代方向**：**ISM 甲基化特徵空間已耗盡**。Paired 改進依賴 Phase 1A multi-bio (+0.0112 F1)；TO 只有 PON-only phasing 能重啟。
- **決策**：2026-04-08 不再探索新 ISM 甲基化距離特徵；C++ 10 新 CSV 欄位永久保留供 characterization。

## References
- MEMORY：`project_fine_pairwise_negative.md`
- Report：`docs/experiments/finalized/2026/04/20260408_fine_pairwise_distance_analysis_01.md`
- Data：`.../20260408_fine_pairwise_analysis/`
- C++：`RegionProcessor.hpp/cpp`（10 新 CSV 欄位）
