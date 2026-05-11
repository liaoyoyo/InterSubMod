---
id: ism-kb-09-conclusions-negative-03-clairs-to-verdict
name: "N03: ClairS-TO Verdict Pilot — F1 增益 NEGATIVE"
description: "2026-04-20 HCC1395 t20_n30 pilot；Verdict 內部校準 POSITIVE（Germline FP 96.1%、Somatic TP 91.8%）但 Verdict_Germline 100% 落 LowQual 從不在 PASS，S1 ΔF1=0；主升級路徑改 Wakhan/SAVANA。"
status: archived
last_verified: 2026-04-20
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_clairs_to_verdict_pilot.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, clairs-to, verdict, external-cn, structural]
canonical_paths: [09_conclusions/negative/N03_clairs_to_verdict.md]
alias_paths: []
---

# N03: ClairS-TO Verdict Characterization Pilot NEGATIVE on F1

## Status
`concluded-negative` (F1 gain)；Verdict annotation 保留供未來 Wakhan/SAVANA 交叉驗證。

## Context
ClairS-TO v0.3.0 的 Verdict 模組（基於 ASCAT + binomial test）宣稱可分類 variant 為 Germline / Somatic / SubclonalSomatic。假說：將 Verdict_Germline 從 PASS 集合移除可提升 F1。

## Decision-tested
1. Verdict 內部校準是否達標（Germline FP rate ≥70%、Somatic TP rate ≥85%）
2. 將 Verdict 整合為 post-hoc filter 是否能提升 F1

## Method
- Sample：HCC1395 subsample t20_n30（purity ≈ 0.40，低純度人工稀釋）
- Step 1：Verdict vs SEQC2 truth 對比（characterization）
- Step 2：Reference F1 計算 — S1 nominal「remove Verdict_Germline from PASS」、S2 aggressive「只保留 PASS ∩ Verdict_Somatic∪Subclonal」

## Result
- **Verdict 內部校準 POSITIVE**：
  - Verdict_Germline FP rate = 96.1%（3,463/3,602）
  - Verdict_Somatic TP rate on PASS = 91.8%（7,785/8,479）
  - Verdict_SubclonalSomatic TP rate = 94.9%
- **F1 NEGATIVE**：
  - Verdict_Germline 4,633 個**全部落在 LowQual**，從不出現在 PASS
  - S1 ΔF1 = +0.0000
  - S2 ΔF1 = −0.2007（recall 0.521→0.210）
  - PASS 95.4% 的 FP（14,892/15,614）落在 PASS ∩ no_Verdict 子集

## Consequences
- **根因（結構性）**：ClairS-TO v0.3.0 的 PASS vs LowQual 判定已**前置嵌入** Verdict 所仰賴的 ASCAT / 二項式資訊，Verdict 與 LowQual 在資訊層面重疊。6 個 in-house 樣本 Verdict flags 全空（5 樣本缺 1000G loci 導致 ASCAT 失敗；HCC1395 purity >0.8 觸發 safety gate 跳過 Verdict）。
- **替代方向**：主升級路徑改為 **Wakhan**（haplotype-specific CN，與 HP-tagged ISM 對應最佳）；SAVANA 作 SV/CNA 互補層。Verdict 標籤保留為 per-variant annotation。
- **決策**：2026-04-20 不補齊 1000G loci、不改 purity safety gate、不整合 Verdict 為 C++ filter。

## References
- MEMORY：`project_clairs_to_verdict_pilot.md`
- Report：`docs/experiments/in_progress/2026/04/20260420_ClairS_TO_Verdict_Characterization_Pilot_01.md`
- Survey：`docs/references/2026/04/20260420_external_CN_tools_survey_01.md`
- 觸發 pilot：`docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md`
