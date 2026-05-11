---
id: ism-kb-09-conclusions-negative-13-option-c-dual-path
name: "N13: Option C 雙路（HP-free clustering）NEGATIVE"
description: "2026-04-07 Python prototype 測試；HP-free ClusterPermanovaF AUC=0.512；5 HP-free features combo 0.564；HP-dependent combo 0.598；全部組合 0.607；C++ 修改取消。"
status: archived
last_verified: 2026-04-07
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_option_c_dual_path_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, clustering, hp-free, option-c, self-phasing]
canonical_paths: [09_conclusions/negative/N13_option_c_dual_path.md]
alias_paths: []
---

# N13: Option C 雙路測試 NEGATIVE

## Status
`concluded-negative`（C++ C-3/C-4 修改取消）。

## Context
Option C 假說：LabelTest 雙路計算（同時以 HP-based 和 HP-free 特徵分析）可在 TO self-phasing 汙染 HP 標籤的情況下，提供乾淨的 methylation-only cluster 品質特徵。用戶指示先做 Python prototype 驗證再修改 C++。

## Decision-tested
HP-free methylation clustering features（ClusterPermanovaF、PairwiseMeanDist 等）在 TO 模式是否有 TP/FP 區分力。

## Method
- Python prototype（TO valid subset, n=90,572）
- HP-free 5 features combo vs HP-dependent 2 features combo vs 全部組合
- 交互作用項測試

## Result
- **架構發現**：`cluster_labels` 已經是 HP-free（來自甲基化 hierarchical clustering via TreeCutter）；ClusterPermanova 被 `passed_gate` 擋住，TO 僅 22% 有效
- **ClusterPermanovaF AUC = 0.512**（純甲基化 cluster 品質完全無區分力）
- **HP-free 5 features combo AUC = 0.564**
- **HP-dependent combo AUC = 0.598**
- **全部組合 AUC = 0.607** — HP-free 只增加 +0.009
- 交互作用項無幫助（0.608 vs 0.607）

## Consequences
- **根因**：純甲基化 clustering 無區分力；germline ASM 創造的 cluster structure 與 somatic heterogeneity 不可區分（identifiability problem，見 N15）。所有信號來自 HP tags，TO 模式被 self-phasing 汙染。
- **替代方向**：正確策略確認為 **PON-only phasing 上游修正**（非 ISM 層級雙路計算）。
- **決策**：2026-04-07 C-3/C-4 取消，不修改 C++ 代碼。

## References
- MEMORY：`project_option_c_dual_path_negative.md`
