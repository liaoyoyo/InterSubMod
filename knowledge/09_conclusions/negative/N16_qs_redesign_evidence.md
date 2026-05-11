---
id: ism-kb-09-conclusions-negative-16-qs-redesign-evidence
name: "N16: Quality Score TO 模式完全失效（雙重驗證）"
description: "2026-04-01/04-05 雙重驗證；舊 TO QS AUC=0.497 + ClairS-TO 矯正後 0.507；TO FP median QS (85) > TP median (75) 方向反轉；LOH penalty 是 TO QS 失效主因；7-Feature LR AUC 0.627；AF 不可硬閾值。"
status: archived
last_verified: 2026-04-05
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_qs_redesign_evidence.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, quality-score, qs, loh-penalty, to-mode]
canonical_paths: [09_conclusions/negative/N16_qs_redesign_evidence.md]
alias_paths: []
---

# N16: Quality Score TO 模式完全失效

## Status
`concluded-negative`（TO QS 現有設計失效；需根本重設計）。

## Context
Quality Score (QS) 是 ISM 核心綜合分數。Paired 模式 AUC 0.754（有效），但 TO 模式疑似失效。2026-03-31 + 2026-04-05 兩輪驗證（後者用矯正後 ClairS-TO VCF）。

## Decision-tested
1. TO QS 是否保留區分力
2. LOH penalty / AF 過濾 / EHR Tier 等子組件是否有效

## Method
- Paired vs TO AUC 對比
- LOH penalty 觸發率與影響量化
- AF 硬閾值跨樣本驗證
- EHR Tier 分佈分析
- 7-Feature LR ensemble（2026-04-05 矯正後）

## Result
### 第一輪（2026-03-31）
- Paired QS AUC = 0.754；TO QS AUC = **0.497**（隨機）
- TO FP median QS (85) > TP median QS (75) — 方向反轉
- TO TP 44.5% 觸發 LOH penalty (mean -11.1)；FP 35.8% (mean -8.9)
- 移除 LOH penalty：TO TP median 75→100, FP 85→100
- TO caller_af AUC = 0.418（劣於隨機）
- TO per-sample TP purity 37-99.9%（HCC1954 災難性），無通用閾值
- TO EHR Tier 各 bin FP rate 平坦 30-37%

### 第二輪（2026-04-05 矯正後）
- **QS AUC = 0.507**（精確隨機）
- TP 91.2% High Tier, FP 90.4% High Tier — 完全無區分
- 7-Feature LR AUC = 0.627（比舊 0.760 大幅下降）
- AlleleDelta TP ≈ FP（germline FP 有真正 allele structure）
- LOH_S AUC 0.895→0.721, LOH_W 0.863→0.739

## Consequences
- **根因**：LOH penalty 在 TO 模式反向（懲罰 TP 多於 FP）；germline FP 在所有維度（QS、AF、Tier、AlleleDelta）與 TP 表現相同。
- **替代方向**：
  1. 立即：移除 TO LOH penalty（或反轉為 bonus）
  2. 短期：TO/Paired 分離 QS 權重配置
  3. 中期：多特徵 ML 模型（但無單一特徵 AUC > 0.6）
  4. AF 只能作 continuous feature，不可硬閾值
- **決策**：QS TO 根本重設計（方向：「消除反向特徵」優於「加入新正向特徵」）。

## References
- MEMORY：`project_qs_redesign_evidence.md`
- 分析：`observation_workspaces/20260331_loh_af_tier_qs_analysis/`
- 審查：`audit_loh_af_hypothesis.md`
- 矯正報告：`20260404_ClairSTO_VCF矯正全面重分析報告_01.md`
