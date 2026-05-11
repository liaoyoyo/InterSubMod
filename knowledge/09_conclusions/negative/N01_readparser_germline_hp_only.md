---
id: ism-kb-09-conclusions-negative-01-readparser-germline-hp-only
name: "N01: ReadParser --germline-hp-only Phase 1 filter 方向 NEGATIVE"
description: "2026-04-21 Phase 0-1 完成；機制正確（12 unit tests pass）但 HCC1395 TO 全量 28,509 TP + 11,606 FP 無 TSV 特徵 AUC delta ≥+0.02；HP 衍生特徵反而下降（HPFineNGroups -0.026）。不進 Phase 2。"
status: archived
last_verified: 2026-04-21
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_readparser_germline_hp_only_phase1_negative.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, germline, hp-tag, readparser, phase1]
canonical_paths: [09_conclusions/negative/N01_readparser_germline_hp_only.md]
alias_paths: []
---

# N01: ReadParser --germline-hp-only Phase 1 filter 方向 CONDITIONAL NEGATIVE

## Status
`concluded-negative` (filter-level)；characterization 方向條件性保留。

## Context
假說：將 somatic HP tag (HP:i:11/21/33) 降級為純 germline HP 分群，可消除 self-phasing 帶來的 noisy HP 衍生特徵，進而提升 TP/FP 區分力。先前 `HPFineNGroups ≥4 + NR≥80` 曾顯示 TP rate 89.1%，但未驗證是否為 somatic tag 人工分組所致。

## Decision-tested
`--germline-hp-only` flag 開啟時，HP 衍生特徵（HPFineNGroups、HPMergedDelta、NHP3 等）的 AUC 是否提升 ≥+0.02，以支持 filter-level F1 改善。

## Method
- Phase 0 smoke：chr19 615 TP regions + 12 unit tests，commit `775036d` on `refactor/phase1-safety`
- Phase 1：HCC1395 TO 全量（28,509 TP + 11,606 FP），flag=on vs flag=off ISM 輸出比對
- Demotion 數學守恆 audit 獨立驗證

## Result
- 無任何 TSV 特徵 AUC delta ≥ +0.02
- HP 衍生特徵反而 AUC 下降：HPFineNGroups -0.026、HPMergedDelta -0.025、NHP3 -0.035
- flag=on 時 NGroups ≥3 regions **完全消失**（0/28,509）
- TP/FP 的 somatic tag 密度近乎相同（TP 27.4/site vs FP 27.7/site）
- Plan baseline HP_Ratio TP median=0.836 與實測 V3-Fixed TO 0.549 不符

## Consequences
- **根因**：somatic HP tag 本身無 TP/FP 訊號；修正只是把 noisy 分群變成無分群，不產生新訊號。HPFineNGroups subclone marker POSITIVE 結論可能為 somatic HP tag 的人工分組 artifact。
- **替代方向**：不進 Phase 2（7 樣本 × 2 模式全量重跑）；修正保留（default=off，研究者可啟用獲乾淨 HP 分群）。HPFineNGroups subclone-related 結論必須在 flag=on 下重驗證才算可信。
- **決策**：2026-04-21 Phase 1 後關閉 filter 路線。

## References
- MEMORY：`project_readparser_germline_hp_only_phase1_negative.md`
- Phase 1 report：`docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md`
- Smoke report：`docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_smoke_01.md`
- Commit：`775036d`
