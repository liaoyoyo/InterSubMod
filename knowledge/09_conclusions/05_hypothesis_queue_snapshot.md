---
id: ism-kb-09-conclusions-hypothesis-queue-snapshot
name: "Hypothesis Queue Snapshot"
description: "research/autoresearch/hypothesis_queue.json 結構與最新快照（2026-06-15 實讀：23 total / 5 pending H019-H023 全 subclonal-axis 對齊 / 7 rejected / 2 adopted / 3 adopted_annotation / 6 closed）。⚠ JSON 為真值 SoT。"
status: active
last_verified: 2026-06-15
content_nature: reference
doc_type: reference
verified_scope: "hypothesis_queue.json status counts read 2026-06-15 + alignment to subclonal axis (CURRENT_FOCUS pinned)"
related_ids:
  - ism-kb-09-conclusions-index
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-10-research-status-active-hypotheses
  - ism-kb-10-research-status-current-focus-snapshot
tags: [conclusions, hypothesis, queue, snapshot, subclonal-reconstruction]
canonical_paths: [09_conclusions/05_hypothesis_queue_snapshot.md]
alias_paths: []
---

> 🔄 **2026-06-15 re-sync**（前版凍在 2026-05-18 Thread-D era）。**真值 SoT = `research/autoresearch/hypothesis_queue.json`**；本檔為結構鏡像。

# Hypothesis Queue Snapshot

> 2026-06-15 實讀 hypothesis_queue.json（非憑記憶）：**23 total**。

## 狀態分布（2026-06-15 實讀）

| status | 數 | 說明 |
|--------|---:|------|
| pending | 5 | H019-H023（全與 subclonal/LOH/ASM/CN 軸對齊，見下）|
| rejected | 7 | 含 H013-H018（filter/caller-F1，全追溯 concluded-DEAD，2026-06 `/pivot-direction` 降權）|
| adopted | 2 | H011 (QS≥50 TO rescue) · H012 (GQ≥3 TO rescue) |
| adopted_annotation | 3 | H008 (PairwiseMedianDist) · H009 (CramersV) · H010 (hp_assign_rate) |
| closed | 6 | 早期已結 |

## 5 pending（與 subclonal 主軸對齊）

| ID | P | 主題 |
|----|--:|------|
| H019 | 82 | longphase-S 在 SEQC2 LOH.bed 標定 LOH 區仍輸出 HP（phasing 行為）|
| H020 | 72 | LOH-pure 區甲基分群分得開且 read 不少 → 驗分離機制（真 ASM？）|
| H021 | 70 | 用 CN 判讀 ploidy 驗/反推 purity |
| H022 | 80 | ASM 訊號漏斗篩選是否夠清楚（審查篩選標準 + 位點列表）|
| H023 | 75 | ASM 漏斗放寬閾值後位點量變化 + TP/FP/FN 量化 |

> 這 5 條皆服務當前主軸（subclonal reconstruction 的 LOH/ASM/CN/甲基 characterization 面）。

## ❌ rejected（勿重開）

H013-H018 = methyl filter / caller-F1-headroom / methylation-as-FP-marker / methylation→TP enrichment — 全追溯 concluded-DEAD（Phase2 Cycle3 P6 NEGATIVE_filter_direction_failed）。

## 相關

- 真值 SoT：`research/autoresearch/hypothesis_queue.json`
- live 主軸：[../10_research_status/01_current_focus_snapshot.md](../10_research_status/01_current_focus_snapshot.md) → `InterSubMod/docs/CURRENT_FOCUS.md`
- 活躍假說：[../10_research_status/02_active_hypotheses.md](../10_research_status/02_active_hypotheses.md)
