---
id: ism-kb-10-research-status-active-hypotheses
name: "Active Hypotheses"
description: "當前 priority=high 活躍假說清單：H011 QS≥50、H012 GQ≥3、H_COMBO 組合 filter 等 adopted 假說。⚠ 2 週有效。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "active hypotheses against research/autoresearch/hypothesis_queue.json"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-09-conclusions-hypothesis-queue-snapshot
  - ism-kb-10-research-status-current-focus-snapshot
tags: [status, hypotheses, active, priority-high]
canonical_paths: [10_research_status/02_active_hypotheses.md]
alias_paths: []
---

# Active Hypotheses

> ⚠️ **此為 2026-04-22 快照**，最新：`research/autoresearch/hypothesis_queue.json`

- 一句結論：當前 adopted 假說 4 個（H011/H012/H_COMBO/H_KDE_001）；rejected 假說不再追
- 適用對象：研究迴圈下一輪選擇、假說優先級決策
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  jq '.hypotheses[] | select(.status == "adopted")' \
    /big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/hypothesis_queue.json
  ```

---

## Adopted 假說（當前 pipeline 內）

### H011：QS≥50 TO rescue
- **Priority**：95
- **Track**：to_pileup
- **Delta F1**：+0.008556
- **狀態**：adopted

### H012：GQ≥3
- **Priority**：~85
- **Delta F1**：+0.009365
- **狀態**：adopted

### H_COMBO：組合 filter
- **狀態**：adopted

### H_KDE_001：KDE audit logging
- **狀態**：adopted（infra）

---

## Rejected 假說（別重跑）

| ID | 名稱 | Priority | 理由 |
|----|------|---------|------|
| H001 | HPP bias | — | rejected |
| H003 | — | 85 | TO haplotagging 不可靠 |
| H005 | VAF/AlleleDelta in TO | — | rejected |
| H007 | — | 85 | 同 H003 |
| H008 | — | — | rejected |
| H009 | CramersV track-dep | — | rejected |

---

## Conditional Positive / Annotation Only

（若有，在此列出；目前主要以 adopted + rejected 為主）

---

## Pending（尚未驗證）

（每 cycle 產生新候選；當前 queue 含 95+ 假說，細節見原始 JSON）

---

## 最近 Cycle 結果（2026-04-11~2026-04-21）

從 `research/autoresearch/evidence_ledger.jsonl`：

| Cycle | 日期 | 主題 | Verdict |
|-------|------|------|---------|
| cycle_011 | 2026-04-11 | PON 驗證 | ✅ |
| cycle_012 | 2026-04-11 | H2009 診斷 | ✅ |
| cycle_013 | 2026-04-13 | Phase B/C/D Dual-BAM | ✅ HCC1395 |
| cycle_014 | 2026-04-14 | LOH × AF | ✅ 7/7 |
| cycle_015 | 2026-04-15 | Per-CpG ASM | ⚠ characterization |
| cycle_016 | 2026-04-18 | HPFineNGroups Part B | ✅ |
| cycle_017 | 2026-04-19-21 | CovM/KDE | ⚠ partial |
| cycle_018 | 2026-04-20 | ClairS-TO Verdict | ❌ |
| cycle_019 | 2026-04-21 | germline-hp-only Phase 1 | ❌ CONDITIONAL |

---

## 下一候選

當前考量：
- 7 樣本 Phase 2A 全量驗證（等 haplotag）
- 2026-04-22 TO 跨樣本 archive Path 2（COLO829+DORADO+H2009 重跑 ISM）

---

## 相關

- Queue snapshot：[../09_conclusions/05_hypothesis_queue_snapshot.md](../09_conclusions/05_hypothesis_queue_snapshot.md)
- Evidence ledger：[03_evidence_ledger_format.md](03_evidence_ledger_format.md)
- 原始檔：`research/autoresearch/hypothesis_queue.json`
