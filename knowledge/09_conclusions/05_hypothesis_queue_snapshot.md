---
id: ism-kb-09-conclusions-hypothesis-queue-snapshot
name: "Hypothesis Queue Snapshot"
description: "research/autoresearch/hypothesis_queue.json 的結構與最新快照；adopted/rejected/conditional 分類。⚠ 每 2 週更新。"
status: active
last_verified: 2026-05-18
content_nature: reference
doc_type: reference
verified_scope: "hypothesis queue against research/autoresearch/hypothesis_queue.json"
related_ids:
  - ism-kb-09-conclusions-index
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-10-research-status-active-hypotheses
tags: [conclusions, hypothesis, queue, snapshot]
canonical_paths: [09_conclusions/05_hypothesis_queue_snapshot.md]
alias_paths: []
---

# Hypothesis Queue Snapshot

> ⚠️ **此為快照，2 週有效**。最新狀態：`research/autoresearch/hypothesis_queue.json`
>
> **📅 2026-05-18 Update Notice**: 新增假說與狀態變更（詳見 `docs/CURRENT_FOCUS.md §2026-05-17`）：
> - 新 category: Pre-registration confirmatory（依 `/scientific-rigor §7.1` 強制 3 欄 H 預測/否證/threshold）
> - HCC1395 chr8 hotspot ⭐3 PARTIAL POSITIVE（2026-05-15 multi-agent fan-out）
> - Q5 biorxiv/ensembl MCP 僵屍誤判 NEGATIVE→修正（已實測非僵屍）
>
> 本快照僅作 2026-04-22 queue 鏡像；下次完整內容更新待後續 session。

- 一句結論：95+ 假說記錄於 queue；狀態分 adopted / rejected / annotation_only / conditional_positive；高優先假說 priority=95
- 適用對象：AI 研究迴圈、假說選擇、證據追蹤
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  jq '.hypotheses | length' /big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/hypothesis_queue.json 2>/dev/null
  ```

---

## Queue 結構

`research/autoresearch/hypothesis_queue.json` 含欄位：

```json
{
  "hypotheses": [
    {
      "id": "H011",
      "name": "QS≥50 TO rescue",
      "priority": 95,
      "status": "adopted",
      "pipeline_track": "to_pileup",
      "created": "2026-03-22",
      "last_updated": "2026-04-15",
      "cycle_ids": ["cycle_001", "cycle_005"],
      "result_summary": "+0.008556 delta F1"
    },
    ...
  ]
}
```

---

## 狀態分類

| 狀態 | 含義 | 範例 ID |
|------|------|---------|
| `adopted` | 驗證通過，納入 pipeline | H011, H012, H_COMBO, H_KDE_001 |
| `rejected` | 驗證失敗，不進一步 | H001-H010（部分）、H003、H007 |
| `conditional_positive` | 部分驗證，需更多證據 | — |
| `annotation_only` | 僅作特徵標註，不 filter | — |
| `active` | 進行中 | — |

---

## 優先級分佈

| Priority | 意義 |
|----------|------|
| 95 | 當前最高優先 |
| 85 | 高優先，但已 rejected（H003, H007） |
| 60-80 | 中優先 |
| <60 | 低優先或 annotation |

---

## 當前 High Priority（priority ≥ 85）

| ID | 名稱 | 狀態 | 備註 |
|----|------|------|------|
| H011 | QS≥50 TO rescue | adopted | +0.008556 ΔF1 |
| H003 | HPP bias | rejected | TO haplotagging 不可靠 |
| H007 | — | rejected | 同上 |

---

## Adopted 假說（進 pipeline）

- **H011**：QS≥50 TO rescue (+0.008556)
- **H012**：GQ≥3 (+0.009365)
- **H_COMBO**：組合 filter
- **H_KDE_001**：KDE audit logging

---

## Cycles 對應

`research/autoresearch/evidence_ledger.jsonl` 記錄每個 cycle：
- 19 cycles (2026-03-22 至 2026-04-21)
- 欄位：`cycle_id`, `hypothesis_id`, `verdict`, `tier_used`, `artifacts_path`

詳見 [../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)

---

## 如何查詢

```bash
# 全 queue 摘要
jq '.hypotheses | group_by(.status) | map({status: .[0].status, count: length})' \
  research/autoresearch/hypothesis_queue.json

# High priority adopted
jq '.hypotheses[] | select(.priority >= 85 and .status == "adopted") | {id, name, priority}' \
  research/autoresearch/hypothesis_queue.json

# 按 cycle 查
jq '.[] | select(.cycle_id == "cycle_019")' \
  research/autoresearch/evidence_ledger.jsonl
```

---

## 相關

- Evidence ledger 格式：[../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)
- Active hypotheses：[../10_research_status/02_active_hypotheses.md](../10_research_status/02_active_hypotheses.md)
- 原始檔：`research/autoresearch/hypothesis_queue.json`
