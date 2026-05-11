---
id: ism-kb-10-research-status-evidence-ledger-format
name: "Evidence Ledger Format"
description: "research/autoresearch/evidence_ledger.jsonl 的欄位規範、cycle 結構、最近 2 週摘要。"
status: active
last_verified: 2026-04-22
content_nature: reference
doc_type: reference
verified_scope: "ledger format against research/autoresearch/evidence_ledger.jsonl"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-active-hypotheses
  - ism-kb-09-conclusions-hypothesis-queue-snapshot
  - ism-kb-06-workflows-cpp-change-pdd
  - ism-kb-06-workflows-pptx-and-weekly-report
tags: [status, ledger, cycle, evidence, format]
canonical_paths: [10_research_status/03_evidence_ledger_format.md]
alias_paths: []
---

# Evidence Ledger Format

- 一句結論：evidence_ledger.jsonl 每行 1 cycle；11 個欄位；19 cycles 累積 2026-03-22 至 2026-04-21
- 適用對象：研究迴圈設計、cycle 追溯
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  head -1 /big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/evidence_ledger.jsonl | jq .
  ```

---

## Ledger 欄位規範（13 欄，v0.5+）

```json
{
  "cycle_id": "cycle_020",
  "hypothesis_id": "H_EXAMPLE_001",
  "hypothesis": "範例：某新假說",
  "pipeline_track": "paired_full",
  "result_f1": 0.9700,
  "delta_f1": 0.0050,
  "verdict": "POSITIVE",
  "human_decision": "adopt",
  "tier_used": "Tier A",
  "timestamp": "2026-04-24T10:00:00Z",
  "artifacts_path": "output/.../cycle_020/",
  "operator": "claude-opus-4.7-agent",
  "reviewer": "liaoyoyo"
}
```

---

## 欄位詳解

| 欄位 | 類型 | v0.5 新增 | 說明 |
|------|------|:---:|------|
| `cycle_id` | string | | 唯一 cycle 識別 |
| `hypothesis_id` | string | | 對應 hypothesis_queue.json 中的 id |
| `hypothesis` | string | | 假說一句話描述 |
| `pipeline_track` | enum | | paired_full / paired_pileup / to / paired_full_to |
| `result_f1` | float | | 絕對 F1 |
| `delta_f1` | float | | vs baseline 的 Δ |
| `verdict` | enum | | POSITIVE / NEGATIVE / CONDITIONAL_NEGATIVE / POSITIVE_CHARACTERIZATION |
| `human_decision` | string | | 用戶決策（adopt / reject / further_investigation） |
| `tier_used` | enum | | Tier A (canonical) / Tier B (synthesis) / Tier C (archive) |
| `timestamp` | ISO 8601 | | 記錄時間 |
| `artifacts_path` | string | | 產出 artifact 路徑 |
| `operator` | string | ✅ | **誰執行此 cycle**（AI agent ID / 研究者名）— ELN accountability |
| `reviewer` | string | ✅ | **誰 review verdict**（通常 PI 或研究者本人）— ELN accountability |

### v0.5 相容性
- **既有 19 cycles（2026-03-22 至 2026-04-21）**：豁免 `operator` / `reviewer` 欄位；不 backfill
- **新 cycles（2026-04-24+）必填兩欄**；`operator` 常見值：
  - `claude-opus-4.7-agent`（AI 主導 cycle）
  - `liaoyoyo`（人類執行）
  - `cpp-change-skill`（自動化 PDD skill 執行）
- `reviewer` 常見值：
  - `liaoyoyo`（PI / 研究者自 review）
  - `pending`（待 review）

**業界對照**：對齊 ELN（Electronic Lab Notebook）與 NIH reproducibility accountability 要求 — 每個實驗紀錄應含「誰執行」+「誰 review」兩欄位。

---

## 最近 2 週摘要（2026-04-11~2026-04-21）

| Cycle | 日期 | 假說 | Verdict | Track |
|-------|------|------|---------|-------|
| cycle_011 | 2026-04-11 | PON 跨樣本驗證 | ✅ POSITIVE (97.77%) | all |
| cycle_012 | 2026-04-11 | H2009 負向根因 | ✅ POSITIVE | paired_full |
| cycle_013 | 2026-04-13 | Phase B/C/D Dual-BAM HCC1395 | ✅ POSITIVE | paired_full |
| cycle_014 | 2026-04-14 | LOH × AF × Methylation | ✅ POSITIVE (7/7) | all |
| cycle_015 | 2026-04-15 | Per-CpG ASM | ⚠ CHARACTERIZATION | paired_full |
| cycle_016 | 2026-04-18 | HPFineNGroups Part B canonical | ✅ POSITIVE (+3.7pp) | paired_full |
| cycle_017 | 2026-04-19-21 | CovM Baseline + KDE | ⚠ PARTIAL | all |
| cycle_018 | 2026-04-20 | ClairS-TO Verdict | ❌ NEGATIVE on F1 | to |
| cycle_019 | 2026-04-21 | --germline-hp-only Phase 1 | ❌ CONDITIONAL_NEGATIVE | to |

---

## 查詢範例

```bash
# 所有 POSITIVE cycle
jq 'select(.verdict == "POSITIVE")' research/autoresearch/evidence_ledger.jsonl

# 最近 5 cycle
tail -5 research/autoresearch/evidence_ledger.jsonl | jq .

# 按 pipeline_track 分類
jq -s 'group_by(.pipeline_track) | map({track: .[0].pipeline_track, count: length})' \
  research/autoresearch/evidence_ledger.jsonl
```

---

## 相關

- Hypothesis queue：[../09_conclusions/05_hypothesis_queue_snapshot.md](../09_conclusions/05_hypothesis_queue_snapshot.md)
- Active 假說：[02_active_hypotheses.md](02_active_hypotheses.md)
- 原始檔：`research/autoresearch/evidence_ledger.jsonl`
