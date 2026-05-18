---
id: ism-kb-09-conclusions-hypothesis-queue-snapshot
name: "Hypothesis Queue Snapshot"
description: "research/autoresearch/hypothesis_queue.json 結構與最新快照；adopted/rejected/conditional/pre-registered 分類；含 2026-05-17 Pre-registration 強制機制。⚠ 每 2 週更新。"
status: active
last_verified: 2026-05-18
content_nature: reference
doc_type: reference
verified_scope: "hypothesis queue + Pre-registration mechanism against docs/CURRENT_FOCUS.md §2026-05-17 + /scientific-rigor §7.1"
related_ids:
  - ism-kb-09-conclusions-index
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-10-research-status-active-hypotheses
  - ism-kb-10-research-status-current-focus-snapshot
tags: [conclusions, hypothesis, queue, snapshot, pre-registration, thread-d, v6-production]
canonical_paths: [09_conclusions/05_hypothesis_queue_snapshot.md]
alias_paths: []
---

# Hypothesis Queue Snapshot

> ⚠️ **此為 2026-05-18 快照（基於 2026-05-17 plan tender-pondering-blossom）**，2 週有效至 2026-06-01。
> 最新狀態：`research/autoresearch/hypothesis_queue.json` + `docs/CURRENT_FOCUS.md §2026-05-17`
>
> **本快照已從 2026-04-22 priority-based queue 鏡像深度更新到 Tier 1-4 雙軌（Thread D paper + V6 production）+ Pre-registration 強制機制階層。**

- 一句結論：95+ 假說記錄於 queue；當前 active focus = H_THREAD_D_MAIN（⭐3 PARTIAL POSITIVE）+ H_V6_PROD（T1.2 Hard Gate）+ 4 baseline adopted；新增 Pre-registration 強制機制（依 /scientific-rigor §7.1）
- 適用對象：AI 研究迴圈、假說選擇、證據追蹤、Pre-registration 條件檢核
- 可直接執行命令：
  ```bash
  jq '.hypotheses | length' /big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/hypothesis_queue.json
  jq '.hypotheses | group_by(.status) | map({status: .[0].status, count: length})' \
    /big7_disk/liaoyoyo2001/InterSubMod/research/autoresearch/hypothesis_queue.json
  ```

---

## §1 Queue 結構

`research/autoresearch/hypothesis_queue.json` 含欄位：

```json
{
  "hypotheses": [
    {
      "id": "H_THREAD_D_MAIN",
      "name": "TP-enriched phasing signatures (LOH × cross_het)",
      "priority": 99,
      "status": "conditional_positive",
      "evidence_tier": "L2",
      "pipeline_track": "thread_d_paper",
      "pre_registration": {
        "H_predict": "cross_het bucket 在 ≥4 樣本同方向 TP-enrichment",
        "falsification": "Wilcoxon p > 0.1 across samples OR Δ < +0.005",
        "decision_threshold": "T2.1 + T2.3 通過 → ⭐4 升級"
      },
      "created": "2026-05-13",
      "last_updated": "2026-05-17",
      "cycle_ids": ["multi_agent_run_20260515"],
      "result_summary": "⭐3 PARTIAL POSITIVE; chr8 H4 LR deviance CN 0.211 > HP 0.063"
    }
  ]
}
```

> **2026-05-17 新增欄位**（依 /scientific-rigor §7.1 + AGENTS.md §3 補強）：
> - `evidence_tier`：L1-L5（依 /scientific-rigor §2 證據分級）
> - `pre_registration`：3 欄強制（H_predict / falsification / decision_threshold）
> - `pipeline_track`：thread_d_paper / v6_production / baseline_supporting

---

## §2 狀態分類（含 2026-05-17 新增）

| 狀態 | 含義 | 範例 ID |
|------|------|---------|
| `adopted` | 驗證通過，納入 pipeline | H011, H012, H_COMBO, H_KDE_001 |
| `conditional_positive` | 部分驗證，需更多證據 | **H_THREAD_D_MAIN（新）**、LOH-constrained phasing NG=2、LOH × AF × Methylation Paired |
| `pre_registered_pending` ⭐ NEW | Pre-registration 註冊但未跑 | **H_V6_PROD**、H_Z_AUTO_RECUR、Tier 2-4 開跑前候選 |
| `rejected` | 驗證失敗，不再追 | H001-H010（部分）、H003、H007 |
| `concluded_negative` | 已完整 postmortem | germline-hp-only Phase 1、TO germline FP、Q5 MCP 僵屍誤判 |
| `annotation_only` | 僅作特徵標註，不 filter | — |
| `active` | 進行中 | — |

---

## §3 當前主軸假說（2026-05-18 Tier 1-4）

### High Priority — Thread D paper Track

| ID | 名稱 | Priority | Status | Evidence Tier |
|----|------|---------|--------|---------------|
| H_THREAD_D_MAIN | TP-enriched phasing signatures | 99 | conditional_positive | L2 (⭐3) |
| H4 | HCC1395 chr8 hotspot CN+AF dominant | 97 | conditional_positive | L2 |
| H_Z_AUTO_RECUR | Z-AUTO KDE cross-sample replication | 96 | pre_registered_pending | L4 |

### High Priority — V6 production Track

| ID | 名稱 | Priority | Status | Gate |
|----|------|---------|--------|------|
| H_V6_PROD | V6 binary production tag finalize | 95 | pre_registered_pending | 🔴 Hard Gate × 2 (T1.2) |

---

## §4 Adopted 假說（baseline supporting layer）

> 自 2026-04-22 以來保留但 priority 已降至 baseline；雙軌主軸下不再 active 推進。

- **H011**：QS≥50 TO rescue (+0.008556)
- **H012**：GQ≥3 (+0.009365)
- **H_COMBO**：組合 filter
- **H_KDE_001**：KDE audit logging

---

## §5 Rejected / Concluded Negative 假說（防重複調查）

> Cross-reference 至 memory「Concluded」區。

| ID | 名稱 | 理由 | Memory |
|----|------|------|--------|
| H001 | HPP bias | rejected | — |
| H003/H007 | (TO haplotagging) | TO haplotagging 不可靠 | — |
| H005 | VAF/AlleleDelta in TO | rejected | — |
| H008/H009 | (CramersV track-dep) | rejected | — |
| — | Read-Level Germline FP Phase 1 | LOSO AUC 0.721 但 FP removal=0% | `project_readparser_germline_hp_only_phase1_negative` |
| — | TO Germline FP | G1-G7 全 AUC<0.64 | `project_germline_fp_identification_nogo` |
| — | Q5 biorxiv/ensembl MCP「僵屍」誤判 | 實測非僵屍 | `feedback_researcher_claim_needs_empirical_verification`（commit `f3611a7`） |

---

## §6 Pre-registration 強制機制（依 /scientific-rigor §7.1，2026-05-17 新增）

### 適用範圍
所有 Tier 2-4 新研究方向開跑前須在 `research/<topic>/00_INDEX.md` 強制 3 欄。

### 3 欄定義
| 欄位 | 內容 | 範例 |
|------|------|------|
| H 預測 | confirmatory 假設 | 「cross_het bucket 在 ≥4 樣本同方向 TP-enrichment」 |
| 否證條件 | NO-GO threshold | 「Wilcoxon p > 0.1 OR Δ < +0.005」 |
| Decision threshold | 升 ⭐4 / 廢棄 / 推進 paper | 「通過 → ⭐4 升級」 |

### Hard Gate 約束
達 NO-GO 條件 → verdict 不可事後改寫（依 AGENTS.md §3 + CLAUDE.md §1 Hard Gate）。

### 範本
`InterSubMod/templates/research_index.md`（commit `a031d21`）

### 區分
| 類型 | 強制 Pre-registration |
|------|---------------------|
| Confirmatory | 強制（達 NO-GO 不可改寫） |
| Exploratory | 建議（不強制） |

---

## §7 優先級分佈

| Priority | 意義 | 當前範例 |
|----------|------|---------|
| 95-99 | 雙軌主軸（最高優先） | H_THREAD_D_MAIN, H4, H_V6_PROD, H_Z_AUTO_RECUR |
| 85-95 | 高優先（baseline / 已 rejected） | H011, H012, H003/H007 (rejected) |
| 60-80 | 中優先 | — |
| <60 | 低優先 / annotation | — |

---

## §8 Cycles 對應

`research/autoresearch/evidence_ledger.jsonl` 記錄每個 cycle + 多 multi-agent runs：

| Cycle / Run | 日期 | 主題 | 對應假說 |
|-------------|------|------|---------|
| 19 cycles | 2026-03-22 ~ 04-21 | PON / KDE / germline-hp-only 等 | H011, H012, baseline |
| multi_agent_A-E + coord | 2026-05-15 | V6 TPFP HP LOH CN characterization | H_THREAD_D_MAIN, H4 |
| T1.1 + T1.3 | 2026-05-16 | Thread D 主軸正名 + scaffolding | infra |
| governance v3 + skill | 2026-05-17 | D2 分流 + /scientific-rigor + 8 cross-ref | infra |
| T1.2 V6 prod tag | ⏳ | 5-day Hard Gate workflow | H_V6_PROD |

詳見 [../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)

---

## §9 如何查詢

```bash
# 全 queue 摘要
jq '.hypotheses | group_by(.status) | map({status: .[0].status, count: length})' \
  research/autoresearch/hypothesis_queue.json

# High priority Tier 1-4
jq '.hypotheses[] | select(.priority >= 95) | {id, name, priority, status, evidence_tier}' \
  research/autoresearch/hypothesis_queue.json

# Pre-registered pending
jq '.hypotheses[] | select(.status == "pre_registered_pending") | {id, pre_registration}' \
  research/autoresearch/hypothesis_queue.json

# 按 cycle 查
jq '.[] | select(.cycle_id == "multi_agent_run_20260515")' \
  research/autoresearch/evidence_ledger.jsonl
```

---

## §10 相關

- Active hypotheses：[../10_research_status/02_active_hypotheses.md](../10_research_status/02_active_hypotheses.md)
- Current focus：[../10_research_status/01_current_focus_snapshot.md](../10_research_status/01_current_focus_snapshot.md)
- Evidence ledger format：[../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)
- 原始檔：`research/autoresearch/hypothesis_queue.json`
- Pre-registration 範本：`InterSubMod/templates/research_index.md`
- /scientific-rigor §7.1：`InterSubMod/.claude/skills/scientific-rigor/SKILL.md`
