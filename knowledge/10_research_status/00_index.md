---
id: ism-kb-10-research-status-index
name: "Research Status 索引"
description: "當前研究狀態（時間敏感）索引：CURRENT_FOCUS 快照、active 假說、evidence ledger 格式、阻塞、里程碑。⚠ 2 週有效。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "research status directory"
related_ids:
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-10-research-status-active-hypotheses
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-10-research-status-blockers-and-risks
  - ism-kb-10-research-status-next-milestones
tags: [status, index, current-focus, active]
canonical_paths: [10_research_status/00_index.md]
alias_paths: []
---

# Research Status 索引

> ⚠️ **本目錄全部文件時間敏感，2 週有效**。最新以 `docs/CURRENT_FOCUS.md` 為準。

- 一句結論：當前進度、活躍假說、阻塞、里程碑；查最新狀態從此開始
- 適用對象：每日站會、決策前、AI agent 回答「現在在做什麼」
- 可直接執行命令（驗證日期：2026-04-22）：`cat /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md | head -30`

---

## 文件列表

| 檔案 | 主題 | 更新頻率 |
|------|------|---------|
| [01_current_focus_snapshot.md](01_current_focus_snapshot.md) | Phase 2 進度快照 | 每週 |
| [02_active_hypotheses.md](02_active_hypotheses.md) | priority=high 假說清單 | 每 2 週 |
| [03_evidence_ledger_format.md](03_evidence_ledger_format.md) | ledger 欄位規範 + 最近 2 週摘要 | 每 2 週 |
| [04_blockers_and_risks.md](04_blockers_and_risks.md) | 當前阻塞與風險 | 每週 |
| [05_next_milestones.md](05_next_milestones.md) | 下一里程碑 | 每 2 週 |

---

## 權威來源優先級

1. **docs/CURRENT_FOCUS.md**（時間最新，隨時更新）
2. **research/autoresearch/evidence_ledger.jsonl**（每 cycle 追加）
3. **research/autoresearch/hypothesis_queue.json**（假說狀態）
4. **本 KB 快照**（每 2 週更新）

**原則**：此 KB 快照是**便利查詢**用途；真正決策應讀最上層權威文件。

---

## 相關

- Project overview：[../01_project_overview/](../01_project_overview/)
- Conclusions：[../09_conclusions/](../09_conclusions/)
- 權威：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
