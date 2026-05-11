---
id: ism-kb-10-research-status-current-focus-snapshot
name: "CURRENT_FOCUS Snapshot"
description: "docs/CURRENT_FOCUS.md 的結構化鏡像；Phase 2 方向 A+D 進行中。⚠ 2 週有效。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "mirror of docs/CURRENT_FOCUS.md on 2026-04-22"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-active-hypotheses
  - ism-kb-10-research-status-blockers-and-risks
  - ism-kb-01-project-overview-current-phase
  - ism-kb-03-pipelines-f1-baseline-canonical
tags: [status, snapshot, current-focus, phase-2]
canonical_paths: [10_research_status/01_current_focus_snapshot.md]
alias_paths: []
---

# CURRENT_FOCUS Snapshot

> ⚠️ **此為 2026-04-22 快照，2 週有效至 2026-05-06**
> 最新：[docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)

- 一句結論：Phase 2 方向 A+D 進行中；HCC1395 pilot POSITIVE (97.3%)；7 樣本 paired_full haplotag 重跑待執行
- 適用對象：決策前快速了解現況
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  cat /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  ```

---

## 快照內容（2026-04-22）

### 當前 Phase
**Phase 2**：方向 A + D
- **A**：Normal Methylation Reference 整合
- **D**：CN/Purity-aware Correction

### 進度
| 項目 | 狀態 |
|------|------|
| Phase A-D Code | ✅ 完成 (2026-04-13) |
| HCC1395 全基因體驗證 | ✅ 完成 (2026-04-13) |
| 7 samples 全量驗證 | ⏳ 待 haplotag + ISM 重跑 |

### Phase 1A Locked（已凍結）
- paired-pure ΔF1 = **+0.0112**（locked）
- TO ΔF1 = **-0.0206**（NEGATIVE locked）
- **完整 provenance（CI、方法、樣本、限制）** → [../03_pipelines/05_f1_baseline_canonical.md](../03_pipelines/05_f1_baseline_canonical.md)（SoT）

---

## 突破策略順序

E → **A+D**（當前）→ B → C

詳見 [../01_project_overview/04_breakthrough_strategy.md](../01_project_overview/04_breakthrough_strategy.md)

---

## 活躍主題

| 主題 | 狀態 |
|------|------|
| Normal BAM 整合 | 🟢 POSITIVE (HCC1395) |
| HPFineNGroups Part B 全量驗證 | ✅ 完成（canonical filter 升級） |
| Zone-Aware Framework | 🟡 characterization only |
| LOH × AF × Methylation | ✅ 7/7 positive |
| KDE Fix 下游量化 | ⚠ PARTIAL (H-CN1 POSITIVE, COLO829 pending) |

---

## 使用情境

若使用者問「現在在做什麼？」：
1. 讀本文件（快速概覽）
2. 若快照日期 > 2 週 → 讀 `docs/CURRENT_FOCUS.md`
3. 需假說細節 → `02_active_hypotheses.md`
4. 需證據細節 → `03_evidence_ledger_format.md` + `research/autoresearch/evidence_ledger.jsonl`

---

## 相關

- 最新權威：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
- 活躍假說：[02_active_hypotheses.md](02_active_hypotheses.md)
- 阻塞：[04_blockers_and_risks.md](04_blockers_and_risks.md)
- Next milestones：[05_next_milestones.md](05_next_milestones.md)
