---
id: ism-kb-01-project-overview-current-phase
name: "Current Phase Snapshot"
description: "當前 Phase 2 方向 A+D（Normal Methylation Reference + CN/Purity-aware）；HCC1395 pilot POSITIVE；主要阻塞為 haplotag + ISM 全量重跑。⚠ 此為快照，2 週有效。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "phase snapshot against docs/CURRENT_FOCUS.md"
related_ids:
  - ism-kb-01-project-overview-index
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-01-project-overview-breakthrough-strategy
tags: [overview, phase, current, snapshot]
canonical_paths: [01_project_overview/03_current_phase.md]
alias_paths: []
---

# Current Phase Snapshot

> ⚠️ **此為快照，有效期 2 週**（至 2026-05-06）。最新狀態請見 [docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)。

- 一句結論：Phase 2 方向 A+D 進行中；HCC1395 pilot POSITIVE (97.3% sig)；阻塞在 haplotag + ISM 全量重跑
- 適用對象：進度追蹤、下一步決策
- 可直接執行命令（驗證日期：2026-04-22）：`cat /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md | head -50`

---

## 當前 Phase：Phase 2 (方向 A+D)

### 方向 A：Normal Methylation Reference 整合
- **進度**：HCC1395 pilot 完成，97.3% significant regions
- **狀態**：🟢 POSITIVE signal
- **下一步**：擴展至 7 樣本

### 方向 D：CN/Purity-aware Correction
- **進度**：Theory framework 定義
- **狀態**：🟡 設計中
- **依賴**：方向 A 的 baseline

---

## 阻塞點

| ID | 阻塞 | 影響 | 解法 |
|----|------|------|------|
| B1 | Haplotag 重跑待執行 | 7 樣本 Phase 2A 無法啟動 | 需 PON-only phasing 執行 |
| B2 | expected_coverage hardcoded 75.0 bug | COLO829 等低覆蓋度樣本失真 | `/cpp-change` 修 KDE + 重跑 |
| B3 | COLO829 TO 無 methylation | TO pipeline 跨樣本不齊 | Accept limitation or skip TO for COLO829 |

---

## 最近 2 週 evidence_ledger 主題（2026-04-11~2026-04-21）

1. ✅ PON 跨樣本驗證（97.77%±1.12%）
2. ✅ H2009 負向根因診斷
3. ❌ 文獻假說 L1-L4 交叉驗證（mostly NEGATIVE）
4. ✅ Phase B/C/D Dual-BAM 驗證（HCC1395）
5. ✅ LOH × AF × 甲基化雙重證據（7/7）
6. ⚠️ Per-CpG ASM 與 Epiallele 指標（characterization positive）
7. ✅ F HPFineNGroups Part B（canonical filter 升級）
8. ⚠️ CovM Baseline 與 KDE 審計（partial/retraction）
9. ❌ ClairS-TO Verdict（F1 不可行；calibration positive）

---

## 活躍假說（priority=high）

| ID | 假說 | 狀態 |
|----|------|------|
| H011 | QS≥50 TO rescue | ✅ adopted |
| H012 | GQ≥3 filter | ✅ adopted |
| H_COMBO | 組合 filter | ✅ adopted |

完整：[../10_research_status/02_active_hypotheses.md](../10_research_status/02_active_hypotheses.md)

---

## 下一里程碑

- **M1**：7 樣本 × paired_full 重跑完成（等 haplotag + ISM）
- **M2**：Phase 2A 正式分析啟動
- **M3**：Normal Reference framework v1 完成

---

## 相關

- 最新權威：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
- 研究狀態：[../10_research_status/](../10_research_status/)
- 突破策略：[04_breakthrough_strategy.md](04_breakthrough_strategy.md)
