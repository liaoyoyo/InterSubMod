---
id: ism-kb-10-research-status-next-milestones
name: "Next Milestones"
description: "Phase 2 下一里程碑：M1 7 樣本 paired_full 重跑、M2 Phase 2A 正式分析、M3 Normal Reference v1。⚠ 2 週有效。"
status: active
last_verified: 2026-05-18
content_nature: runtime-fact
doc_type: reference
verified_scope: "milestones against docs/CURRENT_FOCUS.md phase 2 roadmap"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-blockers-and-risks
  - ism-kb-01-project-overview-breakthrough-strategy
tags: [status, milestones, roadmap, phase-2]
canonical_paths: [10_research_status/05_next_milestones.md]
alias_paths: []
---

# Next Milestones

> ⚠️ **此為 2026-04-22 roadmap**；每 2 週檢視
>
> **📅 2026-05-18 Update Notice**: roadmap 已演進（詳見 `docs/CURRENT_FOCUS.md §2026-05-17 Tier 1-4`）：
> - **Tier 1 (W3 2026-05-15~22)**: T1.1 ✅ / T1.2 🔴 / T1.3 ✅
> - **Tier 2 (W3-W4)**: Z-AUTO KDE 4 樣本擴展 + HCC1395 primary 章節 + 6 樣本 replication
> - **Tier 3 (W4-W6)**: Paper outline + GitHub + Docker + Benchmark suite
> - **Tier 4 (W6+)**: Phase 2A Normal BAM cross-sample, PI Report errata
>
> 本快照僅作 2026-04-22 milestone 鏡像；下次完整內容更新待後續 session。

- 一句結論：3 個短期里程碑（M1-M3）聚焦 Phase 2A 7 樣本全量驗證
- 適用對象：中期規劃、進度預估
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  tail -30 /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  ```

---

## 短期里程碑

### M1：7 樣本 paired_full 重跑完成
- **預估**：2026-05 初
- **依賴**：PON-only phasing + haplotag 執行
- **驗收**：
  - 7 樣本 canonical run 齊全
  - significance_summary 59 欄完整
  - F1 與 baseline 差異 < 0.01
- **優先級**：🔴 高

### M2：Phase 2A 正式分析啟動
- **預估**：M1 完成後 2 週
- **目標**：Normal BAM + ISM 整合跨 7 樣本驗證
- **驗收**：
  - HCC1395 pilot 97.3% sig 是否複製到其他 6 樣本
  - Cross-sample direction consistency ≥6/7

### M3：Normal Methylation Reference Framework v1
- **預估**：2026-06
- **目標**：Phase A+D 整合文件
- **驗收**：
  - Framework 理論 + 實作完成
  - 7 樣本實證
  - Characterization + 論文素材

---

## 中期里程碑

### M4：Zone-Aware 全量驗證
- **預估**：2026-05-06（Phase B 並行）

### M5：HPFineNGroups Part B flag=on 重驗
- **預估**：2026-05（下一 cycle）

### M6：論文 draft v0
- **預估**：2026-07（取決於 M3 完成度）

---

## 已完成里程碑（參考）

- ✅ 2026-04-17：Phase 1A locked (+0.0112 paired_full)
- ✅ 2026-04-18：HPFineNGroups Part B canonical filter 升級
- ✅ 2026-04-13：HCC1395 Phase 2 pilot 97.3% sig
- ✅ 2026-04-11：PON 97.77% 跨樣本驗證

---

## 時間軸（視覺化）

```
2026-04 ─┬─ Phase 1A locked         ✅
         ├─ HPFineNGroups Part B    ✅
         ├─ Phase 2 pilot HCC1395   ✅
         │
2026-05 ─┼─ M1: 7 樣本 paired_full  ⏳ (blocked on haplotag)
         ├─ M2: Phase 2A 正式分析   ⏳
         ├─ M4: Zone-Aware 驗證     ⏳
         ├─ M5: HPFineNGroups flag  ⏳
         │
2026-06 ─┼─ M3: Normal Ref v1        ⏳
         │
2026-07 ─┼─ M6: 論文 draft v0        ⏳
```

---

## 阻塞解除優先級

若資源有限：
1. 先解 B1（haplotag 重跑）→ 解鎖 M1/M2/M4
2. 並行解 B2（expected_coverage bug）→ 品質確保
3. B3（COLO829 TO）可 accept limitation 先不處理

---

## 相關

- Current Focus：[01_current_focus_snapshot.md](01_current_focus_snapshot.md)
- Blockers：[04_blockers_and_risks.md](04_blockers_and_risks.md)
- Breakthrough Strategy：[../01_project_overview/04_breakthrough_strategy.md](../01_project_overview/04_breakthrough_strategy.md)
