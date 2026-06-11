---
id: ism-kb-10-research-status-next-milestones
name: "Next Milestones"
description: "Tier 1-4 序列化里程碑：T1.2 V6 production tag Hard Gate（W3）→ T2 證據強化（W3-W4）→ T3 paper draft（W4-W6）→ T4 reactive 擴展（W6+）。⚠ 2 週有效。"
status: active
last_verified: 2026-05-18
content_nature: runtime-fact
doc_type: reference
verified_scope: "milestones against docs/CURRENT_FOCUS.md §2026-05-17 Tier 1-4 plan tender-pondering-blossom"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-10-research-status-blockers-and-risks
  - ism-kb-01-project-overview-breakthrough-strategy
tags: [status, milestones, roadmap, tier-1-4, thread-d, v6-production]
canonical_paths: [10_research_status/05_next_milestones.md]
alias_paths: []
---

<!-- STALE-REDIRECT-BANNER (scripts/stale_redirect_banner.sh) -->
> ⚠ **此檔為 2026-05-18 前後快照，可能已過時** — 現役主軸/狀態以 `InterSubMod/docs/CURRENT_FOCUS.md` 為準（主軸已於 2026-06-11 pivot 至 Subclonal reconstruction（取代 G6））。本檔僅供歷史對照，勿據此判斷現況。


# Next Milestones

> ⚠️ **此為 2026-05-18 快照（基於 2026-05-17 plan tender-pondering-blossom Tier 1-4 序列化）**
> 每 2 週檢視；最新權威：`docs/CURRENT_FOCUS.md §2026-05-17`
>
> **本快照已從 2026-04-22 M1/M2/M3 鏡像深度更新到 Tier 1-4 雙軌（Thread D paper + V6 production）里程碑階層。**

- 一句結論：Tier 1 (W3) → Tier 2 (W3-W4) → Tier 3 (W4-W6) → Tier 4 (W6+) 序列化；當前最高優先 = T1.2 V6 production tag Hard Gate
- 適用對象：中期規劃、進度預估、resource allocation
- 可直接執行命令：
  ```bash
  sed -n '15,80p' /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  ```

---

## §1 Tier 1 (W3 2026-05-15~22) — 必須前置 🔴

### T1.1：Thread D 主軸正名 ✅ DONE 2026-05-16
- **產出**：`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_*.md` 338→381 行
- **變更**：加 banner + §2.5 paradigm reframe + 主標題改名「TP-enriched phasing signatures (LOH × cross_het)」
- **狀態**：✅ 完成

### T1.2：V6 Production Tag Finalize 🔴 Hard Gate 待執行
- **預估**：W3 結束（2026-05-22）
- **5-day Workflow**：
  - Day 1-2：COLO829 V6 ISM 補完（Archive TO rerun + KDE-corrected） 🟢
  - Day 3：7-sample marker coverage + caller F1 比較（V3F vs V5 vs V6 vs SEQC2） 🟢
  - Day 4：Binary commit hash 寫 `manifest.yaml` 🟡 + `git tag v6-prod-{YYYYMMDD}` 🔴 **Hard Gate**
  - Day 5：PI errata 5 條 + V6 sign-off written email draft 🟡 + User review → send 🔴 **Hard Gate**
- **解鎖**：T2.x Archive TO 7-sample rerun + T4.3 PI errata package
- **計劃**：`InterSubMod/research/selfphasing_v6_production/00_PLAN.md`

### T1.3：init-research scaffolding ✅ DONE 2026-05-16
- **產出**：
  - `InterSubMod/research/thread_d_paper/`（165 行 00_PLAN + 108 行 manifest）
  - `InterSubMod/research/selfphasing_v6_production/`（154 行 00_PLAN + 122 行 manifest）
- **狀態**：✅ 完成

---

## §2 Tier 2 (W3-W4) — 證據強化 🟠

### T2.1：Z-AUTO KDE 跨 4 樣本擴展 ⏳
- **預估**：W4 中
- **依賴**：T1.2 V6 binary tag 完成
- **目標**：4 樣本 (H1437/H2009/HCC1954/HCC1937) KDE 各自做
- **驗收**：⭐3 → ⭐4 升級必要條件 — Z-AUTO mechanism 在 ≥3/4 樣本 recur

### T2.2：HCC1395 primary discovery 章節骨架 ⏳
- **預估**：W4
- **目標**：Strategy A §3「primary discovery」章節（HCC1395 chr8 hotspot + paradigm reframe + V5 over-promote 證據）
- **依賴**：T2.1 部分證據可並行起草

### T2.3：6-sample replication cohort 章節骨架 ⏳
- **預估**：W4
- **目標**：Strategy A §4「replication cohort」章節（HCC1954 / HCC1937 / H1437 / H2009 / COLO829 / DORADO）
- **依賴**：T2.1

---

## §3 Tier 3 (W4-W6) — Paper Draft + 工程平行 🟡

### T3.1：Paper Full Outline ⏳
- **預估**：W4-W5
- **格式**：Tool paper（Bioinformatics / NAR GB）
- **內容**：Abstract + 6 章 + 6 主圖

### T3.2：GitHub Repo 整理 + Public 化 ⏳
- **預估**：W5
- **目標**：Reproducible release（決策 #9）

### T3.3：Docker Image Build + Tutorial ⏳
- **預估**：W5
- **目標**：1-hour install + run

### T3.4：Benchmark Suite 公開化 ⏳
- **預估**：W5-W6
- **目標**：HG002 subset 或 HCC1395 SEQC2 公開部分

### T3.5：Discussion 章節 ⏳
- **預估**：W6
- **內容**：63% framework gap + cancer-only acceptance + Normal BAM future direction

---

## §4 Tier 4 (W6+) — Reactive 擴展 🟢

### T4.1：Phase 2A Normal BAM cross-sample ⏳
- **觸發條件**：G4 characterization 升級 45% → 70%
- **目標**：跨 7 樣本 Normal BAM + ISM 整合驗證

### T4.2：GC/mappability/repeat 新軸 pilot ⏳
- **觸發條件**：Reviewer 質疑 framework 63% gap
- **目標**：補強 framework coverage

### T4.3：PI Report 4-29 errata + V6 sign-off ⏳
- **依賴**：T1.2 完成後一併打包
- **內容**：5 條 errata（PI Report 4-29 對應條目）+ V6 sign-off email

### T4.4：HG002 non-cancer pilot ⏳
- **觸發條件**：Reviewer 強質疑 cancer-only limitation
- **目標**：non-cancer benchmark validation

---

## §5 時間軸（視覺化）

```
2026-05 W3 ─┬─ T1.1 主軸正名             ✅ 2026-05-16
            ├─ T1.3 scaffolding          ✅ 2026-05-16
            ├─ T1.2 V6 prod tag          🔴 W3 結束（Hard Gate × 2）
            │
2026-05 W4 ─┼─ T2.1 Z-AUTO KDE 4 樣本     ⏳ 升 ⭐4 必要
            ├─ T2.2 HCC1395 primary 章節 ⏳
            ├─ T2.3 6-sample replication ⏳
            │
2026-05 W5 ─┼─ T3.1 Paper outline         ⏳
2026-06 W5 ─┤
            ├─ T3.2 GitHub public         ⏳
            ├─ T3.3 Docker tutorial       ⏳
            ├─ T3.4 Benchmark suite       ⏳
            │
2026-06 W6 ─┼─ T3.5 Discussion 章節        ⏳
            │
2026-06+   ─┼─ T4.1 Phase 2A Normal BAM    ⏳ G4 上升才觸發
            ├─ T4.2 新軸 pilot              ⏳ reactive
            ├─ T4.3 PI errata package      ⏳ T1.2 完成同期
            └─ T4.4 HG002 pilot             ⏳ reactive
```

---

## §6 阻塞解除優先級

若資源有限，依以下順序：
1. **T1.2 V6 production tag**（解鎖 T2 + T4.3）
2. **T2.1 Z-AUTO KDE 4-sample**（解鎖 ⭐4 升級）
3. **T2.2/T2.3 章節骨架**（解鎖 T3.1 paper draft）
4. **T3.2-T3.4 工程平行**（可與 T3.1 重疊）
5. **T3.5 Discussion**（依 T3.1 大局已定）
6. **T4.x reactive**（依 reviewer feedback 觸發）

---

## §7 已完成里程碑（歷史參考）

- ✅ 2026-04-13：HCC1395 Phase 2 pilot 97.3% sig
- ✅ 2026-04-17：Phase 1A locked (+0.0112 paired_full)
- ✅ 2026-04-18：HPFineNGroups Part B canonical filter 升級
- ✅ 2026-04-11：PON 97.77% 跨樣本驗證
- ✅ 2026-05-15：V6 ISM ⭐3 PARTIAL POSITIVE（multi-agent fan-out）
- ✅ 2026-05-16：T1.1 + T1.3 完成
- ✅ 2026-05-17：governance v3 D2 分流 + /scientific-rigor skill + 8 cross-ref + Q5 erratum + 5 snapshot notice
- ✅ 2026-05-18：本次 5 snapshot 深度更新到 Tier 1-4 + Thread D paper + V6 production tag

---

## §8 Pre-registration 適用範圍（依 /scientific-rigor §7.1）

所有 Tier 2-4 新研究方向開跑前須在 `research/<topic>/00_INDEX.md` 強制 3 欄：
- H 預測 / 否證條件 / Decision threshold

例外：T1.x 已 DONE 不需追溯 / T4.x reactive 觸發時逐項補

---

## §9 相關

- Current Focus：[01_current_focus_snapshot.md](01_current_focus_snapshot.md)
- Active Hypotheses：[02_active_hypotheses.md](02_active_hypotheses.md)
- Blockers：[04_blockers_and_risks.md](04_blockers_and_risks.md)
- Breakthrough Strategy：[../01_project_overview/04_breakthrough_strategy.md](../01_project_overview/04_breakthrough_strategy.md)
- Plan：`~/.claude/plans/tender-pondering-blossom.md`
- Thread D paper plan：`InterSubMod/research/thread_d_paper/00_PLAN.md`
- V6 production plan：`InterSubMod/research/selfphasing_v6_production/00_PLAN.md`
