---
id: ism-kb-10-research-status-current-focus-snapshot
name: "CURRENT_FOCUS Snapshot"
description: "docs/CURRENT_FOCUS.md 結構化鏡像；主軸 = Subclonal reconstruction (somatic haplotag + methylation, ONT)，2026-06-11 取代 G6。⚠ 2 週有效；live SoT 永遠以 CURRENT_FOCUS.md 為準。"
status: active
last_verified: 2026-06-15
content_nature: runtime-fact
doc_type: reference
verified_scope: "mirror of docs/CURRENT_FOCUS.md pinned 焦點區塊 (2026-06-12 更新) + 2026-06-14 external-validation section"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-active-hypotheses
  - ism-kb-10-research-status-blockers-and-risks
  - ism-kb-10-research-status-next-milestones
  - ism-kb-01-project-overview-current-phase
  - ism-kb-03-pipelines-f1-baseline-canonical
  - ism-kb-09-conclusions-hypothesis-queue-snapshot
tags: [status, snapshot, current-focus, subclonal-reconstruction, somatic-haplotag, methylation]
canonical_paths: [10_research_status/01_current_focus_snapshot.md]
alias_paths: []
---

> 🔄 **2026-06-15 re-sync 至當前主軸**（前版凍在 2026-05-18 Thread-D era，已 2 代過時）。**live SoT 永遠 = `InterSubMod/docs/CURRENT_FOCUS.md`**；本檔為快速定向鏡像。

# CURRENT_FOCUS Snapshot

> 一句結論：主軸 = **Subclonal reconstruction using somatic haplotagging and methylation profiles (Nanopore)**；6 樣本×3 癌種資產齊 + normal 甲基 5/6 ready，⭐3→⭐4 候選；卡關在 **HD-1 gate**（用戶決定）。

- 適用：決策前快速了解現況、判斷下一步與阻塞。
- 深度需求 → 讀 `InterSubMod/docs/CURRENT_FOCUS.md`（pinned 焦點區塊 = 最新）。

## 當前主軸（2026-06-11 pivot，用戶確認）

**Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing** — 取代 G6 phasing；G6 LOH-phasing / G1 ZAR1L/BRCA2 ASM **降為支撐材料**。

**兩個互補 SoT 面**：
1. 甲基-phasing-assist foundation：`InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`（V1-V12 + 6 樣本資產 + G-A~E gap）。
2. ASM-characterization + 四道 NEGATIVE + LOH-phasing 脊柱 + HD-1 gate：`InterSubMod/docs/concepts/2026/06/20260611_Subclonal_Reconstruction_Paper_Focus_整合篇章_01.md`。

| 主軸 | tier | 狀態 |
|------|------|------|
| 🆕 Subclonal reconstruction（新論文主軸）| ⭐3→⭐4 候選 | 兩面 foundation 已立；6 樣本資產 + normal 甲基 5/6 ready；🔴 gate=HD-1 + G-B |
| ~~G6 LOH-constrained phasing~~（降支撐）| ⭐3 | park（=新主軸 phasing 脊柱）|
| ~~G1 ZAR1L/BRCA2 ASM~~（降支撐）| ⭐3 | park（=新主軸 ASM characterization 層）|

## 資產（緩解單樣本限制）

6 cell line × 3 癌種（HCC1395/1937/1954 乳腺 · H1437/H2009 肺 · COLO829 黑色素瘤）皆有 somatic-haplotag BAM + somatic VCF + ISM TP/FP；**matched-normal 甲基 5/6 有**（只 COLO829 缺）→ V10 跨 5 樣本可直接跑衝 ⭐4。契約 `InterSubMod/docs/data_specs/20260612_external_data_dependencies_01.md`（6 normal 全 zhenyu112 帳號 = SPOF）。

## 誠實天花板（對外必守）

- 甲基 = germline-haplotype 層級；非 reconstruction 驅動、非 variant filter（DEAD）。
- T3 存在性窄翻、可用性 NEGATIVE；T2 只證 1-1/2-1 可分（非歸 H3）。
- 外部文獻驗證庫 59 源、0 真 CONFLICT、稽核 CLEAN（`project_external_validation_library`）。

## ❌ DEAD（勿再開）

甲基化當 FP filter（⭐2 L4）；T3 subclone「可用性」（救 ambiguous read 偏向反了）；pure-methylation TP/FP 判別 / TO germline-FP filter / ASM-as-discriminator / LR methyl-augmented filter（全 concluded dead guardrail）。

## 相關

- live 權威：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
- 活躍假說：[02_active_hypotheses.md](02_active_hypotheses.md) · 阻塞：[04_blockers_and_risks.md](04_blockers_and_risks.md) · 里程碑：[05_next_milestones.md](05_next_milestones.md)
- Hypothesis queue：[../09_conclusions/05_hypothesis_queue_snapshot.md](../09_conclusions/05_hypothesis_queue_snapshot.md)
- 機械態 cycle：`/cycle-state`（注意：多週主軸權威 SoT = CURRENT_FOCUS，非 active.json — 見 CLAUDE.md §6 雙層 SoT）
