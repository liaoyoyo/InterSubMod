---
id: ism-kb-10-research-status-active-hypotheses
name: "Active Hypotheses"
description: "當前 active 假說：Subclonal reconstruction 主軸下 G-A 跨樣本重現 (V10, 資料 ready)、G-B subclone 甲基 somatic-specific vs germline-allelic 對照、四道 NEGATIVE methods (HD-1-independent)。⚠ 2 週有效；queue SoT = hypothesis_queue.json。"
status: active
last_verified: 2026-06-15
content_nature: runtime-fact
doc_type: reference
verified_scope: "active hypotheses against docs/CURRENT_FOCUS.md pinned 焦點區塊 (2026-06-12) subclonal axis"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-09-conclusions-hypothesis-queue-snapshot
tags: [status, hypotheses, active, subclonal-reconstruction, somatic-haplotag, methylation]
canonical_paths: [10_research_status/02_active_hypotheses.md]
alias_paths: []
---

> 🔄 **2026-06-15 re-sync 至當前主軸**（前版凍在 2026-05-18 Thread-D/V6 era）。**queue 真值 SoT = `research/autoresearch/hypothesis_queue.json`**；本檔為快速定向鏡像。

# Active Hypotheses

> 主軸：Subclonal reconstruction (somatic haplotag + methylation, ONT)。以下為當前 open 工作項（grounded 自 CURRENT_FOCUS pinned 區塊）。

## 當前 open 工作項（subclonal 主軸）

| ID/面 | 假說/工作 | 狀態 | 資料就緒 |
|------|----------|------|---------|
| **G-A** | V10 跨 5 樣本（乳腺3+肺2）甲基 subclone 重現 → 衝 ⭐4 | 🟢 可直接跑 | tagged BAM+VCF+ISM 齊；normal 甲基 5/6 |
| **G-B** | subclone 甲基是 somatic-specific 還是 germline-allelic（對照）| 🔴 gate（必做才能 claim subclone-甲基）| 待設計 |
| **四道 NEGATIVE methods** | 甲基→TP/FP filter 四道 NEGATIVE 寫成 methods（**HD-1 獨立、今天可寫、防彈**）| 🟢 可寫 | 已 concluded |
| **COLO829 normal 甲基** | 補第 6 樣本 normal 甲基（ONT_PAO 查 or re-basecall）| ⏳ | 缺 |

## 🔴 HD-1 gate（用戶決定，hold）

phasing by-construction 循環依賴 → 跑 **R-SELFREF（~25-50hr C++ 全 7-sample flag-on 負控）** 升 Grade A，**或**降為 characterization。NEGATIVE methods 與 HD-1 **獨立**（不等 HD-1 即可寫）。

## 誠實定論（甲基救 unphase / tag 矯正線 V1-V12，已整理）

真值 SoT = `InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md`。
- T1 unphase 救援 SUPPORTED+caveats（headline 0.885；僅 ~6% unphase 可嘗試）。
- T2 OVERSTATED（只證 1-1/2-1 可分，非歸 H3）。
- **T3 存在性窄翻案 + 可用性 NEGATIVE**（救 ambiguous read lean<0.5）。
- 甲基 = germline-haplotype 層級。memory `project_methyl_phasing_assist_line`。

## ❌ DEAD（勿註冊）

甲基化當 FP filter（⭐2 L4）；T3 subclone 可用性；H013-H018（filter/caller-F1，已降權/rejected — 全追溯 concluded-DEAD）。

## 相關

- live 主軸：[01_current_focus_snapshot.md](01_current_focus_snapshot.md) → [../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
- queue SoT：`research/autoresearch/hypothesis_queue.json`
- 阻塞：[04_blockers_and_risks.md](04_blockers_and_risks.md) · queue 快照：[../09_conclusions/05_hypothesis_queue_snapshot.md](../09_conclusions/05_hypothesis_queue_snapshot.md)
