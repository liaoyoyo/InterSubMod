---
id: ism-kb-10-research-status-blockers-and-risks
name: "Blockers & Risks"
description: "當前阻塞：HD-1 gate (R-SELFREF or 降 characterization)、G-B subclone 甲基歸因對照、COLO829 normal 甲基缺、6 normal BAM 單帳號 SPOF；對外誠實天花板風險。⚠ 2 週有效。"
status: active
last_verified: 2026-06-15
content_nature: runtime-fact
doc_type: reference
verified_scope: "blockers against docs/CURRENT_FOCUS.md pinned 焦點區塊 (2026-06-12) subclonal axis + 2026-06-14 external-validation"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-10-research-status-next-milestones
  - ism-kb-07-derived-features-coverage-multiple
  - ism-kb-05-data-formats-merged-dataset-pitfalls
tags: [status, blockers, risks, hard-gate, subclonal-reconstruction, hd-1-gate]
canonical_paths: [10_research_status/04_blockers_and_risks.md]
alias_paths: []
---

> 🔄 **2026-06-15 re-sync 至當前主軸**（前版凍在 2026-05-18 Thread-D/V6 era）。live SoT = `InterSubMod/docs/CURRENT_FOCUS.md`。

# Blockers & Risks

> 主軸：Subclonal reconstruction (somatic haplotag + methylation, ONT)。

## 🔴 阻塞（gate）

| # | 阻塞 | 性質 | 解法 |
|---|------|------|------|
| **HD-1** | phasing by-construction 循環依賴 | 用戶決定（hold）| 跑 R-SELFREF ~25-50hr C++（7-sample flag-on 負控）升 Grade A，**或**降為 characterization。⚠ NEGATIVE methods 與此**獨立**，可先寫。 |
| **G-B** | subclone 甲基 = somatic-specific 還是 germline-allelic 未分離 | 科學 gate | 設計 germline-het null + somatic-controlled 對照（必做才能 claim subclone-甲基）|
| **COLO829 normal 甲基** | 第 6 樣本 normal 甲基缺 | 資料 | 查 ONT_PAO or re-basecall（缺它仍可用 5 樣本跑 G-A）|

## 🟡 風險

- **6 normal BAM 單帳號 SPOF**：全部 zhenyu112 帳號持有 → 帳號失效則 normal 甲基不可得。契約：`InterSubMod/docs/data_specs/20260612_external_data_dependencies_01.md`。
- **單樣本 / 單 pipeline 天花板**：所有 ASM 單樣本 HCC1395 ⭐3；跨 6 樣本 6/6 是 phenomenon-level（excess-over-null）；private 0/38 underpowered → single-pipeline 封頂 ⭐3，6 樣本資產才推 ⭐4。
- **對外 over-claim 風險（必守口徑）**：① 甲基非 reconstruction 驅動 ② BRCA2 ≠ 乾淨 cis 錨點（subclone/copy-confounded；乾淨 cis 改 chr17/TBC1D16）③ ASM real but non-directional + non-discriminative（勿用「方向 POSITIVE」over-claim）④ 「Grade A」實為 Grade B+（R-SELFREF 對照未跑）。
- **ISM 創新點口徑**（投稿）：= 無監督 read×read 距離矩陣 PERMANOVA + normal-baseline cis-test + somatic-subclone 目標；🔴 禁用「對手二代定序/缺顯著性檢定」當差異（cvlr/ASMS/MethylBERT 都 ONT-capable + 有 randomization 檢定）。

## ❌ 已關閉（非阻塞，勿重開）

甲基→FP filter（⭐2 L4 DEAD）；neuter bug（早已修，燈#2 GREEN，勿當 live）；T3 可用性。

## 相關

- live 主軸：[01_current_focus_snapshot.md](01_current_focus_snapshot.md) → [../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
- 里程碑：[05_next_milestones.md](05_next_milestones.md) · 外部依賴契約：`InterSubMod/docs/data_specs/20260612_external_data_dependencies_01.md`
