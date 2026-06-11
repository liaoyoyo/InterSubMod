---
title: "Handoff 給論文 session — BRCA2「copy」→「subclone/copy」精確修改字串 + Martin-Trujillo category-error 修正"
date: 2026-06-07
status: in_progress
type: handoff (paper-session edit list; I do NOT edit validated files — governance per paper_readiness HD-3)
cycle_id: 20260607_capstone_verification_reconciliation
source: InterSubMod/docs/experiments/in_progress/2026/06/20260607_methylation_clustering_capstone_synthesis_01.md (ledger 94)
affects: docs/CURRENT_FOCUS.md, knowledge/11_external_literature/07_asm_cis_cancer_impact.md, knowledge/11_external_literature/10_paper_readiness_convergence.md, docs/reports/validated/2026/06/20260602_..._phasing_GradeA與ASM收斂/master_draft.md, docs/reports/in_progress/2026/06/20260608_G6_methods_negative_backbone_draft_01.md
---

# Handoff：BRCA2「copy」→「subclone/copy」+ % 不 robust + Martin-Trujillo 修正

> **背景**：論文 session 已做 06-08 amendment（BRCA2 copy-predominant / chr17 clean / dosage refuted）——**方向正確**。本 handoff 是 capstone（4-agent 獨立驗證 wf_700971bc-309）後的 **3 個精煉**，非推翻。**改 validated 文件是 governance 決定，由你/論文 session 核准套用，我不默默改。**

## 3 個精煉（為什麼要改）

| # | 精煉 | 證據 |
|---|------|------|
| **R1** | **「copy」→「subclone/copy」** | longphase-S `HaplotagStrategy.cpp:505-516`：HP1-1 = germline-H1 + 帶 somatic ALT = **somatic SUBCLONE tag，非 copy**（源碼 + 文獻 pubmed 25066126 subclone-specific methylation 三證；SOMATIC_H4 是 dead code）|
| **R2** | **「~80% copy / ~20% cis」精確 % → 質性** | 分解 `d_HP≈d_copy+d_within` **非閉合**（殘差 BRCA2 9% / chr17 40%）→ **% split 不 robust**。保留 d_within=−0.023(marginal) 當 focal 殘餘，但**勿寫「80%/20%」** |
| **R3** | **Martin-Trujillo (CN-dosage) 不 corroborate** | 我們 FAST **決定性 REFUTED copy-DOSAGE**（MW p=0.6183, signed ρ=−0.083）；BRCA2 confound 是 **subclone 非 CN-dosage** → 引 CN-dosage 文獻當 corroborate 是 **category error** |

## 精確修改字串（OLD verbatim → NEW）

### 檔 1 — `InterSubMod/docs/CURRENT_FOCUS.md` line 24

**OLD**：
`3. **BRCA2 ≠ 乾淨 cis 錨點**（06-07 重分析 = ~80% copy-artifact）；乾淨 cis 改 **chr17/TBC1D16**；`

**NEW**：
`3. **BRCA2 ≠ 乾淨 cis 錨點**（06-07 重分析 = **subclone/copy-confounded 主導**，HP1-1 是 somatic subclone tag 非 copy；focal cis d_within=−0.023 marginal·單樣本不可分離；**% split 不 robust**）；乾淨 cis 改 **chr17/TBC1D16**；`

### 檔 2 — `knowledge/11_external_literature/07_asm_cis_cancer_impact.md` line 23

**OLD**（節錄關鍵）：`列為 **predominantly copy-artifact**（HP-axis Δβ=−0.122 ≈ d_copy −0.11 + d_within −0.023，**~80% copy**）；... 誠實口徑 = 「**~80% copy + 邊際 ~20% 真 cis**」`

**NEW**：`列為 **predominantly subclone/copy-confounded**（HP1-1 = longphase-S 的 somatic SUBCLONE tag〔germline-H1 + somatic-ALT, HaplotagStrategy.cpp:505-516〕非 copy；HP-axis Δβ=−0.122 ≈ d_copy −0.11〔subclone vs germline 同 REF〕+ d_within −0.023〔focal allele within subclone, marginal perm p=0.022〕）；... 誠實口徑 = 「**subclone/copy 主導；focal 突變 cis 在單樣本與 subclone 背景不可乾淨分離（CAMDAC 同此限制）**」——**勿寫精確 %（split 不 robust，分解殘差達 40%）**，勿寫「BRCA2 真 cis-driven」亦勿寫「純 copy-artifact」`

### 檔 3 — `07_asm_cis_cancer_impact.md` line 26（Martin-Trujillo，**最重要修正**）

**OLD**：`> - **Martin-Trujillo (CN 解釋 82-92%) 現在反而 corroborate copy-artifact reclassification**（「吃掉 BRCA2 的 confond」），比 ASM×CN ρ=−0.055 pilot 更直接。`

**NEW**：`> - ⚠ **Martin-Trujillo (CN-DOSAGE 解釋 82-92%) 不 corroborate BRCA2 reclassification** —— 我們 06-07 FAST 分析**決定性 REFUTED copy-DOSAGE**（MW p=0.6183, signed ρ=−0.083 反向）；BRCA2 的 confound 是 **subclone（細胞群差異）非 CN-dosage**。引 CN-dosage 文獻當 corroborate 是 category error。正確對照 = **subclone-specific methylation**（pubmed 25066126/24356097, prostate/CLL/myeloma/breast）+ **CAMDAC**（單樣本 focal-cis vs subclone 不可分離, biorxiv 2020.11.03.366252）。`

### 檔 4 — `07_asm_cis_cancer_impact.md` line 9（frontmatter verified_scope）

**OLD 片段**：`(BRCA2=copy-artifact, chr17/TBC1D16 lone clean cis)`
**NEW**：`(BRCA2=subclone/copy-confounded〔HP1-1=subclone tag〕, % split not-robust, chr17/TBC1D16 lone clean cis; per capstone wf_700971bc-309 / ledger 94)`

### 檔 5 — `knowledge/11_external_literature/10_paper_readiness_convergence.md`

- **line 41** OLD `**BRCA2 不再是 cis 錨點 = copy-artifact（06-07）**` → NEW `**BRCA2 不再是 cis 錨點 = subclone/copy-confounded（06-07；HP1-1 是 subclone tag 非 copy，% 不 robust）**`
- **line 52** OLD `copy-partition 分解` → NEW `subclone/copy-partition 分解`；OLD `**BRCA2 + chr5 = copy-artifact**` → NEW `**BRCA2 + chr5 = subclone/copy-confounded**`
- **line 98 (HD-3)** OLD `BRCA2 = copy-artifact；只 chr17 clean` + `~80% copy + 邊際未確立 ~20% 真 cis` → NEW `BRCA2 = subclone/copy-confounded（HP1-1=subclone tag）；只 chr17 clean` + `subclone/copy 主導 + focal d_within=−0.023 marginal（**% split 不 robust，勿寫 80/20**）`；並把句尾 Martin-Trujillo 一句改為檔 3 NEW 口徑。
- **line 133 (FT-2)** OLD `cis-anchor vs copy-artifact（BRCA2 ~80% copy；...）` → NEW `cis-anchor vs subclone/copy-confound（BRCA2 subclone/copy 主導，HP1-1=subclone tag，% 不 robust；hypo≠canonical hyper...）`

### 檔 6 — `master_draft.md`（20260602 phasing 週報）lines 120 + 124

- **line 120 amendment** OLD（verbatim，含 markdown）`**predominantly copy-artifact**（~80% copy）` → NEW `**predominantly subclone/copy-confounded**（HP1-1=somatic subclone tag 非 copy；**% split 不 robust**）`；保留「copy-DOSAGE 決定性非 driver」一句（正確）。
- **line 124 (c)** 已標 SUPERSEDED — 補一句指向檔 1 NEW 口徑即可。

### 檔 7 — `docs/reports/in_progress/2026/06/20260608_G6_methods_negative_backbone_draft_01.md` line 50

**OLD**：`copy-**partition**（非 dosage）仍在最強訊號上 confound HP-axis magnitude` ... `**BRCA2 (chr13) + chr5 = copy-artifact**（BRCA2 HP-axis −0.122 ≈ d_copy −0.11 + 邊際 d_within −0.023, perm p=0.022, n=197）`

**NEW**：`**subclone/copy-partition**（非 dosage；HP1-1 是 somatic subclone tag 非 copy, HaplotagStrategy.cpp:505-516）仍在最強訊號上 confound HP-axis magnitude` ... `**BRCA2 (chr13) + chr5 = subclone/copy-confounded**（BRCA2 HP-axis −0.122 = subclone vs germline 差 d_copy −0.11 + 邊際 focal d_within −0.023, perm p=0.022, n=197；**% split 不 robust**）`

## 一句話對外口徑（所有文件統一）

> BRCA2 的 HP-axis 訊號主要是 **somatic-subclone vs germline 的甲基差**（HP1-1 是 longphase-S 的 subclone tag，非 copy；非 CN-dosage——dosage 已 REFUTED）；**焦點突變的 cis 在單樣本與 subclone 不可乾淨分離（CAMDAC 同此限制）**。**chr17/TBC1D16 是唯一乾淨 cis exemplar。勿寫精確 % split。**

## Provenance

capstone `InterSubMod/docs/experiments/in_progress/2026/06/20260607_methylation_clustering_capstone_synthesis_01.md` · 驗證 workflow `wf_700971bc-309` · ledger 80→94 · scripts 37/62/63/64/**65** · `dwithin_validity.json`（可重跑）。
