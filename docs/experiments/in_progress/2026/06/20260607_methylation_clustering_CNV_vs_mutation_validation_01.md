---
title: "甲基分群驅動來源完整驗證 — CNV/copy-partition vs somatic-mutation-cis"
date: 2026-06-07
status: in_progress
tier: "⭐3 (single-sample HCC1395, multi-channel confound-controlled characterization)"
task_type: B_validation
cycle_id: 20260607_fast_cnv_validation_results
parent_cycles: 20260602_asm_cn_confound_pilot, 20260603_copy_partition_3test_quantified
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/fast_cnv_validation.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/alignment_channel_test.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/copy_partition_confirm.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/cis_scan_full.json
scripts: research/tsg_promoter_asm_reviewer/scripts/62_fast_cnv_validation_batch.py,research/tsg_promoter_asm_reviewer/scripts/63_alignment_channel_test.py,research/tsg_promoter_asm_reviewer/scripts/39_cgi_annotation.py
sample: HCC1395 paired_full (single sample)
partial_flag: "subset — read-level copy-partition adjudication capped at 7/60 leaky-tag loci; alignment test on 24 loci; single sample"
---

# 甲基分群驅動來源完整驗證

> ## ⚠ 修正 banner（2026-06-07，longphase-S 源碼審計後）
> 本報告「**copy-PARTITION**」通道應讀作「**subclone/copy partition**」。源碼審計（`HaplotagStrategy.cpp:505-516`）確認 **HP1-1 = germline 單倍型1 + somatic ALT = somatic SUBCLONE tag，非 copy；非 CN-dosage**（copy-DOSAGE 已決定性 REFUTED, MW p=0.6183）。d_within 經三重證據驗為**有效 subclone-control**（HP1-1/ref 19/19 帶其他 somatic ALT），故 **質性成立**：BRCA2 = **subclone/copy 主導**。⚠ **2026-06-07 capstone（ledger 94）supersede**：原寫「~80% copy →『~81% subclone/copy』量化成立、僅 relabel」的**精確 % 不 robust**（加性分解跨位點不閉合：BRCA2 殘差 9%、chr17 達 40%）→ **改質性，勿引用精確 %（81%/64%）**。raw 焦點對比 `d_focal_CLEAN` 10/11 可算（最強位點 untestable 是焦點突變在 subclone 內 clonal，非覆蓋偏差）。完整：`InterSubMod/docs/experiments/in_progress/2026/06/20260607_HP11_tag_definition_correction_record_01.md`（ledger 93→94）。

## 問題

HCC1395 的「甲基分群 / Δβ」是 **CNV / copy-partition / ploidy 驅動**，還是 **somatic-mutation-cis 驅動**？跟 HP-tag + 突變一致，會不會只是 **CNV 或 alignment 的副作用**？

## 方法：四個 confound 通道逐一測

| 通道 | 假設（若為真）| 測量 | 來源 cycle |
|------|------|------|-----------|
| **copy-DOSAGE** | \|Δβ\| 隨 CN 數量上升 | partial ρ(\|Δβ\|,CN) + CN-state class-contrast | 6-02 + FAST |
| **subclone/copy-PARTITION**（原「copy-partition」）| HP1 vs HP1-1 = somatic-subclone vs germline 甲基差（HP1-1=subclone tag）| d_copy vs d_within（within-tag） | 6-03 + 6-07 audit |
| **alignment** | 擴增/重複區 mismapping 造假分群 | soft-clip / secondary / NM-per-kb / HP1-1-vs-HP1 對齊 delta | FAST ③ |
| **mechanical** | 突變造/毀 CpG → 機械差異 | ref-sequence CpG gain/loss | 6-03 |

## 結果（每通道）

### ① copy-DOSAGE — 決定性 REFUTED

CN-state class-contrast（HP-axis sig \|Δβ\|，`master_o1_cn.tsv`）：

| 地層 | n | median \|Δβ\| |
|------|--:|--:|
| neutral nonLOH (CN=2 平衡) | 76 | 0.0811 |
| gain nonLOH | 2542 | 0.0737 |
| cnLOH (CN=2+LOH) | 913 | 0.0821 |
| gain LOH | 1580 | 0.0701 |
| loss | 31 | 0.0824 |

- 跨 CN-state 全部 **0.070–0.082**；neutral vs gain **MW p=0.62 無法區分**。
- Falsifier：signed ρ(\|Δβ\|,CN) = **−0.083**（p=2.6e-9）；dosage artifact 應為**正** ρ → **反向 = REFUTED**。
- cnLOH (0.082) ≈ 平衡 neutral (0.081) → copy 結構改變不改分群量級。
- **結論：甲基分群量級不是 CN 數量驅動。** 對齊 6-02 H-A (ρ=−0.055)，本輪是更強的直接 class-contrast。

### ② copy-PARTITION — SUPPORTED（主導量級）

6-03 的 within-tag 拆解（7 個可算位點）：d_copy（HP1/ref vs HP1-1/ref，純 copy）撐起 d_HP 量級**中位 ~64%**（BRCA2 90%）。即 HP-axis 報出的「分群量級」**多數來自 copy-partition**，非 somatic allele。⚠ 只在 leaky-tag（somatic tag 有 ALT+REF）位點可算 → **7/60**，最強 6 個純 ALT tag 測不了。

### ③ alignment — REFUTED（5/6；chr11 旗標）

24 位點 BAM 對齊測試（6 untestable + 3 control + 15 neutral baseline），closes 6-02 MAPQ 盲點：

| group | softclip | NM/kb | sec_rate | supp_rate | HP1-1−HP1 softclip |
|-------|--:|--:|--:|--:|--:|
| untestable (6) | 0.027 | 39.3 | 0.0009 | 0.028 | **−0.004** |
| control (3) | 0.021 | 28.2 | 0.0 | 0.042 | +0.030 |
| neutral baseline (15) | 0.025 | 30.7 | 0.39* | 0.053 | −0.018 |

- untestable 的 **soft-clip 不高、secondary≈0、supplementary 正常、HP1-1 不較差對齊（delta −0.004）** → **分群不是 alignment-driven**。
- 唯一紅旗：**chr11:64557316 NM/kb=87.8**（2-3× 偏高）→ divergent 區，demote。其餘 5 個（含最強 chr4 d_cis=0.706）對齊正常。
- *neutral sec_rate 高是少數 baseline 位點落在 repetitive 區（chr21:9924901 sec=3.6）；**最強訊號反而不在 multi-mapping 區**。

### ④ mechanical — 81/816 排除

cis_scan：mechanical_cis CREATES_CpG 52 + DESTROYS_CpG 29 = 81 位點帶機械假象（如 chr19:14434617 DESTROYS_CpG），permutation 排不掉，標 artifact。

### CGI/gene context（調控 vs 結構）

- chr17:79991120 (**TBC1D16**)：差異 CpG **富集在 CGI**（in 0.365 vs out 0.129，MW **p=1.8e-7**）→ regulatory-cis 像。
- BRCA2 (ZAR1L)：差異**避開 CGI**（in 0.017 vs out 0.232，p=9.6e-11）→ 結構/copy。
- 最強 6 個：全 **intergenic + CGI 沙漠**（離 CGI 35kb–716kb）→ 結構 prior。

## 全 816 位點決策表（見樹也見林）

| 類別 | 數量 |
|------|--:|
| not-cis (T0/T2) | 682 |
| mechanical-artifact (CpG 造/毀) | 81 |
| T3-nominal (過不了 Bonferroni 0.05/816) | 44 |
| **Bonferroni-clean T3** | **10** |

**10 個 Bonferroni-clean T3 逐位點**：
- **chr17:79991120 (TBC1D16)** — ✅ 唯一乾淨真 cis（allele 主導 + CGI 富集 + NEUTRAL + alignment 乾淨）
- chr13:32315128 (BRCA2/ZAR1L)、chr5:6201328 — subclone/copy-confounded（HP1-1=somatic subclone tag 非 copy、非 CN-dosage；焦點 cis d_within 邊際，% 不 robust）
- chr19:14434617 — mechanical-artifact (DESTROYS_CpG)
- chr4:133868344 (d_cis=0.706)、chr7:79941963、chr20:59439285/61415690/61564264、chr11:64557316 — **6 個 untestable**（pure-ALT tag，CGI 沙漠，非 dosage/alignment 假象，但 copy-vs-cis 單樣本測不了；chr11 額外帶高 NM/kb）

## 結論：用戶假設三層判定

| 子主張 | 判定 |
|------|------|
| 甲基分群 = copy-**DOSAGE**（CNV 數量）副作用 | **❌ REFUTED**（量級跨 CN-state 平、signed ρ 反向）|
| 甲基分群 = **alignment** 副作用 | **❌ REFUTED**（5/6 對齊正常；chr11 旗標）|
| 分群**量級** = **subclone/copy-PARTITION**（HP1-1=subclone tag，非 copy、非 CN-dosage）| **✅ SUPPORTED 主導**（質性；**精確 % 不 robust**，分解跨位點不閉合，勿引用 ~64%/~81%）|
| 完全 artifact、無真訊號 | **❌ REFUTED**（within-tag permutation 4/4 顯著；chr17 乾淨）|

**綜合**：甲基分群**存在性與量級的主因是 subclone/copy partition（somatic-subclone vs germline 的結構性甲基差；HP1-1 是 subclone tag，非 copy），不是 CN 數量、不是 alignment、不是純 artifact。** 跟 HP-tag 一致是因 HP1-1 = somatic subclone；跟突變一致是因突變定義那個 subclone。**焦點突變的乾淨 cis 在單樣本不可分離；全 816 只有 1 個（chr17/TBC1D16，d_within=0.142 allele 主導）是 subclone-control 後存活的乾淨真 cis。**

## 🔴 硬天花板（誠實）

單樣本 HCC1395 內，**對最強訊號集 copy 與 mutation 近乎不可識別**（純 ALT tag 測不了 = 7/60 才可算；CN=2 乾淨地層只 76 位點且 75/76 是 TP）。要裁決那 6 個最強訊號**必須換 design**（COLO829 跨樣本 / 實驗擾動）——但 COLO829 無 SEQC2-style CN truth，CN-disentanglement 那支綁死 HCC1395。

## 限制

- 單樣本；read-level copy-partition 只 7/60 可算；alignment 測 24 位點；6 untestable 的 copy-vs-cis 未定（CGI-沙漠結構 prior 非證實）。
- 「64%」與「chr17 唯一乾淨」是 7-leaky / 可測子集數字，**引用須 scope-lock**。
- 業界對照：CN→甲基通道真實（dosage compensation，CAMDAC）但 CGI/promoter-localized + directional，與我們 CGI-沙漠最強訊號不符 → desert 訊號非 copy-dosage 簽名；field standard 去 confound = CAMDAC（ASCN+purity），HP-axis 是其 read-level 類比。

## 後續

- chr17:79991120 (TBC1D16) 個案 deep-dive + COLO829 複現（升 ⭐4 必要）。
- 6 untestable：需 per-read copy-phasing 或 COLO829 才能裁決。
- chr11:64557316 額外排除 divergent-region 假象。
