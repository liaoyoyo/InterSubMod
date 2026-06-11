---
title: "Capstone — 甲基分群驅動來源：完整整理 + 獨立驗證 + 背景佐證"
date: 2026-06-07
status: in_progress
tier: "⭐3 (single-sample HCC1395; multi-channel confound-controlled characterization; independently re-verified)"
task_type: B_validation
cycle_id: 20260607_capstone_verification_reconciliation
parent_cycles: 20260603_copy_partition_3test_quantified, 20260607_fast_cnv_validation_results, 20260607_hp11_tag_definition_audit
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/copy_partition_confirm.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/fast_cnv_validation.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/alignment_channel_test.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/subclone_tag_rederivation.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/dwithin_validity.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/cis_scan_full.json
scripts: research/tsg_promoter_asm_reviewer/scripts/37_copy_partition_confirm.py,research/tsg_promoter_asm_reviewer/scripts/62_fast_cnv_validation_batch.py,research/tsg_promoter_asm_reviewer/scripts/63_alignment_channel_test.py,research/tsg_promoter_asm_reviewer/scripts/64_subclone_tag_rederivation.py,research/tsg_promoter_asm_reviewer/scripts/65_dwithin_validity_reproducible.py
verification: workflow wf_700971bc-309 (4-angle independent verification — numbers/source/literature/adversary)
sample: HCC1395 paired_full (single sample)
partial_flag: "subset — read-level subclone-control capped at leaky-tag loci; alignment test 24 loci; single sample; % quantification NOT robust"
---

# Capstone：甲基分群驅動來源 — 完整整理 + 獨立驗證

> 本文整合 5 個 cycle（copy-partition 發現 → CNV 三層驗證 → HP1-1 subclone 審計 → capstone reconcile），所有數字經 **4-agent 獨立驗證 workflow（wf_700971bc-309）從 raw cache 重算確認**。

## 0. 一句話結論（reconciled）

甲基分群**不是 CN 數量(dosage)、不是 alignment、不是純 artifact**；其主因是 **somatic-subclone vs germline 的甲基差異**（HP1-1 是 longphase-S 的 **somatic subclone tag**，非 copy）。**焦點突變的 cis 效應在單樣本與 subclone 背景不可乾淨分離**（文獻 CAMDAC 同此限制）。質性上：**BRCA2 = subclone 主導**（焦點 cis 小）、**chr17/TBC1D16 = allele 主導**（唯一乾淨 cis 候選）。⚠ **精確 % 量化不 robust，勿引用。**

## 1. 調查全弧（時間線）

| 日期 | cycle | 結論 | ledger |
|------|-------|------|--------|
| 06-03 | copy_partition 三項量化 | HP-axis 量級多由「copy」撐起；chr17 vs BRCA2 差異 | 80→83 |
| 06-07 | CNV 三層驗證（FAST + alignment）| dosage REFUTED / alignment REFUTED / 決策表 816→10 | 84,86,87,88 |
| 06-07 | HP1-1 tag 審計（用戶質疑）| HP1-1 = somatic subclone tag，「copy」改「subclone/copy」 | 93 |
| 06-07 | **capstone 驗證 + reconcile** | 全數字 reproducible；% 不 robust → 質性框架 | **94** |

## 2. 完整數據表（每數字經獨立重算 CONFIRMED）

### 2a. 四通道 confound 判定（CNV-vs-mutation 核心問題）

| 通道 | 判定 | 關鍵數字（獨立重算 ✓）|
|------|------|------|
| copy-**DOSAGE**（CN 數量）| **❌ REFUTED** | neutral-nonLOH \|Δβ\|=0.0811(n76) vs gain 0.0737(n2542) **MW p=0.6183**；signed ρ=**−0.0829**(p=2.6e-9) 反向 |
| **alignment** | **❌ REFUTED（5/6）** | untestable softclip 0.027 / sec≈0 / HP1-1−HP1 softclip **−0.004**；唯 chr11 NM/kb=87.8 |
| **subclone/copy-PARTITION** | **✅ 主導量級** | d_copy（subclone vs germline 同 REF）撐起 d_HP；質性主導 |
| **mechanical**（CpG 造/毀）| 81/816 排除 | CREATES 52 + DESTROYS 29 |

### 2b. 焦點位點四量測（bit-for-bit 重算確認）

| 位點 | d_HP | d_copy | d_within | d_focal_CLEAN | 質性判定 |
|------|-----:|-------:|---------:|--------------:|------|
| **chr17 (TBC1D16)** | 0.122 | 0.029 | **0.142** | 0.178 | **allele 主導 = 乾淨 cis** |
| **BRCA2/ZAR1L** | −0.122 | −0.110 | **−0.023** | −0.099 | **subclone 主導** |
| chr5 | 0.129 | 0.055 | 0.059 | 0.115 | 混合 |

> ⚠ `d_HP ≈ d_copy + d_within` **僅近似**（殘差 BRCA2 9% / chr17 40%；chr17 d_within>d_HP）→ **精確 % 不 robust**。質性方向（chr17 allele / BRCA2 subclone）robust。

### 2c. 全 816 決策表

not-cis(T0/T2) 682 · mechanical-artifact 81 · T3-nominal(過不了 Bonferroni) 44 · **Bonferroni-clean T3 10**（1 乾淨 cis chr17 / 2 copy-art BRCA2,chr5 / 1 mechanical chr19 / 6 untestable pure-ALT CGI-沙漠）。

## 3. HP1-1 機制（longphase-S 源碼，獨立複核 CONFIRMED）

`HaplotagStrategy.cpp:505-516`：read 判 **H1_1 ⇔ maxNormalHP==GERMLINE_H1 AND maxTumorHP==SOMATIC_H3（tumorMaxHPcount≠0）**。

- **HP1-1 = germline 單倍型1 + 帶 ≥1 somatic ALT = somatic SUBCLONE tag，非 copy。**
- **NUANCE（源碼新發現）**：`SOMATIC_H4`(REF phase) 是 **dead code**（hpCount[4] 從不 ++；hpCount[3]++ 僅在 `base==TumorAltBase`）→ tag 純粹「帶 somatic ALT」，無「ALT vs REF phase」之分。
- 另有 gate：norHPsimilarity≥0.6 + tumHPsimilarity≥0.6 + 單 PS block（否則降 unTag/H3）。
- **HP1-1/ref reads 必帶其他 somatic ALT**（code-mandated：hpCount[3]>0）。
- 經驗證實：ALT reads **只**出現在 HP1-1 tag（HP1/HP2 alt=0）→ HP1-1 確是 somatic-ALT subclone tag。

## 4. d_within 有效性（可重跑 `scripts/65`）

BRCA2 HP1-1/ref reads 驗證（`dwithin_validity.json`）：

| 檢查 | 結果 |
|------|------|
| n | 19 |
| focal base | 全 **G(REF)**，median BQ=26.0 |
| 帶其他 somatic ALT | **19/19**（32317522/32324831/32339132，平均 1.53）|

→ **d_within 質性有效**（HP1-1/ref 是真 subclone reads；subclone-controlled 焦點 allele 估計）。⚠ **但精確 % 不 robust**（分解非閉合）。

## 5. 背景文獻佐證（敘述合理性，獨立查證 CONFIRMED）

| 主張 | 文獻 | 判定 |
|------|------|------|
| longphase 把 somatic ALT reads 依 germline 相位分到 HP1-1 | longphase-S preprint (biorxiv 2025.11.20.689492) | ✅ |
| 同 germline haplotype 的腫瘤 subclone 有不同甲基化 | subclone-specific methylation (pubmed 25066126 / 24356097；prostate/CLL/myeloma/breast) | ✅ 真實生物 |
| 單一 somatic 點突變的 cis 甲基效應小且難在 bulk 分離 | cis-mQTL/ASM modest (PMC2893533/PMC7679001) | ✅ |
| somatic-ALT reads = subclone reads → 單樣本 focal-cis vs subclone 不可分離 | CAMDAC (biorxiv 2020.11.03.366252) + bulk ASE confounder；解法 = matched-normal diffASE | ✅ |

→ **修正後敘述與領域常識一致，無文獻反對。**

## 6. 對抗驗證後的最終敘述（什麼成立 / 什麼修正 / reviewer 還會挑什麼）

### 成立（reproducible + literature-backed）
1. dosage NO（MW p=0.62 + signed ρ 反向）、alignment NO（5/6）。
2. 全域 NEGATIVE（甲基不能 filter TP/FP）。
3. 質性：**chr17 allele 主導（乾淨 cis）> BRCA2 subclone 主導**——d_within 質性有效。
4. HP1-1 = somatic subclone tag（源碼 + 文獻 + 經驗三證）。

### 已修正（capstone 對抗抓出）
- 「copy-partition」→「**subclone/copy partition**」（label）。
- **% 量化（81%/64%）NOT robust** → 改質性框架（分解非閉合殘差 40%）。
- d_within 有效性證據從 hard-coded 字串 → **可重跑 `scripts/65`**。

### reviewer 仍會挑（誠實未解）
- **subclone vs copy（同細胞不同 amplicon copy）本身單樣本仍不可分**——「subclone/copy」合併標誠實但未解。
- **chr17 的 cis 宣稱仍受 subclone confound**（d_within 大可能是 chr17 subclone 甲基差大，非真 focal cis）——需 COLO829 / matched-normal diffASE。
- 單樣本天花板：最強 6 個 untestable 位點 copy-vs-cis 未定。

## 7. 對論文（CURRENT_FOCUS + master_draft §2.2）

- **chr17/TBC1D16 為最乾淨 cis 候選 — 不變且強化（d_within 質性有效）。**
- BRCA2「~80% copy-artifact」→「**subclone/copy 主導（HP1-1 是 somatic subclone tag，非 copy）；焦點突變 cis 在單樣本不可分離**」——**不寫精確 %**。
- 「copy-partition」用語全面改「subclone/copy partition」。
- （由 paper session 落地；本 capstone 的 `affects:` 標記）

## 8. Artifacts inventory

- scripts：37(copy-partition) / 62(FAST) / 63(alignment) / 64(rederivation) / **65(d_within validity, reproducible)**
- JSON：copy_partition_confirm / fast_cnv_validation / alignment_channel_test / subclone_tag_rederivation / **dwithin_validity** / cis_scan_full
- 報告：20260603_copy_partition / 20260607_CNV_validation / 20260607_correction_record（皆加修正 banner）
- ledger：80→94（append-only）；驗證 workflow wf_700971bc-309
