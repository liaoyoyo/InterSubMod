---
title: "修正紀錄 — HP1-1 tag 是 somatic subclone tag（非 copy）；copy-partition 改 subclone/copy"
date: 2026-06-07
status: in_progress
tier: "correction record (source audit + clean re-derivation + d_within validity)"
cycle_id: 20260607_hp11_tag_definition_audit
corrects: 20260603_copy_partition_3test_quantified, 20260607_fast_cnv_validation_results
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/subclone_tag_rederivation.json
scripts: research/tsg_promoter_asm_reviewer/scripts/64_subclone_tag_rederivation.py
affects: 20260603_ISM_copy_partition_confound_verification_01, 20260607_methylation_clustering_CNV_vs_mutation_validation_01, master_draft(paper §2.2 BRCA2), CURRENT_FOCUS(BRCA2 ~80% copy-artifact)
---

# 修正紀錄：HP1-1 tag 定義審計

## 觸發

用戶嚴格質疑：copy-PARTITION 判定靠 `d_copy = HP1/ref vs HP1-1/ref`，但 **HP1-1 是 longphase-S 依 caller 的 sSNV ALT 指派的**。那「HP1-1/ref」是不是 tagging artifact，整個 copy-partition 判定是否有問題？

## 查證（source + 實證）

### 1. 源碼確認 — HP1-1 = somatic subclone tag

`longphase-s/src/haplotag/HaplotagStrategy.cpp:505-516`：read 被判 **H1_1 的條件 = (a) germline het SNP 支持 GERMLINE_H1 AND (b) `tumorMaxHPcount≠0`（帶 somatic 變異支持 SOMATIC_H3）**。即：

> **HP1-1 = germline 單倍型1 + 帶 somatic ALT = 一個 somatic SUBCLONE tag，不是「一個 copy」。**

BAM 編碼 `HP:Z:1` / `HP:Z:1-1` / `HP:Z:2` / `HP:Z:2-1`（HaplotagType.h:334 enum→字串）。**用戶的觀察正確。**

### 2. d_within 有效性 — 三重證據 = VALID（subclone-control）

驗 BRCA2 的 HP1-1/ref reads（19 條）到底是真 subclone reads 還是 leak/error：

| 檢查 | 結果 |
|------|------|
| focal base | 全 **G（REF）**，median BQ=26.5（confident REF，非 ALT mis-call）|
| 帶其他 somatic ALT？ | **19/19** 都帶（在 chr13:32317522/32324831/32339132，平均 1.53 個）|
| 甲基歸哪群 | **subclone**（d_copy=−0.110 大遠離 HP1/ref；d_within=−0.023 小近 HP1-1/alt）|

→ **HP1-1/ref = 真 somatic-subclone reads，只在焦點 REF。所以 `d_within`（HP1-1 內 alt vs ref）是有效的 subclone-controlled 焦點 allele 估計；分解 `d_HP ≈ d_copy + d_within` 成立。**

### 3. 乾淨重算（繞過 subclone tag）— `d_focal_CLEAN`

用 germline 相位（HP∈{1,1-1}）+ read 實際焦點鹼基（ref/alt）算焦點對比，不靠「-1」tag：

| 位點 | d_HP | d_focal_CLEAN | d_within | subclone/copy frac |
|------|-----:|--------------:|---------:|-------------------:|
| chr17 (TBC1D16) | 0.122 | 0.178 | **0.142** | **−0.16（allele 主導=乾淨 cis）** |
| BRCA2/ZAR1L | −0.122 | −0.099 | −0.023 | **0.81（subclone/copy 主導）** |
| chr5 | 0.129 | 0.115 | 0.059 | 0.54 |
| chr4（之前測不了）| 0.549 | **0.615** | None | —（focal clonal，無 subclone control）|
| chr20:59439285 | 0.307 | 0.272 | None | — |

- **d_focal_CLEAN 10/11 可算**；d_within（subclone control）只 4/11。
- 「最強 6 個測不了」= **焦點突變在其 subclone 內 clonal**（無 within-subclone REF）→ 無 subclone control，**不是覆蓋偏差 artifact**；raw 對比仍可算且大。

## 結論：什麼要修 / 什麼成立

### 要修（label + 機制，非撤回數字）

| 原 | 改 |
|----|----|
| 「copy-partition」 | **「subclone/copy partition」**（HP1-1 是 subclone tag；d_copy = somatic-subclone vs germline 甲基差，已證非 CN-dosage → copy-identity 或 subclone 表觀）|
| BRCA2「~80% copy-artifact」 | 質性：**「subclone/copy 主導，非焦點突變 cis」**（⚠ **精確 % 不 robust**，見下方 capstone reconcile；勿引用「81%」）|
| 「7/60 測不了 = 覆蓋偏差」 | 「subclone CONTROL 限 leaky-tag；**raw 焦點對比 d_focal_CLEAN 到處可算**；untestable = 焦點突變 clonal」|

### 成立（核心未倒）

1. **「是 CNV/dosage 嗎」→ NO**（FAST CN-state contrast 不依賴 tag）。
2. **全域 NEGATIVE**（甲基不能 filter TP/FP）。
3. **chr17 > BRCA2**：chr17 d_within=0.142（allele 主導，cis 乾淨）vs BRCA2 −0.023（subclone 主導）。**強化**。
4. **質性結論 robust**：d_within 質性有效（subclone-controlled，source + `scripts/65` 19/19 重現確認）→ BRCA2 subclone-dominated（d_within 小 −0.023）、chr17 allele-dominated（d_within 大 0.142）。⚠ **但精確 % 量化（64% / 81%）NOT robust**（分解非閉合，chr17 殘差達 40%）——2026-06-07 capstone 對抗驗證定案，**勿引用精確 %，用質性框架**（reconcile ledger 94，supersede ledger 93「全撤回」與「保留 81%」雙向過度）。
5. **誠實核心**：somatic-ALT reads 本質就是 somatic subclone reads → 焦點-allele 甲基對比本質被 subclone/copy confound；單樣本無法分離焦點突變 cis vs subclone copy/表觀狀態。

## 對論文（master_draft §2.2 + CURRENT_FOCUS）

- **chr17/TBC1D16 為最乾淨 cis 候選 — 不變且強化。**
- BRCA2「~80% copy-artifact」→ 改「**subclone/copy 主導（HP1-1 是 somatic subclone tag，非 copy）；焦點突變 cis 在單樣本不可乾淨分離**」——**不寫精確 %**（% 量化不 robust，capstone 定案）。
- 「copy-partition」用語全面改「subclone/copy partition」。

> evidence_ledger `20260607_hp11_tag_definition_audit`（entry 93）。報告 `20260603_*` + `20260607_methylation_clustering_*` 已加修正 banner 指向本紀錄。
