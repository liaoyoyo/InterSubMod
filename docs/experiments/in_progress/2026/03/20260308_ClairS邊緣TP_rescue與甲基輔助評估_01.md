<!--
建立時間: 2026-03-08 12:00
狀態: in_progress
目標: 驗證 ClairS 邊緣 TP rescue 與甲基輔助評估在 paired 與 TO 情境下的效果
處理範圍:
  - final VCF candidate pool
  - caller-first / methylation-support / veto 比較
  - HCC1395 5kHz paired / TO 與 HCC1395_DORADO paired
關聯檔案:
  - docs/reports/validated/2026/03/20260308_HCC1395_5kHz_TO_borderline_rescue特徵證據鏈整理_01.md
  - docs/reports/validated/2026/03/20260309_甲基rescue是否穩定有效的跨樣本判讀_01.md
-->

# 2026-03-08 ClairS 邊緣 TP rescue 與甲基輔助評估

## 研究目標

驗證以下假設：

1. 是否能從 `ClairS` / `ClairS-TO` 的 **final VCF** 找出「信心度介於邊緣、可能被後段流程判成 FN」的候選位點。
2. 是否能在 `caller-first` 的前提下，讓甲基訊號作為第二層 `support` 或 `veto`，在 **少引入 FP** 的條件下救回 TP。

本輪固定限制：

- 只用 **final VCF**，不使用 pileup 中間 candidate / tensor。
- paired 使用 `ClairS final output.vcf.gz`。
- tumor-only 使用 `ClairS-TO final snv.vcf.gz`，且這份 VCF 來自 **pileup 路線**，不是 full model。
- 評估重點是 `HCC1395 5kHz paired`、`HCC1395 5kHz TO`，再用 `HCC1395_DORADO paired` 做第一個交叉驗證。

## 本輪實作

### 新增腳本

1. [extract_borderline_rescue_candidates.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/extract_borderline_rescue_candidates.py)
2. [evaluate_rescue_with_methylation.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/evaluate_rescue_with_methylation.py)
3. [compare_rescue_strategies.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/compare_rescue_strategies.py)

### 本輪修正

1. `evaluate_rescue_with_methylation.py` 改為只保留 `caller_lost_tp / caller_removed_fp` 再 join。
2. `rescue_joined_features.tsv` 改為只輸出：
   - `candidate_eligible`
   - 或已有甲基 summary overlap 的列
3. `HCC1395 TO` 首次 full join 會產出過大無效中間檔，因此本輪改以 **caller-only streaming 統計** 保留有效結果，並明確標記 `TO methylation overlap = 0/0`。

## 使用資料與來源

### Paired 5kHz

- Caller VCF: [output.vcf.gz](/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz)
- Truth VCF: [high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz](/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz)
- Truth BED: [High-Confidence_Regions_v1.2.bed](/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed)
- LongPhase-S baseline TP: [filtered_snv_tp.vcf.gz](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395/20260211/longphase_s/filtered_snv_tp.vcf.gz)
- LongPhase-S baseline FP: [filtered_snv_fp.vcf.gz](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395/20260211/longphase_s/filtered_snv_fp.vcf.gz)
- candidate-specific InterSubMod summary:
  - [lost TP significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_methylation_hcc1395/intersubmod_tp/significance_summary.csv)
  - [removed FP significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_methylation_hcc1395/intersubmod_fp/significance_summary.csv)

### Paired DORADO

- Caller VCF: [output.vcf.gz](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz)
- LongPhase-S baseline TP: [filtered_snv_tp.vcf.gz](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395_DORADO/20260211/longphase_s/filtered_snv_tp.vcf.gz)
- LongPhase-S baseline FP: [filtered_snv_fp.vcf.gz](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395_DORADO/20260211/longphase_s/filtered_snv_fp.vcf.gz)

### TO 5kHz

- Caller VCF: [snv.vcf.gz](/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz)
- LongPhase-TO baseline TP: [tp.vcf](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step04_benchmark_longphase_to/tp.vcf)
- LongPhase-TO baseline FP: [fp.vcf](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step04_benchmark_longphase_to/fp.vcf)
- 現有 baseline InterSubMod summary 與 downstream lost/removed 候選的 overlap：
  - lost TP: `0`
  - removed FP: `0`

## 主要輸出

### Candidate pool

1. [HCC1395 paired candidate_group_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_paired/extract/candidate_group_summary.tsv)
2. [HCC1395 TO candidate_group_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_to/extract/candidate_group_summary.tsv)
3. [HCC1395_DORADO paired candidate_group_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_DORADO_paired/extract/candidate_group_summary.tsv)

### Rescue rule comparison

1. [HCC1395 paired rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_paired/eval/rescue_rule_comparison.tsv)
2. [HCC1395 TO rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_to/eval/rescue_rule_comparison.tsv)
3. [HCC1395_DORADO paired rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_DORADO_paired/eval/rescue_rule_comparison.tsv)
4. [strategy_best_by_mode.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/strategy_compare/strategy_best_by_mode.tsv)

## 結果

### 1. HCC1395 5kHz paired

candidate pool：

- `caller_lost_tp = 1052`，其中 `candidate_eligible = 920`
- `caller_removed_fp = 12974`，其中 `candidate_eligible = 5266`
- `candidate_eligible` 主要 FILTER：
  - `caller_lost_tp`: `LowQual=815`, `PASS=105`
  - `caller_removed_fp`: `LowQual=4729`, `PASS=537`

caller-only：

- 依計劃排序的最佳 safe rule：
  - `caller_candidate_gq_ge_15`
  - `TP rescued = 106`
  - `FP reintroduced = 75`
  - `F1 delta = +0.000825`
- 若以 **最佳 delta F1** 看：
  - `caller_any_gq_ge_20`
  - `TP rescued = 59`
  - `FP reintroduced = 8`
  - `F1 delta = +0.000871`

caller + methylation：

- 目前 **沒有**超過 caller-only。
- `methylation veto` 在這一版 final VCF 邊緣候選上沒有額外增益。
- `methylation support` 目前是負增益：
  - `support_pairwise_ge_020`: `50 TP / 130 FP`, `F1 -0.000763`
  - `support_strong_or_subclone`: `25 TP / 106 FP`, `F1 -0.000881`

重要解讀：

- 這說明在 paired 5kHz 的 final VCF 邊緣救援裡，**caller GQ 仍是主訊號**。
- 甲基資訊目前還沒有穩定提供「救 TP 又不拉回太多 FP」的第二層增益。
- 既有甲基 summary 只覆蓋：
  - `111 / 1052` 個 lost TP
  - `802 / 12974` 個 removed FP
- 也就是說，本輪 paired methylation rescue 只能視為 **部分覆蓋的 first pass**，還不是完整答案。

### 2. HCC1395_DORADO paired

candidate pool：

- `caller_lost_tp = 1489`，其中 `candidate_eligible = 1122`
- `caller_removed_fp = 2533`，其中 `candidate_eligible = 1658`
- `candidate_eligible` 主要 FILTER：
  - `caller_lost_tp`: `LowQual=1029`, `PASS=93`
  - `caller_removed_fp`: `LowQual=1332`, `PASS=326`

caller-only：

- 依計劃排序的最佳 safe rule：
  - `caller_candidate_gq_ge_15`
  - `TP rescued = 97`
  - `FP reintroduced = 88`
  - `F1 delta = +0.000502`
- 若以 **最佳 delta F1** 看：
  - `caller_any_gq_ge_20`
  - `TP rescued = 53`
  - `FP reintroduced = 7`
  - `F1 delta = +0.000782`

caller + methylation：

- 本輪沒有 `DORADO paired` 的 candidate-specific `lost_tp / removed_fp` InterSubMod summary。
- 因此這裡的 methylation 欄位不能當成有效 negative result，只能記成：
  - `caller-only 已可提升`
  - `methylation rescue 尚未被正確測到`

### 3. HCC1395 5kHz TO

candidate pool：

- `caller_lost_tp = 2108`，其中 `candidate_eligible = 773`
- `caller_removed_fp = 2720008`，其中 `candidate_eligible = 298`
- `caller_removed_fp` 的 `pon_hit_fraction = 0.999019`
- `candidate_eligible` 主要 FILTER：
  - `caller_lost_tp`: `PASS=675`, `LowQual=76`, `LowQual;NoAncestry=15`, `LowQual;MultiHap=7`
  - `caller_removed_fp`: `PASS=213`, `LowQual=72`, `LowQual;MultiHap=8`, `LowQual;NoAncestry=5`

caller-only：

- 本輪最佳 safe rule：
  - `caller_candidate_qual_ge_10_or_gq_ge_20`
  - `TP rescued = 491`
  - `FP reintroduced = 118`
  - `F1 delta = +0.006824`
- 第二名：
  - `caller_candidate_gq_ge_15`
  - `TP rescued = 365`
  - `FP reintroduced = 71`
  - `F1 delta = +0.005233`
- `caller_candidate_gq_ge_20`
  - `TP rescued = 143`
  - `FP reintroduced = 38`
  - `F1 delta = +0.001966`
- `Verdict support`
  - `0 TP / 0 FP`
  - 在這一版候選池中沒有提供額外價值

重要解讀：

- TO final VCF 的救援空間 **明顯大於 paired**。
- 而且 TO 的 candidate-eligible downstream pool 主要是 `PASS` + 一小批 `LowQual`，不是大量 PoN/NonSomatic 雜訊。
- 這支持一個很重要的方向：
  - **先從 final VCF 的 borderline PASS / LowQual pool 做 caller-side rescue，可能比直接從甲基端硬救更有效。**

TO 的甲基結論：

- 本輪不能解讀。
- 因為 `caller_lost_tp / caller_removed_fp` 與現有 baseline InterSubMod summary 的 overlap 為：
  - `lost TP = 0`
  - `removed FP = 0`
- 所以若要回答「TO 下甲基資料能否幫助 borderline TP rescue」，下一步必須：
  - 先用這批 TO candidates 單獨跑 InterSubMod
  - 再做 `caller + methylation support/veto`

## 總結

### 已確認

1. 只看 **final VCF** 也能建立有用的邊緣 rescue 候選池。
2. 在 `5kHz paired` 與 `DORADO paired`，caller-only 的 `GQ` 型規則已可穩定救回一批 TP。
3. 在 `5kHz TO`，caller-only rescue 的效果更強，最佳 safe rule 可達 `F1 +0.006824`。
4. 目前最穩的 rescue 主訊號仍在 caller 端，不在甲基端。

### 尚未確認

1. 甲基資料是否能在 **完整 final VCF borderline pool** 上提供真正的第二層 rescue 增益。
2. `low VAF + high AlleleDelta + low CramersV` 這類甲基 artifact veto，是否適合轉用在 **TP rescue 的 veto**，而不是只做 FP triage。
3. TO 下的甲基 rescue 是否能在 candidate-specific run 後優於 caller-only。

## 建議下一步

1. 第一優先：
   - 對 `HCC1395 5kHz TO` 的 `caller_lost_tp / caller_removed_fp` candidate pool 單獨跑 InterSubMod。
2. 第二優先：
   - 在 paired 5kHz 上補齊 **final VCF 全候選池** 的 candidate-specific InterSubMod coverage，不只沿用舊的 PASS 子集。
3. 第三優先：
   - 對 `DORADO paired` 跑同型 candidate-specific InterSubMod，確認甲基 rescue 是否真為平台不適用。
4. 若要先進 pipeline：
   - 先加 `RescueCandidate` annotation，不先加 hard filter / hard keep。
