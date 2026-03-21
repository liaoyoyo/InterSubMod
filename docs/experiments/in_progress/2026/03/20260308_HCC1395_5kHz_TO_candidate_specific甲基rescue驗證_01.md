<!--
建立時間: 2026-03-08 04:20
目標: 對 HCC1395 5kHz tumor-only 的 caller_lost_tp / caller_removed_fp candidate pool 單獨跑 InterSubMod，正式回答 TO 下甲基資料能否幫 borderline TP rescue
處理範圍:
  - HCC1395 5kHz TO candidate pool export
  - candidate-specific InterSubMod
  - label/cluster agreement 驗證
  - caller-only vs methylation rescue 規則比較
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_candidate_pool_vcfs.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/evaluate_rescue_with_methylation.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/
-->

# HCC1395 5kHz TO candidate-specific 甲基 rescue 驗證

## 1. 研究問題

本輪要正式回答：

1. 在 `HCC1395 5kHz tumor-only` 下，若先從 `ClairS-TO final snv.vcf.gz` 的 borderline pool 挑出 `caller_lost_tp / caller_removed_fp` 候選，甲基資料是否能幫忙做 TP rescue。
2. `label-first / cluster-first` 的訊號是否只是存在，還是能進一步提升 borderline rescue 的品質。

固定前提：

1. 使用的是 [snv.vcf.gz](/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz)。
2. 這份 VCF 來自 `ClairS-TO pileup` 路線，不是 full model。
3. 本輪不再用 baseline kept-set summary 硬 join，而是實際對 candidate pool 單獨跑 InterSubMod。

## 2. 本輪實作

### 2.1 新增腳本

1. [export_candidate_pool_vcfs.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_candidate_pool_vcfs.py)

### 2.2 實際執行流程

1. 從 [borderline_candidate_pool.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_to/extract/borderline_candidate_pool.tsv) 匯出 candidate-specific VCF
2. 用既有 [tumor_tagged.bam](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam) 單獨跑 InterSubMod
3. 用 [validate_method_design.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py) 產出 `label_cluster_agreement.tsv`
4. 用 [evaluate_rescue_with_methylation.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/evaluate_rescue_with_methylation.py) 重算 TO rescue 規則

## 3. 主要輸出

### 3.1 Candidate export

1. [candidate_vcf_export_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_vcf_export_summary.tsv)
2. [candidate_lost_tp.vcf](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_lost_tp.vcf)
3. [candidate_removed_fp.vcf](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_removed_fp.vcf)

### 3.2 Candidate-specific InterSubMod

1. [TP significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_tp/significance_summary.csv)
2. [FP significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_fp/significance_summary.csv)
3. [TP label_first_metrics.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_tp/label_first_metrics.tsv)
4. [FP label_first_metrics.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_fp/label_first_metrics.tsv)

### 3.3 Method design / rescue eval

1. [label_cluster_agreement.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/design_validation/label_cluster_agreement.tsv)
2. [method_design_validation.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/design_validation/method_design_validation.tsv)
3. [rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_rule_comparison.tsv)
4. [rescue_rule_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_rule_summary.md)
5. [rescue_joined_features.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv)

## 4. 實際完成與正確性確認

### 4.1 candidate VCF 匯出成功

根據 [candidate_vcf_export_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_vcf_export_summary.tsv)：

1. `caller_lost_tp` target `773`，exported `773`，missing `0`
2. `caller_removed_fp` target `298`，exported `298`，missing `0`

### 4.2 candidate-specific InterSubMod 已實跑完成

實際可分析 region：

1. `caller_lost_tp`: `675 / 773 = 87.32%`
2. `caller_removed_fp`: `213 / 298 = 71.48%`

這代表本輪 `TO methylation rescue` 不再是上次的 `0/0 overlap`，而是正式可解讀的 candidate-specific 結果。

### 4.3 空間與輸出量

1. 本輪新 round 輸出總量約 `321M`
2. 截至完成後磁碟仍約：
   - `/big8_disk`: `226G`
   - `/bip8_disk`: `196G`
3. 本輪沒有新增新的 `200G+` 級 BAM，不需要你手動清理

## 5. 主要結果

### 5.1 最佳整體 rescue 仍是 caller-only

根據 [rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_rule_comparison.tsv)：

最佳 safe rule 仍是 caller-only：

1. `caller_candidate_qual_ge_10_or_gq_ge_20`
2. `TP rescued = 491`
3. `FP reintroduced = 118`
4. `F1 delta = +0.006824`

這表示在 `5kHz TO` 的 borderline rescue，主訊號仍然在 caller 端。

### 5.2 甲基資料不是全無效，確實可以幫忙 rescue

最好的甲基 support 規則是：

1. `caller_candidate_gq_ge_10__support_pairwise_ge_020`
2. `TP rescued = 300`
3. `FP reintroduced = 68`
4. `F1 delta = +0.004219`

這個結果比 baseline 明顯好，說明：

1. `TO` 下甲基資料**可以**幫忙 borderline TP rescue
2. 但目前還沒有超過最佳 caller-only 規則

### 5.3 agreement-based support 比單純 Strong/Subclone 更乾淨

agreement 納入後：

1. `caller_candidate_gq_ge_10__support_agreement_positive`
2. `TP rescued = 148`
3. `FP reintroduced = 25`
4. `F1 delta = +0.002163`

對照不看 agreement 的：

1. `caller_candidate_gq_ge_10__support_strong_or_subclone`
2. `TP rescued = 149`
3. `FP reintroduced = 30`
4. `F1 delta = +0.002134`

也就是說：

1. 把 `label_upgrade / consistent_strong / consistent_subclone` 當成正向 support
2. 可以在幾乎不損 TP 的前提下，少帶回 `5` 個 FP

這表示 `label-first / cluster-first` 的交叉訊號不是只有描述性，對 rescue precision 確實有一點幫助。

### 5.4 artifact veto 對 TO rescue 幫助有限

本輪測到：

1. `caller_candidate_gq_ge_15__veto_lowvaf_highadelta_lowcv`
   - 與 `caller_candidate_gq_ge_15` 完全相同：`365 TP / 71 FP`
2. `caller_candidate_gq_ge_15__veto_combined_artifact`
   - `308 TP / 65 FP`
   - 雖然少了 `6` 個 FP，但同時少了 `57` 個 TP

這說明兩件事：

1. 我們先前對 TO artifact 很有效的 `low VAF + high AlleleDelta + low CramersV`
   - 不太落在這批 `caller_lost_tp` 候選裡
2. `artifact veto` 比較適合 FP triage，不適合直接拿來做 TO borderline TP rescue

### 5.5 label/cluster 訊號在 candidate-specific pool 不是空的

根據 [label_cluster_agreement.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/design_validation/label_cluster_agreement.tsv)：

`lost_tp` 可分析 `675` 筆中：

1. `consistent_weak_or_noise = 260`
2. `label_upgrade = 126`
3. `cluster_plus_weak_label = 106`
4. `consistent_strong = 86`
5. `conflict = 83`

`removed_fp` 可分析 `213` 筆中：

1. `consistent_weak_or_noise = 98`
2. `consistent_strong = 32`
3. `label_upgrade = 30`
4. `conflict = 24`
5. `cluster_plus_weak_label = 23`

這表示：

1. `lost_tp` 中確實存在一批被 label/cluster 正向支持的候選
2. 但 `removed_fp` 也不是完全沒有同型訊號
3. 所以目前最合理的定位仍是 `support / ranking / annotation`，不是 hard keep

## 6. 結論

### 6.1 正式答案

對「`TO` 下甲基資料能不能幫 borderline TP rescue」這題，現在可以正式回答：

1. **可以**
2. 但目前是 **次要輔助訊號**
3. 還沒有超過最佳 caller-only rescue

### 6.2 本輪最重要的研究結論

1. `TO` candidate-specific methylation rescue 已經被正確量測，不再是 `0/0 overlap`
2. 甲基 support 規則有實際增益：
   - 最佳純甲基 support：`300 TP / 68 FP`, `F1 +0.004219`
3. `label-first + cluster-first` 的 agreement 可提升 support precision：
   - `agreement_positive`: `148 TP / 25 FP`
   - 比單純 `Strong/Subclone` 更乾淨
4. 最佳整體 rescue 仍是 caller-only：
   - `491 TP / 118 FP`, `F1 +0.006824`
5. TO artifact veto 目前不適合拿來救 TP：
   - 它比較像 FP triage 工具，不是 TP rescue 主工具

### 6.3 方法學定位

目前最合理的流程定位是：

1. `caller-first`
2. `methylation support second`
3. `artifact veto` 作為獨立 FP triage 線，不和 TP rescue 混在一起

## 7. 下一步建議

1. 第一優先：
   - 對 `HCC1395_DORADO TO` 跑同型 candidate-specific InterSubMod，確認這個結果是否能跨平台重現
2. 第二優先：
   - 對 `agreement_positive rescued TP` 與 `support_pairwise_ge_020` 的 top TP / FP 做 diagnostics，分辨哪些是可信 support、哪些仍是弱支持
3. 第三優先：
   - 將 `TO rescue` 的甲基結果先寫成 annotation / ranking，不先改成 hard keep 規則
