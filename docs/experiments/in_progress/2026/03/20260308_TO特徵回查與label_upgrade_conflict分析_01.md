<!--
建立時間: 2026-03-08 01:26
目標: 回查 HCC1395 5kHz tumor-only LongPhase-TO + InterSubMod 中的 label_upgrade/conflict、Weak->Strong / Noise->Strong，以及 low VAF + high AlleleDelta 在 TP/FP 下的細分狀態
處理範圍:
  - HCC1395 5kHz TO run: 20260307_hcc1395_to_pilot_1
  - label_cluster_agreement.tsv
  - significance_summary.csv
  - benchmark split TP/FP VCF
  - TO 特徵回查腳本與輸出
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_to_feature_recall.py
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/label_cluster_agreement.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step05_intersubmod/intersubmod_tp/significance_summary.csv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step05_intersubmod/intersubmod_fp/significance_summary.csv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_summary.md
-->

# TO 特徵回查與 label_upgrade / conflict 分析

## 1. 本次要回答的問題

針對 `HCC1395 5kHz tumor-only LongPhase-TO + InterSubMod`，本次回查聚焦 3 件事：

1. `label_upgrade=2756` 與 `conflict=3227` 中，哪些主要是 `TP`，哪些可能富集 `FP`
2. `Weak->Strong / Noise->Strong` 在 TO 下是否仍像 paired pure 一樣值得警戒
3. `low VAF + high AlleleDelta` 在 TO 下是否仍有效，以及它主要落在哪些 label/shift 類型

## 2. 本次新增的可重跑工具

- [analyze_to_feature_recall.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_to_feature_recall.py)

輸出位置：

- [to_feature_summary.md](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_summary.md)
- [to_feature_focus_matrix.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_focus_matrix.tsv)
- [to_feature_shift_breakdown.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_shift_breakdown.tsv)
- [to_feature_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_rule_comparison.tsv)
- [to_feature_top_fp_candidates.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_top_fp_candidates.tsv)
- [to_feature_top_tp_controls.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_top_tp_controls.tsv)

## 3. 主要結果

### 3.1 `label_upgrade` 與 `Weak/Noise->Strong` 在 TO 下大多是 TP

根據 [to_feature_focus_matrix.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_focus_matrix.tsv)：

1. `agreement_label_upgrade`
   - `TP=2526`
   - `FP=230`
   - `FP rate=8.35%`
2. `shift_Weak_to_Strong`
   - `TP=1849`
   - `FP=176`
   - `FP rate=8.69%`
3. `shift_Noise_to_Strong`
   - `TP=207`
   - `FP=16`
   - `FP rate=7.17%`

這代表在 `5kHz TO` 場景，`label_upgrade`、`Weak->Strong`、`Noise->Strong` **不能直接視為主要 artifact 池**。和 paired pure 的直覺不同，這些集合在 TO 下以 `TP` 為主。

### 3.2 `conflict` 也不是好的直接刪除目標

`agreement_conflict`：

- `TP=2465`
- `FP=762`
- `FP rate=23.61%`

雖然 `FP rate` 高於 `label_upgrade`，但仍以 `TP` 為主。若直接把整包 `conflict` 當作刪除規則，會嚴重傷害 recall。

更細的 breakdown 見 [to_feature_shift_breakdown.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_shift_breakdown.tsv)：

1. `conflict + Noise->Weak`
   - `TP=2397`
   - `FP=751`
   - `FP rate=23.86%`
2. `conflict + Strong->Subclone`
   - `TP=59`
   - `FP=10`
3. `conflict + Subclone->Strong`
   - `TP=9`
   - `FP=1`

也就是說，`conflict` 的主體其實是 `Noise->Weak`，而不是 `Strong` 類升級。

### 3.3 真正有用的 TO artifact 特徵仍是 `low VAF + high AlleleDelta`

`rule_low_vaf_high_ad`：

- `TP=3`
- `FP=36`
- `FP rate=92.31%`
- `median VAF=0.1043`
- `median AlleleDelta=0.2754`
- `median GQ=7`
- `median QUAL=7.1712`

加上 `CramersV < 0.05` 後：

- `rule_low_vaf_high_ad_cv`
- `TP=2`
- `FP=36`
- `FP rate=94.74%`

這是本輪最關鍵的新結論：在 TO 下，**有用的 FP triage 特徵仍然是低 VAF + 高 AlleleDelta**，而且比整包 `label_upgrade` / `conflict` 更乾淨。

### 3.4 但這個 TO artifact 特徵幾乎不落在 `label_upgrade/conflict`

這也是本次最值得記住的地方。

`combo_label_upgrade_low_vaf_high_ad`：

- `TP=0`
- `FP=1`

`combo_conflict_low_vaf_high_ad`：

- `TP=0`
- `FP=0`

`combo_Weak_to_Strong_low_vaf_high_ad`：

- `TP=0`
- `FP=1`

`combo_Noise_to_Strong_low_vaf_high_ad`：

- `TP=0`
- `FP=0`

也就是說，**TO 下的 low VAF + high AlleleDelta 幾乎是另一個正交特徵池，不是由 label_upgrade/conflict 在撐住**。

### 3.5 TO 下 low VAF + high AlleleDelta 主要落在 `cluster_only / Strong->Weak`

在 `rule_low_vaf_high_ad` 中，top class shift 為：

- `Strong->Weak:27`
- `unchanged:9`
- `Strong->Noise:2`

top agreement type 為：

- `cluster_only:27`
- `consistent_strong:7`
- `consistent_weak_or_noise:2`

所以在 TO 下，真正更像 artifact 的區域是：

1. `cluster-first` 本來叫它 `Strong`
2. 但 `label-first` 並沒有跟上，導致 `Strong->Weak`
3. 同時 caller 特徵又呈現 `低 VAF + 高 AlleleDelta + 低 GQ/QUAL`

這一型比 `Weak->Strong / Noise->Strong` 更值得優先追查。

## 4. 規則對 F1 的影響

根據 [to_feature_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_rule_comparison.tsv)：

1. `rule_low_vaf_high_ad_cv`
   - 移除 `2 TP`
   - 移除 `36 FP`
   - `F1 0.713000 -> 0.713261`
   - `delta=+0.000261`
2. `rule_low_vaf_high_ad`
   - 移除 `3 TP`
   - 移除 `36 FP`
   - `delta=+0.000245`
3. `rule_legacy_combo`
   - 移除 `3 TP`
   - 移除 `24 FP`
   - `delta=+0.000137`

反過來看：

1. `agreement_label_upgrade` 若整包移除
   - `TP removed=2526`
   - `FP removed=230`
   - `F1 delta=-0.040177`
2. `agreement_conflict` 若整包移除
   - `TP removed=2465`
   - `FP removed=762`
   - `F1 delta=-0.034434`
3. `shift_Weak_to_Strong` 若整包移除
   - `TP removed=1849`
   - `FP removed=176`
   - `F1 delta=-0.029070`
4. `shift_Noise_to_Strong` 若整包移除
   - `TP removed=207`
   - `FP removed=16`
   - `F1 delta=-0.003240`

所以這輪 TO 回查的結論非常明確：

- `low VAF + high AlleleDelta` 仍然有用
- `label_upgrade / conflict / Weak->Strong / Noise->Strong` **不適合直接當刪除規則**

## 5. 代表性候選位點

根據 [to_feature_top_fp_candidates.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_top_fp_candidates.tsv)：

目前最值得先回查的 `FP` 有：

1. `chr1:168472226:T:A`
   - `label_upgrade`
   - `Weak->Strong`
   - `VAF=0.1043`
   - `QUAL=6.0834`
   - `GQ=6`
   - `AlleleDelta=0.1688`
   - 是本輪唯一同時命中 `label_upgrade + Weak->Strong + lowVAF/highAD` 的 `FP`
2. 多數其他 `lowVAF/highAD` 的 `FP`
   - 並不屬於 `label_upgrade`
   - 反而多屬於 `cluster_only + Strong->Weak`

對照 [to_feature_top_tp_controls.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_top_tp_controls.tsv)，也確實存在少數 `TP` 會命中同類規則，例如：

1. `chr15:97413083:C:T`
2. `chr8:120817039:A:G`
3. `chr7:153260385:C:T`

所以後續若要把這個規則往前推，仍需要：

1. 加入更細的 caller 條件
2. 或補 `HP balance / NA HP ratio / read-level snapshot`
3. 避免把這 `2~3` 個 TP 一起刪掉

## 6. 本次結論

1. `TO` 下的 `label_upgrade` 和 `Weak/Noise->Strong` **大多是 TP**
2. `conflict` 也不是主要 artifact 集合，主體是 `Noise->Weak`
3. 真正持續有效的 TO artifact 特徵仍是：
   - `low VAF + high AlleleDelta`
   - 若再加 `low CramersV` 稍微更乾淨
4. 但這個特徵池在 TO 下 **幾乎不靠 `label_upgrade/conflict` 撐住**
5. 它主要落在：
   - `cluster_only`
   - `Strong->Weak`
   - `低 GQ / 低 QUAL`

## 7. 下一步建議

1. 下一輪 TO diagnostics 應優先看：
   - `cluster_only + Strong->Weak + low VAF + high AlleleDelta`
   - 而不是先看整包 `Weak->Strong / Noise->Strong`
2. 對 `chr1:168472226:T:A` 與 top `cluster_only + Strong->Weak` FP 做 `samtools` / matrix diagnostics
3. 將同一支 [analyze_to_feature_recall.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_to_feature_recall.py) 直接套到未來的 `HCC1395_DORADO TO`
4. 若 `DORADO TO` 也出現同樣模式，才考慮把：
   - `low VAF + high AlleleDelta (+ low CramersV)`
   正式升級成 tumor-only 的 `ArtifactSuspect` 流程候選
