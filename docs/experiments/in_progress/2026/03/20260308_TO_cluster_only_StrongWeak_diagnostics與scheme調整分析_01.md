<!--
建立時間: 2026-03-08 11:10
目標: 針對 HCC1395 5kHz tumor-only 中原先定義為 cluster_only + Strong->Weak + lowVAF/highAD 的 top FP 做 read-level diagnostics，並檢查 label-first / cluster-first 分類語意調整是否更合理
處理範圍:
  - HCC1395 5kHz TO run: 20260307_hcc1395_to_pilot_1
  - top FP / TP matrix diagnostics
  - samtools snapshot
  - SignificanceAnalyzer / LabelTest / validate_method_design.py
  - verification scheme adjustment analysis
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_region_diagnostics.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/region_samtools_snapshot.sh
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/summarize_samtools_snapshots.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_to_verification_scheme_adjustments.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_clusteronly_strongweak_fp_diagnostics/
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_verification_scheme_adjustments_hcc1395_5khz/
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_validate_method_design_refined_hcc1395_5khz/
-->

# TO `cluster_only + Strong->Weak` diagnostics 與 scheme 調整分析

## 1. 這輪要回答的問題

針對 `HCC1395 5kHz tumor-only LongPhase-TO + InterSubMod`，本輪聚焦 3 個問題：

1. 原先在 TO 特徵回查中被描述為 `cluster_only + Strong->Weak + lowVAF/highAD` 的 top FP，是否有共通的 read-level 異常型態。
2. 目前 core 與分析腳本對 `label-first` / `cluster-first` 的語意是否有混淆，尤其是 `Strong` 的定義是否過鬆。
3. 若調整分類語意或驗證 scheme，是否能比現行 `low VAF + high AlleleDelta + low CramersV` 規則拿到更好的結果。

## 2. 重要前提：本輪 TO 輸入是 `ClairS-TO pileup`，不是 full model

根據 [run.log](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/run.log) 與 [run_longphase_to_intersubmod_pilot.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_longphase_to_intersubmod_pilot.sh)：

1. `Step 01` 使用的是 `ClairS-TO` 的 `ont_r10_dorado_sup_5khz_ssrs` pileup 路線。
2. `LongPhase-TO phase` 吃的是這個 pileup 產生的 [snv.vcf.gz](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step01_clairs_to/snv.vcf.gz)。
3. `haplotagged BAM` 來自原始 `5kHz MM/ML tumor BAM` 經 `LongPhase-TO haplotag` 處理後的 [tumor_tagged.bam](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam)。

因此，本輪 diagnostics 與結論先只對：

- `TO`
- `5kHz`
- `ClairS-TO pileup`
- `LongPhase-TO`
- `InterSubMod`

這條路線成立，不能直接外推到 full model 或 paired pure。

## 3. 實際執行內容

### 3.1 選取的候選位點

原始候選來自 [to_feature_per_region.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_feature_recall_hcc1395_5khz/to_feature_per_region.tsv)，條件為：

- `agreement_type = cluster_only`
- `class_shift = Strong->Weak`
- `low_vaf_high_ad = True`

本輪實際匯出的 5 個 FP 與 2 個 TP 對照如下：

| region | scope | VAF | AlleleDelta | QUAL | GQ | PairwiseMedianDist |
| --- | --- | --- | --- | --- | --- | --- |
| `chr6:165980348:T:C` | FP | `0.0926` | `0.7003` | `10.7685` | `10` | `0.1555` |
| `chr11:18904680:T:G` | FP | `0.0893` | `0.4751` | `5.1500` | `5` | `0.1538` |
| `chr11:18905130:G:C` | FP | `0.0862` | `0.4216` | `6.1602` | `6` | `0.1638` |
| `chr11:18905913:T:C` | FP | `0.1111` | `0.4178` | `11.5394` | `11` | `0.1781` |
| `chr11:18904511:G:A` | FP | `0.1017` | `0.4043` | `8.8101` | `8` | `0.1624` |
| `chr7:153260385:C:T` | TP | `0.1333` | `0.1582` | `8.0416` | `8` | `0.2876` |
| `chr15:97413083:C:T` | TP | `0.1429` | `0.1646` | `15.3639` | `15` | `0.2525` |

### 3.2 新增與修正的工具

1. 新增 [analyze_to_verification_scheme_adjustments.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_to_verification_scheme_adjustments.py)
   - 用既有 `significance_summary.csv + to_feature_per_region.tsv + benchmark_comparison.tsv` 重算多個 verification scheme
2. 修正 [summarize_samtools_snapshots.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/summarize_samtools_snapshots.py)
   - 原本只抓區間第一列，無法對準目標位點
   - 現已改成解析 `region_id`，回到目標 SNV 位置計算：
     - `target_depth`
     - `target_ref_count`
     - `target_alt_count`
     - `target_alt_fraction`
   - 也修正了遞迴掃描 `fp/tp` 分層目錄的問題
3. 修正 [validate_method_design.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py)
   - 將原先過於粗糙的 `cluster_only` 進一步拆出 `cluster_plus_weak_label`
   - 讓「完全沒有 label」與「只有弱 allele label」不再混在同一類

## 4. read-level diagnostics 的主要觀察

輸出主目錄：

- [20260308_to_clusteronly_strongweak_fp_diagnostics](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_clusteronly_strongweak_fp_diagnostics)
- [snapshot_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_clusteronly_strongweak_fp_diagnostics/snapshot_summary/snapshot_summary.tsv)

### 4.1 4 個 top FP 集中在 chr11 同一個局部區段

下列 4 個 FP：

- `chr11:18904511:G:A`
- `chr11:18904680:T:G`
- `chr11:18905130:G:C`
- `chr11:18905913:T:C`

全部落在 `chr11:18904511-18905913`，總 span 只有 `1402 bp`。

這表示它們不像彼此獨立的隨機 noise，而更像同一個局部 artifact block 或局部 read-selection / haplotype 特性造成的連續假訊號。

### 4.2 `HP skew` 很明顯，但單獨看並不能區分 FP / TP

根據 [snapshot_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_clusteronly_strongweak_fp_diagnostics/snapshot_summary/snapshot_summary.tsv)：

1. 5 個 FP 都呈現：
   - `single_hp_dominant`
   - `haplotype_skewed`
2. 但 2 個 TP 對照同樣也呈現：
   - 極度單一 HP 優勢
   - 很低 `NA HP fraction`
   - 高 `collapsed_hp_balance_delta`

例如：

| region | scope | target_alt_fraction | hp1 | hp2 | na_hp_fraction | hp_balance_delta |
| --- | --- | --- | --- | --- | --- | --- |
| `chr6:165980348:T:C` | FP | `0.0862` | `62` | `0` | `0.1014` | `1.0000` |
| `chr11:18904680:T:G` | FP | `0.0833` | `7` | `69` | `0.0706` | `0.8158` |
| `chr15:97413083:C:T` | TP | `0.1586` | `1` | `118` | `0.0083` | `0.9832` |
| `chr7:153260385:C:T` | TP | `0.1267` | `112` | `0` | `0.0667` | `1.0000` |

所以本輪一個很重要的結論是：

- `HP skew` 不能單獨當 artifact 規則
- 它最多只能作為輔助解釋訊號

### 4.3 真正比較有辨識力的是 `高 AlleleDelta + 低 PairwiseMedianDist + 很低 alt fraction`

FP 與 TP 對照的差異比較明顯的是：

1. FP 的 `AlleleDelta` 顯著更高
   - FP: `0.4043 ~ 0.7003`
   - TP: `0.1582 ~ 0.1646`
2. FP 的 `PairwiseMedianDist` 明顯更低
   - FP: `0.1538 ~ 0.1781`
   - TP: `0.2525 ~ 0.2876`
3. FP 的 `target_alt_fraction` 更低
   - FP: `0.0794 ~ 0.1071`
   - TP: `0.1267 ~ 0.1586`
4. FP 的 `QUAL/GQ` 整體更低
   - FP: `QUAL 5.15 ~ 11.54`, `GQ 5 ~ 11`
   - TP: `QUAL 8.04 ~ 15.36`, `GQ 8 ~ 15`

這表示 TO 下比較可疑的不是「單一 HP」，而是：

1. allele 標籤切得很開
2. 但 methylation distance 並沒有真的拉開到同等程度
3. caller 自己的支持度又偏弱

## 5. 程式碼與方法學檢查

### 5.1 現行 core 的 `Strong` 定義確實偏鬆

根據 [SignificanceAnalyzer.cpp](/big8_disk/liaoyoyo2001/InterSubMod/src/core/SignificanceAnalyzer.cpp)：

1. `cluster_significant` 目前是：
   - `passed_gate`
   - 且 `global_alt/global_hp fisher p <= 0.05`
2. `Strong` 目前只要求：
   - `cluster_significant`
   - 且 `label_significant = any_significant`

也就是說，**只要有任何 label signal，就可能進到 `Strong`**，並沒有區分：

- `HP + Allele` 都支持
- 只有 `HP`
- 只有 `Allele`

### 5.2 本輪 top FP 與 TP 對照，全部都是「只有 allele sig，沒有 HP sig」

對 7 個 diagnostics 候選回查 [significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step05_intersubmod/intersubmod_fp/significance_summary.csv) / [TP significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step05_intersubmod/intersubmod_tp/significance_summary.csv) 後，這 7 個位點有共同特徵：

1. `ClusterPermanovaValid = false`
2. `ClusterPermanovaP = 1.0`
3. `HPMergedSig = false`
4. `HPFineSig = false`
5. `LabelHPPermanovaValid = false`
6. `AlleleSig = true`
7. `LabelAllelePermanovaP = 0.001`
8. 但 `VerificationClass = Strong`

也就是說，這些位點被 core 叫做 `Strong`，主要是因為：

- global / local cluster 還有訊號
- allele label 有顯著
- 但 HP 結構其實完全沒有支持

這正是 `Strong` 被混進「只有弱 label」案例的直接證據。

## 6. 多次 scheme 嘗試結果

輸出主目錄：

- [20260308_to_verification_scheme_adjustments_hcc1395_5khz](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_verification_scheme_adjustments_hcc1395_5khz)
- [verification_scheme_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_verification_scheme_adjustments_hcc1395_5khz/verification_scheme_rule_comparison.tsv)
- [verification_scheme_class_counts.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_verification_scheme_adjustments_hcc1395_5khz/verification_scheme_class_counts.tsv)
- [verification_scheme_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_verification_scheme_adjustments_hcc1395_5khz/verification_scheme_summary.md)

### 6.1 嘗試的 4 個 scheme

1. `current`
   - 直接使用現行 `VerificationClass`
2. `label_tiered`
   - `Strong` 只保留給 `hp_sig + allele_sig`
   - 只有 `allele_sig` 時降為 `Weak`
3. `local_hybrid`
   - 在 `label_tiered` 基礎上用 `GlobalP + LocalBestP`
4. `permanova_strict`
   - 只讓 `ClusterPermanova` 成功且顯著的群進入 cluster support
5. `hybrid`
   - `GlobalP + (LocalBestP or ClusterPermanova)` 混合條件

### 6.2 主要結果

根據 [verification_scheme_class_counts.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_verification_scheme_adjustments_hcc1395_5khz/verification_scheme_class_counts.tsv)：

1. 現行 `current Strong`
   - `TP=6442`
   - `FP=2244`
2. 改成 `label_tiered` 後，`Strong` 大幅收斂
   - `TP=3345`
   - `FP=385`
3. 也就是有大量 `current Strong` 被降級

但這不代表可以直接刪，因為：

- 被降級的 `current_strong_allele_only`
  - `TP=5104`
  - `FP=2042`

這一整包其實 **TP 仍然遠多於 FP**，不能整包當 artifact。

### 6.3 對 F1 的實際影響

以本輪 `InterSubMod TO baseline`：

- `TP=28393`
- `FP=11807`
- `FN=11054`
- `F1=0.712971`

做規則試算後：

| rule | TP removed | FP removed | F1 delta |
| --- | --- | --- | --- |
| `current_low_vaf_high_ad_cv` | `2` | `36` | `+0.000290` |
| `current_strong_allele_only_low_cv` | `1` | `27` | `+0.000226` |
| `hybrid_weak_low_cv` | `1` | `29` | `+0.000244` |
| `tiered_weak_low` | `2` | `29` | `+0.000227` |

因此：

1. 調整 scheme **有助於語意與標記更合理**
2. 但在 `5kHz TO` 上，**還沒有超過目前最簡單有效的 `low VAF + high AlleleDelta + low CramersV`**
3. 目前最佳的實際 F1 仍是：
   - `current_low_vaf_high_ad_cv`
   - `F1 +0.000290`

## 7. `cluster_only` 語意已修正為 `cluster_plus_weak_label`

這輪另做了一個必要的分析層修正：

- 原先的 `cluster_only`
  - 會把「完全沒有 label」和「只有弱 allele label」混在一起
- 現在在 [validate_method_design.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py) 中新增：
  - `cluster_plus_weak_label`

重跑後的輸出在：

- [20260308_to_validate_method_design_refined_hcc1395_5khz](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_validate_method_design_refined_hcc1395_5khz)
- [label_cluster_agreement.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_to_validate_method_design_refined_hcc1395_5khz/label_cluster_agreement.tsv)

新的 agreement type 分布：

1. `consistent_weak_or_noise = 23104`
2. `cluster_plus_weak_label = 7428`
3. `conflict = 3227`
4. `label_upgrade = 2756`
5. `cluster_only = 1976`
6. `consistent_strong = 1471`

而本輪原先要追的可疑子集，現在更精確地落在：

- `cluster_plus_weak_label + Strong->Weak + lowVAF/highAD`

其中：

- `lowVAF/highAD`：`25 FP / 2 TP`
- `lowVAF/highAD + low CramersV`：`25 FP / 1 TP`

這比舊的 `cluster_only` 語意更準，因為它真的把「弱 allele label」單獨拉出來了。

## 8. 本輪結論

1. `HCC1395 5kHz TO` 的 top 可疑 FP 確實有共通型態：
   - 極低 `alt fraction`
   - 高 `AlleleDelta`
   - 低 `PairwiseMedianDist`
   - 低 `QUAL/GQ`
   - 且 4/5 聚集在 chr11 的 `1402 bp` 區段
2. `HP skew` 本身不夠區分 FP/TP，因為 TP 對照也可能極度偏單一 HP
3. 現行 core 的 `Strong` 定義偏鬆，因為「只有 allele sig、沒有 HP sig」也可能進到 `Strong`
4. 但把這一整包 `current Strong allele-only` 全部視為 artifact 是錯的，因為：
   - `5104 TP / 2042 FP`
5. 因此本輪最合理的結論是：
   - **調整 scheme 可提升可解釋性**
   - **真正可帶來 F1 改善的仍是特徵子集規則，不是整包類別重標**
6. 在 `5kHz TO` 上，暫時最穩的 artifact 規則仍是：
   - `low VAF + high AlleleDelta + low CramersV`
   - `36 FP / 2 TP`
   - `F1 +0.000290`

## 9. 下一步建議

1. 將後續 TO diagnostics 的主目標從舊語意：
   - `cluster_only + Strong->Weak`
   修正為：
   - `cluster_plus_weak_label + Strong->Weak`
2. 優先把同一套 diagnostics / scheme evaluator 套到 `HCC1395_DORADO TO`
   - 檢查這個 `chr11-like local block + weak allele-only label` 是否能跨平台重現
3. 若未來要動 core C++ 分類，建議方向不是直接刪除 `Strong`
   - 而是新增一個 annotation，例如：
   - `Strong_AlleleOnly`
   - `ArtifactSuspect_WeakLabel`
4. 真正要升級成流程規則前，仍要先在 `DORADO TO` 驗證：
   - `高 AlleleDelta + 低 PairwiseMedianDist + 低 alt fraction + 低 GQ/QUAL`
   是否仍富集 FP
