<!--
建立時間: 2026-03-08 23:20
目標: 將 HCC1395 5kHz TO borderline TP rescue 的關鍵特徵、規則定義、數據結果與推論意義整理成可直接查閱的證據鏈文件
處理範圍:
  - HCC1395 5kHz TO candidate-specific rescue
  - caller-only / methylation-support / agreement-support / artifact-veto
  - 特徵意義、判斷標準與結果比較
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/evaluate_rescue_with_methylation.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/
-->

# HCC1395 5kHz TO borderline rescue 特徵證據鏈整理

## 1. 文件目的

這份文件專門整理 `HCC1395 5kHz tumor-only` 的 borderline TP rescue 問題，回答 4 件事：

1. 本輪到底用了哪些資料與過程。
2. 目前觀察到哪些特徵有用。
3. 每個特徵的意義、判斷標準與風險是什麼。
4. 根據現有證據，這些特徵在流程上應如何定位。

這份文件不是單純重貼結果，而是把「規則定義」和「數據結果」連成可追溯的證據鏈。

## 2. 問題定義與範圍

### 2.1 本次問題

要回答的是：

**在 `ClairS-TO final VCF` 中，那些位於信心度邊緣、被 downstream 流程丟掉的真 TP 候選，是否可以藉由甲基資料補強判定，提升 recall，同時不要重新引入太多 FP。**

### 2.2 固定範圍

1. 樣本：`HCC1395 5kHz tumor-only`
2. Caller 輸入：`ClairS-TO final snv.vcf.gz`
3. Caller 路線：`pileup`
4. 本輪不使用：
   - pileup 中間 candidate / tensor
   - full model VCF
5. 本輪只處理：
   - `caller_lost_tp`
   - `caller_removed_fp`

### 2.3 最重要的限制

這份結論只對：

1. `HCC1395 5kHz`
2. `tumor-only`
3. `ClairS-TO pileup final VCF`

直接成立。

不能直接外推到：

1. `DORADO`
2. `paired`
3. `ClairS-TO full model`

## 3. 資料與證據鏈

### 3.1 原始資料

1. Caller VCF：[snv.vcf.gz](/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz)
2. Candidate pool：[borderline_candidate_pool.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_to/extract/borderline_candidate_pool.tsv)
3. TO tagged BAM：[tumor_tagged.bam](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam)

### 3.2 中間處理

1. 從 candidate pool 匯出：
   - [candidate_lost_tp.vcf](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_lost_tp.vcf)
   - [candidate_removed_fp.vcf](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_removed_fp.vcf)
2. 匯出摘要：[candidate_vcf_export_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_vcf_export_summary.tsv)
3. Candidate-specific InterSubMod：
   - [TP significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_tp/significance_summary.csv)
   - [FP significance_summary.csv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_fp/significance_summary.csv)
4. label/cluster 關係表：
   - [label_cluster_agreement.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/design_validation/label_cluster_agreement.tsv)
5. 最終 rescue 規則比較：
   - [rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_rule_comparison.tsv)

### 3.3 證據鏈順序

本輪證據鏈是：

1. `ClairS-TO final VCF`
2. `extract_borderline_rescue_candidates.py` 先定義 `caller_lost_tp / caller_removed_fp`
3. `export_candidate_pool_vcfs.py` 匯出 candidate-specific VCF
4. `02_intersubmod.sh` + `inter_sub_mod` 對 candidate-specific pool 重跑 InterSubMod
5. `validate_method_design.py` 產出 `agreement_type`
6. `evaluate_rescue_with_methylation.py` 以完整 lost/removed 母數比較各種 rescue 規則

重點是：

**這次不是拿 baseline kept-set 的 InterSubMod summary 硬 join，而是真正對 candidate pool 單獨重跑。**

## 4. 資料量與可分析性

### 4.1 Candidate export

根據 [candidate_vcf_export_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/export/candidate_vcf_export_summary.tsv)：

| 候選類型 | target | exported | missing |
| --- | ---: | ---: | ---: |
| `caller_lost_tp` | 773 | 773 | 0 |
| `caller_removed_fp` | 298 | 298 | 0 |

### 4.2 Candidate-specific InterSubMod 覆蓋率

| 候選類型 | total | 有 InterSubMod summary | 覆蓋率 |
| --- | ---: | ---: | ---: |
| `caller_lost_tp` | 773 | 675 | 87.32% |
| `caller_removed_fp` | 298 | 213 | 71.48% |

這表示：

1. 這輪 `TO methylation rescue` 已經是正式可解讀結果，不再是舊版 `0/0 overlap`
2. 但仍有一部分候選沒有進入可分析 InterSubMod summary，解讀時要記得這不是 100% 全覆蓋

## 5. 主要特徵、意義與判斷標準

### 5.1 Caller-side feature

| 特徵 | 來源 | 判斷標準 | 直觀意義 | 主要風險 |
| --- | --- | --- | --- | --- |
| `GQ` | ClairS-TO final VCF | `gq>=15`、`gq>=20` | caller 對 genotype/variant 的信心 | 高 GQ 不保證一定是真 somatic |
| `QUAL` | ClairS-TO final VCF | `qual>=10` | caller 綜合信心 | 不同 caller / mode 間不一定可直接比 |
| `candidate_eligible` | candidate pool | 已通過 rescue 候選條件 | 先限制在「邊緣但仍值得看」的集合 | 決定候選池本身會影響後續結論 |

### 5.2 Methylation support feature

| 特徵 | 來源 | 判斷標準 | 直觀意義 | 主要風險 |
| --- | --- | --- | --- | --- |
| `PairwiseMedianDist` | InterSubMod `significance_summary.csv` | `PairwiseMedianDist >= 0.20` 且 `gq>=10` | read-level 甲基距離夠高，表示 region 內存在較明顯甲基結構差異 | 距離高不一定代表真實 somatic，也可能來自平台雜訊或局部 block artifact |
| `VerificationClass` | InterSubMod `significance_summary.csv` | `VerificationClass in {Strong, Subclone}` 且 `gq>=10` | InterSubMod 將該 region 判為較強的正向甲基訊號 | 單看 `Strong/Subclone` 會混入不少 FP |

### 5.3 Label / Cluster 交叉特徵

`agreement_positive` 不是原始欄位，而是由 [validate_method_design.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py) 根據 `cluster_class` 和 `label_class` 推得：

| `agreement_type` | 定義 | 解讀 |
| --- | --- | --- |
| `consistent_strong` | `cluster=Strong` 且 `label=Strong` | cluster-first 與 label-first 同時支持 strong |
| `consistent_subclone` | `cluster=Subclone` 且 `label=Subclone` | 兩套方法都支持 subclone |
| `label_upgrade` | `cluster in {Weak, Noise}` 且 `label in {Strong, Subclone}` | label-first 提供比 cluster-first 更強的正向支持 |

因此本輪的 `agreement_positive` 規則，實際標準是：

1. `candidate_eligible`
2. `gq >= 10`
3. `agreement_type in {label_upgrade, consistent_strong, consistent_subclone}`

它的意義是：

**不是只看甲基強不強，而是要求 cluster-first 與 label-first 之間至少出現合理的正向一致或升級。**

### 5.4 Artifact veto feature

在 [evaluate_rescue_with_methylation.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/evaluate_rescue_with_methylation.py) 中，artifact 類規則主要有兩層：

1. `low VAF + high AlleleDelta + low CramersV`
   - `af < 0.24`
   - `AlleleDelta > 0.25`
   - `CramersV < 0.05`
2. `combined_artifact_veto`
   - 命中上面規則
   - 或 `cluster=Strong` 但 `label=Weak`
   - 或 `PairwiseMedianDist < 0.12` 且 `AlleleDelta > 0.15`
   - 或 `hp_assign_rate < 0.50` 且 `cluster=Strong`

這類特徵的意義是：

**更像是「可疑 artifact」的訊號，不是正向救 TP 的訊號。**

## 6. 關鍵規則結果

以下表格是目前最重要的規則比較：

| 規則 | 類型 | TP rescued | FP reintroduced | 觸發集合 precision | F1 delta |
| --- | --- | ---: | ---: | ---: | ---: |
| `qual>=10 or gq>=20` | caller-only | 491 | 118 | 80.62% | +0.006824 |
| `gq>=15` | caller-only | 365 | 71 | 83.72% | +0.005233 |
| `gq>=10 + PairwiseMedianDist>=0.20` | caller+methylation support | 300 | 68 | 81.52% | +0.004219 |
| `gq>=10 + Strong/Subclone` | caller+methylation support | 149 | 30 | 83.24% | +0.002134 |
| `gq>=10 + agreement_positive` | caller+label/cluster support | 148 | 25 | 85.55% | +0.002163 |
| `gq>=15 + combined_artifact_veto` | caller+methylation veto | 308 | 65 | 82.57% | +0.004374 |

資料來源：

1. [rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_rule_comparison.tsv)
2. [label_cluster_agreement.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/design_validation/label_cluster_agreement.tsv)

## 7. 逐條推論與意義

### 7.1 `PairwiseMedianDist >= 0.20` 為什麼值得記住

證據：

1. `300 TP / 68 FP`
2. `F1 +0.004219`
3. 觸發集合 precision `81.52%`

推論：

1. 單看 caller 邊緣信號之外，read-level 甲基距離確實能抓到一批「不像隨機 noise」的候選 TP
2. 這批規則雖然沒有超過最佳 caller-only，但增益已足夠說明甲基資料不是無效訊號
3. 它更像是「中等覆蓋、可實用的 support feature」

限制：

1. 高距離不一定等於真 TP
2. 它仍會帶回 `68` 個 FP
3. 因此不適合單獨升級成 hard keep

### 7.2 `agreement_positive` 為什麼比 `Strong/Subclone` 更值得注意

證據：

1. `Strong/Subclone`: `149 TP / 30 FP`
2. `agreement_positive`: `148 TP / 25 FP`
3. `agreement_positive` 幾乎不損 TP，但少帶回 `5` 個 FP

推論：

1. 單純把 `VerificationClass=Strong/Subclone` 當 support，還不夠精細
2. 加入 `label-first / cluster-first` 的交叉關係後，可以稍微去掉一部分不夠一致的 FP
3. 這表示 `agreement_type` 是一個有方法學意義的 refinement feature，不只是描述欄位

更重要的含義是：

**對 borderline rescue 而言，甲基訊號本身不是只有強弱，訊號之間是否一致也有資訊量。**

### 7.3 為什麼 `artifact veto` 在這裡表現不好

證據：

1. `lowVAF/highAlleleDelta/lowCramersV` 對 `gq>=15` 這批 rescue pool 幾乎沒有額外效果
2. `combined_artifact_veto` 雖少 `6 FP`，卻同時少 `57 TP`

推論：

1. 先前很有效的 TO artifact 特徵，主要是拿來做 FP triage
2. 但在這個問題裡，我們要救的是 `caller_lost_tp`
3. 這兩個集合不是同一種分布

因此：

**artifact-veto 和 TP-rescue 不應混成同一條規則。**

### 7.4 為什麼最佳整體仍是 caller-only

證據：

1. `qual>=10 or gq>=20` 仍有全場最佳 `F1 +0.006824`
2. `gq>=15` 也仍高於所有甲基 support 規則

推論：

1. caller 數值仍是主訊號
2. 甲基資料現階段提供的是輔助支持，不是主要判決來源
3. 這與先前 paired 與 LongPhase-S rescue 的觀察一致：caller 仍是第一層，甲基更適合第二層

## 8. label/cluster 訊號本身的證據

在可分析的 candidate-specific pool 中：

### 8.1 `lost_tp` 主要 agreement 類型

| agreement_type | 數量 |
| --- | ---: |
| `consistent_weak_or_noise` | 260 |
| `label_upgrade` | 126 |
| `cluster_plus_weak_label` | 106 |
| `consistent_strong` | 86 |
| `conflict` | 83 |

### 8.2 `removed_fp` 主要 agreement 類型

| agreement_type | 數量 |
| --- | ---: |
| `consistent_weak_or_noise` | 98 |
| `consistent_strong` | 32 |
| `label_upgrade` | 30 |
| `conflict` | 24 |
| `cluster_plus_weak_label` | 23 |

推論：

1. `agreement_positive` 在 TP 和 FP 中都存在，不能單靠 agreement 類別直接做 hard keep
2. 但配合 `gq>=10` 之後，`agreement_positive` 的 precision 會提升到 `85.55%`
3. 所以 agreement 比較適合做第二層 support，而不是第一層獨立決策

## 9. 現階段最合理的流程定位

根據這輪證據，最合理的分工是：

1. **第一層：caller-first**
   - 先用 `GQ / QUAL / candidate_eligible` 定義可救回的邊緣集合
2. **第二層：methylation-support**
   - `PairwiseMedianDist >= 0.20`
   - `agreement_positive`
3. **獨立旁路：artifact triage**
   - `low VAF + high AlleleDelta + low CramersV`
   - `combined_artifact_veto`

不要做的事：

1. 不要把 artifact-veto 直接拿去做 TP rescue 主規則
2. 不要把 `Strong/Subclone` 直接當 hard keep
3. 不要把這輪 `5kHz pileup` 結論直接當成跨平台規則

## 10. 最終總結

本輪可以明確記住的結論是：

1. `TO` 下甲基資料**不是無效**，它確實能幫 borderline TP rescue
2. 目前最有用的 support 特徵是：
   - `PairwiseMedianDist >= 0.20`
   - `agreement_positive`
3. `agreement_positive` 比單純 `Strong/Subclone` 更乾淨，證明 label-first / cluster-first 的交叉訊號有方法學價值
4. 但最佳整體 rescue 仍是 caller-only，代表甲基訊號目前應定位為第二層 support
5. 先前有效的 artifact 特徵，在這題上不應直接用來救 TP；它比較適合獨立的 FP triage

## 11. 建議下一步

1. 用同一份結構整理 `HCC1395_DORADO TO`
2. 對 `agreement_positive` 與 `PairwiseMedianDist>=0.20` rescue 到的 top TP / FP 補 diagnostics
3. 若未來要進流程，先做 annotation / ranking，不先做 hard keep
