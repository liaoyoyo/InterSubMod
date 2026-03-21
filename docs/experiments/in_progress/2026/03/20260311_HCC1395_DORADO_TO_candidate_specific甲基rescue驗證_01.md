<!--
建立時間: 2026-03-11 00:58
目標: 完成 HCC1395_DORADO tumor-only 的 candidate-specific InterSubMod，並用與 20260311 GQ/甲基矩陣相同的方法正式驗證 TO 下是否仍是 caller-first 有效、甲基只作 support
處理範圍:
  - HCC1395_DORADO TO baseline benchmark
  - candidate-specific pool 萃取與匯出
  - LongPhase-TO phase
  - candidate-window subset BAM 抽取與 haplotag
  - candidate-specific InterSubMod
  - design validation / rescue evaluation
  - 納入 4-dataset GQ / 甲基 rescue 矩陣
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_candidate_pool_vcfs.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/evaluate_rescue_with_methylation.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_gq_methylation_rescue_matrix.py
-->

# HCC1395_DORADO TO candidate-specific 甲基 rescue 驗證

## 1. 破題結論

本輪已正式完成 `HCC1395_DORADO TO` 的 candidate-specific InterSubMod 驗證，而且結果已能回答原本卡住的核心問題：**在 TO 模式下，`GQ` 仍是最穩定的 caller-first 主訊號；甲基資訊不是負效果，也確實能救回 TP，但在相同 caller gate 下，目前仍只適合作為第二層 support，而沒有超過 `GQ only`。**

更具體地說：

1. `DORADO TO` 的 `candidate-specific` 覆蓋已達 `100%`，不再有 `0/0 overlap` 或 coverage 不足的藉口。
2. `GQ only` 在 `DORADO TO` 是正增益：
   - `gq>=10`: `40 TP / 11 FP`, `F1 +0.000540`
   - `gq>=15`: `33 TP / 9 FP`, `F1 +0.000446`
3. 甲基特徵單獨也不是負的：
   - `Quality_Score>=60`: `208 TP / 77 FP`, `F1 +0.002620`
   - `PairwiseMedianDist<=0.20`: `173 TP / 60 FP`, `F1 +0.002217`
   - `hp_assign_rate>=0.99`: `122 TP / 41 FP`, `F1 +0.001577`
4. 但在固定 `gq>=10` gate 後，甲基 support 全部都沒有再超過 `gq>=10 only`：
   - `gq>=10 + Quality_Score>=60`: `delta F1 vs gate = -0.000064`
   - `gq>=10 + PairwiseMedianDist<=0.20`: `-0.000142`
   - `gq>=10 + hp_assign_rate>=0.99`: `-0.000257`
   - `gq>=10 + Strong/Subclone`: `-0.000366`
   - `gq>=10 + agreement_positive`: `-0.000430`

因此，若把本輪結果與 [20260311_GQ與甲基rescue系統性驗證_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md) 的 `5kHz TO` 對照，**TO 下跨樣本都支持同一個高層結論：`caller-first` 穩定成立，甲基資料可作 support，但目前不是主 rescue 規則。**

---

## 2. 研究問題

本輪要正式回答 3 個問題：

1. `HCC1395_DORADO TO` 能否用與 `HCC1395 5kHz TO` 完全相同的 candidate-specific 流程，正式驗證 TO 下的甲基 rescue。
2. 在 `DORADO TO` 下，甲基資料是否真的有效救回 TP，且不會在相同標準下救回太多 FP。
3. 將 `DORADO TO` 納入同一套 `GQ / methylation rescue matrix` 後，能否正式確認：
   - TO 下仍是 `GQ` 有效
   - 甲基仍偏向第二層 support
   - `low VAF + high AlleleDelta` 仍不適合升級成 TO 主 rescue 規則

---

## 3. 輸入資料與輸出位置

### 3.1 輸入資料

1. caller final VCF：
   - [snv.vcf.gz](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/snv.vcf.gz)
2. tumor BAM：
   - [HCC1395.bam](/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam)
3. truth VCF：
   - [high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz](/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz)
4. truth BED：
   - [High-Confidence_Regions_v1.2.bed](/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed)

### 3.2 本輪 round 根目錄

- [20260311_hcc1395_dorado_to_candidate_rescue](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue)

### 3.3 最終矩陣輸出

- [20260311_gq_methylation_rescue_matrix_with_dorado_to](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to)

---

## 4. 空間策略與執行設計

### 4.1 為什麼這輪不能直接生成 full tagged BAM

本輪啟動前的實測空間為：

| 掛載點 | 可用空間 |
| --- | ---: |
| `/` | `289G` |
| `/big8_disk` | `226G` |
| `/bip8_disk` | `191G` |

原始 `DORADO tumor BAM` 的實體大小為：

| 檔案 | 體積 |
| --- | ---: |
| [HCC1395.bam](/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam) | `233G` |

若直接生成 full `tagged BAM`，最合理落點只剩 `/home` 或 `/`；但目前 `/home` 已有 [20260309_hcc1395_dorado_paired_candidate_rescue](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue) 這個 `223G` 級 round。直接再生成一顆 `200G+` full tagged BAM，風險過高。

### 4.2 本輪採用的低空間策略

本輪改用：

1. 先從 candidate pool 匯出小型 candidate VCF
2. 將 candidate SNV 擴成 `±1000 bp` 視窗
3. 用視窗 BED 從原始 tumor BAM 抽出 coordinate-sorted subset BAM
4. 只對這顆 subset BAM 做 `LongPhase-TO haplotag`
5. 再對 candidate-specific VCF 跑 InterSubMod

這個策略的關鍵意義是：

- 不改變研究問題本身，因為本輪就是要回答 candidate pool 上的甲基 rescue
- 大幅降低空間成本，避免為了 341 個 candidate 視窗生成整顆 full tagged BAM

實際結果：

| 產物 | 體積 |
| --- | ---: |
| [tumor_candidate_windows.bam](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows.bam) | `482M` |
| [tumor_candidate_windows_tagged.bam](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam) | `469M` |

本輪 round 總體積最後約為：

- [20260311_hcc1395_dorado_to_candidate_rescue](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue): `3.1G`

---

## 5. 實際執行流程與正確性確認

### 5.1 Step 1: baseline benchmark

使用 [benchmark_split_snv_vcf.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/pipeline/utils/benchmark_split_snv_vcf.sh) 對現有 `ClairS-TO final VCF` 做 baseline benchmark。

輸出：

- [variant_counts.txt](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step02_benchmark_clairs_to/variant_counts.txt)

結果：

| 方法 | PASS_TOTAL | TP | FP | FN | F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `ClairS-TO baseline` | `40437` | `28861` | `11576` | `10586` | `0.722573` |

這與先前 [20260307_LongPhaseTO_空間需求與中間產物確認_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260307_LongPhaseTO_空間需求與中間產物確認_01.md) 內的 `DORADO TO baseline` 一致，表示本輪 baseline 口徑正確。

### 5.2 Step 2: candidate pool 萃取

輸出：

1. [borderline_candidate_pool.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/extract/borderline_candidate_pool.tsv)
2. [candidate_group_summary.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/extract/candidate_group_summary.tsv)
3. [candidate_vcf_export_summary.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/export/candidate_vcf_export_summary.tsv)

重點數據：

| 類別 | 總數 | candidate_eligible |
| --- | ---: | ---: |
| `caller_lost_tp` | `1721` | `247` |
| `caller_removed_fp` | `2713328` | `94` |

這裡有兩個關鍵意義：

1. `DORADO TO` 的 raw `removed_fp` 母體仍然非常大，與 `5kHz TO` 一樣，不適合直接全域放寬 caller 閾值。
2. 但真正進入 rescue-eligible pool 的候選只剩 `341` 筆，這正是 candidate-specific 分析適合存在的原因。

### 5.3 Step 3: LongPhase-TO phase

輸出：

1. [tumor_phased.vcf](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_phased.vcf)
2. [tumor_phased_LOH.bed](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_phased_LOH.bed)
3. [longphase_to_phase.log](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/longphase_to_phase.log)

執行結果：

| 指標 | 數值 |
| --- | ---: |
| wall time | `365s` |
| purity | `1` |

### 5.4 Step 4: candidate-window subset BAM 與 haplotag

candidate 視窗 BED：

- [candidate_windows_1kb_merged.bed](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/export/candidate_windows_1kb_merged.bed)

視窗摘要：

| 指標 | 數值 |
| --- | ---: |
| merged regions | `336` |
| merged bp | `674,733` |

haplotag 輸出：

1. [tumor_candidate_windows_tagged.bam](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam)
2. [tumor_candidate_windows_tagged.bam.bai](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam.bai)
3. [longphase_to_haplotag.log](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/longphase_to_haplotag.log)

haplotag 結果：

| 指標 | 數值 |
| --- | ---: |
| total alignment | `58,766` |
| tagged alignment | `47,784` |
| subset HP assign rate | `0.813123` |

### 5.5 Step 5: phased VCF benchmark

輸出：

- [step04_benchmark_longphase_to/variant_counts.txt](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step04_benchmark_longphase_to/variant_counts.txt)

結果：

| 方法 | TP | FP | FN | F1 |
| --- | ---: | ---: | ---: | ---: |
| `LongPhase-TO phased VCF` | `28861` | `11576` | `10586` | `0.722573` |

這與 `ClairS-TO baseline` 完全相同，表示：

> 在 `HCC1395_DORADO TO` 上，`LongPhase-TO phase` 本身沒有改變 PASS call set；本輪真正要驗證的增益仍是 candidate-specific 甲基層。

### 5.6 Step 6: candidate-specific InterSubMod

輸出：

1. [TP significance_summary.csv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/intersubmod/intersubmod_tp/significance_summary.csv)
2. [FP significance_summary.csv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/intersubmod/intersubmod_fp/significance_summary.csv)
3. [label_cluster_agreement.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/design_validation/label_cluster_agreement.tsv)
4. [rescue_rule_comparison.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_rule_comparison.tsv)
5. [rescue_joined_features.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv)

candidate-specific coverage：

| 類別 | 可分析數 | coverage |
| --- | ---: | ---: |
| `caller_lost_tp` | `247 / 247` | `100%` |
| `caller_removed_fp` | `94 / 94` | `100%` |

這是本輪最重要的正確性確認之一，因為它代表：

> `DORADO TO` 的 candidate-specific run 已不再受 coverage 不足干擾，可以直接解讀規則效果。

---

## 6. 方法學補充：為什麼本輪需要 analysis-only 的 PASS override

本輪第一次直接用 [candidate_lost_tp.vcf](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/export/candidate_lost_tp.vcf) 跑 InterSubMod 時，出現：

- `Loaded: 0, Skipped: 247`

原因不是 BAM 或 phased VCF 問題，而是目前 `inter_sub_mod` 的 VCF loader 只讀 `FILTER=PASS` 的 SNV。這一點可從 [src/core/SomaticSnv.cpp](/big8_disk/liaoyoyo2001/InterSubMod/src/core/SomaticSnv.cpp) 看出，程式會用 `bcf_has_filter(..., "PASS")` 決定是否載入。

因此本輪對 [export_candidate_pool_vcfs.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_candidate_pool_vcfs.py) 補上 `--force-pass-filter`，讓 **analysis-only 的 candidate VCF** 在匯出時把 `FILTER` 欄改寫成 `PASS`，原始 caller final VCF 則完全不變。

這個調整的意義是：

1. 只影響 candidate-specific InterSubMod 的分析輸入
2. 不改 benchmark 計數、不改 caller 原始結果
3. 讓本輪能真正量到 `DORADO TO` 非 PASS borderline 候選的甲基訊號

---

## 7. 主要結果

### 7.1 caller-first：`GQ` 在 DORADO TO 仍有效

資料來源：

- [gq_threshold_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/gq_threshold_sweep.tsv)

| 規則 | TP rescued | FP reintroduced | F1 | delta F1 vs baseline |
| --- | ---: | ---: | ---: | ---: |
| `gq>=10` | `40` | `11` | `0.723113` | `+0.000540` |
| `gq>=15` | `33` | `9` | `0.723019` | `+0.000446` |
| `gq>=20` | `16` | `3` | `0.722801` | `+0.000229` |
| `qual>=10 or gq>=20` | `16` | `3` | `0.722801` | `+0.000229` |

解讀：

1. `DORADO TO` 和 `5kHz TO` 一樣，`GQ` 在 candidate pool 上是有效的 rescue 訊號。
2. 但最佳區間比 `5kHz TO` 更保守，`gq>=10` 與 `gq>=15` 都只帶來小幅增益。
3. 這支持「TO 下仍要 caller-first」，但不支持把某個單一 `GQ` 閾值直接升級成所有樣本的全域標準。

### 7.2 甲基單特徵：不是負效果，而且有明顯正增益

資料來源：

- [methylation_only_rule_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/methylation_only_rule_sweep.tsv)

| 規則 | TP rescued | FP reintroduced | F1 | delta F1 vs baseline | 意義 |
| --- | ---: | ---: | ---: | ---: | --- |
| `Quality_Score>=60` | `208` | `77` | `0.725193` | `+0.002620` | 整體綜合結構品質較高的 candidate pool 子集 |
| `PairwiseMedianDist<=0.20` | `173` | `60` | `0.724790` | `+0.002217` | 在 `DORADO TO`，較低 pairwise 反而更偏向可救回 TP |
| `hp_assign_rate>=0.99` | `122` | `41` | `0.724150` | `+0.001577` | HP 指派很完整的候選較可信 |
| `Strong/Subclone` | `88` | `30` | `0.723707` | `+0.001134` | cluster-first 主類別仍有正訊號 |
| `agreement_positive` | `59` | `27` | `0.723271` | `+0.000698` | label/cluster 交叉正訊號有用，但不算強 |
| `PairwiseMedianDist>=0.20` | `74` | `34` | `0.723447` | `+0.000875` | 高 pairwise 也不是全負，但方向弱於低 pairwise |

解讀：

1. `DORADO TO` 證明甲基特徵不是全負。
2. 但它與 `5kHz TO` 的方向不完全相同：
   - `5kHz TO` 偏向 `高 PairwiseMedianDist`
   - `DORADO TO` 偏向 `低 PairwiseMedianDist`
3. 這代表 `PairwiseMedianDist` 本身是有訊號的，但**方向具樣本依賴性**，目前不適合直接升級成跨 TO 樣本的單一硬閾值。

### 7.3 `GQ + 甲基`：在相同 caller gate 下仍未超過 `GQ only`

資料來源：

- [gq_plus_methylation_rule_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/gq_plus_methylation_rule_sweep.tsv)

固定比較 `gq>=10` gate：

| 規則 | TP rescued | FP reintroduced | delta F1 vs baseline | delta F1 vs gate |
| --- | ---: | ---: | ---: | ---: |
| `gq>=10 + Quality_Score>=60` | `36` | `11` | `+0.000476` | `-0.000064` |
| `gq>=10 + PairwiseMedianDist<=0.20` | `30` | `9` | `+0.000398` | `-0.000142` |
| `gq>=10 + hp_assign_rate>=0.99` | `20` | `4` | `+0.000284` | `-0.000257` |
| `gq>=10 + Strong/Subclone` | `12` | `2` | `+0.000174` | `-0.000366` |
| `gq>=10 + PairwiseMedianDist>=0.20` | `10` | `2` | `+0.000142` | `-0.000398` |
| `gq>=10 + agreement_positive` | `8` | `2` | `+0.000110` | `-0.000430` |

這張表的意義很關鍵：

1. 這些規則相對 baseline 仍然是正的。
2. 但相對同一個 `gq>=10` caller gate，全部都變成負的。
3. 也就是說，甲基訊號在 `DORADO TO` 不是沒用，而是**一旦 caller 已經先用 `GQ` 選出邊緣候選，甲基 support 目前還不能進一步把 F1 再推高。**

這與 `5kHz TO` 的高層結論一致。

### 7.4 label-first / cluster-first 的交叉訊號在 DORADO TO 仍然存在，但不夠強

資料來源：

- [label_cluster_agreement.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/design_validation/label_cluster_agreement.tsv)

`caller_lost_tp` 中：

| agreement_type | count |
| --- | ---: |
| `consistent_weak_or_noise` | `108` |
| `cluster_plus_weak_label` | `51` |
| `consistent_strong` | `30` |
| `label_upgrade` | `29` |
| `conflict` | `24` |
| `cluster_only` | `5` |

`caller_removed_fp` 中：

| agreement_type | count |
| --- | ---: |
| `consistent_weak_or_noise` | `42` |
| `label_upgrade` | `15` |
| `cluster_plus_weak_label` | `14` |
| `consistent_strong` | `12` |
| `conflict` | `9` |
| `cluster_only` | `2` |

解讀：

1. `agreement_positive` 在 `DORADO TO` 不是空的。
2. 但它不像 `5kHz TO` 那樣特別乾淨，因為 `removed_fp` 裡也有 `27` 筆正向 agreement 訊號。
3. 這就是為什麼 `agreement_positive` 在 `DORADO TO` 仍為正增益，但只到 `+0.000698`，且在 `gq>=10` gate 後不再加分。

### 7.5 artifact 規則：仍不適合拿來做 TO 主 rescue

資料來源：

- [artifact_role_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/artifact_role_comparison.tsv)

| 規則 | 角色 | 結果 | 解讀 |
| --- | --- | --- | --- |
| `low VAF + high AlleleDelta` | `early_selection` | `0 TP / 2 FP` | 幾乎沒有 TP rescue 能力 |
| `low VAF + high AlleleDelta (+ low CramersV)` | `rescue_veto @ gq>=10/15/20` | `delta F1 vs gate = 0` | 對 DORADO TO 幾乎沒有作用 |
| `combined_artifact_veto` | `rescue_veto @ gq>=10` | `delta F1 vs gate = -0.000094` | 仍是負向 |

結論：

> `low VAF + high AlleleDelta` 這條線在 `DORADO TO` 仍較適合當後段 artifact triage，不適合提前升級成主 rescue 規則。

---

## 8. 與 5kHz TO 的直接對照

資料來源：

- [20260311_gq_methylation_rescue_matrix_with_dorado_to](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to)

| 類型 | `HCC1395 5kHz TO` | `HCC1395 DORADO TO` | 解讀 |
| --- | --- | --- | --- |
| `gq>=10 only` | `499 TP / 119 FP`, `+0.006943` | `40 TP / 11 FP`, `+0.000540` | 兩者都支持 caller-first；5kHz 增量較大 |
| `Quality_Score>=60 only` | `556 TP / 162 FP`, `+0.007466` | `208 TP / 77 FP`, `+0.002620` | 兩者都為正 |
| `PairwiseMedianDist>=0.20 only` | `408 TP / 132 FP`, `+0.005374` | `74 TP / 34 FP`, `+0.000875` | 高 pairwise 只在 5kHz TO 明顯較強 |
| `PairwiseMedianDist<=0.20 only` | `268 TP / 81 FP`, `+0.003589` | `173 TP / 60 FP`, `+0.002217` | DORADO TO 反而更偏向低 pairwise |
| `gq>=10 + Quality>=60` vs gate | `-0.001391` | `-0.000064` | 兩者都未超過 `gq>=10 only` |
| `gq>=10 + agreement_positive` vs gate | `-0.004780` | `-0.000430` | 兩者都未超過 `gq>=10 only` |

這張對照表支持兩個層次的結論：

1. **高層結論已跨 TO 樣本成立**：
   - `GQ` 有效
   - 甲基單特徵不是負的
   - 但甲基在 caller gate 後仍只作 support，不是主規則
2. **低層特徵方向尚未跨樣本穩定**：
   - `PairwiseMedianDist` 在 `5kHz TO` 與 `DORADO TO` 的最佳方向不同
   - 因此暫時不能把某個 `pairwise` 閾值升級成所有 TO 樣本的固定規則

---

## 9. 結論

本輪正式完成後，現在可把口徑更新為：

1. `HCC1395_DORADO TO candidate-specific InterSubMod` 已完成，`TO` 下是否仍是 `GQ` 有效、甲基只作 support，現在已可直接回答。
2. `TO` 模式下，跨 `5kHz` 與 `DORADO` 兩個樣本都支持：
   - `caller-first` 穩定成立
   - 甲基資料不是無效，也不是全負
   - 但在相同 caller gate 下，甲基 support 目前尚未超過最佳 `GQ only`
3. `low VAF + high AlleleDelta (+ low CramersV)` 在 `DORADO TO` 仍不適合作為主 rescue 規則。
4. `PairwiseMedianDist` 是有訊號的，但方向具樣本依賴；目前不能直接升級成跨 TO 全域規則。
5. 若把這輪與 `DORADO paired` 一起看，現在最穩定的跨樣本結論是：
   - **paired 與 TO 都是 caller-first 穩定**
   - **甲基層可作 support / ranking / annotation**
   - **artifact-veto 保持獨立 triage 支線**

---

## 10. 目前尚缺與下一步

1. `DORADO TO` 雖然已補齊，但目前仍是 candidate-specific 分析；若之後要做完整 TO round，仍可再補 full kept-set 的 InterSubMod 與 dashboard。
2. `PairwiseMedianDist` 的方向在 `5kHz TO` 與 `DORADO TO` 相反，下一步應優先看：
   - [feature_interval_enrichment.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/feature_interval_enrichment.tsv)
   - [pairwise_feature_grid.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_gq_methylation_rescue_matrix_with_dorado_to/pairwise_feature_grid.tsv)
3. `gq>=10` 後的最佳 support 在 `DORADO TO` 已接近持平但仍為負；若要再往前推，最有價值的方向是：
   - `Quality_Score + hp_assign_rate`
   - `Pairwise <= 某區間 + Quality_Score`
   - 而不是直接回頭用 artifact veto

