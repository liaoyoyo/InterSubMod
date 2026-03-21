<!--
建立時間: 2026-03-09 18:40
目標: 補齊 HCC1395_DORADO paired 的 candidate-specific InterSubMod，驗證 5kHz TO 已觀察到的甲基 rescue 規則是否可跨樣本成立
處理範圍: HCC1395_DORADO paired pure sample、ClairS final VCF 邊緣 candidate pool、LongPhase-S tagged BAM 重建、candidate-specific InterSubMod、rescue 規則重算
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣TP_rescue與甲基輔助評估_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/extract_borderline_rescue_candidates.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_candidate_pool_vcfs.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/evaluate_rescue_with_methylation.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py
-->

# HCC1395_DORADO paired candidate-specific 甲基 rescue 驗證

## 1. 研究問題

本輪要回答的是：

1. `HCC1395 5kHz TO` 已成立的甲基 rescue support 規則，是否也能在 `HCC1395_DORADO paired` 成立。
2. 同一套規則在 `DORADO paired` 下，是否真的能救回 TP，且不會在相同標準下救回太多 FP。
3. 若不成立，失敗主因是 coverage 不足、模式差異，還是規則本身不具跨樣本穩定性。

本輪固定比較 3 條甲基 support 規則：

1. `gq>=10 + PairwiseMedianDist>=0.20`
2. `gq>=10 + agreement_positive`
3. `gq>=10 + Strong/Subclone`

caller-first baseline 則維持：

1. `candidate_gq_ge_15`
2. `candidate_gq_ge_20`

## 2. 資料與輸入定義

### 2.1 caller 與 baseline 定義

- caller final VCF: [output.vcf.gz](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz)
- truth VCF: [high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz](/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz)
- truth BED: [High-Confidence_Regions_v1.2.bed](/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed)
- baseline kept TP: [filtered_snv_tp.vcf.gz](/bip8_disk/liaoyoyo2001/InterSubMod_out/output/s-pure/HCC1395_DORADO/20260211/longphase_s/filtered_snv_tp.vcf.gz)
- baseline kept FP: [filtered_snv_fp.vcf.gz](/bip8_disk/liaoyoyo2001/InterSubMod_out/output/s-pure/HCC1395_DORADO/20260211/longphase_s/filtered_snv_fp.vcf.gz)

baseline count：

- baseline TP = `29889`
- baseline FP = `240`
- truth total = `39447`
- baseline F1 = `0.859176`

### 2.2 candidate pool

本輪重新產生 candidate pool，不沿用舊的 `0/0 overlap` 推估。

- candidate summary: [candidate_group_summary.md](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/extract/candidate_group_summary.md)
- candidate pool: [borderline_candidate_pool.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/extract/borderline_candidate_pool.tsv)
- candidate export summary: [candidate_vcf_export_summary.md](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/export/candidate_vcf_export_summary.md)

結果：

- `caller_lost_tp = 1489`
- `candidate_eligible lost_tp = 1122`
- `caller_removed_fp = 2533`
- `candidate_eligible removed_fp = 1658`
- candidate VCF export missing 全部為 `0`

## 3. 實際執行過程

### 3.1 空間與 tagged BAM 重建

原本 `/big8_disk` 與 `/bip8_disk` 可用空間不足以安全重建 `DORADO paired tagged BAM`，因此本輪改寫到 `/home` scratch。

- 重建目錄: [rebuild/longphase_s](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/rebuild/longphase_s)
- LongPhase-S log: [longphase_s.log](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/rebuild/longphase_s/longphase_s.log)
- tagged BAM: [HCC1395_DORADO_tagged.bam](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/rebuild/longphase_s/HCC1395_DORADO_tagged.bam)
- tagged BAM index: [HCC1395_DORADO_tagged.bam.bai](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/rebuild/longphase_s/HCC1395_DORADO_tagged.bam.bai)

重建結果：

- BAM size = `223G`
- LongPhase-S wall time = `18m39s`
- total tagged alignments = `18,794,025 / 40,033,094`
- tagged fraction = `0.469462`
- HP1 = `8,688,357`
- HP2 = `8,951,897`
- HP1-1 = `483,755`
- HP2-1 = `493,850`
- HP3 = `176,166`

註記：

- `haplotag_qc.py` 在 `bam.bai` 產出後被手動中止，因為它只是在完整重掃 `223G` BAM。
- 本輪後續 `validate_method_design` 主要依賴 candidate-specific region 內的 `reads.tsv`，不依賴 sample-level `haplotag_qc.tsv`。

### 3.2 candidate-specific InterSubMod

- TP summary: [significance_summary.csv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/intersubmod/intersubmod_tp/significance_summary.csv)
- FP summary: [significance_summary.csv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/intersubmod/intersubmod_fp/significance_summary.csv)
- design validation: [label_cluster_agreement.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/design_validation/label_cluster_agreement.tsv)
- rescue join: [rescue_joined_features.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_joined_features.tsv)
- rescue rules: [rescue_rule_comparison.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_rule_comparison.tsv)
- strategy summary: [strategy_summary.md](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/strategy_compare/strategy_summary.md)

candidate-specific 可分析率：

- `caller_lost_tp`: `93 / 1122 = 8.29%`
- `caller_removed_fp`: `326 / 1658 = 19.66%`

這表示 `DORADO paired` 的 candidate-specific 甲基 coverage 明顯低於 `5kHz TO`，是解讀結果時必須先記住的第一層限制。

## 4. 規則結果

### 4.1 caller-only baseline

| 規則 | TP rescued | FP reintroduced | F1 delta |
| --- | ---: | ---: | ---: |
| `candidate_gq_ge_15` | 97 | 88 | +0.000502 |
| `candidate_gq_ge_20` | 51 | 7 | +0.000749 |
| `caller_any_gq_ge_20` | 53 | 7 | +0.000782 |

解讀：

- `caller-first` 在 `DORADO paired` 仍穩定成立。
- 若看 safety 排序，`gq>=15` 是最佳 safe rule。
- 若只看 `delta F1`，`gq>=20` 系列更好。

### 4.2 甲基 support 規則

| 規則 | TP rescued | FP reintroduced | F1 delta | 判讀 |
| --- | ---: | ---: | ---: | --- |
| `gq>=10 + PairwiseMedianDist>=0.20` | 10 | 36 | -0.000280 | 無效 |
| `gq>=10 + agreement_positive` | 36 | 72 | -0.000298 | 無效 |
| `gq>=10 + Strong/Subclone` | 25 | 38 | -0.000059 | 接近持平但仍為負 |

直接結論：

1. `DORADO paired` 不支持把這 3 條甲基 support 規則當成 TP rescue 正規則。
2. 這不是單純「提升小」，而是 **同一標準下救回的 FP 明顯多於 TP**。
3. 因此目前不能把 `5kHz TO` 的甲基 rescue 結論升級成「跨樣本成立」。

## 5. 為什麼 DORADO paired 不成立

### 5.1 可分析 coverage 低

這是最先要看的限制：

- `lost TP` 只覆蓋 `8.29%`
- `removed FP` 只覆蓋 `19.66%`

代表大部分 `DORADO paired` 邊緣 caller 候選，最後根本沒有足夠的甲基 summary 可以判。

### 5.2 在可分析子集內，正向甲基訊號也不偏向 TP

在已分析子集內：

`caller_lost_tp` 的 `VerificationClass`：

- `Noise = 36`
- `Weak = 28`
- `Strong = 22`
- `Subclone = 7`

`caller_removed_fp` 的 `VerificationClass`：

- `Noise = 232`
- `Weak = 38`
- `Strong = 41`
- `Subclone = 15`

也就是說：

- `Strong/Subclone` 在 TP pool 內不是沒有。
- 但在 FP pool 內也不少，而且絕對數還更高。

這正是 `gq>=10 + Strong/Subclone` 最後變成 `25 TP / 38 FP` 的根本原因。

### 5.3 agreement_positive 也沒有變乾淨

`agreement_positive = label_upgrade + consistent_strong + consistent_subclone`

在已分析子集內：

- `lost TP`：
  - `label_upgrade = 24`
  - `consistent_strong = 13`
  - `consistent_subclone = 0`
- `removed FP`：
  - `label_upgrade = 59`
  - `consistent_strong = 36`
  - `consistent_subclone = 4`

因此：

- `agreement_positive` 在 `DORADO paired` 並沒有像 `5kHz TO` 一樣變成比較乾淨的 support。
- 它在這個樣本/模式下反而更容易把 removed FP 一起帶回來。

### 5.4 PairwiseMedianDist>=0.20 也偏向 FP

整個 candidate-specific 可分析子集裡：

- `Pairwise>=0.20` 命中 TP = `13`
- `Pairwise>=0.20` 命中 FP = `58`

套上 `gq>=10` 之後仍然是：

- `10 TP / 36 FP`

這代表在 `DORADO paired` 中，單靠較高的 pairwise distance 並不能安全代表「更像真 TP」。

## 6. 與 5kHz TO 的直接對照

| 模式 | 規則 | TP rescued | FP reintroduced | F1 delta |
| --- | --- | ---: | ---: | ---: |
| `HCC1395 5kHz TO` | `gq>=10 + Pairwise>=0.20` | 300 | 68 | +0.004219 |
| `HCC1395_DORADO paired` | `gq>=10 + Pairwise>=0.20` | 10 | 36 | -0.000280 |
| `HCC1395 5kHz TO` | `gq>=10 + agreement_positive` | 148 | 25 | +0.002163 |
| `HCC1395_DORADO paired` | `gq>=10 + agreement_positive` | 36 | 72 | -0.000298 |
| `HCC1395 5kHz TO` | `gq>=10 + Strong/Subclone` | 149 | 30 | +0.002134 |
| `HCC1395_DORADO paired` | `gq>=10 + Strong/Subclone` | 25 | 38 | -0.000059 |

這張表已經足夠支持目前最保守、也最正確的判讀：

1. `5kHz TO` 的甲基 rescue 是真實存在的。
2. 但它 **不是穩定跨樣本 / 跨模式規則**。
3. 真正穩定成立的仍是 `caller-first`，不是 `methylation-first`。

## 7. 本輪結論

### 7.1 已成立的結論

1. `HCC1395_DORADO paired` 的 candidate-specific InterSubMod 已正式補齊，不再是舊版缺 summary 狀態。
2. `caller-only` 在 `DORADO paired` 穩定成立，最佳 safe rule 仍是 `gq>=15`。
3. `gq>=10 + Pairwise>=0.20`、`gq>=10 + agreement_positive`、`gq>=10 + Strong/Subclone` 在 `DORADO paired` 全部不成立。

### 7.2 目前最合理的總結口徑

根據目前證據，最合理的分工仍然是：

1. 第一層：`caller-first`
2. 第二層：`methylation-support`
3. 但「甲基 support 可穩定提升 TP rescue」**目前只在 `HCC1395 5kHz TO` 成立**
4. `DORADO paired` 的新證據反而支持：**甲基 rescue 尚未穩定跨樣本成立**

### 7.3 對研究方向的意義

這輪結果不代表甲基資訊無用，而是代表：

1. 甲基資訊的效益高度依賴樣本、平台、caller 模式與 candidate coverage。
2. 若要宣稱「穩定成立」，至少還需要 `HCC1395_DORADO TO` 的同型 candidate-specific 驗證。
3. 在此之前，不能把 `Pairwise>=0.20` 或 `agreement_positive` 升級成跨樣本 rescue 規則。

## 8. 下一步

1. 若要直接檢驗「TO 下是否可跨平台成立」，下一步應做 `HCC1395_DORADO TO candidate-specific InterSubMod`。
2. 但這需要新的 `LongPhase-TO tagged BAM`，而目前本機可寫大空間只剩 root filesystem，且本輪已被 `223G` paired tagged BAM 佔用一大塊。
3. 因此下一輪若要做 `DORADO TO`，需先明確規劃 scratch 空間與中間產物保留策略。

