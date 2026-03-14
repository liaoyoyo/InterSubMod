<!--
建立時間: 2026-03-11 03:10
目標: 完成 phase 2：paired raw pileup vs full model 直接 benchmark 對照，並補 finer feature interval / orthogonality 驗證
處理範圍:
  - HCC1395 5kHz paired raw ClairS pileup/full/merged benchmark
  - HCC1395_DORADO paired raw ClairS pileup/full/merged benchmark
  - HCC1395 5kHz / HCC1395_DORADO × paired / TO 的 finer feature threshold / interval / orthogonality
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase2_paired_model_feature_analysis.py
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/phase2_summary.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md
-->

# Phase 2：paired raw pileup/full 與 finer feature interval / orthogonality 分析

## 1. 破題結論

這輪 phase 2 有兩個核心結論。

第一，`paired raw pileup vs full model` 的直接 benchmark 已經清楚說明：在 `HCC1395 5kHz paired` 與 `HCC1395_DORADO paired` 中，raw `pileup_filter/full_alignment_filter` 雖然都比 merged output 有更高 recall，但會同時引入大量 FP，導致 F1 低於 merged output，更遠低於 `LongPhase-S` 與 `InterSubMod`。因此 paired 場景下，raw pileup/full 的主要價值是說明「caller 原始召回空間確實存在」，不是可以直接取代 merge 後最終 call set。

第二，finer feature 分析支持先前的高層分工沒有改變：`caller-first` 仍是主規則，甲基特徵在 TO 下有明顯 support 訊號，但大多屬於區間性、資料集依賴（dataset-dependent），而不是全域 hard rule。`PairwiseMedianDist` 的方向差異在 `5kHz TO` 與 `DORADO TO` 仍然成立；paired 象限的甲基正交增益依然很弱，主因仍是 candidate-specific coverage ceiling。

---

## 2. 本輪輸出

### 2.1 執行腳本

- [run_phase2_paired_model_feature_analysis.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase2_paired_model_feature_analysis.py)

### 2.2 主輸出目錄

- [phase 2 output 目錄](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis)
- [phase2_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/phase2_summary.md)

### 2.3 主表

- [paired_raw_model_benchmark.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/paired_raw_model_benchmark.tsv)
- [fine_feature_threshold_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_threshold_sweep.tsv)
- [fine_feature_interval_detail.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_interval_detail.tsv)
- [fine_feature_interval_top_bins.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_interval_top_bins.tsv)
- [feature_orthogonality_detail.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/feature_orthogonality_detail.tsv)
- [feature_orthogonality_top_pairs.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/feature_orthogonality_top_pairs.tsv)

---

## 3. 方法與資料來源

## 3.1 paired raw pileup vs full model benchmark

paired raw model 的比較對象固定為：

- `pileup_filter.vcf`
- `full_alignment_filter.vcf`
- `merged output.vcf.gz`

資料來源：

- `HCC1395 5kHz paired`
  - [pileup_filter.vcf](/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/tmp/vcf_output/pileup_filter.vcf)
  - [full_alignment_filter.vcf](/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/tmp/vcf_output/full_alignment_filter.vcf)
  - [output.vcf.gz](/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz)
- `HCC1395_DORADO paired`
  - [pileup_filter.vcf](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/tmp/vcf_output/pileup_filter.vcf)
  - [full_alignment_filter.vcf](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/tmp/vcf_output/full_alignment_filter.vcf)
  - [output.vcf.gz](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz)

truth set 固定為：

- [SEQC2 SNV truth VCF](/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz)
- [SEQC2 HC BED](/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed)

benchmark 口徑：

1. 先把 raw VCF 清成只保留前 8 欄的合法 SNV VCF。
2. 再用 [benchmark_split_snv_vcf.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/pipeline/utils/benchmark_split_snv_vcf.sh) 做 `PASS + SNV + HC regions` benchmark。
3. 最後與既有的 `ClairS / LongPhase-S / InterSubMod` benchmark 並排。

這裡特別要記住：paired raw VCF 的 sample FORMAT 欄位中存在 `GQ='.'` 之類的格式問題，因此這輪是先轉成「僅前 8 欄」的合法 SNV VCF 再 benchmark。這不改變位點集合，但避免 `bcftools` 在 FORMAT parsing 上直接失敗。

## 3.2 finer feature interval / orthogonality

這輪沿用既有 candidate-specific rescue join，分析四個 dataset：

- `HCC1395 5kHz TO`
- `HCC1395 5kHz paired`
- `HCC1395_DORADO paired`
- `HCC1395_DORADO TO`

join 來源維持與上一輪相同的 `rescue_joined_features.tsv`，只是在本輪新增：

1. 更細的 threshold sweep
2. 更細的 numeric interval bins
3. feature union / overlap / orthogonality

重要限制：

- `5kHz TO` 與 `DORADO TO` 的 read-level snapshot scope 不同，這輪**沒有**重新把 read-level 絕對值拿來硬比。
- paired dataset 的甲基 coverage 仍偏低，因此 paired 的 interval / orthogonality 解讀要更保守。
- `delta_f1_vs_best_single` 的意義是：某兩個特徵做聯集（union）後，是否比這兩者中較好的單獨規則再更好；這個值為正，才比較接近真正的「正交補強」。

---

## 4. paired raw pileup vs full model benchmark

| Dataset | Source | TP | FP | FN | Precision | Recall | F1 | 相對 merged output 的意義 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `HCC1395 5kHz paired` | `pileup_filter` | `30595` | `7450` | `8852` | `0.8042` | `0.7756` | `0.7896` | 比 merged 多 `730 TP`，但多 `6020 FP` |
| `HCC1395 5kHz paired` | `full_alignment_filter` | `30759` | `12684` | `8688` | `0.7080` | `0.7798` | `0.7422` | 比 merged 多 `894 TP`，但多 `11254 FP` |
| `HCC1395 5kHz paired` | `merged_output` | `29865` | `1430` | `9582` | `0.9543` | `0.7571` | `0.8443` | raw caller 最終輸出 |
| `HCC1395 5kHz paired` | `LongPhase-S` | `29754` | `627` | `9693` | `0.9794` | `0.7543` | `0.8522` | 比 merged 少 `111 TP`，但少 `803 FP` |
| `HCC1395 5kHz paired` | `InterSubMod` | `29752` | `544` | `9695` | `0.9820` | `0.7542` | `0.8532` | 比 merged 少 `113 TP`，但少 `886 FP` |
| `HCC1395_DORADO paired` | `pileup_filter` | `30872` | `2232` | `8575` | `0.9326` | `0.7826` | `0.8510` | 比 merged 多 `886 TP`，但多 `1644 FP` |
| `HCC1395_DORADO paired` | `full_alignment_filter` | `31332` | `2559` | `8115` | `0.9245` | `0.7943` | `0.8545` | 比 merged 多 `1346 TP`，但多 `1971 FP` |
| `HCC1395_DORADO paired` | `merged_output` | `29986` | `588` | `9461` | `0.9808` | `0.7602` | `0.8565` | raw caller 最終輸出 |
| `HCC1395_DORADO paired` | `LongPhase-S` | `29889` | `240` | `9558` | `0.9920` | `0.7577` | `0.8592` | 比 merged 少 `97 TP`，但少 `348 FP` |
| `HCC1395_DORADO paired` | `InterSubMod` | `29877` | `238` | `9570` | `0.9921` | `0.7574` | `0.8590` | 比 merged 少 `109 TP`，但少 `350 FP` |

### 4.1 這張表代表什麼

這張表最重要的意義不是「raw pileup/full 是否更會叫到 TP」，而是要量化：

1. raw model 確實保留了多少額外召回空間。
2. 這個額外召回空間要付出多少 FP 代價。
3. merged caller、LongPhase-S、InterSubMod 最後是在做哪一種 precision / recall tradeoff。

### 4.2 明確推論

1. `5kHz paired` 中，raw `full_alignment_filter` 雖然比 raw `pileup_filter` 再多 `164 TP`，但同時多 `5234 FP`，因此 raw full 在這個 dataset 明顯不值得直接放寬。
2. `DORADO paired` 中，raw `full_alignment_filter` 相對 raw `pileup_filter` 的 tradeoff 比 `5kHz` 好很多，因為它多 `460 TP`、只多 `327 FP`；但即使如此，merged output 仍有更高 F1。
3. 兩個 paired dataset 都支持同一件事：**raw pileup/full 可以證明 caller 還有召回空間，但不能直接拿來當更好的最終 benchmark 輸出。**
4. `LongPhase-S` 與 `InterSubMod` 的 paired 增益，仍主要來自**大量削減 FP**，而不是增加 TP。

---

## 5. finer feature interval 分析

## 5.1 TO：區間性明確存在

### `HCC1395 5kHz TO`

從 [fine_feature_interval_top_bins.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_interval_top_bins.tsv) 可看到：

- `GQ` 最佳區間是 `[18,20)`：`85 TP / 12 FP`
- `Quality_Score` 最佳區間是 `[55,60)`：`33 TP / 5 FP`
- `PairwiseMedianDist` 最佳區間是 `[0.18,0.20)`：`60 TP / 12 FP`
- `hp_assign_rate` 最佳區間是 `[0.5,0.7)`：`35 TP / 8 FP`

這代表：

1. `5kHz TO` 的 support 特徵確實比較像**區間訊號**，不是「越高越好」。
2. `PairwiseMedianDist` 在 `5kHz TO` 不是單純極高最好，而是中高區間 `0.18~0.20` 最集中。
3. `Quality_Score` 也不是硬要 `>=80` 才有用，反而 `55~60` 這個帶最富集 TP。

### `HCC1395_DORADO TO`

- `GQ` 最佳區間同樣落在 `[20,25)`：`15 TP / 3 FP`
- `Quality_Score` 最佳區間是 `[55,60)`：`20 TP / 5 FP`
- `PairwiseMedianDist` 最佳區間是 `[0.12,0.15)`：`31 TP / 5 FP`
- `hp_assign_rate` 最佳區間是 `[0.85,0.9)`：`10 TP / 2 FP`

這代表：

1. `DORADO TO` 對 `Quality_Score` 的 support 帶與 `5kHz TO` 接近，都偏在 `55~60`。
2. 但 `PairwiseMedianDist` 的最佳帶往較低區間移動到 `0.12~0.15`，這再次支持它是 dataset-dependent。
3. `hp_assign_rate` 也更像中高區間訊號，而不是非得 `1.0` 才最好。

## 5.2 paired：最穩的 interval 仍是 `GQ`

### `HCC1395 5kHz paired`

- `GQ [20,25)`：`35 TP / 0 FP`

### `HCC1395_DORADO paired`

- `GQ [20,25)`：`38 TP / 3 FP`

paired 的關鍵訊息是：

1. 真正穩定的正向區間仍然是 caller `GQ`。
2. 甲基特徵雖然某些 bins 的 `enrichment_ratio` 可大於 `1`，但因 coverage ceiling 與母數限制，還不足以升級成 paired 主規則。
3. `hp_assign_rate` 在 `DORADO paired` 的 `[0.99,1.01)` 有正向富集，這比較支持它做 `phase/QC annotation`，而不是 truth 主規則。

---

## 6. orthogonality 分析

## 6.1 paired：目前沒有清楚的正交補強

在 [feature_orthogonality_top_pairs.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/feature_orthogonality_top_pairs.tsv) 中，經過排除互補型閾值後，paired dataset 沒有留下明確正值的 top orthogonality pair。

這代表：

1. paired 目前的正向增益幾乎全部還是 caller `GQ` 主導。
2. 甲基特徵與 `GQ` 在 paired 目前沒有形成穩定的正交增益。
3. 這個結果要與低 coverage ceiling 一起解讀，不能直接寫成「paired 甲基一定無效」。

## 6.2 `5kHz TO`：有小幅正交，但仍屬 support

較有代表性的正交補強組合有：

| Pair | Jaccard | Union TP | Union FP | `delta_f1_vs_best_single` | 解讀 |
| --- | ---: | ---: | ---: | ---: | --- |
| `gq_ge_15 + hp_assign_ge_095` | `0.3629` | `561` | `150` | `+0.002174` | 中度重疊，union 有小幅補強 |
| `gq_ge_15 + pairwise_ge_020` | `0.3483` | `551` | `164` | `+0.001994` | `GQ` 與高 pairwise 有部分獨立訊號 |
| `pairwise_le_020 + strong_subclone` | `0.0529` | `450` | `136` | `+0.002420` | 幾乎低重疊，屬另一種 support 帶 |

這些數據表示：

1. `5kHz TO` 中，`GQ` 與某些甲基 / phase / label 特徵之間確實存在有限度的正交性。
2. 但這些正交增益仍屬小幅補強，不足以把甲基升級為第一層主規則。
3. `low VAF + high AlleleDelta + low CramersV` 在 TO 下依然比較像 artifact 旁路，不是主要的 TP support 引擎。

## 6.3 `DORADO TO`：也有正交性，但幅度更小

代表性組合：

| Pair | Jaccard | Union TP | Union FP | `delta_f1_vs_best_single` | 解讀 |
| --- | ---: | ---: | ---: | ---: | --- |
| `pairwise_le_015 + strong_subclone` | `0.0781` | `203` | `66` | `+0.000920` | 低 pairwise 與 label 結果有弱正交 |
| `pairwise_le_020 + strong_subclone` | `0.1286` | `231` | `80` | `+0.000742` | 方向一致但補強幅度較 `5kHz TO` 小 |
| `hp_assign_ge_099 + strong_subclone` | `0.2271` | `176` | `53` | `+0.000752` | `phase/QC` 類訊號可提供弱補強 |

這表示：

1. `DORADO TO` 也不是完全沒有甲基 / phase 的交叉訊號。
2. 但它們的補強幅度明顯小於 `5kHz TO`。
3. 所以跨 TO 樣本目前最合理的結論仍是：
   - 第一層：`caller-first`
   - 第二層：`methylation-support`
   - 不是 `methylation-first`

---

## 7. 現階段最合理的高層結論

1. paired raw `pileup/full` 的直接 benchmark 已正式驗證：raw model 的額外召回空間是真實存在的，但 F1 交換成本太高，因此 paired 主線仍應以 merged output、LongPhase-S、InterSubMod 為主。
2. `GQ` 仍是 paired 與 TO 都最穩的主規則。
3. TO 的 `Quality_Score`、`PairwiseMedianDist`、`hp_assign_rate`、`agreement_positive`、`Strong/Subclone` 確實存在 support 訊號，但方向與最佳區間會隨 dataset 改變。
4. `PairwiseMedianDist` 的 dataset-dependent 現象經過 finer interval 分析後更清楚，現在更不適合升級成全域硬規則。
5. paired 的甲基 rescue 目前仍主要受 coverage ceiling 限制，因此 phase 2 之後最合理的不是立刻做 paired 新 hard rule，而是：
   - 先把有效 support 升級成 annotation / ranking
   - 再視需要擴大 paired candidate-specific coverage

---

## 8. 可複查路徑

- [phase 2 執行腳本](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase2_paired_model_feature_analysis.py)
- [phase 2 summary](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/phase2_summary.md)
- [paired raw benchmark 表](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/paired_raw_model_benchmark.tsv)
- [threshold sweep 表](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_threshold_sweep.tsv)
- [interval 細表](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_interval_detail.tsv)
- [interval top bins](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_interval_top_bins.tsv)
- [orthogonality 細表](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/feature_orthogonality_detail.tsv)
- [orthogonality top pairs](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/feature_orthogonality_top_pairs.tsv)
- [上一輪四象限整合報告](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md)

