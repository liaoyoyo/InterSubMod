<!--
建立時間: 2026-03-07 20:10
目標: 確認 LongPhase-TO + InterSubMod 是否必須產生 haplotagged BAM、目前磁碟空間是否足夠，以及用既有 tumor-only 素材先判斷 LongPhase-TO 純樣本主線的可行性
處理範圍:
  - LongPhase-TO / InterSubMod 接線需求
  - /big8_disk 與 /bip8_disk 可用空間
  - HCC1395 ONT / HCC1395_DORADO tumor-only 現有素材
  - 既有 phased 中間產物是否可直接拿來 benchmark
關聯檔案:
  - /big8_disk/liaoyoyo2001/knowledge/05_tools/longphase_to.md
  - /big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing_workflow.md
  - /big8_disk/liaoyoyo2001/knowledge/03_file_formats/vcf_clairs_to.md
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/pipeline/steps/02_intersubmod.sh
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/
-->

# LongPhase-TO 空間需求與中間產物確認

## 1. 本次要回答的問題

1. 跑 `LongPhase-TO + InterSubMod` 時，是否一定要先產生 `LongPhase-TO haplotagged BAM`。
2. 目前 `/big8_disk` 與 `/bip8_disk` 的剩餘空間，是否足以支撐這條流程。
3. 在不先產生 200GB 級 tagged BAM 的情況下，能否先延續前次觀察，判斷 `LongPhase-TO` 在 pure tumor 上是否更值得投資。

## 2. 直接依據

1. 根據 [LongPhase-TO 文件](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase_to.md)，標準 tumor-only 流程是：
   `phase -> haplotag -> tagged BAM`。
2. 根據 [Phasing workflow](/big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing_workflow.md)，tumor-only 標準輸出為：
   `phased.vcf + tagged.bam + LOH.bed`。
3. 根據 [02_intersubmod.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/pipeline/steps/02_intersubmod.sh)，目前 repo pipeline 入口直接要求：
   `--tagged-bam`、`--tp-vcf`、`--fp-vcf`。

## 3. 本次實際執行

### 3.1 空間與 BAM 體積

2026-03-07 20:06 CST 實測：

| 項目 | 結果 |
|---|---:|
| `/big8_disk` 可用空間 | `487G` |
| `/bip8_disk` 可用空間 | `481G` |
| 現有 `HCC1395 5kHz` LongPhase-S tagged BAM | `263.63 GiB` |
| 現有 `HCC1395_DORADO` LongPhase-S tagged BAM | `225.97 GiB` |
| `HCC1395 ONT/HCC1395.bam` | `223G` |
| `HCC1395 ONT_Dorado/HCC1395.bam` | `233G` |
| `HCC1395 ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam` | `272G` |

實際路徑：

- [HCC1395 ONT tumor BAM](/big8_disk/data/HCC1395/ONT/HCC1395.bam)
- [HCC1395 ONT normal BAM](/big8_disk/data/HCC1395/ONT/HCC1395BL.bam)
- [HCC1395 5kHz tumor BAM](/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam)
- [HCC1395 5kHz normal BAM](/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam)
- [HCC1395_DORADO tumor BAM symlink](/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam)
- [HCC1395_DORADO normal BAM symlink](/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam)
- [HCC1395 5kHz tagged BAM](/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam)
- [HCC1395_DORADO tagged BAM](/big8_disk/mingen112/test_data/Dorado_HCC1395/ONT/somatic_tag_result/tumor/HCC1395_Tumor_ONT.GRCh38.sorted_Tmode_tagged_ClairS_pileup_v040_woTumVAF.bam)

### 3.2 TO 現有素材盤點

#### HCC1395 ONT

- 有 [ClairS-TO baseline VCF](/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz)
- 有 [tmp/phasing_output](/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/tmp/phasing_output/)
- 有既有 `benchmark-test/`，但原目錄下 `tp/fp/fn.vcf` 都是 `0 byte`
- 本次重新跑了 baseline benchmark，產出可用 TP/FP/FN：
  - [ONT baseline TP](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_ONT/benchmark_baseline_clairs_to/tp.vcf)
  - [ONT baseline FP](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_ONT/benchmark_baseline_clairs_to/fp.vcf)
  - [ONT baseline FN](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_ONT/benchmark_baseline_clairs_to/fn.vcf)

#### HCC1395_DORADO

- 有 [ClairS-TO baseline VCF](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/snv.vcf.gz)
- 有 per-chromosome [phased_vcf_output](/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/tmp/phasing_output/phased_vcf_output/)
- 原目錄下沒有現成 `benchmark-test/tp.vcf / fp.vcf / fn.vcf`
- 本次已補跑 baseline benchmark，產出可用 TP/FP/FN：
  - [DORADO baseline TP](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_DORADO/benchmark_baseline_clairs_to/tp.vcf)
  - [DORADO baseline FP](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_DORADO/benchmark_baseline_clairs_to/fp.vcf)
  - [DORADO baseline FN](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_DORADO/benchmark_baseline_clairs_to/fn.vcf)

### 3.3 既有 phased 中間產物檢查

本次將 [HCC1395 ONT merged.vcf.gz](/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/tmp/phasing_output/merged.vcf.gz) 抽成單一 sample 檔：

- [longphase_to_phased_single_sample.vcf.gz](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_ONT/longphase_to_phased_single_sample.vcf.gz)
- [longphase_to_phased_single_sample_nonmissing.vcf.gz](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/bip8_disk_output/research_rounds/20260307_longphase_to_space_eval/HCC1395_ONT/longphase_to_phased_single_sample_nonmissing.vcf.gz)

檢查結果：

1. `merged.vcf.gz` header 有 `24` 個 sample columns，不是單純一個 sample 的最終 callset。
2. 第一個 sample 有效記錄只有 `80,556` 筆，缺值記錄有 `1,990,672` 筆。
3. 因此 `merged.vcf.gz` 不能直接視為 LongPhase-TO 的最終 somatic VCF。

## 4. 本次 benchmark 結果

### 4.1 HCC1395 ONT ClairS-TO baseline

實際 benchmark：

- `TP=28509`
- `FP=11606`
- `FN=10938`
- `F1=0.7166`

這與前面研究紀錄一致，可當作 ONT TO baseline。

### 4.2 HCC1395_DORADO ClairS-TO baseline

實際 benchmark：

- `TP=28861`
- `FP=11576`
- `FN=10586`
- `F1=0.7226`

這是本次新補齊的 methylation-capable TO baseline。

### 4.3 HCC1395 ONT `merged.vcf.gz` 非缺值子集

將非缺值子集直接拿去 benchmark，結果為：

- `TP=1135`
- `FP=74503`
- `FN=38312`
- `F1=0.0197`

這不是 LongPhase-TO 很差，而是證明：

> **`tmp/phasing_output/merged.vcf.gz` 不是可直接拿來宣稱 LongPhase-TO 效果的最終 somatic callset。**

## 5. 結論

### 5.1 跑 `LongPhase-TO + InterSubMod` 是否一定要 tagged BAM

分兩種情況：

1. **若只是跑 InterSubMod 的 cluster-first 甲基分析**
   - 技術上可以直接吃有 `MM/ML` 的 tumor BAM。
   - 但這不叫做 `LongPhase-TO + InterSubMod`，因為沒有用到 TO 的 HP 標記。

2. **若要做真正的 `LongPhase-TO + InterSubMod` 驗證**
   - 目前答案是 **要**。
   - 原因是你要驗證的是 `TO phasing / haplotagging` 和 `InterSubMod label-first` 是否一起提供新訊息。
   - 目前 repo pipeline 也直接要求 `--tagged-bam`。

### 5.2 目前空間夠不夠

結論是：

1. **一次只跑一個 TO haplotagged BAM，空間夠。**
2. **不要同時啟動兩個 200GB 級 tagged BAM run。**
3. **5kHz 可跑，但最緊。**

保守估計：

| 目標樣本 | 預估 tagged BAM 體積 | 判斷 |
|---|---:|---|
| HCC1395 ONT | `223-240 GiB` | 可跑，但不能當 InterSubMod 主驗證 |
| HCC1395_DORADO | `233-245 GiB` | 可跑，且是目前最實際的 TO+甲基候選 |
| HCC1395 5kHz | `272-285 GiB` | 可跑，但空間與上游準備壓力最大 |

**建議門檻**：

- 啟動 TO haplotag 前，目標掛載點至少保留 `>300 GiB`
- 若還要保留其他大型中間檔，建議保留 `>350 GiB`

### 5.3 現在對「LongPhase-TO 純樣本是否更有效」的判斷

目前只能得到這個層級的結論：

1. `ClairS-TO` baseline 已確認：
   - `HCC1395 ONT`：`F1=0.7166`
   - `HCC1395_DORADO`：`F1=0.7226`
2. 既有 `tmp/phasing_output/merged.vcf.gz` **不能**直接拿來宣稱 LongPhase-TO 是否更有效。
3. 因此，**現階段還不能下結論說 LongPhase-TO pure sample 比 paired pure 更有效，或比 baseline TO 更有效。**
4. 但現在已明確知道：
   - `HCC1395 ONT` 只能先做 TO caller/phasing 邏輯研究，不適合作 InterSubMod 主驗證，因為 [HCC1395 樣本文件](/big8_disk/liaoyoyo2001/knowledge/02_samples/HCC1395.md) 已明確說明 `ONT/HCC1395.bam` **沒有 MM/ML**
   - `HCC1395_DORADO` 是目前最值得先做的 methylation-capable TO pilot
   - `HCC1395 5kHz` 仍是主實驗樣本，但若要做 TO，還要先補 upstream TO caller 結果，且空間風險更高

## 6. 對既有觀察與結論的承接

這次結果應該記住三件事，之後研究 TO 都不要重複踩坑：

1. `低 VAF + 高 AlleleDelta`、`GQ >= 20 rescue`、`Weak->Strong / Noise->Strong` 的觀察，要留到 **真正的 TO tagged BAM + TP/FP VCF** 出來後再驗證。
2. **不要再把 `tmp/phasing_output/merged.vcf.gz` 直接當 final LongPhase-TO 結果。**
3. 先做 `HCC1395_DORADO` 的 TO pilot，比直接硬上 `5kHz TO` 更合理，因為：
   - 已有 TO caller VCF
   - 已補齊 baseline benchmark split
   - BAM 體積比 5kHz 小
   - 有 `MM/ML`

## 7. 下一步建議

1. 以 `HCC1395_DORADO` 做第一個 `LongPhase-TO + InterSubMod` 正式 pilot。
2. 先只跑一個樣本、一個輸出目錄，不並行第二個 tagged BAM 任務。
3. 產出真正的 `tagged BAM` 後，再生成 TO 版 `TP/FP/FN` VCF，接續：
   - `label-first / cluster-first`
   - `低 VAF + 高 AlleleDelta`
   - `GQ rescue`
   - `Weak->Strong / Noise->Strong`
4. `HCC1395 5kHz TO` 保留為下一階段主驗證，等 DORADO pilot 接線成功後再做。
