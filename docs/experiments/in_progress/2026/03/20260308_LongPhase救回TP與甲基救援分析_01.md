<!--
建立時間: 2026-03-08 03:55
目標: 反向分析 LongPhase-S 是否誤刪高可信 TP，並評估 caller 與 methylation 特徵能否用於 TP rescue；同時盤點未來 LongPhase-TO + InterSubMod 的接線條件
處理範圍:
  - HCC1395 5kHz paired pure
  - ClairS caller input TP/FP vs LongPhase-S kept TP/FP
  - InterSubMod 對 caller_lost_tp / caller_removed_fp 的甲基分析
  - tumor-only 後續接線盤點
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_hcc1395/longphase_rescue_group_summary.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_hcc1395/longphase_rescue_rule_comparison.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_methylation_hcc1395/intersubmod_tp/significance_summary.csv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_methylation_hcc1395/intersubmod_fp/significance_summary.csv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_hcc1395_methylation/longphase_rescue_methylation_rule_comparison.tsv
-->

# LongPhase 救回 TP 與甲基救援分析

## 1. 本次目標

這輪改從反方向思考：

1. `LongPhase-S` 幫我們移掉很多 FP，但是否也誤刪了一批其實很可信的 TP
2. 若有，最有用的救援訊號是在 caller 端，還是在 InterSubMod 的甲基/HP 端
3. 如果這套方法有用，未來能否直接套到 `LongPhase-TO`

根據 Knowledge：

1. [LongPhase-S 文件](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase_s.md) 說明它的核心是 paired tumor-normal somatic haplotagging 與 purity estimation。
2. [LongPhase-TO 文件](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase_to.md) 說明 TO 可接受 `ClairS-TO / DeepSomatic-TO` 輸出，之後可接 `phase` 與 `haplotag`。
3. [Benchmark workflow](/big8_disk/liaoyoyo2001/knowledge/06_workflows/benchmark_workflow.md) 要求固定 truth scope、PASS filter 與 TP/FP/FN 口徑。

## 2. 本次新增腳本

1. [analyze_longphase_rescue.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_longphase_rescue.py)
   - 比較 caller TP/FP 與 LongPhase kept TP/FP
   - 輸出 `caller_lost_tp` 與 `caller_removed_fp` 的 VCF/TSV
   - 直接量化 caller rescue 規則對 F1 的影響
2. [analyze_longphase_rescue_with_methylation.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_longphase_rescue_with_methylation.py)
   - 將 InterSubMod 的 `VerificationClass / PairwiseMedianDist / AlleleDelta` 與 caller 特徵合併
   - 比較 caller-only 與 caller+methylation rescue 規則

## 3. 分析方法

### 3.1 資料來源

使用 `HCC1395 5kHz` paired pure：

1. caller input TP：`/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/benchmark-test/tp.vcf`
2. caller input FP：`/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/benchmark-test/fp.vcf`
3. LongPhase-S kept TP：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395/20260211/longphase_s/filtered_snv_tp.vcf.gz`
4. LongPhase-S kept FP：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure/HCC1395/20260211/longphase_s/filtered_snv_fp.vcf.gz`
5. tagged BAM：`/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`

### 3.2 定義

1. `caller_lost_tp`
   - 在 caller benchmark TP 中
   - 但不在 LongPhase-S kept TP 中
2. `caller_removed_fp`
   - 在 caller benchmark FP 中
   - 但不在 LongPhase-S kept FP 中

### 3.3 注意事項

本次集合比對採用 **精確 key**：

`Chr:Pos:Ref:Alt`

因此和整體 benchmark aggregate count 有 `2 TP / 2 FP` 的差異：

1. benchmark 表中的 `LongPhase-S TP=29754, FP=627`
2. 精確 key 交集得到：
   - `caller_kept_tp=29752`
   - `caller_lost_tp=113`
   - `caller_kept_fp=625`
   - `caller_removed_fp=805`

這表示仍有少量表示法或 key 不完全一致的位點。  
本輪 rescue 分析只基於可精確對齊的集合。

## 4. caller 端 rescue 結果

主表在 [longphase_rescue_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_hcc1395/longphase_rescue_rule_comparison.tsv)。

### 4.1 群組摘要

來自 [longphase_rescue_group_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_hcc1395/longphase_rescue_group_summary.tsv)：

| group | count | median AF | median DP | median GQ | median ALT count | H fraction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `caller_kept_tp` | 29752 | 0.4328 | 65 | 21 | 27 | 0.5737 |
| `caller_lost_tp` | 113 | 0.5405 | 52 | 17 | 29 | 0.4336 |
| `caller_kept_fp` | 625 | 0.1562 | 70 | 11 | 11 | 0.4464 |
| `caller_removed_fp` | 805 | 0.1429 | 43 | 10 | 6 | 0.4571 |

解讀：

1. `caller_lost_tp` 的 AF、ALT count 明顯不低
2. 這批 lost TP 並不是單純低支持雜訊
3. 最明顯的分界是 `GQ`
   - lost TP median GQ = `17`
   - removed FP median GQ = `10`

### 4.2 最好的 caller rescue rule

| rule | TP rescued | FP reintroduced | F1 | delta vs LongPhase-S |
| --- | ---: | ---: | ---: | ---: |
| `gq_ge_20` | 45 | 0 | 0.852947 | **+0.000739** |
| `af_ge_0_20_and_gq_ge_20` | 45 | 0 | 0.852947 | **+0.000739** |
| `af_ge_0_15_and_alt_ge_10` | 91 | 107 | 0.852398 | +0.000189 |

最重要結論：

1. **單純 `GQ >= 20` 已經能有效救回一批被 LongPhase-S 誤刪的 TP**
2. 而且在這批精確對齊集合中，`0` 個 removed FP 被放回來
3. 代表現在最有力的 TP rescue 訊號其實在 **caller 端**

## 5. methylation / label-first 是否能幫忙 rescue

我另外把：

1. `caller_lost_tp.vcf`
2. `caller_removed_fp.vcf`

重新丟給 InterSubMod，輸出在：

1. [lost TP significance summary](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_methylation_hcc1395/intersubmod_tp/significance_summary.csv)
2. [removed FP significance summary](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_methylation_hcc1395/intersubmod_fp/significance_summary.csv)

### 5.1 類別分布

`caller_lost_tp`：

| VerificationClass | count |
| --- | ---: |
| Noise | 48 |
| Weak | 35 |
| Strong | 22 |
| Subclone | 6 |

`caller_removed_fp`：

| VerificationClass | count |
| --- | ---: |
| Noise | 479 |
| Strong | 138 |
| Weak | 103 |
| Subclone | 82 |

解讀：

1. `lost TP` 裡確實有一部分會被 InterSubMod 判成 `Strong/Subclone`
2. 但 `removed FP` 裡同樣也有很多 `Strong/Subclone`
3. 這表示目前的 label-first class **不適合直接當 LongPhase rescue 規則**

### 5.2 caller + methylation 聯合 rescue

主表在 [longphase_rescue_methylation_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_longphase_rescue_hcc1395_methylation/longphase_rescue_methylation_rule_comparison.tsv)。

| rule | TP rescued | FP reintroduced | delta vs LongPhase-S |
| --- | ---: | ---: | ---: |
| `gq_ge_20` | 45 | 0 | **+0.000739** |
| `gq20_and_pairwise_ge_020` | 23 | 0 | +0.000378 |
| `gq20_and_strong_or_subclone` | 14 | 0 | +0.000230 |
| `strong_or_subclone_and_af_ge_020` | 25 | 29 | +0.000057 |
| `pairwise_ge_020` | 59 | 292 | -0.002581 |

最重要結論：

1. 目前 **最好的 rescue 規則仍然是 caller `GQ >= 20`**
2. 加入目前的 methylation label 或 pairwise distance 後，**沒有超過 `GQ` 單獨使用**
3. `Strong/Subclone` 作為 rescue 條件太寬，會把太多 removed FP 放回來
4. 所以現階段的甲基特徵較適合：
   - 幫 InterSubMod 做 FP triage
   - 幫 `Strong` 細分類
   - **不適合直接當 LongPhase-S 的 TP rescue 主規則**

## 6. 對 LongPhase-S / LongPhase-TO 的推論

### 6.1 LongPhase-S

目前最合理的推論是：

1. `LongPhase-S` 很可能對一部分 `caller 高 GQ TP` 過於保守
2. 若要救回 TP，第一優先應該先檢查：
   - `GQ`
   - `AF`
   - `ALT count`
3. InterSubMod 的甲基訊號目前比較像第二層輔助，而不是第一層 rescue 規則

### 6.2 LongPhase-TO

這輪還沒有實際跑 `LongPhase-TO + InterSubMod`，但材料與流程已更清楚：

1. `HCC1395 ONT ClairS-TO` 已有：
   - `snv.vcf.gz`
   - `benchmark-test/tp.vcf`
   - `benchmark-test/fp.vcf`
   - `benchmark-test/fn.vcf`
2. `HCC1395_DORADO ClairS-TO` 目前只有 `snv.vcf.gz`，還缺現成 split TP/FP benchmark 檔
3. 也就是說：
   - **ONT 版本可先做 TO pilot**
   - `DORADO TO` 需先補 benchmark split

此外，本輪新增的兩個 rescue 腳本本質上是 **LongPhase-S/TO 共用框架**：

1. 只要有 caller tp/fp 與 LongPhase kept tp/fp
2. 就能直接重用到 `LongPhase-TO`

## 7. 本輪結論

1. 反向看 `LongPhase-S` 是值得的，因為它確實有誤刪一批高可信 TP
2. 在 `HCC1395 5kHz`，最佳 TP rescue 訊號目前是 caller `GQ >= 20`
3. 單靠目前的 InterSubMod class 或 pairwise 特徵，**還不能比 `GQ` 更好地救回 TP**
4. 這表示目前甲基訊號對 `LongPhase-S` 的最強價值，仍偏向：
   - 區分/標記 artifact-prone strong
   - 幫助 FP triage
   - 而不是直接接管 TP rescue
5. `LongPhase-TO` 的 pilot 已具體化，但先以 `HCC1395 ONT` 為起點較合理

## 8. 下一步建議

1. 先回頭檢查 `LongPhase-S` 對 `GQ >= 20` TP 的丟棄邏輯，確認是哪些 internal condition 在誤刪
2. 補做：
   - `caller_lost_tp` 的 samtools snapshot
   - 尤其是 `GQ >= 20` 但被 LongPhase-S 丟掉的那 `45` 筆
3. 在 `LongPhase-S` 研究線上新增一個候選想法：
   - `caller high-GQ rescue gate`
4. 後續有資源時，優先啟動：
   - `HCC1395 ONT ClairS-TO -> LongPhase-TO -> InterSubMod` 的 tumor-only pilot
