<!--
建立時間: 2026-03-07 23:35
目標: 延續 2026-03-07 pure paired 主線，驗證 `Strong` 類別是否需要細分、`低 VAF + 高 AlleleDelta` 是否適合加入原始流程判斷，以及用 haplotagged BAM 補做 read-level samtools 驗證
處理範圍:
  - 7 個 pure paired sample bundle
  - HCC1395 5kHz / HCC1395_DORADO refined label 規則比較
  - `Strong` suspect subgroup 的 TP/FP 效果
  - 代表性 region 的 haplotagged samtools snapshot
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260307_低VAF高AlleleDelta與shift跨樣本分析_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_summary.md
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_rule_comparison.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_per_region.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/samtools_snapshot_summary_tagged/snapshot_summary.tsv
-->

# Strong 細分類與 samtools 驗證分析

## 1. 本次目標

延續 pure paired 主線，這輪不是直接改核心判斷，而是先回答 4 個更關鍵的問題：

1. `低 VAF + 高 AlleleDelta` 是否適合直接加入原始流程判斷
2. `Strong` 類別是否太粗，應先拆成較可信與較可疑子群
3. `Weak->Strong / Noise->Strong` 中，哪些真的是 artifact，哪些其實是有效訊號
4. 用 haplotagged BAM 做 read-level snapshot 後，是否支持上述細分類思路

## 2. 本次新增腳本與輸出

### 2.1 新增腳本

1. [refine_strong_labels.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/refine_strong_labels.py)
   - 將 `label-first Strong` 進一步細分
   - 比較多種 refined rule 對 `TP/FP/F1` 的影響
2. [summarize_samtools_snapshots.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/summarize_samtools_snapshots.py)
   - 匯總 `region_samtools_snapshot.sh` 的結果
   - 自動輸出 `snapshot_summary.tsv` / `.md`

### 2.2 本次主要輸出

1. [refined_label_summary.md](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_summary.md)
2. [refined_label_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_rule_comparison.tsv)
3. [refined_label_best_by_sample.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_best_by_sample.tsv)
4. [refined_label_per_region.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_per_region.tsv)
5. [tagged snapshot summary](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/samtools_snapshot_summary_tagged/snapshot_summary.tsv)

## 3. refined label 的核心結果

### 3.1 主要結論

這輪結果支持：

1. `低 VAF + 高 AlleleDelta` 不適合直接升級成全域原始流程規則
2. 但它非常適合先作為 `Strong` 類別內的 suspect 子群標記
3. `Noise->Strong` 本身不能直接當過濾規則，因為它會刪掉大量 TP
4. 真正有用的是把 `Strong` 拆開，而不是把所有 upgrade 類直接刪掉

### 3.2 HCC1395 5kHz 最重要結果

在 [HCC1395 refined rule comparison](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_rule_comparison.tsv) 中：

| rule | trigger_count | TP removed | FP removed | F1 | delta vs InterSubMod |
| --- | ---: | ---: | ---: | ---: | ---: |
| `strong_low_vaf_high_ad` | 69 | 1 | 68 | 0.854006 | **+0.000806** |
| `strong_upgrade_low_vaf_high_ad` | 6 | 0 | 6 | 0.853263 | +0.000063 |
| `strong_noise_upgrade` | 1710 | 1659 | 51 | 0.825864 | -0.027336 |
| `strong_any_upgrade` | 11171 | 10948 | 223 | 0.642082 | -0.211118 |

解讀：

1. 最有用的不是「所有升級型 Strong」都刪
2. 也不是「Noise->Strong」直接刪
3. 而是 `Strong` 裡面那一小群 `低 VAF + 高 AlleleDelta` 個案
4. 這群在 `HCC1395 5kHz` 幾乎只打到 FP

### 3.3 HCC1395_DORADO 與其他樣本

在其他樣本中，這個 refined rule 沒有跨樣本穩定性：

1. `HCC1395_DORADO`
   - `strong_low_vaf_high_ad`：移掉 `9` 個 TP、`0` 個 FP
   - `F1 -0.000144`
2. `COLO829 / H1437 / H2009 / HCC1937 / HCC1954`
   - 多數樣本只有 `0~1` 個觸發或直接傷 TP
   - 沒有形成可泛化的正增益

因此目前最合理定位是：

1. `HCC1395 5kHz` 特化候選規則
2. 或 `ArtifactSuspect` 註記
3. **不是**全域原始流程預設規則

## 4. Strong 類別細分後的觀察

### 4.1 HCC1395 5kHz

在 [refined_label_per_region.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/refined_label_per_region.tsv) 中，`HCC1395 5kHz` 的 `Strong` 可分成：

| source_scope | refined_label | count |
| --- | --- | ---: |
| `tp` | `Strong-Upgraded` | 9289 |
| `tp` | `Strong-Trusted` | 6636 |
| `tp` | `Strong-NoiseUpgrade` | 1659 |
| `tp` | `Strong-SuspectLowVAF` | 1 |
| `fp` | `Strong-Trusted` | 186 |
| `fp` | `Strong-Upgraded` | 166 |
| `fp` | `Strong-SuspectLowVAF` | 62 |
| `fp` | `Strong-NoiseUpgrade` | 51 |
| `fp` | `Strong-UpgradeSuspectLowVAF` | 6 |

這代表：

1. `Strong-SuspectLowVAF` 在 `5kHz` 明顯偏向 FP
2. 但 `Strong-NoiseUpgrade` 並不乾淨，因為 TP 還很多
3. 所以 `Noise->Strong` 適合作為警示，不適合作為直接刪除條件

### 4.2 HCC1395_DORADO

`DORADO` 中：

1. `Strong-SuspectLowVAF` 幾乎全落在 TP
2. 這直接反證「低 VAF + 高 AlleleDelta」不能作全域規則

## 5. haplotagged samtools snapshot

### 5.1 驗證材料修正

本輪先用 raw BAM 跑 snapshot，發現 `HP=NA` 幾乎全滿，不能用來驗證 `label-first`。

因此後續改用 **實際參與 InterSubMod 分析的 haplotagged BAM**：

1. `HCC1395 5kHz`
   - `/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`
2. `HCC1395_DORADO`
   - `/big8_disk/mingen112/test_data/Dorado_HCC1395/ONT/somatic_tag_result/tumor/HCC1395_Tumor_ONT.GRCh38.sorted_Tmode_tagged_ClairS_pileup_v040_woTumVAF.bam`

之後的 read-level 結論都以 [tagged snapshot summary](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_refined_label_analysis/samtools_snapshot_summary_tagged/snapshot_summary.tsv) 為準。

本輪另外新增一個簡單但有效的 read-level 摘要：

1. `hp1_collapsed_count`：將 `1` 與 `1-1` 合併視為 haplotype 1
2. `hp2_collapsed_count`：將 `2` 與 `2-1` 合併視為 haplotype 2
3. `collapsed_hp_balance_delta = |hp1 - hp2| / (hp1 + hp2)`
4. `na_hp_fraction`：未分派 HP 比例

### 5.2 HCC1395 5kHz suspect FP

抽查的 5 個 `Strong-SuspectLowVAF` FP：

1. `chr8:82037650:A:G`
2. `chr8:82037629:T:A`
3. `chr8:82863512:C:T`
4. `chr8:81431203:A:G`
5. `chr8:93469697:C:T`

觀察：

1. 前 4 個 FP suspect 在 `80` 條 reads 的 HP 分布中都呈現很強的單側偏倚
   - 例如 `chr8:82037650:A:G`：`2-1=41`, `2=34`, `1-1=4`, `3=1`
   - 例如 `chr8:82863512:C:T`：`1=40`, `1-1=35`, `2-1=4`, `3=1`
2. `chr8:93469697:C:T` 則是高度 `NA` (`60/80 = 75%`)
3. 這說明 suspect FP 並不是單一型態：
   - 有些是「強單側 HP 偏倚」
   - 有些是「高未分派 HP」
4. 用 `collapsed_hp_balance_delta` 量化後：
   - `chr8:81431203:A:G`：`0.938`
   - `chr8:82037629:T:A`：`0.899`
   - `chr8:82037650:A:G`：`0.899`
   - `chr8:82863512:C:T`：`0.899`
   - `chr8:93469697:C:T`：`1.000`，且 `na_hp_fraction=0.75`
5. 也就是說，這批 `5kHz` suspect FP 是「高度 haplotype skew / 高 NA」的混合群

### 5.3 HCC1395 5kHz suspect TP

唯一被 `strong_low_vaf_high_ad` 誤傷的 TP 是：

1. `chr6:61876832:C:T`

其 HP 分布：

1. `1=38`
2. `2=27`
3. `2-1=13`
4. `1-1=2`
5. `NA=0`

和前述 FP suspect 相比，這個 TP 顯示：

1. 兩個主 haplotype 都有穩定支持
2. 不是單側偏倚
3. `collapsed_hp_balance_delta = 0.000`
3. 因此後續很值得把「低 VAF + 高 AlleleDelta」再和 HP balance 結合，避免這種少數 TP 誤傷

### 5.4 HCC1395_DORADO suspect TP

抽查 3 個 `DORADO` 的 `Strong-SuspectLowVAF` TP：

1. `chr8:111074299:T:A`
   - `1=64`, `1-1=16`
2. `chr6:61876832:C:T`
   - `1=29`, `2=25`, `2-1=9`, `NA=17`
3. `chr15:97413083:C:T`
   - `1=11`, `2-1=12`, `NA=57`

解讀：

1. `DORADO` 中就算落入 `低 VAF + 高 AlleleDelta`，也不代表是 FP artifact
2. 同一條規則在不同平台下對應到的 read-level 結構並不一致
3. 這再次支持它不能直接升級成全域原始流程條件
4. `collapsed_hp_balance_delta` 在 `DORADO` 的 3 個 TP suspect 為：
   - `chr15:97413083:C:T`：`0.043`，但 `na_hp_fraction=0.7125`
   - `chr6:61876832:C:T`：`0.079`，且 `na_hp_fraction=0.2125`
   - `chr8:111074299:T:A`：`1.000`，但仍是 TP
5. 因此 `HP balance` 很可能只在 `HCC1395 5kHz` 這個平台/資料組合下有額外辨識力，不能直接視為跨平台規則

## 6. 對原始流程的判斷

目前不建議：

1. 直接把 `低 VAF + 高 AlleleDelta` 寫成全域硬過濾
2. 直接把 `Noise->Strong` 全部刪除
3. 直接用 `upgrade` 作為刪除邏輯

目前較合理的做法是：

1. 在研究輸出中先新增 `ArtifactSuspect` 概念
2. 預設條件可先用：
   - `label_class == Strong`
   - `VAF < 0.15`
   - `AlleleDelta > 0.15`
3. 再視情況補上：
   - `HP balance`（優先在 `5kHz` 檢查）
   - `NA HP ratio`
   - `class_shift`

換句話說，下一步應優先做的是 **annotation/refinement**，不是直接改掉 InterSubMod 核心決策。

## 7. tumor-only 目前判斷

目前仍不能直接回答這個規則在 tumor-only 是否有效，原因不是概念不足，而是：

1. 目前只有 tumor-only caller 輸出目錄
2. 尚未建立 standardized 的 tumor-only InterSubMod bundle
3. 因此還沒有和 paired pure 同口徑的：
   - `label_cluster_agreement.tsv`
   - `refined_label_per_region.tsv`
   - `benchmark_comparison.tsv`

不過這輪結果已經提供了清楚的 tumor-only pilot 入口：

1. 先以 `HCC1395 5kHz` 為第一 tumor-only pilot
2. 將 TO caller 的 TP/FP benchmark VCF 串到 5kHz tagged BAM
3. 重跑一個小規模 TO label bundle
4. 驗證 `Strong-SuspectLowVAF` 是否仍富集 FP

## 8. 本輪結論

1. `低 VAF + 高 AlleleDelta` 的確隱含可辨識的 FP 訊號
2. 但這個訊號目前只在 `HCC1395 5kHz` 的 `Strong` 子群明顯有效
3. `Noise->Strong` 不等於 artifact，不能直接刪
4. `Strong` 類別值得正式再細分
5. read-level snapshot 支持「suspect Strong 並非單一型態」，至少包含：
   - 單側 HP 偏倚型
   - 高未分派 HP 型
6. 下一步最佳投資不是改核心閾值，而是：
   - 建立 `ArtifactSuspect` 註記
   - 補上 HP balance / NA HP ratio
   - 做 tumor-only pilot

## 9. 下一步建議

1. 在 round 分析層加入 `ArtifactSuspect` 註記，不先改核心過濾
2. 對 `HCC1395 5kHz` 的 `Strong-SuspectLowVAF` 再批次補更多 snapshot
3. 補做 `HP balance` 與 `NA HP ratio` 的定量欄位，驗證是否能比單純 `AD+VAF` 更穩
4. 建立 `HCC1395 5kHz` tumor-only pilot bundle，檢查這組 refined rule 是否仍富集 FP
