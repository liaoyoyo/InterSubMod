<!--
建立時間: 2026-03-11 22:10
目標: 對 HCC1395 5kHz TO 與 HCC1395_DORADO TO 的 Quality_Score、PairwiseMedianDist、hp_assign_rate 代表性 TP/FP 補做 read-level diagnostics，確認其對應現象與在流程中的合理定位
處理範圍:
  - candidate-specific representative region selection
  - matrix diagnostics
  - samtools snapshot
  - feature-level interpretation
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_GQ與甲基rescue系統性驗證_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_HCC1395_DORADO_TO_candidate_specific甲基rescue驗證_01.md
-->

# TO support 特徵 read-level diagnostics

## 1. 破題結論

本輪直接對 `HCC1395 5kHz TO` 與 `HCC1395_DORADO TO` 的代表性 `TP/FP` 位點補做矩陣與 `samtools` diagnostics 後，可以把 3 個特徵的定位講得更清楚：

1. `Quality_Score` 不是假的訊號，但**不能直接當 hard keep**。  
   - 在兩個 TO 資料集中，`Quality_Score` 高的 `TP` 通常有較完整的 alt 支持與較穩定的甲基/標籤結構。  
   - 但 `FP` 也可以同樣有很高的 `Quality_Score`，尤其當該位點同時出現極端 alt fraction、單一 haplotype 主導或明顯偏斜時。  
   - 因此 `Quality_Score` 最合理的定位是：**第一層 support / ranking 訊號**，不是 hard keep 規則。

2. `PairwiseMedianDist` 在 TO 下**不是全域規則**，而是**樣本依賴的 support / annotation**。  
   - `5kHz TO` 的正向 support 偏向 **較高 pairwise**。  
   - `DORADO TO` 的正向 support 偏向 **較低 pairwise**。  
   - 這代表 `PairwiseMedianDist` 有訊號，但方向受樣本與平台影響；目前不應升級成全域固定閾值。

3. `hp_assign_rate` 單獨看時更像 **phase/QC annotation**，不是真正的 truth-discriminative support。  
   - 高 `hp_assign_rate` 的 `TP` 與 `FP` 都很多。  
   - 特別是 `FP` 常同時伴隨極高 alt fraction 與極端 haplotype skew。  
   - 因此 `hp_assign_rate` 值得保留，但應定位為：**phase completeness / QC annotation**，不宜直接作為 rescue 主規則。

4. 兩個 TO 資料集目前共同支持的高層結論沒有變：  
   - 第一層：`caller-first`
   - 第二層：`methylation-support`
   - 獨立旁路：`artifact triage`

---

## 2. 本輪研究問題

本輪要補回答 3 個實作層問題：

1. `Quality_Score`、`PairwiseMedianDist`、`hp_assign_rate` 在代表性 `TP/FP` 的 read-level 層面各自長什麼樣。
2. 哪些特徵是「真的有 support 意義」，哪些只是「技術品質 / phase 完整度 proxy」。
3. 哪些特徵值得升成 annotation，哪些不應升成全域規則。

---

## 3. 資料來源與執行方式

### 3.1 來源資料集

1. `HCC1395 5kHz TO`
   - candidate feature 表：
     [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv)
   - snapshot 用 BAM：
     [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam)
   - 注意：這顆 BAM 來自 `ClairS-TO pileup` 主 pilot 的 full tagged BAM。

2. `HCC1395_DORADO TO`
   - candidate feature 表：
     [/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv)
   - snapshot 用 BAM：
     [/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam)
   - 注意：這顆 BAM 是 **candidate-window subset tagged BAM**，不是 full tagged BAM；優點是省空間，缺點是只覆蓋本輪 candidate windows。

### 3.2 執行腳本與主輸出

- 執行腳本：
  [/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py)
- 本輪輸出根目錄：
  [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics)
- 候選清單：
  [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/candidate_manifest.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/candidate_manifest.tsv)
- 去重後代表性位點：
  [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/selected_regions.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/selected_regions.tsv)
- read-level 彙整：
  [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/diagnostic_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/diagnostic_summary.tsv)

### 3.3 候選選取方式

每個 dataset 都固定做：

1. 只看 `candidate_eligible = True`
2. 分 `caller_lost_tp` 與 `caller_removed_fp`
3. 對每類各取：
   - `Quality_Score` top 3
   - `PairwiseMedianDist` top 3  
     - `5kHz TO` 取較高 pairwise
     - `DORADO TO` 取較低 pairwise
   - `hp_assign_rate` top 3
4. 去重後實際跑 diagnostics

最終代表性位點數量：

| dataset | `caller_lost_tp` | `caller_removed_fp` | 合計 |
| --- | ---: | ---: | ---: |
| `HCC1395 5kHz TO` | 9 | 9 | 18 |
| `HCC1395_DORADO TO` | 8 | 8 | 16 |

---

## 4. 本輪 diagnostics 看的指標代表什麼

### 4.1 `target_alt_fraction`

- 來源：`samtools mpileup`
- 定義：在 target SNV 位置，`ALT` base 佔所有 pileup base 的比例
- 單位：0 到 1 的比例
- 意義：
  - 越高代表該位置的 read-level alt 支持越集中
  - 若非常高但同時是 `FP`，常暗示局部 artifact / mapping 偏差 / 非預期訊號

### 4.2 `na_hp_fraction`

- 來源：snapshot 讀段中的 `HP` tag 統計
- 定義：沒有 `HP` tag 的 reads 比例
- 單位：0 到 1 的比例
- 意義：
  - 越高代表 phase / haplotag 完整度越差
  - 高 `NA HP` 常削弱 `label-first` 的可解釋性

### 4.3 `collapsed_hp_balance_delta`

- 定義：`HP1` 與 `HP2` 的不平衡程度，`0` 表示平衡，`1` 表示完全偏單側
- 單位：0 到 1 的比例
- 意義：
  - 高值通常代表單一 haplotype 主導
  - 但它不是單獨的真偽判別訊號，因為 `TP` 與 `FP` 都可能極度偏側

### 4.4 `Quality_Score`

- 來源：InterSubMod summary
- 意義：甲基與統計訊號的綜合品質指標
- 本輪重點：
  - 看它是否對應到更穩定的 alt 支持、較低 `NA HP` 或較好的 cluster/label 結構

### 4.5 `PairwiseMedianDist`

- 來源：read-level 距離矩陣的 median
- 意義：
  - 越高代表 read 間甲基差異越大
  - 但是否是好事，取決於 dataset / mode

### 4.6 `hp_assign_rate`

- 來源：InterSubMod summary
- 意義：可分派到 HP 的 reads 比例
- 本輪重點：
  - 看高 `hp_assign_rate` 是否真能對應更可信的 `TP`
  - 或只是代表 phase 完整度較高

---

## 5. 主要結果

### 5.1 `HCC1395 5kHz TO`

| 特徵 | 代表性 TP 現象 | 代表性 FP 現象 | 推論 |
| --- | --- | --- | --- |
| `Quality_Score` 高 | 多數 TP 有中等 alt fraction（約 `0.10` 到 `0.43`），且有 `Weak->Strong` / `Weak->Subclone` 類型 | FP 也可同樣有 `Quality=100`，其中一部分 alt fraction 很高（約 `0.48` 到 `0.55`），也可帶有 `Weak->Strong` | `Quality_Score` 有 support 訊號，但不具單獨判別力 |
| `PairwiseMedianDist` 高 | TP 常出現在 `Strong/Subclone` 或 `cluster_plus_weak_label`，但其中也混有高 `NA HP` 或極端單側 HP | FP 也能出現很高 pairwise，且常伴隨較低 alt fraction 或高 `NA HP` | 高 pairwise 在 `5kHz TO` 是真訊號，但仍需要 caller 或其他條件約束 |
| `hp_assign_rate` 高 | 能撈到 TP，但其中不少只是 `Noise/Weak` 且 alt fraction 很高 | FP 也大量命中，且常呈現極高 alt fraction + 極端 haplotype skew | `hp_assign_rate` 更像 phase completeness / QC，不適合作為 hard keep |

#### 5.1.1 5kHz TO 的代表性數值

| 特徵 | TP 平均 `target_alt_fraction` | FP 平均 `target_alt_fraction` | TP 平均 `na_hp_fraction` | FP 平均 `na_hp_fraction` | TP 平均 `hp skew` | FP 平均 `hp skew` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Quality_Score` | `0.2427` | `0.3865` | `0.1056` | `0.0861` | `0.4899` | `0.6197` |
| `PairwiseMedianDist` | `0.2371` | `0.2022` | `0.2583` | `0.3444` | `0.4082` | `0.5358` |
| `hp_assign_rate` | `0.7257` | `0.7531` | `0.0778` | `0.2417` | `0.7206` | `0.9436` |

這張表的意思是：

1. `Quality_Score` 高的 `FP` 在 alt fraction 與 HP skew 上反而更極端，這解釋了為什麼它不能直接當 keep 規則。
2. `PairwiseMedianDist` 高的 `TP` 與 `FP` 差異主要不在「pairwise 本身」，而在 `NA HP` 與 alt fraction 的附帶型態。
3. `hp_assign_rate` 高本身幾乎不區分真偽，反而需要搭配其他 read-level 訊號一起看。

### 5.2 `HCC1395_DORADO TO`

| 特徵 | 代表性 TP 現象 | 代表性 FP 現象 | 推論 |
| --- | --- | --- | --- |
| `Quality_Score` 高 | TP 多為中等 alt fraction（約 `0.24` 到 `0.29`），`NA HP` 低，HP balance 相對較正常 | FP 中有多個極端 alt fraction（約 `0.95` 到 `1.00`）且極度偏單一 HP | 在 `DORADO TO`，`Quality_Score` 的正向 support 更清楚，但仍不能直接 hard keep |
| `PairwiseMedianDist` 低 | TP 多為 `Noise/Weak`，alt fraction 偏低到中等，部分伴隨高 `NA HP` | FP 同樣常是 `Noise`，且 alt fraction 更低或同樣偏單側 HP | 低 pairwise 是 `DORADO TO` 的 dataset-specific support，不可當全域規則 |
| `hp_assign_rate` 高 | TP 與 FP 都很容易出現 `1.0`，且很多都呈現單一 HP 主導 | 對 FP 的富集尤其明顯，常伴隨極高 alt fraction | `hp_assign_rate` 應保留為 phase/QC annotation，而不是 support 主規則 |

#### 5.2.1 DORADO TO 的代表性數值

| 特徵 | TP 平均 `target_alt_fraction` | FP 平均 `target_alt_fraction` | TP 平均 `na_hp_fraction` | FP 平均 `na_hp_fraction` | TP 平均 `hp skew` | FP 平均 `hp skew` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Quality_Score` | `0.2666` | `0.7066` | `0.0111` | `0.0600` | `0.2287` | `0.7558` |
| `PairwiseMedianDist` | `0.1464` | `0.1346` | `0.3454` | `0.1972` | `0.8043` | `0.9884` |
| `hp_assign_rate` | `0.3503` | `0.7636` | `0.0139` | `0.0167` | `0.6861` | `1.0000` |

這張表的意義更直接：

1. `Quality_Score` 高的 `TP` 與 `FP` 差在 alt fraction 與 HP skew，這說明 `Quality_Score` 本身是 useful support，但需要配合 read-level artifact 訊號。
2. `PairwiseMedianDist` 低在 `DORADO TO` 可以同時出現在 `TP` 與 `FP`，因此它更適合作 dataset-aware annotation，而不是統一的 hard rule。
3. `hp_assign_rate` 高的 `FP` 比 `TP` 更極端，支持它是 phase completeness 指標，不是 truth 指標。

---

## 6. 代表性圖像與解讀

### 6.1 `5kHz TO`：高 pairwise 的 TP 代表例

代表位點：
- [chr20:28587067:G:A matrix notes](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/hcc1395_5khz_to/diagnostics/caller_lost_tp/chr20_28587067_G_A/matrix_diagnostics/region_notes.md)

![5kHz TO 高 pairwise TP 的甲基矩陣](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/hcc1395_5khz_to/diagnostics/caller_lost_tp/chr20_28587067_G_A/matrix_diagnostics/heatmap_methylation.png)

![5kHz TO 高 pairwise TP 的距離矩陣](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/hcc1395_5khz_to/diagnostics/caller_lost_tp/chr20_28587067_G_A/matrix_diagnostics/heatmap_distance.png)

圖像解釋：
- 上圖甲基矩陣：
  - X 軸：CpG 位點索引，單位為 `CpG sites`
  - Y 軸：reads 索引，單位為 `reads`
  - 顏色：甲基值，範圍約 `0` 到 `1`
- 下圖距離矩陣：
  - X / Y 軸：reads 索引
  - 顏色：read-read pairwise distance

這張圖搭配數值 `pairwise_median_dist=0.497648`、`target_alt_fraction=0.252525`、`na_hp_fraction=0.591667` 的重點是：
- 在 `5kHz TO`，高 pairwise 確實可能對應到較明顯的 read 間異質性
- 但同時也可能伴隨高 `NA HP`
- 所以它是 support，不是單獨 keep 規則

### 6.2 `DORADO TO`：低 pairwise 的 FP 代表例

代表位點：
- [chr5:78245602:A:T matrix notes](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/hcc1395_dorado_to/diagnostics/caller_removed_fp/chr5_78245602_A_T/matrix_diagnostics/region_notes.md)

![DORADO TO 低 pairwise FP 的甲基矩陣](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/hcc1395_dorado_to/diagnostics/caller_removed_fp/chr5_78245602_A_T/matrix_diagnostics/heatmap_methylation.png)

![DORADO TO 低 pairwise FP 的距離矩陣](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/hcc1395_dorado_to/diagnostics/caller_removed_fp/chr5_78245602_A_T/matrix_diagnostics/heatmap_distance.png)

圖像解釋：
- 上圖甲基矩陣：
  - X 軸：CpG 位點索引，單位為 `CpG sites`
  - Y 軸：reads 索引，單位為 `reads`
- 下圖距離矩陣：
  - X / Y 軸：reads 索引
  - 顏色：pairwise distance

這張圖搭配數值 `pairwise_median_dist=0.003273`、`target_alt_fraction=0.079365`、`na_hp_fraction=0.556962` 的重點是：
- 在 `DORADO TO`，低 pairwise 不代表真變異
- 它也可能只是讀段之間都很一致地不支持該候選
- 所以「低 pairwise 是 support」這件事只在 dataset-specific 的 aggregate 層面有微弱正訊號，不能直接升級成單位點 hard rule

---

## 7. 特徵定位建議

### 7.1 `Quality_Score`

**目前建議定位：**
- `soft support`
- `ranking`
- `annotation`

**不建議：**
- `hard keep`

**理由：**
- 在 `5kHz TO` 與 `DORADO TO` 都能救到 TP
- 但代表性 `FP` 也可以有非常高的 `Quality_Score`
- 比較合理的做法是把它保留為 `MethylSupport_QualityHigh` 類 annotation，再與 caller 指標共同判讀

### 7.2 `PairwiseMedianDist`

**目前建議定位：**
- `dataset-aware annotation`
- `support (analysis-layer only)`

**不建議：**
- 全域 `hard threshold`

**理由：**
- `5kHz TO` 與 `DORADO TO` 的方向相反
- 它有訊號，但不是跨樣本穩定訊號
- 現階段更合理的做法是：
  - 保留數值本身
  - 在報告或後續模型裡做 dataset-aware 使用

### 7.3 `hp_assign_rate`

**目前建議定位：**
- `phase completeness / QC annotation`

**不建議：**
- 當作 truth-level support 主規則

**理由：**
- 高 `hp_assign_rate` 的 `TP` 與 `FP` 都很多
- `FP` 還常伴隨極端 alt fraction 與單一 HP 主導
- 它比較像「資料是否容易被 phase」，不是「該位點是否更真」

---

## 8. 結論

本輪補完 read-level diagnostics 後，原本 aggregate 層面的推論變得更穩：

1. `Quality_Score` 值得保留在 support 層，但還不夠乾淨，不能直接 hard keep。
2. `PairwiseMedianDist` 確實有訊號，但它目前是**樣本依賴的 support / annotation**，不能全域化。
3. `hp_assign_rate` 的主要意義是 **phase/QC 完整度**，不是真偽判別主特徵。
4. 因此目前最合理的流程分工仍是：
   - 第一層：`caller-first`
   - 第二層：`Quality_Score` 等甲基 support
   - annotation 層：
     - `PairwiseMedianDist`
     - `hp_assign_rate`
     - read-level `NA HP / haplotype skew / alt fraction`

---

## 9. 後續建議

1. 對 `Quality_Score` 補做更系統的區間分析，確認是否存在比 `>=60` 更穩的 support 帶。
2. 對 `PairwiseMedianDist` 只保留數值與 dataset-aware 解讀，不先升級成規則。
3. 將 `hp_assign_rate` 明確從「support 候選」降為「phase/QC annotation」。
4. 下一輪若要改流程，優先新增 annotation 欄位，而不是直接改 hard filter / hard keep。

