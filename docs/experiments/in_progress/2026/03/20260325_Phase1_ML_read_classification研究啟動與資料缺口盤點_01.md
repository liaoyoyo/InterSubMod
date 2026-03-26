<!--
建立時間: 2026-03-25 14:35
目標: 依據 2026-03-24 收斂出的突破方向，正式啟動 Phase 1（ML read classification）研究，先盤點既有 read-level 資產、定義資料單位與確認主要缺口
處理範圍:
  - Phase 1 啟動問題定義
  - HCC1395 5kHz / DORADO paired/TO 資料角色
  - read-level / region-level 既有資產盤點
  - training schema 與 exporter 缺口
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_甲基候選研究框架與觀察資產盤點_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase1_training_manifest.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_read_training_table.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_region_diagnostics.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py
-->

# Phase 1 ML read classification 研究啟動與資料缺口盤點

## 1. 破題結論

目前真正的任務不是直接訓練一個模型，而是先把 `Phase 1: ML read classification` 的資料層補齊。

這輪研究已確認 3 件事：

1. 主線已正式切到 `方向 E`，目標是用 read-level pattern modeling 突破既有 unsupervised clustering 的辨識天花板，而不是再做舊規則微調。
2. repo 內其實已經有大量可重用的 read-level 與 region-level 資產：
   - `rescue_joined_features.tsv`
   - `reads.tsv`
   - `methylation.csv`
   - `distance/.../matrix.csv`
   - 既有 `TO support` 與 `paired old output` diagnostics
3. 真正的阻塞點不是「沒有資料」，而是「還沒有統一的 per-read training/export layer」。

因此本輪的近期目標是：

1. 定義 `Phase 1` 的資料單位與 label schema
2. 明確指定 discovery / validation dataset
3. 盤點哪些欄位已存在、哪些欄位必須由 exporter 補齊
4. 將下一步任務收斂成可實作的 exporter 與 pilot manifest

---

## 2. 背景摘要

### 2.1 為什麼現在要做 Phase 1

根據 [20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md)，ISM 現階段最大的瓶頸不是缺少局部 read-level 現象，而是這些現象尚未被提升成比 `clustering + significance` 更穩健的 read classifier。

同一份全域分析已把突破方向優先序定為：

1. `Phase 1`：`ML read classification`
2. `Phase 2`：`normal methylation reference + CN/Purity-aware`
3. `Phase 3`：`gene-level evidence integration`
4. `Phase 4`：`CpG 功能分層與智慧選點`

所以現在的研究任務，已經不是重新證明 `label-first / cluster-first` 合不合理，而是把既有 read-level 觀察真正轉成可 head-to-head 的新方法原型。

### 2.2 為什麼可以從 HCC1395 系列開始

根據 Knowledge [02_samples/HCC1395.md](/big8_disk/liaoyoyo2001/knowledge/02_samples/HCC1395.md)：

1. `HCC1395 5kHz` 有 `5mCG + 5hmCG` MM/ML 標籤，適合做 discovery 與 stress-test。
2. `HCC1395_DORADO` 是同 biological sample 的 cross-platform 對照，但 BAM 只有 `5mCG`，適合做 validation，不適合直接把數值分布和 `5kHz` 當同母體看待。
3. `HCC1395/HCC1395BL` 仍是目前最穩定的 benchmark 樣本組，搭配 SEQC2 high-confidence truth 可維持固定 benchmark 口徑；這點也符合 Knowledge [06_workflows/benchmark_workflow.md](/big8_disk/liaoyoyo2001/knowledge/06_workflows/benchmark_workflow.md) 的既定流程。

### 2.3 既有研究已經提供哪些 Phase 1 前提

根據 [20260315_甲基候選研究框架與觀察資產盤點_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_甲基候選研究框架與觀察資產盤點_01.md) 與 [20260311_TO_support特徵_readlevel_diagnostics_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md)：

1. `paired old output` 與 `TO support diagnostics` 已經提供足夠的 per-region matrix 與 read-level 觀察。
2. 甲基訊號目前最合理的高層定位仍是：
   - 第一層：`caller-first`
   - 第二層：`methylation-support`
   - 第三層：`annotation / QC / artifact triage`
3. 真正值得突破的方向，是 read-level pattern recognition，不是再把某個單欄位升成 hard rule。

---

## 3. 重要術語表

| 術語 | 中文詞義 | 在本輪的角色 |
| --- | --- | --- |
| `Phase 1` | 第一階段突破方向 | 目前主研究任務，目標是建立 ML read classifier 原型 |
| `read-level training unit` | 單一 read 在單一 SNV 視窗中的訓練單位 | 本輪最先要定義的資料單位 |
| `region-level label` | 一個 SNV 視窗對應的真值或流程狀態標籤 | 目前最容易從既有資料表繼承的 supervision 來源 |
| `rescue_joined_features.tsv` | 候選位點與 InterSubMod 摘要欄位整合表 | 可作為 region-level baseline feature 與 label 對照來源 |
| `reads.tsv` | 區域內每條 read 的基本資訊表 | 可作為 per-read exporter 的主鍵與 read metadata 來源 |
| `methylation.csv` | 每條 read 在各 CpG 上的甲基化矩陣 | 可作為 ML 模型的主要 pattern 輸入來源 |
| `head-to-head` | 與既有 baseline 同表比較 | Phase 1 必須同時比較 `VerificationClass / QualityScore / PairwiseMedianDist` |
| `discovery dataset` | 第一輪找方法與調 schema 的資料集 | 以 `HCC1395 5kHz paired/TO` 為主 |
| `validation dataset` | 驗證跨平台穩定性的資料集 | 以 `HCC1395_DORADO paired/TO` 為主 |

---

## 4. 已知前提與適用邊界

### 4.1 已驗證成立的前提

1. 舊的 `paired pure` 主線與 TO baseline 已足夠作 Phase 1 的固定對照，不需要再先重跑一次 baseline。
2. `rescue_joined_features.tsv` 已提供 region-level 的 baseline 特徵與狀態欄位，可直接承接：
   - `truth_status`
   - `downstream_status`
   - `VerificationClass`
   - `PairwiseMedianDist`
   - `Quality_Score`
   - `agreement_type`
   - `hp_assign_rate`
3. 既有 region output 已暴露：
   - `reads/reads.tsv`
   - `methylation/methylation.csv`
   - `distance/<metric>/matrix.csv`
   - `metadata.txt`

### 4.2 必須保守處理的邊界

1. `5kHz TO` 的 read-level diagnostics 目前來自 full tagged BAM，而 `DORADO TO` 多數是 candidate-window subset tagged BAM，因此 read-level 絕對值不可直接橫比，只能做方向比較。
2. `HCC1395 5kHz` 含 `5mCG + 5hmCG`，`DORADO` 目前只有 `5mCG`。若直接把 matrix 數值當同分布輸入，會混入 chemistry / basecalling 差異。
3. `paired` caller 邊界的很多位點缺少可用甲基摘要，不能把 paired rescue 的限制誤解成「read-level ML 一定無效」；更準確的說法是目前 exporter 與 coverage 層還沒補齊。
4. `VerificationClass`、`agreement_type`、`PairwiseMedianDist` 是 baseline 或弱 supervision 候選，不是真正的 gold truth。

### 4.3 本輪明確不做

1. 不重新打開 `H006 window_bp=1000`
2. 不重做 `AlleleDelta standalone filter`
3. 不把 `TO support` 的 sample-specific 規則直接升成新的全域策略
4. 不先跳到 `Phase 2` 的 normal baseline / CN-aware 建模

---

## 5. 本輪研究問題、假設與驗收條件

### 5.1 研究問題

1. 既有 region outputs 與 joined feature tables，能否被重組成統一的 `per-read training table`？
2. 現階段最合理的 supervision，應該是 `region truth / downstream status`、`read alt/ref label`，還是 `VerificationClass` 類的弱 supervision？
3. 在不重跑大型 pipeline 的前提下，能否先用 `HCC1395 5kHz paired/TO` 與 `HCC1395_DORADO paired/TO` 建出第一個 Phase 1 pilot manifest？

### 5.2 研究假設

1. `reads.tsv + methylation.csv + metadata.txt + rescue_joined_features.tsv` 已足夠組成第一版 `read × region` 訓練單位。
2. 本輪最大的缺口是 exporter 與 schema，不是原始觀察資產不足。
3. 若先把 `5kHz paired/TO` 當 discovery、`DORADO paired/TO` 當 validation，可以在不過度放大單樣本結論的前提下啟動 Phase 1。

### 5.3 成功條件

1. 明確定義第一版 training unit、主鍵、欄位與 label schema。
2. 明確指定 discovery / validation dataset 與不該混用的資料邊界。
3. 產出一份 exporter 導向的近期任務清單，能直接接到下一輪腳本實作。

### 5.4 失敗條件

1. 若無法從現有資料穩定對應 `region_key -> region_dir -> read rows`，則 Phase 1 需先補資料索引層。
2. 若 read-level 單位無法定義清楚 supervision，則不能直接宣稱可開始模型訓練。
3. 若 `5kHz` 與 `DORADO` 的輸入結構差異大到無法 harmonize，則須先拆成 dataset-specific pilot，不可硬做 joint training。

### 5.5 本輪評估指標

本輪先不評估模型準確度，而先評估研究啟動是否成功：

1. `資料可對應率`：多少 region 可從 feature table 對應回完整 region outputs
2. `欄位完備度`：training schema 中多少欄位可由既有資產直接提供
3. `樣本覆蓋度`：四個固定 baseline dataset 中，哪些已可直接進 pilot manifest
4. `邊界清晰度`：哪些資料可 cross-platform 比較、哪些只能 within-dataset 使用

---

## 6. 既有資產盤點

### 6.1 region-level baseline feature table 已存在

從 `HCC1395 5kHz TO` 的 `rescue_joined_features.tsv` 可直接讀到下列欄位：

- `sample`
- `platform`
- `caller`
- `mode`
- `region_key`
- `truth_status`
- `downstream_status`
- `qual`
- `gq`
- `dp`
- `af`
- `ad_ref`
- `ad_alt`
- `candidate_eligible`
- `PairwiseMeanDist`
- `PairwiseMedianDist`
- `AlleleDelta`
- `CramersV`
- `VerificationClass`
- `DominantLabel`
- `PassedGating`
- `GlobalP`
- `Quality_Score`
- `agreement_type`
- `class_shift`
- `hp_assign_rate`
- `allele_assign_rate`
- `cluster_class`
- `label_class`

這代表 baseline comparison 所需的 region-level 摘要欄位，其實已經齊全。

### 6.2 read-level 主表已存在

現有 region output 的 `reads.tsv` 目前已提供：

- `read_id`
- `read_name`
- `chr`
- `start`
- `end`
- `mapq`
- `hp`
- `alt_support`
- `is_tumor`
- `strand`

這些欄位足以作為：

1. read 主鍵與 region 內排序依據
2. `ALT / REF / UNKNOWN` 的 read-level allele 標記
3. `HP`、`strand`、`mapq` 等非甲基 covariates

### 6.3 甲基矩陣已存在

現有 `methylation.csv` 的結構是：

1. 第一欄是 `read_id`
2. 後續每欄是一個 CpG genomic position
3. 值域為 `0~1` 的 methylation probability，缺值為 `NA`

這意味著現有輸出已經足以支持兩種路線：

1. `tabular summary route`
   - 先把每條 read 壓成 `mean / std / coverage / high-methyl fraction`
2. `sequence / matrix route`
   - 直接把 per-read CpG vector 當模型輸入

### 6.4 既有腳本可重用，但還不是 exporter

目前最接近 Phase 1 資料層的既有腳本有兩類：

1. [export_region_diagnostics.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_region_diagnostics.py)
   - 已能從既有 region output 匯出：
     - `region_summary.tsv`
     - `heatmap_methylation.png`
     - `heatmap_distance.png`
     - `label_vs_cluster.png`
   - 但它是 diagnostics exporter，不是 training-table exporter
2. [run_to_support_feature_diagnostics.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_to_support_feature_diagnostics.py)
   - 已能從 `rescue_joined_features.tsv` 選 representative regions，並回接 diagnostics
   - 但它做的是代表位點挑選，不是把所有 read flatten 成訓練表

---

## 7. 目前最明確的缺口

### 7.1 缺少統一的 per-read exporter

目前沒有一個腳本能直接輸出：

- 一列 = 一個 `read × region`
- 同時帶有 `region-level baseline features`
- 同時帶有 `read-level metadata`
- 同時帶有 `methylation summary / CpG vector`
- 同時標示 `dataset_id / mode / truth supervision`

這是 Phase 1 目前最核心的工程缺口。

### 7.2 缺少明確的 label schema

目前至少有 4 種可能 supervision：

1. `region truth_status`
2. `region downstream_status`
3. `read alt_support`
4. `VerificationClass / agreement_type / class_shift`

截至本輪 sample80 驗證後，這個問題已不再只是理論推論，而有明確收斂：

1. `Phase 1A` 可立即定義為 `within-tumor alt-support read classification`
   - read-level 主 label：`ALT / REF`
   - region context：`truth_status / VerificationClass / Quality_Score / PairwiseMedianDist`
2. `Phase 1B` 目前不能直接定義為 `tumor vs normal read classification`
   - 原因不是欄位沒匯出，而是目前 read exports 內的 `is_tumor` 在 sampled shard `14,165` 筆 reads 中全部為 `1`
   - 這表示現有輸出 universe 仍是 tumor-only，若要做 Phase 1B，必須先補 normal-read export layer
3. 因此第一版最穩妥的做法，已從「multi-label export 後再決定」收斂成：
   - 先把 `Phase 1A` 當立即可執行任務
   - 將 `Phase 1B` 明確標記為資料層待補的後續任務

正式規格已另寫於：

- [20260325_Phase1_label_schema與harmonization規格_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1_label_schema與harmonization規格_01.md)

### 7.3 缺少資料切分與 harmonization 規則

截至 sample80 驗證後，這部分也已有可執行結論，而不是待討論問題：

1. `5kHz` 與 `DORADO` 第一版不能直接 joint pool
   - sampled shard 顯示 `num_cpg_observed / methyl_mean / methyl_std` 在 `platform|mode` 間有明顯 shift
2. `paired` 與 `TO` 也不能先假設為同分布
3. 第一版至少要保留：
   - `dataset_role`
   - `harmonization_group = platform|mode`
4. 第一版 split 應收斂為：
   - discovery：`HCC1395 5kHz paired/TO`
   - external validation：`HCC1395_DORADO paired/TO`
5. 若未來要做 joint train，必須至少補：
   - platform-aware feature handling
   - mode-aware feature handling
   - group-aware normalization

---

## 8. 第一版資料 schema 建議

### 8.1 建議訓練單位

第一版先定義為：

`一列 = 一個 read 在一個 region_key 下的觀測`

主鍵建議為：

- `dataset_id`
- `mode`
- `region_key`
- `read_name`

### 8.2 必備欄位

#### A. dataset / region metadata

- `dataset_id`
- `sample`
- `platform`
- `mode`
- `region_key`
- `truth_status`
- `downstream_status`

#### B. read metadata

- `read_id`
- `read_name`
- `mapq`
- `hp`
- `alt_support`
- `strand`
- `is_tumor`

#### C. region baseline features

- `qual`
- `gq`
- `af`
- `ad_ref`
- `ad_alt`
- `VerificationClass`
- `PairwiseMedianDist`
- `AlleleDelta`
- `CramersV`
- `Quality_Score`
- `agreement_type`
- `class_shift`
- `hp_assign_rate`
- `allele_assign_rate`

#### D. read-level methylation summary

這部分目前 repo 內尚未直接輸出成 training table，建議 exporter 現場補齊：

- `num_cpg_observed`
- `methyl_mean`
- `methyl_std`
- `methyl_median`
- `methyl_high_fraction`
- `methyl_low_fraction`
- `methyl_na_fraction`

#### E. optional raw vector pointer

- `methyl_vector_json` 或 `methyl_vector_compact`
- `cpg_position_list`

### 8.3 第一版 supervision 建議

本輪推論：

1. 若做 read classifier prototype，最穩妥的是先預測 `read alt_support` 或 `region truth/downstream family` 的可分性，而不是直接宣稱在做 tumor/normal deconvolution。
2. 若要更接近 MethylBERT-style 問題，可以把 `paired normal` 放到第二輪，將 `is_tumor` 與 `paired normal methylation baseline` 一起納入。

也就是說：

- `Phase 1A`：先做 `tumor read within-region discrimination`
- `Phase 1B`：再做 `tumor vs normal read classification`

這是根據現有資產做出的保守切法，不是文獻直接定論。

---

## 9. 立即任務與執行順序

### 9.1 現在必做

1. 定義 `phase1 exporter contract`
2. 建立四個 baseline dataset 的 `region manifest`
3. 實作第一版 `per-read training table exporter`
4. 先在 `HCC1395 5kHz TO` 做小型 pilot
5. 再擴到 `HCC1395 5kHz paired`

### 9.2 接著要做

6. 把 `HCC1395_DORADO TO/paired` 納入 validation manifest
7. 補 `5kHz vs DORADO` harmonization 策略
8. 定義第一版 head-to-head baseline table

### 9.3 暫時後做

9. 直接做 Transformer / sequence model
10. 與 `Phase 2` 的 normal baseline / CN-aware 校正合併
11. 一開始就做四個 dataset joint training

### 9.4 本輪已落地的第一版 exporter

本輪已新增：

- [export_phase1_read_training_table.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_read_training_table.py)
- [build_phase1_training_manifest.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase1_training_manifest.py)
- [export_phase1_manifest_shard.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_manifest_shard.py)

用途：

1. 從 candidate-specific `rescue_joined_features.tsv` 選出 region
2. 對應回既有 `reads.tsv` 與 `methylation.csv`
3. 展平成 `read × region` 的 pilot training table
4. 由四個 baseline dataset 的 `significance_summary.csv` 先建立 `summary-first` 的完整 Phase 1 training manifest

目前已完成最小驗證：

- pilot 輸出根目錄：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_read_training_exporter_pilot](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_read_training_exporter_pilot)
- 本輪實際測試：
  - dataset：`HCC1395 5kHz TO`
  - `selected_regions=2`
  - `regions_found=2`
  - `missing_regions=0`
  - `read_rows_exported=77`

這表示 Phase 1 現在已從「只有 schema 設想」推進到「已有可重跑的 exporter 與實際 pilot rows」。

### 9.5 本輪已落地的 baseline manifest v1

為了讓四個 baseline dataset 都能納入同一個 Phase 1 底座，本輪另外補了 `summary-first manifest` 路線，而不是強依賴已漂移的 candidate-specific DORADO 路徑。

主輸出目錄：

- 完整 manifest：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1)
- 四象限 smoke test：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_smoke4](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_smoke4)

核心結果：

1. 完整 manifest v1 已成功建立，總共 `141,014` 個 baseline regions：
   - `HCC1395 5kHz paired`：`30,367`
   - `HCC1395 5kHz TO`：`40,096`
   - `HCC1395_DORADO paired`：`30,123`
   - `HCC1395_DORADO TO`：`40,428`
2. 這份 manifest 以 `significance_summary.csv` 為主軸，先把四象限都納進統一欄位：
   - `region_key`
   - `truth_status`
   - `VerificationClass`
   - `PairwiseMeanDist / PairwiseMedianDist`
   - `AlleleDelta`
   - `Quality_Score`
   - `Potential_LOH / Coverage_Category / LOH_Subtype`
3. 四象限 smoke test 已完成：
   - 每個 dataset 的 `TP / FP` 各取 `1` 個 region 做 read-level resolve
   - `8/8` 全部成功 resolve
   - `missing_regions = 0`
   - 共匯出 `1,701` 筆 `read × region` rows
4. manifest-driven shard exporter 已完成：
   - 可直接從 full manifest 做 selected-region read export
   - sample80 驗證已完成：
     - 每個 dataset 的 `TP / FP` 各取 `10` 個 region
     - `80/80` 全部成功 resolve
     - 共匯出 `14,165` 筆 reads
5. exporter schema 已補入任務導向欄位：
   - `dataset_role`
   - `harmonization_group`
   - `phase1a_task`
   - `phase1a_region_label`
   - `phase1a_read_label`
   - `phase1b_ready`
   - `phase1b_blocker`
6. `Phase 1A split manifest` 已建立：
   - [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_split_manifest_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_split_manifest_v1)
   - discovery：`70,463` regions
   - external validation：`70,551` regions

這代表目前已可把 `Phase 1` 切成兩層穩定資產：

1. `baseline full manifest`
2. `candidate-specific / selected-region read export`

### 9.6 本輪新增的研究與驗證分支

本輪除了完成主目標，也新增了幾個已確認值得繼續做的分支：

1. `summary-first manifest + on-demand resolve`
   - 已證明比「先全量掃 metadata 再建 manifest」更實際
2. `baseline full outputs` 與 `candidate-specific rescue outputs` 雙路並存
   - baseline 路線適合完整覆蓋
   - candidate-specific 路線適合 borderline / rescue 問題
3. `DORADO candidate-specific 路徑漂移` 被正式辨識為一個非 blocking 問題
   - 不再阻塞 Phase 1 baseline manifest
   - 但之後若要做 candidate-specific head-to-head，仍值得補做 provenance 回收
4. `Phase 1A / 1B` 的任務邊界已由資料驗證收斂
   - `Phase 1A` 可立即開始
   - `Phase 1B` 需先補 normal-read export layer

---

## 10. 本輪閱讀建議

1. [CURRENT_FOCUS.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
2. [20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md)
3. [20260307_純樣本甲基研究執行計畫_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md)
4. [20260315_甲基候選研究框架與觀察資產盤點_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_甲基候選研究框架與觀察資產盤點_01.md)
5. [20260311_TO_support特徵_readlevel_diagnostics_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md)
6. [20260325_Phase1_label_schema與harmonization規格_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1_label_schema與harmonization規格_01.md)
7. Knowledge [02_samples/HCC1395.md](/big8_disk/liaoyoyo2001/knowledge/02_samples/HCC1395.md)
8. Knowledge [06_workflows/benchmark_workflow.md](/big8_disk/liaoyoyo2001/knowledge/06_workflows/benchmark_workflow.md)

---

## 11. 本輪暫定結論

目前的任務與目標已經可以講得很具體：

1. **現在任務**：啟動 `Phase 1: ML read classification`
2. **現在目標**：先建立 `per-read training/export layer`，讓 read-level matrix 能被拿來做 head-to-head 方法比較
3. **現在不做**：直接宣稱模型優於 baseline，或跳去 normal baseline / CN-aware
4. **目前已完成的第一步**：`HCC1395 5kHz TO` 的 exporter pilot 已成功輸出 `77` 筆 `read × region` rows
5. **目前已完成的第二步**：四個 baseline dataset 的完整 manifest v1 已建立，且四象限 smoke resolve `8/8` 成功
6. **目前已完成的第三步**：`Phase 1A split manifest` 已建立，discovery / external validation pool 已固定
7. **目前已完成的第四步**：`Phase 1A` 的 head-to-head baseline table 與第一版 read-classifier benchmark 已建立
8. **下一個最合理的實作點**：
   - 做 `Phase 1A` 的 feature ablation 與 `to-pure` failure diagnosis
   - 規劃 normal-read export layer，作為 `Phase 1B` 的必要前置
   - 視需要補 `DORADO candidate-specific` provenance 回收
