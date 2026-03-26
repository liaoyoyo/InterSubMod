<!--
建立時間: 2026-03-25 18:55
目標: 根據 2026-03-25 Phase 1 exporter / manifest / shard 驗證結果，正式收斂 Phase 1A/1B label schema 與 5kHz/DORADO harmonization 規格
處理範圍:
  - Phase 1A 立即可執行任務定義
  - Phase 1B 目前阻塞條件
  - 5kHz / DORADO / paired / TO harmonization 邊界
  - 下一步任務清單
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1_ML_read_classification研究啟動與資料缺口盤點_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase1_training_manifest.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_manifest_shard.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_read_training_table.py
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_smoke4
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80
-->

# Phase 1 label schema 與 harmonization 規格

## 1. 先講結論

截至 `2026-03-25`，Phase 1 已可正式拆成兩條線：

1. `Phase 1A` 可以立即開始：
   - 任務定義：`within-tumor alt-support read classification`
   - 訓練單位：`一列 = 一個 tumor read 在一個 region 下的觀測`
   - 主要 read label：`ALT / REF`
   - 主要 region context：`truth_status + VerificationClass + Quality_Score + PairwiseMedianDist`
2. `Phase 1B` 目前不能直接開始：
   - 原因不是 exporter 壞掉，而是當前 read exports 裡沒有 normal reads
   - sampled shard 觀察到 `is_tumor = 1` 對所有 `14,165` 筆 read 都成立
   - 這表示「tumor vs normal read classification」需要另外補 normal-read export layer
3. `5kHz` 與 `DORADO` 不能先當成同一母體直接 joint train：
   - 即使是同 biological sample，read-level `num_cpg_observed / methyl_mean / methyl_std` 分布仍有明顯 platform + mode shift
   - 第一版必須保留 `harmonization_group = platform|mode`

因此，本輪真正完成的不是模型，而是 **Phase 1 的資料任務定義已收斂，且 exporter schema 已直接反映這個結論**。

---

## 2. 本輪證據底座

### 2.1 region-level manifest 已完成

完整 manifest：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1)

核心數字：

- 四個 baseline dataset 共 `141,014` 個 regions
- `TP = 116,974`
- `FP = 24,040`

分 dataset：

- `HCC1395 5kHz paired`: `30,367`
- `HCC1395 5kHz TO`: `40,096`
- `HCC1395_DORADO paired`: `30,123`
- `HCC1395_DORADO TO`: `40,428`

### 2.2 manifest-driven shard export 已完成

四象限 smoke：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_smoke4](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_smoke4)
- `8/8` selected regions 成功 resolve
- `missing = 0`
- `read_rows = 1,701`

sampled shard：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80)
- 每個 dataset 的 `TP / FP` 各取 `10` 個 region
- 共 `80` 個 region
- `80/80` 成功 resolve
- `read_rows = 14,165`

### 2.3 exporter schema 已落地 Phase 1 任務欄位

目前 read-level 輸出已新增：

- `dataset_role`
- `harmonization_group`
- `phase1a_task`
- `phase1a_region_label`
- `phase1a_read_label`
- `phase1b_ready`
- `phase1b_blocker`

這代表後續不需要再靠人工讀文件判斷資料屬性。

---

## 3. Phase 1A 正式定義

### 3.1 任務

`Phase 1A = within-tumor alt-support read classification`

### 3.2 為什麼先做這條

1. 目前所有已打通的 read exports 都穩定提供：
   - `alt_support`
   - `hp`
   - `mapq`
   - `strand`
   - per-read methylation summary
   - region-level baseline context
2. `phase1a_read_label` 在 sampled shard 中已可直接使用：
   - `REF = 8,071`
   - `ALT = 6,094`
3. 這條任務不需要假設目前已有 normal reads，也不需要先引入 Phase 2 的 normal baseline 建模。

### 3.3 第一版 supervision

read-level 主 label：

- `phase1a_read_label`
  - `ALT`
  - `REF`
  - 未來若來源出現未知值則記為 `UNKNOWN`

region-level context：

- `phase1a_region_label = truth_status`
- `VerificationClass`
- `Quality_Score`
- `PairwiseMedianDist`
- `PassedGating`

資料角色：

- `dataset_role = discovery`
  - `HCC1395 5kHz paired`
  - `HCC1395 5kHz TO`
- `dataset_role = validation`
  - `HCC1395_DORADO paired`
  - `HCC1395_DORADO TO`

---

## 4. Phase 1B 目前不能直接開始

### 4.1 阻塞條件

`Phase 1B` 原本暫定為 `tumor vs normal read classification`。

但 sampled shard 已明確顯示：

- `phase1b_ready = False` 對全部 `14,165` 筆 reads 都成立
- `is_tumor = 1` 對全部 `14,165` 筆 reads 都成立

這不是單純的匯出欄位缺漏，而是 **目前輸出的 read universe 本身就是 tumor-only**。

### 4.2 來源驗證

paired canonical output 的原始 `reads.tsv` 範例也符合這個現象：

- [reads.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp/chr7/chr7_153978134/chr7_153973134_153983134/reads/reads.tsv)

至少在目前抽查的 paired region 中，`is_tumor` 只有 `1`。

更進一步的 code-level 檢查也支持同一個結論：

1. [RegionProcessor.cpp](/big7_disk/liaoyoyo2001/InterSubMod/src/core/RegionProcessor.cpp) 會初始化：
   - `tumor_reader`
   - `normal_reader`
2. 但 `process_single_region(...)` 目前只接收 `tumor_reader`
3. `normal_reader` 在 [RegionProcessor.cpp](/big7_disk/liaoyoyo2001/InterSubMod/src/core/RegionProcessor.cpp) 中沒有後續使用點
4. 對 paired canonical outputs 執行 `ALT|REF|UNKNOWN\t0\t` 搜尋也沒有命中

因此目前 Phase 1B 的阻塞點已不是單純資料整理問題，而是 **region export pipeline 尚未把 normal reads 納入 read-level 輸出**。

### 4.3 對研究規劃的意義

因此 Phase 1B 不應再被寫成「下一步直接接著做」。

更準確的說法是：

1. `Phase 1A` 先做 tumor-read 內部 discrimination
2. `Phase 1B` 需要先補：
   - normal-read export layer
   - tumor/normal 對齊規則
   - normal baseline join 規格

換句話說，`Phase 1B` 與 `Phase 2 normal baseline` 已開始有依賴關係，但還不需要立刻合併成同一階段。

---

## 5. harmonization 規格

### 5.1 不可直接 joint pool 的證據

sampled shard 顯示不同 `harmonization_group` 的 read-level分布有可見落差。

以 `num_cpg_observed` 與 `methyl_mean` 為例：

- `ONT_5kHz|paired-pure, FP`
  - `num_cpg_observed mean = 81.8548`
  - `methyl_mean mean = 0.7346`
- `ONT_Dorado|paired-pure, FP`
  - `num_cpg_observed mean = 58.6169`
  - `methyl_mean mean = 0.6690`
- `ONT_5kHz|to-pure, TP`
  - `num_cpg_observed mean = 92.3815`
  - `methyl_mean mean = 0.4572`
- `ONT_Dorado|to-pure, TP`
  - `num_cpg_observed mean = 99.3457`
  - `methyl_mean mean = 0.4976`

這些差異足以說明：

1. `5kHz vs DORADO` 有 platform shift
2. `paired vs TO` 也有 mode shift
3. 兩個 shift 同時存在時，不應先做單一 joint distribution 假設

### 5.2 第一版規則

第一版 harmonization 規格收斂為：

1. 永遠保留 `harmonization_group = platform|mode`
2. 第一版 benchmark 先做：
   - `5kHz discovery`
   - `DORADO external validation`
3. 若要 joint train，至少要同時滿足：
   - 顯式帶入 `platform`
   - 顯式帶入 `mode`
   - 對 read-level methylation summary 做 group-aware normalization
4. 在還沒做 normalization 以前：
   - 允許 cross-platform compare
   - 不允許把 raw read features 當同分布直接混訓

### 5.3 split 建議

第一版先採最保守切法：

1. discovery:
   - `HCC1395 5kHz paired`
   - `HCC1395 5kHz TO`
2. external validation:
   - `HCC1395_DORADO paired`
   - `HCC1395_DORADO TO`
3. 報告時至少分三層呈現：
   - overall
   - by `harmonization_group`
   - by `truth_status`

---

## 6. 已完成清單

1. candidate-specific pilot exporter 已完成並重跑：
   - [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_read_training_exporter_pilot](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_read_training_exporter_pilot)
2. baseline full manifest v1 已完成：
   - [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_v1)
3. manifest-driven shard exporter 已完成：
   - [export_phase1_manifest_shard.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_manifest_shard.py)
4. smoke4 與 sample80 均驗證成功：
   - smoke4：`8/8` resolve，`1,701` reads
   - sample80：`80/80` resolve，`14,165` reads
5. exporter schema 已內建 `Phase 1A/1B` 任務欄位。
6. `Phase 1A split manifest` 已完成：
   - [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_split_manifest_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_split_manifest_v1)
   - discovery regions：`70,463`
   - external validation regions：`70,551`
7. `Phase 1A` 的 head-to-head baseline 與 read-classifier benchmark 已完成：
   - baseline table：
     [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_head_to_head_baseline_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_head_to_head_baseline_v1)
   - benchmark round：
     [20260325_Phase1A_read_classifier_benchmark_round1_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1A_read_classifier_benchmark_round1_01.md)
   - sample200 最穩定結果：
     - discovery holdout `F1=0.7882`
     - external validation `F1=0.8908`

---

## 7. 下一步任務

### 7.1 可以直接繼續做

1. 做 `Phase 1A` 的 feature ablation
2. 補 `Phase 1A` 的訓練/驗證 protocol
3. 擴大 sampled shard 或改成分批全量 read export
4. 補 `Phase 1A` 的失敗診斷輸出：
   - by `harmonization_group`
   - by `truth_status`
   - by `VerificationClass`

### 7.2 需要另外補資料層才可做

5. `Phase 1B tumor-vs-normal`
6. normal methylation reference join
7. `5kHz / DORADO` joint-train normalization

### 7.3 已定位但尚未動手的工程缺口

8. paired mode 的 `normal_reader` 目前未進入 `process_single_region`
9. 若要啟動 `Phase 1B`，需要在「改主流程納入 normal reads」與「新增 sidecar normal-read exporter」之間選擇一條實作路線

---

## 8. 本輪結案口徑

本輪真正完成的不是模型訓練，而是：

1. `Phase 1` 的 baseline manifest 已建立
2. shard export 已從 full manifest 成功補 resolve
3. `Phase 1A` 與 `Phase 1B` 的邊界已經不是概念，而是被資料驗證過的可執行規格
4. 後續若再寫「直接做 tumor-vs-normal read classification」，應視為與目前證據不一致
