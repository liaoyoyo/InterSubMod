<!--
建立時間: 2026-03-11 03:45
目標: 將 phase 2 finer interval / orthogonality 結果回接到 analysis-layer annotation，並驗證 annotation tier 是否能改善 caller-first rescue
處理範圍:
  - HCC1395 5kHz TO
  - HCC1395 5kHz paired
  - HCC1395_DORADO TO
  - HCC1395_DORADO paired
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase2_annotation_layer.py
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase2_paired_model_feature_analysis.py
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_layer_summary.md
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/phase2_summary.md
-->

# Phase 2：finer interval 回接 annotation layer 驗證

## 1. 破題結論

這輪已把 `phase 2` 的 finer interval 結果正式回接到 `analysis-layer annotation`，並在四個主要 dataset 上完成實測。結論很清楚：**annotation layer 是有價值的，但目前最合理的定位仍是「分層註記、排序與 QC」，不是直接取代 caller-first 的新主規則。**

更具體地說：

1. `HCC1395 5kHz TO` 與 `HCC1395_DORADO TO` 的 `Quality_Score`、`PairwiseMedianDist` top-bin 確實都富集 `TP`，適合升級成第二層 `support annotation`。
2. `hp_assign_rate` 即使在某些區間有正向 enrichment，也仍較適合保留為 `phase/QC annotation`，不應升級成 truth-level keep 規則。
3. 將這些 annotation 回接成 tier/score 後，**沒有一個 policy 超過原本最佳 caller-first gate**：
   - `HCC1395 5kHz TO`：最佳仍是 `gq>=10`
   - `HCC1395 5kHz paired`：最佳仍是 `gq>=15`
   - `HCC1395_DORADO TO`：最佳仍是 `gq>=10`
   - `HCC1395_DORADO paired`：最佳反而是更嚴格的 `gq>=20`
4. 因此這輪最合理的流程調整不是改 hard filter，而是新增 annotation 欄位，把已驗證有效的 support / QC / artifact 訊號明確寫出來，供後續判讀與排序使用。

---

## 2. 本輪目標

本輪要回答 3 個問題：

1. `phase 2` 的 finer interval 結果，能不能被制度化成一層穩定、可重複使用的 annotation。
2. 這些 annotation 若組成 tier 或 score，是否能真的比既有 `caller-first` 更好。
3. 哪些 annotation 應升級成 `support`，哪些應保留為 `QC` 或 `artifact triage`。

---

## 3. 方法與資料來源

### 3.1 執行腳本

- [build_phase2_annotation_layer.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase2_annotation_layer.py)

### 3.2 輸入資料

1. phase 2 finer interval 結果：
   - [fine_feature_interval_top_bins.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/fine_feature_interval_top_bins.tsv)
   - [feature_orthogonality_top_pairs.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_paired_model_feature_analysis/feature_orthogonality_top_pairs.tsv)
2. 四個 dataset 的 candidate-specific rescue join：
   - [HCC1395 5kHz TO rescue_joined_features.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv)
   - [HCC1395 5kHz paired rescue_joined_features.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260308_borderline_rescue_analysis/HCC1395_paired/eval/rescue_joined_features.tsv)
   - [HCC1395_DORADO paired rescue_joined_features.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_joined_features.tsv)
   - [HCC1395_DORADO TO rescue_joined_features.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv)

### 3.3 主要輸出

- [annotation layer output 目錄](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer)
- [annotation_layer_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_layer_summary.md)
- [annotation_layer_config.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_layer_config.tsv)
- [annotation_presence_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_presence_summary.tsv)
- [annotation_policy_evaluation.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_policy_evaluation.tsv)
- [annotation_overlap_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_overlap_summary.tsv)
- [annotation_best_policy.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_best_policy.tsv)

### 3.4 Annotation 設計原則

本輪固定把 annotation 分成三層：

1. `support`
   - 允許作為第二層 support 訊號
2. `qc`
   - 只表示 phase / coverage / completeness 狀態
3. `artifact`
   - 只作後段 triage，不直接作 TP rescue 主規則

---

## 4. Annotation layer 設計

### 4.1 Dataset-specific 設定

| Dataset | Caller primary | Caller strict | Support annotation | QC annotation | 說明 |
| --- | --- | --- | --- | --- | --- |
| `HCC1395 5kHz TO` | `gq>=10` | `gq>=15` | `Quality top-bin`、`Pairwise top-bin`、`agreement_positive`、`Strong/Subclone` | `hp_assign>=0.95`、`hp_assign>=0.99`、`hp_assign top-bin` | TO 主線，support 比較完整 |
| `HCC1395_DORADO TO` | `gq>=10` | `gq>=15` | `Quality top-bin`、`Pairwise top-bin`、`Strong/Subclone` | `hp_assign>=0.95`、`hp_assign>=0.99`、`hp_assign top-bin` | 不把 `agreement_positive` 升級成 support |
| `HCC1395 5kHz paired` | `gq>=15` | `gq>=20` | 暫無 | `Quality top-bin`、`Pairwise top-bin`、`hp_assign top-bin` | paired 目前只保留 annotation/QC |
| `HCC1395_DORADO paired` | `gq>=15` | `gq>=20` | `Strong/Subclone`（弱 support） | `Quality top-bin`、`Pairwise top-bin`、`hp_assign top-bin` | paired 仍以 caller-first 為主 |

對應設定檔可直接查：
- [annotation_layer_config.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_layer_config.tsv)

### 4.2 top-bin 來源

這一輪沒有重新發明閾值，而是直接沿用 phase 2 驗證出的最佳區間：

| Dataset | `Quality_Score` top-bin | `PairwiseMedianDist` top-bin | `hp_assign_rate` top-bin |
| --- | --- | --- | --- |
| `HCC1395 5kHz TO` | `[55,60)` | `[0.18,0.20)` | `[0.5,0.7)` |
| `HCC1395_DORADO TO` | `[55,60)` | `[0.12,0.15)` | `[0.85,0.9)` |
| `HCC1395 5kHz paired` | `[80,inf)` | `[0.3,inf)` | `[0,0.5)` |
| `HCC1395_DORADO paired` | `[80,inf)` | `[0.1,0.12)` | `[0.99,1.01)` |

這張表的意義是：annotation layer 這次不是憑印象寫規則，而是直接把已驗證的區間回接到資料層。

---

## 5. Annotation enrichment 結果

### 5.1 `HCC1395 5kHz TO`

| Annotation | 角色 | TP | FP | `TP/FP` | `enrichment_ratio` | 解讀 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `quality_topbin` | support | 33 | 5 | 6.60 | 2.0827 | 小而乾淨，適合 support annotation |
| `pairwise_topbin` | support | 60 | 12 | 5.00 | 1.5778 | 有訊號，但不是單獨 keep 規則 |
| `agreement_positive` | support | 213 | 62 | 3.44 | 1.0841 | 較廣的 support |
| `strong_subclone` | support | 208 | 60 | 3.47 | 1.0939 | 與 agreement 高重疊，但仍可留作弱 support |
| `hp_assign_topbin` | qc | 35 | 8 | 4.38 | 1.3806 | 有 enrichment，但仍應視為 QC |
| `artifact_lowvaf_highad_lowcv` | artifact | 2 | 1 | 2.00 | 0.6311 | 在這個 candidate-specific pool 中數量太少，不適合直接主導 |

這組數據表示：`5kHz TO` 的 `Quality` 與 `Pairwise` 確實值得回接成 support annotation，但 `hp_assign` 仍不該誤認成 truth support。

### 5.2 `HCC1395_DORADO TO`

| Annotation | 角色 | TP | FP | `TP/FP` | `enrichment_ratio` | 解讀 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `quality_topbin` | support | 20 | 5 | 4.00 | 1.5223 | support 成立 |
| `pairwise_topbin` | support | 31 | 5 | 6.20 | 2.3595 | dataset-aware support 成立，但方向與 5kHz 不同 |
| `strong_subclone` | support | 88 | 30 | 2.93 | 1.1163 | 可作弱 support |
| `agreement_positive` | annotation | 59 | 27 | 2.19 | 0.8316 | 不建議升成 support |
| `hp_assign_topbin` | qc | 10 | 2 | 5.00 | 1.9028 | 高 enrichment，但仍應保留為 QC，不宜直接 keep |

這裡再次證明：`PairwiseMedianDist` 的方向明顯 dataset-dependent，但 `Quality_Score` 作為 support annotation 比較穩。

### 5.3 paired

`paired` 的結果更保守：

| Dataset | 最值得保留的 paired annotation | 解讀 |
| --- | --- | --- |
| `HCC1395 5kHz paired` | `Quality top-bin`, `Pairwise top-bin`, `artifact_lowvaf_highad_lowcv` | 目前都只適合 annotation / artifact，不支持升成 support |
| `HCC1395_DORADO paired` | `hp_assign top-bin`, `Quality top-bin`, `Strong/Subclone` | 有弱 support / QC 訊號，但 caller-first 仍遠強於它們 |

---

## 6. Annotation policy 實測

### 6.1 四個 dataset 的最佳 policy

| Dataset | `best_policy_by_baseline` | `delta_f1_vs_baseline` | `best_policy_vs_primary_gate` | `delta_f1_vs_primary_gate` |
| --- | --- | ---: | --- | ---: |
| `HCC1395 5kHz TO` | `caller_primary_only` | `+0.006754` | `caller_primary_only` | `0.000000` |
| `HCC1395 5kHz paired` | `caller_primary_only` | `+0.000873` | `caller_primary_only` | `0.000000` |
| `HCC1395_DORADO paired` | `caller_strict_only` | `+0.000635` | `caller_strict_only` | `+0.000405` |
| `HCC1395_DORADO TO` | `caller_primary_only` | `+0.000540` | `caller_primary_only` | `0.000000` |

來源：
- [annotation_best_policy.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotation_best_policy.tsv)

### 6.2 `HCC1395 5kHz TO`：annotation 有價值，但仍輸給 caller-primary

| Policy | TP rescued | FP reintroduced | F1 | `delta_f1_vs_primary_gate` |
| --- | ---: | ---: | ---: | ---: |
| `caller_primary_only (gq>=10)` | 486 | 117 | 0.719451 | `0.000000` |
| `strict_or_support_any` | 436 | 95 | 0.718850 | `-0.000601` |
| `strict_or_support_primary` | 423 | 88 | 0.718705 | `-0.000746` |
| `strict_or_support_score_ge2` | 384 | 78 | 0.718171 | `-0.001280` |

意義：

1. annotation tier 確實讓 FP 更少，但同步也失去太多 TP。
2. 所以在 `5kHz TO`，annotation layer 比較適合當排序、提示和解讀層，而不是直接縮成新 hard keep。

### 6.3 `HCC1395_DORADO TO`：同樣是 support / annotation，不是主規則

| Policy | TP rescued | FP reintroduced | F1 | `delta_f1_vs_primary_gate` |
| --- | ---: | ---: | ---: | ---: |
| `caller_primary_only (gq>=10)` | 40 | 11 | 0.723113 | `0.000000` |
| `strict_or_support_any` | 37 | 9 | 0.723083 | `-0.000030` |
| `strict_or_support_primary` | 35 | 9 | 0.723051 | `-0.000062` |

意義：

1. annotation layer 幾乎追平 caller-primary，但還是沒有超過。
2. 這表示跨兩個 TO dataset 都支持同一個高層結論：
   - 第一層：caller
   - 第二層：support annotation

### 6.4 paired：更支持 caller-first

`paired` 的結果更明確：

- `HCC1395 5kHz paired`：`gq>=15` 已是最佳，annotation 不提供額外增量。
- `HCC1395_DORADO paired`：更嚴格的 `gq>=20` 甚至優於 `gq>=15`；annotation 只能帶來非常小的補充，且仍未超過 `gq>=20`。

---

## 7. overlap / orthogonality 補充判讀

### 7.1 `HCC1395 5kHz TO`

較值得記住的 overlap：

| Feature A | Feature B | Jaccard | 解讀 |
| --- | --- | ---: | --- |
| `agreement_positive` | `strong_subclone` | `0.280660` | 兩者高度相關，但不是完全等價 |
| `agreement_positive` | `hp_assign_high_095` | `0.265432` | label support 與高 HP completeness 有部分交集 |
| `quality_topbin` | `pairwise_topbin` | `0.000000` | 幾乎完全正交，代表 quality 與 pairwise 在 5kHz TO 確實抓不同子群 |

這個結果支持把 `Quality` 與 `Pairwise` 都保留在 annotation layer，因為它們不是重複訊號。

### 7.2 `HCC1395_DORADO TO`

`DORADO TO` 的 annotation overlap 較少，但 `quality_topbin` 與 `pairwise_topbin` 同樣不是完全重疊。這支持它們作為 dataset-aware support annotation，而不是硬合成單一 feature。

---

## 8. 這輪最合理的流程調整

### 8.1 應做

1. 在 analysis-layer 增加 annotation 欄位：
   - `CallerPrimaryTier`
   - `CallerStrictTier`
   - `MethylSupport_QualityTopBin`
   - `MethylSupport_PairwiseTopBin`
   - `MethylSupport_AgreementPositive`
   - `PhaseQC_HPAssignHigh`
   - `Artifact_LowVAF_HighAlleleDelta_LowCV`
2. 在報告與 TSV 層明確區分：
   - `support`
   - `qc`
   - `artifact`
3. 將 `Quality_Score`、dataset-aware `PairwiseMedianDist` 回接到排序與人工判讀層。

### 8.2 先不要做

1. 不要把 `Quality top-bin` 直接當 keep 規則。
2. 不要把 `Pairwise top-bin` 直接做跨樣本全域硬規則。
3. 不要把 `hp_assign_rate` 誤升級成 truth-level rescue 主特徵。
4. 不要把這輪 annotation tier 直接改寫成 pipeline 預設 hard filter。

---

## 9. 驗證狀態

1. `python3 -m py_compile [build_phase2_annotation_layer.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_phase2_annotation_layer.py)` 通過。
2. annotation builder 已實際跑完四個 dataset。
3. 已產出四份 annotated candidates：
   - [hcc1395_5khz_to_annotated_candidates.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotated_candidates/hcc1395_5khz_to_annotated_candidates.tsv)
   - [hcc1395_5khz_paired_annotated_candidates.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotated_candidates/hcc1395_5khz_paired_annotated_candidates.tsv)
   - [hcc1395_dorado_to_annotated_candidates.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotated_candidates/hcc1395_dorado_to_annotated_candidates.tsv)
   - [hcc1395_dorado_paired_annotated_candidates.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_phase2_annotation_layer/annotated_candidates/hcc1395_dorado_paired_annotated_candidates.tsv)

---

## 10. 結論

這輪已經把 `phase 2` 的 finer interval 結論正式回接到 annotation layer，而且結果足夠清楚：

1. `Quality_Score` 值得升級成 `support annotation`。
2. `PairwiseMedianDist` 值得升級成 `dataset-aware support annotation`，但不能全域化。
3. `hp_assign_rate` 應明確降階為 `phase/QC annotation`。
4. `artifact_lowvaf_highad_lowcv` 仍保留在 triage 支線。
5. 目前所有 annotation tier 測試都沒有超過最佳 caller-first gate，因此下個階段最合理的工作是：
   - 先把 annotation 層接到輸出與報告
   - 再觀察它在排序、人工審閱與 cross-dataset 解讀上的價值
   - 而不是直接改 pipeline hard rule

---

## 11. 下一步

1. 把這輪 annotation 欄位接到整合矩陣與週報主報告。
2. 若要再往前走，優先做：
   - annotation score 的 ranking 實驗
   - paired coverage ceiling 的補強
   - 把 annotation 欄位接到 candidate diagnostics / dashboard
3. 若未來要正式改 pipeline，優先做 `annotation export`，不是 `hard keep / hard veto`。
