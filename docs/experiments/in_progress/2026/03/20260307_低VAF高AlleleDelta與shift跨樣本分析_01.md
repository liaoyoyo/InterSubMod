<!--
建立時間: 2026-03-07 19:10
目標: 延續 pure paired 主線，分析低 VAF + 高 AlleleDelta 候選特徵、Weak->Strong / Noise->Strong shift、read-level snapshot 與 tumor-only 後續接線
處理範圍:
  - HCC1395 5kHz / HCC1395_DORADO 新版正式 rerun
  - 其他 5 個 pure paired 既有 run（COLO829, H1437, H2009, HCC1937, HCC1954）
  - 候選規則跨樣本比較
  - shift TP/FP 分層
  - samtools region snapshot
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260307_純樣本round執行與進度報告_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_candidate_rule_analysis/
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_shift_analysis_all_pure/
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_other_pure_existing_rounds/
-->

# 低 VAF + 高 AlleleDelta 與 shift 跨樣本分析

## 1. 本次目標

延續 `pure paired` 主線，確認下列問題：

1. `低 VAF + 偏高 AlleleDelta` 是否值得加入原始流程判斷
2. 這類特徵是否只對 `HCC1395 5kHz` 有效，或在其他樣本也有用
3. `Weak->Strong / Noise->Strong` 是否偏向 TP 或 FP
4. `label-first` 的最終標籤是否合理，是否需要更細的分類
5. tumor-only 下一步可接哪些資料，而不是停留在概念

本次背景知識確認：

1. 根據 Knowledge `02_samples/HCC1395.md`，`HCC1395 5kHz` 為含 `MM/ML` 的主要甲基主線資料，`HCC1395_DORADO` 為同樣本不同定序流程交叉驗證。
2. 根據 Knowledge `06_workflows/benchmark_workflow.md`，benchmark 必須固定 truth scope、固定 `PASS` 比較口徑，且報告 `TP/FP/FN/truth_total`。
3. 根據 Knowledge `05_tools/InterSubMod.md`，InterSubMod 建議搭配 haplotagged BAM 進行更精確的 read-level 甲基分析。

## 2. 本次新增腳本與自動化

### 2.1 新增腳本

1. [analyze_candidate_rules.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_candidate_rules.py)
   - 跨樣本比較候選規則的 `TP/FP/F1 delta`
   - 自動 fallback 到前一版 `longphase_s` TP/FP VCF，支援 `--skip-longphase` run
2. [analyze_shift_patterns.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_shift_patterns.py)
   - 將 `Weak->Strong / Noise->Strong / Weak->Subclone / Noise->Subclone` 依 `TP/FP` 分層
   - 輸出 `shift_tp_fp_summary.tsv` 與 `shift_top_regions.tsv`
3. [region_samtools_snapshot.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/region_samtools_snapshot.sh)
   - 對指定 region 輸出 `depth.tsv`、`mpileup.txt`、`reads.sam`、`hp_tag_counts.tsv`
   - 用於快速做 read-level 驗證

### 2.2 本次 round / analysis 輸出

1. [candidate_rule_summary.md](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_candidate_rule_analysis/candidate_rule_summary.md)
2. [candidate_rule_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_candidate_rule_analysis/candidate_rule_comparison.tsv)
3. [best_rule_by_sample.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_candidate_rule_analysis/best_rule_by_sample.tsv)
4. [shift_pattern_summary.md](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_shift_analysis_all_pure/shift_pattern_summary.md)
5. [shift_tp_fp_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_shift_analysis_all_pure/shift_tp_fp_summary.tsv)
6. [shift_top_regions.tsv](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_shift_analysis_all_pure/shift_top_regions.tsv)
7. [other pure round README](/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/20260307_other_pure_existing_rounds/README.md)

## 3. 候選規則跨樣本結果

候選規則主要比較：

```text
1. 舊 combo：
   (QUAL < 0.75) OR (AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24)

2. 目前 core：
   AlleleDelta > 0.15 AND VAF < 0.15
```

### 3.1 每個樣本最佳規則

| 樣本 | 平台 | baseline F1 | 最佳規則 | best F1 | delta |
|---|---|---:|---|---:|---:|
| HCC1395 | ONT_5kHz | 0.8522 | `current_core_ad015_vaf015` | 0.853190 | +0.000990 |
| HCC1395_DORADO | ONT_Dorado | 0.8592 | `baseline` | 0.859200 | +0.000000 |
| COLO829 | ONT | 0.8921 | `baseline` | 0.892100 | +0.000000 |
| H1437 | ONT | 0.8562 | `baseline` | 0.856200 | +0.000000 |
| H2009 | ONT | 0.8816 | `baseline` | 0.881600 | +0.000000 |
| HCC1937 | ONT | 0.3382 | `old_written_rule_raw_qual` | 0.338219 | +0.000019 |
| HCC1954 | ONT | 0.8048 | `baseline` | 0.804800 | +0.000000 |

### 3.2 主要結論

1. `HCC1395 5kHz` 是唯一明顯正增益樣本
   - `current_core_ad015_vaf015`：`TP removed=2`，`FP removed=83`，`F1 +0.000990`
2. `HCC1395_DORADO` 與其他大多數樣本都不是正增益
   - `HCC1395_DORADO`：最佳仍是 `baseline`
   - `COLO829 / H1437 / H2009 / HCC1954`：所有候選規則都不比 baseline 好
3. `HCC1937` 只有極小的 `+0.000019`，可視為接近無效
4. `CramersV < 0.05` 在 `HCC1395 5kHz` 上幾乎沒有額外篩選力
   - `ad015_vaf015`
   - `ad015_cv005_vaf015`
   兩者在 5kHz 上給出相同最佳結果

### 3.3 對原始流程的判斷

目前不建議把這組規則直接納入全域原始流程。

原因：

1. 它只在 `HCC1395 5kHz` 顯示穩定正增益
2. 在 `DORADO` 與其他樣本上多數是負增益或接近無效
3. 這更像「平台/樣本特化的診斷規則」，而不是跨樣本通用過濾規則

目前較合理的定位：

1. 先作為 `5kHz` 特化的候選 post-filter
2. 或作為 `ArtifactSuspect` 標記訊號，而不是直接刪除

## 4. `Weak->Strong / Noise->Strong` 的 TP/FP 分層

### 4.1 `Weak->Strong`

FP scope rate（`FP 中有多少比例變成 Weak->Strong`）：

| 樣本 | FP rate | FP count / total |
|---|---:|---:|
| HCC1395 | 0.274322 | 172 / 627 |
| HCC1954 | 0.241379 | 7 / 29 |
| COLO829 | 0.157754 | 354 / 2244 |
| HCC1395_DORADO | 0.116667 | 28 / 240 |
| H2009 | 0.093023 | 8 / 86 |
| HCC1937 | 0.061538 | 12 / 195 |
| H1437 | 0.000000 | 0 / 8 |

解讀：

1. `Weak->Strong` 在所有樣本都主要還是 TP 為主
2. 但 `HCC1395 5kHz` 的 FP rate 最高，說明它的升級型 `Strong` 確實更 noisy
3. `HCC1954` 的比例也高，但母數只有 `29` 個 FP，不能和 5kHz 直接等量觀察

### 4.2 `Noise->Strong`

FP scope rate：

| 樣本 | FP rate | FP count / total |
|---|---:|---:|
| HCC1395 | 0.081340 | 51 / 627 |
| HCC1395_DORADO | 0.058333 | 14 / 240 |
| COLO829 | 0.000000 | 0 / 2244 |
| H1437 | 0.000000 | 0 / 8 |
| H2009 | 0.000000 | 0 / 86 |
| HCC1937 | 0.000000 | 0 / 195 |
| HCC1954 | 0.000000 | 0 / 29 |

解讀：

1. `Noise->Strong` 幾乎只出現在 `HCC1395` 的兩種平台
2. `5kHz` 比 `DORADO` 更高
3. 這是本輪最值得注意的「平台/資料集特異」現象之一

### 4.3 對標籤合理性的判斷

目前 `label-first` 並不是無效，但 `Strong` 這個終態可能還太粗。

因為：

1. `Weak->Strong` 大部分是 TP，但在 `5kHz` 有相對高比例的 FP
2. `Noise->Strong` 只有 `HCC1395` 系列大量出現，暗示 `Strong` 內部可能混有 artifact-prone 子群
3. 因此後續值得考慮在 `Strong` 內再切一層，例如：
   - `Strong-Trusted`
   - `Strong-SuspectLowVAF`
   - `Strong-NoiseProne`

## 5. 候選規則與升級型 shift 的交叉觀察

在 `HCC1395 5kHz` 中，`current_core_ad015_vaf015` 觸發集和 `Weak->Strong / Noise->Strong` 的重疊很有限：

1. TP 的升級型 shift 與規則觸發集幾乎沒有重疊
2. FP 的 `Weak->Strong / Noise->Strong` 一共 `223` 筆，其中只和規則觸發集重疊 `6` 筆
3. 且這 `6` 筆都屬於 `Weak->Strong`

解讀：

1. 這個規則抓到的不是「所有會被升級成 Strong 的 noisy 位點」
2. 它抓到的是更窄的一種 5kHz artifact 類型
3. 因此不應把 `Weak->Strong` 全部視為可過濾，也不應把規則當成 `label_upgrade` 的替代

## 6. `samtools` read-level 初步驗證

### 6.1 代表性 FP

使用 [region_samtools_snapshot.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/region_samtools_snapshot.sh) 對 `HCC1395 5kHz` 的代表性 FP `chr8:94170573:G:T` 匯出 snapshot：

1. [FP snapshot dir](/bip8_disk/liaoyoyo2001/InterSubMod_out/output/s-pure/HCC1395/20260307/samtools_snapshots/chr8_94170573_G_T)
2. `read_count=83`
3. 前 `60` 條 reads 的 HP tag 分布：
   - `1 = 28`
   - `1-1 = 11`
   - `2-1 = 4`
   - `3 = 4`
   - `NA = 13`
4. `mpileup` 顯示該位點深度 `66`，存在混合 `G/g/T`

### 6.2 代表性 TP

同樣對規則誤傷的 TP `chr6:61876832:C:T` 匯出 snapshot：

1. [TP snapshot dir](/bip8_disk/liaoyoyo2001/InterSubMod_out/output/s-pure/HCC1395/20260307/samtools_snapshots/chr6_61876832_C_T)
2. `read_count=137`
3. 前 `60` 條 reads 的 HP tag 分布：
   - `1 = 29`
   - `1-1 = 2`
   - `2 = 20`
   - `2-1 = 9`
4. `mpileup` 顯示該位點深度 `132`，混合 `C/c/T`

### 6.3 初步解讀

1. 這兩個 snapshot 都顯示 read-level 支持不是單一 outlier 就能解釋
2. 但 FP 個案的 HP 結構較偏單側與未分派混雜
3. TP 個案在兩個主 haplotype 都有較完整覆蓋
4. 目前這只是一個小型 read-level 檢查，不能單靠 2 個 region 定論，但已證明 `samtools snapshot` 值得納入後續驗證流程

## 7. tumor-only 目前可接資料

目前雖然還沒有 standardized 的 tumor-only InterSubMod round，但資料端已有可接資源：

1. `/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0`
2. `/big8_disk/data/HCC1395/ONT/ClairS_TO_ss_v0_3_0`
3. `/big8_disk/data/HCC1395/ONT/DeepSomatic_TO_v1_8_0`
4. `/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0`
5. `/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_ss_v0_3_0`
6. `/big8_disk/data/HCC1395/ONT_Dorado/DeepSomatic_TO_v1_8_0`

目前缺的是：

1. tumor-only 專用的 InterSubMod benchmark bundle
2. 可直接和 paired pure 同口徑比較的 `label_cluster_agreement.tsv`
3. 可直接驗證 `低 VAF + 高 AlleleDelta` 是否也能抓 tumor-only artifact 的標準流程

因此 tumor-only 現在應列為下一輪可立即接線的方向，但不應和本輪 pure paired 主結論混在一起。

## 8. 本輪結論

1. `低 VAF + 偏高 AlleleDelta` 不是通用規則
2. 它目前只在 `HCC1395 5kHz` 顯示明顯正增益
3. `Weak->Strong` 大多數仍是 TP，但 `HCC1395 5kHz` 的 FP rate 顯著偏高
4. `Noise->Strong` 幾乎只在 `HCC1395` 系列出現，尤其 `5kHz` 最值得警戒
5. `label-first` 不是沒用，但 `Strong` 類別很可能需要再細分
6. `samtools snapshot` 已可作為後續 read-level 驗證工具

## 9. 下一步建議

1. 不要把 `low VAF + high AlleleDelta` 直接升級成全域原始流程規則
2. 先把它做成：
   - `5kHz 特化候選規則`
   - 或 `ArtifactSuspect` 標記欄位
3. 下一輪優先做：
   - 對 `HCC1395 5kHz` 規則觸發集補更多 `region_samtools_snapshot`
   - 對 `Noise->Strong` / `Weak->Strong` 的 FP top regions 做 diagnostics 圖與 read-level 比對
   - 設計更細的 label 類別，測試 `Strong` 再切分是否能改善 F1
4. tumor-only 方向可接：
   - `HCC1395 / HCC1395_DORADO` 的 `ClairS-TO / DeepSomatic-TO`
   - 但應獨立做成新的 benchmark + label bundle，不要直接沿用 paired pure 結論
