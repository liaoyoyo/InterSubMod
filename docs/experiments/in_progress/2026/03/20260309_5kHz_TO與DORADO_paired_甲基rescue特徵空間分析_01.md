<!--
建立時間: 2026-03-09 22:50
目標: 詳細分析 HCC1395 5kHz TO 與 HCC1395_DORADO paired 的 candidate-specific rescue feature space，確認是否存在其他可救回更多 TP 的甲基或標籤特徵，並區分單特徵與多特徵組合的效果
處理範圍:
  - HCC1395 5kHz TO candidate-specific rescue_joined_features
  - HCC1395_DORADO paired candidate-specific rescue_joined_features
  - 單特徵分佈、單特徵 rescue、雙特徵組合 rescue、相對同一 caller gate 的增量
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_methylation_rescue_feature_space.py
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260309_HCC1395_DORADO_paired_candidate_specific甲基rescue驗證_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260308_HCC1395_5kHz_TO_borderline_rescue特徵證據鏈整理_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/01/20260119_甲基化顯著性深入診斷與多層驗證策略報告_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/01/20260107_output_data_format_analysis_01.md
-->

# 5kHz TO 與 DORADO paired 的甲基 rescue 特徵空間分析

## 1. 研究問題

這輪要回答 4 個更細的問題：

1. 在已完成的兩個 candidate-specific `rescue_joined_features.tsv` 上，是否還有其他甲基或 phase/label 特徵能救回更多 TP。
2. 這些特徵若單獨使用，是否真的比既有 `PairwiseMedianDist>=0.20` / `agreement_positive` 更好。
3. 這些特徵若與其他特徵共同使用，是否能在同一 caller gate 下帶來額外增量，而不是只把候選池變窄。
4. 這些現象在 `HCC1395 5kHz TO` 與 `HCC1395_DORADO paired` 之間是否方向一致，是否足以支持跨樣本規則。

## 2. 輸入資料與輸出

### 2.1 兩個主要輸入檔

1. `HCC1395 5kHz TO`
   - [rescue_joined_features.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv)
2. `HCC1395_DORADO paired`
   - [rescue_joined_features.tsv](/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_joined_features.tsv)

### 2.2 本輪新增腳本

1. [analyze_methylation_rescue_feature_space.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_methylation_rescue_feature_space.py)

### 2.3 本輪主要輸出

1. [dataset_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/dataset_summary.tsv)
2. [numeric_feature_distribution.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/numeric_feature_distribution.tsv)
3. [categorical_feature_distribution.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/categorical_feature_distribution.tsv)
4. [gate_baselines.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/gate_baselines.tsv)
5. [single_feature_rule_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/single_feature_rule_sweep.tsv)
6. [combination_rule_sweep.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/combination_rule_sweep.tsv)
7. [top_rules_by_gate.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/top_rules_by_gate.tsv)
8. [feature_space_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260309_methylation_rescue_feature_space/feature_space_summary.md)

## 3. 分析方法

### 3.1 固定原則

1. caller gate 與甲基 support 分開計算。
2. 所有甲基 / label / phase support 規則都另外回報：
   - `delta_f1_vs_baseline`
   - `delta_f1_vs_gate`
3. 主判讀以 `delta_f1_vs_gate` 為準。

原因是：

1. `delta_f1_vs_baseline` 只代表「比最原始 kept set 好多少」。
2. 但若要判定某個甲基特徵是否真的有額外資訊，必須和**同一 caller gate**比較。
3. 否則很容易把「把 gate 變窄」誤解成「甲基真的有幫助」。

### 3.2 本輪固定 caller gate

1. `candidate_only`
2. `gq>=10`
3. `gq>=15`
4. `gq>=20`
5. `qual>=10 or gq>=20`

### 3.3 本輪觀察的甲基 / phase / label 特徵

1. `PairwiseMedianDist`
2. `PairwiseMeanDist`
3. `AlleleDelta`
4. `CramersV`
5. `GlobalP`
6. `Quality_Score`
7. `hp_assign_rate`
8. `allele_assign_rate`
9. `VerificationClass`
10. `agreement_positive`

特徵意義：

1. `GlobalP`
   - 全域 Fisher 檢定 p-value，見 [20260107_output_data_format_analysis_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/01/20260107_output_data_format_analysis_01.md)
2. `Quality_Score`
   - 多層驗證後的整體信心分數，`>=70` 視為高信心，見 [20260119_甲基化顯著性深入診斷與多層驗證策略報告_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/01/20260119_甲基化顯著性深入診斷與多層驗證策略報告_01.md)
3. `agreement_positive`
   - 由 [validate_method_design.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_method_design.py) 推得，定義為 `label_upgrade / consistent_strong / consistent_subclone`
4. `hp_assign_rate`
   - region reads 中帶有有效 HP tag 的比例
5. `PairwiseMedianDist / PairwiseMeanDist`
   - read-level 甲基距離摘要

## 4. 先看資料可分析範圍

| dataset | mode | candidate rows | analyzed rows | lost TP coverage | removed FP coverage | baseline F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `HCC1395 5kHz TO` | `to-pure` | 1071 | 888 | 87.32% | 71.48% | 0.712697 |
| `HCC1395_DORADO paired` | `paired-pure` | 2780 | 419 | 8.29% | 19.66% | 0.859176 |

關鍵差異：

1. `5kHz TO` 的 candidate-specific coverage 很高，足以對甲基特徵做實質比較。
2. `DORADO paired` 的 coverage 明顯不足，所以任何跨樣本判讀都必須保守。

## 5. 各特徵分佈觀察

### 5.1 5kHz TO 的中位數分佈

| feature | TP median | FP median | 觀察 |
| --- | ---: | ---: | --- |
| `gq` | 15.0 | 11.0 | caller 端仍有明顯差異 |
| `qual` | 15.3536 | 11.0829 | caller 端仍有明顯差異 |
| `af` | 0.2179 | 0.1964 | 差異小 |
| `PairwiseMedianDist` | 0.2272 | 0.2303 | 幾乎無分離，FP 甚至略高 |
| `PairwiseMeanDist` | 0.2411 | 0.2449 | 幾乎無分離 |
| `AlleleDelta` | 0.0146 | 0.0095 | TP 略高，但差異不大 |
| `GlobalP` | 0.4590 | 0.5680 | TP 較低，但不是強分離 |
| `Quality_Score` | 75.0 | 75.0 | 中位數完全相同 |
| `hp_assign_rate` | 0.977778 | 0.975610 | 幾乎無差 |
| `allele_assign_rate` | 1.0 | 1.0 | 完全不能當判別特徵 |

推論：

1. `5kHz TO` 的甲基特徵不是「整體分佈明顯分開」的型態。
2. 這也是為什麼 `Pairwise>=0.20` 雖然可以 rescue 一批 TP，但沒有辦法在同一 caller gate 上額外加分。
3. `allele_assign_rate` 與 `CramersV<=0.05` 看似數值很好，但其實更像技術 proxy，不應升級為正式 rescue 特徵。

### 5.2 DORADO paired 的中位數分佈

| feature | TP median | FP median | 觀察 |
| --- | ---: | ---: | --- |
| `gq` | 18.0 | 11.0 | caller 端仍有明顯差異 |
| `qual` | 18.2620 | 11.9390 | caller 端仍有明顯差異 |
| `af` | 0.6543 | 0.9770 | TP 與 FP 的 allele 組成差異大 |
| `PairwiseMedianDist` | 0.1232 | 0.1323 | **方向與 5kHz 不同，TP 較低** |
| `PairwiseMeanDist` | 0.1347 | 0.14465 | **TP 較低** |
| `AlleleDelta` | 0.0074 | 0.0 | TP 較高 |
| `GlobalP` | 0.4525 | 1.0 | TP 較低 |
| `Quality_Score` | 75.0 | 60.0 | TP 較高 |
| `hp_assign_rate` | 1.0 | 0.90625 | TP 較高 |
| `allele_assign_rate` | 1.0 | 1.0 | 同樣不是有效判別特徵 |

推論：

1. `DORADO paired` 最有價值的方向不是「高 Pairwise」，而是：
   - `高 hp_assign_rate`
   - `較低 Pairwise`
   - `較低 GlobalP`
   - `較高 Quality_Score`
2. 這代表 `DORADO paired` 的正向 support 更像「結構穩定且 phase 品質好」，而不是「距離越高越像 TP」。

## 6. 單特徵 rescue：哪些真的有用

### 6.1 5kHz TO

固定拿最重要的 caller gate `gq>=10` 比較：

| 規則 | TP rescued | FP reintroduced | F1 delta vs baseline | F1 delta vs gate | 解讀 |
| --- | ---: | ---: | ---: | ---: | --- |
| `gq>=10` | 499 | 119 | +0.006943 | 0 | caller gate baseline |
| `gq>=10 + Quality_Score>=60` | 395 | 89 | +0.005551 | -0.001392 | 可用，但只是較窄的 support |
| `gq>=10 + PairwiseMedianDist>=0.15` | 390 | 92 | +0.005445 | -0.001498 | 有效，但沒有額外增量 |
| `gq>=10 + PairwiseMedianDist>=0.20` | 300 | 68 | +0.004219 | -0.002724 | 先前主規則成立，但仍低於 gate |
| `gq>=10 + agreement_positive` | 148 | 25 | +0.002163 | -0.004780 | 很乾淨，但非常窄 |
| `gq>=10 + Strong/Subclone` | 149 | 30 | +0.002134 | -0.004809 | 比 `agreement_positive` 稍髒 |

重點：

1. `5kHz TO` 沒有找到任何甲基或 label 特徵，能在 `gq>=10` 之上繼續提高 F1。
2. 但有幾個「可解釋且可用」的 support 特徵：
   - `Quality_Score>=60`
   - `PairwiseMedianDist>=0.15`
   - `PairwiseMedianDist>=0.20`
   - `agreement_positive`
3. 這些規則都能從 baseline 救回 TP，但都**沒有超過同一 caller gate**。

進一步看保留比例：

1. `Pairwise>=0.20`
   - 保留 `300/499 = 60.1%` 的 gate TP
   - 同時保留 `68/119 = 57.1%` 的 gate FP
   - 代表它不是強分離，而是中等強度的 support
2. `agreement_positive`
   - 保留 `148/499 = 29.7%` 的 gate TP
   - 同時保留 `25/119 = 21.0%` 的 gate FP
   - 代表它較窄，但 precision 明顯較好

### 6.2 DORADO paired

固定拿目前 paired 最合理的 caller gate `gq>=15` 比較：

| 規則 | TP rescued | FP reintroduced | F1 delta vs baseline | F1 delta vs gate | 解讀 |
| --- | ---: | ---: | ---: | ---: | --- |
| `gq>=15` | 97 | 88 | +0.000502 | 0 | caller gate baseline |
| `gq>=15 + hp_assign_rate>=0.99` | 50 | 15 | +0.000634 | +0.000132 | 本輪最佳單特徵增量 |
| `gq>=15 + GlobalP<=0.50` | 35 | 20 | +0.000327 | -0.000176 | 方向正確但不夠好 |
| `gq>=15 + PairwiseMedianDist<=0.20` | 60 | 59 | +0.000255 | -0.000248 | 單看低 Pairwise 不夠 |
| `gq>=15 + Quality_Score>=60` | 49 | 47 | +0.000223 | -0.000280 | 單看 Quality 不夠 |

重點：

1. `DORADO paired` 唯一明確單特徵正增量，是 `hp_assign_rate>=0.99`。
2. `Pairwise` 在這裡方向反轉，應該看 **較低的 Pairwise**，不是 `>=0.20`。
3. `Quality_Score` 和 `GlobalP` 單獨使用都還不夠，需要和 `hp_assign_rate` 共同分析。

## 7. 雙特徵組合：哪些有加成

### 7.1 5kHz TO

這輪沒有找到任何雙特徵組合能超過 `gq>=10` gate。

較好的組合仍然只是較窄的 support，不是額外加分：

| 規則 | TP rescued | FP reintroduced | F1 delta vs baseline | F1 delta vs gate |
| --- | ---: | ---: | ---: | ---: |
| `gq>=10 + Pairwise>=0.20 + Quality>=60` | 249 | 56 | +0.003508 | -0.003435 |
| `gq>=10 + Pairwise>=0.20 + Quality>=75` | 227 | 45 | +0.003254 | -0.003689 |
| `gq>=10 + GlobalP<=0.50 + Quality>=60` | 217 | 43 | +0.003111 | -0.003832 |

結論：

1. `5kHz TO` 的甲基 support 目前以**單一第二層 support**最合理。
2. 把多個條件硬交集只會再縮小 TP。

### 7.2 DORADO paired

`DORADO paired` 反而出現了幾個值得記住的雙特徵組合：

| 規則 | TP rescued | FP reintroduced | F1 delta vs baseline | F1 delta vs gate | 解讀 |
| --- | ---: | ---: | ---: | ---: | --- |
| `gq>=15 + Pairwise<=0.20 + hp_assign>=0.99` | 46 | 12 | +0.000606 | +0.000103 | 最佳雙特徵增量 |
| `gq>=15 + Quality>=60 + hp_assign>=0.99` | 36 | 6 | +0.000516 | +0.000013 | 很乾淨，但增量極小 |
| `gq>=15 + Quality>=75 + hp_assign>=0.99` | 33 | 5 | +0.000479 | -0.000023 | 過窄，開始傷 TP |

保留比例：

1. `gq>=15 + hp_assign>=0.99`
   - 保留 `50/97 = 51.5%` 的 gate TP
   - 只保留 `15/88 = 17.0%` 的 gate FP
2. `gq>=15 + Pairwise<=0.20 + hp_assign>=0.99`
   - 保留 `46/97 = 47.4%` 的 gate TP
   - 只保留 `12/88 = 13.6%` 的 gate FP

這代表：

1. `DORADO paired` 的 support 要成立，關鍵不是高甲基距離，而是**高 HP 可用性**。
2. 低 Pairwise 在這裡比較像穩定結構，而不是弱訊號。

## 8. 哪些特徵不能當成正式結論

### 8.1 `allele_assign_rate>=0.99`

這個規則在數值排序上很高，但**不應解讀為甲基 rescue 特徵**。

原因：

1. `5kHz TO` 與 `DORADO paired` 的 TP/FP 中位數都幾乎是 `1.0`
2. 它更像「該 region 讀段成功被 allele 標記」的技術 proxy
3. 在 `5kHz TO` 中，`candidate_only + allele_assign>=0.99` 幾乎等於整個 analyzed subset，不具真正判別力

### 8.2 `CramersV<=0.05`

這個規則也不能升級。

原因：

1. 兩個 dataset 的 TP/FP 中位數都幾乎是 `0`
2. 它更像目前 candidate-specific pool 的常態，而不是可分辨 TP/FP 的特徵

## 9. 證據鍊與整體判讀

### 9.1 5kHz TO

1. `Pairwise>=0.20`、`Quality>=60`、`agreement_positive` 都能從 baseline 救回 TP。
2. 但它們全部都沒有超過 `gq>=10` 這個 caller gate。
3. 所以 `5kHz TO` 的甲基特徵目前定位應該是：
   - `support`
   - `ranking`
   - `annotation`
4. 不應改寫成「比 caller 更強的主 rescue 規則」。

### 9.2 DORADO paired

1. 高 Pairwise 並不成立。
2. 真正有價值的方向是：
   - `hp_assign_rate>=0.99`
   - `Pairwise<=0.20`
   - `Quality_Score>=60`
3. 這些特徵只有在 `gq>=15` 這種 caller gate 之後，才會出現小幅正增量。
4. 即使如此，它們仍然沒有超過 `gq>=20` 的 caller-only 效果。

### 9.3 跨樣本結論

目前不能說：

1. `PairwiseMedianDist>=0.20` 是跨樣本 rescue 規則
2. `agreement_positive` 是跨樣本 rescue 規則
3. 甲基資料本身已可穩定替代 caller

目前可以說的是：

1. `5kHz TO` 的甲基 rescue 是真的，但以 `support` 為主
2. `DORADO paired` 也有 support 訊號，但方向不同，而且更依賴 `hp_assign_rate`
3. 尚未找到可直接共用的單一全域規則

## 10. 結論

### 10.1 正式答案

1. **有其他可用的甲基/phase 特徵，但它們不是全域規則。**
2. `5kHz TO`：
   - 最值得記住的是 `Quality_Score>=60`、`PairwiseMedianDist>=0.15/0.20`、`agreement_positive`
   - 但它們都沒有超過 `gq>=10`
3. `DORADO paired`：
   - 最值得記住的是 `hp_assign_rate>=0.99`
   - 次佳是 `Pairwise<=0.20 + hp_assign>=0.99`
   - 表示這裡要看的是「高 phase 可用性 + 穩定低距離」而不是高距離

### 10.2 最合理的流程定位

目前最合理的 rescue 分工是：

1. 第一層：`caller-first`
2. 第二層：
   - `5kHz TO` 用 `Quality / Pairwise / agreement` 做 support 或排序
   - `DORADO paired` 用 `hp_assign + 低 Pairwise/高 Quality` 做 support
3. 不要做的事：
   - 不要把 `Pairwise>=0.20` 升級為跨樣本規則
   - 不要把 `allele_assign_rate` 或 `CramersV<=0.05` 當成正式甲基 rescue 特徵

## 11. 下一步建議

1. 若要繼續 `5kHz TO`：
   - 對 `Quality>=60`、`Pairwise>=0.15`、`agreement_positive` rescue 到的 top TP / FP 做 diagnostics
2. 若要繼續 `DORADO`：
   - 優先做 `HCC1395_DORADO TO candidate-specific InterSubMod`
   - 驗證 `hp_assign>=0.99` 與 `低 Pairwise + 高 hp_assign` 是否能在 TO 重現
3. 若要進入流程化：
   - 先做 annotation / ranking，不做 hard keep
