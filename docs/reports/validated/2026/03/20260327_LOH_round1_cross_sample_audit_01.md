<!--
建立時間: 2026-03-27
目標: 完成 InterSubMod LOH Round 1 cross-sample audit 正式報告
處理範圍: 7 paired + 7 TO 具甲基資料樣本，LOH-like / HP imbalance / same-locus paired-vs-TO compare / case panel
關聯檔案:
  - docs/plans/2026/03/20260326_LOH盤點執行規格_01.md
  - docs/experiments/in_progress/2026/03/20260327_LOH_round1_cross_sample_audit啟動與決策紀錄_01.md
  - output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/
  - scripts/analysis/build_loh_round1_cross_sample_audit.py
-->

# LOH Round 1 Cross-Sample Audit 正式報告

> 生成時間：2026-03-27
> 性質：observation-first / diagnostics-first，非 intervention round

> ⚠️ **勘誤聲明（2026-03-30）**：本報告中所有 **TO 端** HP 相關統計（LOH enrichment、effective_hp 分佈、eff_hp<10 比例）均基於 HP integer tag bug（ReadParser.cpp HP:i:11/21/33 mapping 錯誤）計算所得，數值無效。**Paired 端結論完全有效，不受影響。**
> 修正後 TO LOH 完整分析見：[20260330_TO_LOH_enrichment_post_hp_fix_01.md](20260330_TO_LOH_enrichment_post_hp_fix_01.md)
> 具體修正數字見各節 `[修正]` 標注。
> 主要輸出 workspace：[20260327_loh_round1_cross_sample_audit](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit)
> 分析腳本：[build_loh_round1_cross_sample_audit.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit.py)

---

## 一、這輪要回答的問題

這輪不是要直接證明 `LOH` 能提升 `F1`，而是先盤清楚三件事：

1. `LOH-like / HP family imbalance` 在 7 個具甲基資料樣本的 paired / TO 中，實際分佈長什麼樣子。
2. 這些 `LOH-like` 現象在 `TP / FP` 之間，是否有穩定且可泛化的差異。
3. paired 與 TO 在同樣本同位點上，究竟是互相支持，還是呈現大量不一致與 mode-specific failure pattern。

這個 round 的核心定位，是替後續 `InterSubMod` 五個研究願景中的：

1. `region-first evidence panel`
2. `tumor-only 獨立主線`
3. `paired 與 TO 的互相驗證`

先建立第一層可追溯的診斷底圖。

---

## 二、固定定義與方法

### 2.1 資料範圍

Round 1 共納入 `14` 個 dataset：

1. `7` 個 canonical paired_full complete_matrix
2. `7` 個 tumor-only pilot / fastresume round

樣本如下：

1. `HCC1395_HKU_5kHz`
2. `HCC1395_DORADO`
3. `COLO829`
4. `H1437`
5. `H2009`
6. `HCC1937`
7. `HCC1954`

主輸入 manifest：

- [input_manifest.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/input_manifest.tsv)

所有 `14/14` dataset 均成功納入，沒有 sample 因缺檔被排除。

### 2.2 本輪主定義

本輪沒有直接使用原始 `HP_Ratio` 當主判斷，而是重算：

1. `effective_hp_reads = HP1FamilyN + HP2FamilyN`
2. `hp_ratio_core = HP1FamilyN / effective_hp_reads`
3. `core_loh_like = (effective_hp_reads > 0) AND (hp_ratio_core < 0.1 OR hp_ratio_core > 0.9)`

保留原始欄位作對照：

1. `Potential_LOH`
2. `LOH_Subtype`
3. `VerificationClass`

### 2.3 PS 的處理方式

1. TO：從 `tumor_phased.vcf` 補查 `PS / GT / GT2 / GT3`
2. paired：目前沒有 variant-level `PS` 輸出，只能在 case panel 用 `tagged.bam` 的 read-level `PS` 輔助觀察
3. 本輪明確不把不同 `PS block` 的 `HP1 / HP2` 當成可跨 block 串聯的親源資訊

### 2.4 本輪實際輸出

Round 1 workspace 主要產物：

1. [all_region_rows.tsv.gz](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz)
2. [loh_enrichment_by_sample_mode.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/loh_enrichment_by_sample_mode.tsv)
3. [hp_coverage_qc_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/hp_coverage_qc_summary.tsv)
4. [verificationclass_by_loh_subtype.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/verificationclass_by_loh_subtype.tsv)
5. [same_locus_compare.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/same_locus_compare.tsv)
6. [cases/case_registry.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/cases/case_registry.tsv)
7. [round_summary.md](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/round_summary.md)

固定圖表：

**Fig01 — LOH-like Fraction Overview（各樣本 paired/TO LOH-like 佔比概覽）**

![Fig01 LOH-like Fraction Overview](../../../../../research/loh_investigation/figures/loh_round1/fig01_loh_like_fraction_overview.png)

**Fig02 — HP Ratio Core Distribution（hp_ratio_core 分佈）**

![Fig02 HP Ratio Core Distribution](../../../../../research/loh_investigation/figures/loh_round1/fig02_hp_ratio_core_distribution.png)

**Fig03 — Effective HP vs HP Ratio Scatter（effective_hp_reads 與 hp_ratio_core 散點圖）**

![Fig03 Effective HP vs HP Ratio Scatter](../../../../../research/loh_investigation/figures/loh_round1/fig03_effective_hp_vs_hp_ratio_scatter.png)

**Fig04 — VerificationClass × LOH Subtype Structure（VerificationClass 與 LOH_Subtype 交叉結構）**

![Fig04 VerificationClass LOH Subtype Structure](../../../../../research/loh_investigation/figures/loh_round1/fig04_verificationclass_lohsubtype_structure.png)

**Fig05 — LOH vs non-LOH Feature Boxplots（LOH-like vs 非 LOH-like 的特徵分佈對比）**

![Fig05 LOH non-LOH Feature Boxplots](../../../../../research/loh_investigation/figures/loh_round1/fig05_loh_nonloh_feature_boxplots.png)

**Fig06 — Sample Bin Heatmap（各樣本 LOH 分類熱圖）**

![Fig06 Sample Bin Heatmap](../../../../../research/loh_investigation/figures/loh_round1/fig06_sample_bin_heatmap.png)

---

## 三、Round 1 的直接事實

### 3.1 盤點規模

本輪最終整合得到：

1. `748,391` 個 region rows
2. `459,782` 個 same-locus paired-vs-TO union rows
3. `5` 個代表案例

這表示 Round 1 已經達成最初要求的：

1. 7 樣本 paired 軸
2. 7 樣本 TO 軸
3. same-locus compare
4. case panel

### 3.2 `HP_Ratio` 本身不能單獨解讀

這輪確認了一個很重要的 schema 事實：

1. `effective_hp_reads = 0` 的 region 總共有 `69,807` 個
2. 這 `69,807` 個 region 的原始 `HP_Ratio` 全部都是 `0.5`
3. 這些 region 的 `Potential_LOH` 全部是 `False`

也就是說：

1. `HP_Ratio = 0.5` 並不一定代表「平衡」
2. 它也可能只是「根本沒有有效 HP1/HP2 read」

因此後續研究若直接拿原始 `HP_Ratio` 當主要特徵，會把：

1. 真正 balanced
2. no effective HP support

混成同一類。

這是本輪最關鍵的 schema 修正結論之一。

### 3.3 paired 與 TO 的 `LOH-like` 整體型態不同

依 [loh_enrichment_by_sample_mode.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/loh_enrichment_by_sample_mode.tsv) 加總：

| mode | TP total | FP total | TP LOH-like fraction | FP LOH-like fraction | FP/TP enrichment |
| --- | --- | --- | --- | --- | --- |
| paired | 325,270 | 3,429 | 29.33% | 35.02% | 1.194× |
| ~~TO~~ | ~~291,310~~ | ~~128,382~~ | ~~60.29%~~ | ~~54.99%~~ | ~~0.912×~~ |

> **[修正 2026-03-30]** TO 行數字無效（HP bug）。修正後正確值：
> | TO | 291,310 | 128,382 | **44.5%** | **35.8%** | **0.805×** |
> 參見：[20260330_TO_LOH_enrichment_post_hp_fix_01.md](20260330_TO_LOH_enrichment_post_hp_fix_01.md)

這代表：

1. paired 整體上有輕度 `FP` enrichment（結論有效）
2. ~~TO 整體上反而不是 `FP` enrichment，而是 `TP` 與 `FP` 都大量 `LOH-like`~~

> **[修正 2026-03-30]** 修正後 TO LOH-like 比例大幅降低（44.5% TP vs 舊版 60.3%），方向確認為 **TP 富集（0.805×）**，即 TO LOH-like 在 TP 比 FP 更常見。

這個結果很重要，因為它直接說明：

1. `LOH-like` 在 paired 可以是偏風險訊號（結論有效）
2. `LOH-like` 在 TO 是 **TP 富集訊號**（不是背景型 annotation，也不是 FP marker）

### 3.4 paired 的 `LOH-like` 具有 sample heterogeneity

paired mode 各樣本 `FP/TP enrichment` 差異很大：

| sample | paired FP/TP enrichment |
| --- | --- |
| HCC1954 | 3.185× |
| H2009 | 2.685× |
| H1437 | 1.795× |
| HCC1937 | 1.505× |
| HCC1395_DORADO | 1.260× |
| COLO829 | 1.155× |
| HCC1395_HKU_5kHz | 1.016× |

這表示 paired 的 `LOH-like` 也不能直接被寫成單一規則；它是：

1. 有些樣本明顯偏 FP
2. 有些樣本幾乎沒有區分力

因此 paired 端比較合理的用法，是把 `LOH-like` 當：

1. sample-aware diagnostics
2. candidate risk feature

而不是全域 hard filter。

### 3.5 TO 幾乎所有樣本都呈現「LOH-like 在 TP 更常見」

~~TO mode 各樣本的 `FP/TP enrichment` 幾乎都接近 `1`~~

> **[修正 2026-03-30]** 舊版 per-sample enrichment 全部無效（HP bug）。修正後正確值：

| sample | TO FP/TP enrichment（修正後）| 說明 |
| --- | --- | --- |
| HCC1954 | **0.852×** | TP 富集（-FP risk）|
| HCC1937 | **0.882×** | TP 富集 |
| HCC1395 5kHz | **0.896×** | TP 富集 |
| HCC1395_DORADO | **0.909×** | TP 富集 |
| H1437 | **0.919×** | TP 富集 |
| H2009 | **0.923×** | TP 富集 |
| COLO829 | **0.956×** | 微弱 TP 富集（p=4.5e-4）|

修正後結論：TO LOH-like **全面呈現 TP 富集**（enrichment < 1），且方向一致。

直接結論修正：

1. `LOH-like` 不適合在 TO 當 FP filter（結論不變，但原因修正）
2. TO LOH-like 是弱 TP 支持訊號（-FP risk feature），而非中性 annotation

### 3.6 TO 的 phase support 問題比 paired 明顯

依 [hp_coverage_qc_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/hp_coverage_qc_summary.tsv)：

~~TO 中 `effective_hp_reads < 10` 的比例很高~~

> **[修正 2026-03-30]** 此節數字全部無效（HP bug）。修正後 TO eff_hp 分佈大幅改善：

| sample | mode | median eff_hp（修正後）| LOH eligible（eff_hp≥30）|
| --- | --- | --- | --- |
| COLO829 | TO | **25**（舊: 5）| 34.7%（舊: 0.7%）|
| H1437 | TO | **66**（舊: 13）| 95.5%（舊: 24.5%）|
| H2009 | TO | **85** | 97.7%（舊: 64.5%）|
| HCC1395 | TO | **61** | 92.2%（舊: 58.5%）|
| HCC1937 | TO | **99** | 97.4%（舊: 86.7%）|
| HCC1954 | TO | **61** | 93.0%（舊: 70.0%）|

原本標注的「COLO829 TO TP eff_hp<10 比例 67.11%」現已作廢。
| COLO829 | TO | FP | 63.52% |
| HCC1395_DORADO | TO | FP | 45.08% |
| H1437 | TO | FP | 45.08% |
| H1437 | TO | TP | 43.82% |

paired 相對乾淨得多，但也有 sample-specific outlier：

| sample | mode | truth | `effective_hp_reads = 0` |
| --- | --- | --- | --- |
| H1437 | paired | FP | 37.50% |
| H1437 | paired | TP | 31.09% |
| HCC1395_DORADO | paired | TP | 11.43% |

這說明：

1. TO 端大量 region 的 HP 支持本身就偏弱
2. paired 端雖然整體較好，但仍有樣本相依的 phase/QC 問題

### 3.7 same-locus paired-vs-TO 的主要不一致，不是 `both_fp`，而是 `TO-only FP`

依 [same_locus_compare.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/same_locus_compare.tsv) 加總：

| concordance | count | fraction |
| --- | --- | --- |
| both_tp | 287,092 | 62.44% |
| to_only_fp | 126,865 | 27.59% |
| paired_only_tp | 38,178 | 8.30% |
| to_only_tp | 4,218 | 0.92% |
| paired_only_fp | 1,912 | 0.42% |
| both_fp | 1,517 | 0.33% |

這裡最重要的結論不是「paired 與 TO 都有一些 FP」，而是：

1. 真正支配 discordance 的，是大量 `TO-only FP`
2. `paired_only FP` 與 `both FP` 都很少
3. `paired_only TP` 也不少，代表 TO 不只 FP 多，還有 sensitivity gap

這個結果直接支援把 TO 視為獨立主線問題，而不是 paired 的簡化版。

### 3.8 `LOH_Subtype` 在 TO 變成更廣泛的背景註解

依 [verificationclass_by_loh_subtype.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/verificationclass_by_loh_subtype.tsv) 聚合：

paired TP：

1. `None = 70.67%`
2. `LOH_Noise = 20.04%`
3. `LOH_Weak = 5.18%`
4. `LOH_Strong = 3.13%`
5. `LOH_Subclone = 0.97%`

TO TP：

1. `None = 39.71%`
2. `LOH_Noise = 30.44%`
3. `LOH_Weak = 17.40%`
4. `LOH_Strong = 10.68%`
5. `LOH_Subclone = 1.77%`

這表示 TO 中 `LOH_*` 子類型佔比大幅提高；因此 `LOH_Subtype` 在 TO 更像：

1. 大背景下的描述標籤
2. failure-mode mapping 的一部分

而不是直接代表高辨識力事件。

### 3.9 TO 的 PS 補查比例高，但 paired 仍缺 variant-level PS

本輪 TO `PS` 非空比例為 `94.6918%`。

這代表：

1. TO 其實有足夠高的 `PS` 可做 summary annotation
2. paired 仍然卡在沒有 variant-level `PS` export

所以目前 phase-aware 解釋仍是：

1. TO 可以做到「這個位點在某個 PS block 內」
2. paired 只能在 case panel 檢查 `tagged BAM` 的 read-level `PS`
3. 兩者都還不能跨 block 講親源故事

---

## 四、Round 1 的推論與判斷

### 4.1 `LOH` 適合當第一步觀察，但不適合當第一個主證明

這輪支持把 `LOH / HP imbalance` 當第一步 observation 的理由是：

1. 不用改 C++ 主流程，就能直接從既有 summary / VCF / tag.bam 得到大量診斷訊號
2. 它很快揭露 paired 與 TO 的 failure-mode 差異
3. 它能立刻指出哪些樣本是 phase/QC outlier

但這輪同時也明確證明：

1. 在 TO，`LOH-like` 不是 TP/FP 的主分離器
2. 在 paired，它也有顯著 sample heterogeneity

所以 `LOH` 最合理的定位是：

1. `evidence panel` 的一層
2. `diagnostics / annotation / risk mapping`

而不是單獨的主決策器。

### 4.2 paired 與 TO 可以互相輔助，但輔助方向不對稱

這輪資料支持的輔助方式是：

1. paired 可以提供較乾淨的 reference，告訴我們哪些 `LOH-like` 在較高 phase support 下仍然存在
2. TO 可以提供大量 `FP` 壓力測試，特別是 `TO-only FP`

也就是：

1. paired 比較適合幫忙定義「可解釋 evidence panel」
2. TO 比較適合幫忙驗證「這組 evidence 是否真能處理最難的 FP 場景」

這與你先前定義的「paired 與 TO 可互相輔助，但 TO 是獨立主線」是一致的。

### 4.3 `region-first` 是對的，但一定要帶上 support quality

本輪沒有直接往 per-CpG 展開，先做 `region-first` 是正確的，因為：

1. 現有 summary / same-locus / case 都是 region-level
2. 先把 `region-first` 的 evidence panel 定義穩定，後續 per-CpG 才有意義

但本輪也指出，`region-first` 不能只看生物標籤，還必須同時帶：

1. `effective_hp_reads`
2. `hp0_ratio`
3. `hp3_ratio`
4. `PS available or not`

否則會把：

1. 真正的 one-sided structure
2. phase support 不足造成的假象

混在一起。

---

## 五、代表案例

本輪自動建立 `5` 個代表案例：

1. [case_01_both_tp_loh_like_HCC1954_chr17_39520424](../output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/cases/case_01_both_tp_loh_like_HCC1954_chr17_39520424)
2. [case_02_both_fp_loh_like_HCC1395_DORADO_chr3_149726729](../output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/cases/case_02_both_fp_loh_like_HCC1395_DORADO_chr3_149726729)
3. [case_03_paired_only_loh_like_H1437_chr5_1715227](../output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/cases/case_03_paired_only_loh_like_H1437_chr5_1715227)
4. [case_04_to_only_loh_like_HCC1954_chr5_1314055](../output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/cases/case_04_to_only_loh_like_HCC1954_chr5_1314055)
5. [case_05_high_hp0_or_hp3_COLO829_chr1_485314](../output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/cases/case_05_high_hp0_or_hp3_COLO829_chr1_485314)

其中最值得優先看的是：

1. `case_01`
   - paired 與 TO 都是 `TP`
   - 兩側都呈現高 `effective_hp_reads`
   - 是「跨 mode 一致的 LOH-like evidence」代表
2. `case_02`
   - paired 與 TO 都是 `FP`
   - 但 paired 端仍有強 one-sided HP 支持，TO 端則退化成 `effective_hp_reads = 0`
   - 這很像「同位點但 mode 間 evidence quality 不對稱」的典型
3. `case_04`
   - `TO-only FP`
   - 正好對應本輪最主要的不一致型態

---

## 六、風險與限制

### 6.1 仍未解的限制

1. paired 沒有 variant-level `PS` export
2. 不同 `PS block` 的 `HP1 / HP2` 不能外推出父母源
3. 本輪只做到 region / locus-level，尚未與 CNV/LOH interval、gene、second-hit 整合
4. `LOH-like` 在 TO 中不是 discriminative signal，而是高背景現象

### 6.2 這輪不應過度宣稱的事情

本輪**不能**直接宣稱：

1. `LOH` 已可提升 TO 的 `F1`
2. `LOH-like` 已可單獨區分 TP/FP
3. TO 的 `PS` 已足夠支撐長程 haplotype 敘事

本輪**可以**穩定宣稱：

1. `LOH / HP imbalance` 已足夠作為第一層 observation workspace
2. `HP_Ratio` 需要和 `effective_hp_reads` 一起解讀
3. `TO-only FP` 是 paired-vs-TO discordance 的主要來源

---

## 七、下一步建議

### 7.1 最適合立即進行的小步驗證

1. 針對 paired 中 `FP/TP enrichment` 高的樣本：
   - `HCC1954`
   - `H2009`
   - `H1437`
   - `HCC1937`
   做第二輪 `LOH-like FP vs LOH-like TP` 深化 case study
2. 對 TO 先做 `support-aware stratification`：
   - `effective_hp_reads`
   - `hp0_ratio`
   - `PS available`
   之後再談 `LOH` 的辨識力
3. 把 `same_locus_compare` 與 `tumor_phased_LOH.bed / CNV` 接起來，從 `locus` 升級到 `interval-aware` 解釋

### 7.2 若要往研究願景繼續推

最合理的順序是：

1. 先把 `region-first evidence panel` 定義穩
2. 再把 paired 的 `PS-aware` 輸出補齊
3. 最後才進入 `事件順序 inference / second-hit` 的正式推理

這樣路線才不會在 phase block 與 evidence quality 還沒釐清前，就過早進入高風險敘事。

---

## 八、最終判斷

Round 1 已經成功回答這輪最重要的問題：

1. `LOH` 值得當第一步觀察軸
2. 但它目前更像 `annotation / diagnostics / evidence panel`，不是主判斷器
3. TO 的主要問題不是「缺 LOH-like 位點」，而是「LOH-like 在 TP/FP 都很多，且 TO-only FP 大量存在」
4. paired 與 TO 可以互相輔助，但角色不同：
   - paired 負責建立較乾淨的可解釋參考
   - TO 負責驗證方法是否真能撐住最難的 FP 場景

因此，這輪結論支持你目前的主軸：

1. 五個目標應作為最終研究願景
2. `LOH` 先作第一步觀察與 evidence panel 建構
3. `tumor-only` 獨立成主線，但持續用 paired 當互相對照與驗證來源
