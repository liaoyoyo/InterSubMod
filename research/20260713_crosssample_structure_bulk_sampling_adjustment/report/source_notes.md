<!--
建立時間: 2026-07-13
目標: 定義全樣本粗結構、bulk-sampling 校正與區域 PyClone-VI 候選標記 HTML 的證據鏈與主張上限
處理範圍: 技術讀者；GRCh38 chr1-22；7 dataset rows；HCC1395 same-cell-line cross-source pair；historical layered-v2
關聯檔案:
  - InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/rerun_a/
  - InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1/
  - InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/topology_pair_analysis.json
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/scripts/topology_pair_analysis.py
  - InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/
狀態: validated-report-input
-->

# HTML report source notes

## Reporting job

- **Task type**：B — comprehensive validation；七個 dataset rows 與 chr1–22 全部納入結構 census，沒有以 subset 取代全量結果。
- **User question**：先看全樣本結構分佈是否暗示相似 clone 演化，再用外部 subclone 軟體檢查 HCC1395 兩技術來源與其他樣本是否一致或類似。
- **Audience**：technical；優先保留分母、資料 gate、模型條件與可反駁的 claim ceiling。
- **Delivery**：canonical `artifact.json` 經 Data Analytics portable renderer 生成單一 self-contained HTML。
- **Historical scope**：InterSubMod 結構層使用 historical layered-v2 VAF-selected rooted-unlabeled mutation-state graph pattern。它是工程 snapshot，不是 clean layered-v3、single-cell truth 或經實驗驗證的 clone tree。
- **External model**：PyClone-VI 0.2.0；外部分群依 allele count/VAF、purity 與 allele-specific CN 推估 cellular prevalence。區域標記是「模型條件下的 possible state」，不是觀察到的 clone identity、clone count、ancestry 或 tree truth。

## Answer-first spine

1. 七個 dataset rows 的 complete-region 五分類是否有共同粗訊號。
2. HCC1395 與 HCC1395_DORADO 在 exact-coordinate matched pre-VAF coarse structure、VAF-selected final shape 與 aggregate composition 三個不同 estimands 下各自呈現什麼。
3. HCC pair 在 raw、Dirichlet、rarefaction、5/10 Mb block bootstrap、EB 與 source-standardization 後的 aggregate 相對距離排名是否穩定。
4. `C_read_groups`、`T`、`Topo`、VAF-selected most-likely result 與 external PyClone cluster 分別回答什麼問題。
5. 七樣本區域 PyClone state 的 full-exact-join coverage、possible subclone signal 與 fail-closed 狀況。
6. HCC pair 的 96.80% same-state 是否由 clonal majority 主導；minor/subclone-focused endpoints 是否仍一致。
7. 哪些數字支持「same-locus coarse 結構訊號跨來源部分再現」，哪些不支持「高度 clone equivalence」或「唯一真實樹 accuracy」。

## Chart map

| Chart | Question | Native family | Denominator | Claim ceiling |
|---|---|---|---|---|
| 7-sample complete five-class | 每個 dataset 的粗結構組成為何 | `stackedBar100` | each sample complete regions | Direct-dominant common output mode；沒有 solver/class-space null，不是 validity evidence 或 clone fraction |
| HCC JSD adjustment rank | 校正 sampling/source 後，HCC pair 是否成為最相似 | `bar` | HCC pair inserted beside 20 cross-biological pairs | robust relative composition rank；不是 biological equivalence |
| PyClone regional coverage | 外部模型在七樣本中哪些區域可判讀、哪些有 possible subclone signal | `stackedBar100` | all primary regions per sample | conditional regional state coverage；不是 clone abundance |

## Deliberate anti-inflation decisions

- 主要 five-class 圖使用 **all complete regions**，不把 unresolved 排除後再正規化；HCC/DORADO 的 unresolved 2.05%/9.90% 必須可見。
- 七樣本 7/7 Direct-dominant 只稱為共同輸出模式；尚無 solver/class-space/eligibility null，不能推論為方法有效性或共同生物演化。
- 不把 resolved-only Direct 73.20%/72.26% 單獨稱為完整 profile 高一致。
- 明確拆開四種 estimand：aggregate complete-five-class、matched pre-VAF coarse topology、matched VAF-selected final shape、external PyClone state/partition。不同 population 與 null 不互借。
- Matched pre-VAF coarse agreement 是 3,969/5,720=69.39%（95% CI 68.18%–70.57%），高於 chromosome-preserving null 39.51%，p=1/5001；這是最強正向的 same-locus coarse reproducibility endpoint，不是 exact-tree accuracy。
- Matched VAF-selected final-shape agreement 是 4,243/5,720=74.18%，composition JSD=0.1967；尚無 matched final-shape permutation null，不能借用 pre-VAF p-value。
- Aggregate HCC rank 9/21 代表 8/20 cross-biological dataset-row comparisons 更近，涵蓋 7 個 unique biological-ID pair combinations；不能只寫「rank 9」而隱去實際比較對象。
- `Topo=1` 共有 21,976 regions，但只有 10,832 是 `T=1`，另有 11,144 是 `T>1` 且共享同一 unlabeled shape；拓撲唯一不等於 labeled tree 唯一。
- HCC external same-state 使用 4,296/4,438=96.80%，但同時拆出 4,267 個 both-single-clonal-like regions；不把 majority-dominated accuracy 當 subclone reproducibility。
- Minor signal 使用 non-vacuous union：40/172=23.26%。Empty-union Jaccard 不設為 1。
- Final semantic-hotfix snapshot 為 32/32 PASS；5,520 個 subclonal union=0 rows 的 Jaccard 明確輸出 NA，partition exact 只在 34 個 `external_partition_informative=True` rows 定義。
- Partition endpoint 只以 both-multicluster 34 regions 作 informative denominator：21/34=61.76%；不使用由 single-cluster regions 主導的 aggregate partition percentage。
- 21/34 是較寬鬆的 **all-joined common-mutation partition endpoint**，可含 incomplete full-exact join；它不是 4,438 個 full-exact state-evaluable regions 的子分母。
- 不使用 `global_vs_regional_cluster_concordance.tsv`：既有 red-team 已發現 duplicate-key multiplicity，且 global/regional cluster label 的數值相等本身沒有跨 fit biological identity。
- `cluster_assignment_prob>=0.80` 只列 sensitivity；若任何 mutation 未達門檻，region fail closed 為 unavailable，不以剩餘高信心 clonal majority 宣稱完美一致。
- COLO829 與 HCC1937 缺少通過 gate 的 allele-specific CN，外部 subclone 層全數 blocked；禁止以 CN=2 代替。
- GENCODE/CGC/DGIdb/CLP enrichment 只用來排序可讀 region examples；DGIdb gene–drug claim 不代表 region variant 可用藥，也不是 topology accuracy evidence。
- FLI1、LRP1B、CBLC 三個醒目 examples 在兩來源皆 `C_read_groups=0`，且至少一側最小 assignment probability <0.80；因此沒有 InterSubMod read-level grouping corroboration。

## Representative rows

HTML 內嵌入 34 個「兩來源皆 multicluster」的 informative HCC partition regions（21 exact、13 different）與 15 個 CGC/approved-antineoplastic-context examples，合計 49 rows，低於 120-row 上限。完整 47,377-row annotation 保留於：

`research/20260713_crosssample_structure_bulk_sampling_adjustment/results/region_possible_clone_annotations_v1/region_possible_clone_annotations.tsv.gz`

## Decision statement

目前最合理結論是：七樣本的 Direct-dominant 只證明這是共同輸出模式，尚不能排除 solver/class-space/eligibility base rate。真正帶有 null 的正向證據是 HCC same-cell-line cross-source pair 在 5,720 個同座標 complete-both regions 的 pre-VAF coarse agreement 69.39%，高於 chromosome-preserving null 39.51%；它支持 coarse same-locus 結構訊號的跨來源部分再現。相對地，aggregate complete-five-class rank 只有 9/21（8/20 cross-biological dataset-row comparisons 更近；source-standardized 為 10/21），external minor subclone region Jaccard 只有 23.26%。因此 internal engineering 判為 **PASS_WITH_CAVEATS**，scientific/external release 在 fresh clean layered-v3 7/7 sample-root `_SUCCESS`、solver/class-space null 與 orthogonal truth 前維持 **NO-GO**；不支持兩來源已被證明具有高度一致的 clone/subclone 演化樹。

## Portable renderer compatibility receipt

- Canonical deliver 的第一次與第二次 browser QA 都在 desktop 1440px 以 `horizontal_overflow` 正確停止；不是 reader timeout、chart error、console error 或 external request。
- Headless-shell probe：`innerWidth=1440`、非 overlay scrollbar 使 `clientWidth=1425`；canonical reader `.analytics-top-bar` 使用 `width:100vw` 並置中，bounding box 為 `left=-7.5/right=1432.5`，造成 `scrollWidth=1433`。
- 包裝器只注入 `.analytics-top-bar{width:100%;margin-left:0;margin-right:0}`；不隱藏內容、不更動 artifact、數據、圖表 spec 或來源。
- Mobile QA 另由 artifact 本身縮短兩張 composition legend labels，並移除 rank chart 的重複 method-color legend；方法名稱與完整數值仍保留於 x-axis/table，沒有隱藏資料。
- 正式發布仍必須通過同一套 canonical `verifyPortableArtifact` 的 desktop 1440 與 mobile 390、source dialog、payload equality、geometry、external request、console error 與 page-level horizontal overflow 檢查；詳細數字寫入 `portable_delivery_receipt.json`。
