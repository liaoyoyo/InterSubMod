<!--
建立時間: 2026-07-13 00:00 +0800
目標: 記錄全 7 dataset rows 的 structure composition、bulk-sampling adjustment 與 HCC1395 technical-pair 相對位置
處理範圍: chr1-22；7 dataset rows；6 biological IDs；historical layered-v2 engineering snapshot
關聯檔案:
  - InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/pre-decision-audit.md
  - InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/scripts/analyze_crosssample_composition.py
  - InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/rerun_a/summary.json
  - InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/results/reproducibility/reproducibility_receipt.json
-->

# 跨樣本結構 composition 與 bulk-sampling adjustment：完整執行紀錄

> **Share with caveats。HCC1395 technical pair 在 complete 五類 primary 並非異常接近：JSD=0.1984，插入 20 個 cross-biological dataset-row pairs 後 rank=9；sample size、5/10 Mb spatial blocks 與 EB shrinkage 都沒有把它變成 rank 1。**

## 1. Round 定義

- Task type：**B Comprehensive validation**；全 7 dataset rows，無 subset。
- 服務目標：**G4 跨樣本一致性／reproducibility**。
- Dataset rows：HCC1395、HCC1395_DORADO、COLO829、H1437、H2009、HCC1937、HCC1954。
- Biological IDs：6；HCC1395與HCC1395_DORADO是同一 biological ID 的 technical pair，不當成兩個 biological replicates。
- Primary estimand：complete五類＝Single、Sister、Direct、Sister+direct、Unresolved。
- Secondary estimand：resolved-only四類；只供和20260712既有報告對照。
- Benchmark TP/FP/FN/precision/recall/F1：不適用；本輪不是caller-vs-truth benchmark，而是composition uncertainty與technical reproducibility分析。

### 研究問題

1. 全7 dataset的類別比例與不確定區間為何？
2. HCC technical pair是否比不同 biological IDs 的 dataset pairs 更接近？
3. 結果是否由樣本量、spatial clustering、resolved closure、selection-source mixture或EB regularization造成？

### 假設與 falsifier

- 假設：若技術來源有強 aggregate composition再現，HCC pair應位於20個cross-biological pair distances的低尾端，且在rarefaction／block bootstrap／EB後保持。
- Falsifier：若有多個cross-biological pairs更近，或加入unresolved／source standardization後rank大幅改變，則只能稱moderate／adjustment-sensitive technical proximity。
- 不測試：不同cell lines是否應有相同演化；不把pattern比例當clone/subclone比例。

## 2. 輸入與資料品質

### 輸入路徑

1. `InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_four_class_summary.tsv`
2. `InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_regions.tsv`
3. `InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_four_class_by_source.tsv`
4. `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv`

### Read-only盤點結果

- 37,039 resolved region rows；sample×region unique。
- 從 complete/coarse表安全補回2,846 unresolved，得到39,885 complete rows；逐dataset與summary完全守恆。
- region均含合法`chr1-22:start-end`；7 datasets各有22 autosomes。
- Complete universe的5 Mb block：全7共同492／union569；共同block保留各sample 88.17%–92.02% complete rows。
- 10 Mb共同267／union294；保留91.71%–96.40%。
- `structural_topo1=21,976`、`vaf_resolved_topogt1=15,063`；逐類／逐sample與resolved summary一致。

## 3. 全7 dataset raw complete五類分布

百分比均以各dataset的complete regions為分母，因此五類相加100%。

| Dataset row | Complete n | Single | Sister | Direct | Sister+direct | Unresolved |
|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | 6,940 | 13.11% | 5.16% | 71.70% | 7.98% | 2.05% |
| HCC1395_DORADO | 6,750 | 19.81% | 2.28% | 65.11% | 2.90% | 9.90% |
| COLO829 | 6,949 | 18.12% | 0.26% | 61.82% | 0.39% | 19.41% |
| H1437 | 6,984 | 12.56% | 1.32% | 77.75% | 3.19% | 5.18% |
| H2009 | 5,882 | 6.90% | 2.02% | 78.17% | 9.15% | 3.76% |
| HCC1937 | 2,557 | 13.45% | 12.98% | 66.25% | 6.73% | 0.59% |
| HCC1954 | 3,823 | 23.38% | 10.10% | 56.34% | 7.85% | 2.33% |

每類的Jeffreys Dirichlet marginal 95% credible interval在`composition_posterior_intervals.tsv`。六個biological-ID composition另列於`biological_id_compositions.tsv`；HCC1395 ID採兩technical sources等權比例平均，避免把它當雙生物重複。

## 4. HCC technical pair：raw distance與相對rank

Rank定義：把HCC pair插入20個`biological_id_a != biological_id_b`的dataset-row distances；rank 1最接近、rank 21最遠。cross-biological pairs不是「應相同」null，只是經驗參照。

| Estimand | TV | JSD base2 | Aitchison | TV rank | JSD rank | Aitchison rank |
|---|---:|---:|---:|---:|---:|---:|
| Complete五類 primary | **0.1455** | **0.1984** | **2.0860** | 8 | 9 | 10 |
| Resolved四類 secondary | 0.0860 | 0.1383 | 1.1388 | 4 | 8 | 6 |

Primary JSD rank=9代表**8個cross-biological pairs比HCC technical pair更近**。所以全域composition只支持「中等technical proximity」，不支持「同cell line兩技術來源高度一致」。

Resolved-only TV看似降到8.60%，主要因先排除unresolved再closure；HCC unresolved為2.05%，DORADO為9.90%，不能忽略，故四類只列secondary。

## 5. Bulk-sampling與posterior sensitivity

固定seed=`20260713`；每個stochastic endpoint 5,000 replicates。

| Method（complete五類 primary） | HCC JSD median | 95% interval | median rank | P(rank=1) |
|---|---:|---:|---:|---:|
| Raw | 0.1984 | — | 9 | 0.000 |
| Dirichlet posterior | 0.1987 | [0.1851, 0.2116] | 9 | 0.000 |
| Equal-n rarefaction，n=2,557 | 0.1991 | [0.1814, 0.2162] | 9 | 0.000 |
| Common 5 Mb block bootstrap | 0.1976 | [0.1858, 0.2095] | 9 | 0.000 |
| Common 10 Mb block bootstrap | 0.1967 | [0.1848, 0.2092] | 9 | 0.000 |
| EB leave-HCC-out，2× concentration | 0.1971 | deterministic sensitivity | 9 | 0.000 |

判讀：raw／posterior／rarefaction／5 Mb／10 Mb／EB的rank均為9；樣本數不同與genomic clustering不是HCC aggregate差異的主要解釋。EB從6個biological compositions或leave-HCC-out 5個估prior，0.5×／1×／2×濃度的全表最大比例位移只有0.204 percentage points。

## 6. Selection-source分層與Simpson-like mixture cancellation

Resolved regions中：

- HCC1395 source weights：structural 3,438/6,798=50.57%，VAF-resolved 49.43%。
- DORADO：structural 2,444/6,082=40.18%，VAF-resolved 59.82%。
- Conditional HCC TV：structural Topo=1 **28.24%**；VAF-resolved Topo>1 **8.17%**。
- Raw resolved mixture TV：**8.60%**。
- 用全7 resolved pooled共同weights（structural=59.33%，VAF-resolved=40.67%）direct-standardize：TV **16.75%**，為raw的1.95×；JSD=0.1867，relative rank=10。

這是明確的**mixture cancellation／Simpson-like sensitivity**：兩來源mix不同，且conditional compositions不同，raw marginal差異被部分抵銷。但source是方法結果／selection state，並非已證明的外生confounder；direct standardization可能over-adjust，所以16.75%不能稱「校正後真差異」。Primary complete三來源共同weight sensitivity TV=15.56%，同樣只作邊界分析。

## 7. 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260713_crosssample_structure_bulk_sampling_adjustment/scripts/analyze_crosssample_composition.py \
  --summary research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_four_class_summary.tsv \
  --regions research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_regions.tsv \
  --by-source research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_four_class_by_source.tsv \
  --coarse-regions research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv \
  --output-dir research/20260713_crosssample_structure_bulk_sampling_adjustment/results/rerun_a \
  --seed 20260713 --replicates 5000
```

第二次完整執行把`--output-dir`改為`results/rerun_b`；其餘完全相同。

### 實際輸出片段

```text
STATUS -> PASS_WITH_CAVEATS (81/81 checks)
HCC primary five-class -> TV=0.145452 JSD=0.198380
HCC resolved source sensitivity -> raw_TV=0.085966 global_standardized_TV=0.167533
```

### Determinism驗證

```text
STATUS -> PASS (18/18 byte-identical files)
```

驗證器逐相對檔名比對SHA-256，涵蓋stochastic summaries、TSV、JSON、report與provenance；receipt在`results/reproducibility/`。

## 8. 輸出索引

- 主機讀摘要：`results/rerun_a/summary.json`
- 人讀分析：`results/rerun_a/analysis_report.md`
- 全7 posterior intervals：`results/rerun_a/composition_posterior_intervals.tsv`
- Raw dataset／biological-ID pairwise：`results/rerun_a/raw_pairwise_distances.tsv`
- Resampling composition／pairwise：`results/rerun_a/resampled_{compositions,pairwise_distances}.tsv`
- HCC rank全敏感度：`results/rerun_a/hcc_pair_relative_rank.tsv`
- Source conditional／standardized：`results/rerun_a/source_{stratified,standardized}_pairwise_distances.tsv`
- EB：`results/rerun_a/empirical_bayes_{priors,shrinkage,pairwise_distances}.tsv`
- Block feasibility：`results/rerun_a/block_bootstrap_feasibility.tsv`
- QA／warnings／provenance：`results/rerun_a/{checks,warnings}.tsv`、`provenance.json`
- 兩次重跑hash receipt：`results/reproducibility/`
- 圖片輸出：無；本子任務指定TSV／JSON／執行紀錄，所有數值關係保留為可重算表格，未另製圖。

## 9. Verdict與claim ceiling

### Overall assessment：Share with caveats

1. **可支持**：四／五類graph-pattern composition可全7樣本重算；HCC technical pair在多種sampling adjustments下呈中等而非例外的composition proximity；source-mixture會明顯改變marginal距離。
2. **不可支持**：兩HCC資料高度一致、同一真實演化樹、方法已由technical pair證實有效、跨cell-line應有相同演化。
3. **無法由本輪單獨區分**：HCC差異究竟來自真region universe差異、basecalling／coverage／resolution mechanism或tree方法本身；需要fresh layered-v3及matched exact-region／orthogonal truth再分解。

> Historical layered-v2 VAF-selected rooted-unlabeled mutation-state graph-pattern composition；不是confirmed clone/subclone truth。fresh layered-v3 7/7 scientific gate仍未完成，因此publication-scope biological claim維持NO-GO。
