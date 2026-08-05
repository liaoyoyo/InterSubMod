<!--
建立時間: 2026-07-13
目標: 全 7 dataset structure composition 的 bulk-sampling／compositional consistency adjustment
處理範圍: chr1-22；7 dataset rows；6 biological IDs；historical layered-v2 engineering snapshot
關聯檔案: InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/
-->

# Cross-sample structure composition adjustment

> **PASS WITH CAVEATS。Primary 是 complete 五類（含 unresolved）；resolved-only 四類只作 secondary。這是 technical compositional proximity，不是 clone-tree truth。**

## 研究問題與設計

- 問題：HCC1395 technical pair 的 pattern composition 是否比不同 biological IDs 更接近；結果對 closure、樣本量、genomic blocks、source mixture 與 EB shrinkage 是否敏感？
- 7 dataset rows／6 biological IDs分開；HCC pair 不算兩個 biological replicates。
- 不同 cell lines 不被強迫相同；20 個 cross-biological dataset-row pairs只作 empirical reference。
- 固定 seed `20260713`；每個 stochastic endpoint `5,000` replicates。

## Raw HCC technical-pair composition distance

| Estimand | TV | JSD (base 2) | Aitchison | JSD rank vs 20 cross-bio pairs | lower-tail percentile |
|---|---:|---:|---:|---:|---:|
| Complete five-class **primary** | 0.1455 | 0.1984 | 2.0860 | 9/21-scale | 40.0% |
| Resolved four-class secondary | 0.0860 | 0.1383 | 1.1388 | 8/21-scale | 35.0% |

Complete五類把 unresolved 2.05% vs 9.90% 保留，所以 HCC TV=14.55%；若先丟掉 unresolved 再 closure，TV降為 8.60%。後者不可當 primary。

Primary JSD的 rank=9 表示把 HCC technical pair插入20個 cross-biological dataset-row distances後，**有8個跨生物 pairs反而更近**；不是「technical pair最接近」。

## Adjustment stability（primary complete五類，JSD）

| Method | HCC JSD median | 95% interval | median rank | P(rank=1) |
|---|---:|---:|---:|---:|
| Raw | 0.1984 | — | 9 | 0.000 |
| Dirichlet posterior | 0.1987 | [0.1851, 0.2116] | 9 | 0.000 |
| Equal-n rarefaction | 0.1991 | [0.1814, 0.2162] | 9 | 0.000 |
| Common-block bootstrap 5 Mb | 0.1976 | [0.1858, 0.2095] | 9 | 0.000 |
| Common-block bootstrap 10 Mb | 0.1967 | [0.1848, 0.2092] | 9 | 0.000 |

**結論：bulk sample size與spatial resampling不是主要解釋。** Raw／posterior／rarefaction／5 Mb／10 Mb的HCC primary JSD與rank幾乎不動；EB最大比例位移也只有 0.204 percentage points。aggregate composition呈現中等、非例外的technical proximity，而不是高度一致。

## Source-mixture sensitivity

- HCC resolved raw TV：8.60%。
- 分 source 看，structural Topo=1 conditional TV=28.24%；VAF-resolved Topo>1 conditional TV=8.17%。
- 用全 7 resolved regions 的共同 source weights（structural=59.33%、VAF-resolved=40.67%）direct-standardize後：TV=16.75%（1.95×）。
- 這是 **Simpson-like mixture cancellation** 診斷：raw source mix部分抵銷 conditional差異；不是「校正後真值」。source是 outcome-dependent selection state，標準化可能 over-adjust。
- Complete三來源共同權重 sensitivity TV=15.56%；同樣不能取代 raw primary。

## Sampling adjustment

- Dirichlet posterior：Jeffreys `alpha=0.5/category`，輸出每類 marginal 95% credible interval與 posterior pair distances。
- Equal-n rarefaction：primary n=2,557；secondary n=2,542，無放回抽樣。
- Spatial bootstrap：region具真實 chr/start/end；共同 5 Mb blocks=492/569 union，另做10 Mb sensitivity。共同 blocks以相同 multiplicity跨7 rows聯合重抽，但不是exact-region配對。
- EB：先把 HCC pair等權合成一個 biological composition，再估 all-6 與 leave-HCC-out-5 priors；0.5×/1×/2× concentration只作 regularization sensitivity。最大比例位移=0.204 percentage points。

## Claim ceiling

Historical layered-v2 VAF-selected rooted-unlabeled mutation-state graph-pattern composition; technical sampling sensitivity only, not clone/subclone truth, biological equivalence, or a validated tree.

fresh layered-v3尚未7/7 scientific gate；任何「相似演化」「同一 clone tree」「方法證明有效」均超出本輪證據。完整 sensitivity數值見 `hcc_pair_relative_rank.tsv`、`source_standardized_pairwise_distances.tsv` 與 `summary.json`。
