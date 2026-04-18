<!--
建立時間: 2026-04-02 18:27
目標: O15 Phase 1 — HCC1395 LOH Zone Metrics Complete Analysis
處理範圍: HCC1395 samples × SEQC2 4-class LOH zone × all ISM metrics
關聯檔案:
  - scripts/analysis/build_observation_O15_loh_zone_metrics_hcc1395.py
  - research/loh_investigation/data/seqc2_vs_to_validation/loh_classified_*.bed
  - research/loh_investigation/figures/o15_p1_fig*.png
  - research/loh_investigation/data/o15_p1_*.tsv
-->

# O15 Phase 1: HCC1395 LOH Zone Metrics Complete Analysis

**Generated**: 2026-04-02 18:27
**Total rows**: 141,014
**Samples**: HCC1395, HCC1395_DORADO
**Modes**: paired, to

## Zone Distribution

```
truth_label         FP     TP
mode   loh_zone              
paired FN          471   1246
       FP            8    626
       TN          150  31438
       TP          238  26313
to     FN          376   1208
       FP          180    566
       TN        13699  29945
       TP         8918  25632
```

## Key Findings

### Top AUC Metrics (TP/FP discrimination)

| Metric | Mode | LOH Group | AUC | n_TP | n_FP |
| --- | --- | --- | --- | --- | --- |
| caller_gq | paired | inside | 0.922 | 26939 | 246 |
| caller_ad_alt | paired | outside | 0.870 | 32684 | 621 |
| caller_gq | paired | outside | 0.860 | 32684 | 621 |
| caller_af | paired | outside | 0.848 | 32684 | 621 |
| caller_ad_ref | to | inside | 0.784 | 26198 | 9098 |
| LabelAllelePermanovaF | to | inside | 0.757 | 26198 | 9098 |
| HPFineF | to | inside | 0.710 | 26198 | 9098 |
| caller_ad_alt | paired | inside | 0.697 | 26939 | 246 |
| AlleleDelta | to | inside | 0.696 | 26198 | 9098 |
| HP_Ratio | paired | outside | 0.656 | 32684 | 621 |

### Top Effect Sizes (Cohen's d, TP vs FP)

| Metric | Mode | LOH Group | Cohen's d | 95% CI |
| --- | --- | --- | --- | --- |
| caller_gq | paired | inside | 2.004 | [1.899, 2.121] |
| caller_gq | paired | outside | 1.520 | [1.450, 1.596] |
| caller_af | to | inside | -1.204 | [-1.224, -1.182] |
| AlleleDelta | paired | outside | -1.183 | [-1.321, -1.033] |
| LabelAllelePermanovaP | to | inside | -1.093 | [-1.112, -1.071] |
| caller_af | paired | outside | 1.064 | [0.997, 1.129] |
| AlleleP | to | inside | -1.029 | [-1.049, -1.008] |
| caller_ad_alt | paired | outside | 0.962 | [0.896, 1.031] |
| caller_ad_alt | to | inside | -0.951 | [-0.977, -0.923] |
| caller_ad_ref | to | inside | 0.936 | [0.920, 0.952] |

## Figures

### Fig01: HP_Ratio by LOH Zone × Mode (Violin)

![Fig01](../figures/o15_p1_fig01_hp_ratio_violin.png)

### Fig02: Extreme HP_Ratio Rate by Zone × Mode × Truth

![Fig02](../figures/o15_p1_fig02_extreme_hp_ratio.png)

### Fig03: ISM Methylation Metrics by Zone × Mode

![Fig03](../figures/o15_p1_fig03_methyl_metrics_violin.png)

### Fig04: AUC Heatmap — TP/FP Discrimination

![Fig04](../figures/o15_p1_fig04_auc_heatmap.png)

### Fig05: VerificationClass by Zone × Truth (TO Only)

![Fig05](../figures/o15_p1_fig05_verification_class.png)

### Fig06: Quality_Score by Zone × Mode × Truth

![Fig06](../figures/o15_p1_fig06_quality_score.png)

### Fig07: AlleleDelta vs caller_af (LOH-inside vs outside, TO)

![Fig07](../figures/o15_p1_fig07_allele_delta_vs_af.png)

### Fig08: Forest Plot — Cohen's d per Metric (TO)

![Fig08](../figures/o15_p1_fig08_forest_plot.png)

### Fig09: PassedGating PASS Rate by Zone × Mode × Truth

![Fig09](../figures/o15_p1_fig09_passed_gating.png)

### Fig10: LOH_Subtype × SEQC2 Zone Cross-Tabulation (TO)

![Fig10](../figures/o15_p1_fig10_loh_subtype_heatmap.png)

### Fig11: PERMANOVA F-stats by Zone × Mode

![Fig11](../figures/o15_p1_fig11_permanova_f_stats.png)

### Fig12: Methylation Significance by Zone × Mode

![Fig12](../figures/o15_p1_fig12_methyl_significance.png)

### Fig13: CpG/Stability/HeuristicScore by Zone × Mode × Truth

![Fig13](../figures/o15_p1_fig13_cpg_stability_heuristic.png)

### Fig14: VCF FILTER & caller_af Distribution

![Fig14](../figures/o15_p1_fig14_vcf_filter_af.png)

### Fig15: Caller Read Counts by Zone × Mode

![Fig15](../figures/o15_p1_fig15_caller_read_counts.png)

### Fig16: Summary Table — All Metrics × Zone × Mode

![Fig16](../figures/o15_p1_fig16_summary_table.png)

## Data Files

- `../data/o15_p1_statistical_summary.tsv` — all metrics × zone × mode: median, IQR, n
- `../data/o15_p1_auc_by_stratum.tsv` — AUC for each metric × zone × mode
- `../data/o15_p1_effect_sizes.tsv` — Cohen's d for TP vs FP within each zone × mode

## Conclusions

1. LOH zones (SEQC2 TP/FP/FN/TN) create distinct metric distributions across all ISM features.
2. The AUC heatmap reveals which metrics retain discriminative power inside vs outside LOH regions.
3. Cohen's d forest plot shows effect size differences are generally attenuated inside LOH zones.
4. VerificationClass distributions shift toward Noise/Weak inside LOH-TP zones for TO mode.
5. Quality_Score is systematically lower inside LOH zones due to LOH penalty in the scoring formula.
