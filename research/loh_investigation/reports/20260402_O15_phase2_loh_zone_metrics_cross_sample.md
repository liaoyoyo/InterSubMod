<!--
建立時間: 2026-04-02 18:50:57
目標: O15 Phase 2 — Cross-sample LOH zone metrics analysis (binary inside/outside)
處理範圍: 7 samples x 2 modes, 748,391 total rows
關聯檔案:
  - scripts/analysis/build_observation_O15b_loh_zone_metrics_cross_sample.py
  - research/loh_investigation/data/o15_p2_cross_sample_auc.tsv
  - research/loh_investigation/data/o15_p2_cross_sample_effects.tsv
  - research/loh_investigation/data/o15_p2_concordance.tsv
-->

# O15 Phase 2: Cross-Sample LOH Zone Metrics Analysis

## Overview

This analysis extends Phase 1 (HCC1395 SEQC2 4-class LOH zones) to ALL 7 samples
using **binary LOH classification** from each sample's LongPhase-TO LOH.bed file.

- **Samples**: HCC1395, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954
- **LOH Source**: Per-sample LongPhase-TO tumor_phased_LOH.bed
- **Classification**: Binary (inside / outside LOH zone)
- **Total rows**: 748,391

## LOH Zone Distribution (TO Mode)

```
loh_zone        inside  outside
sample                         
COLO829          10582    40035
H1437            15531    43384
H2009            34141   103554
HCC1395          17555    22541
HCC1395_DORADO   17785    22643
HCC1937          12419    12236
HCC1954           4205    63081
```

## Key Findings

### 1. Cross-Sample Generalizability of Top Metrics (TO, Inside LOH)

Metrics with AUC > 0.6 across samples:

| Metric | # Samples with AUC > 0.6 |
|--------|--------------------------|
| AlleleDelta | 6 |
| LabelAllelePermanovaF | 6 |
| caller_ad_ref | 6 |
| HPFineF | 5 |
| HeuristicScore | 4 |

### 2. Phase 1 vs Phase 2 Concordance (HCC1395)

  - inside: Pearson R = 0.993 (n=30)
  - outside: Pearson R = 0.997 (n=30)


### 3. LOH Coverage Variation

| Sample | LOH Coverage (Mb) |
|--------|-------------------|
| HCC1395 | 1632 |
| HCC1395_DORADO | 1634 |
| COLO829 | 927 |
| H1437 | 1156 |
| H2009 | 1175 |
| HCC1937 | 1840 |
| HCC1954 | 390 |

## Figures

### Fig01: HP_Ratio Violin (7x2)

![Fig01](../figures/o15_p2_fig01_hp_ratio_violin.png)

### Fig02: Extreme HP Rate Heatmap

![Fig02](../figures/o15_p2_fig02_extreme_hp_heatmap.png)

### Fig03: Methylation Metrics Violin (TO)

![Fig03](../figures/o15_p2_fig03_methyl_metrics_violin.png)

### Fig04: AUC Heatmap (TO-inside)

![Fig04](../figures/o15_p2_fig04_auc_heatmap.png)

### Fig05: VerificationClass Bar

![Fig05](../figures/o15_p2_fig05_verification_class_bar.png)

### Fig06: Cohen's d Heatmap

![Fig06](../figures/o15_p2_fig06_cohend_heatmap.png)

### Fig07: AlleleDelta vs caller_af Scatter

![Fig07](../figures/o15_p2_fig07_allele_delta_vs_af.png)

### Fig08: PassedGating Rate Bar

![Fig08](../figures/o15_p2_fig08_passed_gating_rate.png)

### Fig09: LOH Coverage vs Delta-AUC

![Fig09](../figures/o15_p2_fig09_loh_coverage_vs_delta_auc.png)

### Fig10: Forest Plot (Top 10)

![Fig10](../figures/o15_p2_fig10_forest_plot.png)

### Fig11: Phase1 vs Phase2 Concordance

![Fig11](../figures/o15_p2_fig11_p1_vs_p2_concordance.png)

### Fig12: Methylation Significance Violin

![Fig12](../figures/o15_p2_fig12_methyl_significance.png)

### Fig13: CpG/Stability Box

![Fig13](../figures/o15_p2_fig13_cpg_stability.png)

### Fig14: VCF Quality

![Fig14](../figures/o15_p2_fig14_vcf_quality.png)

### Fig15: Allele Balance Heatmap

![Fig15](../figures/o15_p2_fig15_allele_balance.png)

### Fig16: Summary Table

![Fig16](../figures/o15_p2_fig16_summary_table.png)

## Conclusions

1. **LOH zone effects are sample-dependent**: The magnitude of LOH impact on ISM metrics
   varies substantially across samples, driven largely by LOH coverage differences.

2. **Phase 1 vs Phase 2 agreement**: Binary LOH.bed classification produces comparable
   AUC rankings to SEQC2 4-class zones for HCC1395, validating the simpler approach.

3. **TO mode is most affected**: LOH zones consistently show larger metric distortions
   in tumor-only mode compared to paired mode, consistent with self-phasing artifacts.

4. **Generalizability assessment**: Metrics that discriminate well inside LOH zones for
   HCC1395 do not necessarily generalize to all samples, particularly those with low
   LOH coverage (HCC1954) or different biological characteristics.
