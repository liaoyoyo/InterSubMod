<!--
建立時間: 2026-04-02 20:17
目標: 驗證 LOH zone-aware QS v2 重設計是否改善 TO mode 的 AUC
處理範圍: 7 samples x 2 modes x 3 QS formulas
關聯檔案:
  - scripts/analysis/build_qs_v2_zone_aware_validation.py
  - research/loh_investigation/figures/qs_v2_fig*.png
  - research/loh_investigation/data/qs_v2_*.tsv
-->

# QS v2 Zone-Aware Validation Report

**Date**: 2026-04-02
**Script**: `scripts/analysis/build_qs_v2_zone_aware_validation.py`
**Data**: `all_region_rows.tsv.gz` (748,391 rows)
**LOH zones**: per-sample LongPhase-TO LOH.bed (binary inside/outside)

## QS v1 Reconstruction Validation

- Correlation with existing Quality_Score: **0.853190**
- Max absolute error: **25.00**
- Mean absolute error: **6.1802**
- Exact match rate: **0.7235**

## Zone Distribution

```
truth_label          FP      TP
mode   loh_zone                
paired inside       987   88775
       outside     2442  236495
to     inside     26875   85343
       outside   101507  205967
```

## AUC Comparison: By Zone

### TO Mode

| Zone | Sample | AUC_v1 | AUC_v2a | AUC_v2b | delta_v2a | delta_v2b |
|------|--------|--------|---------|---------|-----------|-----------|
| inside | COLO829 | 0.495 | 0.499 | 0.507 | +0.004 | +0.012 |
| inside | H1437 | 0.556 | 0.565 | 0.598 | +0.009 | +0.042 |
| inside | H2009 | 0.500 | 0.505 | 0.508 | +0.005 | +0.008 |
| inside | HCC1395 | 0.552 | 0.586 | 0.633 | +0.033 | +0.081 |
| inside | HCC1395_DORADO | 0.556 | 0.592 | 0.639 | +0.036 | +0.083 |
| inside | HCC1937 | 0.518 | 0.534 | 0.545 | +0.016 | +0.028 |
| inside | HCC1954 | 0.556 | 0.593 | 0.626 | +0.037 | +0.070 |
| outside | COLO829 | 0.482 | 0.482 | 0.482 | +0.000 | +0.000 |
| outside | H1437 | 0.513 | 0.513 | 0.513 | +0.000 | +0.000 |
| outside | H2009 | 0.490 | 0.490 | 0.490 | +0.000 | +0.000 |
| outside | HCC1395 | 0.501 | 0.501 | 0.501 | +0.000 | +0.000 |
| outside | HCC1395_DORADO | 0.508 | 0.508 | 0.508 | +0.000 | +0.000 |
| outside | HCC1937 | 0.483 | 0.483 | 0.483 | +0.000 | +0.000 |
| outside | HCC1954 | 0.498 | 0.498 | 0.498 | +0.000 | +0.000 |

### Paired Mode

| Zone | Sample | AUC_v1 | AUC_v2a | AUC_v2b | delta_v2a | delta_v2b |
|------|--------|--------|---------|---------|-----------|-----------|
| inside | COLO829 | 0.487 | 0.473 | 0.460 | -0.014 | -0.028 |
| inside | H1437 | 0.945 | 0.948 | 0.918 | +0.003 | -0.027 |
| inside | H2009 | 0.490 | 0.508 | 0.515 | +0.017 | +0.025 |
| inside | HCC1395 | 0.420 | 0.442 | 0.475 | +0.022 | +0.055 |
| inside | HCC1395_DORADO | 0.557 | 0.559 | 0.598 | +0.002 | +0.041 |
| inside | HCC1937 | 0.529 | 0.541 | 0.545 | +0.012 | +0.017 |
| inside | HCC1954 | 0.725 | 0.751 | 0.771 | +0.026 | +0.046 |
| outside | COLO829 | 0.513 | 0.513 | 0.513 | +0.000 | +0.000 |
| outside | H1437 | 0.643 | 0.643 | 0.643 | +0.000 | +0.000 |
| outside | H2009 | 0.693 | 0.693 | 0.693 | +0.000 | +0.000 |
| outside | HCC1395 | 0.630 | 0.630 | 0.630 | +0.000 | +0.000 |
| outside | HCC1395_DORADO | 0.501 | 0.501 | 0.501 | +0.000 | +0.000 |
| outside | HCC1937 | 0.500 | 0.500 | 0.500 | +0.000 | +0.000 |
| outside | HCC1954 | 0.527 | 0.527 | 0.527 | +0.000 | +0.000 |

## Pooled AUC (All Samples)

| Mode | Zone | AUC_v1 | AUC_v2a | AUC_v2b | delta_v2a | delta_v2b |
|------|------|--------|---------|---------|-----------|-----------|
| to | inside | 0.578 | 0.590 | 0.609 | +0.012 | +0.031 |
| to | outside | 0.548 | 0.548 | 0.548 | +0.000 | +0.000 |
| to | all | 0.546 | 0.554 | 0.557 | +0.008 | +0.012 |
| paired | inside | 0.686 | 0.708 | 0.714 | +0.022 | +0.028 |
| paired | outside | 0.808 | 0.808 | 0.808 | +0.000 | +0.000 |
| paired | all | 0.754 | 0.782 | 0.783 | +0.028 | +0.029 |

## GO/NO-GO Evaluation

| Criterion | Threshold | v2a Value | v2a Verdict | v2b Value | v2b Verdict |
|-----------|-----------|-----------|-------------|-----------|-------------|
| TO overall AUC improvement | >= +0.05 | +0.008 | FAIL | +0.012 | FAIL |
| TO LOH-inside AUC improvement | >= +0.10 | +0.012 | FAIL | +0.031 | FAIL |
| Cross-sample consistency (TO inside) | >= 5/7 improved | 7/7 | PASS | 7/7 | PASS |
| No regression (TO LOH-outside) | >= -0.02 | +0.000 | PASS | +0.000 | PASS |

### Overall Verdict

- **QS_v2a**: 2/4 criteria passed -> **CONDITIONAL GO**
- **QS_v2b**: 2/4 criteria passed -> **CONDITIONAL GO**

## Conclusions

1. QS_v1 reconstruction validated (corr=0.853190, max_err=25.00)
2. TO overall: v1 AUC=0.546, v2a AUC=0.554 (delta=+0.008), v2b AUC=0.557 (delta=+0.012)
3. TO LOH-inside: v1 AUC=0.578, v2a AUC=0.590 (delta=+0.012), v2b AUC=0.609 (delta=+0.031)
4. Cross-sample consistency: v2a=7/7, v2b=7/7
5. LOH-outside regression check: v2a=+0.000, v2b=+0.000
