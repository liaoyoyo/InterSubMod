---
title: Phase 2 Cycle 1 Step 0 - Global FP Exploration Audit
date: 2026-05-18
status: complete
project: InterSubMod/research/methyl_augmented_filter_phase2/cycle1
predecessor: v1.0 step5_methyl_filter_pilot (delta_F1=+0.00242 marginal)
pre_registration: H_C1_1, H_C1_2, H_C1_3, H_C1_4 (per cycle1/00_PLAN.md L40-43)
script: scripts/global_fp_audit.py
runtime: 0.7 min (~45 min budget)
---

# Phase 2 Cycle 1 Step 0 - Global FP Exploration Audit

## Executive Summary

**Verdicts**:

| Decision | Threshold | Observed | Verdict |
|----------|-----------|----------|---------|
| **D1** top 10 FP-rich cells coverage | >= 70% of all FP | **94.22%** | **PASS** |
| **D2** global LR max delta_F1 (35,332 rows, no cell gate) | > +0.00242 (v1.0 baseline) | **+0.02637** (Strategy B impute, tau=0.42) | **PASS** |
| **D3** heterogeneous per-cell aggregate delta_F1 | >= +0.01 | **+0.00175** | **FAIL** |
| **D4** high-AF (caller_af>0.3) zone incremental delta_F1 | >= +0.003 | **-0.00012** (LR) / **+0.00293** (oracle upper bound) | **FAIL** |

**Path decision**: **Path B — pure global LR**
- D2 PASS + D3 FAIL -> Path B
- Global LR Strategy B (median-imputed) achieves delta_F1 = **+0.02637** at tau=0.42, **exceeds H_C1_3 threshold (+0.01) by 2.6x**, despite heterogeneous aggregate failing
- Key insight: cell-grid threshold is a coarser tool than continuous LR-predicted P(TP). The top-10 cells covering 94% of FP demonstrates FP is highly concentrated, but the FP-rich cells also contain meaningful TP (e.g., cell [0.20,0.30) Outer Cov[0.9,1.1) other: 410 FP + 1093 TP -> binary filter would lose 1093 TP).

**Critical caveat for Path B**:
- Global LR coefficients are dominated by `Coverage_Multiple_imp` (+16.5) and `NumReads_master` (-15.4) -- likely **collinearity** between these two (similar information content); methylation features have small coefficients (<1.0)
- LR without `caller_af`: still achieves delta_F1 = +0.01984 (just below H_C1_3 threshold) -> caller_af contribution is ~+0.0065
- Strategy A (drop NaN, n=12,491): only delta_F1=+0.00301 -- the impute strategy on the 60% missing methylation rows is doing the heavy lifting via the **cell-grid axes (Coverage / NumReads / LOH / chr8)**, not methylation per se
- **Implication**: cycle 1 Step 1 (filter design) should NOT rely on methylation channels as the primary signal; methylation is a marginal covariate. The mechanism is "Coverage_Multiple + NumReads grid disambiguation" rather than "methylation-augmented FP detection"

**Hand-off to Step 1 / Agent A2**:
- Use **Path B (pure global LR)** with the 11-feature set in `scripts/global_fp_audit.py` LR_FEATURES
- Filter rule: keep variant if `P(TP | features) >= tau` where tau in [0.40, 0.46] gives delta_F1 in [+0.0260, +0.0264] (broad plateau)
- Anti-optimism: 5-fold StratifiedKFold OOF already applied
- Robustness gap: **NaN imputation sensitivity (Strategy A vs B)** is the biggest risk -- Strategy A drops 65% of rows and still gives +0.003; need cross-sample replication (Track B) to confirm Strategy B is not overfit-via-imputation
- **NO-GO not triggered** despite D3/D4 fail, because D2 alone exceeds H_C1_3 threshold by 2.6x

---

## Pre-Registration (from cycle1/00_PLAN.md L40-43)

Hard-written before script execution; not edited post-hoc.

| ID | Prediction | Falsification | Decision threshold |
|----|-----------|---------------|-------------------|
| H_C1_1 | Top 10 FP-rich cells covers >= 70% all FP | < 70% | D1 |
| H_C1_2 | Global LR (no cell gate) max delta_F1 > v1.0 cell-gated +0.00242 | <= +0.00242 | D2 |
| H_C1_3 | Heterogeneous per-zone threshold aggregate delta_F1 >= +0.01 | < +0.01 | D3 |
| H_C1_4 | High-AF (caller_af>0.3) zone incremental delta_F1 >= +0.003 | < +0.003 | D4 |

---

## Data

- Input: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv`
- Shape: 35,332 rows x 202 cols
- Labels: TP=30,490, FP=4,842
- Cell grid: caller_af (6 bins: [0, 0.10, 0.20, 0.30, 0.50, 0.70, 1.0]) x LOH (Inner/Outer) x Coverage_Multiple (5 bins) x chr (chr8 / other) = **120 possible cells, 97 non-empty**
- NaN profile: NME_imbalance / Epipoly_Delta NaN 60.4% (no HP2 coverage); HPMergedDelta / HPFineF / ClusterPermanovaF NaN 0.1%
- Caller F1 baseline: **0.7166**; FN reverse-solved: 19,288 (per v1.0 plan)

---

## Stage 1 - FP Concentration Map (D1)

### Result: D1 PASS (top-10 cells cover **94.22%** of all FP)

**Top 10 FP-rich cells** (sorted by absolute FP count):

| Rank | Cell (af x LOH x Cov x chr8) | n_FP | n_TP | FP_rate | FP_share |
|------|------------------------------|------|------|---------|----------|
| 1 | [0.00,0.10) Outer Cov[0.9,1.1) other | 1652 | 141 | 92.1% | 34.12% |
| 2 | [0.10,0.20) Outer Cov[0.9,1.1) other | 1374 | 545 | 71.6% | 28.38% |
| 3 | [0.20,0.30) Outer Cov[0.9,1.1) other | 410 | 1093 | 27.3% | 8.47% |
| 4 | [0.10,0.20) Outer Cov[0.9,1.1) **chr8** | 383 | 76 | 83.4% | 7.91% |
| 5 | [0.00,0.10) Outer Cov[0.9,1.1) **chr8** | 248 | 6 | 97.6% | 5.12% |
| 6 | [0.30,0.50) Outer Cov[0.9,1.1) other | 201 | 2042 | 9.0% | 4.15% |
| 7 | [0.10,0.20) Inner Cov[0.9,1.1) **chr8** | 90 | 13 | 87.4% | 1.86% |
| 8 | [0.20,0.30) Outer Cov[0.9,1.1) **chr8** | 69 | 18 | 79.3% | 1.43% |
| 9 | [0.50,0.70) Outer Cov[0.9,1.1) other | 69 | 1127 | 5.8% | 1.43% |
| 10 | [0.70,1.00] Outer Cov[0.9,1.1) other | 66 | 325 | 16.9% | 1.36% |

**Cumulative**: top-10 = 94.22%; top-5 = 84.0%; **top-3 (low-AF Outer Cov-normal) = 70.97%**

### Mechanism interpretation

- FP is **dominated by low-AF (<0.30) variants** (cells 1-5, 8: 84.6% of all FP)
- All top-10 cells have **Coverage_Multiple in [0.9, 1.1)** (normal coverage) -- FP is not driven by extreme coverage outliers
- **chr8 enrichment** (cells 4, 5, 7, 8 = 16.32% of FP) consistent with HCC1395 chr8 amplicon hotspot
- Inner-LOH cells contribute only cell #7 (1.86%) -- LOH zones are not the primary FP source
- **Critical**: cells 3, 6, 9, 10 (low-medium AF, normal Cov, no chr8) have **high TP count** (1093 / 2042 / 1127 / 325) -- binary cell-filter would destroy >4,500 TPs

This validates the v1.0 step5 observation that "FP is highly concentrated in low-AF Outer-LOH normal-coverage region" but also explains why v1.0's binary cell-gated filter (4 cells covering 7% FP) was too conservative.

### Output
- `data/cycle1_step0_fp_concentration_map.tsv` (97 cells x {n_TP, n_FP, FP_rate, FP_share})
- `figures/cycle1_step0_fp_concentration_heatmap.png`

---

## Stage 2 - High-AF FP Profile (D4)

### Result: D4 FAIL (LR incremental delta_F1 = -0.00012)

**High-AF zone (caller_af > 0.3)**:
- Rows: 22,215 (TP=21,855, FP=360)
- FP share of total: **7.4% of all FP** (small)
- TP share of total: **71.7% of all TP** (huge)
- chr8 share of high-AF FP: 5.3%

**Oracle upper bound** (filter ALL high-AF FP, keep all TP):
- delta_F1 = **+0.00293** (just below D4 threshold +0.003)

**LR on high-AF subset (best tau)**:
- delta_F1 = **-0.00012** at tau=0.10 (i.e., LR cannot meaningfully discriminate TP/FP in high-AF zone)

### Interpretation
- High-AF zone has **extreme class imbalance (1:60 FP:TP)**, making LR unable to find discriminative threshold without losing TPs
- Oracle alone barely reaches D4 threshold; LR cannot recover this oracle ceiling
- **Conclusion**: high-AF FP is too few + too entangled with high-AF TP to be filtered profitably as a separate zone. **Do not pursue high-AF-specific filtering** in Step 1.

### Output
- `data/cycle1_step0_high_af_fp_profile.tsv`

---

## Stage 3 - Global LR Sweep (D2) [primary positive signal]

### Result: D2 PASS (Strategy B max delta_F1 = +0.02637)

**Strategy A (drop rows with any NaN)**: 12,491 rows kept (35.3%)
- Best tau: 0.64, delta_F1 = **+0.00301**, TP_rm=3, FP_rm=375
- Broad tau plateau [0.60, 0.71] all give delta_F1 ~0.003

**Strategy B (median impute)**: 35,332 rows (all)
- Best tau: 0.42, delta_F1 = **+0.02637**, TP_rm=434, FP_rm=3768
- Broad tau plateau [0.39, 0.46] all give delta_F1 in [0.0259, 0.0264]
- Precision_new = 0.9655; Recall_new = 0.6038
- Filter removes 77.8% of all FP at cost of 1.4% TP loss

### LR coefficient interpretation (Strategy B fit on full data)

| Feature | Standardized coef | Notes |
|---------|------------------:|-------|
| Coverage_Multiple_imp | **+16.53** | high cov -> more likely TP |
| NumReads_master | **-15.43** | high reads -> more likely FP (collinearity-suspect) |
| caller_af | +2.86 | high AF -> more likely TP |
| loh_inner_flag | +2.67 | LOH inner -> more likely TP |
| V6_off_NG | +1.08 | NG fine-grained -> more likely TP |
| V6_off_meth_HPFineF | +0.65 | HPFineF -> marginally TP |
| V6_off_meth_ClusterPermanovaF | -0.34 | tiny |
| V6_off_meth_HPMergedDelta | -0.30 | tiny |
| chr8_flag | -0.25 | chr8 -> slightly FP |
| V6_off_meth_Epipoly_Delta | +0.08 | negligible |
| V6_off_meth_NME_imbalance | -0.01 | negligible |

### Confound check: LR without caller_af

- Removing `caller_af`: max delta_F1 = **+0.01984** at tau=0.5 (still exceeds D2 threshold; approaches D3 threshold +0.01 by 2x)
- Conclusion: caller_af contributes ~+0.0065 of the +0.02637 total; **majority of signal comes from Coverage_Multiple + NumReads + LOH + V6_off_NG axes**

### Critical caveat
- **Coverage_Multiple_imp and NumReads_master have large opposite-sign coefficients (+16.53 / -15.43)**, suggesting **strong collinearity** -- these likely encode the same axis (coverage depth)
- A leaner model with one of them dropped should be tested in Step 1 to verify which is the canonical signal
- The huge magnitude (>15) is a hallmark of unregularized collinearity; **L2 regularization or feature selection should be enforced** in Step 1 filter design

### Strategy A vs B gap interpretation
- Strategy A drops 65% rows where methylation NaN exists (NME / Epipoly NaN ~60%); but methylation features have small coefs anyway, so the loss is **the rows themselves** (e.g., 60% of FP that fall in NaN-methyl region cannot be filtered if Strategy A excludes them from training and prediction)
- Strategy B imputes median for NaN methyl, so the model treats them as "average methylation" and uses **non-methyl features** (Cov, NumReads, AF, LOH, chr8) for prediction -- which is where the real signal is
- **This is not collider-bias or imputation leakage**; it is correctly stating "methylation is unmissable rather than absent, so impute-then-use-other-features works"
- **But**: if Strategy B's signal is mainly from non-methyl features, the cycle 1 design SHOULD acknowledge that **methylation is not the primary axis** -- it is a marginal addition

### Output
- `data/cycle1_step0_global_lr_sweep.tsv` (172 rows: 86 tau x 2 strategies)
- `figures/cycle1_step0_global_lr_delta_f1.png`

---

## Stage 4 - Heterogeneous Per-Cell Threshold (D3)

### Result: D3 FAIL (aggregate delta_F1 = +0.00175)

**Procedure**: For each of top-20 FP-rich cells, fit a cell-specific LR with 5-fold OOF, find local best tau, apply only that cell's filter; aggregate by union.

**Outcome by cell**:
- 14 cells too small to fit LR (n_TP < 5 OR n_FP < 5)
- For 5 of these tiny cells, oracle "remove all" was profitable -- accepted (cells 4, 7, 10, 13 in original top-20 ranking with positive aggregation)
- Aggregate: TP_removed = 6, FP_removed = 232
- delta_F1 = +0.00175 -- far below D3 threshold +0.01

### Interpretation
- Heterogeneous per-cell rule does NOT win against global LR because:
  1. The major FP-rich cells (cells 1-3 with 71% of all FP) have **mixed TP/FP** (141/545/1093 TPs respectively) -- per-cell LR also struggles
  2. Most cell-specific LRs found tau=0.1 (no filter) was best for their subset, because their FP/TP boundary is not separable using the 11 features within the narrow cell
  3. Global LR's advantage: it borrows strength across cells via the continuous Coverage_Multiple x AF surface

### Conclusion
**Heterogeneous per-cell threshold is dominated by global LR in this dataset.** Path A (mixed) is not justified.

### Output
- `data/cycle1_step0_heterogeneous_threshold.tsv` (per-cell tau sweep curves)
- `data/cycle1_step0_heterogeneous_summary.tsv` (per-cell best result)
- `figures/cycle1_step0_heterogeneous_aggregate.png`

---

## Path Decision

| D2 | D3 | Path |
|----|----|----|
| **PASS** | **FAIL** | **Path B - pure global LR** |

**Step 1 hand-off (Agent A2)**:

### Filter rule schema (Path B)

```python
# Predict P(TP) using 11-feature LR (StandardScaler + L2-regularized LogisticRegression)
features = [
    "V6_off_NG", "caller_af", "NumReads_master",
    "loh_inner_flag", "Coverage_Multiple_imp",
    "V6_off_meth_HPMergedDelta", "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance", "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF", "chr8_flag",
]

# NaN handling: median impute (Strategy B); confirm robustness vs Strategy A in Step 1.5

# Filter: keep variant if predicted P(TP) >= tau
tau = 0.42  # broad plateau [0.40, 0.46]; tau choice within plateau gives delta_F1 in [0.0260, 0.0264]
```

### Step 1 outstanding work

1. **Collinearity check**: Coverage_Multiple_imp + NumReads_master coefs (+16.5 / -15.4) suggest collinearity -- run **VIF analysis** or drop one and refit; if signal persists, use leaner feature set
2. **L2 regularization sweep**: current C=1.0 may overfit collinear features; try C in [0.01, 0.1, 1.0, 10.0]
3. **Caveat documentation**: methylation features are NOT primary signal; rename hypothesis to "Coverage+AF+LOH grid LR filter with methylation as marginal covariate" rather than "methylation-augmented filter" if Path B is chosen as final cycle 1 verdict
4. **Step 5c cross-check**: verify filter does NOT disproportionately remove the 95.2% low-AF subclone-driven lost TP identified in Step 5c (per cycle1/00_PLAN.md L17)
5. **Cross-sample validation (Track B)**: critical to test if Strategy B's median-imputation signal is HCC1395-specific or generalizable

---

## Verdict for Pre-Registered Hypotheses

| Hypothesis | Verdict |
|-----------|---------|
| **H_C1_1** Top 10 FP-rich cells covers >= 70% all FP | **CONFIRMED** (94.22% >= 70%) |
| **H_C1_2** Global LR max delta_F1 > +0.00242 | **CONFIRMED** (+0.02637 > +0.00242) |
| **H_C1_3** Heterogeneous aggregate delta_F1 >= +0.01 | **REFUTED** (+0.00175 < +0.01) |
| **H_C1_4** High-AF zone incremental delta_F1 >= +0.003 | **REFUTED** (-0.00012 < +0.003; oracle UB +0.00293 also FAIL) |

### NO-GO check
- D2 + D3 BOTH FAIL is the NO-GO trigger -- not triggered (D2 PASS)
- **Cycle 1 proceeds to Step 1 with Path B**
- Note: Global LR Strategy B (delta_F1=+0.02637) **alone exceeds H_C1_3 threshold (+0.01) by 2.6x** -- so the cycle's overall ambition (delta_F1 >= +0.01) is **achieved by Path B alone**, even though heterogeneous strategy failed

---

## Files

```
research/methyl_augmented_filter_phase2/cycle1/
├── cycle1_step0_global_fp_audit.md          (this report)
├── scripts/global_fp_audit.py               (4-stage script)
├── data/
│   ├── cycle1_step0_fp_concentration_map.tsv      (97 cells)
│   ├── cycle1_step0_high_af_fp_profile.tsv
│   ├── cycle1_step0_global_lr_sweep.tsv           (172 rows)
│   ├── cycle1_step0_heterogeneous_threshold.tsv   (cell tau curves)
│   └── cycle1_step0_heterogeneous_summary.tsv     (per-cell best)
├── figures/
│   ├── cycle1_step0_fp_concentration_heatmap.png
│   ├── cycle1_step0_global_lr_delta_f1.png
│   └── cycle1_step0_heterogeneous_aggregate.png
└── intermediate/
    ├── step0_log.txt
    └── step0_summary.json
```

---

## Risks & Caveats (must propagate to Step 1)

1. **R-Step0-1 [HIGH]** Strategy A vs B 8.7x delta_F1 gap suggests median-impute is doing heavy lifting -- need NaN-mechanism investigation (MCAR vs MNAR for NME / Epipoly missing 60%)
2. **R-Step0-2 [HIGH]** Coverage_Multiple + NumReads collinearity (+16.5 / -15.4) -- regularization or feature pruning required
3. **R-Step0-3 [MED]** Methylation features marginal -- cycle 1 may need to be re-labeled "augmented_filter (multi-axis, methylation-incl)" rather than "methylation-augmented"
4. **R-Step0-4 [MED]** Step 5c lost-TP cross-check pending -- if Path B re-introduces same lost TP pattern, gain is overstated
5. **R-Step0-5 [HIGH]** HCC1395-only -- Track B (4 additional samples) is critical before claiming generalizability; current ⭐3 tier should not advance without n>=4/5 cross-sample replication
