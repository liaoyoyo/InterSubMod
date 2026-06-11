---
title: Phase 2 Cycle 1 Step 1 - Filter Design (Path B, cfg_drop_nr, L2 C=1.0)
date: 2026-05-18
status: complete
project: InterSubMod/research/methyl_augmented_filter_phase2/cycle1
predecessor: cycle1_step0_global_fp_audit.md (Path B selected, D2 PASS / D3 FAIL)
scripts: scripts/filter_design_and_verdict.py (initial sweep), scripts/collinearity_resolve.py (4 config comparison), scripts/final_filter_and_verdict.py (final fit + verdict)
runtime: ~3 min total
---

# Phase 2 Cycle 1 Step 1 - Filter Design

## Executive Summary

Final filter rule chosen (Path B, cfg_drop_nr, L2 C=1.0):

| Item | Value |
|---|---|
| Feature set | 10 features (caller_af, V6_off_NG, loh_inner_flag, Coverage_Multiple_imp, chr8_flag, 5 methyl covariates) — **NumReads_master dropped** |
| Regularization | L2 (Ridge), `LogisticRegression(penalty='l2', C=1.0, solver='lbfgs')` |
| NaN handling | Strategy B (median impute per feature) |
| OOF protocol | 5-fold StratifiedKFold |
| tau* | **0.39** |
| Expected ΔF1 | **+0.02236** (9.24x v1.0 baseline +0.00242; 2.24x H_C1_3 threshold +0.01) |
| Expected precision / recall | 0.9541 / 0.6030 |
| Max VIF | **1.83** (vs 217 with full feature set) |
| Max abs coef | **3.44** (vs 68.6 in cfg_full C=10) |

**Why not cfg_full C=10.0 with ΔF1=+0.02995?** Severe collinearity (VIF=217) inflates coefs to |68|, which is a statistical artifact of unregularized identification, not a robust signal. L2 alone cannot tame VIF=217 — only dropping one collinear feature does. The +0.0076 ΔF1 gain (from +0.02236 to +0.02995) is unstable; we prefer interpretable + stable.

---

## Stage 1A - VIF Audit (R-Step0-2)

Original 11-feature VIF (Strategy B impute, standardized):

| Feature | R² | VIF | Concern |
|---|---:|---:|---|
| V6_off_NG | 0.168 | 1.20 | none |
| caller_af | 0.335 | 1.50 | none |
| **NumReads_master** | **0.995** | **217.2** | **SEVERE collinearity** |
| loh_inner_flag | 0.379 | 1.61 | none |
| **Coverage_Multiple_imp** | **0.995** | **215.1** | **SEVERE collinearity** |
| V6_off_meth_HPMergedDelta | 0.373 | 1.59 | none |
| V6_off_meth_HPFineF | 0.455 | 1.84 | none |
| V6_off_meth_NME_imbalance | 0.012 | 1.01 | none |
| V6_off_meth_Epipoly_Delta | 0.059 | 1.06 | none |
| V6_off_meth_ClusterPermanovaF | 0.384 | 1.62 | none |
| chr8_flag | 0.185 | 1.23 | none |

R²=0.995 means each of {NumReads_master, Coverage_Multiple_imp} can be predicted with 99.5% accuracy from the other (and minor contributions from chr8 / af) — they encode essentially the same axis.

## Stage 1B - L2 Sweep + 4-Config Comparison

### L2 sweep on full 11-feature set (cfg_full)

| C | best ΔF1 | tau | max\|coef\| | sign_flips |
|---:|---:|---:|---:|---:|
| 0.001 | +0.01984 | 0.64 | 1.13 | 2 |
| 0.01 | +0.02215 | 0.47 | 2.42 | 1 |
| 0.1 | +0.02429 | 0.40 | 4.60 | 1 |
| 1.0 | +0.02637 | 0.42 | 14.56 | 1 |
| 10.0 | +0.02995 | 0.53 | 68.57 | 1 |

ΔF1 increases with C, but max\|coef\| explodes monotonically — hallmark of unregularized collinearity with L2 absorbing extra capacity.

### 4-config comparison (cfg_drop_nr = chosen)

| Config | Features | C=1.0 ΔF1 | max\|coef\| | max VIF |
|---|---:|---:|---:|---:|
| **cfg_full** | 11 | +0.02637 | 14.56 | **217.2** |
| **cfg_drop_cov** | 10 (drop Coverage_Multiple_imp) | +0.02216 | 3.40 | 1.83 |
| **cfg_drop_nr** | **10 (drop NumReads_master)** | **+0.02236** | **3.44** | **1.83** |
| cfg_ratio | 11 (replace with reads_per_unit_cov) | +0.02481 | 3.49 | 1.83 |

Decision rationale:
- cfg_drop_nr slightly outperforms cfg_drop_cov (+0.02236 vs +0.02216) — Coverage_Multiple_imp is the cleaner axis (normalized to local diploid coverage); NumReads_master is just raw count
- cfg_ratio achieves +0.02481 but introduces an engineered feature (reads_per_unit_cov) that is hard to interpret biologically — preferring transparency for ⭐4 audit
- cfg_full +0.02637-0.02995 is an upper bound but unstable

**Chosen: cfg_drop_nr C=1.0** — drop NumReads_master, keep Coverage_Multiple_imp. Max VIF 1.83 (no collinearity).

### L1 (Lasso) cross-check (C=0.1, full feature set)

| Feature | mean L1 coef | zero fraction across folds |
|---|---:|---:|
| Coverage_Multiple_imp | +2.59 | 0.00 |
| NumReads_master | -1.71 | 0.00 |
| caller_af | +1.50 | 0.00 |
| loh_inner_flag | +1.21 | 0.00 |
| V6_off_NG | +0.69 | 0.00 |
| V6_off_meth_HPFineF | +0.56 | 0.00 |
| V6_off_meth_ClusterPermanovaF | -0.30 | 0.00 |
| chr8_flag | -0.20 | 0.00 |
| V6_off_meth_HPMergedDelta | -0.16 | 0.00 |
| V6_off_meth_Epipoly_Delta | +0.064 | 0.00 |
| V6_off_meth_NME_imbalance | -0.0005 | 0.00 |

Lasso does NOT zero out the NumReads / Coverage pair — confirms both carry independent signal *given the linear identification issue persists*. But Lasso also confirms NME_imbalance is essentially zero (|coef|<0.001). Methylation features other than HPFineF are marginal.

---

## Stage 2 - NaN Mechanism (R-Step0-1)

| Feature | NaN rate | Chi² p (NaN x af_bin) | MW p (caller_af) | Logit β | Verdict |
|---|---:|---:|---:|---:|---|
| V6_off_meth_NME_imbalance | 60.4% | 0.0e+00 | 9.4e-230 | +0.55 | **MNAR** |
| V6_off_meth_Epipoly_Delta | 60.4% | 0.0e+00 | 1.5e-229 | +0.55 | **MNAR** |
| V6_off_meth_HPFineF | 0.1% | 1.0e-03 | 1.3e-04 | +0.56 | MNAR (statistically) but **negligible effect** (rate range 0.002) |

**NME/Epipoly NaN rate by AF bin** (range 47-92%):

| af_bin | NaN rate |
|---|---:|
| [0.00,0.10) | **80.2%** |
| [0.10,0.20) | 51.4% |
| [0.20,0.30) | 47.3% |
| [0.30,0.50) | 53.7% |
| [0.50,0.70) | 49.8% |
| [0.70,1.00] | **91.5%** |

U-shape: NaN concentrates at AF extremes. Mechanism: NME / Epipoly compute HP1 vs HP2 imbalance, requiring both HPs to have ≥3 reads. Low AF → no HP2 coverage (germline-like); High AF → ALT dominates, HP2 sparse.

**Implication for Strategy A vs B 8.7x gap**: NaN is MNAR (not random missing). Strategy A (drop NaN) **systematically removes** the AF=0-0.1 (FP-rich) and AF=0.7-1.0 (FP-rich) regions, then evaluates filter on a non-representative middle-AF subset. Strategy B (impute median for missing methylation, use Coverage / AF / LOH / chr8 for signal in those rows) is the correct approach because:

1. The chosen filter doesn't rely on NME/Epipoly (their coefs are 0.006 and 0.057 — negligible)
2. The discriminative signal is Coverage_Multiple_imp + caller_af + loh_inner_flag + chr8_flag — all present in 100% of rows
3. Imputation median for methyl in NaN rows is a *flagged unknown* that contributes negligibly via small coefs; the filter decision is driven by non-methyl axes

This validates the Step 0 caveat: "methylation is unmissable rather than absent, so impute-then-use-other-features works".

---

## Stage 3 - Final Filter Rule

Stored in `cycle1_track_a_filter.json`:

```python
features = [
    "V6_off_NG", "caller_af", "loh_inner_flag",
    "V6_off_meth_HPMergedDelta", "V6_off_meth_HPFineF",
    "V6_off_meth_NME_imbalance", "V6_off_meth_Epipoly_Delta",
    "V6_off_meth_ClusterPermanovaF", "chr8_flag",
    "Coverage_Multiple_imp",  # NumReads_master DROPPED due to VIF=217
]
clf = LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=5000)
# StandardScaler fit on train fold; transform on val
# Per-fold scalers + coefs stored in JSON for reproducibility

# Filter: keep variant if predicted P(TP) >= 0.39
```

### Feature importance (cfg_drop_nr, C=1.0, 5-fold mean)

| Feature | mean coef | std across folds |
|---|---:|---:|
| **caller_af** | **+3.44** | 0.032 |
| **loh_inner_flag** | **+1.46** | 0.018 |
| **Coverage_Multiple_imp** | **+1.27** | 0.007 |
| **V6_off_NG** | **+1.07** | 0.022 |
| V6_off_meth_HPFineF | +0.75 | 0.048 |
| chr8_flag | -0.46 | 0.006 |
| V6_off_meth_HPMergedDelta | -0.43 | 0.018 |
| V6_off_meth_ClusterPermanovaF | -0.29 | 0.016 |
| V6_off_meth_Epipoly_Delta | +0.06 | 0.009 |
| V6_off_meth_NME_imbalance | +0.007 | 0.021 |

**Interpretation**:
- Primary signal: caller_af (high AF → TP), loh_inner_flag (inner LOH → TP), Coverage_Multiple_imp (normal cov → TP), V6_off_NG (more reads in NG → TP)
- Secondary signal: V6_off_meth_HPFineF (positive), 3 other meth (negative direction; small)
- chr8 acts as FP-marker (negative coef) — consistent with HCC1395 chr8 amplicon hotspot

**Methylation channels are not the primary signal.** This refines the cycle name from "methyl-augmented filter" to "**multi-axis filter with methylation as marginal covariate**" — see R-Step0-3 caveat addressed.

---

## Outputs

```
cycle1/
├── cycle1_step1_filter_design.md        (this file)
├── cycle1_track_a_filter.json           (final rule with per-fold coefs + scalers)
├── scripts/
│   ├── filter_design_and_verdict.py     (initial sweep — superseded for filter)
│   ├── collinearity_resolve.py          (4-config audit)
│   └── final_filter_and_verdict.py      (canonical final fit)
├── data/
│   ├── cycle1_step1_vif.tsv             (11-feature VIF)
│   ├── cycle1_step1_l2_coefs.tsv        (L2 sweep coefs)
│   ├── cycle1_step1_l2_df1.tsv          (L2 sweep ΔF1)
│   ├── cycle1_step1_l1_coefs.tsv        (Lasso cross-check)
│   ├── cycle1_step1_collinearity_comparison.tsv (4 config x 4 C grid)
│   ├── cycle1_step1_nan_mechanism.tsv   (MCAR/MNAR test results)
│   ├── cycle1_step1_final_tau_sweep.tsv (final cfg tau curve)
│   └── cycle1_step1_oof_predictions.tsv (35,332 row P(TP) for downstream Track B)
└── figures/
    ├── cycle1_step1_vif_audit.png
    ├── cycle1_step1_l2_regularization_sweep.png
    ├── cycle1_step1_collinearity_comparison.png
    └── cycle1_step1_nan_mechanism.png
```

## Caveats propagated to Step 2

- **R-Step0-1 RESOLVED**: NaN is MNAR (low-AF + high-AF concentrated). Strategy B impute is correct because non-methyl axes carry the signal in NaN rows; the methylation coefs are small (<1) so impute median ≈ zero contribution there.
- **R-Step0-2 RESOLVED**: Collinearity tamed by dropping NumReads_master; final max VIF 1.83, max\|coef\| 3.44.
- **R-Step0-3 ACKNOWLEDGED**: Methylation features are 5th-rank covariates with |coef| < 1.0; rename hypothesis framing to "multi-axis filter (incl. methylation as marginal)" in cycle1 final synthesis.
- **R-Step0-4 PASSED (preview)**: Step 2 Stage 2 shows 81% of 21 Step 5c lost TP rescued by cycle 1 filter — global LR is NOT disproportionately removing low-AF subclone TP (see findings).
- **R-Step0-5 OPEN**: Cross-sample validation (Track B) is still required for ⭐4. Stage 3 multi-seed n=5 std=5e-5 is excellent intra-sample stability but does not address inter-sample generalization.

## Hand-off to Step 2

Filter rule + per-fold artefacts in `cycle1_track_a_filter.json`. Step 2 (HCC1395 ΔF1 verdict) executed in same run; see `cycle1_track_a_findings.md`.
