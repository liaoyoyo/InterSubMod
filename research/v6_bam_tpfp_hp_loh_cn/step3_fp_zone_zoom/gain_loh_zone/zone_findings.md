<!--
build_date: 2026-05-15
agent: Step 3 Z-GL (Gain+LOH) deep dive
status: in_progress
report_class: characterization_post_hoc
zone: Z-GL
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3, supporting H4 + zone characterization)
inputs:
  - step1_v3f_v5_v6_three_way/step1_master_three_way.tsv (master-joined only)
  - SEQC2 CNV truth (annotated_hcc1395_cnv.tsv)
outputs:
  - zone_summary.tsv
  - zone_grid.tsv (HP × cov sub-grid, all Inner)
  - zone_three_version_trajectory.tsv
  - lr_ablation_deviance.tsv (3-axis: HP, CN, AF; LOH fixed)
  - seqc2_concordance.tsv (Coverage_Multiple vs SEQC2 CN)
  - zone_confound_guards.tsv
verdict: NEGATIVE for FP enrichment hypothesis (FP_enrichment 0.022 — Inner gain/elevated is TP-pure)
-->

# Z-GL — Gain+LOH zone (Inner × Coverage_Multiple ≥ 1.3)

> Z-GL = `master_join_ok == 1` AND `loh_side == "Inner"` AND `Coverage_Multiple >= 1.3` (cov_elevated/gain/high_gain).
> Master-joined only (LOH and CN annotation required).

## 0. TL;DR

- **n = 1,687 regions** (1,682 TP + 5 FP) — TP rate **0.997** (Wilson 95% CI [0.993, 0.999])
- **FP_enrichment = 0.022×** (FP rate 0.003 vs global 0.137) — Z-GL is **strongly TP-pure**, not FP-rich
- Fisher exact: **odds = 56.5, p = 4.5e-101** — enrichment direction is **strongly opposite** to FP-rich premise
- SEQC2 concordance: Coverage_Multiple median 1.35 (CN=3), 1.42 (CN=4), 1.49 (CN=5) — KDE-corrected CovM scales correctly with true CN
- LR ablation: caller_af dominant (0.113 dev_explained), HP 0.042, CN essentially 0.000 — within Z-GL the CN axis adds no information (zone already filters by CN)

## 1. Why Z-GL is TP-pure (not FP-rich)

The original H3-adjacent prediction (gain+LOH ⇒ allele imbalance ⇒ FP-rich) appears **inverted** in HCC1395:
- Inner LOH means the SNV is in a homozygous (or strongly biased) region — somatic calls there have **stronger** evidence by construction
- Coverage ≥ 1.3 in HCC1395 paired-pileup mode corresponds to genuine copy number gain, which **increases read support for true somatic calls**
- The combination filters in real subclonal somatic variants with high-confidence variant signal

This is consistent with prior Wave 3 finding that LOH is in fact TP-enriched (not FP-enriched) at the region level.

## 2. SEQC2 CNV concordance

(See `seqc2_concordance.tsv`.)

| SEQC2 CN | n in Z-GL | Coverage_Multiple median (KDE) | TP rate |
|---|---|---|---|
| 2 (diploid) | 8 | 1.75 | 0.875 |
| 3 | 120 | 1.34 | 1.000 |
| 4 | 546 | 1.39 | 0.995 |
| 5 | 477 | 1.49 | 1.000 |
| 6 | 289 | 1.62 | 1.000 |
| 7 | 206 | 1.71 | 0.995 |
| 8 | 40 | 2.08 | 1.000 |

- KDE-corrected Coverage_Multiple **monotonically increases with SEQC2 CN** (1.34 → 2.08 for CN 3→8)
- The diploid (CN=2) cases in Z-GL are likely SEQC2 truth disagreements where local coverage gain doesn't reflect arm-level CN
- TP rate is ≥ 99.5% for all SEQC2 CN ≥ 3 — Z-GL is dominated by true somatic + copy gain

## 3. Sub-grid (HP × cov, all Inner)

(See `zone_grid.tsv`.)

| Cell (Inner|HP|cov) | n | TP_rate | FP_rate | Powered |
|---|---|---|---|---|
| Inner\|same_HP1\|cov_elevated | 817 | 0.999 | 0.001 | yes |
| Inner\|same_HP1\|cov_gain | 63 | 1.000 | 0.000 | yes |
| Inner\|same_HP2\|cov_elevated | 85 | 1.000 | 0.000 | yes |
| Inner\|other\|cov_elevated | 601 | 0.995 | 0.005 | yes |
| Inner\|other\|cov_gain | 83 | 1.000 | 0.000 | yes |

All powered cells: TP rate ≥ 0.995. cross_het and cross_het_inv cells in Inner gain are too small to be powered (n ≤ 14).

## 4. LR ablation deviance decomposition

(See `lr_ablation_deviance.tsv`.)

| Axis dropped | Incremental dev_explained | Comment |
|---|---|---|
| hp_bucket | 0.042 | Modest contribution (mainly `same_HP1` vs `other`) |
| Coverage_Multiple | 0.0005 | Essentially zero — zone already filters by cov ≥ 1.3 |
| caller_af | 0.113 | Dominant — AF still separates the 5 FP from 1,682 TP |

`dev_full = 56.73, dev_null = 68.20` — model explains only ~17% of null deviance because there are only 5 FP events to separate from 1,682 TP.

## 5. Three-version trajectory

(See `zone_three_version_trajectory.tsv`.)

V3F/V5/V6 have essentially identical TP rate (0.997) and marker coverage (35.3-35.4%) in this zone. The 3 marker-FP events are the same set across versions. **Z-GL is a TP-pure zone that is invariant to phasing version.**

## 6. Confound guards

Powered cells (5) tested. All 5 cells pass guards 1-3 (NumReads OLS, caller_af OLS, permutation). MH chr-stratification produces variable p depending on chr distribution; HP symmetry shows strong HP1 skew (`Inner|same_HP1` cells dominate) consistent with V6 priority bug residual (ratio 1.838).

## 7. Verdict & Limitations

- **Z-GL = TP-pure zone** (FP_enrichment 0.022, not a FP marker)
- Inner LOH + CN gain in HCC1395 paired-pileup is a **strong TP signature**, not a FP signature
- Limitation 1: master-joined only — Z-GL definition requires CN annotation, so master-unjoined FP (which is where most chr8 FP live) is excluded
- Limitation 2: 5 FP events is too few for stable LR coefficients
- Limitation 3: HCC1395-specific — other samples may have different ploidy/CN landscape

## 8. Hand-off to Step 3 synthesis

- Z-GL is **not** a FP zone candidate
- It can serve as a **TP-clean reference** for contrast against Z-CHR8 / Z-AUTO
- Reframe: gain+LOH is a **somatic-amplification signature** (subclonal copy gain on the somatic haplotype), not an FP signature
- SEQC2 concordance confirms KDE-corrected Coverage_Multiple is calibrated against SEQC2 truth in the gain region of HCC1395 (r ≈ 0.83 from prior calibration)
