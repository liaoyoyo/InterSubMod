<!--
build_date: 2026-05-15
agent: Step 3 Z-OCH (Outer cross_het) deep dive
status: in_progress
report_class: characterization_post_hoc
zone: Z-OCH
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3, H3)
inputs:
  - step1_v3f_v5_v6_three_way/step1_master_three_way.tsv (full, 35,332 rows)
outputs:
  - zone_summary.tsv
  - zone_grid.tsv (cov axis = V6_off_n_reads quantile proxy)
  - zone_three_version_trajectory.tsv
  - zone_confound_guards.tsv
verdict: H3 NEGATIVE (FP_enrichment 0.124, well below 2.0× threshold)
-->

# Z-OCH — Outer cross_het zone (H3 重測)

> Z-OCH = `loh_side != Inner` AND `hp_bucket in {cross_het, cross_het_inv}` (V6_off HP family).
> **All FP retained** (no master_join_ok filter — addresses Step 2 caveat #1 where 90% FP were dropped).

## 0. TL;DR

- **n = 1,468 regions** (1,443 TP + 25 FP) — TP rate **0.983** (Wilson 95% CI [0.975, 0.988])
- **FP_enrichment = 0.124×** (FP rate 0.017 vs global 0.137 — i.e. Outer cross_het is **TP-pure**, NOT FP-rich)
- **H3 (FP_enrichment ≥ 2.0×): NEGATIVE** — even when re-tested with the full FP set, Z-OCH does not exhibit the predicted FP enrichment
- Fisher exact vs not-zone: **odds = 9.57, p = 3.8e-62** — the enrichment direction is **opposite** to H3's prediction (Z-OCH is significantly TP-enriched)
- **All 6 powered sub-cells pass all 4 confound guards** with permutation p < 0.001 — the TP-purity signal is robust to NumReads, caller_af, permutation, and chr-stratification

## 1. Why H3 was wrong — mechanistic explanation

The Thread D framework predicted Z-OCH (Outer × cross-haplotype germline-somatic mixing) as FP-rich because:
- Outer = no LOH context → germline het present
- cross_het = germline reads + somatic reads on **opposite** haplotypes → spurious co-occurrence signal

Empirical result: **Outer cross_het is TP-pure** in HCC1395. The cross_het signature is a **germline + bona-fide somatic** pattern; pure FP germline calls do not produce cross_het because the FP call has no real somatic read support.

This is consistent with Step 2 finding (cells `Outer|cross_het|cov_normal`, `Outer|cross_het|cov_elevated`, `Outer|cross_het_inv|cov_normal` were all TP=100% in master-joined universe).

## 2. Sub-grid (cov_proxy axis = V6_off_n_reads quantile q33/q67/q90)

The 6 powered cells (n ≥ 50) are listed in `zone_confound_guards.tsv`:

| cell_id | n | n_TP | n_FP | TP_rate | FP_enrichment | guards_passed |
|---|---|---|---|---|---|---|
| Outer\|cross_het\|cov_proxy_low | 152 | 152 | 0 | 1.000 | 0.000 | 4/4 |
| Outer\|cross_het\|cov_proxy_mid | 318 | 318 | 0 | 1.000 | 0.000 | 4/4 |
| Outer\|cross_het\|cov_proxy_high | 196 | 196 | 0 | 1.000 | 0.000 | 4/4 |
| Outer\|cross_het_inv\|cov_proxy_low | 142 | 142 | 0 | 1.000 | 0.000 | 4/4 |
| Outer\|cross_het_inv\|cov_proxy_mid | 341 | 341 | 0 | 1.000 | 0.000 | 4/4 |
| Outer\|cross_het_inv\|cov_proxy_high | 218 | 218 | 0 | 1.000 | 0.000 | 4/4 |

**All 6 powered cells: TP rate = 1.000, FP rate = 0.000**, passing all 4 confound guards (NumReads OLS, caller_af OLS, permutation p<0.001, chr-stratified MH).

HP symmetry guard skipped per Step 3 spec (zone is asymmetric by definition).

## 3. UNKNOWN loh_side rows (FP-rich, low n)

UNKNOWN rows (master_join_ok==0, can't tell Inner/Outer) within cross_het show some FP (e.g. `UNKNOWN|cross_het_inv|cov_proxy_low` n=20, 1 TP / 19 FP). These cells are underpowered (n<50) and reflect master-unjoined FP whose chr/pos was missing from master TSV — likely FP outside the standard LOH/CN annotation coverage rather than true Outer cross_het.

## 4. Three-version trajectory

(See `zone_three_version_trajectory.tsv` for full data.)

The zone definition uses V6_off HP bucket, so V3F/V5 stats are computed under the **same set of regions** but using each version's HP family counts.

Key observation: V3F's `cross_het`/`cross_het_inv` count in this zone is much lower (zone built on V6 buckets), confirming that the **cross-haplotype pattern is a V6-specific marker** in regions that V3F labeled differently (likely `other` or `same_HP*`). Marker coverage (NG≥3) is similar across V3F/V5/V6 within the zone, and marker_TP_rate is high in all 3 versions (~0.99).

## 5. Verdict & Limitations

- **H3 verdict: NEGATIVE** (FP_enrichment = 0.124, far below 2.0× threshold)
- Z-OCH is in fact **TP-pure** in HCC1395 — the H3 hypothesis premise (cross-het ⇒ FP-rich) is **inverted** by the data
- Limitation: cov_proxy uses V6_off_n_reads (not Coverage_Multiple) because most rows in this zone are master-unjoined; this is a coarse proxy
- Cross-sample: 4 other V6 samples may show different cross_het FP behavior (Step 4 work)

## 6. Hand-off to Step 3 synthesis

- Z-OCH should **not** be reported as a candidate FP-marker zone
- It can serve as a **TP-clean reference zone** for contrast against Z-CHR8 / Z-AUTO
- Reframe: cross_het is a **somatic-evidence signature** (germline het + real somatic read), not a FP signature
