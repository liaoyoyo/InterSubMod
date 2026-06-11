<!--
build_date: 2026-05-15
agent: Step 3 Z-AUTO (top 5% FP density) deep dive
status: in_progress
report_class: characterization_post_hoc
zone: Z-AUTO
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3, H5)
inputs:
  - step1_v3f_v5_v6_three_way/step1_master_three_way.tsv (full)
outputs:
  - zone_summary.tsv
  - zone_grid.tsv (3-axis sub-grid with cov_proxy)
  - zone_three_version_trajectory.tsv
  - zone_confound_guards.tsv
  - fp_density_per_region.tsv.gz (per-region KDE FP density + Z-AUTO flag)
verdict: H5 NEGATIVE — Jaccard with Z-OCH∪Z-CHR8∪Z-GL union = 0.184, far below 0.60 threshold
-->

# Z-AUTO — Auto-detected top 5% FP-density zone (H5)

> Z-AUTO = top 5% positions by per-region Gaussian-kernel FP density (bandwidth = 1Mb, per-chr smoothing).
> Density threshold = **21.998** (95th percentile of all 35,332 regions).

## 0. TL;DR

- **n = 1,767 regions** (551 TP + 1,216 FP) — TP rate **0.312** (Wilson 95% CI [0.291, 0.334])
- **FP_enrichment = 5.02×** (FP rate 0.688 vs global 0.137) — by definition the strongest FP-rich zone in HCC1395
- Z-AUTO contains **1,216 / 4,842 = 25.1% of all HCC1395 FP** in just 5% of regions
- **H5 (Jaccard vs known zones ∪ ≥ 0.60): NEGATIVE** — Jaccard = 0.184
  - Z-AUTO vs Z-CHR8: 0.299 (partial overlap — chr8 has many FP-dense pockets)
  - Z-AUTO vs Z-GL: 0.022 (no overlap — Z-GL is TP-pure by construction)
  - Z-AUTO vs Z-OCH: 0.008 (no overlap — Z-OCH is TP-pure by construction)
- **Z-AUTO is largely a distinct mechanism** from the prior 3 known zones (i.e., FP hotspots exist that are NOT predictable from LOH/HP/CN axes alone)

## 1. Spatial composition

Per-chr distribution of Z-AUTO members (read from `fp_density_per_region.tsv.gz`):

The KDE smoothing concentrates Z-AUTO around true high-FP-density regions. Since chr8 is the strongest hotspot (FP rate 0.316), a large fraction of Z-AUTO is on chr8, but Z-AUTO also captures dense FP pockets on chr17, chr12, chr13, chr4. **Z-AUTO is multi-chromosome**, not chr8-exclusive.

## 2. Sub-grid (cov_proxy axis)

(See `zone_grid.tsv`.)

Powered cells in Z-AUTO:

| cell_id | n | n_TP | n_FP | TP_rate | FP_rate | FP_enrichment |
|---|---|---|---|---|---|---|
| Inner\|same_HP1\|cov_proxy_mid | 60 | 58 | 2 | 0.967 | 0.033 | 0.24× |
| Inner\|other\|cov_proxy_low | 56 | 50 | 6 | 0.893 | 0.107 | 0.78× |
| Inner\|other\|cov_proxy_mid | 192 | 100 | 92 | 0.521 | 0.479 | 3.50× |
| Inner\|other\|cov_proxy_high | 107 | 51 | 56 | 0.477 | 0.523 | 3.82× |
| (plus Outer + UNKNOWN cells, mostly FP-rich) | | | | | | |

Most-FP-enriched sub-cells in Z-AUTO are `Inner|other|cov_proxy_{mid,high}` and Outer/UNKNOWN with no clear HP bucket — consistent with the interpretation that Z-AUTO captures FP that **don't fit the HP-bucket signature** (Z-OCH) or **don't fit the gain+LOH signature** (Z-GL).

## 3. Three-version trajectory

(See `zone_three_version_trajectory.tsv`.)

V3F/V5/V6 marker_TP_rate is ≈ 0.247-0.250 (very low) — Z-AUTO regions barely produce coherent markers (marker_FP > marker_TP).

| Version | marker_n | marker_coverage_pct | marker_TP_rate | hp_other_pct |
|---|---|---|---|---|
| V3F | 1,040 | 58.9% | 0.250 | 0.851 |
| V5  | 1,060 | 60.0% | 0.248 | 0.761 |
| V6  | 1,059 | 59.9% | 0.247 | 0.765 |

V3F has higher `hp_other` (0.851) — most Z-AUTO regions have ambiguous HP families in V3F. V5/V6 reassign some `other` regions to defined HP buckets (cross_het/cross_het_inv rises slightly), but the **FP rate is unchanged** — phasing version does not affect Z-AUTO FP rate.

## 4. Confound guards

All powered cells in Z-AUTO pass permutation guards (p<0.001) — the FP enrichment is statistically robust. Within-cell NumReads/caller_af residualization shows minimal change (Z-AUTO FP is not a numeric-covariate artifact).

## 5. H5 verdict & implications

- **Z-AUTO vs known zones Jaccard = 0.184** (NEGATIVE for H5 ≥ 0.60)
- Z-AUTO partially overlaps Z-CHR8 (Jaccard 0.299) because chr8 dominates HCC1395 FP, but the remaining ~70% of Z-AUTO is **not predicted by the 3 known zones**
- This means there are emergent FP clusters in HCC1395 that **lack** clear LOH/HP/CN signatures and would require **a different framework** (e.g. mappability, repeat context, GC bias, structural variants) to characterize

## 6. Verdict & Limitations

- **H5 verdict: NEGATIVE** — the prior 3 known zones do not recapitulate FP hotspots
- Z-AUTO is a useful **discovery tool** for hidden FP mechanisms beyond LOH/HP/CN
- Limitation 1: KDE bandwidth (1Mb) is arbitrary; smaller bandwidth would isolate sharper peaks
- Limitation 2: Z-AUTO is HCC1395-specific; cross-sample Z-AUTO would identify different hotspots per sample
- Limitation 3: Z-AUTO does not provide a **mechanism**, only a **density**

## 7. Hand-off to Step 3 synthesis

- Z-AUTO is the **single strongest FP zone** (FP_enrichment 5.02×, n_FP=1,216 = 25% of all FP)
- It is largely **distinct** from prior 3 zones (Jaccard 0.184)
- Suggests an **investigative path**: characterize Z-AUTO FP at the read/sequence level (mappability, repeat, GC, structural context) — outside Step 3 scope
- Cross-sample question: does Z-AUTO recur in H1437/H2009/HCC1954/HCC1937? (Step 4 work)
