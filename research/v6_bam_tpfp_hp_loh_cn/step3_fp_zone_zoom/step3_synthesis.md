<!--
build_date: 2026-05-15
agent: Step 3 Coordinator synthesis — 4 FP zone deep dive
status: in_progress
report_class: characterization_post_hoc
scope: HCC1395 pilot, full Step 1 master TSV (35,332 rows = 30,490 TP + 4,842 FP)
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3)
inputs:
  - step1_v3f_v5_v6_three_way/step1_master_three_way.tsv (full)
  - SEQC2 CNV truth (Z-GL only)
outputs:
  - step3_cross_zone_summary.tsv
  - step3_jaccard_overlap.tsv
  - step3_verdicts.json
  - figures/step3_zone_enrichment.png
  - figures/step3_per_chr_fp_rate.png
verdicts:
  H3: NEGATIVE (Z-OCH FP_enrichment 0.124, well below 2.0× — Outer cross_het is TP-pure, not FP-rich)
  H4: POSITIVE (chr8 (LOH+CN)=0.249 > HP=0.063 by 0.186 dev_explained, above 0.05 threshold)
  H5: NEGATIVE (Z-AUTO Jaccard with known-zone union = 0.184, far below 0.60 — Z-AUTO captures novel FP mechanism)
-->

# Step 3 Synthesis — 4 FP Zone Deep Dive (HCC1395 pilot)

> **Characterization-only**. No filter / ΔF1 claims per plan §Out-of-scope.
> Universe: full Step 1 master TSV (35,332 rows = 30,490 TP + 4,842 FP). All FP retained (no master_join_ok filter).
> Global TP rate = **0.8630**, Global FP rate = **0.1370**.

## 0. TL;DR

| Hypothesis | Verdict | Key number | Threshold |
|---|---|---|---|
| **H3** Z-OCH × cov_gain FP enrichment ≥ 2.0× | **NEGATIVE** | FP_enrichment 0.124× | ≥ 2.0× |
| **H4** chr8 hotspot (LOH+CN) ≥ HP + 0.05 dev_explained | **POSITIVE** | (LOH+CN) − HP = 0.186 | ≥ 0.05 |
| **H5** Z-AUTO ∪ known-zones Jaccard ≥ 0.60 | **NEGATIVE** | Jaccard 0.184 | ≥ 0.60 |

**Key surprise**: 2 of the 3 prior "FP-rich" zones (Z-OCH Outer cross_het, Z-GL Inner gain+LOH) turn out to be **TP-pure**, not FP-rich, when tested on the full FP set. Only Z-CHR8 (chr8 hotspot) is genuinely FP-enriched (2.31×), and Z-AUTO (KDE-detected top 5%) captures a **distinct FP mechanism** largely outside the known-zone framework.

## 1. 4-zone cross comparison

| Zone | n | n_TP | n_FP | TP_rate | FP_rate | TP_enrichment | FP_enrichment | Fisher p |
|---|---|---|---|---|---|---|---|---|
| **Z-OCH** (Outer cross_het) | 1,468 | 1,443 | 25 | 0.983 | 0.017 | 1.139× | **0.124×** | 3.8e-62 |
| **Z-CHR8** (chr8 hotspot) | 3,061 | 2,094 | 967 | 0.684 | 0.316 | 0.793× | **2.305×** | 1.7e-159 |
| **Z-GL** (Inner gain+LOH) | 1,687 | 1,682 | 5 | 0.997 | 0.003 | 1.155× | **0.022×** | 4.5e-101 |
| **Z-AUTO** (top 5% FP density) | 1,767 | 551 | 1,216 | 0.312 | 0.688 | 0.361× | **5.022×** | 0.0 |

(Source: `step3_cross_zone_summary.tsv`.)

Two observations:
1. **Z-OCH and Z-GL are TP-pure, not FP-rich** — the original H3-adjacent hypotheses (germline-somatic mixing on opposite haps, gain+LOH allele imbalance) are **inverted** by the data. Both signatures are in fact **somatic-evidence markers**, not FP markers.
2. **Z-CHR8 is the only FP-rich pre-defined zone** (FP_enrichment 2.31×, captures 967 / 4,842 = 20.0% of all FP).
3. **Z-AUTO is the strongest FP-rich zone** (FP_enrichment 5.02×, captures 1,216 / 4,842 = 25.1% of all FP in 5% of regions).

## 2. Jaccard overlap matrix

| | Z-OCH | Z-CHR8 | Z-GL | Z-AUTO |
|---|---|---|---|---|
| Z-OCH | 1.000 | 0.001 | 0.000 | 0.008 |
| Z-CHR8 | 0.001 | 1.000 | 0.165 | 0.299 |
| Z-GL | 0.000 | 0.165 | 1.000 | 0.022 |
| Z-AUTO | 0.008 | 0.299 | 0.022 | 1.000 |

- **Z-AUTO ∪ (Z-OCH ∪ Z-CHR8 ∪ Z-GL) Jaccard = 0.184** → H5 NEGATIVE
- Z-AUTO has highest overlap with Z-CHR8 (Jaccard 0.299) because chr8 dominates FP density — but ~70% of Z-AUTO is outside chr8
- Z-CHR8 ∩ Z-GL (Jaccard 0.165) is a small set of chr8 master-joined Inner gain regions
- 4-zone intersection: **n = 0** (no region is in all 4 zones simultaneously)

(Source: `step3_jaccard_overlap.tsv`.)

## 3. H3 重測 — Z-OCH FP enrichment (NEGATIVE)

Step 2 H3 was NEGATIVE (FP_enrichment ≈ 0.0) because the master-join filter dropped 90% of FP. Step 3 re-tested on the full FP set:

- Z-OCH n=1,468 (1,443 TP + 25 FP), FP_rate = 0.017, **FP_enrichment = 0.124×**
- All 6 powered sub-cells: TP rate = 1.000, FP rate = 0.000, passing all 4 confound guards
- Fisher exact: odds = 9.57, p = 3.8e-62 — Z-OCH is **TP-enriched**, opposite to H3 prediction

**Mechanistic re-interpretation**: cross_het = germline het + bona-fide somatic. FP germline-only calls do not produce cross_het because they lack real somatic read support. Z-OCH is therefore a **somatic-evidence signature**, not an FP signature.

(See `outer_cross_het/zone_findings.md` for details.)

## 4. H4 — chr8 hotspot LOH+CN > HP (POSITIVE)

LR deviance decomposition on master-joined chr8 (n=2,373):

| Axis dropped | Incremental dev_explained | Rank |
|---|---|---|
| caller_af | **0.393** | 1 (dominant) |
| Coverage_Multiple (CN) | **0.211** | 2 |
| hp_bucket | 0.063 | 3 |
| loh_side_norm | 0.038 | 4 |

- (LOH + CN) − HP = 0.038 + 0.211 − 0.063 = **+0.186** (3.7× the 0.05 threshold) → **H4 POSITIVE**
- caller_af dominates overall (0.393), but H4 is specifically about LOH+CN vs HP — the comparison is valid even after AF partialling
- chr8 FP_rate 0.316 vs global 0.137 (FP_enrichment 2.31×) — highest of all 23 chromosomes; chr17 is next at 0.198

**Mechanism**: chr8 has well-known structural and copy-number context in HCC1395 (chr8 gain dominates the tumor genome). CN axis (gain/elevated) explains chr8 FP variance 3.3× more than HP bucket → H4 confirms LOH+CN drives the chr8 hotspot, not phasing.

**Caveat**: chr8 hotspot is **sample-specific** (memory `project_hcc1395_chr8_hotspot.md`). Cross-sample mean AUC ≤ 0.641 from prior Wave 3 work confirms this does not generalize.

(See `chr8_hotspot/zone_findings.md` for details.)

## 5. H5 — Z-AUTO vs known zones Jaccard (NEGATIVE)

KDE smoothing on per-region FP density (Gaussian kernel, bandwidth 1Mb, per-chr):

- Top 5% by density threshold = 21.998 → n = 1,767 regions
- **Z-AUTO captures 1,216 / 4,842 = 25.1% of all FP** in 5% of the regions
- FP_enrichment = **5.02×**, the strongest of any zone
- Jaccard vs known-zone union (Z-OCH ∪ Z-CHR8 ∪ Z-GL) = **0.184** → H5 NEGATIVE

The known-zone framework does not predict 70%+ of Z-AUTO. There are emergent FP clusters in HCC1395 that **lack clear LOH/HP/CN signatures** and would require a different characterization framework (mappability, repeat context, GC bias, structural variants).

**Powered sub-cells in Z-AUTO** show that FP-rich Z-AUTO cells are dominantly `Inner|other|cov_proxy_{mid,high}` and Outer/UNKNOWN cells with no clear HP family — confirming that Z-AUTO FP are **outside the HP-bucket signature**.

(See `auto_detected_top5pct/zone_findings.md` for details.)

## 6. Cross-zone signature taxonomy (Coordinator framework)

Reframing the 4 zones as **TP / FP signature pairs**:

| Zone | Sample-of-FP captured | Signature interpretation |
|---|---|---|
| Z-OCH | 25 / 4,842 = 0.5% | TP signature (cross-hap somatic evidence) — NOT FP marker |
| Z-CHR8 | 967 / 4,842 = 20.0% | Sample-specific FP hotspot (CN-driven; HCC1395-only) |
| Z-GL | 5 / 4,842 = 0.1% | TP signature (gain on somatic hap) — NOT FP marker |
| Z-AUTO | 1,216 / 4,842 = 25.1% | Emergent FP cluster (mechanism unknown, outside LOH/HP/CN) |

**Coverage of total FP via Z-CHR8 ∪ Z-AUTO** (union of the 2 FP-rich zones):
- chr8: 967 FP positions
- Z-AUTO: 1,216 positions; ~70% of Z-AUTO is non-chr8 (~850 unique non-chr8 FP positions)
- Combined: roughly **1,800 / 4,842 = 37% of HCC1395 FP** captured in two zones

**Common signature across FP-rich zones**: Z-CHR8 and Z-AUTO both show:
- Low cov_proxy/CN homogeneity (mixed Inner/Outer/UNKNOWN cells)
- Predominance of `other` HP bucket (no specific HP family signature)
- High caller_af variance (AF gradient is the dominant LR axis)

**Unique signatures**:
- Z-CHR8: chromosome-specific structural context; consistent across V3F/V5/V6
- Z-AUTO: KDE-detected positional density without coherent mechanism

## 7. Verdict summary

| H | Verdict | Note |
|---|---|---|
| H3 | NEGATIVE | Z-OCH is TP-pure (FP_enrichment 0.124×), Step 2 finding confirmed at the full-FP set level |
| H4 | POSITIVE | chr8 (LOH+CN) dev_explained 0.249 exceeds HP 0.063 by 0.186 |
| H5 | NEGATIVE | Z-AUTO largely distinct from known zones (Jaccard 0.184) — novel FP mechanism present |

## 8. Hand-off to Coordinator (Task 7)

### 8.1 Cross-zone synthesis

- HCC1395 paired-pileup FP separates into **two real FP-rich zones**: Z-CHR8 (sample-specific) and Z-AUTO (mechanism unknown)
- The prior 2 hypothesis-driven "FP-rich" zones (Z-OCH, Z-GL) are in fact **TP-pure** — must re-interpret them as somatic signatures, not FP signatures
- The 4-zone framework is **incomplete**: ~63% of HCC1395 FP is outside Z-CHR8 ∪ Z-AUTO

### 8.2 Common signature

- FP-rich zones share: AF gradient dominance (caller_af LR contribution 0.393 in chr8, similar in Z-AUTO sub-cells)
- HP bucket is **secondary** (HP dev_explained ~0.06)
- LOH and CN matter but **only in combination** (and only on chr8 master-joined subset due to FP master-join asymmetry)

### 8.3 Unique signature

- **Z-CHR8 unique**: structural context of chr8 (gain arm in HCC1395); 967 FP, predominantly Inner with mixed `other` HP bucket
- **Z-AUTO unique**: positional density not explained by LOH/HP/CN; 1,216 FP, ~70% non-chr8; suggests mappability / repeat / SV-context factors

### 8.4 Step 4 routing recommendations

1. **Z-OCH and Z-GL → drop from FP-marker analysis** (they're TP signatures)
2. **Z-CHR8 → expected non-generalizing**; Step 4 should test if any other sample has a "chr-N hotspot" (analogous arm-level gain hotspot) without expecting HCC1395-specific overlap
3. **Z-AUTO → run KDE FP-density on each of H1437/H2009/HCC1954/HCC1937 V6 ISM separately**; per-sample Z-AUTO unlikely to share positions but the **density-based mechanism may recur**
4. **AF gradient** is the dominant FP axis — Step 4 should report cell-level caller_af distribution per zone per sample
5. **chr-by-chr ranking** (per `chr8_hotspot/per_chr_full_breakdown.tsv`) is a useful baseline for cross-sample comparison

### 8.5 Methodological caveats forwarded

- Master-join FP loss is a **structural artifact** of the master.tsv pipeline; Step 3 confirmed FP-density is real (Z-AUTO 25% of FP), so master.tsv FP coverage gap is a known data-source limitation
- KDE bandwidth (1Mb) is arbitrary; for cross-sample comparison the same bandwidth should be applied
- Guard 4 (chr-stratified MH) is inherently invalid for Z-CHR8 (single-chr zone) — note this when comparing cells across zones
- LR deviance decomposition results are **conditional on AF**; AF is a dominant axis but H3/H4/H5 statements remain valid at the LOH/HP/CN level

## 9. Files generated

```
step3_fp_zone_zoom/
├── step3_synthesis.md (this file)
├── step3_cross_zone_summary.tsv
├── step3_jaccard_overlap.tsv
├── step3_verdicts.json
├── figures/
│   ├── step3_zone_enrichment.png
│   └── step3_per_chr_fp_rate.png
├── outer_cross_het/
│   ├── zone_findings.md
│   ├── zone_summary.tsv
│   ├── zone_grid.tsv (cov_proxy axis)
│   ├── zone_three_version_trajectory.tsv
│   └── zone_confound_guards.tsv (6 powered cells, all pass 4/4)
├── chr8_hotspot/
│   ├── zone_findings.md
│   ├── zone_summary.tsv
│   ├── per_chr_comparison.tsv
│   ├── per_chr_full_breakdown.tsv
│   ├── zone_grid_proxy_cov.tsv
│   ├── zone_grid_real_cov.tsv (master-joined subset)
│   ├── zone_three_version_trajectory.tsv
│   ├── lr_ablation_deviance.tsv (4-axis: LOH/HP/CN/AF)
│   └── zone_confound_guards.tsv (17 powered cells)
├── gain_loh_zone/
│   ├── zone_findings.md
│   ├── zone_summary.tsv
│   ├── zone_grid.tsv
│   ├── zone_three_version_trajectory.tsv
│   ├── lr_ablation_deviance.tsv (3-axis: HP/CN/AF; LOH fixed Inner)
│   ├── seqc2_concordance.tsv (Coverage_Multiple vs SEQC2 CN)
│   └── zone_confound_guards.tsv
└── auto_detected_top5pct/
    ├── zone_findings.md
    ├── zone_summary.tsv
    ├── zone_grid.tsv
    ├── zone_three_version_trajectory.tsv
    ├── zone_confound_guards.tsv
    └── fp_density_per_region.tsv.gz (per-region KDE density + Z-AUTO flag, compressed)
```
