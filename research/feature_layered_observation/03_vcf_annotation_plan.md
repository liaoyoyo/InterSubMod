---
title: Phase B1 VCF Annotation Plan
date: 2026-04-23
status: completed
phase: B1
owner: feature_layered_observation
related:
  - research/feature_layered_observation/scripts/01_vcf_annotate.py
  - research/feature_layered_observation/data/merged_with_vcf.tsv.gz
  - research/feature_layered_observation/logs/01_vcf_parse_stats.tsv
  - research/feature_layered_observation/logs/01_vcf_join_rate.tsv
---

# Phase B1 — VCF Annotation Plan

## Summary

Annotated the Thread B merged master (`merged_7samples_paired_full_plus_7to_full_20260423.tsv.gz`, 748,676 rows) with caller-level signals from the corresponding ClairS / ClairS-TO VCFs for 7 samples × {paired_full, to_pileup}. All 14 (sample, mode) cells joined at 100% coverage; no row was dropped; 32 `vcf_*` columns were appended.

Output: `research/feature_layered_observation/data/merged_with_vcf.tsv.gz` (748,676 rows × 60 cols).

## VCF path table (14 cells)

| Sample | Mode | Caller / schema | Path (relative to repo root) |
|---|---|---|---|
| HCC1395 | paired_full | ClairS v0.4.0 | `output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/longphase_s/` |
| HCC1395_DORADO | paired_full | ClairS v0.4.0 | `output/canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full/longphase_s/` |
| COLO829 | paired_full | ClairS v0.4.0 | `output/canonical/COLO829/paired_full/20260421_COLO829_paired_full_full/longphase_s/` |
| H2009 | paired_full | ClairS v0.4.0 | `output/canonical/H2009/paired_full/20260421_H2009_paired_full_full/longphase_s/` |
| H1437 | paired_full | ClairS v0.4.0 | `output/canonical/H1437/paired_full/20260421_H1437_paired_full_full/longphase_s/` |
| HCC1937 | paired_full | ClairS v0.4.0 | `output/canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full/longphase_s/` |
| HCC1954 | paired_full | ClairS v0.4.0 | `output/canonical/HCC1954/paired_full/20260421_HCC1954_paired_full_full/longphase_s/` |
| HCC1395 | to_pileup | ClairS-TO v0.3.0 | `output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step04_benchmark_longphase_to/` |
| HCC1395_DORADO | to_pileup | ClairS-TO v0.3.0 | `output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step04_benchmark_longphase_to/` |
| COLO829 | to_pileup | ClairS-TO v0.3.0 | `output/synthesis/research_rounds/20260423_colo829_to_pilot/step04_benchmark_longphase_to/` |
| H2009 | to_pileup | ClairS-TO v0.3.0 | `output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step04_benchmark_longphase_to/` |
| H1437 | to_pileup | ClairS-TO v0.3.0 | `output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step04_benchmark_longphase_to/` |
| HCC1937 | to_pileup | ClairS-TO v0.3.0 | `output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step04_benchmark_longphase_to/` |
| HCC1954 | to_pileup | ClairS-TO v0.3.0 | `output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step04_benchmark_longphase_to/` |

Each directory holds `filtered_snv_tp.vcf.gz` and `filtered_snv_fp.vcf.gz`. The note below on TO archive vs hp_fix_reanalysis is in section "Problems & resolutions".

## Extracted INFO / FORMAT tags

### ClairS paired_full (v0.4.0) — superset actually present

| Tag | Source | Type | Notes |
|---|---|---|---|
| QUAL | record | float | filtered to PASS in TP/FP VCFs |
| FILTER | record | string | "PASS" in TP; FP VCFs contain labels such as Germline, LowQual |
| H | INFO flag | 0/1 | single-haplotype variant |
| FAU, FCU, FGU, FTU | INFO int | count | forward-strand tumor allele counts |
| RAU, RCU, RGU, RTU | INFO int | count | reverse-strand tumor allele counts |
| GT, GQ, DP, AF | FORMAT | mixed | per-sample genotype quality + tumor VAF |
| NAF, NDP | FORMAT | float / int | normal BAM VAF and depth |
| AD_ref, AD_alt | FORMAT tuple split | int | from AD |
| NAD_ref, NAD_alt | FORMAT tuple split | int | from NAD |

### ClairS-TO to_pileup (v0.3.0) — superset actually present

| Tag | Source | Type | Notes |
|---|---|---|---|
| QUAL | record | float | |
| FILTER | record | string | includes PASS, NonSomatic, LowQual, LowAltBQ, LowAltMQ, ReadStartEnd, VariantCluster, NoAncestry, MultiHap, StrandBias, LowSeqEntropy, RefCall |
| H | INFO flag | 0/1 | |
| FAU..RTU | INFO int | count | same meaning as paired |
| SB | INFO float | p-value | Fisher strand bias (only TO has this in INFO) |
| PoN_1..PoN_4 | INFO flag | 0/1 | gnomAD, dbSNP, 1000G, CoLoRSdb PoN hits |
| Verdict_Germline, Verdict_Somatic, Verdict_SubclonalSomatic | INFO flag | 0/1 | ClairS-TO post-hoc classification |
| GT, GQ, DP, AF | FORMAT | mixed | (no NAF/NDP/NAD – TO has no normal BAM) |
| AD_ref, AD_alt | FORMAT tuple split | int | |

All tags stored with prefix `vcf_` in the merged dataframe. Flags become 0/1 integers; missing tags are NaN.

## Join-rate results (per sample × mode)

| sample | mode | master_rows | vcf_qual_hits | join_rate |
|---|---|---:|---:|---:|
| COLO829 | paired_full | 37,458 | 37,458 | 1.00 |
| COLO829 | to_pileup | 50,617 | 50,617 | 1.00 |
| H1437 | paired_full | 67,476 | 67,476 | 1.00 |
| H1437 | to_pileup | 58,915 | 58,915 | 1.00 |
| H2009 | paired_full | 132,995 | 132,995 | 1.00 |
| H2009 | to_pileup | 137,695 | 137,695 | 1.00 |
| HCC1395 | paired_full | 30,381 | 30,381 | 1.00 |
| HCC1395 | to_pileup | 40,115 | 40,115 | 1.00 |
| HCC1395_DORADO | paired_full | 30,129 | 30,129 | 1.00 |
| HCC1395_DORADO | to_pileup | 40,428 | 40,428 | 1.00 |
| HCC1937 | paired_full | 12,588 | 12,588 | 1.00 |
| HCC1937 | to_pileup | 24,655 | 24,655 | 1.00 |
| HCC1954 | paired_full | 17,938 | 17,938 | 1.00 |
| HCC1954 | to_pileup | 67,286 | 67,286 | 1.00 |

Total: 748,676 / 748,676 = **100.0%** join coverage. No (sample, mode) had missing VCF.

No duplicate (sample, mode, Chr, Pos) keys within any VCF — the TP/FP files partition the truth space cleanly.

## caller_AF vs ISM AF validation

User requested Pearson r > 0.95 between `vcf_AF` (caller) and `AF` (master).

**Observation (HCC1395 paired_full, tp_label=1, n=29,754):**
- `Pearson(vcf_AF, AF)` = **0.6721**
- `Pearson(AlleleDelta, vcf_AF)` = 0.6891
- `Pearson(AF, |vcf_AF − 0.5|)` = 0.5000

**Root cause (not a join bug):** the `AF` column in the merged master is **not** the caller VAF. Sample of first 15 TP rows (HCC1395 paired_full):

```
Chr    Pos        AF       vcf_AF  AlleleDelta
chr1   877772     0.2038   0.8235  0.2038
chr1   1004726    0.0000   0.5600  -0.0153
chr1   1049980    0.0351   0.8864  0.0351
chr1   1212740    0.2110   0.9778  0.2110
```

`AF` mirrors `|AlleleDelta|` (ISM read-level allele-imbalance feature derived from aggregated reads in the region), not the caller single-site VAF. Mean caller `vcf_AF` ≈ 0.50 across TP rows (as expected for a somatic call set with clonal + subclonal variants), but mean master `AF` ≈ 0.058 because it is a derived imbalance score, not a frequency.

**Conclusion:** 100% row-level join is correct. The low Pearson is expected given the semantic mismatch. Downstream analyses should treat `vcf_AF` as *the* caller VAF and `AF` / `AlleleDelta` as ISM aggregate features. Recommend renaming `AF` → `AlleleDelta_abs` in a future master rebuild to avoid confusion; not done in this Phase to keep row-stable compatibility.

## Problems & resolutions

1. **archive TO rounds (Mar 2026) vs hp_fix_reanalysis (Mar 30)** — `202603_early_pilots/*/step04_benchmark_longphase_to/` retained the raw ClairS-TO + LongPhase output; `20260330_to_hp_fix_reanalysis/*/step04_*` only has downstream TSVs (no re-generated VCFs). VCF path therefore uses the March-17/18 fastresume rounds, which pre-date HP integer tag fix but share the same caller output (HP fix affected ISM read parsing only, not VCF content).
2. **COLO829 to_pileup** — the March archive `20260317_colo829_to_pilot/step04_benchmark_longphase_to/` is empty; switched to the new `20260423_colo829_to_pilot` pilot, which has the full TP/FP VCF set.
3. **HCC1395 to_pileup** — user suggested either canonical HCC1395/paired_full (wrong mode) or a non-existent `20260421_ReadParser_GermlineHPOnly_HCC1395` directory. Fell back to `20260315_hcc1395_to_pilot/step04_benchmark_longphase_to/` which has the matching TP/FP VCF pair.
4. **HCC1395_DORADO paired_full FP = 240 rows** but master lists 30,129 rows for this cell; join still 100% because FP VCFs only hold non-PASS variants not retained in the master (master already filtered by ISM). 100% join implies every master row sits in the TP VCF.
5. **FP VCF interpretation** — `filtered_snv_fp.vcf.gz` contains records that truth-matching classified as FP; the master already contains `tp_label` derived from the same truth comparison, so `vcf_tp_label_from_vcf` serves as a redundancy/consistency check (not a new source of truth).

## Verification commands

```bash
# Row-stability check
python3 -c "import pandas as pd; d=pd.read_csv('research/feature_layered_observation/data/merged_with_vcf.tsv.gz', sep='\t', low_memory=False); print(len(d), d.columns.size)"
# Expected: 748676 60

# Join rate by cell
cat research/feature_layered_observation/logs/01_vcf_join_rate.tsv

# Parse stats
cat research/feature_layered_observation/logs/01_vcf_parse_stats.tsv
```

## Files produced

- `research/feature_layered_observation/scripts/01_vcf_annotate.py` (parser)
- `research/feature_layered_observation/data/merged_with_vcf.tsv.gz` (output, 748,676 × 60)
- `research/feature_layered_observation/logs/01_vcf_parse_stats.tsv` (per-cell parse telemetry)
- `research/feature_layered_observation/logs/01_vcf_join_rate.tsv` (per-cell join coverage)
- `research/feature_layered_observation/logs/01_vcf_annotate.log` (full run log)

## Next steps

- Phase C G3 (VCF caller signals): stratify `vcf_FILTER`, `vcf_H`, `vcf_SB`, `vcf_Verdict_*`, `vcf_QUAL` by LOH × AF × CN groups.
- Consider rebuilding master with renamed `AF` → `AlleleDelta_abs` in a future cycle.
