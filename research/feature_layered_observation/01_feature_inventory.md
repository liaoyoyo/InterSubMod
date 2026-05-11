# 01 · Feature Inventory (Phase A Output + Phase E Verdict Update)

**Date**: 2026-04-23
**Status**: Phase A registry + Phase E verdict overlay
**Registry source**: `research/feature_layered_observation/data/feature_registry.tsv` (137 features, 10 groups)
**Downstream**: Phase C (G1–G10 agents) consume via `scripts/00_feature_registry.py`
**Phase E verdict 彙整**：見 `00_main_observation.md` §9，以及每群組 `.md` 的 §10。關鍵更新：
- `NumReads` → **CONDITIONAL_POSITIVE** (resid 0.777)
- `HPFineN_HP1S` / `same_hap_marker` → **CONDITIONAL_POSITIVE** / POSITIVE_FINGERPRINT (G6 核心)
- `SampleASM_Delta` / `NormalBaseline_Mean` / `ClusterPermanovaF` / `AlleleDelta` / `LabelAllelePermanovaF` → **CONFOUND_COLLAPSED**
- `Combined_HP_Signed_Delta` / `Normal_HP_*` / `HP_Signed_Residual` → **DATA_GAP** (0% populated canonical)
- `Stability` → **NOT_IMPLEMENTED** (全 0，疑為 bug)
- `Coverage_Multiple` / `NumCpGs` / `Quality_Score` / `HeuristicScore` → **NEGATIVE** (resid <0.53)

---

## Purpose

Central catalog for every feature that enters the feature-layered-observation framework. Each entry lists:

- **C++ / VCF / BAM source** — file:line or pysam extractor (BAM QC deferred to Phase B2)
- **One-sentence computation / formula**
- **Upstream deps** — which features must be computed first
- **Filter condition** — when the value is valid
- **Data type & range** — continuous / categorical / binary / ordinal
- **Prior conclusion ID** — if a previous experiment already characterized this feature

Features with `INFERRED_FROM_NAME` tag require source verification in a later pass (see §Summary).

---

## Summary statistics

| Metric | Count | Note |
| --- | --- | --- |
| Total features | **137** | 117 ISM + 11 VCF (G3) + 9 BAM QC (G4) |
| Continuous | 72 | Main statistical targets |
| Categorical | 14 | Includes LOH_Subtype, Quality_Tier, VerificationClass, DominantLabel |
| Binary (true/false) | 23 | Significance / validity / overlap flags |
| Ordinal (counts) | 28 | Read / CpG / group counts |
| With `INFERRED_FROM_NAME` | **10** | 1 ISM (`Combined_HP_Signed_Delta`) + 9 BAM QC (G4, extractor TBD) |

### Per-group distribution

| Group | Name | n | Continuous | Categorical | Binary | Ordinal | Inferred |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| G1 | Read Counting & Coverage | 12 | 4 | 4 | 0 | 4 | 0 |
| G2 | LOH & Subclone Annotation | 7 | 1 | 4 | 2 | 0 | 0 |
| G3 | Caller VCF Signals | 11 | 4 | 1 | 4 | 0 | 0 |
| G4 | BAM Read-Level QC | 9 | 9 | 0 | 0 | 0 | 9 |
| G5 | HP Merged (2-bucket) | 7 | 4 | 0 | 1 | 2 | 0 |
| G6 | HP Fine (4-bucket LOH-phasing) | 17 | 10 | 0 | 1 | 5 | 0 |
| G7 | Cluster & Global Methylation | 11 | 7 | 1 | 2 | 0 | 0 |
| G8 | Entropy / Epipoly / PerCpgASM | 11 | 8 | 0 | 1 | 2 | 0 |
| G9 | ASM (Allele + Sample + Residual) | 25 | 14 | 0 | 5 | 6 | 1 |
| G10 | Quality Summary & Verification | 27 | 11 | 4 | 7 | 5 | 0 |
| **Total** | | **137** | **72** | **14** | **23** | **25** | **10** |

### Features flagged `INFERRED_FROM_NAME` (require Phase B2 source lookup)

All 9 G4 BAM QC features are INFERRED because their extractor script (`scripts/analysis/bam_readqc_per_region.py`) is not yet authored — Phase B2 will define and verify during BAM sweep. Plus 1 ISM:

- `Combined_HP_Signed_Delta` (G9) — writer at `RegionProcessor.cpp:1246`, but the aggregation weighting between tumor and normal signed delta is not obvious from the call site alone. The G9 Phase C agent should re-read `RegionProcessor.cpp` (search for `combined_hp_signed_delta =`) to confirm.

---

## Feature tables (by group)

> Each row below is a compact projection of the registry TSV. Full fields (upstream_deps, filter_condition, range, prior_conclusion) remain in the TSV — reproduce them in Phase C group `.md` files as needed. "Used in B/D" flags Thread B (LOH×AF×CN stratification) or Thread D (LOH-constrained phasing) prior use.

### G1 · Read Counting & Coverage (12)

| Feature | Source | Type | Computation (1-liner) | Used in B/D |
| --- | --- | --- | --- | --- |
| RegionID | RegionProcessor.cpp:1078 | ordinal | Region index | — |
| Chr | RegionProcessor.cpp:1157 | categorical | Chromosome name | B |
| Pos | RegionProcessor.cpp:1157 | ordinal | Genomic position | B |
| Ref | RegionProcessor.cpp:1157 | categorical | Reference allele | — |
| Alt | RegionProcessor.cpp:1157 | categorical | Alternate allele | — |
| NumReads | RegionProcessor.cpp:1157 / MatrixBuilder.hpp:80 | continuous | Read count in region | B confound |
| NumCpGs | MatrixBuilder.hpp:83 | continuous | CpG count | B |
| Coverage_Multiple | RegionProcessor.cpp:63 `compute_coverage_multiple` | continuous | NumReads / diploid_coverage | B (CovM axis) |
| Diploid_Coverage_Used | RegionProcessor.cpp:697 | continuous | Per-sample diploid cov | KDE fix pending |
| Coverage_Category | RegionProcessor.cpp | categorical | Binned Low/Normal/High/Amp | B |
| NTumorReads | RegionProcessor.cpp num_tumor_reads | ordinal | Tumor read count | — |
| NNormalReads | RegionProcessor.cpp num_normal_reads | ordinal | Normal read count | — |

### G2 · LOH & Subclone Annotation (7)

| Feature | Source | Type | Computation | Used in B/D |
| --- | --- | --- | --- | --- |
| HP_Ratio | RegionProcessor.cpp:39 `compute_hp_ratio` | continuous | HP1 / (HP1 + HP2) | B (LOH axis) |
| Potential_LOH | RegionProcessor.cpp | binary | HP_Ratio < 0.1 OR > 0.9 | B primary LOH |
| LOH_Subtype | RegionProcessor.cpp:188 `determine_loh_subtype` | categorical | None / Noise / Weak / Strong / Subclone | B 5-class |
| LOH_Bed_Overlap | LohBedAnnotator.cpp `classify_loh_source` | binary | SNV in external LOH BED | B (TO 16-52x mismatch) |
| LOH_Source | LohBedAnnotator.cpp `LohSource` enum | categorical | NONE / BED_ONLY / RATIO_ONLY / BOTH | B |
| LOH_Bed_Annotation | LohBedAnnotator.cpp | categorical | BED name column | — |
| Subclone_ID | SubcloneAnalyzer.cpp:100 `assign_subclones` | categorical | -1 / 0..K-1 | — |

### G3 · Caller VCF Signals (11, from ClairS-TO v0.3.0)

| Feature | Source | Type | Computation | Used in B/D |
| --- | --- | --- | --- | --- |
| VCF_QUAL | VCF col 6 | continuous | Phred-like variant quality | — |
| VCF_FILTER | VCF col 7 | categorical | PASS / LowQual / NonSomatic / LowAltBQ / LowAltMQ / ReadStartEnd / VariantCluster / NoAncestry / MultiHap / StrandBias / LowSeqEntropy / Realignment / RefCall | — |
| VCF_AF | FORMAT/AF | continuous | Tumor allele frequency | B (AF axis, AUC 0.654) |
| VCF_DP | FORMAT/DP | continuous | Depth at site (caller) | — |
| VCF_AD | FORMAT/AD | categorical tuple | [ref_count, alt_count] | — |
| VCF_SB | INFO/SB | continuous | Fisher exact p for strand bias | — |
| VCF_PoN | INFO/PoN_1..4 flags | categorical | gnomAD / dbSNP / 1000g / CoLoRSdb | — |
| VCF_Verdict_Germline | INFO/Verdict_Germline flag | binary | Caller Germline verdict | ClairS-TO Pilot NEG on F1 |
| VCF_Verdict_Somatic | INFO/Verdict_Somatic flag | binary | Caller Somatic verdict | — |
| VCF_Verdict_Subclonal | INFO/Verdict_SubclonalSomatic flag | binary | Caller Subclonal verdict | — |
| VCF_H_OneHap | INFO/H flag | binary | Variant only in one haplotype | — |

**Note**: G3 requires Phase B1 VCF-to-TSV join (script `scripts/01_vcf_annotate.py`, TBD). All 7 sample × 2 mode (paired_full + TO) filtered VCFs need to be parsed and joined on (Chr, Pos, Ref, Alt) to extend the master TSV.

### G4 · BAM Read-Level QC (9, extractor TBD — all INFERRED_FROM_NAME)

| Feature | Extractor (planned) | Type | Computation |
| --- | --- | --- | --- |
| BAM_MapQ_mean | pysam fetch + mean(read.mapping_quality) | continuous | Mean MAPQ |
| BAM_MapQ_median | pysam fetch + median | continuous | Median MAPQ |
| BAM_MapQ_p10 | pysam fetch + np.percentile(10) | continuous | 10th pct MAPQ |
| BAM_MapQ_p90 | pysam fetch + np.percentile(90) | continuous | 90th pct MAPQ |
| BAM_NM_mean | read.get_tag('NM') | continuous | Mean edit distance |
| BAM_Softclip_Frac | CIGAR S ops / total reads | continuous | Softclip read fraction |
| BAM_Strand_Bias | read.is_reverse counts | continuous | Fisher p OR |fwd-rev|/(fwd+rev) (TBD) |
| BAM_Read_Length_mean | read.query_length | continuous | Mean read length |
| BAM_LowMQ20_Frac | mapping_quality < 20 | continuous | Fraction MAPQ < 20 |

**Note**: Phase B2 background sweep pending (~6-8 hr). G4 agent blocked until extractor and results available.

### G5 · HP Merged (2-bucket) (7)

| Feature | Source | Type | Computation |
| --- | --- | --- | --- |
| GlobalP_HPFamily | GlobalTest.cpp:252 `test_hp_family` | continuous | Fisher FFH p on cluster × HP-family |
| CramersV_HPFamily | GlobalTest.cpp | continuous | Cramér's V |
| HPMergedDelta | RegionProcessor.cpp:1621 `hp_ms.merged.delta` | continuous | |mean_HP1fam − mean_HP2fam| |
| HPMergedP | LabelTest.cpp | continuous | Permutation p |
| HPMergedSig | LabelTest.cpp | binary | Significance flag |
| HP1FamilyN | RegionProcessor.cpp hp_merged_n_hp1 | ordinal | Count HP1 ∪ HP1-1 |
| HP2FamilyN | RegionProcessor.cpp hp_merged_n_hp2 | ordinal | Count HP2 ∪ HP2-1 |

### G6 · HP Fine (4-bucket LOH-phasing) (17, Thread D core)

| Feature | Source | Type | Computation |
| --- | --- | --- | --- |
| GlobalP_HPFine | GlobalTest.cpp `test_hp_fine` | continuous | Fisher FFH p on cluster × HP_fine (4-bucket) |
| CramersV_HPFine | GlobalTest.cpp | continuous | Cramér's V |
| HPFine_NGroups_CF | RegionProcessor.cpp:1604 | ordinal | Cluster-First n_valid_groups |
| HPFineF | LabelTest.cpp | continuous | PERMANOVA pseudo-F on HP_fine labels |
| HPFineP | LabelTest.cpp | continuous | Permutation p |
| HPFineSig | LabelTest.cpp | binary | Significance flag |
| HPFineNGroups | LabelTest.cpp:265 `hp_to_fine_labels` | ordinal | N non-empty among {HP1, HP1-1, HP2, HP2-1}, HP0/HP3 excluded |
| HPFineN_HP1 | RegionProcessor.cpp:1186 | ordinal | Count HP1 |
| HPFineN_HP1S | RegionProcessor.cpp:1186 | ordinal | Count HP1-1 |
| HPFineN_HP2 | RegionProcessor.cpp:1187 | ordinal | Count HP2 |
| HPFineN_HP2S | RegionProcessor.cpp:1187 | ordinal | Count HP2-1 |
| HPFineD_HP1_HP1S | RegionProcessor.cpp:1189 | continuous | Mean distance HP1 vs HP1-1 |
| HPFineD_HP1_HP2 | RegionProcessor.cpp:1189 | continuous | Mean distance HP1 vs HP2 |
| HPFineD_HP1_HP2S | RegionProcessor.cpp:1190 | continuous | Mean distance HP1 vs HP2-1 |
| HPFineD_HP1S_HP2 | RegionProcessor.cpp:1190 | continuous | Mean distance HP1-1 vs HP2 |
| HPFineD_HP1S_HP2S | RegionProcessor.cpp:1191 | continuous | Mean distance HP1-1 vs HP2-1 |
| HPFineD_HP2_HP2S | RegionProcessor.cpp:1191 | continuous | Mean distance HP2 vs HP2-1 |

### G7 · Cluster & Global Methylation Structure (11)

| Feature | Source | Type | Computation |
| --- | --- | --- | --- |
| GlobalP | GlobalTest.cpp `test_clusters` | continuous | Fisher FFH p on cluster × methylation |
| CramersV | MathUtils.hpp:98 | continuous | Effect size; sparse in TO |
| PairwiseMeanDist | DistanceMatrix.cpp | continuous | Mean pairwise distance |
| PairwiseMedianDist | DistanceMatrix.cpp | continuous | Median pairwise distance |
| ClusterPermanovaF | SignificanceAnalyzer.cpp | continuous | PERMANOVA pseudo-F (cluster) |
| ClusterPermanovaP | SignificanceAnalyzer.cpp | continuous | Permutation p |
| ClusterPermanovaValid | SignificanceAnalyzer.cpp | binary | Group-size validity |
| ClusterDispersionP | SignificanceAnalyzer.cpp | continuous | PERMDISP p |
| ClusterDispersionWarn | SignificanceAnalyzer.cpp | binary | Dispersion warning |
| LocalBestCluster | LocalTest.cpp | categorical | Best local-match cluster ID |
| LocalBestP | LocalTest.cpp | continuous | Best local p |

### G8 · Entropy / Epipolymorphism / PerCpgASM (11)

| Feature | Source | Type | Computation |
| --- | --- | --- | --- |
| PerCpgASM_Valid | PerCpgAsm.cpp | binary | Coverage gate |
| Fisher_N_Sig | PerCpgAsm.cpp:336 | ordinal | Count CpGs with BH-FDR sig |
| Fisher_Frac_Sig | PerCpgAsm.cpp:336 | continuous | Fraction significant (AUC 0.726 paired) |
| Fisher_N_Tested | PerCpgAsm.cpp | ordinal | Tested CpG count |
| Fisher_MaxNegLogFDR | PerCpgAsm.cpp:337 | continuous | max(-log10 BH-FDR) |
| NME_HP1 | PerCpgAsm.cpp:356 `compute_nme(hp1,2)` | continuous | Normalized methylation entropy HP1 |
| NME_HP2 | PerCpgAsm.cpp:357 `compute_nme(hp2,2)` | continuous | NME HP2 |
| Entropy_Imbalance | PerCpgAsm.cpp:360 | continuous | |NME_HP1 − NME_HP2| |
| Epipoly_HP1 | PerCpgAsm.cpp:367 `compute_epipolymorphism(hp1,4)` | continuous | Scherer epipoly HP1 |
| Epipoly_HP2 | PerCpgAsm.cpp:368 | continuous | Epipoly HP2 |
| Epipoly_Delta | PerCpgAsm.cpp:371 | continuous | |Epipoly_HP1 − Epipoly_HP2| |

### G9 · Allele-Specific Methylation (25)

| Feature | Source | Type | Computation |
| --- | --- | --- | --- |
| AlleleDelta | RegionProcessor.cpp:1648 label_allele.delta | continuous | |mean_ALT − mean_REF| — **known L2 collider w/ AF** |
| AlleleP | RegionProcessor.cpp:1649 | continuous | Permutation p |
| AlleleSig | RegionProcessor.cpp:1650 | binary | Significance flag |
| SampleASM_Delta | SampleAsm.cpp | continuous | d_between(T,N) − d_within — **R1-Global NEGATIVE residualized AUC 0.527** |
| SampleASM_P | SampleAsm.cpp | continuous | Permutation p |
| SampleASM_Sig | SampleAsm.cpp | binary | Significance |
| SampleASM_NTumor | SampleAsm.cpp | ordinal | Tumor n used |
| SampleASM_NNormal | SampleAsm.cpp | ordinal | Normal n used |
| NormalBaseline_Mean | NormalBaseline.cpp / RegionProcessor.cpp:894 | continuous | Mean methylation (all normal reads) |
| NormalBaseline_Coverage | NormalBaseline.cpp | continuous | Normal read coverage |
| HP_Residual_Delta | RegionProcessor.cpp | continuous | Tumor_HP_Delta − Normal_HP_Delta |
| HP_Residual_P | RegionProcessor.cpp | continuous | Permutation p |
| HP_Residual_Sig | RegionProcessor.cpp | binary | Significance |
| Tumor_HP_Delta | RegionProcessor.cpp | continuous | Unsigned tumor HP1-HP2 delta |
| Tumor_HP_Valid | RegionProcessor.cpp | binary | Min-count gate |
| Normal_HP_Delta | NormalBaseline.cpp | continuous | Unsigned normal HP1-HP2 delta |
| Normal_HP_Valid | NormalBaseline.cpp | binary | Min-count gate |
| Tumor_HP1 | RegionProcessor.cpp:943 | ordinal | Tumor HP1 count |
| Tumor_HP2 | RegionProcessor.cpp | ordinal | Tumor HP2 count |
| Normal_HP1 | NormalBaseline.cpp | ordinal | Normal HP1 count |
| Normal_HP2 | NormalBaseline.cpp | ordinal | Normal HP2 count |
| Tumor_HP_Signed_Delta | RegionProcessor.cpp:1009 `compute_signed_hp_delta` | continuous | Signed tumor delta (HP1 − HP2) |
| Normal_HP_Signed_Delta | RegionProcessor.cpp | continuous | Signed normal delta |
| HP_Signed_Residual | RegionProcessor.cpp:1021 | continuous | Tumor − Normal signed |
| Combined_HP_Signed_Delta | RegionProcessor.cpp | continuous | Aggregate signed delta — **INFERRED_FROM_NAME** |

### G10 · Quality Summary & Verification (27)

| Feature | Source | Type | Computation |
| --- | --- | --- | --- |
| HeuristicScore | SignificanceAnalyzer.cpp:378 | continuous | Composite score [0-1] |
| PassedGating | RegionProcessor.cpp:1608 | binary | global_p ≤ 0.1 |
| LabelHPPermanovaF | SignificanceAnalyzer.cpp | continuous | PERMANOVA F on HP-family |
| LabelHPPermanovaP | SignificanceAnalyzer.cpp | continuous | Permutation p |
| LabelHPPermanovaValid | SignificanceAnalyzer.cpp | binary | Validity |
| LabelHPDispersionP | SignificanceAnalyzer.cpp | continuous | PERMDISP p |
| LabelHPDispersionWarn | SignificanceAnalyzer.cpp | binary | Dispersion warning |
| LabelAllelePermanovaF | SignificanceAnalyzer.cpp | continuous | PERMANOVA F on allele — **AF proxy?** |
| LabelAllelePermanovaP | SignificanceAnalyzer.cpp | continuous | Permutation p |
| LabelAllelePermanovaValid | SignificanceAnalyzer.cpp | binary | Validity |
| LabelAlleleDispersionP | SignificanceAnalyzer.cpp | continuous | PERMDISP p |
| LabelAlleleDispersionWarn | SignificanceAnalyzer.cpp | binary | Dispersion warning |
| UnassignedAffinity | LabelTest.cpp:68 `test_unassigned_affinity` | continuous | HP3/HP0 → HP1 vs HP2 affinity |
| UnassignedAffinityP | LabelTest.cpp | continuous | Permutation p |
| UnassignedDir | LabelTest.cpp | categorical | HP1 / HP2 / None |
| NHP3 | LabelTest.cpp | ordinal | Count HP3-tagged reads |
| NHP0 | LabelTest.cpp | ordinal | Count HP0-tagged reads |
| DominantLabel | SignificanceAnalyzer.cpp:281 | categorical | hp / allele / cluster / none |
| Stability | SignificanceAnalyzer.cpp cluster_stability | continuous | Bootstrap stability score |
| VerificationClass | SignificanceAnalyzer.cpp:330-339 | categorical | Strong / Subclone / Weak / Noise |
| Quality_Score | RegionProcessor.cpp:241 `compute_quality_score` | continuous | Composite [0-100] |
| Quality_Tier | RegionProcessor.cpp:296 `determine_quality_tier` | categorical | High / Medium / Low / VeryLow |
| Significant | RegionProcessor.cpp:1143 | binary | Gate ∧ p≤0.05 ∧ V≥0.1 ∧ reads≥20 |
| SuggestFilter | RegionProcessor.cpp:1149 | binary | label_delta > 0.3 |
| NHP_Somatic11 | ReadParser.cpp | ordinal | Raw HP tag 11 count (self-phasing audit) |
| NHP_Somatic21 | ReadParser.cpp | ordinal | Raw HP tag 21 count |
| NHP_Somatic33 | ReadParser.cpp | ordinal | Raw HP tag 33 count |

---

## Pitfalls & notes to propagate to Phase C agents

1. **`HPFine_NGroups_CF` vs `HPFineNGroups`**: both are computed from `hp_to_fine_labels`; the CF (Cluster-First) column writes `global_hp_fine_n_groups` from the Stage-1 cluster-first path, while `HPFineNGroups` comes from the Stage-2 LabelTest path. They are usually equal but can diverge when cluster-first validity gates fail. Treat as near-duplicates; G6 agent must test collinearity.
2. **`AlleleDelta` collider**: Hard-coded memory — L2 collider bias w/ AF; any AUC≥0.58 must pass `/auc-confound-guard` with within-AF-bin residualization.
3. **`Coverage_Multiple` hardcoded 75.0 default**: `compute_coverage_multiple` defaults `expected_coverage=75.0` (RegionProcessor.cpp:522); per-sample diploid_coverage auto-estimation may not be active for all runs — check `Diploid_Coverage_Used` column.
4. **`HPMergedDelta` / `Tumor_HP_Delta`**: unsigned absolute values; signed versions (`Tumor_HP_Signed_Delta`, `HP_Signed_Residual`, `Combined_HP_Signed_Delta`) preserve direction for downstream ASM interpretation.
5. **`NHP_Somatic11/21/33`**: captured **regardless** of `--germline-hp-only` flag (audit trail). Do not confuse with `HPFineN_HP1S` which respects the flag.
6. **`Fisher_Frac_Sig`**: highest known paired-track AUC (~0.726), but memory note: `characterization-only; CI crosses random for F1 filter gain`.
7. **G3 VCF features**: currently NOT joined in master TSV — Phase B1 script needed before G3 agent can run.
8. **G4 BAM features**: all INFERRED from naming convention — definitions must be pinned when `scripts/analysis/bam_readqc_per_region.py` is authored in Phase B2.
9. **`Combined_HP_Signed_Delta`**: exact aggregation weighting (simple mean? coverage-weighted?) is not obvious from writer alone; the G9 agent should re-grep `combined_hp_signed_delta =` in `RegionProcessor.cpp` to confirm before publishing formula.
10. **`VerificationClass`** and **`LOH_Subtype`** are derived — do not treat as independent signals. `LOH_Subtype = f(Potential_LOH, VerificationClass)` per `RegionProcessor.cpp:188`.

---

## Downstream usage

Phase C group agents (G1..G10) should:

```python
from research.feature_layered_observation.scripts.feature_registry import (
    get_feature_metadata, list_features_by_group
)

features = list_features_by_group("G6")  # ['GlobalP_HPFine', 'HPFineNGroups', ...]
for fid in features:
    md = get_feature_metadata(fid)
    print(md["computation"], md["source_file"], md["filter_condition"])
```

The registry TSV is the single source of truth; do **not** hardcode feature lists in analysis scripts.
