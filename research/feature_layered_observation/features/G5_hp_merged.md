# G5 · HP Merged (2-bucket) — Feature Layered Observation

**Date**: 2026-04-23
**Scope**: 7 samples × 2 modes (paired_full, to_pileup); 748,676 rows in master, 697,078 (93.1%) joined with G5 columns. COLO829 to_pileup (50,617 rows) has no archive — G5 NaN in that cell.
**Registry**: `research/feature_layered_observation/data/feature_registry.tsv` group G5 (7 features)
**Methodology**: `research/feature_layered_observation/02_methodology.md`
**Data source**: `data/G5/master_g5.tsv.gz` (join of `merged_with_vcf.tsv.gz` + per-sample `significance_summary.csv`)

---

## 1. Feature definition & source

| Feature | C++ source | Type | Definition |
| --- | --- | --- | --- |
| HP_Ratio | RegionProcessor.cpp:39 `compute_hp_ratio` | continuous [0,1] | `HP1_count / (HP1_count + HP2_count)` — no HP-family merge |
| HPMergedDelta | RegionProcessor.cpp:1621 `hp_ms.merged.delta` | continuous ≥0 | `|mean_HP1family − mean_HP2family|` on HP-family buckets {HP1 ∪ HP1-1} vs {HP2 ∪ HP2-1} |
| HPMergedP | LabelTest.cpp | continuous [0,1] | Permutation p for merged-family PERMANOVA |
| HPMergedSig | LabelTest.cpp | binary | `HPMergedP ≤ 0.05 ∧ validity` |
| HP1FamilyN | RegionProcessor.cpp hp_merged_n_hp1 | ordinal | Count HP1 ∪ HP1-1 reads |
| HP2FamilyN | RegionProcessor.cpp hp_merged_n_hp2 | ordinal | Count HP2 ∪ HP2-1 reads |
| GlobalP_HPFamily | GlobalTest.cpp:252 | continuous | Fisher FFH p on cluster × HP-family |
| CramersV_HPFamily | GlobalTest.cpp | continuous | Cramér's V effect size |

**Derived (for this study)**
- `HP_Ratio_norm = |HP_Ratio − 0.5|`   — symmetric extremeness, 0 = balanced, 0.5 = full LOH
- `HP_FamilyN_sum = HP1FamilyN + HP2FamilyN`  — should ≤ NumReads
- `HP_FamilyN_frac = HP_FamilyN_sum / NumReads`   — coverage of HP-family assignment

---

## 2. Observation goals

Answer the three spec questions plus the feature group's Verdict under `02_methodology.md` rubric.

1. Does `min(HP_Ratio, 1-HP_Ratio)` (equivalently `|HP_Ratio−0.5|`) give tighter TP/FP separation than raw HP_Ratio?
2. Is `HPMergedDelta` essentially `|AlleleDelta|` (AF echo) → L2 collider suspect?
3. Is `HP1FamilyN + HP2FamilyN` a proxy for `NumReads`?

---

## 3. Global distribution (Step 1)

See `figures/G5_hp_merged/fig01_global_distribution.png`.

| Feature | AUC (Wilson 95% CI) | Cohen d | mean_TP | mean_FP | Verdict | 
| --- | --- | --- | --- | --- | --- |
| HP_Ratio | 0.519 [0.517, 0.521] | +0.053 | 0.538 | 0.521 | NEGATIVE |
| **HP_Ratio_norm** | **0.466 [0.464, 0.467]** | −0.113 | 0.264 | 0.285 | **FLIPPED (below random)** |
| HPMergedDelta | 0.481 [0.479, 0.483] | −0.052 | 0.014 | 0.016 | NEGATIVE |
| HPMergedP | 0.516 [0.514, 0.518] | +0.058 | 0.561 | 0.533 | NEGATIVE |
| HP1FamilyN | 0.530 [0.528, 0.532] | +0.074 | 38.6 | 36.0 | NEGATIVE (near NumReads effect) |
| HP2FamilyN | 0.507 [0.505, 0.509] | +0.029 | 34.5 | 33.4 | NEGATIVE |
| **HP_FamilyN_sum** | **0.567 [0.565, 0.569]** | +0.187 | 69.1 | 60.1 | BELOW 0.58 (NumReads proxy) |
| HP_FamilyN_frac | 0.387 [0.385, 0.389] | −0.269 | 0.631 | 0.725 | FLIPPED (TP has less HP coverage) |
| **GlobalP_HPFamily** | **0.577 [0.575, 0.579]** | +0.361 | 0.839 | 0.719 | BELOW 0.58 threshold |
| CramersV_HPFamily | 0.493 [0.491, 0.495] | −0.109 | 0.015 | 0.027 | NEGATIVE |
| HPMergedSig | 0.487 [0.485, 0.489] | — | rate_TP=0.313 | rate_FP=0.338 | NEGATIVE |

**No G5 feature crosses the 0.58 Beyond-AUC ceiling globally.**

Notable directional findings:
- `HP_Ratio_norm` is **below 0.5** globally — TP has **lower** extremeness than FP. This is counterintuitive but explained by LOH-conditional direction flipping (see §4).
- `HP_FamilyN_frac` < 1 for TP (0.63) more than for FP (0.73) — TP tends to have more HP-unassigned reads (HP0/HP3) in families.

---

## 4. LOH × AF × CN stratified observation (Step 2–4)

See `figures/G5_hp_merged/fig04_HP_Ratio_stratified_auc.png`, `fig04_HP_Ratio_norm_stratified_auc.png`.

### HP_Ratio stratified AUCs

| Stratum | AUC | n_pos | n_neg | Interpretation |
| --- | --- | --- | --- | --- |
| global | 0.519 | 583k | 114k | near random |
| LOH_Weak | **0.620** | 68,756 | 6,447 | mild elevation — but TP rate in this stratum is 0.914 (base rate inflation) |
| LOH_Noise | 0.459 | 86,665 | 24,817 | reversed |
| LOH_Strong | **0.628** | 35,438 | 2,980 | strong direction (see below) |
| LOH_Subclone | 0.402 | 7,318 | 925 | reversed |
| AF_Extreme | 0.521 | 522,914 | 106,832 | near random |
| AF_Near-half | 0.517 | 6,343 | 170 | low-FP; unstable |
| AF_Intermediate | 0.517 | 53,772 | 7,047 | near random |
| mode=paired_full | 0.534 | 325,507 | 3,458 | near random (extreme imbalance) |
| mode=to_pileup | 0.552 | 257,522 | 110,591 | slight |

The per-LOH direction flip (Weak/Strong > 0.6 but Noise/Subclone < 0.5) confirms Simpson's paradox pattern.

### HP_Ratio_norm stratified AUCs — the "hidden signal" puzzle

| Stratum | AUC | n_pos | n_neg |
| --- | --- | --- | --- |
| global | 0.466 | 583k | 114k |
| AF_Extreme | 0.441 | 522,914 | 106,832 |
| **AF_Near-half** | **0.782** | 6,343 | 170 |
| **AF_Intermediate** | **0.682** | 53,772 | 7,047 |
| CN_Diploid | 0.431 | 238,466 | 40,764 |

AF_Near-half AUC=0.782 looks very strong — but the per-sample breakdown destroys it.

### HP_Ratio_norm per-sample per-AF (n≥100)

```
AF_class          Extreme  Intermediate  Near-half
HCC1395        paired_full   0.422    0.644         NaN    (n_neg too small)
HCC1395        to_pileup     0.522    0.536         0.563
HCC1395_DORADO paired_full   0.469    0.452         0.402
HCC1395_DORADO to_pileup     0.516    0.621         NaN
HCC1937        paired_full   0.392    0.464         NaN
HCC1937        to_pileup     0.535    0.681         NaN
HCC1954        paired_full   0.310    NaN           NaN
HCC1954        to_pileup     0.519    0.519         0.610
H2009          paired_full   0.219    0.267         NaN
H2009          to_pileup     0.520    0.570         0.520
H1437          paired_full   0.397    NaN           NaN
H1437          to_pileup     0.518    0.546         NaN
COLO829        paired_full   0.450    0.515         NaN
```

- **Near-half per-sample**: 1/7 > 0.58, 0/7 > 0.65. Median AUC 0.542.
- **Intermediate per-sample**: 3/13 > 0.58, 1/13 > 0.65. Median AUC 0.536.

**Conclusion**: the pooled Near-half AUC=0.782 is a pooling artifact, driven by the extreme class imbalance (170 FP vs 6,343 TP) and a handful of samples (HCC1937 TO AUC=0.681, HCC1395_DORADO TO Intermediate AUC=0.621). Not a cross-sample robust signal.

### Near-half AF decomposed by LOH

```
AF=Near-half  n=6,513  tp_rate=0.974
  LOH=0  n=1,130  tp_rate=0.881  median_norm_TP=0.106  FP=0.165
  LOH=1  n=5,383  tp_rate=0.993  median_norm_TP=0.500  FP=0.500
```

Inside LOH=1 rows, HP_Ratio is clamped at 0 or 1 (norm=0.5) for both TP and FP — no discrimination at all. Discrimination only exists in the LOH=0 subset, where TP has lower extremeness (norm 0.106) than FP (0.165). This is **NumReads × AF ceiling effect**, not a methylation signal: FP with Near-half AF and no LOH are usually germline/low-evidence sites whose few reads skew to one side.

### HPMergedDelta stratified

Uniformly ≤ 0.51, with AF_Near-half AUC = 0.184 (strongly flipped). No discriminative power.

---

## 5. Cross-sample consistency (Step 3)

See `fig03_*_per_sample.png`. Spearman concordance matrix (LOH × AF cell TP rate across samples) shows moderate positive agreement (most cells 0.5–0.9), confirming the LOH × AF stratification itself is a stable TP-rate structure, but this is a property of LOH/AF, not G5 features.

Per-sample HP_Ratio_norm AUC is inconsistent: direction flips between samples (see §4 table).

---

## 6. Stratified AUC (Step 4)

See `fig04_*.png`. Full table in `data/G5/G5_auc_table.tsv`. Summary:

- No feature has AUC > 0.58 in any layer that generalizes across samples.
- HP_Ratio LOH_Weak/LOH_Strong AUC ≈ 0.62 looks like a stratum-specific lift, but the direction (more LOH = TP) is base-rate driven (TP rate in LOH_Weak is 0.91, LOH_Strong 0.92 vs global 0.84).

---

## 7. Confound check (Step 5)

See `fig05_*_confound.png` and `data/G5/G5_confound.tsv`.

| Feature | Raw AUC | Residualized AUC (on AlleleDelta, NumReads, vcf_AF, Coverage_Multiple) | AF_Extreme raw | AF_Near-half raw | AF_Intermediate raw |
| --- | --- | --- | --- | --- | --- |
| HP_Ratio | 0.519 | 0.546 | 0.521 | 0.517 | 0.517 |
| HP_Ratio_norm | 0.466 | **0.532** | 0.441 | 0.782 | 0.682 |
| HPMergedDelta | 0.481 | 0.436 (dropped) | 0.495 | 0.184 | 0.365 |

- `HP_Ratio_norm` residualized AUC (0.532) is still sub-threshold.
- `HPMergedDelta` residualization **decreased** AUC — confirms its residual weight beyond AlleleDelta/NumReads/vcf_AF is anti-correlated with TP after removing AF component.

### Q2 specific check — HPMergedDelta vs AlleleDelta collinearity

```
Pearson(HPMergedDelta, AlleleDelta) = 0.209
Spearman(HPMergedDelta, AlleleDelta) = 0.225
Pearson(HPMergedDelta, vcf_AF) = −0.146
Spearman(HPMergedDelta, vcf_AF) = −0.194
```

Correlation is **low** (< 0.25), so HPMergedDelta is *not* a simple AF echo. But it is mildly anti-correlated with vcf_AF, which combined with the AF_Near-half AUC=0.184 flip shows HPMergedDelta is confounded with caller AF structure. Not a collider in the strong L2 sense, but residualized AUC (0.436) below raw (0.481) indicates the covariates absorb even the weak signal it carries.

### Q3 specific check — HP1FamilyN+HP2FamilyN vs NumReads

```
Pearson(HP_FamilyN_sum, NumReads) = 0.677
Spearman(HP_FamilyN_sum, NumReads) = 0.522
HP_FamilyN_frac median: 0.612
  fraction frac > 0.95: 0.188
  fraction frac < 0.20: 0.141
```

Correlation 0.68 is substantial but not saturating. Only 19% of regions have ≥95% HP-family coverage; 14% have <20% (likely high HP0/HP3 regions in non-LOH territory). HP_FamilyN_sum carries NumReads signal (raw AUC 0.567 ≈ NumReads AUC per G1) but should be dropped in favor of NumReads itself.

---

## 8. Spatial auto-correlation (Step 6)

See `fig06_*_spatial.png` and `data/G5/G5_spatial.tsv`.

- HP_Ratio, HP_Ratio_norm, HPMergedDelta: 566 5Mb bins tested, no artifact flag (|ΔAUC hi_tp vs lo_tp| < 0.08).

No spatial red flag.

---

## 9. Prior conclusions & knowledge-base context

- `HP_Ratio` extreme behavior in LOH (≈ 0 or 1) is exactly the design of `Potential_LOH = (HP_Ratio < 0.1 OR > 0.9)` gate in `RegionProcessor.cpp:188`. G5 `HP_Ratio` therefore cannot provide **independent** LOH signal — it is the primary LOH detector input, and `Potential_LOH` is derived from it.
- `HPMergedDelta` parallels `Combined_HP_Signed_Delta` (G9) but using HP-family bucket means. G9 has been shown to provide no discrimination on paired after AF residualization.
- `HPMergedSig` is currently **characterization-only**: paired_full TP rate = 0.989 vs 0.988 (nearly identical), in TO the class TP rate actually flips (sig=1 TP rate = 0.662 < sig=0 TP rate = 0.716).

Prior memory flags relevant here:
- `feedback_L2_collider_bias` — residualize on AF before trusting AUC > 0.58. Followed.
- `project_beyond_auc_exhaustion_confirmed` — pure methylation features ≤ 0.58. Consistent with G5.
- `feedback_feature_name_vs_definition_rule` (2026-04-22) — I confirmed HP_Ratio is raw HP1/(HP1+HP2), not signed; HPMergedDelta uses unsigned abs of family means.

### 9.1 Knowledge base references (Phase D)

Queries: `haplotype phasing` (top_score 126.5, full), `HP tag merged` (top_score 33.6, partial). All high confidence.

| kb_path | kb_title | Relevance to G5 |
|---|---|---|
| `03_file_formats/bam-format.md` | BAM 格式文件 | Defines HP tag storage (`HP:i:1` / `HP:i:2` / `HP:i:3`); direct spec for how `HP_Ratio` numerator/denominator is read from BAM |
| `05_tools/longphase.md` | LongPhase (germline) | Germline phasing emits the baseline HP1/HP2 tags; `HP_Ratio_norm` assumes this symmetric baseline — if germline phasing is incomplete, HP_Ratio_norm biases |
| `03_file_formats/vcf-longphase.md` | Phased VCF 規格 | GT pipe / PS phase set defines how phasing blocks segment; explains why `HPMergedDelta` can be NaN when a region spans a phase-set break |
| `06_workflows/phasing-workflow.md` | Phasing 工作流程 | Merged HP tag is paired-mode flavor (LongPhase-S HP 1/2/1-1/2-1/3 collapsed to HP1 vs HP2 family); canonical reference for the "family merge" step |
| `05_tools/longphase-to.md` | LongPhase-TO | TO mode phased VCF outputs HP1/HP2/HP2-1 — the upstream of TO-track HPMergedDelta; explains TO-specific TP rate inversion (sig=1 < sig=0) as a phased-TO artefact |

**HP_Merged**: 2/2 topic queries high-confidence. Gap: KB does not document the G5 `HPMergedDelta = |mean(HP1_family_methyl) − mean(HP2_family_methyl)|` aggregation rule, or `Combined_HP_Signed_Delta` sign convention. Recommend Phase F KB authoring for the HP-family merge semantics (especially the somatic/germline HP merge under `--germline-hp-only` flag).

### 9.2 External literature (Phase D)

**Haplotype-aware somatic calling references:**

1. **Shafin, K., Pesout, T., Chang, P.-C. et al. (2021).** "Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads." *Nature Methods* 18, 1322–1332. DOI: **10.1038/s41592-021-01299-w** — establishes that **HP-aware** deep-learning SNV calling outperforms HP-agnostic calling, especially in segmental duplications and low-mappability regions. Relevance to G5: confirms that HP tags carry real per-read information that ClairS/ClairS-TO already exploit internally — which is precisely why aggregated HP family means (HPMergedDelta) leave **no independent signal** once caller QS is controlled for. **Supports** our NEGATIVE verdict on HPMergedDelta as a filter.

2. **Lin, J.-H., Chen, L.-C., Yu, S.-C. & Huang, Y.-T. (2022).** "LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants." *Bioinformatics* 38(7), 1816–1822. DOI: **10.1093/bioinformatics/btac058** — LongPhase is the upstream tool producing the HP1/HP2 tags that G5 aggregates. Co-phasing SNPs+SVs yields N50 = 25 Mb haplotype blocks, ~10× faster than WhatsHap/HapCUT2/Margin. Relevance to G5: `HP_Ratio_norm` directly depends on LongPhase block symmetry — **phase-set break at a region boundary** is a known failure mode producing bimodal HP_Ratio. **Supports** our interpretation that HP_Ratio extremes are LOH-equivalent (and thus not independent of Potential_LOH).

3. **Fu, Y., Aganezov, S., Mahmoud, M. et al. (2024).** "MethPhaser: methylation-based long-read haplotype phasing of human genomes." *Nature Communications* 15, 4859. DOI: **10.1038/s41467-024-49588-0** — extends SNV phasing by methylation signal, increasing phase length N50 by 78–151% at 83–99% accuracy. Relevance to G5: demonstrates that **methylation + HP tag is a meaningful joint signal** at the read level — but only across phase-set breaks, not within blocks. **Challenges** our pooled G5 analysis by suggesting that any residual HPMergedDelta signal should concentrate at phase-set boundaries; future G5 reanalysis should stratify by distance to phase-set break.

**Direction summary**: Papers 1 and 2 *support* our NEGATIVE verdict (HP information is caller-internal and LOH-entangled). Paper 3 *challenges* the pooled aggregation by pointing at phase-set boundaries as the remaining potential value.

---

## 10. Verdict, challenges, follow-ups

### Verdict (per `02_methodology.md` rubric)

**NEGATIVE** — no G5 feature crosses 0.58 global AUC. Best pooled feature (`GlobalP_HPFamily`, AUC 0.577) sits below threshold and is a p-value (TP has larger p for HP-family structure, i.e. less significant HP stratification — backwards from the intuition that TP = structured methylation).

### Three challenges (質疑)

1. **The "huge" Near-half AUC=0.782 for HP_Ratio_norm is a sampling artifact.** Class imbalance (6,343 TP vs 170 FP) and 1/7 per-sample AUC > 0.58 mean the pooled number is driven by a minority. Reject as signal.
2. **HP_Ratio is not an independent G5 feature.** It is the input to `Potential_LOH`; using it as a TP/FP discriminator downstream leaks the LOH stratification dimension. The apparent LOH_Weak/Strong AUC ≈ 0.62 is recoverable from LOH_Subtype alone (per-stratum TP rate differences).
3. **HPMergedDelta has non-zero but low correlation with AlleleDelta (r≈0.21).** It is not a pure L2 collider in the "AF echo" sense we feared, but its residualized AUC drops below raw (0.436 vs 0.481), so the covariate set already captures its weak TP signal. No independent contribution.

### Logic chain

Global null (AUC ~0.5 for every feature) → stratified lifts (LOH 0.62, AF_Near-half 0.78 for HP_Ratio_norm) → per-sample null (median ~0.54 within AF stratum) → residualized null → spatial null → **no independent TP/FP discriminative power**. HP-family features restate information already in Potential_LOH + NumReads + AF.

### Follow-ups / recommendations

1. **Drop G5 continuous features from filter candidate pool.** Keep them as *characterization* axes (family-bucket p-value maps for the paper's ISM figure).
2. **HP_Ratio normalization suggestion**: do not use `|HP_Ratio−0.5|` as a TP/FP score. It only works in a narrow AF+LOH sub-cell that is already captured by `LOH_Subtype` + `AF_class` jointly. The true utility of HP_Ratio is for the binary LOH gate, not a continuous TP score. For visualization / annotation purposes, stratify by `LOH_Subtype` first, then show HP_Ratio histograms within each LOH class.
3. **HPMergedSig reversal in TO mode** deserves a narrow follow-up: sig=1 TP rate is **lower** than sig=0 in to_pileup (0.662 vs 0.716). This is the "TO germline → HP-structure → low TP rate" pattern already characterized by ClairS-TO self-phasing audits. Cross-reference with `project_self_phasing_causal_chain_confirmed`.
4. **GlobalP_HPFamily** shows a weakly below-threshold signal (AUC 0.577, d=0.361). Because it is a p-value and direction is "higher p → more TP", it is likely a NumReads / NumCpGs power artifact (higher n → smaller p only if signal exists; low-power regions fail to reject). Consider dropping or restating with `CramersV_HPFamily` instead.
5. Next: **G6 (HP Fine 4-bucket)** — this is the Thread D core, where LOH-constrained phasing lives. HP-merged (G5) is a lossy aggregation; expect G6 to recover the structure that G5 collapses.

### Notes / caveats

- COLO829 to_pileup is missing (archive step05 empty); G5 coverage 93.1% not 100%.
- HCC1395 to_pileup came from `bip8_output_archive/20260307_hcc1395_to_pilot_1/step05_intersubmod/` — lacks `GlobalP_HPFamily` and `CramersV_HPFamily` columns (older ISM build), hence 87.9% coverage for those two.
- All AUC CIs are Hanley-McNeil (conservative). Wilson CI was used for proportions.

---

## Deliverables

- **Figures**: `research/feature_layered_observation/figures/G5_hp_merged/fig00..fig08` (23 PNG)
- **Tables**: `research/feature_layered_observation/data/G5/G5_*.tsv`
  - `G5_global_stats.tsv` — Step 1 pooled stats
  - `G5_relations.tsv` — Step 0 sanity correlations
  - `G5_auc_table.tsv` — Step 4 stratified AUCs
  - `G5_cell_delta.tsv` — Step 2 cell-level TP-FP feature deltas
  - `G5_confound.tsv` — Step 5 residualization
  - `G5_spatial.tsv` — Step 6 spatial screen
  - `G5_hp_ratio_normalization.tsv` — Step 7 normalization benchmark
  - `G5_HPMergedSig_per_class.tsv` — Step 8 per-class TP rate
  - `G5_HP_Ratio_norm_persample_AF.tsv` — per-sample AUC by AF_class
- **Scripts**: `scripts/g5_build_extended_master.py`, `scripts/observe_G5_hp_merged.py`, `scripts/g5_persample_ratio_norm_af.py`
- **Master**: `data/G5/master_g5.tsv.gz` (748,676 rows × 68 cols, 697,078 G5-joined)
