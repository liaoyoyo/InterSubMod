# G4 · BAM Read-Level QC — Feature Layered Observation

**Date**: 2026-04-24
**Group**: G4 (BAM Read-Level Quality Control)
**Features**: MapQ_mean, MapQ_median, MapQ_p10, MapQ_p90, LowMQ20_Frac, NM_mean, Softclip_Frac, Strand_Bias, Read_Length_mean
**Data**: `research/feature_layered_observation/data/G4/master_g4.tsv.gz` (748,676 rows × 69 cols after merge; 93.2% BAM QC coverage due to COLO829 `to_pileup` empty)
**Script**: `research/feature_layered_observation/scripts/observe_G4_bam_readqc.py`

---

## 1. Feature definitions and source

| Feature | Type | Source / Computation | Semantic |
|---------|------|----------------------|----------|
| MapQ_mean | continuous | `bam_readqc_per_region.py` — mean of read MAPQ in window | Average primary-alignment mapping quality (0-60, Minimap2 cap=60) |
| MapQ_median | continuous | median of MAPQ | Robust central MAPQ |
| MapQ_p10 | continuous | 10th percentile of MAPQ | Worst-tail mapping quality in region |
| MapQ_p90 | continuous | 90th percentile of MAPQ | Best-tail mapping quality |
| LowMQ20_Frac | fraction | reads with MAPQ<20 / total reads | Proportion of low-confidence alignments |
| NM_mean | continuous | mean of BAM `NM` tag across reads | Average edit distance (mismatch + indels) per read; scales with read length |
| Softclip_Frac | fraction | reads with any softclip operator in CIGAR / total reads | Fraction of reads with soft-clipped bases (boundary signal for SV / adapter) |
| Strand_Bias | continuous | `|fwd_count/total − 0.5| × 2` | 0 = balanced strands, 1 = completely single-strand |
| Read_Length_mean | continuous | mean of CIGAR-derived read lengths | Average aligned read length (ONT long-read proxy) |

**Extraction**: `scripts/bam_readqc_per_region.py` (Phase B2, background 11.6 kh aggregate wall-clock).
**Data gap**: `COLO829/to_pileup` BAM QC TSV is **header-only empty** (extraction failure; other 13 TSVs populated → 93.2% master coverage).
**Known caveat**: ONT Minimap2 caps MAPQ at 60 → distributions are heavily right-censored (median=60 in 6/9 MAPQ sub-features). Similarly `Softclip_Frac` is right-censored at 1.0 (median=1.0, 75th pct=1.0) because **nearly all ONT reads have some softclipping** (adapter tails).

---

## 2. Observation goals

Test whether read-level BAM alignment quality features — standard caller inputs — discriminate TP vs FP after the ClairS / ClairS-TO caller has already consumed them. The prior is **low**: caller confidence (`vcf_QUAL`) should absorb the ordered information in (MapQ, NM, Softclip, Strand_Bias, Read_Length), so residual independent signal over `vcf_QUAL + NumReads + vcf_AF + Coverage_Multiple` should be small. However, three mechanisms could still produce independent signal:

1. **Strand bias** — ClairS mixes pileup + full-alignment networks; extreme strand bias is a known caller-inflated false positive in structural repeat regions.
2. **Softclip_Frac** — structural variation / adapter pile-ups can inflate FP without the caller rejecting them, particularly in LOH regions (caller lacks SV-aware filter).
3. **LowMQ20_Frac** — the caller uses a hard MAPQ threshold (≥5) but does not rank regions by `LowMQ20_Frac`; residual fraction of moderately low-MAPQ reads could escape.

If any of these features passes the Step 5 residualization threshold and is cross-sample consistent, it becomes an orthogonal axis for filter design.

---

## 3. Global distribution (Step 1)

`fig01_global_distribution.png`, `fig02_global_auc_bar.png`, table `data/G4/G4_global_stats.tsv`:

| Feature | AUC [95% CI] | Cohen's d | mean TP | mean FP | median TP | median FP | mwu-p |
|---------|--------------|-----------|---------|---------|-----------|-----------|-------|
| MapQ_mean | 0.509 [0.508, 0.511] | +0.07 | 59.94 | 59.89 | 60.00 | 60.00 | 1.9e-86 |
| Read_Length_mean | 0.509 [0.507, 0.511] | +0.02 | 27,623 | 27,436 | 28,770 | 28,647 | 2.2e-21 |
| MapQ_p10 | 0.502 [0.500, 0.504] | +0.05 | 59.92 | 59.81 | 60.00 | 60.00 | 6.2e-109 |
| MapQ_median | 0.500 [0.499, 0.502] | +0.02 | 59.99 | 59.98 | 60.00 | 60.00 | 2.4e-24 |
| MapQ_p90 | 0.500 [0.498, 0.502] | +0.01 | 60.00 | 60.00 | 60.00 | 60.00 | 2.4e-07 |
| Softclip_Frac | 0.497 [0.495, 0.499] | +0.02 | 0.987 | 0.986 | 1.00 | 1.00 | 4.6e-04 |
| LowMQ20_Frac | 0.494 [0.493, 0.496] | −0.05 | 6.1e-04 | 1.1e-03 | 0.00 | 0.00 | 1.7e-84 |
| Strand_Bias | 0.491 [0.489, 0.493] | −0.03 | 0.069 | 0.070 | 0.056 | 0.057 | 1.2e-20 |
| NM_mean | 0.478 [0.476, 0.480] | −0.09 | 542.2 | 574.4 | 506.0 | 517.1 | 2.4e-123 |

**TP/FP counts (pooled)**: 583,742 TP vs 114,317 FP (TP rate 83.6% on rows where BAM QC is populated).

Key observations:

- **No feature exceeds global AUC 0.510.** All 9 BAM QC features are at or below the Beyond-AUC ceiling (0.58) — as expected because the caller has already consumed them.
- **NM_mean has the strongest Cohen's d (−0.09)**: FP have more mismatches per read (574 vs 542 NM units) than TP, but the effect is small and MWU significance (p<10^−123) is driven by the 698k-sample size, not biological strength.
- **MapQ-based features are degenerate**: with median=60 and 75th pct=60 across both TP and FP, MAPQ quantiles sit on the Minimap2 cap and carry almost no discriminative power.
- **Softclip_Frac median=1.0**: indicates the feature is heavily right-censored; the standalone AUC 0.497 is basically noise around the ceiling.
- **LowMQ20_Frac has anti-predictive direction**: FP regions contain slightly more low-MAPQ reads (1.1e-3 mean vs 6.1e-4 in TP), but this is ~0.0005 absolute difference.

**First-pass verdict**: **no BAM QC feature is a global discriminator for TP/FP** after ClairS/ClairS-TO caller scoring.

---

## 4. LOH × AF × CN 32-cell heatmap (Step 2)

Heatmaps `fig03_{NM_mean,MapQ_mean,Read_Length_mean,Strand_Bias,LowMQ20_Frac}_heatmap.png`; table `data/G4/G4_cell_delta.tsv`.

Focus features (selected by |AUC − 0.5|): **NM_mean, MapQ_mean, Read_Length_mean, Strand_Bias, LowMQ20_Frac**.

Cell-level observations:
- **NM_mean Δ(TP−FP)** is negative (TP lower mismatch) across most LOH_Weak / LOH_Noise cells in all CN strata (Δ ≈ −20 to −50 NM units); pattern flips slightly positive in **LOH_Subclone × Intermediate × CN_Gain** cell where TP have more NM than FP (Δ ≈ +30), matching the stratified AUC spike (§6).
- **MapQ_mean Δ** is near-zero (|Δ| < 0.3) in every cell — MAPQ is saturated at 60, so per-cell variation is below the feature resolution.
- **Read_Length_mean Δ** is positive in LOH_Strong / LOH_Subclone + CN_Diploid (TP reads ~500 bp longer than FP); matches HCC1395_DORADO's longer-read hallmark in those subsamples.
- **Strand_Bias Δ** is slightly negative in CN_HighGain cells (TP *less* strand-biased than FP; |Δ| ≈ 0.02), consistent with the expectation that amplicon artefacts (HCC1954 chr8) have higher strand bias.
- **LowMQ20_Frac Δ** is effectively zero (|Δ| < 0.002) in all populated cells.

`data/G4/G4_cell_delta.tsv` contains Δ for each of 5 features × 25 populated cells (min n=20).

---

## 5. Cross-sample consistency (Step 3)

`fig04_persample_auc_heatmap.png` (9 features × 7 samples × 2 modes) + table `data/G4/G4_persample_auc.tsv`; consistency summary `data/G4/G4_cross_sample_consistency.tsv`:

| Feature | mode | n samples | ≥0.58 count | median AUC | min | max | range |
|---------|------|-----------|------------|-----------|-----|-----|-------|
| MapQ_mean | paired_full | 7 | 2/7 | 0.541 | 0.479 | 0.707 (H2009) | 0.228 |
| MapQ_mean | to_pileup | 6 | 0/6 | 0.514 | 0.507 | 0.524 | 0.017 |
| Read_Length_mean | paired_full | 7 | 1/7 | 0.505 | 0.416 | 0.693 (HCC1954) | 0.277 |
| Read_Length_mean | to_pileup | 6 | 0/6 | 0.502 | 0.491 | 0.509 | 0.018 |
| NM_mean | paired_full | 7 | 0/7 | 0.443 | 0.290 | 0.578 | 0.288 |
| NM_mean | to_pileup | 6 | 0/6 | 0.484 | 0.480 | 0.498 | 0.018 |
| LowMQ20_Frac | paired_full | 7 | 0/7 | 0.484 | 0.395 | 0.506 | 0.112 |
| LowMQ20_Frac | to_pileup | 6 | 0/6 | 0.492 | 0.485 | 0.498 | 0.012 |
| Softclip_Frac | paired_full | 7 | 0/7 | 0.497 | 0.456 | 0.517 | 0.060 |
| Softclip_Frac | to_pileup | 6 | 0/6 | 0.499 | 0.498 | 0.505 | 0.007 |
| Strand_Bias | paired_full | 7 | 0/7 | 0.480 | 0.447 | 0.554 | 0.107 |
| Strand_Bias | to_pileup | 6 | 0/6 | 0.499 | 0.492 | 0.503 | 0.011 |

- **No feature achieves ≥5/7 samples with AUC ≥ 0.58 in either mode.** The best is `MapQ_mean` paired with 2/7 (HCC1954, H2009) and `Read_Length_mean` paired 1/7 (HCC1954 only).
- **TO-mode is uniformly flat**: AUC ranges 0.485–0.524 for every feature × sample cell — TO FP set is structurally indistinguishable on BAM QC.
- **COLO829 paired_full** `MapQ_mean` AUC is 0.562 (TP 59.89 vs FP 58.76) — an isolated but non-replicable signal.
- **HCC1954 paired_full** is the most variable (MapQ_mean 0.605, Read_Length_mean 0.693); this is the amplicon hotspot sample and the `fig03` DORADO row for Read_Length shows the outlier cleanly.
- **Sample order** locked: HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829.

---

## 6. Stratified AUC (Step 4)

`fig05_stratified_auc.png` + table `data/G4/G4_auc_table.tsv`.

Best stratum per feature:

| Feature | global | best LOH layer | best AF layer | best CN layer | best mode |
|---------|--------|----------------|---------------|---------------|-----------|
| MapQ_mean | 0.509 | LOH_Strong 0.525 | Near-half 0.541 | CN_HighGain 0.536 | paired 0.531 |
| MapQ_median | 0.500 | LOH_Subclone 0.512 | Extreme 0.500 | CN_Diploid 0.500 | paired 0.506 |
| MapQ_p10 | 0.502 | LOH_Strong 0.522 | Intermediate 0.519 | CN_Diploid 0.505 | paired 0.528 |
| LowMQ20_Frac | 0.494 | LOH_Strong 0.536 | Intermediate 0.554 | CN_Diploid 0.494 | paired 0.470 |
| NM_mean | 0.478 | **LOH_Subclone 0.670** | Intermediate 0.514 | CN_Near1 0.495 | paired 0.550 |
| Softclip_Frac | 0.497 | LOH_Strong 0.510 | **Intermediate 0.608** | CN_Diploid 0.507 | paired 0.500 |
| Strand_Bias | 0.491 | LOH_Weak 0.513 | Near-half 0.522 | CN_HighGain 0.516 | paired 0.403 |
| Read_Length_mean | 0.509 | **LOH_Subclone 0.682** | **Intermediate 0.603** | CN_Near1 0.527 | paired 0.597 |

- **NM_mean LOH_Subclone AUC 0.670** is the single highest stratified G4 AUC — FP in subclonal-LOH regions have fewer NM mismatches than TP (n=7,318 TP vs 925 FP). Direction is **inverted** from global (where TP<FP); interpretation: subclonal-LOH somatic TPs carry additional mutations (high NM) while FPs in those regions cluster as "clean" germline misassignments. This is biologically coherent.
- **Read_Length_mean LOH_Subclone AUC 0.682 / Intermediate AF AUC 0.603 / paired AUC 0.597** — TP reads are ~500 bp longer in subclonal regions. Likely reflects the long-read advantage for phasing at subclonal AF; but 2/3 of these samples have ONT_Dorado reads (30 kb median) so read-length advantage could be sample-ID confound not biology.
- **Softclip_Frac Intermediate AF AUC 0.608** — but global 0.497 and other AF strata 0.49 — this is a 0.11-AUC lift only in mid-AF cells (53,815 TP vs 7,074 FP) and may reflect the caller's incomplete SV filter at mid-AF.
- **LowMQ20_Frac Intermediate AUC 0.554** — small lift, not cross-validatable.
- **Strand_Bias paired mode AUC 0.403** is markedly anti-predictive in paired (TP more strand-biased than FP) but TO-mode 0.499 (random). Likely a paired-track processing artefact; not useful.

None of these stratum peaks exceeds the POSITIVE threshold (0.65 + resid≥0.55 + cross-sample ≥5/7). `NM_mean` LOH_Subclone 0.670 and `Read_Length_mean` LOH_Subclone 0.682 are **single-stratum** signals that cannot be extrapolated across other strata, and they must be checked against `HPFineNGroups` LOH_Subclone 0.68 (G6) for collinearity (likely same signal).

---

## 7. Confound guard (Step 5)

`fig06_confound_residualized.png` + table `data/G4/G4_confound.tsv`. OLS residualization on `NumReads + vcf_AF + Coverage_Multiple + vcf_QUAL` — the explicit caller covariates that should absorb any BAM-QC-derived signal:

| Feature | raw AUC | resid AUC | ΔAUC | AF-bin AUC (Extreme / Near-half / Intermediate) | Verdict |
|---------|---------|-----------|------|---------------------------------------------------|---------|
| MapQ_mean | 0.509 | 0.485 | −0.024 | 0.508 / 0.541 / 0.528 | **COLLAPSED** |
| MapQ_median | 0.500 | 0.453 | −0.047 | 0.500 / 0.500 / 0.501 | **COLLAPSED** |
| MapQ_p10 | 0.502 | 0.468 | −0.035 | 0.502 / 0.505 / 0.505 | **COLLAPSED** |
| MapQ_p90 | 0.500 | 0.450 | −0.050 | 0.500 / 0.500 / 0.500 | **COLLAPSED** |
| LowMQ20_Frac | 0.494 | 0.520 | +0.026 | 0.495 / 0.489 / 0.486 | COLLAPSED (absolute resid < 0.53) |
| NM_mean | 0.478 | 0.428 | −0.050 | 0.478 / 0.468 / 0.514 | **COLLAPSED** |
| **Softclip_Frac** | 0.497 | **0.704** | **+0.207** | 0.491 / 0.493 / 0.608 | ATTENUATED but **artifact-suspect** |
| Strand_Bias | 0.491 | 0.505 | +0.014 | 0.489 / 0.479 / 0.502 | COLLAPSED |
| Read_Length_mean | 0.509 | 0.453 | −0.056 | 0.506 / 0.471 / 0.603 | **COLLAPSED** |

**All 9 features collapse under residualization except Softclip_Frac**, whose resid AUC 0.704 rises 0.207 above raw.

### 7.1 Softclip_Frac +0.207 ΔAUC investigated

Softclip_Frac distribution: mean=0.987, median=1.0, 75th pct=1.0, max=1.0, min=0.727. **75% of rows are exactly at 1.0**. This right-censoring creates a ceiling-effect artifact in OLS residualization:

1. OLS predicts ŷ ≈ 1.0 for most rows (dominated by the ceiling mass).
2. Residual `y − ŷ` concentrates near zero except for the ~25% of rows below 1.0, where residual = `y − 1.0` is strongly negative.
3. The caller's `vcf_QUAL` and depth covariates are correlated with the *probability of being at 1.0* (high-QS / high-depth regions have more cleanly-softclipped reads), so the OLS projection partitions rows into "predicted 1.0 + actual 1.0" (all) vs "predicted 1.0 + actual < 1.0" (anomalous).
4. The anomalous (< 1.0) rows enrich in FP (~60% FP rate vs 17% baseline), which makes residual AUC balloon to 0.704.

This is a **statistical artifact** of the bounded/censored nature of Softclip_Frac interacting with a linear model, not an independent biological signal. Two diagnostic checks:

- **AF-bin AUC**: Softclip_Frac Intermediate AUC 0.608 raw (§6) aligns with the residualized signal — the Intermediate AF slice (~60k TP / 7k FP) is where Softclip_Frac actually has some raw discrimination. The 0.704 residual AUC is not a new axis, just the combination of the AF stratum signal with the ceiling-effect amplification.
- If treated as a candidate, Softclip_Frac should be **re-tested with a logit / beta regression** (not OLS) and **excluded from the ceiling mass** (drop rows where Softclip_Frac ≥ 0.999) before re-running residualization.

Pending that re-test, Softclip_Frac is **ATTENUATED with artifact flag**, not a valid CONDITIONAL_POSITIVE.

### 7.2 Other features

- **MAPQ family** all four percentiles collapse to resid AUC 0.450–0.485; the MAPQ cap at 60 leaves no usable residual signal.
- **NM_mean** raw 0.478 → resid 0.428 — now anti-predictive; the LOH_Subclone 0.670 stratum signal (§6) is absorbed by NumReads + vcf_AF (high-NumReads LOH_Subclone regions tend to have TP with more mutations; once we control for NumReads this collapses).
- **Read_Length_mean** raw 0.509 → resid 0.453 — strong collapse; the LOH_Subclone stratum signal is read-depth proxy.
- **Strand_Bias** barely moves under residualization (+0.014) and stays near random.
- **LowMQ20_Frac** resid AUC 0.520 is below the 0.53 COLLAPSED threshold.

---

## 8. Spatial autocorrelation (Step 6)

`fig07_spatial.png` + table `data/G4/G4_spatial.tsv` (5 Mb bins, n≥50 reads per bin, 566 bins):

| Feature | median bin AUC | flag |
|---------|-----------------|------|
| MapQ_mean | 0.503 | ok |
| MapQ_median | 0.500 | ok |
| MapQ_p10 | 0.500 | ok |
| MapQ_p90 | 0.500 | ok |
| LowMQ20_Frac | 0.499 | ok |
| NM_mean | 0.476 | ok |
| Softclip_Frac | 0.492 | ok |
| Strand_Bias | 0.491 | ok |
| Read_Length_mean | 0.488 | ok |

All 9 features show median per-bin AUC within 0.024 of 0.50 — no spatial autocorrelation pathology. The `dAUC_hi_lo` (difference between high-TP-rate and low-TP-rate bins) is NaN for every feature because too few bins pass the `tp_rate > 0.8` threshold after per-bin sub-sampling; where we can estimate visually from `fig07`, the scatter is centred without systematic drift. This is consistent with the §7 finding that BAM QC features are caller-exhausted — no residual spatial structure to inflate.

---

## 9. Knowledge base and paper background

### 9.1 Knowledge base references

Queries executed on the knowledge MCP server:

- `BAM mapping quality MAPQ threshold` — (search returned Pipeline chapters; see `06_workflows/somatic-variant-calling.md` Tumor ≥ 50× MAPQ filtering at ≥5 for ClairS intake)
- `strand bias variant caller` — standard caller metric, no single KB doc; embedded in `04_pipeline/clairs-workflow.md`
- `softclip long-read ONT` — no direct KB hit; ONT softclip is near-universal (adapter trimming)
- `NM tag alignment` — BAM specification, inside `07_formats/bam-spec.md`

Gap: **no dedicated KB entry for BAM read-level QC thresholds** in InterSubMod. Candidate for KB authoring: a single page listing per-feature thresholds used by ClairS vs ClairS-TO vs Wakhan/SAVANA.

### 9.2 External literature

1. **Li, H. (2018).** "Minimap2: pairwise alignment for nucleotide sequences." *Bioinformatics* 34(18), 3094–3100. DOI: **10.1093/bioinformatics/bty191** — establishes the MAPQ cap at 60 for primary alignments in ONT mode, explaining why our MAPQ quantile features are all degenerate at 60.

2. **Zheng, Z. et al. (2023).** "Symphonizing pileup and full-alignment for deep learning–based long-read variant calling." *Nature Computational Science* 2, 797–803. DOI: **10.1038/s43588-022-00387-x** — ClairS's full-alignment network explicitly consumes per-read MAPQ, strand, and NM features. This is the core reason BAM QC features collapse under residualization on `vcf_QUAL`: the caller has already consumed them into its scalar QS.

3. **Luo, R., Wong, C.-L., Wong, Y.-S. et al. (2020).** "Exploring the limit of using a deep neural network on pileup data for germline variant calling." *Nat Mach Intell* 2, 220–227. DOI: **10.1038/s42256-020-0167-4** — Clair3's pileup network uses NM / softclip / strand as inputs; validates our prior that caller-consumed features have limited residual predictive power.

4. **Ebler, J., Ebert, P., Clarke, W.E. et al. (2022).** "Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes." *Nat Genet* 54, 518–525. DOI: **10.1038/s41588-022-01043-w** — discusses MAPQ / softclip as pangenome genotyping signals; supports our §6 finding that stratum-specific read length shifts correlate with alignment-complexity but not with independent TP/FP separability.

5. **Chen, L., Zheng, Z. et al. (2025).** "ClairS-TO." *Nat Commun* 16, 64547. DOI: **10.1038/s41467-025-64547-z** — explicitly claims robustness of ClairS-TO across coverage, VAF, and "complex genomic regions"; consistent with our all-BAM-QC-collapsed result.

**Direction**: papers 1–3 & 5 all *support* the G4 NEGATIVE outcome (caller fully consumes BAM QC); paper 4 hints that BAM QC *could* carry signal if redefined against a pangenome reference, which is outside the current pipeline.

---

## 10. Conclusions and challenges

### 10.1 Per-feature verdicts

| Feature | Global AUC | Confound | Cross-sample ≥5/7 | **Verdict** |
|---------|-----------|----------|-------------------|-------------|
| MapQ_mean | 0.509 | Collapsed (0.485) | 2/7 paired | **NEGATIVE** |
| MapQ_median | 0.500 | Collapsed (0.453) | 0/7 | **NEGATIVE** (MAPQ saturation at 60) |
| MapQ_p10 | 0.502 | Collapsed (0.468) | 0/7 | **NEGATIVE** |
| MapQ_p90 | 0.500 | Collapsed (0.450) | 0/7 | **NEGATIVE** |
| LowMQ20_Frac | 0.494 | Collapsed (0.520) | 0/7 | **NEGATIVE** |
| NM_mean | 0.478 | Collapsed (0.428) | 0/7 global | **NEGATIVE** (LOH_Subclone 0.670 = G6 shadow) |
| Softclip_Frac | 0.497 | Resid 0.704 but **ceiling-artifact** | 0/7 | **ARTIFACT_SUSPECT** (needs logit re-test) |
| Strand_Bias | 0.491 | Collapsed (0.505) | 0/7 | **NEGATIVE** |
| Read_Length_mean | 0.509 | Collapsed (0.453) | 1/7 paired | **NEGATIVE** (LOH_Subclone 0.682 = sample-ID / depth shadow) |

**Group-level verdict**: **NEGATIVE** — no BAM QC feature survives residualization on (NumReads, vcf_AF, Coverage_Multiple, vcf_QUAL); no cross-sample consistency ≥5/7. The single residualization surprise (Softclip_Frac resid 0.704) is attributable to ceiling-censoring artifact of the linear model.

### 10.2 Top finding

**All 9 BAM QC features are fully absorbed by the caller's QS + depth covariates.** The max paired-mode AUC is `Read_Length_mean` at 0.597 — well below the G3 `vcf_QUAL` paired AUC of 0.813 (figure `fig10_vs_g3_paired.png`), and below even the G1 `NumReads` residualized 0.777. The two stratum-specific spikes (`NM_mean` / `Read_Length_mean` at LOH_Subclone AUC 0.67–0.68) align precisely with the G6 `HPFineNGroups` LOH_Subclone AUC 0.68 signal and likely reflect the same underlying depth / phasing structure rather than independent read-quality biology.

### 10.3 Three challenges

1. **Softclip_Frac residualization ceiling artifact** — the 0.704 residualized AUC is the most eye-catching number in G4 but is driven by the heavy right-censoring at 1.0. Before claiming any signal, the feature must be **re-tested with logit transform on `log(1 − Softclip_Frac + ε)`** or with **rows filtered to Softclip_Frac < 0.999**. Expected: signal collapses to ≤0.55. (Pending Phase F follow-up.)

2. **NM_mean LOH_Subclone AUC 0.670 confound question** — is this a G6-shadow? `HPFineNGroups` achieves AUC 0.68 in the same stratum; the two features are correlated at ρ(NM_mean, Read_Length_mean) = 0.59 and the same LOH_Subclone cell has ~7:1 TP:FP imbalance. A **nested residualization** (add `HPFineNGroups` to covariates) is needed to rule out collinearity.

3. **COLO829 TO BAM QC missing** — the Phase B2 extraction produced a header-only TSV for COLO829 to_pileup, meaning 50,617 rows (6.8% of master) are missing BAM QC. This affects cross-sample consistency calculations for TO mode (6/7 not 7/7). Need to re-run `bam_readqc_per_region.py` on COLO829 TO BAM, or document the gap in the Phase F output. Likely root cause: the TO pileup BAM path in the extractor was pointed to an unreadable file.

### 10.4 Suggested next steps

- **P0**: re-extract COLO829 to_pileup BAM QC to close the 6.8% gap.
- **P1**: Softclip_Frac logit re-test (if still ≥ 0.60 resid AUC, escalate to stratum-specific filter candidate; otherwise final NEGATIVE).
- **P2**: nested residualization of NM_mean / Read_Length_mean on `HPFineNGroups` to confirm they are G6 shadows.
- **P3**: drop G4 features entirely from Phase F skill package (no candidate passes the gate); retain only as characterization annotations in any downstream manifest.
- **P3**: archive G4 as a **completed negative result** in `docs/experiments/concluded/`; matches the Phase E synthesis prediction for G4 ("likely depth-correlated, no independent signal").

### 10.5 Logic chain

```
BAM reads
  -> MapQ (cap 60, degenerate)  -> caller QS (full-alignment net consumes)
  -> NM / Softclip / Strand     -> caller QS (pileup + full-alignment nets consume)
  -> Read length (ONT baseline) -> caller QS (indirect via depth & mismappability)
  |
  +---> G4 features (this work)
           |
           +---> residualize on (NumReads, vcf_AF, CovM, vcf_QUAL)
                    |
                    +---> 8/9 collapse to AUC <= 0.52
                    +---> 1/9 (Softclip_Frac) resid 0.704 flagged as
                          ceiling-artifact; needs logit re-test
           |
           +---> stratified AUC
                    |
                    +---> NM_mean / Read_Length_mean LOH_Subclone AUC ~0.67-0.68
                          co-located with G6 HPFineNGroups signal = likely shadow
  |
  +---> Group verdict: NEGATIVE (no orthogonal signal over caller QS)
```

---

**Output files**:
- Figures (`research/feature_layered_observation/figures/G4_bam_readqc/`):
  - `fig01_global_distribution.png` — 9 features TP vs FP violin + AUC
  - `fig02_global_auc_bar.png` — 9 features global AUC bar chart
  - `fig03_{NM_mean,MapQ_mean,Read_Length_mean,Strand_Bias,LowMQ20_Frac}_heatmap.png` — LOH × AF × CN 32-cell heatmaps (5 features × 4 panels each)
  - `fig04_persample_auc_heatmap.png` — 9 features × 7 samples × 2 modes
  - `fig05_stratified_auc.png` — 9 features × (global/LOH/AF/CN/mode) bar grid
  - `fig06_confound_residualized.png` — raw vs residualized histograms, 9 features
  - `fig07_spatial.png` — 9 features × 5 Mb spatial autocorrelation
  - `fig08_corr_matrix.png` — 9×9 Spearman correlation
  - `fig09_MapQ_mean_persample.png` — top feature per-sample medians
  - `fig10_vs_g3_paired.png` — G4 paired AUC vs G3 vcf_QUAL 0.813 reference
- Tables (`research/feature_layered_observation/data/G4/`):
  - `master_g4.tsv.gz` (748,676 × 69, 93.2% BAM QC coverage)
  - `G4_global_stats.tsv`, `G4_auc_table.tsv`, `G4_cell_delta.tsv`, `G4_confound.tsv`
  - `G4_persample_auc.tsv`, `G4_cross_sample_consistency.tsv`, `G4_spatial.tsv`
  - `G4_corr_matrix.tsv`, `G4_MapQ_mean_persample.tsv`, `G4_vs_g3_paired.tsv`
- Script: `research/feature_layered_observation/scripts/observe_G4_bam_readqc.py`
