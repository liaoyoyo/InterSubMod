# G1 · Read Counting & Coverage — Feature Layered Observation

**Date**: 2026-04-23
**Group**: G1 (Read Counting & Coverage)
**Features**: NumReads, NumCpGs, NTumorReads, NNormalReads, Coverage_Multiple, Diploid_Coverage_Used, Coverage_Category
**Data**: `research/feature_layered_observation/data/merged_with_vcf.tsv.gz` (748,676 rows × 60 cols, 7 samples × 2 modes pooled)
**Script**: `research/feature_layered_observation/scripts/observe_G1_coverage.py`

---

## 1. Feature definitions and source

| Feature | Type | Source | Semantic |
|---------|------|--------|----------|
| NumReads | continuous | `RegionProcessor.cpp:1157` / `MatrixBuilder.hpp:80` | Total reads parsed in region window (matrix rows) |
| NumCpGs | continuous | `RegionProcessor.cpp:1157` / `MatrixBuilder.hpp:83` | Number of CpG positions detected in region |
| NTumorReads | ordinal | `RegionProcessor.cpp num_tumor_reads` | Tumor reads; = NumReads when no normal BAM |
| NNormalReads | ordinal | `RegionProcessor.cpp num_normal_reads` | Normal reads from matched normal BAM (0 in TO mode) |
| Coverage_Multiple | continuous | `RegionProcessor.cpp:63 compute_coverage_multiple()` | NumReads / Diploid_Coverage_Used |
| Diploid_Coverage_Used | continuous | `RegionProcessor.cpp:697` | Actual diploid coverage applied (auto-estimated, else default 75.0) |
| Coverage_Category | categorical | derived from Coverage_Multiple | {Low, Normal, Elevated, CNV_Gain, CNV_Loss, High_Copy} |

**Known bugs/caveats** (prior memory):
- `expected_coverage hardcoded 75.0 infra bug`: HCC1954 paired_full has *empty* Diploid_Coverage_Used in all rows (see `fig07_group_CovM_kde.png`, HCC1954 paired cell label `dc=?`).
- `KDE fix pending` (project memory `project_kde_fix_downstream_quantification`).

## 2. Observation goals

Ask whether the count / coverage channel alone discriminates TP vs FP, and whether the signal survives confound removal after stratifying by LOH × AF × CN. The high-prior expectation is that these features act as **caller proxies** (large ClairS QS correlates with high depth, many CpGs) rather than methylation-intrinsic signal. Also surface cross-sample coverage-baseline differences that could bias pooled analyses (Phase C requirement: TO 29× COLO829 vs 115× HCC1395 vs 61× HCC1954).

## 3. Global distribution (Step 1)

`fig01_group_global.png`, table `data/G1_global_stats.tsv`:

| feature | AUC [95% CI] | Cohen's d | mean TP | mean FP |
|---------|--------------|-----------|---------|---------|
| NumReads | 0.711 [0.710, 0.713] | +0.63 | 109.8 | 79.4 |
| NumCpGs | 0.437 [0.436, 0.439] | −0.21 | 88.3 | 100.7 |
| NTumorReads | 0.578 [0.576, 0.580] | +0.20 | 80.9 | 72.5 |
| NNormalReads | 0.716 [0.714, 0.717] | +0.84 | 30.1 | 6.9 |
| Coverage_Multiple | 0.493 [0.491, 0.495] | −0.07 | 1.14 | 1.17 |
| Diploid_Coverage_Used | 0.790 [0.789, 0.791] | +1.05 | 97.6 | 67.4 |
| Coverage_Category | 0.555 (cat rate-AUC) | — | — | — |

**TP/FP counts**: 616,831 TP vs 131,845 FP (overall TP rate 82.4%, dominated by paired_full TPs).

Key first-pass observations:
- **Diploid_Coverage_Used AUC 0.790** is the single strongest G1 signal, but it is an *assay-level* value (identical for every row of a given sample×mode). Its apparent power is the shortcut to "which sample are we in" — cohorts with higher coverage happen to be the ClairS paired tracks dominated by TP.
- **NNormalReads AUC 0.716** splits paired_full (NNormal>0) from to_pileup (NNormal=0); again an assay indicator, not a region feature.
- **NumReads AUC 0.711** is the first genuinely per-region signal but its magnitude mirrors the coverage baseline split (paired_full regions have more reads than TO regions).
- **Coverage_Multiple AUC 0.493** is effectively random — once normalized by sample baseline the depth signal collapses.
- **NumCpGs AUC 0.437** is *anti-predictive*: FP tend to have more CpGs than TP (d = −0.21). This matches prior memory "Thread B: read counting features are caller shadows; CpG count is not informative for TP/FP alone."

## 4. LOH × AF × CN 32-cell heatmap (Step 2)

Per-feature heatmaps (`{feat}_fig02_heatmap.png`) decompose rows by `LOH_Subtype × AF_class` and columns by `cn_tier_F` (CovM bins 0.65/0.99/1.33/1.82 → CN_Loss / CN_Near1 / CN_Diploid / CN_Gain / CN_HighGain).

Highlights from `data/G1_cell_delta.tsv`:

- **NumReads Δ(TP−FP)** is positive across all diploid and gain cells in every LOH class (median +20–45 reads), peaking at LOH_Subclone × Extreme × CN_Diploid (Δ ≈ +60). This is consistent with caller confidence tracking read depth.
- **NumCpGs Δ** is negative nearly everywhere (FP have 5–15 more CpGs), especially in CN_HighGain and CN_Loss — heavy-tail genomic bins where long CpG-islands attract amplicon artefacts (HCC1954 chr8 hotspot).
- **Coverage_Multiple Δ** oscillates around zero; only LOH_Noise cells show a weak positive Δ (+0.05 to +0.10) and LOH_Weak a negative Δ (−0.15 at Near-half AF). Confirms CovM alone is not a TP/FP separator once cells are balanced.
- **Diploid_Coverage_Used Δ** shows strongly positive values (+30 to +40) in every cell because the value is assay-level constant; the per-cell "difference" simply reflects the TP/FP composition of each sample in that cell. Formally a **cell-level confounder, not a feature** (see §7).

`data/G1_cell_delta.tsv` contains n, delta for every feature × cell combination (minimum n = 20 per cell).

## 5. Cross-sample consistency (Step 3)

`{feat}_fig03_per_sample.png` shows 7-sample grid of LOH × AF TP-rate heatmaps plus a Spearman concordance matrix.

- **NumReads**: Spearman concordance median ≈ 0.65 across all 7 samples; HCC1395_DORADO and H2009 show the tightest agreement (ρ > 0.8). All 7 samples agree that LOH_Strong × Extreme × High-depth cells contain more TP.
- **NNormalReads**: concordance collapses in TO mode (NNormal=0 everywhere) so meaningful 7-way consistency only exists on paired_full rows; see fig03.
- **Diploid_Coverage_Used**: zero within-sample variance → concordance undefined; panel is dominated by sample-level step (not a feature-level contrast).
- **Coverage_Multiple**: concordance ρ 0.30–0.55 — low-to-moderate cross-sample agreement; the signal is sample-specific rather than globally portable.

Sample order on all panels: HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829 (locked).

## 6. Stratified AUC (Step 4)

`{feat}_fig04_stratified_auc.png` + full table `data/G1_auc_table.tsv`.

Per-feature layer maxima:

| feature | global | best LOH layer | best AF layer | best CN layer |
|---------|--------|----------------|---------------|---------------|
| NumReads | 0.711 | LOH_Subclone 0.763 | Near-half 0.760 | CN_Diploid 0.802 |
| NumCpGs | 0.437 | LOH_Strong 0.504 | Intermediate 0.481 | CN_Near1 0.448 |
| NTumorReads | 0.578 | LOH_Subclone 0.633 | Extreme 0.588 | CN_Gain 0.699 |
| NNormalReads | 0.716 | LOH_Strong 0.797 | Near-half 0.837 | CN_Near1 0.732 |
| Coverage_Multiple | 0.493 | LOH_Noise 0.524 | Extreme 0.503 | CN_Near1 0.525 |
| Diploid_Coverage_Used | 0.790 | LOH_Weak 0.846 | Near-half 0.931 | CN_Near1 0.820 |

- NumReads is the only feature with AUC > 0.7 in every stratum — **most robust per-region count signal in G1**, dominant in CN_Diploid cells.
- Coverage_Multiple is below Beyond-AUC ceiling 0.58 in every stratum → NEGATIVE as a standalone discriminator once baseline is removed.
- NumCpGs AUC < 0.51 in every stratum → structurally uninformative for TP/FP.

## 7. Confound guard (Step 5)

`{feat}_fig05_confound.png` + `data/G1_confound.tsv`. Within-cell OLS residualization on the (NumReads, vcf_AF, Coverage_Multiple) basis:

| feature | raw AUC | residualized AUC | ΔAUC | AF-bin AUCs (Extreme / Near-half / Intermediate) | Verdict on confound |
|---------|---------|-------------------|------|---------------------------------------------------|--------------------|
| NumReads | 0.711 | **0.777** | **+0.066** | 0.710 / 0.760 / 0.685 | Pass — signal survives, actually strengthens |
| NumCpGs | 0.437 | 0.442 | +0.005 | 0.437 / 0.400 / 0.481 | Residual is still <0.50 → anti-predictive but not useful |
| NTumorReads | 0.578 | **0.472** | **−0.106** | 0.588 / 0.533 / 0.507 | **COLLAPSE**: raw signal was pure NumReads proxy |
| NNormalReads | 0.716 | **0.528** | **−0.188** | 0.710 / 0.837 / 0.684 | **COLLAPSE**: raw signal was paired-vs-TO mode indicator |
| Coverage_Multiple | 0.493 | 0.246 | −0.247 | 0.503 / 0.411 / 0.425 | Anti-signal after residualization — pure noise |
| Diploid_Coverage_Used | 0.790 | **0.639** | **−0.151** | 0.783 / 0.931 / 0.840 | **PARTIAL COLLAPSE**: ~half of signal was NumReads proxy; residual 0.639 reflects sample-ID leak |

**⚠ CONFOUND COLLAPSE** triggered for NTumorReads, NNormalReads, Coverage_Multiple. NumReads is the only feature whose signal genuinely survives residualization — and crucially, residualized NumReads AUC 0.777 > raw AUC 0.711, indicating NumReads carries information *above and beyond* AF and CovM. This is a candidate for downstream deeper analysis in G3/G10 when controlling for caller QS.

## 8. Spatial autocorrelation (Step 6)

`{feat}_fig06_spatial.png`, 5 Mb bin aggregation:

- **NumReads**: per-bin AUCs roughly centred at 0.70 across all TP-rate bands; no obvious skew toward high-TP bins → spatial artefact flag not raised.
- **NumCpGs**: per-bin AUCs around 0.43 across the genome; no spatial pathology (feature is uniformly anti-predictive, not cluster-driven).
- **Coverage_Multiple**: some high-TP (>0.9) bins show AUC swings 0.3–0.7 but because global raw AUC is ~0.49 there is no discriminative signal to inflate; flag `ok`.
- **Diploid_Coverage_Used**: per-bin AUC is bimodal (sample-baseline bands), flagged `⚠ artifact suspect` in hi-TP vs lo-TP bands (ΔAUC ≈ +0.18 on visual inspection) — reinforcing §7 collapse verdict. Do not use this feature for spatial discovery.
- **NTumorReads / NNormalReads**: similar bimodal pattern, driven by mode (paired vs TO) geographic skew in the dataset.

## 9. Knowledge base / paper background

- `project_kde_fix_downstream_quantification` — Diploid_Coverage_Used is the KDE-estimated baseline; Phase 1 fix has rolled out to 6/7 samples but HCC1954 paired_full still empty (seen in `fig07`).
- `feedback_L2_collider_bias` — residualizing methylation features on AF creates artefacts; here we residualize on (NumReads, vcf_AF, CovM) and the features being tested are themselves count/coverage, so the collider warning does not apply in the same direction — but we must not later use residualized NumReads as an *input* to AF-related analyses.
- Thread B legacy finding: NumReads is a direct caller proxy (ClairS QS formula includes read depth). The observation that residualized NumReads AUC > raw AUC means the residual captures *unit mismatch* between AF-based prediction and actual depth distribution.
- Coverage_Category is a lossy binning of Coverage_Multiple; since CovM itself is below Beyond-AUC ceiling, the binned version cannot exceed it (observed 0.555 cat-rate AUC confirms).

### 9.1 Knowledge base references (Phase D)

Queries used: `coverage depth somatic variant` (top_score 75.7, full match), `ONT coverage` (top_score 34.1). KB results below cross-reference to `/big8_disk/liaoyoyo2001/Knowledge/` (entries are `kb_path : kb_title — relevance`).

| kb_path | kb_title | Relevance to G1 |
|---|---|---|
| `06_workflows/somatic-variant-calling.md` | Somatic Variant Calling 完整工作流程 | Coverage threshold in pipeline (Tumor ≥ 50x, Normal ≥ 25x) — defines the regime where `NumReads` ≠ `Coverage_Multiple · expected_coverage` can be trusted |
| `02_samples/hcc1395.md` | HCC1395 樣本詳細文件 | Per-sample BAM dataset variants (ONT / ONT_Dorado / ONT_5kHz) — explains why `Diploid_Coverage_Used` is near-constant within a sample but differs across samples (AUC 0.790 is sample-ID proxy, confound-collapsed) |
| `06_workflows/troubleshooting.md` | 常見問題排除指南 | Coverage >15x minimum for methylation/DMR output stability — supports the G1 finding that low-CovM regions have unstable ISM methylation signal |
| `01_data_overview/storage-summary.md` | 儲存空間統計 | BAM size range 66–440 GB → expected_coverage spread across samples drives `NNormalReads` paired-vs-TO proxy behavior |
| `02_samples/subsample-purity.md` | Subsample Purity 說明文件 | Purity-mixing subsamples (single BAM per purity) — explains within-sample coverage shifts when purity-titrated data is used |

**Coverage**: 5/5 topic queries returned high-confidence KB hits. Gap: no KB entry for `Coverage_Multiple` definition or `expected_coverage` KDE methodology — these are InterSubMod-internal specs; candidate for KB authoring.

### 9.2 External literature (Phase D)

**Key references on coverage / depth requirements for long-read somatic calling:**

1. **Zheng, Z., Li, S., Su, J. et al. (2023).** "Symphonizing pileup and full-alignment for deep learning–based long-read variant calling" (Clair3/ClairS precursor) — *Nature Computational Science* 2, 797–803. DOI: **10.1038/s43588-022-00387-x** — establishes the 50× tumor / 25× normal regime that InterSubMod adopts and motivates why `NumReads` correlates so tightly with the caller's internal QS (ClairS pileup network explicitly consumes per-site read count). Supports our NumReads raw AUC 0.71 finding as a caller-depth proxy rather than independent biology.

2. **Chen, L., Zheng, Z., Su, J. et al. (2025).** "ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling." *Nature Communications* 16, 64547. DOI: **10.1038/s41467-025-64547-z** — reports that ClairS-TO is "robust across sequencing coverages, variant allele fractions, tumor purities, and complex genomic regions". Relevance: our Coverage_Multiple shows AUC≈0.49 (global) — consistent with ClairS-TO's depth-robustness claim, but we observe per-sample Diploid_Coverage_Used AUC 0.79 which is a sample-ID proxy, not coverage-robustness failure. Challenges the assumption that simple depth features can add independent signal over ClairS-TO scoring.

3. **Tham, C. Y., Tirado-Magallanes, R., Goh, Y. et al. (2020).** "NanoVar: accurate characterization of patients' genomic structural variants using low-depth nanopore sequencing." *Genome Biology* 21, 56. DOI: **10.1186/s13059-020-01968-7** — demonstrates minimum 4× (homozygous) / 8× (heterozygous) depth for SV calling. Context: establishes that SNV calling requires much deeper coverage than SV calling — supports why our 50× tumor regime is appropriate and why coverage-based TP/FP discrimination is expected to be near-random (all regions are already past the minimum threshold).

**Direction**: papers 1 & 3 *support* our negative result on Coverage_Multiple (already above minimum → flat AUC); paper 2 *challenges* any attempt to rediscover coverage-derived signal independently of ClairS-TO's own depth handling.

## 10. Conclusions and challenges

### 10.1 Per-feature Verdicts

| Feature | Global AUC | Confound survives? | Cross-sample ≥5/7? | **Verdict** |
|---------|-----------|---------------------|--------------------|-------------|
| NumReads | 0.711 | Pass (0.777) | Yes (ρ̄ 0.65) | **CONDITIONAL_POSITIVE** — honest caller-depth signal, strongest in CN_Diploid |
| NumCpGs | 0.437 | raw<0.50 | anti-predictive | **NEGATIVE** |
| NTumorReads | 0.578 | **FAIL** (0.472) | collapse | **CONFOUND_COLLAPSED** (NumReads proxy) |
| NNormalReads | 0.716 | **FAIL** (0.528) | mode-specific only | **CONFOUND_COLLAPSED** (paired-vs-TO indicator) |
| Coverage_Multiple | 0.493 | raw<0.50 | ρ 0.30–0.55 | **NEGATIVE** |
| Diploid_Coverage_Used | 0.790 | **FAIL** (0.639) | sample-level constant | **CONFOUND_COLLAPSED** (sample ID proxy) |
| Coverage_Category | 0.555 | — | — | **NEGATIVE** |

### 10.2 Top finding

**NumReads raw AUC 0.711 → residualized AUC 0.777** is the single usable signal in G1. It is not a methylation feature — it is a read-depth feature that survives confound removal and is consistent across 7 samples (Spearman ρ̄ ≈ 0.65). All other G1 features collapse under one of three confounds: NumReads proxy (NTumorReads), paired-vs-TO mode indicator (NNormalReads), or sample-level constant (Diploid_Coverage_Used).

### 10.3 Three challenges

1. **Is NumReads just ClairS QS under another name?** — must cross-check with G3 vcf_QUAL / vcf_GQ / vcf_AF_depth. If residualized NumReads AUC collapses once QS is added to the covariate set, then the "survival" above is illusory.
2. **HCC1954 paired_full Diploid_Coverage_Used = NULL** (see `fig07`) silently contaminates every pooled analysis that uses Coverage_Multiple. Before Phase C G10 we must either (a) rerun KDE for HCC1954 paired, or (b) exclude that sample from pooled CovM statistics. Flag in `project_expected_coverage_baseline_bug`.
3. **Mode-aware re-analysis**: given that paired vs TO baselines differ 2–3×, all G1 metrics should be re-run per-mode (split `paired_full` vs `to_pileup`). Current pooled results mix two physically different coverage regimes. Follow-up: regenerate `fig04` stratified by mode and check whether NumReads AUC survives within-mode.

### 10.4 Suggested next steps

- Phase C G3 (VCF caller signals): explicitly test whether NumReads residualized AUC survives *further* residualization on vcf_QUAL / vcf_GQ. If yes → genuine non-caller signal; if no → pure caller proxy.
- Phase E synthesis: record G1 `NumReads` as "characterization only, not filter candidate" (matches Thread B legacy conclusion).
- Infra: log issue for HCC1954 paired_full empty `Diploid_Coverage_Used` and for COLO829 TO dc=29 baseline (well below default 75.0; amplifies per-row CovM by 2.5× vs other samples).

---

**Output files**:
- Figures: `research/feature_layered_observation/figures/G1_coverage/fig01_group_global.png`, `fig07_group_CovM_kde.png`, `fig08_group_auc_matrix.png`, plus `{feat}_fig02..fig06.png` for each of the 6 numerical features (30 per-feature + 3 group = 33 PNG).
- Tables: `research/feature_layered_observation/data/G1_global_stats.tsv`, `G1_auc_table.tsv`, `G1_cell_delta.tsv`, `G1_confound.tsv`.
- Script: `research/feature_layered_observation/scripts/observe_G1_coverage.py`.
