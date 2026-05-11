# G3 · VCF Caller Signals — Feature Layered Observation

**Date**: 2026-04-23
**Group**: G3 (Caller VCF fields)
**Features (paired_full + to_pileup)**: `vcf_QUAL`, `vcf_AF`, `vcf_DP`, `vcf_AD_{ref,alt}`, `vcf_H`
**Paired-only**: `vcf_NAF`, `vcf_NDP`, `vcf_NAD_{ref,alt}`, `vcf_GQ`
**TO-only**: `vcf_SB`, `vcf_PoN_1..4`, `vcf_Verdict_{Germline,Somatic,SubclonalSomatic}`
**Strand counts (paired)**: `vcf_FCU`, `vcf_RCU` (ref/alt × strand)
**Data**: `research/feature_layered_observation/data/merged_with_vcf.tsv.gz` (748,676 rows · PASS-only)
**Script**: `research/feature_layered_observation/scripts/g3_vcf_caller_analysis.py`
**Figures**: `research/feature_layered_observation/figures/G3_vcf_caller/` (~50 PNG)

> **Critical pre-condition**: master dataset is pre-filtered to `vcf_FILTER == PASS` for 100% of rows. All Verdict analyses in the PASS subset therefore see `Verdict_Germline = Verdict_Somatic = Verdict_SubclonalSomatic = 0` for every row (these flags only exist in the LowQual stratum removed upstream). The ClairS-TO Pilot (`project_clairs_to_verdict_pilot`) already confirmed this on a smaller subset; G3 now validates it on all 7 samples.

---

## 1. Feature definitions

- `vcf_QUAL` — ClairS/ClairS-TO aggregate QS (Phred-like, higher = higher confidence)
- `vcf_AF` — caller VAF (alt / total depth from the pileup network)
- `vcf_DP` — total depth at site in tumor BAM
- `vcf_AD_{ref,alt}` — per-allele depth at site
- `vcf_NAF / NDP / NAD_{ref,alt}` — paired-only: same quantities from **matched normal BAM**
- `vcf_GQ` — genotype quality (paired only)
- `vcf_SB` — strand bias score (TO only)
- `vcf_FCU / RCU` — forward/reverse strand C count (reference-strand); `FAU/FGU/FTU/RAU/RGU/RTU` are the other nucleotides. Most relevant as **strand imbalance diagnostic**.
- `vcf_PoN_1..4` — Panel-of-Normals hit flags (TO only, 4 tiers)
- `vcf_Verdict_{Germline,Somatic,SubclonalSomatic}` — ClairS-TO post-hoc verdict labels (TO only, LowQual-gated — see §2)

## 2. Observation goals

Establish whether the caller's own posterior signals (QUAL, GQ, AF, DP) separate TP from FP beyond the label they already drive (PASS/LowQual) — and whether `vcf_NAF` (normal-BAM alt frequency) provides an independent germline-contamination filter in paired mode. Re-confirm on all 7 samples the Pilot finding that Verdict flags live exclusively in LowQual rows.

## 3. Step 1 — global distribution (`fig01_global_*.png`)

Global AUC (|AUC|, ≥0.5 by reflection):

| feature | mode | raw AUC [95% CI] | direction | n non-nan |
|---|---|---|---|---|
| **vcf_QUAL** | paired_full | **0.813** [0.807, 0.818] | TP>FP | 328,965 |
| vcf_DP | paired_full | 0.786 [0.780, 0.791] | TP>FP | 328,965 |
| **vcf_NAF** | paired_full | 0.561 [0.552, 0.571] | **FP>TP** | 328,965 |
| vcf_FCU | paired_full | 0.577 | TP>FP | 328,965 |
| vcf_RCU | paired_full | 0.577 | TP>FP | 328,965 |
| vcf_AF | paired_full | 0.508 | TP>FP | 328,965 |
| vcf_AF | to_pileup | 0.582 | FP>TP | 419,711 |
| vcf_DP | to_pileup | 0.564 | TP>FP | 419,711 |
| vcf_QUAL | to_pileup | 0.531 | FP>TP | 419,711 |
| vcf_SB | to_pileup | 0.507 | — | 419,711 |

Paired `vcf_QUAL` and `vcf_DP` carry the strongest G3 signal. **TO mode has no QUAL/DP discrimination** — `vcf_QUAL` AUC 0.531 is effectively random, confirming the PASS stratum is a single low-separability mass for the caller.

## 4. Step 2 — LOH × AF × CN stratification (`fig02_*_strata.png`)

Per-feature 4-row heatmap (TP rate / mean TP / mean FP / Δ) for every numerical feature above. Highlights (paired):

- `vcf_QUAL` Δ(TP−FP) uniformly positive and largest in LOH_Strong × Near-half × CN_Normal (Δ ≈ +20 QS units). LOH_Subclone cells still show +Δ but n shrinks to <4.5k rows.
- `vcf_NAF` Δ is **negative in every cell** (FP have higher NAF) — consistent with germline-contamination hypothesis (see §6).
- `vcf_DP` Δ positive everywhere, largest in LOH_Weak/Strong.
- `vcf_FCU/RCU` Δ symmetric forward vs reverse → no strand-bias signal.

## 5. Step 3 — cross-sample consistency (`fig03_*_per_sample.png`)

Per-sample AUC for top signals (paired_full):

| feature | HCC1395 | HCC1395_D | HCC1937 | HCC1954 | H2009 | H1437 | COLO829 |
|---|---|---|---|---|---|---|---|
| **vcf_QUAL** | 0.898 | 0.865 | 0.825 | 0.755 | 0.943 | 0.951 | 0.788 |
| **vcf_NAF** | 0.502 | 0.763 | **0.990** | 0.948 | 0.907 | 0.784 | 0.525 |
| vcf_DP | 0.539 | 0.563 | 0.636 | 0.588 | 0.636 | 0.577 | 0.508 |
| vcf_AF | 0.925 | 0.522 | 0.906 | 0.848 | 0.727 | 0.750 | 0.532 |

`vcf_QUAL` ≥ 0.75 in **7/7** samples → strongest consistent caller signal in G3. `vcf_NAF` splits into two regimes: HCC1395 and COLO829 near-random (NAF uninformative), the other 5 samples AUC 0.78–0.99 — reflects how much the matched normal sees the somatic site (see §10.3). Cell-wise Spearman concordance of TP rates across samples is low (median ρ ≈ 0.16) because TP-rate cells are saturated near 1.0 in paired — variance for ordering is tiny.

## 6. Step 4 — stratified AUC (`fig04_*_stratified_auc.png`)

`vcf_QUAL` (paired_full) |AUC| by stratum:

| stratum | label | n | AUC |
|---|---|---|---|
| LOH | LOH_Strong | 22,973 | **0.926** |
| LOH | LOH_Weak | 46,233 | 0.866 |
| LOH | LOH_Subclone | 4,550 | 0.828 |
| LOH | LOH_Noise | 23,010 | 0.813 |
| LOH | None | 232,199 | 0.793 |
| AF | Near-half | 5,652 | **0.962** |
| AF | Intermediate | 42,385 | 0.905 |
| AF | Extreme | 280,928 | 0.796 |
| CN | Loss / Normal / Elevated / Gain / Low | all | 0.82–0.85 |

`vcf_QUAL` AUC > 0.79 in every stratum — **portable caller signal**, strongest in LOH_Strong and Near-half AF cells where caller uncertainty is naturally highest. No CN tier collapses the signal (all 0.82–0.85).

## 7. Step 5 — confound guard (`fig05_*_confound.png`)

OLS residualization on `(NumReads, Coverage_Multiple)` basis (paired_full):

| feature | raw AUC | resid AUC | ΔAUC | interpretation |
|---|---|---|---|---|
| vcf_QUAL | 0.813 | **0.797** | −0.016 | **Survives** — genuine caller signal beyond depth |
| vcf_DP | 0.786 | 0.745 | −0.041 | Partial survival — ~60% of signal is independent of NumReads |
| **vcf_NAF** | 0.561 | **0.729** | **+0.168** | **Signal amplified post-residualization** — NAF is masked by depth confound |
| vcf_AF | 0.508 | 0.524 | +0.016 | Flat — caller AF already too close to master AF |
| vcf_FCU/RCU | 0.577 | 0.53 | −0.05 | Weak, partially NumReads-proxy |

**Key finding**: `vcf_NAF` residualized AUC **0.729** (up from raw 0.561). The raw signal is suppressed because high-DP rows mix paired-TP (NAF=0) with paired-FP (NAF=small but depth-proportional). Once depth is removed, the direction `FP>TP` dominates and separation jumps. Per-AF-class validation (Step 5 sub-analysis): NAF AUC = Extreme 0.534 / Near-half **0.940** / Intermediate 0.690 — the Near-half stratum is where NAF most cleanly discriminates (these are exactly the cells where germline contamination is plausible).

## 8. Step 6 — spatial autocorrelation (`fig06_*_spatial.png`)

Per 5 Mb bin AUC:
- `vcf_QUAL` AUC-per-bin spread 0.75–0.95; no high-TP-bin inflation → no spatial artefact.
- `vcf_NAF` per-bin AUC highly bimodal (Extreme 0.5 vs Near-half 0.94) — matches §7 AF-class split, **not** a spatial artefact.
- `vcf_DP` per-bin AUC ≈ 0.70–0.80; stable across bins.
- `vcf_SB` (TO) AUC per-bin flat around 0.50 across every chromosome.

No G3 feature triggers the `project_spatial_autocorrelation_confound` flag.

## 9. vcf_FILTER and vcf_Verdict_* (critical negative)

- **PASS TP rate by mode**:
  - paired_full PASS: 98.9% (325,507 TP / 3,458 FP out of 328,965)
  - to_pileup PASS: 69.4% (291,324 TP / 128,387 FP out of 419,711)
  - TO PASS carries the entire FP burden of the dataset (128k / 131k total FP).

- **vcf_Verdict_* in PASS**: all three Verdict flags are **0 for 100% of PASS rows** in every sample. The Verdict annotations only fire on LowQual rows, which are absent from `merged_with_vcf.tsv.gz`. This reproduces the 2026-04-20 ClairS-TO Pilot finding on all 7 samples.

- **vcf_FILTER PASS vs LowQual AUC**: undefined in G3 (dataset is PASS-only). To evaluate LowQual-gated Verdict flags, a re-ingestion that preserves LowQual rows is required (see `docs/experiments/in_progress/2026/04/20260420_ClairS_TO_Verdict_Characterization_Pilot_01.md`).

- **TO QUAL threshold scan** (PASS subset): TP rate is flat at 0.694 for QUAL ≥ 4 up to ≥ 12; drops slightly to 0.671 only at QUAL ≥ 20 (but n falls from 420k to 183k). **No useful PASS-subset QUAL threshold** for TO mode FP reduction.

### Knowledge references (Phase D)

- **ClairS-TO paper** (Chen et al. 2025, *Nat. Commun.* 10.1038/s41467-025-64547-z) — documents Verdict categories (Germline, Somatic, SubclonalSomatic, Refcall) and explicitly states that Verdict tags primarily annotate rejected variants; our finding that PASS-subset Verdicts = 0 is consistent.
- **Clair3/ClairS** pileup-network depth dependence (Zheng et al. 2023, *Nat. Comput. Sci.* 10.1038/s43588-022-00387-x) — supports why `vcf_QUAL` partially co-moves with `NumReads` / `vcf_DP` (§7 residualization Δ −0.016 is small but the correlation is real).
- Internal `project_clairs_to_verdict_pilot` — prior warning that Verdict-on-PASS carries no signal.
- Internal `project_read_level_germline_fp` — hypothesised that matched-normal features could filter paired germline FPs; G3 §7 shows `vcf_NAF` resid AUC 0.729 is the best direct evidence for this, but note FP count in paired is already small (3,458 total).

## 10. Conclusions and challenges

### 10.1 Per-feature verdicts

| feature | mode | raw AUC | resid AUC | cross-sample ≥5/7? | Verdict |
|---|---|---|---|---|---|
| **vcf_QUAL** | paired_full | 0.813 | 0.797 | **7/7 ≥ 0.75** | **POSITIVE — dominant G3 signal** |
| vcf_DP | paired_full | 0.786 | 0.745 | 0/7 >0.65 | CONDITIONAL (partially NumReads-proxy) |
| **vcf_NAF** | paired_full | 0.561 | **0.729** | 5/7 ≥ 0.76 | **CONDITIONAL_POSITIVE — germline filter candidate** |
| vcf_AF | paired_full | 0.508 | 0.524 | mixed | NEGATIVE (already correlated with master AF) |
| vcf_AF | to_pileup | 0.582 | 0.579 | mixed | WEAK (matches Beyond-AUC ceiling) |
| vcf_QUAL | to_pileup | 0.531 | 0.528 | 7/7 ≤0.60 | NEGATIVE in PASS subset |
| vcf_DP | to_pileup | 0.564 | 0.602 | 7/7 ≤0.52 | NEGATIVE |
| vcf_FCU/RCU | paired_full | 0.577 | 0.53 | weak | WEAK (strand-symmetric) |
| vcf_SB | to_pileup | 0.507 | 0.501 | — | NEGATIVE |
| vcf_Verdict_* | to_pileup | — | — | all-zero in PASS | **STRUCTURALLY ABSENT** |

### 10.2 Key answers (requested)

- **`vcf_QUAL` paired global AUC = 0.813** (resid 0.797). Strongest per-region caller signal in G3.
- **`vcf_NAF` paired global AUC = 0.561 raw but 0.729 residualized** — germline contamination in normal BAM. NAF > 0.1 in 8.6% of FP vs 0.0% of TP.
- **`vcf_Verdict_Germline` in PASS subset = 0 everywhere.** FP rate cannot be computed from Verdict because no PASS row carries a non-zero Verdict tag. S1 ΔF1 = 0 confirmed on all 7 samples.
- **`vcf_FILTER` PASS-only dataset**: TP rate paired 98.9% / TO 69.4%. No LowQual data in master to compute contrast; Verdict / LowQual analysis must be rerun on the full (pre-PASS-filter) VCFs.

### 10.3 Three challenges

1. **vcf_QUAL vs NumReads collinearity** — QUAL resid AUC 0.797 vs NumReads resid AUC 0.777 (G1). G3 and G1 may be double-counting the same underlying depth-confidence axis. Action: regress `tp_label` on `(vcf_QUAL, NumReads)` jointly and report marginal contribution.
2. **vcf_NAF bimodality** — AUC ≈ 0.5 in HCC1395 & COLO829 but ≥ 0.76 in the other five samples. Correlates with sample purity (HCC1395 purity 99%; COLO829 known low-purity subclonal). Need purity-aware interpretation before proposing `vcf_NAF > 0.05` as a global filter.
3. **Verdict dark matter** — all Verdict_* are 0 in PASS; re-run the ingestion preserving LowQual rows to quantify Verdict_Germline precision. The 0420 Pilot gave ΔF1 = 0 because the flag never fires in PASS; the real test is whether Verdict-recovered LowQual rows rescue FN (not filter FP).

### 10.4 Next steps

- G10 integration: join `(vcf_QUAL, NumReads, vcf_NAF)` as the three-feature "caller baseline" against which every ISM methylation feature must prove marginal gain.
- Ingestion fix: write a new `merged_with_vcf_including_lowqual.tsv.gz` that retains all FILTER strata so Verdict flags and `vcf_FILTER` contrast can be analysed.
- Per-sample NAF thresholding: stratify `vcf_NAF > 0.05` filter by sample, report paired FP reduction per sample (HCC1937, HCC1954, H2009, H1437 likely high-yield; HCC1395 and COLO829 null).

---

**Output files**
- Figures: `research/feature_layered_observation/figures/G3_vcf_caller/fig01..06_*_{paired_full,to_pileup}_*.png` (~50 files — global / strata / per_sample / stratified_auc / confound / spatial for each feature×mode)
- Tables: produced by script in-memory; persisted tables follow G1/G2/G7 convention under `research/feature_layered_observation/data/G3/` (if materialised by script run). Rerun `python3 research/feature_layered_observation/scripts/g3_vcf_caller_analysis.py` to regenerate.
- Script: `research/feature_layered_observation/scripts/g3_vcf_caller_analysis.py`
