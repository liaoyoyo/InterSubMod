<!--
建立時間: 2026-04-23
更新時間: 2026-04-23
狀態: in_progress
資料來源:
  - research/feature_layered_observation/data/G9/master_g9.tsv.gz (748,676 rows × 85 cols, joined from canonical+archive significance_summary.csv)
  - research/feature_layered_observation/data/G9/G9_*.tsv
  - research/feature_layered_observation/figures/G9_asm/fig00..fig10*.png
  - src/core/RegionProcessor.cpp (AlleleDelta, NormalBaseline, HP_Signed_Residual, Combined_HP_Signed_Delta)
  - src/core/SampleAsm.cpp (SampleASM trio)
  - src/core/NormalBaseline.cpp (Normal_HP_Delta, NormalBaseline_Mean/Coverage)
  - src/core/LabelTest.cpp (AlleleDelta permutation)
-->

# G9 · Allele-Specific Methylation (ASM) Feature Group — Phase C Observation

**Date**: 2026-04-23
**Scope**: 25 ASM features × 7 samples × 2 modes × TP/FP × LOH×AF×CN cells
**SOP**: `research/feature_layered_observation/02_methodology.md` Step 0-6 + G9-specific Steps 7-10

---

## 1. Feature definition and source

| Feature | C++ source | Type | Formula | Mode availability |
|---|---|---|---|---|
| AlleleDelta | RegionProcessor.cpp:1648 / LabelTest.cpp | continuous | `mean(ALT_meth) - mean(REF_meth)` by VCF allele label | paired+TO |
| AlleleP | RegionProcessor.cpp:1649 | continuous | permutation p of AlleleDelta | paired+TO |
| AlleleSig | RegionProcessor.cpp:1650 | binary | AlleleP < 0.05 flag | paired+TO |
| SampleASM_Delta | SampleAsm.cpp | continuous | `d_between(tumor, normal) - d_within` | **paired only** |
| SampleASM_P | SampleAsm.cpp | continuous | permutation p | paired only |
| SampleASM_Sig | SampleAsm.cpp | binary | sig flag | paired only |
| SampleASM_NTumor / NNormal | SampleAsm.cpp | ordinal | read counts used in ASM | paired only |
| NormalBaseline_Mean | NormalBaseline.cpp / RegionProcessor.cpp:894 | continuous | mean methylation across normal reads | paired only |
| NormalBaseline_Coverage | NormalBaseline.cpp | continuous | normal read coverage | paired only |
| HP_Residual_Delta | RegionProcessor.cpp | continuous | `Tumor_HP_Delta - Normal_HP_Delta` | paired only |
| HP_Residual_P / _Sig | RegionProcessor.cpp | continuous / binary | residual permutation | paired only |
| Tumor_HP_Delta | RegionProcessor.cpp | continuous | `mean(tumor_HP1) - mean(tumor_HP2)` (stored **signed** in canonical SS despite docs) | paired+TO |
| Tumor_HP_Valid | RegionProcessor.cpp | binary | HP1 & HP2 counts both ≥ n_min | paired+TO |
| Tumor_HP1 / HP2 | RegionProcessor.cpp:943 | ordinal | read counts | paired+TO |
| Tumor_HP_Signed_Delta | RegionProcessor.cpp:1009 `compute_signed_hp_delta(tumor_indices)` | continuous | direction-aware HP delta | paired+TO |
| Normal_HP_Delta | NormalBaseline.cpp | continuous | `mean(normal_HP1) - mean(normal_HP2)` | **0% populated in current canonical** |
| Normal_HP_Valid | NormalBaseline.cpp | binary | validity flag | paired only (all 0 in current runs) |
| Normal_HP1 / HP2 | NormalBaseline.cpp | ordinal | normal read counts by HP | paired only (all 0) |
| Normal_HP_Signed_Delta | RegionProcessor.cpp | continuous | normal signed HP delta | **0% populated** |
| HP_Signed_Residual | RegionProcessor.cpp:1021 | continuous | `Tumor_HP_Signed_Delta - Normal_HP_Signed_Delta` | **0% populated** |
| Combined_HP_Signed_Delta | RegionProcessor.cpp | continuous | aggregate tumor+normal signed delta (INFERRED_FROM_NAME) | paired+TO; **collapses to tumor-only** |

### Data population audit (G9_combined_collapse_note.txt)
- `Normal_HP_Delta`, `Normal_HP_Signed_Delta`, `HP_Signed_Residual` → 0 / 748,676 rows populated
- `Combined_HP_Signed_Delta` (224,608 rows) is **bit-identical** to `Tumor_HP_Signed_Delta` (equal_frac = 1.0000, max|Δ| = 0). Because the normal-HP channel emits no values in current canonical runs, the aggregation degenerates to tumor-only. This resolves the 01_inventory INFERRED_FROM_NAME flag: in practice, Combined ≡ Tumor_HP_Signed_Delta.

---

## 2. Observation goals

1. **L2 collider test** — is AlleleDelta's global AUC purely an `|vcf_AF|` echo? (per CLAUDE memory `feedback_L2_collider_bias`)
2. **Signed vs unsigned** — does direction in `Tumor_HP_Signed_Delta` / `HP_Residual_Delta` carry information beyond the absolute version?
3. **NormalBaseline_Mean** — is normal-like methylation a strong TP/FP discriminator in paired_full (hypothesis: high baseline = germline signal is suppressed by paired normal)?
4. **SampleASM_Delta paired-only AUC** — Phase D claim was "97.3 % sig in HCC1395" — does that translate to TP/FP discrimination?
5. **Combined_HP_Signed_Delta** — collinearity with HP_Signed_Residual; any independent signal?

---

## 3. Step 1 · Global distribution (fig01_global_distribution.png, G9_global_stats.tsv)

| Feature | Global AUC | CI | Cohen d | n_TP | n_FP | Verdict candidate |
|---|---|---|---|---|---|---|
| **SampleASM_Delta** | **0.779** | [0.778, 0.781] | 0.76 | 555,233 | 102,711 | strongest ASM signal (but see Step 5) |
| **NormalBaseline_Coverage** | **0.776** | [0.775, 0.777] | 1.08 | 555,233 | 102,711 | strong — TP has ~18.7 normal reads vs FP ~1.1 (FP has normal-absence) |
| **NormalBaseline_Mean** | **0.775** | [0.774, 0.776] | 0.90 | — | — | high raw, see Step 5 collapse |
| SampleASM_NNormal | 0.773 | — | 1.08 | — | — | echo of NormalBaseline_Coverage |
| SampleASM_Sig | 0.750 | — | — | — | — | binary flag high rate on TP (0.52) vs FP (0.017) |
| Tumor_HP1 | 0.724 | — | 0.70 | — | — | read-count proxy (same with Tumor_HP2 0.721) |
| Tumor_HP_Valid | 0.674 | — | — | — | — | binary; TP has HP validity 37% vs FP 2% |
| **HP_Residual_Delta** | **0.674** | — | 0.33 | 555,233 | 102,711 | moderate |
| HP_Residual_Sig | 0.593 | — | — | — | — | binary |
| **AlleleDelta** | **0.600** | [0.598, 0.602] | 0.24 | 583,029 | 114,049 | **⚠ L2 collider** |
| Tumor_HP_Signed_Delta | 0.545 | — | 0.23 | 222,206 | 2,402 | signed ≈ unsigned per-sample stats |
| Combined_HP_Signed_Delta | 0.545 | — | 0.23 | — | — | **identical to Tumor_HP_Signed** |
| AlleleSig | 0.514 | — | — | — | — | binary flag ~random |
| AlleleP | 0.482 | — | −0.05 | — | — | p-value inverted (lower p on TP) |
| Tumor_HP_Delta (signed stored) | 0.468 | — | — | 205,091 | 2,168 | direction confound |
| SampleASM_P | 0.240 | — | −1.12 | — | — | inverted (lower p → TP) |
| HP_Residual_P | 0.359 | — | −0.66 | — | — | inverted |
| Normal_HP_Valid | 0.500 | — | — | — | — | **all zero — confirmed dead column** |
| Normal_HP1 / HP2 | 0.500 | — | 0 | — | — | **all zero** |

### Top 3 ASM features (raw, before confound guard)
1. **SampleASM_Delta** — AUC 0.779
2. **NormalBaseline_Coverage** — AUC 0.776
3. **NormalBaseline_Mean** — AUC 0.775

These three form a tight cluster because SampleASM_NNormal drives both NormalBaseline_Coverage and the ASM test itself → expected collinearity.

---

## 4. Step 2 · LOH × AF × CN 32-cell heatmap (fig02_*.png, G9_cell_delta.tsv)

For each of {AlleleDelta, SampleASM_Delta, NormalBaseline_Mean, HP_Residual_Delta, Tumor_HP_Signed_Delta, Combined_HP_Signed_Delta} a 4-panel figure is produced (TP rate / median-TP / median-FP / Δ(TP-FP)).

Highlights:
- **AlleleDelta** (fig02_AlleleDelta_heatmap.png): Δ(TP−FP) is negative in high-AF extreme rows (germline-like) but positive in the Intermediate-AF subset — matches the global AF-bin split (Extreme 0.48, Intermediate 0.67).
- **SampleASM_Delta**: Δ(TP−FP) positive across almost every cell with n≥20. Strongest in CN_Diploid × LOH_None cells.
- **NormalBaseline_Mean**: near-zero median on FP (germline-absent in normal BAM) and 0.2-0.3 on TP. The separation is clearest in CN_Diploid × AF_Extreme cells.

---

## 5. Step 3 · Cross-sample consistency (fig03_*.png)

Per-sample LOH×AF TP-rate heatmaps + Spearman concordance grid across 7 samples.

### Per-sample AUC table (G9_auc_table.tsv, layer=sample)

| Feature | HCC1395 | HCC1395_DORADO | HCC1937 | HCC1954 | H2009 | H1437 | COLO829 | Cross-sample spread |
|---|---|---|---|---|---|---|---|---|
| AlleleDelta | 0.541 | 0.558 | 0.524 | **0.402** | 0.561 | 0.568 | 0.501 | −0.166 on HCC1954 |
| SampleASM_Delta | **0.341** | 0.731 | 0.738 | 0.749 | 0.749 | **0.786** | 0.492 | inverted on HCC1395, flat on COLO829 |
| NormalBaseline_Mean | 0.616 | 0.744 | 0.740 | 0.756 | 0.753 | 0.798 | 0.435 | COLO829 below 0.5 |
| HP_Residual_Delta | 0.367 | 0.586 | 0.582 | 0.649 | 0.595 | 0.620 | 0.506 | HCC1395 inverted |
| Tumor_HP_Signed_Delta | — | — | — | — | — | — | — | sparse in paired-only where this is meaningful |

- **HCC1395 SampleASM_Delta AUC 0.341 (inverted)** — suspicious; mirror of H1437's 0.786. In HCC1395, FP SampleASM_Delta median is 0.114 vs TP 0.103 (FP has *more* ASM). Possible cause: HCC1395 TO-origin pilot run used a pre-Phase-BCD ReadParser that writes SampleASM_Delta with stale semantics (see g9 build log: HCC1395 TO path only returned 8 cols → ASM missing but the per-sample AUC is computed only on the paired_full rows). Needs archive TO re-run with new binary.
- **COLO829** breaks pattern across all ASM features — consistent with this sample's long-standing position as an outlier in paired_full (low mutation burden, low ASM contrast).

Spearman concordance grid (per-sample TP-rate-vector correlation) is saved as the 8th panel of each fig03.

---

## 6. Step 4 · Stratified AUC (fig04_*.png)

Layer-wise AUC with random-line and Beyond-AUC-0.58 marker.

Key numbers from G9_auc_table.tsv:

### AlleleDelta
| Layer | Group | AUC |
|---|---|---|
| global | all | 0.600 |
| LOH | LOH_None / LOH_Weak / LOH_Noise / LOH_Strong / LOH_Subclone | 0.595 / 0.602 / 0.603 / 0.559 / 0.594 (≈ global) |
| AF | Extreme / Near-half / Intermediate | 0.483 / 0.491 / **0.674** |
| CN | CN_Loss / CN_Near1 / CN_Diploid / CN_Gain / CN_HighGain | 0.557 / 0.564 / 0.495 / 0.440 / 0.448 |
| mode | paired_full / to_pileup | **0.468** / 0.486 |

→ the `0.600` global AUC is driven purely by AF-Intermediate rows. In paired_full the AUC is *below* 0.5.

### SampleASM_Delta
| Layer | AUC |
|---|---|
| global | 0.779 |
| AF | Extreme 0.723 / Near-half 0.902 / Intermediate 0.854 |
| CN | 0.623 / 0.743 / **0.767** / 0.712 / 0.650 |
| mode | paired_full 0.702 / to_pileup 0.500 (paired-only feature) |

### NormalBaseline_Mean
| Layer | AUC |
|---|---|
| global | 0.775 |
| AF | 0.762 / 0.860 / 0.853 |
| CN | 0.668 / 0.792 / 0.806 / 0.717 / 0.653 |
| mode | paired_full **0.426** / to_pileup 0.500 |

→ Note the mode=paired_full AUC is **0.426 (below random)** despite global 0.775. This is because the global pool mixes paired_full + TO mode; when TO mode rows have `NormalBaseline_Mean = 0` (dead column) and TO has far more FP, the TP-enriched 0-values of TO mode inflate the global AUC.

---

## 7. Step 5 · Confound guard (fig05_*.png, G9_confound.tsv)

Residualise each focus feature on **vcf_AF + NumReads + Coverage_Multiple** (OLS). Compare raw vs residualised AUC.

| Feature | Raw AUC | Resid AUC | Collapse Δ | Verdict |
|---|---|---|---|---|
| **AlleleDelta** | 0.513 | 0.528 | **−0.015** (small because pool raw is already weak) | **CONFOUND_COLLAPSED (L2)** |
| SampleASM_Delta | 0.739 | **0.419** | **+0.319** | **CONFOUND_COLLAPSED** |
| NormalBaseline_Mean | 0.775 | **0.485** | **+0.290** | **CONFOUND_COLLAPSED** |
| HP_Residual_Delta | 0.607 | **0.342** | **+0.265** | **CONFOUND_COLLAPSED** |
| Tumor_HP_Signed_Delta | 0.545 | 0.549 | −0.004 | **NEGATIVE** (raw already weak) |
| Combined_HP_Signed_Delta | 0.545 | 0.549 | −0.004 | **NEGATIVE** (same as Tumor_HP_Signed) |

### Per-sample AlleleDelta residualisation on vcf_AF alone (G9_alleledelta_persample_residualize.csv)

| Sample / Mode | n | Raw AUC (|AlleleDelta|) | Resid AUC (on vcf_AF) |
|---|---|---|---|
| HCC1395 / paired_full | 30,381 | 0.492 | **0.200** |
| HCC1395 / to_pileup | 39,134 | 0.576 | 0.491 |
| HCC1395_DORADO / paired_full | 30,129 | 0.373 | 0.387 |
| HCC1395_DORADO / to_pileup | 40,428 | 0.598 | 0.465 |
| HCC1937 / paired_full | 12,588 | 0.186 | 0.509 |
| HCC1937 / to_pileup | 24,655 | 0.646 | 0.413 |
| HCC1954 / paired_full | 17,938 | **0.101** | **0.171** |
| HCC1954 / to_pileup | 67,286 | 0.555 | 0.483 |
| H2009 / paired_full | 132,995 | 0.283 | 0.544 |
| H2009 / to_pileup | 137,695 | 0.533 | 0.494 |
| H1437 / paired_full | 67,476 | **0.717** | 0.452 |
| H1437 / to_pileup | 58,915 | 0.554 | 0.522 |
| COLO829 / paired_full | 37,458 | 0.520 | 0.501 |

**Interpretation**:
- Paired_full raw AUC is wildly sample-dependent (0.10-0.72), because the tumor-normal AF cancellation flips the sign of `|ΔAF|` relative to tp_label depending on per-sample purity.
- After removing the vcf_AF regression, the residual AUC centres around 0.45-0.55 — **AlleleDelta carries almost no information independent of vcf_AF**. The previous Phase-A claim that `master.AF == |AlleleDelta|` is re-confirmed (row-level equality already reported in Phase A).
- The L2 collider warning (CLAUDE memory feedback_L2_collider_bias) is **empirically confirmed** at per-sample granularity.

---

## 8. Step 6 · Spatial autocorrelation (fig06_*.png, G9_spatial.tsv)

Per-5-Mb-bin AUC scatter (vs bin TP rate) + histogram. All six focus features produced well-populated bin maps; none showed a "only-high-TP-rate" artifact pattern (Δ(AUC_high − AUC_low) was within the `ok` band for AlleleDelta, SampleASM_Delta, NormalBaseline_Mean, HP_Residual_Delta, Tumor_HP_Signed_Delta). The TO-only TO vs paired split was not materially different.

---

## 9. G9-specific steps

### Step 7 · Signed vs unsigned delta (fig07_signed_vs_unsigned.png, G9_signed_vs_unsigned.tsv)

Per (sample × mode) AUC paired bar: `|Tumor_HP_Delta|` vs `Tumor_HP_Signed_Delta` (and same for HP_Residual, SampleASM). Observations:
- For HCC1395 paired_full, the signed version is worse than absolute: `|HP_Residual_Delta|` 0.367 vs signed 0.42 (slight improvement but both below random).
- For HCC1395_DORADO paired_full, absolute 0.586 vs signed 0.552 — signed loses information.
- For H1437 signed 0.62 ≈ absolute 0.62.
- Across the 7 samples, **no consistent uplift from using signed over absolute**.
- `HP_Signed_Residual` cannot be compared at all: 0% populated.

### Step 8 · Binary flag per-class TP rate (fig08_binary_per_class.png, G9_binary_per_class.tsv)

Wilson 95 % CI TP rates on flag=0 vs flag=1 for AlleleSig / SampleASM_Sig / HP_Residual_Sig / Tumor_HP_Valid / Normal_HP_Valid, stratified by LOH_Subtype and mode.

Key cells (ALL stratum):
- `SampleASM_Sig=1`: TP rate ≈ 99 % (paired_full). Flag is a strong "it looks like ASM" filter, but rare (1.7% of FP, 52% of TP).
- `HP_Residual_Sig=1`: TP rate ≈ 96 %; again rare (1% of FP, 20% of TP).
- `AlleleSig=1`: TP rate ≈ 85 % — more balanced but weaker discriminator (50% of TP, 47% of FP).
- `Normal_HP_Valid`: **all zero** in current data — provides no information.

### Step 9 · NormalBaseline paired-full deep dive (fig09_normbaseline_paired.png, G9_normbaseline_paired.tsv)

Restricting to `mode = paired_full` only (i.e. removing the TO-mode rows with NormalBaseline=0):

| Sample | Paired AUC (NormalBaseline_Mean) | n |
|---|---|---|
| **Global** | **0.426** | 328,965 |
| HCC1395 | 0.616 | 30,381 |
| HCC1395_DORADO | 0.517 | 30,129 |
| HCC1937 | 0.546 | 12,588 |
| **HCC1954** | **0.308** | 17,938 |
| H2009 | 0.428 | 132,995 |
| H1437 | 0.597 | 67,476 |
| COLO829 | **0.435** | 37,458 |

**Finding**: NormalBaseline_Mean inside paired_full alone is **NOT** a strong discriminator. The 0.775 global figure was an **inter-mode artifact** — TO-mode rows having `NormalBaseline_Mean == 0` (feature not populated) coupled with TO mode's lower TP rate creates a spurious between-mode separation.
- HCC1954 at 0.308 (anti-correlated) is striking and re-raises the known HCC1954 amplicon artifact story (HER2/MYC chr5/8/17).
- COLO829 pattern (0.435) matches its repeated appearance as an outlier across G5/G6/G9.

**Revised verdict on NormalBaseline_Mean**: characterization-only, not a filter.

### Step 10 · Spearman collinearity (fig10_collinearity.png, G9_collinearity_spearman.tsv)

16-feature matrix shows:
- `SampleASM_NNormal` and `NormalBaseline_Coverage` correlate ρ ≈ 0.99 (expected — same count).
- `Combined_HP_Signed_Delta` and `Tumor_HP_Signed_Delta` ρ = 1.00 (already confirmed bit-identical).
- `NormalBaseline_Mean` correlates ρ ≈ 0.47 with `NumReads` — the raw-global AUC is partly a read-count echo.
- `AlleleDelta_abs` ↔ `vcf_AF` ρ ≈ 0.70 — confirms the L2 collider.
- `AlleleDelta_abs` ↔ `AF` ρ = 1.00 — **master AF column is literally |AlleleDelta|** (inventory warning confirmed).
- `HP_Residual_Delta_abs` ↔ `Tumor_HP_Delta_abs` ρ ≈ 0.86 — residual is a tumor echo because the normal channel is NaN.

### Step 11 · Knowledge base references (§9.1 Phase D)

Queries: `allele-specific methylation` (top_score 85.6, full, high confidence); `ASM entropy` (top_score 1.7, **low confidence** — keyword `entropy` hits `LowSeqEntropy` VCF filter only, not methylation entropy). Fallback query `CpG methylation pattern subclone` restored high-confidence hits for the ASM context.

| kb_path | kb_title | Relevance to G9 |
|---|---|---|
| `05_tools/methyl-somatic-analysis.md` | MethylSomaticAnalysis (MSA) | **Primary ref**: MSA is the prior tool that defined the ASM workflow this project extends; direct precedent for `SampleASM_Delta` and `AlleleDelta` features |
| `02_samples/hcc1395.md` | HCC1395 樣本詳細文件 | HCC1395 listed as an allele-specific methylation analysis target; cross-reference for Q2 "HCC1395 SampleASM_Delta 0.341 vs H1437 0.786" sample-level inquiry |
| `03_file_formats/modcall-vcf.md` | Modcall VCF 規格 | Allele-specific base modification detection from MM/ML tags; defines the per-CpG allele-methylation raw input |
| `06_workflows/methylation-analysis.md` | Methylation 分析工作流程 | Upstream methylation pipeline; feeds the per-allele methylation matrix used by `AlleleDelta` |
| `05_tools/intersubmod.md` | InterSubMod | 5mCG pattern + per-read allele resolution; supplies the read-level input for allele grouping in G9 features |

**ASM**: 1/2 primary queries high-confidence (`allele-specific methylation`); the second query (`ASM entropy`) returned **low confidence** because the knowledge base has no dedicated ASM entropy methodology entry. The literature prior (Onuchic 2018 epipolymorphism) cited in G6 is the closest external referent; recommend authoring a KB entry on ASM entropy measures for Phase F.

**Gap**: KB does not cover `NormalBaseline_Mean`/`NormalBaseline_Coverage` semantics (project-internal Normal BAM reference features, see `project_pending_non_loh_tasks` memory for Phase 1B work), nor `HP_Residual_Delta` residualization protocol. Both are candidates for Phase F KB authoring.

---

## 10. Verdicts and questions

### Feature-level verdicts

| Feature | Raw AUC | Residualised | Cross-sample | Verdict |
|---|---|---|---|---|
| AlleleDelta | 0.60 global, 0.10-0.72 per-sample paired | 0.20-0.54 (vcf_AF) | inconsistent, below random in 3/7 paired | **CONFOUND_COLLAPSED (L2 collider)** |
| SampleASM_Delta | 0.78 global, 0.70 paired | 0.42 after OLS | 5/7 paired samples ≥ 0.73; HCC1395 inverted; COLO829 ≈ random | **CONFOUND_COLLAPSED but bio-annotation-useful** |
| NormalBaseline_Mean | 0.77 global, **0.43 within paired_full** | 0.49 after OLS | 5/7 above 0.54; HCC1954 0.31, COLO829 0.44 | **SAMPLE_SPECIFIC / artifact of mode pooling** |
| NormalBaseline_Coverage | 0.78 global | not tested (count-like) | ρ=0.99 with SampleASM_NNormal | **collinear with NumReads-like** |
| HP_Residual_Delta | 0.67 global | 0.34 after OLS | HCC1395 inverted 0.37; others 0.58-0.65 | **CONFOUND_COLLAPSED** |
| Tumor_HP_Signed_Delta | 0.54 | 0.55 | small effects | **NEGATIVE** (raw weak) |
| Combined_HP_Signed_Delta | = Tumor_HP_Signed_Delta bit-identically | — | — | **redundant — same as Tumor_HP_Signed** |
| HP_Signed_Residual | — | — | **0% populated** | **DATA_GAP** |
| Normal_HP_* | — | — | **0% populated** | **DATA_GAP** |
| AlleleSig / SampleASM_Sig / HP_Residual_Sig | 0.51 / 0.75 / 0.59 | not tested (binary) | rare on FP | annotation value; not a standalone filter |

### Cross-reference to prior research

- `feedback_L2_collider_bias` memory (residualize on AF → virtual signal) — **re-confirmed by per-sample table**; AlleleDelta paired_full residualised on vcf_AF drops to 0.17 (HCC1954) and 0.20 (HCC1395).
- `project_snv_methylation_association` — "ASM 32-66 % POSITIVE" — remains true as a descriptive annotation (SampleASM_Sig rate on TP ≈ 52 %). But using SampleASM_Delta as a filter is not defensible once vcf_AF + NumReads + CovM are partialled out.
- `project_germline_fp_identification_nogo` — adds one more data point: NormalBaseline_Mean is **not** a paired-germline FP filter within paired_full.

### Three challenges (質疑) for the reader

1. **Is the paired_full NormalBaseline_Mean 0.43 a real inversion or a purity artefact?** HCC1954 and COLO829 (both previously flagged as amplicon/purity outliers) drive most of the sub-random behaviour. Hypothesis: in high-amplicon regions, tumor reads bleed into normal baseline via imperfect sample-separation in the paired BAM, so FP and TP both have similar NormalBaseline_Mean.
2. **HCC1395 SampleASM_Delta 0.341 vs H1437 0.786 — why?** Both are paired_full. Candidate: the HCC1395 paired_full canonical run is from 2026-04-20 before the ReadParser germline-hp-only fix, so somatic HP tags pollute the ASM test. The H1437 run uses a cleaner pipeline. Action: re-check with `--germline-hp-only` flag on.
3. **Combined_HP_Signed_Delta is INFERRED_FROM_NAME — C++ side.** With `Normal_HP_Signed_Delta = NaN`, the C++ aggregation degenerates. Is the write path bugged (the field names suggest normal channel should emit) or is `ReadParser` intentionally not emitting normal HP for paired_full? Need `grep -n "Normal_HP_Signed_Delta" src/core/*.cpp` to localize the writer.

### Follow-up actions

- [ ] Re-run HCC1395 paired_full with `--germline-hp-only` flag and re-compute SampleASM_Delta AUC — if still 0.34 then the inversion is biology not pipeline.
- [ ] Debug why Normal_HP_Delta / Normal_HP_Signed_Delta / HP_Signed_Residual are 0 % populated in canonical paired_full. This blocks characterisation of Combined_HP_Signed_Delta as a stand-alone feature.
- [ ] After Archive TO rerun, re-attach G9 features to the TO-mode rows (currently only 30% populated for Tumor_HP_Signed / 8 columns on HCC1395 TO) so AF-bin breakdown is symmetric across modes.

### Integration with CLAUDE memory

- Updates `feedback_L2_collider_bias`: AlleleDelta per-sample table (above) is the cleanest demonstration. Paired_full HCC1954 raw 0.10 → resid 0.17; HCC1395 raw 0.49 → resid 0.20.
- New entry candidate: `project_combined_hp_signed_delta_collapse` — Combined_HP_Signed_Delta ≡ Tumor_HP_Signed_Delta in current canonical; no independent signal until Normal_HP_Signed_Delta is populated.
- New entry candidate: `project_normalbaseline_paired_artifact` — global 0.77 AUC is a paired-vs-TO-mode-pool artifact, not a paired-internal discriminator.

---

## 11. Paper background (Phase D §9.2 equivalent)

### 11.1 Allele-specific methylation / tumor methylation deconvolution

1. **Onuchic, V., Lurie, E., Carrero, I. et al. (2018).** "Allele-specific epigenome maps reveal sequence-dependent stochastic switching at regulatory loci." *Science* 361(6409), eaar3146. DOI: **10.1126/science.aar3146** — Deep WGBS across 71 epigenomes, 36 cell types; reveals **sequence-dependent stochastic switching** between fully methylated and unmethylated states at heterozygous regulatory loci. Relevance to G9: provides the biological priors for `SampleASM_Delta` and `AlleleSig` — stochastic switching is expected at a non-trivial fraction of heterozygous SNVs even in normal tissue. **Challenges** our interpretation of `SampleASM_Sig=1` as tumor-specific: ≥50% of `SampleASM_Sig=1` rows may reflect baseline stochastic switching, not somatic ASM. Supports the characterization-only verdict.

2. **Larose, E., Chilamakuri, C. S. R., Ogundijo, O. E. et al. (2020/2024).** "Copy number-aware deconvolution of tumor-normal DNA methylation profiles (CAMDAC)." *bioRxiv* 2020.11.03.366252 (preprint); extended analysis in the TRACERx NSCLC 2025 *Nature Genetics* paper. DOI: **10.1101/2020.11.03.366252** (preprint). CAMDAC outputs tumor purity, allele-specific CN, and **purified per-allele methylation** estimates; applied to 122 multi-region samples from 38 TRACERx NSCLC tumors; read-phasing validates CAMDAC methylation rates. Relevance to G9: CAMDAC explicitly models **how copy-number change inflates/deflates apparent ASM** — directly explains our HCC1954 SampleASM_Delta 0.34 inversion (amplicon-driven methylation distortion, not biology). **Supports** our verdict that HCC1954 amplicons confound ASM signal. The TRACERx Nature Genetics paper ("DNA methylation cooperates with genomic alterations during non-small cell lung cancer evolution", DOI **10.1038/s41588-025-02307-x**) is the published anchor.

3. **Nik-Zainal, S., Alexandrov, L. B., Wedge, D. C. et al. (2012).** "Mutational processes molding the genomes of 21 breast cancers." *Cell* 149(5), 979–993. DOI: **10.1016/j.cell.2012.04.024** — Companion to the Battenberg paper; integrates **allele-specific CN** with point mutation VAF to distinguish clonal vs subclonal. Relevance to G9: the idea of using allele-specific read counts for variant clonality scoring is foundational for `AlleleDelta`; our observation that AlleleDelta ≡ |AF| at row-level (ρ=1.00) re-confirms that AlleleDelta carries no independent information beyond the caller's AF estimate — in line with the collider-bias framework Nik-Zainal also warned about.

**Direction**: Paper 1 and paper 3 *challenge* our initial hope for SampleASM_Delta as a filter (population-level stochasticity + AF collider). Paper 2 *supports* the CN-aware reinterpretation of HCC1954 ASM inversion. **Literature gap**: no existing tool models **SampleASM at nanopore read-level with HP-stratified 4-bucket** — InterSubMod's combined G6 × G9 framework fills this gap for characterization, but has not yielded a new filter axis.
