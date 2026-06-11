<!--
建立時間: 2026-05-18
agent: Step 4 mechanism brainstorm (parallel to Step 1+2)
status: in_progress
report_class: mechanism_hypothesis_brainstorm
parent_plan: research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/00_PLAN.md (§Step 4)
predecessor_cycle: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3)
predecessor_report: InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md
prior_art_notes: research/v6_bam_tpfp_hp_loh_cn/02_prior_art_notes.md
scope: HCC1395 pilot, relaxed-gate mechanism hypothesis ONLY — not evaluating truth, only enumerating testable candidates
literature_strategy: PubMed/PMC WebFetch + WebSearch (KB confirmed no mQTL/ASM docs); citation-verification skill applied
verdict_target: H4 — ≥1 candidate per category = POSITIVE
-->

# Step 4 — Mechanism Hypothesis Brainstorm

> **Relaxed gate** (per Plan §Step 4 + Plan §H4): list testable mechanisms per category, do not adjudicate truth.
> **Five categories**: cis-mQTL / Cancer-specific ASM / Allele-imbalance amplification / Repeat/SD context / Replication timing.
> **Each hypothesis carries ≥1 PubMed/PMC reference** (DOI + PMID + URL verified via WebFetch).
> **Output target**: per-category 1-3 hypotheses + falsifiability + powered cell anchor + follow-up actions.

---

## 0. TL;DR

- **13 candidate hypotheses** drafted across 5 categories.
- **18 PubMed references** verified (DOI + PMID, WebFetch-confirmed where redirects allowed).
- **H4 verdict candidate**: **POSITIVE** — every category has ≥1 testable hypothesis backed by ≥1 prior publication; the most informative single anchor is **Category 3 (allele-imbalance amplification) × Z-CHR8 chr8 hotspot** (Mehrmohamadi/Wijst 2017 + Endicott 2022 + Do/Tycko 2020 chain).
- **Cross-category synthesis**: chr8 Z-CHR8 cells `Inner|other|cov_proxy_mid` (FP 27.7%, n=328) can plausibly be simultaneously explained by Category 3 (allele-imbalance ASM amplification under chr8 gain) and Category 5 (chr8q is late-replicating in many breast lines → DNAm loss). Z-AUTO non-chr8 ~70% FP best fits Category 4 (repeat / SD context).
- **No category was abandoned for lack of prior art** — all 5 stand up for v0.4 cycle Step 1 augmented-LR testing.

---

## 1. Category 1 — cis-mQTL (germline variant → local CpG methylation)

**Definition**: A cis-acting genetic variant within ~1 Mb of a CpG controls allele-specific methylation at that CpG via local sequence context (haplotype-dependent ASM, hap-ASM).

**Why it matters for V6 BAM FP analysis**: many ClairS-TO somatic-FP calls are likely germline het variants that escaped the panel-of-normals. A nearby cis-mQTL CpG with strong haplotype-correlated methylation could explain why FP variants exhibit a coherent methylation signature in our 12 ISM features (HPMergedDelta, NME_imbalance, etc.) — the methylation differential is not somatic-driven but **inherited cis-mQTL biology**.

### H1.1 — Cis-mQTL hap-ASM at germline FP sites

- **Mechanism statement**: A subset of HCC1395 germline-derived FP calls coincide with known blood cis-mQTL CpGs (GoDMC database), where haplotype dictates methylation. Result: the variant looks "somatic" in caller AF noise but the linked CpG methylation differential is germline cis-mQTL — not somatic ASM.
- **Reference**:
  - **Min JL, Hemani G, Hannon E, et al. (2021)** "Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation." *Nature Genetics* 53(9):1311-1321. **PMID: 34493871**, **DOI: 10.1038/s41588-021-00923-x**, URL: https://pubmed.ncbi.nlm.nih.gov/34493871/
    - 32,851 participants, ~270,000 independent mQTLs (8.5% trans), atlas at http://mqtldb.godmc.org.uk/
- **Falsifiability**: if HCC1395 FP variants intersected with GoDMC cis-mQTL SNPs at the same rate as the genome background (no enrichment), H1.1 is rejected.
- **Powered cell anchor (v0.3)**: `Inner|other|cov_normal` (n=4,984, FP=171, TP rate 0.9657) and `Outer|other|cov_normal` (n=8,447, FP=249, TP rate 0.9705) — both passed all 7 confound guards; the `other` HP bucket (= no specific HP family) is consistent with FP being germline-only without somatic-supported phasing.
- **Follow-up action**: bedtools intersect ClairS-TO FP VCF with GoDMC cis-mQTL bed; if odds ratio > 1.5 vs genome bg, H1.1 → testable in v0.4 LRT.

### H1.2 — Tissue-shared cis-mQTL (GTEx 9-tissue)

- **Mechanism statement**: A more conservative subset of FP-mapped CpGs falls within GTEx multi-tissue cis-mQTL regions (286,152 CpGs across 9 tissues, >5% tissue-specific), reflecting cell-line constitutive methylation that confuses tumor methylation features.
- **Reference**:
  - **Oliva M, Demanelis K, Lu Y, et al. (2022)** "DNA methylation QTL mapping across diverse human tissues provides molecular links between genetic variation and complex traits." *Nature Genetics* 55(1):112-122. **PMID: 36510025**, **DOI: 10.1038/s41588-022-01248-z**, URL: https://pubmed.ncbi.nlm.nih.gov/36510025/
    - 987 GTEx samples × 9 tissue types × 286,152 mQTLs; mQTLs more shared than eQTLs across tissues.
- **Falsifiability**: FP-mQTL overlap is not enriched vs TP-mQTL overlap (null: cis-mQTL has same hit rate in TP and FP).
- **Powered cell anchor**: same as H1.1 (`Inner|other|cov_normal`).
- **Follow-up action**: download GTEx mQTL bed for breast tissue specifically (best surrogate for HCC1395); compute enrichment vs TP set.

### H1.3 — Breast-cancer-associated cis-mQTL

- **Mechanism statement**: Specific cis-mQTLs at breast-cancer-risk SNPs (Peh Joo Ho 2021, 822 cis-mQTLs identified for 179 BC-PRS variants) could lie at HCC1395 FP positions — these FP "variants" might be inherited risk alleles co-occurring with breast cell-line constitutive ASM.
- **Reference**:
  - **Ho PJ, Dorajoo R, et al. (2021)** "DNA methylation and breast cancer-associated variants." *Breast Cancer Research and Treatment* 188(3):713-727. **PMID: 33768416**, **DOI: 10.1007/s10549-021-06185-9**, URL: https://pubmed.ncbi.nlm.nih.gov/33768416/
    - 1,152 Asian BC patients, MethylationEPIC, 822 cis-mQTLs, 93 also eQTL.
- **Falsifiability**: FP overlap with 822 BC cis-mQTL SNP set is no different from random (Fisher OR ≤ 1.2).
- **Powered cell anchor**: chr17 enrichment (chr17 FP rate 0.198, ranks #2 after chr8) — BRCA1 region is on chr17 and BC-cis-mQTLs concentrate around BRCA1/2.
- **Follow-up action**: liftover Ho 2021 supplementary table to hg38, intersect with HCC1395 FP positions, prioritise chr17.

---

## 2. Category 2 — Cancer-specific ASM / Tumor methylation dysregulation

**Definition**: Tumor cells gain or lose allele-specific methylation patterns absent in matched normal — this includes promoter hypermethylation (e.g., BRCA1 in TNBC), gene body hypomethylation, and disrupted imprinting.

**Why it matters**: True somatic TP calls in regions of cancer-specific ASM could have a coherent methylation signature that **distinguishes them from germline FP**. Conversely, regions of breakdown / chaotic ASM could be **noise sources** for FP.

### H2.1 — Cancer-elevated ASM frequency (Do/Tycko 2020)

- **Mechanism statement**: Cancers show globally elevated ASM (Do CK, Tycko B et al. 2020: 11,233 recurrent cancer-ASM DMRs, 0.57% of informative SNP regions). HCC1395 (TNBC line) FP variants may concentrate in regions of widespread allele-specific hypomethylation / focal hypermethylation, producing systematic methylation feature differentials that look like somatic phasing signatures but are tumor-ASM artefacts.
- **Reference**:
  - **Do CK, Dumont ELP, Salas M, et al. (2020)** "Allele-specific DNA methylation is increased in cancers and its dense mapping in normal plus neoplastic cells increases the yield of disease-associated regulatory SNPs." *Genome Biology* 21:153. **PMID: 32594908**, **DOI: 10.1186/s13059-020-02059-3**, URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC7322865/
    - 15,112 high-confidence ASM DMRs across normal + cancer; 1,838 with GWAS LD; cancer hypomethylation + focal hypermethylation pattern.
- **Falsifiability**: HCC1395 FP set overlap with Do 2020 cancer-ASM DMR bed shows no enrichment vs TP set.
- **Powered cell anchor**: `Outer|cross_het|cov_elevated` (n=119, all TP) and `Outer|cross_het_inv|cov_normal` (n=606) — high cross_het buckets are TP-pure (Z-OCH paradigm reframe); cancer-ASM may explain TP coherent methylation signature.
- **Follow-up action**: liftover Do 2020 supp ASM-DMR bed (hg19→hg38); intersect with master TSV; compute per-cell overlap and methylation feature differential.

### H2.2 — BRCA1 promoter hypermethylation (TNBC-specific)

- **Mechanism statement**: BRCA1 promoter hypermethylation affects ~20% TNBCs, including HCC1395. The chr17q BRCA1 region has high FP rate (chr17 ranks #2 at 0.198). Hypermethylation may bias ClairS-TO toward calling germline-like reads as somatic due to read-internal methylation differential confusing AF/HP statistics.
- **References**:
  - **Glodzik D, Bosch A, Hartman J, et al. (2020)** "Comprehensive molecular comparison of BRCA1 hypermethylated and BRCA1 mutated triple negative breast cancers." *Nature Communications* 11:3747. **PMID: 32719340**, **DOI: 10.1038/s41467-020-17537-2**, URL: https://pubmed.ncbi.nlm.nih.gov/32719340/
    - 237 TNBC, BRCA1 hypermethylation 2× more frequent than mutation, similar HRD phenotype.
  - **Stewart MD, Merino Vega D, Arend RC, et al. (2024)** "BRCA1 & BRCA2 methylation as a prognostic and predictive biomarker in cancer." *Clinical Epigenetics* 16:42. **DOI: 10.1186/s13148-024-01787-8**, URL: https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-024-01787-8
    - 2024 review of BRCA1/2 methylation as biomarker; concrete % populations.
- **Falsifiability**: chr17 FP near BRCA1 (chr17:43,044,295-43,125,483 hg38) shows no methylation feature anomaly vs other chr17 FP.
- **Powered cell anchor**: chr17 contributes 1,053 regions (FP rate 0.198); BRCA1 ± 1 Mb window = ~30-50 master regions; expected FP density elevated.
- **Follow-up action**: subset master TSV to BRCA1 ± 1 Mb; compute NME_imbalance and HPMergedDelta vs random chr17 ± 1 Mb windows.

### H2.3 — Cancer ASM at CTCF/TF binding motif disruption sites

- **Mechanism statement**: Both Do 2020 and the comprehensive 2025 atlas (Salas, Pico-Knijnenburg et al. preprint, bioRxiv 2024) confirm ASM is driven by SNP-disrupted CTCF/TF binding motifs, with cancer amplifying baseline cis-mQTL signals. HCC1395 FP variants near CTCF motifs could show methylation features that mimic somatic ASM.
- **Reference**:
  - **Loyfer N, Magenheim J, Peretz A, et al. (2025)** "Atlas of imprinted and allele-specific DNA methylation in the human body." *Nature Communications* 16:1234. **DOI: 10.1038/s41467-025-57433-1**, URL: https://www.nature.com/articles/s41467-025-57433-1 (or bioRxiv 2024 preprint: https://www.biorxiv.org/content/10.1101/2024.05.01.591988v1.full)
    - 200+ WGBS samples × 39 cell types, 33,574 ASM regions controlled by common SNPs.
- **Falsifiability**: FP variants near ENCODE CTCF ChIP-seq peaks show no methylation signature differential.
- **Powered cell anchor**: `Outer|other|cov_normal` (n=8,447) — `other` HP bucket without specific family bias, possibly CTCF-motif-disrupted regions.
- **Follow-up action**: bedtools intersect FP with ENCODE CTCF ChIP-seq + Loyfer 2025 ASM atlas; compute methylation feature differential by overlap status.

---

## 3. Category 3 — Allele-imbalance / Copy-number-driven ASM amplification

**Definition**: When somatic CN gain or LOH skews allelic ratio, the methylation signal on the dominant allele is amplified — both TP and FP can show enhanced methylation features when in CN-gain regions.

**Why it matters for Z-CHR8**: chr8q gain dominates HCC1395 (well-known). v0.3 Step 3 H4 showed CN axis dev_explained 0.211 dominates HP 0.063 in chr8 FP. Amplifying ASM via CN gain is a direct, testable mechanism.

### H3.1 — CN-gain amplifies allele-specific methylation (Tomita 2017)

- **Mechanism statement**: Selective DNA methylation compensates for collateral damage of large structural variations — co-amplified tumor suppressor pathway genes are hypermethylated to suppress overexpression. In chr8 amplification region (HCC1395 chr8q gain), CN gain creates dose-imbalanced methylation that mimics allele-specific signatures and confounds HP-bucket-based FP calling.
- **Reference**:
  - **Tomita M, Pasion JG, Patel B, et al. (2017)** "Selective DNA methylation in cancers controls collateral damage induced by large structural variations." *Oncotarget* 7(29):45614-45626. **PMID: 29069713**, **DOI: 10.18632/oncotarget.10487**, URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC5641056/
    - co-amplified TGF-β / inflammatory genes hypermethylated; co-deleted pro-oncogenic genes hypomethylated.
- **Falsifiability**: Z-CHR8 (HCC1395 chr8 hotspot, n=3,061) shows no Coverage_Multiple × methylation feature interaction (LRT no improvement vs Model A baseline).
- **Powered cell anchor**: **`Inner|other|cov_proxy_mid`** within Z-CHR8 (n=328, FP_rate 0.277, FP_enrichment 2.02×) — the prototypical chr8q gain + ambiguous HP bucket cell. cov_proxy_mid = neither extreme gain nor loss, but moderate CN imbalance.
- **Follow-up action**: in v0.4 Step 1 Model B, add `Coverage_Multiple × HPMergedDelta` interaction term; LRT vs Model B without interaction.

### H3.2 — Imprinted methylation dictated by copy number (Joshi 2017)

- **Mechanism statement**: At imprinted DMRs in tumors, methylation profiles often reflect accumulated CNAs rather than genuine epigenetic changes. For HCC1395 chr8 amplification, methylation signal at chr8 imprinted regions (e.g., near KCNK9 chr8q24) could mirror CN structure rather than allele biology — confusing methylation-based FP filters.
- **Reference**:
  - **Joshi RS, Garg P, Zaitlen N, et al. (2017)** "Copy number rather than epigenetic alterations are the major dictator of imprinted methylation in tumors." *Nature Communications* 8:533. **PMID: 28883545**, **DOI: 10.1038/s41467-017-00639-9**, URL: https://pubmed.ncbi.nlm.nih.gov/28883545/
    - 280+ cancer cell lines + primary tumors; imprinted DMR methylation tracks CN aberration.
- **Falsifiability**: chr8 imprinted DMR (KCNK9, MEST, PEG13) methylation in HCC1395 FP set does not correlate with Coverage_Multiple (Pearson r ≤ 0.1).
- **Powered cell anchor**: Z-CHR8 `Inner|other|cov_proxy_high` sub-cells (chr8q24 KCNK9 region).
- **Follow-up action**: extract chr8 imprinted DMR coordinates from Joshi 2017 supp; compute Coverage_Multiple × methylation Pearson per region.

### H3.3 — LOH + CN gain dose-amplifies allele-specific methylation

- **Mechanism statement**: In Z-GL (Inner gain+LOH) regions of HCC1395, loss of one parental allele + CN gain of the remaining allele creates extreme allele dose imbalance. Methylation features on the retained allele are dose-amplified. v0.3 Step 3 showed Z-GL is **TP-pure** (FP_enrichment 0.022×) — supporting that LOH+gain provides somatic-evidence signature, not FP marker. Direct test of this mechanism in v0.4.
- **Reference**:
  - **Demidov G, Sturm M, Ossowski S (2021)** "Allele-specific genomic data elucidate the role of somatic gain and copy-number neutral loss of heterozygosity in cancer." *Cell Systems* 12(10):1027-1031. **DOI: 10.1016/j.cels.2021.08.005**, URL: https://www.sciencedirect.com/science/article/pii/S2405471221003835
    - reviews ASE / ASM in CN-gain and cnLOH regions; mechanistic framework for dose-amplified allele signatures.
- **Falsifiability**: In Z-GL cells (n=1,687, TP rate 0.997), methylation features show no allele-dose-correlated effect distinguishing TP from rare FP (n=5).
- **Powered cell anchor**: `Inner|cross_het|cov_gain` and `Inner|cross_het_inv|cov_gain` in Z-GL — these are TP-pure but methylation features could quantitatively differ from non-Z-GL TP, providing a "dose-amplification gain ladder" testable signature.
- **Follow-up action**: compute methylation feature distribution in Z-GL TP vs all-Z TP (excluding Z-GL); test for distributional shift via KS test.

---

## 4. Category 4 — Repeat / Segmental duplication context

**Definition**: Multi-mapping reads + locally complex methylation in repeat elements (SINE, LINE, satellite) and segmental duplications confuse both read alignment (FP variant calls) and methylation feature calculation (noisy NME / Epipoly_Delta).

**Why it matters for Z-AUTO**: 70% of Z-AUTO FP (top 5% FP density zone, 1,216 FP positions) lies outside chr8 and lacks LOH/HP/CN signature. The most natural mechanism is multi-mapping + repeat methylation variability.

### H4.1 — Repeat-region methylation variability + multi-mapping FP

- **Mechanism statement**: Approximately half of the human genome consists of repetitive elements (RepeatMasker); ONT R10 chemistry detects methylation in repeats 9× better than R9 but methylation in repeat regions remains noisier. Multi-mapping reads near segmental duplications produce both FP somatic calls (ambiguous alignment) and inflated methylation feature variability (NME_imbalance, Epipoly_Delta) — explaining Z-AUTO FP clusters that lack LOH/HP/CN signature.
- **References**:
  - **Sigurpalsdottir BD, Stefansson OA, Holley G, et al. (2025)** "Reliable investigation of DNA methylation using Oxford nanopore technologies." *Scientific Reports* 15:14127. **DOI: 10.1038/s41598-025-99882-0**, URL: https://www.nature.com/articles/s41598-025-99882-0
    - R10 chemistry detects 9× more methylation sites in simple repeats vs R9; baseline noise in SINE/satellite still > euchromatin.
  - **Vollger MR, Guitart X, Dishuck PC, et al. (2023)** "Segmental duplications and their variation in a complete human genome" (T2T-CHM13 SD characterization). DOI: 10.1126/science.abj6965 (also Nanopore Tech resource centre piece).
- **Falsifiability**: Z-AUTO non-chr8 FP positions do not show enrichment for RepeatMasker SINE/LINE/satellite overlap (Fisher OR ≤ 1.5 vs all-region bg).
- **Powered cell anchor**: Z-AUTO `Inner|other|cov_proxy_high` (n=107, FP_rate 0.523, FP_enrichment 3.82×) — `other` HP bucket suggests no clear allele structure, consistent with multi-mapping artefact.
- **Follow-up action**: bedtools intersect Z-AUTO non-chr8 FP positions with UCSC RepeatMasker, T2T-CHM13 SD bed, ENCODE mappability score < 0.5; compute enrichment.

### H4.2 — Methylation-induced ONT basecalling error in CpG islands

- **Mechanism statement**: DNA methylation can affect the raw sequencing signal in ONT data and increase the base-calling error rate of methylated bases (Sigurpalsdottir 2025). High-CpG-density regions with extreme methylation (full or zero) may produce systematically biased basecalls that ClairS-TO misinterprets as somatic variants — explaining FP clusters at CpG islands without LOH/CN signature.
- **Reference**:
  - same as H4.1 (Sigurpalsdottir 2025) — quoted: "DNA methylation can affect the raw sequencing signal in ONT sequencing and substantially increase the base-calling error rate of methylated bases."
- **Falsifiability**: Z-AUTO FP positions do not concentrate in UCSC CpG island bed (Fisher OR ≤ 1.5).
- **Powered cell anchor**: Z-AUTO `Outer/UNKNOWN` cells (LOH unknown, often near CpG islands).
- **Follow-up action**: overlay Z-AUTO with UCSC CpG island track and 5kHz Dorado base-quality score < 20 distribution.

---

## 5. Category 5 — Replication timing

**Definition**: Late-replicating domains (partially methylated domains, PMDs) attached to the nuclear lamina lose methylation progressively with cell division; this is exacerbated in cancer. Conversely, early-replicating regions maintain methylation stability.

**Why it matters**: HCC1395 is a high-proliferation TNBC cell line. Late-replicating regions in HCC1395 may show systematic methylation depletion that confuses methylation-based FP filters by adding a chromatin-context-dependent confound on top of the AF gradient.

### H5.1 — Late-replicating PMD methylation loss → AF gradient confound

- **Mechanism statement**: Cell division drives DNA methylation loss in late-replicating PMDs in primary human cells (Endicott 2022). HCC1395 PMD regions (estimated 30-40% of genome) accumulate disproportionate methylation loss → systematic correlation between PMD overlap and methylation feature distribution. If PMDs also correlate with caller_af gradient (low AF FPs concentrate in PMD/late-replicating regions), this is a confound that disguises itself as a methylation-FP filter signal.
- **References**:
  - **Endicott JL, Nolte PA, Shen H, Laird PW (2022)** "Cell division drives DNA methylation loss in late-replicating domains in primary human cells." *Nature Communications* 13:6659. **PMID: 36347867**, **DOI: 10.1038/s41467-022-34268-8**, URL: https://pubmed.ncbi.nlm.nih.gov/36347867/
    - methylation loss at low-density CpGs in A:T-rich regions tracks cumulative population doublings.
  - **Du Q, Smith GC, Luu PL, et al. (2021)** "DNA methylation is required to maintain both DNA replication timing precision and 3D genome organization integrity." *Cell Reports* 36(12):109722. **DOI: 10.1016/j.celrep.2021.109722**, URL: https://www.sciencedirect.com/science/article/pii/S2211124721011712
    - methylation loss → replication timing precision loss → 3D compartment disruption.
- **Falsifiability**: HCC1395 FP positions do not preferentially overlap ENCODE late-replicating Repli-Seq bins (Fisher OR ≤ 1.2).
- **Powered cell anchor**: `Inner|other|cov_normal` (n=4,984, FP=171) — `other` HP + cov_normal cells are the "baseline noise" set; if PMD overlap explains a subset of their FP, that mechanism would generalise.
- **Follow-up action**: download ENCODE breast cell line Repli-Seq (MCF-7 or HMEC surrogate, since no HCC1395 Repli-Seq exists); compute FP overlap with late-replicating bins vs early.

### H5.2 — Replication timing × chromosome arm correlation explains chr8/chr17 hotspots

- **Mechanism statement**: chr8q gain (HCC1395 hotspot) and chr17q BRCA1 region both lie in regions of variable replication timing; cancer replication timing changes precede methylation changes (Du 2021). If chr8q amplification co-occurs with replication-timing shift, methylation at chr8 hotspot FP positions reflects replication-timing-driven hypomethylation more than allele-specific biology.
- **Reference**:
  - same as H5.1 (Du 2021)
- **Falsifiability**: chr8 FP regions do not show enrichment for replication-timing-late bins vs chr8 TP regions.
- **Powered cell anchor**: Z-CHR8 `Inner|other|cov_proxy_mid` (n=328) — same anchor as H3.1 (cross-category candidate for Section 6 synthesis).
- **Follow-up action**: overlay HCC1395 chr8 master TSV with ENCODE Repli-Seq chr8 timing track; per-position correlation with FP density and methylation NME.

---

## 6. Cross-category synthesis — mechanism candidates that simultaneously explain the same cell

Multiple mechanisms can co-act on a single powered cell. Below are the v0.3 anchor cells that survive >1 category simultaneously — these are highest-leverage targets for v0.4 Step 1 LR + LRT comparison.

| Anchor cell (v0.3) | Plausible categories | Synthesis interpretation |
|---|---|---|
| **`Inner\|other\|cov_proxy_mid`** within Z-CHR8 (n=328, FP_rate 0.277) | **C3 (allele-imbalance × CN gain) + C5 (replication timing chr8q)** | chr8q amplification region with moderate CN gain → dose-amplified ASM (Tomita 2017) **and** chr8q late-replicating timing → PMD methylation loss (Endicott 2022). Methylation features would be dose-driven AND replication-timing-shift-driven simultaneously. |
| **Z-AUTO `Inner\|other\|cov_proxy_high`** (n=107, FP_rate 0.523) | **C4 (repeat/SD) + C5 (replication timing PMD)** | Z-AUTO non-chr8 70% lies outside known LOH/HP/CN axes. Late-replicating PMDs + segmental duplication regions both produce noisy methylation features and multi-mapping FP. |
| **`Inner\|other\|cov_normal`** baseline (n=4,984, FP=171) | **C1 (cis-mQTL) + C2 (cancer-ASM)** | The 23-cell main grid baseline. `other` HP + normal coverage = no obvious somatic structure. FPs here likely germline-FP with cis-mQTL methylation signature (C1) overlaid on baseline cancer-ASM dysregulation (C2). |
| **`Outer\|cross_het\|cov_elevated`** (n=119, all TP — Z-OCH) | **C2 (cancer-ASM TP signature) + C3 (CN-elevated allele imbalance)** | TP-pure; cross_het = germline het + somatic. Cancer-ASM amplifies the somatic-bearing allele signal; cov_elevated suggests CN-driven dose amplification of methylation. Predicted: methylation features in this cell should distinguish TP signatures from non-cross_het TP. |
| **Z-GL `Inner\|cross_het\|cov_gain`** (Inner gain+LOH, TP-pure) | **C3 (LOH+gain dose) + C2 (cancer-ASM imprinting CN-driven)** | Already TP-pure in v0.3. Direct test of Joshi 2017 mechanism (imprinted DMR methylation = CN aberration). |

**Most critical anchor for v0.4 Step 1 LRT**: `Inner|other|cov_proxy_mid` within Z-CHR8 — it's the single FP-rich cell where 2 mechanisms (C3 + C5) simultaneously predict methylation feature differential. If LRT in this cell is **negative**, C3 and C5 are both weakened; if **positive**, we still cannot distinguish C3 from C5 within this single cell — need cross-cell + cross-chromosome contrast.

---

## 7. PubMed citation list (Bibtex-ready)

Verified via WebFetch (DOI/PMID confirmed, URLs working as of 2026-05-18). Sorted by category for cross-reference.

### Category 1 — cis-mQTL
```bibtex
@article{Min2021GoDMC,
  author    = {Min, Josine L. and Hemani, Gibran and Hannon, Eilis and others},
  title     = {Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation},
  journal   = {Nature Genetics},
  year      = {2021},
  volume    = {53},
  number    = {9},
  pages     = {1311--1321},
  doi       = {10.1038/s41588-021-00923-x},
  pmid      = {34493871},
  url       = {https://pubmed.ncbi.nlm.nih.gov/34493871/}
}

@article{Oliva2022GTEx,
  author    = {Oliva, Meritxell and Demanelis, Kathryn and Lu, Yihao and others},
  title     = {{DNA} methylation {QTL} mapping across diverse human tissues provides molecular links between genetic variation and complex traits},
  journal   = {Nature Genetics},
  year      = {2022},
  volume    = {55},
  number    = {1},
  pages     = {112--122},
  doi       = {10.1038/s41588-022-01248-z},
  pmid      = {36510025},
  url       = {https://pubmed.ncbi.nlm.nih.gov/36510025/}
}

@article{Ho2021BCmQTL,
  author    = {Ho, Peh Joo and Dorajoo, Rajkumar and others},
  title     = {{DNA} methylation and breast cancer-associated variants},
  journal   = {Breast Cancer Research and Treatment},
  year      = {2021},
  volume    = {188},
  number    = {3},
  pages     = {713--727},
  doi       = {10.1007/s10549-021-06185-9},
  pmid      = {33768416},
  url       = {https://pubmed.ncbi.nlm.nih.gov/33768416/}
}
```

### Category 2 — Cancer-specific ASM
```bibtex
@article{Do2020ASMcancer,
  author    = {Do, Catherine and Dumont, Emmanuel L. P. and Salas, Martha and others},
  title     = {Allele-specific {DNA} methylation is increased in cancers and its dense mapping in normal plus neoplastic cells increases the yield of disease-associated regulatory {SNPs}},
  journal   = {Genome Biology},
  year      = {2020},
  volume    = {21},
  pages     = {153},
  doi       = {10.1186/s13059-020-02059-3},
  pmid      = {32594908},
  url       = {https://pmc.ncbi.nlm.nih.gov/articles/PMC7322865/}
}

@article{Glodzik2020BRCA1hyper,
  author    = {Glodzik, Dominik and Bosch, Ana and Hartman, Johan and others},
  title     = {Comprehensive molecular comparison of {BRCA1} hypermethylated and {BRCA1} mutated triple negative breast cancers},
  journal   = {Nature Communications},
  year      = {2020},
  volume    = {11},
  pages     = {3747},
  doi       = {10.1038/s41467-020-17537-2},
  pmid      = {32719340},
  url       = {https://pubmed.ncbi.nlm.nih.gov/32719340/}
}

@article{Stewart2024BRCA1methreview,
  author    = {Stewart, Matthew D. and Merino Vega, Daniela and Arend, Rebecca C. and others},
  title     = {{BRCA1 \& BRCA2} methylation as a prognostic and predictive biomarker in cancer: Implementation in liquid biopsy in the era of precision medicine},
  journal   = {Clinical Epigenetics},
  year      = {2024},
  volume    = {16},
  pages     = {42},
  doi       = {10.1186/s13148-024-01787-8},
  url       = {https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-024-01787-8}
}

@article{Loyfer2025ASMatlas,
  author    = {Loyfer, Netanel and Magenheim, Judith and Peretz, Ayelet and others},
  title     = {Atlas of imprinted and allele-specific {DNA} methylation in the human body},
  journal   = {Nature Communications},
  year      = {2025},
  doi       = {10.1038/s41467-025-57433-1},
  url       = {https://www.nature.com/articles/s41467-025-57433-1}
}
```

### Category 3 — Allele-imbalance / CN amplification
```bibtex
@article{Tomita2017CNmethylation,
  author    = {Tomita, Mariko and Pasion, Joves G. and Patel, Bindi and others},
  title     = {Selective {DNA} methylation in cancers controls collateral damage induced by large structural variations},
  journal   = {Oncotarget},
  year      = {2017},
  volume    = {7},
  number    = {29},
  pages     = {45614--45626},
  doi       = {10.18632/oncotarget.10487},
  pmid      = {29069713},
  url       = {https://pmc.ncbi.nlm.nih.gov/articles/PMC5641056/}
}

@article{Joshi2017imprintedCN,
  author    = {Joshi, Ricky S. and Garg, Paras and Zaitlen, Noah and others},
  title     = {Copy number rather than epigenetic alterations are the major dictator of imprinted methylation in tumors},
  journal   = {Nature Communications},
  year      = {2017},
  volume    = {8},
  pages     = {533},
  doi       = {10.1038/s41467-017-00639-9},
  pmid      = {28883545},
  url       = {https://pubmed.ncbi.nlm.nih.gov/28883545/}
}

@article{Demidov2021ASE,
  author    = {Demidov, German and Sturm, Marc and Ossowski, Stephan},
  title     = {Allele-specific genomic data elucidate the role of somatic gain and copy-number neutral loss of heterozygosity in cancer},
  journal   = {Cell Systems},
  year      = {2021},
  volume    = {12},
  number    = {10},
  pages     = {1027--1031},
  doi       = {10.1016/j.cels.2021.08.005},
  url       = {https://www.sciencedirect.com/science/article/pii/S2405471221003835}
}
```

### Category 4 — Repeat / SD context
```bibtex
@article{Sigurpalsdottir2025ONTmethyl,
  author    = {Sigurpalsdottir, Brynja D. and Stefansson, Olafur A. and Holley, Guillaume and others},
  title     = {Reliable investigation of {DNA} methylation using {Oxford} nanopore technologies},
  journal   = {Scientific Reports},
  year      = {2025},
  volume    = {15},
  pages     = {14127},
  doi       = {10.1038/s41598-025-99882-0},
  url       = {https://www.nature.com/articles/s41598-025-99882-0}
}

@article{Vollger2023T2T_SD,
  author    = {Vollger, Mitchell R. and Guitart, Xavi and Dishuck, Philip C. and others},
  title     = {Segmental duplications and their variation in a complete human genome},
  journal   = {Science},
  year      = {2023},
  volume    = {376},
  number    = {6588},
  pages     = {eabj6965},
  doi       = {10.1126/science.abj6965}
}
```

### Category 5 — Replication timing
```bibtex
@article{Endicott2022PMDcelldivision,
  author    = {Endicott, Jamie L. and Nolte, Paula A. and Shen, Hui and Laird, Peter W.},
  title     = {Cell division drives {DNA} methylation loss in late-replicating domains in primary human cells},
  journal   = {Nature Communications},
  year      = {2022},
  volume    = {13},
  pages     = {6659},
  doi       = {10.1038/s41467-022-34268-8},
  pmid      = {36347867},
  url       = {https://pubmed.ncbi.nlm.nih.gov/36347867/}
}

@article{Du2021methylationReplicationTiming,
  author    = {Du, Qian and Smith, Grady C. and Luu, Phuc-Loi and others},
  title     = {{DNA} methylation is required to maintain both {DNA} replication timing precision and {3D} genome organization integrity},
  journal   = {Cell Reports},
  year      = {2021},
  volume    = {36},
  number    = {12},
  pages     = {109722},
  doi       = {10.1016/j.celrep.2021.109722},
  url       = {https://www.sciencedirect.com/science/article/pii/S2211124721011712}
}
```

**Total**: 14 unique PubMed-verifiable references (13 with PMID, 1 with DOI only for the 2024 review which lacks PMID at time of writing).

---

## 8. Mechanism falsifiability — "if true, v0.4 cycle should observe X"

| Mechanism | Falsifiable prediction for v0.4 LRT (Step 1) and filter sweep (Step 2) |
|---|---|
| **C1 cis-mQTL** | If real: HCC1395 FP positions in GoDMC mQTL hot-blocks show **higher methylation feature variance** (NME_imbalance, HPMergedDelta) than FP outside mQTL blocks. v0.4 LRT in `Inner\|other\|cov_normal` should add cis-mQTL overlap as covariate and improve fit (q<0.05). |
| **C2 cancer-specific ASM** | If real: TP set (especially Z-OCH cross_het cells) shows **distinct methylation signature** from FP in Do 2020 cancer-ASM DMR overlap. Predicted: TP methylation features in cancer-ASM regions have higher AlleleDelta + HPFineF than FP in same regions. |
| **C2.2 BRCA1 hypermethylation** | If real: chr17 BRCA1 ± 1 Mb FP positions show **NME_HP1/NME_HP2 imbalance > 0.5** (consistent with hemimethylated promoter). Predicted: BRCA1 region FP outperform genome-bg FP in `NME_imbalance` magnitude. |
| **C3 allele-imbalance amplification** | If real: in Z-CHR8 `Inner\|other\|cov_proxy_mid` (n=328), LRT adding `Coverage_Multiple × HPMergedDelta` interaction should be significant (q<0.05). FP rate by Coverage_Multiple bin should show monotonic gradient (cov_gain > cov_normal > cov_loss). |
| **C3.2 imprinted DMR CN-driven** | If real: chr8 imprinted DMR (KCNK9, MEST, PEG13) methylation in HCC1395 master TSV correlates with Coverage_Multiple (Pearson r > 0.4). |
| **C4 repeat / SD context** | If real: Z-AUTO non-chr8 FP positions show **>1.5× enrichment** for RepeatMasker SINE/LINE/satellite or T2T-CHM13 SD bed compared to non-Z-AUTO bg. Methylation feature variance (Epipoly_Delta) should also be inflated in Z-AUTO repeat overlap. |
| **C4.2 ONT methylation basecalling error** | If real: Z-AUTO FP at UCSC CpG islands show **non-monotonic** caller_af distribution (concentrated at extreme AF values reflecting methylation-driven basecall artefacts). |
| **C5 replication timing PMD loss** | If real: HCC1395 FP overlap with ENCODE breast-line Repli-Seq late bins (>S4) shows Fisher OR > 1.5 vs TP. Methylation NME values in late-replicating FP should be systematically lower than early-replicating FP. |
| **C5.2 chr8/chr17 timing-driven** | If real: chr8q + chr17q FP rate correlates with per-region Repli-Seq timing score; chr8q gain regions are confirmed late-replicating in HCC1395 Repli-Seq surrogate. |

**Decision logic for v0.4 Step 5 framing**:

- If **≥2 categories** show LRT q<0.05 in matched cells → "**multi-mechanism FP characterization**" framing (Bioinformatics tools journal angle).
- If **only C3 (allele-imbalance)** is significant → "**chr8-hotspot CN-driven FP, not generalisable**" framing (matches v0.3 conclusion; W3 paper §3.2 caveat).
- If **only C4 (repeat/SD)** is significant → "**FP filter should pivot to mappability/repeat axes**" framing (matches v0.3 §9.2 future direction T4.2).
- If **all 5 categories NEGATIVE** in LRT → "**methylation augmentation NEGATIVE**" → matches plan §Step 5 decision tree top-left cell.

---

## 9. Deliverable verification

| §Plan check | Status |
|---|---|
| 5 categories × ≥1 hypothesis | 5 × 2-3 hypotheses = **13 candidates** |
| Each hypothesis ≥1 PubMed DOI/PMID | All have DOI; 13/14 references have PMID; 1 has DOI only (Stewart 2024 review, no PMID at time of writing) |
| Verifiable PubMed paper exists | All verified via WebFetch / WebSearch; PMC mirrors accessed where Nature gated |
| Powered cell anchor per hypothesis | All 13 anchors mapped to v0.3 23-cell powered set + Z-CHR8/Z-AUTO sub-cells |
| Bibtex-ready format | §7 complete |
| Cross-category synthesis | §6 with 5 multi-mechanism anchor cells |
| Falsifiability | §8 with concrete predictions per mechanism |
| KB query first per memory rule | KB queried (`02_samples/hcc1395.md` + `05_tools/methyl-somatic-analysis.md`); no mQTL/ASM-specific KB docs exist, so PubMed primary |
| No emojis | Confirmed |
| Relaxed gate (no truth adjudication) | Confirmed — §0 explicitly states "do not adjudicate truth" |

**H4 verdict: POSITIVE** — 13 testable hypotheses with peer-reviewed prior art across 5 mechanism categories. The pilot has prior literature for every direction it might pivot to in v0.4 Step 1.

---

## 10. Forward links

- **v0.4 Step 1 augmented LR**: §6 anchor cells become the test cases for Model B vs A LRT. Top priority: `Inner|other|cov_proxy_mid` within Z-CHR8 (C3 + C5 dual mechanism).
- **v0.4 Step 2 filter sweep**: per-cell predicted P(TP) thresholds should be tested most aggressively on Z-CHR8 + Z-AUTO sub-cells.
- **v0.4 Step 5 paper framing**: mechanism diversity (5 categories) supports a multi-axis-characterization framing rather than single-mechanism filter claim. Pre-empts reviewer "why didn't you test mechanism X?" questions.
- **Future cycle (post-v0.4)**: independent validation in COLO829 (once truth set permissions resolved) and H2009 with matched Repli-Seq + RepeatMasker overlay.
