<!--
建立時間: 2026-08-12
目標: fresh 重算 chr2:18M 個案中的 somatic linkage、LOH、methylation 與可辨識性邊界
處理範圍: HCC1395 單樣本、HKU/DORADO 技術路徑、chr2:18,066,481–18,110,828；不作全樣本生物外推
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/chr2_18M_independent_audit_fresh.json
  - InterSubMod/docs/explain/04_subclone-reconstruction-chr2-18M.standalone.html
  - InterSubMod/docs/explain/05_subclone-correction-audit-chr2-18M.standalone.html
-->

# Independent chr2:18M Subclone Audit

## Coordinate correction

- Figure `(4.1)` is `chr2:18,096,341` (CpG G coordinate; C coordinate `chr2:18,096,340`).
- Transcribed `chr2:18,096,041` is an `A`, not a CpG.

## SEQC2 truth, HC, and LOH status

- Truth VCF: `/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- HC BED: `/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`
- LOH BED: `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`
- Overlapping LOH intervals: `[{'chrom': 'chr2', 'start0': 16146119, 'end0': 22100000, 'extra': ['loh']}]`

| Site | Position | Status | In HC | Truth records |
|---|---:|---|---|---|
| 1 | 18068480 | `truth_confirmed` | `True` | `[{'ref': 'C', 'alts': ['G'], 'filters': ['PASS', 'HighConf']}]` |
| 2 | 18072546 | `truth_confirmed` | `True` | `[{'ref': 'G', 'alts': ['C'], 'filters': ['PASS', 'HighConf']}]` |
| 3 | 18086020 | `out_of_hc_truth_unevaluable` | `False` | `[]` |
| 4 | 18096269 | `truth_confirmed` | `True` | `[{'ref': 'C', 'alts': ['G'], 'filters': ['PASS', 'MedConf']}]` |
| 5 | 18099697 | `truth_confirmed` | `True` | `[{'ref': 'G', 'alts': ['C'], 'filters': ['PASS', 'HighConf']}]` |
| 6 | 18108828 | `truth_confirmed` | `True` | `[{'ref': 'C', 'alts': ['G'], 'filters': ['PASS', 'HighConf']}]` |

## HKU_T

- Role: `tumor`
- Records / unique / duplicate extras: `340` / `280` / `60`
- HP parent counts: `{'2': 232, '1': 4, 'None': 44}`

### BQ20 allele counts

| Site | Position | REF | Counts |
|---|---:|---|---|
| 1 | 18068480 | C | `{'C': 22, 'G': 16}` |
| 2 | 18072546 | G | `{'C': 21, 'G': 25}` |
| 3 | 18086020 | G | `{'A': 29, 'G': 30}` |
| 4 | 18096269 | C | `{'T': 13, 'C': 21, 'DEL': 10, 'G': 4}` |
| 5 | 18099697 | G | `{'G': 23, 'C': 11, 'DEL': 1}` |
| 6 | 18108828 | C | `{'C': 20, 'G': 9}` |

### Key linkage tables

`00/10/01/11` indicate absence/presence of the two events among reads that cover both sites.

| Pair | n | 00 | 10 | 01 | 11 |
|---|---:|---:|---:|---:|---:|
| E1_vs_E2 | 31 | 18 | 0 | 0 | 13 |
| E1_vs_E3 | 13 | 3 | 3 | 7 | 0 |
| E2_vs_E3 | 17 | 3 | 5 | 9 | 0 |
| E1_vs_E6 | 1 | 1 | 0 | 0 | 0 |
| E2_vs_E6 | 2 | 2 | 0 | 0 | 0 |
| E3_vs_E5 | 13 | 8 | 1 | 0 | 4 |
| E3_vs_E4ANY | 22 | 0 | 9 | 13 | 0 |
| E3_vs_E6 | 6 | 0 | 3 | 3 | 0 |
| E4ANY_vs_E5 | 32 | 5 | 18 | 9 | 0 |
| E4ANY_vs_E6 | 19 | 9 | 0 | 0 | 10 |
| E4G_vs_E6 | 19 | 9 | 0 | 8 | 2 |
| E4T_vs_E6 | 19 | 9 | 0 | 5 | 5 |
| E4D_vs_E6 | 19 | 9 | 0 | 7 | 3 |
| E5_vs_E6 | 21 | 3 | 7 | 11 | 0 |

- Lineage anchors: `{'alpha_pos3_A': 30, 'pos3_ref_G': 30, 'beta_strict_E1_or_E2_or_E6': 49, 'alpha_and_beta_strict_violation': 0}`
- All-REF coverage check: `{'2': {'checked': 76, 'all_ref': 17}, '3': {'checked': 35, 'all_ref': 4}, '4': {'checked': 9, 'all_ref': 0}, '5': {'checked': 4, 'all_ref': 0}, '6': {'checked': 1, 'all_ref': 0}}`
- Position 4 subtype counts: `{'G': {'n': 4, 'strand_counts': {'-': 4}, 'with_E6_G': 2}, 'T': {'n': 14, 'strand_counts': {'+': 9, '-': 5}, 'with_E6_G': 5}, 'DEL': {'n': 10, 'strand_counts': {'-': 7, '+': 3}, 'with_E6_G': 3}}`

### Methylation: pos3 A vs pos3 reference G

| CpG | A n/mean | G n/mean | delta | MW FDR | Fisher FDR |
|---|---|---|---:|---:|---:|
| 2.1 | 9/0.751 | 8/0.483 | 0.267 | 0.176 | 0.335 |
| 2.2 | 14/0.901 | 9/0.224 | 0.677 | 0.00318 | 0.00149 |
| 3.1 | 20/0.222 | 19/0.917 | -0.695 | 0.000327 | 7.37e-06 |
| 3.2 | 19/0.97 | 17/0.144 | 0.826 | 4.14e-05 | 2.44e-07 |
| 3.3 | 14/0.416 | 16/0.844 | -0.428 | 0.025 | 0.0208 |
| 3.4 | 14/0.0448 | 15/0.874 | -0.829 | 4.14e-05 | 7.37e-06 |
| 3.5 | 9/0.00959 | 13/0.865 | -0.856 | 0.00173 | 5.03e-05 |
| 4.1 | 9/0.98 | 12/0.143 | 0.837 | 0.00273 | 6.8e-05 |
| 5.1 | 5/0.0298 | 8/0.95 | -0.92 | 0.00599 | 0.0013 |
| 5.2 | 5/0.029 | 7/0.994 | -0.965 | 0.00599 | 0.00158 |

## HKU_N

- Role: `normal`
- Records / unique / duplicate extras: `390` / `195` / `195`
- HP parent counts: `{'1': 66, '2': 76, 'None': 53}`

### BQ20 allele counts

| Site | Position | REF | Counts |
|---|---:|---|---|
| 1 | 18068480 | C | `{'C': 18, 'DEL': 1}` |
| 2 | 18072546 | G | `{'G': 21}` |
| 3 | 18086020 | G | `{'G': 18}` |
| 4 | 18096269 | C | `{'C': 21, 'T': 1}` |
| 5 | 18099697 | G | `{'G': 17}` |
| 6 | 18108828 | C | `{'C': 19}` |

### Matched-normal methylation: HP1 vs HP2

| CpG | HP1 n/mean | HP2 n/mean | delta | MW FDR | Fisher FDR |
|---|---|---|---:|---:|---:|
| 2.1 | 10/0.0275 | 9/0.0375 | -0.01 | 0.412 | 1 |
| 2.2 | 14/0.0448 | 9/0.00566 | 0.0392 | 0.0549 | 1 |
| 3.1 | 11/0.0538 | 12/0.0373 | 0.0166 | 0.412 | 1 |
| 3.2 | 12/0.0454 | 11/0.12 | -0.0747 | 0.496 | 0.957 |
| 3.3 | 11/0.814 | 11/0.046 | 0.768 | 0.00703 | 0.000737 |
| 3.4 | 10/0.0329 | 10/0.85 | -0.817 | 0.00703 | 0.000595 |
| 3.5 | 10/0.583 | 8/0.0118 | 0.571 | 0.0496 | 0.0321 |
| 4.1 | 10/0.887 | 9/0.078 | 0.809 | 0.00703 | 0.000595 |
| 5.1 | 8/0.0397 | 7/0.0325 | 0.00721 | 0.299 | 1 |
| 5.2 | 8/0.0255 | 6/0.0301 | -0.00458 | 0.412 | 1 |

## DORADO_TAGGED_T

- Role: `tumor`
- Records / unique / duplicate extras: `341` / `341` / `0`
- HP parent counts: `{'2': 279, 'None': 61, '1': 1}`

### BQ20 allele counts

| Site | Position | REF | Counts |
|---|---:|---|---|
| 1 | 18068480 | C | `{'G': 18, 'C': 22}` |
| 2 | 18072546 | G | `{'DEL': 1, 'G': 27, 'C': 18}` |
| 3 | 18086020 | G | `{'G': 27, 'A': 21}` |
| 4 | 18096269 | C | `{'T': 15, 'C': 29, 'DEL': 7, 'G': 6}` |
| 5 | 18099697 | G | `{'C': 18, 'G': 22, 'DEL': 1}` |
| 6 | 18108828 | C | `{'C': 43, 'G': 10, 'DEL': 3}` |

### Key linkage tables

`00/10/01/11` indicate absence/presence of the two events among reads that cover both sites.

| Pair | n | 00 | 10 | 01 | 11 |
|---|---:|---:|---:|---:|---:|
| E1_vs_E2 | 18 | 11 | 1 | 0 | 6 |
| E1_vs_E3 | 0 | 0 | 0 | 0 | 0 |
| E2_vs_E3 | 1 | 0 | 1 | 0 | 0 |
| E1_vs_E6 | 0 | 0 | 0 | 0 | 0 |
| E2_vs_E6 | 0 | 0 | 0 | 0 | 0 |
| E3_vs_E5 | 5 | 1 | 2 | 0 | 2 |
| E3_vs_E4ANY | 13 | 0 | 7 | 5 | 1 |
| E3_vs_E6 | 0 | 0 | 0 | 0 | 0 |
| E4ANY_vs_E5 | 25 | 6 | 10 | 9 | 0 |
| E4ANY_vs_E6 | 2 | 1 | 0 | 0 | 1 |
| E4G_vs_E6 | 2 | 1 | 0 | 1 | 0 |
| E4T_vs_E6 | 2 | 1 | 0 | 1 | 0 |
| E4D_vs_E6 | 2 | 1 | 0 | 0 | 1 |
| E5_vs_E6 | 8 | 1 | 7 | 0 | 0 |

- Lineage anchors: `{'alpha_pos3_A': 28, 'pos3_ref_G': 27, 'beta_strict_E1_or_E2_or_E6': 44, 'alpha_and_beta_strict_violation': 0}`
- All-REF coverage check: `{'2': {'checked': 60, 'all_ref': 15}, '3': {'checked': 6, 'all_ref': 0}, '4': {'checked': 0, 'all_ref': 0}, '5': {'checked': 0, 'all_ref': 0}, '6': {'checked': 0, 'all_ref': 0}}`
- Position 4 subtype counts: `{'G': {'n': 6, 'strand_counts': {'+': 2, '-': 4}, 'with_E6_G': 0}, 'T': {'n': 16, 'strand_counts': {'-': 12, '+': 4}, 'with_E6_G': 0}, 'DEL': {'n': 7, 'strand_counts': {'-': 4, '+': 3}, 'with_E6_G': 1}}`

### Methylation: pos3 A vs pos3 reference G

| CpG | A n/mean | G n/mean | delta | MW FDR | Fisher FDR |
|---|---|---|---:|---:|---:|
| 2.1 | 1/0.0118 | 1/0.0157 | -0.00392 | NA | NA |
| 2.2 | 5/0.802 | 2/0.00784 | 0.794 | 0.135 | 0.143 |
| 3.1 | 15/0.0742 | 13/0.998 | -0.924 | 2.41e-05 | 2.68e-06 |
| 3.2 | 15/1 | 9/0.00566 | 0.994 | 2.41e-05 | 2.68e-06 |
| 3.3 | 12/0.257 | 9/0.993 | -0.736 | 0.00609 | 0.00261 |
| 3.4 | 8/0.256 | 7/0.858 | -0.602 | 0.0601 | 0.0568 |
| 3.5 | 8/0.38 | 6/0.835 | -0.455 | 0.084 | 0.143 |
| 4.1 | 8/0.971 | 5/0.403 | 0.567 | 0.0839 | 0.0568 |
| 5.1 | 2/0.502 | 1/1 | -0.498 | NA | NA |
| 5.2 | 2/0.0118 | 1/1 | -0.988 | NA | NA |

## DORADO_RAW_T

- Role: `tumor`
- Records / unique / duplicate extras: `341` / `341` / `0`
- HP parent counts: `{'None': 341}`

### BQ20 allele counts

| Site | Position | REF | Counts |
|---|---:|---|---|
| 1 | 18068480 | C | `{'G': 18, 'C': 22}` |
| 2 | 18072546 | G | `{'DEL': 1, 'G': 27, 'C': 18}` |
| 3 | 18086020 | G | `{'G': 27, 'A': 21}` |
| 4 | 18096269 | C | `{'T': 15, 'C': 29, 'DEL': 7, 'G': 6}` |
| 5 | 18099697 | G | `{'C': 18, 'G': 22, 'DEL': 1}` |
| 6 | 18108828 | C | `{'C': 43, 'G': 10, 'DEL': 3}` |

### Key linkage tables

`00/10/01/11` indicate absence/presence of the two events among reads that cover both sites.

| Pair | n | 00 | 10 | 01 | 11 |
|---|---:|---:|---:|---:|---:|
| E1_vs_E2 | 18 | 11 | 1 | 0 | 6 |
| E1_vs_E3 | 0 | 0 | 0 | 0 | 0 |
| E2_vs_E3 | 1 | 0 | 1 | 0 | 0 |
| E1_vs_E6 | 0 | 0 | 0 | 0 | 0 |
| E2_vs_E6 | 0 | 0 | 0 | 0 | 0 |
| E3_vs_E5 | 5 | 1 | 2 | 0 | 2 |
| E3_vs_E4ANY | 13 | 0 | 7 | 5 | 1 |
| E3_vs_E6 | 0 | 0 | 0 | 0 | 0 |
| E4ANY_vs_E5 | 25 | 6 | 10 | 9 | 0 |
| E4ANY_vs_E6 | 2 | 1 | 0 | 0 | 1 |
| E4G_vs_E6 | 2 | 1 | 0 | 1 | 0 |
| E4T_vs_E6 | 2 | 1 | 0 | 1 | 0 |
| E4D_vs_E6 | 2 | 1 | 0 | 0 | 1 |
| E5_vs_E6 | 8 | 1 | 7 | 0 | 0 |

- Lineage anchors: `{'alpha_pos3_A': 28, 'pos3_ref_G': 27, 'beta_strict_E1_or_E2_or_E6': 44, 'alpha_and_beta_strict_violation': 0}`
- All-REF coverage check: `{'2': {'checked': 60, 'all_ref': 15}, '3': {'checked': 6, 'all_ref': 0}, '4': {'checked': 0, 'all_ref': 0}, '5': {'checked': 0, 'all_ref': 0}, '6': {'checked': 0, 'all_ref': 0}}`
- Position 4 subtype counts: `{'G': {'n': 6, 'strand_counts': {'+': 2, '-': 4}, 'with_E6_G': 0}, 'T': {'n': 16, 'strand_counts': {'-': 12, '+': 4}, 'with_E6_G': 0}, 'DEL': {'n': 7, 'strand_counts': {'-': 4, '+': 3}, 'with_E6_G': 1}}`

### Methylation: pos3 A vs pos3 reference G

| CpG | A n/mean | G n/mean | delta | MW FDR | Fisher FDR |
|---|---|---|---:|---:|---:|
| 2.1 | 1/0.0118 | 1/0.0157 | -0.00392 | NA | NA |
| 2.2 | 5/0.802 | 2/0.00784 | 0.794 | 0.135 | 0.143 |
| 3.1 | 15/0.0742 | 13/0.998 | -0.924 | 2.41e-05 | 2.68e-06 |
| 3.2 | 15/1 | 9/0.00566 | 0.994 | 2.41e-05 | 2.68e-06 |
| 3.3 | 12/0.257 | 9/0.993 | -0.736 | 0.00609 | 0.00261 |
| 3.4 | 8/0.256 | 7/0.858 | -0.602 | 0.0601 | 0.0568 |
| 3.5 | 8/0.38 | 6/0.835 | -0.455 | 0.084 | 0.143 |
| 4.1 | 8/0.971 | 5/0.403 | 0.567 | 0.0839 | 0.0568 |
| 5.1 | 2/0.502 | 1/1 | -0.498 | NA | NA |
| 5.2 | 2/0.0118 | 1/1 | -0.988 | NA | NA |

## DORADO_N

- Role: `normal`
- Records / unique / duplicate extras: `167` / `167` / `0`
- HP parent counts: `{'None': 167}`

### BQ20 allele counts

| Site | Position | REF | Counts |
|---|---:|---|---|
| 1 | 18068480 | C | `{'C': 21}` |
| 2 | 18072546 | G | `{'G': 22, 'DEL': 1}` |
| 3 | 18086020 | G | `{'G': 27}` |
| 4 | 18096269 | C | `{'C': 27, 'DEL': 1}` |
| 5 | 18099697 | G | `{'G': 24}` |
| 6 | 18108828 | C | `{'C': 16}` |

### Matched-normal methylation: HP1 vs HP2

| CpG | HP1 n/mean | HP2 n/mean | delta | MW FDR | Fisher FDR |
|---|---|---|---:|---:|---:|
| 2.1 | 0/NA | 0/NA | NA | NA | NA |
| 2.2 | 0/NA | 0/NA | NA | NA | NA |
| 3.1 | 0/NA | 0/NA | NA | NA | NA |
| 3.2 | 0/NA | 0/NA | NA | NA | NA |
| 3.3 | 0/NA | 0/NA | NA | NA | NA |
| 3.4 | 0/NA | 0/NA | NA | NA | NA |
| 3.5 | 0/NA | 0/NA | NA | NA | NA |
| 4.1 | 0/NA | 0/NA | NA | NA | NA |
| 5.1 | 0/NA | 0/NA | NA | NA | NA |
| 5.2 | 0/NA | 0/NA | NA | NA | NA |

## Cross-sample methylation: tumor versus matched normal

These are all-read tumor/normal comparisons. They establish sample-level differential methylation, but do not by themselves separate haplotype ASM from tumor-acquired methylation.

### HKU_T_vs_HKU_N_all_reads

| CpG | Tumor n/mean | Normal n/mean | delta | MW FDR | Fisher FDR |
|---|---|---|---:|---:|---:|
| 2.1 | 50/0.525 | 22/0.0299 | 0.495 | 0.000771 | 9.02e-06 |
| 2.2 | 50/0.585 | 23/0.0295 | 0.555 | 0.000429 | 2.36e-06 |
| 3.1 | 50/0.561 | 23/0.0452 | 0.516 | 0.00119 | 4.57e-06 |
| 3.2 | 49/0.568 | 23/0.0812 | 0.487 | 0.00605 | 1.75e-05 |
| 3.3 | 50/0.594 | 22/0.43 | 0.164 | 0.488 | 0.249 |
| 3.4 | 54/0.467 | 20/0.442 | 0.025 | 0.945 | 1 |
| 3.5 | 45/0.508 | 21/0.284 | 0.224 | 0.527 | 0.0995 |
| 4.1 | 49/0.552 | 23/0.507 | 0.0445 | 0.957 | 0.891 |
| 5.1 | 45/0.56 | 21/0.0332 | 0.526 | 0.000771 | 1.08e-05 |
| 5.2 | 38/0.51 | 20/0.0298 | 0.48 | 0.00119 | 9.85e-05 |

### DORADO_TAGGED_T_vs_DORADO_N_all_reads

| CpG | Tumor n/mean | Normal n/mean | delta | MW FDR | Fisher FDR |
|---|---|---|---:|---:|---:|
| 2.1 | 45/0.497 | 34/0.016 | 0.481 | 0.000353 | 1.41e-06 |
| 2.2 | 46/0.68 | 28/0.099 | 0.581 | 0.000353 | 1.41e-06 |
| 3.1 | 58/0.476 | 23/0.0171 | 0.459 | 0.000997 | 3.87e-05 |
| 3.2 | 62/0.551 | 27/0.0954 | 0.456 | 0.000353 | 3.87e-05 |
| 3.3 | 59/0.561 | 25/0.283 | 0.279 | 0.0324 | 0.0382 |
| 3.4 | 57/0.531 | 24/0.546 | -0.015 | 0.837 | 1 |
| 3.5 | 58/0.505 | 28/0.292 | 0.212 | 0.0324 | 0.0754 |
| 4.1 | 57/0.758 | 28/0.436 | 0.322 | 0.00609 | 0.00597 |
| 5.1 | 47/0.429 | 28/0.00602 | 0.423 | 0.00609 | 3.87e-05 |
| 5.2 | 45/0.352 | 26/0.0486 | 0.304 | 0.0324 | 0.00584 |
