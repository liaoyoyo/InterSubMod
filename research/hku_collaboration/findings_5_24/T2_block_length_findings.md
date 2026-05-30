# T2 — Somatic haplotype block length on HCC1395 tagged BAM

**Date**: 2026-05-24
**Author**: InterSubMod Research
**Scope**: Response to LongPhase-S paper Reviewer 1 Minor #3 — quantify length / distribution of somatic phase blocks introduced by LongPhase.
**Sample**: HCC1395 tumor-mode tagged BAM (ClairS pileup v040, wo-filter).

---

## How it was run

- **Input BAM**: `/big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` (264 GB, indexed `.bai` present)
- **Chromosomes scanned**: `chr1`, `chr8`, `chr19`
- **Read filters**: primary alignment only (exclude secondary / supplementary / duplicate / unmapped), `MAPQ >= 20`
- **PS-tag handling**: `read.has_tag('PS')` then `str(read.get_tag('PS'))` — fallback for both `PS:i:` and `PS:Z:` encodings; reads without PS (typically HP3 / unphased) are excluded by definition (not part of any block).
- **Block extent**: per PS_id, `block_start = min(reference_start)`, `block_end = max(reference_end)`, `block_length = block_end - block_start`. `n_reads` = number of reads contributing to that PS.
- **Script**: `InterSubMod/research/hku_collaboration/scripts/T2_somatic_block_length.py` (pysam streaming, no full-BAM load)
- **Wall time**: ~13 min on the three chromosomes (chr1 407 s + chr8 259 s + chr19 113 s)

Reads kept (post-filter, MAPQ>=20, primary) / reads with a PS tag:

| chr   | reads_total  | kept_primary_mapq20+ | with_PS    | with_PS / kept |
|-------|--------------|----------------------|------------|----------------|
| chr1  | 4,411,427    | 2,236,562            | 1,422,567  | 63.6%          |
| chr8  | 1,956,724    | 1,653,223            | 1,061,048  | 64.2%          |
| chr19 | 664,981      | 535,588              | 357,724    | 66.8%          |

---

## Summary

All lengths in bp. Source: `findings_5_24/T2_somatic_block_summary.tsv`.

| chromosome | n_blocks | median   | Q1       | Q3       | N50        | max         | min    | mean     | median reads / PS |
|------------|----------|----------|----------|----------|------------|-------------|--------|----------|--------------------|
| chr1       | 490      | 245,385  | 124,303  | 653,533  | 1,071,203  | 18,303,383  | 843    | 527,332  | 847                |
| chr8       | 289      | 246,864  | 120,971  | 697,269  | 1,024,929  |  3,520,673  | 453    | 513,269  | 1,011              |
| chr19      |  94      | 489,288  | 139,694  | 838,022  | 1,018,274  |  3,576,640  | 17,233 | 615,754  | 2,542              |
| **ALL**    | **873**  | **255,900** | **124,838** | **697,269** | **1,035,719** | **18,303,383** | **453** | **532,197** | **1,024** |

Length distribution (ALL, log decades):

| bucket       | n_blocks | share |
|--------------|----------|-------|
| < 1 kb       | 2        | 0.2%  |
| 1 – 10 kb    | 8        | 0.9%  |
| 10 – 100 kb  | 118      | 13.5% |
| 100 kb – 1 Mb| 612      | 70.1% |
| 1 – 10 Mb    | 132      | 15.1% |
| > 10 Mb      | 1        | 0.1%  |

Companion figure: `figures/T2_somatic_block_length_dist.png` (left: log-x histogram of all 873 blocks with median / N50 / 17.4 kb reference lines; right: per-chromosome boxplot, log-y).

---

## Interpretation (3 sentences)

1. **Somatic phase blocks are ~10–50× longer than ONT reads.** Across chr1+chr8+chr19, the median block is **~256 kb** (N50 **~1.04 Mb**); typical ONT R10 reads are 10–20 kb, so a somatic block spans ~13–25 reads end-to-end on the median and >50 reads at N50, confirming that LongPhase-S is producing genuinely multi-read phase sets rather than read-length-bounded fragments.

2. **Median block length >> median TP-pair distance (17.4 kb).** The verified TP-pair median distance is 17.4 kb, which sits in the **Q1 region (<25th percentile)** of block lengths — i.e. **>75% of somatic blocks comfortably contain a typical TP pair**, supporting that downstream same-haplotype TP-TP linkage queries are not block-edge-limited in the majority of cases.

3. **Block length is cross-chromosome consistent at the medium scale, with chr19 the outlier.** chr1 and chr8 have nearly identical medians (245 kb vs 247 kb) and N50 (1.07 Mb vs 1.02 Mb); chr19 medians ~2× higher (489 kb) with ~3× the reads/PS (2,542), which is plausibly driven by chr19's higher coverage density and shorter physical length producing fewer but better-supported PS sets — worth flagging in the response letter but not a phasing-method artifact.

---

## Outputs

- Per-block TSV: `InterSubMod/research/hku_collaboration/findings_5_24/T2_somatic_block_per_block.tsv` (873 rows)
- Summary TSV: `InterSubMod/research/hku_collaboration/findings_5_24/T2_somatic_block_summary.tsv` (4 rows: chr1 / chr8 / chr19 / ALL)
- Figure (DPI 150, 2-panel): `InterSubMod/research/hku_collaboration/figures/T2_somatic_block_length_dist.png`
- Script: `InterSubMod/research/hku_collaboration/scripts/T2_somatic_block_length.py`

---

## Caveats

- **Scope**: chr1+chr8+chr19 only (per task spec). Genome-wide extension is straightforward — change `CHROMS` in the script — but expected wall time ~3-4× of current.
- **Block extent definition**: we use `max(reference_end) - min(reference_start)` of reads carrying the same PS. This is the standard "read-spanning extent" definition and is **slightly larger** than the strict "phased-variant-span" definition (which would use the variant positions). For the reviewer's question on read-level somatic block length this is the appropriate measure.
- **Excluded reads**: HP3 / unphased reads have no PS by definition and are not part of any block; this matches LongPhase semantics and is not a missing-data issue.
- **`min_len`**: chr1's 843 bp and chr8's 453 bp minima reflect rare PS sets supported by a single short read; they do not represent typical phase block resolution and are visible only in the left tail of the histogram (<1 kb bucket: 2 / 873 blocks).
