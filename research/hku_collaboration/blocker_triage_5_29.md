# HKU External-Delivery Blocker Triage — 2026-05-29

> Read-only investigation of 3 HKU standalone-HTML delivery blockers. No heavy pipeline / ISM / BAM job run (disk 97% full, ASM batch active).
> HKU standalone HTML: `/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/in_progress/2026/05/20260523_HKU_collab_tag_proposal_01.standalone.html`
> Date: 2026-05-29

---

## TASK 1 — Verify F4 HP 6-value counts (HCC1395 baseline)

**Source of truth**: `/big7_disk/liaoyoyo2001/InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A1_HP_6values_5sample.tsv`
Row: `HCC1395 / baseline / chr8+chr19`
BAM: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam`

| Quoted in HTML/docs | Source column | Exact source number | Verdict |
|---|---|---|---|
| HP1 = 863K | HP_1_count | **863065** | **MATCH** (863K) |
| HP2 = 266K | HP_2_count | **266191** | **MATCH** (266K) |
| HP1-1 = 79K | HP_11_count | **79118** | **MATCH** (79K) |
| HP2-1 = 23K | HP_21_count | **23007** | **MATCH** (23K) |
| HP3 = 630 | HP_33_count | **630** | **MATCH** (exact) |

**VERDICT: ALL 5 VALUES MATCH.** Quoted figures are correct rounding of the real source (HP3=630 is exact). Supporting source fields: total=2240163, HP1_HP2_ratio=3.2423, HP33_pct=0.0281, HP_untagged_count=1008152, elapsed_sec=1002.9.

Note: the figures README (`research/hku_collaboration/figures/README.md`) flagged these counts as "記憶引用 (post-haplotag-fix)，正式報告前須 cross-check evidence_ledger". This triage cross-checks them against the A1_HP_6values_5sample.tsv SoT directly — confirmed accurate. Placeholder caveat in README can be cleared.

---

## TASK 2 — base64 inline feasibility (external `<img>` in HKU standalone HTML)

**Total external `<img>` referenced: 18.** All resolve under `research/hku_collaboration/figures/` (HTML relative `../../../../../research/hku_collaboration/figures/<name>`). **All 18 exist.**

| # | Figure file | exists | bytes |
|---|---|---|---|
| 1 | F1_pipeline_flow.png | yes | 164602 |
| 2 | A2_1_ps_set_tp_count_dist.png | yes | 235541 |
| 3 | A2_2_common_reads_xy.png | yes | 296842 |
| 4 | A2_3_loh_hp_missing_vs_seqc2.png | yes | 165894 |
| 5 | F2v2_6state_with_tree_inset.png | yes | 216466 |
| 6 | F3v2_x_imputation_flow.png | yes | 415539 |
| 7 | F4v2_evolution_tree_4scenarios.png | yes | 279032 |
| 8 | F5_ism_readcpg_heatmap.png | yes | 105284 |
| 9 | T1_HP3_TP_rate_per_chr.png | yes | 66228 |
| 10 | F7_chr2_region_evidence.png | yes | 320802 |
| 11 | F6_caller_window_vs_phasing.png | yes | 198667 |
| 12 | A6_loh_cnv_stratified_hp_distribution.png | yes | 201130 |
| 13 | A6_subclone_ratio_violin.png | yes | 309883 |
| 14 | A7_phase_block_size_distribution.png | yes | 275500 |
| 15 | A7_n50_per_chr_bar.png | yes | 102805 |
| 16 | A8_chr_ideogram_hp_tag.png | yes | 206687 |
| 17 | A8_per_chr_hp_summary_table.png | yes | 58207 |
| 18 | A8_chr_loh_cnv_coverage.png | yes | 70431 |

**Total raw PNG bytes: 3,689,540 (3.51 MB).**
base64 encoding inflates ~33% → inlined HTML weight ≈ **4.9 MB** of data-URI text (plus existing HTML/JS).

**RECOMMENDATION: Inlining all 18 as base64 is FEASIBLE.** A single self-contained HTML of ~5 MB opens fine in modern browsers (Chrome/Firefox/Safari handle multi-MB data-URI documents without issue). This makes the standalone HTML truly portable for HKU external delivery (no broken `../` relative paths on their end). Caveats: (a) the document becomes non-diff-friendly and slower to load on very weak connections if served remotely — for an emailed/handed-off single file this is irrelevant; (b) re-encode whenever a figure changes. No blocker. (This task is report-only; HTML was NOT modified.)

---

## TASK 3 — F5 real-ISM availability

### F5 is confirmed SYNTHETIC
- Generator: `research/hku_collaboration/scripts/gen_concept_figures.py`, fn `fig5_ism_readcpg_heatmap()` (lines 501-605). Matrix built with `rng.beta(...)` (lines 508-515) — 10 read × 20 CpG random data.
- `research/hku_collaboration/figures/README.md` placeholder table explicitly: *"F5 內 10×20 矩陣 ... 示意數據 rng.beta · 非真實 HCC1395 ISM 輸出"*.

### Real ISM region outputs DO exist on disk — replacement is possible WITHOUT any new run

Real HCC1395 ISM per-region outputs (full schema: `methylation.csv` read×CpG matrix + `reads.tsv` + `distance/NHD/matrix.csv` + `clustering/tree.nwk` + `significance.json` + `metadata.txt`) exist in bulk under:
`/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/` (and `.../phaseC_genome_three_way/`).
Provenance (per that dir's README.md): **HCC1395, ClairS-TO ssrs + LongPhase-TO TO-mode, V3F/V5/V6/baseline binaries** — genuine somatic ISM output, NOT synthetic. (Dir name says "paired" for historical reasons but content is tumor-only.)

**Closest real region to target chr2:18,086,020:**

| Region dir (V3F_on_tp) | SNV pos | distance to target | matrix |
|---|---|---|---|
| `.../V3F_on_tp/filtered_snv_tp/chr2/chr2_18072546/chr2_18067546_18077546/` | chr2:18,072,546 | **~13.5 kb** | **43 reads × 53 CpG** |
| `.../chr2/chr2_18068480/...` | chr2:18,068,480 | ~17.5 kb | (full ISM dir present) |
| `.../chr2/chr2_18099697/...` | chr2:18,099,697 | ~13.7 kb | (full ISM dir present) |

**Recommended replacement source for F5:**
`/big7_disk/liaoyoyo2001/InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/V3F_on_tp/filtered_snv_tp/chr2/chr2_18072546/chr2_18067546_18077546/`

Verified contents (metadata.txt): real TP SNV chr2:18,072,546 G→C, SNV quality 0.9887, 43 reads (25 fwd / 18 rev), 53 CpG sites, matrix 43×53.
- `methylation/methylation.csv` (14157 B) — real read×CpG methylation-probability matrix, header = genomic CpG coordinates, NA for uncovered sites, values 0.00-1.00.
- `reads/reads.tsv` (3326 B) — read_id, read_name, chr, start, end, mapq, **hp**, alt_support, is_tumor, strand (real ONT read IDs + HP tags + ALT/REF support — supports the HP annotation bar that F5 currently fakes).
- `distance/NHD/matrix.csv` + `clustering/tree.nwk` / `linkage_matrix.csv` / `leaf_order.txt` — real hierarchical clustering already computed (the dendrogram F5 needs).

**No new ISM run is required to replace F5.** A pure plotting script (seaborn clustermap over the real `methylation.csv` + HP annotation from `reads.tsv` + dendrogram from `clustering/linkage_matrix.csv`) reproduces a genuine ISM read×CpG heatmap from real HCC1395 data — no BAM, no binary, no disk-heavy step. This sidesteps the disk-full / ASM-batch constraint entirely.

### If a region exactly centered on chr2:18,086,020 is later required (NOT currently on disk)
No region dir centered exactly on 18,086,020 exists (nearest is 18,072,546, ~13.5 kb away; the F6/F7 figures use 18,086,020 only as an illustrative caller-target label, not a core claim). A minimal ISM run to generate the exact region (deferred until disk frees up):
- **Binary**: `/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod`
- **Input BAM**: HCC1395 TO tagged BAM — `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam` (baseline) or the V3F equivalent matching the phaseC outputs.
- **Region**: chr2:18,081,020-18,091,020 (10 kb window centered on 18,086,020, matching the existing ±5 kb region convention).
- Single-region run — small footprint, but still deferred per the no-heavy-job constraint.

**RECOMMENDATION for F5 blocker: USE EXISTING real region chr2_18072546 (plot-only, no run needed).** Document the 13.5 kb offset from the 18,086,020 label, or relabel F5 to its true coordinate 18,072,546. Generating the exact 18,086,020 region is optional and can wait for disk headroom.

---

## Summary verdicts

| Task | Verdict |
|---|---|
| **F4 HP counts** | **ALL 5 MATCH** (863065 / 266191 / 79118 / 23007 / 630) vs A1_HP_6values_5sample.tsv |
| **base64 inline** | **18 imgs, all exist, 3,689,540 B (3.51 MB) raw → ~4.9 MB base64. FEASIBLE.** |
| **F5 real ISM** | **Real region AVAILABLE — no run needed.** chr2_18072546 (43×53 matrix, ~13.5 kb from target). Plot-only replacement. |

Sources cross-checked: A1_HP_6values_5sample.tsv · figures/README.md · gen_concept_figures.py · paired_priority_bug_audit/README.md · region metadata.txt/methylation.csv/reads.tsv. evidence_ledger.jsonl had no chr2:18072546 / phaseC entry (region outputs are intermediate artifacts, not ledger-tracked).
