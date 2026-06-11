# HD-4: AF→NGroups Attribution — Methylation Biology vs Phasing Artifact

> **Date:** 2026-06-08 · **Cohort:** paired mode, LOH-bed TP, CN1 (matches step2 H1p), n=36,854
> **All numbers:** `research/loh_subclone_af_paired/data/hd4_attribution.json` (+ `hd4_supp.json`)
> **Scripts:** `scripts/hd4_attribution.py`

## Verdict: **PHASING-DRIVEN (mechanism b)** — a documented haplotype-occupancy confound, NOT a methylation-biology finding.

The "intermediate-AF → higher NGroups" association is real and reproduces (r=0.656), but
NGroups mechanically counts populated **phasing sub-families**, and AF mechanically controls
how many sub-families can be populated. It does not measure methylation diversity.

---

## 1. What NGroups actually measures (this alone largely answers HD-4)

`HPFineNGroups` = `LabelTest::hp_to_fine_labels()` → **count of distinct populated sub-families
among {HP1, HP1-1, HP2, HP2-1}** (max 4), with a `min_reads_per_group` threshold.
- Source: `src/core/LabelTest.cpp:265-305` (label assignment) and `:607-633` (`fine_n_groups = unique_groups.size()`).
- Each read's group comes **only from its HP tag string** (`ReadParser.cpp:128-138`: longphase-TO integer 1/2/11/21/33 → "1"/"2"/"1-1"/"2-1"/"3").
- **The methylation distance matrix is used ONLY for the separate PERMANOVA F-stat (`HPFineF`), NOT for the count.** No methylation enters NGroups.
- `ReadParser.cpp:145-148` code comment (project's own): *"LongPhase-TO marks somatic-only phase blocks with HP:i:11/21/33; these reflect the somatic variant phasing itself (circular dependency)..."* → the `-1` sub-families ARE the somatic-variant phasing.

**Conclusion:** NGroups is a haplotype/phasing-tag occupancy count. In an LOH region (single
germline haplotype), going from 1→2 populated sub-families requires the **somatic ALT
sub-family** to appear — which is an arithmetic function of allele frequency.

## 2. Quantitative confirmation (landed data)

| Test | Result | File |
|---|---|---|
| Baseline NGroups~AF-centrality (Spearman) | **r=0.656** | hd4_attribution.json |
| NGroups is a capped count (0–4) | 99.5% are ≤2; max=4 | hd4_attribution.json |
| η² of 3-level AF class on NGroups | **0.455** (near step-function) | hd4_attribution.json |
| Within-class continuous AF→NGroups | Extreme r=0.13, Inter r=0.21, Near-half **r=−0.01** | hd4_attribution.json |
| Entire effect = binary "2nd sub-family appears" (NGroups≥2) vs AF-centrality | **point-biserial r=0.80** | hd4_attribution.json |
| % with NGroups≥2: Extreme→Inter | **2.4% → 76.6%** | step2_paired_ngroups_by_af_class.tsv |

The association collapses to a single binary event (a 2nd haplotype sub-family becoming
populated), which is exactly what AF controls. There is no graded, continuous methylation
signal underneath.

## 3. Why the partial-correlation controls are not decisive on landed data

- **Controlling for methylation features** (HPFineF, AlleleDelta, PairwiseMeanDist, CramersV,
  ClusterPermanovaF) drops r 0.656→0.322. **This is NOT evidence for methylation** — those
  features are *downstream mediators/colliders* of having ≥2 groups (e.g. HPFineF can only be
  nonzero when ≥2 groups exist). Controlling for them is mechanistically the wrong direction.
- **Controlling for the available phasing proxies** (hp_balance, hp_minority_frac, hp0/hp3
  ratios) barely moves it (0.656→0.640). This is **not** because the effect is methylation —
  it is because the master TSV only exports germline-family read counts (`HP1FamilyN`=HP1+HP1-1,
  `HP2FamilyN`=HP2+HP2-1) and **does NOT split out the somatic sub-families (HP1-1/HP2-1)** that
  actually constitute the NGroups increment. The right occupancy variable is not in the data.

## 4. Computable vs needs-a-run

**Computable on landed data (done):** NGroups definition trace, baseline corr, capped-count
distribution, η², within-class corr, binary-event point-biserial, by-class occupancy %.

**Needs a new ISM run to fully close (both already computed in C++, just not exported to TSV):**
1. Export `fine_group_counts` = per-variant {HP1, HP1-1, HP2, HP2-1} read counts (`Stats.hpp:372`).
   Then partial-out `has_somatic_subfamily` exactly; NGroups~AF | sub-family-presence should → ~0,
   mechanically nailing (b).
2. Export `n_clusters` (methylation-ONLY hierarchical/silhouette cluster count,
   `HierarchicalClustering.hpp` / `RegionProcessor.cpp:1689`). Test whether a *true methylation-derived*
   group count tracks AF at all. If it doesn't, that is the clean negative control for (a).

## 5. Consistency with prior record

- `docs/experiments/INDEX.md`: HPFineNGroups documented as "somatic heterogeneity marker"; a
  separate analysis already found it "68% NumReads confound (raw 0.495→resid 0.160)".
- MEMORY `project_hpfinengroups_subclone_marker.md`: ⭐4→⭐3 downgrade, mechanism reinterpreted as a
  **phasing signature** (pipeline-dependent). HD-4 supplies the mechanistic + quantitative basis for that.

## Bottom line for the subclone thread
The AF→NGroups "multi-sample positive" should be reported as a **phasing/haplotype-occupancy
signature** (genetic, AF-mechanical), **not** as an epigenetic/methylation finding. It does not
support a claim that intermediate-AF subclones carry distinct methylation patterns. Any
methylation claim must come from a methylation-only group count or the HPFineF-style distance
test, evaluated independently of the sub-family count.
