---
title: Formal Wilcoxon signed-rank test — LOH-constrained phasing TP gap (NG=2 Inner same-hap vs Outer cross-het)
task: phasing-formal-wilcoxon-tpgap
date: 2026-05-29
author: InterSubMod Research (AI-assisted, read-only investigation)
status: COMPUTED — data existed, test re-verified independently from raw TSV
mode: read-only (disk 97% full on /big7_disk; no pipeline/BAM/ISM run executed)
---

# Formal Wilcoxon signed-rank test — LOH-constrained phasing TP gap

## TL;DR (BLUF)

The underlying per-sample data **already existed** on disk. I located it, parsed the
raw values, and **independently re-ran `scipy.stats.wilcoxon`** (scipy 1.13.1) from the
source TSV — no number is fabricated; every value below is parsed from a file.

> **Wilcoxon signed-rank, one-sided (Inner same-hap TP > Outer cross-het TP):**
> **W = 21.0, p = 0.015625, n = 6, median gap = +0.365 (all 6/6 positive).**
> **Sign test (robustness): 6/6 positive, p = 0.015625.**

p = 0.015625 = 1/64 is the **smallest achievable p-value for n=6** in an exact
one-sided signed-rank test (all six paired differences ranked positive). The
descriptive "6/6 samples, TP gap +0.37" is therefore confirmed as statistically
significant at α=0.05.

⚠ **Caveat (do not over-claim):** n=6 is small; this is a per-sample paired test on
6 cell-line samples, not a per-variant test. The exact test floors at 1/64, so
p=0.0156 reflects perfect rank consistency, not a large effective sample. The
two-sided exact p = 0.03125. For publication, report n, the exact method, and the
sign-test concordance alongside W/p.

---

## Source provenance (all paths absolute)

| Artifact | Path |
|----------|------|
| Raw per-sample composition (TP rate × side × combo) | `/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv` |
| Inner/Outer same-hap proportions (93-99% claim) | `/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_proportion.tsv` |
| Pre-existing Wilcoxon JSON (2026-04-23, reproduced below) | `/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data/obs18_wilcoxon_B1.json` |
| Generating script | `/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/scripts/obs18_TO_NG2_composition.py` |
| Discovery write-up | `/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/09_TO_sample_af_lohside_ng.md` |
| Memory pointer | `project_loh_constrained_phasing_discovery.md` (2026-04-22) |

The pre-existing `obs18_wilcoxon_B1.json` (dated 2026-04-23) already contained
W=21.0 / p=0.015625. Per project rule "既有 artifact 必先驗證才能引用", I did **not**
trust it blindly — I re-parsed the raw `obs18_NG2_composition_by_sample.tsv` and
recomputed with scipy; the result matches to the bit.

---

## The 6 paired values (parsed from raw TSV)

Pairing: **Inner × NG=2 `same_HP1 (HP1 + HP1-1)` TP rate** vs
**Outer × NG=2 `cross_het (HP1 + HP2-1)` TP rate**, per sample.

| Sample | Inner same_HP1 TP | (n) | Outer cross_het TP | (n) | gap |
|--------|:-----------------:|:---:|:------------------:|:---:|:---:|
| HCC1395        | 0.958904 | 219   | 0.500000 | 2    | **+0.45890** |
| HCC1395_DORADO | 0.938504 | 7220  | 0.553101 | 3305 | **+0.38540** |
| HCC1937        | 0.758586 | 6901  | 0.236364 | 880  | **+0.52222** |
| HCC1954        | 0.428865 | 3402  | 0.084198 | 4050 | **+0.34467** |
| H2009          | 0.931847 | 27086 | 0.882431 | 3785 | **+0.04942** (baseline-saturated) |
| H1437          | 0.919593 | 4614  | 0.688443 | 4439 | **+0.23115** |

- n paired = 6
- gaps = [0.458904, 0.385403, 0.522222, 0.344668, 0.049416, 0.231149]
- median gap = **+0.365035** (rounds to the documented +0.37)
- mean gap = +0.331960; min = +0.049416 (H2009); max = +0.522221 (HCC1937)
- direction: **6/6 positive** (Inner same-hap TP > Outer cross-het TP for every sample)

> Note on HCC1937 vs memory table: the 2026-04-22 memory note rendered HCC1937 as
> "0.76 / 0.24 / +0.52" — consistent with the raw 0.7586 / 0.2364 / +0.5222 here.
> The memory's column header "Inner same_HP1 TP" is what the test actually uses.

### Same-hap proportion claim (93-99%) — verified

From `obs18_NG2_composition_proportion.tsv`, Inner `same_pct`:
HCC1395 93.2%, HCC1395_DORADO 99.0%, HCC1937 98.8%, HCC1954 96.5%,
H2009 98.3%, H1437 97.0% → **93.2%–99.0%**, reproducing the "93-99% same-hap" claim.

---

## Test results (recomputed 2026-05-29, scipy 1.13.1)

### Primary: Wilcoxon signed-rank
```
wilcoxon(inner, outer, alternative='greater', zero_method='wilcox', method='exact')
  -> W = 21.0,  p = 0.015625      (one-sided, Inner > Outer)
wilcoxon(inner, outer, alternative='two-sided', method='exact')
  -> W = 0.0,   p = 0.031250
```
- W=21.0 = sum of positive ranks = 1+2+3+4+5+6 (all six differences positive) =
  maximum possible for n=6 → exact one-sided p = 1/2^6 = **0.015625** (the floor).

### Robustness #1: Sign test (binomial, exact)
```
binomtest(6, 6, 0.5, alternative='greater')  -> p = 0.015625
binomtest(6, 6, 0.5, alternative='two-sided') -> p = 0.031250
```
6/6 concordant positive — agrees with the signed-rank test.

### Robustness #2: alternate same-hap reference (same_HP2 instead of same_HP1)
Using Inner `same_HP2 (HP2 + HP2-1)` TP rate as the reference (instead of same_HP1):
gaps = [0.4081, 0.1976, 0.3267, 0.2083, 0.0340, 0.0982], all 6/6 positive,
**W = 21.0, p = 0.015625** → conclusion is invariant to the HP1/HP2 reference choice.

**Verdict: POSITIVE.** The Inner-same-hap vs Outer-cross-het TP gap is significant
(Wilcoxon one-sided p = 0.0156, sign test p = 0.0156, 6/6 directionally consistent,
median +0.365), robust to the same-hap reference allele.

---

## Reproduce / re-run skeleton (read-only, no heavy compute)

Ready to run as-is — only reads the existing TSV, ~1 s, no pipeline:

```python
import csv
import numpy as np
from scipy.stats import wilcoxon, binomtest

SRC = ("/big7_disk/liaoyoyo2001/InterSubMod/research/"
       "tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv")

rows = list(csv.DictReader(open(SRC), delimiter='\t'))

def tp_rate(sample, side, combo_prefix):
    for r in rows:
        if (r['sample'] == sample and r['loh_side'] == side
                and r['combo'].startswith(combo_prefix)):
            return float(r['tp_rate'])
    return None

SAMPLES = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954", "H2009", "H1437"]
inner = [tp_rate(s, 'Inner', 'same_HP1') for s in SAMPLES]
outer = [tp_rate(s, 'Outer', 'cross_het (HP1') for s in SAMPLES]
gaps = np.array(inner) - np.array(outer)

w = wilcoxon(inner, outer, alternative='greater', zero_method='wilcox', method='exact')
b = binomtest(int((gaps > 0).sum()), int((gaps != 0).sum()), 0.5, alternative='greater')

print("n =", len(gaps))
print("gaps =", [round(x, 6) for x in gaps])
print("median gap =", round(float(np.median(gaps)), 6))
print(f"Wilcoxon (one-sided): W={w.statistic}, p={w.pvalue}")
print(f"Sign test: {(gaps>0).sum()}/{(gaps!=0).sum()} pos, p={b.pvalue}")
# Expected: W=21.0, p=0.015625; sign test 6/6, p=0.015625
```

---

## Notes for the publication / HKU handoff

1. **Report all three together**: Wilcoxon W=21.0 / one-sided p=0.0156 (exact);
   sign-test 6/6 p=0.0156; median paired gap +0.365 (range +0.049 … +0.522).
   State n=6 explicitly and that p hits the exact-test floor (perfect rank order).
2. **Effect size**: report the per-sample paired difference (median +0.365) and its
   bootstrap 95% CI. The pre-existing JSON reports a 1000-resample median-gap 95% CI
   of [0.140, 0.491] (seed 20260423); this is a small-n percentile bootstrap — treat
   as indicative, re-derive with BCa if a referee asks.
3. **H2009 is near-saturated** (Outer cross-het TP already 0.88), so its gap (+0.049)
   is small but still positive — it does not break the 6/6 monotone direction.
4. **Negative control not yet run**: the discovery's planned `--germline-hp-only`
   negative control (expected to collapse the gap) is **not** part of this test and
   would require a pipeline run — deferred (disk 97% full; out of scope for this
   read-only task).
5. **Scope**: this is a 6-sample TO-mode test (the LOH-constrained phasing discovery
   is TO-mode-specific). The companion paired-mode test `B3_wilcoxon_gap_stats.json`
   shows the gap is essentially **absent in paired mode** (paired-gap median ≈ -0.0002,
   Wilcoxon p = 0.578), which is the intended contrast: the signature is a TO-mode
   phenomenon. Worth citing as a built-in specificity control.

## Environment / provenance

- Computed: 2026-05-29 on /big7_disk, read-only; PID 411460 ASM batch not impacted
  (no pipeline/BAM/ISM executed; only TSV reads + in-memory stats).
- Tooling: Python3 / scipy 1.13.1 / numpy 1.23.5.
- Repo commit at time of writing: `274f152` (branch refactor/phase1-safety).
- All numbers parsed from `obs18_NG2_composition_by_sample.tsv` /
  `obs18_NG2_composition_proportion.tsv`; matches pre-existing `obs18_wilcoxon_B1.json`.
