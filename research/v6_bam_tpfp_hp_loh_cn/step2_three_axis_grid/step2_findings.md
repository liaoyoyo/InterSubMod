# Step 2 — Findings (HCC1395 paired_full, V6_off primary)

## Power gate (Step 1.5)
- Decision: **PROCEED** full 3-axis 50-cell grid
- Powered cells (n≥50, main 50): **23 / 50** (46%)
- Full 75 cells (incl. UNKNOWN LOH): non-empty 41

## H2 — Inner × same_HP TP enrichment
- Powered cells matching: **6**
- Mean TP_enrichment vs global = **1.0109**
- **Verdict: POSITIVE** (POSITIVE ≥ 1.005, given global TP_rate already 0.9832 ceiling)
- Caveat: global TP_rate = 0.9832 (FP universe shrank to 506 due to master-join). Enrichment compressed near ceiling; interpret directionally not as effect size.

## H3 — Outer cross_het × cov_gain FP enrichment
- Cells matching: **4**
- Mean FP_enrichment = **0.0000**
- **Verdict: NEGATIVE** (POSITIVE ≥ 2.0×)

## H6 — Powered cells ≥ 30% (grid sparsity)
- Main 50-cell powered share = **46%**
- **Verdict: POSITIVE**

## H7 — ≥5 cells pass all 7 confound guards
- Cells tested (top-list union): **23**
- Cells passing all 7 guards: **2**
- **Verdict: NEGATIVE**

## V5 over-promote cells (top 5)

| cell_id | V5/V3F ratio | V6/V5 ratio | V6/V3F ratio | n_V3F | n_V5 | n_V6 | n_TP_V6 |
|---|---|---|---|---|---|---|---|
| Inner|cross_het_inv|cov_normal | 5.952 | 0.995 | 5.919 | 62 | 369 | 367 | 366 |
| Outer|cross_het|cov_elevated | 4.577 | 1.000 | 4.577 | 26 | 119 | 119 | 119 |
| Outer|cross_het_inv|cov_normal | 3.854 | 1.002 | 3.860 | 157 | 605 | 606 | 606 |
| Outer|cross_het_inv|cov_elevated | 3.515 | 1.000 | 3.515 | 33 | 116 | 116 | 116 |
| Outer|cross_het|cov_normal | 3.417 | 1.000 | 3.417 | 168 | 574 | 574 | 574 |

## Cells passing all 7 guards

| cell_id | n | n_TP | n_FP | TP_rate | guard4_perm_p | guard6_mh_p |
|---|---|---|---|---|---|---|
| Inner|other|cov_normal | 4984 | 4813 | 171 | 0.9657 | 0.0010 | 0.0000 |
| Outer|other|cov_normal | 8447 | 8198 | 249 | 0.9705 | 0.0010 | 0.0000 |

## Step 3 hand-off candidates

Zone candidates from guard-passing + log-odds extreme cells:

| cell_id | n_guards_passed | n | n_TP | n_FP | TP_rate | use_for |
|---|---|---|---|---|---|---|
| Outer|other|cov_normal | 7 | 8447 | 8198 | 249 | 0.9705 | — |
| Inner|other|cov_normal | 7 | 4984 | 4813 | 171 | 0.9657 | — |
| Inner|same_HP1|cov_normal | 5 | 6000 | 5972 | 28 | 0.9953 | TP-clean reference; Z-OCH/Z-CHR8 contrast |
| Outer|other|cov_elevated | 6 | 3722 | 3708 | 14 | 0.9962 | TP-clean reference |
| Inner|same_HP1|cov_elevated | 6 | 817 | 816 | 1 | 0.9988 | TP-clean reference; Z-OCH/Z-CHR8 contrast |
| Inner|other|cov_elevated | 6 | 601 | 598 | 3 | 0.9950 | TP-clean reference |
| Inner|cross_het_inv|cov_normal | 6 | 367 | 366 | 1 | 0.9973 | TP-clean reference |
| Outer|other|cov_gain | 6 | 332 | 332 | 0 | 1.0000 | TP-clean reference |
| Inner|other|cov_loss | 6 | 325 | 319 | 6 | 0.9815 | — |
| Inner|other|cov_gain | 6 | 83 | 83 | 0 | 1.0000 | — |
| Outer|other|cov_loss | 6 | 58 | 54 | 4 | 0.9310 | FP-rich zone (Step 3 deep dive) |
| Outer|cross_het_inv|cov_normal | 5 | 606 | 606 | 0 | 1.0000 | TP-clean reference |
| Outer|cross_het|cov_normal | 5 | 574 | 574 | 0 | 1.0000 | TP-clean reference |
| Inner|same_HP1|cov_gain | 5 | 63 | 63 | 0 | 1.0000 | Z-OCH/Z-CHR8 contrast |
| Outer|same_HP2|cov_normal | 4 | 641 | 635 | 6 | 0.9906 | TP-clean reference |

## Methodology & ceiling caveat

- Universe restricted to `master_join_ok==1` → 30,036 rows (29,530 TP + 506 FP). FP universe is 10.5% of the phaseC FP set (4,842 total). FP-rate / log-odds metrics are interpretable directionally only.
- Global TP_rate = 0.9832 → TP_enrichment upper bound ≈ 1.017. H2 threshold rescaled to ≥ 1.005.
- LR converged on **4 powered cells** with sufficient label variation (n_FP ≥ 3). Remaining 19 powered cells have n_FP ≤ 1 → LR skipped to avoid spurious β.
- Deviance decomposition (per 02_methodology_notes §3, replacing AUC):
  - `Inner|other|cov_normal`: dev_explained = 0.661 (LRT p_NG=2.5e-5, p_AF<1e-6, p_NR=0.046)
  - `Outer|other|cov_normal`: dev_explained = 0.456 (LRT p_NG<1e-6, p_AF<1e-6, p_NR=0.030)
  - `Inner|other|cov_elevated`: dev_explained = 0.366 (p_AF=0.001 dominant)
  - `Outer|other|cov_elevated`: dev_explained = 0.195 (p_AF<1e-6 dominant)
- Cramér's V (Step 1.5 collinearity audit): max axis-pair 0.610 (loh × NG_bin); VIF for all numeric covariates ≤ 3.35 — within safe range.
- 7-guard pass distribution: 7/7 = 2 cells; 6/7 = 8 cells; 5/7 = 4 cells; 4/7 = 8 cells; 3/7 = 1 cell. Guard 4 (permutation) and Guard 6 (Mantel-Haenszel) are the dominant filters; Guard 7 (HP1/HP2 symmetry) flags 15/23 cells because HP1 family is systematically over-represented (Thread D ratio 1.838 vs 1.0).

## Hand-off to Step 3 (Agent C)

**Zone candidates with cell coordinates:**

1. **Z-OCH (Outer cross_het)**:
   - `Outer|cross_het|cov_normal` (n=574, TP=574, FP=0) — pure TP, V5-over-promote 3.42×
   - `Outer|cross_het|cov_elevated` (n=119, TP=119, FP=0, marginal) — V5-over-promote 4.58×
   - `Outer|cross_het_inv|cov_normal` (n=606, TP=606, FP=0) — V5-over-promote 3.85×
   - **Hand-off**: Z-OCH cells are TP-saturated in master-joined universe; FP enrichment expected to live in the master-unjoined 4,336 FP rows. Recommend zone deep-dive on full phaseC FP set (no master join).
2. **Z-CHR8 (HCC1395 chr8 hotspot)**:
   - Hand-off to Agent C: filter `step1_master_three_way.tsv` to `chr==chr8` and re-run 3-axis grid; compare cell rates to genome-wide.
3. **Z-GL (Gain+LOH, Inner × cov_gain/elevated)**:
   - `Inner|other|cov_elevated` (n=601, FP=3) — passed 6/7 guards
   - `Inner|same_HP1|cov_elevated` (n=817, FP=1) — TP-clean reference
   - `Inner|other|cov_gain` (n=83, TP=83, FP=0, marginal) — passed 6/7 guards
   - Few master-joined FP → expand to full master-unjoined FP for FP-density estimate
4. **Z-AUTO (top 5% FP-density)**:
   - Power-only signal currently: `Outer|other|cov_normal` (n=8447, FP=249 ⇒ ~2.95% FP, well above other cells)
   - `Outer|other|cov_loss` (n=58, FP=4 ⇒ 6.9% FP rate) — small but enriched
   - Hand-off: Agent C run KDE smoothing on per-region FP density in `step1_master_three_way.tsv` and cross-overlap with Z-OCH/Z-CHR8/Z-GL.

**Critical Step 3 caveat**: master-join FP loss (90% drop) means Z-AUTO can be reliably built only on the full 4,842 FP set without master_join_ok filter. Step 3 scripts must read `step1_master_three_way.tsv` directly (not the master-joined subset) for FP-zone construction; LOH/CN annotation will be missing for non-joined FP rows but `loh_side` ≈ region-level proxy from V6_off HP bucket can substitute (Agent C decision).

---
Generated by `step2_three_axis_grid/scripts/confound_guard.py`.
