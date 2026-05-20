<!--
建立時間: 2026-05-19
agent: Coordinator subagent (Phase 2 Cycle 2 Step B3' + B4')
status: complete
report_class: cycle step verdict (Cycle 2 cross-sample n=5)
audience: PI / lab member / Coordinator hand-off (Step #33)
scope: 5 samples {HCC1395 + H1437 + H2009 + HCC1954 + HCC1937} V6_off cross-sample apply
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md v2.0 §Step B3'/B4'
predecessor_cycle: Cycle 2 Step C1 (H_C1_6 PASS), Cycle 1 ⭐3 strong (HCC1395 ΔF1 +0.02236)
verdict: H_C1_5 NEGATIVE (transfer) / MIXED-positive direction (refit) — see §4
last_verified: 2026-05-19
-->

# Phase 2 Cycle 2 Step B3' + B4' — Cross-sample H_C1_5 Verdict (n=5)

> **H_C1_5 verdict** (per pre-registered rules):
> - **Transfer-fit mode** (cycle 1 V6-trained coefs, τ=0.39 fixed): **DIRECTION_NEGATIVE** — 1/5 ΔF1 > 0, 4/5 ΔF1 < 0; Wilcoxon p = 0.1875 (sign-imbalanced but small n=5)
> - **Re-fit mode** (per-sample 5-fold StratKFold OOF, swept τ): **MIXED** — 3/5 strict positive, 0/5 strict negative, 2/5 near-zero (≈1e-7); positive direction trend but only 3/5 same direction
>
> **Bottom line**: The cycle 1 HCC1395-trained filter does NOT transfer to 4 new samples (caller_af weight overfits HCC1395-specific FP/TP distribution). Per-sample re-fit recovers ΔF1≈0 (no harm but also no uplift) on 4 new samples — caller F1 baselines are already 0.84–0.89 (vs HCC1395 0.7166), leaving little headroom for a global LR filter.

---

## 0. TL;DR

| Mode | HCC1395 | H1437 | H2009 | HCC1954 | HCC1937 | Direction | Wilcoxon p | Verdict |
|------|--------:|------:|------:|--------:|--------:|----------:|-----------:|---------|
| Transfer (τ=0.39 fixed) | **+0.02232** | -0.03597 | -0.00450 | **-0.37744** | -0.07068 | 1+ / 4- | 0.1875 | **DIRECTION_NEGATIVE** |
| Transfer (swept τ best) | +0.02246 | -0.00744 | -0.00085 | -0.16972 | -0.02859 | 1+ / 4- | 0.3125 | **DIRECTION_NEGATIVE** |
| Re-fit (per-sample 5-fold CV swept τ) | +0.02236 | -0.0000002 | +0.0000052 | +0.00095 | +0.00761 | 3+ / 0- / 2≈0 | 0.125 | **MIXED** (positive trend) |

**Core observations:**

1. **Transfer-fit catastrophically overfilters on 4 new samples** — HCC1954 ΔF1 = **-0.377** because cycle 1's positive caller_af coefficient (+3.44) makes the LR aggressively keep high-AF variants, but HCC1954 has a different AF distribution (caller F1 = 0.84, mostly high-AF TP) so τ=0.39 drops 11,480 TP (59% of TP) → recall collapses 0.74 → 0.30.
2. **Per-sample re-fit recovers near-zero ΔF1** — Best τ degenerates to 0.10 (H1437/H2009) keeping almost everything, because the 4 new samples have **high caller F1 (0.37–0.89)** and so few TP/FP confusable variants that an LR cannot find a useful τ threshold.
3. **HCC1937 (caller F1=0.37, BRCA1 outlier) re-fit ΔF1 = +0.00761** — Only sample where re-fit yields nominal uplift > +0.005; this is the only sample with enough FP (2,697) and low enough caller F1 for the filter to extract value.
4. **H_C1_5 NEGATIVE** with caveats: the failure mode is **not** that methylation features are unhelpful — it's that the 4 new samples' caller F1 is already near ceiling, leaving no headroom for **any** LR filter (not just one with methylation). Filter generalization is **caller-F1-headroom-bounded**, not methylation-bounded.

**Decision for cycle 2 hand-off**:
- Reject "cycle 1 filter as released production filter" — fails ΔF1 ≥ 0 on 4/5 samples in transfer mode.
- Down-tier cycle 1 result from "⭐3 strong" → "**⭐3 marginal, HCC1395-only**" (still HCC1395-internal-valid, not cross-generalizable).
- Open cycle 3 candidate: **caller-F1-headroom-aware re-design** — only apply filter when caller F1 < 0.80, or build per-sample re-fit infrastructure as production model.

---

## 1. Context

### 1.1 Pre-registered H_C1_5 (cycle 1 deferred to cycle 2)

> "Apply cycle 1 filter (10-feature L2 LR, τ*=0.39) to 4 new samples {H1437, H2009, HCC1954, HCC1937} V6 BAM with significance.csv. Wilcoxon paired signed-rank vs 0 on per-sample ΔF1.
> - POSITIVE: ≥4/5 ΔF1 > 0 AND Wilcoxon p < 0.0625 (requires 5/5 same direction at n=5 exact)
> - DIRECTION: ≥4/5 same direction but p ≥ 0.0625 (n=5 power-limited)
> - NEGATIVE: ≤2/5 same direction"

n=5 Wilcoxon exact min p-value is **0.0625** — achievable only if **all 5 same direction**.

### 1.2 Cycle 1 ⭐3 strong (single-sample HCC1395)

- Cycle 1 trained 10-feature global LR on HCC1395 V6 BAM, 35,332 region-rows
- Pre-reg ΔF1 ≥ +0.01 (Cohen small) → achieved **+0.02236** = 2.24× threshold
- caller_af (+3.44) / LOH_inner (+1.46) / Cov_imp (+1.27) / NG (+1.07) / HPFineF (+0.75) top-5 by |coef|
- Methylation features = 5th-rank covariate (HPFineF +0.75) not dominant

### 1.3 4 new samples — ISM with significance.csv

Step B1' (predecessor task) re-ran V6 ISM with significance on 4 samples:

| Sample | V6_off TP rows | V6_off FP rows | Caller F1 baseline | Truth set |
|--------|---:|---:|---:|---|
| H1437 | 70,134 | 713 | 0.8670 | Orthogonal |
| H2009 | 135,353 | 1,334 | 0.8863 | Orthogonal |
| HCC1954 | 19,446 | 685 | 0.8385 | Orthogonal |
| HCC1937 | 13,906 | 2,475 | 0.3692 | Orthogonal (BRCA1 outlier) |

Per-sample F1 baselines: HCC1395 F1 = 0.7166 from **ClairS-TO v0.3.0 (ssrs model) @ 0.93 purity vs SEQC2 truth** (cycle 1 constant; source: `InterSubMod/research/paired_priority_bug_audit/09_V6_caller_F1_verification.md` Table L70-72 — TP=28,509 / FP=11,606 / FN=10,938；canonical pipeline = ClairS-TO ssrs → LongPhase-TO V3F/V5/V6 → ISM; tumor-only mode 全程無 normal BAM). 4 新樣本 caller F1 from orthogonal truth-set benchmarks (per-sample ClairS-TO ssrs baseline)。**重要**：先前 cycle 2 引用「20260216 paired pipeline baseline」是 framing 錯誤（已 2026-05-20 patch；參見 cycle2_NEGATIVE_postmortem.md §5 "Caller Mode Framing Inheritance Error"）。ClairS-TO 沒有 pileup/full_alignment model 區分（只 ss/ssrs），先前命名「to_pileup」雙重誤命，已修為 canonical 「clairs_to_ssrs」對齊 ClairS-TO 官方 model name。

### 1.4 Step B2' — Augmented master TSV per sample

`build_4sample_master.py` joined V6 significance.csv methylation features (13 + 2 derived) onto `step4_cross_sample_extension/per_sample_master/{sample}_v6_master.tsv` chassis on (chr, pos). Output: `cycle2/data/{sample}_master_augmented.tsv` × 4 files.

| Sample | Rows | n_TP | n_FP | Cycle1 essential min non-null % | Methylation min non-null % |
|---|---:|---:|---:|---:|---:|
| H1437 | 70,964 | 70,191 | 773 | 95.08% | 79.47% |
| H2009 | 136,701 | 135,359 | 1,342 | 97.29% | 79.02% |
| HCC1954 | 20,136 | 19,449 | 687 | 89.08% | 50.93% |
| HCC1937 | 16,607 | 13,910 | 2,697 | 75.80% | 34.33% |

NaN methylation features → impute by cycle 1 medians (Strategy B, MNAR-justified) for transfer mode, by per-sample medians for re-fit mode.

---

## 2. Methods — Step B3' Cross-Sample Apply

### 2.1 Two modes per sample (n=5)

| Mode | Coefficients | Scaler | Imputation | τ |
|---|---|---|---|---|
| **transfer_fixed** | cycle 1 per-fold coefs (5-fold mean ensemble) | cycle 1 per-fold scalers | cycle 1 stored medians (V6) | τ*=0.39 (cycle 1 optimum) |
| **transfer_sweep** | cycle 1 per-fold coefs ensemble | cycle 1 scalers | cycle 1 medians | τ swept 0.10–0.95 (per-sample best ΔF1) |
| **refit** | re-fit L2 LR (C=1.0, 5-fold StratKFold, lbfgs, seed=42) | per-sample fold scalers | per-sample medians | τ swept (per-sample best ΔF1) |

### 2.2 Metric calculation (region-level ΔF1)

`compute_metrics()` uses **region-level** TP/FP counts as TP_TOTAL/FP_TOTAL (because the LR predicts region-level), but caller_F1 baseline + FN_caller are at **variant-level** (consistent with cycle 1's `CALLER_F1_BASELINE = 0.7166` + `FN_CALLER = 19288` constants).

For HCC1395 we hard-code cycle 1 region-level counts (TP=30,490, FP=4,842, FN=19,288) to bit-exactly reproduce cycle 1's published ΔF1 = +0.02236. For the 4 new samples:

- TP_total / FP_total = region-level master row counts
- FN_caller = reverse-solved from variant-level F1: `FN = 2*TP_variant/F1 − 2*TP_variant − FP_variant`

Reverse-solved FN_caller:

| Sample | TP_variant | FP_variant | F1 | FN_caller (solved) |
|---|---:|---:|---:|---:|
| H1437 | 70,191 | 773 | 0.8670 | 20,762 |
| H2009 | 135,359 | 1,342 | 0.8863 | 33,387 |
| HCC1954 | 19,449 | 687 | 0.8385 | 6,805 |
| HCC1937 | 13,910 | 2,697 | 0.3692 | 44,835 |

ΔF1 = F1_post − caller_F1, where F1_post is computed from filter-kept region-level TP/FP plus variant-level FN_caller.

### 2.3 Output

- `cycle2/data/cycle2_cross_sample_delta_f1.tsv` (5 samples × 3 modes × 18 cols)
- `cycle2/figures/cycle2_b3_per_sample_delta_f1.png`
- `cycle2/intermediate/cycle2_b3_summary.json` (full per-fold coefs / feature importance)
- `cycle2/intermediate/cross_sample_apply.log`

---

## 3. Results — Step B3' Cross-Sample

### 3.1 Transfer mode (cycle 1 V6-trained coefs, fixed τ=0.39)

| Sample | τ | TP_kept | FP_kept | TP_removed | FP_removed | Precision_post | Recall_post | F1_post | ΔF1 |
|---|--:|---:|---:|---:|---:|--:|--:|---:|---:|
| HCC1395 | 0.39 | 30,015 | 1,447 | 475 | 3,395 | 0.9540 | 0.6030 | 0.7389 | **+0.02232** |
| H1437 | 0.39 | 64,931 | 383 | 5,260 | 390 | 0.9941 | 0.7139 | 0.8310 | **-0.03597** |
| H2009 | 0.39 | 133,959 | 1,125 | 1,400 | 217 | 0.9917 | 0.7938 | 0.8818 | **-0.00450** |
| HCC1954 | 0.39 | 7,969 | 345 | **11,480** | 342 | 0.9585 | 0.3035 | 0.4611 | **-0.37744** |
| HCC1937 | 0.39 | 10,717 | 2,340 | 3,193 | 357 | 0.8208 | 0.1824 | 0.2985 | **-0.07068** |

**Transfer mode catastrophe on HCC1954**: τ=0.39 drops 11,480/19,449 = 59% of TP. Mechanism: cycle 1's caller_af coefficient is +3.44 trained on HCC1395 (TP centered at AF≈0.45, FP at AF≈0.10). HCC1954 has TP centered at a different AF point; standardization mean ≈0.45 from cycle 1 scaler shifts HCC1954 TP into the "filter out" half-plane.

### 3.2 Transfer mode with swept τ (best per-sample threshold)

| Sample | best τ | ΔF1 (swept best) |
|---|--:|--:|
| HCC1395 | 0.38 | +0.02246 |
| H1437 | 0.10 | -0.00744 (still negative even at most permissive τ) |
| H2009 | 0.10 | -0.00085 |
| HCC1954 | 0.10 | -0.16972 (still catastrophic) |
| HCC1937 | 0.10 | -0.02859 |

Even at τ=0.10 (most permissive), HCC1954 still loses 6,029 TP — the transfer coefs themselves over-confidently flag many HCC1954 TP as FP-like.

### 3.3 Re-fit mode (per-sample 5-fold StratKFold OOF, swept τ)

| Sample | best τ | TP_kept | FP_kept | TP_removed | FP_removed | F1_post | ΔF1 |
|---|--:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 0.39 | 30,015 | 1,443 | 475 | 3,399 | 0.7390 | **+0.02236** |
| H1437 | 0.10 | 70,191 | 773 | 0 | 0 | 0.8670 | **-2.4e-7** (~0) |
| H2009 | 0.61 | 135,349 | 1,328 | 10 | 14 | 0.8863 | **+5.2e-6** (~0) |
| HCC1954 | 0.71 | 19,401 | 568 | 48 | 119 | 0.8395 | **+0.00095** |
| HCC1937 | 0.24 | 13,830 | 830 | 80 | 1,867 | 0.3768 | **+0.00761** |

**Re-fit mode observation**: For H1437, the optimal τ degenerates to 0.10 (keep everything) because per-sample LR cannot find a useful threshold — at any threshold > 0.10, more TP than FP get removed (caller F1 = 0.87, FP/TP ratio = 1.1%, FP density too sparse to gain from filtering).

H2009 best τ=0.61 keeps 99.99% of variants → essentially no filtering.

Only HCC1937 yields meaningful re-fit uplift (+0.00761), driven by its low caller F1 (0.37) — its FP rate (16% FP / 84% TP) is high enough for an LR to find FP-rich regions to drop.

---

## 4. Step B4' — H_C1_5 Wilcoxon Verdict (n=5)

### 4.1 Per-mode Wilcoxon

| Mode | n_pos | n_neg | n_zero | Direction | Wilcoxon p (ΔF1) | Wilcoxon p (ΔP) | Wilcoxon p (ΔR) | Verdict |
|------|------:|------:|-------:|:---------:|-----:|-----:|-----:|---------|
| transfer_fixed | 1 | 4 | 0 | 1+ / 4- | 0.1875 | 1.0 | **0.0625** | **DIRECTION_NEGATIVE** |
| transfer_sweep | 1 | 4 | 0 | 1+ / 4- | 0.3125 | 0.125 | **0.0625** | **DIRECTION_NEGATIVE** |
| refit | 3 | 0 | 2 | 3+ / 0- / 2≈0 | 0.125 | 0.0679 | 0.0679 | **MIXED** (positive trend) |

(near-zero threshold = |ΔF1| < 1e-5; refit H1437 = -2.4e-7 and H2009 = +5.2e-6 fall in "zero" bucket.)

**Direction concordance interpretation:**
- Transfer modes: 4/5 NEGATIVE direction (only HCC1395 positive, the training sample). Sign-imbalance enough to reach DIRECTION verdict but n=5 lacks power for p < 0.0625 strict majority (would need 5/5 same direction).
- Re-fit mode: 3/5 strict POSITIVE, 2/5 near-zero (no negative). Positive trend but cannot reach DIRECTION_POSITIVE rule (needs ≥4 strict positives).

### 4.2 Wilcoxon p-value scope notes

- n=5 exact Wilcoxon (with no zeros) min p = 1/2⁵ × 2 = 0.0625 (two-sided), achievable only if all 5 same direction.
- Refit mode has 2 zeros (H1437 / H2009); SciPy switches to normal approximation → p = 0.125 (but n=5 too small for asymptotic accuracy).
- For DIRECTION_NEGATIVE verdicts the Wilcoxon p on **ΔRecall** = 0.0625 (5/5 ΔRecall < 0 in transfer modes), confirming that the filter uniformly reduces recall across all 5 samples — which is the proximal mechanism of the F1 collapse.

### 4.3 H_C1_5 final verdict

**Transfer-fit mode**: **DIRECTION_NEGATIVE** (1/5 ΔF1 > 0, 4/5 ΔF1 < 0, Wilcoxon p=0.1875 not significant at α=0.0625 due to power limit but direction unambiguous). The cycle 1 HCC1395-trained filter is **NOT** cross-sample generalizable as a transferable rule.

**Re-fit mode**: **MIXED** (3/5 strict POS, 2/5 near-zero). Positive trend exists but ΔF1 effect sizes on the 4 new samples are negligible (3 of 4 < +0.005 = below Cohen small). Only HCC1937 (caller F1=0.37) yields meaningful uplift (+0.00761).

---

## 5. Mechanism — Why does the filter fail to generalize?

### 5.1 Caller F1 headroom hypothesis

The cycle 1 filter was trained on HCC1395 with caller F1 = **0.7166**. The 4 new samples have caller F1:

| Sample | Caller F1 | Headroom (1 − F1) | TP density | FP density | TP/FP ratio |
|---|--:|--:|--:|--:|--:|
| HCC1395 | 0.7166 | 0.2834 | 86% | 14% | 6.3 |
| H1437 | 0.8670 | 0.1330 | 98.9% | 1.1% | 90.8 |
| H2009 | 0.8863 | 0.1137 | 99.0% | 1.0% | 100.9 |
| HCC1954 | 0.8385 | 0.1615 | 96.6% | 3.4% | 28.3 |
| HCC1937 | 0.3692 | 0.6308 | 83.8% | 16.2% | 5.2 |

For samples with caller F1 > 0.83, **FP density is < 4%** — the filter's job (find FP-rich regions to drop) becomes nearly trivial: there's almost no FP signal to filter. The "best" filter just keeps everything (τ ≈ 0.10 → no rows dropped).

For HCC1937 (caller F1=0.37, FP density 16%), the filter has work to do — and re-fit ΔF1 = +0.00761, the only meaningful uplift.

For HCC1395 (caller F1=0.72, FP density 14%, TP/FP=6.3), the filter has the most balanced training signal — and yields the strongest uplift (+0.02236).

### 5.2 Transfer mode failure root cause

The cycle 1 LR coefficients are baked to HCC1395's caller_af distribution:

- HCC1395 TP caller_af mean ≈ 0.45, FP caller_af mean ≈ 0.10
- LR learns `caller_af` coef = +3.44 → "high AF → more likely TP"
- Applied to HCC1954 with mean caller_af ≈ 0.35 (different distribution), the standardized scaler subtracts 0.45 mean → 0.35-0.45 = -0.1/0.28 = -0.36 standard units → many TP get z-score < 0 → P(TP) < 0.5 → flagged as FP-like at τ=0.39.

This is classic **distributional shift** failure mode — the LR is not robust to the new samples' caller_af mean.

### 5.3 What the failure does NOT mean

- It does NOT mean methylation features are unhelpful per se. The methylation coef (HPFineF = +0.75) is the 5th-rank weight; the failure mechanism is dominated by caller_af + LOH + Coverage (rank 1-3) overfitting.
- It does NOT mean cycle 1 is wrong on HCC1395 — H_C1_2/H_C1_3 still pass on HCC1395.
- It does NOT mean the 13 methylation features have no signal — they're masked here by the high TP density on 3/4 new samples (no room for any filter to improve).

---

## 6. HCC1937 Outlier Analysis

HCC1937 is the only sample with caller F1 = 0.37 — a known BRCA1 driver outlier. Step B2' showed marker rate 0.817 (edge case).

- Re-fit best τ=0.24, F1_post=0.3768 vs caller F1=0.3692 → ΔF1 = **+0.00761**
- TP_removed=80 (0.58% of 13,910 TP), FP_removed=1,867 (69% of 2,697 FP)
- Mechanism: HCC1937 has 2,697 FP vs only 13,910 TP — FP density 16% (similar to HCC1395), so the filter can identify FP-rich regions

HCC1937 also passes methylation non-null sanity at lower rates (34% NME_imbalance) but cycle 1 essential features stay at 75%+ non-null. The HCC1937 result is **consistent with cycle 1 mechanism** — methylation augmentation contributes when FP density is non-trivial.

---

## 7. HCC1954 Caller Ceiling Outlier Analysis

HCC1954 has caller F1 = 0.8385 (third-highest among 5) but its TP rate at the region-level is **0.084 (Outer cross_het, Thread D §B2 finding)** — a known **caller FP background ceiling** outlier (see `InterSubMod/docs/reports/validated/2026/04/20260426_HCC1954_Outlier_CallerCeiling_CasePanel_01.md`).

Transfer-fit ΔF1 = -0.377 (catastrophic): not because methylation features are bad on HCC1954, but because **cycle 1's high-AF prior** (caller_af coef +3.44) collapses HCC1954's broader TP AF distribution.

Re-fit ΔF1 = +0.00095: the re-fit per-sample LR finds best τ=0.71 (very conservative — keep nearly everything because filtering hurts). HCC1954 is **bounded by caller ceiling**, not by methylation features.

---

## 8. Files Inventory

### 8.1 Scripts

- `cycle2/scripts/build_4sample_master.py` — Step B2' per-sample augmented master TSV builder
- `cycle2/scripts/cross_sample_apply.py` — Step B3' 5-sample × 3 modes filter apply
- `cycle2/scripts/wilcoxon_verdict.py` — Step B4' Wilcoxon paired signed-rank n=5 verdict

### 8.2 Data

- `cycle2/data/{H1437,H2009,HCC1954,HCC1937}_master_augmented.tsv` — Step B2' output (4 files, 5–40 MB each)
- `cycle2/data/cycle2_cross_sample_delta_f1.tsv` — Step B3' 5 samples × 3 modes × 18 cols
- `cycle2/data/cycle2_b4_wilcoxon_summary.tsv` — Step B4' per-mode Wilcoxon stats + verdict

### 8.3 Figures

- `cycle2/figures/cycle2_b3_per_sample_delta_f1.png` — 5-sample × 3-mode grouped bar chart with H_C1_3 +0.01 and cycle 1 +0.02236 reference lines
- `cycle2/figures/cycle2_b4_wilcoxon.png` — 3-panel per-mode ΔF1 bar with per-sample annotation + Wilcoxon p + verdict box

### 8.4 Intermediate

- `cycle2/intermediate/build_4sample_master.log`
- `cycle2/intermediate/step_b2_prime_sanity.md`
- `cycle2/intermediate/cross_sample_apply.log`
- `cycle2/intermediate/cycle2_b3_summary.json`
- `cycle2/intermediate/cycle2_b4_summary.json`

---

## 9. Coordinator hand-off (Step #33)

### 9.1 Tier proposal for cycle 1

- **Original cycle 1 proposal: ⭐3 strong** (HCC1395 single-sample verified)
- **Cycle 2 B3'/B4' result**: cycle 1 filter does NOT cross-sample generalize (transfer DIRECTION_NEGATIVE)
- **Revised proposal: ⭐3 marginal (HCC1395-only)** — cycle 1 result remains internally valid for HCC1395 (H_C1_2 + H_C1_3 pass, +0.02236 = 2.24× Cohen small), but H_C1_5 verdict NEGATIVE breaks the ⭐4 cross-sample upgrade path.

### 9.2 Cycle 3 candidates (if pursued)

1. **Caller-F1-headroom-aware filter**: Only apply LR filter when caller F1 < 0.80 (HCC1395, HCC1937 qualify); else skip. Then re-run H_C1_5 on the qualifying subset (n=2 — too small for Wilcoxon).
2. **Per-sample re-fit as production model**: Provide template LR but require re-fit per sample. Pre-reg: "per-sample re-fit mean ΔF1 > 0 across n=5" → currently mean=+0.00619 but median=+0.00095. Effect size weak.
3. **Add caller_af distributional adaptation**: Robust scaler / IQR-based, or per-sample mean/scale matching before applying transfer coefs. Test on HCC1395 hold-out + 4 new samples.
4. **HCC1937-style low-F1 outlier deep-dive**: HCC1937 ΔF1 = +0.00761 in re-fit suggests low-F1 samples are the natural target. Search for more low-F1 samples (panel) to reach n=4+ for Wilcoxon.

### 9.3 Decision required from Coordinator

- ☐ Accept ⭐3 marginal HCC1395-only down-tier for cycle 1
- ☐ Decide cycle 3 direction (one of 4 above, or abandon methylation filter direction)
- ☐ Update `evidence_ledger.jsonl` with cycle 2 NEGATIVE entry referencing this findings

---

## 10. Reproducibility

Reproduce all 3 steps from scratch:

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/methyl_augmented_filter_phase2/cycle2/scripts/build_4sample_master.py    # ~18 s
python3 research/methyl_augmented_filter_phase2/cycle2/scripts/cross_sample_apply.py     # ~11 s
python3 research/methyl_augmented_filter_phase2/cycle2/scripts/wilcoxon_verdict.py       # ~2 s
```

Wall clock: ~30 s total. Cycle 1 V6 re-fit ΔF1 reproduces exactly +0.02236 (bit-exact, matches `cycle1_track_a_filter.json.expected_delta_F1`).

Seed: PRIMARY_SEED = 42 (matches cycle 1). sklearn `lbfgs` solver deterministic.
