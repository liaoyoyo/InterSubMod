<!--
建立時間: 2026-05-20
agent: main session
status: complete
report_class: cycle 4 H_NEW_2 + H_NEW_4 LOSO synthesis
audience: PI / user / cycle 5+ planner
scope: 2-feature LR + 9-feature drop-caller_af LR, LOSO sample-level CV (5 samples)
parent_plan: bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md v5.0
predecessor: 20260520_loso_sample_level_circularity_revealed
verdict: H_NEW_2 FAIL + H_NEW_4 SANITY VIOLATED (unexpected positive HCC1395 +0.00699)
last_verified: 2026-05-20
-->

# Cycle 4 H_NEW_2 + H_NEW_4 LOSO Findings

> **Two-line bottom line**:
> - **H_NEW_2 FAIL**: 2-feature LR (loh_inner_flag + HPFineF) LOSO mean = -0.00002, 0/5 above +0.002 threshold — robust feature 在 LR framework 不夠強到對 cross-sample 有用
> - **H_NEW_4 SANITY VIOLATED (unexpected positive)**: 9-feature LR (drop caller_af) LOSO HCC1395 = **+0.00699** (vs baseline 10-feature LOSO -0.00012). Drop caller_af 後 4-sample-trained LR 對 HCC1395 反而找到 marginal cross-sample signal — caller_af 確實是 LOSO 災難主因，但同時也 mask 了 small cross-sample signal.

---

## 0. TL;DR — 3-way LOSO Comparison

| Sample | Baseline LOSO (10 feat) | H_NEW_2 LOSO (2 feat) | H_NEW_4 LOSO (9 feat, drop caller_af) |
|---|--:|--:|--:|
| HCC1395 | -0.00012 | -0.00012 | **+0.00699** ★ |
| HCC1937 | +0.00000 | +0.00000 | +0.00000 |
| HCC1954 | -0.00008 | -0.00000 | -0.00008 |
| H1437 | -0.00001 | -0.00000 | -0.00001 |
| H2009 | -0.00001 | +0.00000 | -0.00001 |
| **Mean** | **-0.00004** | **-0.00002** | **+0.00138** |
| n_pos > +0.002 | 0/5 | 0/5 | **1/5** |
| Verdict | DIRECTION_NEGATIVE | **FAIL** (pre-reg) | **SANITY VIOLATED** + marginal positive |

---

## 1. H_NEW_2 結果 (Pre-reg FAIL)

### Pre-reg
- **PASS**: ≥2/5 LOSO ΔF1 > +0.002 AND mean ΔF1 > 0
- **FAIL**: 0/5 above +0.002 OR mean ≤ 0

### Observed
- HCC1395 / HCC1954 / H1437 / H2009 / HCC1937 全部 ΔF1 ≈ 0 (mean = -0.00002)
- 全部 best τ = 0.10 (= keep everything, LR 找不到 useful threshold)
- 0/5 above +0.002

### Pre-reg 預測對照
- 事前 prior 25% PASS — 實際 0% → 符合 prior 偏 conservative 方向

### 含義
即使移除 caller_af 並只保留 cross-sample observation-robust features (loh_inner_flag + HPFineF)，**LR 在 sample-level 仍找不到 useful filter threshold**。Observation-driven 找到的 robust signal 在 LR threshold-based filtering 框架內**不夠強**。

---

## 2. H_NEW_4 結果 (Pre-reg SANITY VIOLATED — Unexpected positive)

### Pre-reg
- 預測 HCC1395 LOSO 仍 ≈ 0 (drop caller_af 只移除 in-distribution overfit 來源)
- Confidence ~80%

### Observed
- HCC1395 LOSO = **+0.00699** (best τ = 0.95 strict filter)
- 其他 4 sample 仍 ≈ 0 (best τ = 0.10)
- Mean = +0.00138 (largely driven by HCC1395 alone)

### Pre-reg 違反 (HARKing 防護)

> **明確聲明**：此結果違反 pre-reg prediction。事前 80% confidence 預期 HCC1395 仍 ≈ 0，實際 +0.00699。Post-hoc 解釋如下，**不是 confirmatory result**。

### Post-hoc 機制解讀

| 觀察 | 解讀 |
|---|---|
| Baseline 10-feature LOSO HCC1395 = -0.00012 | caller_af 在 4-sample combined train 中 coef 異常 (HCC1937 反向 d=-1.41 vs H1437 d=+0.37) → LR 學到「混亂的 caller_af signal」→ 對 HCC1395 (caller_af d=+1.60) 套用無效 |
| H_NEW_4 LOSO HCC1395 = +0.00699 | drop caller_af 後 LR 純粹用 9 features (LOH + Coverage + NG + 5 methyl + chr8) train → 4 sample 的「共同 signal」(主要 LOH + 部分 methylation) 對 HCC1395 some marginal signal |
| Best τ = 0.95 (strict filter) | 與其他 4 sample τ=0.10 完全不同 → 4-sample-trained LR 對 HCC1395 學到「高 confidence 才 keep」rule，犧牲 recall 換 precision marginal gain |

**含義**: caller_af 在 cross-sample LR train 中是「**confusing signal**」(direction-inconsistent)，drop 後其他 features 反而能形成 weak coherent signal — 但只對 HCC1395 (caller F1 = 0.72 + FP density 14%) effective，其他 4 sample 仍 caller-F1-ceiling 卡住。

### 為何不是 cycle 5+ 直接 go-ahead？

| Caveat | 細節 |
|---|---|
| 仍只 1/5 sample positive | HCC1395 alone — 不是 cross-sample robust signal |
| ΔF1 +0.00699 小於 Cohen +0.01 | small effect, not "meaningful improvement" by ClairS-TO paper standard |
| Post-hoc finding | 違反 pre-reg, 未在 confirmatory framework 內 |
| HCC1395 仍是 unique sample | F1 0.72 + FP density 14% — 與其他 sample 性質不同；可能仍是 sample-specific |

---

## 3. 3-Way LOSO Comparison

![3-way comparison](figures/loso_3way_comparison.png)

| Metric | Baseline | H_NEW_2 | H_NEW_4 |
|---|--:|--:|--:|
| Features | 10 (incl. caller_af) | 2 (loh + HPFineF) | 9 (drop caller_af) |
| Mean LOSO ΔF1 | -0.00004 | -0.00002 | **+0.00138** |
| HCC1395 LOSO | -0.00012 | -0.00012 | **+0.00699** |
| n_pos > +0.002 | 0/5 | 0/5 | 1/5 |
| Verdict | DIRECTION_NEGATIVE | FAIL | SANITY VIOLATED + marginal positive on HCC1395 |

---

## 4. Decision Rule Applied

依 plan v5.0 §Decision Rule：

> "**H_NEW_2 FAIL + H_NEW_4 confirms caller_af is root cause** → 確認 LR framework 在 sample-level 完全失敗；pivot 必要"

但這次結果是 **partial confirmation**:
- ✅ caller_af 確實是 LOSO 災難主因 (drop 後 HCC1395 從 -0.00012 → +0.00699)
- ❌ "LR framework 在 sample-level 完全失敗" — partial false (HCC1395 有 marginal cross-sample signal +0.00699)
- ❌ "其他 4 sample 仍 ≈ 0" — TRUE (caller-F1-headroom 仍是 binding constraint)

→ **Revised verdict**: LR framework 在 **caller-F1 < 0.80 + FP density > 10%** samples 上有 *very weak* cross-sample signal (HCC1395 ~+0.007)，**not deployable as universal production filter**.

---

## 5. Implications for Cycle 5+

### 5.1 確認的事

| 事 | 證據 |
|---|---|
| caller_af 在 LR 中是 direction-inconsistent signal | v6_observation AUC 0.20-0.92 / cohen -1.41 ~ +1.60 |
| Drop caller_af 後 LR 找到 marginal HCC1395 signal | H_NEW_4 +0.00699 |
| 4 高 caller F1 samples (H1437/H2009/HCC1954/HCC1937) 受 ceiling 限制 | 全 3 個 LOSO variant ΔF1 ≈ 0 |
| loh_inner_flag + HPFineF 2-feature LR 不夠強 | H_NEW_2 全部 ≈ 0 |

### 5.2 仍 unknown 的事

| Unknown | 後續驗證 |
|---|---|
| HCC1395 +0.00699 是真 cross-sample signal 還是 single-sample artifact? | 需找更多 caller F1 < 0.80 + FP density > 10% samples (cycle 3 Step 2 panel survey) |
| 9-feature LR drop caller_af 在其他 low-F1 sample 上的行為? | 找 panel candidate 後驗 |
| Per-zone LR / RF / interaction terms 是否能改善 H_NEW_2 result? | Cycle 4 Trial A/B/C 待跑 — 但 prior 已下降 |
| Non-LR framework (e.g. zone-stratified rules / boolean filter) 是否有 sample-level signal? | Out-of-scope for cycle 4 |

### 5.3 建議 Cycle 5 入口

依 user 5/20 偏好「先觀察找合適特徵與機制」：

1. **Path A (推薦)**: 接受 LR-based filter direction 是 final-NEGATIVE-for-universal-production；**pivot phase_block_3d 或 thread_d** (per CURRENT_FOCUS plan v2)。LR 結論清楚：受 caller-F1-headroom 限制 + caller_af direction-inconsistent。
2. **Path B**: 若仍要繼續 methyl_filter，建議 H_NEW_3 (chr8-specific gate) 而非 H_NEW_1/Trial A (LR variants) — chr8 zone-specific 是 unique angle 不依賴 LR threshold-finding 能力。
3. **Path C**: 跑 Cycle 3 Step 2 panel survey 找 low-F1 candidates — 驗 HCC1395 +0.00699 是否 generalize 到其他 low-F1 sample。

---

## 6. Files

```
cycle4/loso_validation/
├── scripts/
│   ├── run_loso_cv.py           (baseline 10-feature)
│   ├── run_loso_hnew2.py        (NEW H_NEW_2 2-feature)
│   ├── run_loso_hnew4.py        (NEW H_NEW_4 9-feature drop caller_af)
│   └── loso_3way_comparison.py  (NEW 3-way bar)
├── data/
│   ├── loso_cv_results.tsv      (baseline)
│   ├── loso_hnew2_results.tsv   (NEW)
│   └── loso_hnew4_results.tsv   (NEW)
├── figures/
│   ├── loso_5sample_dF1.png     (baseline)
│   ├── loso_hnew2_5sample.png   (NEW)
│   ├── loso_hnew4_5sample.png   (NEW)
│   └── loso_3way_comparison.png (NEW 3-way bar)
├── intermediate/
│   ├── loso_hnew2_summary.json
│   └── loso_hnew4_summary.json
├── loso_findings.md             (baseline LOSO core)
└── loso_hnew_findings.md        ← this report
```

---

## 7. Reproducibility

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_hnew2.py
python3 research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_hnew4.py
python3 research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/loso_3way_comparison.py
```

Wall clock: ~15 sec each. Deterministic — PRIMARY_SEED=42, lbfgs solver, identical LOSO protocol to baseline (cycle 2 cross_sample_apply.py functions reused).

---

## 8. Pre-Reg Prediction vs Observed (HARKing 防護)

| Hypothesis | Prior PASS | Observed | Match? |
|---|---:|---|---|
| H_NEW_2 ≥2/5 ΔF1 > +0.002 | 25% | 0/5 above threshold | ✅ matches conservative prior |
| H_NEW_4 HCC1395 LOSO ≈0 | 80% | HCC1395 = +0.00699 | ❌ **VIOLATED** |

**H_NEW_4 violation post-hoc analysis** is explicitly labeled NOT confirmatory — needs Path C panel survey to confirm whether +0.00699 generalizes.

---

**End of H_NEW_2 + H_NEW_4 LOSO Findings**
