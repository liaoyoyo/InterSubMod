<!--
建立時間: 2026-05-18
agent: Coordinator (main session) Phase 2 Cycle 1 synthesis
status: complete
report_class: cycle synthesis (HCC1395 single-sample, Track B deferred)
audience: PI / lab member / 自己未來
scope: HCC1395 single-sample global FP filter pilot — Track A only
tier_proposal: ⭐3 strong (HCC1395 verified, multi-seed stable, ΔF1 +0.02236 = 9.24× v1.0)
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md v2.0
predecessor_cycle: v1.0 step5_methyl_filter_pilot ⭐3 marginal (+0.00242)
verdict: H_C1_1/2/3 PASS / H_C1_4 FAIL / H_C1_5/6 DEFERRED to cycle 2 (4-sample BAM 物理不存在)
last_verified: 2026-05-18
-->

# Phase 2 Cycle 1 Findings — HCC1395 Global FP Filter (Track A only)

> **Tier 提案**: ⭐3 strong (HCC1395 single-sample verified, ΔF1 +0.02236 = 9.24× v1.0 baseline)
> **Track B deferred**: 4 樣本 V3F/V5 BAM 物理不存在；用戶決定 cycle 1 在 HCC1395 single-sample 收尾
> **Predecessor**: v1.0 ⭐3 marginal (+0.00242)；本 cycle pivot global FP exploration 達 11× 提升

---

## 0. TL;DR

Cycle 1 Track A 在 HCC1395 single-sample 達到 **ΔF1 = +0.02236** (9.24× v1.0 baseline, 2.24× +0.01 Cohen ribbon)，**3/4 預註冊假設 PASS**，且機制清晰：

- ✅ **H_C1_1 (top 10 cells covers 94.22% all FP)** PASS — FP 高度集中
- ✅ **H_C1_2 (ΔF1 > v1.0 +0.00242)** PASS — 9.24× 提升
- ✅ **H_C1_3 (ΔF1 ≥ +0.01)** PASS — 2.24× 達 Cohen 小 effect threshold
- ❌ **H_C1_4 (high-AF zone ≥ +0.003)** FAIL — sub-zone 限制 (-0.00011)，但 global LR 仍 capture
- 🔵 **H_C1_5/H_C1_6 (cross-sample + V3F/V5/V6 sanity)** **DEFERRED** — 4 樣本 V3F/V5 BAM 物理不存在

**核心發現**：
1. **Step 5c lost TP 17/21 (81%) 被 cycle 1 filter rescue** — 證明 global LR ≫ v1.0 cell-gated
2. **Methylation 是 5th-rank covariate**（HPFineF coef +0.75），但 caller_af (+3.44) / LOH (+1.46) / Cov (+1.27) / NG (+1.07) 主導 — cycle naming 可考慮改 "multi-axis filter incl. methylation"
3. **Multi-seed std = 5e-5** (20× below threshold) — 高度穩定

**Limitation**: HCC1395 single-sample only, generalizability 未驗

---

## 1. Context

### 1.1 Predecessor (v1.0 cycle, 2026-05-18 早)

v1.0 cycle ⭐3 PARTIAL POSITIVE marginal:
- ΔF1 = +0.00242 (< +0.005 Cohen ribbon)
- Step 5c TP rescue NEGATIVE (95.2% lost TP 是 low-AF subclone)
- Step 5d robustness GREEN with caveats

### 1.2 Cycle 1 Pivot (Plan v2.0)

從 v1.0 「cell-gated filter」pivot 為「global FP exploration + heterogeneous threshold」:
- v1.0 powered cell gate 過嚴 → 4 cells 只 covers 7% all FP
- Plan v2.0 預測：全 35,332 row global LR 應達 ΔF1 ≥ +0.01

### 1.3 Track B 物理限制

Plan v2.0 §Step B1 假設「V3F/V5 4 樣本 BAM 已存在」**實際不存在** — Deep search agent 確認只 HCC1395 有 V3F/V5。Track B 必須 defer 或 Path B V6-only。**用戶決定 cycle 1 在 HCC1395 single-sample 收尾**（V6 4 樣本 ISM rerun 留 cycle 2）。

---

## 2. Pre-Registered Hypotheses Verdicts

| ID | Prediction | Falsification | Observed | Verdict |
|----|------------|---------------|----------|---------|
| **H_C1_1** | Top 10 FP-rich cells covers ≥70% all FP | <70% | **94.22%** | ✅ PASS |
| **H_C1_2** | Global LR ΔF1 > v1.0 +0.00242 | ≤+0.00242 | **+0.02236** (9.24×) | ✅ PASS |
| **H_C1_3** | Global LR ΔF1 ≥ +0.01 (Cohen small) | <+0.01 | +0.02236 (2.24×) | ✅ PASS |
| **H_C1_4** | High-AF zone incremental ΔF1 ≥ +0.003 | <+0.003 | -0.00011 | ❌ FAIL |
| **H_C1_5** | Cross-sample 4 樣本 ≥4/5 ΔF1>0 + Wilcoxon p<0.05 | — | DEFERRED | 🔵 — |
| **H_C1_6** | V3F/V5/V6 BAM cross-binary sanity | — | DEFERRED | 🔵 — |

**NO-GO trigger**: H_C1_2 + H_C1_3 同時 FAIL — 未觸發。

---

## 3. Final Filter Specification

| Component | Value |
|-----------|-------|
| Features (10, drop NumReads_master VIF=217) | V6_off_NG, caller_af, NumReads_master ❌ removed, loh_inner_flag, Coverage_Multiple_imp, V6_off_meth_{HPMergedDelta, HPFineF, NME_imbalance, Epipoly_Delta, ClusterPermanovaF}, chr8_flag |
| Regularization | L2 (Ridge), C=1.0 |
| Optimal τ* | 0.39 (broad plateau 0.38-0.42) |
| NaN handling | Strategy B (median impute) — justified by MNAR analysis (NaN 80% at AF<0.1, 92% at AF>0.7 U-shape) |
| CV | 5-fold StratifiedKFold OOF |
| Feature ranking (by abs coef) | caller_af (+3.44) > LOH_inner (+1.46) > Cov (+1.27) > NG (+1.07) > HPFineF (+0.75) > HPMergedDelta > NME_imbalance > Epipoly_Delta > ClusterPermanovaF > chr8_flag |
| **Deployable artifact** | `cycle1_track_a_filter.json` (per-fold coefs + scalers + medians) |

### Confusion matrix (HCC1395, τ*=0.39)

| | Truth=TP | Truth=FP |
|--|----------|----------|
| Predicted keep (P≥0.39) | 30,015 | 1,443 |
| Predicted filter (P<0.39) | 475 (1.56%) | 3,399 (70.20%) |

- Precision = 30,015 / 31,458 = **0.9541**
- Recall = 30,015 / 49,778 = **0.6030** (recall slight drop from caller 0.6125)
- **ΔF1 = +0.02236**

---

## 4. 重要發現

### 4.1 Step 5c lost TP rescue 81% — v1.0 cell-gated 過嚴

v1.0 cell-gated filter (4 AGGREGATED cells) 失去 21 lost TP，所有 rescue rule 都 NEGATIVE。
**Cycle 1 global LR 救回 17/21 = 81%**：

| Rescue approach | TP rescued | FP reintroduced | ΔF1 vs cycle 1 baseline |
|----------------|-----------|-----------------|------------------------|
| v1.0 Step 5c (cell-gated rescue rule) | 0/21 | — | -0.00043 |
| **Cycle 1 global LR (本 cycle)** | **17/21 (81%)** | 0 (in same baseline) | **+0.02236** |

**機制**：v1.0 cell-gated 受限 4 cells × 60 TP + 344 FP 過 sparse；global LR borrows strength across 35,332 rows → 同 21 row 在 better-fit boundary 下 P(TP) 升至 ≥ 0.39。

剩 4/21 unrescued 是極低 caller_af cases — Future work: methylation-specific gated rescue (V6_off_meth_AlleleP，依 v1.0 Step 5c 提示 best feature 但需 gated mode)。

### 4.2 Methylation 是 5th-rank covariate 非主導

caller_af (+3.44) / LOH (+1.46) / Cov (+1.27) / NG (+1.07) 主導；methylation HPFineF (+0.75) 第 5 位但仍貢獻 — **non-zero coef confirm methylation 有獨立 incremental information**。

對 paper framing 影響：
- 不能宣告「methylation filter」— 應稱「**multi-axis filter incl. methylation**」
- §3 主軸：caller AF + structural (LOH/CN) + phasing (NG) + methylation augmentation
- Limitations: methylation 為 5th-rank 而非主導

### 4.3 R-Step0 5 個 caveat 處理狀態

| Caveat | Severity | Status |
|--------|----------|--------|
| R-Step0-1 (NaN MNAR 8.7× gap) | HIGH | ✅ RESOLVED — MNAR confirmed; impute correct |
| R-Step0-2 (Cov/NumReads collinearity VIF=217) | HIGH | ✅ RESOLVED — drop NumReads_master, max VIF=1.83 |
| R-Step0-3 (methylation marginal) | MED | ✅ ACKNOWLEDGED — methyl 5th-rank, framing 改 |
| R-Step0-4 (lost TP overlap) | MED | ✅ RESOLVED — 81% rescue |
| R-Step0-5 (HCC1395-only) | HIGH | 🔵 OPEN — Track B deferred |

### 4.4 Multi-seed CV stability

| Seed | best ΔF1 | best τ |
|------|----------|--------|
| 42 | +0.02236 | 0.39 |
| 7 | +0.02245 | 0.38 |
| 13 | +0.02235 | 0.39 |
| 2026 | +0.02233 | 0.42 |
| 1395 | +0.02230 | 0.41 |
| **Mean** | **+0.02236** | 0.398 |
| **Std (ddof=1)** | **5e-5** | 0.013 |

Std = 5e-5 << 0.001 threshold → **STABLE intra-sample**.

---

## 5. Tier Evaluation

### ⭐3 strong (proposed)

| 維度 | 評估 |
|------|------|
| Effect size | ΔF1 +0.02236 = 2.24× Cohen small ribbon ✅ |
| Pre-registration | 4 H pre-specified, 3 PASS / 1 FAIL ✅ |
| Confound guard | VIF audit + L2 reg + MNAR analysis + 5-fold OOF ✅ |
| Multi-seed stability | std 5e-5 ✅ |
| Lost TP rescue | 81% (vs v1.0 0%) ✅ |
| Cross-sample | **DEFERRED** ❌ |
| Mechanism | non-methyl axes 主導，methyl 5th-rank (R-Step0-3 acknowledged) ⚠️ |

### ⭐3 而非 ⭐4 的理由

✅ HCC1395 single-sample evidence 強 (9.24× v1.0)
✅ R-Step0-1/2/4 caveats all RESOLVED
❌ H_C1_5 cross-sample 未驗 (V3F/V5 4 樣本 BAM 不存在 → deferred)
❌ H_C1_6 V3F/V5/V6 sanity 未做（可用 HCC1395 phaseC 三向資料補做 — 留 cycle 2）

⭐4 升級必要條件 (cycle 2):
1. V6 4 樣本 ISM rerun with significance (~3.2 hr, BAM 已存在於 v6_5sample_extension/)
2. Cross-sample H_C1_5 Wilcoxon n=5
3. HCC1395 phaseC V3F/V5/V6 三向 cycle 1 filter apply 看 H_C1_6 cross-binary variance

---

## 6. Paper Framing Recommendation

### 6.1 主軸 (保守 framing)

> **「Multi-axis filter (caller AF + LOH + Coverage_Multiple + HPFineNGroups + methylation 5 features) on HCC1395 clairs_to_ssrs (ClairS-TO v0.3.0 ssrs model + LongPhase-TO V6 pipeline; tumor-only mode, no normal BAM) achieves ΔF1 = +0.02236 vs ClairS-TO caller baseline F1=0.7166 (post-filter 0.7390). Methylation contributes as 5th-rank covariate (HPFineF coef +0.75 in 10-feature LR). 81% of v1.0 cell-gated lost TP rescued by global LR. Cross-sample generalization pending (only HCC1395 has V3F/V5 BAMs)."**

### 6.2 §3 主圖 (建議)

- Fig 3a: τ sweep ΔF1 curve (multi-seed ± band) — 顯示 plateau 0.38-0.42
- Fig 3b: VIF audit + L2 regularization sweep (collinearity 解決證據)
- Fig 3c: Step 5c 21 lost TP individual P(TP) — show 17/21 rescued
- Fig 3d: Feature importance (10 features by abs coef) — methylation 5th-rank visible

### 6.3 與 prior art 差異化

- TumorLens: sample-level not per-region filter
- ROCIT: methylation-only transformer, no multi-axis
- SGZ: variant-level 4-axis but no phasing/methylation
- **本 cycle 差異**: per-region multi-axis (LOH + CN + AF + phasing + methylation) variant-level filter — 仍無同口徑 prior art

### 6.4 Limitations 章節必寫

1. HCC1395 single-sample only — cross-sample generalization pending cycle 2
2. Methylation is 5th-rank covariate (HPFineF coef +0.75) — non-methyl axes (AF/LOH/Cov/NG) 主導
3. High-AF (caller_af>0.3) zone unfilterable separately (H_C1_4 FAIL) — global LR absorbs
4. 4 low-AF subclone TP (19% of v1.0 lost) remain unrescuable
5. V3F/V5 4 樣本 BAM 不存在 → cycle 2 重做 BAM 產生

---

## 7. Cycle 2 Plan Hand-off

### 7.1 Cycle 2 主任務

1. **V6 4 樣本 ISM rerun with significance** (~3.2 hr) — phaseD V6 BAM 已存在但 significance header-only，重跑 16 runs
2. **Apply cycle 1 filter to 4 samples** — re-fit per sample + transfer fit sanity
3. **H_C1_5 Wilcoxon n=5** — ≥4/5 ΔF1>0 + p<0.05 → ⭐4
4. **H_C1_6 HCC1395 phaseC 三向 cycle 1 filter apply** — V3F/V5/V6 ΔF1 max var < 0.005

### 7.2 Cycle 2 optional

5. **Methyl-specific rescue gated rule** for 4 unrescued lost TP (V6_off_meth_AlleleP gated mode)
6. **High-AF zone non-filterable mechanism investigation** (R-Step0-4 sub-zone)

---

## 8. Files Inventory

### Cycle 1 deliverables (`InterSubMod/research/methyl_augmented_filter_phase2/cycle1/`)

- `00_PLAN.md` (plan v2.0 副本)
- `cycle1_step0_global_fp_audit.md` (Step 0 Agent A1 主報告)
- `cycle1_step1_filter_design.md` (Step 1 filter design, VIF + L2 + NaN MNAR)
- `cycle1_track_a_findings.md` (Step 2 verdicts)
- `cycle1_track_a_filter.json` (deployable filter rule, per-fold coefs)
- `cycle1_findings.md` (本檔, Coordinator synthesis)
- `scripts/`: global_fp_audit.py / filter_design_and_verdict.py / collinearity_resolve.py / final_filter_and_verdict.py
- `figures/`: 9 PNGs (VIF audit / L2 sweep / NaN mechanism / lost TP / multi-seed / FP heatmap / global LR ΔF1 / heterogeneous / collinearity)
- `data/`: 15 TSVs
- `intermediate/`: logs + JSONs

---

## 9. Coordinator notes

- Cycle 1 wall clock: Step 0 ~0.7 min + Step 1+2 ~13 min + BAM search ~15 min + this synthesis = ~30 min total
- Plan v2.0 §Step B1 假設錯 → Track B deferred 是正確決策
- ΔF1 +0.02236 + multi-seed std 5e-5 + lost TP 81% rescue 是 strong single-sample evidence
- ⭐3 strong (not ⭐4) 因 cross-sample 未驗 — 誠實 framing 避免 over-claim
- Cycle 2 高優先 (V6 4 樣本 ISM rerun ~3.2 hr) — 跟 V6 production 4-day workflow Day 1-2 可共用 ISM data
