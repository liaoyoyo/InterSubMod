<!--
建立時間: 2026-05-18
更新時間: 2026-05-18
agent: Coordinator (main session) Phase 2 Cycle 1 synthesis
status: in_progress
report_class: cycle synthesis (HCC1395 single-sample, Track B deferred)
audience: PI / lab member / 自己未來
scope: HCC1395 clairs_to_ssrs single-sample global FP filter (V0.3 → V1.0 → Phase 2 Cycle 1 lineage)
tier: ⭐3 strong (HCC1395-internal validated, multi-seed stable, 9.24× v1.0 ΔF1)
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md v2.0
predecessor_cycles:
  - v0.3 characterization (2026-05-15, ⭐3)
  - v1.0 step5_methyl_filter_pilot (2026-05-18 早, ⭐3 marginal +0.00242)
inputs:
  - v1.0 master TSV: research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv (35,332 × 202)
  - SEQC2 truth: research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv
outputs:
  - 本檔（主報告）
  - research/methyl_augmented_filter_phase2/cycle1/cycle1_findings.md (詳細 synthesis)
  - research/methyl_augmented_filter_phase2/cycle1/cycle1_track_a_filter.json (deployable filter)
verdict: H_C1_1/2/3 PASS / H_C1_4 FAIL / H_C1_5/6 DEFERRED → ⭐3 strong; cycle 2 升 ⭐4 條件已明確
last_verified: 2026-05-18
report_template: cycle-pilot v1.0
-->

# Phase 2 Cycle 1 — Global FP Filter (HCC1395) — ⭐3 strong

> **Tier**: ⭐3 strong (HCC1395 single-sample verified, ΔF1 +0.02236 = 9.24× v1.0 baseline)
> **Predecessor**: v1.0 step5_methyl_filter_pilot ⭐3 marginal (+0.00242)
> **Goal achieved**: ΔF1 從 v1.0 marginal +0.00242 → 提升至 +0.02236 (Cohen 小 effect ribbon 2.24×)

---

## 0. TL;DR

從 v1.0 「cell-gated filter」pivot 為「global FP exploration + multi-axis LR」，HCC1395 single-sample 達 **ΔF1 = +0.02236** (9.24× v1.0 baseline, 2.24× +0.01 Cohen ribbon)，**3/4 預註冊假設 PASS**：

- ✅ **H_C1_1** Top 10 cells covers 94.22% all FP (vs 70% threshold)
- ✅ **H_C1_2** Global LR ΔF1 +0.02236 vs v1.0 +0.00242 (9.24×)
- ✅ **H_C1_3** ΔF1 ≥ +0.01 Cohen ribbon (2.24×)
- ❌ **H_C1_4** High-AF zone FAIL (-0.00011, sub-zone limitation 但 global LR 仍 capture)
- 🔵 **H_C1_5/6** DEFERRED to cycle 2 (V3F/V5 4 樣本 BAM 物理不存在)

**核心發現**：
1. **Step 5c lost TP 17/21 (81%) 被 cycle 1 filter rescue** — 證明 v1.0 cell-gated 過嚴；global LR 不需 rescue rule 就能救
2. **Methylation 是 5th-rank covariate (HPFineF coef +0.75)** — 主導為 caller_af + LOH + Cov + NG；filter naming 改 "multi-axis filter incl. methylation"
3. **Multi-seed std = 5e-5** (20× below threshold) — 高度穩定 intra-sample

---

## 1. Context — 為什麼這個 cycle 存在

### 1.1 Predecessor (v1.0 cycle, 2026-05-18 早晨)

V1.0 step5_methyl_filter_pilot ⭐3 marginal:
- ΔF1 +0.00242 < +0.005 Cohen ribbon
- v1.0 powered cell gate (n≥50 + Model B converged + n_TP≥5) 過嚴 → 4 cells 只 covers **7% all FP**
- Step 5c TP rescue NEGATIVE (95% lost TP 是 low-AF subclone)

### 1.2 用戶 pivot direction (2026-05-18 下午)

「在某些 FP 比例高的組合區域用較嚴格條件去除 FP，少量 TP 影響，達到更好 F1 提升超過 1%？」

→ Plan v2.0 pivot strategy: **global LR (no cell gate) + heterogeneous threshold + 4 zone audit**

### 1.3 Track B 物理限制（cycle 中段發現）

Plan v2.0 §Step B1 假設「V3F/V5 4 樣本 BAM 已存在」**錯**：
- Deep search agent 確認只 HCC1395 有 V3F/V5
- 4 樣本 (H1437/H2009/HCC1954/HCC1937) 只 V6 在 `v6_5sample_extension/`
- 用戶決定 cycle 1 在 HCC1395 single-sample 收尾（cross-sample 留 cycle 2）

---

## 2. Pre-Registration (Hard Gate per scientific-rigor §7.1)

| ID | Prediction | Falsification | Decision threshold |
|----|-----------|---------------|-------------------|
| H_C1_1 | Top 10 FP-rich cells covers ≥70% all FP | <70% | concentration map D1 |
| H_C1_2 | Global LR ΔF1 > v1.0 +0.00242 | ≤+0.00242 | global LR D2 |
| H_C1_3 | Heterogeneous OR global aggregate ΔF1 ≥ +0.01 Cohen | <+0.01 | aggregate D3 |
| H_C1_4 | High-AF zone incremental ΔF1 ≥ +0.003 | <+0.003 | high-AF D4 |
| H_C1_5 | Cross-sample n=5 ≥4/5 ΔF1>0 + Wilcoxon p<0.05 | ≤2/5 | cross-sample (DEFERRED) |
| H_C1_6 | V3F/V5/V6 BAM max var <0.005 | >0.01 | sanity (DEFERRED) |

**NO-GO trigger**: H_C1_2 + H_C1_3 同時 FAIL — **未觸發**

---

## 3. Methods

### 3.1 Workflow (Multi-agent fan-out)

```
Agent A1 (Step 0 audit, foreground ~0.7 min):
  Stage 1 — FP concentration map (120-cell)
  Stage 2 — High-AF FP profile
  Stage 3 — Global LR sweep [τ 0.10-0.95]
  Stage 4 — Heterogeneous top 20 cells
  → D1/D2/D3/D4 verdicts → Path selection

Agent B1 (Step B1 pre-flight, ~5 min):
  Pre-flight check → discover V3F/V5 4-sample BAMs NOT EXIST → halt

Agent A2 (Step 1+2 Track A, ~13 min):
  Step 1 — VIF audit + L2 sweep + NaN MNAR + filter design
  Step 2 — HCC1395 ΔF1 verdict + Step 5c lost TP cross-check + multi-seed

Deep BAM search Agent (~15 min):
  Confirm V3F/V5 4-sample BAMs物理不存在

Coordinator (main session):
  Track B deferred decision + cycle 1 synthesis
```

### 3.2 Final Filter Specification

| Component | Value |
|-----------|-------|
| Features (10) | V6_off_NG, caller_af, loh_inner_flag, Coverage_Multiple_imp, V6_off_meth_{HPMergedDelta, HPFineF, NME_imbalance, Epipoly_Delta, ClusterPermanovaF}, chr8_flag |
| Excluded (VIF=217) | NumReads_master |
| Regularization | L2 (Ridge), C=1.0 |
| Optimal τ* | **0.39** (broad plateau 0.38-0.42) |
| NaN handling | Strategy B (median impute), MNAR-justified |
| CV | 5-fold StratifiedKFold OOF |
| Deployable artifact | `InterSubMod/research/methyl_augmented_filter_phase2/cycle1/cycle1_track_a_filter.json` |

### 3.3 嚴謹度防護

- Pre-registration 6 H pre-specified (manifest.yaml + plan v2.0)
- VIF audit + L2 regularization sweep (R-Step0-2 resolved)
- NaN MNAR test (R-Step0-1 resolved: NaN 80% at AF<0.1, 92% at AF>0.7 U-shape)
- 5-fold StratifiedKFold OOF (anti-optimism)
- 5 random seeds × 5-fold = 25 CV runs (multi-seed stability)
- BH-FDR q<0.05 全局校正

---

## 4. Results

### 4.1 Step 0 Global FP Audit

| Decision | Threshold | Observed | Verdict |
|----------|-----------|----------|---------|
| **D1** Top 10 cells coverage | ≥70% | **94.22%** | ✅ PASS |
| **D2** Global LR ΔF1 vs v1.0 baseline | > +0.00242 | **+0.02637** (Strategy B impute, τ=0.42) | ✅ PASS |
| **D3** Heterogeneous aggregate ΔF1 | ≥ +0.01 | +0.00175 | ❌ FAIL |
| **D4** High-AF zone incremental | ≥ +0.003 | -0.00012 | ❌ FAIL |

→ **Path B selected** (pure global LR; D2 PASS + D3 FAIL)

**FP concentration (top 5)**:
1. [0.00,0.10) Outer Cov[0.9,1.1) other — **34.12%**
2. [0.10,0.20) Outer Cov[0.9,1.1) other — 28.38%
3. [0.20,0.30) Outer Cov[0.9,1.1) other — 8.47%
4. [0.10,0.20) Outer Cov[0.9,1.1) **chr8** — 7.91%
5. [0.00,0.10) Outer Cov[0.9,1.1) chr8 — 5.12%

Top-3 (low-AF Outer normal-Cov) 涵蓋 **71% all FP**。

### 4.2 Step 1 Filter Design

**VIF audit**: NumReads_master VIF=**217.17**, Coverage_Multiple_imp VIF=**215.06** (severe collinearity)
**L2 sweep**: C=0.001→10.0 ΔF1 +0.0198→+0.0299，但 max|coef| 爆 1.13→68.6 (unregularized saddle)
**Resolution**: drop NumReads_master + L2 C=1.0 → max VIF=1.83, max|coef|=3.44

**NaN MNAR test**: NME_imbalance / Epipoly_Delta NaN 60.4% are MNAR (chi-square p=0, MW p<1e-229, U-shape across AF bins)
**Strategy B (impute) 正確** — chosen filter methyl coefs <1.0, non-methyl axes 主導

### 4.3 Step 2 HCC1395 ΔF1 Verdict

| Metric | Value |
|--------|-------|
| Filter config | cfg_drop_nr (10 features), L2 C=1.0 |
| τ* | 0.39 |
| TP kept | 30,015 / 30,490 (TP loss 1.56%) |
| FP removed | 3,399 / 4,842 (**70.20%**) |
| Precision | 0.9541 |
| Recall | 0.6030 |
| **ΔF1 vs caller F1=0.7166** | **+0.02236** (9.24× v1.0) |

**Step 5c lost TP cross-check**: **17/21 (81%) rescued** by cycle 1 filter (vs v1.0 Step 5c 0%)

**Multi-seed CV (5 seeds)**:
- Mean ΔF1 +0.02236, **std = 5e-5** (20× below 0.001 threshold) → **STABLE intra-sample**
- τ* plateau [0.38-0.42] across seeds

### 4.4 Feature Importance (10-feature LR, L2 C=1.0)

| Rank | Feature | Coef | Type |
|------|---------|------|------|
| 1 | caller_af | +3.44 | Caller |
| 2 | loh_inner_flag | +1.46 | Structural |
| 3 | Coverage_Multiple_imp | +1.27 | Structural |
| 4 | V6_off_NG | +1.07 | Phasing |
| **5** | **V6_off_meth_HPFineF** | **+0.75** | **Methylation** |
| 6 | V6_off_meth_HPMergedDelta | +0.42 | Methylation |
| 7 | V6_off_meth_NME_imbalance | +0.31 | Methylation |
| 8 | V6_off_meth_Epipoly_Delta | +0.28 | Methylation |
| 9 | V6_off_meth_ClusterPermanovaF | +0.18 | Methylation |
| 10 | chr8_flag | +0.15 | Sample-specific |

**Methylation 為 5th-rank covariate** — 非主導但有 incremental information。

---

## 5. Track B Deferred (V3F/V5 4-sample BAM 物理不存在)

Deep BAM search agent 確認：

| Sample | V3F BAM | V5 BAM |
|--------|---------|--------|
| HCC1395 | ✅ `/big7_disk/.../pononly_v3_fixed/tumor_tagged.bam` | ✅ `/big7_disk/.../pononly_v5_somatic_fallback/tumor_tagged.bam` |
| H1437 | ❌ NOT_FOUND | ❌ NOT_FOUND |
| H2009 | ❌ NOT_FOUND | ❌ NOT_FOUND |
| HCC1954 | ❌ NOT_FOUND | ❌ NOT_FOUND |
| HCC1937 | ❌ NOT_FOUND | ❌ NOT_FOUND |

→ Plan v2.0 §Step B1 「V3F+V5 16 runs」**物理不可能**
→ 用戶決定 cycle 1 在 HCC1395 single-sample 收尾，**cross-sample 留 cycle 2 (Path B V6-only)**

---

## 6. R-Step0 5 個 Caveats 處理狀態

| ID | Severity | Status | Resolution |
|----|----------|--------|-----------|
| R-Step0-1 NaN gap 8.7× | HIGH | ✅ RESOLVED | MNAR confirmed, impute correct |
| R-Step0-2 Collinearity VIF=217 | HIGH | ✅ RESOLVED | Drop NumReads_master, max VIF 1.83 |
| R-Step0-3 Methylation marginal | MED | ✅ ACKNOWLEDGED | Methyl 5th-rank, framing 改 |
| R-Step0-4 Step 5c lost TP overlap | MED | ✅ RESOLVED | 81% rescue rate |
| R-Step0-5 HCC1395-only | HIGH | 🔵 OPEN | Track B deferred to cycle 2 |

---

## 7. Limitations

1. **HCC1395 single-sample only** — generalizability 未驗
2. **Methylation 是 5th-rank** — non-methyl axes (AF/LOH/Cov/NG) 主導；不能宣告「methylation filter」
3. **High-AF zone unfilterable separately** (H_C1_4 FAIL)
4. **4 low-AF subclone TP (19% of v1.0 lost) remain unrescuable** — future work: methyl-specific gated rescue
5. **V3F/V5 4 樣本 BAM 不存在** — cycle 2 用 Path B V6-only OR 重做 BAM 產生

---

## 8. Tier Evaluation → **⭐3 strong**

| 維度 | 評估 |
|------|------|
| Effect size | ΔF1 +0.02236 = 2.24× Cohen small ribbon ✅ |
| Pre-registration | 4 H pre-specified, 3 PASS / 1 FAIL ✅ |
| Confound guard | VIF + L2 + MNAR + 5-fold OOF + multi-seed ✅ |
| Stability | std 5e-5 ✅ |
| Lost TP rescue | 81% (vs v1.0 0%) ✅ |
| Cross-sample | DEFERRED ❌ |
| Mechanism | non-methyl axes 主導 (acknowledged) ⚠️ |

⭐3 而非 ⭐4 的理由: H_C1_5/6 DEFERRED → cross-sample 未驗

### ⭐4 升級必要條件 (cycle 2)
1. V6 4 樣本 ISM rerun with significance (~3.2 hr, BAM 已存在 v6_5sample_extension/)
2. Apply cycle 1 filter → Wilcoxon n=5
3. HCC1395 phaseC V3F/V5/V6 三向 cycle 1 filter 驗 H_C1_6

---

## 9. Paper Framing (保守 framing)

> **"Multi-axis filter (caller AF + LOH inner + Coverage_Multiple + HPFineNGroups + 5 methylation covariates) on HCC1395 clairs_to_ssrs achieves ΔF1 = +0.02236 vs ClairS-TO caller baseline F1=0.7166 (post-filter F1=0.7390). Methylation contributes as 5th-rank covariate (HPFineF coef +0.75). 81% of v1.0 cell-gated lost TP rescued by global LR. Cross-sample generalization pending (only HCC1395 has V3F/V5 BAMs available)."**

### §3 主圖建議
- Fig 3a: τ sweep ΔF1 curve (multi-seed ± band) → plateau 0.38-0.42
- Fig 3b: VIF audit + L2 regularization sweep
- Fig 3c: Step 5c 21 lost TP individual P(TP) → 17/21 rescued
- Fig 3d: Feature importance ranking (methyl 5th-rank)

### Prior art 差異化
- TumorLens: sample-level, not per-region
- ROCIT: methylation-only transformer
- SGZ: variant-level 4-axis no phasing/methyl
- **本 cycle**: per-region multi-axis variant-level filter — **無同口徑 prior art**

---

## 10. Cycle 2 Plan Preview

1. V6 4 樣本 ISM rerun with significance (~3.2 hr background)
2. Per-sample augmented master TSV (4 files)
3. Apply cycle 1 filter cross-sample + Wilcoxon n=5
4. HCC1395 phaseC V3F/V5/V6 三向 cycle 1 filter apply (H_C1_6)
5. (Optional) Methyl-specific gated rescue for 4 unrescuable TP
6. 與 V6 production 4-day workflow 共用 ISM data（Day 2 6-sample marker coverage）

---

## 11. Reproducibility — File Inventory

### Cycle 1 deliverables (`InterSubMod/research/methyl_augmented_filter_phase2/cycle1/`)
- `00_PLAN.md` (plan v2.0 副本)
- `cycle1_step0_global_fp_audit.md` (Step 0 主報告)
- `cycle1_step1_filter_design.md` (Step 1 filter design)
- `cycle1_track_a_findings.md` (Step 2 verdicts)
- `cycle1_findings.md` (Coordinator synthesis, 詳細版)
- `cycle1_track_a_filter.json` (deployable filter rule)
- `scripts/`: global_fp_audit.py / filter_design.py / collinearity_resolve.py / final_filter_and_verdict.py
- `figures/`: 9 PNGs
- `data/`: 15 TSVs

### Predecessor reference
- v1.0 master TSV: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv`
- v1.0 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260518_V6_Methyl_Filter_Pilot_01.md`
- v0.3 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md`
- Phase 2 project plan: `InterSubMod/research/methyl_augmented_filter_phase2/00_PLAN.md`

### Plan
- Cycle 1 plan v2.0: `/bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md`

---

## 12. Decision Audit Trail

| Step | Decision | 來源 |
|------|----------|------|
| Plan v2.0 設計 (3 user decisions) | 1) Path B; 2) PubMed mechanism; 3) cross-sample 留 | 2026-05-18 AskUserQuestion |
| Step -1 amendment | Add phaseC rerun (v1.0 Step -1 同等發現) | 2026-05-18 Agent A1 discover |
| Track B deferred decision | V3F/V5 4 樣本 BAM 不存在 → HCC1395 only | 2026-05-18 用戶 「只先驗證 HCC1395」 |
| Manual SoT update | INDEX + CURRENT_FOCUS + ledger + memory | 本 commit |
