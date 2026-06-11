<!--
建立時間: 2026-05-21
更新時間: 2026-05-21
作者: Phase A4-F1Audit sub-agent (Opus 4.7)
類別: Methodology completeness audit (F1 6-step audit trail)
狀態: in_progress
標籤: F1_audit, LOSO, data_leakage, overfit, multi_algorithm, PI_4goal
證據鏈: cycle 1 +0.02236 (in-dist) → cycle 4 LOSO -0.00012 (LR) → A4 multi-algo LOSO (DT/RF/XGB)
       → 本 audit (6-step input/output trail + leakage + overfit + cross-algo sanity check)
-->

# Phase A4-F1Audit — F1 6-Step Audit Trail

> **PI 目標 4 要求**：「方法定義與各步驟的輸出與結果狀況，還有 F1 值得狀況，數據要清楚驗證與審查，確認沒有錯誤與誤用」。

---

## §1 Audit Motivation

PI 4-goal Report Goal 4 (methodology completeness) 要求對每一步驟從 caller VCF 到最終 ΔF1 的 input/output 數字逐步驗證、確認 data leakage 不存在、確認 overfit 與 underfit 不掩蓋真實效果，並 cross-algorithm 確認 LR 結論不是 model-specific 偶然。

本 audit 不重跑 LOSO；它只對既有 cycle 1 (in-distribution +0.02236)、cycle 4 (LOSO LR baseline)、A4 (multi-algorithm LOSO) 既有 artifact 做數字逐步對齊（pull source path + line 對齊）。

**Audit scope**：6 個步驟（caller VCF → V6 BAM tagging → ISM region call → LR predict → τ sweep → F1 計算 → LOSO held-out）+ 3 個 leakage check + overfit symptom check + cross-algorithm sanity check。

---

## §2 6-Step Input/Output Trail

### Step 1 — Caller VCF Input (ClairS-TO ssrs paired-pileup)

| Sample | Phased VCF PASS | Total | SNV TP | SNV FP | caller F1 | Source |
|---|---:|---:|---:|---:|---:|---|
| HCC1395 | 47,798 | 3,187,275 | 30,490 | 4,842 | **0.7166** | `20260511_V6_binary_complete_documentation_01.md:587,597` |
| HCC1937 | — | — | 13,910 | 2,697 | **0.3692** | `cycle4/loso_validation/data/loso_cv_results.tsv:3` |
| HCC1954 | — | — | 19,449 | 687 | **0.8385** | `cycle4/loso_validation/data/loso_cv_results.tsv:4` |
| H1437 | — | — | 70,191 | 773 | **0.8670** | `cycle4/loso_validation/data/loso_cv_results.tsv:5` |
| H2009 | — | — | 135,359 | 1,342 | **0.8863** | `cycle4/loso_validation/data/loso_cv_results.tsv:6` |

**Sanity check ✓**：HCC1395 PASS=47,798 與 V3F/V5/V6 三向 phased VCF **file-identity** 完全一致（V6 patch 不改 phasing layer，重用 V5 phased VCF — `20260511 §:597`）。caller F1 0.7166 / 0.6273 cross-binary invariant 已 A0 audit verified。

**Caveat (A0 audit C2.25)**：errata §3a.5b.12 寫 "47,838 SNV PASS" 與 main doc "47,798" 差 40 variants（SNV-only filter 差異，<0.1%，非影響 F1 結論）。

### Step 2 — V6 BAM Tagging (LongPhase-TO V6 patch)

| Sample | HP1 reads | HP2 reads | HP ambig (hp=3) reads | HP1:HP2 ratio | baseline ratio |
|---|---:|---:|---:|---:|---:|
| HCC1395 | 87.8M | 80.7M | **138,317** | **1.609** | 17.3 (4.19:1 4-Layer) |
| 5-sample | — | — | — | **0.611-1.840** (中位數 ≈0.96) | 17.3 |

**Sanity check ✓**：
- 跨樣本 HP1:HP2 ratio 中位數 ≈0.96（vs baseline 17.3:1）→ V6 priority bug 完全修正
- hp_ambig (hp=3) reads ×13.2 over baseline 10,440 → V6 保守化策略
- NG_on=2 rate 5/5 ≥ 0.83 / marker rate 4/5 ≥0.85 gate

**Source**：`20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md:86-88,98`；`20260511_V6_binary_complete_documentation_01.md:36`。

### Step 3 — ISM Region Call (NG ≥ 3 marker filter)

| Pipeline | Marker count | Marker rate (TP purity) | TP markers | FP markers |
|---|---:|---:|---:|---:|
| baseline | 15,738 | 0.8967 | 14,114 | 1,624 |
| V3F | 21,997 | **0.9175** | 20,183 | 1,814 |
| V5 | 18,382 | 0.8937 | — | — |
| **V6** | **23,980** | **0.9093** | **21,806** | **2,174** |

**Δ V6 vs baseline**：marker count **+52.4%** (+8,242) / marker rate **+1.26pp** / hp_ambig reads **×13.2** / LOH-positive marker rate 0.9801 (4 BAM 最高 vs baseline 0.9710)。

**Source**：`20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md:128-130`；`20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md:88,98`。

### Step 4 — Cycle 1 LR Predict (in-distribution 5-fold OOF)

**Input**：HCC1395 35,332 rows × **10 features** (drop NumReads_master VIF=217)
- `V6_off_NG`, `caller_af`, `loh_inner_flag`, `Coverage_Multiple_imp`, `chr8_flag`, plus 5 V6_off_meth_* features (HPMergedDelta / HPFineF / NME_imbalance / Epipoly_Delta / ClusterPermanovaF)

**Model**：L2 Ridge LR (C=1.0, lbfgs) + StandardScaler + 5-fold StratifiedKFold OOF + Strategy B median impute (NaN MNAR)

**Output**：OOF p̂(TP) array length 35,332；feature ranking by |coef|:
- caller_af (+3.44) > LOH_inner (+1.46) > Cov (+1.27) > NG (+1.07) > HPFineF (**+0.75**, methylation 5th-rank) > HPMergedDelta > NME_imbalance > Epipoly_Delta > ClusterPermanovaF > chr8_flag

**Multi-seed stability**：5 seed mean ΔF1=+0.02236 / std=**5e-5**（20× below threshold 0.001）。

**Source**：`cycle1/cycle1_findings.md:82-88,140-150`。

### Step 5 — τ Sweep → Best ΔF1 (in-distribution)

**τ grid**：[0.10, 0.95] step 0.01。

**Best τ\***：0.39（broad plateau 0.38-0.42，multi-seed best-τ std=0.013）。

**Confusion matrix (HCC1395, τ\*=0.39)**：

| | Truth=TP | Truth=FP |
|--|---:|---:|
| Predicted keep (P≥0.39) | **30,015** | 1,443 |
| Predicted filter (P<0.39) | 475 (1.56%) | **3,399 (70.20%)** |

- caller F1 = 0.7166
- post-filter F1 = **0.7390**
- Precision = 30,015 / 31,458 = 0.9541
- Recall = 30,015 / 49,778 = 0.6030 (caller recall 0.6125 slight drop)
- **ΔF1 = +0.02236** = 9.24× v1.0 baseline (+0.00242) = 2.24× Cohen small ribbon (+0.01)

**Source**：`cycle1/cycle1_findings.md:90-99`。

### Step 6 — LOSO Held-out F1 (sample-level cross-validation)

**Protocol**：每次 hold out 1 sample；剩 4 sample concat 為 training pool；per-train median impute + StandardScaler fit on train only；apply 同 medians+scaler 到 held-out test；τ sweep 取 best τ\*。

| Test sample | train rows | test rows | caller F1 | LOSO ΔF1 | best τ\* | TP kept | FP kept |
|---|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 244,408 | 35,332 | 0.7166 | **-0.00012** | 0.10 | 30,490 | 4,842 |
| HCC1937 | 263,133 | 16,607 | 0.3692 | **+0.00000** | 0.10 | 13,910 | 2,697 |
| HCC1954 | 259,604 | 20,136 | 0.8385 | **-0.00008** | 0.10 | 19,449 | 687 |
| H1437 | 208,776 | 70,964 | 0.8670 | **-0.00001** | 0.10 | 70,191 | 773 |
| H2009 | 143,039 | 136,701 | 0.8863 | **-0.00001** | 0.10 | 135,359 | 1,342 |
| **Mean** | — | — | — | **-0.00004** | — | — | — |

- Wilcoxon paired vs 0 (n=5): p=0.125
- Sign: 0 positive / 4 negative / 1 zero → **DIRECTION_NEGATIVE**
- 全 5 best τ → 0.10 = keep everything → filter 在 sample-level **找不到 useful threshold**

**Circularity gap quantification**：

| Metric | Cycle 1 (in-dist 5-fold OOF) | Cycle 4 LOSO (sample-level CV) | Δ |
|---|---:|---:|---:|
| HCC1395 ΔF1 | **+0.02236** | **-0.00012** | **+0.02248** |

**100% effect from sample-level circularity** — 移除 HCC1395 from training set → filter effect 完全消失。

**Source**：`cycle4/loso_validation/data/loso_cv_results.tsv:2-6`；`loso_findings.md:31-39,84-99`。

---

## §3 Data Leakage Audit

### 3.1 Median imputation leakage check ✓ NO LEAKAGE

```python
# run_loso_cv.py:103-104
# Impute using train medians (NOT test sample medians, to avoid leakage)
train_imp, train_medians = impute_with_medians(train_combined, CYCLE1_FEATURES, medians=None)

# run_loso_cv.py:118-120
# Apply to test sample (impute using TRAIN medians to avoid leakage)
test_imp, _ = impute_with_medians(test_df, CYCLE1_FEATURES, medians=train_medians)
```

**Verified**：
- Line 104 `medians=None` → compute new medians from `train_combined` only
- Line 120 `medians=train_medians` → reuse train medians on test → **no test-sample median leak**

**Source**：`research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_cv.py:103-120`。

### 3.2 τ sweep oracle status ⚠ BEST-CASE (over-optimistic)

τ\* 是 **per-fold 在 held-out test 上找的 best ΔF1**（line 129 `sweep_tau(p_test, y_test, ...)`）。

- 這個 τ\* 是 **oracle**（test ground truth 已知才能 sweep）
- 真實 production deployment 必須在 train fold 找 τ 然後 apply（會更差）
- 即使如此 best-case：5 sample LOSO mean=-0.00004 / 全 τ\*=0.10 → filter 在最寬鬆條件仍無效

**Conclusion**：LOSO ΔF1 是 **upper bound（over-optimistic）**；真實 production ΔF1 ≤ -0.00004。即使在這個寬鬆條件下，filter generalization 仍為 NEGATIVE → 結論 robust。

### 3.3 Feature selection per-fold check ✓ NO LEAKAGE

- Cycle 1 step 1 決定 `drop NumReads_master VIF=217` → **feature set 固定為 10 個**，跨 5 LOSO fold 全部一致
- 沒有 per-fold 重選 feature → 沒有 sample-conditioned feature selection leakage
- Multi-seed LR (seed 42/7/100/999/2026) 跑出 std=0 (lbfgs deterministic)

**Source**：`cycle1/cycle1_findings.md:82`；A4 TSV 25 LR rows seed 不影響 ΔF1。

---

## §4 Overfit Symptom Audit

### 4.1 Per-algorithm train-test gap (overall mean)

| Algorithm | overall gap mean | HCC1395 | HCC1937 | HCC1954 | H1437 | H2009 |
|---|---:|---:|---:|---:|---:|---:|
| LR  | +0.0576 | +0.0946 | **+0.4841** | -0.0422 | -0.0893 | -0.1593 |
| DT  | +0.0544 | +0.0616 | **+0.4854** | -0.0390 | -0.0839 | -0.1523 |
| RF  | +0.0551 | +0.0647 | **+0.4845** | -0.0419 | -0.0816 | -0.1501 |
| XGB | +0.0629 | +0.0977 | **+0.4852** | -0.0363 | -0.0822 | -0.1501 |

### 4.2 Key diagnosis: gap 不是 model overfit，是 baseline F1 artifact

**HCC1937 gap +0.48 across 4 algo 完全一致** (差 < 0.002)。如果是 model-specific overfit，DT/RF/XGB（capacity 高）應 > LR；實測四者一致 → **artifact source 是 caller F1 baseline distribution shift, 非 model variance**：
- HCC1937 caller F1 = 0.3692（5 sample 最低）
- train pool 4 sample 平均 caller F1 ≈ 0.83
- train_F1 ≈ 0.85，test_F1 ≈ 0.37 → gap +0.48 自然產生

H1437/H2009 gap negative (-0.08 ~ -0.16) 反向但同源：兩樣本 caller F1 > train pool mean → train_F1 反而低。

### 4.3 排除 LR underfit hypothesis

LR overall gap +0.0576 ≈ DT +0.0544 ≈ RF +0.0551（差 < 0.005）→ 三者 train-side 等價 → **LR 並未 underfit**。

LR 失敗純粹是 **held-out (out-of-distribution) generalization** 問題，與 capacity / linearity 無關。

**Source**：`20260521_A4_multi_algorithm_LOSO_methodology_completeness_01.md:164-185`；A4 TSV 100 rows train_test_gap column。

---

## §5 Cross-Algorithm Sanity Check (4 algo × 5 sample × 5 seed)

### 5.1 Per-algorithm overall ΔF1 (cross-sample, cross-seed)

| Algorithm | mean ΔF1 ± std | n_positive / n_total | Cohen ribbon verdict |
|---|---|---:|---|
| LR  | -0.00004 ± 0.00005 | 0 / 25 | ribbon-null |
| DT  | **+0.00255** ± 0.00573 | 5 / 25 | ribbon-null (< +0.005 small) |
| RF  | **+0.00267** ± 0.00512 | 15 / 25 | ribbon-null (< +0.005 small) |
| XGB | +0.00102 ± 0.00235 | 5 / 25 | ribbon-null |

**沒有任何 algorithm 達到 cross-sample mean > Cohen +0.005 small ribbon** → **跨樣本 stable boost 不存在**。

### 5.2 Per-algorithm × per-sample ΔF1

| Sample | LR | DT | RF | XGB | caller F1 baseline |
|---|---:|---:|---:|---:|---:|
| HCC1395 | -0.00012 | **+0.01378** ✓ | **+0.01269** ✓ | -0.00016 | 0.7166 (中) |
| HCC1937 | +0.00000 | +0.00000 | +0.00000 | **+0.00562** ✓ | 0.3692 (極低) |
| HCC1954 | -0.00008 | -0.00030 | -0.00000 | -0.00020 | 0.8385 (高) |
| H1437 | -0.00001 | -0.00030 | +0.00006 | -0.00009 | 0.8670 (高) |
| H2009 | -0.00001 | -0.00045 | +0.00061 | -0.00008 | 0.8863 (高) |

✓ = ΔF1 > Cohen +0.005 medium ribbon。

**關鍵 cross-algo 觀察**：
- **rescue 機會只出現在 caller F1 < 0.75 的中低 baseline sample** (HCC1395 + HCC1937)
- 高 baseline sample (>0.83) 三條 ML 跑線皆 ribbon-null
- **DT/RF rescue HCC1395 +0.0127~+0.0138** 證明 LR 在該樣本確實 underfit（非線性 interaction 可被 tree capture），但這是 **HCC1395-specific in-distribution capacity**，跨樣本 mean 仍 ribbon-null
- **XGB rescue HCC1937 +0.0056** 是 baseline F1 0.37 低 → ΔF1 上限大導致

**Verdict (final, ⭐3)**: **PARTIAL H_method on HCC1395 only; H_circularity 大方向確認** — sample-level distribution shift 主導，非 linearity 主導。

**Source**：`20260521_A4_multi_algorithm_LOSO_methodology_completeness_01.md:112-150`；`A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv` 100 rows。

---

## §6 PI Q&A 預備

### Q1: 「為何不用 DT / RF？」

**A**: A4 已測 LR / DT / RF / XGB 4 algorithm × 5 sample × 5 seed = 100 LOSO fold (`A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv`)。

- DT/RF 在 HCC1395 rescue +0.0127~+0.0138（超 Cohen medium ribbon），證實 LR 在該樣本 underfit
- **但 cross-sample mean ΔF1 全 algo < Cohen +0.005 small ribbon** (LR -0.00004, DT +0.00255, RF +0.00267, XGB +0.00102)
- **無 algorithm 達 5/5 positive** → 非 model class 限制，是 distribution shift 限制

→ 「用 DT/RF 也救不了 cross-sample generalization」是實證結論，不是猜測。

### Q2: 「Cycle 1 +0.02236 為何不能 production?」

**A**: Sample-level circularity gap = **+0.02236 (in-dist) − (-0.00012 LOSO) = +0.02248**。100% effect 來自 HCC1395 自身 distribution overfit。

- in-dist 5-fold OOF: train 80% HCC1395 → predict 20% HCC1395 → row-level 不重疊，但 **sample-level 完全 circular**
- LOSO: 4 sample train → predict held-out → **5 best τ\* 全部 →0.10 (=keep everything) → filter 無作用**
- production deployment 用新 sample → 等價 LOSO 條件 → ΔF1 ≈ 0

→ 「Cycle 1 +0.02236 是 in-distribution case study, 不是 production filter」。

### Q3: 「會不會是 hyperparameter 不夠？」

**A**: A4 設計避免 hyperparameter leakage:
- 5 seed × 4 algo 全部 mean 落在 [-0.00004, +0.00267] → **seed 對結果幾乎無影響** (LR std=5e-5, RF std=5.1e-3)
- 若用 per-sample tuned hyperparameter → 引入 hyperparameter leakage（held-out τ 已是 oracle，再加 hyperparam oracle 是雙重 leakage）
- DT/RF/XGB 已用「文獻常用 hyperparameter」（max_depth, n_est, lr 等），不是 default-only

→ 「hyperparameter 換不變」（A4 §6 caveat #2）。

### Q4: 「F1 還能怎麼提升？」

**A** (依 A4-Ext 6 算法盤點結論, `20260521_A4_Ext_other_algorithms_inventory_01.md`)：

- **在 LR/DT/RF/XGB framework 內**：基本耗盡（LOSO ≈ 0）
- **Pivot Path 1 - Read-level features**: 從 region-level (NG/AF/LOH) 改 read-level（per-read methylation pattern）→ 完全不同 feature space，cycle 5+ 探索
- **Pivot Path 2 - thread_d Tier 2 already-known signature**: 不依賴 LR filter framework
- **Pivot Path 3 - 加新樣本** (HCC1395_5kHz / COLO829 / DORADO / H2009): cycle 6+ 增 sample variability
- **Honest path - Phase 2 NEGATIVE archive**: accept LOSO 結論，移轉資源到 phase_block_3d

→ 「在 region-level + 4 model class + LOSO 設計下，F1 ceiling 已撞」。

**Source**：`20260521_A4_Ext_other_algorithms_inventory_01.tsv`；`loso_findings.md:200-215`。

---

## §7 結論

### 7.1 Audit verdict

| 維度 | 結論 |
|---|---|
| 6-step input/output trail | **PASS** — 每數字皆有 source path:line；caller F1 0.7166 三向 file-identity，marker count +52.4% 5 處對齊，OOF +0.02236 multi-seed std=5e-5 |
| Data leakage | **NO LEAKAGE** — line 104/120 verified per-train median + train scaler；τ sweep 是 oracle（已知 over-optimistic）；feature set 固定（drop NumReads VIF=217）跨 fold 一致 |
| Overfit symptom | **NOT model-specific overfit** — HCC1937 gap +0.48 across 4 algo 完全一致 → caller F1 baseline artifact；LR 沒 underfit（gap ≈ DT ≈ RF） |
| Cross-algorithm sanity | **H_circularity 確認 + H_method partial** — 4 algo cross-sample mean 全 < Cohen +0.005；DT/RF rescue HCC1395 是 in-distribution capacity；XGB rescue HCC1937 是 baseline F1 低 artifact |
| PI Q&A 預備 | Q1-Q4 預備完整、皆有 A4 實證支撐 |

### 7.2 PI report 引用建議

- 引用 **Step 1 (caller F1 0.7166)** + **Step 6 (LOSO -0.00012)** 並列，明確標 「100% sample-level circularity gap = +0.02248」
- 不單獨引用 Cycle 1 +0.02236 而不附 LOSO caveat
- A4 multi-algo TSV 作為「不是 LR model 限制」的關鍵證據（PI Q1 必問）
- 引用 4 algo HCC1937 gap +0.48 一致 作為「gap 是 baseline artifact 非 overfit」實證

### 7.3 Audit follow-up actions

- 本 audit 屬 read-only methodology audit，**無需更新 evidence_ledger**
- TSV `A4_F1_step_audit_trail.tsv` 可作 PI report 目標 4 §「F1 步驟驗證表」附件
- PI Q1-Q4 預備答案可整合到 PI signoff email draft (`20260521_PI_V6_signoff_email_draft_5goal.md`) 的 「教授可能問題」section

---

## Provenance

- **Audit base scripts**:
  - cycle 1 LR: `research/methyl_augmented_filter_phase2/cycle1/scripts/final_filter_and_verdict.py`
  - cycle 4 LOSO LR: `research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_cv.py:103-129`
  - A4 multi-algo: `research/methyl_augmented_filter_phase2/phase2_completeness_audit/scripts/A4_multi_algo_LOSO.py`
- **Source artifacts**:
  - `20260511_V6_binary_complete_documentation_01.md` — caller F1 0.7166 / 47,798 PASS / V6 file-identity
  - `20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md` — marker +52.4% / hp_ambig ×13.2 / HP ratio 1.609
  - `cycle1/cycle1_findings.md` — in-dist +0.02236 / τ\*=0.39 / coef ranking
  - `cycle4/loso_validation/data/loso_cv_results.tsv` — 5 sample LOSO LR results
  - `A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv` — 4 algo × 5 sample × 5 seed
  - `20260521_A4_multi_algorithm_LOSO_methodology_completeness_01.md` — per-algo overall + gap analysis
- **Audit output**:
  - 本檔 `docs/experiments/in_progress/2026/05/20260521_A4_F1_step_by_step_audit_01.md`
  - TSV `research/methyl_augmented_filter_phase2/phase2_completeness_audit/A4_F1_step_audit_trail.tsv`
- **Wall clock**: ~1.5 hr (read-only audit, no re-execution)
