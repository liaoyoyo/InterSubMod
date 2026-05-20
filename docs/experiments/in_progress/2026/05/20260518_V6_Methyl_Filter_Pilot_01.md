<!--
建立時間: 2026-05-18
更新時間: 2026-05-18
agent: Coordinator synthesis (main session) of v1.0 Methylation-Augmented FP Filter Pilot
status: in_progress
report_class: filter-pilot-evaluation (cross v0.3 characterization → filter F1 boundary)
audience: PI / lab member / 自己未來
scope: HCC1395 clairs_to_ssrs single-sample pilot, V3F+V5+V6 三 BAM × off/on flag × 4 powered-gate FP-rich cells
tier: ⭐3 PARTIAL POSITIVE (5/5 H POSITIVE but ΔF1 marginal +0.00242)
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md (v1.0 + Step -1 amendment)
predecessor_cycle: 2026-05-15 V3F/V5/V6 ISM 三向 × LOH × HP × CN characterization ⭐3 (InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md)
inputs:
  - research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/ (NEW 2026-05-18, 12 ISM runs with significance enabled)
  - research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_three_way.tsv (v0.3 base)
  - research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv (SEQC2 truth ref)
outputs:
  - 本檔（主報告）
  - research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step{0-5}_*.md + figures + scripts
upstream_reports:
  - 2026-04-26 Thread D 主軸 (TP-enriched phasing signatures, LOH × cross_het)
  - 2026-05-15 V0.3 characterization ⭐3
  - paired_priority_bug_audit/09_V6_caller_F1_verification.md
verdict: ALL H1-H5 POSITIVE; ΔF1 = +0.00242 (POSITIVE-but-MARGINAL, < +0.005 Cohen ribbon); ⭐3 candidate eligible for cross-sample phase 2
decision: Option A (依 Coordinator recommendation + 用戶 confirm 2026-05-18) — 進入 cross-sample phase 2 pilot
last_verified: 2026-05-18
report_template: filter-pilot v1.0
-->

# Methylation-Augmented FP Filter Pilot (HCC1395 single-sample, V3F+V5+V6)

> **Verdict**: ⭐3 PARTIAL POSITIVE — ALL 5 hypotheses POSITIVE, but **ΔF1 = +0.00242 is MARGINAL** (< +0.005 Cohen ribbon)
> **Decision (用戶 confirmed 2026-05-18)**: Option A — 進入 cross-sample phase 2 pilot
> **Predecessor**: v0.3 characterization cycle ⭐3 PARTIAL POSITIVE (2026-05-15)

---

## 0. TL;DR

本 cycle 在 v0.3 characterization 之上跨越到 **filter F1 evaluation** (objective per resolved decision: 評估 FP_removal_pct > TP_loss_pct AND ΔF1 > 0 vs ClairS-TO caller baseline F1=0.7166)。

**Workflow**：phaseC ISM 12 runs 重跑加 significance computation (Step -1) → master TSV 加 13 methylation features × 3 BAM × 2 flag (Step 0) → 138 augmented LR + LRT (Step 1) → FP-rich cells τ sweep (Step 2) → ΔF1 vs caller (Step 3) → mechanism brainstorm + PubMed (Step 4) → decision synthesis (Step 5)。

**5 hypothesis verdicts**:

| H | Verdict | Effect | 證據 |
|---|---------|--------|------|
| **H1** ≥1 methylation covariate LRT q<0.05 | **POSITIVE** | 16/30 testable cells | top `Outer\|other\|cov_normal` p=1.8e-58 |
| **H2** max ΔF1 > 0 vs caller 0.7166 | **POSITIVE-marginal** | **+0.00242** | post-filter F1=0.71902, < +0.005 Cohen ribbon |
| **H3** FP_removal > TP_loss at τ* | **POSITIVE** | 98.26% > 35.00% | filter_signal +0.633 |
| **H4** ≥1 mechanism candidate | **POSITIVE** | 13 hypotheses × 14 PubMed refs | 5 categories (cis-mQTL/cancer ASM/allele-imbalance/repeat/replication timing) |
| **H5** V5 ≈ V6 LR β sanity | **POSITIVE** | Median Δβ = 1.87e-5 << 0.005 | V6 重用 V5 phased VCF 預期內 |

**Decision tree**: ALL H1-H4 POSITIVE + H5 sanity OK → **「Methylation-augmented FP filter ⭐3 candidate」** (依 plan v1.0 §Step 5 decision tree)

**但 ΔF1 marginal 是重要 caveat** — 必須在 paper 與 lab meeting 中誠實標示「需 cross-sample 才能升 ⭐4」。

---

## 1. Context — 為什麼這個 cycle 存在

### 1.1 v0.3 留下的問題（precursor）

v0.3 characterization cycle 找到 3 個核心 finding：
1. H4 chr8 hotspot LR deviance: CN 0.211 > HP 0.063（(LOH+CN)-HP = +0.186）
2. Paradigm reframe: Z-OCH/Z-GL 是 TP signatures 而非 FP markers
3. **Framework gap**: 4 zone framework 只 cover ~37% FP，63% FP unexplained

用戶於 2026-05-18 提出新問題：**「63% unexplained FP 用 methylation augmentation 能否轉成 actionable filter，使 FP removal % > TP loss % 且 ΔF1 > 0 (相對於 ClairS-TO caller-only F1=0.7166)？」**

### 1.2 與歷史 filter NO-GO 的差異化（已記錄於 known-pitfalls）

| 歷史 NO-GO | 本 cycle 差異 |
|------------|-------------|
| LOH binary filter 10/10 (2026-04-06) | 不用 binary；用 cell-level multi-axis LR-predicted threshold |
| CN zone-aware filter 跨樣本崩 (2026-04-10) | HCC1395 pilot 先；cross-sample 留 phase 2 |
| Pure methylation clustering AUC 0.512 (Option C) | 不用 pure；當 5th-9th covariate **on top of** 4 axis |
| O9-O13 甲基化 NEGATIVE | 加 within-cell within-AF OLS confound guard 防 collider bias |

### 1.3 Caller F1 baseline (0.7166)

- 來源: `InterSubMod/research/paired_priority_bug_audit/06_V3F_vs_V5_evaluation.md` §1 + `09_V6_caller_F1_verification.md`
- HCC1395 ClairS-TO vs SEQC2 truth set
- V3F = V5 = V6 三版相同（重用 phased VCF，caller F1 不變）
- Reverse-solve: TP=30,490 / FP=4,842 / **FN=19,288** (反解) / N_truth=49,778

---

## 2. Pre-Registration (依 scientific-rigor §7.1, 不可事後改寫)

| ID | Prediction | Falsification | Decision threshold |
|----|-----------|---------------|-------------------|
| H1 | ≥1 methylation covariate LRT q<0.05 | 0 cells reach q<0.05 | BH-FDR q<0.05 |
| H2 | HCC1395 ΔF1 > 0 vs caller 0.7166 | ΔF1 ≤ 0 | post-filter F1 - 0.7166 > 0 |
| H3 | FP removal % > TP loss % at τ* | FP removal ≤ TP loss | region-level confusion matrix |
| H4 | ≥1 candidate mechanism | 0 候選 | relaxed gate per resolved decision |
| H5 | V5 vs V6 LR β diff ≤ 0.003 | > 0.005 | sanity check |

**NO-GO 條件**: H1 OR H2 OR H3 任一 NEGATIVE → 寫 evidence_ledger 不可事後改寫

**Pre-reg 結果**: 全 H POSITIVE → **不觸發 NO-GO** → 進 paper framing

---

## 3. Methods

### 3.1 Workflow chain

```
Step -1 (新增): rerun phaseC ISM 12 runs with significance
   │ (Removed --no-distance-matrix flag; ~80 min, V3F+V5+V6 × on/off × tp/fp)
   ▼
Step 0: build augmented master TSV (12.4 sec)
   │ Join 13 methylation features × 3 BAM × 2 flag = 78 cols → master 88 + 114 = 202 cols
   ▼
Step 1: 3 BAM × Augmented LR + LRT (~16 min Agent G)
   │ 138 fits Model A (baseline 3 covariates) vs Model B (+5 methyl covariates)
   │ Within-cell within-AF-bin OLS + 5-fold CV
   │ BH-FDR q<0.05 全局校正
   ▼
Step 2: FP-rich cells filter threshold sweep (~16 min Agent G)
   │ 12 FP-rich cells × τ 0.5-0.95 by 0.01 = 540 evaluations
   │ TP_loss_pct vs FP_removal_pct
   ▼
Step 3: ΔF1 vs caller baseline (Agent H ~6 min)
   │ CV out-of-fold predictions (anti-optimism)
   │ τ 擴展 [0.10, 0.95] (Step 2 floor recovery)
   │ Strategy 1: max ΔF1 / Strategy 2: sanity check
   ▼
Step 4: Mechanism brainstorm + PubMed (Agent I ~8 min, 平行 with Step 1-3)
   │ 5 categories × 1-3 hypotheses × ≥1 PubMed DOI/PMID
   ▼
Step 5: Coordinator decision synthesis (main session)
   │ H1-H5 integrated verdict + paper framing
```

### 3.2 Data sources

- **3 BAM (V3F+V5+V6)**: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/{pononly_v3_fixed, threshold_compare/v5_flag, v6_germline_absent_revert}/tumor_tagged.bam` (each ~287 GB)
- **VCF (TP/FP)**: `/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_{tp,fp}.vcf.gz`
- **LOH BED**: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_phased_LOH.bed`
- **SEQC2 truth ref**: `InterSubMod/research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv`
- **NEW 12 ISM runs**: `InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/` (~25 GB)

### 3.3 12 methylation features used

從 ISM significance_summary.csv 撈：
1. HPMergedDelta, HPMergedP, HPMergedSig
2. HPFineF, HPFineP, HPFineSig, HPFineNGroups
3. ClusterPermanovaF, ClusterPermanovaP
4. AlleleDelta, AlleleP
5. Entropy_Imbalance
6. Derived: NME_imbalance = |NME_HP1 - NME_HP2|
7. Derived: Epipoly_imbalance = |Epipoly_HP1 - Epipoly_HP2|

Augmented LR Model B 主要用 5 個 (HPMergedDelta, HPFineF, NME_imbalance, Epipoly_Delta, ClusterPermanovaF)

### 3.4 嚴謹度防護 (scientific-rigor §rigor)

| 防護 | 落實狀態 |
|------|---------|
| Pre-registration 5 H 事先寫死 | ✅ 寫於 plan v1.0 §Pre-Registration |
| NO-GO 條件不可事後改 | ✅ ALL POSITIVE，未觸發 |
| Within-cell within-AF OLS | ✅ Step 1 confound guard |
| 5-fold CV out-of-fold | ✅ Step 1+3 用 CV preds |
| BH-FDR q<0.05 全局 | ✅ 138 cells 統一校正 |
| V5 vs V6 sanity (H5) | ✅ Median Δβ = 1.87e-5 |
| 12 methylation features 全 join | ✅ Step 0 不 cherry-pick |
| Mechanism relaxed gate | ✅ Step 4 列 13 候選 + PubMed |
| Effect size calibration | ⚠️ ΔF1 +0.00242 < +0.005 marginal — 必標示 |

---

## 4. Results

### 4.1 Step -1: phaseC rerun (Hard Gate amendment 2026-05-18)

**問題發現**：原 phaseC ISM 12 runs 用 `--no-distance-matrix` flag → silently 跳過 significance computation → 12 features 不存在
**解決**：移除 flag 重跑 12 runs

| BAM | flag | label | rows | time |
|-----|------|-------|------|------|
| V3F | off | tp/fp | 30,476 / 4,821 | ~12 min each |
| V3F | on | tp/fp | 30,476 / 4,821 | ~12 min each |
| V5 | off | tp/fp | 30,476 / 4,821 | ~12 min each |
| V5 | on | tp/fp | 30,476 / 4,821 | ~12 min each |
| V6 | off | tp/fp | 30,476 / 4,821 | ~12 min each |
| V6 | on | tp/fp | 30,476 / 4,821 | ~12 min each |

**12 runs total = 80 min** (vs estimated 6-10 hr — 加 significance 只多 ~2 min/run，因 PERMANOVA threading 高效)
**Success rate**: 99.95% TP / 99.57% FP
**13 methylation features 全 confirmed** in header

### 4.2 Step 0: Augmented master TSV

- Shape: **35,332 rows × 202 cols** (88 base + 114 methyl)
- Non-null:
  - HP/Cluster/Allele/HPFine* 99.90%
  - NME/Epipoly 19-40%（生物學上 sparse — intra-HP entropy 需 HP 內有足夠 reads）
- Build time: **12.4 sec** (pandas merge 高效)
- File: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv` (25.3 MB)

### 4.3 Step 1: Augmented LR + LRT (H1, H5)

#### H1 verdict (POSITIVE)
- Model A: 138/138 converged
- Model B: **30/138 converged** (22%) — 108 small-FP cells (n_FP<5) saturated（expected behavior）
- **LRT q<0.05 BH-FDR**: 16/30 cells POSITIVE
- Top LRT signals:
  - `Outer|other|cov_normal` × 3 BAMs × 2 flags: LRT p down to **1.8e-58**, ΔAUC +0.02
  - `Outer|other|cov_elevated`
  - `Inner|other|cov_normal`
  - `Inner|other|cov_elevated`

#### H5 verdict (POSITIVE)
- Raw mean abs Δβ (V5 vs V6) = 5.81 (inflated by singular fits in tiny cells)
- **Robust median abs Δβ** (converged + n_eff_B≥100 + methyl-only, n=30): **1.87e-5** << 0.005 threshold
- Confirms Step 0 sanity (V5 vs V6 methyl Pearson r = 1.0000)

#### In-sample optimism warning
- In-sample AUC 0.97 vs CV AUC 0.90 (~7 AUC points)
- Step 3 必用 CV out-of-fold preds 計算 ΔF1

### 4.4 Step 2: FP-rich cells filter sweep (H3 direction)

- 從 v0.3 step3 抓 12 FP-rich cells (Z-CHR8 sub-cells + Z-AUTO sub-cells)
- τ sweep: 0.5-0.95 by 0.01 = 45 values × 12 cells = 540 evaluations
- **9/11 FP-rich cells 有 τ 使 FP_removal_pct > TP_loss_pct**（2 cells trivial UNKNOWN）
- **Best filter signal cell**: `auto|Outer|other|cov_proxy_mid` filter_signal=**+0.80** (FP_remov 99.2%, TP_loss 19.2%)
- **Caveat**: 所有 best τ 撞 grid floor 0.50 → in-sample bias，Step 3 擴展 τ 範圍

### 4.5 Step 3: ΔF1 (H2, H3)

#### Caller baseline (反解)
```
TP_caller = 30,490 / FP_caller = 4,842
P = 0.8630
F1 = 0.7166 → R = 0.6125 → FN = 19,288 → N_truth = 49,778
F1 sanity = 0.71648 ✓
```

#### Per-cell ΔF1 (Strategy 1: max ΔF1, CV preds)

| Scope | τ* | ΔF1 | TP_loss% | FP_remov% | F1_post | CV AUC |
|-------|----|----|---------|----------|---------|--------|
| **AGGREGATED** (4 cells merged) | **0.52** | **+0.00242** | 35.0% | 98.3% | **0.71902** | n/a |
| chr8|Outer|other|cov_normal | 0.52 | +0.00136 | 42.3% | 98.5% | 0.71796 | **0.928** |
| chr8|Inner|other|cov_normal | 0.18 | +0.00089 | 23.5% | 94.1% | 0.71749 | 0.849 |
| auto|Outer|other|cov_proxy_mid | 0.79 | +0.00071 | 38.5% | 98.3% | 0.71731 | 0.853 |
| auto|Inner|other|cov_proxy_mid | 0.22 | +0.00045 | 29.4% | 91.6% | 0.71705 | 0.796 |

#### H2 verdict (POSITIVE-marginal)
- max ΔF1 = **+0.00242** @ τ*=0.52 (well off Step 2 floor 0.50)
- post F1 = 0.71902 vs caller 0.7166
- **POSITIVE direction** but **< +0.005 Cohen ribbon "marginal"** → MARGINAL effect size

#### H3 verdict (POSITIVE)
- At τ*=0.52: FP_removal_pct = **98.26%**, TP_loss_pct = **35.00%**
- filter_signal = +0.633 (FP_removal − TP_loss > 0)

### 4.6 Step 4: Mechanism + PubMed (H4)

**13 hypotheses × 14 PubMed refs** (full table in `step4_mechanism_candidates.md`):

| Category | Hypotheses | Key references |
|----------|-----------|---------------|
| **C1 cis-mQTL** | 3 (GoDMC hap-ASM / GTEx multi-tissue / BC-PRS) | Min 2021 PMID 34493871, Oliva 2022 PMID 36510025, Ho 2021 PMID 33768416 |
| **C2 Cancer ASM** | 3 (Do/Tycko / BRCA1 hyper TNBC / CTCF motif) | Do 2020 PMID 32594908, Glodzik 2020 PMID 32719340, Loyfer 2025 NComms |
| **C3 Allele-imbalance** | 3 (CN-collateral / imprinted CN-driven / LOH+gain dose) | Tomita 2017 PMID 29069713, Joshi 2017 PMID 28883545 |
| **C4 Repeat/SD** | 2 (RepeatMasker multi-mapping / ONT basecalling) | Sigurpalsdottir 2025, Vollger 2023 T2T SD |
| **C5 Replication timing** | 2 (PMD methyl loss / chr8/17 timing-driven) | Endicott 2022 PMID 36347867, Du 2021 |

**5 multi-mechanism anchor cells**:
- `Inner|other|cov_proxy_mid` (Z-CHR8 n=328 FP_rate 0.277): C3 + C5
- Z-AUTO `Inner|other|cov_proxy_high` (n=107 FP_rate 0.523): C4 + C5
- `Inner|other|cov_normal` (n=4,984 FP=171): C1 + C2
- `Outer|cross_het|cov_elevated` (Z-OCH n=119 all TP): C2 + C3
- Z-GL `Inner|cross_het|cov_gain`: C3 + C2.2

**H4 verdict (POSITIVE)** — 5 categories 全 ≥1 候選 + ≥1 peer-reviewed prior literature；relaxed gate 通過

---

## 5. Critical Caveats

### 5.1 ΔF1 marginal (highest priority)
- **+0.00242 < +0.005** Cohen ribbon "marginal"
- 歷史對照 Phase 1A paired-pure methyl+context = +0.0112 (5 樣本 CI [+0.003, +0.020])
- **不可宣告 production-ready filter** — 需 cross-sample 才能 claim

### 5.2 Sample shrinkage
- Model B convergence 22% (30/138) — 108 cells saturated (n_FP<5)
- Step 2 gate: 4/11 cells 通過 (Model B converged + n_fit≥100 + n_TP≥5)
- 4 cells 不足以宣告 framework completeness

### 5.3 In-sample optimism (已 mitigated)
- In-sample AUC 0.97 vs CV AUC 0.90 (~7 點差)
- Step 3 用 CV out-of-fold → ΔF1 +0.00242 是 CV-adjusted

### 5.4 Counts 雙口徑警示
- **Plan-spec (post-ISM region-level)**: TP=30,490 / FP=4,842
- **Caller variant-level (09_V6_caller_F1_verification)**: TP=28,509 / FP=11,606 / FN=10,938
- 兩組數字不一致！Step 3 用 plan-spec，但 paper 必須誠實標 framing

### 5.5 Single-sample only
- HCC1395 clairs_to_ssrs only
- 4 cells × 1 sample = 4 cell-sample pairs，**不能宣告 generalizability**

### 5.6 Mechanism unverified
- H4 relaxed gate — 13 hypotheses 都未 cycle-specific 驗證
- 哪個 mechanism 真正驅動 LR β？未測

---

## 6. Tier Evaluation → **⭐3 PARTIAL POSITIVE**

| 維度 | 評估 |
|------|------|
| 樣本數 | 1 (HCC1395, single-sample pilot) |
| 跨樣本 | 未做（留 phase 2） |
| Pre-reg falsifiability | ✅ 5 H 全 pre-specified |
| Confound guard | ✅ within-cell within-AF OLS + 5-fold CV |
| Effect size | ⚠️ ΔF1 marginal (+0.00242 < +0.005) |
| Mechanism | ✅ 13 候選 + 14 PubMed refs (relaxed gate) |
| Prior art differentiation | ✅ TumorLens/ROCIT/SGZ/Wakhan/SAVANA 無同口徑 |

**⭐3 而非 ⭐4 的理由**：
- ❌ Single-sample 不可 generalize
- ❌ Effect size marginal < +0.005
- ❌ Mechanism 候選但未驗證
- ✅ 5 H 全 POSITIVE + 4 cells CV AUC 0.80-0.93 + V5/V6 一致
- ✅ Filter direction signal +0.633 確認

升 ⭐4 必要條件 (見 §7):
1. Cross-sample validation (H1437/H2009/HCC1954/HCC1937 V3F+V5+V6 ISM)
2. Mechanism cycle-specific 驗證
3. ΔF1 提升到 ≥ +0.005 marginal threshold

---

## 7. Future Direction — Phase 2 Plan

### 7.1 Phase 2 cross-sample pilot (升 ⭐4 必經)
- 對 4 樣本 (H1437/H2009/HCC1954/HCC1937) phaseD V6 + 補 V3F/V5 ISM 12 runs each (~80 min × 4 = ~5 hr 平行)
- 在 cross-sample master TSV 上重做 Step 1-3
- Decision rule: ≥4/5 樣本 ΔF1 > 0 + Wilcoxon p<0.05 → ⭐4 升級
- COLO829 deferred (truth set 0600 權限)

### 7.2 Mechanism cycle-specific 驗證
- 對 H4 5 categories 至少 1 個做 dedicated cell-level mechanism analysis
- 例：chr8|Outer|other|cov_normal × cis-mQTL GoDMC 對照（用 mcp__knowledge / WebFetch 跑）
- 對 chr8 region in HCC1395 cohort，比對 Min 2021 GoDMC mQTL hotspot positions

### 7.3 與 CURRENT_FOCUS T2.1 整合
- T2.1 Z-AUTO KDE 跨 4 樣本擴展 — 與本 cycle phase 2 共用 V3F+V5+V6 cross-sample ISM
- 一次跑 cross-sample ISM → 同時餵 T2.1 (Z-AUTO recur) + phase 2 (filter generalize)

### 7.4 Reactive — 若 phase 2 NEGATIVE
- 啟動 T4.2 GC/mappability/repeat 新軸 pilot（依 v0.3 framework 63% gap)
- 不投入更多 methylation 路徑

---

## 8. Paper Framing Recommendation

### 8.1 主軸保守 framing

> **"Read-level methylation features add LRT-significant discriminative information (q<0.05 in 16/30 cell-combos) on top of LOH × HP × CN framework. In an aggressive cell-level filter pilot on HCC1395 clairs_to_ssrs, this yielded ΔF1 = +0.0024 vs ClairS-TO caller baseline F1=0.7166 — directionally POSITIVE but effect size below clinical relevance threshold (+0.01). Cross-sample validation required before production claim."**

### 8.2 §3 主圖建議

- Fig 3a: H1 LRT q heatmap (cells × BAM × covariate)，highlight 16/30 q<0.05
- Fig 3b: ΔF1 vs τ curve (per cell + aggregated)，標 τ*=0.52
- Fig 3c: FP_removal vs TP_loss scatter，diagonal y=x reference (show 4 cells 在 y>x 區)
- Fig 3d: H4 mechanism tree (5 categories × 13 hypotheses) + anchor cell labels

### 8.3 §5 Limitations 必寫
1. Single-sample HCC1395 clairs_to_ssrs only
2. ΔF1 +0.0024 < +0.005 marginal Cohen ribbon
3. Model B convergence 22% (small-FP cell saturation)
4. Counts 雙口徑 (post-ISM region-level vs caller variant-level)
5. H4 mechanism 候選但未 cycle-specific 驗證

### 8.4 與 prior art 差異化（依 02_prior_art_notes.md 已驗證）
- TumorLens: sample-level，不做 per-region filter
- ROCIT: per-read methylation transformer，不做 variant-level filter
- SGZ: 4-axis variant-level filter，無 methylation
- **本 cycle 差異化**: first read-level methylation-augmented multi-axis (LOH × HP × CN × AF + 5 methyl) variant-level filter — prior art 無同口徑

---

## 9. Reproducibility — File Inventory

### 9.1 計畫與方法
- Plan: `/bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md` (v1.0 + Step -1 amendment)
- Plan copy in research dir: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/00_PLAN.md`

### 9.2 Step deliverables
所有 deliverable 於 `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/`:
- `step0_schema.md` (35,332 × 202 cols full schema)
- `step5_master_augmented.tsv` (25.3 MB)
- `step1_lrt_per_cell.tsv` (138 LR fits, BH-FDR q-value)
- `step1_findings.md` (H1 + H5 verdicts)
- `step2_filter_sweep.tsv` (540 evals)
- `step2_findings.md` (H3 direction POSITIVE)
- `step3_delta_f1.tsv` (430 rows × ΔF1 per τ)
- `step3_optimal_tau_summary.md` (Strategy 1+2 winner)
- `step3_findings.md` (H2 + H3 verdicts)
- `step4_mechanism_candidates.md` (5 categories × 13 hypotheses × 14 PubMed refs)
- `step5_findings.md` (Coordinator synthesis, 完整版)
- `figures/`: step1_lrt_heatmap, step2_roc_per_cell, step2_filter_signal, step3_deltaf1_vs_tau, step3_filter_signal_curve
- `scripts/`: build_augmented_master, augmented_lr, filter_sweep, delta_f1, _common_step1

### 9.3 NEW ISM data
- `InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/` (12 ISM runs, ~25 GB)
- Script: `InterSubMod/research/paired_priority_bug_audit/scripts/run_phaseC_genome_three_way_with_significance.sh`

### 9.4 Upstream reference
- v0.3 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md`
- v0.3 standalone HTML: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.standalone.html`
- 09_V6_caller_F1_verification.md: `InterSubMod/research/paired_priority_bug_audit/09_V6_caller_F1_verification.md`
- prior art: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/02_prior_art_notes.md`

---

## 10. Decision Audit Trail

| Step | Decision | 來源 |
|------|----------|------|
| Plan v1.0 design | 4 user decisions (τ 0.5-0.95 max ΔF1, +PubMed, 3 BAM 完整, manual SoT) | 2026-05-18 AskUserQuestion |
| Step -1 amendment | Add phaseC rerun after blocker discovery | 2026-05-18 Agent F2 escalation + user Option A confirm |
| 12 runs setting | 12 全跑 (V3F+V5+V6 × on/off × tp/fp) | 2026-05-18 user choice |
| Final decision Option A | ⭐3 candidate + cross-sample phase 2 | 2026-05-18 user confirm post step5_findings.md review |
| Manual SoT update | INDEX + CURRENT_FOCUS + evidence_ledger | 本 commit |

---

## 11. Coordinator Notes

- Multi-agent fan-out (Agent F, F2, G, H, I + Coordinator) 合計 ~30 min ISM rerun + ~30 min 後分析 = ~1 hr wall clock (主 cycle，不含 Step -1 80 min)
- Plan v1.0 全 Pre-reg 落地，無 NO-GO 觸發
- 上游 phaseC blocker (--no-distance-matrix) 由 plan v0.3 設計時未發現 → 加入 known-pitfalls 追加 P-15 候選
- ΔF1 +0.00242 是 single-sample marginal — 對 PI lab meeting framing 必須誠實標示「需 cross-sample 才能升 ⭐4」
