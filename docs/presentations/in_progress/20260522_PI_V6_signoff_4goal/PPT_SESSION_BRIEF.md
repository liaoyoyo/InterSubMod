# PPT Session Brief — 對話重點 + 圖片路徑 + slide-ready 敘述

> **建立日期**: 2026-05-21
> **用途**: 改 PPT session 用素材庫
> **基礎來源**: 5 個 HTML（index + 4 target）+ master_draft.md + 對話教學整理

---

## A. 4 目標主軸 (Slide-ready 敘述)

### Target 1: V6 vs baseline tag.bam 層

**標題建議**: V6 priority bug fix — BAM-level 跨 5 樣本一致確認

**核心敘述（≤ 50 字 / slide）**:
> V6 修補 longphase-to priority bug；HCC1395 baseline ALT-only HP1:HP2=4.41 → V6=0.43（10× HP2 reverse），5 V6 樣本全 < 1 確認跨樣本通用。

**證據卡 (講稿)**:
- baseline ALT-only 4.41 ≈ PI Errata 全 victim 4.19（口徑對齊）
- V6 fix 後 5 sample ALT-only range 0.37-0.79 全 < 1（5/5 sign-consistent）
- All-reads ratio baseline 3.44 → V6 3.36（不變 -0.08）= V6 只動 ALT-in-bug-trap 不污染 REF
- HP=33 (Layer 1.5 ambig) baseline 0.028% → V6 1.68%（60× firing）；H2009 V6 max 4.27%

**Tier**: ⭐⭐⭐⭐ L2

---

### Target 2: ISM region 層 V6 vs baseline

**標題建議**: V6 ISM 量質同進 — NG≥3 marker +52.4% + rate +1.26pp

**核心敘述**:
> HCC1395 4-way ISM (baseline/V3F/V5/V6) 4 軸全勝；V6 NG≥3 marker count 15,738 → 23,980 (+52.4%) 且 rate 0.8967 → 0.9093 (+1.26pp) 量質同步上升，net TP:FP=14:1。

**證據卡**:
- hp=3 ambig reads 10,440 → 138,317（×13.2）= V6 Layer 1.5 開始正確標 ambig
- NG_on=2 purity 0.8570 → 0.8285（-0.0285 trade-off, net 仍 14:1 正向）
- HP1:HP2 全 BAM 1.696 → 1.609（-0.087 offset）
- Caller F1 baseline = V6 = 0.7166（不傷 caller-level）
- master TSV 112 cols × 35,332 regions

**Counterexample (要主動講, transparency)**:
- chr12: ΔF1 -0.0158
- AF<0.1: -0.0304
- 350 lost TP: 100% 在 LOH=NA boundary

**跨樣本 caveat**: HCC1395 paired only（4 樣本 baseline BAM 不存在，Track B deferred）；cross-sample 用 V6-only HPFineN ρ=0.845

**Tier**: HCC1395 paired ⭐⭐⭐⭐ L2 / cross-sample ⭐⭐⭐ L3

---

### Target 3: 5 features 獨立+組合分析

**標題建議**: loh_inner_flag 是最強 cross-sample anchor

**核心敘述**:
> 9 個 features × 5 樣本驗證；5/9 features sign-consistent ≥4/5；**loh_inner_flag 5/5 sign-consistent**（Gap 8-43pp）為最強 anchor。caller_af direction sample-dependent（purity-driven, HCC1395 d=+1.33 vs HCC1937 d=-1.46）。

**loh_inner_flag 量化（實測）**:

| Sample | TP_inner% | FP_inner% | Gap (pp) | AUC | Cohen's d | 等級 |
|---|---:|---:|---:|---:|---:|---|
| HCC1395 | **45.97%** | 4.63% | **+41.34** | 0.71 | 0.88 | ⭐⭐⭐ |
| HCC1937 | **49.49%** | 6.04% | **+43.45** | 0.72 | 0.93 | ⭐⭐⭐ |
| H2009 | 28.07% | 4.92% | +23.15 | 0.62 | 0.52 | ⭐⭐ |
| H1437 | 20.09% | 0.39% | +19.70 | 0.60 | 0.49 | ⭐⭐ |
| HCC1954 | 9.97% | 1.46% | +8.51 | 0.54 | 0.29 | ⭐ |

→ 5/5 sign-consistent ✓；3 strong/moderate（HCC1395+HCC1937+H2009）；2 weak 因 FP base rate 太低（HCC1954 FP n=687, H1437 n=773）

**Boxplot caveat（主動講）**: loh_inner_flag 是 binary {0,1}，boxplot median 全壓 0 一條線並非 feature 失效；正確視覺化 = grouped bar of inner% + AUC/d + 9-cell heatmap

**其他 features 摘要**:
- caller_af direction 反復（5 sample 中 2+/3-）= cross-sample LOSO 失敗主因
- NG (NGroups) 4/5 sign-consistent
- HPFineF / NME_imbalance / ClusterPermanovaF 4/5 sign-consistent

**9-cell heatmap finding**: HCC1395 outer × NG2 × AF_L 是最強 FP 集中區（n=3,034, TP-rate=0.22）

**Tier**: 整體 ⭐⭐ L4 descriptive；loh_inner_flag 5/5 anchor ⭐⭐⭐ L3

---

### Target 4: TO 純樣本能否超 caller F1?

**標題建議**: Production filter NOT deployable — distribution shift 非 model class

**核心敘述**:
> 4 algorithm (LR/DT/RF/XGB) LOSO 100-fold benchmark；4 algo cross-sample mean ΔF1 全 < Cohen +0.005 ribbon。Cycle 1 +0.02236 = 100% sample-level circularity (drop +0.02248)。root cause 在 distribution shift 不在 model class。

**核心數字表**:

| Algorithm | cross-sample mean | best τ | n > +0.005 |
|---|---:|---:|---:|
| LR | -0.00004 | 0.10 | 0/5 |
| DT | +0.00255 | varies | 1/5 (HCC1395 +0.013) |
| RF | +0.00267 | varies | 1/5 (HCC1395 +0.013) |
| XGB | +0.00102 | varies | 1/5 (HCC1937 +0.006) |

**Nuanced verdict**:
- DT/RF HCC1395 +0.013 是 in-distribution capacity（LR 該樣本 underfit）, **不是 generalize ability**
- 4/5 sample 仍 LOSO ≈ 0 → filter 不可 production deploy

**LOSO vs in-distribution reframe**:

| HCC1395 | Cycle 1 (5-fold OOF, 同 sample 內 split rows) | LOSO (4 train + 1 held-out) | Drop |
|---|---:|---:|---:|
| ΔF1 | +0.02236 ⭐⭐ L4 | -0.00012 ⭐⭐⭐⭐ L2 | **+0.02248 (100% circularity)** |

**Train-test gap audit**: 4 algo HCC1937 gap 全 +0.48 一致 = caller F1 baseline shift artifact 非 model overfit

**A4-Ext 額外 6 算法盤點**（為何不試）:
- Ensemble (LR+RF voting): 0+0 voting lower bound ≈ 0
- Calibration (Platt/isotonic): F1 monotone re-mapping invariant
- Per-AF-bin LR: caller_af 是 collider bias
- Boolean rule: cycle 3 已試 +0.01499 marginal
- Per-zone LR / Threshold-only: vote 3/5 仍 cycle 5 optional

**Pivot direction**:
- Path A (推薦): phase_block_3d / thread_d
- Path B: H_NEW_3 chr8 zone gate
- Path C: low-F1 panel 驗 HCC1395 H_NEW_4 +0.00699
- Path D: Read-level epigenetic context

**Tier**: ⭐⭐⭐⭐ L2

---

## B. 教學 takeaway 6 點（slide candidate / Q&A 預備）

| # | Takeaway | 1-句版本 |
|---|---|---|
| 1 | HP tag 6 值意義 | HP=0 unphased / 1,2 paired-germline / 11,21 self-phasing somatic / 33 V6 Layer 1.5 ambig |
| 2 | Priority bug root cause | Read 上無 germline het 可投票時 hardcode default=HP=1（不是 HP=0 unphased）|
| 3 | ALT-only ratio 為何是 bug 指紋 | 95% REF 不受影響 + 5% ALT 受影響 → ALT-only ratio amplify bug 100×；4.41 vs all-reads 3.44 差 +0.97 = 指紋 |
| 4 | V6 fix 跨樣本通用證據 | 5 sample V6 ALT-only 0.37-0.79 全 < 1 + range 接近 0.5 + all-reads 不變 = 點到 bug trap 不亂改 |
| 5 | NG 分群 + marker 量質同進 | priority bug fix → HP=33 暴露 → cluster 多樣 → NG≥3 marker +52.4% 且 rate +1.26pp 量質同步 |
| 6 | Filter NOT deployable + pivot | Cycle 1 +0.02236 = 100% circularity + 4 algo LOSO 全 < Cohen ribbon → pivot read-level / phasing-signature / low-F1 panel |

---

## C. 新圖片路徑（按 4 target 分類）

### Target 1: V6 vs baseline tag.bam 層

| Figure | 路徑 | 用途 |
|---|---|---|
| HP 6 值分佈 (5 BAM stacked bar) | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A1_HP_distribution.png` | 1A slide |
| ALT-only vs all-reads ratio (5 BAM grouped bar) | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A2_ALT_only_vs_all_reads_ratio.png` | 1C slide (核心) |
| (既有可用) 17.3:1 rootcause 4-layer diagram | 從 `InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_17to1_rootcause_explained_01.standalone.html` 截圖 | 1C support |
| (既有) V6 vs baseline summary | 從 `InterSubMod/docs/reports/validated/2026/05/20260514_V6_vs_baseline_HTML_summary_01.html` 截圖 | 1D background |

### Target 2: ISM region 層

| Figure | 路徑 | 用途 |
|---|---|---|
| (既有 4-way 主圖) HCC1395 ISM characterization | 從 `InterSubMod/docs/experiments/in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.standalone.html` 截圖 | 2A slide |
| (既有) V6 TPFP HP-LOH-CN heatmap | 從 `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.standalone.html` 截圖 | 2B-2D |

### Target 3: 5 features 獨立+組合分析

| Figure | 路徑 | 用途 |
|---|---|---|
| **★ NEW** loh_inner_flag grouped bar (TP vs FP inner%) | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A3_loh_inner_flag_grouped_bar_5sample.png` | 3A 核心（替代 boxplot 一條線誤導） |
| 5 features boxplot per sample | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A3_5features_boxplot_per_sample.png` | 3A 連續 features part |
| Methylation 5 sub boxplot | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A3_methyl_5sub_boxplot.png` | 3B |
| 9-cell heatmap (LOH × NG × AF, 5 sample) | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A3_9cell_heatmap_5sample.png` | 3C 核心 |
| AUC + Cohen's d cross-sample bar | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A3_AUC_cohend_bar_5sample.png` | 3D |
| Direction consistency Spearman ρ heatmap | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A3_direction_consistency_heatmap.png` | 3E |

### Target 4: TO 純樣本 F1 方法論

| Figure | 路徑 | 用途 |
|---|---|---|
| 4 algorithm LOSO ΔF1 (LR/DT/RF/XGB × 5 sample) | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A4_4algorithm_LOSO_dF1.png` | 4A 核心 |
| Train-test F1 gap (overfit symptom) | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A4_train_test_gap.png` | 4B sanity |
| A4-Ext 6 algorithm vote ranking | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A4_Ext_algorithm_vote_ranking.png` | 4D |
| (既有可用) phase2 PI trust framework figures | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_pi_verification/figures/` | 4E support |

---

## D. 精煉金句（PPT headline / TL;DR 候選）

按 slide 重要度排序：

| 句 | 用途 |
|---|---|
| **「V6 priority bug fix 跨 5 樣本一致 — ALT-only ratio 全 < 1 ✓」** | Target 1 headline |
| **「Marker +52.4% × Rate +1.26pp = 量質同進，net TP:FP=14:1」** | Target 2 headline |
| **「loh_inner_flag 5/5 sign-consistent 是最強 cross-sample anchor」** | Target 3 headline |
| **「Cycle 1 +0.02236 = 100% sample-level circularity → filter NOT deployable」** | Target 4 headline |
| **「換 DT/RF/XGB 不能解 — root cause 在 distribution shift 不在 model class」** | Target 4 PI Q1 |
| **「V6 修對 bug 不是亂改：ALT-only 翻轉但 all-reads 不變」** | Target 1 deep |
| **「Pivot read-level / phasing-signature / low-F1 panel — algorithm-space 已 exhausted」** | Decisive Next Step |
| **「A0 audit 273 claims / 254 verified / 0 fabrication」** | Trust framing |

---

## E. PI 必問 4-7 個 + 標準答案

| Q | 標準答案 |
|---|---|
| **Q1**: 為何不用 DT/RF/XGBoost？ | A4 100-fold 4 algo benchmark cross-sample mean 全 < Cohen +0.005 ribbon；DT/RF HCC1395 +0.013 是 in-distribution capacity 非 generalize；root cause 在 distribution shift 不在 model class |
| **Q2**: Filter production-deployable? | No. Cycle 1 +0.02236 vs LOSO -0.00012 drop +0.02248 = 100% sample-level circularity; in-distribution case study only |
| **Q3**: hyperparameter 沒調好? | 5 random seed × 4 algo mean 不變 (LR deterministic / DT/XGB deterministic given seed / RF std=0.00041)；換 hyperparam 不能翻轉 sample-level circularity |
| **Q4**: F1 還能怎麼提升? | 4 path: A) phase_block_3d / thread_d, B) chr8 zone gate, C) low-F1 panel 驗 HCC1395 H_NEW_4 +0.00699, D) Read-level epigenetic context |
| **Q5**: 17.3:1 跟 1.696 怎麼差這麼多? | 不同口徑: 17.3:1 是 SP1/2/3 region ALT-only cherry-pick (IGV 113:0 / 109:1 / 108:0)；1.696 是全 BAM all-reads；4.41 (本實驗 chr8+chr19 ALT-only) 跟 4.19 (PI Errata 全 victim) 對齊 |
| **Q6**: loh_inner_flag boxplot 一條線是 bug? | 不是，是 binary feature 預期行為（median=Q1=Q3=0）；正確視覺化是 grouped bar of inner% (Gap 8-43pp 都正向) + AUC 0.54-0.72 + 9-cell heatmap |
| **Q7**: 為何 4 樣本 baseline BAM 不存在? | Track B cross-sample DEFERRED (MEMORY project_phase2_cycle1_global_fp_filter)；建議下週初補跑或暫接受 V6-only cross-sample (ρ=0.845) |

---

## F. PPT 結構建議（20 slide / 20-30 min PI 1-on-1）

| Slide # | 標題 | 重點 figure |
|---|---|---|
| 1 | Title + agenda | — |
| 2 | This Week's Verdict (≤ 50 字) | text only |
| 3 | 4 目標 grid + completion% | 4 card |
| 4 | Target 1 headline | text |
| 5 | Target 1: HP 6 值分佈 | A1 figure |
| 6 | Target 1: ALT-only ratio reverse | A2 figure |
| 7 | Target 1: priority bug 4-layer diagram | 17to1 截圖 |
| 8 | Target 2 headline | text |
| 9 | Target 2: marker +52% × rate +1.26 | 20260519 截圖 |
| 10 | Target 2: counterexample (chr12 / AF<0.1) | 20260519 截圖 |
| 11 | Target 3 headline | text |
| 12 | Target 3: loh_inner_flag grouped bar (NEW) | **A3 NEW figure** |
| 13 | Target 3: 9-cell heatmap | A3 heatmap |
| 14 | Target 4 headline | text |
| 15 | Target 4: 4 algorithm LOSO | A4 figure |
| 16 | Target 4: Cycle 1 vs LOSO reframe | reframe diagram |
| 17 | Target 4: 6 algorithm vote ranking | A4-Ext figure |
| 18 | Decisive Next Step + Pivot 4 paths | text + table |
| 19 | A0 audit transparency | 9-box stat |
| 20 | PI Q&A 預備（7 Q standby） | text |

---

## G. 引用 ground rules（嚴守, 不可違反）

| 應引用 | 不應引用（A0 caveat）|
|---|---|
| ✅ 20260521 PI signoff email (22/22 verified) | ❌ 20260514 priority bug eng .md §8.2 (row-mislabel) |
| ✅ 20260519 HCC1395 4-way HTML (110/110) | ❌ 17to1 HTML §4 "~95% GT2" framing |
| ✅ 20260515 V6 TPFP characterization (14/14) | ❌ 4-Layer "Phase N50 ~20,000+" (用 verified 8,109) |
| ✅ phase2_pi_trust_framework HTML | ❌ "+36% pipeline speed" (unverifiable) |
| ✅ Source TSV (A1/A2/A3/A4) row-by-row | ❌ "V5 BAM = V6 BAM 0 diff" (claim 沒 backup TSV) |
| ✅ master TSV col 112 | ❌ master TSV col 116 (誤) |

---

## H. 完整檔案總覽

### HTML（5 個）
- `InterSubMod/docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/index.html` (491 行)
- `InterSubMod/docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/target1_V6_vs_baseline_tagbam.html` (863 行)
- `InterSubMod/docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/target2_ISM_V6_vs_baseline.html` (690 行)
- `InterSubMod/docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/target3_5features_distributions.html` (1248 行, 含 §3A.1 NEW loh_inner_flag grouped bar)
- `InterSubMod/docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/target4_TO_pure_F1_methodology.html` (1096 行)

### Master Draft + Decisions Log
- `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/master_draft.md` (569 行, 41 KB)
- `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/auto_decisions_for_user_review.md` (238 行)
- `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/next_week_plan.md` (69 行)
- `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/AUTO_DECISIONS_LOG.md` (D001-D058 全紀錄)

### Phase A 補測 .md (6 個)
- `InterSubMod/docs/experiments/in_progress/2026/05/20260521_A1_baseline_vs_V6_HP_6values_scan_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/05/20260521_A2_ALT_only_HP_ratio_口徑對齊_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/05/20260521_A3_5features_distribution_combinatorial_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/05/20260521_A4_multi_algorithm_LOSO_methodology_completeness_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/05/20260521_A4_F1_step_by_step_audit_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/05/20260521_A4_Ext_other_algorithms_inventory_01.md`

### Audit 兩欄表
- `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_existing_artifacts_verification.md`
