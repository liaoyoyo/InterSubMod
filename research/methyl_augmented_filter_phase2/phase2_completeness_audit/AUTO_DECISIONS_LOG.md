# 自動決策紀錄 (AUTO_DECISIONS_LOG)

> 用戶 2026-05-21 指令：「持續自動完成直到所有都有結論與結果，有任何要決策的都自動選用預設方式，並紀錄註記讓我最後可以知道理解」
> 
> 本檔記錄整個 Phase A0 → Phase D 全自動執行過程中**所有自動決策**，供用戶最終 audit 理解。

## 執行模式: 全自動 (auto)

遵循 `InterSubMod/.claude/CLAUDE.md §1`:
- C0/C2/C3 由 AI 預設快速通過
- C1/C4 (weekly-report) 仍必停
- Hard Gate (刪檔/C++ commit/NO-GO 判定/覆寫 evidence_ledger) **永遠必停**

## 工作流程

| Phase | Status | 預估時間 |
|---|---|---:|
| Phase A0 | in_progress | 1-2 hr |
| Phase A1 | pending | 30-60 min |
| Phase A2 | pending | 2 hr |
| Phase A3 | pending | 2-3 hr |
| Phase A4 | pending | 3-4 hr |
| Phase A4-Ext | pending | 1-2 hr |
| Phase A4-F1Audit | pending | 1 hr |
| Phase B | pending | 2-3 hr |
| Phase C | pending | 4-6 hr |
| Phase D | pending | 1-2 hr |

## Auto-Decisions Log (chronological)

### 2026-05-21 — Phase A0 啟動

**D001**: Multi-agent 分批策略
- 11 artifacts → 6 個 Explore agents 並行 (而非 11 個)
- 理由: 控制 main agent context 壓力，每 agent 處理 1-2 個相關 artifact

**D002**: Artifact 分組策略
- Agent 1: 20260511_V6_binary_complete_documentation_01.md + 20260513_V6_Attribution_Errata_01.md (V6 binary 主體文件)
- Agent 2: 20260514_HP_tag_priority_bug_engineering_report_01.md + 20260521_PI_V6_signoff_email_draft_5goal.md (priority bug + PI signoff)
- Agent 3: 20260514_HP_tag_17to1_rootcause_explained_01.standalone.html + 20260514_HP_tag_priority_bug_engineering_report_01.standalone.html (priority bug HTML 2)
- Agent 4: 20260514_Self_Phasing_4Layer_Synthesis_Engineering_Report_01.html + 20260514_V6_vs_baseline_HTML_summary_01.html (4-layer + V6 summary)
- Agent 5: 20260515_V6_TPFP_HP_LOH_CN_Characterization_01.standalone.html + 20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.standalone.html (ISM characterization 2)
- Agent 6: phase2_engineering_completeness.standalone.html + phase2_pi_verification/phase2_pi_trust_framework.standalone.html (phase2 verification 2)

**D003**: Audit 輸出格式
- 每 Agent 產出 1 個 audit report .md 到 `A0_audits/<agent_id>_audit.md`
- 統一 columns: claim_id / numerical_value / source_path / source_line / cross_ref_ledger / verdict (verified / suspect / unverifiable)
- main agent 合併到 1 個 A0_existing_artifacts_verification.md

**D004**: Audit verdict 預設規則
- **verified**: 數字 grep 對到 source TSV/CSV 原檔
- **suspect**: 數字只在 .md/HTML 內出現，無對應原始 data file
- **unverifiable**: 數字源 file 已被刪除/搬移/找不到

**D005**: Audit 完成後流程
- Main agent 合併產出 A0_existing_artifacts_verification.md (verified / 待查驗 兩欄表)
- 自動進入 Phase A1-A4 並行（task dependencies 已設）

---

### 2026-05-21 — Phase A0 啟動後 BAM inventory

**D006**: BAM 檔案位置盤點（讀取 disk 後）
- HCC1395 baseline tag.bam: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_tagged.bam` (278 GB) ✓
- HCC1395 V6 tag.bam: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_germline_absent_revert/tumor_tagged.bam` ✓
- 4 樣本 V6 tag.bam: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v6_5sample_extension/{HCC1937,HCC1954,H1437,H2009}/tumor_tagged.bam` ✓
- **4 樣本 baseline tag.bam: 不存在** ✗ (MEMORY 中 [project-phase2-cycle1-global-fp-filter] 已記載 "Track B cross-sample DEFERRED")

**D007**: Phase A1/A2 scope adjust（受 D006 限制）
- 原計畫: baseline 4 樣本 + 5 樣本 HP scan/ratio
- **adjust 為**:
  - A1: HCC1395 baseline vs V6 完整 HP 6 值對比 (1 sample) + V6 5 樣本 HP 6 值分佈 (cross-sample stability，無 baseline 對比)
  - A2: HCC1395 baseline vs V6 ALT-only ratio (1 sample) + V6 5 樣本 ALT-only ratio (cross-sample) 
- 預期影響: 「baseline 4 樣本」相關 claim 在最終 HTML 標為「baseline BAM 4 樣本不存在 — known gap」，不可宣稱比較完成

**D008**: Phase A3 5 features 可用 data source
- 已有 `cycle2/data/{HCC1395,HCC1937,HCC1954,H1437,H2009}_master_augmented.tsv` (5 樣本 per-sample master TSV)
- 已有 `cycle2/data/per_bam_master_V6.tsv` (V6 整合)
- 已有 `cycle2/data/per_bam_master_V5.tsv` / V3F (對照組)
- A3 boxplot + 9-cell heatmap 直接讀這些 TSV，**不需重跑 BAM analysis**

**D009**: Phase A4 LOSO 可用 framework
- 已有 `cycle4/loso_validation/scripts/run_loso_cv.py` (LR baseline)
- 已有 `cycle4/loso_validation/scripts/run_loso_hnew{2,4}.py` (feature subset variants)
- 已有 `cycle4/loso_validation/data/loso_cv_results.tsv` (LR baseline 已跑)
- A4 補 DT / RF / XGBoost 三算法（複用 LOSO framework + 替換 model）

**D010**: Audit agents 啟動完成
- 6 個 general-purpose agents 同時 background 跑
- agentId: a4e372afea6384d4c / aad6be39d72cec79d / a51af57877544ea41 / a26d1090d8e5c1194 / a632836ab194bffbc / ac594472c6471b5ae
- 等待 completion notification

---

### 2026-05-21 — Audit Agent 2 完成（首個回報）

**D011**: Agent 2 audit result (priority_bug + PI signoff)
- 50 claims audited / 45 verified / 4 suspect / 1 unverifiable
- **Critical issue** (4 suspect): `20260514_HP_tag_priority_bug_engineering_report_01.md §8.2` 4-sample ratio table **row-mislabel**：
  - H1437 Doc=0.611 vs Source TSV=1.243 ❌
  - H2009 Doc=1.243 vs Source TSV=0.901 ❌
  - HCC1954 Doc=1.014 vs Source TSV=0.958 ❌
  - HCC1937 Doc=0.978 vs Source TSV=0.611 ❌
  - Aggregate range "0.611-1.243" 正確，但 per-row pair 錯
  - Source: `InterSubMod/research/paired_priority_bug_audit/phaseD_v6_5sample/v6_cross_sample_summary.tsv`
- **20260521_PI_V6_signoff_email_draft_5goal.md (Doc 2): 22 claims 全 verified** ✓ — 可直接給 PI

**D012**: §8.2 row-mislabel 處理決策
- **不直接修改 20260514 validated 報告**（避免動 validated 報告）
- 在 PI 4-goal HTML 中用「**Source TSV 重新引用**」，並加 footnote 註記 Doc 1 §8.2 mislabel
- 後續 weekly-report 母稿 §17 教授可能問答需準備此 mislabel 解釋

**D013**: Agent 6 audit result (phase2_engineering_completeness + phase2_pi_trust_framework HTML)
- **Verdict: PASS** — both HTMLs PI-distributable as-is
- LOSO numbers verified row-by-row vs `cycle4/loso_validation/data/loso_cv_results.tsv` ✓
- Cycle 1 +0.02236 / 9.24× / τ=0.39 verified vs `cycle1/cycle1_findings.md` ✓
- 14 hypothesis verdicts cross-reference correctly to ledger ✓
- H_NEW_2 / H_NEW_4 numbers verified vs `loso_hnew{2,4}_results.tsv` ✓
- minor: drop value +0.02248 vs raw +0.02247 (1e-5 rounding, immaterial)
- **No fabricated numbers / no tier inflation / LOSO reframe consistent**

**D014**: Agent 5 audit result (20260515 V6 TPFP Characterization + 20260519 V6 vs baseline HCC1395)
- **Verdict: BOTH PASS** — 100% metrics verified, no retraction risk
- 20260519 ~110/110 metrics PASS (marker count/rate/HP ratio/hp=3/NG_on=2 purity/F1 sweep/counterexample/LR substitution/per-CpG cross-tab)
- 20260515 ~14/14 H1a-H7 metrics PASS (NEGATIVE Δgap/Z-OCH 0.124×/Z-zone/cross-sample 5/5)
- **校正**: master TSV 實測 **112 cols** (plan 與 prompt 寫 116 cols 是誤；HTML 內 112 正確)
- 核心數字 100% 對齊：+52.4% / +1.26pp / ×13.2 / 1.696→1.609 / hp=3 ×13.2 / NG_on=2 0.8285 全 verified
- 13/13 source artifacts (master_four_way.tsv 11MB / master_three_way.tsv 9.3MB / step1-4 JSON / figures dir) 全存在
- 兩份 HTML 主動揭露 limitations + counterexample (chr12 / AF<0.1 / LOH=NA boundary 350 lost TP) — transparency 極高

**D015**: Plan 修訂備註（不重寫 plan file，留紀錄）
- Plan 寫 "HCC1395 4-way ISM 116 cols × 35,332 regions" → 校正為 **112 cols × 35,332 regions**
- 後續 HTML target2 / weekly-report 母稿 採 112 cols

**D016**: NO-GO 方向確認（Research Direction Guard 提示）
- Guard 提示 TO Germline FP 鑑別 NO-GO (G1-G7) 已關閉
- 本次工作為 **V6 BAM tagging 改進 + ISM characterization + LR cross-sample 方法論 audit**
- **與 G1-G7 NO-GO 方向不同**：不是試圖用 ML 區分 germline FP；而是 V6 修 priority bug 改 HP tag → 影響下游 ISM 計算
- 繼續執行

**D017**: Agent 3 audit result (20260514_HP_tag_17to1_rootcause + priority_bug_engineering_report HTML)
- **Doc 1 (17to1 rootcause HTML)**: 21/22 verified, 1 SUSPECT
  - §4 TL;DR "VCF GT2 約 95% 偏 `1|.`" 把上游 L1 比例 (1.77:1 = 63.9%) 跟下游 BAM HP1 占比 (94.6%) 混淆
  - 屬 **pre-erratum framing**（已被 4-Layer Synthesis Engineering Report §1 line 488 釐清為 26,436/14,931 = 1.77:1）
  - 所有 headline 數字 (17.3:1, 94.6%, 752, 34,855) 都對；只有「~95% GT2 偏 1|.」這句技術上不準確
- **Doc 2 (priority_bug engineering HTML)**: 22/22 verified, SUPERSEDED-banner present (被 4-Layer Synthesis 取代)
  - 4-Layer 分解 17.3 = 1.77 × 9.8 嚴格 source-traced
  - 唯一 minor: "V5 baseline ~1.86:1" 是 1.863/1.838 取平均 acceptable
- **No fabrication detected**, Cross-doc consistency 6 shared claims 全對, Risk LOW
- 引用要點: PROV-V5-002 / -003 / -SYNTH-001 / -PAIRED-001/002 / -V6-001 各 ledger entry 全有

**D018**: Doc 1 §4 TL;DR suspect 處理決策
- **不直接修改 validated HTML**（保 validated 狀態）
- PI 4-goal HTML target1 採 4-Layer Synthesis §1 line 488 (1.77:1 = 26,436/14,931 = 63.9% L1 ratio) 作為「VCF GT2 分佈」描述源
- 文字避免「~95% GT2 偏 1|.」這句；改用「L1 assignment 1.77:1 (~64% : 36%)」清楚分層敘述

**D019**: Agent 1 audit result (V6 binary main docs: 20260511 complete documentation + 20260513 attribution errata)
- 73 claims: **72 verified / 1 suspect / 0 unverifiable / 0 fabrication**
- C2.25 suspect (minor): errata §3a.5b.12 baseline phased PASS = 47,838 vs main doc §8.6.2 + `09_V6_caller_F1_verification.md` = 47,798 (差 40, < 0.1%, 可能 SNV-only filter 差異)
- PI-critical 全 verified:
  - F1 = 0.7166 / 0.6273 跨 baseline/V3F/V5/V6 → 47,798 PASS triple match ✓
  - Phase D 4 樣本 ratio 0.611-1.243 + marker_off_rate 0.817-0.993 + NG_on=2 rate 0.904-0.992 全對 row-by-row to `phaseD_v6_5sample/v6_cross_sample_summary.tsv` ✓
  - V6 vs V3F marker +9.0% = 23,980/21,997 ✓
  - priority bug amplification 9.4× = 17.3/1.85 (5/15 amend) ✓
  - 17,404 unique victim reads ✓
- 5 V6 cycles ledger stability ⭐3-5 全 consistent，無 tier inflation

**D020**: Doc 1 §7 整合 matrix + Doc 2 §3a.5b.7-9 final attribution **可直接引用到 PI HTML target1**
- Phase D 4 樣本數字 **優先從 source TSV 重新引用**（而非從 Doc 1 §8.2 — D012 已知 row-mislabel）

**D021**: Agent 4 audit result (4-Layer Synthesis + V6 vs baseline summary HTML)
- 25 claims / 20 verified / 3 suspect / 2 unverifiable / 0 fabrication
- Suspect:
  - "17,404 V3F victim subset" (Summary §5.2) undocumented
  - "Phase N50 ~20,000+" (4Layer §4.1) vs verified **8,109** internal inconsistency
  - "hp=33" dual-population usage: 571 (HAP3 altHaplotype) vs 138,317 (BAM HP:i:33 reads) 同 label 不同 population
- Unverifiable: "+36% pipeline speed", "V5 BAM = V6 BAM 0 diff" (claim 沒 backup TSV)
- Verified 核心: 17.3:1 / 1.84:1 / 1.77:1 (= 26,436/14,931) / 9.8×/9.4× / F1=0.7166 invariant / +99.7% N50 / +9% marker / Phase D 4-sample 0.611-1.243

**D022**: Phase A0 audit 全完成 6/6
- 總 ~273 claims, 254 verified (93%), 17 suspect (6%), 2 unverifiable (1%), **0 fabrication**
- A0_existing_artifacts_verification.md unified report 已寫
- 所有 PI-critical 數字 verified
- Mark Phase A0 task #1 completed

**D023**: Phase A1-A4 並行啟動 (2026-05-21)
- 5 個 general-purpose agents fan-out background
- A1: HCC1395 baseline vs V6 HP 6 值 + V6 5 sample HP 分佈
- A2: HCC1395 baseline vs V6 ALT-only HP ratio + V6 5 sample ALT-only ratio
- A3: 5 features boxplot + 9-cell heatmap (cycle2 master TSV 已有 5 sample data)
- A4: LR/DT/RF/XGBoost LOSO benchmark (cycle4 framework 重用)
- A4-Ext: 其他算法盤點 + 預期結果（不跑只盤點）

**D024**: BAM 處理策略
- 278GB BAM 直接 scan 慢 → 允許 agent 自行 sampling strategy
- Recommended: pysam sample by region (chr19 / chr8 / chr1 subset) 或 samtools view -s 0.1 sampling
- 全 BAM HP 計數可用 `samtools view -h | awk` pipe (不用全 load memory)

**D025**: Phase A4-Ext 完成 (首個回報)
- 6 算法 vote ranking:
  - Top 1 (3/5): Threshold-only τ-sweep + Per-zone LR (HCC1937) — 合計 < 3 hr cycle 5 補測
  - 中段 (2/5): Boolean rule / Per-AF-bin LR — defer
  - 最低 (1/5): Ensemble + Calibration — 理論失敗 (0+0=0 / F1 monotone-invariant)
- **根本結論**: 6 算法都無法翻轉 LOSO sample-level ≈ 0 finding；root cause 在 distribution shift (sample-level circularity)，不在 model bias
- 真正 pivot 在 read-level (cycle 5+) 或 phasing-signature track，不在 algorithm-space
- 為 PI 報告 §4「為何不試 X」實證 anchor:
  - Ensemble: 兩個 0 base voting lower bound ≈ 0
  - Calibration: F1 monotone re-mapping invariant
  - Per-AF-bin: caller_af 是 collider (L2 collider bias)
  - Boolean rule: cycle 3 已試 +0.01499 marginal
- Artifacts: InterSubMod/docs/experiments/in_progress/2026/05/20260521_A4_Ext_other_algorithms_inventory_01.md + A4_Ext_algorithm_inventory.tsv + vote_ranking.png

### 2026-05-21 — Phase A3 完成

**D026**: Phase A3 result (5 features TP/FP boxplot + 9-cell heatmap + direction consistency)
- 5 features × 5 sample analysis complete (45 row TSV)
- **AUC range 0.010-0.923 (median 0.533)** across 9 features × 5 samples
- **Sign-consistent ≥4/5 → 5/9 features**:
  - `loh_inner_flag` 5/5 ✓ 最強 anchor
  - `NG` (NGroups)
  - `HPFineF` (methylation sub)
  - `NME_imbalance` (methylation sub)
  - `ClusterPermanovaF` (methylation sub)
- **Pairwise Spearman ρ on |Cohen's d| ranking**: -0.07 to 0.97 (median 0.44)
  - H1437 × H2009 ρ=0.97 (同 NCI-H lineage)
  - HCC1954 outlier ρ≈0
- **9-cell TP-rate range 0.22-1.00**, 最強 FP 集中 cell = HCC1395 outer × NG2 × AF_L (TP-rate=0.22, n=3,034)
- **caller_af direction sample-dependent**: HCC1395 AUC=0.92 vs HCC1937 AUC=0.20 (purity-driven)
- **H1437 Coverage_Multiple AUC=0.01** 真實反向 (FP n=8 集中 amplicon high-cov)
- Tier 自我約束: ⭐2 descriptive; 升 ⭐3 需 within-group OLS + AUC confound-guard

**D027**: Direction Guard 確認 (Fine-Pairwise NEGATIVE)
- A3 不使用 fine-pairwise distance；是 5 features TP/FP boxplot + 組合圖
- 不衝突繼續

**D028**: 跨樣本一致性 strongest features
- `loh_inner_flag` 5/5 sign-consistent 是 PI 4-goal HTML target3 核心 anchor
- caller_af direction inconsistent **正式 evidence** for cycle 4 LOSO conclusion

**D029**: Artifacts (A3)
- InterSubMod/docs/experiments/in_progress/2026/05/20260521_A3_5features_distribution_combinatorial_01.md
- TSV: A3_per_feature_AUC_cohend_5sample.tsv / A3_9cell_heatmap_data_5sample.tsv / A3_direction_consistency_spearman.tsv
- 5 figures (boxplot / methyl sub / AUC bar / 9-cell heatmap / direction consistency heatmap)

---

### 2026-05-21 — Phase A4 grid 完成 + nuanced verdict

**D030**: A4 grid 跑完 100 folds (4 algo × 5 sample × 5 seed)
- TSV: A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv (100 rows)
- Summary JSON: data/A4_summary.json
- Figures: A4_4algorithm_LOSO_dF1.png + A4_train_test_gap.png

**D031**: **Nuanced verdict** (修訂 A4-Ext 的 H_circularity 全成立結論)

```
Per algorithm overall mean (cross 5 sample, 5 seed):
  LR:  -0.00004 (0/5 > +0.005)         [原 H_circularity confirmed]
  DT:  +0.00255 (1/5: HCC1395 +0.01378) [新 H_method partial]
  RF:  +0.00267 (1/5: HCC1395 +0.01269) [新 H_method partial]
  XGB: +0.00102 (1/5: HCC1937 +0.00562) [新 H_method partial]
```

**新發現**:
- **LR 跨 algo 不是 100% model failure** — DT/RF 確實在 HCC1395 hold-out 比 LR 多 capture +0.013 F1
- **但僅 1/5 sample 超 Cohen ribbon** — 4/5 sample 仍 LOSO ≈ 0
- 非線性 model 可在「中低 caller-F1 sample」(HCC1395 0.72, HCC1937 0.37) 找 small boost
- 高 caller-F1 sample (0.83-0.89) 仍全 algo ceiling-bound
- **Production filter 仍 NOT deployable** — 無 algorithm 5/5 positive

**D032**: A4-Ext (D025) verdict 部分修訂
- 原 D025 寫「6 算法均不能翻轉 LOSO sample-level ≈ 0」過於 absolute
- 修訂：「6 algo 不能 5/5 positive；DT/RF/XGB 在 single-sample 找 small boost (HCC1395 / HCC1937)，但不 generalize」
- root cause 仍是 distribution shift (sample-level circularity) + caller-F1 ceiling 雙因素

**D033**: 啟動 A4 markdown 補完 agent (a2f959dc8f342d5bf)
- 補 §3 / §4 / §5 pending 部分
- 反映 nuanced verdict
- 為 PI 報告 target4 提供準確結論基礎

**D034**: 啟動 A4-F1Audit agent (a136a3590a1385b0a)
- F1 step-by-step audit trail
- data leakage / overfit / cross-algo sanity check
- PI Q&A 預備 5-7 個

**D035**: Phase C 專案資料夾結構已建 (預備, 不依賴 A 完成)
- InterSubMod/docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/
- assets/{figures,data,js,css}/ + README.md 已建

---

### 2026-05-21 — Phase A4-F1Audit 完成

**D036**: Phase A4-F1Audit F1 6-step audit trail
- Output: InterSubMod/docs/experiments/in_progress/2026/05/20260521_A4_F1_step_by_step_audit_01.md (7 sections, 390 lines)
- TSV: A4_F1_step_audit_trail.tsv (22 audit rows)

**6-step audit verdict: PASS**
- Step 1 Caller VCF 47,798 PASS triple-file-identity verified
- Step 2 V6 BAM tagging 5 sample TP/FP from LOSO TSV verified
- Step 3 ISM region: V6 marker +52.4%, hp_ambig ×13.2 對齊 20260519 ✓
- Step 4 Cycle 1 in-dist 10 feature L2 OOF +0.02236, multi-seed std=5e-5 (deterministic)
- Step 5 τ*=0.39 confusion matrix self-consistent (precision 0.9541, recall 0.6030)
- Step 6 LOSO HCC1395 -0.00012, best τ→0.10 → circularity gap **+0.02248 = 100% sample-level effect**

**Leakage audit: NO LEAKAGE ✓**
- run_loso_cv.py:104 train medians compute, :120 apply to test → no leak
- τ sweep on test fold = oracle (over-optimistic upper bound; 即使如此 mean -0.00004 → robust)
- Feature set drop NumReads VIF=217 跨 5 fold 固定 → no per-fold leakage

**Overfit audit: NOT model-specific**
- 4 algo HCC1937 train-test gap 全 +0.48 一致 (LR 0.4841 ≈ DT 0.4854 ≈ RF 0.4845 ≈ XGB 0.4852)
- 是 caller F1 baseline (0.37 vs train 0.83) artifact, 非 model overfit
- LR overall gap +0.0576 ≈ DT/RF → 排除 LR underfit hypothesis

**Cross-algo: H_circularity 確認 + H_method partial**
- 4 algo cross-sample mean ΔF1 全 < Cohen +0.005
- DT/RF rescue HCC1395 +0.0127~+0.0138 是 in-distribution capacity (LR 該樣本 underfit)
- XGB rescue HCC1937 +0.0056 是 baseline F1 低 artifact

**PI Q1-Q4 預備完整**:
- Q1 (DT/RF?) — A4 實證跨 algo 都失敗 5/5 positive
- Q2 (production?) — circularity gap +0.02248 量化
- Q3 (hyperparam?) — 5 seed × 4 algo mean 不變
- Q4 (F1 提升?) — read-level / thread_d Tier 2 / 加樣本 / NEGATIVE archive 4 path 盤點

**D037**: Phase A 進度
- Completed: A0 / A3 / A4 (含 markdown) / A4-Ext / A4-F1Audit (5/7)
- In_progress: A1 (HP 6 值, 5 BAM scan) / A2 (ALT-only ratio, 5 BAM scan)
- 預估 A1/A2 剩 ~60-80 min

---

### 2026-05-21 — Phase A1 完成

**D038**: Phase A1 result (HCC1395 baseline vs V6 + 4 樣本 V6 HP 6 值)

| Sample | Ver | HP1:HP2 | HP33% | HP=11 | HP=21 | HP=33 |
|---|---|---|---|---|---|---|
| HCC1395 | baseline | 3.24 | 0.03% | 79,118 | 23,007 | 630 |
| HCC1395 | V6 | 3.63 | 1.68% | 98,449 | 25,075 | 37,656 |
| HCC1937 | V6 | 2.16 | 0.10% | 1,065 | 962 | 1,630 |
| HCC1954 | V6 | 2.73 | 0.05% | 357 | 593 | 559 |
| H1437 | V6 | 1.08 | 0.53% | 4,862 | 5,637 | 4,138 |
| H2009 | V6 | 1.09 | 4.27% | 45,312 | 41,923 | 49,360 |

**Key findings**:
1. **HP=33 (Layer 1.5 ambiguous handling)** 跨樣本 firing 不一致:
   - HCC1395 baseline 0.028% → V6 1.68% (**60× firing**)
   - H2009 V6 4.27% (最大 firing 跨樣本)
   - HCC1937/HCC1954/H1437 firing 較弱 (0.05-0.53%)
2. **HP1:HP2 ratio (chr8+chr19 subset)**:
   - V6 跨樣本 range 1.08-3.63 (高 heterogeneity)
   - HCC1395 chr8+chr19 V6 比 baseline **更偏 HP1** (3.24→3.63)
   - **與全 BAM 1.696→1.609 反向**！表示 chr8+chr19 (hotspot) 的 V6 fix 行為與全 BAM 不一致
3. **HP=11 (paired-aware somatic)**:
   - HCC1395 baseline 79,118 → V6 98,449 (+24.5%)
   - 全樣本 V6 都 fire (1k-100k+ range)
4. **HP=21 (somatic alt)**: 跨樣本變化大

**Artifacts**:
- TSV: InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A1_HP_6values_5sample.tsv
- Figure: phase2_completeness_audit/figures/A1_HP_distribution.png
- Raw JSON: A1_HP_6values_raw.json

**D039**: A1 finding 對 PI HTML target1 影響
- HCC1395 chr8+chr19 V6 reverse direction 需 explicit caveat："chr8+chr19 hotspot subset behavior ≠ 全 BAM behavior"
- HP=33 firing 跨樣本 unequal 是新 finding，與既有 20260519 HCC1395 only report 互補

---

### 2026-05-21 — Phase A2 完成 (含 fallback)

**D040**: Phase A2 result (ALT-only HP1:HP2 ratio 5 sample)

| Sample | Ver | scope | ALT-only ratio | All-reads ratio |
|---|---|---|--:|--:|
| HCC1395 | baseline | chr8+chr19 | **4.41** | 3.44 |
| HCC1395 | V6 | chr8+chr19 | **0.43** | 3.36 |
| HCC1937 | V6 | chr8+chr19 | 0.43 | 2.03 |
| HCC1954 | V6 | chr8+chr19 | 0.37 | 2.90 |
| H1437 | V6 | chr8+chr19 | 0.72 | 1.06 |
| H2009 | V6 | **chr19** | 0.79 | 1.07 |

**Key findings**:
1. **PI Errata 17.3:1 對齊**: baseline HCC1395 chr8+chr19 ALT-only = 4.41 ≈ PI Errata 4.19，**不到 17.3** → 17.3:1 應是特定 region/AF filter，與全 chr8+chr19 不同口徑
2. **V6 priority bug fix 跨樣本一致**: 5 V6 sample ALT-only ratio 0.37-0.79 **全部 < 1 (HP2-leaning)**
3. **HCC1395 V6 翻轉**: baseline 4.41 → V6 0.43 (10× HP2 reverse) — V6 修補機制確認
4. **all-reads vs ALT-only 差異**: HCC1395 baseline 3.44 (all) vs 4.41 (ALT-only) — 顯示 ALT-only 偏 HP1 更嚴重，符合 priority bug 對 ALT reads 特別有 effect

**D041**: H2009 scope fallback decision
- 原計畫 chr8+chr19 (12,495 SNV → 預估 85+ min, 超 Monitor 60 min timeout)
- Fallback: H2009 改 chr19 only (2,963 SNV → 18 min)
- TSV scope column 標 "chr19" 區分
- PI HTML target1 H2009 row 加 footnote "scope=chr19 only (timeout dodge); 全 chr8+chr19 待後續補測"

**Artifacts**:
- TSV: InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A2_ALT_only_HP_ratio_5sample.tsv (6 rows)
- Figure: phase2_completeness_audit/figures/A2_ALT_only_vs_all_reads_ratio.png
- JSON: A2_ALT_only_HP_ratio_5sample.json

**D042**: Phase A 全完成 (7/7 tasks)
- ✅ A0 (既有素材 audit)
- ✅ A1 (HP 6 值分佈)
- ✅ A2 (ALT-only ratio with fallback)
- ✅ A3 (5 features TP/FP)
- ✅ A4 (4 algo LOSO)
- ✅ A4-Ext (其他算法盤點)
- ✅ A4-F1Audit (F1 step audit)

**下一步**: 啟動 Phase B (weekly-report W1-W7)

---

### 2026-05-22 — Phase B (weekly-report W1-W7) 完成

**D043**: 載入 weekly-report skill v2 (44 個 skill 列表中)，全自動 mode 對齊 SKILL.md §「執行模式」— C0/C2/C3 用 AI 預設快速通過; C1/C4 標記必停但因用戶明示「持續自動」先產 draft 後 user review

**D044**: Raw data 收集 — 11 個 source artifact (per user instruction)，含 CURRENT_FOCUS.md / AUTO_DECISIONS_LOG.md / A0 audit unified / Phase A 5 task .md / A* TSV / evidence_ledger entries #42-51 / PI V6 signoff email draft

**D045**: 主線類型 (C1) = `progress:problem` 混合
- Thread A (主, progress): V6 BAM-level 進展 (目標 1+2)
- Thread B (sub, problem): LR 方法論 reframe (目標 4)
- Thread C (sub, observation): 5 features characterization (目標 3)

**D046**: 內容 4 層分類 (C2) — [F] 17 / [O] 6 / [I] 4 / [U] 5 共 32 筆素材；4 桶分流 PPT Tier 1=8 / 講稿 Tier 2=11 / 備註 Tier 3=9 / 暫存=4

**D047**: 紅旗檢查 (C3 W5) — 4 紅旗修正
- 過度宣稱 ("V6 解決 priority bug") → 加 specific 量化 "5/5 same direction, range 0.37-0.79"
- 顯著無 CI ("DT/RF 顯著 rescue") → 加 Cohen ribbon + 1/5 sample qualifier
- 流水帳 (A0/A1/.../A4-F1Audit 7 件平列) → 改 4 narrative cluster
- 教授視角缺 → §16 加 pivot 三條路 decision branching point

**D048**: 教授問答預測 (C3 W6) — 7 個追問 (Must-Answer 3 + May-Ask 4)；每個都已預備 evidence-anchor 回答 (含 source path/line)

**D049**: 母稿產出 (C4 W7) — 17 段 Layer 0-4 完整；frontmatter 含 main_thesis / report_type / audience_scenario / source_artifacts list

**D050**: Number-source-grep audit (per SKILL.md §10.2 升級) — 母稿 24 個 `[source: <path>]` 標註；user 可 `grep -E "\[source: " master_draft.md` 確認 coverage

**D051**: scientific-rigor 繼承 §2-§7 (高影響 PI report) — §2 證據分級 / §3 Cohen ribbon / §4 DAG / §7 Pre-reg alignment 全跑

**D052**: handoff 4 選 (C4 後) = **D** (default) — 母稿留檔 + next_week_plan.md 加寫；不主動 trigger pptx-build / html-report-build，等 user review

**D053**: 不重複既有 .md 內容 (per user instruction) — 母稿引用 source path + 本週 delta 補測 / audit / reframe + 整合判斷三層；不複製大段 source 內容

**D054**: status_summary 段 (< 500 字) 寫入 `auto_decisions_for_user_review.md` 最後 §STATUS_SUMMARY

**D055**: Phase B 產出 3 個檔案
- `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/master_draft.md` (569 行)
- `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/auto_decisions_for_user_review.md` (238 行)
- `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/next_week_plan.md` (69 行)

**下一步**: Phase B 完成 → 等 user review (handoff option D); 後續可選 `/html-report-build standalone` 接 PI 4-goal HTML build 或 `/pptx-build --from-draft master_draft.md` 接 PPTX

---

### 2026-05-21 — Phase D fresh-context evaluator 完成

**D056**: Phase D evaluator verdict
- **PASS with 4 minor recommendations** (PI-distributable as-is)
- **Tier: ⭐⭐⭐⭐ L2** (multi-sample cross-validated, 0 fabrication, MEMORY-consistent)
- 上限非 ⭐⭐⭐⭐⭐ 因 4 樣本 baseline BAM 缺 (Track B DEFERRED); only HCC1395 paired

**10/10 number spot-check 全 match source TSV**:
- A2 HCC1395 baseline ALT 4.4078 ↔ "4.41" ✓
- A2 V6 HCC1395 0.4252 ↔ "0.43" ✓
- A2 V6 cross-sample 0.372-0.793 ↔ "0.37-0.79" ✓
- A1 HCC1395 baseline HP1:HP2 3.2423 ↔ "3.24" ✓
- A1 H2009 HP=33% 4.2746 ↔ "4.27%" ✓
- A4 LR/DT/RF/XGB overall mean -0.00004/+0.00255/+0.00267/+0.00102 全對齊
- target2 marker +52.4% / +1.26pp / ×13.2 全 verified

**A0 caveat compliance 100%**:
- 不引 20260514 priority bug eng .md §8.2 row-mislabel ✓
- 不引 17to1 §4 "~95% GT2" framing ✓
- N50 用 verified 8,109 ✓
- master TSV col=112 ✓

**4 minor recommendations** (非阻擋, polish 為 lab meeting 級):
1. index.html:363「證實 capacity 存在」過強 → 改「指向 capacity 在 HCC1395 in-distribution 可被 tree-based capture」
2. target2:229「V6 嚴格優於 baseline」需弱化措辭
3. 完備度 85% 9 處 self-declared 缺計算公式 → 加 ✱footnote 註腳
4. target4 §4A.1 加 95% CI column (5 seed bootstrap)

**D057**: 啟動 polish agent 修 4 個 minor (background)
- agentId: 待啟動完成 notification

**D058**: Phase D 主流程結論
- PASS → 主 agent 可推進 PI sign-off email send + Phase 2 cycle close + push v6-prod-20260520 tag (需用戶明示授權)
- 4 minor recommendation 為 polish, 不阻擋 PI 1-on-1 上報

---

## 全自動執行終結

**任務完成度**: 10/10 tasks completed
1. ✅ Phase A0 (既有素材 audit - 11 artifacts, 273 claims, 254 verified, 0 fabrication)
2. ✅ Phase A1 (HP 6 值分佈 - 6 BAM, A1 TSV + figure)
3. ✅ Phase A2 (ALT-only ratio - 6 BAM with H2009 chr19 fallback, A2 TSV + figure)
4. ✅ Phase A3 (5 features TP/FP - AUC + 9-cell heatmap + Spearman ρ, 3 TSV + 5 figures)
5. ✅ Phase A4 (4 algo LOSO - 100 fold, A4 TSV + 2 figures + summary JSON)
6. ✅ Phase A4-Ext (6 算法盤點 - vote ranking)
7. ✅ Phase A4-F1Audit (F1 step trail)
8. ✅ Phase B (weekly-report W1-W7 - master_draft.md 41KB + auto_decisions + next_week_plan)
9. ✅ Phase C (Hybrid HTML 1 index + 4 子 HTML - 29-63KB each, 全 standalone)
10. ✅ Phase D (evaluator fresh-context PASS ⭐⭐⭐⭐ L2)

**Total Auto-Decisions**: D001-D058 全紀錄

**用戶後續需要做的事**:
1. Review `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/master_draft.md` 確認週報母稿
2. Review `auto_decisions_for_user_review.md` 確認 default 選項
3. Review 5 個 HTML (`InterSubMod/docs/presentations/in_progress/20260522_PI_V6_signoff_4goal/index.html` 起)
4. (Optional) 4 minor polish 完成後 lab meeting 級
5. (Hard Gate) 親自 send V6 sign-off email + git push v6-prod-20260520 tag (需用戶明示授權)

---

## (此 log 將隨執行過程持續追加 — END)

---

## (此 log 將隨執行過程持續追加)

---

## (此 log 將隨執行過程持續追加)

---

## (此 log 將隨執行過程持續追加)

---

## (此 log 將隨執行過程持續追加)

---

## (此 log 將隨執行過程持續追加)

---

## (此 log 將隨執行過程持續追加)
