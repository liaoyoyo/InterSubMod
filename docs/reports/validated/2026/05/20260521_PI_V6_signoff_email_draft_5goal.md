<!--
build_date: 2026-05-20
agent: V6 production tag PI sign-off email draft, 5-Goal lens
status: draft (awaiting user review + send)
report_class: pi-email-draft
audience: PI (advisor) — V6 production tag sign-off
parent_workflow: InterSubMod/research/selfphasing_v6_production/4day_compressed_workflow.md Day 5
inputs:
  - InterSubMod/docs/experiments/in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md (5/19 主報告 + 5/20 §13)
  - InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md (V6 binary 完整說明)
  - InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md (errata 5 條 + 5/18 V6 補強骨架)
verdict: V6 production tag = 🟢 GO with calibrated caveats; honest 5-Goal lens assessment
last_verified: 2026-05-20
-->

# Email Draft: V6 Production Binary Sign-Off — 5-Goal Validation Summary

**To**: PI  
**Subject**: V6 Production Binary Sign-Off Request — 5-Goal Validation Complete (HCC1395 single-sample + 4 Phase D)  
**Tag**: `v6-prod-20260520` → commit `8a90532` on `fix/pon-only-phasing` (longphase-to-mod, **local only**, awaiting push authorization)  
**Estimated reading time**: 5 min full / 1 min Executive summary only  
**Attachments**: 5 reports + 13 figures + 4 deployable scripts (paths in §6)

---

## 1. Executive Summary (3 sentences)

V6 binary patch (移除 V5 Layer 1.5 over-promote, 重用 V5 phasing layer) 在 HCC1395 全基因組 5 大目標驗證後**主要在 Goal 2 marker filter + Goal 5 hard NG threshold 嚴格優於 baseline**（marker count +52.4%, marker rate +1.26pp, F1 hard ΔF1 +0.001~+0.197 跨 T=1-5）；對 Goal 5 multi-axis LR (+0.023 ΔF1) 為 BAM-independent (baseline LR +0.02302 ≈ V6 LR +0.02236)，Goal 3+4 unchanged。**建議 ✅ GO V6 production tag**，但 sign-off email 明示**不應宣稱 V6 改進 LR-based F1 或 Goal 3+4** — 升級價值在 marker filter downstream 與 hard threshold F1 robustness。

## 2. What changed (4-BAM 演進 + Day 2-5 work)

| BAM | Binary commit | 設計變更 |
|---|---|---|
| baseline | LongPhase-TO upstream | priority bug 17.3:1 偏 HP1 (PI 4-29 確認) |
| V3F | `41ff147` 2026-04-10 | two-layer getVote 修補 priority bug；最 balanced |
| V5 | `938f0df` 2026-04-30 | + Layer 1.5 somatic fallback；but germline-absent 區域繼承 4.19:1 priority bug |
| **V6** | `8a90532` 2026-05-20, tag `v6-prod-20260520` | V5 phasing + V3F-style germline-absent hp=33 conservative |

**Day 2-5 work timeline**:
- 5/19 Day 2: baseline vs V6 head-to-head（filter + F1 sweep + counter-example + LR @ baseline）
- 5/20 Day 3: 5-Goal validation (LR ablation methylation + per-CpG × HP × ALT + HPFineN cross-sample)
- 5/20 Day 4 (1 day early): ✅ commit `8a90532` + tag `v6-prod-20260520` (local only, push pending user authorization)
- 5/21 Day 5: PI email send (本檔，🔴 Hard Gate awaiting user copy to mail client)

## 3. 4-BAM × 5-Goal Verdict Table（簡版）

| Goal | Best BAM | V6 升級價值 |
|---|---|---|
| 1 per-CpG × HP × ALT | V3F (imbalance 0.275) | 🟡 partial (V6 0.377 比 baseline 0.446 好) |
| 2 marker filter | **V6** (+52% count) | ✅ Strong dominance |
| 2 cross-sample (V6×5) | V6 (ρ 0.845) | ✅ Valid (HCC1937 outlier documented) |
| 3 二次打擊 | blocked | ❌ HP tag 改動不解決 (需 per-read redesign) |
| 4 TO normal 補強 | blocked | ❌ caller_af 0.654 已主導 |
| 5 F1 hard T=1 | V5/V6 tie | ✅ marginal +0.0004 (Cohen ribbon 內) |
| 5 F1 hard T=2-3 | **V6** | ✅ ΔF1 V6-baseline +0.03~+0.16 |
| 5 F1 LR full | tie (baseline+0.001) | ⚠️ BAM-independent |
| 5 F1 LR no-meth | tie | ⚠️ methylation contrib 僅 +0.0005 |

**Sharp findings (5/19-20)**:
1. baseline ≈ V5 在 HP × ALT confound (Cramer V 0.1068/0.1069 identical) — V5 Layer 1.5 沒實質修補
2. baseline LR ΔF1=+0.02302 略勝 V6 LR +0.02236 — Cycle 1 +0.02236 NOT V6-specific
3. ISM methylation features 對 ΔF1 貢獻僅 +0.0005-0.0007（5th-rank coef，去掉後幾乎不變）
4. ISM 對 F1 提升的核心 leverage 是 caller_af (+3.44 coef) — caller 已知特徵
5. **Source code dual bug 鐵證** (5/20 對話澄清): baseline 同時有兩個獨立 source-code-level bug — Bug 1 `getVote` priority order (line 510-513, somatic 倒置在 germline 之前 → HP1:HP2 = 17.3:1 偏移)；Bug 2 `judgeHaplotype` dead code branch (line 697-701, enum value 3/4 vs hpResult mapped integer 0/1/2/11/21/33 永遠不等 → low-confidence fallback 9 年 dead code)。V3F commit `41ff147` two-layer rewrite + integer literal fallback **同時修補兩者**；V5/V6 繼承此修補。baseline 10,440 hp=3 reads 是 getVote 第 2 priority pair (HAPLOTYPE3 winning) 的小邊界 case (0.42% reads)，**不是 Bug 範圍**。Util.h:19-26 enum + ISM ReadParser.cpp:130-141 mapping (HP:i:33 BAM → reads.tsv "3") 鐵證已寫入主報告 §6.4 + errata §3.4 / §5.4a。

## 4. V6 Production Tag Recommendation: **🟢 GO with calibrated caveats**

| 升級理由 | 充分性 |
|---|---|
| Marker filter +52.4% count | ✅ Strong |
| F1 hard T=2-3 strict dominance | ✅ Strong |
| Cross-sample V6 一致性 ρ 0.845 | ✅ Valid |
| Caller F1 不變 (五階段相同 0.7166) | ✅ Verified |
| V6 = V5 phasing + V3F-style hybrid 證實 | ✅ Confirmed |

| 不應 over-claim | 證據 |
|---|---|
| V6 改進 LR F1 | ❌ BAM-independent (5/20 step 1+2) |
| V6 改進 Goal 3+4 | ❌ HP tag 改動不解決 |
| V6 ≥ V3F 跨所有 use case | ❌ Goal 1 per-CpG V3F 仍 best |

**推薦行動**：
1. ✅ 執行 `git tag v6-prod-20260521`（marker filter + cross-sample ISM downstream production use）
2. ⚠️ Goal 1 per-CpG 純研究 publication 仍考慮用 V3F BAM（lower HP × ALT confound）
3. ⚠️ Goal 5 cross-sample F1 deployment 需 caller-F1-headroom-aware (Phase 2 Cycle 3, 後續工作)
4. ⚠️ Goal 3+4 unblock 需 ISM 框架 redesign（per-read epigenotype + methyl-phasing），不在 V6 升級範圍

## 5. Open Questions + Next Phase

| Open Q | 規劃 |
|---|---|
| Goal 5 跨樣本 LR 失敗 mechanism (Cycle 2 caller-F1-headroom-bounded) | Cycle 3 redesign：caller F1 < 0.80 + FP density > 0.10 gate；plan v2 R-MENTAL-DRIFT 紀律後啟動 |
| Goal 1 per-CpG × HP × ALT 在 V3F vs V6 trade-off | 後續：V3F dedicated 跑 per-CpG 純研究 (~2 hr ISM rerun) |
| baseline LR > V6 LR 原因 (NG distribution low + better feature decoupling?) | 後續：VIF 對比 baseline vs V6 LR feature space |
| ALT-only HP1:HP2 ratio 重現 PI 4-29 17.3:1 (本實驗 all-reads 1.696) | 後續：~2 hr 重跑 ISM with alt_support filter |
| COLO829 V6 (本輪 skip due to truth-set permission) | Phase 2 V6_phase2: chmod 660 or 替代 truth set |

## 6. Attachments (引用 artifacts)

### 主報告
- `InterSubMod/docs/experiments/in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md` — Day 2-3 完整報告 (13 sections, 含 §13 4-BAM × 5-Goal cross-tab)
- `InterSubMod/docs/experiments/in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.standalone.html` — PI-friendly HTML 80 KB (13 sections, 13 figures, sticky TOC)

### V6 binary 完整說明
- `InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md` — V6 binary + 4-phase 驗證 (chr19+全基因組+4 樣本 cross-sample)

### Errata 補強
- `InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md` — PI Report 4-29 errata 5 條 + 5/18 V6 補強骨架 (§2.4 / §4.4 / §5.6)

### Data + Scripts (deployable)
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv` (35,332 × 112)
- `InterSubMod/research/methyl_augmented_filter_phase2/cycle1/cycle1_track_a_filter.json` (Cycle 1 LR deployable filter)
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/scripts/{baseline_vs_v6_analysis,plot_baseline_vs_v6,baseline_vs_v6_counterexample,baseline_vs_v6_f1_analysis,plot_f1_analysis,apply_lr_to_baseline,lr_ablation_methylation,per_cpg_hp_alt_4bam,hpfine_cross_sample_consistency}.py` (9 deployable scripts)

### Figures (13 total in `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/figures/baseline_vs_v6/`)
- F1 global summary 4-panel
- F2 per-chr HP ratio
- F3 marker rate by LOH zone
- F4 HP family distribution
- F5 narrative before/after
- F6 F1 vs NG threshold
- F7 precision-recall curve
- F8 F1 by LOH zone
- F9 ΔF1 V6 vs others
- F10 LR applied to baseline
- F11 LR ablation methylation
- F12 per-CpG HP × ALT 4-BAM
- F13 HPFineN cross-sample

---

## Sign-off Request

Per **4day_compressed_workflow.md Day 4-5 Hard Gates**:

✅ **Day 4 (2026-05-20, 1 day early)**: `git commit 8a90532` + `git tag -a v6-prod-20260520` 完成（local only on `fix/pon-only-phasing`）  
🔴 **Day 5 (next)**: 本 email send (user 親自 copy 到 mail client) + `git push origin fix/pon-only-phasing v6-prod-20260520`（push 需用戶明示授權）

**Anything to clarify or push back before send?**

— Generated 2026-05-20 by InterSubMod Day 3 5-Goal validation workflow.
