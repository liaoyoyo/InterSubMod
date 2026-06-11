# A0 — 既有素材完整性查驗總彙總

> **執行日期**: 2026-05-21
> **執行方式**: 6 個並行 audit agents (fan-out → main agent funnel merge)
> **總 audit 範圍**: 11 個 artifact (4 個 .md + 7 個 standalone.html)
> **總 numerical claim 數**: ~273 (跨 6 agents)
> **總 verification rate**: 254/273 verified (~93%)
> **Fabrication detected**: **0** (無 LLM 內插假數字)

---

## 1. 高階結論

| Verdict | 數量 | 含義 |
|---|---:|---|
| **PASS (PI-distributable as-is)** | 6/11 artifacts | 數字 100% 對 source TSV/CSV/原始 .md ✓ |
| **PASS with minor caveats** | 4/11 artifacts | 核心數字對齊，1-4 個 minor inconsistency 待 footnote 處理 |
| **PASS with row-mislabel** | 1/11 artifacts | 20260514 priority bug engineering .md §8.2 4-sample ratio row 錯位 |
| **FAIL / retraction** | 0/11 artifacts | 無 |

**核心 finding**：所有 PI-critical 數字 (F1=0.7166/0.6273, marker +52.4%, HP ratio 1.696→1.609, hp=3 ×13.2, 17.3:1, LOSO -0.00012, Cycle 1 +0.02236) **全部 verified against source TSV / JSON / commit hash**。**無 2026-05-20 +0.057 等級 fabrication**。

---

## 2. Verified 可直接引用 artifact 清單

| # | Artifact 路徑 | Verdict | Claims | Verified | 引用建議 |
|---|---|---|---:|---:|---|
| 1 | `InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md` | PASS minor | 43 | 42 | Doc 1 §7 整合 matrix 可直接引用；§7 整體 PI-safe |
| 2 | `InterSubMod/docs/reports/validated/2026/05/20260513_V6_Attribution_Errata_01.md` | PASS minor | 30 | 29 | §3a.5b.7-9 final attribution 可直接引用；C2.25 (47,838 vs 47,798) 為 SNV-only filter 差異不影響結論 |
| 3 | `InterSubMod/docs/reports/validated/2026/05/20260521_PI_V6_signoff_email_draft_5goal.md` | **PASS 100%** | 22 | 22 | 全 22 個 claim verified 可直接給 PI |
| 4 | `InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_priority_bug_engineering_report_01.standalone.html` | PASS minor | 22 | 22 | SUPERSEDED-banner 顯示已被 4-Layer Synthesis 取代；數字仍全對 |
| 5 | `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.standalone.html` | **PASS 100%** | ~14 | 14 | G1/G2/G3 主力素材；可直接重用 |
| 6 | `InterSubMod/docs/experiments/in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.standalone.html` | **PASS 100%** | ~110 | 110 | 目標 2 主力重用素材；可直接重用 |
| 7 | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/phase2_engineering_completeness.standalone.html` | **PASS 100%** | — | — | 目標 4 方法論主力 |
| 8 | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_pi_verification/phase2_pi_trust_framework.standalone.html` | **PASS 100%** | — | — | LOSO reframe 主力；可直接重用 |

---

## 3. 待謹慎引用 artifact (有 row-mislabel 或 framing 不準)

### 3.1 `InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_priority_bug_engineering_report_01.md` §8.2

**問題**: §8.2 4-sample cross-sample ratio table **row-mislabel**：

| Sample | Doc §8.2 寫的 ratio ❌ | Source TSV 實際 ratio ✓ |
|---|---|---|
| H1437 | 0.611 ❌ | 1.243 ✓ |
| H2009 | 1.243 ❌ | 0.901 ✓ |
| HCC1954 | 1.014 ❌ | 0.958 ✓ |
| HCC1937 | 0.978 ❌ | 0.611 ✓ |

- Aggregate range "0.611-1.243" 正確；但 per-row pair 錯
- Marker rate column **未錯位**，只 ratio column 錯
- 文檔 narrative outside §8.2 不受影響
- **Source TSV**: `InterSubMod/research/paired_priority_bug_audit/phaseD_v6_5sample/v6_cross_sample_summary.tsv`

**處理決策 (D012)**: 不直接修改 validated 報告；PI 4-goal HTML target1 從 source TSV **重新引用**並加 footnote 註記 §8.2 mislabel

### 3.2 `InterSubMod/docs/reports/validated/2026/05/20260514_HP_tag_17to1_rootcause_explained_01.standalone.html` §4 TL;DR

**問題**: §4 "VCF GT2 約 95% 偏 `1|.`" 把上游 L1 ratio (1.77:1 = 63.9%) 與下游 BAM HP1 占比 (94.6%) 混淆

- 屬 pre-erratum framing (已被 4-Layer Synthesis §1 line 488 釐清為 26,436/14,931 = 1.77:1)
- 所有 headline 數字 (17.3:1, 94.6%, 752, 34,855) 正確
- 只有「~95% GT2 偏 1|.」這句技術上不準確

**處理決策 (D018)**: 不直接修改 validated HTML；PI 4-goal HTML 採 4-Layer Synthesis §1 line 488 (1.77:1 = 26,436/14,931 = 63.9% L1 ratio) 作為「VCF GT2 分佈」描述源；避免「~95% GT2 偏 1|.」這句

### 3.3 `InterSubMod/docs/reports/validated/2026/05/20260514_Self_Phasing_4Layer_Synthesis_Engineering_Report_01.html`

**3 個 suspect**:
1. "17,404 V3F victim subset" (§5.2) — undocumented in audit dirs
2. "Phase N50 ~20,000+" (§4.1) vs verified `8,109` — **internal inconsistency**
3. "hp=33" dual-population usage: `571` (HAP3 altHaplotype assignment in §1.3) vs `138,317` (BAM HP:i:33 reads in §3.4) — same label different populations

**2 個 unverifiable**:
- "+36% pipeline speed" (§3.1) — 無 benchmark TSV
- "V5 BAM = V6 BAM 0 diff read-by-read" (§5.2) — claim 沒 backup TSV

**處理決策**: PI 4-goal HTML target1 引用此 HTML 時 — 避開 hp=33 dual-usage 寫法（改 "altHaplotype assignment HAP3 = 571" 與 "BAM HP:i:33 reads = 138,317" 兩個分開敘述）；不引用 +36% speed claim；N50 統一用 verified 8,109

---

## 4. PI 4-goal HTML 引用 ground rules（基於 audit）

| 數字 | 安全引用源 (verified) | 避免引用源 |
|---|---|---|
| F1 = 0.7166 / 0.6273 (caller invariant) | `09_V6_caller_F1_verification.md` 或 phase2_pi_trust_framework HTML | — |
| Marker +52.4% / +1.26pp / ×13.2 | 20260519 standalone HTML 或 `step1_baseline_vs_v6_findings.json` | — |
| HP ratio 17.3:1 baseline / 1.84:1 V6 / 1.85:1 V3F-only | 20260514_HP_tag_priority_bug_engineering_report HTML (verified) 或 4-Layer §S5 | 17to1 standalone HTML §4 TL;DR 那句 |
| Phase D 4 樣本 ratio (0.611-1.243) | **`phaseD_v6_5sample/v6_cross_sample_summary.tsv` row-by-row** | 20260514 priority bug .md §8.2 (row-mislabel) |
| F1 LOSO -0.00012 / mean -0.00004 / Wilcoxon p=0.125 | `cycle4/loso_validation/data/loso_cv_results.tsv` | — |
| Cycle 1 +0.02236 / 9.24× / τ=0.39 | `cycle1/cycle1_findings.md` 或 `cycle1_track_a_filter.json` | — |
| 1.77:1 L1 ratio (= 26,436/14,931) | 4-Layer Synthesis §1 line 488 | 17to1 standalone HTML "~95% GT2" framing |
| master TSV col count | **112 cols** (Plan 與 D015 校正) | 116 cols (誤) |
| N50 (V6 phased) | 8,109 (verified) | "~20,000+" (4-Layer §4.1 hand-wave) |

---

## 5. 6 個 Agent audit report deep-dive paths

| Agent | 範圍 | Audit report path |
|---|---|---|
| 1 | V6 binary main docs (20260511 + 20260513) | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_audits/agent1_v6_binary_main_docs_audit.md` |
| 2 | priority_bug eng .md + PI signoff email | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_audits/agent2_priority_bug_PI_signoff_audit.md` |
| 3 | 17to1 + priority_bug eng HTML 2 個 | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_audits/agent3_priority_bug_HTML_audit.md` |
| 4 | 4-Layer Synthesis HTML + V6 vs baseline HTML summary | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_audits/agent4_4layer_summary_HTML_audit.md` |
| 5 | 20260515 V6 TPFP characterization + 20260519 V6 vs baseline HCC1395 HTML | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_audits/agent5_ISM_TPFP_HTML_audit.md` |
| 6 | phase2_engineering_completeness + phase2_pi_trust_framework HTML | `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_audits/agent6_phase2_verification_HTML_audit.md` |

---

## 6. Verdict 統計

```
Total artifacts: 11
- 4 .md (V6 binary doc / errata / priority bug eng / PI signoff email)
- 7 .standalone.html (17to1 / priority bug eng / 4-layer / V6 summary / V6 TPFP / V6 vs baseline / phase2 engineering + pi trust)

Total numerical claims audited: ~273
- verified: ~254 (93%)
- suspect: ~17 (6%) - all minor framing / row-mislabel / internal inconsistency
- unverifiable: ~2 (1%) - +36% pipeline speed / V5≡V6 0-diff (claim 沒 backup TSV)
- fabrication: 0 ✓

Verdict by artifact:
- PASS 100%: 5 artifacts (Doc 3, Doc 4, Doc 5, Doc 6, Doc 7, Doc 8 in §2)
- PASS minor: 4 artifacts (Doc 1, Doc 2 in §2; Doc 4 priority bug eng .md, Doc 17to1 HTML — 已列 §3.1, §3.2)
- PASS with caveats: 1 artifact (4-Layer Synthesis — 已列 §3.3)
- FAIL: 0
```

---

## 7. 下一步

✅ Phase A0 audit 完成
→ 啟動 **Phase A1-A4 並行**（不被 audit 結果阻擋的數據補測）

Phase A1: HCC1395 baseline vs V6 HP 6 值對比 + V6 5 sample HP 分佈
Phase A2: HCC1395 baseline vs V6 ALT-only ratio + V6 5 sample ALT-only ratio  
Phase A3: 5 features TP/FP boxplot + 9-cell heatmap (cycle2 master TSV 已有 5 sample data)
Phase A4: LR/DT/RF/XGBoost LOSO benchmark (cycle4 framework 重用)
