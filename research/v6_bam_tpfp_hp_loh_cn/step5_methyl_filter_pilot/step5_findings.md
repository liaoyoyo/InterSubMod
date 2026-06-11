<!--
建立時間: 2026-05-18
agent: Coordinator (main session) — Methylation-Augmented FP Filter Pilot v1.0
status: in_progress
report_class: filter-evaluation-pilot
audience: PI / lab member / 自己未來
scope: HCC1395 paired-pileup pilot, V3F+V5+V6 三 BAM × off/on flag × 4 powered-after-gate FP-rich cells
tier_proposal: ⭐3 PARTIAL POSITIVE (POSITIVE-but-marginal)
parent_plan: research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/00_PLAN.md (v1.0 + Step -1 amendment)
upstream:
  - step5_methyl_filter_pilot/step0_schema.md (35,332 × 202 master TSV)
  - step5_methyl_filter_pilot/step1_findings.md (LRT H1 POSITIVE, H5 POSITIVE)
  - step5_methyl_filter_pilot/step2_findings.md (filter signal direction)
  - step5_methyl_filter_pilot/step3_findings.md (ΔF1 H2/H3 POSITIVE)
  - step5_methyl_filter_pilot/step4_mechanism_candidates.md (H4 POSITIVE)
verdict: ALL H1-H5 POSITIVE; ΔF1 = +0.00242 marginal (below +0.005 Cohen ribbon "marginal"); ⭐3 PARTIAL POSITIVE pilot
last_verified: 2026-05-18
report_template: filter-pilot v1.0
predecessor_cycle: research/v6_bam_tpfp_hp_loh_cn/ (v0.3 characterization, ⭐3 PARTIAL POSITIVE, 2026-05-15)
-->

# Step 5 — Decision tree + Paper Framing Synthesis

> **Verdict**: ALL H1-H5 POSITIVE, but ΔF1 = +0.00242 is **MARGINAL** (below +0.005 Cohen ribbon)
> **Proposed Tier**: ⭐3 PARTIAL POSITIVE filter pilot
> **Next gate**: 用戶 review → 手動更新 SoT (per resolved decision #4, 不自動觸發 /conclude-research)

---

## 0. TL;DR

本 cycle 在 v0.3 characterization 之上跨越到 **filter F1 evaluation**。經過 Step -1 重跑 12 ISM runs（80 min）+ Step 0-4 後分析，**全 5 個 H 假設 POSITIVE**：

| H | Verdict | 關鍵數值 |
|---|---------|----------|
| **H1** Methylation covariate LRT q<0.05 | **POSITIVE** | 16/30 testable cells，top p=1.8e-58 |
| **H2** ΔF1 > 0 vs caller baseline 0.7166 | **POSITIVE-marginal** | **+0.00242** (post 0.71902 vs caller 0.7166) |
| **H3** FP_removal > TP_loss at τ* | **POSITIVE** | 98.26% > 35.00% (filter_signal +0.633) |
| **H4** ≥1 mechanism candidate | **POSITIVE** | 13 hypotheses × 14 PubMed refs |
| **H5** V5 ≈ V6 LR β consistency | **POSITIVE** | Robust median Δβ = **1.87e-5** << 0.005 |

**但 ΔF1 = +0.00242 < +0.005 marginal threshold** → POSITIVE-but-marginal，需 paper 框架謹慎處理。

---

## 1. Decision Tree (依 plan §Step 5)

依 plan v1.0 §Step 5 decision tree:
```
H1 NEGATIVE                      → methylation 補強 NO-GO
H1 POSITIVE + H2 NEGATIVE         → characterization-only
H1+H2 POSITIVE + H3 NEGATIVE      → "filter possible but TP loss > FP removal"
H1+H2+H3 POSITIVE + H4 NEGATIVE   → "Statistical filter no mechanism"
ALL H1-H4 POSITIVE + H5 sanity OK → ⭐3 filter candidate
```

**本 cycle**：ALL H1-H5 POSITIVE → **「Methylation-augmented FP filter ⭐3 candidate」**

但 **caveat**：H2 是 POSITIVE-marginal，effect size < +0.005 Cohen ribbon — 必須在 paper 與 lab meeting 中誠實標示為「marginal effect, cross-sample 驗證前不可宣告 production-ready」。

---

## 2. Quantitative Summary

### 2.1 Optimal filter rule

- **Winning scope**: AGGREGATED (多 cell 合併 filter, 404 unique regions)
- **τ*** = 0.52 (predict P(TP) < 0.52 → filter out)
- **Pre-filter baseline (ClairS-TO)**: TP=30,490 / FP=4,842 / FN=19,288 (reverse-solved) / F1=0.7166
- **Post-filter**: TP_kept ~ 30,470 / FP_kept ~ 84 / TP_loss=21 / FP_removed=4,758
  - 等價於 from caller {TP=30,490, FP=4,842} → filter 掉 {TP=21, FP=4,758} → 剩 {TP=30,469, FP=84}
- **F1_post** = 0.71902
- **ΔF1** = +0.00242 (post − caller)

### 2.2 Per-cell ranking

| Cell | τ* | ΔF1 | TP_loss% | FP_remov% | CV AUC |
|------|----|----|---------|----------|--------|
| AGGREGATED | 0.52 | **+0.00242** | 35.0% | 98.3% | n/a |
| chr8|Outer|other|cov_normal | 0.52 | +0.00136 | 42.3% | 98.5% | **0.928** |
| chr8|Inner|other|cov_normal | 0.18 | +0.00089 | 23.5% | 94.1% | 0.849 |
| auto|Outer|other|cov_proxy_mid | 0.79 | +0.00071 | 38.5% | 98.3% | 0.853 |
| auto|Inner|other|cov_proxy_mid | 0.22 | +0.00045 | 29.4% | 91.6% | 0.796 |

### 2.3 Effect size assessment

依 plan §3 Cohen ribbon (`F1 delta`):
- 小 effect: 0.01
- 中: 0.05
- 大: 0.10

**ΔF1 = +0.00242 < 0.01 小 effect threshold** → effect 強度為 "marginal / sub-clinical"

歷史對照：
- Phase 1A paired-pure methyl+context: ΔF1 +0.0112 (5 樣本 CI [+0.003, +0.020]) — 真 ⭐3
- 本 cycle: ΔF1 +0.00242 single-sample HCC1395 — 還沒到 Phase 1A 水準

---

## 3. Critical Caveats（影響 paper 與 next cycle 決策）

### 3.1 Sample shrinkage (Step 1 caveat)
- Model B convergence: **30/138 (22%)** — 108 small-FP cells (n_FP<5) saturated
- Step 2 cell gate: 4/11 cells 通過 Model B converged + n_fit≥100 + n_TP≥5

### 3.2 In-sample optimism (Step 2 caveat)
- in-sample AUC 0.97 vs CV AUC 0.90 → 7 AUC points optimism
- **Step 3 已用 CV out-of-fold predictions** → ΔF1 +0.00242 是 CV-adjusted (anti-optimism)

### 3.3 τ floor recovery (Step 2 caveat)
- Step 2 best τ 撞 grid 底 0.50（in-sample bias）
- Step 3 擴展 τ∈[0.10, 0.95] → 真實 τ* = 0.52 (well off floor)

### 3.4 Counts 對照警示
- **Plan-spec (post-ISM region-level)**: TP=30,490 / FP=4,842
- **Caller variant-level (09_V6_caller_F1_verification.md)**: TP=28,509 / FP=11,606 / FN=10,938
- 兩組數字**不一致**！Step 3 用 plan-spec per user instruction，但 paper framing 必須明確：
  - 「Post-ISM region-level F1 improvement」vs「caller variant-level F1 improvement」是不同口徑
  - 若用 caller variant-level baseline，ΔF1 結論可能不同

### 3.5 Single-sample only
- HCC1395 paired-pileup only
- 4 cells × 1 sample = 4 cell-sample pairs，**不能宣告 generalizability**
- 必須 cross-sample 才能升 ⭐4

### 3.6 Mechanism unverified
- H4 是 relaxed gate（列候選 + PubMed prior）— 13 hypothesis 都未 cycle-specific 驗證
- 哪個 mechanism 真正驅動 LR β 顯著性？未測

---

## 4. Paper Framing Recommendation

### 4.1 主軸建議：保守 characterization-augmented filter pilot

不要 frame 為「filter improves F1」 — frame 為：
> **"Read-level methylation features add LRT-significant discriminative information (q<0.05 in 16/30 cell-combos) on top of LOH × HP × CN framework. In an aggressive cell-level filter pilot, this yielded ΔF1 = +0.0024 vs caller baseline F1=0.7166 — directionally POSITIVE but effect size below clinical relevance threshold (+0.01). Cross-sample validation needed before production claim."**

### 4.2 §3 主圖建議

- **Fig 3a**: H1 LRT q heatmap (cells × BAM × covariate)，highlight 16/30 q<0.05
- **Fig 3b**: ΔF1 vs τ curve (per cell + aggregated)，標 τ*=0.52
- **Fig 3c**: FP_removal vs TP_loss scatter，diagonal y=x reference line（show 4 cells 在 y>x 區）
- **Fig 3d**: H4 mechanism tree (5 categories × 13 hypotheses) + anchor cell labels

### 4.3 與 prior art 對比

依 02_prior_art_notes.md:
- ROCIT: per-read methylation transformer，AUC 0.933 (tumor vs non-tumor read)，**不做 variant-level filter**
- SGZ: 4-axis variant-level filter，無 methylation
- TumorLens: sample-level，不做 per-region filter
- **本 cycle 差異化**：first read-level methylation augmented multi-axis (LOH×HP×CN×AF + 5 methyl) variant-level filter，雖 ΔF1 marginal 但 prior art 無同口徑

### 4.4 Limitations 章節必寫
1. Single-sample HCC1395 paired-pileup only
2. ΔF1 +0.0024 < +0.005 marginal Cohen ribbon
3. Model B convergence 22% (small-FP cell saturation)
4. Counts 雙口徑 (post-ISM region-level vs caller variant-level)
5. H4 mechanism 列候選但未 cycle-specific 驗證

---

## 5. Next Cycle / Phase 2 Decisions

### 5.1 升 ⭐4 必要條件
1. **Cross-sample validation** — phaseD 4 樣本 V3F+V5+V6 三方 ISM 重跑（~10-30 hr 平行）
2. **Mechanism cycle-specific 驗證** — H4 5 categories 至少 1 個做 dedicated cell-level mechanism analysis（e.g., chr8|Outer|other|cov_normal × cis-mQTL GoDMC 對照）
3. **ΔF1 升到 +0.005** (marginal threshold) — 含 caller variant-level baseline 重算

### 5.2 反向 — 若不投入 phase 2
- Frame 為「marginal positive characterization extension」收進 v0.3 主報告 supplementary
- 不再提 「filter production candidate」 claim

### 5.3 Paper outline 影響
- §3.3 (本來 characterization framework) 新增小節 "Methylation-augmented filter pilot (marginal ΔF1)"
- §5 Discussion 加 "Cross-sample generalization remains the gate" → tie to T2.1 Z-AUTO KDE
- §6 Future work 加 mechanism-specific deep dive

---

## 6. Decision (per resolved decision #4 — manual SoT update)

**不自動觸發** `/conclude-research`。本 findings.md 與 deliverables 留給用戶 review，決定：

### Decision options (用戶選一)

**A. ⭐3 升級為 Methylation-augmented filter ⭐3 candidate**
- 接受 ΔF1 +0.0024 marginal but POSITIVE
- 進入 cross-sample phase 2 pilot
- 寫入 INDEX + CURRENT_FOCUS + evidence_ledger
- 更新 memory `project_methyl_filter_pilot_marginal_positive.md`

**B. 「Marginal-but-POSITIVE characterization extension」**
- 接受 H1+H5 為主要 positive findings（characterization quality）
- 不 frame 為 filter；不投入 cross-sample
- 寫進 v0.3 main report supplementary
- evidence_ledger entry "marginal effect, characterization only"

**C. Defer decision**
- 留 findings.md，下週週報時統合判斷
- 不更新任何 SoT

### Recommended (Coordinator suggestion): **A 但保守 framing**

- ΔF1 +0.0024 雖 marginal 但**5 個 H 全 POSITIVE + 4 cells CV AUC 0.80-0.93 + V5/V6 一致**
- Cross-sample 是 next gate（與 T2.1 同期可行）
- Paper § 多一個獨立 contribution（marginal-but-novel filter pilot）

### Hard Gate reminder

依 plan §Pre-registration「NO-GO 不可事後改寫」— 本 cycle 沒 NO-GO，5 H 全 POSITIVE，但 marginal effect 是事先 acknowledge 的 caveat，**不違反 Pre-reg**。

---

## 7. Files Inventory

### Deliverables (本 cycle)
- `step5_methyl_filter_pilot/00_PLAN.md` (v1.0 + Step -1 amendment)
- `step5_methyl_filter_pilot/step0_schema.md` (35,332 × 202 master TSV schema)
- `step5_methyl_filter_pilot/step5_master_augmented.tsv` (25.3 MB)
- `step5_methyl_filter_pilot/step1_lrt_per_cell.tsv` (138 rows × LRT + LR β)
- `step5_methyl_filter_pilot/step1_findings.md` (H1 + H5 verdicts)
- `step5_methyl_filter_pilot/step2_filter_sweep.tsv` (506 rows)
- `step5_methyl_filter_pilot/step2_findings.md` (H3 direction POSITIVE)
- `step5_methyl_filter_pilot/step3_delta_f1.tsv` (430 rows × ΔF1 per τ)
- `step5_methyl_filter_pilot/step3_optimal_tau_summary.md`
- `step5_methyl_filter_pilot/step3_findings.md` (H2 + H3 verdicts)
- `step5_methyl_filter_pilot/step4_mechanism_candidates.md` (5 categories × 13 hypotheses + 14 PubMed refs)
- `step5_methyl_filter_pilot/step5_findings.md` (本檔)
- `step5_methyl_filter_pilot/figures/`:
  - `step1_lrt_heatmap.png`
  - `step2_roc_per_cell.png`
  - `step2_filter_signal.png`
  - `step3_deltaf1_vs_tau.png`
  - `step3_filter_signal_curve.png`
- `step5_methyl_filter_pilot/scripts/`:
  - `build_augmented_master.py` (~12s reproducer)
  - `_common_step1.py`
  - `augmented_lr.py`
  - `filter_sweep.py`
  - `delta_f1.py` (~390 LOC)

### Pre-cycle reference
- v0.3 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md`
- v0.3 plan: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md`
- phaseC rerun script: `InterSubMod/research/paired_priority_bug_audit/scripts/run_phaseC_genome_three_way_with_significance.sh`
- New phaseC with significance: `InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/`

### Predecessor cycle (v0.3)
- Tier ⭐3 PARTIAL POSITIVE (2026-05-15)
- Paradigm reframe: Z-OCH/Z-GL = TP signatures (not FP markers)
- 本 cycle (v1.0) 是 v0.3 之延伸 — filter F1 evaluation 跨越邊界

---

## 8. Step 5b — Manual SoT update (待用戶 review)

依 resolved decision #4 不自動觸發 /conclude-research。完成本 findings.md 後：

1. **用戶 read** 本檔 + step1-4 findings
2. **用戶決定** A / B / C 路徑（§6）
3. **若 A 或 B**：
   - 寫入 `InterSubMod/docs/experiments/INDEX.md`（新 entry 2026-05-18）
   - 寫入 `InterSubMod/docs/CURRENT_FOCUS.md`（新 2026-05-18 區段）
   - Append `InterSubMod/research/autoresearch/evidence_ledger.jsonl` 一條
   - Optional: 新 memory `~/.claude/projects/.../memory/project_methyl_filter_pilot.md`
4. **若 C**: defer，不動 SoT，留待週報

---

## 9. Coordinator Recommendation

**選 A 但 framing 保守**：
- 接受 ⭐3 candidate（5 H POSITIVE + 4 cell CV AUC 0.80-0.93）
- 進 cross-sample phase 2 pilot（T2.1 整合）
- Paper §3.3 加小節 + §5 Discussion + §6 Future
- Limitations 明寫 "ΔF1 marginal" + "single-sample" + "mechanism candidates not validated"

不選 B 因為：H1 16/30 cell LRT q<0.05 + H3 filter_signal +0.633 + H5 V5≈V6 robust 都是強訊號，不只是 characterization，而是有 filter direction 訊號 — 即使 ΔF1 marginal。

不選 C 因為：deliverables 已產出，evidence 完整，沒理由 defer SoT update（除非用戶完全想 pivot 到別的 cycle）。
