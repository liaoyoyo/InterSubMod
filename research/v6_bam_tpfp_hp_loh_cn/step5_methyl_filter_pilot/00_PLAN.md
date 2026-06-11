# Methylation-Augmented FP Filter Pilot (HCC1395)

> Plan v1.0 — 2026-05-18
> **Predecessor**: v0.3 characterization cycle ⭐3 PARTIAL POSITIVE (completed 2026-05-15)
>   - Plan: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md`
>   - 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md`
> **Different task**: 跨越 characterization → filter F1 evaluation
> **Scope**: HCC1395 single-sample pilot (cross-sample 留 phase 2)
> **Mechanism gate**: relaxed — 列候選機制供後續驗證即可，不強制 literature prior

---

## Context — 為什麼需要這個新 cycle

### v0.3 留下的問題

v0.3 是嚴格 characterization-only（plan §Out-of-scope 第 1 條），找到 3 個關鍵 finding：
1. **H4 chr8 hotspot LR deviance: CN 0.211 > HP 0.063**（(LOH+CN)-HP = +0.186 POSITIVE）
2. **Paradigm reframe**: Z-OCH/Z-GL 是 TP signatures 而非 FP markers
3. **Framework gap**: 4 zone framework 只 cover ~37% FP，63% FP unexplained

用戶問題：**「63% unexplained FP 用 methylation augmentation 能否轉成 actionable filter，使 FP removal % > TP loss % 且 ΔF1 > 0 (相對於 ClairS-TO caller-only F1=0.7166)？」**

### 與歷史 filter NO-GO 的差異化（防重蹈覆轍）

| 歷史 NO-GO | 本 cycle 差異 |
|------------|-------------|
| LOH binary filter 10/10 失敗 (2026-04-06) | 不用 binary；用 cell-level multi-axis LR-predicted threshold |
| CN zone-aware filter 跨樣本崩 (2026-04-10) | HCC1395 pilot 先；cross-sample generalization 留 phase 2 |
| Pure methylation clustering AUC 0.512 (Option C, 2026-04-07) | 不用 pure methylation；當作 5th-9th covariate 補強現有 4 軸 |
| O9-O13 甲基化系列 NEGATIVE | 加 within-cell within-AF OLS confound guard 防 collider bias |

### 與 ClairS-TO caller baseline 的比較
- **F1 baseline** = 0.7166 (ClairS-TO HCC1395 vs SEQC2 truth)
- 來源: `InterSubMod/research/paired_priority_bug_audit/06_V3F_vs_V5_evaluation.md` §1 + `09_V6_caller_F1_verification.md`
- **Filter target**: ΔF1 > 0 表示「post-filter F1 > caller-only F1」
- **本 cycle 是 pilot** — POSITIVE 也只是 first signal，需 cross-sample 才能升 ⭐4

---

## Pre-Registration (Hard, 不可事後改寫)

| ID | Prediction | Falsification | Decision threshold |
|----|-----------|---------------|-------------------|
| **H1** | ≥1 methylation covariate 在 powered cells 內 LRT q<0.05 (vs baseline LR) | 0 cells reach q<0.05 | BH-FDR q<0.05 |
| **H2** | FP-rich cells LR-threshold filter 使 HCC1395 ΔF1 > 0 vs caller F1=0.7166 | ΔF1 ≤ 0 | post-filter F1 - 0.7166 > 0 |
| **H3** | FP removal % > TP loss % (post-filter) | FP removal ≤ TP loss | 用 region-level confusion matrix 計算 |
| **H4** | 至少提出 1 個候選生物學機制（cis-mQTL / cancer-specific ASM / allele-imbalance / repeat-context）對應 powered cell | 0 候選機制 | relaxed — 不要求 literature prior，列候選供 follow-up |
| **H5** | V5 vs V6 BAM 對應 ΔF1 差異 ≤ 0.003（V6 重用 V5 phased VCF 預期 methylation 一致） | ΔF1 V5 vs V6 > 0.005 | sanity check |

**NO-GO 條件**：H1 OR H2 OR H3 任一 NEGATIVE → 整 cycle NEGATIVE，寫入 evidence_ledger 不可事後改寫

---

## Plan Structure — 5 Steps

### Step 0 — Build augmented master TSV (BLOCKING prerequisite)

**問題**：v0.3 的 `step1_master_three_way.tsv` 沒有 methylation features（Agent A 沒 join）。

**任務**：從 `research/paired_priority_bug_audit/phaseC_genome_three_way/{V5,V6}_{off,on}_{tp,fp}/filtered_snv_{tp,fp}/{chr}/{region_id}/significance_summary.csv` 撈 12 個 methylation 欄位 join 進 master TSV。

**12 features 完整 list**：
1. HPMergedDelta（HP 間甲基化差 mean delta）
2. HPFineF（HP fine bucket PERMANOVA F-stat）
3. HPMergedSig（HP merge significance）
4. NME_HP1, NME_HP2 → 衍生 `NME_imbalance = |NME_HP1 - NME_HP2|`
5. Epipoly_HP1, Epipoly_HP2, Epipoly_Delta
6. ClusterPermanovaF, ClusterPermanovaP
7. AlleleDelta, AlleleP
8. Entropy_Imbalance

**Schema 設計**：
- 既有 step1 master 64 cols + 新增 V5/V6 各 12 features = 88 cols
- V5 vs V6 methylation 預期完全一致（重用 phased VCF）— H5 sanity check 用

**輸出**：
- `step5_methyl_filter_pilot/step5_master_augmented.tsv`
- `scripts/build_augmented_master.py`

**估時**：~30 min (region-level 數千個 CSV file I/O；可平行處理 by chr)

### Step 1 — Augmented LR + LRT comparison (V3F + V5 + V6 三 BAM 對照)

**3 BAM 對照設計**（依 resolved decision #3）：對每 BAM 各自跑 Model A / Model B，看 V3F → V5 → V6 trajectory of LRT 顯著性與 LR β 變化。H5 sanity check 比 V5 vs V6 LR β 一致性（預期 ≤ 0.003 差異）。

**模型對照**：
```
Model A (baseline, v0.3): logit(P(TP)) = β0 + β1·NG + β2·caller_af + β3·NumReads
Model B (augmented):     + β4·HPMergedDelta + β5·HPFineF + β6·NME_imbalance + β7·Epipoly_Delta + β8·ClusterPermanovaF
```

對每 BAM × 每 powered cell 各跑一次 → 共 (3 × 23 powered cells) × 2 models = 138 LR fits

**對每 powered cell (n≥50)** 跑：
1. Fit Model A、Model B（statsmodels logit）
2. LRT: H_0 (Model A) vs H_1 (Model B), df=5
3. BH-FDR q-value (全 cells 校正)
4. 報每個 covariate 的 β + LRT p（drop-one-axis ablation）

**Confound guard**（依 v0.3 ⁒2 §4 + memory `feedback_L2_collider_bias.md`）：
- Within-cell within-AF-bin OLS（避免 L2 collider）
- 5-fold cross-validation hold-out

**輸出**：`step5_methyl_filter_pilot/step1_lrt_per_cell.tsv` (cells × LRT_p / q / 5 methyl betas / dev_explained delta)
**H1 verdict**：是否 ≥1 cell q<0.05

### Step 2 — FP-rich cells filter threshold sweep

**範圍**：只在 FP_rate > global × 1.5 的 cells 內做（v0.3 已標記出 12 cells；Z-CHR8 sub-cells + Z-AUTO sub-cells）

**流程**：
1. 對每 FP-rich cell 用 Model B (有顯著 LRT) 或 Model A (LRT NEGATIVE) predict P(TP) per region
2. Threshold τ sweep: 0.5 → 0.95 by 0.01 (45 個 τ values)
3. 對每 τ 算:
   - TP_kept / TP_total (recall after filter)
   - FP_removed / FP_total (precision improvement)
   - TP_loss_pct = 1 - TP_kept/TP_total
   - FP_removal_pct = FP_removed/FP_total
4. 找 **「FP_removal_pct > TP_loss_pct」** 的 τ 範圍（H3 評估）

**輸出**：
- `step5_methyl_filter_pilot/step2_filter_sweep.tsv` (per cell × τ × TP/FP metrics)
- `step5_methyl_filter_pilot/figures/step2_roc_per_cell.png` (precision-recall curve per cell)

### Step 3 — ΔF1 vs caller-only baseline

**全域 ΔF1 計算**：
```
Caller-only baseline:
  TP_caller = 30,490 (HCC1395 ClairS-TO post-PON)
  FP_caller = 4,842
  TN, FN 從 SEQC2 truth 算 (需先重建)
  P_caller = TP_caller / (TP_caller + FP_caller)
  R_caller = TP_caller / (TP_caller + FN_caller)
  F1_caller = 2·P·R / (P+R) = 0.7166

Post-filter (for selected τ*):
  TP_post = TP_caller - TP_lost
  FP_post = FP_caller - FP_removed
  FN_post = FN_caller + TP_lost (filtered TP 變成 FN)
  P_post, R_post, F1_post
  ΔF1 = F1_post - F1_caller
```

**選 τ* 策略**：max ΔF1 across all FP-rich cells 合併 filter rule

**輸出**：
- `step5_methyl_filter_pilot/step3_delta_f1.tsv` (per τ × ΔF1)
- `step5_methyl_filter_pilot/step3_optimal_tau_summary.md` (max ΔF1 + 對應 cells + 對應 τ)

**H2 verdict**：max ΔF1 > 0 → POSITIVE pilot
**H3 verdict**：在 optimal τ* 是否 FP_removal_pct > TP_loss_pct

### Step 4 — Mechanism hypothesis brainstorm (relaxed gate)

對 H1 LRT 顯著的 cells，每 cell 列 1-3 候選機制（**不要求 prior literature**）：

**候選機制 categories**：
- **cis-mQTL**：germline variant 影響附近 CpG methylation → 解釋為什麼 germline FP 在某 methylation pattern 富集
- **Cancer-specific ASM**：腫瘤特異甲基化失調（如 promoter hypermethylation）→ 解釋 TP 的 methylation signature
- **Allele-imbalance**：當 LOH + CN gain 時，allele-specific methylation 可能被放大
- **Repeat / SD context**：重複序列附近 methylation 變異性大 → confound 訊號
- **Replication timing**：early replicating region methylation 較穩定

**每候選機制需提供**：
- 一句話 mechanism statement
- **PubMed literature search** (依 resolved decision #2)：每候選機制扎 1-2 篇 reference（cis-mQTL / ASM / cancer methylation atlas）
- 可驗證 follow-up（例：「比對 GTEx mQTL 資料庫」、「比對 TCGA methylation atlas」）
- 對應 region 範例（chr + pos + 該 cell 內 top TP/FP examples）

**Literature search 規範**（用 `/citation-verification` skill 確認 paper 存在）：
- PubMed query 範本：`(cis-mQTL OR allele-specific methylation) AND (BRCA1 OR TP53 OR chr8 amplicon OR triple-negative breast cancer)`
- 每 candidate 至少 1 paper DOI / PMID
- 引用品質：偏好 ≥2024 review papers + 領域 canonical papers (Tycko, Greally, Pidsley, etc.)

**輸出**：`step5_methyl_filter_pilot/step4_mechanism_candidates.md`

**H4 verdict**：至少 1 個候選即 POSITIVE（不評估 mechanism 正確性，只評估「我們有沒有想出可測試的假設」）

### Step 5 — Decision tree + Paper framing

```
H1 NEGATIVE (LRT 全不顯著)        → methylation 補強 NO-GO confirmed (v0.3 已知 ceiling)
H1 POSITIVE + H2 NEGATIVE (ΔF1≤0) → "methylation characterization-only", 補強 v0.3 結論
H1+H2 POSITIVE + H3 NEGATIVE       → "filter possible but TP loss > FP removal" — 暫不投入 production
H1+H2+H3 POSITIVE + H4 NEGATIVE    → "Statistical filter no mechanism" — 高風險，需 phase 2 機制驗
ALL H1-H4 POSITIVE + H5 sanity OK  → "Methylation-augmented FP filter ⭐3 candidate" — 進 cross-sample pilot
```

**輸出**：`step5_methyl_filter_pilot/step5_findings.md`（含 verdict + decision + paper framing 建議）

### Step 5b — Manual SoT update (依 resolved decision #4)

**不自動觸發** `/conclude-research`。Step 5 完成後產出 deliverables，用戶 review：
1. 讀 `step5_methyl_filter_pilot/step5_findings.md`
2. 確認 H1-H5 verdicts 站得住
3. **用戶決定**：是否寫入 SoT
4. 若 OK → 手動執行：
   - 更新 `InterSubMod/docs/experiments/INDEX.md` 加新 entry
   - 更新 `InterSubMod/docs/CURRENT_FOCUS.md` 加 2026-05-XX entry
   - append `InterSubMod/research/autoresearch/evidence_ledger.jsonl` 一條
   - 若新增 memory file，寫 `~/.claude/projects/.../memory/project_methyl_filter_pilot.md`

---

## 關鍵檔案路徑

### 輸入（read-only）
| 用途 | 路徑 |
|------|------|
| v0.3 step1 master TSV (base) | `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_three_way.tsv` |
| Region-level methylation source | `InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way/{V5,V6}_{off,on}_{tp,fp}/filtered_snv_{tp,fp}/{chr}/{region_id}/significance_summary.csv` |
| Caller F1 baseline 數值 | `InterSubMod/research/paired_priority_bug_audit/06_V3F_vs_V5_evaluation.md` §1 |
| SEQC2 truth set (for FN 計算) | `InterSubMod/research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv` |
| v0.3 powered cells list | `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step2_three_axis_grid/step2_grid_3d.tsv` |
| v0.3 FP-rich cells | `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step3_fp_zone_zoom/step3_cross_zone_summary.tsv` |

### 輸出（新 cycle 寫到）
新目錄：`InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/`
- `00_PLAN.md` (本計畫副本)
- `step5_master_augmented.tsv`
- `step1_lrt_per_cell.tsv`
- `step2_filter_sweep.tsv`
- `step3_delta_f1.tsv`
- `step3_optimal_tau_summary.md`
- `step4_mechanism_candidates.md`
- `step5_findings.md`
- `scripts/{build_augmented_master,augmented_lr,filter_sweep,delta_f1,mechanism_brainstorm}.py`
- `figures/step2_roc_per_cell.png`, `figures/step3_deltaf1_vs_tau.png`

---

## Reuse 既有 functions / data (避免重寫)

| 既有資源 | 路徑 | 怎麼 reuse |
|---------|------|-----------|
| HP bucket 計算邏輯 | `step1_v3f_v5_v6_three_way/scripts/build_three_way_master.py` | 直接 import / 改寫 join logic |
| Power gate + Wilson CI | `step2_three_axis_grid/scripts/build_grid.py` | reuse cell-level n/CI 計算 |
| Confound guard 5+2 道 | `step2_three_axis_grid/scripts/confound_guard.py` | 加入 Model B 版本 LRT |
| LR per-cell pattern | `step3_fp_zone_zoom/scripts/analyze_zones.py` (chr8 4-axis LR) | 加 5 methylation features |
| 12 methylation features 來源 | `phaseC_genome_three_way/V6_off_tp/filtered_snv_tp/chr*/region/significance_summary.csv` | direct read csv header → join |

---

## 驗證方式 (end-to-end)

### 階段 verification
1. **Step 0 通過條件**: augmented master TSV 含全 12 methylation features × V5/V6 = 24 個新 cols；non-null rate > 95%
2. **Step 1 通過條件**: 至少 1 cell q<0.05 (H1) OR 全 NEGATIVE 也算 valid result (寫進 evidence_ledger)
3. **Step 2 通過條件**: per cell ROC curve 產出；FP-rich cells (12 個) 全 sweep 完成
4. **Step 3 通過條件**: ΔF1 vs τ curve；optimal τ* 標出
5. **Step 4 通過條件**: 至少 1 候選機制提案（不評正確性）
6. **Step 5 通過條件**: decision tree verdict 寫入 evidence_ledger + memory

### 命令模板
```bash
# Step 0
python step5_methyl_filter_pilot/scripts/build_augmented_master.py \
    --base step1_v3f_v5_v6_three_way/step1_master_three_way.tsv \
    --methyl-source research/paired_priority_bug_audit/phaseC_genome_three_way \
    --bams V5,V6 \
    --output step5_methyl_filter_pilot/step5_master_augmented.tsv

# Step 1
python step5_methyl_filter_pilot/scripts/augmented_lr.py \
    --master step5_methyl_filter_pilot/step5_master_augmented.tsv \
    --powered-cells step2_three_axis_grid/step2_grid_3d.tsv \
    --output step5_methyl_filter_pilot/step1_lrt_per_cell.tsv

# Step 2
python step5_methyl_filter_pilot/scripts/filter_sweep.py \
    --lr-models step1_lrt_per_cell.tsv \
    --fp-rich-cells step3_fp_zone_zoom/step3_cross_zone_summary.tsv \
    --tau 0.5,0.95,0.01 \
    --output step5_methyl_filter_pilot/step2_filter_sweep.tsv

# Step 3
python step5_methyl_filter_pilot/scripts/delta_f1.py \
    --sweep step2_filter_sweep.tsv \
    --baseline-f1 0.7166 \
    --output step5_methyl_filter_pilot/step3_delta_f1.tsv
```

### Final deliverable verification
- 寫入 `InterSubMod/docs/experiments/in_progress/2026/05/20260518_V6_Methyl_Filter_Pilot_01.md` (doc-standards frontmatter)
- 同步 `InterSubMod/docs/experiments/INDEX.md`
- 寫入 evidence_ledger.jsonl 一條 entry（不可事後改）

---

## Out-of-scope (避免 scope creep)

- ❌ **Cross-sample generalization**（user 明示「先 HCC1395 為先」；H1437/H2009/HCC1954/HCC1937 留 phase 2）
- ❌ **Mechanism literature validation**（H4 只要求列候選，不要求事先有 PubMed paper）
- ❌ **Repeat / mappability / SV proximity 新軸**（留 phase 3，依本 cycle 結果）
- ❌ **Read-level methylation embedding (類 ROCIT transformer)**（大工程，留 future work）
- ❌ **C++ pipeline 修改**（純 post-hoc TSV 分析）
- ❌ **重跑 BAM / ISM**（methylation features 已存在於 phaseC region-level CSV）
- ❌ **V3F filter 評估**（用戶要 V5 + V6，V3F 已是 baseline pre-Layer-1.5 對照）

---

## 預期時程（依 4 條 resolved decisions 更新）

| Step | 估時 | 依賴 | 變動 |
|------|------|------|------|
| Step 0 (build augmented master, **3 BAM = V3F+V5+V6**) | ~40 min (+10min, 3 BAM) | phaseC region-level CSV | +V3F |
| Step 1 (3 BAM × Augmented LR + LRT) | ~1.3 hr (+0.3 hr, 3 BAM) | Step 0 | +V3F |
| Step 2 (filter sweep, τ 0.5-0.95 by 0.01) | ~45 min | Step 1 | confirmed |
| Step 3 (ΔF1 vs caller F1=0.7166, τ* = max ΔF1) | ~30 min | Step 2 | confirmed |
| Step 4 (mechanism + **PubMed literature search**) | ~3.5 hr (+2 hr, citation verify) | Step 1 (parallel) | +literature |
| Step 5 (decision, 手動 review 後寫 SoT) | ~30 min (+ user review) | Step 3+4 | manual SoT |
| **HCC1395 pilot 合計** | **~6.5 hr** (+2 hr 文獻 + V3F) | ~1 工作日 | 4.5 → 6.5 hr |

**Multi-agent 設計**：
- Agent F (Step 0): build augmented master TSV
- Agent G (Step 1+2 平行): LR + filter sweep
- Agent H (Step 3): ΔF1 計算
- Agent I (Step 4 平行): mechanism brainstorm + KB query
- Coordinator (主 session): Step 5 decision + paper framing

---

## 嚴謹度防護（scientific-rigor §rigor）

| 防護 | 章節 |
|------|------|
| Pre-registration 5 H 事先寫死 | §7.1 |
| NO-GO 條件不可事後改 | §7.1 Hard Gate |
| Confound guard within-cell within-AF OLS | §4 DAG + memory `feedback_L2_collider_bias.md` |
| 5-fold cross-validation hold-out | §5 |
| BH-FDR q<0.05 multi-cell 校正 | §3 Effect size |
| V5 vs V6 sanity check (H5) | §5 對照 |
| 12 methylation features 全 join (not cherry-pick) | §6 消融 |
| Mechanism hypothesis relaxed (用戶接受) | §4 mechanism — but lowered bar |
| Effect size calibration (ΔF1 +0.005 = marginal) | §3 Cohen ribbon |

---

## Resolved Decisions (2026-05-18 用戶確認 4 條)

1. ✅ **τ sweep**: 0.5-0.95 by 0.01 (45 個 τ values)，τ* = **max ΔF1** (最 aggressive)
2. ✅ **Step 4 mechanism**: 含 **PubMed literature search** (cis-mQTL / ASM，每候選機制 1-2 references) — +2 hr 總時程
3. ✅ **Step 1 LR**: 3 BAM (V3F + V5 + V6) 完整對照 — schema augmented master 100 cols 而非 88 cols；計算量 +30%
4. ✅ **Step 5 NEGATIVE verdict**: **手動更新 SoT** — 不自動觸發 `/conclude-research`；用戶 review 後才寫 INDEX / MEMORY / CURRENT_FOCUS / evidence_ledger

---

## Predecessor reference (v0.3 cycle, completed)

- v0.3 plan: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md` (還在原位)
- v0.3 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md`
- v0.3 HTML preview: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.standalone.html`
- v0.3 INDEX entry: `InterSubMod/docs/experiments/INDEX.md` (2026-05-15 line)
- v0.3 CURRENT_FOCUS entry: `InterSubMod/docs/CURRENT_FOCUS.md` (2026-05-15 區段)
- v0.3 figures: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/figures/synthesis/fig{1..5}_*.png`

---

## 與 CURRENT_FOCUS 2026-05-17 plan 接口

本 cycle 可同步運作於：
- **T2.1 Z-AUTO KDE 跨 4 樣本擴展**（⏳ ⭐3 → ⭐4 升級必要條件）：兩 cycle 共用 V5/V6 master TSV，可平行
- **T1.2 V6 production tag finalize (Hard Gate)**: 不阻擋 — 本 cycle 不修 V6 binary
- **T3.1 Paper outline (Bioinformatics / NAR GB)**: 本 cycle POSITIVE → paper §3 主軸從「characterization framework」擴成「characterization + augmented filter pilot」
- **T4.2 GC/mappability/repeat 新軸 pilot**: 本 cycle 結果決定 T4.2 是否啟動（NEGATIVE → 啟動 T4.2；POSITIVE → T4.2 降權）

---

## Verification & Decision Hooks

執行完 cycle 後，必做（依 scientific-rigor §11 協作圖）：
1. `/auc-confound-guard` 三關 (Gate 1-3) 過 LRT 顯著結果
2. `/validation-protocol` L1-L4 對應 verdict
3. `/run-evaluator` 升 tier 評估（若 ⭐4 候選）
4. `/conclude-research` 寫 evidence_ledger + memory + INDEX
5. Spaced recall: 30 天後檢核此 cycle 結論是否仍 hold
