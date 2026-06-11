# Phase 2 Cycle 1 — Global FP Filter Redesign + Cross-Sample Validation

> Plan v2.0 — 2026-05-18
> **Predecessor**: v1.0 step5_methyl_filter_pilot ⭐3 marginal + Step 5c rescue NEGATIVE + Step 5d audit GREEN
> **Project**: `InterSubMod/research/methyl_augmented_filter_phase2/` (init-research 2026-05-18)
> **Different task**: pivot from "TP rescue" to "global FP exploration + heterogeneous threshold filter"
> **Goal**: 達 ΔF1 ≥ +0.01 (Cohen 小 effect) + cross-sample ≥4/5 same direction

---

## Context

### v1.0 cycle 結果摘要

- ALL 5 H POSITIVE 但 ΔF1 = **+0.00242 marginal** (< +0.005 Cohen ribbon)
- Step 5c rescue NEGATIVE: 95.2% lost TP 是 low-AF subclone + caller_af 與 FP 嚴重重疊，methylation 訊號其實是 caller_af proxy → rescue rule 必然 reimport FP
- Step 5d robustness GREEN with caveats: ΔF1 stable (std 2e-5), 4 unique LRT cells, NaN median 0.503 borderline
- v1.0 Step 2/3 「powered cell gate」過嚴 → 4 cells covers **只 7% all FP** (350/4,842)

### 用戶 pivot direction (2026-05-18)

「在某些 FP 比例高的組合區域用較嚴格條件去除 FP，少量 TP 影響，達到更好 F1 提升超過 1%？」

→ Pivot strategy: **多個 FP-rich zone heterogeneous threshold filter** (each zone 獨立 τ) + global LR (不限 powered cell gate)

### 與歷史 NEGATIVE 方向差異化 (Research Direction Guard)

| 已關閉 | 本 cycle 差異化 |
|--------|--------------|
| LOH binary filter 10/10 (2026-04-06) | 不用 binary；用 per-zone LR-predicted threshold (heterogeneous) |
| TO Germline FP G1-G7 (AUC<0.64) | 在 paired-pileup + multi-axis LR + methyl covariate |
| Pure methylation NO-GO (Beyond-AUC 2026-04-09) | Methyl 作 5th-9th covariate，不獨用 |
| Step 5c TP rescue NEGATIVE (2026-05-18) | 不嘗試 rescue lost TP；改攻 high-AF FP zone |

---

## Pre-Registration (Hard Gate per scientific-rigor §7.1)

| ID | Prediction | Falsification | Decision threshold |
|----|-----------|---------------|-------------------|
| **H_C1_1** | Top 10 FP-rich cells covers ≥70% all FP | < 70% | concentration map D1 |
| **H_C1_2** | Global LR (35,332 row, no cell gate) max ΔF1 > v1.0 cell-gated +0.00242 | ΔF1 ≤ +0.00242 | global LR sweep D2 |
| **H_C1_3** | Heterogeneous per-zone threshold aggregate ΔF1 ≥ +0.01 | ΔF1 < +0.01 | aggregate D3 |
| **H_C1_4** | High-AF (caller_af>0.3) zone incremental ΔF1 ≥ +0.003 | < +0.003 | high-AF D4 |
| **H_C1_5** | Cross-sample n=5: ≥4/5 ΔF1 > 0 + Wilcoxon paired p<0.05 | ≤2/5 same direction OR p ≥ 0.05 | cross-sample test |
| **H_C1_6** | V3F/V5/V6 BAM 在 cross-sample ΔF1 一致 (max var < 0.005) | max var > 0.01 | sanity |

**NO-GO 條件**:
- H_C1_2 + H_C1_3 同時 FAIL → pivot to T4.2 GC/mappability/repeat 新軸 cycle (per CURRENT_FOCUS)
- H_C1_5 FAIL → cycle 結論「HCC1395-specific, not generalizable」，⭐3 保持

---

## Plan Structure — Track A 序列 → Track B

### Track A: HCC1395 重設計 filter (預期 ~3-4 hr active)

#### Step 0 — Global FP exploration audit (原 cycle 1.5)

**4 Stages in 1 script** (`scripts/global_fp_audit.py`):

| Stage | Goal | ETA |
|-------|------|-----|
| Stage 1 — Concentration map | (caller_af bin × LOH inner/outer × Coverage_Multiple bin × chr8 vs other) 120-cell FP heatmap | 10 min |
| Stage 2 — High-AF FP profile | caller_af > 0.3 FP 數量 + cells 分布 + 是否 chr8 amplicon 主導 | 5 min |
| Stage 3 — Global LR sweep | 全 35,332 row 跑 LR (no cell gate) + τ sweep [0.10, 0.95] | 20 min |
| Stage 4 — Heterogeneous threshold | Top 20 FP-rich cells 各 cell-specific τ + aggregate | 10 min |

**Pre-registered D1-D4 decision criteria**:
- D1: top 10 FP-rich cells covers ≥ 70% all FP → H_C1_1
- D2: global LR max ΔF1 > +0.00242 → H_C1_2
- D3: heterogeneous aggregate ΔF1 ≥ +0.01 → H_C1_3
- D4: high-AF zone incremental ΔF1 ≥ +0.003 → H_C1_4

**輸出**: `methyl_augmented_filter_phase2/cycle1_step0_global_fp_audit.md` + figures/tsv

#### Step 1 — Filter design (依 D1-D4 結果決定 path)

| D2 | D3 | Path |
|----|-----|------|
| PASS | PASS | **Path A** — global LR + heterogeneous threshold 混合 |
| PASS | FAIL | **Path B** — pure global LR (簡單) |
| FAIL | PASS | **Path C** — pure heterogeneous (FP-rich zone 集中) |
| FAIL | FAIL | **NO-GO** — pivot T4.2 GC/mappability cycle |

依 path 寫 `scripts/cycle1_filter_design.py`，輸出 final filter rule (τ + cell list + covariate set)

#### Step 2 — HCC1395 ΔF1 verdict (Track A 收尾)

- Apply final filter to HCC1395 master TSV
- 5-fold CV out-of-fold ΔF1 (anti-optimism per Step 5d)
- 對 Step 5c S1-S4 4 sources check rescue rule 是否 helpful (negative control)
- 輸出: `cycle1_track_a_findings.md` (H_C1_2/3/4 verdicts)

### Track B: Cross-sample ISM rerun (background parallel from Step 0)

#### Step B1 — 4 樣本 V3F+V5 ISM 補跑 (~10 hr background)

**任務**: phaseD_v6_5sample/ 已有 4 樣本 V6 ISM 但無 V3F+V5
- 補 V3F + V5 ISM with significance (per Step -1 amendment script `run_phaseC_genome_three_way_with_significance.sh`)
- 8 ISM runs per sample × 4 樣本 = 32 runs
- 但 V6 已存在所以只需 V3F + V5 = 4 runs/sample × 4 samples = **16 runs**
- 預期 ~12 min/run × 16 = ~3.2 hr 序列 OR ~1 hr parallel 4 concurrent
- 寫到 `phaseC_genome_three_way_with_significance_4sample/{sample}/{BAM}_{flag}_{label}/`

**TMPDIR safety**: export TMPDIR=/big7_disk/tmp (per memory feedback_tmp_disk_full_pipeline_pitfall.md)

#### Step B2 — Per-sample augmented master TSV

對 4 樣本各自 build augmented master (per Step 0 of v1.0)：
- 12 ISM runs per sample (V3F+V5+V6 × on/off × tp/fp)
- Join methylation features per sample
- 輸出: `{sample}_master_augmented.tsv` 4 files

#### Step B3 — Cross-sample apply Track A filter

- 對 4 樣本 + HCC1395 = 5 samples 套用 Track A 設計的 filter
- 計算每樣本 ΔF1 vs caller F1
- Wilcoxon paired vs 0 (n=5, min p=0.0625)

#### Step B4 — Verdict H_C1_5 + H_C1_6

- ≥4/5 ΔF1 > 0 + Wilcoxon p<0.05 → POSITIVE 升 ⭐4
- ≤2/5 same direction → HCC1395-specific 保持 ⭐3

---

## Cycle 1 Decision Tree

```
Step 0 (audit) → D1/D2/D3/D4 outcomes
   │
   ├─ D2+D3 FAIL → NO-GO pivot T4.2 (cycle 結束)
   │
   ├─ D2 OR D3 PASS → Step 1+2 Track A filter design + ΔF1 verdict
   │   │
   │   ├─ H_C1_2/3/4 ALL PASS → 進 Track B cross-sample
   │   │
   │   ├─ H_C1_3 FAIL (ΔF1 < +0.01) → 保持 ⭐3 marginal，Track B optional
   │   │
   │   └─ H_C1_2 PASS only (+0.00242 ~ +0.01) → Track B 弱優先
   │
   └─ Track B → H_C1_5 verdict
       │
       ├─ ≥4/5 ΔF1>0 + p<0.05 → ⭐4 cross-sample validated
       │
       └─ ≤2/5 → HCC1395-specific 保持 ⭐3
```

---

## 關鍵檔案路徑

### Input (read-only)
| 用途 | 路徑 |
|------|------|
| v1.0 master TSV | `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv` (35,332 × 202) |
| Phase 2 project plan | `InterSubMod/research/methyl_augmented_filter_phase2/00_PLAN.md` |
| Step 5c TP rescue | `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5c_tp_rescue_analysis.md` |
| Step 5d audit | `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5d_robustness_audit.md` |
| HCC1395 phaseC ISM 12 runs | `InterSubMod/research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/` |
| 4 樣本 V6 ISM (phaseD) | `InterSubMod/research/paired_priority_bug_audit/phaseD_v6_5sample/` |
| V3F + V5 BAM (4 樣本 待跑) | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/{pononly_v3_fixed,threshold_compare/v5_flag}/...tagged.bam` (per-sample) |
| Rerun script | `InterSubMod/research/paired_priority_bug_audit/scripts/run_phaseC_genome_three_way_with_significance.sh` (modify per-sample BAM paths) |
| SEQC2 truth | `InterSubMod/research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv` |

### Output (write to)
| 用途 | 路徑 |
|------|------|
| Cycle 1 主目錄 | `InterSubMod/research/methyl_augmented_filter_phase2/cycle1/` |
| Audit output | `methyl_augmented_filter_phase2/cycle1/cycle1_step0_global_fp_audit.md` + data + figures |
| Track A findings | `methyl_augmented_filter_phase2/cycle1/cycle1_track_a_findings.md` |
| Cross-sample ISM | `phaseC_genome_three_way_with_significance_4sample/` |
| Track B findings | `methyl_augmented_filter_phase2/cycle1/cycle1_track_b_findings.md` |
| Cycle 1 synthesis | `methyl_augmented_filter_phase2/cycle1/cycle1_findings.md` |
| Final report | `InterSubMod/docs/experiments/in_progress/2026/05/20260519_Phase2_Cycle1_Global_FP_Filter_01.md` |

### Reuse 既有 scripts
| 既有 | 怎麼 reuse |
|------|-----------|
| `step5_methyl_filter_pilot/scripts/augmented_lr.py` | Step 1 Path A/B/C 共用 LR fit |
| `step5_methyl_filter_pilot/scripts/filter_sweep.py` | Step 1 τ sweep |
| `step5_methyl_filter_pilot/scripts/delta_f1.py` | Step 2 + Step B3 ΔF1 計算 |
| `step5_methyl_filter_pilot/scripts/build_augmented_master.py` | Step B2 per-sample master TSV build |

---

## 嚴謹度防護

- Pre-registration 6 H hard-write 於本 plan + manifest.yaml
- Within-cell within-AF OLS confound guard (依 step 5d C.7 提示 NaN handling Strategy A+B 雙跑)
- 5-fold StratifiedKFold OOF (anti-optimism)
- BH-FDR q<0.05 全局校正
- Wilcoxon paired n=5 exact min p=0.0625 (cross-sample)
- Counts 雙口徑警示: paper-spec (30490/4842 post-ISM) vs caller variant-level (28509/11606)
- TMPDIR=/big7_disk/tmp safety (Step B1 rerun)

---

## 預期時程

| Step | 估時 | 依賴 |
|------|------|------|
| Step 0 Audit (4 stages) | ~45 min | v1.0 master TSV |
| Step 1 Filter design (依 D2/D3 path) | ~1.5 hr | Step 0 |
| Step 2 HCC1395 ΔF1 verdict | ~30 min | Step 1 |
| Step B1 ISM rerun (16 runs background, 4 concurrent) | ~1 hr | independent (parallel) |
| Step B2 4-sample augmented master | ~30 min | Step B1 |
| Step B3 Cross-sample filter apply | ~45 min | Step 1 + B2 |
| Step B4 Wilcoxon n=5 verdict | ~15 min | Step B3 |
| **Track A subtotal** | **~3 hr** | — |
| **Track B subtotal** | **~2.5 hr** (1 hr ISM background + 1.5 hr active) | parallel from Step 0 |
| **Cycle 1 total** | **~4-5 hr active** | (Track B background overlap reduces wall clock) |

---

## Multi-agent fan-out 設計

| Agent | 任務 | Background? |
|-------|------|------------|
| Agent A1 | Step 0 audit (4 stages) | foreground (need result for Step 1) |
| Agent A2 | Step 1+2 filter design + ΔF1 (Track A) | foreground after Step 0 |
| Agent B1 | Step B1 4-sample V3F+V5 ISM rerun | **background** (independent, ~1 hr parallel) |
| Agent B2 | Step B2 build 4-sample master TSV | foreground after Agent B1 |
| Agent B3 | Step B3+B4 cross-sample ΔF1 + Wilcoxon | foreground after Agent A2 + B2 |
| Coordinator (main session) | Cycle 1 synthesis + paper framing | — |

---

## Out-of-scope (避免 scope creep)

- ❌ H_PHASE2_3 mechanism cycle-specific 驗證 (GoDMC + Loyfer atlas) → cycle 2
- ❌ H_PHASE2_4 5-goal impact assessment → cycle 3+
- ❌ V3F/V5/V6 三 BAM 對 4 樣本完整比較 (本 cycle 只用 V6 + V3F + V5)
- ❌ COLO829 (truth set 0600 deferred)
- ❌ HG002 non-cancer pilot (T4.4 reactive only)
- ❌ Read-level methylation embedding (類 ROCIT)
- ❌ C++ pipeline 修改

---

## Verification (end-to-end)

### 階段 verification
1. Step 0 通過: 4 stages 全完成 + D1-D4 verdicts 明確
2. Step 1 通過: filter rule 寫好 + path 明確 (A/B/C/NO-GO)
3. Step 2 通過: H_C1_2/3/4 verdicts 寫入 ledger
4. Step B1 通過: 16 ISM runs 全成功 row count ≈ 30,490 TP / 4,842 FP per run
5. Step B2 通過: 4 augmented master TSV 完整 (~178 cols each)
6. Step B3 通過: 5-sample ΔF1 計算完成
7. Step B4 通過: H_C1_5/6 verdicts 寫入 ledger

### 命令模板
```bash
# Step 0
python methyl_augmented_filter_phase2/cycle1/scripts/global_fp_audit.py \
    --master step5_methyl_filter_pilot/step5_master_augmented.tsv \
    --output methyl_augmented_filter_phase2/cycle1/

# Step B1 (4 samples × V3F+V5 = 16 runs background)
for sample in H1437 H2009 HCC1954 HCC1937; do
  for bam in V3F V5; do
    bash phaseC_with_significance_sample.sh $sample $bam off tp &
    bash phaseC_with_significance_sample.sh $sample $bam off fp &
  done
done; wait

# Step B3 cross-sample filter apply
python methyl_augmented_filter_phase2/cycle1/scripts/cross_sample_apply.py \
    --filter-rule cycle1_track_a_filter.json \
    --samples HCC1395,H1437,H2009,HCC1954,HCC1937
```

---

## Decision rule (post-cycle 1)

1. **H_C1_3 PASS + H_C1_5 PASS** → ⭐4 cross-sample validated, paper §3 主軸升級
2. **H_C1_3 PASS + H_C1_5 FAIL** → "HCC1395-specific filter, marginal cross-sample"，paper case study
3. **H_C1_3 FAIL + H_C1_2 PASS** → ⭐3 stays, ΔF1 +0.005-0.01 incremental，cross-sample optional
4. **H_C1_2/3 BOTH FAIL** → NO-GO pivot T4.2 GC/mappability/repeat cycle (per CURRENT_FOCUS)

依 plan v1.0 resolved decision #4: **手動 SoT update** — 用戶 review cycle1_findings.md 後決定 INDEX / CURRENT_FOCUS / evidence_ledger 更新
