<!--
建立時間: 2026-05-18
agent: init-research skill
status: initiated
report_class: multi-week project plan
audience: PI / lab member / 自己未來
scope: HCC1395 TP rescue + 4 樣本 cross-sample validation + 5 goal impact assessment
parent_cycle: v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot (2026-05-18 ⭐3)
upstream_reports:
  - InterSubMod/docs/experiments/in_progress/2026/05/20260518_V6_Methyl_Filter_Pilot_01.md
  - InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_findings.md
verdict: TBD — 預期 cross-sample phase 2 升 ⭐4 OR 保持 ⭐3 marginal
last_verified: 2026-05-18
predecessor: v6_bam_tpfp_hp_loh_cn (v0.3 characterization ⭐3 + v1.0 methyl pilot ⭐3-marginal)
-->

# Methylation-Augmented FP Filter — Phase 2 Project Plan

> **Multi-week project**，4 hypotheses，預期 3-5 cycles。
> **Predecessor**: `step5_methyl_filter_pilot` (v1.0) 完成 ⭐3 PARTIAL POSITIVE (2026-05-18)
> **Main goal**: 升 ⭐4 cross-sample validated filter + 評估對 G1-G5 五大目標的 impact

---

## 1. 背景與動機

### 1.1 Predecessor 留下的問題

`step5_methyl_filter_pilot` (v1.0, 2026-05-18) 結果：
- ALL 5 H POSITIVE: H1 LRT 16/30 cells q<0.05 / H2 ΔF1 +0.00242 / H3 FP_rem 98.3% > TP_loss 35% / H4 13 mechanism + 14 PubMed / H5 V5≈V6 Δβ=1.87e-5
- **但 ΔF1 marginal** (+0.00242 < +0.005 Cohen ribbon)
- **35% TP loss** 在 production-readiness 是 unacceptable
- Single-sample HCC1395 only — generalizability 未驗

3 個關鍵 follow-up 需要 multi-cycle 才能完成：
1. **TP rescue analysis**: 35% lost TP 是否有共同特徵可救回 → 升 ΔF1 至 ≥ +0.01
2. **Cross-sample validation**: 4 樣本 ΔF1 一致性 → 升 ⭐4
3. **Mechanism cycle-specific 驗證**: 13 hypothesis 哪個真實主導

### 1.2 為什麼開新 project（不繼續在 v6 cycle 內做）

- 任務跨 ≥3 hypothesis × ≥3 cycles
- Cross-sample phase 2 需要重跑 4 樣本 V3F+V5 ISM (~10 hr 平行)
- Mechanism 驗證需 external DB query (GoDMC / Loyfer atlas)
- 累積 evidence ledger 多 entries → project-level scaffolding 比 single cycle 更合適

### 1.3 與 CURRENT_FOCUS 2026-05-17 plan 對接

本 project 對應：
- **T2.1 Z-AUTO KDE 跨 4 樣本擴展** — 與 H_PHASE2_2 共用 V3F+V5+V6 cross-sample ISM
- **T2.2/T2.3 paper 章節骨架** — 直接餵 paper
- **T3.1 Tool paper outline** — Methylation filter 章節
- **T4.2 GC/mappability/repeat 新軸 pilot** — reactive if H_PHASE2_1/2 NEGATIVE

---

## 2. 假說

### H_PHASE2_1：TP rescue analysis (HCC1395 internal)

**陳述**：在 4 FP-rich aggregated cells (τ*=0.52) 內，35% lost TP 有可識別共同特徵 (S1 caller_af 邊緣 / S2 methylation 弱訊號 / S3 sub-clone / S4 chr-specific)，可設計二級 whitelist rule 使 ΔF1 從 +0.00242 升至 ≥ +0.01。

**前提條件**：
- Step 5c TP rescue analysis 完成（in-flight at v6 cycle）
- Lost TP / Kept TP / Removed FP 四群有顯著特徵差異

**已知 Confound**：
- Multi-test inflation across feature candidates → BH-FDR q<0.05
- Cross-validation overfit → 5-fold CV out-of-fold rescue rule
- Sample-size: lost TP ~21 — Wilson CI 寬

**驗證標準**：
- **Positive**：Rescue rule 在 5-fold CV 平均使 ΔF1 ≥ +0.01 + 不引入 >10% new FP
- **Negative**：Rescue rule 在 hold-out 無 ΔF1 改善 OR 新 FP > 10%

### H_PHASE2_2：Cross-sample validation (4 樣本)

**陳述**：Methylation-augmented filter (Model B + rescue rule) 在 H1437/H2009/HCC1954/HCC1937 上 ≥4/5 ΔF1 > 0 (Wilcoxon paired p<0.05) → ⭐4 升級。

**前提條件**：
- H_PHASE2_1 POSITIVE OR NEUTRAL（rescue 不一定 required）
- phaseC_genome_three_way_with_significance 補跑 4 樣本 V3F+V5+V6 ISM (~10 hr parallel)

**已知 Confound**：
- Sample-specific (歷史 HCC1395 chr8 hotspot 跨樣本崩 mean AUC ≤0.641)
- HCC1937 BRCA1 outlier marker rate 0.817 邊緣
- HCC1954 caller ceiling outlier

**驗證標準**：
- **Positive**：≥4/5 ΔF1 > 0 + Wilcoxon paired p<0.05 → ⭐4 升級
- **Negative**：≤2/5 same direction → 保持 ⭐3，paper §Limitations 標 single-sample

### H_PHASE2_3：Mechanism cycle-specific verification

**陳述**：13 個 mechanism hypotheses 中至少 2 個 (cis-mQTL + cancer ASM) 可 cycle-specific 驗證 — top 5 lost-TP / removed-FP region 與 GoDMC mQTL hotspots ≥30% overlap，與 Loyfer 2025 atlas ≥40%。

**驗證標準**：
- **Positive**：≥1 category 達 overlap threshold (Fisher p<0.05 vs random)
- **Negative**：All categories overlap ≤ random expectation

### H_PHASE2_4：Five-goal impact assessment

**陳述**：Methyl-augmented filter 對 5 goals (G1 ASM / G2 subclone / G3 second-hit / G4 normal BAM / G5 evidence panel F1) ≥1 個 incremental signal 貢獻。

**Constraint**: G1/G2/G4 characterization-only / G3 暫緩 / G5 ceiling +0.0112

---

## 3. 方法

### 數據來源

| 數據集 | 路徑 | 描述 |
|--------|------|------|
| Predecessor master TSV | `research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_master_augmented.tsv` | 35,332 × 202 |
| HCC1395 ISM 12 runs | `research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/` | ~25 GB |
| 4 樣本 V6 ISM | `research/paired_priority_bug_audit/phaseD_v6_5sample/` | 待補 V3F+V5 |
| SEQC2 truth | `research/seqc2_cnv_stratification/data/annotated_hcc1395_cnv.tsv` | 31,640 rows |
| GoDMC / Loyfer atlas | **待 cycle 4 WebFetch** | external DB |

### 分析步驟（3-5 cycles）

```
Cycle 1 (H_PHASE2_1): TP rescue analysis
  → Mann-Whitney lost vs kept TP → BH-FDR → source taxonomy → rescue rule → CV ΔF1
  → 驗證: 5-fold CV ΔF1 ≥ +0.01

Cycle 2 (prerequisite): 4 樣本 V3F+V5 ISM rerun
  → 10 hr parallel
  → 產 4 augmented master TSV per sample
  → 驗證: 4 樣本 schema 一致 + V5 vs V6 Pearson r > 0.95

Cycle 3 (H_PHASE2_2): Cross-sample validation
  → Apply Cycle 1 rescue rule + augmented LR
  → Wilcoxon n=5
  → 驗證: ≥4/5 ΔF1 > 0 + p<0.05

Cycle 4 (條件，依 H_PHASE2_2): Mechanism verification
  → bedtools intersect lost-TP × GoDMC / Loyfer
  → Random shuffle null
  → 驗證: ≥1 category Fisher p<0.05 vs random

Cycle 5 (條件): G1-G5 impact assessment
  → 對每 G 跑 incremental analysis
  → 驗證: ≥1 G 目標 measurable incremental signal
```

### 統計方法

- Mann-Whitney U (lost vs kept TP)
- BH-FDR q<0.05
- Logistic regression with statsmodels
- 5-fold cross-validation (sklearn KFold)
- Wilcoxon paired signed-rank (cross-sample)
- Fisher exact for overlap (mechanism)
- Random shuffle permutation null (1000 iter)

---

## 4. 可行性評估

| 因素 | 評估 |
|------|------|
| 數據可用性 | △ HCC1395 ISM 已有；4 樣本 V3F+V5 ISM **待補 ~10 hr 平行** |
| 計算資源 | △ /big7_disk 92% used，剩 3.6T；4 樣本 ISM ~ +25 GB OK |
| 與已有結論衝突 | ✓ 不衝突 — 已與 4 個 NEGATIVE 方向 (pure methyl / G1-G7 / O12 / Option C) 明確差異化 |
| Sample size adequacy | △ Lost TP ~21 → MW 有效但 power 弱；4+1=5 樣本 Wilcoxon exact min p=0.0625 |
| Pre-registration | ✓ 4 H 已 hard-write 於 manifest.yaml + 本 plan |

---

## 5. 已知風險

1. **Lost TP sample size ~21**: 多重比較 inflate FDR
   - **緩解**: BH-FDR q<0.05 + Bonferroni sanity + permutation test
2. **Cross-sample HCC1395-specific signal 跨樣本崩** (歷史 mean AUC ≤0.641)
   - **緩解**: Pre-reg 4/5 同方向 decision threshold；HCC1954/HCC1937 outlier 容忍
3. **GoDMC database tissue mismatch (blood vs breast)**
   - **緩解**: 用 Loyfer 2025 含 cancer tissue + null shuffle
4. **重跑 4 樣本 V3F+V5 ISM ~10 hr + /big7_disk 92% full**
   - **緩解**: 共用 v6_methyl pilot 釋出空間 + TMPDIR=/big7_disk/tmp
5. **Effect size marginal 即使 rescue 救回也可能 < +0.005**
   - **緩解**: Pre-reg threshold +0.01 為 H_PHASE2_1 POSITIVE 條件

---

## 6. 與已關閉方向的差異化（依 Research Direction Guard）

| 已關閉方向 | 本 project 差異 |
|-----------|---------------|
| Pure methylation NO-GO (Beyond-AUC 2026-04-09) | 用 methyl 作 5th-9th covariate **on top of** 4 axis framework，不獨用 |
| TO Germline FP G1-G7 NO-GO | Paired-pileup 模式，多 axis LR-predicted threshold + rescue rule |
| O12 LOH methylation NEGATIVE (collider) | within-cell within-AF OLS confound guard (predecessor 已驗) |
| Option C dual-path NEGATIVE | HP bucket 是 grid 主軸，methyl 補強而非取代 |

---

## 7. 相關檔案

- Predecessor cycle 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260518_V6_Methyl_Filter_Pilot_01.md`
- v0.3 characterization 主報告: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md`
- step5_methyl_filter_pilot/00_PLAN.md
- step5_methyl_filter_pilot/step4_mechanism_candidates.md (13 mechanism hypotheses)
- prior art: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/02_prior_art_notes.md`
- CURRENT_FOCUS 2026-05-17 plan tender-pondering-blossom T1-T4

---

## 8. Pre-registration commitment (per scientific-rigor §7.1, Hard Gate)

- 任 1 個 hypothesis NEGATIVE → 寫 evidence_ledger 不可事後改
- H_PHASE2_2 NEGATIVE → project verdict "marginal pilot, no cross-sample evidence"
- 全 4 hypothesis NEGATIVE → conclude-research as NEGATIVE project

**Pre-reg commit hash**: TBD (本 plan 寫完 commit 後填入)
