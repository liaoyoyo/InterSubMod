<!--
建立時間: 2026-05-20
agent: main session Coordinator
status: complete
report_class: cycle pre-exploration priority synthesis
audience: PI / user 決定 Trial scope
scope: Cycle 4 Stage 1+2 結束 — 13 candidates ROI 排序 + Trial scope 建議
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md v3.1 §Cycle 4 Stage 1+2
predecessor: 20260520_cycle3_step1_5_ism_ablation_vestigial
verdict: Stage 1+2 mapping complete — 4 highest-ROI candidates 建議進 Trial A
last_verified: 2026-05-20
-->

# Cycle 4 Stage 1+2 — 13 Candidates ROI 排序 + Decision Matrix

> **本檔職責**: 整合 Stage 1 mapping (12 features + 10 algorithms + 13 observations inventory) → Stage 2 ROI 排序 → 建議 user 看完 HTML 後決定哪些進 Trial A/B/C

---

## 0. Stage 1 統計總覽

| 維度 | Tried | Not tried | Reserve / low ROI | 總 |
|------|---:|---:|---:|---:|
| Methylation features | 5 (F_A1-F_A5) | 7 (F_B1-F_B7) | — | 12 |
| Algorithms | 3 (A_T1-A_T3) | 7 (A_U1-A_U8) | — | 10 |
| Observations | 5 (B_T1-B_T5) | 4 priority (O1-O4) | 4 reserve (O5-O8) | 13 |
| Python 新 features | 0 | 5 (F1-F5) | — | 5 |

**Coverage**: 三個維度合計 5/12 + 3/10 + 5/13 = 30% 試過 → **70% 未試**。Cycle 3 結論「ISM vestigial」僅 narrow window，可探索空間大。

---

## 1. 13 Candidates ROI 排序

依「對 H_C4 4 個 pre-reg hypothesis 直接挑戰能力 × 計算成本 × 預期回報」三向評分（1-5 各維度）。

| Rank | ID | 名稱 | 對應 H | 挑戰力 | 成本 (days) | 預期回報 | ROI score | 進 Trial? |
|------|----|------|-------|:---:|:---:|:---:|:---:|---|
| 1 | O3 | Interaction LR + LRT | H_C4_INT | 5 | 1.0 | 4 | **20** | ✅ Trial A |
| 2 | F1 | Methylation entropy per region | H_C4_NEW_FEAT | 4 | 0.5 | 4 | **32** | ✅ Trial A |
| 3 | F2 | Per-CpG × HP χ² Fisher aggregation | H_C4_NEW_FEAT | 4 | 0.5 | 3 | **24** | ✅ Trial A |
| 4 | O1 | Methylation entropy histogram TP vs FP + AUC | H_C4_NEW_FEAT | 3 | 0.5 | 3 | **18** | ✅ Trial A |
| 5 | F4 | Methylation variance per region | H_C4_NEW_FEAT | 3 | 0.5 | 2 | 12 | ✅ Trial A (cheap) |
| 6 | A_U1 | Random Forest | H_C4_NL | 4 | 1.5 | 3 | **8** | ✅ Trial B |
| 7 | A_U2 | XGBoost / LightGBM | H_C4_NL | 4 | 1.5 | 3 | **8** | ✅ Trial B |
| 8 | O4/A_U5 | Per-zone heterogeneous LR | H_C4_ZONE | 4 | 2.0 | 3 | **6** | ✅ Trial C |
| 9 | F3 | CpG density per region | H_C4_NEW_FEAT | 2 | 0.25 | 2 | **16** | ✅ Trial A (cheap) |
| 10 | F5 | HP1/HP2 methylation correlation | H_C4_NEW_FEAT | 2 | 0.5 | 2 | 8 | ✅ Trial A (cheap) |
| 11 | O2 | Per-CpG × HP × ALT-base χ² | H_C4_NEW_FEAT | 3 | 1.0 | 2 | 6 | Trial A (optional) |
| 12 | A_U3 | Polynomial / Spline LR | H_C4_NL | 3 | 1.5 | 2 | 4 | Trial B (optional) |
| 13 | A_U8 | Elastic Net | refining LR | 2 | 1.0 | 2 | 4 | reserve |

ROI score = 挑戰力 × 預期回報 / 成本（注意：F1-F5 cheap so ROI 偏高）。

### Reserve (Trial A FAIL 後不執行)

| ID | 名稱 | 成本 |
|----|------|---:|
| O5 | Distance matrix eigenvector clustering | 1.5 |
| O6 | Strand symmetry test | 1.0 |
| O7 | Outlier read detection | 1.5 |
| O8 | Per-CpG × caller_af correlation | 0.5 (已知 caller_af proxy 風險) |
| A_U6 | Neural network MLP | 2.0 |
| A_U7 | Bayesian shrinkage | 2.0 |

---

## 2. Trial Scope 建議

### Trial A — Interaction LR + 5 Python 新 features（建議優先 ~1.5 day）

| 內容 | 對應 |
|------|------|
| Baseline LR (cycle 1 10-feature) | reference |
| Interaction LR: + caller_af × {HPMergedDelta, HPFineF, NME_imbalance, Epipoly_Delta, ClusterPermanovaF} | O3 (H_C4_INT) |
| 5 Python 新 features 加入 LR (F1-F5) | F1/F2/F3/F4/F5 (H_C4_NEW_FEAT) |
| LRT vs baseline | overall LR significance |
| Per-feature AUC + Cohen's d | individual contribution |

**Decision after Trial A**:
- H_C4_INT PASS + H_C4_NEW_FEAT PASS → 進 Trial B 驗 non-linear capture
- H_C4_INT PASS only → 進 Trial B optional
- H_C4_NEW_FEAT PASS only → 確認新 feature 哪一個有效；可能跳 Trial B 直接做 stratified analysis
- 全 FAIL → cycle 4 結論 final NEGATIVE，archive methyl_filter_phase2 軌

### Trial B — Random Forest + XGBoost + SHAP（conditional, ~2 day）

| 內容 | 對應 |
|------|------|
| sklearn RandomForestClassifier (n_est=500) | A_U1 |
| xgboost.XGBClassifier (learning_rate=0.1) | A_U2 |
| 5-fold StratifiedKFold OOF | anti-optimism |
| SHAP per-feature contribution | mechanism check |

### Trial C — Per-zone heterogeneous LR（conditional, ~2 day）

| 內容 | 對應 |
|------|------|
| 5 zones split (Z-OCH / Z-GL / Z-AUTO / Z-CHR8 / Z-other) | O4/A_U5 |
| Per-zone L2 LR + swept τ | per-zone optimization |
| Aggregate ΔF1 + zone-stratified Wilcoxon | H_C4_ZONE |

---

## 3. 與 Cycle 3 結論的關係

| Cycle 3 結論 | Cycle 4 挑戰路徑 |
|---|---|
| ISM 5 features 在 cycle 1 LR vestigial | Trial A 加 interaction + 5 new features 看是否 conditional 有效 |
| caller_af 是 HCC1395 主導 (74.6%) | Trial B 看 RF/XGB 是否 capture 非 caller_af 訊號 |
| 全局 LR ⭐3 strict subset | Trial C 看 per-zone heterogeneous 是否 zone-specific signal |
| Filter 4-feature production | 若 Trial A/B/C 全 FAIL → 確認 4-feature 已 saturate |

---

## 4. Out-of-scope (per plan v3.1)

- ❌ BAM-level re-extraction (5mC/5hmC, strand-specific)
- ❌ Cross-sample n=5 Wilcoxon per trial（集中 HCC1395 + HCC1937）
- ❌ Read-level embedding
- ❌ Subclonal deconvolution
- ❌ C++ pipeline 修改

---

## 5. Files Inventory

```
cycle4/
├── stage1_mapping/
│   ├── feature_inventory.tsv      (12 features × 9 cols)
│   ├── algorithm_inventory.tsv    (10 algorithms × 8 cols)
│   └── observation_inventory.tsv  (13 observations × 8 cols)
├── stage2_priority/
│   ├── cycle4_priority_synthesis.md     ← this report
│   └── cycle4_priority.standalone.html  ← (next step: PI HTML viz)
├── figures/
│   ├── stage1_gap_matrix.png      (3-panel feature×algo×obs status)
│   └── stage1_features_boxplot.png (12 methyl features × 5 samples)
└── scripts/
    └── stage1_inventory_figures.py
```

---

## 6. Pre-reg Predictions (for HARKing防護)

事前主觀估計 H_C4 PASS 機率（per plan v3.1）：

| H | Prior | Reason |
|---|---|---|
| H_C4_INT | ~40% | Interaction effect 可能 capture 但 effect 小 |
| H_C4_NL | ~30% | RF/XGB 通常 > 線性 +0.005 但 small data 過擬合 |
| H_C4_ZONE | ~35% | Zone-specific signal 可能在 chr8/LOH zones |
| H_C4_NEW_FEAT | ~25% | 新 feature 也可能 caller_af proxy |
| 至少 1 H PASS | ~75% | 組合機率 (1 - 0.6 × 0.7 × 0.65 × 0.75) |

如最終 cycle 4 verdict 與 prior 不同 → 寫 "expected vs observed" 對照防 confirmation bias。

---

**End of Cycle 4 Stage 1+2 Priority Synthesis**
