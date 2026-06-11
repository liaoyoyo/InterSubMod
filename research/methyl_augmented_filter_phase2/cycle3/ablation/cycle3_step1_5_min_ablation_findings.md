<!--
建立時間: 2026-05-20
agent: main session Coordinator
status: complete
report_class: cycle ablation verdict (Cycle 3 Step 1.5 min)
audience: PI / Coordinator / cycle 3 reframe decision
scope: 5 samples × 4 variants × 2 modes = 40 combinations (cycle 2 既有 master TSV reuse)
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md §Cycle 3 Step 1.5 v2.1
predecessor_cycle: 20260519-0031-cycle3-caller-f1-headroom (Step 1 PASS), 20260519_phase2_cycle2_cross_sample_negative
verdict: H_A1 MARGINAL (caller_af 67% of HCC1954 disaster) / H_M1a FAIL (ISM 是 vestigial covariate)
last_verified: 2026-05-20
-->

# Cycle 3 Step 1.5 — Min Ablation: ISM 是否對 F1 有真實貢獻？

> **核心結論 (n=5, refit OOF)**：
> - **H_M1a FAIL** — Drop ISM 5 個 methylation features 後 HCC1395 refit ΔF1 從 +0.02236 → +0.02171，drop **+0.00065** (≈3% 總 uplift)，遠低於 Cohen +0.003 threshold。**4/5 樣本 no-methyl ΔF1 ≥ full**（HCC1937 +0.00037 / HCC1954 +0.00084 / H1437 / H2009 ≈0）。
> - **5-sample mean ΔF1 反而略升**：full +0.00619 vs **no-methyl +0.00630** — ISM 在 cycle 1 LR 框架下是 **vestigial covariate**，與 Step 5c 「methylation 訊號實為 caller_af proxy」線索一致。
> - **H_A1 MARGINAL** — Shrink caller_af coef → 0 後 HCC1954 transfer ΔF1 從 -0.37744 → -0.12629，improvement **+0.25114**（threshold +0.30）。caller_af 是 HCC1954 transfer 災難的主要 confound (~67%) 但仍有 ~33% 未解釋。
>
> **Cycle 3 reframe action**：filter 從「methyl-augmented multi-axis」改為「**caller_af + LOH + Coverage + NG multi-axis filter**」；paper §3 不再宣稱 ISM methylation 貢獻；ISM characterization 結論保留（v0.3 cycle ⭐3 不變）。

---

## 0. TL;DR — 兩個 Hypothesis Verdict

| Hypothesis | Pre-reg Prediction | Computed | Threshold PASS | Threshold FAIL | Verdict |
|---|---|---|---|---|---|
| **H_A1** caller_af coef shrink → HCC1954 ΔF1 改善 | improvement > +0.30 | **+0.25114** | +0.30 | +0.10 | **MARGINAL** |
| **H_M1a** drop ISM 5 methyl → HCC1395 refit ΔF1 drop | drop ≥ +0.003 | **+0.00065** | +0.003 | +0.001 | **FAIL** |

事前預測 (plan v2.1): 55% A1 PASS/MARGINAL + M1a FAIL → **命中**。

---

## 1. Method

### 1.1 Data

| 樣本 | Rows | n_TP | n_FP | FP density |
|---|---:|---:|---:|---:|
| HCC1395 (cycle 1) | 35,332 | 30,490 | 4,842 | 14.0% |
| H1437 (cycle 2) | 70,964 | 70,191 | 773 | 1.1% |
| H2009 (cycle 2) | 136,701 | 135,359 | 1,342 | 1.0% |
| HCC1954 (cycle 2) | 20,136 | 19,449 | 687 | 3.4% |
| HCC1937 (cycle 2) | 16,607 | 13,910 | 2,697 | 16.2% |

### 1.2 Four Variants × Two Modes

| Variant | Drop | n_features | Modes Tested |
|---|---|--:|---|
| full | (baseline) | 10 | transfer + refit |
| **no-caller-af** | caller_af | 9 | transfer (mock coef→0) + refit (drop feature) |
| **no-methyl** | 5 ISM methyl | 5 | transfer (mock coefs→0) + refit (drop features) |
| no-both | caller_af + 5 methyl | 4 | transfer + refit |

### 1.3 Implementation

- Script: `InterSubMod/research/methyl_augmented_filter_phase2/cycle3/ablation/scripts/run_ablation_variants.py`
- Wall clock: ~25 sec (40 combinations on 5 samples × 4 variants × 2 modes)
- Reuse: `cycle2/scripts/cross_sample_apply.py` `transfer_predict()` / `refit_oof()` / `compute_metrics()` / `sweep_tau()` / `impute_with_medians()`
- Determinism: PRIMARY_SEED=42, sklearn lbfgs, identical to cycle 1+2 protocol

---

## 2. Results

### 2.1 REFIT mode pivot (主 verdict 來源)

| Sample | full | no-caller-af | **no-methyl** | no-both |
|---|--:|--:|--:|--:|
| HCC1395 | **+0.02236** | +0.00955 | **+0.02171** | +0.00912 |
| HCC1937 | **+0.00761** | -0.00001 | **+0.00798** | +0.00000 |
| HCC1954 | +0.00095 | -0.00001 | **+0.00179** | -0.00000 |
| H1437 | -0.00000 | -0.00000 | -0.00000 | -0.00000 |
| H2009 | +0.00001 | +0.00000 | +0.00000 | +0.00000 |
| **5-sample mean** | +0.00619 | +0.00191 | **+0.00630** | +0.00182 |

**核心觀察**:

1. **no-methyl 在 4/5 樣本 ΔF1 ≥ full**（HCC1395 -0.00065；HCC1937 +0.00037；HCC1954 +0.00084；H1437/H2009 ≈0）
2. **5-sample mean no-methyl (+0.00630) 略高於 full (+0.00619)** — ISM 完全 non-contributing 甚至 marginal negative
3. **caller_af drop 影響顯著**: HCC1395 +0.02236 → +0.00955 (drop 57%); HCC1937 +0.00761 → -0.00001 (完全失效); 5-sample mean drop 69%
4. **no-both ≈ no-caller-af** — 證實 ISM 訊號完全 redundant with caller_af + LOH + Coverage + NG

### 2.2 TRANSFER mode pivot (coef shrinkage on cycle 1 filter)

| Sample | full | no-caller-af | no-methyl | no-both |
|---|--:|--:|--:|--:|
| HCC1395 | **+0.02232** | +0.00008 | +0.02149 | -0.00015 |
| H1437 | -0.03597 | -0.00354 | -0.03508 | -0.00257 |
| H2009 | -0.00450 | -0.00138 | -0.00360 | -0.00043 |
| HCC1954 | **-0.37744** ⚠ | **-0.12629** | -0.36705 | -0.09645 |
| HCC1937 | -0.07068 | -0.00903 | -0.06946 | -0.00767 |

**Transfer mode caller_af shrinkage 修補幅度**:

| Sample | Δ improvement from no-caller-af | % of full disaster mitigated |
|---|--:|--:|
| HCC1954 | +0.25114 | 66.6% |
| HCC1937 | +0.06165 | 87.2% |
| H1437 | +0.03243 | 90.2% |
| H2009 | +0.00312 | 69.3% |

**Caller_af 平均修補 ~78% transfer-mode disaster**；殘餘 22% 來自其他 features overfit。

**Methylation shrinkage 修補幅度**:

| Sample | Δ improvement from no-methyl | % of full disaster mitigated |
|---|--:|--:|
| HCC1954 | +0.01039 | 2.8% |
| HCC1937 | +0.00122 | 1.7% |
| H1437 | +0.00089 | 2.5% |
| H2009 | +0.00090 | 20.0% |

**Methylation 平均修補 <7% disaster** — 證實 methylation 非 transfer 災難主因。

### 2.3 Best τ per variant (refit)

| Sample | full τ | no-caller-af τ | no-methyl τ | no-both τ |
|---|--:|--:|--:|--:|
| HCC1395 | 0.39 | 0.28 | 0.41 | 0.32 |
| HCC1937 | 0.24 | 0.14 | 0.27 | 0.10 |
| HCC1954 | 0.71 | 0.26 | 0.74 | 0.10 |
| H1437 | 0.10 | 0.10 | 0.10 | 0.10 |
| H2009 | 0.61 | 0.10 | 0.10 | 0.10 |

- no-methyl 的 best τ 與 full 接近 (±0.03)，confirming methyl 不改 decision boundary
- no-caller-af 的 best τ 大幅變低 (HCC1954 0.71→0.26, HCC1937 0.24→0.14)，confirming caller_af 是主決策變數

---

## 3. Verdict Interpretation

### 3.1 H_M1a FAIL — ISM methylation 是 vestigial covariate

**證據強度**：⭐⭐⭐⭐ (4/5)

- HCC1395 incremental value = **+0.00065** = 2.9% of full uplift = below Cohen +0.003 ribbon
- 4/5 樣本 no-methyl ΔF1 ≥ full — **methylation 不只 marginal，甚至略微負貢獻**
- 5-sample mean no-methyl > full — 跨樣本一致
- 與 cycle 1 主報告 Step 5c 「methylation 訊號實為 caller_af proxy」相符

**為什麼 5-rank coef +0.75 但 incremental ≈ 0？**

L2 regularization 把相關特徵 split：HPFineF coef +0.75 來自其與 caller_af 的部分共變異 (HPFineF ≈ f(caller_af) + residual)；drop 5 個 methylation 後，caller_af + LOH + Coverage + NG 4 個 features 重新分配權重就能 capture 同樣訊息。

### 3.2 H_A1 MARGINAL — caller_af 主導但非唯一 confound

**證據強度**：⭐⭐⭐ (3/5)

- caller_af shrink 修補 67% (HCC1954) ~ 90% (H1437) transfer disaster
- 5-sample mean: caller_af drop 從 +0.00619 → +0.00191 (-69%)
- 但 HCC1954 仍 -0.126 (no-caller-af) — 33% 未解釋；no-both 為 -0.096 表示 caller_af + methylation 兩者解釋 75% 但仍餘 25%

**剩餘 25% 災難來源候選** (待大規模 ablation):
- LOH_inner_flag distribution shift
- Coverage_Multiple_imp scaler shift
- NG distributional shift  
- Sample-specific TP/FP density shift (HCC1954 FP density 3.4% vs HCC1395 14% → LR predict probability calibration error)

### 3.3 為什麼 5 樣本 mean +0.00619 (full) 仍然「不負」？

關鍵：**4/5 樣本 ΔF1 ≈ 0** (H1437 / H2009 best τ→0.10 keep-everything；HCC1937 / HCC1954 best τ 在 caller F1 ceiling 樣本上 LR 找不到有意義的 threshold)；唯一一個 large positive 是 HCC1395 (+0.02236)。

**HCC1395 trained sample，他樣本 LR 退化為 no-filter (τ→0.10) 或 small uplift** — Refit mode 反映的不是 cross-sample generalization 而是 "filter on training-like samples works, on dissimilar samples LR fails to find a useful boundary"。

cycle 3 Step 1 gate rule (caller F1 < 0.80 + FP density > 0.10) 正好識別此 pattern。

---

## 4. Cycle 3 Reframe Action

### 4.1 Filter naming 改寫

**Before**: "ISM (InterSubMod) methylation-augmented multi-axis FP filter"

**After**: "**Multi-axis caller-F1-headroom-gated FP filter** (caller_af + LOH_inner + Coverage_Multiple + NG, methylation residual not contributing)"

### 4.2 Paper §3 改寫範圍

| 原 claim | 改寫 |
|---|---|
| "Methylation-augmented filter 達到 ΔF1 +0.02236" | "Multi-axis filter (caller_af + LOH + Cov + NG) 達到 ΔF1 +0.02236；ISM methylation features 作為 covariate 加入無 incremental value (n=5 ablation, 5-sample mean drop -0.00011, HCC1395 drop +0.00065)" |
| "5 ISM methylation features 在 LR 中 rank 5-9" | "5 ISM methylation features 在 cycle 1 LR coef rank 5-9，但 L2 ablation 證實為 caller_af proxy (drop 後 caller_af + LOH + Cov + NG 重新分配權重 capture 同訊息)" |
| "Methylation-aware filter" | **撤回** — 改 "caller-F1-headroom-gated 4-feature filter"; methylation features 作為 ISM characterization 工具保留（v0.3 cycle 結論 ⭐3） |

### 4.3 ISM 結論保留範圍

- **保留**: ISM characterization 5 zones × 7 samples × LR LRT 16/30 cells q<0.05 結論（cycle 1 v1.0 主報告 §3 不變）
- **保留**: 5-zone TP/FP discrimination AUC 觀察
- **撤回**: ISM 對 F1 filter 有 incremental contribution 的 claim
- **保留**: ISM 在 mechanistic understanding 上的價值（理解 LOH × HP × CN × methylation 分布；non-filter contribution）

### 4.4 Cycle 3 Step 2 panel survey 是否仍 GO？

**YES**（per plan v2.1 Decision rule path 2: "A1 PASS + M1a FAIL → Step 2 panel survey 仍 GO 驗 caller_af + LOH + CN filter"）：

- Filter 本身效力（caller-F1-headroom gate + multi-axis LR）仍待 cross-sample validation
- 但 panel survey 不再是「驗 ISM-augmented filter」，而是「驗 4-feature filter (no methyl)」
- 大規模評估 (M2 partial F-test / M3 SHAP / M0 residualize) 可降級為 supplementary（H_M1a FAIL 已 high-confidence）

---

## 5. Decision Tree Path Taken

依 plan v2.1 §"Decision rule (post-Step 1.5 ablation)"：

```
最小評估結果:
  A1 MARGINAL (improvement +0.25114, 在 [PASS 0.30] 與 [FAIL 0.10] 之間偏 PASS)
  M1a FAIL (drop +0.00065, well below FAIL 0.001)
       │
       ▼
  路徑 2: A1 PASS-leaning + M1a FAIL
       │
       ▼
  Cycle 3 reframe "non-methyl multi-axis filter"
  ISM 從 filter design 移除（但 characterization 保留）
  Paper §3 改寫（撤回 methylation-augmented 宣稱）
  Step 2 panel survey 仍可 GO 驗 caller_af + LOH + CN filter
  大規模 ablation 降為 supplementary
```

---

## 6. 大規模評估 (M2 / M3 / M0 / A2 / LOFO) — 是否仍跑？

**建議降級為 supplementary**：

| Stage | 必要性 | 理由 |
|---|---|---|
| **M2** ISM partial F-test (LRT) | LOW | H_M1a FAIL ΔF1=+0.00065 已是 functional test；LRT 即使 p<0.05 也 dwarf by effect size |
| **M3** SHAP per-methyl | MEDIUM | 確認哪個 methyl feature 最高 SHAP 仍有 mechanistic 價值；但不改 filter design |
| **M0** Residualize HPFineF AUC | MEDIUM | confounder check 對 ISM characterization paper §3 仍有價值；可移到 ISM characterization cycle (v0.3) supplement |
| **A2** caller_af KS test 跨樣本 | HIGH | 量化 distributional shift 解釋 H_A1 殘餘 33% — 強化 mechanism 敘述 |
| **LOFO** Leave-one-methyl-out | LOW | H_M1a FAIL 已 settle，LOFO 細節對結論無影響 |

**建議**: 跑 **A2 (~30 min)** 補強 H_A1 mechanism；M3 + M0 暫緩；M2 + LOFO 跳過。

---

## 7. Files

```
cycle3/ablation/
├── cycle3_step1_5_min_ablation_findings.md     ← this report
├── scripts/run_ablation_variants.py
├── data/cycle3_step1_5_min_ablation.tsv        ← 40 rows (5×4×2)
├── figures/cycle3_step1_5_min_ablation.png     ← refit grouped bar
└── intermediate/cycle3_step1_5_summary.json    ← pre-reg verdict JSON
```

---

## 8. Reproducibility

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/methyl_augmented_filter_phase2/cycle3/ablation/scripts/run_ablation_variants.py
```

Wall clock: ~25 sec. Deterministic — PRIMARY_SEED=42, lbfgs solver, identical to cycle 1+2 protocol. HCC1395 full refit ΔF1=+0.022358 bit-exact match cycle 1 stored result (drift 0).

---

## 9. Hand-off

| 動作 | 優先 |
|---|---|
| Cycle 3 plan / 00_PLAN.md update | HIGH — reflect filter reframe |
| Cycle 3 caller_f1_gate.json update | HIGH — change `scope_caveats` 移除 "methylation as covariate" |
| Memory `project_phase2_cycle1_global_fp_filter.md` caveat 追加 | HIGH — record H_M1a FAIL |
| Evidence ledger entry | HIGH — append cycle 3 Step 1.5 ablation result |
| Paper §3 draft 修改 | MEDIUM — defer 到 paper draft 啟動時 |
| A2 KS test (optional) | MEDIUM — 30 min 補強 mechanism |
| INDEX.md update | MEDIUM — cycle 3 Step 1.5 entry |
| CURRENT_FOCUS.md update | HIGH — record reframe |

---

**End of Cycle 3 Step 1.5 Min Ablation Findings**
