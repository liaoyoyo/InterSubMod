<!--
建立時間: 2026-05-21
agent: Phase A4-Ext subagent (inventory only, no execution)
status: inventory_complete
report_class: A4-Ext algorithm scan / PI Goal 4 due-diligence companion
audience: PI / Phase 2 cycle 4 planner / paper §3 framing
scope: 盤點除 LR/DT/RF/XGBoost (A4 in-flight) 外 6 種其他算法的合理性 + 預期 LOSO 結果 + 補測投票
parent: research/methyl_augmented_filter_phase2/cycle4/loso_validation/loso_findings.md
predecessor: research/methyl_augmented_filter_phase2/cycle3/cycle3_caller_f1_gate.json
verdict: 6 算法均不預期翻轉 LOSO sample-level 0 ΔF1; top vote = Per-zone LR (3/5) + Threshold-only (3/5); 其餘 1-2/5
last_verified: 2026-05-21
-->

# A4-Ext: 其他算法方法 F1 提升合理性盤點

> **任務定位**: PI 4-goal 報告目標 4 due-diligence companion。盤點不是執行 — 只給 6 種「除 LR/DT/RF/XGBoost 外」算法的定義、預期 LOSO ΔF1、補測投票。**全部 6 種預期都不能翻轉 LOSO sample-level ≈ 0 的核心 finding**，但盤點本身為「為何不試 X」提供實證 anchor。

---

## §1 背景與盤點動機

### 1.1 為何要盤點

用戶 2026-05-21 提問：「目標 4 confirm 還有沒有其他也合理適合提升 F1 的算法方法，幫我確認是否合理可用於驗證」。

A4 (cycle 4 Trial B) 已在跑 RF / XGBoost / GBM 三個非線性 model。但「ML 演算法宇宙」遠不止 LR/DT/RF/XGB — PI 報告必須對「為何不試 X」給實證 anchor，否則被質疑「沒窮舉」時無依據。

### 1.2 三條 anchor 約束所有候選

| Anchor | 來源 | 量化 |
|---|---|---|
| **A1 LOSO 已證 LR 5 樣本全 ≈ 0 ΔF1** | `cycle4/loso_validation/loso_findings.md` | mean ΔF1 = -0.00004, best τ → 0.10 全退化 |
| **A2 caller_af direction inconsistent** | `loso_cv_results.tsv` coefs_json | HCC1395 -1.20, HCC1937 +1.14, HCC1954 +0.13, H1437 -0.09, H2009 -0.15 (5 樣本反向) |
| **A3 F1 ceiling-bounded 4 樣本** | cycle 2 + LOSO | caller F1 ≥ 0.7: HCC1395 0.72, HCC1954 0.84, H1437 0.87, H2009 0.89; 只剩 HCC1937 (0.37) 有 headroom |

**A1+A2+A3 推論**: 任何 ML 算法「在 5 features 框架」+「在 4 個 high-F1 樣本上」**不可能** improvement > +0.005 (Cohen ribbon)。剩 HCC1937 唯一 headroom 但 single-sample n=1 無法升 tier。

### 1.3 盤點 6 算法清單

A4 已涵蓋 (不盤點): LR baseline / Decision Tree / Random Forest / XGBoost

本文盤點 (A4 不涵蓋):
1. Boolean rule / zone-specific gate
2. Threshold-only optimization (per-feature)
3. Ensemble (LR + RF voting)
4. Calibration (Platt / isotonic)
5. Per-zone LR (LOH-in vs LOH-out)
6. Per-AF-bin LR (low/mid/high caller_af)

---

## §2 6 個算法 Panel

### Panel 1: Boolean rule / zone-specific gate

| Dim | Detail |
|---|---|
| **定義** | 完全不依賴 ML，per-LOH-zone 寫死 hard rule。Example: `keep if (loh_inner_flag=1) OR (caller_af > 0.30 AND V6_off_NG ≥ 2)`。Hyperparam: 無（rule 本身就是參數）。適合場景：當 ML threshold 不穩 + decision boundary 由領域知識主導 |
| **預期 LOSO ΔF1 範圍** | **-0.005 ~ +0.003** (4 high-F1 樣本); **+0.005 ~ +0.020 single-sample (HCC1937 only, n=1)** |
| **為何合理試** | (1) cycle 3 caller-F1 gate 已驗證 "qualifying subset mean +0.01499" — 已是 partial Boolean rule (per-sample circular caveat); (2) Boolean rule 不受 caller_af direction inconsistent 影響（不擬合 coef）; (3) 業界 reference: VEP `IMPACT=HIGH` 類 hard-coded variant filtering 標準作法 |
| **風險** | (a) Hard rule 在 4 high-F1 樣本同樣受 F1-ceiling 限制；(b) per-LOH-zone × per-feature × per-threshold 組合爆炸（5^3 = 125 combinations），手動 sweep 易 multiple-testing inflate FP; (c) cycle 3 step 1 已試類似邏輯 mean +0.01499 仍 sample-level circular（per-sample re-fit） |
| **vote** | **2 / 5** (低 ROI) — cycle 3 已試類似邏輯結論 marginal；同 anchor A3 ceiling 限制 |
| **推薦 followup** | 暫緩。若 HCC1937 single-sample case study 仍想追，建議直接寫 paper case-study section 而非 cycle 5 全跑 |

---

### Panel 2: Threshold-only optimization (per-feature)

| Dim | Detail |
|---|---|
| **定義** | 對單一 feature 直接 LOSO τ-sweep，不訓練 model。Example: `caller_af > τ*` 找最佳 τ on 4 train samples, apply 在 1 test sample。Hyperparam: 1 個 τ per feature。適合場景: 當 single feature 已知是 dominant predictor + ML 反而引入 noise |
| **預期 LOSO ΔF1 範圍** | **-0.002 ~ +0.005** (best feature = caller_af, 但 LOSO best τ → grid 端點 = keep everything; H_new4 drop caller_af 已實證 HCC1395 +0.00699 / HCC1937 +0.00000) |
| **為何合理試** | (1) LR 在 LOSO 退化 best τ → 0.10 = keep everything → 暗示問題出在「需要 cross-sample transferable τ」不是「需要更複雜 model」; (2) Single-feature τ 對 caller_af direction flip 仍 robust (因為 unidirectional sweep); (3) 業界 reference: bcftools filter `INFO/AF > 0.05`、GATK VariantFiltration 都是 single-threshold rule |
| **風險** | (a) 5 features 各跑 → multiple-testing inflate；(b) cycle 1+2+3 已隱含驗證: cycle 1 caller_af LR coef +3.44 dominant 但 cycle 2 cross-sample 4/5 NEGATIVE → 單 feature τ 同樣 sample-specific; (c) cycle 4 LOSO H_new4 (drop caller_af) 已測得 mean +0.001 / Wilcoxon p=0.31 → threshold-only on remaining 4 features 不會更高 |
| **vote** | **3 / 5** (值得試 — 但有 conditional) — 已有 H_new4 部分證據; conditional = "只跑 caller_af + V6_off_NG 兩個最強 feature, ~1 hr, 用作 'ML negative 但 hard threshold 也 negative' 的 robustness check" |
| **推薦 followup** | **Phase A4 內補跑 quick τ-sweep (2 features × 5 LOSO folds = 10 runs)**, 純 robustness check 用; 不期待翻轉但補上「窮舉 evidence」 |

---

### Panel 3: Ensemble (LR + RF voting)

| Dim | Detail |
|---|---|
| **定義** | sklearn `VotingClassifier(voting='soft')` 平均 LR 與 RF 的 predict_proba。Hyperparam: weights (default 1:1) + 個別 model hyperparam。適合場景: 當 LR 與 RF 互補 (linear + nonlinear) 且兩者各自 marginal POSITIVE |
| **預期 LOSO ΔF1 範圍** | **-0.001 ~ +0.002** — LR LOSO 0, RF LOSO 預期也 ~0 (A4 正在跑驗證), ensemble 不創造新 signal |
| **為何合理試** | (1) 業界 trick: Kaggle 常用 stacking/voting 提升 0.001-0.01；(2) 在 calibration 不一致時，ensemble averaging 可降 prediction variance；(3) sklearn 一行 import 成本極低 |
| **風險** | (a) 數學上, 兩個各自 ≈ 0 ΔF1 的 model voting 仍 ≈ 0 (期望值線性);  (b) ensemble 無法克服「sample-level circularity」這個根本問題 — root cause 在 distribution shift 不在 model bias; (c) Kaggle 0.001-0.01 提升假設前提是「兩個 base model 已 > random + 互補」, LOSO 已證 LR base ≈ 0 → ensemble lower bound ≈ 0 |
| **vote** | **1 / 5** (不值得) — 理論上知道失敗 (兩個 0 + 0 = 0); 唯一例外: A4 RF 跑出 POSITIVE 才 worth 加 LR 平均 |
| **推薦 followup** | 條件式 — 等 A4 RF/XGB 結果, 若 RF/XGB LOSO ΔF1 > +0.005 (任一 sample) 再考慮 ensemble 是否能 stabilize |

---

### Panel 4: Calibration (Platt / isotonic)

| Dim | Detail |
|---|---|
| **定義** | sklearn `CalibratedClassifierCV(method='sigmoid' 或 'isotonic')`. 修正 p_pred 讓 p=0.5 真實對應 50% TP 機率。Hyperparam: method (sigmoid=Platt / isotonic) + cv folds。適合場景: 當 model AUC 良好但 threshold 不易選 (rank correct, scale wrong) |
| **預期 LOSO ΔF1 範圍** | **-0.001 ~ +0.001** — 純機率重映射, F1 在 monotone-aware τ-sweep 下不變 (re-mapping 不改 ranking) |
| **為何合理試** | (1) scientific-rigor §3 推薦對 PROBABILITY claim 校準; (2) 若 PI 想用 p_pred 報告 (e.g. "high-confidence FP > 0.9")，calibration 是 prerequisite; (3) Platt scaling 在 LR 本身已內建（sigmoid output），所以 isotonic 才有意義 — isotonic 對 non-linear miscalibration robust |
| **風險** | (a) F1 metric 本質上是 ranking-based + threshold-sweep，calibration 在 sweep 下數學上 invariant (任何 monotone transformation g(p) sweep over g(τ) 等價於 sweep over τ); (b) 若 best τ 在 grid 端點 (= 0.10 / 0.95), calibration 後 best τ 也 shift 到對應端點，ΔF1 不變 |
| **vote** | **1 / 5** (不值得 for F1 metric) — 數學上知道 invariant; **但** 若 PI 報告需要 p_pred quote → 升 3/5 |
| **推薦 followup** | 暫不為 F1 跑; 若 paper §3 要 quote "high-confidence FP threshold p > 0.9" → 補做 1 個 isotonic for documentation |

---

### Panel 5: Per-zone LR (LOH-in vs LOH-out)

| Dim | Detail |
|---|---|
| **定義** | 把 5 features 拆成 2 個 LR — `LR_LOH_in` (train on `loh_inner_flag=1` 子集) + `LR_LOH_out` (train on `loh_inner_flag=0` 子集), 各自找 best τ。Hyperparam: 同 cycle 1 LR (C=1.0 L2). 適合場景: 當 LOH zone 內外 feature distribution 顯著不同 (interaction strong) |
| **預期 LOSO ΔF1 範圍** | **-0.003 ~ +0.008** — 在 HCC1937 (low caller F1, LOH heterogeneous) 可能 +0.005~+0.015, 其餘 4 高 F1 樣本仍 ≈ 0 |
| **為何合理試** | (1) cycle 1 LR `loh_inner_flag` coef +1.95 = 最強單 feature (僅次 caller_af) → LOH 是 known interaction axis; (2) cycle 3 caller-F1 gate logic 內已 implicit per-zone — formalize 成 per-zone LR 是 natural extension; (3) 業界 reference: Mixture-of-experts / stratified modeling 在 heterogeneous population 常見; (4) **A4 已涵蓋 'interaction terms in single LR'** (Trial A); per-zone = explicit hard split, 數學上比 interaction 更 expressive (allow 完全不同 coef sign) |
| **風險** | (a) 樣本量切半 → 每個 LR variance ↑ (HCC1937 LOH_in subset 可能 < 5k rows 不夠 train); (b) per-zone 仍受 sample-level circularity 限制 (LOSO 4-sample-trained 對 1 held-out 仍 distribution shift); (c) cycle 3 step 1 partial gate 已實證 per-zone-like rule 在 LOSO 後仍 ≈ 0 |
| **vote** | **3 / 5** (值得試 — 有 conditional 機會) — A4 interaction 已 cover 一半, 但 explicit per-zone 在 HCC1937 single-sample 仍可能有 0.005-0.015 visible effect; conditional = "primarily as HCC1937 case study supplement" |
| **推薦 followup** | **建議 cycle 5 補測 1 hr quick test**: 只跑 LOSO HCC1937-held-out (因為其他 4 樣本 F1-ceiling-bounded 不會有 signal); 結果用作 single-sample case-study evidence (不升 tier) |

---

### Panel 6: Per-AF-bin LR (low/mid/high caller_af)

| Dim | Detail |
|---|---|
| **定義** | 把 caller_af 切 3 bins (例 [0,0.1) / [0.1,0.3) / [0.3,1.0]), 對每個 bin 訓 1 個 LR。Hyperparam: bin 邊界 + 同 cycle 1 LR C=1.0 L2. 適合場景: 當 caller_af 與其他 feature 有 strong interaction (例如 low-AF read-level signal 與 high-AF callset-level signal 完全不同 mechanism) |
| **預期 LOSO ΔF1 範圍** | **-0.005 ~ +0.005** — caller_af direction inconsistent (A2 anchor) 強烈暗示「bin 內也不一致」, 反而 bin 切割可能放大 variance |
| **為何合理試** | (1) caller_af direction inconsistent 是 Phase 2 核心未解 confound; (2) 業界假說: low-AF (clonal-emerging) variants 由 read-level methylation 主導，high-AF (heterozygous-like) 由 mappability 主導 — 若真，per-bin LR 可分離兩 regime; (3) cycle 3 Step 1.5 ablation 顯示 HPFineF coef = caller_af L2 ridge-split (not independent) — per-AF-bin 可解耦這個 collinearity |
| **風險** | (a) **數學上**: 若 direction flip 來自 sample-level distribution 不是 AF-level interaction, 切 AF bin 無法解; (b) bin 邊界 sensitivity (0.1 vs 0.15 vs 0.2 可能完全不同結果) → 引入 researcher DoF (Pre-reg 風險); (c) 樣本量切 3 倍稀疏, HCC1937 low-bin 可能 < 2k rows; (d) **L2 collider bias 警告 (MEMORY.md)**: caller_af 是 outcome causal pathway 上的 collider — 切 bin 後 conditioning on caller_af 可能引入虛假 association |
| **vote** | **2 / 5** (低 ROI) — 風險 d (collider bias) 是 deal-breaker; 若想驗證「low-AF read-level signal」更好用 read-level model (cycle 5 paper §5 候選) 而非 per-bin LR |
| **推薦 followup** | 暫緩。若 paper §5 想討論 low-AF regime, 用 ReadParser read-level 框架而非 per-bin variant-level LR |

---

## §3 補測投票排序 (highest first)

| Rank | Algorithm | Vote | Recommended action | Expected effort | Expected ΔF1 (LOSO mean) |
|---:|---|:--:|---|---|---|
| 1 | **Threshold-only (caller_af + V6_off_NG only)** | **3/5** | Quick LOSO τ-sweep 補入 A4 robustness check | 1 hr | -0.001 ~ +0.002 (mean ≈ 0) |
| 1 | **Per-zone LR (LOH-in / LOH-out)** | **3/5** | HCC1937 single-sample LOSO test for case-study supplement | 1.5 hr | HCC1937 +0.005~+0.015 / 4 others ≈ 0 |
| 3 | Boolean rule / zone-gate | 2/5 | 暫緩 (cycle 3 已 partial cover) | (skip) | -0.005 ~ +0.003 |
| 3 | Per-AF-bin LR | 2/5 | 暫緩 (collider bias 風險) | (skip) | -0.005 ~ +0.005 |
| 5 | Ensemble (LR + RF voting) | 1/5 | Conditional: only if A4 RF LOSO > +0.005 | (gated) | -0.001 ~ +0.002 |
| 5 | Calibration (Platt / isotonic) | 1/5 | F1-invariant 不為 F1 跑; 若 PI 要 p_pred quote 補做 isotonic | (optional) | ~0 (math invariant) |

### Top 1 推薦 (二者並列)

**Top 1a: Threshold-only τ-sweep on (caller_af + V6_off_NG)** — vote 3/5, effort 1 hr, 為 PI 報告補上「窮舉 evidence」: 「LR LOSO ≈ 0 + threshold-only LOSO ≈ 0 → 兩條獨立路徑同向, 不是 model 複雜度問題」。
**Top 1b: Per-zone LR HCC1937-held-out** — vote 3/5, effort 1.5 hr, 為 paper §3 補 HCC1937 single-sample case study (不升 tier, 僅 anecdotal evidence on caller-F1-headroom mechanism)。

兩者合計 < 3 hr, 推薦 cycle 4 完成後 1-day 補測。

---

## §4 PI 報告中「為何不試 X」實證對應段落

> 以下段落供 PI 報告 §4.3「Algorithm-space due diligence」直接引用。

### 4.4.1 為何不試 Boolean rule / zone-gate

> 「Cycle 3 caller-F1-headroom gate (qualifying subset mean +0.01499) 已為 Boolean rule 的代表; LOSO 後因 per-sample circular re-fit 降至 ⭐1 L4。Boolean rule 整體 vote 2/5 — 與 cycle 3 已驗結論同源，補測 ROI 低。」

### 4.4.2 為何不試 Threshold-only optimization

> 「Cycle 4 LOSO 已隱含驗證: 5 樣本 best τ 全退化為 0.10 (grid 端點 = keep everything) → 即使 single-feature τ optimization 也找不到 transferable threshold。H_new4 (drop caller_af) sensitivity 補測顯示 4-feature τ-sweep 對 HCC1395 +0.00699 / HCC1937 +0.00000，與 LR 同向 ≈ 0。Vote 3/5 推薦 1 hr quick robustness check 補入完整 evidence。」

### 4.4.3 為何不試 Ensemble (LR + RF voting)

> 「Ensemble 在 base model LOSO ≈ 0 時數學上 lower bound ≈ 0 (期望值線性)。A4 RF/XGBoost 正驗證非線性 base model；只有當 RF LOSO ΔF1 > +0.005 才考慮 ensemble。Vote 1/5 conditional。」

### 4.4.4 為何不試 Calibration (Platt / isotonic)

> 「F1 metric 在 monotone re-mapping 下 invariant — calibration 不影響 ranking, 對 F1 數學上 0 改善。Vote 1/5 for F1. 若 PI 報告要 quote 'high-confidence FP threshold p > 0.9' → 補做 1 個 isotonic for documentation (3/5 in that context)。」

### 4.4.5 為何不試 Per-zone LR

> 「A4 Trial A interaction terms 已部分涵蓋 per-zone effect (math: interaction = soft per-zone); explicit per-zone LR 在 HCC1937 single-sample 仍有 +0.005~+0.015 潛在 visible effect (vote 3/5)，但其他 4 高 F1 樣本 F1-ceiling-bounded 不會貢獻。建議 cycle 5 補 HCC1937 LOSO 作 case-study supplement (不升 tier)。」

### 4.4.6 為何不試 Per-AF-bin LR

> 「caller_af 是 outcome 因果路徑上的 collider (MEMORY.md L2 collider bias warning); 切 AF bin 後 conditioning on caller_af 可能引入虛假 association。若想討論 low-AF regime 應用 read-level 模型 (cycle 5 paper §5 候選), 不應用 per-AF-bin variant-level LR。Vote 2/5, 暫緩。」

### 4.4.7 整體 sanity check

> 「6 個非 A4 算法 vote 平均 = 12/30 = **40th percentile**, 對應 'low expected ROI / sample-level signal absent' 結論。**沒有一個算法理論上能翻轉 LOSO sample-level ≈ 0 finding**，因 root cause 是 distribution shift (sample-level circularity) 不是 model bias。Phase 2 真正 pivot direction 在 read-level 框架 (cycle 5+) 或 phasing-signature track (LOH-constrained phasing discovery, paper main thesis 候選), 不在 algorithm-space 內。」

---

## §5 Files

```
docs/experiments/in_progress/2026/05/20260521_A4_Ext_other_algorithms_inventory_01.md   ← this report
research/methyl_augmented_filter_phase2/phase2_completeness_audit/A4_Ext_algorithm_inventory.tsv
research/methyl_augmented_filter_phase2/phase2_completeness_audit/figures/A4_Ext_algorithm_vote_ranking.png
```

---

## §6 Caveats & Limits of this inventory

1. **不執行**: 本文 vote 基於 mechanism reasoning + cycle 1-4 既有實證, 非真跑;
2. **6 算法非窮舉**: 仍有 SVM / k-NN / Neural Net / Bayesian filter 等; 但這 4 個與 RF/XGB 同樣 fundamental limited by sample-level circularity, 不另列;
3. **LOSO sample-level circularity 是 root cause**: 任何 variant-level model 都受限; 真正 pivot 需 read-level 或 phasing track;
4. **F1 ceiling-bounded 4 樣本**: caller F1 ≥ 0.7 樣本 (HCC1395 / HCC1954 / H1437 / H2009) 即使 model 完美, F1 上限受 caller recall 與 TP/FP imbalance 限制; 任何 +0.005 以上 effect 都需要在 ceiling-bounded regime 內找 — 結構性困難;
5. **HCC1937 唯一 headroom 是 n=1**: 5 樣本中 single-sample evidence 升 tier 受 multi-sample-consistency skill 限制 (⭐3+ 需 4+ sample 一致), HCC1937 alone 永遠 ≤ ⭐2 L4。

---

**End of A4-Ext inventory — supplements PI Goal 4 due-diligence**
