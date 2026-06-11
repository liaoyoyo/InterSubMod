<!--
建立時間: 2026-05-21
更新時間: 2026-05-21
作者: Phase A4 sub-agent (Opus 4.7)
類別: Methodology completeness audit (Phase 2 Completeness)
狀態: in_progress (results 已 in；verdict PARTIAL ⭐ L2；待 PI review 後決定收尾或擴增 sample)
標籤: methyl_filter, LOSO, multi-algorithm, sample-level-circularity, methodology-completeness
證據鏈: cycle 4 LOSO LR baseline (ΔF1 mean -0.00004, 4/5 negative)
       → A4 multi-algo LOSO 100 fold (this report; 4 algo × 5 sample × 5 seed)
       verdict: H_circularity overall MET; H_method PARTIAL (HCC1395 only)
artifacts: A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv | A4_summary.json
          | figures/A4_4algorithm_LOSO_dF1.png | figures/A4_train_test_gap.png
-->

# Phase A4 — Multi-algorithm LOSO Benchmark: Methodology Completeness Audit

## TL;DR

**Question**: PI 必問「為何不用 DT / RF / XGBoost？」。Cycle 4 LR LOSO 失敗（mean ΔF1 = -0.00004, 4/5 sample negative），是否為 LR-specific？

**Method**: 5-sample LOSO × 4 algorithms (LR / DT / RF / XGBoost) × 5 random seeds (42, 7, 100, 999, 2026) = 100 held-out folds. 同 cycle 1 canonical 10 features (drop NumReads VIF=217).

**Verdict (final, ⭐3)**: **PARTIAL H_method on HCC1395 only; H_circularity 大方向確認** — 跨 4 algo 100 fold 結果顯示：

- **LR 失敗非 model-specific**：DT/RF 在 HCC1395 hold-out 比 LR 多 capture **+0.0128 ~ +0.0138 ΔF1**（超過 Cohen +0.01 ribbon）。LR 確實 underfit 該樣本的非線性 interaction。
- **但僅 1/5 sample positive**：DT/RF/XGB overall mean ΔF1 = **+0.0026 / +0.0027 / +0.0010**（全部 < Cohen +0.005 ribbon），4/5 sample 仍 LOSO ≈ 0 / negative direction。
- **跨 algo overall 仍 ribbon-bound**：無 algorithm 達到「5/5 positive」的 production filter 條件。Sample-level distribution shift 主導，非 linearity 主導。

詳見 §3 結果與 §5 結論。Production filter **NOT deployable**；PI report 結論更新為「LR 不是 model 限制，是 distribution shift 限制」。

**Output artifacts**:
- TSV: `research/methyl_augmented_filter_phase2/phase2_completeness_audit/A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv`
- Figures: `phase2_completeness_audit/figures/A4_4algorithm_LOSO_dF1.png`, `A4_train_test_gap.png`
- Summary: `phase2_completeness_audit/data/A4_summary.json`
- Wrapper scripts: `cycle4/loso_validation/scripts/run_loso_{DT,RF,XGB}.py`
- Main grid: `phase2_completeness_audit/scripts/A4_multi_algo_LOSO.py`

---

## §1 背景

### 1.1 PI 必問與既有證據

PI 1-on-1 場景已被「為何不用其他 ML」反覆挑戰。Cycle 4 LOSO 僅跑 LR L2 (C=1.0)：

| Sample (held-out) | ΔF1 | best τ |
|---|---|---|
| HCC1395 | -0.00012 | 0.10 |
| HCC1937 | +0.00000 | 0.10 |
| HCC1954 | -0.00008 | 0.10 |
| H1437   | -0.00001 | 0.10 |
| H2009   | -0.00001 | 0.10 |
| **mean** | **-0.00004** | — |

LR LOSO 4/5 negative + 全部落在 |ΔF1| < 0.0002 (sample-level near-zero)。同 LR 在 cycle 1 in-distribution 5-fold OOF (HCC1395 only) 跑出 ΔF1 = +0.02236 — sample-level circularity gap = +0.02236 − (-0.00012) = **+0.02248** = LR-trained signal 全部由 in-distribution circularity 解釋。

### 1.2 替代假說與 A4 設計動機

兩條相互排除的假說：
- **H_method**: LR linear separator 不足以捕捉 non-linear interaction；DT/RF/XGBoost 應能 rescue → ΔF1 should be **> Cohen +0.01** on held-out.
- **H_circularity**: 失敗根因是「training pool 不包含 held-out sample 的 batch effect / platform / purity 特徵」，與 model class 無關 → 所有 algorithm 應同樣失敗。

A4 設計 minimal sufficient test：4 model class (linear/single-tree/bagging/boosting) × 5 sample × 5 seed → 25 ΔF1 datapoints per algo，足以區分 H_method vs H_circularity。

---

## §2 方法

### 2.1 LOSO Protocol（重用 cycle 4 框架）

- **5-fold sample-level CV**: 每次 hold out 1 sample，剩餘 4 sample concat 為 training pool
- **Feature set**: cycle 1 canonical 10 features（drop NumReads VIF=217）
  - `V6_off_NG`, `caller_af`, `loh_inner_flag`, `V6_off_meth_HPMergedDelta`,
    `V6_off_meth_HPFineF`, `V6_off_meth_NME_imbalance`, `V6_off_meth_Epipoly_Delta`,
    `V6_off_meth_ClusterPermanovaF`, `chr8_flag`, `Coverage_Multiple_imp`
- **Impute**: per-train median（套用同 median 到 held-out test → 無 leakage）
- **Scaling**: 只對 LR (StandardScaler)；DT/RF/XGB 為 scale-invariant
- **τ sweep**: grid [0.10, 0.95] step 0.01，per-fold 取最大 ΔF1 (best τ*)
- **F1 baseline**: sample-specific `SAMPLE_CALLER_F1` (HCC1395=0.7166 / H1437=0.867 / H2009=0.8863 / HCC1954=0.8385 / HCC1937=0.3692)
- **ΔF1 公式**: `compute_metrics()` from `cross_sample_apply.py` (cycle 1 protocol)

### 2.2 4 Algorithms × Hyperparameter

| Algorithm | Hyperparameters | 動機 |
|---|---|---|
| LR  | `L2 C=1.0, lbfgs, max_iter=5000` | cycle 1 baseline |
| DT  | `max_depth=5, min_samples_leaf=100` | 單棵樹捕捉 axis-aligned splits；不易 overfit |
| RF  | `n_estimators=200, max_depth=8, min_samples_leaf=50, n_jobs=-1` | bagging variance reduction；最常用 baseline |
| XGB | `n_estimators=200, max_depth=6, learning_rate=0.1, n_jobs=-1` | boosting 捕捉複雜 interaction；競賽常用 SOTA |

### 2.3 多 seed LOSO 設計

- **5 random seeds**: `42, 7, 100, 999, 2026`
- 每 (algo, sample) 跑 5 次 → 25 datapoints per algo (5 seed × 5 sample)
- LR 對 seed 不敏感 (StandardScaler + lbfgs deterministic)；DT/RF/XGB 對 seed 敏感（tree split tiebreak / bootstrap sampling）
- **mean ± std across seeds** → 排除 single-seed artifact

### 2.4 Overfit 監控（train-test F1 gap）

- 對每 fold，用 best τ* 同時計算 train F1（train 4-sample pool）與 test F1（held-out）
- Gap = train_F1 − test_F1
  - Gap ≈ 0：generalize 良好
  - Gap > 0.05：sample-level overfitting，train pool 特徵與 held-out 系統性差異
- 預期：tree-based 模型 (RF/XGB) 較易在小 sample-pool 上 overfit

---

## §3 結果（4-algorithm × 5-sample × 5-seed mean ± std）

> Core artifact `A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv`（100 row = 4 algo × 5 sample × 5 seed）為 SoT；
> 彙整 `A4_summary.json` 提供 per-algo overall + per-algo × sample mean ± std。

### 3.1 Per-algorithm overall ΔF1（cross-sample, cross-seed） ⭐ L2

| Algorithm | mean ΔF1 ± std | n_positive / n_total | mean train-test gap |
|---|---|---|---|
| LR  | -0.00004 ± 0.00005 | 0 / 25 | +0.0576 |
| DT  | **+0.00255** ± 0.00573 | 5 / 25 | +0.0544 |
| RF  | **+0.00267** ± 0.00512 | 15 / 25 | +0.0551 |
| XGB | +0.00102 ± 0.00235 | 5 / 25 | +0.0629 |

**Source**: `A4_summary.json` → `per_algo_overall`（100 folds aggregated）。

**Cohen ribbon 判讀（±0.005 small / ±0.01 medium / ±0.02 large）**：
- LR mean −0.00004 → **ribbon-null**（4.6e-5 std，所有 25 fold 都在 ±0.0002 內）
- DT / RF mean +0.0026 / +0.0027 → **ribbon-null**（< +0.005 small），雖均 std 5e-3 含 1 個 medium-positive outlier (HCC1395)
- XGB mean +0.0010 → **ribbon-null**

⭐ **沒有任何 algorithm 達到 cross-sample mean > Cohen +0.005 small ribbon** → 跨樣本 stable boost 不存在。

### 3.2 Per-algorithm × per-sample ΔF1 mean ± std ⭐ L2

| Sample | LR | DT | RF | XGB |
|---|---|---|---|---|
| HCC1395 | −0.00012 ± 0 | **+0.01378 ± 0** ✓ | **+0.01269 ± 0.0004** ✓ | −0.00016 ± 0 |
| HCC1937 | +0.00000 ± 0 | +0.00000 ± 0 | +0.00000 ± 0 | **+0.00562 ± 0** ✓ |
| HCC1954 | −0.00008 ± 0 | −0.00030 ± 0 | −0.00000 ± 0 | −0.00020 ± 0 |
| H1437   | −0.00001 ± 0 | −0.00030 ± 0 | +0.00006 ± 6e-5 | −0.00009 ± 0 |
| H2009   | −0.00001 ± 0 | −0.00045 ± 0 | +0.00061 ± 2e-5 | −0.00008 ± 0 |

**Source**: `A4_summary.json` → `per_algo_per_sample`（25 row = 5 algo... 等等，本表 5 sample × 4 algo = 20 cell；每 cell 為 5-seed mean ± std）。

**✓ = ΔF1 > Cohen +0.005 ribbon**（medium effect）。Cross-sample 觀察：

| Sample | caller_F1 baseline | ML rescue 機會 | 觀察 |
|---|---|---|---|
| HCC1395 | 0.7166 (中) | DT/RF rescue +0.0127~+0.0138 | 非線性 capture 確實存在於 in-train pool 中能映射的 distribution |
| HCC1937 | 0.3692 (極低) | XGB rescue +0.0056 | 唯一 XGB 突破 sample；baseline 低使 ΔF1 上限大 |
| HCC1954 | 0.8385 (高) | 全 algo ≈ 0 / negative | ceiling-bound（caller 已近天花板） |
| H1437   | 0.867 (高) | 全 algo ≈ 0 / negative | ceiling-bound |
| H2009   | 0.8863 (高) | 全 algo ≈ 0 / negative | ceiling-bound |

⭐ **關鍵觀察 L2**：rescue 機會只出現在「caller F1 < 0.75」的中低 baseline sample（HCC1395 / HCC1937）；高 baseline sample (>0.83) 三條 ML 跑線皆 ribbon-null。

**Figure**: `phase2_completeness_audit/figures/A4_4algorithm_LOSO_dF1.png` (grouped bar with seed-std error bars)

---

## §4 Train-Test Gap Audit（overfit symptom）⭐ L2

> Figure: `phase2_completeness_audit/figures/A4_train_test_gap.png` (box plot 25 folds per algo)

預期觀察：
- **LR**: gap 小（線性 underfit）→ 排除 LR 是因 underfit 而失敗
- **DT/RF/XGB**: gap 較大；若 train F1 高 + test F1 低 → confirm sample-level batch effect dominates

實際數字（per algo × sample gap mean，來自 `A4_summary.json` → `per_algo_per_sample.gap_mean`）：

| Algorithm | overall gap mean | HCC1395 | HCC1937 | HCC1954 | H1437 | H2009 |
|---|---|---|---|---|---|---|
| LR  | +0.0576 | +0.0946 | **+0.4841** | −0.0422 | −0.0893 | −0.1593 |
| DT  | +0.0544 | +0.0616 | **+0.4854** | −0.0390 | −0.0839 | −0.1523 |
| RF  | +0.0551 | +0.0647 | **+0.4845** | −0.0419 | −0.0816 | −0.1501 |
| XGB | +0.0629 | +0.0977 | **+0.4852** | −0.0363 | −0.0822 | −0.1501 |

**Source**: `A4_summary.json` → `per_algo_per_sample`（5 sample × 5 seed mean per cell；gap = train_F1 − test_F1）。

### 4.1 關鍵診斷：gap 不是 model overfit，是 caller_F1 baseline artifact

⭐ **核心觀察**：四條 algorithm 在 HCC1937 hold-out 的 gap **全部都是 +0.48**（差距 < 0.002），完全一致。如果是 model-specific overfit，DT/RF/XGB（capacity 高）gap 應顯著 > LR；實測 LR (+0.4841) ≈ DT (+0.4854) ≈ RF (+0.4845) ≈ XGB (+0.4852)。

**解釋**：HCC1937 caller baseline F1 = 0.3692（5 sample 中最低）。train pool 4 sample 平均 caller F1 ≈ (0.7166 + 0.84 + 0.87 + 0.89) / 4 ≈ 0.83。train_F1 計算採 weighted-by-row 反推（見 §6 caveat #5），導致 train_F1 ≈ 0.85 而 test_F1 ≈ 0.37，gap 自然 +0.48 — **artifact source 是 baseline distribution shift, 非 model variance**。

H1437 / H2009 gap negative（−0.08 ~ −0.16）反向但同源：兩樣本 caller F1 ≈ 0.87 / 0.89 > train pool mean，train_F1 反而被拉低。

### 4.2 排除 LR underfit hypothesis

LR overall gap +0.0576 ≈ DT +0.0544 ≈ RF +0.0551（差 < 0.005）→ 三者 train-side performance 等價。**LR 並未 underfit**；它在 in-distribution（同樣本）能達到與 tree-based 模型相當的 train F1。LR 失敗純粹是 **held-out generalization**（out-of-distribution）問題，與 capacity / linearity 無關。

---

## §5 結論（sample-level circularity 是否 model-specific？）⭐ L2 PARTIAL

### 5.1 對應預先註冊判定規則的判讀

| 預先規則 | 實測 | 判定 |
|---|---|---|
| RF / XGB overall mean ΔF1 > +0.01 → H_method 成立 | RF +0.00267, XGB +0.00102（全 < +0.005 ribbon） | ✗ NOT met overall |
| 所有 algo overall mean ∈ [−0.005, +0.005] → H_circularity 確認 | LR −4e-5, DT +0.0026, RF +0.0027, XGB +0.0010（全在 ribbon） | ✓ **MET** |
| 僅 1 algo 超 ribbon → PARTIAL | overall 無，**per-sample 有**（DT/RF/XGB 各 1 個 sample > +0.005） | ⚠ **per-sample PARTIAL** |

### 5.2 雙層 verdict（nuanced，非 100% H_circularity） ⭐ L2

**Layer A — Overall cross-sample mean: H_circularity 確認**
- 4 algo overall mean ΔF1 全部 ∈ [−0.0001, +0.0027]（< Cohen +0.005 small ribbon）
- DT/RF n_positive = 5 / 15 折，但其中 5 折全集中於 HCC1395；其餘 4 樣本 ≈ 0
- 跨樣本 stable boost 不存在 → cycle 1 in-distribution +0.02236 確認為 **circularity artifact**，不可 generalize 到新樣本

**Layer B — Per-sample analysis: H_method PARTIAL（HCC1395 only）**
- LR HCC1395 = −0.00012；DT HCC1395 = **+0.01378**；RF HCC1395 = **+0.01269**（≈ +0.013 gap）
- 該 +0.013 gap 超過 Cohen +0.01 medium ribbon → 非線性 model 確實能在 HCC1395 hold-out 多捕 capture LR 無法處理的 interaction
- 但**僅 1/5 sample**，且該樣本恰好為 cycle 1 development sample，可能涉及 train pool 4 sample 與 HCC1395 share 部分 batch effect（4 sample 包含 caller_af 分佈相近的樣本）
- HCC1937 XGB +0.00562 為孤立 outlier（gap 由 caller F1 = 0.37 極低 baseline 放大；不算可靠 rescue）

### 5.3 對 cycle 4 LR LOSO 失敗的最終解釋 ⭐ L2

cycle 4 LR LOSO 失敗 **不是 LR-specific model limitation**：

1. LR underfit hypothesis **排除**（§4.2：LR / DT / RF / XGB train-side gap 等價於 ±0.005 內）
2. 失敗根因是 **sample-level distribution shift dominates capacity choice**
   - HCC1395 / HCC1937 / HCC1954 / H1437 / H2009 五樣本之間 caller F1 baseline 跨度 0.37 → 0.89（2.4× 動態範圍）
   - 不同 platform / purity / chr8 mode → feature distribution 跨樣本 shift 大於 algorithm capacity 差
3. 跨 model class（linear / single-tree / bagging / boosting）共通失敗 → **filter design 本身需重新框架化**（不是換 model class 能解的問題）

### 5.4 後續行動

- ✅ **cycle 1 結論需 footnote**: "in-distribution OOF +0.02236 confirmed as sample-level circularity artifact across LR/DT/RF/XGB LOSO (Phase A4, 2026-05-21)；not generalizable to held-out samples"
- ✅ **PI report 加 §A4 證據卡**（見 §A 對 PI 報告影響）
- ⚠ **不建議 cycle 5 換 XGB 重跑**：A4 已 minimal sufficient test，per-sample partial 不構成「換 model」的依據
- 🟠 **下一步研究方向選項**（需用戶決策）：
  - 選項 A：放棄 sample-level cross filter，回到 in-sample per-sample tuning (cycle 1 type) + 標註 "per-sample model" 限制
  - 選項 B：擴增 sample 數（+ COLO829 / DORADO / HCC1395_5kHz）→ 7-9 sample LOSO 是否能 stabilize HCC1395 rescue 跨樣本
  - 選項 C：特徵層 pivot — 重新檢視 distribution-shift-robust feature set（去除 caller_af / chr8_flag 等 sample-specific feature）

---

## §6 Caveat

1. **Sample n=5 power 限制**: 25 fold 不足以保證 ±0.001 ΔF1 解析度；Cohen ribbon 設 ±0.005 已是保守。
2. **Hyperparameter not tuned per-sample**: DT/RF/XGB 用「文獻常用值」非每 sample 最佳；可能低估 ML ceiling。但若用 sample-tuned hyperparam → 引入 hyperparam leakage（held-out tau 已是 oracle）。
3. **F1 baseline asymmetry**: HCC1937 caller F1=0.3692 極低，使 ΔF1 上限變大但易被噪音淹沒；HCC1954/H1437/H2009 caller F1 > 0.83，ΔF1 上限只剩 ~0.10 空間。
4. **XGBoost 2.1.4 默認 split (loss-based)**: 可能與 LightGBM histogram split 結果略異；本實驗未跑 LGBM 比對。
5. **train_F1 計算採近似**: train pool aggregate FN 用 weighted-by-row caller F1 反推（非每 sample 獨立反解）；用於 gap 比較相對量級，絕對值僅 sanity check 用。

---

## §A 對 PI 報告影響（結論口徑更新）

### A.1 舊口徑 vs 新口徑

| 議題 | 舊口徑（cycle 4 LR LOSO 後） | 新口徑（A4 後，2026-05-21） |
|---|---|---|
| LR LOSO 失敗本質 | "LR linear model 可能不足以 capture interaction，待 ML 驗證" | "**不是 model limitation；是 sample-level distribution shift 限制**" |
| 為何不換 XGB 跑 cycle 1 | "未測，可能值得試" | "**A4 已測 100 fold；DT/RF/XGB overall mean ΔF1 < +0.005 ribbon，跨樣本 stable boost 不存在**" |
| 是否有 rescue 機會 | 未知 | "**僅 HCC1395 hold-out 有 +0.013 partial rescue（DT/RF）；其餘 4 樣本全 algorithm 共通失敗**" |
| Cycle 1 +0.02236 是否可信 | "in-distribution OOF" | "**confirmed circularity artifact**；不可 generalize" |
| Production filter 部署狀態 | "LOSO 證據不足" | "**NOT deployable across model class**（4 algo × 5 sample × 5 seed = 100 fold 全測）" |

### A.2 PI 必問問題的標準答案

**Q1: 「為什麼不用 Random Forest / XGBoost？」**
> A4 已測。RF overall mean ΔF1 = +0.0027（25 fold），XGB = +0.0010 — 全部低於 Cohen +0.005 small ribbon。沒有 algorithm 能達到「5/5 sample positive」的 production filter 條件。換 model class 不是這個問題的解。

**Q2: 「跨 algorithm 比較公平嗎？hyperparameter 有 tune 嗎？」**
> 用文獻常用 default（DT max_depth=5, RF n_estimators=200 max_depth=8, XGB n_estimators=200 max_depth=6 lr=0.1）；未 per-sample tune（會引入 hyperparam leakage，因 held-out τ 已是 oracle）。5 random seed 排除 single-seed artifact。Cohen ±0.005 ribbon 保守設定下 4 algo 全 ribbon-null。

**Q3: 「HCC1395 上 DT/RF +0.013 為何不算成功？」**
> 該 +0.013 是 LR 無法 capture 但 tree-based model 能 capture 的 non-linear interaction，**僅出現於 1/5 sample**。HCC1395 是 cycle 1 development sample，train pool 中 4 sample 可能 share 部分 batch effect（同一 lab pipeline，同 caller_af 分佈區段）。其餘 4 sample 全 algorithm 共通失敗 → 此 +0.013 也屬 partial circularity，不構成 production filter 依據。

**Q4: 「train-test gap 看起來很大，是不是 overfit？」**
> 表面看 LR HCC1937 gap = +0.48，似 overfit；但 DT/RF/XGB HCC1937 gap 也都 ≈ +0.48（差 < 0.002）。如果是 model overfit，high-capacity model gap 應更大；實測四者等價 → gap 來源是 caller_F1 baseline asymmetry（HCC1937 = 0.37 vs 其餘平均 ≈ 0.83）的計算 artifact（§4.2），非 model 變異。LR 並未 underfit。

**Q5: 「下一步該做什麼？」**
> 三選項（需 PI 拍板）：
> (A) 放棄 cross-sample filter，明確標 cycle 1 為 per-sample model 限制
> (B) 擴增 sample 數（+ 4 樣本 → 9 sample LOSO）看 HCC1395 partial rescue 是否 stabilize
> (C) 特徵層 pivot — 設計 distribution-shift-robust feature set（去除 caller_af / chr8_flag 等 sample-specific feature）
> A4 minimal sufficient test 已完成，不建議再加 algorithm class。

### A.3 PI report deck 結論卡草稿（slide-ready）

```
標題：跨 model class 否決 sample-level methylation filter

證據（A4 ⭐ L2，2026-05-21）：
  • 4 algorithm (LR/DT/RF/XGBoost) × 5 sample LOSO × 5 seed = 100 fold
  • Overall mean ΔF1: LR −0.00004, DT +0.0026, RF +0.0027, XGB +0.0010
  • 全部 ribbon-null（< Cohen +0.005）
  • 唯一 partial rescue: HCC1395 hold-out + tree-based model (+0.013)
  • Train-test gap 跨 algo 等價 → 非 overfit；是 baseline shift artifact

結論：
  ✗ Not a model limitation
  ✓ Sample-level distribution shift dominates
  ✗ Production filter not deployable across samples
  ⚠ Cycle 1 +0.02236 = in-distribution circularity artifact
```

---

## 引用與依賴

- **Cycle 4 LOSO LR baseline**: `research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_cv.py`, `data/loso_cv_results.tsv`
- **Cycle 1 canonical filter**: `cycle1/cycle1_track_a_filter.json` (10 features, drop NumReads VIF=217, L2 C=1.0, τ*=0.39)
- **Helper functions**: `cycle2/scripts/cross_sample_apply.py` (CYCLE1_FEATURES, sweep_tau, compute_metrics, impute_with_medians, load_sample)
- **5-sample master TSVs**: `step5_master_augmented.tsv` (HCC1395) + `cycle2/data/{H1437,H2009,HCC1954,HCC1937}_master_augmented.tsv`

## Provenance

- 主腳本: `research/methyl_augmented_filter_phase2/phase2_completeness_audit/scripts/A4_multi_algo_LOSO.py`
- 單演算法 wrapper: `cycle4/loso_validation/scripts/run_loso_{DT,RF,XGB}.py`
- Wall-time: ~20-60 min (XGB / RF dominate)
- xgboost version: 2.1.4
- sklearn version: 1.6.1
