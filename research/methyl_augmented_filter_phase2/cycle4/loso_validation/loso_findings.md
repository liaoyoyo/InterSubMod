<!--
建立時間: 2026-05-20
agent: main session Coordinator
status: complete
report_class: cycle critical validation finding
audience: PI / user / cycle 1 framing 修正
scope: Phase 2 sample-level LOSO cross-validation
parent_plan: /bip7_disk/liaoyoyo2001/.claude/plans/v6-optimized-wadler.md v4.0
predecessor: 20260520_cycle3_step1_5_ism_ablation_vestigial
trigger: user 2026-05-20 query "LR filter 用 HCC1395 數據訓練,又用 HCC1395 數據驗證,這敘述合理嗎"
verdict: Cycle 1 +0.02236 **100% by sample-level circularity** — filter has no cross-sample generalization signal
last_verified: 2026-05-20
-->

# LOSO Cross-Validation — Phase 2 核心質疑 直接驗證

> **核心發現**: Cycle 1 LR filter HCC1395 +0.02236 effect size **100% 來自 sample-level circularity bias**。當 HCC1395 從 training set 移除（LOSO）後，filter 對 HCC1395 ΔF1 = **-0.00012**（drop +0.02248 = 完全失效）。
>
> **5 樣本 LOSO ΔF1 全部 ≈ 0** (-0.00012 ~ +0.00000), Wilcoxon p=0.125, 4/5 marginally negative, 1/5 zero — **DIRECTION_NEGATIVE**。
>
> 所有 best τ 退化為 0.10 (=keep everything) → **4-sample-trained LR 找不到 useful filter threshold for any held-out sample**。

---

## 0. TL;DR — 直接回答用戶 2026-05-20 質疑

**用戶質疑**: 「LR filter 用 HCC1395 數據訓練,又用 HCC1395 數據驗證,這敘述合理嗎」

**答**: **質疑完全成立**。LOSO 直接證實：

| 驗證模式 | HCC1395 ΔF1 | 解讀 |
|---|--:|---|
| Cycle 1 5-fold OOF (row-level, sample circular) | **+0.02236** | row-level 不 overlap 但 sample-level 完全 circular |
| **LOSO (sample-level, HCC1395 held out)** | **-0.00012** | filter 失效 (best τ→0.10 keep everything) |
| **效果完全來自 circularity bias** | **+0.02248** | 100% effect size 流失 |

→ Cycle 1 +0.02236 **不是真正的 filter effect**，是 **HCC1395 自身 distribution 的 in-sample fit**。

---

## 1. Method

### 1.1 LOSO Protocol

```
For each test_sample in [HCC1395, HCC1937, HCC1954, H1437, H2009]:
    train_samples = other 4 samples
    train_data = concat(master_TSV[train_samples])  # 共 ~143k - 263k rows
    
    # Per-train median impute (避免 test leakage)
    train_medians = compute_medians(train_data, 10 features)
    
    # StandardScaler fit on train only
    scaler = StandardScaler().fit(X_train)
    
    # L2 LR (C=1.0, lbfgs, seed=42)
    clf = LogisticRegression().fit(scaler.transform(X_train), y_train)
    
    # Apply to held-out test (use train medians + train scaler)
    X_test = scaler.transform(impute(test_data, train_medians))
    p_test = clf.predict_proba(X_test)[:, 1]
    
    # Sweep τ ∈ [0.10, 0.95] → best ΔF1
    ΔF1[test_sample] = max(τ ∈ TAU_GRID): compute_ΔF1(p_test, τ)
```

### 1.2 與 cycle 1+2+3 的差別

| Validation | Train sample(s) | Test sample | Sample leakage? |
|---|---|---|---|
| Cycle 1 5-fold OOF | HCC1395 (80% folds) | HCC1395 (20% fold) | ✅ row-level OK / ❌ **sample-level circular** |
| Cycle 2 transfer | HCC1395 (cycle 1 coef) | 4 新 samples | ✅ true cross-sample (但 single train sample) |
| Cycle 2 refit | sample s itself | sample s (per-sample OOF) | ❌ **sample-level circular per sample** |
| Cycle 3 Step 1 gate | per-sample re-fit (n=2) | same n=2 | ❌ **circular within qualifying subset** |
| **LOSO (this)** | **4 samples combined** | **1 held-out sample** | ✅ **真正 sample-level CV** |

LOSO 是 Phase 2 第一次跑「真正 sample-level cross-validation」。

---

## 2. Results

### 2.1 Per-sample LOSO ΔF1

| Test sample | Caller F1 | LOSO ΔF1 | Best τ | TP kept | FP kept | Verdict |
|---|--:|--:|--:|---:|---:|---|
| HCC1395 | 0.7166 | **-0.00012** | 0.10 | 30,490 | 4,841 | filter 無效 (keep all) |
| HCC1937 | 0.3692 | **+0.00000** | 0.10 | 13,910 | 2,697 | filter 無效 |
| HCC1954 | 0.8385 | **-0.00008** | 0.10 | 19,449 | 687 | filter 無效 |
| H1437 | 0.8670 | **-0.00001** | 0.10 | 70,191 | 773 | filter 無效 |
| H2009 | 0.8863 | **-0.00001** | 0.10 | 135,359 | 1,342 | filter 無效 |
| **Mean** | — | **-0.00004** | — | — | — | — |
| **Median** | — | **-0.00001** | — | — | — | — |

### 2.2 統計檢定

- Wilcoxon paired vs 0 (n=5): **p = 0.125** (n=5 exact min p=0.0625 requires 5/5 same sign; 0/5 positive)
- Sign: **0 positive / 4 negative / 1 zero**
- **Verdict: DIRECTION_NEGATIVE**

### 2.3 Best τ 全部退化為 0.10

5 樣本 best τ 全部 = **0.10** (lowest in grid)。
- τ = 0.10 means **保留所有 P(TP) ≥ 0.10 的 row**
- 與「不 filter / keep everything」實質相同
- 表示 **LR fitted on 4 samples 對 held-out sample 沒有任何 useful discriminative threshold**

---

## 3. Mechanism 解讀

### 3.1 為何 LOSO 完全失效

**Cycle 1 +0.02236 = HCC1395 自身 distribution overfit**

| Cycle 1 LR (HCC1395-trained) | LOSO LR (4-sample combined trained) |
|---|---|
| caller_af coef +3.44 (HCC1395 AF distribution) | caller_af coef 變化 (4-sample average AF distribution) |
| 完美匹配 HCC1395 自身 TP/FP boundary | 對任何單一 sample 都 sub-optimal |
| τ*=0.39 = HCC1395-specific optimum | 4 sample 平均 τ optimum 漂移 |

cycle 2 transfer mode (HCC1395 → 4 新 samples) NEGATIVE 已暗示這個問題；LOSO 從另一個方向確認 — **即使 4 樣本一起 train，也找不到對任何 sample useful 的 filter**。

### 3.2 為何 cycle 2 transfer mode 仍有 HCC1395 +0.02232？

Cycle 2 transfer mode 是 **HCC1395 自己 train + HCC1395 自己 test** → +0.02232 是 HCC1395 in-distribution 結果。
而 cycle 2 對 4 新 samples 的 transfer 是 NEGATIVE (4/5)。

LOSO HCC1395 -0.00012 ≠ cycle 2 transfer HCC1395 +0.02232，因為：
- Cycle 2 用 cycle 1 HCC1395-trained coef apply 在 HCC1395 → circular
- LOSO 用 4-sample-trained coef apply 在 HCC1395 → 真正 held-out
- **差距 +0.02244 = sample-level circularity 量化**

---

## 4. 與 cycle 1-3 結論的 reconciliation

### 4.1 Cycle 1 結論需修正

| 之前 (cycle 1 + Phase 2 PI Trust HTML) | LOSO 之後 修正 |
|---|---|
| "HCC1395 internal validated ⭐3 strong" | "HCC1395 **in-distribution 5-fold OOF +0.02236**, sample-level 不 valid (LOSO -0.00012)" |
| "Phase 2 KPI1 = 224% 達成" | "KPI1 (HCC1395 in-distribution) 224% / **KPI1 (sample-level LOSO) ~0%**" |
| Tier ⭐⭐⭐⭐ L2 | **降為 ⭐⭐ L4 anecdotal** (single-sample in-distribution only) |

### 4.2 Cycle 3 Step 1.5 ablation 結論不變

ISM vestigial in LR 結論仍成立，**但 scope 更窄**:
- 之前: "ISM vestigial in cycle 1 LR framework"
- LOSO 後: "ISM vestigial **even within HCC1395 in-distribution circular validation** — true cross-sample LR 本身就無 signal, ISM contribution discussion 變得 moot"

### 4.3 Cycle 2 transfer NEGATIVE 結論強化

Cycle 2 transfer mode 4/5 NEGATIVE 已 indicate sample-level 失效；LOSO 從不同方向 (4-sample combined train) 完全 confirm — **sample-level LR filter signal 不存在**，無論 train sample 多少。

### 4.4 Cycle 3 gate rule 結論需重評

| 之前 | LOSO 之後 |
|---|---|
| "qualifying subset (HCC1395+HCC1937) mean +0.01499 PASS" | "qualifying subset PASS 是 per-sample circular refit；LOSO 後 HCC1937 = +0.00000 證實 sample-level 也無 signal" |
| "Caller-F1-headroom mechanism" | mechanism 部分仍真實 (HCC1937 caller F1 + FP density 確實 qualifying)，但 **filter 本身在 sample-level 無效**，gate rule 是「在無效 filter 上加 gate」，整體仍 ⭐⭐ marginal |

---

## 5. Phase 2 真正結論（LOSO 後）

| Claim | 修正前 tier | 修正後 tier |
|---|---|---|
| HCC1395 in-distribution ΔF1 +0.02236 | ⭐⭐⭐⭐ L2 | **⭐⭐ L4** (single-sample, no sample-level cross-val) |
| BAM-invariant V3F/V5/V6 | ⭐⭐⭐⭐ L2 | **⭐⭐⭐⭐ L2 unchanged** (BAM 跨 BAM 是不同問題) |
| Cross-sample transfer NEGATIVE | ⭐⭐⭐ L3 | **⭐⭐⭐⭐ L2** (LOSO 確認) |
| ISM vestigial in LR | ⭐⭐⭐⭐ L2 | **⭐⭐⭐ L3** (scope 變窄 — LR 本身就無 sample-level signal) |
| Caller-F1-headroom mechanism | ⭐⭐⭐ L3 | **⭐⭐ L4** (mechanism 仍真實但 filter 整體無效讓機制變 moot) |
| Cycle 3 qualifying mean +0.01499 | ⭐⭐ | **⭐ L4** (per-sample circular per cycle 3) |
| **LOSO sample-level negative ΔF1** | (new) | **⭐⭐⭐⭐ L2** (本 LOSO 驗證的核心新結論) |

**整體 Phase 2 verdict**: 從「⭐3 internal valid with cross-sample caveat」**降為「sample-level FILTER FAILED, in-distribution case study only」**。

---

## 6. SoT Update Action Items

| Artifact | Action |
|---|---|
| **PI Trust HTML Section 2 dashboard** | KPI1 加 "LOSO-corrected" row 顯示 -0.00012 / KPI2 升為 main; 整體 banner 改 caution → critical |
| **PI Trust HTML Section 3 Evidence Ladder** | HCC1395 +0.02236 ⭐⭐⭐⭐→⭐⭐ ; 加新 LOSO claim ⭐⭐⭐⭐ |
| **PI Trust HTML Section 7 Uncertainty** | 移除「caller-F1-headroom 是真實機制 L3」(降為 L4) ; 加 LOSO mechanism 確認 |
| **Cycle 1 主報告** (`InterSubMod/docs/experiments/in_progress/2026/05/20260518_Phase2_Cycle1_Global_FP_Filter_01.md`) | 加 head banner "LOSO 2026-05-20 reframe: HCC1395 in-distribution only, sample-level filter failed" |
| **Memory** `project_phase2_cycle1_global_fp_filter.md` | description + body LOSO caveat |
| **evidence_ledger** | append `20260520_loso_sample_level_circularity_revealed` |
| **INDEX.md** | Cycle 1 entry status 改為 ⭐2 (in-distribution case study) + Cycle 4 LOSO entry |
| **CURRENT_FOCUS** | prepend 2026-05-20 LOSO 結論 section |
| **Paper §3 framing** | 完全撤回 "ISM-augmented filter" 宣稱; 改為 "ISM characterization study + LR filter sample-level negative finding" |

---

## 7. 下一步建議

### Cycle 4 Trial A/B/C 是否仍 GO？

| Trial | LOSO 後重評 | 建議 |
|---|---|---|
| Trial A (Interaction LR + 5 Python new features) | LOSO 已證實 LR 本身在 sample-level 無 signal；interaction 與 new features 在 LR framework 內**很難**翻轉這個結論 | **降級** — 可選跑但 prior PASS 從 75% 降至 ~15% |
| Trial B (RF/XGBoost) | 非線性 model 仍受 sample-level circularity 同等限制 | 降級 prior 從 30% 至 ~10% |
| Trial C (Per-zone LR) | per-zone heterogeneous 不改變 sample-level 問題 | 降級 prior 從 35% 至 ~10% |

### Pivot 建議（強化 ROI 評估）

| Path | LOSO 後 ROI |
|---|---|
| **Pivot phase_block_3d** | **顯著上升** — 在不同 framework 重新尋找 cross-sample signal |
| **Pivot thread_d Tier 2** | **顯著上升** — 已知 signature track 不依賴 LR filter |
| Cycle 4 Trial A/B/C | 顯著下降 (LOSO 已 implicitly 跑了一個 generalized 4-sample-trained LR test) |
| BAM-level new feature | 仍未測，但需要 cycle 5+ if 真要繼續 |
| Phase 2 final NEGATIVE archive | **最 honest** — accept LOSO 結論 |

---

## 8. Reproducibility

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_cv.py
```

Wall clock: ~30 sec. Deterministic — PRIMARY_SEED=42, lbfgs solver, identical to cycle 1+2 protocol. Reuses cycle 2 cross_sample_apply.py functions (CYCLE1_FEATURES, compute_metrics, impute_with_medians, sweep_tau, load_sample, reverse_solve_fn).

---

## 9. Files

```
cycle4/loso_validation/
├── loso_findings.md                    ← this report
├── scripts/run_loso_cv.py              (~10 KB, reuse cycle 2 functions)
├── data/loso_cv_results.tsv            (5 rows × 22 cols)
├── figures/loso_5sample_dF1.png        (5-sample bar with annotation)
└── intermediate/loso_summary.json      (machine-readable verdict)
```

---

## 10. 答覆用戶 2026-05-20 質疑

> "LR filter 用 HCC1395 數據訓練,又用 HCC1395 數據驗證,這敘述合理嗎"

**這個質疑直接導致 Phase 2 整體結論 reframe**：

1. **Row-level 5-fold OOF 合理** — 避免 train-test 同行
2. **Sample-level circular 不合理** — LOSO 直接證實
3. **HCC1395 +0.02236 = 100% sample-level circularity bias** — 移除後 -0.00012
4. **整個 LR filter framework 在 sample-level 無 signal** — 5 samples LOSO 全 ≈ 0, best τ → 0.10
5. **ISM characterization 結論不變** (v0.3 cycle ⭐3 仍保留)
6. **Production filter direction 失敗** — Phase 2 後續 archive 或 pivot

---

**End of LOSO Findings — Phase 2 核心 reframe trigger**
