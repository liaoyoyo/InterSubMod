---
created: 2026-04-20
status: complete
scope: KDE expected_coverage C++ 修正後下游量化 cycle（Step 0-6）
---

# KDE Fix Downstream Validation

## 目的

C++ KDE expected_coverage 修復（commits `374fad4` + `12d9b3e`）後，量化四大下游影響：

1. **跨樣本 bias 真實收斂**（原 ±39% 僅 HCC1395，全 7 樣本 × 2 mode 跨度？）
2. **H-CN1 verdict 是否翻轉**（原 stale 75× 下 Gain recall 14.6% 的 artifact）
3. **Z3 amplicon blacklist ΔF1**（CovM 閾值下游效應）
4. **COLO829 既有結論**（ratio 2.59 最極端樣本）

**不回跑 C++ binary**（已經 Step 1-4 驗證通過）；**不改既有 master dataset**（只對比）。

## 目錄結構

```
research/kde_fix_validation/
├── scripts/
│   ├── build_kde_rerun_master.py      # Step 0: 14 combos → 單一 TP-only master
│   ├── step1_hcc1395_seqc2_kde.py     # Step 1: HCC1395 vs SEQC2 CN truth
│   ├── step2_z3_amplicon_kde.py       # Step 2: Z3 blacklist S1/S2/S3 ΔF1
│   └── step3_quantile_drift.py        # Step 3: 跨樣本 CovM 分位數漂移
├── outputs/
│   ├── step1_hcc1395_seqc2/
│   │   ├── paired_pileup/{per_cn_bin_metrics,confusion_matrix}.tsv + fig
│   │   └── paired_full/...
│   ├── step2_z3_amplicon/
│   │   ├── step2_summary.tsv          # 7 samples × S1/S2/S3(stale)/S3(fixed)
│   │   └── step2_compare.md
│   ├── step3_quantile_drift/
│   │   ├── quantile_drift_per_sample.tsv
│   │   └── category_migration_per_sample.tsv
│   └── step4_colo829_audit/
│       └── impacted_conclusions.md    # 9 筆結論：4 CRITICAL + 3 HIGH + 2 LOW
└── README.md (本檔)
```

**圖輸出**：`docs/experiments/in_progress/2026/04/figures/20260420_kde_fix_acceptance/`
- `fig5_quantile_drift_per_sample.png`（14 panels）
- `fig6_category_migration_per_sample.png`

## 4-步驟量化結果摘要

### Step 1 — HCC1395 SEQC2 Gain/Loss recall 重算

| 指標 | Stale 75× | Fixed KDE | 提升 |
|------|----------:|----------:|-----:|
| paired_pileup Gain recall | 14.6% | **41.87%** | ×2.87 |
| paired_full Gain recall | — | **45.78%** | ×3.14 |
| Loss recall | — | 65-68% | — |
| Spearman(CovM, SEQC2 CN) | — | **0.845** | — |

→ **H-CN1 = 🟢 PARTIAL POSITIVE**（方向性 proxy 成立，定量仍需輔助特徵）

### Step 2 — Z3 amplicon blacklist ΔF1

| 策略 | Stale Σ\|ΔF1\| | Fixed Σ\|ΔF1\| | 差異 |
|------|--------------:|--------------:|-----:|
| S1 文獻 blacklist（coordinate-based） | — | — | **0.0000**（sanity pass） |
| S2 whole-chr | — | — | **0.0000**（sanity pass） |
| S3 CovM 95%ile | — | — | **0.0000** |

→ S3 strategy **scale-invariant**（per-sample 95%ile 閾值隨 baseline 等比縮放，mask 成員不變）。原 Z3 pilot CONDITIONAL NEGATIVE verdict **不受 KDE 修正影響**。

### Step 3 — Cross-sample CovM quantile drift

| Sample | Stale/KDE ratio | Δp50(CovM) |
|--------|----------------:|-----------:|
| COLO829 | 2.59 | **+0.613** |
| HCC1395_DORADO | 1.42 | +0.364 |
| HCC1395 | 1.36 | +0.343 |
| HCC1954 | 1.23 | +0.207 |
| H1437 | 1.09 | +0.081 |
| H2009 | 0.95 | −0.060 |
| HCC1937 | 0.82 | **−0.241** |

→ Δp50 方向與 `ratio − 1` 完全一致；COLO829 右移最劇、HCC1937 反向漂移（stale 75× 高估 true 91×）。

### Step 4 — COLO829 結論審計

9 筆影響：
- 🔴 **CRITICAL ×4**：O1-O10 Fig 11 QS 解釋、TO LOH M2 mask 91.7%、ISM-for-melanoma 假說、O11-O13 類跨樣本 QS 觀察
- 🟠 **HIGH ×3**：LOH 內外 ISM 區分力、PCA 跨樣本分離、LOH enrichment post-HP-fix
- 🟡 **LOW ×2**：planning doc（trash）、Z3 cycle（本輪已連動）

→ 詳見 `outputs/step4_colo829_audit/impacted_conclusions.md`

## 下一 cycle 必要 follow-up（本 cycle 不做）

1. COLO829 QS 重跑（依賴 CNV_Loss penalty 的 QS 可能大幅回升）
2. M2/M4 mask 重算（91.7% → 預期 <30%）
3. PCA 跨樣本分離重跑
4. 以新 QS 重測「ISM 對黑色素瘤無效」假說

## 對照文件

- **主結論**：`docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md`
- **Baseline 方法驗證**：`docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md` §4.4
- **C++ 修正歷程**：`docs/methodology/20260419_KDE_expected_coverage_audit_01.md`
