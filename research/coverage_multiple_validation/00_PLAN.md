# Coverage_Multiple 作為 CN proxy 的可靠性驗證

> **狀態**：blocked — 等待以現行 binary 重跑 master dataset（2026-04-19 更新）
> **建立日期**：2026-04-19
> **專案目錄**：`research/coverage_multiple_validation/`

## 2026-04-19 狀態更新

Step 1/2 以 `/big7_disk/.../20260330_.../all_region_rows.tsv.gz` 為資料來源，該 master dataset
由 2026-03-30 時的 binary 產出；但 C++ KDE per-sample diploid 估計於 **2026-04-13** 才經 commit
`8d0a0c8` 啟用，master dataset 因此全部 rows 的 expected_coverage=75.0（stale-binary artifact）。

2026-04-19 C++ 再次強化（commits `374fad4` + `12d9b3e`）：
- KDE fallback 失效時升為 WARN 並帶失效原因
- TSV 新增 `Diploid_Coverage_Used` 欄位供後續直接讀取
- 回傳型別改為 `DiploidEstimate{value, used_fallback}` 避免 value==75.0 sentinel 誤判

→ Step 1/2 的 per-CN-bin recall 與 per-sample bias 量化**需以現行 binary 重跑 master dataset**
才能確立最終結果；目前報告的數值仍具研究價值（證實 baseline 失準確實發生），但 root cause
已從「CLI fallback 未啟用 KDE」更正為「stale binary」。

## 背景與動機

Z3 amplicon blacklist pilot (2026-04-19) 顯示 HCC1954 whole-chr5/8/17 黑名單能帶來 ΔF1=+0.0065（local 改善），但其他 6 樣本 5/6 hurt（mean ΔF1=−0.0044）→ 無法作為 global canonical filter。

診斷發現 HCC1954 Z3 FP Coverage_Multiple=0.73 vs TP=0.49（比值 ≈ 1.5×），這與「HER2 amplicon + 8p loss + 17p TP53 LOH」的 CN 架構一致。然而跨樣本失敗提示：**Coverage_Multiple (KDE-estimated diploid baseline) 本身可能是失準的 CN proxy**，導致其他樣本誤判。

若 CovM 是可靠 proxy，Z3 amplicon blacklist 應可泛化；若 CovM 失準，需以更準確的 CN 度量（per-sample re-centered baseline 或 SEQC2 oracle CN）重驗才能定案。

## 假說

### H1：Coverage_Multiple 可作為跨樣本可靠的 CN proxy

**陳述**：KDE 估計的 diploid baseline 在所有 7 樣本上都能反映真實 CN（≥0.8 recall）。

**前提條件**：
- HCC1395 有 SEQC2 CNV gold standard 可供 ground truth 驗證
- 其他 6 樣本可用 pan-cancer proxy-neutral 區間（chr1p/2q/3p/9q）近似中性基線

**已知 Confound**：
- expected_coverage 的 KDE 估計本身可能被高 ploidy 區域污染而上抬
- 樣本間絕對測序深度差異
- Proxy-neutral 區間在特定癌症類型仍可能異常

**驗證標準**：
- **Positive**：HCC1395 Gain recall ≥ 0.80 **且** 跨樣本 baseline bias |Δ| < 10%
- **Negative**：Gain recall < 0.60 **或** baseline bias > 20% **或** expected_coverage 為硬編碼值

**現階段跡象**（Step 1 已完成）：HCC1395 Gain recall = 14.6% → 強烈 leaning rejected。

### H2：Per-sample re-centered baseline 或 oracle CN 能修正 Z3 amplicon blacklist 的跨樣本 collateral damage

**陳述**：若 H1 rejected，採用：
- **Variant A (Re-center)**：`CovM_v2 = NumReads / neutral_median_NumReads`（per-sample）
- **Variant B (Oracle CN)**：以 SEQC2 truth CN 直接取代 CovM（僅 HCC1395 可得）

**前提條件**：
- H1 確認原始 CovM 有偏誤
- Step 2 已量化各樣本 baseline shift 方向/幅度

**已知 Confound**：
- Re-center 所用 proxy-neutral 區間本身可能含 CNV（H2009/H1437 特別）
- Oracle CN 僅 HCC1395 可得，無法作為 pipeline 通用解
- 修正後 CN 範圍可能與 S2 blacklist 閾值不相容，需重新校準

**驗證標準**：
- **Positive**：Variant A/B 使 HCC1954 ΔF1 ≥ +0.005 **且** 其他 6 樣本 mean ΔF1 ≥ 0
- **Negative**：HCC1954 改善消失，或其他樣本仍 mean ΔF1 < 0

## 方法

### 數據來源

| 數據集 | 路徑 | 描述 |
|--------|------|------|
| Master region | `/big7_disk/.../20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz` | 7 樣本 × paired/to 全量 region-level |
| SEQC2 CNV truth | `/big8_disk/data/HCC1395/SEQC2/CNV/` | HCC1395 gold standard gain/loss/LOH + CN 值 + exclusion |

### 分析步驟

```
Step 1: HCC1395 per-region × SEQC2 truth 對齊，計算 confusion matrix + CN bin metrics
        → 驗證: Gain recall, Loss recall, CovM vs truth CN correlation（Pearson/Spearman）
        → 產出: data/step1_{confusion_matrix,per_cn_bin_metrics,correlation_per_label}.tsv

Step 2: 跨 7 樣本估計 expected_coverage baseline bias
        → HCC1395 用 SEQC2 Neutral；其他 6 樣本用 pan-cancer proxy-neutral 區間
        → 驗證: bias = (est_expected - true_neutral_median) / true_neutral_median
        → 產出: data/step2_expected_coverage_per_sample.tsv

Step 3: [跳過] 原規劃的 bias 修正敏感性分析併入 Step 4 Variant A

Step 4: Oracle CN pilot — Variant A (re-center) + Variant B (HCC1395 oracle)
        → 驗證: 各 variant 在 HCC1954 ΔF1 與其他 6 樣本 ΔF1 的分佈
        → 產出: data/step4_{covm_v2_per_sample,oracle_pilot_results}.tsv
```

### 統計方法

- **相關性**：Pearson + Spearman（per label + global）
- **Recall**：per CN bin（0-1, 1-2, 2, 3, 4, ≥5）precision/recall
- **ΔF1**：per sample × strategy，與原始 Z3 pilot ceiling 對比

## 可行性評估

| 因素 | 評估 |
|------|------|
| 數據可用性 | ✓ Master + SEQC2 CNV BED 就位 |
| 計算資源 | ✓ 單機 pandas 可處理 |
| 與已有結論衝突 | ✓ 補強 Z3 pilot CONDITIONAL 結論的根因分析 |

## 已知風險

1. **Proxy-neutral 區間本身含 CNV**：H2009/H1437 等小細胞肺癌常見 chr3p loss → chr3p 不適合作為 proxy neutral。緩解：若 Step 2 發現該區間偏離嚴重，改用 chr11/13 其他 region。
2. **Oracle CN 僅 HCC1395**：無法泛化至其他樣本；H2 Variant B 只能作為 upper bound 論證，不構成 pipeline 通用解。
3. **S2 閾值與 CovM_v2 scale 不相容**：re-center 後的閾值需重新校準；若 Step 4 發現 ΔF1 退步僅因 threshold 未調，需 parameter sweep。

## 相關檔案

- 觸發來源：`docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md`
- 上游依賴：`research/z3_internal_feature_exploration/`（S2 原始 pilot）
- Master dataset：observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix
- 報告輸出位置：`reports/`（待 Step 4 完成後撰寫總報告）
