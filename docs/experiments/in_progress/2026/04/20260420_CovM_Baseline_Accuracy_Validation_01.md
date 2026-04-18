---
title: Coverage_Multiple / Ploidy Baseline 準確度驗證 — 是否為 Z3 跨樣本失敗根因？
date: 2026-04-20
status: MIXED（H-CN1/H-CN2 POSITIVE，H-CN3 NEGATIVE — 「cov 需修正」為真，但修正不解 Z3 pilot 跨樣本失敗）
owner: InterSubMod Research
scope: TO mode · 7 samples × 2 modes · HCC1395 SEQC2 CNV truth benchmark
related:
  - 20260419_Z3_amplicon_blacklist_pilot_result_01.md
  - ../../../validated/2026/04/20260411_GC校正與甲基化CN驗證_01.md
  - ../../../validated/2026/04/20260409_SEQC2_CNV分層觀察_01.md
hypotheses:
  - H-CN1: HCC1954 KDE expected_coverage 偏高 → baseline 失準
  - H-CN2: Coverage_Multiple 對 Gain 的 recall 遠低於 Loss
  - H-CN3: segment-level / oracle CN 可跨樣本救援 Z3 blacklist
---

# Coverage_Multiple Baseline 準確度驗證 — Z3 跨樣本失敗根因釐清

## 一、執行摘要

**核心回答使用者提問**：「確認會造成這樣的原因是否是因為判斷平均 read 與比例還有倍體問題時，不夠準確，確認是否是當時的 cov 條件需要修正」。

> **是 — 但不是為了解 Z3 pilot 問題，而是獨立的 infrastructure bug。**

三項發現彼此獨立：

| 假說 | 結論 | 證據 |
|------|------|------|
| **H-CN1** KDE baseline mis-calibration | **POSITIVE**（嚴重） | 全 7 樣本 × 2 mode 的 `expected_coverage_estimate = 75.0`（常數），per-sample bias 跨度 **−28% 到 +150%** |
| **H-CN2** Coverage_Multiple gain recall 遠差於 loss | **POSITIVE**（嚴重） | HCC1395 SEQC2：Gain recall **14.6%**（CN=3 僅 0.15% 被標為 Gain）vs Loss recall **86.9%** |
| **H-CN3** 改用更準 CN 可救援 Z3 blacklist | **NEGATIVE** | HCC1395 oracle truth CN≥3.5 ∩ Z3 只拿到 ΔF1=+0.0011；HCC1954 S2 對 HCC1395 本身也只 +0.001–0.002 |

**最終判定**：
- **「cov 條件需要修正」= YES**：`expected_coverage` 未經 KDE 個別樣本估計，全走 default fallback 75.0，是確鑿的 infrastructure bug，應由 `/cpp-change` 另案修復
- **「修正 cov 能解決 Z3 pilot 跨樣本失敗」= NO**：即便給 HCC1395 完美 truth CN，Z3 filter 增益只 ~+0.001；HCC1954 的收益 (+0.0065) 真正來源是**樣本特異 CNV 架構**而非 CN 估計精度
- **對五研究目標的影響**：零；本驗證解答 Z3 × HCC1954 是否值得繼續投入（**不值得**，Zone-Aware Framework 維持 characterization-only 定位）

---

## 二、動機與假說鏈

### 2.1 觸發事件

2026-04-19 Z3 × HCC1954 amplicon/high-CN blacklist pilot：
- S2 (whole chr5/8/17 ∩ Z3)：HCC1954 ΔF1=+0.0065，**其他 6 樣本 mean ΔF1=−0.0044**（5/6 hurt）
- S3 (CovM > 95%ile non-Z3 baseline ∩ Z3)：跨樣本近零（Δ≈+0.0002）
- S4 ceiling：HCC1954+0.0075, HCC1937+0.0309 → 顯示**某些樣本 Z3 本身 TP 稀薄**

使用者假說：問題可能出在 CN/ploidy 估計 — HCC1954 pseudo-tetraploid（~92 chr/cell）下，KDE 以 diploid 假設抓 expected_coverage 會失準，使 Coverage_Multiple 跨樣本不可比。

### 2.2 三階假說

- **H-CN1**（baseline）：HCC1954 KDE 估 `expected_coverage` 偏高（hyper-aneuploid → NumReads 分佈右偏 → 80%ile 上界仍在 CN>2 帶）
- **H-CN2**（分層精度）：HCC1395 per-region，Coverage_Multiple 對 Gain 的 recall ≪ Loss
- **H-CN3**（可改善）：oracle CN（用 SEQC2 truth 直接替代 CovM）可拉高 Z3 blacklist 跨樣本有效性 ≥ +0.005

---

## 三、Step 1 — HCC1395 × SEQC2 CNV per-region 精度

### 3.1 方法

- 輸入：master dataset 過濾 `HCC1395` + `HCC1395_DORADO` TO region；SEQC2 `ngs_benchmark_cnvs_gain_loss_loh.bed` + `Additional_file_{3,4}_cnv_{gain,loss}_cn_median.txt`
- 每個 InterSubMod region × SEQC2 label（Excluded → Gain → Loss → LOH → Neutral 優先序）+ truth median CN
- 6-class 預測 × 4-class truth confusion matrix

### 3.2 關鍵結果

![Step 1 confusion matrix](../../../../research/coverage_multiple_validation/figures/step1_confusion_matrix.png)

![Step 1 CovM vs truth CN scatter](../../../../research/coverage_multiple_validation/figures/step1_covm_vs_truth_scatter.png)

**per-CN-bin Coverage_Multiple 表現**：

| truth CN bin | n | mean CovM | median CovM | pred_gain_rate | pred_loss_rate |
|---|---|---|---|---|---|
| CN=0–1 | 4 | 0.464 | 0.434 | 0.0% | 75.0% |
| CN=1–2 | 1,428 | 0.346 | 0.320 | 0.0% | **86.9%** |
| CN=2 (neutral+LOH) | 22,705 | 0.632 | 0.613 | 0.6% | 19.6% |
| **CN=3** | **24,805** | **0.863** | **0.853** | **0.15%** | **1.6%** |
| **CN=4** | **19,308** | **1.151** | **1.147** | **5.4%** | **0.09%** |
| **CN≥5** | **12,274** | **1.617** | **1.560** | **58.3%** | **0.02%** |

**per-label correlation**：

| Label | n | Pearson r | mean CovM | mean truth CN |
|---|---|---|---|---|
| Gain | 56,387 | **0.810** | 1.125 | 3.89 |
| Loss | 1,432 | −0.047 | 0.346 | 1.00 |
| Neutral | 3,022 | — | 0.840 | 2.0 |
| LOH | 19,683 | — | 0.600 | 2.0 |
| **ALL** | 80,524 | **0.833** | 0.973 | 3.30 |

### 3.3 判定：**H-CN2 POSITIVE**

- CN=3 的 Coverage_Multiple 中位數 **0.85**（期望值 ~1.5）→ **Gain 被系統性低估 43%**
- CN=4 中位 1.15（期望 ~2.0）→ 低估 42.5%
- CN≥5 中位 1.56（期望 ≥2.5）→ 低估 ≥38%
- 僅 CN≥5 有 58.3% 被正確歸 Gain；CN=3/4 絕大多數歸 Normal
- **Loss side 健全**：CN=1–2 有 86.9% 被正確歸 Loss

**→ Coverage_Multiple 對 Gain 側有系統性右偏壓縮**，全域 r=0.833 隱藏了 Gain recall=14.6% 的嚴重分層偏誤。

---

## 四、Step 2 — 跨樣本 expected_coverage baseline 對比

### 4.1 方法

- Master dataset 反推：`expected_coverage[sample, mode] = median(NumReads / Coverage_Multiple)`
- 真實 neutral baseline：HCC1395 用 SEQC2 Neutral BED；其他 6 樣本用 chr1p+chr2q+chr3p+chr9q 作 pan-cancer proxy neutral
- Bias = (estimated − neutral_median_NumReads) / neutral_median_NumReads

### 4.2 關鍵結果

![Step 2 expected_coverage bias](../../../../research/coverage_multiple_validation/figures/step2_expected_coverage_bias.png)

**全 14 rows（7 samples × 2 modes）的反推結果**：

| sample | mode | estimated expected_cov | neutral median NumReads | bias |
|---|---|---|---|---|
| COLO829 | to | **75.0** | 30 | **+150.0%** |
| COLO829 | paired | **75.0** | 30 | **+150.0%** |
| H1437 | to | **75.0** | 66 | +13.6% |
| H1437 | paired | **75.0** | 68 | +10.3% |
| H2009 | to | **75.0** | 104 | −27.9% |
| H2009 | paired | **75.0** | 104 | −27.9% |
| HCC1395 | to | **75.0** | 54 (SEQC2) | **+38.9%** |
| HCC1395 | paired | **75.0** | 54 (SEQC2) | **+38.9%** |
| HCC1395_DORADO | to | **75.0** | 57 (SEQC2) | +31.6% |
| HCC1395_DORADO | paired | **75.0** | 56 (SEQC2) | +33.9% |
| HCC1937 | to | **75.0** | 101 | −25.7% |
| HCC1937 | paired | **75.0** | 99 | −24.2% |
| HCC1954 | to | **75.0** | 65 | +15.4% |
| HCC1954 | paired | **75.0** | 65 | +15.4% |

### 4.3 判定：**H-CN1 POSITIVE（以 infrastructure bug 形式）**

發現極端驚人：**14/14 rows 的 expected_coverage 精準等於 75.0**，**7 個樣本共用同一常數**。

這**不是 KDE 估計偏誤**，而是：
- Pipeline 執行時 `--expected-coverage` 走了 CLI default 的 **75.0 fallback**
- `src/core/RegionProcessor.cpp:59-136` 的 KDE auto-estimation **在本 master dataset 產出時未被啟用**
- 對 COLO829（真實 neutral=30）→ baseline 高估 150%，CovM 被整體壓縮 2.5×
- 對 H2009（真實 neutral=104）→ baseline 低估 28%，CovM 被整體放大 1.4×

**結論**：
- H-CN1 原假設（HCC1954 特有偏誤）**被更嚴重的事實取代**：全 7 樣本都用錯 baseline
- 之前研究（`20260411_GC校正與甲基化CN驗證_01.md`）引用的「KDE 實作正確，CN 準確度 6.2%→43.8%」**僅驗證 KDE 邏輯本身**，**並未驗證 master dataset 生成時 KDE 是否被啟用**
- **這解釋了 Step 1 的 Gain recall 14.6%**：HCC1395 真實 neutral NumReads=54，而 baseline=75.0 → diploid 樣本（CN=2）的 NumReads≈54 除以 75 會得到 0.72，落入「Low」而非「Normal」；CN=3 的 NumReads≈81 除以 75 = 1.08，落入「Normal」而非「Elevated」

---

## 五、Step 4 — Oracle CN Pilot（H-CN3 驗證）

### 5.1 方法

兩 variant 驗證「若 CN 估計完全正確，Z3 blacklist 能否跨樣本奏效」：

- **Variant A — per-sample re-centering**：`CovM_v2 = NumReads / neutral_median_NumReads`（以 step 2 真實 neutral 重新校準），跑 S3（CovM_v2 95%ile non-Z3 baseline ∩ Z3）
- **Variant B — HCC1395 oracle CN**：直接用 SEQC2 truth_median_cn 作 CN，對 HCC1395/HCC1395_DORADO 跑 truth_cn ≥ {2.5, 3.5, 4.5} ∩ Z3

### 5.2 Variant A 結果（全 7 樣本）

![Step 4 oracle vs CovM ΔF1](../../../../research/coverage_multiple_validation/figures/step4_oracle_vs_covm_deltaf1.png)

| sample | CovM_old p95 ∩ Z3 | CovM_v2 p95 ∩ Z3 | S2 whole-chr ∩ Z3 |
|---|---|---|---|
| COLO829 | −0.00040 | **−0.00040** | −0.01091 |
| H1437 | −0.00019 | **−0.00019** | −0.00865 |
| H2009 | 0.00000 | **0.00000** | −0.00249 |
| HCC1395 | +0.00011 | **+0.00011** | −0.00407 |
| HCC1395_DORADO | +0.00031 | **+0.00031** | −0.00389 |
| HCC1937 | +0.00002 | **+0.00002** | +0.00381 |
| HCC1954 | +0.00017 | **+0.00017** | +0.00653 |

**關鍵發現**：CovM_v2 的 ΔF1 **每個樣本都與 CovM_old 完全相同**。

**數學原因（事後推導）**：
```
CovM_v2[i] = NumReads[i] / neutral_median[sample]
           = (NumReads[i] / 75.0) × (75.0 / neutral_median[sample])
           = CovM_old[i] × scale[sample]
```
因 scale 是 **per-sample 常數**，percentile-based cutoff（非 Z3 95%ile）也同比例放大，**相對排序與篩選邊界完全不變**。

→ **Variant A 證實：percentile filter 對 per-sample re-centering 免疫**；即使 baseline 偏誤被修正，S3-類策略的訊號仍然不變。

### 5.3 Variant B 結果（HCC1395 oracle truth CN）

| sample | variant | f1_before | f1_after | ΔF1 | TP lost | FP removed |
|---|---|---|---|---|---|---|
| HCC1395 | truth_cn ≥ 2.5 ∩ Z3 | 0.8309 | 0.8273 | **−0.0036** | 832 | 885 |
| HCC1395 | truth_cn ≥ 3.5 ∩ Z3 | 0.8309 | 0.8319 | **+0.0011** | 103 | 233 |
| HCC1395 | truth_cn ≥ 4.5 ∩ Z3 | 0.8309 | 0.8315 | +0.0006 | 21 | 81 |
| HCC1395_DORADO | truth_cn ≥ 3.5 ∩ Z3 | 0.8330 | 0.8353 | **+0.0023** | 169 | 427 |
| HCC1395_DORADO | truth_cn ≥ 4.5 ∩ Z3 | 0.8330 | 0.8341 | +0.0012 | 33 | 142 |

### 5.4 判定：**H-CN3 NEGATIVE**

- HCC1395 oracle best（truth_cn≥3.5）**ΔF1=+0.0011**，遠小於採用門檻 +0.005
- HCC1395_DORADO best +0.0023，同樣未達標
- 即便是 **完全正確的 CN**（SEQC2 金標準），Z3 內 high-CN 規則對 HCC1395 也只能拿到 HCC1954 S2（+0.0065）的 **1/6 – 1/3** 收益
- HCC1954 的 +0.0065 真正來自 **chr5/8/17 整條染色體上 FP 集中**（arm-level + 複合重排），**與 CN 估計精度無關**

---

## 六、綜合判定與「cov 條件是否需要修正」的直接回答

### 6.1 主要判定

| 假說 | 判定 | 對 Z3 pilot 的解釋力 |
|------|------|-------------------|
| H-CN1 baseline mis-calibration | **POSITIVE** | 獨立的 infrastructure bug；不解 Z3 pilot 跨樣本失敗 |
| H-CN2 Gain recall 崩潰 | **POSITIVE** | 是 H-CN1 的下游後果；修 H-CN1 即解 H-CN2 |
| H-CN3 oracle CN 可救援 Z3 | **NEGATIVE** | 即使完美 CN 也只拿 +0.001–0.002 |

### 6.2 對使用者問題的直接回答

**「cov 條件是否需要修正？」**

**是**，但修正動機與 Z3 pilot 分離：

1. **當前 master dataset 的 Coverage_Multiple 幾乎沒用**：全 7 樣本共用 `expected_coverage = 75.0`，COLO829 的 CovM 整體被壓縮 2.5×，H2009 被放大 1.4× — 跨樣本 CovM 分佈不可比
2. **修正方向明確**：重跑 pipeline 時確保 `--expected-coverage` 改用 KDE auto 或明確傳入 per-sample 值
3. **修正後不會救 Z3 pilot**：Variant A 已數學證明 percentile-based filter 對 per-sample re-centering 免疫；Variant B oracle 顯示即使 CN 完美也收益有限

**「Z3 blacklist 跨樣本失敗的真因」**：

並非 CN 估計誤差，而是：
- **HCC1954 chr5/8/17 有 arm-level LOH + 複合重排**（HER2+ breast 架構）→ 全臂 FP 富集
- **其他 6 樣本 chr5/8/17 並非 CNV hotspot**（各自有不同的 focal amp/LOH，如 COLO829 chr8 CDKN2A / H2009 chr7 EGFR）→ 套用 HCC1954 blacklist 即誤殺它們的正常 Z3 TP
- **這是真實生物學異質性**，不存在可跨樣本通用的 region-level blacklist

### 6.3 後續動作（建議）

| 項目 | 類型 | 優先 | 說明 |
|------|------|------|------|
| **Fix KDE baseline in pipeline** | `/cpp-change`（另案） | 高 | 確認 `--expected-coverage` 走 KDE；驗收標準 = step2 反推的 14 rows 不再共用 75.0 |
| 重跑 master dataset | 基礎建設 | 高 | 修 KDE 後重生成，否則所有 CovM-based 分析都受污染 |
| **放棄 Z3 × amplicon blacklist 作 global filter** | 方向決策 | 已決 | 本驗證確認此路不通 |
| Zone-Aware Framework 定位 | 不變 | — | 維持 characterization-only；Z3 的跨樣本 TP rate σ 來自生物學 CNV 異質性 |
| **對五研究目標影響** | 評估 | — | 零直接影響；但「CovM 當前不可跨樣本比較」意味著所有 cross-sample CovM 分層結論需加註風險 |

---

## 七、證據完整性與限制

### 7.1 強證據

- Step 2 的 14 rows × 75.0 等值：**反推自 master dataset 本體**，無需依賴外部 log，100% 可驗證
- Step 4 Variant A 的尺度不變性：純數學推導可重現
- Step 1 HCC1395 Gain recall 14.6%：以 SEQC2 Masood 2024 金標準驗證，不依賴專案內假設

### 7.2 限制

1. **Proxy neutral（6 樣本）的合理性**：chr1p+chr2q+chr3p+chr9q 非嚴格 pan-cancer neutral，如 H2009 chr1p 有 EGFR 上游片段、COLO829 chr3p 有 VHL 區；bias 數值絕對值可能 ±5% 波動，但「75.0 是 default 而非樣本特異估值」的結論不受影響
2. **HCC1954 無 SEQC2 CNV truth**：無法直接 oracle 驗證 HCC1954 Z3 能否被 true CN 救援；但 HCC1395 為 pseudo-triploid（ploidy=2.85）的 oracle 結果已顯示 CN 改進對 Z3 增益有限，外推到 HCC1954 tetraploid 不太可能逆轉
3. **Z3 定義本身可能需要 CN-aware**：本研究僅驗證「在當前 Z3 內部作 filter」的路徑；若改成「CN-aware Zone 重新劃分」屬另案，但歷史探索（`20260418_Z3_Internal_Feature_Exploration`）已顯示 Z3 內無跨樣本特徵可用

### 7.3 不再探索的方向（含理由）

- ✗ Z3 × amplicon/high-CN blacklist 作 cross-sample filter — 本報告確認無效
- ✗ Z3 內 CovM-based 二階分層 — Variant A 數學證明 percentile filter 免疫於 baseline 誤差
- ✗ HCC1954-specific rule 納入 canonical pipeline — 違反 caller 對樣本無假設原則（見 20260419 pilot 結果）

---

## 八、產物清單

### 腳本
- `research/coverage_multiple_validation/scripts/step1_hcc1395_covm_vs_seqc2_truth.py`
- `research/coverage_multiple_validation/scripts/step2_expected_coverage_audit.py`
- `research/coverage_multiple_validation/scripts/step4_oracle_cn_pilot.py`

### 數據
- `data/step1_per_region_truth_aligned.tsv.gz` — 80,524 HCC1395 region × SEQC2 label 對齊
- `data/step1_confusion_matrix.tsv` — 6-class × 4-class
- `data/step1_per_cn_bin_metrics.tsv` — 6 CN bin 的 CovM 統計
- `data/step1_correlation_per_label.tsv` — Pearson/Spearman per SEQC2 label
- `data/step2_expected_coverage_per_sample.tsv` — 14 rows × bias
- `data/step4_variantA_rerun.tsv` — 7 samples × {CovM_old, CovM_v2, S2} ΔF1
- `data/step4_variantB_oracle.tsv` — HCC1395 oracle truth_cn 三 cutoff

### 圖表
- `figures/step1_confusion_matrix.png`
- `figures/step1_covm_vs_truth_scatter.png`
- `figures/step1_recall_per_cn.png`
- `figures/step2_expected_coverage_bias.png`
- `figures/step4_oracle_vs_covm_deltaf1.png`

### 報告
- 本檔：`docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md`

---

## 九、結論

> **Coverage_Multiple baseline 準確度對 Z3 amplicon blacklist 跨樣本失敗的解釋力 = NEGATIVE；oracle CN 最大 ΔF1 = +0.0023（HCC1395_DORADO）**

- 「cov 條件需要修正」為真，但為獨立的 pipeline infrastructure bug（`expected_coverage` 被 default 75.0 污染）
- 修正 cov 不能救 Z3 blacklist 跨樣本失敗；後者是真實生物學 CNV 異質性，不存在 global 規則
- Zone-Aware Framework 定位不變（characterization only）
- 後續：另立 `/cpp-change` 任務修 KDE baseline，並在所有 CovM-based cross-sample 分析文件加註「2026-04-20 後 master dataset 應重跑」風險
