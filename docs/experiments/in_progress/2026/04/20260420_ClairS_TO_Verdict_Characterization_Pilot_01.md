<!--
建立時間: 2026-04-20
更新時間: 2026-04-20
狀態: validated
判定: NEGATIVE（F1 gain path）/ POSITIVE（Verdict 內部校準）
scope: HCC1395 subsample t20_n30 · Characterization only · 不改 C++ · 不整合 filter
資料來源:
  - research/clairs_to_verdict_pilot/00_INDEX.md
  - /big8_disk/data/HCC1395/ONT/subsample/t20_n30/ClairS_TO_v0_3_0/snv.vcf.gz
  - /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
  - /big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed
  - docs/references/2026/04/20260420_external_CN_tools_survey_01.md
related:
  - 20260419_Z3_amplicon_blacklist_pilot_result_01.md
  - 20260420_CovM_Baseline_Accuracy_Validation_01.md
-->

# ClairS-TO Verdict 模組 Characterization Pilot — 結果報告

## 一、執行摘要

**最終判定（兩層次）**：

- **Verdict 內部校準 POSITIVE**：在 HCC1395 t20_n30（purity ≈ 0.40）資料分佈下，Verdict 標籤與 SEQC2 truth 的符合度非常好：
  - H-V1 Verdict_Germline FP rate = **96.1%**（3,463/3,602, Wilson 95% CI 0.955–0.967） — 遠超 70% gate
  - H-V2 Verdict_Somatic (on PASS) TP rate = **91.8%**（7,785/8,479） — 超過 85% gate
  - Verdict_SubclonalSomatic TP rate = **94.9%**（517/545）
- **F1 直接增益路徑 NEGATIVE**：
  - Verdict_Germline **全部** 4,633 個標籤只落在 ClairS-TO `LowQual` 變異上，PASS 變異中無 Verdict_Germline
  - 本計劃的「remove Verdict_Germline from PASS」假設性 filter → **ΔF1 = +0.0000**
  - 最激進的「只保留 PASS ∩ Verdict_Somatic∪Subclonal」 → **ΔF1 = −0.2007**（recall 從 0.521 崩到 0.210）
  - 結論：啟用 Verdict 不會改善本案主 pipeline 的 F1，因 Verdict 與 LowQual 決策共用同一套 ASCAT/binomial 資訊
- **對五研究目標影響**：
  - 目標 5（F1）— NEGATIVE；無新增增益路徑
  - 目標 1/2/3/4 — 零衝擊；Verdict 標籤可保留作 per-variant CN-aware annotation 供後續觀察

**後續動作**：依計劃進入 **Wakhan / SAVANA external CN tools pilot**。詳見 `docs/references/2026/04/20260420_external_CN_tools_survey_01.md`。

---

## 二、動機與假說

### 2.1 觸發鏈

- **2026-04-19**：Z3 × HCC1954 amplicon blacklist pilot 確認 region-level filter 跨樣本不可行（HCC1954-only CONDITIONAL）
- **文獻調查**：ClairS-TO v0.2.0+ 自帶 **Verdict 模組**（ASCAT CN segmentation → purity/ploidy 估計 → 二項式檢定 → per-variant germline/somatic/subclonal 分類）— Nat Commun 2025 · [HKU-BAL/ClairS-TO verdict.md](https://github.com/HKU-BAL/ClairS-TO/blob/main/docs/verdict.md)
- **現況盤點**：本專案 v0.3.0 主 pipeline 輸出 VCF 中 Verdict flags 全為空：
  - 6 樣本（COLO829/H1437/H2009/HCC1937/HCC1954）：Docker 容器缺 `clairs-to_cna_data` 1000G loci → ASCAT 失敗
  - HCC1395：ASCAT 成功（purity = 0.91）但 v0.3.0 在 purity > 0.8 時跳過 Verdict tagging
  - **唯獨 HCC1395 t20_n30 subsample**（purity = 0.40）實際產出 14,875 tags → 現成測試材料

### 2.2 假說（Plan 宣告）

| ID | 假說 | POSITIVE 門檻 | Null 門檻 |
|----|------|-------------|----------|
| H-V1 | Verdict_Germline 富集於 SEQC2 FP | FP rate ≥ 70% | < 55% |
| H-V2 | Verdict_Somatic 與 SEQC2 TP 一致率 | TP rate ≥ 85% | FP rate > 20% |
| H-V3 | 若 H-V1 POSITIVE 可推廣至其他 6 樣本 | — | — |

### 2.3 GO 條件

H-V1 + H-V2 同時 POSITIVE **且** 假設性 filter ΔF1 ≥ +0.005 → 進入「7 樣本 loci 補齊 + safety gate 繞開」任務。
否則判定 NEGATIVE，升級 Wakhan / SAVANA external CN 路徑並保留 Verdict 作 annotation。

---

## 三、方法

### 3.1 資料

- **ClairS-TO VCF**：`/big8_disk/data/HCC1395/ONT/subsample/t20_n30/ClairS_TO_v0_3_0/snv.vcf.gz`（88 MB，ClairS-TO v0.3.0，平台 `ont_r10_dorado_sup_5khz_ssrs`）
- **SEQC2 sSNV truth**：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`（39,447 sSNV）
- **SEQC2 HighConf BED**：`/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`（22 chroms, 709,209 intervals）

### 3.2 腳本

兩個獨立純 Python 分析腳本（pysam + pandas + matplotlib）：

- `research/clairs_to_verdict_pilot/scripts/step1_verdict_vs_seqc2.py`：每個 sSNV 記錄 × {chrom, pos, ref, alt, FILTER, Verdict 類別, in_highconf_BED, SEQC2 label (TP/FP), AF, DP}；輸出 confusion matrix × Wilson CI
- `research/clairs_to_verdict_pilot/scripts/step2_reference_f1.py`：列舉 4 個假設性 filter 場景（S0 baseline / S1 remove Verdict_Germline from PASS / S2 只保留 PASS ∩ Verdict_Somatic∪Subclonal / S3 加 LowQual ∩ Verdict_Somatic）並計算 F1 + ΔF1

### 3.3 口徑

- 僅分析 ClairS-TO VCF 中落在 **SEQC2 HighConf BED** 的 sSNV 記錄（truth 集只在此範圍完整）
- 僅考慮 single-base ref/alt（InDel 不在 SEQC2 sSNV truth 中）
- TP/FP 透過 (chrom, pos, ref, alt) key 在 SEQC2 truth 中 membership test 判定
- SEQC2 FILTER 接受 `PASS` 和 `HighConf`（共 39,447 筆）

---

## 四、Step 1 — Verdict × SEQC2 交叉分析結果

### 4.1 預 flight：Verdict × ClairS-TO FILTER 分佈

在解讀結果前先註記一個**對後續判定至關重要**的觀察：

| Verdict class | 數量 | ClairS-TO FILTER 分佈 |
|---------------|------|----------------------|
| Verdict_Germline | 4,633 | **100% LowQual**（無 PASS） |
| Verdict_Somatic | 9,602 | 100% PASS |
| Verdict_SubclonalSomatic | 640 | 100% PASS |

→ 這表明：ClairS-TO v0.3.0 的 Verdict 判定**已經前置嵌入**在 `PASS` vs `LowQual` 的判定裡。Verdict_Germline 永遠不會出現在 PASS 集合上。本 pilot 的「移除 Verdict_Germline from PASS」filter 在 subsample t20_n30 的預期 ΔF1 = 0。

### 4.2 主表：Verdict × SEQC2 confusion matrix（in HighConf BED）

![Verdict FP rate per class](../../../../../research/clairs_to_verdict_pilot/figures/step1_verdict_fp_rate_per_class.png)

#### 全變異（ALL：PASS + LowQual + 其他 filter）

| Verdict | n | TP | FP | TP rate | FP rate | FP rate 95% Wilson CI |
|---------|---|----|----|---------|---------|----------------------|
| Verdict_Germline | 3,602 | 139 | 3,463 | 0.039 | **0.961** | 0.955 – 0.967 |
| Verdict_Somatic | 8,479 | 7,785 | 694 | **0.918** | 0.082 | 0.076 – 0.088 |
| Verdict_SubclonalSomatic | 545 | 517 | 28 | **0.949** | 0.051 | 0.036 – 0.073 |
| None (no Verdict) | 3,329,765 | 14,439 | 3,315,326 | 0.004 | 0.996 | 0.996 – 0.996 |

#### PASS only

| Verdict | n | TP | FP | TP rate | FP rate |
|---------|---|----|----|---------|---------|
| Verdict_Somatic | 8,479 | 7,785 | 694 | **0.918** | 0.082 |
| Verdict_SubclonalSomatic | 545 | 517 | 28 | **0.949** | 0.051 |
| None (no Verdict on PASS) | 27,125 | 12,233 | 14,892 | 0.451 | **0.549** |

#### LowQual only

| Verdict | n | TP | FP | TP rate | FP rate |
|---------|---|----|----|---------|---------|
| Verdict_Germline | 3,602 | 139 | 3,463 | 0.039 | 0.961 |
| None (no Verdict on LowQual) | 32,576 | 1,540 | 31,036 | 0.047 | 0.953 |

### 4.3 假說判定

| 假說 | 觀察值 | 門檻 | 判定 |
|------|-------|------|------|
| **H-V1** Verdict_Germline FP rate | 0.961 | ≥ 0.70 | ✅ **POSITIVE** |
| **H-V2** Verdict_Somatic TP rate (PASS) | 0.918 | ≥ 0.85 | ✅ **POSITIVE** |
| 額外觀察：Verdict_SubclonalSomatic TP rate | 0.949 | — | 額外 POSITIVE |

**結論**：**Verdict 模組在此資料分佈下的內部校準優秀** — 對 purity = 0.40 的低純度 subsample，Verdict 可以正確區分 germline 殘留 vs 體細胞變異。

### 4.4 關鍵觀察 — Verdict 與 LowQual 資訊重疊

對比 LowQual only subset：

- LowQual ∩ Verdict_Germline FP rate = **0.961** （n = 3,602）
- LowQual ∩ no_Verdict FP rate = **0.953** （n = 32,576）
- Δ ≈ 0.008 絕對差 → Verdict_Germline **僅比 LowQual 本身的 FP 排除力多 0.8%**

這指出 Verdict 與 LowQual 在**資訊層面高度重疊**：Verdict 的校準正確性並不能單獨為 PASS 集合帶來新的 FP 排除力。

---

## 五、Step 2 — 假設性 F1 計算

### 5.1 Scenario 表

![Reference F1 barchart](../../../../../research/clairs_to_verdict_pilot/figures/step2_reference_f1_barchart.png)

| Scenario | n | TP | FP | Precision | Recall | F1 | ΔF1 vs baseline |
|----------|---|----|----|-----------|--------|------|----------------|
| S0 baseline (PASS as-is) | 36,149 | 20,535 | 15,614 | 0.568 | 0.521 | **0.5433** | 0 |
| S1 PASS − Verdict_Germline | 36,149 | 20,535 | 15,614 | 0.568 | 0.521 | 0.5433 | **+0.0000** |
| S2 PASS ∩ (Somatic ∪ Subclonal) | 9,024 | 8,302 | 722 | 0.920 | 0.210 | 0.3426 | **−0.2007** |
| S3 PASS ∪ (LowQual ∩ Verdict_Somatic) | 36,149 | 20,535 | 15,614 | 0.568 | 0.521 | 0.5433 | +0.0000 |

**ΔF1 gate**：

- Plan 門檻：ΔF1 ≥ +0.01 強 POSITIVE；+0.005 ≤ Δ < +0.01 CONDITIONAL；< +0.005 NEGATIVE
- S1（nominal）= 0 → **NEGATIVE**
- S2（aggressive）= −0.2 → 因 recall 崩塌，無法採用

### 5.2 Verdict 為什麼救不了 PASS 的 FP

PASS 集合中 FP 的主體是 **PASS ∩ no_Verdict** 子集：

- n = 27,125 records, FP = 14,892（PASS 總 FP 15,614 的 **95.4%**）
- Verdict 根本沒有在這 27,125 筆上產生標籤 → Verdict 對 ClairS-TO PASS 主要 FP 問題**沒有決策權**

這是結構性限制：Verdict 的 binomial 檢定只在 ASCAT segment call 且 CN estimate 可信的變異上活躍。PASS 內**沒被 Verdict 挑中的**這 27,125 筆才是 FP 大戶，Verdict 對它們不提供任何訊號。

---

## 六、Step 3 — ISM region zone 交叉（skipped）

Plan 規定 Step 3 僅在 Step 1/2 至少一項 POSITIVE 時執行。本 pilot Step 1 H-V1/V2 POSITIVE，因此 Step 3 條件成立。

**實際狀態**：**skipped**。

- 檢查 `/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/HCC1395/purity_t20_n30_20260213_011454/` 目錄為空，**t20_n30 subsample 未產出對應 ISM region-level 分析**
- 主 canonical ISM baseline（`output/canonical/HCC1395/*`）使用的是原始高純度 BAM，不是 t20_n30 subsample
- 若強行改用原始 HCC1395 的 ISM，又會失去 Verdict tagging（因高純度下 ASCAT 被跳過）→ 口徑不一致

**Limitation 欄**：本 pilot 無法在 t20_n30 資料下建立「Verdict × Zone × HPFineNGroups」三維交叉。

此 limitation **不改變主結論**：即使 Step 3 能執行，Step 2 已證 Verdict 對 PASS FP 無決策權，region-level zone 交叉只會重述這一事實。

---

## 七、最終判定

### 7.1 主結論

**F1 路徑 NEGATIVE · Verdict 校準 POSITIVE**

- 對現行 ClairS-TO v0.3.0 pipeline 的 PASS 輸出，**啟用 Verdict 不會改善本案 somatic calling F1**（ΔF1 = 0）
- Verdict 的內部校準是正確的：在 purity ≈ 0.40 的資料分佈下，Verdict_Germline 96.1% 是真 germline、Verdict_Somatic 91.8% 是真 TP
- 但這份校準資訊已經嵌入在 LowQual / PASS 的決策裡，沒有可以從 PASS 再榨出的增益

### 7.2 對 5 研究目標的影響

| 目標 | 影響 | 說明 |
|------|------|------|
| 1. per-CpG ASM | 零 | 本 pilot 不涉及甲基化層 |
| 2. Clone 結構 | 零 | Verdict 不提供 subclonal 結構信息 |
| 3. 二次打擊順序 | 零 | 無新訊號 |
| 4. Biomarker | 零 | Verdict 僅是 per-variant 分類，不生成 biomarker 特徵 |
| 5. F1 優化 | **NEGATIVE** | 無路徑；改以外部 CN 工具補正 |

### 7.3 對 Zone-Aware Framework 定位的影響

**不變**。Verdict 未提供跨樣本可用的新機制，Zone-Aware Framework 定位維持 characterization annotation only。

### 7.4 Verdict 標籤保留計劃

依計劃，Verdict 標籤作 **per-variant annotation 永久保留**：

- subsample t20_n30 的 14,875 tags 已存在 VCF INFO 欄位，未來可隨 VCF 被 ISM 讀入
- 當 Wakhan/SAVANA external CN pilot 啟動後，可作為 **independent CN caller 交叉驗證錨點**（Verdict-Somatic ∩ Wakhan-Neutral 區域是最高信心子集）
- 不整合為 filter；不在主 pipeline 加 ClairS-TO `--enable_verdict`（v0.3.0 預設已啟用，本 pilot 不改 Docker 調用）

---

## 八、後續行動（與 Plan 的跨文件銜接）

### 8.1 升級路徑：Wakhan / SAVANA External CN Pilot

詳見 `docs/references/2026/04/20260420_external_CN_tools_survey_01.md`。

建議優先順序：

1. **Wakhan**（haplotype-specific CN，與 HP-tagged ISM 架構對應最佳）單樣本 HCC1395 pilot（2-3 天）
2. 若 Wakhan 在 HCC1395 通過 → 橫擴 6 in-house 樣本
3. SAVANA 作為 SV 互補層（可選 parallel）
4. 整合策略：外部 CN 以 BED 形式進入 RegionProcessor external annotation；不改 ClairS-TO Verdict path

### 8.2 不做事項

- **不補齊 ClairS-TO Docker 的 1000G loci**：本 pilot Step 2 已量化此路徑無 F1 增益
- **不改 ClairS-TO v0.3.0 的 purity > 0.8 Verdict safety gate**：改了也不會救 PASS 的 no_Verdict FP 子集（95.4% 的 FP 就在這裡）
- **不把 Verdict 當 filter 納入 ClairS-TO 後處理**：為 0 的 ΔF1 不值得工程成本

### 8.3 記憶條目

- 已新增：`memory/project_clairs_to_verdict_pilot.md`（本 pilot 主要發現 + 判定）
- 已更新：`docs/CURRENT_FOCUS.md` Zone-Aware Framework 區塊
- 已新增：`docs/references/2026/04/20260420_external_CN_tools_survey_01.md`

---

## 九、產物清單

### 腳本
- `research/clairs_to_verdict_pilot/scripts/step1_verdict_vs_seqc2.py`
- `research/clairs_to_verdict_pilot/scripts/step2_reference_f1.py`

### 數據
- `research/clairs_to_verdict_pilot/data/step1_verdict_vs_seqc2.tsv`（3,831,681 records）
- `research/clairs_to_verdict_pilot/data/step1_confusion_matrix.tsv`（3 subsets × 4 Verdict classes）
- `research/clairs_to_verdict_pilot/data/step2_reference_f1.tsv`（4 scenarios）

### 圖表
- `research/clairs_to_verdict_pilot/figures/step1_verdict_fp_rate_per_class.png`
- `research/clairs_to_verdict_pilot/figures/step2_reference_f1_barchart.png`

### 文件
- `research/clairs_to_verdict_pilot/00_INDEX.md`
- 本報告（`docs/experiments/in_progress/2026/04/20260420_ClairS_TO_Verdict_Characterization_Pilot_01.md`）
- `docs/references/2026/04/20260420_external_CN_tools_survey_01.md`

---

## 十、結論速報

```
ClairS-TO Verdict 對本案 TO somatic calling F1 的增益解釋力 = NEGATIVE (F1 path)
                                           + POSITIVE (internal calibration)
Verdict_Germline FP rate = 96.1%（但 100% 落在 LowQual，不在 PASS）
參考 ΔF1 = +0.0000（S1 nominal filter）/ −0.2007（S2 aggressive）
→ 後續動作：升級 Wakhan/SAVANA external CN 路徑
              Verdict 標籤永久保留作 ISM region annotation
```
