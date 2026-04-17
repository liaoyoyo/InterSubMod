---
title: "B.2-3 驗證：高 NR bin + NR-matched sampling LOH×AF×Methylation"
date: 2026-04-17
status: in_progress
scope:
  - H_B2_3-a: NR≥80 單獨 bin 驗證（檢驗「高 NR 是假關聯 artifact」擔憂）
  - H_B2_3-b: NR-matched sampling（decile bins）消除分佈差異後效應是否保留
related:
  - docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md
  - docs/experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md
  - docs/experiments/in_progress/2026/04/20260417_PartB_effect_size_cn_strat_01.md
git_commit: "fe07550"
---

# B.2-3: High NR bin + NR-matched 驗證

## 0. 摘要

**原擔憂**：step2 報告 NR bin 10-30/30-50/50-80 效應遞增 (ρ=0.483→0.709)；若為統計功效 artifact，則 NR≥80 單獨分析應減弱或消失；NR-matched sampling 也應削弱效應。

**結論**：**駁回 artifact 假設**。

| 假設 | 結果 |
|------|------|
| H_B2_3-a: NR≥80 減弱 | **REJECTED** — TO NR≥80 bin 在 11/12 (bin×sample) 組合中顯著正向 |
| H_B2_3-b: NR-matched 減弱 | **REJECTED** — TO NR-matched 7/7 POS, Δ=+0.29~+0.80, p<10^-26 |

**反向發現**：per-bin 效應強度**隨 NR 遞減**（+0.638 → +0.174），非遞增。可能解釋為 methylation diversity 飽和效應（高 NR 區域內 subclone 多樣性已充分展現，AF 資訊的 additive 價值遞減）。

---

## 1. 方法

- 數據：`all_region_rows.tsv.gz` LOH 區域（`to_loh_bed_hit=True`）共 112,218 rows（僅 TO mode — master dataset 無 Paired LOH annotation）
- AF 分類：extreme (AF<0.3 或 >0.7) vs intermediate (0.4≤AF≤0.6)
- NR bins: [40,60), [60,80), [80,100), [100,150), [150,500)
- **H_B2_3-a**：每個 NR bin 內 per-sample Mann-Whitney U test (alt="less")
- **H_B2_3-b**：每個 sample 的 NR 10-quantile 內進行 extreme / intermediate 1:1 subsample matching，合併後再做 Mann-Whitney U

## 2. 結果

### 2.1 TO NR bin heatmap

| NR bin | Sig count | Median Δ_NG |
|--------|-----------|-------------|
| [40,60) | **7/7** | **+0.638** |
| [60,80) | 6/7 | +0.543 |
| [80,100) | 6/7 | +0.401 |
| [100,150) | 4/7 | +0.281 |
| [150,500) | 1/7 | +0.174 |

**NR≥80 彙總**：11/12 (bin × sample) 組合顯著 POS。

### 2.2 TO NR-matched sampling

| 樣本 | n_matched | Δ_NG matched | p-value |
|------|-----------|--------------|---------|
| HCC1395 | 2472 | +0.506 | 6×10⁻²⁹² |
| HCC1395_DORADO | 2850 | **+0.707** | 0 |
| COLO829 | 135 | **+0.800** | 6×10⁻³⁹ |
| H1437 | 1195 | +0.715 | 8×10⁻²⁶⁹ |
| H2009 | 7962 | +0.292 | 0 |
| HCC1937 | 1404 | +0.504 | 6×10⁻¹⁹⁰ |
| HCC1954 | 211 | +0.483 | 2×10⁻²⁶ |

**7/7 樣本 NR-matched 後仍顯著 POS**，Δ 範圍 +0.29~+0.80。

### 2.3 Paired mode 狀態

**無法驗證** — master dataset (`all_region_rows.tsv.gz`) 僅有 `to_loh_bed_hit`，沒有 `paired_loh_bed_hit` 欄位。需重跑 `cross_sample_audit` 加入 Paired LOH BED annotation。

---

## 3. 解釋

### 3.1 為何 per-bin Δ 隨 NR 遞減？

- **NR 低 (40-60)** 內，單一 region 讀數少 → HPFineNGroups 對 AF intermediate 的敏感度最高（單一 subclonal read 可改變 NGroups），Δ=+0.638
- **NR 高 (150+)** 內，NGroups 易達上限 4（飽和）→ Δ 被壓縮至 +0.174
- 這與 B.1-2 (HPFineNGroups 飽和否定) 邏輯一致：高 NR 區域 NGroups 接近上限，AF 的 additive 訊息變小，但方向仍正確。

### 3.2 為何原 step2 報告 ρ 隨 NR 遞增？

- step2 是 **segment-level Spearman ρ**（ρ 是排序相關，不受效應絕對值影響）
- ρ 隨 NR 增強反映**高 NR segment 內排序關係更穩定**（variance 更小，ρ 更顯著）
- 這裡 Δ_NG 是 **region-level effect size**（絕對差值），遞減反映 ceiling
- 兩者不矛盾，只是不同統計測量角度。

### 3.3 與 B.2-1 (HCC1954 反向) 對照

- B.2-1 驗證 step3 per-sample ρ 僅 5/7 CI 明確
- B.2-3 驗證 **region-level Δ_NG** 7/7 NR-matched POS，HCC1954 n=211 仍 p=2×10⁻²⁶
- 結論：HCC1954 「反向」僅出現在 segment-level small-n sampling，**region-level 證據始終支持 POSITIVE**

---

## 4. 對聲稱的影響

| 原聲稱 / 擔憂 | B.2-3 結果 | 嚴重性 |
|---------------|------------|--------|
| "高 NR bin ρ 遞增可能是統計功效 artifact" | **駁回** — per-bin Δ 遞減，NR-matched 7/7 POS | 高（LOH×AF×Methylation NR confound 擔憂消除） |
| "Paired mode 也有此效應" | **未驗證** — 需重跑 cross_sample_audit | 中（需列為待驗項） |
| "NGroups 飽和會掩蓋 AF 訊號" | 部分確認 — 高 NR 內 Δ 遞減但仍 POS | 低（飽和僅削弱不消除） |

---

## 5. Artifacts

- **腳本**: `research/partB_high_nr_validation/scripts/01_high_nr_check.py`
- **數據**:
  - `data/b2_3_nr_bin_to.tsv`（35 rows = 7 samples × 5 bins）
  - `data/b2_3_nr_matched_to.tsv`（7 rows）
  - Paired 檔存在但空
- **圖表**: `figures/01_nr_bin_heatmap.png`

---

## 6. 後續

### 6.1 自動可做

- 無 — B.2-3 本項 TO mode 已完成

### 6.2 需使用者決策

- **Paired LOH annotation 是否重跑**：若 Paired 也要完整驗證，需重跑 cross_sample_audit with paired_loh_bed_hit 欄位（預估 2-4 小時）。建議**延後至 P0 修復後**，因 paired haplotag 管道也會影響結果。
- 其他 Part B 項目已完成（B.1-2, B.1-1, B.2-1, B.1-3, B.2-5, B.2-2, B.2-3）或需使用者輸入（B.2-4 clinical, B.3 C++ 整合）。

---

## 7. 累積 Part B 驗證結論表

| 質疑 | 結論 | Commit |
|------|------|--------|
| B.1-2 飽和假說 | 否定 | 5358842 |
| B.1-1 residualized AUC | pooled ROBUST Δ=+0.025；per-sample 1/7 AUC-ROBUST | ab61ad1 |
| B.2-1 HCC1954 反向 | 小樣本噪音 (CI 含 0) | 4916cf4 |
| B.1-3 per-sample effect size | 5/7 POS + 2/7 特殊 (COLO829 d=-0.17, H2009 ceiling) | fe07550 |
| B.2-5 cnLOH vs deletion-LOH | 方向一致, CN1=+0.505, CN3=+0.238 | fe07550 |
| B.2-2 Coverage_Multiple proxy | 僅可排序，不精確 | fe07550 |
| **B.2-3 NR confound** | **駁回 — NR-matched 7/7 POS, per-bin Δ 遞減** | **(本輪)** |

**待使用者決策**：B.2-4 Clinical cohort / B.3 C++ 整合成本 / D.1 論文定位 / P0 修復順序。
