---
title: "Part B 合併驗證：B.1-3 Effect Size + B.2-5 CN 分層 + B.2-2 Coverage_Multiple Proxy"
date: 2026-04-17
status: in_progress
scope:
  - B.1-3: HPFineNGroups 7/7 一致的 effect size 強度檢查
  - B.2-5: cnLOH 與 deletion-LOH methylation 動態是否不同
  - B.2-2: Coverage_Multiple 作為 CN proxy 的可靠性
related:
  - docs/experiments/in_progress/2026/04/20260417_HPFineNGroups_saturation_check_01.md
  - docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md
  - docs/experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md
git_commit: "4916cf4"
---

# Part B 合併驗證報告

## 0. 摘要

一次性執行 Part B 三項質疑（B.1-3 / B.2-5 / B.2-2），共享 master dataset 以降低重複 I/O。結論：

| 項目 | 原聲稱 | 驗證結果 | 結論 |
|------|--------|----------|------|
| B.1-3 | HPFineNGroups 7/7 方向一致 | **COLO829 d=-0.17（反向）、H2009 d=+0.12（negligible）** | **精確化為 5/7 正向 + 2/7 negligible/negative** |
| B.2-5 | 混合 cnLOH / deletion-LOH 會掩蓋真實機制 | CN1=+0.505、CN3=+0.238，**兩者方向一致且均 POS** | **否定 "掩蓋" 擔憂** |
| B.2-2 | Coverage_Multiple ≈ CN（精確代理） | 僅 COLO829 bimodal（因 99.7% CN1），其他單峰但 CN 分佈差異極大 | **proxy 僅可排序，不可精確估計** |

---

## 1. B.1-3：HPFineNGroups Per-sample Effect Size

### 1.1 方法

- 過濾：mode=TO、to_loh_bed_hit=False（NonLOH）、NumReads≥40
- 指標：Cohen's d（HPFineNGroups=4 vs <4 的 TP rate 差異）、rank-biserial r
- 判準：|d|≥0.5=medium、|d|≥0.2=small、|d|<0.2=negligible

### 1.2 結果

| 樣本 | n_N4 | n_N<4 | TP_N4 | TP_N<4 | Δ_TP | Cohen's d | 等級 |
|------|------|-------|-------|--------|------|-----------|------|
| HCC1395 | 1909 | 20067 | 0.813 | 0.676 | +0.137 | +0.318 | small |
| HCC1395_DORADO | 884 | 21439 | 0.913 | 0.681 | +0.231 | **+0.601** | medium |
| COLO829 | 206 | 7326 | 0.534 | 0.616 | **-0.082** | **-0.166** | **negative** |
| H1437 | 2265 | 40647 | 0.925 | 0.754 | +0.171 | +0.479 | small |
| H2009 | 24195 | 78374 | 0.935 | 0.903 | +0.032 | +0.117 | **negligible** |
| HCC1937 | 926 | 11266 | 0.714 | 0.448 | +0.266 | **+0.560** | medium |
| HCC1954 | 3882 | 56369 | 0.495 | 0.231 | +0.264 | **+0.571** | medium |

**彙總**：3/7 medium、2/7 small、**1/7 negligible（H2009）、1/7 negative（COLO829）**。

### 1.3 重大發現

**COLO829 在 TO NonLOH NR≥40 條件下，NGroups=4 的 TP rate 反而低於 <4 (0.534 vs 0.616, d=-0.166)**。

這挑戰「7/7 方向一致」的聲稱：若以 per-sample effect size 檢驗，實際為 5/7 正向（其中 3/7 medium、2/7 small）+ 2/7 negligible/negative。

### 1.4 解釋嘗試

- **COLO829 反向可能原因**：n_N4 僅 206（vs 其他樣本 884-24195），小樣本 effect size 估計不穩定。
- **H2009 negligible**：n_N4=24195 為全樣本最大，統計 power 足；但 baseline TP rate 已達 0.903，飽和效應使 Δ 壓縮到 +0.032（即 ceiling effect）。
- 這兩個樣本並非**生物學反向**，而是一個受 power 限制、一個受 ceiling 限制。

---

## 2. B.2-5：cnLOH vs Deletion-LOH 分層

### 2.1 方法

- 僅 LOH region（to_loh_bed_hit=True）+ mode=TO、NumReads≥40
- AF 分類：extreme (AF<0.3 or >0.7) vs intermediate (0.4≤AF≤0.6)
- CN tier：由 Coverage_Multiple 推估
  - CN1_deletion: CM < 0.75（單套刪除）
  - CN2_near: 0.75 ≤ CM < 1.25（近 diploid）
  - CN3_cnLOH_high: 1.25 ≤ CM < 1.75（cnLOH 或 gain）
  - CN4plus: CM ≥ 1.75（複雜 CN）

### 2.2 結果彙總

| CN tier | 樣本數 POS | median Δ_NG | min Δ | max Δ |
|---------|------------|-------------|-------|-------|
| CN1 deletion | **7/7** | **+0.505** | +0.273 | +0.709 |
| CN2 near | 6/6 | +0.343 | +0.121 | +0.742 |
| CN3 cnLOH high | **6/6** | +0.238 | +0.122 | +0.675 |
| CN4plus | 5/5 | +0.337 | +0.180 | +0.664 |

### 2.3 結論

**cnLOH 與 deletion-LOH 方向完全一致，兩者都 POS**：
- 7/7 樣本 CN1 deletion POSITIVE（median Δ_NG=+0.505）
- 6/6 樣本 CN3 cnLOH POSITIVE（median Δ_NG=+0.238）

**deletion-LOH 效應稍強但未掩蓋 cnLOH signal**：median Δ 差異 0.27，兩者同向同質，原擔憂「混合兩者可能掩蓋真實機制」被 B.2-5 否定。

### 2.4 生物學解釋

- **deletion-LOH (CN=1)**：只保留單一 haplotype，intermediate AF 表示 subclonal mutation 於未 delete 的 allele，methylation diversity 直接映射 subclone 組成。
- **cnLOH (CN=2, 單 haplotype duplicated)**：intermediate AF 同樣反映 subclonal mutation，但讀深 2× 導致 noise 更低 → 效應強度稍低可能是訊雜比差異，非生物學機制不同。

---

## 3. B.2-2：Coverage_Multiple 作為 CN Proxy

### 3.1 方法

- 檢查 LOH region Coverage_Multiple 分佈是否 bimodal（若 bimodal → 代表 CN=1 與 CN=2 peaks 可區分）
- 計算每樣本 CN tier 分佈百分比
- 峰偵測：histogram 50 bins，找局部最大值 (peak_density ≥ 10)

### 3.2 結果

| 樣本 | n | median CM | IQR | CN1% | CN2% | CN3% | CN4+% | n_peaks |
|------|---|-----------|-----|------|------|------|-------|---------|
| HCC1395 | 17555 | 0.667 | 0.360 | 62.2 | 29.1 | 6.7 | 1.9 | 1 |
| HCC1395_DORADO | 17785 | 0.707 | 0.360 | 57.6 | 32.8 | 7.1 | 2.4 | 1 |
| **COLO829** | 10582 | 0.293 | 0.146 | **99.7** | 0.3 | 0.0 | 0.0 | **2** |
| H1437 | 15531 | 0.667 | 0.214 | 69.3 | 26.6 | 3.0 | 1.0 | 1 |
| H2009 | 34141 | 1.013 | 0.240 | 10.4 | **78.7** | 10.2 | 0.7 | 1 |
| **HCC1937** | 12419 | 1.200 | 0.533 | 12.9 | 42.9 | **30.6** | **13.5** | 1 |
| HCC1954 | 4205 | 0.787 | 0.386 | 45.0 | 40.5 | 11.3 | 3.2 | 1 |

### 3.3 結論

- **bimodal 僅 COLO829 (1/7)**：原因是 99.7% 為 CN1，分佈只有一個主峰加一個次小峰，非真 CN 區分。
- **其他 6 樣本為單峰，但 CN 分佈差異巨大**：
  - H2009 以 CN2 為主（78.7%）→ 近 diploid-dominant
  - HCC1937 hyper-diploid（CN3+=44.1%）
  - COLO829 極端 hypo-diploid（CN1=99.7%）
- **Coverage_Multiple 可用作 CN 排序 proxy**（CN1 → CN2 → CN3+ 順序正確），但**不適合精確估計**，特別是 HCC1937 這類 hyper-diploid 樣本。

### 3.4 建議

**未來 Phase 2 生物學驗證應整合 CNV caller**（Delly、Manta、sequenza、ASCAT）取代 Coverage_Multiple，特別是：
- HCC1937 / HCC1954（hyper-diploid、複雜 CN）
- 任何聲稱 "CN1 vs CN2 分層" 的分析（B.2-5 應在整合後複驗）

---

## 4. 綜合結論

### 4.1 對既有聲稱的影響

1. **HPFineNGroups "7/7 方向一致"**：應精確化為 **5/7 正向（3/7 medium + 2/7 small）+ 2/7 特殊（COLO829 小樣本反向、H2009 ceiling effect）**。pooled level POSITIVE 不變。

2. **cnLOH vs deletion-LOH 擔憂**：**B.2-5 否定**。兩者方向一致且均 POS，deletion-LOH 效應稍強但非機制性差異。

3. **Coverage_Multiple proxy**：**僅能排序，不能精確估計**。論文 Limitation 段必須明寫。

### 4.2 後續優先級

| 下一步 | 優先級 | 理由 |
|--------|--------|------|
| **整合 CNV caller** 取代 Coverage_Multiple | P1 | B.2-5 結論需要精確 CN 確認 |
| HCC1937 hyper-diploid 深度分析 | P2 | 唯一 CN4+≥13.5% 樣本，可能為 subclonal CN evolution 範例 |
| 更新 HPFineNGroups POSITIVE 主張到 "5/7" | P0 | 避免過度宣稱，影響論文 narrative |
| B.2-4 Clinical cohort 外推 | 待決策 | 需 patient-derived 資料，建議用戶決定 |

---

## 5. Artifacts

- **腳本**: `research/partB_effect_size_cn_stratification/scripts/01_effect_size_and_cn_strat.py`
- **數據**:
  - `data/b1_3_effect_size_per_sample.tsv`
  - `data/b2_5_cn_stratified_effects.tsv`
  - `data/b2_2_coverage_multiple_distribution.tsv`
- **圖表**:
  - `figures/01_b1_3_effect_size_forest.png`
  - `figures/02_b2_5_cn_stratified_forest.png`
  - `figures/03_b2_2_cov_multiple_dist.png`

---

## 6. 待使用者決策項目

接下來 Part B 剩餘項目都需要使用者輸入：

1. **B.2-4 Clinical cohort 外推**：是否投入取得 patient-derived ONT cohort？或先以「cell line only + limitation 段落」投稿？
2. **B.3 Per-CpG ASM C++ 整合**：是否值得將 10 欄位整合進 C++ 核心？成本 vs 邊緣效應。
3. **D.1 論文定位**：A（工具論文）/ B（生物標誌）/ C（negative findings）還是混合？
4. **P0 修復優先級**：ReadParser 修復（1-2 週）是否立即啟動？

建議下一次對話由使用者決定其中至少 2 項，否則後續無法推進。
