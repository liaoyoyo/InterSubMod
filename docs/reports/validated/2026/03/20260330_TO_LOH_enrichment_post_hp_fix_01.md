<!--
建立時間: 2026-03-30 07:00
目標: HP integer tag 修正後，重新評估 TO LOH enrichment 特性（COLO829、H1437、H2009 為重點）
處理範圍: 7個 TO 樣本，修正後 significance_summary.csv（2026-03-30），LOH eligibility + enrichment 分析
關聯檔案:
  - src/core/ReadParser.cpp (HP integer tag fix)
  - scripts/analysis/build_to_loh_enrichment_post_hp_fix.py
  - docs/experiments/in_progress/2026/03/20260330_TO_HP_integer_tag_fix與全樣本重跑結果確認_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_enrichment_post_hp_fix/
-->

# TO LOH Enrichment 重新評估報告（HP Integer Tag Fix 後）

**分析日期**：2026-03-30
**數據版本**：HP fix 後重跑（所有 TO significance_summary.csv 更新時間：2026-03-30 03:53 ~ 06:21）
**前提條件**：2026-03-30 之前的所有 TO HP/LOH 統計數據均無效（HP:i:11/21/33 mapping bug）
**分析腳本**：`scripts/analysis/build_to_loh_enrichment_post_hp_fix.py`
**輸出目錄**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_to_loh_enrichment_post_hp_fix/`

---

## 執行摘要

HP integer tag 修正後，TO LOH 分析首次可以基於完整、正確的 read phasing 數據。本報告重點回答三個問題：

1. **LOH eligibility 修復程度**：修正後各樣本有多少 region 真正具備足夠的 phasing 支持（eff_hp ≥ 30）？
2. **TO LOH enrichment 的真實方向**：修正後 LOH 是 FP 指標還是 TP 指標？
3. **COLO829/H1437/H2009 的特殊觀察**：三個受影響最嚴重的樣本有何新發現？

**核心結論**：修正後 TO LOH enrichment **全面呈現 TP 富集方向**（enrichment < 1），與 paired 端方向相反。TO LOH 不是 FP filter，而是弱 TP 支持訊號（0.77～0.96×）。

---

## 1. LOH Eligibility 修復（eff_hp ≥ 30）

### Table 1. LOH Eligible Region 修正前後比較（TP）

| 樣本 | 修正前 LOH Elig% | 修正後 LOH Elig% | Δ (pp) | fold change |
|------|----------------|----------------|--------|-------------|
| H1437 | 24.5% | **95.5%** | +71.0pp | +3.9× |
| COLO829 | 0.7% | **34.7%** | +34.0pp | +49.5× |
| HCC1395_DORADO | 46.9% | **93.1%** | +46.2pp | +2.0× |
| H2009 | 64.5% | **97.7%** | +33.2pp | +1.5× |
| HCC1395 | 58.5% | **92.2%** | +33.7pp | +1.6× |
| HCC1954 | 70.0% | **93.0%** | +23.0pp | +1.3× |
| HCC1937 | 86.7% | **97.4%** | +10.7pp | +1.1× |

**關鍵觀察**：
- H1437 最顯著：69% 的 reads 原本不可見，修正後 95.5% regions 達到 LOH eligible
- COLO829 雖然 fold change 最大（49.5×），但最終只達 34.7%（低於其他樣本的 92-98%）
  - 原因：COLO829 的 median eff_hp 仍只有 25（Tier B），說明 COLO829 本身 phasing coverage 較低（可能為生物特性）
- 修正後所有樣本 LOH eligible 均提升，確認 bug 已完全解決

---

## 2. TO LOH Enrichment：全域聚合（7個樣本）

### Table 2. TO LOH Enrichment by eff_hp Tier（全樣本聚合，修正後）

| Tier | eff_hp 範圍 | TP total | TP LOH% | FP total | FP LOH% | FP/TP Enrich | 解讀 |
|------|-----------|---------|---------|---------|---------|------------|------|
| C0 | =0 | 732 | 0.0% | 338 | 0.0% | N/A | 無 phasing |
| C | <10 | 4,084 | 76.1% | 2,027 | 77.3% | 1.015× | 無鑑別力 |
| B | 10-29 | 27,422 | 55.8% | 14,506 | 52.4% | **0.938×** | 微弱 TP 富集 |
| **A** | **30-49** | **44,177** | **61.4%** | **25,255** | **43.3%** | **0.706×** | **TP 富集（最強）** |
| A+ | ≥50 | 214,895 | 39.1% | 86,256 | 30.0% | **0.766×** | TP 富集 |
| A≥30 | ≥30 | 259,072 | 42.9% | 111,511 | 33.0% | **0.769×** | TP 富集（p≈0） |
| All | 全部 | 291,310 | 44.5% | 128,382 | 35.8% | **0.805×** | TP 富集 |

**重要發現**：

1. **TO LOH 全面呈現 TP 富集**（enrichment < 1）：與 paired 端的 FP 富集方向完全相反
2. **Tier A（30-49）是 TP 富集最強的 tier**（0.706×）：eff_hp 30-49 範圍的 LOH-like region 明顯更傾向是 TP
3. **Tier C（<10 eff_hp）基本無鑑別力**（1.015×）：phasing 不足時 LOH 訊號是噪音
4. **A≥30 全域統計**：0.769×，Fisher p ≈ 0，統計極顯著

---

## 3. Per-Sample 分析（重點樣本）

### Table 3. TO LOH Enrichment by Sample × Key Tiers

#### All Tier（無門檻）

| 樣本 | TP LOH% | FP LOH% | Enrich | p-value | median eff_hp (TP) |
|------|---------|---------|--------|---------|-------------------|
| HCC1395 5kHz | 60.6% | 54.3% | 0.896× | p≈0 | 61 |
| HCC1395 DORADO | 61.6% | 56.0% | 0.909× | p≈0 | 63 |
| COLO829 | 35.4% | 33.8% | 0.956× | 4.5e-4 | **25** |
| H1437 | 41.8% | 38.4% | 0.919× | p≈0 | 66 |
| H2009 | 40.9% | 37.8% | 0.923× | p≈0 | 85 |
| HCC1937 | 64.8% | 57.2% | 0.882× | p≈0 | 99 |
| HCC1954 | 25.0% | 21.3% | 0.852× | p≈0 | 61 |

#### Tier A（eff_hp 30-49）

| 樣本 | TP n | TP LOH% | FP n | FP LOH% | Enrich | p-value |
|------|------|---------|------|---------|--------|---------|
| HCC1395 5kHz | 7,277 | 84.4% | 2,867 | 82.9% | 0.982× | 0.062 |
| HCC1395 DORADO | 6,784 | 86.7% | 2,680 | 84.9% | 0.979× | 0.022 |
| **COLO829** | **10,881** | **15.6%** | **6,066** | **14.8%** | **0.951×** | **0.189（不顯著）** |
| H1437 | 8,635 | 79.8% | 2,586 | 75.4% | 0.945× | 3.0e-6 |
| H2009 | 6,431 | 76.6% | 545 | 74.5% | 0.973× | 0.270（不顯著） |
| HCC1937 | 645 | 91.3% | 629 | 85.2% | 0.933× | 8.8e-4 |
| HCC1954 | 3,524 | 28.0% | 9,882 | 25.3% | 0.904× | 0.002 |

#### Tier A+（eff_hp ≥50）

| 樣本 | TP n | TP LOH% | FP n | FP LOH% | Enrich | p-value |
|------|------|---------|------|---------|--------|---------|
| HCC1395 5kHz | 18,999 | 48.0% | 7,820 | 39.9% | 0.832× | p≈0 |
| HCC1395 DORADO | 20,073 | 50.8% | 8,371 | 45.4% | 0.892× | p≈0 |
| COLO829 | **594** | 15.3% | **510** | 16.3% | **1.062×** | 0.679（不顯著，樣本極少）|
| H1437 | 34,805 | 30.4% | 10,154 | 26.5% | 0.872× | p≈0 |
| H2009 | 116,427 | 38.0% | 11,189 | 35.0% | 0.922× | p≈0 |
| HCC1937 | 11,645 | 62.8% | 11,010 | 54.6% | 0.868× | p≈0 |
| HCC1954 | 12,352 | 20.7% | 37,202 | 16.8% | 0.808× | p≈0 |

---

## 4. COLO829、H1437、H2009 深度觀察

### 4.1 COLO829

**修正前**：99.3% TP regions 的 effective_hp = 0（完全不可見），LOH 統計完全無效。
**修正後**：

- LOH eligible（eff_hp ≥ 30）：34.7%（仍低於其他樣本的 92-98%）
- Median eff_hp（TP）= 25，落在 Tier B（10-29）
- 整體 enrichment = 0.956×（微弱 TP 富集，p=4.5e-4）
- Tier A（30-49）enrichment = 0.951×，**p=0.189（不顯著）**
- Tier A+ enrichment = 1.062×，**p=0.679（不顯著，樣本只有 594+510 個）**

**結論**：COLO829 TO LOH 訊號微弱，統計力不足。兩個主要原因：
1. Phasing coverage 本身較低（median eff_hp = 25，多數 region 在 Tier B）
2. COLO829 的 LOH-like 比例偏低（35%），可能是樣本生物特性（Melanoma 基因組異質性較高）

**行動**：COLO829 的 TO LOH feature 在 Phase 1 應標記為「低信心」，暫不作為 feature engineering 主材料。

### 4.2 H1437

**修正前**：69% reads 不可見，LOH eligible 僅 24.5%，分析完全扭曲。
**修正後**：

- LOH eligible 飆升至 **95.5%**（最大改善幅度）
- Median eff_hp（TP）= 66，主要落在 Tier A+
- 整體 enrichment = 0.919×（TP 富集，p≈0）
- Tier A（30-49）enrichment = **0.945×**，p=3e-6
- Tier A+ enrichment = **0.872×**，p≈0

**結論**：H1437 TO LOH 訊號清晰且統計顯著。LOH 在 TP 中更常見（尤其是 A+ tier），可以作為可信的 TP 富集指標。這是修正後最明確的新發現之一。

### 4.3 H2009

**修正前**：51% reads 不可見，LOH eligible 64.5%，已經偏低。
**修正後**：

- LOH eligible 達到 **97.7%**（接近天花板）
- Median eff_hp（TP）= 85，主要 Tier A+
- 整體 enrichment = 0.923×（TP 富集，p≈0）
- Tier A（30-49）enrichment = 0.973×，**p=0.270（不顯著，因為 FP 數量很少 n=545）**
- Tier A+ enrichment = **0.922×**，p≈0

**結論**：H2009 是高 eff_hp 樣本，TO LOH 訊號穩定（0.92×），但統計力在 Tier A 受限於 FP 數量少（只有 545 個）。整體方向與 H1437 一致，都是 TP 富集。

---

## 5. 與 Paired LOH Enrichment 的比較

| 樣本 | TO All（新）| TO All（舊 Round 1）| Paired（Round 1）|
|------|------------|---------------------|------------------|
| COLO829 | 0.956× | 1.010× | **1.155×** |
| H1437 | 0.919× | 1.017× | 1.795×（FP 過少，僅 8 個）|
| H2009 | 0.923× | N/A | **2.685×** |
| HCC1395 | 0.896× | 0.912×（Round 1 global）| 1.194×（global）|
| HCC1954 | 0.852× | N/A | N/A |

**關鍵對比**：

1. **TO 和 Paired 的 LOH enrichment 方向相反**：
   - Paired：LOH 傾向 FP enrichment（≥1），特別是 H2009（2.685×）
   - TO：LOH 傾向 TP enrichment（<1），所有樣本一致

2. **這個差異有生物學意義**：
   - Paired 模式中，高 eff_hp LOH region 對應「somatic caller 在 LOH block 中誤認 germline SNP」的 FP 機制
   - TO 模式中，LOH region 的 TP 富集反映「tumor-only 呼叫的真實 somatic SNV 更集中在 allele-imbalanced 區域」

3. **Round 1 TO 舊數據（1.010×）是無效的**：
   - 舊數據 median eff_hp_tp = 5.0（COLO829），幾乎全部是 C0/C tier
   - 舊版 LOH-like 比例 67.9%（COLO829）是基於錯誤 HP_Ratio 計算的假訊號
   - 修正後 LOH-like 比例 35.4%，更接近生物學真實值

---

## 6. Phase 1 Feature Engineering 含義

### 修正後 TO LOH Feature 設計建議

| Feature | 計算方式 | 方向（FP risk）| 強度 | 適用樣本 | 備注 |
|---------|---------|--------------|------|---------|------|
| `to_tier_a_loh_like` | eff_hp 30-49 AND Potential_LOH | **-FP risk（TP 富集）** | 0.706× (agg) | H1437 \| HCC1937 \| HCC1395 | Tier A 最強 |
| `to_tier_a_plus_loh_like` | eff_hp ≥50 AND Potential_LOH | **-FP risk（TP 富集）** | 0.766-0.832× | H1437 \| H2009 \| HCC1395 \| HCC1937 \| HCC1954 | 樣本多，統計最穩 |
| `to_loh_eligible_confidence` | eff_hp ≥30（LOH eligible flag）| 純 region quality QC | — | 全樣本 | 不直接作為 enrichment feature |
| ~~`to_loh_fp_marker`~~ | ~~（不存在）~~| ~~（TO LOH 不是 FP marker）~~ | — | — | Round 1 前的假設已推翻 |

### 重要修正

原本 Round 1 的 TO LOH 結論（0.912× global，以 buggy 數據計算）方向雖偶然接近正確（<1），但數值不可信：
- 舊版 COLO829 TO：1.010×（完全無鑑別力）→ 修正後 0.956×（微弱，不顯著）
- 舊版的「全域 0.912×」主要由 HCC1395 帶動，不代表 COLO829/H1437 的真實情況

**結論不變**：TO LOH 不能用作直接 FP filter
**方向澄清**：TO LOH 是 **負方向 FP risk**（= TP 富集），這一方向修正後更清晰、更一致

---

## 7. 下一步

1. **Phase 1 feature table 更新**：
   - 加入 `to_tier_a_loh_like`（0.706×，-FP risk）和 `to_tier_a_plus_loh_like`（0.766×）
   - COLO829 標記為低信心（LOH feature 不穩定）

2. **COLO829 TO phasing 根因調查**（可選）：
   - COLO829 eff_hp median = 25，為何低於其他樣本？
   - 可能是 tumor purity、haplotype block length、或 melanoma 特異的 genome structure

3. **TO vs Paired LOH 機制差異**（Phase 2 方向）：
   - TO LOH = TP 富集：tumor-only somatic SNV 傾向落在 allele-imbalanced region（正常）
   - Paired LOH = FP 富集：germline SNP 在 LOH block 中被誤認（異常）
   - 這個差異可用於 mode-specific evidence weighting

4. **此報告已完成的分析**：
   - Round 1 TO LOH（舊數據）已正式作廢，本報告為正確基準
   - 所有 Phase 1 引用 TO LOH 數字應更新為本報告結論

---

## 附錄：輸出檔案清單

```
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260330_to_loh_enrichment_post_hp_fix/
├── loh_eligibility_before_after.tsv   # 修正前後 LOH eligible 比較（全樣本 TP/FP）
├── loh_enrichment_by_tier.tsv         # 全域聚合 enrichment by Tier（C0/C/B/A/A+/A≥30）
├── loh_enrichment_by_sample_tier.tsv  # Per-sample × Tier 完整 enrichment 表
├── cross_mode_comparison.tsv          # TO new vs TO old vs paired 對比
├── summary.json                       # 結構化摘要
└── figures/
    ├── fig1_loh_eligibility_before_after.png   # LOH eligible 修正前後柱狀圖
    ├── fig2_loh_enrichment_by_tier.png         # Tier enrichment 折線/柱狀複合圖
    ├── fig3_per_sample_enrichment_heatmap.png  # Per-sample × Tier enrichment heatmap
    └── fig4_eff_hp_distribution_focus.png      # COLO829/H1437/H2009 eff_hp 分佈圖
```
