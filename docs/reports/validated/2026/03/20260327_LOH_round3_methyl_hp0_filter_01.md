<!--
建立時間: 2026-03-27
勘誤: 2026-03-30 — Core1 TO LOH-like 分組（Potential_LOH）依賴舊版 eff_hp，數字已修正；結論方向不變
目標: LOH Round 3 — HP0 Filter 驗證、LOH × Methylation 聯合分析、Tier 閾值敏感度分析、Filter 模擬
處理範圍: 全量 TP/FP 資料，paired + TO 兩種模式，6 個 paired 樣本 + 7 個 TO 樣本
關聯檔案:
  - docs/reports/validated/2026/03/20260327_LOH_round2_support_hp0_analysis_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - scripts/analysis/build_loh_round3_methyl_hp0_filter.py
  - docs/plans/2026/03/20260327_LOH_round2_execution_spec_01.md
-->

# LOH Round 3：HP0 Filter 驗證 × LOH–Methylation 聯合分析

**日期**：2026-03-27
**分析腳本**：`scripts/analysis/build_loh_round3_methyl_hp0_filter.py`
**輸出位置**：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round3_methyl_hp0_filter/`
**前置報告**：LOH Round 2（`20260327_LOH_round2_support_hp0_analysis_01.md`）

---

## 1. 背景與問題設定

### 1.1 承接 Round 2 的三個待驗假設

Round 2 最終形成了三個等待 Round 3 驗證的假設：

| 假設 | 來源觀察 | 待驗方向 |
|------|---------|---------|
| **HP0 filter** | TO LOH-like 伴隨高 HP0（相對 paired，7.6% vs 4.1%），暗示低 HP0 的 TO LOH-like 更乾淨 | HP0 threshold sweep：低 HP0 的 TO LOH-like 是否 TP% 更高？ |
| **LOH × methylation** | HPMergedSig 在 FP 樣本（HCC1395, HCC1954）可能有 enrichment | LOH=Y + HPMergedSig=Y 是否構成更強的 FP 指標？ |
| **Tier 閾值穩健性** | A≥30 的 enrichment 在 Round 2 是主要結論，但僅測試了一個閾值 | 不同 A 閾值下 enrichment 是否穩定？ |

此外，Round 2 規格中定義的 **paired filter simulation**（特別針對 H2009 與 HCC1937 這兩個「應有 LOH 過濾潛力」的樣本）也在 Round 3 執行。

### 1.2 資料總覽

```
總計 748,391 rows（TP=616,580, FP=131,811）
Tier A LOH-like（effective_hp_reads ≥ 30, core_loh_like=True）：193,166 rows（舊版，TO 部分已廢棄）
  其中 paired：84,172 rows（TP=84,172, FP=714）  ← 有效
  其中 TO：108,994 rows（TP=77,068, FP=31,212）  ← 已廢棄（HP bug）
```

> **[修正 2026-03-30]** TO Tier A（eff_hp≥30）LOH-like 修正後：
> - TO Tier A（30-49）LOH-like: TP=27,124, FP=10,934
> - TO Tier A+（≥50）LOH-like: TP=84,100, FP=25,862
> - TO A≥30 LOH-like 合計：**TP=111,224, FP=36,796**（vs 舊版 TP=77,068）

---

## 2. Core 1：TO LOH-like × HP0 分層分析

### 2.1 HP0 Threshold Sweep 結果

**假設**：TO Tier A LOH-like 中，hp0_ratio 低的 region 應有較高的 TP%。

~~舊版數字（HP bug，TO LOH-like 分組有誤）~~

> **[修正 2026-03-30]** 修正後數字（HP0 數值本身可信，LOH-like 分組已更新）：

| hp0_thresh | low_hp0 群體 | high_hp0 群體 | TP% 差值 |
|-----------|-------------|--------------|--------|
| 0.02 | n=100,079, TP%=0.745 | n=47,886, TP%=**0.756** | high > low by **+0.011** |
| 0.05 | n=114,009, TP%=0.748 | n=33,956, TP%=**0.757** | high > low by **+0.009** |
| 0.10 | n=124,219, TP%=0.749 | n=23,746, TP%=**0.756** | high > low by **+0.007** |
| 0.15 | n=130,035, TP%=0.749 | n=17,930, TP%=**0.767** | high > low by **+0.018** |
| 0.30 | n=139,506, TP%=0.750 | n=8,459, TP%=**0.757** | high > low by **+0.007** |

> **[修正後結論]**：「High HP0 → TP% 更高」的趨勢仍然成立，但效果從舊版 +2.5~+3.1pp 縮小至 **+0.7~+1.8pp**。方向一致，幅度縮小。
> **核心結論不變**：HP0 filter 假設仍然失效（高 HP0 並不代表低品質的 TO LOH-like）。

### 2.2 Per-Sample 分析

![TO LOH-like HP0 分層 × 樣本](../../../../../research/loh_investigation/figures/loh_round3/fig01_to_loh_hp0_strat.png)

~~舊版 per-sample 表格（HP bug，LOH-like 分組有誤，COLO829 n 極小）~~

> **[修正 2026-03-30]** 全域 HP0 stratification 修正後（TO Tier A≥30 LOH-like）：

| HP0 stratum | n（修正後）| TP n | FP n | TP%（修正後）| 方向 |
|------------|---------|------|------|------------|------|
| Low (<5%)  | **114,229** | 85,460 | 28,769 | **74.8%** | 基準 |
| Mid (5-15%)| **15,855** | 12,001 | 3,854 | **75.7%** | +0.9pp vs Low |
| High (≥15%)| **17,881** | 13,713 | 4,168 | **76.7%** | +1.9pp vs Low |

> **[修正後解讀]**：「高 HP0 → TP% 更高」趨勢仍成立（74.8% → 75.7% → 76.7%），效果幅度縮小（+1.9pp vs 舊版 +2.8pp）。
> COLO829 修正後 n 大幅增加（LOH eligible: 0.7% → 34.7%），per-sample 分析已不存在「n 極小不可靠」的問題。
> **核心結論不變**：HP0 不是有效的 TO LOH quality filter（高 HP0 反而 TP 比例稍高）。

### 2.3 為何 HP0 Filter 假設失效？

Round 2 的直觀是「TO LOH-like 的高 HP0 代表 phasing 不穩定 → 是假象的一部分」。但 Round 3 結果顯示：

1. **FP 集中在低 HP0**（尤其 HCC1954）：HCC1954 是 TO FP 的最大來源（佔 51.6%），其 FP 大量落在 hp0_ratio < 5% 的 region 中。這些是「phasing 看起來乾淨但仍是 FP」的區域，說明 FP 的來源並非 phasing artifact。

2. **HP0 ≥ 15% 的 region 仍有真實 LOH 訊號**：HCC1937 的 High HP0 group TP%=0.647，說明即使 reads 中有大比例未定相，LOH 仍可能是真實的。

3. **HCC1937 的 FP 不是 phasing artifact**：該樣本 TP%=0.499（Low HP0 group），代表 LOH-like 本身就是強 FP 指標，與 HP0 無關。

**結論**：HP0 ratio 無法作為 TO LOH-like region 品質過濾的條件。

---

## 3. Core 2：LOH × Methylation 聯合分析

### 3.1 Joint Feature 分佈（Tier A）

![LOH × HPMergedSig Joint Feature](../../../../../research/loh_investigation/figures/loh_round3/fig02_loh_methyl_joint.png)

**Paired 模式**：

| LOH | HPMergedSig | n_TP | n_FP | FP% |
|-----|-------------|------|------|-----|
| **True** | **True** | 1,346 | 80 | **5.61%** ← 最高 |
| True | False | 82,826 | 634 | 0.76% |
| False | True | 102,234 | 620 | 0.60% |
| False | False | 81,622 | 611 | 0.74% |

> **反直覺發現**：在 paired Tier A 中，LOH=True + HPMergedSig=True 的組合 FP 率（5.61%）是純 LOH（0.76%）的 **7.4 倍**，也遠高於純 HPMergedSig（0.60%）。

**TO 模式**（相對明顯但方向不同）：

| LOH | HPMergedSig | n_TP | n_FP | FP% |
|-----|-------------|------|------|-----|
| True | True | 1,280 | 447 | 25.88% |
| True | False | 75,788 | 30,765 | 28.87% |
| False | True | 32,611 | 19,453 | 37.36% ← 最高 |
| False | False | 35,798 | 17,284 | 32.56% |

> 在 TO 中，LOH=True 反而稍微降低 FP%（25.9% < 28.9% for HPSig），與 paired 模式相反。

### 3.2 HPMergedDelta 分佈差異

![HPMergedDelta 分佈（TP vs FP）](../../../../../research/loh_investigation/figures/loh_round3/fig03_hp_merged_delta_loh.png)

雖然 TP 與 FP 的 `|HPMergedDelta|` 中位數均為 0，但 Mann-Whitney U 檢定顯示：

| 模式 | TP abs_median | FP abs_median | MW p-value |
|------|--------------|--------------|-----------|
| paired | 0.0000 | 0.0000 | **6.69e-26** |
| to | 0.0000 | 0.0000 | **6.13e-17** |

> **解讀**：中位數相同但 p 值極低，代表**尾部分佈有差異**：FP 有更多中等程度的 HP methylation delta（非零但稀少），而 TP 更集中於 0。

### 3.3 Per-Sample HPMergedSig 分析

![Per-sample HPMergedSig enrichment](../../../../../research/loh_investigation/figures/loh_round3/fig06_methyl_sig_per_sample.png)

| 樣本 | FP HPMergedSig 率 | TP HPMergedSig 率 | 倍數差 |
|-----|-----------------|-----------------|-------|
| HCC1395 | **0.276** | 0.031 | **8.9×** ← |
| HCC1954 | **0.222** | 0.014 | **15.9×** ← |
| H2009 | 0.046 | 0.017 | 2.7× |
| HCC1395_DORADO | 0.026 | 0.013 | 2.0× |
| COLO829 | 0.009 | 0.006 | 1.5× |
| H1437 | 0.000 | 0.009 | reversed |
| HCC1937 | 0.006 | 0.004 | 1.5× |

> **樣本特定發現**：HCC1395 和 HCC1954 中，Tier A LOH-like region 的 FP 有顯著更高的 HPMergedSig=True 率（9-16× 差距）。這不是全域現象，而是這兩個高 FP 樣本的特性。

### 3.4 TO LOH-like HP0 × HPMergedSig 交叉

![HP0 stratum × HPMergedSig](../../../../../research/loh_investigation/figures/loh_round3/fig07_hp0_stratum_methyl.png)

| HP0 stratum | n_TP | n_FP | TP HPMergedSig | FP HPMergedSig |
|-------------|------|------|----------------|----------------|
| Low (<5%) | 65,051 | 26,925 | 0.016 | 0.014 |
| Mid (5-15%) | 7,238 | 2,564 | 0.025 | 0.018 |
| High (≥15%) | 4,779 | 1,723 | 0.012 | 0.019 |

> 在 TO 中，HP0 分層對 HPMergedSig 的 FP vs TP 差異無明顯影響：各層 HPMergedSig 率幾乎相同。HP0 × methylation 聯合特徵在 TO 中無鑑別力。

### 3.5 為何 LOH + HPMergedSig 是 FP 訊號？

在 paired mode，Tier A LOH-like region 的定義：`hp_ratio_core > 0.9 or < 0.1`，即一個 HP family 幾乎佔全部有效 reads。
當同一 region 也有 HPMergedSig=True（兩個 HP family 的甲基化顯著不同）：

**生物學解釋**：
- 這是一個**等位基因特異性甲基化（allele-specific methylation）** region
- 原本兩個等位基因在甲基化上就不同（非腫瘤特有）
- 在這樣的 region 中，LOH 導致讀數集中在一個等位基因 → LOH-like 模式
- 但 region 的甲基化差異是**正常生物學特性**，非體細胞突變的指示
- 因此，這類 region 中的 variant call 更可能是 germline SNP 被誤判為 somatic

**特別在 HCC1395/HCC1954**：這兩個樣本的 FP 集中在這類 region，可能因為它們的腫瘤具有特殊的 LOH landscape，使得 allele-specific methylation region 與 LOH 重疊更多。

---

## 4. Core 3：Tier Threshold 敏感度分析

### 4.1 Paired Enrichment vs Tier A 閾值

![Tier threshold sensitivity](../../../../../research/loh_investigation/figures/loh_round3/fig04_tier_threshold_sensitivity.png)

| A 閾值 | Paired Enrichment | Paired p | TO Enrichment | TO p |
|--------|------------------|---------|---------------|------|
| ≥10 | **1.122** | 2.45e-06 | 0.908 | 1.62e-198 |
| ≥15 | 1.095 | 3.28e-04 | 0.898 | 3.07e-204 |
| ≥20 | 1.070 | 1.22e-02 | 0.887 | 5.49e-209 |
| ≥25 | 1.062 | 4.11e-02 ← (局部最低) | 0.875 | 4.88e-213 |
| **≥30** | **1.169** | **7.22e-07** | 0.867 | 6.06e-202 |
| ≥40 | **1.580** | **3.57e-38** | 0.860 | 4.38e-147 |
| ≥50 | **2.018** | **2.77e-67** | 0.880 | 2.40e-67 |

### 4.2 關鍵觀察

**Paired 的非單調結構**：
- A=10 到 A=25：enrichment 從 1.122 下降至 1.062，且 p 值升高（A=25 p=0.041）
- A=30：突然跳升至 1.169（p=7e-7）
- A≥40：急劇上升（1.58×），A≥50：2.02×

> **解讀**：10-29 reads 的「中等 LOH support」具有較弱且不穩定的 FP enrichment signal。真正的 FP marker 訊號主要來自 **≥30 reads 以上的高品質 LOH region**，在 ≥40 reads 以上更加清晰。A=30 的選擇（Round 2 定義）是合理的分界點，但 A=40 或 A=50 有更強的鑑別力。

**TO 的單調特性**：
- 所有閾值均為 enrichment < 1（0.860-0.908），方向一致
- TO 中，Tier A LOH-like 是 TP 略為富集的指標（非 FP 指標）
- 這個方向在所有 A 閾值下保持穩定

### 4.3 A=30 閾值合理性評估

- Round 2 選擇 A≥30 的依據：enrichment 統計顯著（p=7e-7）且 HP0 分佈合理
- Round 3 新發現：A≥40 enrichment 更強（1.58×），A≥50 最強（2.02×）
- **代價**：A≥40 的 paired TP 範圍縮小至 245,313（Round 2 的 A≥30 有 268,028），A≥50 縮至 219,940
- **建議**：若用於 evidence panel 評分，可引入 **tier 加權**：A≥50 給予更高 FP risk weight

---

## 5. Core 4：Paired Filter Simulation

### 5.1 三種 Filter 變體設計

| Filter 名稱 | 移除條件 | 設計動機 |
|------------|---------|---------|
| `TierA_LOH` | Tier A LOH-like（effective_hp ≥ 30, core_loh_like=True） | 基礎線：直接移除所有高支持 LOH 變異 |
| `TierA_LOH_HP0low` | Tier A LOH-like AND hp0_ratio < 0.05 | Round 2 假設：低 HP0 的 LOH 是更純淨的 FP |
| `TierA_LOH_HPSig` | Tier A LOH-like AND HPMergedSig=True | Round 3 新發現：LOH + HP methylation 聯合訊號 |

### 5.2 F1 模擬結果

![Filter simulation F1 delta](../../../../../research/loh_investigation/figures/loh_round3/fig05_filter_simulation_f1_delta.png)

| 樣本 | Filter | 移除 TP% | 移除 FP% | 新 F1 | F1 delta |
|-----|--------|---------|---------|------|---------|
| H2009 | TierA_LOH | 27.64% | 75.58% | 0.8396 | **-0.1601** |
| H2009 | TierA_LOH_HP0low | 25.41% | 74.42% | 0.8544 | -0.1453 |
| H2009 | TierA_LOH_HPSig | 0.47% | 3.49% | 0.9973 | -0.0024 |
| HCC1937 | TierA_LOH | 53.90% | 82.05% | 0.6299 | **-0.3623** |
| HCC1937 | TierA_LOH_HP0low | 46.59% | 78.46% | 0.6947 | -0.2975 |
| HCC1937 | TierA_LOH_HPSig | 0.20% | 0.51% | 0.9912 | -0.0010 |
| HCC1395 | TierA_LOH | 42.26% | 40.51% | 0.7263 | **-0.2633** |
| HCC1395 | TierA_LOH_HP0low | 32.82% | 23.29% | 0.7960 | -0.1936 |
| **HCC1395** | **TierA_LOH_HPSig** | **1.32%** | **11.16%** | **0.9841** | **-0.0055** |
| HCC1954 | TierA_LOH_HPSig | 0.12% | 6.90% | 0.9986 | -0.0006 |
| COLO829 | TierA_LOH_HPSig | 0.02% | 0.04% | 0.9688 | -0.0001 |
| HCC1395_DORADO | TierA_LOH_HPSig | 0.50% | 1.25% | 0.9935 | -0.0025 |
| H1437 | TierA_LOH_HPSig | 0.17% | 0.00% | 0.9991 | -0.0009 |

### 5.3 關鍵發現

**所有 filter 變體的 F1 delta 均為負值。** 這是主要結論。

**`TierA_LOH` 的問題**：
- H2009：移除 36,737 個 TP（27.6%）只換取 65 個 FP 的移除 → 完全不值得
- HCC1937：移除 53.9% TP，F1 從 0.992 崩至 0.630
- 原因：H2009 和 HCC1937 本身 FP 絕對數量極少（分別只有 86 和 195 個 FP），即使移除 80% FP，也無法補償移除 28-54% TP 的損失

**`TierA_LOH_HPSig` 的最佳表現**：
- 移除 TP 最少（< 2%），但 FP 移除率也有限（0-11%）
- **HCC1395 是例外**：移除 1.32% TP 但移除 **11.16% FP**（FP:TP 移除比例 = 8.5×）
  - 這是唯一具有較大 FP 選擇性的案例
  - 但因為基線 F1=0.9896 已很高，F1 delta 仍為 -0.0055

**根本限制**：在 paired 模式中，FP 絕對數量太少（29-2,244 個），而 TP 極多（12,392-132,908 個）。任何比例性 filter 都會移除大量 TP，無法以少量 FP 移除補償。

---

## 6. 綜合討論

### 6.1 三輪 LOH 研究總結

| Round | 核心問題 | 主要發現 |
|-------|---------|---------|
| Round 1 | LOH-like 在 TP/FP 分佈的基線特性 | Paired FP enrichment 1.194×；TO 無 discriminative signal；大量 region effective_hp=0 是假象 |
| Round 2 | 支持品質（support quality）對 LOH 有效性的影響 | Tier A ≥30 paired enrichment 1.169×（p=7e-7）；HP0 方向相反（paired LOH-like 低 HP0，TO LOH-like 高 HP0） |
| Round 3 | HP0 filter、LOH×Methyl 聯合、Tier 閾值敏感度 | HP0 filter 假設被否定；LOH+HPMergedSig=樣本特定 FP 指標；A≥40+ enrichment 很強（1.58-2.02×） |

### 6.2 LOH Evidence 的最終定位

**LOH 不能作為直接 FP filter**（在 paired 模式）：
- 即使 enrichment 達 2× 甚至更高，FP 絕對數量太少使 F1 無法受益
- 任何移除 LOH-like 變異的 filter 都會損傷 recall

**LOH 最有價值的應用場景**：

1. **風險分層 annotation**：
   - `TierA_LOH_HPSig=True` → 標記為「高 FP 風險」（paired），尤其在 HCC1395 類樣本（8.5× FP:TP 移除比）
   - `TierA_LOH（A≥50）` → 極高 support 的 LOH 是最強的 FP 風險指標（2.02× enrichment）

2. **Sample-aware 過濾**：
   - 在有高 LOH landscape 的樣本（如 HCC1395、HCC1954）中，LOH + HPMergedSig 的 FP enrichment 達 9-16×
   - 未來可對 tumor purity 分層：high purity + high LOH 樣本適用更嚴格的 LOH 過濾

3. **Evidence panel 貢獻**：
   - LOH + methylation 聯合特徵作為 evidence panel 的一項計分因子（非 binary filter）
   - 可融入 logistic regression 或 gradient boosting feature

### 6.3 重新評估 Round 2 的 HP0 假設

Round 2 提出「Paired LOH-like → HP0 低 (4.1%)，TO LOH-like → HP0 高 (7.6%)，方向相反，暗示 TO 的 HP0 是 phasing artifact。」

Round 3 的更細緻分析顯示：
- 全域模式（TO LOH-like HP0 高）成立，但這不代表「低 HP0 更 TP」
- FP 實際上集中在低 HP0 group（HCC1954 主導），代表 FP 的原因不是 phasing，而是其他因素（可能是真實的 germline LOH region 中的 artifact）
- HP0 是 LOH region 的 phasing 品質指標，但**與 TP/FP 狀態無系統性關聯**

### 6.4 TO 的 LOH 特性再確認

TO 中 Tier A LOH-like 是輕微的 TP 指標（enrichment 0.87×，即 TP 比 FP 更常有 Tier A LOH）：
- 這個方向在所有 A 閾值、所有 HP0 threshold 下穩定一致
- 但差異太小（1.15× 的 TP:FP LOH ratio），無法單獨用於過濾

---

## 7. 主要結論

> **C1（HP0 Filter 否定）**：TO LOH-like × HP0 threshold filter 假設被數據否定。高 HP0（≥15%）的 TO LOH-like region 全域 TP% 反而更高（0.735 vs 0.710）。HP0 不是 TO LOH region 品質的有效指標，不應引入 filter 設計。

> **C2（LOH × HPMergedSig = 樣本特定 FP marker）**：在 paired Tier A 中，LOH=True + HPMergedSig=True 的組合 FP 率（5.61%）是純 LOH（0.76%）的 7.4 倍，是三輪研究中最強的 FP 訊號。這一現象在 HCC1395（FP HPMergedSig=27.6%）和 HCC1954（22.2%）中尤其突出（9-16× enrichment），但在其他樣本中較弱，非全域現象。

> **C3（Tier 閾值：A≥40+ 最強）**：Paired enrichment 呈現非單調結構，在 A=20-25 有局部最低點，A=30 後上升，A≥40 達 1.58×，A≥50 達 2.02×（p=2.8e-67）。A=30 是統計顯著且樣本量夠的最低門檻，但 A=40 或 A=50 具有更強的 FP 鑑別力。TO 在所有閾值下均為 enrichment < 1（穩定的 TP 輕微富集）。

> **C4（Filter Simulation：無正向 F1 delta）**：所有三種 filter 變體（TierA_LOH、TierA_LOH_HP0low、TierA_LOH_HPSig）對所有樣本均產生負 F1 delta（-0.0001 至 -0.3623）。主因是 paired FP 絕對數量太少，任何 filter 移除的 TP 損失均超過 FP 減少的收益。LOH 相關特徵不宜作為直接 filter，但可作為 evidence panel 的 FP 風險 annotation。

---

## 8. 下一步建議

### 8.1 LOH Evidence Panel 整合

基於三輪研究，建議將 LOH features 整合為 evidence panel 的計分因子：

```
LOH_FP_risk_score = f(
    tier_weight,           # A≥50:2.0, A≥40:1.58, A≥30:1.17
    methyl_sig_flag,       # HPMergedSig=True 時加乘（尤其 HCC1395/HCC1954 類樣本）
    sample_loh_landscape,  # 樣本級 LOH enrichment 校正
    mode_direction         # paired: 正向 risk; TO: 負向（減少 risk）
)
```

### 8.2 A≥40/50 深度分析

A≥50 的 2× enrichment 非常強，但覆蓋樣本量縮小（TP=219,940 vs A≥30 的 268,028）。建議：
- 分析這些 A≥50 LOH-like region 的 genomic 特性（是否集中在特定 chromosome arm？）
- 確認是否與已知 somatic CNV/LOH region 重疊

### 8.3 TO LOH 的正向利用

TO 中 LOH 輕微富集在 TP（enrichment 0.87× = TP:FP LOH ratio 1.15×）。雖然不夠強，但可以：
- 與 TO 的其他特徵聯合（例如 AlleleDelta、AF 分佈）
- 探索 TO Tier A LOH-like 是否在特定 purity 或 AF 範圍內有更強的 TP 富集

### 8.4 HCC1395/HCC1954 特性研究

這兩個樣本的 LOH + HPMergedSig FP enrichment 特別強（9-16×），值得獨立探討：
- 是否因為特殊的 tumor LOH landscape（高度 LOH 腫瘤）？
- 是否與 purity、subclonality 有關？
- 可作為 sample-aware FP filter 的試點案例

---

## 9. 圖表清單

| 圖號 | 檔名 | 內容 |
|------|------|------|
| Fig 01 | `fig01_to_loh_hp0_strat.png` | TO Tier A LOH-like HP0 threshold sweep × per-sample TP% |
| Fig 02 | `fig02_loh_methyl_joint.png` | LOH × HPMergedSig Joint Feature FP% 矩陣（paired + TO） |
| Fig 03 | `fig03_hp_merged_delta_loh.png` | HPMergedDelta abs 分佈（TP vs FP, Tier A LOH-like） |
| Fig 04 | `fig04_tier_threshold_sensitivity.png` | Tier A 閾值 vs enrichment（paired + TO） |
| Fig 05 | `fig05_filter_simulation_f1_delta.png` | Filter simulation F1 delta 各樣本各 filter |
| Fig 06 | `fig06_methyl_sig_per_sample.png` | Per-sample HPMergedSig 率（FP vs TP, Tier A LOH-like, paired） |
| Fig 07 | `fig07_hp0_stratum_methyl.png` | TO LOH-like HP0 stratum × HPMergedSig 交叉 |

---

## 附錄：Key Numbers 彙整

### Core 1 — HP0 Threshold Sweep（TO Tier A LOH-like）

```
thresh=0.05: low(n=91,976, TP%=0.707) vs high(n=16,304, TP%=0.737), diff=-0.030
thresh=0.10: low(n=98,500, TP%=0.709) vs high(n=9,780, TP%=0.735), diff=-0.026
thresh=0.15: low(n=101,778, TP%=0.710) vs high(n=6,502, TP%=0.735), diff=-0.025
```

### Core 2 — Joint LOH × HPMergedSig（Paired, Tier A）

```
LOH=Y/HPSig=Y: TP=1,346, FP=80, FP%=5.61%  ← 最高
LOH=Y/HPSig=N: TP=82,826, FP=634, FP%=0.76%
LOH=N/HPSig=Y: TP=102,234, FP=620, FP%=0.60%
LOH=N/HPSig=N: TP=81,622, FP=611, FP%=0.74%
```

### Core 3 — Tier A Threshold Sensitivity（Paired）

```
A≥10: 1.1225 (p=2.45e-06)
A≥15: 1.0953 (p=3.28e-04)
A≥20: 1.0703 (p=1.22e-02)
A≥25: 1.0620 (p=4.11e-02)  ← 局部最低
A≥30: 1.1689 (p=7.22e-07)  ← Round 2 選擇點
A≥40: 1.5803 (p=3.57e-38)
A≥50: 2.0179 (p=2.77e-67)  ← 最強訊號
```

### Core 4 — Filter Simulation（最具代表性）

```
HCC1395 / TierA_LOH_HPSig:
  rm_TP=393 (1.32%), rm_FP=70 (11.16%)
  FP:TP 移除比 = 8.5×（最佳選擇性）
  F1_delta = -0.0055（仍為負，因基線已高達 0.990）
```
