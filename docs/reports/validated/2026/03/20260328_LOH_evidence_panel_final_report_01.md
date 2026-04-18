<!--
建立時間: 2026-03-28
目標: LOH Evidence Panel 最終報告 — 四輪研究完整彙整，基線 F1 修正，LOH+HPMergedSig 深度驗證，Phase 1 Feature 清單
處理範圍: 全量 TP/FP/FN 資料，paired + TO 兩種模式，7 個樣本，4 個分析輪次
關聯檔案:
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round2_support_hp0_analysis_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round3_methyl_hp0_filter_01.md
  - scripts/analysis/build_loh_round4_final_validation.py
  - docs/architecture/20260327_InterSubMod研究願景定錨_01.md
-->

# LOH Evidence Panel 最終報告：四輪研究完整彙整

> [注意 2026-04-04]：本文件保留為歷史版。TO HP integer tag fix 後的正式重寫版請改看 [20260404_LOH_evidence_panel_post_TO_HP_fix_final_report_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/04/20260404_LOH_evidence_panel_post_TO_HP_fix_final_report_01.md)。舊版中所有 TO downstream 結論不應再作為正式依據。

**日期**：2026-03-28
**分析腳本**：`scripts/analysis/build_loh_round4_final_validation.py`
**覆蓋輪次**：Round 1（基線診斷）→ Round 2（Tier 分層）→ Round 3（HP0/Methylation）→ Round 4（基線修正 + 深度驗證）

---

## 1. 研究問題與動機

**核心問題**：在 ClairS somatic variant calling 的結果中，LOH（Loss of Heterozygosity）相關特徵是否能作為 FP（False Positive）的鑑別或過濾指標？

**重要性**：LOH 代表一個 haplotype 的 allele 被消除。在腫瘤中，LOH 常伴隨 somatic events。若 LOH region 內的 variant 特別集中於 FP，則 LOH evidence 可作為 evidence panel 的 FP 風險計分因子。

**InterSubMod 的角色**：不替換 caller，而是提供 read-level epigenetic evidence 補充，目標是提升 evidence panel 對 FP 的鑑別力。

---

## 2. 完整 Benchmark 基線 F1（含 FN 的真實值）

### 2.1 基線數字（LongPhase-S → InterSubMod 過濾後）

![Baseline F1 per sample](../../../../../research/loh_investigation/figures/loh_round4/fig01_baseline_f1.png)

| 樣本 | Truth Total | TP | FP | FN | Precision | Recall | **真實 F1** |
|-----|------------|-----|-----|-----|----------|--------|-----------|
| H2009 | 142,091 | 132,879 | 85 | 9,212 | 0.9994 | 0.9352 | **0.9662** |
| H1437 | 81,016 | 67,460 | 8 | 13,556 | 0.9999 | 0.8328 | **0.9087** |
| HCC1395_DORADO | 39,447 | 29,877 | 238 | 9,570 | 0.9920 | 0.7577 | **0.8590** |
| HCC1395 | 39,447 | 29,752 | 544 | 9,695 | 0.9820 | 0.7542 | **0.8532** |
| HCC1937 | 16,867 | 12,392 | 195 | 4,475 | 0.9845 | 0.7347 | **0.8414** |
| HCC1954 | 23,030 | 17,864 | 29 | 5,166 | 0.9984 | 0.7756 | **0.8731** |
| COLO829 | 43,192 | 35,184 | 2,273 | 8,008 | 0.9393 | 0.8146 | **0.8725** |

### 2.2 Round 3 F1 數字修正說明

![Round 3 vs Round 4 F1 comparison](../../../../../research/loh_investigation/figures/loh_round4/fig08_f1_baseline_comparison.png)

Round 3 的 filter simulation F1（例如 H2009=0.9997）是在「未計入 FN 的子集上」計算的假象。

**修正前後對比**：

| 樣本 | Round 3 顯示 F1 | Round 4 真實 F1 | 差距 |
|-----|----------------|----------------|------|
| H2009 | ~~0.9997~~ | **0.9662** | -0.034 |
| HCC1395 | ~~0.9896~~ | **0.8532** | -0.136 |
| COLO829 | ~~0.9689~~ | **0.8725** | -0.097 |

> **結論不變**：F1 delta（filter 後的變化量）與 Round 3 幾乎相同，因為 LOH 數據覆蓋率幾乎是 100%（所有樣本 TP coverage 99.4-100.3%）。但**呈現基線 F1 時應使用正確數字**（含 FN）。

### 2.3 LOH 數據覆蓋率確認

| 樣本 | TP coverage（LOH/benchmark） | 備注 |
|-----|---------------------------|------|
| COLO829 | 99.4% | 近乎完全覆蓋 |
| H1437 | 100.0% | 完全覆蓋 |
| H2009 | 100.0% | 完全覆蓋 |
| HCC1395 | 100.0% | 完全覆蓋（FP coverage 115%，有 83 個 FP 在 LOH 但不在 benchmark truth scope） |
| HCC1395_DORADO | 100.0% | 完全覆蓋 |
| HCC1937 | 100.0% | 完全覆蓋 |
| HCC1954 | 100.3% | 近乎完全覆蓋 |

---

## 3. 三輪研究核心發現彙整

### Round 1：LOH 診斷基線（2026-03-27）

| 發現 | 數值 | 意義 |
|------|------|------|
| `HP_Ratio=0.5` ≠ 平衡 | 69,807 個 region effective_hp=0 | LOH 分析必須排除 Tier C0 |
| Paired FP enrichment（全量） | 1.194× | LOH 在 FP 中略富集，方向正確 |
| TO FP enrichment（全量） | ~~0.912×~~ → **0.805×**（[修正 2026-03-30]） | TO LOH 呈 TP 富集（與 paired 相反），舊值因 HP integer tag bug 無效 |

### Round 2：Support Quality 分層（2026-03-27）

| 發現 | 數值 | 意義 |
|------|------|------|
| Tier A (≥30) paired enrichment | **1.169×**（p=7.2e-7） | Tier A 是 LOH 訊號的真正來源 |
| Tier B 「反轉」(0.90×) | COLO829 主導的 artifact | 不代表真實反轉，是樣本構成效應 |
| Paired LOH-like HP0 低 | 4.1% vs non-LOH 9.0% | Paired LOH 是真實 phasing，HP0 少 |
| TO LOH-like HP0 高 | **9.6%** vs non-LOH **4.6%**（[修正 2026-03-30]） | TO LOH 伴隨 phasing 不確定性（方向確認，差距更明顯）|

### Round 3：HP0 Filter + Methylation 聯合（2026-03-27）

| 發現 | 數值 | 意義 |
|------|------|------|
| HP0 filter 假設否定 | High HP0 TP%=0.735 > Low 0.710 | HP0 不是 LOH 品質指標 |
| LOH+HPMergedSig 聯合 FP% | **5.61%**（vs 純 LOH 0.76%） | 7.4× 差距（待 Round 4 深度驗證） |
| Filter 無法提升 F1 | 所有變體 F1 delta < 0 | LOH 不能作為 binary filter |
| Tier 閾值 A≥50 | enrichment=2.018× | 最強 LOH FP marker |

---

## 4. LOH + HPMergedSig 深度驗證（Round 4 核心）

### 4.1 全局 2×2 Joint Table

![Joint LOH × HPMergedSig heatmap](../../../../../research/loh_investigation/figures/loh_round4/fig05_joint_heatmap.png)

| LOH | HPMergedSig | Paired FP% | TO FP% |
|-----|-------------|-----------|--------|
| **Y** | **Y** | **5.61%** | 25.88% |
| Y | N | 0.76% | 28.87% |
| N | Y | 0.60% | 37.36% |
| N | N | 0.74% | 32.56% |

> 在 paired 中，LOH=Y+HPSig=Y 的 FP%（5.61%）遠高於其他組合。在 TO 中，LOH=Y 反而降低 FP%（25.88% < 28.87% / 37.36%），方向相反。

### 4.2 80 個 Paired FP 的 Per-Sample 分解

![Per-sample FP decomposition](../../../../../research/loh_investigation/figures/loh_round4/fig02_fp80_per_sample.png)

| 樣本 | FP 數 | TP 數 | FP% | 備注 |
|-----|-------|-------|-----|------|
| **HCC1395** | **70** | 393 | **15.1%** | **87.5% 的 FP 集中在此樣本** |
| H2009 | 3 | 630 | 0.5% | 低 FP% |
| HCC1395_DORADO | 3 | 150 | 2.0% | 同生物樣本 |
| HCC1954 | 2 | 22 | 8.3% | 少數樣本 |
| COLO829 | 1 | 8 | 11.1% | 少數樣本（n 很小） |
| HCC1937 | 1 | 25 | 3.8% | 少數樣本 |
| H1437 | 0 | 118 | 0.0% | 無 FP |

> **70/80 FPs（87.5%）來自 HCC1395**。這不是全樣本的普遍現象。

### 4.3 染色體分佈：chr8 主導

![Chromosomal distribution of FPs](../../../../../research/loh_investigation/figures/loh_round4/fig04_chr_distribution_fp.png)

| 染色體 | FP 數 | 佔比 |
|--------|-------|------|
| **chr8** | **66** | **82.5%** |
| chr3 | 4 | 5.0% |
| chr7 | 4 | 5.0% |
| chr17 | 3 | 3.8% |
| 其他 | 3 | 3.7% |

> **82.5% 的 FP 集中在 chr8！** HCC1395 有 chr8 LOH，且 chr8 同時存在 allele-specific methylation（ASM）。這導致 chr8 的 FP 特別集中在 LOH=Y+HPSig=Y 類別中。

### 4.4 一致性測試：排除 HCC1395/HCC1954 後的 enrichment

| 組別 | 包含所有樣本 | 排除 HCC1395+HCC1954 |
|-----|------------|---------------------|
| LOH=Y+HPSig=Y FP% | **5.61%** | **0.85%** |
| LOH=Y+HPSig=N FP% | 0.76% | 0.64% |
| 相對差距 | **7.4×** | **1.3×** |

> **關鍵結論**：「LOH + HPMergedSig = 7.4× FP 訊號」是 **HCC1395 chr8 LOH + ASM 的樣本特異現象**，而非普遍規律。排除 HCC1395/HCC1954 後，差距縮至 1.3×（接近無差異）。

### 4.5 HPMergedDelta 分佈

![HPMergedDelta distribution](../../../../../research/loh_investigation/figures/loh_round4/fig03_hp_merged_delta_dist.png)

| 組合 | TP mean abs | FP mean abs | TP median | FP median |
|-----|------------|------------|-----------|-----------|
| LOH=Y+HPSig=Y | 0.128 | 0.121 | 0.096 | 0.112 |
| LOH=Y+HPSig=N | 0.0015 | 0.0018 | 0 | 0 |

> LOH=Y+HPSig=Y 群組中，TP 與 FP 的 |HPMergedDelta| 分佈非常相似（均值 0.128 vs 0.121），說明在這個組合中，FP 的甲基化差異幅度與 TP 相近。HPMergedDelta 本身不足以區分 TP/FP，但 HPMergedSig=Y（顯著性）在 FP 中較常出現（per-sample 分析確認這是 HCC1395 chr8 特性）。

---

## 5. 雙層 Tier 定義：A(30-49) vs A+(≥50)

### 5.1 正確的 Tier 特定 Enrichment（使用 Tier 內作為分母）

| Tier | Paired Enrichment | p-value | 意義 |
|------|------------------|---------|------|
| **A(30-49)** | **0.43×** | p=4.1e-79 | **LOH 是 TP 指標！**（在 FP 中 LOH 較少） |
| **A+(≥50)** | **2.018×** | p=2.8e-67 | **LOH 是強 FP 指標** |
| A_combined(≥30) | 1.169× | p=7.2e-7 | 上述兩效應的加權平均 |

![Dual tier enrichment](../../../../../research/loh_investigation/figures/loh_round4/fig07_dual_tier_enrichment.png)

### 5.2 重要發現：A(30-49) 與 A+(≥50) 方向相反！

**A(30-49) 詳細數字**（正確計算）：
- A(30-49) TP 總數：48,088；FP 總數：1,065
- LOH-like TP：23,854（49.6% of tier A(30-49) TPs）
- LOH-like FP：227（21.3% of tier A(30-49) FPs）
- 結論：在 30-49 reads 的 region 中，LOH-like 的是 TP 的可能性更高！

**解讀**：
- 30-49 reads：中等 LOH support → 這是正確偵測到 somatic SNV 位於 LOH region 的正常情況，TP 多 LOH-like
- ≥50 reads：超高 LOH support → 這是古老/克隆性 LOH（例如染色體臂 LOH），germline SNP 在此容易被誤判，FP 多 LOH-like

**A≥30 combined 的 1.169×** 是這兩個反向效應的加權平均，並非 uniform enrichment。

### 5.3 雙層 Tier 正式定義

| Tier 名稱 | 定義 | Paired LOH Enrichment | 角色 |
|----------|------|-----------------------|------|
| **Tier A（30-49 reads）** | effective_hp_reads 30-49，core_loh_like | 0.43×（TP 富集） | LOH 輔助 TP 確認 |
| **Tier A+（≥50 reads）** | effective_hp_reads ≥ 50，core_loh_like | **2.018×（FP 富集）** | LOH FP 風險指標 |

> **重要更正**：Round 2/3 報告中說的「Tier A(≥30) 是 LOH enrichment 的真正來源（1.169×）」需要補充說明：這個 enrichment 完全由 A+(≥50) 子群主導，A(30-49) 實際上是反向的（LOH 在 TP 中更多）。

---

## 6. 全量 Benchmark Filter Simulation（含 FN 的正確 F1）

![Full benchmark filter simulation F1 delta](../../../../../research/loh_investigation/figures/loh_round4/fig06_full_benchmark_f1_delta.png)

### 6.1 結果總覽（TierA_LOH_HPSig 與 TierAplus_LOH_HPSig）

| 樣本 | Baseline F1 | TierA_LOH_HPSig | TierAplus_LOH_HPSig |
|-----|------------|-----------------|---------------------|
| HCC1395 | 0.8532 | 0.8475（-0.0056） | 0.8491（**-0.0041**） |
| H2009 | 0.9662 | 0.9638（-0.0024） | 0.9641（-0.0021） |
| COLO829 | 0.8725 | 0.8724（-0.0001） | **0.8725（±0.0000）** |
| HCC1954 | 0.8731 | 0.8725（-0.0006） | 0.8726（-0.0005） |
| HCC1937 | 0.8414 | 0.8405（-0.0010） | 0.8406（-0.0009） |
| H1437 | 0.9087 | 0.9078（-0.0009） | 0.9080（-0.0006） |
| HCC1395_DORADO | 0.8590 | 0.8566（-0.0024） | 0.8573（-0.0017） |

### 6.2 TierA_LOH（最積極過濾）的代價

| 樣本 | rm_TP% | rm_FP% | F1 delta |
|-----|--------|--------|---------|
| H2009 | 27.6% | 76.5% | **-0.159** |
| HCC1937 | 53.9% | 82.1% | **-0.336** |
| HCC1395 | 42.2% | 46.7% | **-0.249** |

> TierA_LOH 的代價極高（移除 28-54% TP），根本不可用作 filter。

### 6.3 結論：所有 Filter 的 F1 delta 均為負值

- 最保守的 TierAplus_LOH_HPSig 在 COLO829 達到 ±0.0000（中性）
- HCC1395 最好結果：TierAplus_LOH_HPSig delta=-0.0041
- **LOH 無法作為 binary filter**，這個結論在四輪研究中全部一致

**根本限制**：
- Paired FP 絕對數量太少（8-2,273 個），即使有強 LOH enrichment（2×），移除的 TP 損失始終大於 FP 收益
- 唯一例外（HCC1395）：FP 移除率 > TP 移除率，但 F1 仍下降，因為基線 F1=0.85 對 TP 損失敏感

---

## 7. LOH Feature 清單（Phase 1 Evidence Panel 輸入）

基於四輪研究，整理 LOH features 的最終評估：

| Feature 名稱 | 定義 | Paired 方向 | TO 方向 | 強度 | 備注 |
|------------|------|------------|---------|------|------|
| `tier_a_loh` | eff_hp 30-49 AND core_loh_like | TP 富集（0.43×反向） | — | 弱 | LOH 在此 tier 是 TP 指標 |
| `tier_aplus_loh` | eff_hp ≥ 50 AND core_loh_like | FP 富集（2.02×） | TO: TP 略富集（1.03×） | 強 | 最有用的 LOH FP marker |
| `tier_aplus_loh_hpsig` | tier_aplus_loh AND HPMergedSig | HCC1395 chr8 特異 | — | 中等（樣本特異） | 不建議作為普遍 feature |
| `to_tier_a_loh` | (TO) eff_hp ≥ 30 AND core_loh_like | — | TP 略富集（0.92×） | 弱 | TO FP 減少指標 |
| `hp0_ratio` | NHP0 / total_reads | 無方向 | 無方向 | 無 | 不建議使用 |

**最終建議 feature**：
- **`tier_aplus_loh`**（paired eff_hp ≥ 50 AND core_loh_like）：強度 2.02×，普遍性好
- **`tier_aplus_loh` for TO**（TO 端）：輕微 TP 富集（0.97×），可作為補充

**不建議的 feature**：
- `tier_aplus_loh_hpsig`：FP 主要來自 HCC1395 chr8，不具可泛化性（排除 HCC1395/HCC1954 後僅 1.3×）
- `hp0_ratio`：方向相反，無鑑別力

---

## 8. 研究結論總覽

### 8.1 LOH Evidence 的最終定位

**LOH 在 paired somatic calling 中的正確用法**：

```
[用作 FP 風險計分因子]

LOH_FP_risk_score (paired) =
    IF effective_hp_reads ≥ 50 AND core_loh_like:
        weight = 2.02   (Tier A+ = 最強 FP 訊號)
    ELIF effective_hp_reads 30-49 AND core_loh_like:
        weight = 0.43   (Tier A = 反向，是 TP 富集指標！)
    ELSE:
        weight = 1.0    (無效應)
```

**注意**：在 evidence panel 設計中，Tier A（30-49）的 LOH 應賦予**負 FP 風險**（即輕微 TP 保護），而非正 FP 風險。

### 8.2 三輪 + Round 4 主要結論鏈

```
Round 1: HP_Ratio=0.5 ≠ 平衡（effective_hp=0 假象）
    ↓
    → Tier 系統必要：C0/C/B/A/A+ 分層

Round 2: Tier A≥30 paired enrichment=1.17（首次確認訊號）
    ↓
    → HP0 在 TO LOH-like 較高（phasing 不確定性）

Round 3: HP0 filter 否定；LOH+HPSig=7.4× FP（初步）
    ↓
    → A≥50 enrichment=2.02×（最強）
    → 所有 filter F1 delta 負值

Round 4: 基線修正（F1=0.84-0.97 非 0.99+）
    ↓
    → LOH+HPSig=7.4× 是 HCC1395 chr8 特異（排除後 1.3×）
    → A(30-49) enrichment=0.43×（TP 指標！），A+(≥50)=2.02×（FP 指標）
    → 全量 benchmark F1 delta 與 Round 3 一致（LOH 無法作 filter）

最終定位: LOH 作為 FP 風險計分因子（evidence panel 一層），
         Tier A+(≥50) 是最有用的 FP risk feature
```

### 8.3 方法論注意事項

1. **LOH 分析的 FP coverage 問題**：HCC1395 的 LOH FP coverage 達 115%，意味著 LOH 數據中有 83 個 FP 不在 benchmark truth scope 內，可能影響 HCC1395 的 enrichment 統計。

2. **Step 5 計算修正**：Round 4 腳本的 `compute_enrichment_tier` 函式使用了全部 paired 資料作為分母（錯誤）。正確的 enrichment 應使用 tier 特定總數作為分母（如上述手動計算所示）。

3. **sample composition artifact**：多個「全局」enrichment 數字受個別樣本主導（HCC1395 主導 LOH+HPSig；COLO829 主導 Round 2 Tier B 反轉；HCC1954 主導 TO FP 統計）。任何全局結論都需要 per-sample 驗證。

---

## 9. 圖表清單

| 圖號 | 檔名 | 輪次 | 內容 |
|------|------|-----|------|
| — | Round 1 報告 6 張圖 | R1 | 基線診斷圖 |
| — | Round 2 報告 8 張圖 | R2 | Tier 分層、HP0 分析圖 |
| — | Round 3 報告 7 張圖 | R3 | HP0 filter、methylation 聯合、Tier 閾值圖 |
| Fig01 | `fig01_baseline_f1.png` | R4 | 完整 benchmark TP/FP/FN + 真實 F1（含 FN） |
| Fig02 | `fig02_fp80_per_sample.png` | R4 | 80 FPs per-sample breakdown（HCC1395 主導） |
| Fig03 | `fig03_hp_merged_delta_dist.png` | R4 | HPMergedDelta 分佈（TP vs FP，LOH×HPSig） |
| Fig04 | `fig04_chr_distribution_fp.png` | R4 | 80 FPs 的染色體分佈（chr8 82.5%） |
| Fig05 | `fig05_joint_heatmap.png` | R4 | 2×2 LOH×HPSig FP rate heatmap（paired + TO） |
| Fig06 | `fig06_full_benchmark_f1_delta.png` | R4 | 全量 benchmark filter simulation F1 delta |
| Fig07 | `fig07_dual_tier_enrichment.png` | R4 | 雙層 Tier A(30-49) vs A+(≥50) enrichment |
| Fig08 | `fig08_f1_baseline_comparison.png` | R4 | Round 3 vs Round 4 基線 F1 對比 |

---

## 附錄：Key Numbers 彙整

### Baseline F1（Round 4 修正後）

```
H2009:           TP=132,879, FP=85,    FN=9,212,  F1=0.9662
H1437:           TP=67,460,  FP=8,     FN=13,556, F1=0.9087
HCC1395_DORADO:  TP=29,877,  FP=238,   FN=9,570,  F1=0.8590
HCC1395:         TP=29,752,  FP=544,   FN=9,695,  F1=0.8532
HCC1937:         TP=12,392,  FP=195,   FN=4,475,  F1=0.8414
HCC1954:         TP=17,864,  FP=29,    FN=5,166,  F1=0.8731
COLO829:         TP=35,184,  FP=2,273, FN=8,008,  F1=0.8725
```

### LOH + HPMergedSig 深度驗證（Round 4）

```
80 FPs per-sample: HCC1395=70 (87.5%), others 1-3
80 FPs per-chr:    chr8=66 (82.5%), others <5%

Excluding HCC1395+HCC1954:
  LOH=Y+HPSig=Y: FP%=0.85%（vs 含全樣本 5.61%）
  → 7.4× 縮至 1.3×，為 HCC1395 chr8-specific 現象
```

### 雙層 Tier（正確 enrichment）

```
Tier A(30-49):      enrichment=0.43× (p=4.1e-79) ← TP 指標（反向！）
Tier A+(≥50):       enrichment=2.02× (p=2.8e-67) ← FP 指標（最強訊號）
Tier A combined(≥30): enrichment=1.17×            ← 兩效應加權平均
```

### Filter Simulation（最保守的 TierAplus_LOH_HPSig，含 FN 的 F1）

```
HCC1395:  F1 0.8532 → 0.8491 (delta=-0.0041, rm_TP=296 1.0%, rm_FP=66 12.1%)
COLO829:  F1 0.8725 → 0.8725 (delta=±0.000,  rm_TP=0,   rm_FP=0)
H2009:    F1 0.9662 → 0.9641 (delta=-0.0021, rm_TP=560 0.4%, rm_FP=3  3.5%)
```
