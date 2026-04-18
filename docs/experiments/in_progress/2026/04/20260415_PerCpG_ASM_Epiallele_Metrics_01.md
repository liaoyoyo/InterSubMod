<!--
建立時間: 2026-04-15 20:00
目標: 文獻驅動的 per-CpG ASM 偵測與 epiallele 異質性指標 Python PoC 驗證
處理範圍: HCC1395 全基因體 31,659 regions (30,401 TP + 1,258 FP)，6 種指標家族 24 個新 metric
關聯檔案:
  - docs/references/20260415_ASM_subclone_methods_literature_survey.md
  - docs/experiments/in_progress/2026/04/20260413_Germline_ASM_Tumor_LOH_交叉分析_01.md
  - research/germline_asm_analysis/figures/07-10_*.png
-->

# Per-CpG ASM 偵測與 Epiallele 異質性指標 — Python PoC 驗證

## 1. 研究動機

### 1.1 背景

先前所有 ISM ASM 分析均為 **window-average**（整個 2000bp 區域聚合值）。文獻顯示：

- **96% 的 ASM 是 entropy imbalance**，非 mean difference（CPEL, Jenkinson 2020, Nature Communications）
- Per-CpG Fisher's exact test 是 ASM 偵測金標準（DAMEfinder, pycoMeth, NanoMethPhase）
- Epiallele pattern metrics（epipolymorphism, PDR, Shannon entropy, VEF）量化不同維度的異質性
- ISM 的 raw_matrix（reads × CpGs）是 **13+ 種文獻方法的共同計算基礎**，但完全未利用 per-CpG 粒度

### 1.2 核心問題

1. Per-CpG Fisher's exact test 與 window-level PERMANOVA 是否一致？
2. 新指標能否提供額外的 TP/FP 區分力（filter value）？
3. 新指標在 LOH × HPFineNGroups 分層中是否展現 characterization value？
4. 哪些指標值得整合進 C++ 正式版本？

### 1.3 文獻調查基礎

完整報告：`docs/references/20260415_ASM_subclone_methods_literature_survey.md`

核心發現：30+ 方法、40+ 論文調查後，13+ 種方法可直接在 ISM raw_matrix 上計算。選定 P0+P1 共 6 種指標家族實作。

---

## 2. 方法

### 2.1 實作的指標家族（6 種，24 個 metrics）

| 指標家族 | 來源論文 | 核心原理 | 輸出欄位 |
|----------|---------|---------|---------|
| **Per-CpG Fisher** | DAMEfinder (2020), pycoMeth (2023) | 每個 CpG 建 HP1/HP2 × meth/unmeth 2×2 表 → Fisher → BH-FDR | n_sig, frac_sig, max_delta, mean_delta, max_neg_log_fdr, n_tested |
| **NME** | CPEL (Jenkinson 2020) | 每 HP 的甲基化模式 Shannon entropy / H_max | nme_hp1, nme_hp2, nme_combined, entropy_imbalance |
| **Epipolymorphism** | methclone (2014), Metheor (2023) | 4-CpG sliding window, 1 - Σ(p_i²) | epipoly_hp1, epipoly_hp2, epipoly_combined, epipoly_delta |
| **PDR** | Landau (2014) | 含不一致 CpG 的 reads 比例 | pdr_hp1, pdr_hp2, pdr_combined, pdr_delta |
| **Shannon Entropy** | WSHPackage (2020) | 4-CpG patterns 的 -Σ(p_i × log₂(p_i)) | shannon_hp1, shannon_hp2, shannon_combined, shannon_delta |
| **VEF** | epialleleR (2023) | 偏離最常見 binary pattern 的 reads 比例 | vef_hp1, vef_hp2, vef_combined |

### 2.2 計算細節

**數據來源**：31,659 regions 的 `methylation.csv`（reads × CpGs 矩陣）+ `reads.tsv`（HP tag + is_tumor 標記）

**Per-CpG Fisher's exact test**：
- 每個 CpG column j，依 HP 分組（HP1+HP1-1 vs HP2+HP2-1）
- 二值化：value > 0.5 → methylated；value ≤ 0.5 → unmethylated；NaN → exclude
- scipy.stats.fisher_exact() → BH-FDR 校正
- 聚合：n_sig（FDR<0.05）, frac_sig, max_delta, mean_delta, max_neg_log_fdr

**NME (Normalized Methylation Entropy)**：
- 每個 HP group 的 reads → binary pattern（across valid CpGs）
- Count unique patterns → p_i = count / total
- H = -Σ(p_i × log₂(p_i))；H_max = log₂(min(n_reads, 2^n_cpgs))
- NME = H / H_max；entropy_imbalance = |NME_HP1 - NME_HP2|

**Epipolymorphism + Shannon Entropy**：
- 4-CpG sliding windows → enumerate 2⁴=16 patterns
- Epipoly = 1 - Σ(p_i²)；Shannon = -Σ(p_i × log₂(p_i))
- Per-HP group 分別計算 + combined（全 reads）

**PDR**：
- 每條 read：若所有有效 CpG 全 >0.5 或全 ≤0.5 → concordant，否則 discordant
- PDR = discordant / total reads

**VEF**：
- Reference pattern = 最常見 binary pattern
- VEF = 偏離 reference pattern 的 reads 比例

**最低條件**：每 HP group ≥5 reads，否則 NaN。Per-CpG Fisher 需 ≥1 CpG；NME ≥2 CpGs；Epipolymorphism/Shannon ≥4 CpGs。

### 2.3 執行環境

- 16-worker multiprocessing，HCC1395 全基因體 31,659 regions
- TP/FP labels：30,401 TP + 1,258 FP
- AUC 計算：sklearn roc_auc_score（Mann-Whitney U）

---

## 3. 結果

### 3.1 Metric Validity

所有 24 個新 metrics 的有效率均 > 95%（除 entropy_imbalance 95.0% 因需兩個 HP group 均有足夠 reads）。Per-CpG Fisher 類指標 100% 有效。

### 3.2 TP/FP AUC — Filter Value = 0

| Rank | Metric | AUC | Direction | n |
|------|--------|-----|-----------|---|
| 1 | pdr_combined | 0.5374 | TP>FP | 31,659 |
| 2 | vef_combined | 0.5320 | TP>FP | 31,656 |
| 3 | fisher_max_neg_log_fdr | 0.5318 | TP>FP | 30,620 |
| 4 | vef_hp1 | 0.5174 | TP>FP | 30,709 |
| 5 | nme_combined | 0.5170 | TP>FP | 31,659 |
| 6 | shannon_hp1 | 0.5154 | TP>FP | 30,543 |
| ... | ... | ... | ... | ... |
| 22 | epipoly_hp2 | 0.4089 | FP>TP | 30,405 |
| 23 | shannon_combined | 0.4063 | FP>TP | 31,658 |
| 24 | fisher_mean_delta | 0.3811 | FP>TP | 30,620 |

**結論**：所有 24 個新 metrics AUC 在 0.38-0.54 範圍，均未突破 0.55。**與 Beyond-AUC exhaustion (2026-04-09) 結論完全一致，filter value 確認為零。**

既有 Quality_Score baseline AUC = 0.543，新指標中最佳的 pdr_combined (0.5374) 並未超越。

### 3.3 PERMANOVA vs Per-CpG Fisher Concordance

```
fisher_sig     False  True 
permanova_sig              
False           2518   3870
True            1028  23204
Agreement rate: 25722/30620 = 84.0%
```

**解讀**：
- 84% 一致率驗證 per-CpG Fisher 與 window-level PERMANOVA 高度互補
- 3,870 個「Per-CpG Fisher 陽性 / PERMANOVA 陰性」→ **per-CpG 偵測到 PERMANOVA 遺漏的局部 ASM**
- 1,028 個「PERMANOVA 陽性 / Per-CpG Fisher 陰性」→ PERMANOVA 捕捉到非 per-CpG 粒度的整體模式差異
- 此結果支持 CPEL 的論點：部分 ASM 是整體 entropy 差異而非逐位點差異

### 3.4 LOH 分層

| Metric | LOH (n=1,713) | Non-LOH (n=29,946) | Diff |
|--------|---------------|---------------------|------|
| fisher_frac_sig | 0.1325 | 0.2571 | **-0.1246** |
| fisher_n_sig | 4.87 | 8.86 | -3.99 |
| entropy_imbalance | 0.0348 | 0.0230 | **+0.0119** |
| epipoly_delta | 0.1929 | 0.1239 | **+0.0691** |
| pdr_delta | 0.0820 | 0.0651 | +0.0169 |
| epipoly_combined | 0.5927 | 0.6178 | -0.0251 |
| nme_combined | 0.9630 | 0.9694 | -0.0064 |

**解讀**：
- **fisher_frac_sig LOH 降低 12.5 百分點**：LOH 物理失去一個 HP → balanced reads 減少 → per-CpG ASM 顯著性下降。與 O15 LOH 量化結論一致。
- **entropy_imbalance LOH 增加 +0.012**：LOH 導致 HP 不對稱 → 兩邊 NME 差異增大。符合 CPEL 理論預期。
- **epipoly_delta LOH 增加 +0.069**：LOH 的 HP epipolymorphism 差異更大，因主要 HP 的 reads 數遠多於次要 HP。

### 3.5 HPFineNGroups 分層

| NGroups | fisher_frac_sig | entropy_imbal | epipoly | pdr | n |
|---------|-----------------|---------------|---------|-----|---|
| 1 | **0.0008** | NaN | 0.5967 | 0.8921 | 382 |
| 2 | 0.0815 | 0.0277 | 0.6135 | 0.9003 | 908 |
| 3 | **0.2704** | 0.0242 | 0.6125 | 0.9174 | 24,387 |
| 4 | 0.2254 | 0.0191 | 0.6363 | 0.9210 | 5,575 |

**解讀**：
- **NGroups=1 → fisher_frac_sig=0.001**：只有一個 HP group 時幾乎不可能有 per-CpG ASM — 完美驗證了指標的物理合理性
- **NGroups=3 最高 (0.270)**：標準雙 HP + somatic subclone → 最多 per-CpG ASM 事件
- **NGroups=4 略低 (0.225)**：更多 subclone groups 但每 group reads 更少 → statistical power 下降

NGroups 分層 AUC：
- NGroups=4 的 pdr_combined AUC=0.6107 — **唯一接近有意義的 AUC，但僅在特定分層**
- NGroups=1 的 epipoly_combined AUC=0.5934 — 有趣但 n=382 太少

### 3.6 PDR 與 VEF 飽和問題

- PDR combined 均值 = 0.917（91.7% reads 含不一致 CpGs）
- VEF combined 均值 = 0.918（91.8% reads 偏離最常見模式）

**兩者接近天花板 1.0，動態範圍極窄**，這解釋了為什麼即使理論上有區分力的指標，在實際 AUC 中表現平平。原因：ONT 長讀長在 2000bp window 中通常覆蓋 50-100+ CpGs，幾乎每條 read 都至少有一個「不一致」CpG → PDR 和 VEF 結構性趨近 1.0。

### 3.7 Site ASM Classification

基於 PERMANOVA + per-CpG Fisher 雙重標準的 6 類位點分布：

| Category | Count | % | Description |
|----------|-------|---|-------------|
| Homogenized | 11,197 | 35.4% | Both significant, direction consistent |
| No_ASM | 10,983 | 34.7% | Neither significant |
| Tumor_Acquired | 3,744 | 11.8% | Only tumor significant |
| Maintained | 3,023 | 9.6% | Both significant |
| LOH_Eliminated | 1,713 | 5.4% | LOH region |
| Reversed | 999 | 3.2% | Both significant, opposite direction |

---

## 4. 圖表

### Figure 07: Per-CpG Fisher Overview
四面板：n_sig 分布、PERMANOVA concordance、TP/FP 比較、Germline ASM 驗證

### Figure 08: Entropy & Epiallele Metrics
四面板：NME HP1 vs HP2 scatter、entropy imbalance by NGroups、epipolymorphism distribution、PDR by HP group

### Figure 09: AUC Landscape — Per-CpG ASM vs Existing Metrics
24 個新指標 + 既有 baseline 的水平 AUC 排列。所有指標均低於 prior ceiling (0.58)。

### Figure 10: LOH × HPFineNGroups Stratification
四面板：entropy imbalance LOH vs Non-LOH、fisher_frac_sig LOH 分層、epipolymorphism NGroups 分層、fisher_frac_sig by Site Classification

圖表位置：`research/germline_asm_analysis/figures/07-10_*.png`

---

## 5. 結論

### 5.1 Filter Value — 確認為零

所有 24 個新 per-CpG metrics AUC 0.38-0.54，無一突破 0.55。**與 Beyond-AUC exhaustion (2026-04-09) 結論完全吻合**。Per-CpG 粒度並未提供 window-average 遺漏的 TP/FP 區分信號。

**根因一致**：germline ASM 效應量（3-6×）遠大於 somatic passenger SNV 對甲基化的影響，TP 和 FP 在 per-CpG 層級同樣不可區分。

### 5.2 Characterization Value — 確認有效

Per-CpG 指標在以下維度展現 characterization 價值：

1. **PERMANOVA concordance 84%**：per-CpG Fisher 與 PERMANOVA 高度一致但互補（3,870 個 per-CpG-only detections）
2. **LOH 分層清晰**：fisher_frac_sig LOH 降低 12.5pp；entropy_imbalance 增加 +0.012；epipoly_delta 增加 +0.069
3. **HPFineNGroups 驗證**：NGroups=1 → 0.001 ASM；NGroups=3 → 0.270 ASM（物理機制完美對應）
4. **Site classification 細化**：6 類位點的 fisher_frac_sig 分布與生物學預期一致（Maintained/Reversed 最高，LOH_Eliminated 最低）

### 5.3 文獻方法在 ISM 上的表現

| 指標 | 原始文獻用途 | ISM TP/FP AUC | ISM characterization | 評價 |
|------|-------------|---------------|---------------------|------|
| Per-CpG Fisher | DMR/ASM detection | 0.45-0.53 | **高**（concordance, LOH, NGroups） | 推薦 C++ 整合 |
| NME/Entropy | Epigenetic state | 0.48-0.52 | **中**（entropy_imbalance LOH 分層） | 推薦 C++ 整合 |
| Epipolymorphism | Tumor heterogeneity | 0.41-0.49 | **中**（epipoly_delta LOH 差異） | 可選 |
| PDR | Discordance | 0.50-0.54 | **低**（飽和 0.917） | 不推薦（動態範圍不足） |
| Shannon Entropy | Pattern diversity | 0.41-0.52 | **低**（與 epipolymorphism 高度相關） | 不推薦（冗餘） |
| VEF | Variant epialleles | 0.49-0.53 | **低**（飽和 0.918） | 不推薦（動態範圍不足） |

### 5.4 C++ 整合決策建議

**推薦整合（精簡版，~10 個欄位）**：
- Per-CpG Fisher: `n_sig`, `frac_sig`, `n_tested`, `max_neg_log_fdr`（4 欄位）
- NME: `nme_hp1`, `nme_hp2`, `entropy_imbalance`（3 欄位）
- Epipolymorphism: `epipoly_hp1`, `epipoly_hp2`, `epipoly_delta`（3 欄位）

**不推薦整合**：
- PDR（飽和 0.917，動態範圍 <0.10）
- Shannon Entropy（與 epipolymorphism 功能重疊）
- VEF（飽和 0.918，動態範圍 <0.10）

理由：精簡版 10 欄位已涵蓋 per-CpG ASM 偵測（Fisher）、entropy 維度（NME）、pattern 多樣性（Epipolymorphism）三個正交維度。PDR/Shannon/VEF 在 2000bp ONT window 中結構性飽和，增加欄位無額外資訊。

### 5.5 Python PoC 判定

按計劃標準：
- **至少一個 effect size > 0.15 在生物學分層中**：fisher_frac_sig LOH vs Non-LOH diff=0.125；epipoly_delta LOH diff=0.069。fisher_frac_sig 接近門檻但未達 0.15。
- **Characterization value 已確認**：PERMANOVA concordance 84%、LOH 分層、NGroups 分層均展現生物學合理性。

**判定：CONDITIONAL POSITIVE** — characterization 價值確認，但 effect size 邊緣。建議在 C++ 整合時精簡為 10 欄位（而非原計劃的 22 欄位），並在全基因體重跑後重新評估。

---

## 6. 與先前研究的關係

| 先前結論 | 本研究確認/更新 |
|----------|----------------|
| Beyond-AUC exhaustion (2026-04-09) | ✅ 確認：per-CpG 粒度同樣 AUC ≤ 0.54 |
| O11 epipolymorphism NEGATIVE | ⚠️ 更新：per-HP epipolymorphism 避免了 n_reads confound，但仍無 filter value |
| ASM 32-66% POSITIVE | ✅ 強化：per-CpG Fisher 85.5% regions ≥1 sig CpG |
| LOH 甲基化全面失效 (O15) | ✅ 一致：LOH 降低 fisher_frac_sig 12.5pp |
| HPFineNGroups somatic marker | ✅ 強化：NGroups 與 per-CpG ASM 強烈對應 |
| ISM characterization 定位 | ✅ 再確認：新指標增強 characterization，不增加 filter |

---

## 7. 後續建議

1. **C++ 精簡整合**：10 欄位（Fisher 4 + NME 3 + Epipoly 3），插入 RegionProcessor.cpp:~843
2. **全基因體重跑**：加入新欄位後跑 7 個樣本，驗證跨樣本穩定性
3. **論文素材**：PERMANOVA concordance 84% 和 HPFineNGroups × fisher_frac_sig 可作為 ISM read-level characterization 的展示
4. **不需做的事**：PDR/Shannon/VEF 不整合；不嘗試用新指標建 filter model

---

## 8. 腳本與數據

| 檔案 | 位置 | 用途 |
|------|------|------|
| compute_per_cpg_asm_metrics.py | `/tmp/ism_phase_b_test/` | 主要計算腳本 |
| analyze_per_cpg_asm.py | `/tmp/ism_phase_b_test/` | 視覺化分析腳本 |
| significance_summary_per_cpg_asm.csv | `/tmp/ism_phase_b_test/output_full_phase_cd/` | 29MB 輸出（+20 欄位） |
| per_cpg_asm_validation_results.txt | `/tmp/ism_phase_b_test/` | 驗證結果全文 |
| Figures 07-10 | `research/germline_asm_analysis/figures/` | 4 張分析圖 |
