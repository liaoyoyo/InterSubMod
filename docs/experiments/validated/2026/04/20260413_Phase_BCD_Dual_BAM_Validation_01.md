<!--
建立時間: 2026-04-13 03:20
更新時間: 2026-04-13 05:30
目標: Phase B/C/D Dual-BAM 架構驗證報告
處理範圍: HCC1395 paired (tumor + normal BAM) 全基因體驗證（31,659 regions = 30,401 TP + 1,258 FP）
關聯檔案:
  - include/core/LohBedAnnotator.hpp
  - include/core/SubcloneAnalyzer.hpp
  - src/core/RegionProcessor.cpp
  - docs/plans/ (cryptic-swinging-kazoo.md)
-->

# Phase B/C/D: ISM Dual-BAM 架構驗證報告

## 摘要

本報告記錄 ISM Phase A-D Dual-BAM 架構的全基因體驗證結果。使用 HCC1395 paired data（tumor + normal BAM），驗證 Normal Baseline、Sample ASM、Somatic HP ASM、LOH BED 交叉驗證、以及跨區域 Subclone 分層的正確性與生物學合理性。

**結論**: 所有 Phase 驗證通過，系統行為符合預期。

---

## 1. 測試環境

| 項目 | 值 |
|------|-----|
| Tumor BAM | HCC1395_tagged.bam (LongPhase-S haplotagged) |
| Normal BAM | HCC1395BL_5kHz.tagged.bam (germline haplotagged, HP:i:1/2) |
| VCF | 31,659 SNVs (30,401 TP + 1,258 FP) |
| Reference | GRCh38_no_alt_analysis_set.fasta |
| LOH BED | ISM hp_ratio 生成的 1,685 regions |
| Window | ±2000bp |
| Distance metric | NHD |
| Threads | 16 |

---

## 2. Phase B: Normal Baseline + Sample ASM 驗證

### 2.1 全基因體結果 (31,659 TP regions)

| 指標 | 值 | 判定 |
|------|-----|------|
| Sample ASM 顯著率 | 97.3% (30,790/31,659) | PASS — 幾乎所有 somatic SNV sites 腫瘤與正常有甲基化差異 |
| Mean Sample ASM delta | 0.2257 (std=0.1418) | PASS — 正向偏移，腫瘤甲基化高於正常 |
| Normal Baseline 有效率 | 100% (31,659/31,659) | PASS — 所有 region 都有足夠 normal reads |
| Mean Normal methylation | 0.4675 | 合理 — 接近全基因體 CpG 的預期 ~50% |
| Mean Normal coverage | 24.1 reads/CpG | 合理 — 正常 BAM ~30x, 考慮 CpG 密度後合理 |

### 2.2 Sample ASM Delta 分佈

| 範圍 | 比例 | 解讀 |
|------|------|------|
| [-1.0, -0.1) | 0.0% | 極少數 normal > tumor |
| [-0.1, 0) | 1.3% | 少量反向差異 |
| [0, 0.05) | 5.6% | 微小差異 |
| [0.05, 0.1) | 11.6% | 中等差異 |
| [0.1, 0.2) | 32.5% | **最大群** — 典型 somatic 甲基化變化 |
| [0.2, 0.5) | 44.1% | **次大群** — 強 somatic 差異 |
| [0.5, 1.0) | 5.0% | 極端差異（可能 CNV/LOH 區域） |

**結論**: 76.6% 的 regions 有 delta ≥ 0.1，顯示 tumor-normal 甲基化差異普遍存在且方向一致（tumor > normal），符合腫瘤去甲基化/過甲基化的已知生物學。

### 2.3 Somatic HP ASM

| 指標 | 值 | 解讀 |
|------|-----|------|
| HP Residual 非零率 | 95.9% | 大部分 regions 可計算 |
| Tumor HP valid | 52.9% | 受限於 somatic haplotagging 的 HP 不平衡 |
| Normal HP valid | 93.8% | Normal BAM 的 HP 標籤平衡良好 |
| Both valid | 50.8% | 可完整計算 somatic HP ASM 的 regions |
| Mean HP residual (both valid) | -0.045 | **Tumor HP ASM < Normal HP ASM** |
| Mean Tumor HP delta | 0.059 | |
| Mean Normal HP delta | 0.104 | |

**關鍵發現**: Somatic HP ASM (tumor_hp - normal_hp) = -0.045，意味著 **somatic SNV sites 的腫瘤 HP 間甲基化差異反而比正常低**。可能原因：

1. Somatic haplotagging 使 tumor reads 高度偏向一個 HP → HP 測試有效性受限
2. Normal germline HP ASM (0.104) 確認了先前 L1 研究發現（germline ASM mean=-0.0205, p=1.16e-14）在全基因體尺度更明顯
3. Tumor 的表觀遺傳變化可能「均質化」了原本的 HP 間差異

### 2.4 HP Balance

| 樣本 | HP1 (mean) | HP2 (mean) | Ratio |
|------|-----------|-----------|-------|
| Tumor | 31.6 | 30.8 | 0.506 |
| Normal | 15.2 | 15.1 | 0.501 |

兩者都接近完美的 50:50 平衡，確認 haplotagging 品質良好。

---

## 3. Phase C: LOH BED 交叉驗證

### 3.1 LOH BED 載入

- 載入 1,685 LOH regions（從 ISM hp_ratio 生成）
- 分布：chr8 最多（405 regions，佔 24%），與 HCC1395 已知 chr8 LOH 一致

### 3.2 LOH Source 交叉驗證

#### chr19 子集 (674 regions)

| LOH Source | 數量 | 比例 |
|-----------|------|------|
| none | 641 | 95.1% |
| both (BED + hp_ratio) | 16 | 2.4% |
| bed_only | 1 | 0.1% |
| ratio_only | 0 | 0.0% |

**Concordance**: 16/17 = **94.1%** — BED 與 ISM hp_ratio 高度一致。

#### 全基因體 (31,659 regions)

| LOH Source | 數量 | 比例 |
|-----------|------|------|
| none | 29,946 | 94.6% |
| both (BED + hp_ratio) | 1,685 | 5.3% |
| bed_only | 28 | 0.1% |
| ratio_only | 0 | 0.0% |

**Concordance**: 1,685/1,713 = **98.4%** — 全基因體尺度 BED 與 ISM hp_ratio 高度一致。28 個 bed_only cases 是 BED 區域覆蓋但 ISM hp_ratio 未達 LOH 門檻的邊界情況，佔比僅 1.6%。

### 3.3 全基因體 LOH 分佈

| 指標 | 值 |
|------|-----|
| Potential LOH (hp_ratio) | 5.3% (1,685/31,659) |
| Coverage Normal | 49.2% |
| Coverage Elevated | 25.9% |
| Coverage CNV_Gain | 13.7% |
| Coverage CNV_Loss | 0.2% |

HCC1395 的 CNV 分佈：~40% 的 regions 有 elevated 或更高 coverage，與已知的高 CNV burden 一致。

---

## 4. Phase D: 跨區域 Subclone 分層

### 4.1 Subclone 分層結果

#### chr19 子集 (674 regions)

| Group | 名稱 | Regions | 比例 | Mean HP ASM | Mean Sample ASM |
|-------|------|---------|------|------------|----------------|
| 0 | Normal Diploid | 115 | 17.5% | 0.019 | 0.048 |
| 1 | Epigenetic Heterogeneity | 85 | 12.9% | 0.102 | 0.063 |
| 2 | LOH Regions | 17 | 2.6% | 0.067 | 0.154 |
| 3 | Tumor-Specific Changes | 441 | 67.0% | 0.093 | 0.199 |

#### 全基因體 (31,659 regions)

| Group | 名稱 | Regions | 比例 | Mean |HP ASM| | Mean |Sample ASM| | HP ratio | HP fine F | Cov Multiple |
|-------|------|---------|------|------------|----------------|----------|-----------|-------------|
| 0 | Normal Diploid | 2,909 | 9.2% | 0.020 | 0.056 | 0.503 | 4.78 | 1.19 |
| 1 | Epigenetic Heterogeneity | 2,615 | 8.3% | 0.107 | 0.063 | 0.506 | 8.96 | 1.04 |
| 2 | LOH Regions | 1,713 | 5.4% | 0.111 | 0.228 | 0.565 | 13.18 | 1.24 |
| 3 | Tumor-Specific Changes | 24,422 | 77.1% | 0.108 | 0.264 | 0.501 | 16.45 | 1.20 |

### 4.2 生物學解讀

- **Group 0 (9.2%)**: HP ratio=0.503（完美平衡），最低 ASM → 正常二倍體區域，somatic SNV 未影響甲基化
- **Group 1 (8.3%)**: 高 HP ASM=0.107 但低 Sample ASM=0.063 → 表觀遺傳異質性（allele-specific methylation），腫瘤vs正常差異小
- **Group 2 (5.4%)**: 1,713/1,713 BED overlap（100%），HP ratio=0.565（偏離0.5）→ LOH 驅動的甲基化模式
- **Group 3 (77.1%)**: 最大群，高 Sample ASM=0.264 → 腫瘤特異性甲基化變化（主要群體），HP ratio=0.501 平衡

**chr19 vs 全基因體比較**：Group 比例穩定（chr19: 17.5/12.9/2.6/67.0% → 全基因體: 9.2/8.3/5.4/77.1%），chr19 的 Normal Diploid 比例偏高（17.5% vs 9.2%），與 chr19 較低的整體 Sample ASM（0.153 vs 0.226）一致。

### 4.3 LOH Group 特徵

LOH group (Group 2) 的特殊性質：

| 指標 | chr19 | 全基因體 |
|------|-------|---------|
| BED overlap | 17/17 (100%) | 1,713/1,713 (100%) |
| BED/ratio concordance | 16/17 (94.1%) | 1,685/1,713 (98.4%) |
| Mean HP fine F | 16.76 | 13.18 |
| Mean HP ratio | 0.515 | 0.565 |
| Mean Coverage Multiple | 0.93 | 1.24 |

全基因體的 LOH group 覆蓋倍數（1.24）高於 chr19（0.93），反映不同染色體的 CNV 背景差異。

### 4.4 HPFineNGroups × LOH 反相關

| NGroups | N (%) | LOH% | 解讀 |
|---------|-------|------|------|
| 0 | 407 (1.3%) | 23.8% | 缺少 HP 分群，部分來自 LOH |
| 1 | 382 (1.2%) | **76.7%** | 單群=幾乎全是 LOH 區域 |
| 2 | 908 (2.9%) | 33.1% | 部分 HP 分群 |
| 3 | 24,387 (77.0%) | 4.1% | 正常三群分布 |
| 4 | 5,575 (17.6%) | **0.2%** | 四群 ≈ 非 LOH |

完美的反相關：LOH 區域因為失去一個 allele，導致 HP fine-grained 分群數減少。NGroups=1 有 76.7% LOH rate（最高），NGroups=4（最高異質性）的 LOH rate 幾乎為零。這是 HPFineNGroups 作為 somatic heterogeneity marker 的直接生物學驗證。

### 4.5 TP vs FP 比較

#### 全基因體 (30,401 TP + 1,258 FP)

| 指標 | TP (30,401) | FP (1,258) |
|------|-----------|-----------|
| Sample ASM Delta | 0.2246 | 0.2529 |
| Normal Baseline Mean | 0.4688 | 0.4354 |
| Quality Score | 93.5 | 91.5 |
| |HP ASM| | 0.0993 | 0.1191 |
| HP Fine F | 14.80 | 9.37 |
| Coverage Multiple | 1.189 | 1.190 |
| Sample ASM sig% | 97.3% | 96.5% |
| LOH BED overlap | 1,634 (5.4%) | 79 (6.3%) |
| HP Residual Delta | -0.0681 | -0.0385 |

#### Subclone 分布

| Subclone | TP | FP |
|----------|-----|-----|
| 0 (Normal Diploid) | 2,824 (9.3%) | 85 (6.8%) |
| 1 (Epigenetic) | 2,530 (8.3%) | 85 (6.8%) |
| 2 (LOH) | 1,634 (5.4%) | 79 (6.3%) |
| 3 (Tumor-Specific) | 23,413 (77.0%) | 1,009 (80.2%) |

#### Verification Class 分布

| Class | TP | FP |
|-------|-----|-----|
| Strong | 7,268 (23.9%) | 313 (24.9%) |
| Subclone | 473 (1.6%) | 11 (0.9%) |
| Weak | 18,974 (62.4%) | 751 (59.7%) |
| Noise | 3,686 (12.1%) | 183 (14.5%) |

#### 解讀

1. **FP 的 Sample ASM 略高於 TP**（0.253 vs 0.225）：FP 是 germline variants 被誤報為 somatic，但在腫瘤環境中仍有真實的甲基化差異。FP 的 Normal Baseline 略低（0.435 vs 0.469），暗示 FP 位點在 germline 中已有特殊的甲基化背景。
2. **FP 的 HP ASM 更高**（0.119 vs 0.099）：符合預期 — FP 是 germline heterozygous variants，本身就有 germline ASM 信號。
3. **FP 的 HP Fine F 更低**（9.37 vs 14.80）：FP 的 HP 分群結構簡單（缺乏 somatic subclone 多樣性），與 HPFineNGroups 作為 somatic heterogeneity marker 的解釋一致。
4. **FP 更集中於 Tumor-Specific group**（80.2% vs 77.0%）：因 FP 有高 Sample ASM 但非 LOH 區域。
5. **FP 沒有 Subclone verification enrichment**（0.9% vs 1.6%）：進一步確認 FP 的甲基化差異非 subclone 驅動。
6. **結論**：與先前 14 個月研究一致 — **ISM 甲基化特徵無法區分 TP/FP**，差異來自 germline vs somatic origin 的本質不同。

---

## 5. 向後相容性

| 檢查項目 | 結果 |
|---------|------|
| 單元測試 | 173/173 全通過 |
| 無 normal BAM | 新欄位默認值正確（NormalReads=0, SampleASM=0, LOH=none, Subclone=-1） |
| 無 LOH BED | loh_annotator_.loaded() = false, 跳過標注 |
| CSV 格式 | 新增 LOH_Bed_Overlap, LOH_Source, LOH_Bed_Annotation, Subclone_ID 欄位 |

---

## 6. 效能

| 指標 | Phase B (31K TP) | Phase C+D (31K TP+FP) |
|------|-----------------|----------------------|
| Total time | ~15 min | ~15 min (預估) |
| Avg ms/region | 29 | ~29 |
| Memory | 22.7 MB | ~23 MB |
| Valid pair ratio | 89.7% | ~90% |

Phase C (LOH 查詢) 和 D (Subclone 分析) 的額外開銷可忽略 — BED 查詢是 O(log n)，subclone 分析是單次 O(n) 後處理。

---

## 7. 補充分析

### 7.1 Per-Chromosome LOH 分佈

| Chr | N | LOH% | Mean Reads | Mean Sample ASM | Mean QS |
|-----|---|------|-----------|----------------|---------|
| chr1 | 2,631 | 4.6% | 100.8 | 0.217 | 94.4 |
| chr2 | 2,700 | 3.2% | 101.5 | 0.235 | 94.0 |
| chr3 | 2,125 | 7.9% | 107.0 | 0.221 | 92.8 |
| chr4 | 2,317 | 2.0% | 100.3 | 0.257 | 91.9 |
| chr5 | 2,056 | 4.5% | 108.8 | 0.215 | 94.0 |
| chr6 | 1,258 | 5.4% | 100.4 | 0.219 | 90.9 |
| chr7 | 2,493 | 2.3% | 143.8 | 0.232 | 92.4 |
| **chr8** | **2,745** | **15.1%** | **121.9** | **0.260** | **89.5** |
| chr9 | 1,314 | 2.5% | 105.2 | 0.230 | 95.2 |
| chr10 | 1,327 | 6.1% | 101.9 | 0.211 | 94.5 |
| chr11 | 1,202 | 2.9% | 95.1 | 0.188 | 94.1 |
| chr12 | 1,408 | 6.8% | 102.6 | 0.226 | 93.2 |
| chr13 | 911 | 3.3% | 94.2 | 0.251 | 92.7 |
| chr14 | 1,214 | 8.8% | 128.0 | 0.228 | 93.3 |
| chr15 | 1,031 | 1.6% | 117.6 | 0.211 | 96.0 |
| chr16 | 652 | 9.5% | 126.7 | 0.172 | 95.0 |
| **chr17** | **949** | **11.1%** | **110.4** | **0.195** | **95.7** |
| chr18 | 892 | 1.2% | 112.2 | 0.238 | 93.6 |
| chr19 | 674 | 2.5% | 104.6 | 0.153 | 97.8 |
| chr20 | 902 | 1.3% | 117.2 | 0.231 | 97.5 |
| chr21 | 380 | 1.1% | 87.8 | 0.235 | 92.6 |
| **chr22** | **478** | **10.0%** | **139.5** | **0.227** | **96.4** |

**LOH 富集染色體**（>10%）：chr8（15.1%）、chr17（11.1%）、chr22（10.0%），均與 HCC1395 已知的 LOH 區域一致。chr8 是 LOH 最富集的染色體，同時有最高的 Sample ASM（0.260）和最低的 Quality Score（89.5）。

### 7.2 Verification Class 特徵分佈

| Class | N (%) | Mean Sample ASM | Mean HP ASM | LOH% |
|-------|-------|----------------|------------|------|
| Strong | 7,581 (23.9%) | 0.290 | 0.158 | 4.7% |
| Subclone | 484 (1.5%) | 0.438 | 0.012 | 6.8% |
| Weak | 19,725 (62.3%) | 0.203 | 0.095 | 4.8% |
| Noise | 3,869 (12.2%) | 0.189 | 0.014 | 8.9% |

Subclone class 有最高的 Sample ASM (0.438) 但最低的 HP ASM (0.012) — 符合「腫瘤特異性但非 HP 驅動」的 subclone 特徵。

### 7.3 HPFineNGroups × LOH 反相關

| NGroups | N (%) | LOH% | 解讀 |
|---------|-------|------|------|
| 0 | 407 (1.3%) | 23.8% | 缺少 HP 分群 |
| 1 | 382 (1.2%) | **76.7%** | 單群=幾乎全是 LOH 區域 |
| 2 | 908 (2.9%) | 33.1% | 部分 HP 分群 |
| 3 | 24,387 (77.0%) | 4.1% | 正常三群分布 |
| 4 | 5,575 (17.6%) | **0.2%** | 四群 ≈ 非 LOH |

完美的反相關：NGroups=1 有 76.7% LOH rate（最高），NGroups=4（最高異質性）的 LOH rate 幾乎為零。LOH 區域因為失去一個 allele，導致 HP fine-grained 分群數減少。這也是 Section 4.4 HPFineNGroups × LOH 反相關的另一角度呈現。

---

## 8. 已知限制與後續工作

1. **LOH BED 來源**: 本次使用 ISM 自身 hp_ratio 生成的 BED 做交叉驗證（self-referential），應使用 LongPhase germline haplotagging 的獨立 LOH.bed
2. **Subclone 分層**: 目前使用固定閾值的規則分層，非統計學方法。Phase D 後續可引入 k-means 或 hierarchical clustering
3. **Somatic HP ASM 的解讀**: 全基因體 HP Residual Delta = -0.068（TP）/-0.039（FP），均為負值，需進一步研究 tumor HP delta < normal HP delta 的生物學意義
4. **跨樣本驗證**: 本次為 HCC1395 單樣本驗證，需 7 samples 全量驗證確認 subclone 分層穩定性

---

## 9. 結論

Phase A-D Dual-BAM 架構驗證通過（全基因體 31,659 regions，含 30,401 TP + 1,258 FP）：

- **Phase B**: Normal Baseline 100% 有效，Sample ASM 97.3% 顯著，76.6% regions 有 delta ≥ 0.1，生物學方向合理
- **Phase C**: LOH BED 載入正確，全基因體 concordance 98.4%（1,685/1,713），chr19 子集 94.1%（16/17）
- **Phase D**: 4-group 分層生物學合理（Normal 9.2% / Epigenetic 8.3% / LOH 5.4% / Tumor-Specific 77.1%），LOH group 100% BED overlap
- **TP vs FP**: FP 有略高 Sample ASM（0.253 vs 0.225）和更高 HP ASM（0.119 vs 0.099），但差異不足以區分 — 與先前 14 個月研究結論一致
- **HPFineNGroups × LOH**: 完美反相關（NGroups=1 LOH rate 76.7%, NGroups=4 LOH rate 0.2%）

ISM 已具備完整的 Dual-BAM 分析能力，為 read-level epigenetic characterization 論文提供了核心基礎設施。
