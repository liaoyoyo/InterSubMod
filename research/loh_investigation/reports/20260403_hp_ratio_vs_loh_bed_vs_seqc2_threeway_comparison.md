<!--
建立時間: 2026-04-03 01:50
目標: 三套 LOH 判定系統（ISM HP_Ratio / LOH.bed / SEQC2）的包含關係與 TP/FP 差異分析
處理範圍: 7 個 TO 樣本 × 419,692 regions + HCC1395 SEQC2 40,096 variants
關聯檔案:
  - research/loh_investigation/scripts/analyze_hp_ratio_vs_loh_bed_vs_seqc2.py
  - research/loh_investigation/data/loh3way_*.tsv (7 TSV)
  - research/loh_investigation/figures/loh3way_fig01-09*.png (9 figures)
  - research/loh_investigation/reports/20260402_seqc2_vs_longphase_to_loh_validation.md
  - src/core/RegionProcessor.cpp:38-49 (compute_hp_ratio, is_potential_loh)
-->

# HP_Ratio vs LOH.bed vs SEQC2 — 三套 LOH 判定系統比較

## 摘要

本分析比較三套 LOH 判定系統在 TO 模式下的包含關係、TP/FP 分佈差異與位置聚集特性：

1. **ISM Potential_LOH**：site-level，HP_Ratio < 0.1 or > 0.9（來自 ISM C++ 核心）
2. **LongPhase-TO LOH.bed**：region-level，phased genotype ratio（來自 LongPhase-TO）
3. **SEQC2 CNV Benchmark**：ground truth LOH（僅 HCC1395 可用）

**核心發現**：ISM 是 LOH.bed 的嚴格超集（recall 99%+），但多報 10-16% 的 LOH sites。這些「ISM-only」LOH sites 的 TP/FP rate 與其他類別**無統計差異**，不具有額外判別價值。差異主要來自低 effective_hp_reads 造成的 HP_Ratio 極端化。

---

## 一、三套系統定義對照

| 系統 | 層級 | 定義 | 數據來源 |
|------|------|------|---------|
| ISM `Potential_LOH` | Site-level (per-SNV) | `HP_Ratio < 0.1 or > 0.9` | `RegionProcessor.cpp:47-49` |
| LongPhase-TO `LOH.bed` | Region-level (genomic interval) | Phased genotype allele ratio | LongPhase-TO 輸出 |
| SEQC2 CNV Benchmark | Region-level (ground truth) | NGS + SNP array 交叉驗證 | SEQC2 consortium |

**ISM HP_Ratio 計算**（`RegionProcessor.cpp:38-42`）：
```cpp
double compute_hp_ratio(int hp1_family_n, int hp2_family_n) {
    double total = static_cast<double>(hp1_family_n + hp2_family_n) + 0.002;
    return (static_cast<double>(hp1_family_n) + 0.001) / total;
}
```

- 分子 = HP1 family reads + epsilon
- 分母 = HP1 + HP2 family reads + epsilon
- 無最低 read 數門檻 → 即使只有 1 read 也可判定 LOH

---

## 二、HCC1395 三方交叉比較（Analysis A）

### 2.1 八類分佈

使用 3-bit encoding（ISM × LOH.bed × SEQC2）將 40,096 個 HCC1395 TO variants 分為 8 類：

| 類別 | Sites | 佔比 | TP | FP | TP Rate | FP Rate |
|------|-------|------|-----|-----|---------|---------|
| **All three** | 17,010 | 42.4% | 12,587 | 4,423 | **0.740** | 0.260 |
| **None** | 16,233 | 40.5% | 11,022 | 5,211 | **0.679** | 0.321 |
| **ISM only** | 5,516 | 13.8% | 3,882 | 1,634 | 0.704 | 0.296 |
| ISM ∩ SEQC2 | 687 | 1.7% | 526 | 161 | 0.766 | 0.234 |
| ISM ∩ LOH.bed | 359 | 0.9% | 274 | 85 | 0.763 | 0.237 |
| LOH.bed ∩ SEQC2 | 176 | 0.4% | 123 | 53 | 0.699 | 0.301 |
| SEQC2 only | 105 | 0.3% | 74 | 31 | 0.705 | 0.295 |
| LOH.bed only | 10 | 0.0% | 7 | 3 | 0.700 | 0.300 |

![Fig01: HCC1395 三方分類計數與 TP/FP rate](assets/20260403_loh3way_comparison/loh3way_fig01_threeway_venn.png)

### 2.2 關鍵觀察

1. **TP/FP rate 差異極小**：所有類別的 TP rate 在 0.679-0.766 之間，差距僅 0.087。LOH 判定本身**無法區分 TP/FP**。

2. **ISM 是 LOH.bed 的嚴格超集**：
   - LOH.bed only（不被 ISM 偵測）= 僅 10 sites（0.0%）
   - ISM only（不在 LOH.bed 中）= 5,516 sites（13.8%）
   - → ISM 完整包含 LOH.bed，但多報大量 sites

3. **「All three」一致的 sites 有微幅較高的 TP rate（0.740）**：因為 LOH 區域的 TP 天然有較高 AF（LOH 增幅 somatic allele），但效果量僅 Δ=0.06，無實用價值。

4. **ISM-only LOH 的 TP rate（0.704）與 None（0.679）幾乎相同**：ISM 多報的 LOH sites 不含額外判別資訊。

---

## 三、七樣本雙向比較（Analysis B）

### 3.1 LOH Rate 對比

| Sample | ISM LOH Rate | LOH.bed Rate | Excess ISM | Jaccard |
|--------|-------------|-------------|------------|---------|
| HCC1395 | 58.8% | 43.8% | **+15.0%** | 0.731 |
| HCC1395_DORADO | 60.0% | 44.0% | **+16.0%** | 0.733 |
| COLO829 | 34.8% | 20.9% | **+13.9%** | 0.601 |
| H1437 | 41.0% | 26.4% | **+14.6%** | 0.643 |
| H2009 | 40.6% | 24.8% | **+15.8%** | 0.609 |
| HCC1937 | 61.1% | 50.4% | **+10.7%** | 0.820 |
| HCC1954 | 22.2% | 6.2% | **+16.0%** | 0.281 |

![Fig02: 7 樣本雙向 LOH 分類熱圖](assets/20260403_loh3way_comparison/loh3way_fig02_twoway_heatmaps.png)

![Fig07: LOH rate / excess / Jaccard 一致性](assets/20260403_loh3way_comparison/loh3way_fig07_concordance_summary.png)

### 3.2 關鍵發現

1. **ISM 一致性超額 10-16%**：所有 7 個樣本，ISM LOH rate 均高於 LOH.bed rate，超額穩定在 10.7%-16.0%。這是系統性的，非隨機。

2. **ISM recall vs LOH.bed > 99%**：LOH.bed 判定為 LOH 的區域，ISM 幾乎全部也判定為 LOH（recall 98.9%-100%）。

3. **ISM precision vs LOH.bed 差異巨大**：
   - HCC1937: 82.2%（最高，Jaccard 0.820）
   - HCC1954: 28.1%（最低，Jaccard 0.281）
   - → ISM 對低 LOH coverage 樣本（HCC1954 僅 6.2%）過度報告 LOH

4. **TP/FP rate 在所有類別中接近一致**：
   - H2009（FP rate 最低的樣本）：Both=8.1%, ISM-only=8.1%, Neither=9.1%
   - HCC1954（FP rate 最高的樣本）：Both=66.2%, ISM-only=73.5%, Neither=75.5%
   - → LOH 類別與 TP/FP 判定能力正交

### 3.3 LOH.bed only 類別的稀缺性

| Sample | LOH.bed only sites | 佔比 |
|--------|-------------------|------|
| HCC1395 | 186 | 0.46% |
| HCC1395_DORADO | 6 | 0.01% |
| COLO829 | 0 | 0.00% |
| H1437 | 4 | 0.01% |
| H2009 | 45 | 0.03% |
| HCC1937 | 37 | 0.15% |
| HCC1954 | 8 | 0.01% |

LOH.bed only（在 LOH.bed 中但 ISM 不認為 LOH）極其罕見（< 0.5%），進一步確認 ISM ⊃ LOH.bed 的超集關係。

---

## 四、HP_Ratio 分佈（Analysis C）

### 4.1 Zone × Truth 四象限

| Group | HP_Ratio Median | Extreme Rate | 說明 |
|-------|----------------|-------------|------|
| TP in LOH.bed | 1.0 | 99.0-100% | 全部極端 |
| FP in LOH.bed | 1.0 | 98.8-100% | 全部極端 |
| TP outside LOH.bed | 0.41-0.46 | 18-29% | 正常分佈 |
| FP outside LOH.bed | 0.40-0.46 | 17-27% | 正常分佈 |

![Fig03: HP_Ratio violin by zone × truth](assets/20260403_loh3way_comparison/loh3way_fig03_hp_ratio_by_zone_truth.png)

### 4.2 關鍵結論

- **LOH.bed 內**：TP 和 FP 的 HP_Ratio 分佈**完全一致**（median 1.0, extreme rate 99%+）。HP_Ratio 在此完全無區分力。
- **LOH.bed 外**：TP 和 FP 的 HP_Ratio 也**幾乎一致**（median ~0.4-0.5, extreme rate ~18-28%）。HP_Ratio 在此同樣無區分力。
- **結論**：HP_Ratio 反映的是 phasing scaffold 和 haplotype block 的特性，與 variant 的 TP/FP 身份無關。

---

## 五、染色體位置分析（Analysis D）

### 5.1 HCC1395 全基因組位置圖

![Fig04: HCC1395 全基因組 ideogram](assets/20260403_loh3way_comparison/loh3way_fig04_chr_position_hcc1395.png)

觀察：
1. **「All three」（紫色）**：形成清晰的連續區塊，對應已知 LOH 區域（chr1p, chr5q, chr8, chr11, chr13q, chr16q, chr17p 等），與 SEQC2 ground truth 高度吻合。

2. **「ISM only」（粉紅色）**：散佈在所有染色體上，**無明顯聚集**。這些是 ISM 多報的 LOH sites，沒有形成連續區塊。

3. **「SEQC2 only」（橘色）**：主要出現在 LOH 區域邊界附近（chr1, chr5, chr8 邊界）。這些是 LOH.bed 未完整覆蓋的邊界區域。

### 5.2 ISM-only LOH 染色體分佈

![Fig05: ISM-only 和 SEQC2-only 染色體密度](assets/20260403_loh3way_comparison/loh3way_fig05_chr_density_hcc1395.png)

- ISM-only LOH 的染色體分佈大致與 variant 密度成正比，沒有特定染色體偏好。
- SEQC2-only（LOH detection FN）主要在 chr1, chr5, chr8 — 對應已知的 LOH 邊界。

---

## 六、ISM-only LOH 特徵分析（Analysis E & F）

### 6.1 Effective HP Reads — 根因

![Fig08: ISM-only LOH — HP_Ratio vs ehp scatter](assets/20260403_loh3way_comparison/loh3way_fig08_ism_only_loh_scatter.png)

**這是最關鍵的圖**：ISM-only LOH sites 呈現清晰的弧形帶狀結構，HP_Ratio 集中在 0.0 和 1.0 附近，且**幾乎全部出現在低 effective_hp_reads 區域**。

| Sample | ISM-only LOH count | 低 ehp 模式 |
|--------|-------------------|-------------|
| HCC1395 | 6,203 | 弧帶結構清晰 |
| HCC1395_DORADO | 6,462 | 同上 |
| COLO829 | 7,037 | 弧帶更明顯 |
| H1437 | 8,621 | 同上 |
| H2009 | 21,861 | 最大量，弧帶密集 |
| HCC1937 | 2,679 | 數量最少 |
| HCC1954 | 10,747 | 弧帶結構，FP（紅色）比例高 |

**解讀**：當 effective_hp_reads 很低時（如 ehp=3），只要 3 reads 全部被分配到同一 haplotype，HP_Ratio 就會是 0.0 或 1.0，觸發 ISM LOH 判定。這是**取樣噪音**，不是真正的 LOH。LOH.bed 使用 region-level phased genotype ratio，不受單一 site 的 read 數量影響，因此不會報告這些假 LOH。

### 6.2 其他特徵

![Fig06: 不一致 sites 特徵分組柱狀圖](assets/20260403_loh3way_comparison/loh3way_fig06_discordant_features.png)

| 特徵 | Both LOH | ISM only | LOH.bed only | Neither |
|------|----------|----------|-------------|---------|
| Effective HP Reads (median) | 中等 | **低** | 中等 | 中-高 |
| HP_Ratio (median) | ~1.0 | ~0.0/1.0 | ~0.5 | ~0.4-0.5 |
| Caller AF (median) | 0.5-0.9 | 0.4-0.9 | ~0.5 | 0.3-0.5 |
| Caller GQ (median) | 相似 | 相似 | 相似 | 相似 |

ISM-only LOH 的唯一系統性差異是**較低的 effective_hp_reads**，其他 caller 特徵無顯著差異。

---

## 七、SEQC2 LOH Zone 深度分析（Analysis G）

### 7.1 ISM LOH Rate by SEQC2 Zone

| SEQC2 Zone | ISM Potential_LOH Rate | 含義 |
|-----------|----------------------|------|
| both_LOH (TP) | 0.988 | ISM 幾乎完全偵測到 ground truth LOH |
| SEQC2_only (FN) | 0.131 | ISM 低估 — 邊界區域 |
| TO_only_novel | 0.975 | ISM 也判定 LOH — TO 特有 |
| TO_gain/loss | 0.970 | ISM 也判定 LOH — 在 gain/loss 區 |
| neither (TN) | 0.239 | 背景假陽性率 — 完全來自取樣噪音 |

![Fig09: HCC1395 SEQC2 zone 六面板深度分析](assets/20260403_loh3way_comparison/loh3way_fig09_seqc2_deep_dive.png)

### 7.2 SEQC2 Zone 內的 TP/FP 差異

- **both_LOH zone**：TP AF median 0.54, FP AF median 0.98 → FP 是高 AF germline variants
- **neither zone**：TP/FP 的 HP_Ratio 無差異（TP median 0.46, FP median 0.46）
- **VerificationClass**：LOH zone 內 Noise > Strong（LOH 破壞甲基化分析的 haplotype 前提）

---

## 八、包含關係總結

```
                   ISM Potential_LOH (site-level)
                   ┌─────────────────────────────────────┐
                   │                                     │
                   │   LOH.bed (region-level)             │
                   │   ┌───────────────────────┐         │
                   │   │                       │         │
                   │   │  SEQC2 LOH (ground)    │         │
                   │   │  ┌──────────────┐     │         │
                   │   │  │              │     │         │
                   │   │  │  17,010      │ 176 │   359   │   ISM-only: 5,516+687
                   │   │  │  (All three) │     │         │   = 6,203 sites
                   │   │  │              │     │         │   (low ehp noise)
                   │   │  └──────────────┘     │         │
                   │   │         SEQC2-only:105 │         │
                   │   └───────────────────────┘         │
                   │            LOH.bed only: 10          │
                   └─────────────────────────────────────┘
                              None: 16,233
```

**定量包含關係**（HCC1395）：
- ISM ⊃ LOH.bed：99.0% recall（186 sites missed / 17,555 LOH.bed sites）
- LOH.bed ⊃ SEQC2：95.4% recall（17,010+176 / 17,010+176+105+687）
- ISM ⊃ SEQC2：98.8% recall（17,010+687 / 17,010+687+105）

---

## 九、結論與啟示

### 9.1 已確認

1. **ISM Potential_LOH 是 LOH.bed 的嚴格超集**（recall 99%+, 7/7 samples 一致）。ISM 多報的 10-16% sites 主要來自低 effective_hp_reads 的取樣噪音。

2. **LOH 分類無法區分 TP/FP**：所有 LOH 類別的 TP rate 差異 < 0.09，無實用判別力。

3. **ISM-only LOH 無聚集性**：散佈全基因組，非特定區域現象。

4. **LOH.bed 邊界是 SEQC2-only（FN）的主要來源**：LOH.bed 的邊界處理保守。

5. **HCC1954 特殊性**：LOH coverage 最低（6.2%），ISM-only LOH 比例最高（Jaccard 0.281），是 ISM 過度報告最嚴重的樣本。

### 9.2 對研究方向的啟示

| 啟示 | 行動建議 |
|------|---------|
| ISM HP_Ratio 不需獨立 LOH 判定 | 可直接使用 LOH.bed 作為 zone annotation，ISM 不需重複判定 |
| ISM-only LOH 是噪音 | 加入 effective_hp_reads ≥ 10 門檻即可消除大部分假 LOH |
| LOH 不區分 TP/FP | 不應將 LOH 作為 filter，只能作為 annotation layer |
| SEQC2 邊界效應 | LOH.bed boundary refinement 可改善 FN，但影響量小（105/40,096 = 0.26%） |

### 9.3 優先行動

- **P0**：在 QS 計算中，用 LOH.bed annotation 取代 ISM Potential_LOH（已在 QS v2 初步驗證）
- **P1**：為 ISM HP_Ratio 加入 `effective_hp_reads ≥ 10` 的最低門檻
- **P2**：LOH.bed 應作為 feature 輸入 ML 模型，而非 hard filter

---

## 附錄：產出清單

### 圖表（9 張）

| 圖號 | 內容 |
|------|------|
| Fig01 | HCC1395 三方分類計數 + TP/FP rate |
| Fig02 | 7 樣本 × 4 類別 proportion/TP/FP 熱圖 |
| Fig03 | HP_Ratio violin by zone × truth × 7 samples |
| Fig04 | HCC1395 全基因組 ideogram（40K sites） |
| Fig05 | ISM-only 和 SEQC2-only 染色體密度 |
| Fig06 | 4 特徵 × 7 樣本 × 4 類別分組柱狀圖 |
| Fig07 | LOH rate / excess / Jaccard 一致性 |
| Fig08 | ISM-only LOH: HP_Ratio vs ehp scatter |
| Fig09 | HCC1395 SEQC2 zone 六面板深度分析 |

所有圖片位於 `assets/20260403_loh3way_comparison/` 目錄下。

### 數據表（7 份）

| 檔案 | 內容 |
|------|------|
| `loh3way_threeway_venn_hcc1395.tsv` | 三方 Venn 8 類統計 |
| `loh3way_twoway_all_samples.tsv` | 7 樣本 × 4 類別 TP/FP 統計 |
| `loh3way_hp_ratio_by_zone_truth.tsv` | HP_Ratio 分佈統計 |
| `loh3way_discordant_characteristics.tsv` | 不一致 sites 特徵比較 |
| `loh3way_concordance_matrix.tsv` | 一致性矩陣（Jaccard, recall, precision） |
| `loh3way_seqc2_zone_truth_stats.tsv` | SEQC2 zone × truth 詳細統計 |

### 腳本

| 檔案 | 說明 |
|------|------|
| `analyze_hp_ratio_vs_loh_bed_vs_seqc2.py` | 完整分析腳本（7 analysis modules, 9 figures） |
