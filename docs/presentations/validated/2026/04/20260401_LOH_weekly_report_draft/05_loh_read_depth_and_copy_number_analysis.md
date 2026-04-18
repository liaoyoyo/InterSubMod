<!--
建立時間: 2026-04-02 00:30
目標: LOH Read Depth、Copy Number、Concordance 四象限深度分析 — 驗證 LOH 判定是否受 read depth 與 copy number 影響
處理範圍: 748,391 regions (all_region_rows) + 288,609 matched loci (same_locus_compare)
關聯檔案:
  - 01_full_narrative_report.md
  - scripts/analysis/build_loh_read_depth_analysis.py
  - scripts/analysis/build_loh_concordance_depth_analysis.py
  - scripts/analysis/build_loh_copy_number_analysis.py
-->

# LOH Read Depth 與 Copy Number 驗證分析

## 問題定義

在我們的 LOH 分析中，有三個需要驗證的基礎假設：

### 假設 H1：LOH 區域的 read depth 是否真的少一半？

**理論預期**：如果 LOH 代表一個 allele 真正被刪除（copy-loss LOH），region 的總 read 數應該約為 non-LOH 區域的 50%。如果 LOH 是 copy-neutral（UPD），read depth 應該不變。如果有 copy gain，read depth 可能增加。

**為什麼重要**：如果 LOH 區域的 reads 沒有減少，表示 LOH 主要是 copy-neutral 或 copy-gain，而非真正的 allelic deletion。這影響我們對 LOH enrichment 翻轉現象的機制解釋。

### 假設 H2：TO-only LOH 是否由 read depth 不足造成？

**理論預期**：如果 TO-only LOH（只有 TO 判為 LOH，Paired 認為 non-LOH）是因為 TO 的 effective_hp_reads 較低，導致 HP_Ratio 因隨機偏移而偏向極端，那麼 TO-only LOH 的讀取深度應該顯著低於 both_LOH。

**替代假說**：如果 self-phasing circular dependency 才是主因，read depth 差異應該很小，但 HP ratio 偏移很大。

### 假設 H3：Copy number 異常是否系統性影響 LOH 判定？

**理論預期**：癌症細胞系有大量 copy number 變異（amplification/deletion），如果特定染色體的 depth 偏低（deletion），該染色體的 LOH rate 會升高。反之，如果 copy gain 造成 LOH，需要額外注意。

---

## 分析方法

### 數據來源
- **Master dataset**: `all_region_rows.tsv.gz` — 748,391 regions × 116 features × 7 samples × 2 modes
- **Concordance dataset**: `same_locus_compare.tsv` — 288,609 個 Paired/TO 同位點配對（TP-only, post HP-fix）

### LOH 判定標準
- **core_loh_like**: HP_Ratio < 0.1 or > 0.9 且 effective_hp_reads ≥ 30
- **effective_hp_reads**: HP1FamilyN + HP2FamilyN（排除 HP0 和 HP3 reads）

### Copy Number 分類方法
以每個 sample 的 genome-wide median read depth 為基線，將每個 region 分類：
- **Deep deletion**: NumReads < sample_median × 0.3
- **Copy-loss**: sample_median × 0.3 ~ × 0.7
- **Copy-neutral**: sample_median × 0.7 ~ × 1.3
- **Copy-gain**: sample_median × 1.3 ~ × 2.0
- **High amplification**: > sample_median × 2.0

### 統計方法
- Wilcoxon rank-sum test（非參數比較）
- Cohen's d（效應量）
- Spearman correlation（染色體 depth deviation vs LOH rate）

---

## 結果一：LOH 區域 Read Depth ≈ 0.73x（不是 0.5x）

### 整體數據

| Mode | Truth | LOH median | non-LOH median | Ratio | p-value |
|------|-------|-----------|----------------|-------|---------|
| Paired | TP | 64 | 79 | **0.81** | ~0 |
| Paired | FP | 85 | 66 | **1.29** | 3.8e-7 |
| TO | TP | 69 | 82 | **0.84** | ~0 |
| TO | FP | 64 | 71 | **0.90** | ~0 |

### Per-Sample TP Read Depth Ratio

| Sample | Paired TP | TO TP | 備註 |
|--------|-----------|-------|------|
| HCC1395 | 0.68 | 0.72 | 最低，接近 copy-loss |
| H1437 | 0.69 | 0.74 | |
| HCC1395_DORADO | 0.71 | 0.74 | |
| HCC1937 | 0.72 | 0.75 | |
| COLO829 | 0.74 | 0.77 | |
| H2009 | 0.77 | 0.82 | |
| HCC1954 | **0.97** | **1.02** | LOH depth ≈ non-LOH，強烈暗示 copy gain |

![Figure 1: LOH vs non-LOH Read Depth Distribution](figures/loh_read_depth_distribution.png)

**解讀**：Violin plot 顯示 LOH 區域（右側）的 NumReads 分佈向下偏移，但幅度不到一半。TP（藍色）在 LOH 區域 depth 降低，FP（紅色）在 Paired 的 LOH 區域 depth 反而偏高。

![Figure 2: Per-Sample LOH Read Depth Ratio](figures/loh_read_depth_ratio_per_sample.png)

**解讀**：紅色虛線標記 0.5（純 copy-loss 理論值）和 1.0（copy-neutral）。所有 TP sample 的 ratio 落在 0.68-1.02 之間，遠高於 0.5。HCC1954 最異常（ratio ≈ 1.0），暗示大量 copy gain 伴隨 LOH。

### HP1 vs HP2 Balance 驗證

![Figure 3: HP1 vs HP2 Read Count Balance](figures/loh_hp_balance_scatter.png)

**解讀**：
- **紅點（LOH）**：聚集在 X 軸和 Y 軸附近 → 一邊 HP 接近 0，確認 LOH 判定正確
- **藍點（non-LOH）**：沿對角線分布 → 兩邊 HP 平衡
- Paired 和 TO 的 pattern 一致，HP imbalance 中位數：LOH = 1.000（完全不平衡），non-LOH = 0.25-0.27

### Per-Chromosome Read Depth Heatmap

![Figure 4: Per-Chromosome LOH Read Depth Ratio](figures/loh_chr_depth_heatmap.png)

**解讀**：熱圖每格 = LOH median depth / non-LOH median depth（per sample × chromosome）。藍色 < 0.5 暗示 copy-loss LOH，紅色 > 1.0 暗示 copy-gain LOH。黑框 = 極端偏差（> 0.5 或 < -0.3）。Paired 和 TO 的染色體 pattern 高度一致，確認 depth 反映真實 copy number 而非 calling artifact。

### H1 結論

> **LOH 區域 read depth 典型為 non-LOH 的 0.73x，不是 0.5x。** 這意味著大部分 LOH 並非純粹的 allele deletion (1+0)，而是伴隨部分 copy gain (2+0 或更高)。HCC1954 是極端案例（ratio ≈ 1.0），與其已知的大量 copy number aberration 一致。Paired FP 在 LOH 區域 depth 反向偏高（ratio = 1.29），暗示高 depth 可能帶來 artifacts。

---

## 結果二：TO-only LOH 由 Self-Phasing Bias 造成（非 Read Depth）

### 四象限 Read Depth 比較

| 象限 | n | TO eff_hp median | Paired eff_hp median |
|------|-------|-----------------|---------------------|
| both_LOH | 88,109 | **60** | **60** |
| TO_only_LOH | 39,978 | **68** | **61** |
| Paired_only_LOH | 1,874 | 65 | 43 |
| neither_LOH | 158,648 | **79** | **79** |

**關鍵發現**：TO-only LOH 的 TO effective_hp_reads 中位數 = **68**，**高於** both_LOH 的 60。Cohen's d = +0.29（小效應，方向相反）。

![Figure 1: Effective HP Reads by Concordance Quadrant](figures/concordance_quadrant_depth.png)

**解讀**：Box plot 清楚顯示 TO-only LOH（橘色）的 TO reads 與 both_LOH（綠色）相當甚至更高。neither_LOH（紅色）的 reads 最高，因為高 depth 區域 HP ratio 更穩定。

### HP Ratio 散點圖 — Self-Phasing 的決定性證據

![Figure 2: HP Ratio Paired vs TO, Coloured by Concordance Quadrant](figures/concordance_hp_ratio_scatter.png)

**解讀**：
- **藍色（neither_LOH）**：聚集在中央（HP ratio 0.2-0.8 for both modes）
- **綠色（both_LOH）**：聚集在四個角落（<0.1 或 >0.9 for both modes）
- **橘色（TO_only_LOH）**：分佈在**上下兩條水平帶** — X 軸（Paired）在 0.2-0.8（平衡），Y 軸（TO）在 <0.1 或 >0.9（極端）

這是 self-phasing circular dependency 的**直接可視化**：在同一基因組位點，Paired 模式看到平衡的 haplotype 分配，但 TO 模式因為 LongPhase-TO 將 somatic variant 自身納入 phasing anchor，將 reads 系統性推向一側。

### 定量統計

| 比較 | 指標 | Cohen's d | p-value | 解讀 |
|------|------|-----------|---------|------|
| TO-only vs both_LOH | TO hp_ratio | **-1.20** | ~0 | **巨大效應** |
| TO-only vs both_LOH | Paired hp_ratio | **0.005** | 0.005 | **無效應** |
| TO-only vs both_LOH | TO eff_hp_reads | **+0.29** | 2e-254 | 小效應，方向反向 |
| TO-only vs both_LOH | Paired eff_hp_reads | +0.04 | 5e-12 | 幾乎無差異 |

### 86.5% 的 TO-only LOH 位點在 Paired 模式下 HP ratio 完全平衡

39,978 個 TO-only LOH 位點中：
- **86.5%** 的 Paired HP_Ratio 在 0.2-0.8（平衡）
- **0%** 的 Paired HP_Ratio 在 <0.1 或 >0.9（極端）
- **100%** 的 TO HP_Ratio 在 <0.1 或 >0.9（by definition）

### Read Depth vs HP Ratio 的 2D 關係

![Figure 3: Read Depth vs HP Ratio (Hexbin)](figures/concordance_depth_vs_hpratio.png)

**解讀**：
- **TO-only LOH（左上）**：即使在高 depth (>100) 區域，TO HP ratio 仍被推到 <0.1 或 >0.9，排除「低 depth 隨機偏移」假說
- **both_LOH（中上）**：HP ratio 極端值集中在 depth 30-100 區域
- **neither（右上）**：HP ratio 集中在 0.3-0.7，不受 depth 影響
- 底排（Paired）的 pattern 與 TO 不同：Paired 的 TO-only 位點（左下）HP ratio 分布在 0.3-0.7，進一步確認 Paired 看不到 LOH

### Read Depth 差異分佈（TO - Paired）

![Figure 5: Read Depth Difference (TO - Paired) by Quadrant](figures/concordance_depth_diff.png)

**解讀**：四象限的 TO - Paired read depth 差異中位數都接近 0（-0.3 ~ +0.5），說明同位點的 Paired 和 TO effective_hp_reads 幾乎相同。TO-only LOH 的差異略偏正（TO reads 稍多），進一步否定「TO reads 少所以 LOH」假說。

### Per-Sample 分佈

![Figure 4: Per-Sample Median Effective HP Reads by Quadrant](figures/concordance_per_sample_depth.png)

**解讀**：所有 7 個 sample 的 TO-only LOH 佔比在 9.7%-17.4%（1.8x 範圍），無特定 sample 主導。這是 **TO self-phasing 的系統性特性**，不是 sample-specific artifact。

### H2 結論

> **Self-phasing circular dependency 是 TO-only LOH 的主因，read depth 不是。** 在 39,978 個 TO-only LOH 位點中，TO effective_hp_reads 中位數（68）反而高於 both_LOH（60）。決定性證據：同一位點 Paired HP_Ratio 中位數 = 0.509（完全平衡），TO HP_Ratio 中位數 = 0.026（極端 LOH），Cohen's d = -1.20。86.5% 的 TO-only LOH 位點在 Paired 下完全平衡。

---

## 結果三：Copy Number 組成分析

### LOH 的 Copy Number 分類

| Mode | Copy-neutral | Copy-loss + Deep del | Copy-gain + High amp |
|------|-------------|---------------------|---------------------|
| Paired | **61.8%** | 29.8% | 8.4% |
| TO | **60.1%** | 24.2% | **15.8%** |

![Figure 3: LOH Proportion by Copy Number Category](figures/cn_category_loh_bar.png)

**解讀**：
- **Copy-neutral**（depth 在 sample median 的 0.7-1.3x）：LOH rate 最高（Paired-TP ~40%），因為包含大量 copy-neutral LOH (UPD)
- **Copy-loss**（depth < 0.7x）：LOH rate 中等（~25-35%），代表真正的 allele deletion
- **Copy-gain**（depth > 1.3x）：LOH rate 偏低（Paired-TP ~15%），但 TO-FP 的 copy-gain LOH 比例明顯偏高
- TO copy-gain LOH（15.8%）是 Paired（8.4%）的 **2 倍**

### Per-Chromosome Depth Deviation 熱圖

![Figure 1: Per-Sample Per-Chromosome Median Read Depth](figures/cn_chr_depth_heatmap.png)

**解讀**：Paired 和 TO 的染色體 depth pattern 高度一致（同一 BAM 的不同 calling mode），確認 depth 反映真實 copy number。黑框標記 deviation > ±50%，這些染色體可能有大段 amplification 或 deletion。

### LOH Rate vs Chromosome Depth Deviation

![Figure 2: LOH Rate vs Chromosome Depth Deviation](figures/cn_loh_rate_vs_depth.png)

**解讀**：
- Overall Spearman rho = **-0.33**（Paired）/ **-0.35**（TO），p < 1e-4
- 中等強度負相關：depth 偏低的染色體 LOH rate 偏高
- 但相關性不強（|rho| < 0.4），多數 LOH 仍發生在 copy-neutral 區域

### Per-Sample LOH 的 Copy Number 組成

![Figure 5: Copy Number Composition of LOH Regions Per Sample](figures/cn_loh_composition_per_sample.png)

**解讀**：
- 所有 sample 的 LOH 都以 **copy-neutral 為主**（灰色，50-70%）
- Copy-loss（紅色）佔 20-35%
- Copy-gain（藍色）佔 5-15%，TO mode 普遍高於 Paired
- HCC1954 的 copy-gain 比例最高，與已知的大量基因組重排一致

### H3 結論

> **LOH 主要是 copy-neutral（~60%），不是 copy-loss。** Copy number 對 LOH 判定有中等程度影響（Spearman rho ≈ -0.35），但不是主要決定因素。TO mode 的 copy-gain LOH 比例是 Paired 的 2 倍（15.8% vs 8.4%），且集中在 FP 中，暗示 TO 在 copy-gain 區域可能過度呼叫 LOH。

---

## 綜合結論

### 1. LOH ≠ Allelic Deletion

LOH 區域 read depth 約為 non-LOH 的 **0.73x**（而非理論的 0.5x）。約 60% 的 LOH 是 **copy-neutral**（可能是 UPD 或 phasing artifact），24-30% 是 **copy-loss**（真正的 allele deletion），8-16% 是 **copy-gain**（一個 allele 丟失但另一個擴增）。

### 2. TO-only LOH = Self-Phasing Artifact（確認）

在 39,978 個 TO-only LOH 位點：
- TO effective_hp_reads **更高**（median 68 vs 60），排除 read depth 假說
- 同位點 Paired HP_Ratio = 0.509（平衡），TO = 0.026（極端），Cohen's d = **-1.20**
- 86.5% 的 TO-only 位點在 Paired 模式下完全平衡
- 效應在所有 7 個 sample 中一致（9.7%-17.4%）
- **結論：LongPhase-TO 的 self-phasing circular dependency 是 TO-only LOH 的主因**

### 3. Copy Number 是 Confounding Factor（中等影響）

- Chromosome depth deviation 與 LOH rate 有 moderate 負相關（rho ≈ -0.35）
- TO mode 在 copy-gain 區域過度呼叫 LOH（15.8% vs Paired 8.4%）
- 建議未來考慮在 QS 中納入 depth deviation 校正

### 4. 對週報結論的影響

本分析**強化**了週報的核心結論：
- ✅ LOH enrichment 翻轉（Paired FP-enriched / TO TP-enriched）**不是 read depth artifact** — 已確認 TO-only LOH 的 depth 與 both_LOH 相當
- ✅ Self-phasing circular dependency **被定量確認** — Cohen's d = -1.20 是巨大效應
- ⚠️ 新發現：LOH 主要是 copy-neutral，暗示大量 LOH 可能是 **phasing artifact 而非真正的 allelic imbalance** — 這對 paired mode 的 LOH 判定也提出疑問
- ⚠️ 新發現：TO 在 copy-gain 區域的 LOH 判定可能有額外的系統性偏差

---

## 輸出清單

### 分析腳本
- `scripts/analysis/build_loh_read_depth_analysis.py`
- `scripts/analysis/build_loh_concordance_depth_analysis.py`
- `scripts/analysis/build_loh_copy_number_analysis.py`

### 圖表（figures/ 目錄下）

| 圖表 | 檔名 | 回答的問題 |
|------|------|-----------|
| Read Depth Distribution | `loh_read_depth_distribution.png` | LOH vs non-LOH reads 差多少？ |
| Per-Sample Depth Ratio | `loh_read_depth_ratio_per_sample.png` | 哪些 sample LOH depth 異常？ |
| HP Balance Scatter | `loh_hp_balance_scatter.png` | LOH 判定是否正確？ |
| Chr Depth Heatmap | `loh_chr_depth_heatmap.png` | 哪些 chr 有 CN 異常？ |
| Quadrant Depth | `concordance_quadrant_depth.png` | TO-only LOH depth 低嗎？ |
| HP Ratio Scatter | `concordance_hp_ratio_scatter.png` | Self-phasing 的直接證據 |
| Depth vs HP Ratio | `concordance_depth_vs_hpratio.png` | 低 depth 造成極端 HP ratio？ |
| Per-Sample Quadrant Depth | `concordance_per_sample_depth.png` | 哪些 sample TO-only 異常？ |
| Depth Difference | `concordance_depth_diff.png` | TO vs Paired reads 差異 |
| CN Chr Heatmap | `cn_chr_depth_heatmap.png` | Per-chr copy number profile |
| CN vs LOH Rate | `cn_loh_rate_vs_depth.png` | Copy number 影響 LOH rate？ |
| CN Category LOH Bar | `cn_category_loh_bar.png` | LOH 在哪個 CN 類別最多？ |
| CN LOH Composition | `cn_loh_composition_per_sample.png` | LOH 主要是 loss/neutral/gain？ |

### 統計表
- `figures/loh_read_depth_stats.tsv`
- `figures/concordance_depth_stats.tsv`
- `figures/cn_analysis_stats.tsv`
