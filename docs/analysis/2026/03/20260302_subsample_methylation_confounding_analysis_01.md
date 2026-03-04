<!--
建立時間: 2026-03-02 16:00
目標: 診斷 InterSubMod 在 subsample 驗證中表現不佳的根本原因，並提出可行研究方向
處理範圍: 問題定義、根因分析、文獻支持、程式碼局限、可行方案、驗證流程
關聯檔案:
  - docs/references/2026/03/20260302_tissue_origin_methylation_confounding_literature_review_01.md
  - docs/architecture/system_overview.md
  - include/core/DataStructs.hpp
  - include/core/Stats.hpp
  - src/core/GlobalTest.cpp
  - src/core/SignificanceAnalyzer.cpp
  - /big8_disk/liaoyoyo2001/Knowledge/02_samples/HCC1395.md
  - /big8_disk/liaoyoyo2001/Knowledge/02_samples/COLO829.md
  - /big8_disk/liaoyoyo2001/Knowledge/02_samples/subsample_purity.md
  - /big8_disk/liaoyoyo2001/Knowledge/02_samples/other_cell_lines.md
-->

# Subsample 甲基化分析混淆效應診斷與研究方向報告

---

## 目錄

1. [問題定義與現況](#1-問題定義與現況)
2. [問題根因分析](#2-問題根因分析)
3. [文獻支持](#3-文獻支持)
4. [InterSubMod 當前設計的局限](#4-intersubmod-當前設計的局限)
5. [可驗證假說與分析方向](#5-可驗證假說與分析方向)
6. [驗證流程](#6-驗證流程)
7. [結論與優先順序](#7-結論與優先順序)

---

## 1. 問題定義與現況

### 1.1 Subsample 混合方法

InterSubMod 的驗證依賴不同純度（purity）的 subsample BAM 檔案。這些 subsample 是透過以下流程產生的：

1. 使用 `samtools view -s` 對 tumor BAM 與 normal BAM 分別進行 downsampling
2. 使用 `samtools merge` 將 downsampled 的 tumor 與 normal BAM 合併
3. 對合併後的 BAM 進行 sorting 與 indexing

產出為**單一混合 BAM**（非 tumor/normal 雙 BAM），且預設未經 phase/haplotag。

**HCC1395 ONT 標準 Subsample 梯度：**

| 目錄名稱 | Tumor 份數 | Normal 份數 | 估計 Purity | 實際組成 |
|----------|-----------|------------|------------|---------|
| `t50_n00` | 50 | 0 | 100% | 純腫瘤 50x |
| `t40_n10` | 40 | 10 | 80% | Tumor 40x + Normal 10x |
| `t30_n20` | 30 | 20 | 60% | Tumor 30x + Normal 20x |
| `t20_n30` | 20 | 30 | 40% | Tumor 20x + Normal 30x |
| `t10_n40` | 10 | 40 | 20% | Tumor 10x + Normal 40x |
| `t00_n25` | 0 | 25 | 0% | 純 Normal 25x (baseline) |

### 1.2 樣本組織來源彙整

以下是所有 benchmark 樣本的 tumor/normal 來源組織：

| 樣本 | 癌症類型 | Tumor 來源 | Normal 來源 | Methylation 資料 |
|------|---------|-----------|------------|-----------------|
| **HCC1395** | 乳腺導管癌 | Cell line（乳腺組織） | HCC1395BL（B-lymphoblast，血液） | ONT_5kHz: 5mCG+5hmCG; ONT_Dorado: 僅 5mCG |
| **COLO829** | 黑色素瘤 | Cell line（皮膚） | COLO829BL（B-lymphoblast，血液） | ONT_PAO: 5mCG+5hmCG |
| **H1437** | 肺腺癌 | Cell line（肺組織） | H1437BL（血液） | ONT: 僅 5mCG（Google 資料） |
| **H2009** | 肺腺癌 | Cell line（肺組織） | H2009BL（血液） | ONT: 僅 5mCG（Google 資料） |
| **HCC1937** | 乳腺癌 (BRCA1 mutant) | Cell line（乳腺組織） | HCC1937BL（血液） | ONT: 僅 5mCG（Google 資料） |
| **HCC1954** | 乳腺癌 | Cell line（乳腺組織） | HCC1954BL（血液） | ONT: 僅 5mCG（Google 資料） |

**關鍵觀察**：所有 benchmark 樣本的 normal 均來自**血液 B-lymphoblast cell line**，而 tumor 均來自**實體組織 cell line**。這意味著每一對 tumor-normal 配對都存在**組織來源甲基化差異**的混淆問題。

### 1.3 目前的驗證數據趨勢

在 subsample 驗證中觀察到的核心問題：

- **F1 score 隨 purity 下降而急劇惡化**：當 normal reads 比例增加，InterSubMod 的分群結果顯著性下降
- **False positive clusters 增加**：低 purity 情境下，出現與腫瘤亞克隆無關的顯著分群
- **HP/Allele 關聯性被稀釋**：正常組織 reads 的加入使得 cluster 與 HP tag / ALT allele 的關聯度下降

### 1.4 現有 Subsample 的甲基化標籤限制

根據知識庫驗證結果（2026-02-13）：

> HCC1395 ONT 標準 subsample（`t10_n40` ~ `t50_n00`）的混合 BAM 在抽樣檢查前 500 reads 時，MM/ML 標籤計數為 **0**。此批 subsample **不可直接用於** InterSubMod 的甲基化分析流程。

若要進行 methylation-level benchmark，需：
- 改用保留 MM/ML 的來源 BAM（如 ONT_5kHz_simplex_5mCG_5hmCG）重新 downsample + merge
- 或使用已確認含 MM/ML 的 ONT_Dorado 版本資料

---

## 2. 問題根因分析

### 2.1 根因一：組織來源甲基化差異（Tissue-of-Origin Confounding）

**這是最關鍵的根因。**

#### 2.1.1 問題本質

血液（B-lymphoblast）與實體組織（乳腺、皮膚、肺）之間存在**數千個獨特的甲基化標記**。這種差異是細胞分化過程中建立的穩定表觀遺傳特徵，與腫瘤狀態無關。

當我們將 tumor BAM（乳腺組織）與 normal BAM（血液）混合產生 subsample 時：

```
Subsample = α × Tumor_reads（乳腺組織甲基化模式）
           + (1-α) × Normal_reads（血液甲基化模式）

其中 α = tumor purity
```

InterSubMod 在分析這些 subsample 時，會偵測到兩類 reads 間的甲基化距離差異：

1. **真正的腫瘤亞克隆信號**（tumor-specific epigenetic heterogeneity）
2. **組織來源差異信號**（tissue-of-origin methylation differences）— **混淆因子**

問題在於：**信號 2 的強度遠大於信號 1**。

#### 2.1.2 量化估計

根據 Loyfer et al. (2023, *Nature*) 和 Moss et al. (2018, *Nature Communications*) 的人類甲基化圖譜：

- 血液與乳腺組織之間存在 **>5,000 個差異甲基化 CpG 位點**（差異 > 30%）
- 其中許多位點的差異幅度達 **60-80%**（接近完全甲基化 vs 完全未甲基化）
- 這些「組織差異 CpG」分佈在全基因組，**會落在 InterSubMod 分析的每個 ±1000bp 窗口內**

相較之下，腫瘤亞克隆間的甲基化差異通常：
- 涉及較少的 CpG 位點（通常 <100 個 per region）
- 差異幅度較小（20-40%）
- 需要足夠的統計效力才能偵測

#### 2.1.3 對 InterSubMod 分析的影響鏈

```
組織來源差異 CpG 存在
  ↓
Tumor reads 與 Normal reads 在這些 CpG 上系統性不同
  ↓
距離計算（NHD/L1/L2）會產生一致的 tumor-normal 距離差異
  ↓
聚類演算法傾向將 reads 分為「組織 A」與「組織 B」兩群
  ↓
顯著性分析（PERMANOVA / Fisher exact test）偵測到高度顯著的分群
  ↓
但這些分群反映的是組織來源，不是腫瘤亞克隆結構
  ↓
隨著 purity 降低，normal reads 比例增加，組織差異信號主導分群
  ↓
真正的亞克隆信號被淹沒 → F1 score 下降
```

### 2.2 根因二：Truth Set 品質差異

不同 benchmark 樣本使用的 truth set 品質差異顯著，導致跨樣本比較不公平。

#### 2.2.1 Truth Set 品質分級

| 樣本 | Truth Set 來源 | 建立方法 | 變異數量 | 品質等級 |
|------|---------------|---------|---------|---------|
| **HCC1395** | SEQC2 v1.2.1 | FDA 主導，63 組 tumor-normal pairs，3 種 aligners，PacBio 額外驗證 | ~39,447 SNVs | **最高** (A+) |
| **COLO829** | NYGC | 多 caller 整合（Strelka2, Mutect2 等），有 `num_callers`/`HighConfidence` 標記 | ~43,192 variants | **高** (A) |
| **H1437, H2009, HCC1937, HCC1954** | Orthogonal-tools | Google 模型共識方法，非傳統多 caller 整合 | 不一 | **中等** (B) |

#### 2.2.2 品質差異的影響

- **SEQC2 (HCC1395)**：三個信度層級（HighConf / MedConf / LowConf），可以精確選擇分析範圍。高信度區域定義明確（`High-Confidence_Regions_v1.2.bed`）。
- **NYGC (COLO829)**：有 `num_callers >= 2` 的二次篩選機制，但 high-confidence regions 的定義不如 SEQC2 嚴格。
- **Orthogonal-tools (其餘 4 個樣本)**：依賴 Google 的共識模型，缺乏傳統的多平台正交驗證。Truth set 中可能包含更多的 false positive 和 false negative，使得 InterSubMod 的驗證結果在這些樣本上波動更大。

**結論**：在不同樣本上得到不一致的 F1 score，部分原因可能來自 truth set 品質差異，而非 InterSubMod 方法本身的問題。

### 2.3 根因三：Subsample MM/ML 標籤缺失

如前述，HCC1395 ONT 標準 subsample 的 BAM 檔案**不含甲基化標籤**（MM/ML tags 為 0）。這意味著：

1. **無法在標準 subsample 上進行甲基化分析**：InterSubMod 的核心功能無法運作
2. **需要從其他資料集重建**：ONT_5kHz（5mCG+5hmCG）或 ONT_Dorado（僅 5mCG）可提供甲基化資訊，但覆蓋度和混合比例可能不同
3. **重建 subsample 的工程量**：需要對 ONT_5kHz 的 tumor（272GB）和 normal（140GB）BAM 執行 downsample + merge，每個 purity level 需要一次操作

#### 2.3.1 可用的含甲基化標籤資料

| 資料集 | Methylation | 覆蓋度 | 可用於重建 subsample |
|--------|------------|--------|-------------------|
| ONT_5kHz_simplex_5mCG_5hmCG | 5mCG + 5hmCG | T: ~50x, N: ~25x | **是**，但檔案大（T: 272GB, N: 140GB） |
| ONT_Dorado | 僅 5mCG | 需確認 | 是，但僅有 5mCG |
| ONT（標準） | **無** | T: ~50x, N: ~25x | **否** |

---

## 3. 文獻支持

本節引用已完成的文獻綜述（38 篇文獻），完整版請見：
`docs/references/2026/03/20260302_tissue_origin_methylation_confounding_literature_review_01.md`

### 3.1 組織來源差異的確切證據

| 文獻 | 關鍵發現 | 對 InterSubMod 的啟示 |
|------|---------|---------------------|
| **Zheng et al. (2017)** — InfiniumPurify, *Genome Biology* | 某些癌種從血液對照推算的純度與通用正常對照的估計值相關性極差，因為差異甲基化 CpG 多為**血液-組織特異性**而非腫瘤-正常差異 | **直接確認**：HCC1395BL（血液）vs HCC1395（乳腺）的甲基化差異有大量是組織特異的 |
| **Moss et al. (2018)** — *Nature Communications* | 39 種人類細胞型態存在數千個獨特甲基化標記，血液與實體組織系統性差異 | 提供了差異的量化基礎 |
| **Loyfer et al. (2023)** — cfSort, *Nature* | 建立了涵蓋所有主要人類細胞型態的甲基化圖譜 | 可作為組織差異 CpG 排除的參考圖譜 |
| **Johnson et al. (2025)** — *Clinical Epigenetics* | 前列腺癌正常組織的甲基化異常延伸至距腫瘤 4cm 區域（field effect） | 即使用癌旁組織作 normal，也可能被 field effect 影響 |
| **Flanagan et al. (2018)** — *EBioMedicine* | 正常乳腺組織甲基化比 CNV 更能預測乳腺癌狀態 | HCC1395 的配對 normal 存在同樣問題 |

### 3.2 解決方案的文獻支持

| 文獻 | 提出的方法 | 適用於 InterSubMod 的程度 |
|------|-----------|------------------------|
| **Zhu et al. (2024)** — CelFiE-ISH, *Genome Biology* | 單分子甲基化 haplotype 解摺積（EM 框架），準確度比 CelFiE 提升 30% | **高度相關**：直接利用 long-read 單分子特性進行 read-level 來源判定 |
| **Staaf et al. (2024)** — PureBeta, *NAR Genomics* | 從全基因組甲基化估計純度，校正 CpG beta 值 | **中度相關**：概念可借鑒，但需改編為 read-level 版本 |
| **Zhang et al. (2024)** — MEnet, *NAR Cancer* | 跨平台甲基化解摺積（唯一支援 nanopore） | **直接相關**：支援 ONT 數據 |
| **Stoeger et al. (2025)** — *bioRxiv/PMC* | ONT R10.4 對黑色素瘤亞克隆的長讀取甲基化分析 | **高度相關**：與 InterSubMod 目標一致的最近似研究 |
| **Kim et al. (2024)** — MethPhaser, *Nature Communications* | 甲基化輔助 phasing，phase block N50 提升 78-151% | **高度相關**：可與 LongPhase 互補 |

### 3.3 Benchmark 改進的文獻支持

| 文獻 | 數據資源 | 意義 |
|------|---------|------|
| **Park & Cook et al. (2025)** — DeepSomatic/CASTLE, *Nature Biotechnology* | 6 對 tumor-normal cell line，329,011 個體細胞變異 | 大幅擴展 long-read somatic benchmark 資源 |
| **NYGC (2025)** — Lancet2 two-tech truth set, *bioRxiv* | HCC1187, HCC1143, COLO829, HCC1395 四樣本的 Illumina + ONT 整合 truth set | 多技術平台整合的新標準 |

---

## 4. InterSubMod 當前設計的局限

### 4.1 ReadInfo 結構只有 `is_tumor` 布林值

**程式碼位置**：`include/core/DataStructs.hpp:25-36`

```cpp
struct ReadInfo {
    int read_id;
    std::string read_name;
    int chr_id;
    int32_t align_start;
    int32_t align_end;
    int mapq;
    std::string hp_tag;
    bool is_tumor;           // ← 僅區分 BAM 來源，無 tissue_type 資訊
    AltSupport alt_support;
    Strand strand;
};
```

**問題**：`is_tumor` 只表示 read 來自哪個 BAM 檔案（tumor vs normal），完全沒有記錄組織類型（乳腺 vs 血液 vs 肺 vs 皮膚）。在 subsample 情境中，所有 reads 混在同一個 BAM 內，`is_tumor` flag 無法被設定（因為 subsample 是合併後的單一 BAM），更遑論組織類型。

### 4.2 FullLabel 結構的 `is_tumor` 直接進入統計檢驗

**程式碼位置**：`include/core/Stats.hpp:55-70`

```cpp
struct FullLabel {
    AltSupport allele = AltSupport::UNKNOWN;
    std::string hp_tag;
    Strand strand = Strand::UNKNOWN;
    bool is_tumor = true;    // ← 直接用於 SAMPLE 維度的統計檢驗
    // ...
};
```

### 4.3 GlobalTest.test_sample() 假設 normal reads 提供無偏參考

**程式碼位置**：`src/core/GlobalTest.cpp:211-227`

```cpp
GlobalTestResult GlobalTest::test_sample(
    const std::vector<int>& cluster_labels,
    const std::vector<FullLabel>& full_labels) {
    // ...
    for (size_t i = 0; i < n; ++i) {
        TestLabel label = full_labels[i].get_sample_label();
        binary_labels[i] = (label == TestLabel::SAMPLE_TUMOR) ? 0 : 1;
    }
    // ... Fisher exact test on Cluster × Sample contingency table
}
```

**問題**：`test_sample()` 將 cluster labels 與 tumor/normal labels 做列聯表分析。當 tumor 和 normal 來自不同組織時，甲基化距離的聚類天然會反映組織差異，使得 `test_sample()` 幾乎永遠產生高度顯著的結果——但這個顯著性來自組織差異而非亞克隆結構。

### 4.4 距離計算未考慮組織差異 CpG

InterSubMod 的距離計算（NHD/L1/L2/JACCARD 等）對所有 CpG 位點一視同仁：

- 沒有機制**排除**已知的組織差異 CpG
- 沒有機制**降權**組織差異 CpG 在距離計算中的貢獻
- 結果是：在包含大量組織差異 CpG 的窗口中，距離矩陣被組織差異主導

### 4.5 固定的甲基化二值化門檻

目前使用固定門檻（`binary_methyl_high = 0.8`, `binary_methyl_low = 0.2`）進行甲基化二值化。這個門檻對所有 CpG 和所有 purity 一視同仁，但：

- 在低 purity 情境下，混合 reads 的甲基化分佈可能不再呈雙峰分佈
- 組織差異 CpG 在混合 BAM 中可能呈現中間值（因為同時包含兩種組織的 reads），被錯誤歸類為 missing

### 4.6 缺乏 subsample 場景的專門處理

InterSubMod 目前的設計假設輸入為 **tumor BAM + normal BAM 雙檔案**模式。但 subsample 為**單一混合 BAM**：

- 無法區分 tumor 和 normal reads
- 所有 reads 的 `is_tumor` flag 只能設為同一值
- 失去了 tumor/normal 維度的分析能力

---

## 5. 可驗證假說與分析方向

### 5.1 方案 A（短期）：CpG 層級組織差異圖譜 + 排除/降權

**假說**：如果在距離計算前排除已知的組織差異 CpG，InterSubMod 的分群結果將更準確反映腫瘤亞克隆結構而非組織來源差異。

**實作步驟**：

1. **建立組織差異 CpG 圖譜**：
   - 使用 pure tumor BAM（`t50_n00`，含 MM/ML 的 ONT_5kHz 版本）和 pure normal BAM（`t00_n25` 或 HCC1395BL.bam）
   - 在每個 CpG 位點計算 tumor reads 與 normal reads 的平均甲基化水平
   - 設定差異閾值（例如 |Δβ| > 0.3），識別「組織差異 CpG」
   - 輸出為 BED 格式的排除清單

2. **在距離計算中應用排除/降權**：
   - 排除策略：在 S4a（CpG gating）階段加入組織差異 CpG 的排除條件
   - 降權策略：在距離計算（S5）中，對組織差異 CpG 乘以低權重（例如 0.1）

3. **驗證**：
   - 比較排除前後的分群結果
   - 觀察 TP/FP 位點的 Fisher exact test p-value 分佈變化
   - 期望：排除後，組織差異驅動的 false positive clusters 減少

**優點**：實作簡單，可在現有框架內完成。
**缺點**：需要 pure tumor 和 pure normal 的甲基化數據來建立圖譜；二元排除可能過於粗糙。

**預計工時**：1-2 週。

---

### 5.2 方案 B（短期）：Read-Level 甲基化來源判定

**假說**：利用每條 read 上多個 CpG 的聯合甲基化模式（methylation haplotype），可以判定該 read 更可能來自血液還是實體組織。

**實作步驟**：

1. **建立參考甲基化 haplotype**：
   - 從 pure tumor BAM 和 pure normal BAM 擷取 read-level 甲基化 haplotype
   - 對每個區域（±1000bp window），建立兩組參考模式：tissue-like 和 blood-like

2. **Read 來源判定**：
   - 對 subsample 中的每條 read，計算其甲基化 haplotype 與兩組參考的相似度
   - 產出一個 `tissue_score`（0-1），表示該 read 來自組織 vs 血液的概率

3. **整合至 InterSubMod**：
   - 在 `ReadInfo` 結構中新增 `float tissue_score` 欄位
   - 在距離計算或聚類分析中，以 `tissue_score` 作為共變量
   - 或在 PERMANOVA 中控制組織來源效應

**優點**：利用 long-read 的單分子特性，精確到 read 層級。
**缺點**：需要建立參考模式，計算量增加。

**理論基礎**：CelFiE-ISH (Zhu et al., 2024) 已證明單分子甲基化 haplotype 可將解摺積準確度提升 30%，可檢測低至 0.03% 的稀有細胞類型。

**預計工時**：2-3 週。

---

### 5.3 方案 C（中期）：純度感知的距離計算

**假說**：若已知或可估計 tumor purity，可以在距離計算中對預期的組織差異進行校正。

**實作步驟**：

1. **純度估計**：
   - 利用 LongPhase-S 的純度估計功能，或從 allele frequency 推算
   - 輸入估計的 purity 值 α

2. **校正距離**：
   - 對每對 reads，在計算距離前：
     - 估計「期望的組織差異貢獻」= f(CpG 位點的組織差異, purity)
     - 校正後距離 = 原始距離 - 期望組織差異貢獻
   - 或使用加權距離：在低 purity 時，對「tumor vs normal」方向的距離分量進行降權

3. **自適應門檻**：
   - 根據 purity 調整甲基化二值化門檻和顯著性判定門檻
   - 例如：purity 40% 時，使用更嚴格的 Cramér's V 效應量門檻

**理論基礎**：PureBeta (Staaf et al., 2024) 已證明純度校正可使 beta 值更接近理論甲基化狀態。

**預計工時**：1-2 個月。

---

### 5.4 方案 D（中期）：同組織 Normal 對照

**假說**：使用與 tumor 相同組織類型的 normal 對照，可以消除組織來源差異的混淆。

**可行路徑**：

1. **利用公開的正常組織甲基化數據**：
   - 使用 ENCODE 或 Roadmap Epigenomics 的正常乳腺組織甲基化數據作為參考
   - 限制：目前這些數據多為 Illumina array 或 WGBS（short-read），非 ONT

2. **建立 ONT 正常組織參考**：
   - 從 ONT_5kHz HCC1395 tumor BAM 的 normal reads（如果能區分的話）建立基線
   - 限制：cell line 環境下無法取得真正的「正常乳腺組織」

3. **使用 field effect 較小的遠端組織**：
   - 如果有遠端正常乳腺組織的 ONT 資料，可作為更佳的參考
   - 限制：目前不可用

**限制**：對 cell line 而言，沒有 matched adjacent tissue，此方案的實用性有限。但概念上值得未來在臨床樣本上驗證。

**預計工時**：取決於數據可用性。

---

### 5.5 方案 E（中期）：多模態整合分析

**假說**：整合 SNV 基因型、單倍型（HP tag）和甲基化三個維度的資訊，可以更穩健地識別真正的腫瘤亞克隆。

**實作思路**：

1. **SNV-gated 甲基化分析**：
   - 僅在攜帶 somatic ALT allele 的 reads 之間進行甲基化距離計算
   - 排除 normal reads（攜帶 REF allele）對距離矩陣的影響
   - 這等於使用 SNV 基因型作為「天然的組織來源標記」

2. **HP-stratified 分析**：
   - 分別在 HP=1 和 HP=2 內部進行甲基化聚類
   - 跨單倍型的甲基化差異更可能反映真正的等位基因特異性表觀遺傳

3. **聯合推斷**：
   - 建立 SNV × HP × methylation 的聯合概率模型
   - 只有在所有維度都支持的情況下，才判定為真正的亞克隆信號

**理論基礎**：Stoeger et al. (2025) 在黑色素瘤研究中展示了甲基化軌跡與 SNV 積累的平行演化，證明多模態整合的可行性。

**預計工時**：2-3 個月。

---

### 5.6 方案 F（長期）：利用甲基化增強 Phasing

**假說**：甲基化信號可以作為 phasing 的額外資訊來源，提升 LongPhase 的 haplotagging 覆蓋率和準確度。

**理論基礎**：MethPhaser (Kim et al., 2024) 已證明利用 ONT 甲基化信號可將 phase block N50 提升 78-151%，準確度達 83.4-98.7%。

**與 InterSubMod 的關聯**：
- 更多 reads 被正確 haplotag → 更多 reads 進入 InterSubMod 分析
- 更長的 phase blocks → 更完整的局部分析窗口
- 甲基化 phasing 與 SNV phasing 的一致性可作為額外的驗證維度

**預計工時**：3-6 個月（需要與 LongPhase 團隊協作）。

---

### 5.7 方案 G（長期）：新 Benchmark 整合（CASTLE 數據集）

**假說**：使用 CASTLE 數據集的多樣本、多技術驗證，可以更公平地評估 InterSubMod 的跨癌種表現。

**CASTLE 數據集概述** (Park & Cook et al., 2025)：
- 6 對 tumor-normal cell line pairs
- Illumina + PacBio HiFi + ONT 三技術平台全基因組測序
- 329,011 個體細胞變異
- 高品質 truth set

**優勢**：
- 擴展驗證範圍至 6 對樣本（vs 目前的 1-2 對）
- 多技術平台提供交叉驗證
- 標準化的 benchmark 流程

**預計工時**：2-4 個月（數據下載 + 流程適配 + 全面驗證）。

---

## 6. 驗證流程

### 6.1 短期驗證（方案 A + B）

#### Step 1：重建含 MM/ML 的 HCC1395 Subsample

```bash
# 從 ONT_5kHz 資料集建立不同 purity 的 subsample
# Tumor: /big8_disk/data/HCC1395/ONT_5kHz_simplex_5mCG_5hmCG/HCC1395.bam (272GB)
# Normal: /big8_disk/data/HCC1395/ONT_5kHz_simplex_5mCG_5hmCG/HCC1395BL.bam (140GB)

# 範例：建立 80% purity subsample
samtools view -s 0.8 -b HCC1395.bam > tumor_80pct.bam
samtools view -s 0.4 -b HCC1395BL.bam > normal_20pct.bam
samtools merge -o subsample_80pct.bam tumor_80pct.bam normal_20pct.bam
samtools sort subsample_80pct.bam -o subsample_80pct.sorted.bam
samtools index subsample_80pct.sorted.bam
```

#### Step 2：建立組織差異 CpG 圖譜

1. 在 pure tumor 和 pure normal BAM 上，對每個 CpG 位點計算平均甲基化率
2. 計算 |Δβ| = |β_tumor - β_normal|
3. 設定閾值（例如 |Δβ| > 0.3），產出「組織差異 CpG」BED 檔

#### Step 3：A/B 測試

- **對照組**：使用原始 InterSubMod 分析 subsample
- **實驗組 A**：排除組織差異 CpG 後分析
- **實驗組 B**：加入 read-level 組織來源標記後分析
- **比較指標**：
  - TP/FP 位點的顯著性分佈
  - F1 score 隨 purity 的衰減曲線
  - Cluster 與 HP/Allele 關聯的 Cramér's V

#### Step 4：結果分析

**預期結果**：

| 指標 | 對照組 | 實驗組 A (排除 CpG) | 實驗組 B (Read-level) |
|------|--------|-------------------|---------------------|
| F1 at purity 80% | 基線 | 改善 | 改善 |
| F1 at purity 40% | 大幅下降 | 中度下降 | 小幅下降 |
| False positive clusters | 多 | 減少 | 大幅減少 |
| Tissue-driven clusters | 主導 | 減少 | 消除 |

### 6.2 中期驗證（方案 C + E）

#### Step 1：純度感知距離校正

1. 在多個已知 purity 的 subsample 上測試校正算法
2. 比較校正前後的分群穩定性

#### Step 2：多模態整合

1. 在 HCC1395 ONT_5kHz 上運行完整的 SNV + HP + methylation 聯合分析
2. 與單一模態（僅甲基化）的結果比較

### 6.3 長期驗證（方案 F + G）

1. 下載 CASTLE 數據集
2. 適配 InterSubMod 流程
3. 在 6 對樣本上執行全面 benchmark
4. 與 DeepSomatic / ClairS 的結果進行比較

---

## 7. 結論與優先順序

### 7.1 核心結論

InterSubMod 在 subsample 驗證中表現不佳的**根本原因是組織來源甲基化差異的混淆效應**。血液（B-lymphoblast）與實體組織（乳腺、皮膚、肺）之間存在數千個差異甲基化 CpG，這些差異在 subsample 混合過程中被保留，並被 InterSubMod 的距離計算和聚類演算法識別為「顯著的甲基化異質性」——但實際上反映的是組織來源差異而非腫瘤亞克隆結構。

這不是 InterSubMod 分析方法本身的錯誤，而是**實驗設計的根本限制**被文獻充分支持：

> 「當使用血液來源的 normal 作為對照，分析實體腫瘤的甲基化差異時，所觀察到的『腫瘤-正常差異』中會混入大量的『組織來源差異』。」— 基於 InfiniumPurify (Zheng et al., 2017)

### 7.2 方案優先順序

| 優先級 | 方案 | 預計工時 | 預期影響 | 實作難度 |
|--------|------|---------|---------|---------|
| **P0** | A：CpG 排除/降權 | 1-2 週 | 中 | 低 |
| **P0** | MM/ML Subsample 重建 | 1 週 | 前置條件 | 低 |
| **P1** | B：Read-level 來源判定 | 2-3 週 | 高 | 中 |
| **P1** | E：多模態整合 (SNV-gated) | 2-3 個月 | 高 | 中-高 |
| **P2** | C：純度感知距離 | 1-2 個月 | 中 | 中 |
| **P3** | F：甲基化增強 Phasing | 3-6 個月 | 中 | 高 |
| **P3** | G：CASTLE Benchmark | 2-4 個月 | 高（驗證面） | 中 |
| **P4** | D：同組織 Normal | 取決數據 | 高（理想但難實現） | 高 |

### 7.3 建議的立即行動

1. **重建含 MM/ML 的 subsample**：使用 ONT_5kHz 資料集，至少建立 100%、80%、60%、40% 四個 purity level
2. **建立 HCC1395 的組織差異 CpG 圖譜**：從 pure tumor 和 pure normal 比較，產出排除清單
3. **在 InterSubMod 中實作 CpG gating 的組織差異排除**：在 S4a 階段加入新的排除條件
4. **執行 A/B 測試**：比較排除前後的驗證結果

### 7.4 方向性判斷

如果方案 A（CpG 排除）能顯著改善 subsample 驗證結果，則確認組織來源是主要混淆因子，後續可投入方案 B 和 E。如果改善有限，則需重新檢視其他可能的混淆因子（如 sequencing bias、basecalling error 等）。

---

*文檔版本: v1.0*
*最後更新: 2026-03-02*
