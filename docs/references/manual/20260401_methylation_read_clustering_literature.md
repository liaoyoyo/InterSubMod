# Read-Level 甲基化模式用於 Variant Classification 文獻調查報告

<!--
建立時間: 2026-04-01 23:30
目標: 調查 read-level 甲基化模式（而非 per-site 統計量）用於 variant classification 的現有方法、工具、與適用性
處理範圍: 學術文獻搜尋與整理，涵蓋 read-level 甲基化分類、甲基化分群、haplotype-specific 甲基化、ONT 整合方法
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
  - docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md
  - docs/references/20260401_ASM偵測方法與預期比例文獻調查_01.md
  - docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
-->

## 1. 搜尋概述

- **核心問題**: 能否利用 **read-level** 甲基化模式（而非 per-site 統計量如 epipolymorphism、mean methylation）來區分 germline FP 和 somatic TP，特別是在 ONT tumor-only 場景下？
- **搜尋時間**: 2026-04-01
- **資料來源數**: 25+ 篇論文/預印本/工具文件
- **關鍵詞**: read-level methylation pattern, variant classification, allele-specific methylation, nanopore methylation read clustering, haplotype methylation tumor, methylation phasing variant filtering, single-read methylation, epigenetic variant classification long reads, tumor read classification methylation
- **背景**: O11-O13 已證明 **per-site** 甲基化統計量（epipolymorphism、cross-region correlation、LOH scenario discrimination）在控制 confounds 後 AUC < 0.60，正式排除。本次搜尋聚焦於利用 **read-level pattern**（整條 read 上多個 CpG 的聯合甲基化狀態）的新方向。

---

## 2. 核心發現：Read-Level 方法分類

### 概念區分：Per-Site vs. Read-Level

| 分析層級 | 定義 | 我們已有的結果 | 本次搜尋焦點 |
|----------|------|----------------|--------------|
| **Per-site** | 單個 CpG 或 region 的聚合統計量（mean、entropy、epipolymorphism） | O11-O13 NEGATIVE（AUC < 0.60） | 不是 |
| **Read-level** | 整條 read 上多個 CpG 的聯合甲基化狀態向量，作為分類特徵 | ISM 已有距離矩陣基礎設施 | **是** |
| **Model-based read classification** | 用 DL/ML 模型對每條 read 進行 tumor/normal 分類 | 尚未嘗試 | **是** |

---

## 3. 方法與工具詳細分析

### 3.1 ROCIT — Read-Level Tumor/Normal Classification（最直接相關）

**論文**: Baker et al., "Genome-wide classification of tumor-derived reads from bulk long-read sequencing," bioRxiv, 2026-03-05. DOI: 10.64898/2026.03.03.709085

**核心思路**:
- Transformer-based model，對每條 long read 進行 tumor vs. non-tumor 二元分類
- **訓練資料生成**: 利用已知 somatic mutations 作為 training labels（攜帶 somatic ALT allele 的 reads 標記為 tumor，其他為 normal）
- **分類特徵**: 每條 read 上的 CpG 甲基化模式（read-level methylation pattern）
- **不需要 matched normal**: 透過 somatic mutations 自動生成 labels，無需正常組織
- **不需要預先定義 DMR**: 全基因組適用，不依賴預先辨識的 tumor DMR

**技術細節**:
- 平台: PacBio HiFi（作者表示 ONT 也適用，需 base calls + methylation）
- 癌種驗證: 前列腺癌、卵巢癌
- 全基因組分類準確度: 高（具體數字需查閱全文 PDF）
- Perturbation analysis: 分析了哪些 CpG 位置對分類最有資訊量

**與本研究的適用性評估**:

| 維度 | 評估 | 說明 |
|------|------|------|
| 概念匹配度 | **高** | 直接證明 read-level methylation patterns 足以區分 tumor/normal reads |
| 目標差異 | **中** | ROCIT 分類 tumor/normal **reads**，我們要分類 germline/somatic **variants** |
| Tumor-only 適用性 | **中-低** | 需要「已知 somatic mutations」做 training labels，但在 TO 場景中 somatic mutations 本身不確定（chicken-and-egg problem） |
| 技術遷移難度 | **中** | 需要 Transformer 訓練基礎設施；ISM 的 read-level 矩陣可作為特徵輸入 |
| ONT 相容性 | **高** | 理論上完全相容，僅需 BAM + MM/ML tags |

**關鍵洞察**: ROCIT 的核心假設是「腫瘤細胞的表觀遺傳景觀（特別是全域低甲基化）使得 tumor reads 與 normal reads 在 CpG 甲基化模式上有系統性差異」。這與我們的假說一致：germline variants 上的 reads 可能有「正常」甲基化模式，somatic variants 上的 reads 可能有「腫瘤特異性」甲基化模式。但直接套用有 chicken-and-egg 問題。

---

### 3.2 MethylBERT — Read-Level Methylation Pattern + Tumor Deconvolution

**論文**: Jeong et al., "MethylBERT enables read-level DNA methylation pattern identification and tumour deconvolution using a Transformer-based model," Nature Communications 16, 2025. DOI: 10.1038/s41467-025-55920-z

**核心思路**:
- 修改版 BERT 模型，以 DNA 3-mer tokens + 甲基化嵌入為輸入
- 對每條 read 進行 tumor/normal cell type 分類
- 輸出後驗機率，用於 Bayesian 推斷 tumor purity

**技術細節**:
- 架構: 12 encoder layers, 768-hidden dimension（500bp reads 使用 6 layers）
- 甲基化編碼: 0 = unmethylated CpG, 1 = methylated CpG, 2 = non-CpG
- 預訓練: Masked Language Modeling on 3-mers
- 微調: 二元分類（tumor/normal）
- 準確度: > 0.95（150bp reads, 低至 10 reads coverage）
- 癌種: DLBCL, CRC, PDAC
- 優於: CancerDetector, DISMIR, Houseman method

**與本研究的適用性評估**:

| 維度 | 評估 | 說明 |
|------|------|------|
| 概念匹配度 | **高** | 證明 read-level methylation + sequence context 足以區分 cell type |
| 目標差異 | **中** | MethylBERT 做 cell-type deconvolution，不直接做 variant classification |
| Tumor-only 適用性 | **中** | 需要 normal cell-type methylation atlas 作為 reference（不需 matched normal） |
| 技術遷移可行性 | **中-高** | 可借鑑 3-mer + methylation 的 encoding 策略；但需要大規模訓練集 |
| Read 長度需求 | **靈活** | 150bp 和 500bp 均有效；ONT long reads 理論上更有優勢 |

**關鍵洞察**: MethylBERT 證明即使在 150bp 短 reads 上，read-level methylation pattern 就足以區分 tumor/normal。ONT long reads（數 kb 到數十 kb）提供的 CpG 數量遠多於此，理論上分類 power 更強。但需要 **normal methylation reference** 作為對比基準。

---

### 3.3 Alpha — Read-Level Methylation Deconvolution

**論文**: "Read-level DNA methylation deconvolution enhances circulating tumor DNA detection," Briefings in Bioinformatics 26(5), 2025. DOI: 10.1093/bib/bbaf551

**核心思路**:
- 結合無偏基因組分割（dynamic programming）與 read-level 甲基化信號偵測
- 對每條 read 計算 alpha 值（聚合鄰近 CpG 的甲基化水平）
- 使用 NNLS deconvolution 估計 tumor cell fraction

**技術細節**:
- 分割策略: Maximum Likelihood 演算法，最小化 segment 內甲基化變異
- 每個 segment 至少 4 個 CpG sites
- Alpha marker regions 平均 ~400bp
- 準確度: MAE=0.004（vs CelFEER 0.014, UXM 0.13）
- 靈敏度: 可偵測低至 0.5% 的 tumor fraction
- 比較: 優於 CelFEER 和 UXM

**與本研究的適用性評估**:

| 維度 | 評估 | 說明 |
|------|------|------|
| 概念匹配度 | **中** | 證明 read-level aggregated methylation 有效，但應用場景為 cfDNA deconvolution |
| 直接適用性 | **低** | 專為 cfDNA 設計，非直接用於 solid tumor variant classification |
| 可借鑑之處 | **中** | alpha 值的計算概念（per-read multi-CpG aggregation）可用於 ISM 特徵工程 |

---

### 3.4 CancerDetector / DISMIR — Per-Read Cancer Detection

**論文**:
- Li et al., "CancerDetector: ultrasensitive and non-invasive cancer detection at the resolution of individual reads," Nucleic Acids Research 46(15), 2018.
- Yin et al., "DISMIR: Deep learning-based noninvasive cancer detection by integrating DNA sequence and methylation information of individual cell-free DNA reads," Briefings in Bioinformatics 22(6), 2021.

**核心思路**:
- **CancerDetector**: 利用鄰近 CpG sites 的甲基化狀態 local correlation 預測每條 read 的來源（tumor vs normal）
- **DISMIR**: 深度學習模型，整合 DNA sequence + methylation state，預測每條 cfDNA read 的來源。引入「switching region」概念定義 cancer-specific DMR

**與本研究的適用性評估**:

| 維度 | 評估 | 說明 |
|------|------|------|
| 概念匹配度 | **中** | 同樣是 per-read classification，但專為 cfDNA 設計 |
| 直接適用性 | **低** | 依賴 bisulfite sequencing 而非 ONT native methylation |
| 可借鑑之處 | **中** | local CpG correlation 概念可借鑑；switching region 概念可能有用 |

---

### 3.5 CPEL — Ising Model for Haplotype-Dependent ASM

**論文**: Jenkinson et al., "Detection of haplotype-dependent allele-specific DNA methylation in WGBS data," Nature Communications 11, 2020. DOI: 10.1038/s41467-020-19077-1

**核心思路**:
- 使用一維 Ising model（統計物理模型）建模多 CpG 位點的聯合甲基化分布
- 計算三個統計量: MML（mean methylation level）、NME（normalized methylation entropy）、JSD（Jensen-Shannon distance between alleles）
- **關鍵發現**: **96.19% 的顯著 haplotype-dependent ASM 僅表現為 entropy imbalance**（mean 相同但 pattern variability 不同）

**技術細節**:
- 需要 phased WGBS data
- 可偵測 mean-based 方法完全遺漏的 entropy ASM
- 整合 haplotype 內多個 CpG 的 joint probability distribution

**與本研究的適用性評估**:

| 維度 | 評估 | 說明 |
|------|------|------|
| 概念匹配度 | **高** | 直接處理 read-level multi-CpG pattern，且發現 entropy > mean |
| 與 ISM 的關係 | **互補** | ISM 的 NHD/distance matrix 已部分捕捉 dispersion 信號，但未明確分離 entropy vs. mean |
| 關鍵啟示 | **非常重要** | 如果 96% 的 ASM 是 entropy imbalance，那麼 per-site mean-based 方法（O11 的 epipolymorphism）失敗是預期的；需要的是 read-level joint distribution 比較 |
| Tumor-only 限制 | **中** | 需要 phased data（LongPhase-TO 可提供），但計算 JSD 需要足夠 reads per allele |

**關鍵洞察**: Jenkinson 2020 的發現直接解釋了我們 O11 為何失敗——epipolymorphism 是 per-site entropy，而真正的 ASM 信號在 **多 CpG 聯合分布** 的 entropy imbalance 中。ISM 的 PERMANOVA on distance matrix 理論上可以部分捕捉這個信號（因為 NHD distance 反映 pattern 差異），但需要更精確的 entropy decomposition（mean vs. dispersion）。

---

### 3.6 甲基化增強 Phasing 工具族

#### 3.6.1 MethPhaser

**論文**: Cheng et al., "MethPhaser: methylation-based long-read haplotype phasing of human genomes," Nature Communications 15, 2024. DOI: 10.1038/s41467-024-49588-0

- 利用 ONT 甲基化信號延伸 SNV-based phasing
- Phase length N50 增加 78-151%，準確度 83.4-98.7%
- 發現 haplotype-specific methylation 在人類基因組中廣泛存在

#### 3.6.2 HapBridge

**論文**: "HapBridge: A Methylation-Guided Approach for Correcting Switch Errors and Bridging Phased Blocks in Long-Read Phasing," bioRxiv, 2025-11. DOI: 10.1101/2025.11.07.687303

- 利用甲基化修正 switch errors 並橋接 phased blocks
- Switch errors 降低 3.07-18.72%，N50 提升 5.84-68.61%
- 優於 MethPhaser

#### 3.6.3 LongHap

**論文**: "Harnessing methylation signals inherent in long-read sequencing data for improved variant phasing," bioRxiv, 2026-03-11.

- 第一個統一整合 sequence + 5mC methylation 的 read-based phasing 方法
- 優於 WhatsHap、HapCUT2、LongPhase、MethPhaser
- 更低 switch error rate，更大 phase block contiguity

**與本研究的適用性評估**:

| 維度 | 評估 | 說明 |
|------|------|------|
| 直接用於 variant classification | **低** | 這些工具改善 phasing，不直接分類 variants |
| 間接價值 | **高** | 更好的 phasing -> 更準確的 HP tags -> ISM 的 HP-based 分群更可靠 |
| Tumor-only 場景 | **中** | 甲基化 phasing 不依賴 matched normal，可在 TO 中使用 |
| 與 LongPhase-TO 的關係 | **互補** | LongPhase-TO 使用 variant-based phasing；加入甲基化可能改善 LOH 區域的 phasing |

---

### 3.7 NanoMethPhase / NANOME — Allele-Specific Methylation Detection

**論文**:
- Akbari et al., "Megabase-scale methylation phasing using nanopore long reads and NanoMethPhase," Genome Biology 22, 2021. DOI: 10.1186/s13059-021-02283-5
- "NANOME: A Nextflow pipeline for haplotype-aware allele-specific consensus DNA methylation detection by nanopore long-read sequencing," bioRxiv, 2025. DOI: 10.1101/2025.06.29.662079

**NanoMethPhase**:
- ONT 專用 allele-specific methylation 偵測
- ~10x coverage 即可偵測全基因組 ASM
- 涵蓋 26.5M autosomal CpGs（95% of human autosomal methylome）

**NANOME**:
- XGBoost ensemble model 整合 Megalodon, Nanopolish, DeepSignal
- Nextflow pipeline，支援 haplotype-aware ASM
- 不需要 parental genotype（單樣本即可）

**與本研究的適用性評估**:

| 維度 | 評估 | 說明 |
|------|------|------|
| ASM 偵測能力 | **高** | 提供全基因組 haplotype-resolved ASM 數據 |
| 直接用於 variant classification | **低** | 偵測 ASM 位置，不直接分類 variants |
| 作為特徵來源 | **中-高** | ASM 結果可作為 ISM 的額外特徵（某位點是否在 normal reference 中已知為 ASM） |

---

### 3.8 pycoMeth — ONT Differential Methylation Testing

**論文**: Leger et al., "pycoMeth: a toolbox for differential methylation testing from Nanopore methylation calls," Genome Biology 24, 2023. DOI: 10.1186/s13059-023-02917-w

- ONT 專用 differential methylation testing
- MetH5 高效儲存格式
- 支援 haplotype-aware, multi-sample consensus segmentation
- 段落可測試 differential methylation 或 allele-specific methylation

**與本研究的適用性**: 提供 ONT 原生的 ASM 測試框架，可作為 ISM 外部驗證工具。

---

### 3.9 modkit — ONT 官方甲基化工具

**工具**: nanoporetech/modkit (GitHub)

- ONT 官方 modified base 處理工具
- `modkit pileup`: 產生 BEDMethyl 格式的甲基化摘要
- **Phased output**: 可分別產生 hp1.bedmethyl 和 hp2.bedmethyl
- `modkit dmr`: Likelihood ratio + HMM 進行 differential methylation 分析
- 支援 haplotype-stratified 甲基化分析

**與本研究的適用性**: 作為 ISM 的上游/驗證工具，提供 per-haplotype 甲基化摘要。

---

### 3.10 LRSomatic — 整合甲基化的 Somatic Calling Pipeline

**論文**: "LRSomatic: a highly scalable and robust pipeline for somatic variant calling in long-read sequencing data," bioRxiv, 2026-02-28. DOI: 10.64898/2026.02.26.707772

- Nextflow-based pipeline，支援 PacBio HiFi + ONT
- 支援 paired tumor-normal 和 **tumor-only** 設計
- 整合 Fiber-seq 甲基化/chromatin accessibility 數據
- 可偵測 haplotype-specific methylation（如 imprinted loci）

**與本研究的適用性**: 證明 somatic calling pipeline 正在向甲基化整合發展，但目前主要用於 annotation 而非 variant filtering。

---

### 3.11 crossNN — 跨平台甲基化分類

**論文**: "crossNN is an explainable framework for cross-platform DNA methylation-based classification of tumors," Nature Cancer 6, 2025. DOI: 10.1038/s43018-025-00976-5

- 可分類 > 170 種腫瘤類型
- 使用 randomly masked 訓練策略處理 sparse methylomes
- 支援不同平台和覆蓋深度

**與本研究的適用性**: 概念上支持甲基化用於分類，但目標是 tumor type classification 而非 variant classification。不直接適用，但 masked training 策略可借鑑。

---

### 3.12 Loyfer 2025 — ASM Atlas

**論文**: Rosenski, Peretz, Magenheim, Loyfer et al., "Atlas of imprinted and allele-specific DNA methylation in the human body," Nature Communications, 2025. DOI: 10.1038/s41467-025-57433-1

- 39 種 normal human cell types 的 deep WGS
- 325K bimodal methylation regions（6% genome, 11% CpGs）
- 34K regions 有 genetic variant-driven ASM
- 460 imprinting regions（含 78 已知印記基因）

**與本研究的適用性**: 提供 **normal ASM reference**。如果某個 variant 位點在 normal atlas 中已知為 ASM，則該 variant 更可能是 germline。這是 Phase 2 Normal Methylation Reference 方向的基礎。

---

## 4. Read-Level vs. Per-Site 方法的技術需求比較

| 需求 | Per-Site（已失敗） | Read-Level Pattern（新方向） |
|------|-------------------|---------------------------|
| BAM MM/ML tags | 需要 | 需要 |
| Reference genome | 需要 | 需要 |
| Phased BAM | 有助但非必要 | 多數方法需要（HP-stratified analysis） |
| Normal reference | 不需要 | ROCIT 不需要；MethylBERT 需要 atlas |
| Training data | 不需要 | ROCIT/MethylBERT 需要；CPEL 不需要 |
| GPU | 不需要 | Transformer 方法需要；CPEL 不需要 |
| 每位點 CpG 數量 | 任意 | 越多越好（> 4-10 個 CpG per read） |
| 每位點 Read 數量 | > 5 | > 10（分群需要統計 power） |
| ISM 已有基礎設施 | 完全重疊 | read-level 矩陣已有，需新增分類器 |

---

## 5. 與 ISM 場景的適用性評估（重點）

### 5.1 直接可行的方向

| 方向 | 方法來源 | ISM 改造需求 | 預期難度 | 預期效果 |
|------|---------|-------------|---------|---------|
| **A: Entropy decomposition** | CPEL (Jenkinson 2020) | 在 ISM 距離矩陣基礎上，分離 MML vs NME 信號 | **低** | **中-高**（96% ASM 是 entropy imbalance，ISM PERMANOVA 已部分捕捉） |
| **B: Per-read alpha score** | Alpha (2025) | 計算每條 read 的 multi-CpG aggregated methylation score | **低** | **中**（可作為新 ISM 欄位） |
| **C: Normal methylation reference lookup** | Loyfer 2025 atlas | 查詢 variant 位點是否在 normal ASM atlas 中 | **低** | **中**（需下載 atlas 數據） |

### 5.2 需要顯著投入的方向

| 方向 | 方法來源 | ISM 改造需求 | 預期難度 | 預期效果 |
|------|---------|-------------|---------|---------|
| **D: Transformer read classifier** | ROCIT/MethylBERT | 建立 DL 訓練基礎設施，收集訓練集 | **高** | **高**（ROCIT 證明可行） |
| **E: Methylation-guided phasing** | MethPhaser/HapBridge/LongHap | 整合到 LongPhase-TO pipeline | **中-高** | **中**（改善 phasing 間接改善 ISM） |

### 5.3 Tumor-Only 場景的特殊限制

| 限制 | 說明 | 影響的方法 |
|------|------|----------|
| **無 matched normal** | 無法直接計算 tumor-specific DMR | ROCIT training（但可用 external atlas） |
| **Somatic variants 不確定** | TO 的 TP/FP 標籤本身不可靠 | ROCIT training labels, evaluation |
| **Germline FP 是罕見 germline** | 逃脫 PoN 的 FP 在所有特徵上與 TP 相似 | 所有方法（根本限制） |
| **LOH 區域 phasing 不完整** | LOH 區域只有一個 allele，HP-based 分群退化 | CPEL entropy, HP-stratified methods |

---

## 6. 衝突觀點分析

### 6.1 共識觀點

| 觀點 | 支持證據 | 信心度 |
|------|---------|--------|
| Read-level methylation patterns 有足夠資訊量區分 cell type | ROCIT, MethylBERT, CancerDetector, DISMIR 均成功 | **高** |
| Germline variants 有更強更穩定的 ASM | mQTL 文獻, Do et al. 2020, 我們的 O11-ASM 分析 | **高** |
| 96% ASM 是 entropy imbalance 而非 mean shift | Jenkinson 2020 CPEL | **高**（單一研究，但方法嚴謹） |
| ONT long reads 理論上比 short reads 更適合 read-level methylation analysis | 更多 CpG per read, native methylation | **高** |

### 6.2 衝突/未驗證觀點

| 觀點 A | 觀點 B | 可能原因 |
|--------|--------|----------|
| ROCIT: 甲基化可區分 tumor/normal reads 全基因組 | 我們 O11-O13: 甲基化無法區分 TP/FP | ROCIT 做的是 cell-type classification（全域模式差異），我們做的是 variant-level classification（局部窗口）。兩者任務不同。 |
| MethylBERT: 150bp 足夠區分 tumor/normal | ISM: 5-10kb 窗口仍無法區分 TP/FP | MethylBERT 使用的是 pre-trained 全基因組模式，ISM 使用的是 per-site 統計量。**分析粒度和方法完全不同。** |
| Do et al.: Germline ASM > Somatic ASM | 我們 O1-O13: 甲基化特徵無法區分 | ISM 測量的是 cluster-based 特徵（距離、顯著性），不是直接的 ASM 效應量。且 FP 中混合了 germline + artifact，稀釋了信號。 |
| Jenkinson 2020: entropy imbalance 是主要 ASM 形式 | ISM PERMANOVA: 已部分捕捉 dispersion | ISM 的 NHD distance 隱含 entropy 資訊，但 PERMANOVA p-value 是聚合統計，可能遺失個別 haplotype 的 entropy 差異。需要明確的 per-haplotype NME 計算。 |

---

## 7. 資料來源評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| ROCIT (Baker et al., 2026) | 預印本 | **中-高** | 未 peer review，但方法新穎；待正式發表確認 |
| MethylBERT (Jeong et al., 2025) | 正式論文 (Nature Communications) | **高** | 已 peer review，有 GitHub 實作 |
| Alpha (2025) | 正式論文 (Briefings in Bioinformatics) | **高** | 已 peer review |
| CPEL (Jenkinson et al., 2020) | 正式論文 (Nature Communications) | **高** | 已 peer review，引用數高 |
| MethPhaser (Cheng et al., 2024) | 正式論文 (Nature Communications) | **高** | 已 peer review |
| HapBridge (2025) | 預印本 | **中** | 未 peer review |
| LongHap (2026) | 預印本 | **中** | 未 peer review |
| NanoMethPhase (Akbari et al., 2021) | 正式論文 (Genome Biology) | **高** | 已 peer review |
| NANOME (2025) | 預印本/正式論文 | **中-高** | 近期發表 |
| LRSomatic (2026) | 預印本 | **中** | 未 peer review |
| Do et al. (2020) | 正式論文 (Genome Biology) | **高** | 已 peer review |
| Loyfer et al. (2025) | 正式論文 (Nature Communications) | **高** | 已 peer review |
| crossNN (2025) | 正式論文 (Nature Cancer) | **高** | 已 peer review |
| Fu, Timp, Sedlazeck (2025) | 綜述 (Nature Reviews Genetics) | **高** | 權威綜述 |

---

## 8. 與已失敗方向的差異化分析

| 已失敗方向 | 為何失敗 | Read-level 方法的不同之處 |
|-----------|---------|-------------------------|
| **O11: epipolymorphism** | per-site entropy 被 n_reads confound；殘差化後 AUC=0.530 | CPEL 的 NME 是 per-haplotype joint entropy（多 CpG），非 per-site；且 96% ASM 在此維度 |
| **O12: LOH scenario** | AlleleDelta 被 AF confound；CramersV 有 collider bias | Read-level classifier（ROCIT/MethylBERT）使用全域模式，不依賴單一 site 的 delta |
| **O13: cross-region correlation** | Shared read count confound | Read-level 方法對每條 read 獨立分類，不依賴跨區域 read sharing |
| **G1-G7: VCF features** | 殘餘 germline FP 在 VCF 特徵上與 TP 本質相似 | 甲基化是 VCF 以外的獨立資訊維度；ROCIT 不使用 VCF 特徵 |
| **ML combination (LR/RF)** | ISM 甲基化特徵 additive value < +0.013 AUC | ISM 使用的是 per-site aggregate features；read-level pattern 是更豐富的特徵空間 |

**核心差異**: 之前的失敗都是在 **per-site aggregate** 層級上操作（mean methylation, entropy, distance）。Read-level 方法直接處理 **每條 read 的完整 CpG 甲基化向量**，保留了更多資訊。ROCIT/MethylBERT 的成功證明這個額外資訊確實存在。

---

## 9. 建議行動（按優先順序）

### P0: 可立即執行（使用 ISM 現有基礎設施）

**A. Entropy Decomposition（MML vs. NME 分離）**

基於 Jenkinson 2020 的發現（96% ASM 是 entropy imbalance），在 ISM 的 HP-stratified methylation 矩陣上：
1. 計算每個 haplotype 的 NME（normalized methylation entropy）
2. 計算兩個 haplotype 之間的 NME 差異（delta NME）
3. 與現有 MML-based 特徵（AlleleDelta）分開評估 AUC
4. 預期: 如果 Jenkinson 的發現成立，delta NME 可能比 AlleleDelta 更有鑑別力

**技術需求**: 僅需修改 ISM 的 Python analysis scripts，不需改 C++ 核心
**風險**: 中。ISM 的 CpG 窗口可能不夠大（~5-10 CpG per region），降低 entropy 估計精度

**B. Per-Read Alpha Score**

借鑑 Alpha 方法，為每條 read 計算 multi-CpG aggregated methylation score：
1. 從 ISM 的 methylation.csv 提取 per-read methylation vector
2. 計算 per-read alpha（CpG 甲基化水平的聚合）
3. 比較 TP/FP regions 的 per-read alpha 分布差異
4. 可與 HP 分層結合

### P1: 短期可執行（需少量新數據/工具）

**C. Normal ASM Reference Lookup**

利用 Loyfer 2025 atlas：
1. 下載 325K bimodal regions + 34K ASM regions 數據
2. 對 ISM 分析的每個 variant 位點，查詢是否落在已知 normal ASM region 內
3. 如果 FP（germline）更常落在 normal ASM regions 中，可作為 annotation feature
4. 不需改 ISM 核心代碼

**風險**: 低。Atlas 基於 array/WGBS，與 ONT 平台的 CpG 覆蓋可能不完全重疊

### P2: 中期探索方向（需要較大投入）

**D. Lightweight Read Methylation Classifier**

受 ROCIT/MethylBERT 啟發，但簡化為 ISM 場景：
1. 不建 Transformer，改用簡單 classifier（logistic regression / random forest）
2. 特徵: 每條 read 上的 CpG methylation vector（ISM 已有）
3. 訓練標籤: 使用 paired 樣本的 truth set（paired TP 位點的 reads = "somatic context", paired FP 位點的 reads = "germline context"）
4. 在 paired 上訓練，測試能否 transfer 到 TO

**風險**: 高。可能重蹈 O11 覆轍（read-level 甲基化本身在局部窗口可能資訊不足）

**E. Methylation-Enhanced Phasing for ISM**

整合 MethPhaser/HapBridge 概念到 LongPhase-TO：
1. 甲基化信號輔助 phasing，改善 LOH 區域的 HP assignment
2. 更好的 HP -> 更好的 ISM HP-stratified 分析
3. 需要與 LongPhase-TO 開發者（同實驗室）協調

---

## 10. 綜述文獻推薦

以下綜述提供了本領域的全面視角：

- **Fu, Timp & Sedlazeck, "Computational analysis of DNA methylation from long-read sequencing," Nature Reviews Genetics 26, 2025** — 最權威的長讀長甲基化計算方法綜述
- **Do et al., "Allele-specific DNA methylation is increased in cancers," Genome Biology 21, 2020** — Cancer ASM 的全面分析
- **Jenkinson et al., "Detection of haplotype-dependent ASM," Nature Communications 11, 2020** — Entropy imbalance ASM 的發現
- **Rosenski/Loyfer et al., "Atlas of imprinted and ASM," Nature Communications, 2025** — Normal human ASM atlas

---

## 11. 結論

### 核心判斷

1. **Read-level 方法與 per-site 方法是本質不同的分析層級**。ROCIT/MethylBERT 的成功不與 O11-O13 的失敗矛盾——前者利用全域 read-level pattern，後者使用 per-site aggregate statistics。

2. **Entropy decomposition 是最有希望的低成本方向**。Jenkinson 2020 發現 96% ASM 是 entropy imbalance，而 ISM 現有的 per-site epipolymorphism（O11）恰好測量的不是這個。在 ISM 框架內加入 per-haplotype NME 計算成本低，值得優先驗證。

3. **Normal methylation reference 是 TO 場景的根本突破口**。無論使用什麼方法，TO 的核心問題是缺少 normal baseline。Loyfer 2025 atlas 提供了外部 normal reference，可能比任何 tumor-internal 特徵更有效。

4. **Transformer-based read classifier 是高投入高潛力方向**。ROCIT 直接證明 read-level methylation pattern 足以區分 tumor/normal reads，但在 TO 場景存在 training label 的 chicken-and-egg 問題。可考慮在 paired 樣本上訓練、TO 上 inference 的 transfer learning 策略。

5. **ISM 已有的 distance matrix 基礎設施是獨特優勢**。文獻中沒有其他工具同時提供: (a) per-SNV anchored methylation matrix, (b) read-read distance matrix, (c) HP/Allele-stratified PERMANOVA。重點是在此基礎上加入 entropy decomposition 和 read-level classification。

### 與研究策略的對齊

本次調查支持 CURRENT_FOCUS 中「Phase 2 Normal Methylation Reference」的策略定位，並提供了具體的技術路徑。最具前景的短期行動是 P0-A（entropy decomposition）和 P1-C（normal ASM reference lookup），兩者均可在現有 ISM 基礎設施上低成本實現。

---

## 12. 參考文獻完整列表

1. Baker TM et al. Genome-wide classification of tumor-derived reads from bulk long-read sequencing. bioRxiv (2026). DOI: 10.64898/2026.03.03.709085 | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.03.709085v1) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12991090/)
2. Jeong S et al. MethylBERT enables read-level DNA methylation pattern identification and tumour deconvolution. Nature Communications 16 (2025). DOI: 10.1038/s41467-025-55920-z | [Nature Comms](https://www.nature.com/articles/s41467-025-55920-z) | [GitHub](https://github.com/CompEpigen/methylbert)
3. Read-level DNA methylation deconvolution enhances circulating tumor DNA detection. Briefings in Bioinformatics 26(5) (2025). DOI: 10.1093/bib/bbaf551 | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12536878/)
4. Jenkinson G et al. Detection of haplotype-dependent allele-specific DNA methylation in WGBS data. Nature Communications 11 (2020). DOI: 10.1038/s41467-020-19077-1 | [Nature Comms](https://www.nature.com/articles/s41467-020-19077-1) | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7567826/)
5. Cheng YY et al. MethPhaser: methylation-based long-read haplotype phasing of human genomes. Nature Communications 15 (2024). DOI: 10.1038/s41467-024-49588-0 | [Nature Comms](https://www.nature.com/articles/s41467-024-49588-0)
6. HapBridge: A Methylation-Guided Approach for Correcting Switch Errors. bioRxiv (2025). DOI: 10.1101/2025.11.07.687303 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.07.687303v3.full)
7. Harnessing methylation signals inherent in long-read sequencing data for improved variant phasing (LongHap). bioRxiv (2026). | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.11.710820v1.full)
8. Akbari V et al. Megabase-scale methylation phasing using nanopore long reads and NanoMethPhase. Genome Biology 22 (2021). DOI: 10.1186/s13059-021-02283-5 | [Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02283-5)
9. NANOME: A Nextflow pipeline for haplotype-aware allele-specific consensus DNA methylation detection. bioRxiv (2025). DOI: 10.1101/2025.06.29.662079 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.29.662079v1)
10. Li W et al. CancerDetector: ultrasensitive and non-invasive cancer detection at the resolution of individual reads. Nucleic Acids Research 46(15) (2018). DOI: 10.1093/nar/gky423 | [NAR](https://academic.oup.com/nar/article/46/15/e89/5036349)
11. Yin Q et al. DISMIR: Deep learning-based noninvasive cancer detection by integrating DNA sequence and methylation information. Briefings in Bioinformatics 22(6) (2021). DOI: 10.1093/bib/bbab250 | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8575022/)
12. Leger A et al. pycoMeth: a toolbox for differential methylation testing from Nanopore methylation calls. Genome Biology 24 (2023). DOI: 10.1186/s13059-023-02917-w | [Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02917-w)
13. LRSomatic: a highly scalable and robust pipeline for somatic variant calling in long-read sequencing data. bioRxiv (2026). DOI: 10.64898/2026.02.26.707772 | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.02.26.707772v1)
14. crossNN: an explainable framework for cross-platform DNA methylation-based classification of tumors. Nature Cancer 6 (2025). DOI: 10.1038/s43018-025-00976-5 | [Nature Cancer](https://www.nature.com/articles/s43018-025-00976-5)
15. Rosenski J, Peretz A, Magenheim J, Loyfer N et al. Atlas of imprinted and allele-specific DNA methylation in the human body. Nature Communications (2025). DOI: 10.1038/s41467-025-57433-1 | [Nature Comms](https://www.nature.com/articles/s41467-025-57433-1)
16. Loyfer N et al. A DNA methylation atlas of normal human cell types. Nature 613 (2023). DOI: 10.1038/s41586-022-05580-6 | [Nature](https://www.nature.com/articles/s41586-022-05580-6)
17. Do C et al. Allele-specific DNA methylation is increased in cancers. Genome Biology 21 (2020). DOI: 10.1186/s13059-020-02059-3 | [Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02059-3)
18. Fu Y, Timp W, Sedlazeck FJ. Computational analysis of DNA methylation from long-read sequencing. Nature Reviews Genetics 26 (2025). DOI: 10.1038/s41576-025-00822-5 | [Nature Reviews Genetics](https://www.nature.com/articles/s41576-025-00822-5)
19. nanoporetech/modkit. GitHub. | [GitHub](https://github.com/nanoporetech/modkit)
20. Ho ME et al. LongPhase-S: purity estimation and variant recalibration with somatic haplotyping for long-read sequencing. bioRxiv (2025). DOI: 10.1101/2025.11.20.689492 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.20.689492v1.full)
21. ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling. Nature Communications (2025). DOI: 10.1038/s41467-025-64547-z | [Nature Comms](https://www.nature.com/articles/s41467-025-64547-z)
22. Evidence of DNA methylation heterogeneity and epipolymorphism in kidney cancer tissue samples. Oncogene (2025). DOI: 10.1038/s41388-024-03270-3 | [Nature/Oncogene](https://www.nature.com/articles/s41388-024-03270-3)
23. epihet: intra-tumoral epigenetic heterogeneity analysis and visualization. Scientific Reports (2021). DOI: 10.1038/s41598-020-79627-x | [Scientific Reports](https://www.nature.com/articles/s41598-020-79627-x)
