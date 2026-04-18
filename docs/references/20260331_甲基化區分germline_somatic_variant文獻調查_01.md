# 甲基化區分 Germline vs Somatic Variant 文獻調查報告

<!--
建立時間: 2026-03-31 23:00
目標: 調查 tumor-only ONT long-read 情境下，能否利用 CpG 甲基化資訊區分 germline variant 和 somatic variant
處理範圍: 學術文獻搜尋與整理，涵蓋生物學基礎、已有工具、mQTL、ASM、tumor-only pipeline 前沿研究
關聯檔案:
  - docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
  - docs/references/external/2026/03/20260302_tissue_origin_methylation_confounding_literature_review_01.md
-->

## 1. 搜尋概述

- **核心問題**: 在 tumor-only（無 matched normal）的 ONT long-read 情境下，能否利用 CpG 甲基化資訊區分 germline variant 和 somatic variant，從而過濾 germline FP？
- **搜尋時間**: 2026-03-31
- **資料來源數**: 30+ 篇論文/預印本/綜述
- **關鍵詞**: germline somatic methylation discrimination, mQTL, allele-specific methylation, tumor-only variant calling, read-level methylation, long-read sequencing, CpG-SNP, ASM cancer

---

## 2. 核心發現

### 2.1 生物學基礎: Germline 與 Somatic Variant 周圍的甲基化差異

#### 2.1.1 Germline Variant 與甲基化的關係

**mQTL (methylation Quantitative Trait Loci) 的廣泛存在**

Germline variant 對周圍 CpG 甲基化有強烈的 cis-regulatory 影響，這是目前最有力的生物學基礎:

- 大規模研究已鑑定出 **4.7 百萬個 cis-mQTL variant-CpG pairs** 和 63 萬個 trans-mQTL，影響超過 12 萬個 CpG 位點 (Oliva et al., Nature Genetics, 2023)
- 約 **93% 的 mQTL 為 cis 作用**（距離 CpG < 1 Mb），且距離越近效應越強 (Gaunt et al., Genome Biology, 2016)
- cis-mQTL 中約 **12% 的 SNP 直接落在 CpG 位點內**（CpG-SNP），徹底破壞該位點的甲基化能力 (Gaunt et al., 2016)
- mQTL 與 2,254 個 GWAS hits 存在 colocalization，涵蓋 83 種 traits (Oliva et al., 2023)

**CpG-SNP 的特殊角色**

- 38-88% 的 ASM 區域依賴於 CpG 二核苷酸中的雜合 SNP，這些 SNP 直接破壞甲基化潛力 (Shoemaker et al., Genome Research, 2010)
- 在任一細胞系中，23-37% 的雜合 SNP 位點表現出 ASM (Shoemaker et al., 2010)
- 超過 50% 的 ASM 區域含有 SNP (Shoemaker et al., 2010)

**關鍵推論**: Germline heterozygous variant 傾向於在兩個 allele 上產生**不對稱但穩定**的甲基化模式（即 ASM）。這是由遺傳決定的、跨組織穩定的表觀遺傳特徵。

#### 2.1.2 Somatic Variant 與甲基化的關係

**癌症中的全域甲基化改變**

- 癌症細胞的典型特徵是: **全基因組低甲基化**（global hypomethylation）+ **CpG island 局部高甲基化**（focal hypermethylation）(Esteller, Human Molecular Genetics, 2007; Feinberg & Vogelstein, 1983)
- Partially Methylated Domains (PMDs) 中的低甲基化與 **somatic mutation 密度正相關** (Zhou et al., Nature Genetics, 2018; Johnstone et al., Nature Communications, 2019)
- PMD 區域覆蓋基因組的 24.3-63.4%，癌症中低甲基化程度遠高於正常組織 (Johnstone et al., 2019)

**Somatic mutation 與甲基化的因果方向**

- 甲基化的 CpG 位點的 C>T 突變率比未甲基化位點高 ~10 倍（deamination 機制）(Cytosine Methylation Affects Mutability, Genetics, 2020)
- 但相反方向也存在: somatic driver mutations（如 IDH1, SETD2, BRAF）可以驅動整個甲基化景觀的改變 (Wen et al., PLoS Computational Biology, 2017)
- DNA methylation 與 genomic alterations 在 NSCLC 進化中**協同作用** (Tarazona et al., Nature Genetics, 2025)

**Somatic mutation 周圍的甲基化特殊性**

- 甲基化 CpG 的鄰近核苷酸（+/-3 bp）突變率反而**降低** (Cytosine Methylation Affects Mutability, Genetics, 2020)
- 此現象在 germline 和 somatic（癌症）基因組中均存在
- 癌症中 de novo ASM（由 somatic mutation 引起）僅佔所有 ASM 增加的 **6-17%**，數量上遠少於 germline SNP 驅動的 ASM (Do et al., Genome Biology, 2020)

#### 2.1.3 關鍵差異總結

| 特徵 | Germline Variant | Somatic Variant |
|------|-----------------|-----------------|
| ASM 驅動比例 | 主要來源（23-37% het SNPs 有 ASM） | 少量（6-17% somatic mutations 產生 de novo ASM） |
| 甲基化穩定性 | 跨組織穩定、遺傳決定 | 腫瘤特異性、動態變化 |
| CpG-SNP 效應 | 38-88% ASM 由 CpG-SNP 驅動 | 不適用（somatic mutation 不直接破壞 CpG 二核苷酸） |
| 全域背景 | 正常甲基化景觀 | 全域低甲基化 + 局部高甲基化 |
| mQTL 關聯 | 強 cis 效應（93% cis） | 無 mQTL 概念（非遺傳的） |

---

### 2.2 已有方法與工具

#### 2.2.1 ROCIT: 最直接相關的方法

**Baker et al., "Genome-wide classification of tumor-derived reads from bulk long-read sequencing", bioRxiv, 2026**
- GitHub: https://github.com/tobybaker/rocit
- **架構**: Transformer-based model，對每個 read 進行 tumor/non-tumor 分類
- **核心創新**: 利用 somatic mutations 做 training labels，用 read-level methylation patterns 進行分類，**不需要 matched normal tissue**
- **癌種**: Prostate cancer, ovarian cancer
- **技術平台**: PacBio HiFi（但作者表示 ONT 也適用，只需 base calls + methylation）
- **應用**: 可改善 somatic variant calling（具體改善幅度未在摘要中公開）
- **Perturbation analysis**: 分析了哪些 CpG 位置最有信息量
- **與本研究的關聯**: ROCIT 證明了 read-level methylation patterns 足以區分 tumor 和 normal reads。如果 tumor reads 攜帶 somatic variants 而 normal reads 攜帶 germline variants，則此方法在概念上支持利用甲基化區分 germline/somatic。

**局限性**:
- ROCIT 分類的是 "tumor read vs normal read"，不是直接分類 "germline vs somatic variant"
- 需要已知 somatic mutations 做 training labels（chicken-and-egg 問題）
- 在 tumor-only 場景中，somatic mutations 本身就不確定

#### 2.2.2 MethylBERT: Read-level 甲基化分類

**Jeong et al., "MethylBERT enables read-level DNA methylation pattern identification and tumour deconvolution using a Transformer-based model", Nature Communications, 2025**
- GitHub: https://github.com/CompEpigen/methylbert
- **架構**: 修改版 BERT（12 layers, 12 heads, 768-dim），以 3-mer + 甲基化狀態為 input
- **功能**: 識別 tumor-derived reads 並估計 tumor cell fraction
- **性能**: 150bp reads 準確率 >95%，優於 CancerDetector、DISMIR、HMM 方法
- **癌種**: DLBCL, CRC, PDAC, 前列腺癌
- **不需要 matched normal**: 可使用 normal cell-type methylation atlas
- **與本研究的關聯**: 證明 read-level methylation 足以區分 cell type。但 MethylBERT 專注於表觀遺傳分類，不直接處理 variant classification。

#### 2.2.3 NanoMethPhase: Allele-specific Methylation 偵測

**Akbari et al., "Megabase-scale methylation phasing using nanopore long reads and NanoMethPhase", Genome Biology, 2021**
- **功能**: 利用 ONT nanopore long reads 進行 megabase-scale methylation phasing
- **成果**: 可在 ~10x coverage 下偵測全基因組 ASM
- **覆蓋**: 26.5M autosomal CpGs（95% of human autosomal methylome）
- **與本研究的關聯**: 提供了在 ONT 平台上偵測 ASM 的技術基礎。若 germline variants 有更穩定的 ASM 模式，此工具可用於特徵提取。

#### 2.2.4 MethPhaser: 甲基化增強 Phasing

**MethPhaser, "Methylation-based long-read haplotype phasing of human genomes", Nature Communications, 2024**
- **功能**: 利用甲基化信號延伸 SNV-based phasing
- **性能**: Phase length N50 增加 78-151%，phasing accuracy 83.4-98.7%
- **與本研究的關聯**: 改善 phasing 可間接幫助 germline/somatic 區分（因為 germline het variants 應在兩個 haplotype 上均勻分佈）

#### 2.2.5 MethylPurify: 無需 Normal 的腫瘤純度估計

**Zheng et al., "MethylPurify: tumor purity deconvolution and differential methylation detection from single tumor DNA methylomes", Genome Biology, 2014**
- **功能**: 從單一腫瘤 bisulfite-seq 樣本估計 purity 並找出 DMR
- **不需要 matched normal**: 完全基於甲基化模式的雙模態分佈
- **與本研究的關聯**: 證明甲基化資訊單獨即可分離 tumor/normal 成分

---

### 2.3 Allele-Specific Methylation (ASM) 在 Germline/Somatic 區分中的潛力

#### ASM 在 Germline 更常見的證據

- **正常組織**: 23-37% 的 germline heterozygous SNPs 表現 ASM (Shoemaker et al., 2010)
- **癌症組織**: ASM 頻率比正常組織高 **5-9 倍** (Do et al., Genome Biology, 2020):
  - Multiple myeloma: 5x
  - B cell lymphomas: 8.5x
  - Glioblastoma: 9x
- **癌症中 ASM 增加的來源**: 72-76% 為 allele-specific loss of methylation，49% 在 GBM (Do et al., 2020)
- **De novo ASM（somatic mutation 驅動）**: 僅佔 6-17% (Do et al., 2020)

#### 利用 ASM 區分 Germline/Somatic 的理論基礎

**支持方向:**
1. Germline het variants 通過 mQTL/CpG-SNP 機制產生穩定的 ASM
2. 這種 ASM 應在 tumor 和 normal cells 中一致
3. Somatic variants 通常不會產生穩定的 ASM（除非恰好影響 TF binding site）

**反對方向（混淆因素）:**
1. 癌症本身大幅增加 ASM（5-9x），可能掩蓋 germline 的 ASM signal
2. LOH 消除了一個 allele，使 germline het variant 表現為 homozygous，ASM 消失
3. Tumor 甲基化的 stochastic 變異很大，降低 ASM 偵測的 signal-to-noise ratio
4. 癌症中 allele switching（ASM 方向翻轉）比正常組織多 3x（43% vs 14%）(Do et al., 2020)

---

### 2.4 Tumor-Only 情境的特殊挑戰與前沿研究

#### 2.4.1 cfDNA / Liquid Biopsy 的經驗

液態活檢也面臨無 matched normal 的情境:
- **CUPiD**: 利用 cfDNA methylation patterns 進行 tissue-of-origin 分類，29 種 tumor classes，sensitivity 84.6% (Tao et al., Nature Communications, 2024)
- **Nanopore cfDNA methylation**: ONT 可在 cfDNA 中同時偵測 CNA + methylation + fragmentation (Katsman et al., Genome Biology, 2022)
- **關鍵洞見**: 在 liquid biopsy 中，germline variants 的 VAF ~50%，而 somatic variants 的 VAF 通常 <10%。但在 solid tumor with high purity 中此區分變弱。

#### 2.4.2 Tumor-Only Calling 的傳統方法與局限

- 傳統方法依賴 **Panel of Normals (PON)** + 人群資料庫（gnomAD, dbSNP）過濾 germline (Chen et al., ClairS-TO, Nature Communications, 2025)
- **Private germline variants**（不在任何資料庫中）是主要的 FP 來源 (Teer et al., BMC Medical Genomics, 2017)
- Ancestry 相關的 germline FP rate 差異顯著 (Teer et al., 2017)
- **UNMASC**: 使用 unmatched normal controls 的方法 (Becker et al., NAR Cancer, 2021)
- **目前無任何工具** 直接利用甲基化作為 tumor-only pipeline 中的 germline filter

#### 2.4.3 ONT Long-read 的獨特優勢

- **同時偵測 base + methylation**: ONT 在定序時天然偵測 5mC/5hmC，無需額外處理 (Nanda et al., Epigenetics & Chromatin, 2024)
- **Read-level phasing + methylation**: 同一條 read 上可同時看到 variant + CpG methylation pattern
- **Long-range 資訊**: ONT reads (10-100kb) 可以涵蓋 variant 周圍大量 CpG 位點
- **Cell-of-origin 資訊**: Read-level methylation pattern 攜帶 cell-of-origin 資訊 (Baker et al., 2026; Jeong et al., 2025)

---

### 2.5 衝突觀點與不確定性

| 觀點 A（支持可行性） | 觀點 B（質疑可行性） | 可能原因 |
|---------------------|---------------------|---------|
| mQTL 普遍存在（4.7M cis-mQTL），germline 對甲基化有強 cis 效應 | Cancer 中甲基化 landscape 劇烈改變（全域低甲基化），可能破壞 mQTL signal | Cancer 的表觀遺傳重編程可能覆蓋 germline 的 cis 效應 |
| ASM 在 germline het SNPs 中很常見（23-37%） | Cancer 中 ASM 增加 5-9x，增加的 ASM 不一定與 germline 相關 | Cancer-specific ASM 增加的是 stochastic 的 allele-specific loss，而非 sequence-dependent ASM |
| ROCIT 證明 read-level methylation 可區分 tumor/normal reads | ROCIT 分類的是 cell origin，不是 variant origin | Tumor cell 可攜帶 germline variant，normal cell 不攜帶 somatic variant，但兩者的甲基化 "context" 不同 |
| CpG-SNP 機制清晰（38-88% ASM） | 只有 ~12% 的 mQTL SNP 直接在 CpG 位點 | 大多數 mQTL 通過 TF binding 間接影響甲基化，signal 可能較弱 |
| ONT 同時偵測 variant + methylation，read-level 整合理論上可行 | 實際操作中 ONT methylation calling accuracy ~95%，加上低 CpG density 區域信息量不足 | 需要足夠的 CpG sites per read 才能有統計力 |

---

## 3. 資料來源評估

| 來源 | 類型 | 可信度 | 年份 | 備註 |
|------|------|--------|------|------|
| Oliva et al., Nature Genetics | 論文 | 高 | 2023 | mQTL across diverse tissues |
| Shoemaker et al., Genome Research | 論文 | 高 | 2010 | ASM 與 CpG-SNP 的奠基研究 |
| Do et al., Genome Biology | 論文 | 高 | 2020 | ASM 在 cancer vs normal 的最大規模比較 |
| Baker et al., bioRxiv (ROCIT) | 預印本 | 中高 | 2026 | 最新最相關，但尚未 peer review |
| Jeong et al., Nature Communications (MethylBERT) | 論文 | 高 | 2025 | Read-level methylation classification |
| Akbari et al., Genome Biology (NanoMethPhase) | 論文 | 高 | 2021 | ONT ASM detection |
| MethPhaser, Nature Communications | 論文 | 高 | 2024 | Methylation-based phasing |
| Katsman et al., Genome Biology | 論文 | 高 | 2022 | ONT cfDNA methylation |
| Tarazona et al., Nature Genetics | 論文 | 高 | 2025 | Methylation + genomic alteration in NSCLC |
| Cytosine Methylation Mutability, Genetics | 論文 | 高 | 2020 | Methylation 影響 neighborhood mutability |
| Johnstone et al., Nature Communications | 論文 | 高 | 2019 | PMD hypomethylation in cancer |
| Teer et al., BMC Medical Genomics | 論文 | 中高 | 2017 | Tumor-only germline FP |
| Fu et al., Nature Reviews Genetics | 綜述 | 高 | 2025 | Long-read methylation analysis review |
| Nanda et al., Epigenetics & Chromatin | 綜述 | 高 | 2024 | ONT methylation clinical review |
| O'Neill et al., Cell Genomics | 論文 | 高 | 2024 | ONT advanced cancer cohort |

---

## 4. 對 InterSubMod 的具體啟示

### 4.1 可行性評估: 利用甲基化區分 Germline/Somatic

**結論: 理論上有生物學基礎，但直接實現面臨重大挑戰。**

**有利因素:**
1. ONT 同時提供 variant + methylation 的 read-level 資訊，是最佳平台
2. mQTL/ASM 的生物學基礎強大，germline variants 確實影響周圍甲基化
3. ROCIT/MethylBERT 證明 read-level methylation patterns 有足夠的信息量
4. InterSubMod 已有 read-level 甲基化解析能力

**不利因素:**
1. **目前沒有任何工具直接利用甲基化區分 germline/somatic variant**
2. Cancer 中甲基化 landscape 劇烈改變，mQTL signal 可能被淹沒
3. Signal 方向不明確: germline variant 的 ASM 在 cancer 中被放大（5-9x），不是減弱
4. CpG density 分佈不均，很多 variant 周圍 CpG 不足
5. 需要大量 training data（已知 germline/somatic label 的 variants with methylation）

### 4.2 可能的實現路徑

**路徑 A: ASM-based Filter（較簡單，較弱）**
- 原理: Germline het variants 更可能表現 stable ASM（兩個 haplotype 甲基化不同但穩定），somatic variants 更可能在 cancer 中表現 disorganized methylation
- 方法: 計算 variant 周圍 haplotype-resolved methylation 的 concordance/variance
- 缺點: 癌症本身增加 ASM，區分力不確定

**路徑 B: Read-level Methylation Context Classification（中等難度）**
- 原理: 仿 ROCIT 思路，利用 variant-bearing reads 的甲基化 context 分類
- 方法: 取每個 variant 的 supporting reads，提取周圍 CpG methylation pattern，用 ML 分類
- 需要: 已知 germline/somatic labels 的 training set（可從 paired-mode 樣本取得）
- 優勢: 直接與 InterSubMod 現有框架整合

**路徑 C: Population-level Methylation Signature（較複雜）**
- 原理: Germline variants 在人群中的 mQTL 效應是已知的，somatic variants 沒有
- 方法: 建立 mQTL database，查詢 variant 是否在已知 mQTL 中，觀察甲基化是否符合預期
- 缺點: 需要外部 mQTL database，且只適用於常見 germline variants

**路徑 D: 甲基化 Heterogeneity Score（最簡單）**
- 原理: Germline variant 的 reads 甲基化模式一致（反映正常 mQTL），somatic variant 的 reads 甲基化更異質（反映 tumor heterogeneity）
- 方法: 計算 ISM 已有的甲基化距離指標，比較 variant-bearing reads 的 within-group methylation variance
- 優勢: 不需要額外 training，直接利用 InterSubMod 的距離度量
- **推薦作為第一個探索方向**

### 4.3 與現有研究主線的關係

根據 `docs/experiments/INDEX.md` 和 `CURRENT_FOCUS.md`:
- TO 場景下所有特徵 AUC < 0.58，甲基化方向在 paired/TO 中有 5/9 反轉
- **甲基化 heterogeneity 作為 germline filter** 是全新方向，與現有 QS 改進正交
- 與 "E: Read-level epigenetic context" 研究方向（2026-Q2 策略）直接相關
- ROCIT 的出現（2026-03-03）是重要的外部驗證，證明 read-level methylation classification 是可行的

---

## 5. 建議行動

### 短期（P0, 1-2 週）
1. **探索性分析**: 取已知 germline/somatic variants（從 paired-mode 樣本已有 truth labels），比較兩組 variant 周圍的甲基化 heterogeneity（using ISM distance metrics）
2. **可視化**: 繪製 germline vs somatic variant 周圍的 per-read methylation heatmap，直觀觀察差異

### 中期（P1, 2-4 週）
3. **Feature engineering**: 設計甲基化相關 features（ASM score, methylation variance, haplotype methylation concordance）
4. **Pilot classifier**: 在 HCC1395 paired data 上訓練 germline/somatic classifier，評估 AUC
5. **CpG density stratification**: 分析 CpG density 對分類效果的影響

### 長期（P2, 1-2 月）
6. **ROCIT replication/adaptation**: 評估 ROCIT 在我們的 ONT data 上的表現
7. **整合進 InterSubMod pipeline**: 若效果良好，將甲基化 germline filter 整合為 QS 的一個新 component
8. **Cross-sample validation**: 在 COLO829, H1437, H2009 上驗證泛化性

---

## 6. 完整參考文獻列表

### 核心文獻

1. **Baker TM et al.** (2026). Genome-wide classification of tumor-derived reads from bulk long-read sequencing. *bioRxiv*. DOI: 10.64898/2026.03.03.709085
   - https://www.biorxiv.org/content/10.64898/2026.03.03.709085v1

2. **Jeong Y et al.** (2025). MethylBERT enables read-level DNA methylation pattern identification and tumour deconvolution using a Transformer-based model. *Nature Communications* 16:849.
   - https://www.nature.com/articles/s41467-025-55920-z

3. **Do C et al.** (2020). Allele-specific DNA methylation is increased in cancers and its dense mapping in normal plus neoplastic cells increases the yield of disease-associated regulatory SNPs. *Genome Biology* 21:153.
   - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02059-3

4. **Shoemaker R et al.** (2010). Allele-specific methylation is prevalent and is contributed by CpG-SNPs in the human genome. *Genome Research* 20:883-889.
   - https://pmc.ncbi.nlm.nih.gov/articles/PMC2892089/

5. **Oliva M et al.** (2023). DNA methylation QTL mapping across diverse human tissues provides molecular links between genetic variation and complex traits. *Nature Genetics* 55:112-122.
   - https://www.nature.com/articles/s41588-022-01248-z

### 技術方法

6. **Akbari V et al.** (2021). Megabase-scale methylation phasing using nanopore long reads and NanoMethPhase. *Genome Biology* 22:68.
   - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02283-5

7. **MethPhaser** (2024). MethPhaser: methylation-based long-read haplotype phasing of human genomes. *Nature Communications* 15:5327.
   - https://www.nature.com/articles/s41467-024-49588-0

8. **Zheng X et al.** (2014). MethylPurify: tumor purity deconvolution and differential methylation detection from single tumor DNA methylomes. *Genome Biology* 15:419.
   - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0419-x

### Cancer Methylation 基礎

9. **Tarazona N et al.** (2025). DNA methylation cooperates with genomic alterations during non-small cell lung cancer evolution. *Nature Genetics* 57:2226-2237.
   - https://www.nature.com/articles/s41588-025-02307-x

10. **Johnstone SE et al.** (2019). Partially methylated domains are hypervariable in breast cancer and fuel widespread CpG island hypermethylation. *Nature Communications* 10:1749.
    - https://www.nature.com/articles/s41467-019-09828-0

11. **Wen B et al.** (2017). Significant associations between driver gene mutations and DNA methylation alterations across many cancer types. *PLoS Computational Biology* 13:e1005840.
    - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005840

12. **Haggerty C et al.** (2020). Cytosine Methylation Affects the Mutability of Neighboring Nucleotides in Germline and Soma. *Genetics* 214:809-821.
    - https://pmc.ncbi.nlm.nih.gov/articles/PMC7153944/

### Tumor-Only Calling

13. **Teer JK et al.** (2017). A method to reduce ancestry related germline false positives in tumor only somatic variant calling. *BMC Medical Genomics* 10:61.
    - https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-017-0296-8

14. **Becker T et al.** (2021). UNMASC: tumor-only variant calling with unmatched normal controls. *NAR Cancer* 3:zcab040.
    - https://academic.oup.com/narcancer/article/3/4/zcab040/6382329

### cfDNA / Liquid Biopsy

15. **Katsman E et al.** (2022). Detecting cell-of-origin and cancer-specific methylation features of cell-free DNA from Nanopore sequencing. *Genome Biology* 23:158.
    - https://link.springer.com/article/10.1186/s13059-022-02710-1

16. **Tao Y et al.** (2024). A cfDNA methylation-based tissue-of-origin classifier for cancers of unknown primary. *Nature Communications* 15:3404.
    - https://www.nature.com/articles/s41467-024-47195-7

### ONT Long-read Cancer 應用

17. **O'Neill K et al.** (2024). Long-read sequencing of an advanced cancer cohort resolves rearrangements, unravels haplotypes, and reveals methylation landscapes. *Cell Genomics* 4:100464.
    - https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00293-3

18. **Nanda AS et al.** (2024). Shedding light on DNA methylation and its clinical implications: the impact of long-read-based nanopore technology. *Epigenetics & Chromatin* 17:40.
    - https://pmc.ncbi.nlm.nih.gov/articles/PMC11684317/

### 綜述

19. **Fu Y, Timp W, Sedlazeck FJ** (2025). Computational analysis of DNA methylation from long-read sequencing. *Nature Reviews Genetics* 26:620-634.
    - https://www.nature.com/articles/s41576-025-00822-5

20. **Nishiyama A, Nakanishi M** (2021). Navigating the DNA methylation landscape of cancer. *Trends in Genetics* 37:1012-1027.
    - https://www.cell.com/trends/genetics/fulltext/S0168-9525(21)00130-X

### mQTL / Germline-Somatic Interaction

21. **Gaunt TR et al.** (2016). Systematic identification of genetic influences on methylation across the human life course. *Genome Biology* 17:61.
    - https://link.springer.com/article/10.1186/s13059-016-0926-z

22. **Vali-Pour M, Lehner B, Supek F** (2022). The impact of rare germline variants on human somatic mutation processes. *Nature Communications* 13:3724.
    - https://www.nature.com/articles/s41467-022-31483-1

23. **Gao T et al.** (2022). Elucidating the genetic architecture of DNA methylation to identify promising molecular mechanisms of disease. *Scientific Reports* 12:19564.
    - https://www.nature.com/articles/s41598-022-24100-0

### Epigenetic Tumor Heterogeneity

24. **Gallegos JE et al.** (2022). Epigenetic tumor heterogeneity in the era of single-cell profiling with nanopore sequencing. *Clinical Epigenetics* 14:107.
    - https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-022-01323-6

25. **Zhou W et al.** (2018). DNA methylation loss in late-replicating domains is linked to mitotic cell division. *Nature Genetics* 50:591-602.
    - https://www.nature.com/articles/s41588-018-0073-4
