# Allele-Specific Methylation (ASM) 與亞克隆甲基化偵測方法學文獻調查

<!--
建立時間: 2026-04-15 14:00
目標: 全面調查 ASM 偵測、read-level 甲基化分類、tumor-normal ASM 比較、亞克隆甲基化推斷、長讀長特有方法、統計模型等方法學，評估與 ISM 的相容性
處理範圍: 6 大方向 × 30+ 工具/方法，含核心原理、輸入/輸出、ISM 相容性
關聯檔案:
  - docs/references/20260401_ASM偵測方法與預期比例文獻調查_01.md
  - docs/references/20260409_beyond_auc_biology_literature_review.md
  - docs/references/20260414_LOH_subclone_AF_methylation_literature.md
  - docs/concepts/2026/04/20260409_競品與學界工具比較_01.md
-->

## 1. 搜尋概述

- **關鍵詞**: allele-specific methylation, ASM, DAMEfinder, CAMDAC, methclone, epipolymorphism, CHALM, epialleleR, MethPhaser, NanoMethPhase, pycoMeth, CPEL, DSS, beta-binomial, BPRMeth, scMET, Methrix, PRISM, MethylBERT, EVOFLUx, modkit, modbam2bed, Shannon entropy, methylation haplotype block, CpG methylation pattern, subclonal methylation, LOH methylation
- **搜尋時間**: 2026-04-15
- **資料來源數**: 40+ 篇論文、工具文件與資料庫
- **搜尋範圍**: PubMed、bioRxiv、GitHub、Google Scholar、Bioconductor

---

## 2. 方法學分類總覽

| 分類 | 工具/方法 | 數量 |
|------|----------|------|
| A. Per-CpG ASM 偵測工具 | DAMEfinder, pycoMeth, NanoMethPhase, modkit, DSS, methylKit, epialleleR, CanASM | 8 |
| B. Read-level 甲基化模式分類 | MethPhaser, methclone, CHALM, epialleleR, Methrix, Metheor, WSHPackage, epihet | 8 |
| C. Tumor vs Normal ASM 比較 | CAMDAC, Do et al. 2020, DSS/DMLtest, MethylBERT, CanASM | 5 |
| D. Subclone 甲基化偵測 | PRISM, EVOFLUx, methclone, MethylBERT, Shannon entropy/epipolymorphism | 5 |
| E. 長讀長特有方法 | MethPhaser, NanoMethPhase, pycoMeth, modkit, modbam2bed, NANOME | 6 |
| F. 統計方法 | CPEL/Ising model, DSS/beta-binomial, BPRMeth, scMET, beta-mixture models, Fisher exact test | 6 |

---

## 3. 詳細方法學整理

### 3.1 Per-CpG ASM 偵測工具

---

#### 3.1.1 DAMEfinder

- **論文**: De Waele L, Fourne L, Lent J, et al. (2020). "DAMEfinder: a method to detect differential allele-specific methylation." *Epigenetics & Chromatin*, 13:25. DOI: 10.1186/s13072-020-00346-8
- **核心原理**: 在每個樣本計算所有 CpG 位點或 CpG pairs 的 ASM score，然後量化兩組條件間 ASM 的變化。不需 SNP 資訊時可純用 reads 量化 ASM，此方法與依賴 SNP 的方法比較表現良好。可偵測 Differentially Allele-specific Methylated regions (DAMEs)。
- **輸入需求**: Bisulfite sequencing (WGBS/RRBS) BAM 檔；可選 SNP 資訊
- **輸出指標**: Per-CpG ASM score；per-region DAME status（比較兩組條件）
- **平台**: Short-read bisulfite 設計；R/Bioconductor
- **ISM 相容性**: **中**。概念可借鑒（per-CpG ASM score），但 DAMEfinder 本身為短讀長 bisulfite 設計。ISM 的 methylation.csv 矩陣可作為輸入，需自行實作 ASM score 計算邏輯。
- **可在 raw_matrix 上計算**: 部分可行。需將 ISM 的 reads x CpGs 矩陣依 allele 分組後計算 per-CpG methylation fraction 差異。

---

#### 3.1.2 pycoMeth

- **論文**: Snajder R, Leger A, Stegle O, Bonder MJ. (2023). "pycoMeth: a toolbox for differential methylation testing from Nanopore methylation calls." *Genome Biology*, 24:83. DOI: 10.1186/s13059-023-02917-w
- **核心原理**: ONT 原生甲基化分析工具箱。包含 MetH5 高效儲存格式、Bayesian segmentation (Meth_Seg) 產生共識分段、Fisher's exact test + IHW (Independent Hypothesis Weighting) 多重校正進行差異甲基化檢測。支援 haplotype-aware ASM testing，接受 phased BAM 以 read groups 標記 haplotypes。
- **輸入需求**: ONT BAM + 甲基化 calls (MetH5 格式)；phased BAM for ASM
- **輸出指標**: Per-segment differential methylation p-values；ASM segments；bedMethyl-like output
- **平台**: ONT 原生設計；Python
- **ISM 相容性**: **高**。pycoMeth 是 ONT 生態系工具，與 ISM 使用相同輸入（ONT haplotagged BAM）。可作為 ISM 下游的獨立 ASM 驗證管道。但 pycoMeth 是全基因組分析，不是 per-SNV window。
- **可在 raw_matrix 上計算**: 不直接。pycoMeth 需要自己的 MetH5 格式輸入，不能直接讀取 ISM 的 methylation.csv。但概念上 Fisher's exact test + IHW 可在 ISM 矩陣上實作。

---

#### 3.1.3 epialleleR

- **論文**: Kiselev V, Leite JV, et al. (2023). "epialleleR: an R/Bioconductor package for sensitive allele-specific methylation analysis in NGS data." *GigaScience*, 12:giad087. DOI: 10.1093/gigascience/giad087
- **核心原理**: 從 BAM 檔直接計算 cytosine 甲基化以及 hypermethylated epiallele 頻率 (Variant Epiallele Frequency, VEF)。VEF 定義為在指定基因組區域中攜帶異常甲基化模式的 reads 比例。可提取甲基化 patterns (epialleles)，並測試甲基化與 SNV 的關聯性。
- **輸入需求**: BAM 檔（bisulfite 或 ONT 均可）；BED 區域
- **輸出指標**: Per-region VEF（variant epiallele frequency）；per-CpG methylation rates；epiallele pattern frequencies
- **平台**: R/Bioconductor；支援多種 NGS 平台
- **ISM 相容性**: **高**。epialleleR 可直接處理 BAM 檔計算 VEF，這個概念可直接整合到 ISM。ISM 的 read x CpG 矩陣本身就是 epiallele 的完整表達，VEF 可在 ISM output 上計算。
- **可在 raw_matrix 上計算**: **是**。將 methylation.csv 中每行（read）視為一個 epiallele，統計各種 pattern 頻率即可得 VEF。

---

#### 3.1.4 modkit (ONT 官方)

- **論文**: Oxford Nanopore Technologies. (2024+). modkit: A bioinformatics tool for working with modified bases. GitHub: https://github.com/nanoporetech/modkit
- **核心原理**: ONT 官方修飾鹼基分析工具。`modkit pileup --phased` 可依 HP tag 分群產生 hp1.bedmethyl + hp2.bedmethyl；`modkit dmr pair` 可比較兩個 haplotype 的差異甲基化區域，使用 Likelihood ratio scoring + Cohen's h effect size + HMM segmentation。
- **輸入需求**: modBAM（帶 MM/ML tags 的 BAM）；phased/haplotagged BAM for ASM
- **輸出指標**: bedMethyl 格式（per-CpG 甲基化率）；DMR regions with p-values and effect sizes
- **平台**: ONT 原生；Rust
- **ISM 相容性**: **非常高**。modkit 直接處理 ISM 使用的相同 haplotagged BAM，可作為 ISM 的正交 ASM 驗證。`modkit dmr pair` 的 hp1 vs hp2 比較結果可與 ISM 的 PERMANOVA 結果交叉驗證。
- **可在 raw_matrix 上計算**: 不直接（需 BAM 輸入），但概念上 modkit 的 likelihood ratio 統計量可在 ISM 矩陣上實作。

---

#### 3.1.5 DSS / DMLtest

- **論文**: Feng H, Conneely KN, Wu H. (2014). "A Bayesian hierarchical model to detect differentially methylated loci." *Nucleic Acids Research*; Park Y, Wu H. (2016). "Differential methylation analysis for BS-seq data under general experimental design." *Bioinformatics*, 32(10):1446-1453. DOI: 10.1093/bioinformatics/btw026
- **核心原理**: Beta-binomial model + Wald test。在每個 CpG 位點使用 beta-binomial 分布建模 count data，考慮生物學變異性（overdispersion），用 shrinkage estimation 穩定化參數估計，再用 Wald test 進行差異甲基化檢測。
- **輸入需求**: Per-CpG count data（methylated/total counts per condition）；至少兩個條件
- **輸出指標**: Per-CpG p-value, FDR；differentially methylated loci (DML)；differentially methylated regions (DMR)
- **平台**: R/Bioconductor
- **ISM 相容性**: **高**。ISM 的 methylation.csv 可依 allele (ALT/REF) 分組後，轉換為 per-CpG methylated/unmethylated counts，直接作為 DSS 的輸入。需將 ISM 的機率值二值化（threshold 0.5 或使用 ML probability 直接）或保留 count 形式。
- **可在 raw_matrix 上計算**: **是**。raw_matrix 依 allele 分組 -> per-CpG (methylated, total) counts -> DMLtest。

---

#### 3.1.6 methylKit

- **論文**: Akalin A, Kormaksson M, Li S, et al. (2012). "methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles." *Genome Biology*, 13:R87.
- **核心原理**: Fisher's exact test 或 logistic regression 進行差異甲基化分析。可整合 DSS 作為後端統計引擎。
- **輸入需求**: Per-CpG coverage + methylation fraction（bisulfite 或其他平台）
- **輸出指標**: Per-CpG / per-region differential methylation；p-values, q-values, methylation difference
- **平台**: R/Bioconductor
- **ISM 相容性**: **中**。需要先將 ISM 的 read-level data 聚合為 per-CpG statistics。
- **可在 raw_matrix 上計算**: 間接。需先將 raw_matrix 依 allele 分組後聚合為 per-CpG 的 coverage 和 methylation 比例。

---

#### 3.1.7 CanASM (資料庫)

- **論文**: (2025). "CanASM: a comprehensive database for genome-wide allele-specific DNA methylation identification and annotation in cancer." *BMC Genomics*. DOI: 10.1186/s12864-025-11849-7
- **核心原理**: 整合 226 個 BS-seq 樣本、31 種癌症類型的 ASM 資料庫，提供多維度調控註解（增強子、super enhancers、ATAC-seq、3D 互作區域、GWAS 關聯）。含三個互動分析工具：TF 結合親和力變化、3D chromatin 互作預測 ASM 靶基因、cell marker 變異分析。
- **輸入需求**: 查詢用（非分析工具）
- **輸出指標**: ASM 區域註解、調控元件重疊、GWAS 關聯
- **平台**: Web 資料庫
- **ISM 相容性**: **中**。作為 ISM 發現的 ASM 位點的下游註解參考。可將 ISM 偵測到的 ASM-positive regions 與 CanASM 交叉比對。
- **可在 raw_matrix 上計算**: 不適用（查詢式資料庫）

---

### 3.2 Read-Level 甲基化模式分類

---

#### 3.2.1 MethPhaser

- **論文**: Caetano SS, Timp W, Treangen TJ, et al. (2024). "MethPhaser: methylation-based long-read haplotype phasing of human genomes." *Nature Communications*, 15:5271. DOI: 10.1038/s41467-024-49588-0
- **核心原理**: 利用 ONT 甲基化信號延伸 SNV-based phasing。三步驟：(1) 從 SNV-haplotagged reads 計算統計上顯著差異的甲基化 CpG 位置；(2) 迭代地用甲基化資訊 haplotag reads；(3) 用新 haplotagged reads 橋接斷開的 phase blocks。ONT R9/R10 data 可增加 phase N50 78%-151%。
- **輸入需求**: ONT BAM + phased VCF (WhatsHap/HapCut2)
- **輸出指標**: 延伸後的 phased VCF（更長 phase blocks）；甲基化 CpG 差異位置
- **平台**: ONT 長讀長；Python
- **ISM 相容性**: **高（上游工具）**。MethPhaser 可作為 ISM 前的 phasing 改進步驟。ISM 目前依賴 LongPhase-S/TO 的 phasing，MethPhaser 可進一步改善 phase block 連續性，減少 unphased regions。ISM 的 raw_matrix 數據可用於評估 MethPhaser 的 ASM CpG 位點選擇。
- **可在 raw_matrix 上計算**: 不直接。MethPhaser 是 phasing 工具非分析工具。但 ISM 的矩陣可用於驗證 MethPhaser 選擇的 ASM CpG 是否與 ISM PERMANOVA 結果一致。

---

#### 3.2.2 methclone

- **論文**: Li S, Garrett-Bakelman FE, Akalin A, et al. (2014). "Dynamic evolution of clonal epialleles revealed by methclone." *Genome Biology*, 15:472. DOI: 10.1186/s13059-014-0472-5
- **核心原理**: 在每個含 4 個相鄰 CpG 的位點 (locus) 定義 2^4 = 16 種可能的 epiallele 模式。計算 Delta entropy（兩個時間點或兩個條件間的 epiallele 組成 Shannon entropy 變化）來評估 clonal dynamics。Epiallele shift 定義為 entropy 顯著改變的位點。
- **輸入需求**: Bisulfite BAM（Bismark 輸出）
- **輸出指標**: Per-locus epiallele composition；Delta entropy；epiallele shift loci
- **平台**: C++；設計用於 short-read bisulfite
- **ISM 相容性**: **高概念相容**。methclone 的 4-CpG epiallele 概念可直接在 ISM 的 methylation.csv 上實現。每 4 個相鄰 CpG 列作為一個 locus，每條 read 的甲基化模式即為一個 epiallele。可按 allele/HP 分組比較 epiallele 組成。
- **可在 raw_matrix 上計算**: **是**。直接取 methylation.csv 中相鄰 4 個 CpG 列，將每條 read 的 4-bit 甲基化模式計數，即可得 epiallele 頻率分布與 Shannon entropy。

---

#### 3.2.3 CHALM (Cell Heterogeneity-Adjusted cLonal Methylation)

- **論文**: Xu J, Shi J, Cui X, et al. (2021). "Cellular Heterogeneity-Adjusted cLonal Methylation (CHALM) improves prediction of gene expression." *Nature Communications*, 12:400. DOI: 10.1038/s41467-020-20492-7
- **核心原理**: 將 read-level 甲基化重新定義為 clonal methylation status。一條 read 的甲基化狀態不再是多個 CpG 的平均值，而是二元分類 -- 只要有 >=1 個 mCpG 就視為「methylated read」。Promoter 甲基化量化為 methylated reads / total reads，此定義更好地預測基因表達。
- **輸入需求**: Bisulfite BAM；gene annotation (promoter regions)
- **輸出指標**: Per-gene CHALM methylation score（methylated reads fraction）
- **平台**: Python；short-read bisulfite 設計
- **ISM 相容性**: **高**。CHALM 概念可直接在 ISM 的 methylation.csv 上實現。每條 read 只要有 >=1 個甲基化 CpG 即標記為 methylated。可計算 per-region 的 CHALM score，並按 allele/HP 分組比較。
- **可在 raw_matrix 上計算**: **是**。在 methylation.csv 中，每行（read）只要有 >=1 個值 > 0.5（或使用連續機率值），即可計算 CHALM score。

---

#### 3.2.4 Metheor

- **論文**: Lee D, Kim S. (2023). "Metheor: Ultrafast DNA methylation heterogeneity calculation from bisulfite alignments." *PLOS Computational Biology*.
- **核心原理**: 超快速（Rust 實作）計算所有六種已知的 DNA 甲基化異質性指標：PDR (Proportion of Discordant Reads)、FDRP (Fraction of Discordant Read Pairs)、qFDRP、Epipolymorphism、Methylation Entropy (ME)、MHL (Methylation Haplotype Load)。計算速度較原始實作快 300 倍，記憶體使用低 60 倍。
- **輸入需求**: Bisulfite BAM
- **輸出指標**: Per-locus PDR, FDRP, qFDRP, Epipolymorphism, ME, MHL
- **平台**: Rust；short-read bisulfite 設計
- **ISM 相容性**: **高概念相容**。六種異質性指標中 PDR、Epipolymorphism、ME 均可在 ISM 矩陣上直接計算。ISM 已有類似概念（距離矩陣本質上量化了 read 間的甲基化不一致性）。
- **可在 raw_matrix 上計算**: **是**。PDR = 含不一致 CpG 的 reads 比例；Epipolymorphism = 1 - sum(p_i^2)；ME = -sum(p_i * log(p_i))。均可從 raw_matrix 的 reads 直接計算。

---

#### 3.2.5 WSHPackage / epihet

- **論文**: Scherer M, et al. (2020). "Quantitative comparison of within-sample heterogeneity scores for DNA methylation data." *Nucleic Acids Research*, 48(8):e46; Zheng H, et al. (2020). "epihet for intra-tumoral epigenetic heterogeneity analysis and visualization." *Scientific Reports*, 10:21462.
- **核心原理**: WSHPackage 整合六種 within-sample heterogeneity scores（FDRP, qFDRP, PDR, Epipolymorphism, Methylation Entropy, MHL）。epihet 提供 intra-tumoral epigenetic heterogeneity 的分析與視覺化，含 PDR、Epipolymorphism、Shannon entropy。
- **輸入需求**: Bisulfite BAM 或 methclone 輸出
- **輸出指標**: Per-region heterogeneity scores
- **平台**: R
- **ISM 相容性**: **高概念相容**（同 Metheor）
- **可在 raw_matrix 上計算**: **是**

---

#### 3.2.6 Methrix

- **論文**: Kloser D, et al. (2020). "Methrix: an R/Bioconductor package for systematic aggregation and analysis of bisulfite sequencing data." *Bioinformatics*, 36(22-23):5524-5525.
- **核心原理**: 高效 bedGraph 讀取與聚合，建立 CpG x Sample 甲基化矩陣。支援批次處理 >100 樣本。可轉換為 bsseq 物件與 DSS/dmrseq 整合。
- **輸入需求**: bedGraph 或 tab-separated text（各種 methylation caller 輸出）
- **輸出指標**: CpG x Sample 矩陣；QC 報告；PCA 分析
- **平台**: R/Bioconductor
- **ISM 相容性**: **低**。Methrix 是 bulk-level 矩陣工具，ISM 是 read-level。但可用於管理 ISM 產出的 per-CpG 聚合統計。
- **可在 raw_matrix 上計算**: 不直接。Methrix 處理的是 site-level 聚合資料，不是 read x CpG 矩陣。

---

### 3.3 Tumor vs Normal ASM 比較

---

#### 3.3.1 CAMDAC (Copy-number Aware Methylation Deconvolution Analysis of Cancers)

- **論文**: Chamberlain EC. (2020/2022). "Copy number-aware deconvolution of tumor-normal DNA methylation profiles." bioRxiv: 10.1101/2020.11.03.366252; 用於 epiTRACERx: Martens JHA, et al. (2025). "DNA methylation cooperates with genomic alterations during non-small cell lung cancer evolution." *Nature Genetics*. DOI: 10.1038/s41588-025-02307-x
- **核心原理**: 形式化甲基化率、copy number 與 tumor purity 的關係，從 bulk tumor + tissue-matched normal bisulfite 測序數據中提取 purified tumor methylome。核心公式考慮每個 allele 的 copy number 和 purity，反向推導 cancer-cell-specific methylation。支援 RRBS 和 WGBS。Deconvolved methylation rates 可進行無偏的 tumor-normal 和 tumor-tumor differential methylation calling，read-phasing 驗證 CAMDAC methylation rates 並直接連結 genotype 和 epitype。
- **輸入需求**: RRBS/WGBS BAM（tumor + matched normal）；SNP array 或 WGS（CN profile）
- **輸出指標**: Purified tumor methylation rates；allele-specific CN；tumor purity；inter-tumor methylation distance (ITMD)；DMRs
- **平台**: R
- **ISM 相容性**: **中-高（概念層面非常相關）**。CAMDAC 的 CN-aware purity deconvolution 概念與 ISM Phase 2 的 Normal Methylation Reference 方向直接相關。ISM 目前已有 Coverage_Multiple 作為 CN proxy (r=0.831)。CAMDAC 的 read-phasing 驗證方法論可借鑒。但 CAMDAC 是 bulk-level 工具，ISM 是 read-level -- ISM 不需做 deconvolution 因為已經有 per-read 資訊。
- **可在 raw_matrix 上計算**: 不直接。CAMDAC 需要 bulk-level pileup data + CN profile。ISM 已有 read-level 解析，概念上更精確。

---

#### 3.3.2 Do et al. 2020 方法

- **論文**: Do C, Duber CB, Tycko B, et al. (2020). "Allele-specific DNA methylation is increased in cancers and its dense mapping in normal plus neoplastic cells increases the yield of disease-associated regulatory SNPs." *Genome Biology*, 21:153. DOI: 10.1186/s13059-020-02059-3
- **核心原理**: 全基因組甲基化定序（WGBS）用於多種正常細胞和組織及三種癌症類型。核心發現：(1) Cancer ASM 是 normal 的 5-9 倍；(2) 49-76% 的 cancer-only ASM 由 allele-specific loss of methylation (LOM) 驅動；(3) 43% cancer allele switching vs 14% in normal；(4) 6-17% somatic mutations 伴隨 ASM gain；(5) 71% de novo ASM TF binding motifs 與 germline-driven ASM 相同。
- **輸入需求**: WGBS（tumor + matched normal）；SNP calling for allele assignment
- **輸出指標**: Per-region ASM status；allele switching rate；LOM/GOM 分類
- **平台**: Custom pipeline
- **ISM 相容性**: **高（理論基礎）**。Do et al. 的 cancer ASM 增幅數據 (5-9x) 與 ISM 觀察到的 FP (germline leak) vs TP (somatic) ASM 差異一致。ISM 的 PERMANOVA 結果可以與 Do et al. 的預期比例交叉驗證。
- **可在 raw_matrix 上計算**: 部分可行。ISM 的 read-level allele assignment + methylation.csv 可計算類似的 ASM prevalence 統計。但完整的 LOM/GOM 分類需要 normal sample 作為 baseline。

---

#### 3.3.3 MethylBERT

- **論文**: Lee D, Kim S, et al. (2025). "MethylBERT enables read-level DNA methylation pattern identification and tumour deconvolution using a Transformer-based model." *Nature Communications*, 16:788. DOI: 10.1038/s41467-025-55920-z
- **核心原理**: 使用 BERT Transformer 模型，在 read-level 識別 tumor-derived 和 normal-derived 序列讀段。Pre-trained BERT 編碼 read-level methylome（序列 + 甲基化模式），分類每條 read 為 tumor 或 normal cell type。用 posterior probability + Bayes' theorem + MLE 估計 tumor purity。
- **輸入需求**: WGBS/RRBS BAM；training data（tumor/normal 標籤）
- **輸出指標**: Per-read tumor/normal classification；tumor cell fraction estimation
- **平台**: Python (PyTorch)
- **ISM 相容性**: **中-高**。MethylBERT 的 read-level classification 概念與 ISM 高度相似，但目的不同（MethylBERT 做 deconvolution，ISM 做 variant characterization）。ISM 的 reads x CpGs 矩陣理論上可作為 MethylBERT 的輸入特徵。若 ISM 有 tumor/normal dual-BAM（Phase 2），可訓練類似分類器。
- **可在 raw_matrix 上計算**: 需重新訓練模型。raw_matrix 格式需轉換為 MethylBERT 的序列 + 甲基化 token 表示。

---

### 3.4 Subclone 甲基化偵測

---

#### 3.4.1 PRISM

- **論文**: Lee D, Lee S, Kim S. (2019). "PRISM: methylation pattern-based, reference-free inference of subclonal makeup." *Bioinformatics*, 35(14):i520-i529. DOI: 10.1093/bioinformatics/btz470
- **核心原理**: Reference-free 的亞克隆組成推斷。使用 DNMT1-like HMM 進行 in silico proofreading，校正 bisulfite 轉化錯誤。在短基因組區域中尋找 dichotomous patterns（可分為完全甲基化和完全非甲基化的模式），其兩種模式的頻率構成亞克隆豐度的充分統計量。用 beta-binomial mixture + EM 算法推斷亞克隆組成。
- **輸入需求**: BS-seq BAM
- **輸出指標**: Subclone proportions；evolutionary tree
- **平台**: Python
- **ISM 相容性**: **高**。PRISM 的 dichotomous pattern 偵測概念可直接在 ISM 的 methylation.csv 上實現。ISM 的 raw_matrix 本質上就是 PRISM 分析所需的 reads x CpGs 矩陣。ISM 已有的 hierarchical clustering 也在做類似的模式分群。但 PRISM 不整合 HP/allele 資訊。
- **可在 raw_matrix 上計算**: **是**。ISM 的 methylation.csv 可直接識別 dichotomous patterns（如果某 region 的 reads 可明確分為 fully methylated 和 fully unmethylated 兩組）。Beta-binomial mixture fitting 可在此基礎上計算。

---

#### 3.4.2 EVOFLUx

- **論文**: Gabbutt C, Duran-Ferrer M, et al. (2025). "Fluctuating DNA methylation tracks cancer evolution at clinical scale." *Nature*. DOI: 10.1038/s41586-025-09374-4
- **核心原理**: 利用天然 DNA 甲基化 barcodes（fluctuating CpGs, fCpGs -- 隨時間開關甲基化的位點）推斷腫瘤演化動力學。只需 bulk 甲基化 profile 即可推斷初始生長速率、惡性腫瘤年齡、epimutation rates。應用於 1,976 個淋巴癌樣本。使用群體遺傳學模型 (population genetics framework) + Maximum Likelihood Estimation。
- **輸入需求**: Bulk methylation array 或 RRBS/WGBS
- **輸出指標**: Tumor growth rate；malignancy age；epimutation rate；subclonal selection events
- **平台**: R
- **ISM 相容性**: **低-中**。EVOFLUx 是 bulk-level 工具，ISM 是 read-level per-SNV 分析。但 fCpG 概念有參考價值 -- ISM 視窗內的 CpG variance 可用類似框架分析。ISM 已在 L4 文獻驗證中測試 fCpG 概念（NEGATIVE 結果：TP/FP CpG variance 完全相同）。
- **可在 raw_matrix 上計算**: 概念上可行但不直接。需要先定義 ISM 視窗內的 fCpG 位點（high variance CpGs），再計算 pattern diversity。

---

#### 3.4.3 Methylation Haplotype Blocks (MHBs)

- **論文**: Guo S, Diep D, Plongthongkum N, et al. (2017). "Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples." *Nature Genetics*, 49:443-449; Lee J, et al. (2025). "Toward the DNA methylation haplotype map of 11 common solid cancers." *Cell Reports*.
- **核心原理**: 相鄰 CpG 位點形成 methylation haplotype blocks (MHBs)，在同一 block 內 CpG 的甲基化狀態高度相關。MHB 可用於組織 deconvolution 和 tumor tissue-of-origin mapping。2025 年研究在 110 個原發腫瘤中定義 81,567 MHBs，具高度癌症類型特異性。
- **輸入需求**: WGBS BAM
- **輸出指標**: MHB 定義（基因組座標）；per-MHB discordance score；cancer-type classification
- **平台**: Python/R
- **ISM 相容性**: **高概念相容**。ISM 的 methylation.csv 天然記錄了 CpG 間的 co-methylation 模式。可在 ISM 視窗內定義 local MHBs 並計算 per-allele discordance。
- **可在 raw_matrix 上計算**: **是**。raw_matrix 中相鄰 CpG 的 co-methylation 可直接計算（pairwise correlation），定義 local MHB 邊界。

---

### 3.5 長讀長 (ONT/PacBio) 特有方法

---

#### 3.5.1 NanoMethPhase

- **論文**: Akbari V, Garant JM, O'Neill K, et al. (2021). "Megabase-scale methylation phasing using nanopore long reads and NanoMethPhase." *Genome Biology*, 22:68. DOI: 10.1186/s13059-021-02283-5
- **核心原理**: 端到端 pipeline 整合 SNV calling + WhatsHap phasing + 甲基化 phasing。將 ONT reads phase 到 haplotype 後轉為 mock-WGBS 格式，再用 DSS R package 做 differential methylation analysis (DMA)。含 SNVoter 子工具改善低覆蓋度區域的 SNV calling 準確度。~10x genome-wide 覆蓋度即可偵測 ASM。
- **輸入需求**: ONT BAM + reference genome；可選 VCF
- **輸出指標**: Haplotype-resolved methylation calls (mock-WGBS format)；DMRs between haplotypes
- **平台**: ONT；Python
- **ISM 相容性**: **中-高**。NanoMethPhase 的 haplotype-resolved 甲基化概念與 ISM 高度相似（都是依 HP 分群甲基化）。但 NanoMethPhase 是全基因組分析，ISM 是 per-SNV window。NanoMethPhase 的 DMA 結果可作為 ISM ASM 分析的正交驗證。
- **可在 raw_matrix 上計算**: 不直接。NanoMethPhase 需要 raw BAM 輸入。但其 DSS-based DMA 步驟可在 ISM 矩陣上實現。

---

#### 3.5.2 NANOME

- **論文**: NANOME consortium. (2025). "NANOME: XGBoost consensus methylation + Clair3 variant calling." *Genome Biology*.
- **核心原理**: 整合 Megalodon、Nanopolish、DeepSignal 三個甲基化 caller 的 XGBoost 共識模型。結合 variant calling (Clair3) 和 long-read phasing 做 haplotype-aware 分析。單鹼基解析度 MSE 改善 +12%，單分子 F1 改善 +3%。
- **輸入需求**: ONT raw signal data (FAST5/POD5)
- **輸出指標**: Consensus methylation calls；haplotype-aware methylation profiles
- **平台**: ONT；Nextflow pipeline
- **ISM 相容性**: **低**。NANOME 是上游 methylation calling 工具。ISM 已使用 Dorado 的 methylation calls (MM/ML tags)，不需要 NANOME 的 consensus 方法。
- **可在 raw_matrix 上計算**: 不適用

---

#### 3.5.3 modbam2bed

- **論文**: Oxford Nanopore Technologies. modbam2bed. GitHub: https://github.com/epi2me-labs/modbam2bed
- **核心原理**: 將 modBAM 中的修飾鹼基資訊轉換為 BED 格式。計算 per-CpG 的 coverage 和 methylation 百分比。可搭配下游分析工具（DSS, methylKit）使用。
- **輸入需求**: modBAM（帶 MM/ML tags）
- **輸出指標**: BED format（per-CpG methylation frequencies）
- **平台**: ONT；C
- **ISM 相容性**: **中**。作為 ISM 的輔助轉換工具，可將 ISM 使用的 BAM 轉為 BED 格式供其他工具使用。但 ISM 已自行從 BAM 提取 per-read 甲基化資訊。
- **可在 raw_matrix 上計算**: 不適用

---

### 3.6 統計方法

---

#### 3.6.1 CPEL / Ising Model

- **論文**: Jenkinson G, Pujadas E, Goutsias J, Feinberg AP. (2018). "An information-theoretic approach to the modeling and analysis of whole-genome bisulfite sequencing data." *BMC Bioinformatics*, 19:129; Jenkinson G, et al. (2020). "Detection of haplotype-dependent allele-specific DNA methylation in WGBS data." *Nature Communications*, 11:5238. DOI: 10.1038/s41467-020-19077-1
- **核心原理**: 使用統計物理學的一維 Ising 模型聯合建模相鄰 CpG 位點的甲基化狀態，捕捉 CpG 間的空間相關性。三個檢驗統計量：T_MML (mean methylation level 差異)、T_NME (normalized methylation entropy 差異)、T_PDM (遺傳-表觀遺傳關聯)。Bootstrap hypothesis testing + BH correction。關鍵發現：**96% 的 ASM 是 entropy imbalance（mean 無差但 pattern 不同）**，只有 CPEL 能偵測此類 ASM。
- **輸入需求**: WGBS BAM + phased VCF (haplotype assignment)
- **輸出指標**: Per-haplotype-block MML, NME, PDM statistics；ASM classification (MML-hap, NME-hap)
- **平台**: MATLAB/Python
- **ISM 相容性**: **非常高（理論意義最大）**。CPEL 的 entropy imbalance 概念與 ISM 已有的 PERMANOVA 高度互補。ISM 的 PERMANOVA 偵測 overall pattern difference，CPEL 的 NME 可進一步拆解為 mean difference vs entropy difference。96% ASM 是 entropy imbalance 這一發現暗示 ISM 的 PERMANOVA（對 pattern variance 敏感）比 per-CpG mean comparison 更有優勢。
- **可在 raw_matrix 上計算**: **是，且非常直接**。ISM 的 methylation.csv 依 allele 分組後，每組的 reads 定義一個甲基化模式分布。NME = normalized entropy of pattern distribution。可直接計算兩組的 NME 差異和 MML 差異。Ising 模型的 coupling parameters 可用 MLE 從 read patterns 估計。
- **計算成本注意**: 原始 CPEL 極耗時（~48h per sample, 20 CPUs）。但在 ISM 的小視窗 (5-50 CpGs) 內，計算量大幅降低，預計可在 per-region 秒級完成。

---

#### 3.6.2 Beta-Binomial Model (DSS 框架)

- **論文**: Wu H, Xu T, Feng H, et al. (2015). "Detection of differentially methylated regions from whole-genome bisulfite sequencing data without replicates." *Nucleic Acids Research*; Park Y, Wu H. (2016). *Bioinformatics*, 32(10):1446.
- **核心原理**: 在每個 CpG 位點，用 beta-binomial 分布建模 count data (methylated reads / total reads)。Beta-binomial 同時捕捉 binomial sampling noise 和 biological overdispersion (between-sample variability)。使用 shrinkage estimation 穩定化散度參數，再進行 Wald test。
- **輸入需求**: Per-CpG (methylated count, total count) per condition
- **輸出指標**: Per-CpG test statistic, p-value, FDR；DML/DMR
- **平台**: R/Bioconductor (DSS package)
- **ISM 相容性**: **高**。ISM 的 methylation.csv 依 allele 分組後，per-CpG 聚合為 (methylated count, total count)，直接套入 DSS 框架。
- **可在 raw_matrix 上計算**: **是**。raw_matrix -> allele split -> per-CpG aggregation -> DSS DMLtest。

---

#### 3.6.3 BPRMeth (Bayesian Probit Regression for Methylation)

- **論文**: Kapourani CA, Sanguinetti G. (2018). "BPRMeth: a flexible Bioconductor package for modelling methylation profiles." *Bioinformatics*, 34(14):2485-2486. DOI: 10.1093/bioinformatics/bty129
- **核心原理**: 用 Binomial Probit Regression 建模甲基化 profiles，提取高階特徵。Variational inference 提供 Bayesian posterior confidence。可用於 methylation profile 的 clustering 和基因表達預測。10 kb 視窗、5000 promoters 約 5 分鐘推論、20 分鐘 clustering。
- **輸入需求**: Per-CpG methylation 值 + 基因組位置；可用於 single cell 和 array
- **輸出指標**: Smoothed methylation profiles；profile-based clustering；gene expression prediction
- **平台**: R/Bioconductor
- **ISM 相容性**: **中**。BPRMeth 的 profile modeling 可用於 ISM 視窗內的甲基化 profile smoothing 和 clustering。但 ISM 已有 distance-based clustering，BPRMeth 提供的是 model-based 替代方案。
- **可在 raw_matrix 上計算**: 間接。需先將 raw_matrix 轉為 per-CpG aggregated methylation 值（by allele 或 by cluster），再進行 profile fitting。

---

#### 3.6.4 scMET (Bayesian Methylation Heterogeneity)

- **論文**: Kapourani CA, Argelaguet R, Sanguinetti G, Vallejos CA. (2021). "scMET: Bayesian modeling of DNA methylation heterogeneity at single-cell resolution." *Genome Biology*, 22:114. DOI: 10.1186/s13059-021-02329-8
- **核心原理**: 階層式 beta-binomial (BB) 模型 + Generalized Linear Model (GLM) 框架。Feature-specific 參數 mu_j 量化 overall DNAm，overdispersion 參數 gamma_j 作為 cell-to-cell 甲基化異質性的代理。可偵測 mean methylation 差異和 variability 差異。用於高度變異特徵 (HVF) 鑑定。
- **輸入需求**: scBS-seq data（single-cell level methylation counts）
- **輸出指標**: Feature-specific mean (mu) and variability (gamma)；HVF identification；differential mean and variability testing
- **平台**: R/Bioconductor
- **ISM 相容性**: **中-高**。scMET 的 beta-binomial 框架可將每條 ONT read 視為一個 "cell"，ISM 的 per-read methylation 天然對應 scMET 的 single-cell data。gamma 參數（variability）概念與 ISM 的 distance matrix 異質性高度相關。特別是 **differential variability testing** 可用於比較 TP vs FP 的甲基化異質性差異。
- **可在 raw_matrix 上計算**: **是，概念上完全對應**。ISM 的每條 read = scMET 的一個 cell；每個 CpG = 一個 feature。raw_matrix 的 binary/probability 值直接對應 scMET 的 methylation count data。需要 R 環境執行。

---

#### 3.6.5 Beta-Mixture Models

- **論文**: Paulo C, et al. (2024). "A novel family of beta mixture models for the differential analysis of DNA methylation data." *PLOS ONE*. DOI: 10.1371/journal.pone.0314014
- **核心原理**: 用 beta distribution mixture 建模甲基化值的多模態分布。Objectively infer methylation state thresholds（不需預設 cutoff），用 model-based clustering 鑑定 differentially methylated CpG sites。
- **輸入需求**: Per-CpG methylation beta values（array 或 sequencing aggregated）
- **輸出指標**: Mixture components；state classification；differential methylation
- **平台**: R
- **ISM 相容性**: **中**。可用於 ISM 視窗內的甲基化分布建模（e.g., 是否為 bimodal 分布指示 ASM）。
- **可在 raw_matrix 上計算**: 間接。需先 per-CpG 聚合為 beta value 分布。

---

#### 3.6.6 Fisher's Exact Test (Per-CpG)

- **論文**: 廣泛使用於多個工具中（DAMEfinder, pycoMeth, methylKit, NanoMethPhase 等）
- **核心原理**: 對每個 CpG 位點建立 2x2 列聯表（allele [ALT/REF] x methylation [methylated/unmethylated]），執行 Fisher's exact test。多重檢驗校正使用 BH-FDR。定義 ASM region: >= 1 significant CpG with |delta_methylation| >= threshold。
- **輸入需求**: Per-CpG per-allele methylated/unmethylated counts
- **輸出指標**: Per-CpG p-value, FDR, odds ratio, delta methylation
- **平台**: 任何統計軟體
- **ISM 相容性**: **非常高且最容易實作**。ISM 的 methylation.csv 依 allele 分組，每個 CpG 直接建立 2x2 表。已在 ISM 的 ASM 分析策略中被推薦為 Phase 1 方法之一。
- **可在 raw_matrix 上計算**: **是，最直接**。raw_matrix + allele label -> per-CpG 2x2 table -> Fisher's exact test -> BH-FDR correction。

---

### 3.7 其他重要參考資源

---

#### 3.7.1 Loyfer et al. 2025 -- ASM Atlas

- **論文**: Rosenski J, Peretz A, Magenheim J, Loyfer N, et al. (2025). "Atlas of imprinted and allele-specific DNA methylation in the human body." *Nature Communications*, 16:2415. DOI: 10.1038/s41467-025-57433-1
- **核心發現**: 39 種正常人類細胞類型的 deep WGS 定義 325K bimodal methylation 區域（6% genome, 11% CpGs）；其中 34K 區域由遺傳變異驅動 ASM；460 區域為 parental imprinting。ASM 經常是 cell-type-specific。
- **ISM 相關性**: **非常高**。提供 ISM 視窗內 CpG 的 expected ASM baseline。可將 ISM 偵測到的 ASM 位點與 Atlas 交叉比對，區分 germline-driven vs tumor-specific ASM。

---

#### 3.7.2 epiTRACERx / CAMDAC in NSCLC

- **論文**: Martens JHA, et al. (2025). "DNA methylation cooperates with genomic alterations during non-small cell lung cancer evolution." *Nature Genetics*. DOI: 10.1038/s41588-025-02307-x
- **核心發現**: 使用 CAMDAC 在 421 NSCLC regions 取得 cancer-cell-specific methylation。發現 epigenomic ITH 在 chromatin-accessible regions 被 specifically reduced（暗示 immune surveillance），全域低甲基化是 subclone-independent（trunk event）。
- **ISM 相關性**: **高**。ISM Phase 2 的 tumor-normal ASM 比較可參考此方法論。

---

#### 3.7.3 DNA Methylation Haplotype Map of Solid Cancers (2025)

- **論文**: Lee J, et al. (2025). "Toward the DNA methylation haplotype map of 11 common solid cancers." *Cell Reports*.
- **核心發現**: 110 原發腫瘤 x 11 癌種定義 81,567 MHBs。MHB discordance 與 driver mutations 和 inflammatory pathways 關聯。co-methylation patterns 優於 mean methylation 做 cancer detection。
- **ISM 相關性**: **高**。ISM 的 read-level 矩陣天然記錄 co-methylation（同一 read 上多個 CpG 的聯合狀態），可直接定義 window-internal MHBs。

---

## 4. ISM 相容性評估矩陣

### 4.1 可直接在 ISM raw_matrix (reads x CpGs) 上計算的方法

| 方法 | 實作難度 | 新增 insight | 優先級 | 說明 |
|------|---------|------------|--------|------|
| **Fisher's exact test per-CpG** | 低 | 中-高 | **P0** | 最直接；依 allele 分組 -> 2x2 table -> p-value |
| **Epipolymorphism** | 低 | 中 | P1 | 4-CpG epiallele pattern frequencies -> 1-sum(p_i^2) |
| **Shannon entropy (ME)** | 低 | 中 | P1 | 4-CpG patterns -> -sum(p_i * log(p_i)) |
| **PDR** | 低 | 中 | P1 | discordant reads / total reads |
| **CHALM score** | 低 | 低-中 | P2 | reads with >=1 mCpG / total reads |
| **MHL** | 中 | 中 | P2 | weighted methylation haplotype load |
| **NME (CPEL entropy)** | 中 | **高** | **P0** | normalized entropy per-allele -> entropy imbalance |
| **methclone epiallele shift** | 中 | 中 | P1 | Delta entropy between two conditions |
| **DSS beta-binomial per-CpG** | 中 | 中-高 | P1 | allele-split -> per-CpG counts -> Wald test |
| **VEF (epialleleR)** | 低 | 中 | P1 | variant epiallele frequency |
| **PRISM dichotomous patterns** | 中 | 中 | P2 | fully meth vs fully unmeth read fraction |
| **scMET variability (gamma)** | 中-高 | **高** | P1 | read-as-cell beta-binomial heterogeneity |
| **Local MHB definition** | 中 | 中 | P2 | CpG pairwise correlation -> block boundaries |

### 4.2 需要額外輸入或獨立管道的方法

| 方法 | 額外需求 | ISM 整合方式 | 優先級 |
|------|---------|-------------|--------|
| **modkit dmr pair** | Haplotagged BAM | 正交 ASM 驗證 | P1 |
| **pycoMeth Meth_Comp** | MetH5 format | 全基因組 ASM 比較 | P2 |
| **NanoMethPhase DMA** | Raw ONT BAM | Haplotype-resolved DMA | P2 |
| **CAMDAC deconvolution** | Normal WGBS + CN profile | Tumor-specific methylation | P2 |
| **MethylBERT** | Training data + GPU | Read-level tumor/normal | P3 |
| **MethPhaser** | Phased VCF | 上游 phasing 改進 | P2 |
| **EVOFLUx** | Bulk methylation array | 演化推斷（不同框架） | P3 |
| **BPRMeth** | Aggregated profiles | Profile-based clustering | P3 |

---

## 5. ISM 架構優勢與獨特定位

基於本次文獻調查，ISM 的核心架構具有以下獨特優勢：

### 5.1 ISM 已有且與文獻方法對齊的能力

| ISM 能力 | 對應文獻方法 | 優勢說明 |
|----------|-------------|---------|
| Per-read methylation matrix | methclone, PRISM, epialleleR | ISM 已有完整 reads x CpGs 矩陣，是最多方法的共同基礎 |
| PERMANOVA on distance matrix | 無直接對應 | **ISM 獨有**：距離矩陣 + permutation test 捕捉 overall pattern 差異 |
| HP tag 四群組分類 | NanoMethPhase (2 group), pycoMeth (2 group) | ISM 的 HP1/HP2 x ALT/REF 四群組分類是文獻中最細緻的 |
| LOH-aware 分析 | CAMDAC (CN-aware), Do et al. | ISM 的 HP ratio LOH detection + LOH-aware QS 是長讀長工具中獨有的 |
| 六種距離度量 | 無直接對應 | 學界工具通常只用 1-2 種統計量 |
| Somatic VCF anchoring | 無直接對應 | **ISM 獨有**：以 somatic variant 為錨點分析甲基化，所有其他工具都是 genome-wide 或 region-based |

### 5.2 ISM 目前缺少且文獻中重要的能力

| 缺少能力 | 對應文獻方法 | 實作難度 | Phase 2 相關性 |
|----------|-------------|---------|---------------|
| Per-CpG ASM 統計 | Fisher's exact test, DSS | **低** | **直接相關** -- Phase 2A Sample ASM |
| Entropy imbalance 量化 | CPEL NME | **中** | 高 -- 96% ASM 是 entropy 差異 |
| Epiallele pattern 計數 | methclone, Metheor | **低** | 中 -- 異質性量化 |
| Normal baseline 比較 | CAMDAC, Do et al. | **中** (Phase 2 已規劃) | **直接相關** -- Phase 2A 核心 |
| CpG 間相關性建模 | CPEL Ising, BPRMeth | **中-高** | 中 -- 改善統計 power |
| Subclone proportion 估計 | PRISM, beta-binomial mixture | **中** | 高 -- Phase 2D 目標 |

---

## 6. 建議行動

### 6.1 短期可實作（在 ISM 現有 Python 分析腳本中加入）

1. **Per-CpG Fisher's exact test + BH-FDR** -- 從 methylation.csv + reads.tsv 的 allele label 計算，標記每個 region 的 ASM-positive CpG 數量和位置。
2. **Epipolymorphism + Shannon Entropy** -- 4-CpG sliding window 計算 pattern diversity，可按 allele/HP 分組比較。
3. **NME (Normalized Methylation Entropy) per-allele** -- CPEL 概念的簡化版，計算每個 allele group 的 pattern entropy，比較 entropy imbalance。
4. **VEF (Variant Epiallele Frequency)** -- 計算 hypermethylated/hypomethylated epiallele 的比例。

### 6.2 中期整合（Phase 2 框架內）

1. **Normal Baseline ASM** -- 使用 Phase 2A 的 Normal BAM integration，計算 tumor vs normal 的 per-CpG ASM 變化（CAMDAC/Do et al. 概念）。
2. **modkit dmr pair 正交驗證** -- 在 haplotagged BAM 上執行 modkit dmr，與 ISM 內部 ASM 結果交叉比較。
3. **scMET-inspired variability testing** -- 將 reads 視為 cells，用 beta-binomial 框架量化 per-region 甲基化異質性。
4. **Local MHB 定義** -- 利用 ISM 的 read-level co-methylation 資訊定義 window 內的 methylation haplotype blocks。

### 6.3 長期研究方向

1. **CPEL Ising Model 的輕量化實作** -- 在 ISM 的小視窗 (5-50 CpGs) 內，Ising 模型計算量可控，可作為 ISM 的新統計模組。
2. **MethylBERT-inspired read classifier** -- 若 Phase 2 dual-BAM 提供足夠 training data，可訓練 read-level tumor/normal 分類器。
3. **PRISM-like subclone proportion estimation** -- 從 ISM 的 per-region methylation patterns 推斷 subclone 組成。

---

## 7. 資料來源與可信度評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| CPEL - Jenkinson et al. 2020, Nature Communications | 同行評審論文 | 高 | 96% ASM 是 entropy imbalance 的關鍵發現 |
| DAMEfinder - De Waele et al. 2020, Epigenetics & Chromatin | 同行評審論文 | 高 | Short-read 設計但概念可遷移 |
| pycoMeth - Snajder et al. 2023, Genome Biology | 同行評審論文 | 高 | ONT 原生工具，最相關 |
| CAMDAC - 2020 bioRxiv + 2025 Nature Genetics (epiTRACERx) | Preprint + 頂級期刊應用 | 高 | CN-aware deconvolution 業界標準 |
| Do et al. 2020, Genome Biology | 同行評審論文 | 高 | Cancer ASM 5-9x 增幅的核心參考 |
| methclone - Li et al. 2014, Genome Biology | 同行評審論文 | 高 | Epiallele 分析經典方法 |
| CHALM - Xu et al. 2021, Nature Communications | 同行評審論文 | 高 | Read-level clonal methylation 定義 |
| epialleleR - Kiselev et al. 2023, GigaScience | 同行評審論文 | 高 | VEF 概念直接可用 |
| MethPhaser - 2024, Nature Communications | 同行評審論文 | 高 | 甲基化 phasing 最新方法 |
| NanoMethPhase - Akbari et al. 2021, Genome Biology | 同行評審論文 | 高 | ONT haplotype-resolved methylation |
| PRISM - Lee et al. 2019, Bioinformatics | 同行評審論文 | 高 | Reference-free subclone inference |
| MethylBERT - Lee et al. 2025, Nature Communications | 同行評審論文 | 高 | Transformer read-level 分類最新方法 |
| EVOFLUx - Gabbutt et al. 2025, Nature | 同行評審論文 | 高 | 甲基化演化推斷最新方法 |
| modkit - ONT official | 工具文件 | 中-高 | 官方工具，持續更新 |
| DSS - Park & Wu 2016, Bioinformatics | 同行評審論文 | 高 | Beta-binomial DML 業界標準 |
| scMET - Kapourani et al. 2021, Genome Biology | 同行評審論文 | 高 | Bayesian heterogeneity 建模 |
| BPRMeth - Kapourani & Sanguinetti 2018, Bioinformatics | 同行評審論文 | 高 | Bayesian profile modeling |
| Metheor - Lee & Kim 2023, PLOS Comp Bio | 同行評審論文 | 高 | 超快速異質性計算 |
| Loyfer/Rosenski et al. 2025, Nature Communications | 同行評審論文 | 高 | ASM Atlas 最全面參考 |
| CanASM - 2025, BMC Genomics | 同行評審論文 | 中-高 | Cancer ASM 資料庫 |
| MHB map - Lee et al. 2025, Cell Reports | 同行評審論文 | 高 | Solid cancer MHB 最新參考 |
| Shoemaker et al. 2010, Genome Research | 同行評審論文 | 高 | ASM 基線數據經典參考 |
| WSHPackage / epihet | 同行評審論文 + R package | 中-高 | 異質性 scores 綜合比較 |
| Methrix - 2020, Bioinformatics | 同行評審論文 | 中 | 矩陣管理工具，非直接相關 |

---

## 8. 總結

本調查涵蓋 30+ 種 ASM 偵測與亞克隆甲基化分析方法。核心發現：

1. **ISM 的 read x CpG 矩陣是大多數方法的共同基礎**，至少 13 種方法可在此矩陣上直接計算新指標。

2. **CPEL 的 entropy imbalance 發現** (96% ASM 是 entropy 差異) 為 ISM 的 PERMANOVA 提供理論支持 -- PERMANOVA 對 pattern variance 敏感，正好捕捉此類 ASM。

3. **Per-CpG Fisher's exact test** 是最簡單且最廣泛使用的 ASM 偵測方法，可在 ISM 現有架構上零成本添加。

4. **Normal baseline 比較** (CAMDAC/Do et al. 框架) 與 ISM Phase 2A 方向完全對齊，是下一步最重要的方法學整合。

5. **ISM 的 SNV-anchored + HP-integrated + LOH-aware 組合在學界是獨有的**，沒有任何現有工具提供相同維度的分析。

6. **Epiallele diversity metrics** (epipolymorphism, Shannon entropy, VEF) 可作為 ISM 的 region-level 特徵，豐富 per-region 的表觀遺傳特徵化。
