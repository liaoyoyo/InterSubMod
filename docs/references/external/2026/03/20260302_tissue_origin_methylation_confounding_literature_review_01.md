<!--
建立時間: 2026-03-02 15:00
目標: 系統性收集與分析「組織來源甲基化差異對腫瘤亞克隆分析的混淆效應」相關文獻
處理範圍: 五大主題 -- (1)組織特異性甲基化差異 (2)Long-read甲基化分析進展 (3)Tissue-of-origin confounding解決方案 (4)Truth set與benchmarking (5)亞克隆甲基化分析工具
關聯檔案:
  - /big8_disk/liaoyoyo2001/Knowledge/02_samples/HCC1395.md
  - /big8_disk/liaoyoyo2001/Knowledge/02_samples/subsample_purity.md
  - /big8_disk/liaoyoyo2001/Knowledge/06_workflows/methylation_analysis.md
  - /big8_disk/liaoyoyo2001/Knowledge/08_references/paper_index.md
-->

# 組織來源甲基化混淆效應 文獻綜述報告

## 1. 搜尋概述

- **研究背景**: 本實驗室使用 Oxford Nanopore 長讀取測序進行腫瘤亞克隆甲基化分析（InterSubMod）。在建立不同純度（purity）的 subsample 時，將 tumor BAM（來自癌旁組織 adjacent tissue / cell line）與 normal BAM（來自血液 B-lymphoblast）混合。由於這兩種組織的甲基化模式存在本質差異（tissue-of-origin methylation differences），導致分析結果被組織來源差異所混淆（confounding），而非真正反映腫瘤亞克隆的表觀遺傳異質性。
- **搜尋時間**: 2026-03-02
- **搜尋範圍**: PubMed, Google Scholar, bioRxiv/medRxiv, Nature系列期刊
- **文獻年代**: 2023-2026（以 2024-2025 為主）
- **關鍵詞組合**:
  - tissue-of-origin methylation differences / tumor normal confounding
  - blood vs tissue methylation / tumor purity adjustment
  - nanopore long-read methylation / subclonal analysis
  - methylation deconvolution / tissue-specific correction
  - SEQC2 / CASTLE / somatic variant benchmark
  - subclonal methylation tools / allele-specific / single molecule
- **資料來源數**: 30+ 篇文獻

---

## 2. 核心發現

### 2.1 主題一：組織特異性甲基化差異對腫瘤分析的影響

#### 2.1.1 血液 vs 實體組織的甲基化差異

| 文獻 | 核心發現 | 與本問題的關聯 |
|------|---------|--------------|
| **Zheng et al. (2017)** - InfiniumPurify, Genome Biology | 某些癌種（如 DLBC, LAML, THYM）從血液對照推算的純度與通用正常對照的估計值相關性極差，因為這些組織的甲基化特徵與血液差異極大，使得差異甲基化 CpG 多為血液組織特異性而非真正的腫瘤-正常差異 | **直接相關**: HCC1395BL 是 B-lymphoblast（血液來源），而 HCC1395 是乳腺導管癌，兩者甲基化基線差異極大 |
| **Moss et al. (2018)** - Nature Communications | 建立了涵蓋 39 種人類細胞型態的甲基化圖譜，發現不同組織間存在數千個獨特的甲基化標記，血液細胞與實體組織的甲基化模式存在系統性差異 | **直接相關**: 提供了血液 vs 乳腺組織差異的量化基礎 |

**關鍵認知**:
> 當使用血液來源的 normal 作為對照，分析實體腫瘤的甲基化差異時，所觀察到的「腫瘤-正常差異」中會混入大量的「組織來源差異」。這不是分析方法的問題，而是實驗設計的根本限制。

#### 2.1.2 癌旁組織（Adjacent Tissue）的 Field Effect

| 文獻 | 核心發現 | 與本問題的關聯 |
|------|---------|--------------|
| **Johnson et al. (2025)** - Clinical Epigenetics | 在前列腺癌中發現，正常組織的差異甲基化模式不僅存在於腫瘤鄰近區域（<1mm），甚至延伸至距腫瘤 4cm 的遠端區域，提供了表觀遺傳 field effect 的證據 | 若使用癌旁組織作為 normal，其甲基化可能已受 field effect 影響而非真正「正常」 |
| **Li et al. (2024)** - JNCI | 在結直腸癌中，正常盲腸組織的甲基化已呈現系統性失調（558 個 unique loci），沿腺瘤-癌序列存在進行性甲基化擾亂 | Field effect 使得「正常對照」的定義更為複雜 |
| **Wang et al. (2025)** - Clinical Epigenetics | DNA 甲基化異質性與肺腺癌的 field cancerization 及預後相關，腫瘤鄰近正常組織的表觀遺傳改變可作為風險標記 | 進一步支持 adjacent tissue 並非理想的 normal 對照 |
| **Flanagan et al. (2018)** - EBioMedicine | 正常乳腺組織的甲基化模式比拷貝數變異更能預測乳腺癌狀態，暗示正常組織已包含早期癌化的表觀遺傳改變 | HCC1395 的 adjacent tissue 同樣可能受此影響 |

#### 2.1.3 Cell Line 配對的特殊考量

| 文獻 / 背景知識 | 核心發現 | 與本問題的關聯 |
|------|---------|--------------|
| **SEQC2 (Fang et al., 2021)** - Nature Biotechnology | HCC1395 為乳腺導管癌 cell line，配對 normal 為 HCC1395BL（同一患者的 B-lymphoblast cell line） | **HCC1395（乳腺）vs HCC1395BL（血液）的甲基化差異是我們問題的核心** |
| 根據 Knowledge/02_samples/HCC1395.md | Subsample 混合使用 tumor BAM + normal BAM，而現有 ONT subsample 缺乏 MM/ML tags | 需從 ONT_5kHz 資料集重建含甲基化標籤的 purity subsample |

---

### 2.2 主題二：Long-read 測序甲基化分析的最新進展

#### 2.2.1 Nanopore 甲基化檢測技術進展

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Boltasseva et al. (2025)** "Reliable investigation of DNA methylation using ONT" | Scientific Reports, 2025 | R10.4.1 流通池較 R9.4.1 的序列準確度顯著提升（97.1% vs 93.1%），研究跨化學版本甲基化數據的一致性與潛在偏差 | InterSubMod 使用 ONT 5kHz simplex 資料，需注意不同 chemistry 版本的甲基化檢測差異 |
| **Guichou et al. (2024)** "DeepMod2" | Nature Communications, 2024 | 基於深度學習的甲基化與 epihaplotype 快速準確檢測方法，適用於 ONT 資料 | 可考慮作為甲基化分析的替代或驗證工具 |
| **Dyshlovoy et al. (2024)** "Applications of Nanopore sequencing in precision cancer medicine" | Int J Cancer, 2024 | 綜述 Nanopore 在精準癌症醫學中的應用，包括腫瘤分類、甲基化分析、SV 檢測 | 提供了 Nanopore 甲基化分析在臨床應用中的完整生態 |
| **Vetter & Aganezov (2024)** "Shedding light on DNA methylation and its clinical implications" | Epigenetics & Chromatin, 2024 | 綜合評述 ONT 技術對甲基化研究的影響，R10 chemistry 已成主流 | 確認我們使用的技術平台符合最新發展趨勢 |
| **Lim & He (2025)** "Computational analysis of DNA methylation from long-read sequencing" | Nature Reviews Genetics, 2025 | 長讀取甲基化分析計算工具的全面綜述，涵蓋信號識別、比較分析、細胞類型多樣性分析 | **重要參考**: 最新且最全面的方法論綜述 |

#### 2.2.2 腫瘤表觀遺傳異質性的 Single-Molecule 分析

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Stoeger et al. (2025)** "Long-read sequencing of melanoma subclones" | bioRxiv/PMC, 2025 | 使用 ONT R10.4 對 23 個黑色素瘤單細胞衍生亞克隆進行長讀取測序，發現品系特異性 DMR 與侵襲性表型相關，甲基化軌跡與 SNV 積累平行演化 | **高度相關**: 直接展示了 Nanopore 在亞克隆甲基化分析中的應用，與 InterSubMod 的目標高度一致 |
| **Gigante et al. (2022)** "Epigenetic tumor heterogeneity in the era of single-cell profiling with nanopore sequencing" | Clinical Epigenetics, 2022 | 表觀遺傳異質性在腫瘤生物學變異和臨床預後中的重要性日益被認識，單分子/單鹼基解析度的評估至關重要 | 提供了 InterSubMod 研究方向的理論支持 |
| **Luo & Zhang (2024)** "Long-read sequencing unveils novel somatic variants and methylation patterns in early lung cancer" | Computers in Biology and Medicine, 2024 | 長讀取測序提供肺癌表觀遺傳與基因組改變的整體視角，整合甲基化與變異資訊揭示癌症相關通路 | 整合分析方法的參考 |

#### 2.2.3 MethPhaser：甲基化輔助 Phasing

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Kim et al. (2024)** "MethPhaser: methylation-based long-read haplotype phasing of human genomes" | Nature Communications, 2024 | 利用 ONT 甲基化信號擴展 SNV-based phasing，對 R9/R10 資料可將 phase block N50 提升 78%-151%（準確度 83.4-98.7%），改善 HLA 等醫學相關基因的 phasing | **高度相關**: InterSubMod 依賴 LongPhase 的 phasing/haplotagging，MethPhaser 的甲基化 phasing 概念可互補 |

---

### 2.3 主題三：解決 Tissue-of-Origin Confounding 的方法

#### 2.3.1 甲基化數據的純度校正

| 文獻 | 期刊/年份 | 方法概述 | 優缺點 | 關聯性 |
|------|----------|---------|--------|--------|
| **Staaf et al. (2024)** "PureBeta" | NAR Genomics and Bioinformatics, 2024 | 單樣本統計框架，從全基因組甲基化數據估計純度，再校正個別 CpG 的 beta 值。三模組架構：參考建立、純度估計、Beta 校正 | 優：不需匹配 normal，與測序純度估計相關性 >0.8。缺：需大型參考隊列，對 cell line 不適用，tissue-specific | **中度相關**: 概念可借鑒，但 PureBeta 設計於 Illumina array，非 long-read |
| **Zheng et al. (2017)** "InfiniumPurify" | Genome Biology, 2016 | 從差異甲基化位點估計腫瘤純度，進行統計校正 | 優：成熟方法。缺：在低純度時高估，需要 tissue-specific 參考 | 奠基性工作，但有 tissue-specific bias 問題 |
| **Aran et al. (2022)** | PLOS ONE, 2022 | 純度校正後的 beta 值改善高維甲基化數據的生物學可解釋性 | 校正後 beta 值呈更二元分布，更接近理論甲基化狀態 | 確認純度校正的必要性 |

#### 2.3.2 甲基化解摺積（Deconvolution）方法

| 文獻 | 期刊/年份 | 方法概述 | 關鍵特點 | 關聯性 |
|------|----------|---------|---------|--------|
| **Zhu et al. (2024)** "CelFiE-ISH" | Genome Biology, 2024 | 利用單分子甲基化 haplotype 進行多細胞類型解摺積的概率模型（EM 框架）。考慮 read 內多個 CpG 的聯合概率而非獨立處理 | 比 CelFiE 準確度提升 30%，可檢測低至 0.03% 的稀有細胞類型 | **高度相關**: 直接利用 long-read 的單分子特性進行甲基化解摺積，可用於區分 tumor/blood 來源的 reads |
| **Zhang et al. (2024)** "MEnet" | NAR Cancer, 2024 | 基於神經網路的甲基化解摺積工具，透過數據增強提升穩健性，**唯一**跨平台支持（nanopore, WGBS, RRBS, methylation arrays） | 跨平台兼容性，適用於不同實驗設計 | **直接相關**: 支援 nanopore 數據的解摺積 |
| **Luo et al. (2024)** "Systematic evaluation of methylation-based cell type deconvolution methods for plasma cfDNA" | Genome Biology, 2024 | 系統性比較 5 種主要 cfDNA 解摺積方法（MethAtlas, cfNOMe toolkit, CelFiE, CelFEER, UXM） | 參考標記選擇、測序深度、參考圖譜完整性均影響性能 | 提供方法選擇的指引 |
| **Loyfer et al. (2023)** "cfSort" | PNAS, 2023 | 基於深度學習的組織解摺積，全面超越既有方法 | 高靈敏度和準確度 | 替代方法參考 |
| **MethylResolver (Titus et al., 2020)** | Communications Biology, 2020 | 基於 Least Trimmed Squares 迴歸的白血球亞群推斷方法，不需要癌症特異性簽名即可解析純度 | 隨未知成分增加而準確度提升 | 方法論參考 |
| **Zhao et al. (2024)** "MetDecode" | Bioinformatics, 2024 | 甲基化解摺積方法，可即時學習圖譜中缺失的成分，考慮每個標記區域的覆蓋度 | 處理不完整參考圖譜的問題 | 穩健性設計值得借鑒 |

#### 2.3.3 組織來源識別分類器

| 文獻 | 期刊/年份 | 方法概述 | 關聯性 |
|------|----------|---------|--------|
| **Loft et al. (2024)** "CUPiD Classifier" | Nature Communications, 2024 | 基於 cfDNA 甲基化的機器學習分類器，跨 29 個腫瘤類型預測組織來源，整體靈敏度 84.6%，組織來源準確度 96.8% | 證明甲基化信號足以區分組織來源 |
| **Nguyen et al. (2024)** "TSMA + GCNN" | J Translational Medicine, 2024 | 結合腫瘤特異性甲基化圖譜（TSMA）與全基因組甲基化密度的圖卷積神經網路方法。TSMA 在腫瘤組織中表現良好但在 cfDNA 中受 WBC DNA 壓制 | 說明了血液成分對甲基化分析的干擾程度 |
| **Li et al. (2024)** "MFCUP" | Clinical Epigenetics, 2024 | 200-CpG 甲基化特徵分類器，25 種癌症類型的驗證準確度 97.2% | 少量 CpG 即可區分組織來源 |

#### 2.3.4 Nanopore 特異的解摺積方法

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Ni et al. (2023)** "Single-molecule methylation profiles of cfDNA in cancer with nanopore sequencing" | Genome Medicine, 2023 | 利用 nanopore cfDNA 的單分子甲基化特徵區分腫瘤來源，不同亞群的腫瘤細胞在甲基化模式上表現顯著異質性 | **直接相關**: nanopore 單分子甲基化用於腫瘤異質性分析的先驅工作 |
| **Han et al. (2022)** "Detecting cell-of-origin and cancer-specific methylation features of cfDNA from Nanopore sequencing" | Genome Biology, 2022 | 從 nanopore cfDNA 數據直接檢測細胞來源和癌症特異性甲基化特徵 | 方法論可借鑒至 InterSubMod 的 read-level 分析 |

---

### 2.4 主題四：Truth Set 與 Benchmarking 的挑戰

#### 2.4.1 SEQC2 的局限性

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Fang et al. (2021)** "SEQC2" | Nature Biotechnology, 2021 | HCC1395-HCC1395BL 是**目前唯一**公開可用的高品質體細胞變異 benchmark。僅包含約 40,000 個 SNV + 2,000 個 Indel，覆蓋約 82% 基因組的高信度區域 | 是我們的主要 benchmark，但代表性有限 |
| 局限性分析 | - | (1) HCC1395 基因組不一定準確代表 TNBC；(2) Cell line 經長期培養可能累積人為突變；(3) Germline variant caller 有 GIAB 7 個參考樣本、數百萬 variant，而 somatic 僅有此一組 | benchmark 結果的泛化能力受限 |

#### 2.4.2 CASTLE 數據集與新 Benchmark

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Park & Cook et al. (2025)** "DeepSomatic" | Nature Biotechnology, 2025 | Google 開放 CASTLE (Cancer Standards Long-read Evaluation) 數據集：6 對匹配 tumor-normal cell line pairs，使用 Illumina + PacBio HiFi + ONT 全基因組測序，329,011 個體細胞變異 | **直接相關**: 大幅擴展了 long-read somatic benchmark 資源 |
| **改善 benchmark 的號召** | Nature Biotechnology, 2025 | 社論呼籲改善體細胞變異 benchmark，指出訓練和測試資源的匱乏是整個領域的瓶頸 | 說明 benchmark 不足是公認問題 |

#### 2.4.3 Two-Tech Truth Set

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Lancet2 (NYGC, 2025)** | bioRxiv, 2025 | 為 HCC1187, HCC1143, COLO829, HCC1395 四個 cell line 建立了 "two-tech" truth set，整合高覆蓋 Illumina + ONT 數據。Lancet2 使用可解釋機器學習模型訓練於 HCC1395 truth set | **重要**: 提供了整合 short-read + long-read 的新 benchmark 標準 |

#### 2.4.4 Severus -- SV Benchmark

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Severus (Kolmogorov et al., 2025)** | Nature Biotechnology, 2025 | 基於斷點圖的體細胞 SV 檢測算法，支援 haplotype-specific 調用，在綜合多技術 cell line panel 上 F1 score 優於其他方法 | 展示了長讀取在複雜變異檢測中的優勢 |

---

### 2.5 主題五：亞克隆甲基化分析工具與方法

#### 2.5.1 單分子甲基化分析工具

| 文獻/工具 | 期刊/年份 | 功能概述 | 關聯性 |
|----------|----------|---------|--------|
| **NanoMethViz (v3.3.3)** | F1000Research, 2024 | Bioconductor 套件，處理 long-read 甲基化數據，支援 modBAM 格式和 ONT modkit 輸出，提供多解析度視覺化 | 可作為 InterSubMod 視覺化的補充工具 |
| **NanoNOMe** | - | 利用長讀取（>10kb）同時測量甲基化和染色質可及性，評估 allele-specific 表觀遺傳狀態 | 方法論參考 |
| **SMRT-Tag & SAMOSA-Tag (2025)** | Gladstone Institutes, 2025 | 單分子分析工具，DNA 需求量降低 90-95%，結合甲基化與染色質可及性分析，已應用於前列腺癌轉移研究 | 前沿技術方向 |

#### 2.5.2 Epiallele 檢測與亞克隆分析

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Landan et al. (2017)** "BADER" | BMC Bioinformatics, 2017 | 基於貝葉斯的 epiallele 檢測方法，可區分具有不同 epialleles 的細胞亞群，識別腫瘤亞克隆 | 與 InterSubMod 的目標高度一致 -- 從甲基化模式識別亞克隆 |
| **Hifimeth (2024)** | bioRxiv, 2024 | PacBio HiFi 數據的高準確度甲基化識別工具，可檢測 allelic DMR（15,908 個，>=20% allelic 差異） | 方法可借鑒至 ONT 平台 |

#### 2.5.3 Allele-Specific 甲基化分析

| 文獻 | 期刊/年份 | 核心發現 | 關聯性 |
|------|----------|---------|--------|
| **Zuo et al. (2024)** | Am J Human Genetics, 2024 | 探討單細胞甲基化數據檢測 allele-specific methylation (ASM) 和 imprinting 的潛力 | ASM 分析方法可整合至 InterSubMod |
| **scDEEP-mC (2025)** | Nature Communications, 2025 | 改進的單細胞全基因組亞硫酸鹽測序方法，可分析 allele-specific methylation、複製動態和 X-inactivation，識別相對均質群體內的亞克隆 | 展示了從甲基化角度發現亞克隆的可行性 |
| **Long-Read POG Cohort (2024)** | Cell Genomics, 2024 | 晚期癌症隊列的 nanopore 測序，揭示 allele-specific methylation 在癌症基因中的角色，支持長讀取在個體化癌症醫學中的應用 | 臨床應用驗證 |

---

## 3. 衝突觀點與討論

### 3.1 Normal 對照的選擇

| 觀點 A：使用匹配血液作為 Normal | 觀點 B：使用匹配組織作為 Normal | 可能原因/權衡 |
|------|------|------|
| 血液易取得，基因型匹配保證（同一患者），是體細胞變異 calling 的標準做法（SEQC2 設計） | 組織匹配可避免 tissue-of-origin 甲基化差異的混淆，但受 field effect 影響可能不「正常」 | **基因變異分析**（SNV/Indel）選血液較佳（因為只看序列差異）；**甲基化分析**時血液對照會引入大量系統性偏差 |

### 3.2 Cell Line vs 臨床樣本

| 觀點 A：Cell Line 是穩定的 Benchmark | 觀點 B：Cell Line 代表性有限 | 可能原因/權衡 |
|------|------|------|
| Cell line 基因型穩定，可重複，SEQC2 提供高信度 truth set | Cell line 經長期培養，基因組/表觀基因組可能偏離原始腫瘤。HCC1395 不一定代表 TNBC。缺乏微環境和免疫浸潤 | 對 variant calling benchmark 而言 cell line 足夠；對甲基化分析而言，cell line 的表觀遺傳狀態可能已偏移，且缺乏腫瘤微環境的複雜性 |

### 3.3 純度校正的適用性

| 觀點 A：全局純度校正（PureBeta 式） | 觀點 B：Read-level 解摺積（CelFiE-ISH 式） | 可能原因/權衡 |
|------|------|------|
| 基於全基因組甲基化估計樣本純度，再統一校正所有 CpG 的 beta 值 | 利用單分子甲基化 haplotype 在 read 層級判定來源 | 全局校正假設純度均勻，忽略區域性異質性；Read-level 方法可捕捉局部異質性，但需要良好的參考甲基化圖譜且計算量大 |

---

## 4. 資料來源評估

### 4.1 高可信度來源（Tier 1 -- 同行評審期刊論文）

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| [PureBeta (Staaf et al., 2024)](https://academic.oup.com/nargab/article/6/4/lqae146/7874836) | 期刊論文 | 高 | NAR Genomics and Bioinformatics, 有 GitHub 可重現 |
| [CelFiE-ISH (Zhu et al., 2024)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03275-x) | 期刊論文 | 高 | Genome Biology, IF > 10 |
| [DeepSomatic (Park et al., 2025)](https://www.nature.com/articles/s41587-025-02839-x) | 期刊論文 | 高 | Nature Biotechnology, Google 團隊 |
| [SEQC2 (Fang et al., 2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/) | 期刊論文 | 高 | Nature Biotechnology, FDA 主導 |
| [MethPhaser (Kim et al., 2024)](https://www.nature.com/articles/s41467-024-49588-0) | 期刊論文 | 高 | Nature Communications |
| [CUPiD (Loft et al., 2024)](https://www.nature.com/articles/s41467-024-47195-7) | 期刊論文 | 高 | Nature Communications |
| [InfiniumPurify (Zheng et al., 2017)](https://link.springer.com/article/10.1186/s13059-016-1143-5) | 期刊論文 | 高 | Genome Biology, 奠基性工作 |
| [DNA methylation atlas (Loyfer et al., 2023)](https://www.nature.com/articles/s41586-022-05580-6) | 期刊論文 | 高 | Nature |
| [Melanoma subclones (Stoeger et al., 2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424993/) | 期刊/預印本 | 高 | PMC 收錄，ONT R10.4 |
| [Shedding light on ONT methylation (Vetter et al., 2024)](https://link.springer.com/article/10.1186/s13072-024-00558-2) | 期刊論文 | 高 | Epigenetics & Chromatin, 綜述 |
| [Computational analysis of DNA methylation from long-read (Lim & He, 2025)](https://www.nature.com/articles/s41576-025-00822-5) | 期刊論文 | 高 | Nature Reviews Genetics, 權威綜述 |
| [Severus (Kolmogorov et al., 2025)](https://www.nature.com/articles/s41587-025-02618-8) | 期刊論文 | 高 | Nature Biotechnology |
| [ClairS-TO (Chen et al., 2025)](https://www.nature.com/articles/s41467-025-64547-z) | 期刊論文 | 高 | Nature Communications |
| [Prostate cancer field effect (Johnson et al., 2025)](https://link.springer.com/article/10.1186/s13148-025-01932-x) | 期刊論文 | 高 | Clinical Epigenetics |
| [Lung cancer field cancerization (Wang et al., 2025)](https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-025-01845-9) | 期刊論文 | 高 | Clinical Epigenetics |

### 4.2 中等可信度來源（Tier 2 -- 預印本或較舊文獻）

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| [Lancet2 two-tech truth set (2025)](https://www.biorxiv.org/content/10.1101/2025.02.18.638852v2.full) | bioRxiv 預印本 | 中-高 | NYGC 團隊，有公開數據 |
| [DeepMod2 (Guichou et al., 2024)](https://www.nature.com/articles/s41467-024-45778-y) | 期刊論文 | 中-高 | Nature Communications |
| [MEnet (Zhang et al., 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11094754/) | 期刊論文 | 中 | NAR Cancer |
| [Hifimeth (2024)](https://www.biorxiv.org/content/10.1101/2024.08.14.607879v1.full) | bioRxiv 預印本 | 中 | 未經同行評審 |
| [cfDNA nanopore methylation (Ni et al., 2023)](https://pubmed.ncbi.nlm.nih.gov/37138315/) | 期刊論文 | 中-高 | Genome Medicine |

---

## 5. 針對 InterSubMod 的具體建議

### 5.1 短期可行方案（1-2 周）

#### 方案 A：Read-Level 組織來源標記
- **概念**: 在分析前，根據每條 read 的甲基化 haplotype 特徵判定其可能來自 tumor（乳腺組織）或 normal（血液），並標記或過濾
- **依據**: CelFiE-ISH 已證明可從單分子甲基化 haplotype 進行細胞類型解摺積，準確度提升 30%
- **實作方式**:
  1. 建立血液 vs 乳腺組織的甲基化參考圖譜（使用公開的 methylation atlas 數據）
  2. 對每條 read 計算其與各組織類型的甲基化相似度
  3. 在 InterSubMod 的顯著性分析中加入組織來源作為共變量

#### 方案 B：使用同組織 Normal 對照
- **概念**: 避免使用血液 normal，改用同組織類型的 normal 作為基線
- **依據**: 根據 Knowledge/02_samples/HCC1395.md，目前使用的 normal 是 HCC1395BL（B-lymphoblast）
- **限制**: 對 cell line 而言沒有 matched adjacent tissue，但可考慮使用公開的正常乳腺組織甲基化數據作為參考

#### 方案 C：甲基化差異的組織來源校正
- **概念**: 在距離計算前，先減去已知的組織來源甲基化差異
- **實作方式**:
  1. 計算 pure tumor BAM 與 pure normal BAM 在每個 CpG 位點的平均甲基化水平
  2. 識別「組織差異 CpG」（在純 tumor 和純 normal 之間差異顯著但與腫瘤狀態無關的位點）
  3. 在 InterSubMod 的距離計算中排除或降權這些 CpG

### 5.2 中期方案（1-3 個月）

#### 方案 D：整合 PureBeta 式的純度感知分析
- **概念**: 在 InterSubMod 中加入純度估計模組，根據已知或估計的純度調整甲基化距離計算
- **依據**: PureBeta 已證明純度校正可使 beta 值更接近理論狀態
- **挑戰**: 需要將 array-based 方法改編為適用於 long-read 的 read-level 版本

#### 方案 E：建立 Nanopore 特異的組織甲基化參考
- **概念**: 使用 pure tumor (t50_n00) 和 pure normal (t00_n25) 的 ONT_5kHz 資料建立 Nanopore 平台特異的甲基化參考圖譜
- **注意**: 根據 Knowledge，現有 subsample BAM 的 MM/ML tags 為 0，需從 ONT_5kHz 原始資料重新建立

### 5.3 長期方向

#### 方向 F：Multi-modal 亞克隆分析
- **概念**: 借鑒 Stoeger et al. (2025) 的黑色素瘤亞克隆研究，整合 SNV + SV + CNA + methylation 的多層次分析
- **依據**: 該研究證明甲基化軌跡與基因組變異平行演化，多模態整合可提供更可靠的亞克隆結構推斷

#### 方向 G：利用甲基化增強 Phasing
- **概念**: 借鑒 MethPhaser 的方法，利用甲基化信號輔助 LongPhase 的 somatic haplotagging
- **潛在效益**: 更多 reads 被正確 haplotag，提升 InterSubMod 的統計效力

---

## 6. 完整文獻清單

### 6.1 組織特異性甲基化與混淆效應

1. **Zheng Q, et al.** (2017) "Estimating and accounting for tumor purity in the analysis of DNA methylation data from cancer studies." *Genome Biology* 17:164. DOI: [10.1186/s13059-016-1143-5](https://link.springer.com/article/10.1186/s13059-016-1143-5)

2. **Moss J, et al.** (2018) "Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease." *Nature Communications* 9:5068. DOI: [10.1038/s41467-018-07466-6](https://www.nature.com/articles/s41467-018-07466-6)

3. **Loyfer N, et al.** (2023) "A DNA methylation atlas of normal human cell types." *Nature* 613:355-364. DOI: [10.1038/s41586-022-05580-6](https://www.nature.com/articles/s41586-022-05580-6)

4. **Johnson KC, et al.** (2025) "DNA methylation in normal-appearing tissue associated with prostate cancer recurrence and metastasis." *Clinical Epigenetics*. DOI: [10.1186/s13148-025-01932-x](https://link.springer.com/article/10.1186/s13148-025-01932-x)

5. **Wang Y, et al.** (2025) "DNA methylation heterogeneity correlates with field cancerization and prognosis in lung adenocarcinoma patients." *Clinical Epigenetics*. DOI: [10.1186/s13148-025-01845-9](https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-025-01845-9)

6. **Li X, et al.** (2024) "DNA-methylation variability in normal mucosa: a field cancerization marker in patients with adenomatous polyps." *JNCI* 116(6):974. DOI: [10.1093/jnci/djae001](https://academic.oup.com/jnci/article/116/6/974/7589935)

7. **Flanagan JM, et al.** (2018) "DNA Methylation Patterns in Normal Tissue Correlate more Strongly with Breast Cancer Status than Copy-Number Variants." *EBioMedicine*. DOI: [10.1016/j.ebiom.2018.05.030](https://www.sciencedirect.com/science/article/pii/S2352396418301531)

### 6.2 Long-Read 甲基化分析進展

8. **Lim F, He C.** (2025) "Computational analysis of DNA methylation from long-read sequencing." *Nature Reviews Genetics*. DOI: [10.1038/s41576-025-00822-5](https://www.nature.com/articles/s41576-025-00822-5)

9. **Boltasseva E, et al.** (2025) "Reliable investigation of DNA methylation using Oxford nanopore technologies." *Scientific Reports*. DOI: [10.1038/s41598-025-99882-0](https://www.nature.com/articles/s41598-025-99882-0)

10. **Guichou JF, et al.** (2024) "A signal processing and deep learning framework for methylation detection using Oxford Nanopore sequencing (DeepMod2)." *Nature Communications* 15:1404. DOI: [10.1038/s41467-024-45778-y](https://www.nature.com/articles/s41467-024-45778-y)

11. **Dyshlovoy SA, et al.** (2024) "Applications of Nanopore sequencing in precision cancer medicine." *Int J Cancer*. DOI: [10.1002/ijc.35100](https://onlinelibrary.wiley.com/doi/10.1002/ijc.35100)

12. **Vetter M, Aganezov S.** (2024) "Shedding light on DNA methylation and its clinical implications: the impact of long-read-based nanopore technology." *Epigenetics & Chromatin*. DOI: [10.1186/s13072-024-00558-2](https://link.springer.com/article/10.1186/s13072-024-00558-2)

13. **Stoeger MK, et al.** (2025) "Long-read sequencing of single cell-derived melanoma subclones reveals divergent and parallel genomic and epigenomic evolutionary trajectories." *bioRxiv/PMC*. [PMC12424993](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424993/)

14. **Gigante S, et al.** (2022) "Epigenetic tumor heterogeneity in the era of single-cell profiling with nanopore sequencing." *Clinical Epigenetics* 14:110. DOI: [10.1186/s13148-022-01323-6](https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-022-01323-6)

15. **Luo W, Zhang S.** (2024) "Long-read sequencing unveils novel somatic variants and methylation patterns in the genetic information system of early lung cancer." *Computers in Biology and Medicine*. DOI: [10.1016/j.compbiomed.2024.108124](https://www.sciencedirect.com/science/article/abs/pii/S0010482524002580)

16. **Kim G, et al.** (2024) "MethPhaser: methylation-based long-read haplotype phasing of human genomes." *Nature Communications* 15:5262. DOI: [10.1038/s41467-024-49588-0](https://www.nature.com/articles/s41467-024-49588-0)

### 6.3 解摺積與純度校正方法

17. **Staaf J, et al.** (2024) "Tumor purity estimated from bulk DNA methylation can be used for adjusting beta values of individual samples to better reflect tumor biology (PureBeta)." *NAR Genomics and Bioinformatics* 6(4):lqae146. DOI: [10.1093/nargab/lqae146](https://academic.oup.com/nargab/article/6/4/lqae146/7874836)

18. **Zhu MC, et al.** (2024) "CelFiE-ISH: a probabilistic model for multi-cell type deconvolution from single-molecule DNA methylation haplotypes." *Genome Biology* 25:151. DOI: [10.1186/s13059-024-03275-x](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03275-x)

19. **Zhang X, et al.** (2024) "Neural-net-based cell deconvolution from DNA methylation reveals tumor microenvironment associated with cancer prognosis (MEnet)." *NAR Cancer* 6(2):zcae022. DOI: [10.1093/narcan/zcae022](https://academic.oup.com/narcancer/article/6/2/zcae022/7673516)

20. **Luo Y, et al.** (2024) "Systematic evaluation of methylation-based cell type deconvolution methods for plasma cell-free DNA." *Genome Biology*. DOI: [10.1186/s13059-024-03456-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03456-8)

21. **Loyfer N, et al.** (2023) "Comprehensive tissue deconvolution of cell-free DNA by deep learning for disease diagnosis and monitoring (cfSort)." *PNAS*. DOI: [10.1073/pnas.2305236120](https://www.pnas.org/doi/10.1073/pnas.2305236120)

22. **Titus AJ, et al.** (2020) "MethylResolver -- a method for deconvoluting bulk DNA methylation profiles into known and unknown cell contents." *Communications Biology* 3:422. DOI: [10.1038/s42003-020-01146-2](https://www.nature.com/articles/s42003-020-01146-2)

23. **Zhao J, et al.** (2024) "MetDecode: methylation-based deconvolution of cell-free DNA for noninvasive multi-cancer typing." *Bioinformatics* 40(9):btae522. DOI: [10.1093/bioinformatics/btae522](https://academic.oup.com/bioinformatics/article/40/9/btae522/7739698)

### 6.4 組織來源分類器

24. **Loft A, et al.** (2024) "A cfDNA methylation-based tissue-of-origin classifier for cancers of unknown primary (CUPiD)." *Nature Communications*. DOI: [10.1038/s41467-024-47195-7](https://www.nature.com/articles/s41467-024-47195-7)

25. **Nguyen HN, et al.** (2024) "Tissue of origin detection for cancer tumor using low-depth cfDNA samples through combination of TSMA and GCNN." *J Translational Medicine*. DOI: [10.1186/s12967-024-05416-z](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-024-05416-z)

26. **Li Y, et al.** (2024) "Tissue of origin prediction for cancer of unknown primary using a targeted methylation sequencing panel (MFCUP)." *Clinical Epigenetics*. DOI: [10.1186/s13148-024-01638-6](https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-024-01638-6)

27. **Han LKM, et al.** (2022) "Detecting cell-of-origin and cancer-specific methylation features of cfDNA from Nanopore sequencing." *Genome Biology*. DOI: [10.1186/s13059-022-02710-1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02710-1)

28. **Ni Y, et al.** (2023) "Single-molecule methylation profiles of cell-free DNA in cancer with nanopore sequencing." *Genome Medicine*. DOI: [10.1186/s13073-023-01178-3](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01178-3)

### 6.5 Benchmark 與 Truth Set

29. **Fang LT, et al.** (2021) "Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing (SEQC2)." *Nature Biotechnology* 39:1151-1160. DOI: [10.1038/s41587-021-00993-6](https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/)

30. **Park J, Cook DE, et al.** (2025) "Accurate somatic small variant discovery for multiple sequencing technologies with DeepSomatic." *Nature Biotechnology*. DOI: [10.1038/s41587-025-02839-x](https://www.nature.com/articles/s41587-025-02839-x)

31. **Lancet2 team (NYGC).** (2025) "Lancet2: Improved and accelerated somatic variant calling with joint multi-sample local assembly graphs." *bioRxiv*. DOI: [10.1101/2025.02.18.638852](https://www.biorxiv.org/content/10.1101/2025.02.18.638852v2.full)

32. **Kolmogorov M, et al.** (2025) "Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing." *Nature Biotechnology*. DOI: [10.1038/s41587-025-02618-8](https://www.nature.com/articles/s41587-025-02618-8)

33. **Chen L, Zheng Z, et al.** (2025) "ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling." *Nature Communications* 16:9630. DOI: [10.1038/s41467-025-64547-z](https://www.nature.com/articles/s41467-025-64547-z)

### 6.6 亞克隆甲基化分析工具

34. **Landan G, et al.** (2017) "Quantification of tumour evolution and heterogeneity via Bayesian epiallele detection (BADER)." *BMC Bioinformatics* 18:354. DOI: [10.1186/s12859-017-1753-2](https://link.springer.com/article/10.1186/s12859-017-1753-2)

35. **Hifimeth (2024)** "High accuracy methylation identification tools on single molecular level for PacBio HiFi data." *bioRxiv*. DOI: [10.1101/2024.08.14.607879](https://www.biorxiv.org/content/10.1101/2024.08.14.607879v1.full)

36. **Zuo C, et al.** (2024) "Investigating the potential of single-cell DNA methylation data to detect allele-specific methylation and imprinting." *Am J Human Genetics*. DOI: [10.1016/j.ajhg.2024.01.014](https://www.cell.com/ajhg/fulltext/S0002-9297(24)00041-7)

37. **scDEEP-mC (2025)** "High-coverage allele-resolved single-cell DNA methylation profiling reveals cell lineage, X-inactivation state, and replication dynamics." *Nature Communications*. DOI: [10.1038/s41467-025-61589-1](https://www.nature.com/articles/s41467-025-61589-1)

38. **Long-Read POG Cohort (2024)** "Long-read sequencing of an advanced cancer cohort resolves rearrangements, unravels haplotypes, and reveals methylation landscapes." *Cell Genomics*. DOI: [10.1016/j.xgen.2024.100653](https://www.sciencedirect.com/science/article/pii/S2666979X24002933)

---

## 7. 總結

### 7.1 核心問題確認

本次文獻搜尋充分確認了以下事實：

1. **血液 vs 實體組織的甲基化差異是已知且被廣泛研究的現象**：不同組織具有數千個獨特的甲基化標記，血液（B-lymphoblast）與乳腺組織的甲基化基線差異極大。
2. **這種差異會系統性地混淆腫瘤-正常甲基化比較**：已有文獻明確指出，使用血液作為 normal 對照時，所得的差異甲基化位點多為組織特異性而非腫瘤特異性。
3. **已有多種方法可部分解決此問題**：從全局純度校正（PureBeta）到 read-level 解摺積（CelFiE-ISH），再到組織來源分類器（CUPiD），方法學工具箱已相當豐富。
4. **Long-read 單分子特性提供獨特優勢**：ONT 的 read-level 甲基化 haplotype 可用於更精確的解摺積和亞克隆分析。

### 7.2 對 InterSubMod 的最重要啟示

- **目前觀察到的甲基化「亞克隆」信號中，有一部分可能來自 tumor/normal 組織來源差異，而非真正的腫瘤內部異質性**。
- **最可行的短期改進**是建立 CpG 層級的組織來源甲基化差異圖譜，並在距離計算中排除或降權這些 CpG。
- **最有潛力的中期方向**是借鑒 CelFiE-ISH 的 read-level 解摺積概念，在 InterSubMod 中加入 read 來源判定功能。
- **Benchmark 方面**，CASTLE 數據集的公開為 long-read somatic variant calling 提供了更豐富的驗證資源。

---

*文檔版本: v1.0*
*最後更新: 2026-03-02*
