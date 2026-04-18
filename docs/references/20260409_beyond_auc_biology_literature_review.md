# ISM 甲基化特徵在 Somatic Variant 判別中的生物學基礎文獻調查

<!--
建立時間: 2026-04-09 00:00
目標: 系統性文獻調查 — 為 ISM 15+ 輪 NEGATIVE 結果（AUC < 0.58）尋找生物學解釋與方法學啟示
處理範圍: mQTL/ASM 生物學、somatic mutation 甲基化效應、ONT 技術精度、subclonal heterogeneity、統計方法學
關聯檔案:
  - docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md
  - docs/references/20260401_ASM偵測方法與預期比例文獻調查_01.md
  - docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
-->

## 1. 搜尋概述

- **關鍵詞**: mQTL, cis-mQTL, allele-specific methylation, somatic mutation methylation, ONT 5mC accuracy, per-read methylation noise, tumor methylation heterogeneity, subclonal evolution, epimutation rate, AUC vs precision-recall, residualization, confound adjustment, EWAS, biomarker evaluation
- **搜尋時間**: 2026-04-09
- **資料來源數**: 45+ 篇論文/預印本/綜述
- **涵蓋範圍**: 2010-2026，重點 2022-2026

---

## 2. 主題一：mQTL 與 ASM 的普遍性

### 2.1 cis-mQTL 的規模與效應範圍

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Oliva et al., Nature Genetics, 2023** | GTEx 9 組織 987 樣本：286,152 CpG 位點有顯著 mQTL；cis 窗口定義為 +/-500 kb；37% mCpGs 跨組織共享，5% 為組織特異；2,254 GWAS hits 共定位 | cis-mQTL 影響範圍（500 kb）遠大於 ISM 分析的 variant 周圍窗口（5 kb），意味著 germline variant 的甲基化效應可能延伸到 ISM 觀察範圍之外 |
| **Gaunt et al., Genome Biology, 2016** | 全生命週期 mQTL 分析：93% mQTL 為 cis 作用（<1 Mb）；12% SNP 直接位於 CpG 二核苷酸內（CpG-SNP）；效應隨距離衰減，近端效應最強 | CpG-SNP 機制（12%）對 ISM 最直接相關 — 當 germline SNP 破壞 CpG 位點時，甲基化變化是確定性的而非機率性的 |
| **Min et al., Nature Communications, 2024** | 歐洲（n=3,701）與東亞（n=2,099）跨群體比較：mQTL 效應量跨群體高度保守；差異主要來自 allele frequency 與 LD 結構差異 | mQTL 效應的跨群體保守性說明其生物學意義穩定，不是隨機噪聲 |
| **Huan et al., Scientific Reports, 2019** | 鑑定 55,000 個可複製 mQTL；超過 34,000 個 CpG sites 在獨立世代間驗證 | 可複製性進一步確認 germline-methylation 關聯的穩固性 |
| **IMAGE method, Genome Biology, 2019** | 整合 mQTL mapping 和 allele-specific analysis 可大幅提升 mQTL 偵測力；利用雜合個體的 ASM 資訊為核心創新 | ISM 的 haplotype-resolved 距離矩陣概念上類似於 IMAGE 的 ASM 整合分析 |

### 2.2 ASM 在正常組織與腫瘤中的比例

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Shoemaker et al., Genome Research, 2010** | 23-37% 雜合 SNP 位點在正常細胞系中表現 ASM；38-88% ASM 由 CpG-SNP 驅動；>50% ASM 區域含 SNP | 這意味著在 ISM 分析的每 3 個 germline het variant 中，就有 1 個有穩定的 allele-specific 甲基化差異 |
| **Do et al., Genome Biology, 2020** | 癌症 ASM 頻率比正常組織高 5-9 倍（MM: 5x, B-cell lymphoma: 8.5x, GBM: 9x）；**但 somatic mutation 驅動的 de novo ASM 僅佔 6-17%**；72-76% 為 allele-specific loss of methylation；allele switching 頻率：正常 14% vs 癌症 43% | **關鍵數字**：癌症中 ASM 增加大部分是全域去甲基化導致的，而非 somatic mutation 特異性效應。6-17% 的 de novo somatic ASM 解釋了為何 ISM 甲基化特徵無法區分 TP/FP |
| **Do et al., Genome Biology, 2020** | ASM delta 閾值：>20% allele 間甲基化差異為嚴格標準；de novo somatic ASM 中 71% 屬於與 germline ASM 相同的 motif class | Germline 與 somatic ASM 使用相同的生物學機制（TF binding disruption），使兩者在甲基化特徵空間中更難區分 |

### 2.3 mQTL 效應量的典型範圍

根據多篇文獻的綜合推斷：
- **強 cis-mQTL**：allele 間甲基化差異 (delta beta) 可達 20-50%，尤其是 CpG-SNP（直接破壞 CpG 位點）
- **中等 cis-mQTL**：delta beta 5-20%，通過 TF binding 間接影響
- **弱 cis-mQTL**：delta beta < 5%，數量最多但信號最弱
- **距離衰減**：效應量隨 SNP-CpG 距離增加而衰減，最強信號在 <10 kb 內

---

## 3. 主題二：Somatic Mutation 對局部甲基化的影響

### 3.1 Somatic SNV 是否改變周圍甲基化？

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Tarazona et al., Nature Genetics, 2025** | NSCLC TRACERx 研究（217 腫瘤區域 + matched normal，59 患者）：開發 ITMD（intratumoral methylation distance）量化腫瘤內甲基化異質性；ITMD 與 mutation heterogeneity 弱相關，但與 SCNA ITH 和 ITED 顯著相關；發現基因高甲基化可以是 ubiquitous 或 subclonal | ITMD 弱相關 mutation heterogeneity 的發現直接支持 ISM 觀察：somatic mutation 本身對局部甲基化的影響弱且不一致 |
| **Wen et al., PLoS Computational Biology, 2017** | 跨癌種分析：driver mutations（IDH1, SETD2, BRAF）可驅動全域甲基化重塑；但 passenger somatic SNV 對局部甲基化無顯著效應 | ISM 分析的 somatic variants 多數為 passenger（非 driver），預期其甲基化效應微弱 |
| **Haggerty et al., Genetics, 2020** | 甲基化 CpG 的 C>T 突變率高 ~10 倍（deamination）；但反過來，**甲基化降低鄰近核苷酸（+/-3 bp）的突變率** | 因果方向是 methylation -> mutation（高甲基化促進 C>T deamination），而非 mutation -> methylation change |
| **Zhou et al., Nature Genetics, 2018** | PMD 區域低甲基化與 somatic mutation 密度正相關；但這是因為 PMD 內複製錯誤率高，而非 mutation 導致低甲基化 | 相關非因果：somatic mutation 富集在已經低甲基化的區域，不是 mutation 本身改變甲基化 |

### 3.2 效應量比較：Germline mQTL vs Somatic Mutation

**綜合文獻證據**：

| 效應類型 | 預期 delta methylation | 比例 | 穩定性 |
|---------|----------------------|------|--------|
| Germline CpG-SNP | 30-50% | 12% of mQTL SNPs | 跨組織穩定 |
| Germline TF-binding mQTL | 5-20% | ~81% of mQTL | 37% 跨組織共享 |
| Somatic driver mutation (e.g., IDH1) | 10-30% (全域) | <1% of somatic mutations | 腫瘤特異 |
| Somatic passenger SNV | **<5%** (局部) | >99% of somatic mutations | 隨機、不一致 |
| Cancer global hypomethylation (background) | 20-40% (PMD) | 24-63% genome | 覆蓋性背景噪聲 |

**結論**：Germline mQTL 的 ASM 效應 >> Somatic passenger mutation 的局部效應，且方向相反（germline 產生穩定 ASM，somatic 被全域甲基化噪聲淹沒）。

### 3.3 是否有研究嘗試用甲基化區分 Somatic vs Germline？

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Baker et al. (ROCIT), bioRxiv, 2026** | Transformer 模型分類單一 read 為 tumor/non-tumor origin；AUC 0.933；使用 cell-type-specific reference methylation + sample-wide methylation distribution + normalized position；PacBio HiFi 平台 | ROCIT 分類的是 **read origin（tumor vs normal cell）**，不是 **variant origin（somatic vs germline）**。Normal cell 也攜帶 germline variant，所以不能直接等同 |
| **Jeong et al. (MethylBERT), Nature Communications, 2025** | Modified BERT 分類 read-level methylation pattern；150 bp reads 準確率 >95%；可估計 tumor cell fraction | 同 ROCIT，是 cell-of-origin 分類，非 variant 分類。但證明 read-level methylation 有足夠信息量 |

**重要澄清**：截至 2026-04，**尚無任何發表工具直接利用甲基化模式區分 somatic 與 germline variant**。ROCIT 和 MethylBERT 解決的是不同問題（cell origin vs variant origin）。

---

## 4. 主題三：長讀長測序中的甲基化分析

### 4.1 ONT 5mC Detection Accuracy

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Comprehensive Benchmarking, bioRxiv, 2024** | Dorado v4 (4kHz): per-read F1=0.93; DeepMod2 (4kHz): F1=0.88; **Dorado 5kHz_v5 優於 4kHz_v4**；per-site concordance with bisulfite: r>0.95；per-site F1 ~0.99 | Per-read F1=0.93 意味著 **每 100 個 CpG call 有 7 個錯誤**。ISM 分析依賴 read-level methylation pattern（非 site-level aggregation），因此 7% per-read 錯誤率是核心限制 |
| **DeepBAM, Briefings in Bioinformatics, 2024** | DeepBAM per-read AUC 平均 0.9847，F1 平均 0.9497；優於 Dorado（AUC 0.8825-0.9852 不等）；在不同數據集間 DeepBAM 更穩定 | DeepBAM 達到更高的單分子精度，但 ISM 目前使用 Dorado |
| **Genner et al., Genome Research, 2025** | R9 vs R10 chemistry 比較：R10 median alignment identity 98.72% vs R9 95.05%；R10 與 Illumina bisulfite 有最強 concordance；**跨 chemistry 甲基化數據可比較但需注意系統性差異** | ISM 數據使用 5kHz（R10 era），精度優於早期 R9 數據 |
| **Fu et al., Nature Reviews Genetics, 2025** | 長讀長甲基化計算分析綜述：ONT/PacBio 可同時偵測 base modification + sequence variant；per-read methylation 是長讀長的獨特優勢；挑戰包括 neighboring modification interference 和 model calibration | 綜述確認 read-level methylation 是長讀長的核心優勢，但也指出 per-read 精度仍有改善空間 |

### 4.2 Single-Read Methylation Pattern 的噪聲水準

**關鍵量化**：

- **Per-read CpG call accuracy**: ~93% (Dorado v4 4kHz), ~95% (DeepBAM)
- **Per-site accuracy** (aggregated): >99%
- **含義**: 如果一個 region 有 10 個 CpG，每個 read 的 methylation pattern 預期有 ~0.7 個 CpG 被錯誤判讀
- **對 ISM 距離計算的影響**: NHD/L1/L2 距離會因為 per-read 噪聲產生基礎 "floor"，降低真實生物信號與技術噪聲的比值

### 4.3 Read-Level Methylation 用於 Somatic Variant Characterization

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Baker et al. (ROCIT), bioRxiv, 2026** | 用 somatic mutation labels 訓練 read-level methylation classifier；perturbation analysis 顯示高 std dev CpGs（跨 cell type 變異大）最有資訊量 | ISM 可以借鑑 ROCIT 的特徵選擇思路：不是所有 CpG 都有相同資訊量 |
| **Long-read lung cancer study, Computers in Biology and Medicine, 2024** | ONT 長讀長同時偵測 somatic variants + methylation patterns in early lung cancer | 概念驗證：同一 read 同時提供 variant + methylation 資訊確實可行 |
| **O'Neill et al., Cell Genomics, 2024** | Advanced cancer cohort ONT 長讀長：resolve rearrangements + unravel haplotypes + reveal methylation landscapes | 臨床場景下的 comprehensive genomic + epigenomic profiling 可行性 |

---

## 5. 主題四：Subclonal Methylation Heterogeneity

### 5.1 腫瘤內甲基化異質性的模式

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Gabbutt & Duran-Ferrer et al. (EVOFLUx), Nature, 2025** | 1,976 淋巴腫瘤樣本：利用 DNA methylation fluctuation 推斷 subclonal architecture 和 evolutionary trajectory；subclonal selection 在 bulk 樣本中不常見（CLL ~30%，DLBCL <5%）；epimutation rate 跨疾病差異達數個數量級；evolutionary history 是獨立預後因子 | EVOFLUx 證明甲基化 heterogeneity 本身有生物學意義（進化歷史追蹤），但其用途是推斷 clonal evolution 而非區分 variant origin |
| **Tarazona et al., Nature Genetics, 2025** | ITMD 量化腫瘤內甲基化異質性；ITMD 與 SCNA ITH 顯著相關但與 mutation ITH 弱相關；發現甲基化可補償 oncogene 共擴增的 essential gene 過表達 | ITMD 與 mutation heterogeneity 弱相關的發現解釋了為何 ISM 的甲基化距離指標無法有效反映 somatic variant 資訊 |
| **Melanoma subclone long-read study, bioRxiv, 2025** | 23 個小鼠黑色素瘤 subclone 長讀長分析：DMR 數量與 SNV 數量正相關（r > 0.5）；lineage-specific methylation trajectory 與 aggressive phenotype 相關 | 在 **subclone-resolved** 分析中，甲基化與突變可以正相關 — 但這需要 single-cell-derived subclone 解析度，非 bulk tumor |
| **Liver cancer single-cell, Scientific Reports, 2025** | 肝癌 multinodular 型通常為 polyclonal；DNA methylation distance 比 gene expression distance 更適合重建腫瘤演化 trajectory | 支持甲基化作為 evolutionary distance metric 的有效性，但同樣不適用於 variant-level classification |

### 5.2 甲基化 Heterogeneity 是否可推斷 Tumor Purity / Subclone Fraction？

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **MethylPurify, Genome Biology, 2014** | 從單一腫瘤 bisulfite-seq 樣本估計 purity；利用 bimodal methylation distribution；正確推斷 purity 並識別 >96% DMR | 概念上與 ISM 的 read-level distance 分析相關，但 MethylPurify 依賴 site-level aggregation 而非 read-level pattern |
| **Jeong et al. (MethylBERT), Nature Communications, 2025** | Read-level methylation pattern 可準確估計 tumor cell fraction；Transformer 架構，150 bp reads 即可工作 | 最接近 ISM 的 read-level 分析概念，但需要已知的 normal methylation reference |
| **MONTE, bioRxiv, 2026** | Cancer-label-free purity estimation：從 bulk methylation 推斷 purity without cancer-specific reference；Bayesian transfer learning 可快速適應新 cancer type | 最新方法，不需要 cancer-specific reference — 概念上可為 ISM 提供 purity-aware 甲基化校正 |
| **PureBeta, NAR Genomics and Bioinformatics, 2024** | 單樣本 purity 估計框架：用全基因組甲基化數據估計 purity 並校正每個 CpG 的 beta value | Purity correction 可能幫助 ISM 減少 tumor-normal 混合比例帶來的甲基化混淆 |

### 5.3 Epigenetic Heterogeneity 與 Clonal Evolution 的關係

**核心洞見**（綜合 EVOFLUx + Tarazona + Melanoma subclone studies）：
- 甲基化異質性反映的是 **clonal age 和 division history**，而非特定 somatic mutation 的效應
- Epimutation rate（隨機甲基化變化率）與細胞分裂次數正相關
- ISM 觀察到的 inter-read methylation distance 主要反映的是 **tumor heterogeneity 背景噪聲**，而非 variant-specific 信號

---

## 6. 主題五：統計方法學

### 6.1 AUC 以外的 Biomarker 評估方法

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Richardson et al., Patterns (Cell), 2024** | **ROC-AUC 在 class imbalance 下是穩定的**（不受 prevalence 影響）；PR-AUC 隨 class imbalance 劇烈變化；PR-AUC 無法簡單 normalize 來校正；推薦 partial ROC-AUC (ROC-AUC_0.1) 關注 high-confidence 區域 | ISM 的 TP/FP 比例不均衡（大部分 variant 是 TP），但根據此研究，ROC-AUC < 0.58 的結論是 **真實的**（不是 class imbalance artifact）。Partial AUC 可能揭示 FP-enriched 子集中的弱信號 |
| **Saito & Rehmsmeier, PLoS One, 2015** | 在強 imbalanced 數據中，precision-recall plot 比 ROC plot 更 informative；PR 基線隨 positive rate 變化 | 早期主張 PR curve 更好，但被 Richardson 2024 挑戰 |
| **Cook & Ramadas, Stata Journal, 2020** | 實務指引：當 positive rate <5% 且主要關心 positive predictive value 時使用 PR curve | ISM 的 FP rate 在 TO 模式約 20-30%（非極度 imbalanced），ROC-AUC 仍是合理指標 |

### 6.2 Partial AUC 與 Subgroup Analysis

**對 ISM 的啟示**：
- **Partial ROC-AUC (FPR < 0.1)**: 可評估 ISM 特徵在 "high-confidence FP removal" 場景下的表現，但 ISM 先前研究（G1-G7）已在 TP loss <= 2% 約束下發現 FP removal = 0%
- **Subgroup analysis by cancer type / AF bin**: ISM 先前研究（O1-O13）已系統性地按 LOH/non-LOH、AF bin、sample 分層，仍無突破
- **Calibration**: 對 ISM 不太適用，因為核心問題是 discrimination 不足而非 calibration 不佳

### 6.3 高維度 -omics 數據中的 Confound Adjustment

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **SVA/SmartSVA, BMC Genomics, 2017** | Surrogate Variable Analysis 用於 EWAS 中的未知 confound 校正；SmartSVA 改進收斂和速度；可自動識別 cell composition 和 batch effect | ISM 已使用 residualization 處理 AF confound（O12 L2 collider bias 發現），但 SVA 框架可能提供更系統性的 confound detection |
| **McCartney et al., Clinical Epigenetics, 2019** | 批次效應偵測與校正（ComBat、FunNorm）；principal component partial R-squared 用於量化系統性變異來源 | ISM 跨 7 個 samples 的分析可能受 sample-level batch effect 影響；但 LOSO cross-validation 已部分控制此問題 |
| **Xu et al., Genome Biology, 2020** | Covariate-adaptive FDR control 在 EWAS 中的應用；利用輔助共變項（如 CpG density、distance to TSS）提升 FDR 偵測力 | 概念上可用於 ISM：利用 genomic context（CpG density, chromatin state）作為 covariate 可能改善特徵的 signal-to-noise |

### 6.4 Multiple Testing Correction 最佳實踐

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **Mansell et al., Epigenomics, 2019** | EWAS multiple testing: Bonferroni 太保守（450K array 有 450K tests）；FDR (BH procedure) 是標準；建議報告 effect size + P value + FDR q-value | ISM 的 per-region 分析（748K regions）面臨嚴重 multiple testing 問題；但核心問題是 effect 本身不存在，而非 power 不足 |

---

## 7. 補充主題：甲基化增強 Phasing 的最新進展

| 文獻 | 關鍵發現 | ISM 關聯性 |
|------|---------|-----------|
| **LongHap, bioRxiv, 2026** | 整合 sequence + methylation signals 進行 read-based phasing；動態識別 informative differentially methylated sites；outperforms WhatsHap, HapCUT2, LongPhase, MethPhaser (lower switch error, greater contiguity) | 直接相關於 ISM 的 self-phasing 問題 — 如果 methylation-informed phasing 能改善 HP tag 品質，可能間接改善 ISM 的 HP-dependent 特徵 |
| **HapBridge, bioRxiv, 2025** | Methylation-guided switch error correction and phase block bridging | 進一步證明甲基化資訊可改善 phasing 品質 |

---

## 8. Synthesis：文獻證據的整合判斷

### Q1: "Germline mQTL 的 ASM 效應 >> Somatic mutation 的 ASM 效應" 是否有文獻支持？

**強烈支持。** 這是本次文獻調查最明確的結論。

1. **量化比較**：
   - Germline: 23-37% het SNP 有 ASM (Shoemaker 2010)，delta >20%
   - Somatic: 僅 6-17% 的 cancer ASM gain 歸因於 somatic mutation (Do 2020)
   - **Germline ASM 效應約為 Somatic 的 3-6 倍**（以比例計）

2. **生物學機制**：
   - Germline: CpG-SNP（12% of mQTL，確定性效應）+ TF binding disruption（穩定、跨組織）
   - Somatic: Passenger SNV 無明確 cis-regulatory 機制；driver mutation 可影響全域但非局部

3. **ISM 觀察的解釋**：
   - ISM 發現 Paired 模式下 Pairwise distance AUC < 0.50（反轉）— 這正是因為 **germline (FP) 的 ASM 比 somatic (TP) 更強**，導致 FP 的 inter-haplotype methylation distance 更大
   - 文獻完美解釋了 ISM Fine-Pairwise distance 分析中 LOH 區域 Paired AUC=0.132（極端反轉）的發現

### Q2: "Read-level methylation pattern 用於 somatic variant classification" 是否有先例？

**無直接先例，但有相關概念驗證。**

1. **Read-level methylation classification 已證明可行**：
   - ROCIT (Baker 2026): AUC 0.933 for tumor/non-tumor read classification
   - MethylBERT (Jeong 2025): >95% accuracy for cell-of-origin classification

2. **但問題定義不同**：
   - ROCIT/MethylBERT 分類的是 **read 來源**（tumor cell vs normal cell），不是 **variant 性質**（somatic vs germline）
   - Normal cell 也攜帶 germline variant，tumor cell 也攜帶 germline variant
   - 直接將 "tumor read" 等同於 "somatic variant read" 是不正確的

3. **ISM 面臨的獨特挑戰**：
   - ISM 嘗試在 **同一 tumor cell population 內** 區分不同 variant 周圍的甲基化模式差異
   - 但 tumor cell 的甲基化背景是共享的（global hypomethylation + focal hypermethylation）
   - 這與 ROCIT 分析 tumor vs normal cells（甲基化背景不同）的場景根本不同

### Q3: "AUC 以外的評估方法可能揭示隱藏信號" 的方法學基礎是什麼？

**根據最新文獻，答案是否定的 — AUC < 0.58 的結論是穩健的。**

1. **Richardson et al. (2024) 明確證明**：ROC-AUC 不受 class imbalance 影響，ISM 的 AUC < 0.58 不是 imbalance artifact

2. **PR-AUC 不會提供額外信息**：
   - PR-AUC 對 class imbalance 敏感，但 ISM 的 TP/FP 比例並非極端（~70:30）
   - 在 ISM 的場景下，ROC-AUC 和 PR-AUC 的結論方向一致

3. **可能有價值的替代分析**：
   - **Partial ROC-AUC (FPR < 0.05)**: 評估極高 specificity 區域，但 ISM G1-G7 已在 TP loss <= 2% 下測試，FP removal = 0%
   - **Subgroup-specific AUC**: ISM 已系統性分層（LOH/non-LOH, AF bins, per-sample），均未突破
   - **Covariate-adaptive analysis**: SVA 框架可能發現未知 confounders，但 ISM 已 residualize AF, NumReads, sample 等主要 confounders

4. **真正的方法學啟示**：問題不在評估方法，而在 **信號本身不存在**（或太弱以至於被 ONT per-read 7% 錯誤率淹沒）

---

## 9. 資料來源評估

| 來源 | 類型 | 可信度 | 年份 | 核心貢獻 |
|------|------|--------|------|---------|
| Oliva et al., Nature Genetics | 論文 | 高 | 2023 | 最大規模多組織 mQTL mapping |
| Gaunt et al., Genome Biology | 論文 | 高 | 2016 | cis-mQTL 奠基研究 |
| Min et al., Nature Communications | 論文 | 高 | 2024 | 跨群體 mQTL 保守性 |
| Shoemaker et al., Genome Research | 論文 | 高 | 2010 | ASM + CpG-SNP 奠基研究 |
| Do et al., Genome Biology | 論文 | 高 | 2020 | 癌症 vs 正常 ASM 最大規模比較 |
| Baker et al. (ROCIT), bioRxiv | 預印本 | 中高 | 2026 | 最相關的 read-level methylation classification |
| Jeong et al. (MethylBERT), Nature Comms | 論文 | 高 | 2025 | Read-level methylation deconvolution |
| Tarazona et al., Nature Genetics | 論文 | 高 | 2025 | NSCLC ITMD + methylation-genomic cooperation |
| Gabbutt & Duran-Ferrer (EVOFLUx), Nature | 論文 | 高 | 2025 | Methylation-based subclonal architecture inference |
| Richardson et al., Patterns (Cell) | 論文 | 高 | 2024 | AUC vs PR-AUC under class imbalance |
| Comprehensive Benchmarking, bioRxiv | 預印本 | 中高 | 2024 | ONT methylation calling tool comparison |
| DeepBAM, Briefings in Bioinformatics | 論文 | 高 | 2024 | Per-read 5mC detection accuracy |
| Genner et al., Genome Research | 論文 | 高 | 2025 | ONT R9 vs R10 methylation comparison |
| Fu et al., Nature Reviews Genetics | 綜述 | 高 | 2025 | Long-read methylation computational review |
| LongHap, bioRxiv | 預印本 | 中高 | 2026 | Methylation-informed phasing |
| MONTE, bioRxiv | 預印本 | 中 | 2026 | Cancer-label-free purity estimation |
| IMAGE, Genome Biology | 論文 | 高 | 2019 | Haplotype-aware mQTL mapping |
| MethylPurify, Genome Biology | 論文 | 高 | 2014 | Single-sample methylation purity estimation |
| Huan et al., Scientific Reports | 論文 | 高 | 2019 | 55K replicated mQTL |
| Melanoma subclone study, bioRxiv | 預印本 | 中高 | 2025 | Single-cell-derived subclone methylation evolution |
| Wen et al., PLoS Comp Bio | 論文 | 高 | 2017 | Driver mutation -> methylation landscape |
| Haggerty et al., Genetics | 論文 | 高 | 2020 | CpG methylation affects neighboring mutability |
| Zhou et al., Nature Genetics | 論文 | 高 | 2018 | PMD hypomethylation + mutation density |
| SVA/SmartSVA, BMC Genomics | 論文 | 高 | 2017 | Reference-free confound adjustment in EWAS |
| Mansell et al., Epigenomics | 論文 | 高 | 2019 | Multiple testing correction in EWAS |
| Xu et al., Genome Biology | 論文 | 高 | 2020 | Covariate-adaptive FDR control |
| PureBeta, NAR Genomics | 論文 | 高 | 2024 | Single-sample purity estimation + beta correction |

---

## 10. 完整參考文獻列表

### mQTL 與 ASM

1. **Oliva M et al.** (2023). DNA methylation QTL mapping across diverse human tissues provides molecular links between genetic variation and complex traits. *Nature Genetics* 55:112-122. DOI: 10.1038/s41588-022-01248-z
2. **Gaunt TR et al.** (2016). Systematic identification of genetic influences on methylation across the human life course. *Genome Biology* 17:61. DOI: 10.1186/s13059-016-0926-z
3. **Min JL et al.** (2024). Genetic control of DNA methylation is largely shared across European and East Asian populations. *Nature Communications* 15:3363. DOI: 10.1038/s41467-024-47005-0
4. **Huan T et al.** (2019). Identification of 55,000 Replicated DNA Methylation QTL. *Scientific Reports* 9:18193. DOI: 10.1038/s41598-018-35871-w
5. **Shoemaker R et al.** (2010). Allele-specific methylation is prevalent and is contributed by CpG-SNPs in the human genome. *Genome Research* 20:883-889. DOI: 10.1101/gr.104695.109
6. **Do C et al.** (2020). Allele-specific DNA methylation is increased in cancers and its dense mapping in normal plus neoplastic cells increases the yield of disease-associated regulatory SNPs. *Genome Biology* 21:153. DOI: 10.1186/s13059-020-02059-3
7. **Zeng X et al.** (2019). IMAGE: high-powered detection of genetic effects on DNA methylation using integrated methylation QTL mapping and allele-specific analysis. *Genome Biology* 20:220. DOI: 10.1186/s13059-019-1813-1

### Somatic Mutation 與甲基化

8. **Tarazona N et al.** (2025). DNA methylation cooperates with genomic alterations during non-small cell lung cancer evolution. *Nature Genetics* 57:2226-2237. DOI: 10.1038/s41588-025-02307-x
9. **Wen B et al.** (2017). Significant associations between driver gene mutations and DNA methylation alterations across many cancer types. *PLoS Computational Biology* 13:e1005840. DOI: 10.1371/journal.pcbi.1005840
10. **Haggerty C et al.** (2020). Cytosine Methylation Affects the Mutability of Neighboring Nucleotides in Germline and Soma. *Genetics* 214:809-821. DOI: 10.1534/genetics.120.303028
11. **Zhou W et al.** (2018). DNA methylation loss in late-replicating domains is linked to mitotic cell division. *Nature Genetics* 50:591-602. DOI: 10.1038/s41588-018-0073-4

### Read-Level Methylation Classification

12. **Baker TM et al.** (2026). Genome-wide classification of tumor-derived reads from bulk long-read sequencing. *bioRxiv*. DOI: 10.64898/2026.03.03.709085
13. **Jeong Y et al.** (2025). MethylBERT enables read-level DNA methylation pattern identification and tumour deconvolution using a Transformer-based model. *Nature Communications* 16:849. DOI: 10.1038/s41467-025-55920-z

### ONT 甲基化精度

14. **Comprehensive Benchmarking** (2024). Comprehensive benchmarking of tools for nanopore-based detection of DNA methylation. *bioRxiv*. DOI: 10.1101/2024.11.09.622763
15. **Zhou Y et al.** (2024). DeepBAM: a high-accuracy single-molecule CpG methylation detection tool for Oxford nanopore sequencing. *Briefings in Bioinformatics* 25:bbae413. DOI: 10.1093/bib/bbae413
16. **Genner R et al.** (2025). Assessing DNA methylation detection for primary human tissue using Nanopore sequencing. *Genome Research* 35:gr.279159.124. DOI: 10.1101/gr.279159.124
17. **Fu Y, Timp W, Sedlazeck FJ** (2025). Computational analysis of DNA methylation from long-read sequencing. *Nature Reviews Genetics* 26:620-634. DOI: 10.1038/s41576-025-00822-5

### Subclonal Methylation Heterogeneity

18. **Gabbutt C, Duran-Ferrer M et al.** (2025). Fluctuating DNA methylation tracks cancer evolution at clinical scale. *Nature*. DOI: 10.1038/s41586-025-09374-4
19. **Melanoma subclone study** (2025). Long-read sequencing of single cell-derived melanoma subclones reveals divergent and parallel genomic and epigenomic evolutionary trajectories. *bioRxiv*. DOI: 10.1101/2025.08.28.672865
20. **Liver cancer single-cell** (2025). Interrogating subclonal heterogeneity of liver cancer with single-cell multi-omics analysis. *Scientific Reports*. DOI: 10.1038/s41598-025-24732-y

### Tumor Purity Estimation

21. **Zheng X et al.** (2014). MethylPurify: tumor purity deconvolution and differential methylation detection from single tumor DNA methylomes. *Genome Biology* 15:419. DOI: 10.1186/s13059-014-0419-x
22. **MONTE** (2026). MONTE: Methylation-based Observation Normalization and Tumor purity Estimation. *bioRxiv*. DOI: 10.64898/2026.01.22.701164
23. **PureBeta** (2024). Tumor purity estimated from bulk DNA methylation can be used for adjusting beta values of individual samples to better reflect tumor biology. *NAR Genomics and Bioinformatics* 6:lqae146. DOI: 10.1093/nargab/lqae146

### 統計方法學

24. **Richardson E et al.** (2024). The receiver operating characteristic curve accurately assesses imbalanced datasets. *Patterns (Cell)* 5:100994. DOI: 10.1016/j.patter.2024.100994
25. **Saito T, Rehmsmeier M** (2015). The Precision-Recall Plot Is More Informative than the ROC Plot When Evaluating Binary Classifiers on Imbalanced Datasets. *PLoS One* 10:e0118432. DOI: 10.1371/journal.pone.0118432
26. **Chen Y et al.** (2017). Fast and robust adjustment of cell mixtures in epigenome-wide association studies with SmartSVA. *BMC Genomics* 18:413. DOI: 10.1186/s12864-017-3808-1
27. **Mansell G et al.** (2019). Correction for Multiple Testing in Candidate-Gene Methylation Studies. *Epigenomics* 11:1009-1021. DOI: 10.2217/epi-2018-0204
28. **Xu Z et al.** (2020). Leveraging biological and statistical covariates improves the detection power in epigenome-wide association testing. *Genome Biology* 21:88. DOI: 10.1186/s13059-020-02001-7

### Methylation-Informed Phasing

29. **LongHap** (2026). Harnessing methylation signals inherent in long-read sequencing data for improved variant phasing. *bioRxiv*. DOI: 10.64898/2026.03.11.710820
30. **MethPhaser** (2024). Methylation-based long-read haplotype phasing of human genomes. *Nature Communications* 15:5327. DOI: 10.1038/s41467-024-49588-0

---

## 11. 對 ISM 後續研究的建議

基於本次文獻調查，以下建議按優先序排列：

### 高優先

1. **Phase 2A Normal Methylation Reference 仍是正確方向**：文獻明確支持 tumor/normal 甲基化背景差異是可區分的（ROCIT AUC 0.933），ISM 的 normal reference baseline 策略與此一致

2. **ISM 定位為 Somatic Heterogeneity Characterization 工具（而非 Variant Filter）是正確的**：EVOFLUx 和 ITMD 的成功案例表明甲基化 heterogeneity 在 evolutionary trajectory inference 中有明確價值

3. **CpG 選點策略值得探索**：ROCIT perturbation analysis 發現高 cross-cell-type std dev CpGs 最有資訊量；ISM 可考慮類似的 informative CpG selection

### 中優先

4. **Purity-aware methylation correction**：MONTE/PureBeta 的 purity correction 概念可整合到 ISM 中，減少 tumor-normal 混合比例的甲基化混淆

5. **Methylation-informed phasing 改善 HP tag 品質**：LongHap 證明甲基化信號可改善 phasing，可能間接改善 ISM 的 HP-dependent 特徵（但 ISM 已確認所有信號來自 HP tags）

### 低優先

6. **替代評估指標**：Richardson 2024 確認 ROC-AUC 是穩健的，ISM 的 AUC < 0.58 結論不需要用其他指標重新評估
