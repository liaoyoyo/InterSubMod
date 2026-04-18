# 08. 文獻調查與理論基礎

<!--
建立時間: 2026-04-01 00:00
目標: 彙整 LOH-甲基化關係與 germline/somatic 區分的文獻調查，連結本週觀察結論
處理範圍: 兩份文獻調查報告的核心發現整理、與本週實驗觀察的對照分析
關聯檔案:
  - docs/references/20260331_LOH區域甲基化模式與事件順序文獻調查_01.md
  - docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md
  - docs/experiments/INDEX.md
-->

---

## 1. LOH 與甲基化的生物學關係

### 1.1 LOH 區域的甲基化特性（文獻共識）

LOH（Loss of Heterozygosity）消除基因組中某一 allele 的貢獻，在甲基化層面有以下已確認的效應：

| 特性 | 文獻依據 | 關鍵數據 |
|------|---------|---------|
| LOH 後甲基化趨向單一 allele pattern | CAMDAC / TRACERx (Tarazona et al., Nat Genet, 2025) | 驗證相關性 Pearson r = 0.97 |
| LOH 消除一個 allele 的表達與甲基化貢獻 | Ha et al. (Genome Res, 2012) | 大部分 monoallelic expression 可由 LOH 解釋 |
| LOH 區域可破壞 imprinted ASM | Zink et al. (Nat Commun, 2025) | 鑑定 460 個 parental allele-specific methylation 區域 |
| LOH 與 promoter methylation 可互補失活 TSG | Jones & Baylin (Cancer Res, 2022) | Knudson 二次打擊的表觀遺傳延伸 |

**核心共識**：LOH 區域的腫瘤甲基化反映單一 allele 的 epigenetic snapshot，這已被 CAMDAC 在 122 個 TRACERx 多區域樣本中驗證（r = 0.97）。

### 1.2 三場景理論

LOH 區域中 somatic SNV 的出現可分為三種事件順序場景，各有不���的甲基化預測：

#### 場景 1: LOH 先發生 -> sSNV 後出現

- **機制**：LOH 消除一個 haplotype，所有剩餘 reads 來自同一 haplotype；後續 sSNV 影響微弱
- **文獻支持**：CAMDAC 驗證 LOH 區域確實趨向單 allele pattern (r = 0.97)；de novo ASM（somatic mutation 引起）僅佔 ASM 增加的 6-17% (Do et al., Genome Biol, 2020)
- **ISM 預測**：HP 之間甲基化差異小、within-HP concordance 高
- **驗證狀態**：假說大致成立

#### 場景 2: sSNV 先發生 -> LOH 後出現

- **機制**：sSNV 先出現在一個 haplotype，LOH 隨後消除另一個 haplotype；結果取決於哪個 allele 被保留
- **文獻支持**：PCAWG timing framework (Gerstung et al., Nature, 2020) 提供事件順序推論方法；duplicated vs non-duplicated mutations 的比例可估計 molecular time
- **ISM 預測**：HP 甲基化差異取決於殘留 allelic ratio，sSNV 本身的甲基化效應通常微弱
- **驗證狀態**：假說部分成立，複雜度高

#### 場景 3: Germline variant 被 TO caller 誤判為 somatic

- **機制**：germline het SNP 在 LOH 區域 VAF 偏離 50%，TO caller 將其誤判為 somatic
- **文獻支持**：Sun et al. (JCO Precis Oncol, 2022) 證實 LOH 區域 germline het VAF 可升高至 60-100%；LOHGIC (Khiabanian et al., J Mol Diagn, 2018) 建立了 VAF 模型但承認 ambiguous cases；76% 攜帶 germline mutation 的腫瘤存在 LOH
- **ISM 預測**：需修正——部分 LOH 時 HP 差異大（反映 ASM），完全 LOH 時差異消失
- **驗證狀態**：假說需要 LOH 完整性校正

### 1.3 事件順序（temporal order）為什麼重要

1. **VAF 分佈不同**：LOH 前 vs 後的 sSNV 有不同的 allelic copies，影響 caller 的 VAF 判讀 (Gerstung et al., 2020)
2. **甲基化足跡不同**：LOH-first 場景的甲基化較均質，sSNV-first 可能保留突變 allele 的 altered pattern
3. **ISM 特徵方向不同**：三場景的 HP methylation divergence 方向和幅度各異，混為一談會抵消信號
4. **甲基化可作為獨立計時器**：EVOFLUx (Gabbutt et al., Nature, 2025) 證明 fluctuating CpGs 可作為分子鐘，與 VAF-based timing 互補

---

## 2. 甲基化區分 Germline vs Somatic 的文獻

### 2.1 mQTL (methylation Quantitative Trait Loci) 效應

mQTL 是甲基化區分 germline/somatic 最有力的生物學基礎：

| 發現 | 來源 | 關鍵數據 |
|------|------|---------|
| 大規模 cis-mQTL 存在 | Oliva et al. (Nat Genet, 2023) | 4.7M cis-mQTL variant-CpG pairs，63 萬 trans-mQTL |
| 93% mQTL 為 cis 作用 | Gaunt et al. (Genome Biol, 2016) | 距離越近效應越強 |
| 12% SNP 直接落在 CpG 位點 | Gaunt et al. (2016) | CpG-SNP 徹底破壞甲基化能力 |
| 38-88% ASM 依賴 CpG-SNP | Shoemaker et al. (Genome Res, 2010) | 23-37% het SNPs 表現 ASM |
| De novo ASM (somatic) 僅 6-17% | Do et al. (Genome Biol, 2020) | 遠少於 germline-driven ASM |

**關鍵推論**：Germline heterozygous variant 傾向產生**不對稱但穩定**的甲基化模式（ASM），這是遺傳決定的、跨組織穩定的表觀遺傳特徵。Somatic variant 產生的 de novo ASM 數量級弱於 germline。

### 2.2 現有工具

| 工具 | 核心方法 | 與本研究關聯 | 局限性 |
|------|---------|-------------|--------|
| **ROCIT** (Baker et al., bioRxiv, 2026) | Transformer 對每條 read 做 tumor/non-tumor 分類，利用 somatic mutations 做 training labels | 證明 read-level methylation 足以區分 cell origin | 分類 cell origin 非 variant origin；需已知 somatic mutations 做 label |
| **MethylBERT** (Jeong et al., Nat Commun, 2025) | 修改版 BERT，150bp reads 準確率 >95% | 證明 read-level methylation 有足夠信息量 | 專注表觀遺傳分類，不直接處理 variant classification |
| **NanoMethPhase** (Akbari et al., Genome Biol, 2021) | ONT long reads megabase-scale methylation phasing | ONT 平台偵測 ASM 的技術基礎 | 非分類工具 |
| **MethPhaser** (Nat Commun, 2024) | 甲基化增強 phasing，N50 +78-151% | 改善 phasing 可間接幫助區分 | 非分類工具 |
| **MethylPurify** (Zheng et al., Genome Biol, 2014) | 從單一腫瘤 methylome 估計 purity | 證明甲基化可獨立分離 tumor/normal | 非 variant-level 工具 |

**重要發現**：目前**沒有任何工具**直接利用甲基化作為 tumor-only pipeline 中的 germline filter。

### 2.3 文獻的主要結論

**支持可行性的方向：**
1. ONT 同時提供 variant + methylation 的 read-level 資訊，是理想平台
2. mQTL/ASM 生物學基礎強大
3. ROCIT/MethylBERT 證明 read-level patterns 有足夠信息量

**質疑可行性的方向：**
1. 癌症甲基化 landscape 劇烈改變（全域低甲基化 + 局部高甲基化），可能淹沒 mQTL signal
2. 癌症中 ASM 增加 5-9 倍 (Do et al., 2020)，增加的 ASM 不一定與 germline 相關
3. CpG density 分佈不均，很多 variant 周圍 CpG 不足
4. 癌症中 allele switching 比正常組織多 3 倍（43% vs 14%）

---

## 3. 文獻如何支持/不支持本週觀察

### 3.1 TO LOH 系統性過判

**本週觀察**：TO 同位點 LOH 判定率為 paired 的 16-52 倍

**文獻對應**：
- Sun et al. (2022)：LOH 區域 germline het VAF 升高至 60-100%，導致 TO caller 誤判
- LOHGIC (2018)：LOH 存在於 76% 攜帶 germline mutation 的腫瘤中
- ClairS-TO (2025)：即使有三層過濾，LOH 區域 germline het variant 仍可能逃過 PON 過濾

**文獻支持度**：**強支持**。TO 場景下 LOH 過判是已知問題，文獻中多次報導。本週觀察的 16-52 倍過判率與文獻描述的機制一致：LongPhase-TO 缺乏 normal BAM anchor，phasing 依賴 germline het SNPs，LOH 使 VAF 偏移導致 phasing artifact。

### 3.2 甲基化區分力弱（TO 所有 AUC < 0.58）

**本週觀察**：TO 模式下所有 ISM 甲基化特徵 AUC < 0.58

**文獻對應**：
- 癌症中 ASM 增加 5-9 倍 (Do et al., 2020)，增加的 stochastic ASM 可掩蓋 germline ASM signal
- PMD 區域的 stochastic drift (Zhou et al., 2018; Johnstone et al., 2019) 增加甲基化雜訊
- 癌症 allele switching 3 倍增加 (Do et al., 2020) 破壞 ASM 方向穩定性

**文獻支持度**：**部分支持**。文獻中的工具（ROCIT, MethylBERT）在特定條件下有效（充足 CpGs、已知 training labels、read-level 而非 region-level），但這些條件在 ISM 的 TO 場景下並不滿足。文獻工具有效的前提包括：
- ROCIT：需要已知 somatic mutations 做 training（TO 場景的 chicken-and-egg 問題）
- MethylBERT：在 150bp reads 上準確率 >95%，但依賴 training atlas

**結論**：ISM 目前的 region-level 甲基化特徵在 TO 下 AUC < 0.58 不令人意外。文獻支持的是 **read-level** 分類可行性，而非 region-level aggregate features 的區分力。這指向 InterSubMod 需要從 region-level 轉向 read-level 特徵工程的研究方向（Direction E）。

### 3.3 O12 LOH 三場景不可區分

**本週觀察**：O12 發現 TO LOH 三場景甲基化區分失敗（AF confound + L2 collider bias）；所有 corrected AUC < 0.58

**文獻對應**：
- 三場景的甲基化差異理論上存在（文獻支持），但混淆因素眾多：
  - PMD stochastic drift 使 LOH-first 和 non-LOH 場景看起來相似
  - 癌症全域低甲基化掩蓋 sSNV 微弱效應
  - 癌症增加的 stochastic ASM 使部分 somatic ASM 也穩定
- **目前沒有文獻直接驗證 ISM-level 三場景可區分性**

**文獻支持度**：**不矛盾但也不支持**。文獻中三場景的甲基化差異是理論推導，尚無實證。O12 的 negative result 是首次在 ISM 框架下的實證檢驗，其 negative 結論與文獻提出的大量混淆因素一致。

### 3.4 O11 甲基化 Heterogeneity 假說否決

**本週觀察**：Within-group heterogeneity 無法區分 TP/FP（epipolymorphism AUC 0.845 -> 0.530 after n_reads correction）

**文獻對應**：
- Do et al. (2020)：癌症 ASM 增加 5-9 倍，主要來源是 allele-specific loss of methylation（72-76%）
- Brocks et al. (2024)：epigenomic heterogeneity 在 enhancer 和 PMD 區域最大，但這不是 germline/somatic-specific 的

**文獻支持度**：**間接支持 negative result**。Within-group heterogeneity 的 high AUC 源自 n_reads confound（read 數越多，距離計算越穩定），而非生物學信號。文獻中的甲基化 heterogeneity 主要反映 tumor heterogeneity 程度，而非 variant origin，這與 O11 的 negative 結論一致。

---

## 4. 關鍵文獻摘要表

| 文獻 | 年份 | 核心貢獻 | 與本週觀察的連結 |
|------|------|---------|---------------|
| CAMDAC / TRACERx (Tarazona et al.) | 2025 | LOH 區域甲基化驗證金標準 (r=0.97) | LOH 甲基化趨同理論基礎 |
| PCAWG timing (Gerstung et al.) | 2020 | 突變與 CN 事件順序推論框架 | 三場景的分子鐘理論 |
| Do et al. | 2020 | 癌症 ASM 5-9x，de novo ASM 僅 6-17% | TO 甲基化區分力弱的解釋 |
| Shoemaker et al. | 2010 | 23-37% het SNPs 表現 ASM | Germline ASM 的基礎數據 |
| Oliva et al. | 2023 | 4.7M cis-mQTL pairs | mQTL 普遍性支持 |
| ROCIT (Baker et al.) | 2026 | Read-level methylation 可區分 tumor/normal reads | Direction E 的外部驗證 |
| MethylBERT (Jeong et al.) | 2025 | Read-level methylation 有足夠信息量 | Read-level 特徵可行性 |
| LOHGIC (Khiabanian et al.) | 2018 | LOH 區域 germline/somatic VAF 模型 | TO LOH 過判的理論解釋 |
| Sun et al. | 2022 | LOH 區域 germline VAF 升高驗證 | TO LOH 過判的實證支持 |
| EVOFLUx (Gabbutt et al.) | 2025 | 甲基化作為分子鐘 | 事件順序推論的新方向 |

---

## 待驗證問題（全為未來研究方向 — Direction E/A）

1. **ROCIT 在 ONT HCC1395 上的表現**：是否能直接應用於 ISM read-level 分類？定位 Direction E P1。
2. **mQTL signal 在高純度腫瘤中的可測性**：需在已知 germline/somatic variants 上實測。定位 Direction A P2。
3. **CpG density 分層效果**：高 CpG 區域甲基化區分力是否提升？可用現有數據分析但優先級低。定位 P2。
4. **MHB concordance 預測力**：Guo et al. 2017 框架在 ISM 中的適用性。定位 Direction E P2。

## 認知門檻補充建議

- **ASM (Allele-Specific Methylation)**：同一位點的兩個 allele 具有不同的甲基化狀態。正常組織中由 germline SNP 通過 mQTL 驅動，癌症中額外增加 stochastic 成分。
- **mQTL (methylation QTL)**：germline genetic variant 對附近 CpG 甲基化的 cis-regulatory 效應。93% 為 cis 作用（< 1 Mb），是 germline variant 產生穩定 ASM 的主要機制。
- **CpG-SNP**：SNP 直接落在 CpG 二核苷酸中，徹底破壞該位點的甲基化能力。佔 mQTL 的 ~12%，但驅動 38-88% 的 ASM。
- **PMD (Partially Methylated Domain)**：癌症中甲基化喪失的大尺度區域（佔基因組 24-63%），甲基化 stochastic drift 的主要來源。
- **CAMDAC**：TRACERx 團隊開發的癌症 allele-specific methylation 解卷積工具，利用 LOH 區域做驗證，是目前 LOH 甲基化分析的標準方法。
- **Direction E**：InterSubMod 2026-Q2 研究策略中的 "Read-level epigenetic context" 方向，受 ROCIT/MethylBERT 啟發，目標是從 region-level aggregate 轉向 read-level 甲基化特徵。
