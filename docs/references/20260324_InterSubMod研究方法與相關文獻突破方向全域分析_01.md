<!--
建立時間: 2026-03-24 15:00
目標: InterSubMod 研究方法、相關文獻、突破方向的全域分析報告
處理範圍: ISM 核心方法 + 12 篇外部研究 + 學界共識/分歧 + 5 條突破方向 + 漸進策略
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
  - docs/architecture/deep-research-report.md
  - docs/methodology/20260324_方法學審查全域結論報告_01.md
-->

# InterSubMod 研究方法、相關文獻、突破方向：全域分析報告

> **建立日期**：2026-03-24
> **分析範圍**：ISM 2025-11 至 2026-03 全實驗歷程 + 12 篇核心外部研究
> **研究者確認策略**：E → A+D → B → C（漸進式方法學升級）

## Context

InterSubMod 經過數月密集實驗（2025-11 至 2026-03），已建立完整的 SNV-anchored 甲基化聚類框架並完成方法學審查（16 份觀察文件）。核心發現是：**甲基化訊號作為 TP/FP 變異過濾器的鑑別力有限**（AUROC~0.5-0.6），但作為 annotation/support/ranking 仍有價值。目前最穩定的架構是 `caller-first → methylation-support → artifact triage`。

本報告旨在：(1) 系統整理 InterSubMod 的研究方法與主軸；(2) 深度對照 8+ 篇核心外部研究；(3) 釐清學界共識/分歧/未驗證之處；(4) 推導 InterSubMod 可行的突破方向。

---

## 一、InterSubMod 的研究方法與主軸

### 1.1 核心假設鏈

InterSubMod 的中心假設是：

> 在每個 somatic SNV 的局部視窗（預設 ±5000bp）內，read-level 甲基化圖樣應能形成穩定群集；若該 SNV 是真實腫瘤事件或與亞克隆/等位機制相關，則群集應與 ALT/REF、HP1/HP2、Tumor/Normal 標籤呈統計顯著關聯。

此假設鏈可拆為三個子問題：
- **變異真偽辨識**（TP/FP discrimination）
- **等位/單倍型特異性甲基化**（ASM / second hit）
- **亞克隆表觀遺傳結構**（subclonal epigenetic architecture）

### 1.2 技術架構

| 層 | 方法 | 實作 |
|---|---|---|
| 甲基化解析 | ONT BAM MM/ML 標籤解碼 + CIGAR 座標校正 | `MethylationParser.cpp` |
| 矩陣建構 | read×CpG 機率矩陣 + 二值矩陣 | `RegionProcessor.cpp` |
| 距離計算 | NHD/L1/L2/Bernoulli/Jaccard/CORR（6 種） | `DistanceCalculator.cpp` |
| 聚類 | UPGMA 階層式聚類 | 內建 |
| 顯著性檢定 | Fisher-Freeman-Halton + Cramér's V（cluster-first）；PERMANOVA（label-first） | `SignificanceAnalyzer` |
| 雙向驗證 | cluster-first ↔ label-first 一致性 → VerificationClass（Strong/Weak/Subclone/Noise） | 多層框架 |
| 品質分層 | HP ratio → LOH detection；coverage multiple → CNV；QualityScore (0-100) | LOH/CNV 分層 |

### 1.3 實驗結論總結（截至 2026-03-24）

**已確認的正面發現：**
- VerificationClass 框架整體合理：Strong Precision 81-100%，Weak 90-99%
- LOH_Subclone 是可靠的 Subclone 子類（HCC1395/HCC1937 佔 87-89%）
- `VC != Noise` 替代 `Significant=True` 可將捕獲率從 1-12% 提升至 33-68%
- TO 下 GQ>=3 可達 delta F1=+0.0094（HCC1395 5kHz pilot）
- QualityScore>=50 可作為 TO rescue support（Precision 76.9%）

**已確認的瓶頸：**
- ISM 對 TP/FP 的鑑別力不足（AUC~0.5-0.62）
- 甲基特徵在固定 caller gate 後**未超過** caller-only 的 GQ 規則
- FP 主要來自 germline ASM，ISM gate 無法區分體細胞 vs 胚系 ASM
- PairwiseMedianDist 方向具樣本依賴性，不可跨樣本全域化
- paired_persistent_final_fp（45 個跨平台共享）確認為 irreducible FP
- ISM 覆蓋 FN pool 僅 7%（773/11051），rescue 上限 F1 delta ≈+0.0097

**已 Rejected 的方向：**
- AlleleDelta 作為 standalone filter（跨樣本 AUROC 0.41-0.76）
- LOH 整合 VerificationClass（TP/FP LOH 率相似）
- HP_Coverage_Symmetry（AUROC ≈ 0.53）
- PERMANOVA 作為 filter（全基因組 AUROC 0.53-0.59，僅作 annotation）

---

## 二、核心外部研究對照分析

### 2.1 Long-Read POG（Cell Genomics 2024）

**論文**：O'Neill et al., "Long-read sequencing of an advanced cancer cohort resolves rearrangements, unravels haplotypes, and reveals methylation landscapes," Cell Genomics 4(11), 2024.

**核心方法**：
- 189 例晚期癌症 + 41 matched normal，ONT PromethION 定序
- 長距離相位（LongPhase / WhatsHap）確認 TSG biallelic inactivation
- NanoMethPhase 偵測 allelically differentially methylated regions (aDMRs)
- MBASED 做 ASE 的 beta-binomial 檢定
- IMPALA 整合 CNV + allelic methylation + somatic calls

**主要發現**：
- 在 copy-number balanced 區域，promoter allelic methylation 比 balanced ASE 更常見（mean 0.16 vs 0.01）
- Allelic promoter methylation 常與 major expressed allele 呈 trans 關係（被甲基化的等位基因表現較低）
- 多數雙突變為 trans，也有 cis；BRCA1/RAD51C promoter methylation + LOH → HRD
- CNV/LOH 是 ASE 的主要驅動因素

**與 InterSubMod 的關聯**：
- **直接可借鑑**：POG 以「基因/機制」為單位整合多種訊號（CNV + aDMR + ASE + somatic），而 ISM 以「SNV 視窗」為單位。ISM 可在上層加 gene-level evidence integrator
- **關鍵差異**：POG 不做 read-level clustering，而是用統計量（MBASED）做等位層級推論
- **啟示**：ISM 需要顯式引入 CNA/LOH 作為解釋變項，而非事後品質扣分

### 2.2 TRACERx NSCLC（Nature Genetics 2025）

**論文**：TRACERx consortium, "DNA methylation cooperates with genomic alterations during non-small cell lung cancer evolution," Nat Genet, 2025.

**核心方法**：
- 217 tumor samples / 59 NSCLC patients，RRBS + CAMDAC deconvolution
- **CAMDAC**：copy-number aware methylation deconvolution，校正 purity + CN 對 bulk methylation 的影響
- **ITMD**（intratumoral methylation distance）：pairwise Pearson distance 量化甲基化異質性
- **PDR 校正**：用 LOH-SNV 做天然 tumor/normal 分離，驗證 CAMDAC PDR（R>0.8）
- **M_R/M_N ratio**：regulatory vs nonregulatory CpG 高甲基化比，類似 dN/dS，偵測正向選汰

**主要發現**：
- 腫瘤 inter-patient 甲基化異質性為 normal 的 25 倍
- ITMD 與 SCNA-ITH 相關（r=0.47-0.66），但與 SNV-ITH 相關性弱
- Promoter 甲基化異質性最低（調控更緊密）
- 61/68 canonical TSGs 出現 CN loss 或 hypermethylation；19 例顯示平行演化
- **AllChAT 機制**：oncogene amplification 導致鄰近 essential gene 被高甲基化以維持蛋白化學計量

**與 InterSubMod 的關聯**：
- **最可直接搬用的驗證框架**：LOH-SNV 驗證 PDR 的方法。ISM 已有 per-read ALT/REF 標記，可直接升級為估計腫瘤/正常甲基化背景的工具
- **M_R/M_N 框架**：ISM 若要從「聚類顯著」走到「驅動機制推論」，需區分 regulatory/nonregulatory CpG
- **CAMDAC 的啟示**：ISM 目前在混合訊號上直接聚類，應先做 purity/CN deconvolution
- **ITMD 可參考**：ISM 的 PairwiseMedianDist 類似 ITMD，但未做 purity 校正

### 2.3 EVOFLUx（Nature 2025）

**論文**：Gabbutt & Duran-Ferrer et al., "Fluctuating DNA methylation tracks cancer evolution at clinical scale," Nature, 2025.

**核心方法**：
- 識別 978 個 fluctuating CpGs (fCpGs)：高異質性、近似中間甲基化、近似獨立波動
- beta-mixture model 將 fCpG methylation 離散化為 0/1/2（0%, 50%, 100%）
- 用分佈形狀推回族群演化史（growth rate, malignancy age, epimutation rate）
- Nanopore 長讀序驗證 fCpG 變動非由底層 somatic mutation 造成
- 應用於 1,976 例淋巴腫瘤

**主要發現**：
- fCpGs 可作為「演化條碼」，low-cost 追蹤腫瘤演化
- 初始腫瘤生長速率、惡性年齡、epimutation rate 跨疾病類型差異達數量級
- Subclonal selection 在 bulk 樣本中並不頻繁
- 更快的初始生長見於更具侵略性的亞型；演化歷史是獨立預後因子

**與 InterSubMod 的關聯**：
- **CpG 選擇策略的啟示**：ISM 目前把視窗內所有 CpG 平等看待。若目標是「亞克隆/演化」，應優先挑選 fCpG-like 的高異質性 CpG 集合
- **互補觀點**：EVOFLUx 用 bulk 甲基化推演化，ISM 用 read-level 甲基化看局部結構。兩者結合可能有價值
- **方法學差異**：EVOFLUx 不需要 somatic variant，純粹依賴甲基化波動；ISM 以 variant 為錨點

### 2.4 DeepSomatic（bioRxiv 2024）

**論文**：Park et al., "DeepSomatic: Accurate somatic small variant discovery for multiple sequencing technologies," bioRxiv, 2024.

**核心方法**：
- 改造 DeepVariant 的三階段架構（make_examples → call_variants → postprocess）
- 6 通道 tensor 表示 tumor-normal 讀序特徵
- 支援 Illumina/PacBio/ONT
- 跨 5 個 matched tumor-normal cell line（含 HCC1395），304,663 高信心 somatic variants

**主要發現**：
- Illumina SNV F1=0.9829，PacBio F1=0.9536，ONT 次之但持續改善
- Multi-cancer model 比 HCC1395-only 減少 38.7% FP/FN（SNV）
- Tumor-only 模式使用 population AF + PON，precision 低於 paired

**與 InterSubMod 的關聯**：
- ISM 的上游 caller（ClairS/ClairS-TO）競爭者。DeepSomatic 不使用甲基化資訊
- ISM 確認的 `caller-first` 結論在 DeepSomatic 語境下同樣成立
- **機會**：ISM 的甲基化 annotation 可作為 DeepSomatic/ClairS 的後處理 supplement

### 2.5 LongPhase-S（bioRxiv 2025）

**論文**：Ho et al., "LongPhase-S: Somatic haplotyping for long-read cancer sequencing," bioRxiv, 2025.

**核心方法**：
- 將每條 somatic read 錨定到 parental germline lineage
- Phase-aware purity estimation
- Variant recalibration 用 purity 資訊
- 處理 tumor heterogeneity / aneuploidy / contamination

**主要發現**：
- ClairS SNV F1 改善 +4.5%，DeepSomatic +1.2%
- Indel F1 改善更大（ClairS +7.1%）

**與 InterSubMod 的關聯**：
- LongPhase-S 是 ISM 目前 pipeline 的上游組件
- ISM 實驗確認 LongPhase-S phasing 沒有改變 TO benchmark call set
- **關鍵啟示**：somatic haplotyping 的 purity estimation 可直接輸入 ISM 作為 purity-aware methylation 的先驗

### 2.6 t-nanoEM（Cell Reports Methods 2025）

**論文**：Kunigo et al., "Targeted long-read methylation analysis using hybridization capture suitable for clinical specimens," Cell Reports Methods, 2025.

**核心方法**：
- Modified EM-seq + hybridization capture + nanopore PromethION
- 8-50ng DNA input，570x mean bait coverage
- 自建 WhatsHap-based phasing pipeline（處理 base-conversion 問題）
- Somatic variant-linked methylation detection

**主要發現**：
- 466-297 allele-biased methylated regions per cell line
- 108 DMRs 與 SNV-bearing alleles 關聯
- 乳癌/肺腺癌多區域微切割顯示 tumor-specific methylation patterns
- Loss of imprinting at MEST/NAP1L5

**與 InterSubMod 的關聯**：
- **臨床落地路線參考**：targeted panel + 高深度是 second hit / promoter methylation / LOH 問題的可行路線
- **方法學差異**：t-nanoEM 用 converted reads，ISM 用 native MM/ML。ISM 的優勢在不需要轉換
- **共同目標**：somatic variant-linked methylation = ISM 的核心概念

### 2.7 NanoMethPhase（Genome Biology 2021）

**論文**：Akbari et al., "Megabase-scale methylation phasing using nanopore long reads and NanoMethPhase," Genome Biology 22, 2021.

**核心方法**：
- SNV detection (Clair) → haplotype phasing (WhatsHap) → methylation phasing (NanoMethPhase)
- SNVoter 改善 nanopore SNV 精確度（precision 74.6% → 89%）
- ~10x 覆蓋即可有效偵測 ASM
- 93% phased/trio DMR 重疊

**與 InterSubMod 的關聯**：
- NanoMethPhase 是 ISM 參考的同領域工具，但定位不同：
  - NanoMethPhase: 全基因組 ASM 偵測（以 haplotype 為軸）
  - ISM: per-SNV 局部甲基化結構（以 variant 為錨）
- **互補空間**：NanoMethPhase 提供全域 ASM landscape，ISM 提供變異位點特異的深度分析

### 2.8 MethSig（Cancer Discovery 2021）

**論文**：Landau et al., "Discovery of candidate DNA methylation cancer driver genes," Cancer Discovery 11(9), 2021.

**核心方法**：
- Beta regression model 校正背景隨機高甲基化率
- Patient-specific p-value → Wilkinson combination → cohort-level significance
- RRBS + Infinium 450K array

**主要發現**：
- CLL 中識別 189 個候選甲基化驅動基因
- CRISPR knockout 驗證 DUSP22/RPRM/SASH1 conferring fitness advantage
- Pan-cancer：SOX17（13 cancer types）、RASSF1（6 types）

**與 InterSubMod 的關聯**：
- MethSig 定義「甲基化驅動」的統計框架，可啟發 ISM 從「顯著性」走向「功能性判讀」
- **機會**：ISM 的 per-SNV 甲基化結構若能對應到 MethSig 已知 driver gene，可增加生物意義

### 2.9 其他重要相關研究

| 研究 | 年份 | 核心貢獻 | 與 ISM 關聯 |
|------|------|---------|-------------|
| **MethPhaser** (Nat Commun 2024) | 2024 | 用甲基化訊號延伸 SNV-based phasing，N50 增加 78-151% | 可作為 ISM 的 phasing 補強 |
| **MethylBERT** (Nat Commun 2025) | 2025 | Transformer-based read-level 甲基化 pattern recognition + tumor deconvolution | 直接競爭 ISM 的 read-level 分析概念 |
| **PRISM** (Bioinformatics 2019) | 2019 | Reference-free subclonal inference from methylation patterns（beta-binomial mixture） | 最接近 ISM subclone 目標的先行研究 |
| **MHB pan-cancer** (Cell Reports 2025) | 2025 | 11 種癌症 81,567 個 MHB，MHB-associated DEGs → oncogenic pathways | ISM 可引入 MHB 做 CpG 功能分層 |
| **NANOME** (PMC 2025) | 2025 | Nextflow pipeline for haplotype-aware ASM detection | ISM 可整合 NANOME 的標準化流程 |
| **CRC TSG methylation** (MDPI 2025) | 2025 | 同時分析 somatic variants + methylation in TSG | 最接近 ISM「variant-methylation integration」目標 |

---

## 三、學界共識、分歧、與未驗證之處

### 3.1 已有共識的事項

1. **Purity/CNV 校正必要性**：癌症甲基化訊號強烈受 purity 與 CNV/LOH 影響（TRACERx CAMDAC、Long-Read POG、多項 deconvolution 研究一致）。不做校正容易把混合訊號誤判成腫瘤特異甲基化。

2. **長讀序的多模態價值**：長讀序不只做 methylation，而是把 phasing + mutation + SV + methylation 放在同一分子上（Long-Read POG、NanoMethPhase、t-nanoEM、LongPhase-S 一致）。

3. **LOH-SNV 作為天然 tumor/normal 分離標籤**：在 LOH 區域，variant-bearing reads 近似純腫瘤來源（TRACERx、Long-Read POG 均使用此策略）。

4. **Promoter methylation + LOH = classical second hit**：BRCA1/RAD51C/CDKN2A 等案例已被多項研究驗證（Long-Read POG、t-nanoEM、CRC TSG 研究）。

5. **ONT 甲基化/錯誤型態的耦合**：Nanopore 在甲基化位點容易出現系統性 SNP 誤呼叫（LongPhase 方法段落直接討論此問題）。

### 3.2 存在分歧或未充分驗證之處

1. **局部甲基化 cluster 的生物意義**：
   - TRACERx 用 ITMD/PDR 量化異質性但在 cohort 層級
   - EVOFLUx 找到 fCpG 可追蹤演化但在 bulk 層級
   - **ISM 在 read-level per-SNV 做 clustering 的生物意義尚缺獨立驗證**
   - 尚無研究直接證明「某 SNV 周圍的甲基化 cluster 結構」能可靠區分 TP/FP

2. **甲基化驅動 vs 乘客的定義多樣**：
   - TRACERx: M_R/M_N ratio + expression correlation → 正向選汰
   - MethSig: beta regression + cohort statistics → driver enrichment
   - EVOFLUx: fCpGs 大多中性 → lineage tracking
   - **三種觀點互補但也意味 ISM 需明確定位在哪個層級**

3. **Read-level methylation clustering 的技術限制**：
   - ONT methylation probability (ML) 在中間值 (p~0.5) 時雜訊大
   - ISM 用 Bernoulli distance 降權中間值，但效果因樣本而異
   - MethylBERT 用 Transformer 學習 pattern，可能比硬規則更 robust
   - **尚無 head-to-head 比較 clustering-based vs ML-based read classification**

4. **Cross-platform/cross-sample 可攜性**：
   - ISM 實驗顯示 PairwiseMedianDist 方向具樣本依賴
   - **尚無研究系統性驗證 read-level methylation features 的跨平台穩定性**

5. **Germline ASM vs Somatic ASM 的區分**：
   - ISM 方法學審查確認 FP 主因是 germline ASM（94% FP 帶 AlleleSig=True）
   - NanoMethPhase 可偵測 ASM 但不區分 germline/somatic
   - **需要 normal sample 甲基化參考才能區分，目前 tumor-only 場景缺乏解法**

---

## 四、InterSubMod 是否已被其他研究覆蓋？

### 4.1 最接近的競爭/先行研究

| 工具/研究 | 相似度 | 關鍵差異 |
|-----------|-------|---------|
| **NanoMethPhase** | 高（同為 ONT + methylation + phasing） | NanoMethPhase: 全基因組 ASM；ISM: per-SNV clustering |
| **PRISM** | 高（methylation pattern → subclone） | PRISM: reference-free bulk deconvolution；ISM: read-level with variant anchor |
| **MethylBERT** | 中高（read-level methylation pattern） | MethylBERT: supervised ML deconvolution；ISM: unsupervised clustering + significance |
| **t-nanoEM** | 中（somatic variant-linked methylation） | t-nanoEM: converted reads + targeted capture；ISM: native MM/ML + WGS |
| **CAMDAC (TRACERx)** | 中（methylation deconvolution） | CAMDAC: bulk CN-aware deconvolution；ISM: read-level per-SNV |

### 4.2 ISM 的獨特定位（尚未被完全覆蓋之處）

**沒有人完全做了 ISM 想做的事。** 具體來說：

1. **「以 somatic SNV 為錨點 + read-level methylation clustering + 雙向顯著性驗證」的完整框架** — 這個組合是 ISM 獨有的
2. **在 variant verification 場景中使用甲基化 pattern** — 大多數研究用甲基化做 deconvolution 或 ASM detection，而非 variant verification
3. **同時輸出 cluster-first 與 label-first 的 VerificationClass** — 這個多層驗證框架無直接對標

**但問題在於**：ISM 的獨特性也是它的瓶頸 — read-level per-SNV methylation clustering 在目前的實驗中未能產生超越 caller-only 的 discriminative power。

---

## 五、突破方向分析與建議

### 方向 A：Copy-number / Purity-aware Read-label Modeling

**靈感來源**：TRACERx CAMDAC + Long-Read POG

- ISM 已有 per-read ALT/REF/UNKNOWN 標記
- 在 LOH 區域，ALT reads ~ tumor-derived，REF reads ~ normal contamination
- 用此標籤估計「局部 tumor purity」與「腫瘤 vs 正常甲基化背景」
- **在校正後的腫瘤甲基化上做聚類**，而非在混合訊號上聚類
- 有 TRACERx LOH-SNV → PDR 驗證（R>0.8）作為外部支撐

### 方向 B：Gene-level / Mechanism-level Evidence Integration

**靈感來源**：Long-Read POG IMPALA + TRACERx TSG analysis

- 加一層 gene-level integrator：把同一基因的多個 SNV window 結果彙總
- 結合 CNV/LOH 資訊，回答「此基因是否有 biallelic inactivation 證據」
- 定義 second-hit archetypes：mutation→LOH, LOH→mutation, promoter_methylation→LOH, LOH→promoter_methylation
- 直接對應 BRCA1/RAD51C/CDKN2A 等已知臨床問題

### 方向 C：CpG 功能分層與智慧選擇

**靈感來源**：EVOFLUx fCpGs + TRACERx M_R/M_N + MHB 研究

- 將 CpG 依功能分為 Regulatory / fCpG-like / MHB-core
- 不同目標選不同 CpG 子集（driver/ASM vs lineage vs structural）
- EVOFLUx 已識別 978 fCpGs 可直接借用

### 方向 D：Normal Sample 甲基化參考

**靈感來源**：方法學審查結論（FP 來自 germline ASM）

- 在 paired 場景，normal BAM 已有 MM/ML 標籤
- 建立 per-region normal methylation baseline
- ISM 聚類結果若與 normal 甲基化 pattern 一致 → downgrade 為 germline ASM

### 方向 E：MethylBERT-style ML Read Classification

**靈感來源**：MethylBERT (2025)

- 用 ML 模型學習 read-level methylation pattern → tumor/normal classification
- 取代或補充 ISM 的 unsupervised clustering
- MethylBERT 已展示 read-level pattern 可超越 mean methylation

---

## 六、研究者確認的漸進策略

> **確認日期**：2026-03-24
> **策略邏輯**：先在現有框架下做方法學修改與測試 → 引入 normal 資料分析 purity/LOH/CNV 影響 → 多點關聯分析 → CpG 功能分層與高關聯位點篩選輸出 → 推斷 second hit 順序、subclone 演化分化、突變時間/區域 → 最終回頭加強 somatic TP 驅動判斷

| 階段 | 方向 | 目標 | 外部先例 | 預期產出 |
|------|------|------|---------|---------|
| **Phase 1** | ML read classification（方向 E） | 在現有框架內改進 read-level pattern recognition，突破 clustering 天花板 | MethylBERT | 改良的 read classification method + head-to-head 比較 |
| **Phase 2** | Normal methylation ref + CN/Purity-aware（方向 A+D） | 引入 normal 甲基化基線 + purity/LOH/CNV 校正，解決 germline ASM 根因 | CAMDAC (TRACERx) | Purity-corrected methylation matrix + germline/somatic ASM 區分 |
| **Phase 3** | Gene-level evidence integration（方向 B） | 多點關聯 → gene-level biallelic event narrative → second hit 順序推斷 | Long-Read POG IMPALA | Gene-level evidence panel + second-hit archetype classification |
| **Phase 4** | CpG 功能分層（方向 C） | 篩選重要甲基位點與高關聯位點 → 輸出供後續學者觀察研究 | EVOFLUx/TRACERx M_R/M_N/MHB | Annotated CpG catalogue + subclone evolution timeline |
| **最終目標** | 整合回饋 | 統整數據 → 加強判斷 somatic 位點與甲基位點是否為癌症驅動 TP | MethSig + all above | Methylation-informed variant driver score |

---

## 七、值得進一步觀察與研究的方向

### 7.1 立即可做的觀察

1. **LOH-SNV PDR 驗證**：在 HCC1395 LOH 區域，用 ALT/REF reads 分別計算 methylation pattern
2. **fCpG overlap 分析**：將 EVOFLUx 的 978 fCpGs 與 ISM 視窗內 CpGs 交叉
3. **Gene-level aggregation pilot**：選 10 個已知 TSGs 彙總多 SNV window 結果
4. **Normal methylation baseline**：對 normal BAM 跑 ISM 同一視窗做比較

### 7.2 中期研究方向

5. **MHB-aware distance**：將 MHB 資訊納入距離計算
6. **Purity-aware clustering**：引入 LongPhase-S 的 purity estimation 作為先驗
7. **Cross-cancer MethSig overlap**：比較 ISM 在 driver gene vs non-driver gene 的表現
8. **MethPhaser integration**：用甲基化延伸 phase block

### 7.3 長期突破性方向

9. **Single-molecule variant-methylation co-evolution model**
10. **Methylation-aware variant recalibration**
11. **Pan-cancer methylation verification benchmark**

---

## 八、論文定位建議

ISM 最強的學術定位不是「variant filter」，而是：

> **「以長讀序 somatic variant 為錨點，整合甲基化 pattern、haplotype、copy-number 資訊，提供 read-level epigenetic context for variant interpretation」**

這個定位與 Long-Read POG、TRACERx、t-nanoEM 互補，且目前無人完整實作。

---

## 附錄：論文出處

1. O'Neill et al., Cell Genomics 4(11), 2024 — Long-Read POG
2. TRACERx consortium, Nat Genet, 2025 — NSCLC methylation heterogeneity
3. Gabbutt & Duran-Ferrer et al., Nature, 2025 — EVOFLUx
4. Park et al., bioRxiv, 2024 — DeepSomatic
5. Ho et al., bioRxiv, 2025 — LongPhase-S
6. Kunigo et al., Cell Reports Methods, 2025 — t-nanoEM
7. Akbari et al., Genome Biology 22, 2021 — NanoMethPhase
8. Landau et al., Cancer Discovery 11(9), 2021 — MethSig
9. Treangen et al., Nat Commun 15:5327, 2024 — MethPhaser
10. Jeong et al., Nat Commun, 2025 — MethylBERT
11. Lee et al., Bioinformatics 35(14), 2019 — PRISM
12. Cell Reports 2025 — MHB pan-cancer map
