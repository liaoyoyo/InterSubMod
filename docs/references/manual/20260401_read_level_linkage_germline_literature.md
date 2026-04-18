# Read-Level Haplotype Linkage 用於 Germline/Somatic 區分：文獻調查報告

<!--
建立時間: 2026-04-01 16:00
目標: 調查是否有已發表方法利用 long-read haplotype linkage / read-level variant co-occurrence 區分 germline vs somatic variants，特別是 tumor-only 場景下利用已知 germline anchor 標記未知位點的方法
處理範圍: 學術文獻搜尋（PubMed, bioRxiv, GitHub, Web），涵蓋 haplotype-aware variant calling, read-level phasing, tumor-only approaches, variant co-occurrence
關聯檔案:
  - docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md
  - docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
  - docs/CURRENT_FOCUS.md
-->

## 1. 搜尋概述

- **核心問題**: 在 tumor-only ONT long-read 情境下，能否利用 read-level variant co-occurrence（同一條 read 同時覆蓋 PON 已標記的 germline 位點 + 未標記的候選位點）來推斷候選位點也是 germline，從而過濾殘餘 FP？
- **搜尋時間**: 2026-04-01
- **資料來源數**: 25+ 篇論文/預印本/工具文檔
- **關鍵詞**: long read phasing somatic germline discrimination, haplotype-aware variant filtering, read-level linkage, variant co-occurrence, ancestral haplotype, tumor-only somatic calling, guilt by association

---

## 2. 核心發現

### 2.1 最直接相關：已有方法利用 Haplotype/Phasing 資訊區分 Germline vs Somatic

#### 2.1.1 smrest — 無 Matched Normal 的 Haplotype-Based Somatic Calling（最相關）

**論文**: Simpson, J.T. (2024). "Detecting Somatic Mutations Without Matched Normal Samples Using Long Reads." bioRxiv, 2024.02.26.582089.

**核心思想**:
- **理論框架**: Somatic mutation 只出現在一個 haplotype 上（混合 REF+ALT alleles），而 germline heterozygous variant 在該 haplotype 上是純 ALT。因此，一旦 reads 被 haplotagged，可透過 Bayesian 分類判斷每個位點是 somatic、germline、還是 sequencing error。
- **方程式**: P(read|somatic) = (1-epsilon) * alpha + epsilon * (1-alpha)，其中 alpha = tumor purity, epsilon = sequencing error rate。
- **實作流程**: genotype known SNPs -> WhatsHap phasing -> haplotagging -> per-haplotype Bayesian classification。

**關鍵結果**:
- 模擬：80x, 1% error rate 下，phased data F1 > 0.9（purity 0.23-0.84），unphased 僅 0.42-0.52
- 真實數據 COLO829：phased regions recall 0.87-0.96，整體 0.61-0.86
- 限制：僅測試高突變率 cell lines；不處理 indels；性能在極端 purity 下退化

**與本研究假說的關係**: smrest 的核心思路與我們的假說高度相似——都是利用 haplotype information 來區分 germline vs somatic。但 smrest 是 **de novo calling**，而我們的假說更窄：利用 **已知 PON-flagged germline anchor** 在同一條 read 上的 co-occurrence 來 **tag 未知位點**。smrest 不直接使用 PON-flagged anchors，而是從 haplotype-level allele mixture pattern 推斷。

#### 2.1.2 ClairS — Ancestral Haplotype Filter（Paired Mode）

**論文**: Zheng et al. (2023). "ClairS: a deep-learning method for long-read somatic small variant calling." bioRxiv, 2023.08.17.553778.

**核心思想**:
- Step 1: 用 Clair3 + LongPhase phasing germline variants + haplotagging reads
- Step 3: **Ancestral haplotype filter** — 排除無法確定單一 ancestral haplotype origin 的 somatic variant candidates。somatic variant 應來自單一 parental haplotype (HP1 或 HP2)；如果 ALT alleles 出現在兩個 haplotype 上，可能是 germline artifact。
- 例外: 高 VAF somatic variants 可能因 CNA/clonal duplication 出現在雙 haplotype。

**與本研究假說的關係**: ClairS 的 ancestral haplotype filter 是目前最接近「利用 haplotype 資訊過濾非 somatic variants」的實作。但它需要 **matched normal**，且判斷邏輯是「somatic should be single-haplotype origin」而非「co-occurrence with known germline = germline」。

#### 2.1.3 LongPhase-S — Somatic Haplotype Recalibration

**論文**: Ho et al. (2025). "LongPhase-S: purity estimation and variant recalibration with somatic haplotagging for long-read sequencing." bioRxiv, 2025.11.20.689492.

**核心思想**:
- 將 somatic reads 錨定到 parental germline haplotypes（HP1-1 / HP2-1 / HP3）
- 重建 somatic haplotype 後，識別與 somatic haplotype 不一致的 false somatic variants，標記為 "LowQual"
- 整合 purity estimation 進行 purity-aware recalibration

**關鍵結果**: ClairS SNV F1 提升 +4.5%, Indel +7.1%

**與本研究假說的關係**: LongPhase-S 利用 germline haplotype 作為 anchor 來 recalibrate somatic variants。conceptually 與我們的 "germline anchor on same read" 相近，但 (1) 需要 matched normal，(2) 是全域 haplotype reconstruction 而非 read-level co-occurrence。

#### 2.1.4 Octopus — Bayesian Haplotype-Aware Caller（Co-Phasing Somatic + Germline）

**論文**: Cooke et al. (2021). "A unified haplotype-based method for accurate and comprehensive variant calling." Nature Biotechnology, 39:1284-1292.

**核心思想**:
- Polymorphic Bayesian genotyping model，在統一框架內同時 call germline + somatic variants
- **關鍵特性**: 能將 somatic mutations 與 germline variants co-phase。在 52,373 somatic mutations 中，13,584 (26%) 被與一個或多個 germline variant phased 在一起
- Germline/somatic VAF 決策邊界會自動根據測序深度調整（不需 preset threshold）
- 部分 high-VAF somatic mutations 因為與 nearby germline variants phased 而被正確分類

**與本研究假說的關係**: Octopus 直接 co-phase somatic 和 germline variants，概念上最接近「利用 phasing linkage 來輔助 classification」。但 Octopus 主要設計給 short reads（Illumina），long-read 支援有限。

### 2.2 間接相關：Read-Level / Linked-Read Variant Co-Occurrence Analysis

#### 2.2.1 SomaticHaplotype — Linked-Read Somatic Mutation Phasing

**論文**: Khandkar et al. (2024). "Somatic mutation phasing and haplotype extension using linked-reads in multiple myeloma." Genome Biology, 25:218.

**核心思想**:
- **Linked alleles 概念**: 在同一 barcode 上 co-occur 的 alleles 很可能來自同一條 HMW DNA 分子
- 利用已知 germline variants 的 haplotype assignment（Long Ranger phased, ~99% accuracy）推斷同一 barcode 上的 somatic mutation 的 haplotype
- **Threshold**: 要求 >= 91% 的 linked alleles 一致指向同一 haplotype（precision 0.997, recall 0.936）
- **Phase block extension**: 利用同一患者不同 time-point 樣本的 germline variants 橋接 phase blocks，延伸 4.6 倍

**關鍵結果**: 79.4% (16,440/20,705) somatic mutations 成功 phased

**與本研究假說的關係**: 這是最接近「guilt by association」概念的實作。SomaticHaplotype 利用同一 barcode（分子）上的 germline variant haplotype 來推斷 somatic mutation 的 haplotype。我們的假說可視為 **逆向版本**：已知 germline anchor 用來推斷同一 read 上的未知 variant 也是 germline。核心差異：SomaticHaplotype 是 **assign haplotype**，我們是 **assign germline identity**。

#### 2.2.2 Lancet — Linked-Read Haplotype-Aware Somatic Caller

**論文**: Wala, Berger et al. (2021). "Somatic variant analysis of linked-reads sequencing data with Lancet." Bioinformatics, 37(9):1310-1316.

**核心思想**:
- 將 barcode + haplotype 資訊整合進 colored De Bruijn graph local-assembly framework
- 計算 barcode-aware coverage，識別與 local haplotype structure 不一致的 variants
- 輸出 VCF 包含 allele-specific barcodes，支援下游識別同一 haplotype 上 co-occurring somatic mutations

**關鍵結果**: 相較原始 Lancet (COLO829)，移除 69% SNV FP 和 77% InDel FP，sensitivity 僅微幅下降

**與本研究假說的關係**: Lancet 展示了利用 haplotype structure 進行 variant filtering 的有效性。但 (1) 需要 matched normal，(2) 基於 linked-reads (10x Genomics)，非 long reads。

### 2.3 Tumor-Only Specific Approaches

#### 2.3.1 ClairS-TO — Long-Read Tumor-Only Caller

**論文**: Chen, Zheng et al. (2025). "ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling." Nature Communications, 16:9630.

**核心思想**:
- 三層過濾: 9 hard filters + 4 PON databases + Verdict 統計分類（基於 purity, ploidy, copy number）
- **Haplotype 使用**: MultiHap flag — ALT allele 出現在多於一個 haplotype 時標記
- 並未實作 read-level co-occurrence 分析

**與本研究假說的關係**: ClairS-TO 是我們使用的 caller，其 MultiHap flag 概念上與 "single haplotype = somatic" 有關，但不使用 PON-flagged anchors 做 read-level linkage。我們的假說是在 ClairS-TO **之後** 的額外 post-hoc filtering。

#### 2.3.2 SGZ — VAF-Based Somatic/Germline Zygosity Classification

**論文**: Sun et al. (2018). "A computational approach to distinguish somatic vs. germline origin of genomic alterations from deep sequencing of cancer specimens without a matched normal." PLoS Computational Biology, 14(2):e1005965.

**核心思想**:
- 利用 VAF 結合 tumor purity, ploidy, local copy number 建立 somatic vs germline 的 statistical model
- 需要 >500x 深度（targeted panel）
- 85% of cases 可做出預測，95-99% accuracy

**與本研究假說的關係**: SGZ 是 site-level VAF 方法，不使用 read-level linkage。我們已在 O1-O13 確認 site-level features 對 TO FP 的 AUC < 0.64，因此 SGZ 類方法的效果上限已知。

#### 2.3.3 LRSomatic — Comprehensive Long-Read Pipeline

**論文**: bioRxiv, 2026.02.26.707772 (Preprint, Feb 2026).

**核心思想**: Nextflow pipeline 整合 SNV/indel/SV/CN calling，支援 PacBio HiFi + ONT，支援 tumor-only 模式。整合 Fiber-seq 做 haplotype-specific chromatin accessibility 分析。

**與本研究假說的關係**: 最新的 long-read somatic pipeline，但未見 read-level germline anchor linkage 的概念。

### 2.4 SV-Level Haplotype-Aware Methods（概念借鑑）

#### 2.4.1 Severus — Haplotype-Aware Somatic SV Calling

**論文**: Kolmogorov et al. (2025). "Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing." Nature Biotechnology.

**核心概念**: Somatic SV 的 supporting reads 應 phased 到單一 haplotype；artifact 則分散在兩個 haplotype。這與 "somatic SNV should be single-haplotype" 原則一致。

#### 2.4.2 SAVANA — Single-Haplotype Resolution SV + CNA

**論文**: Espejo Valle-Inclan et al. (2025). "SAVANA: reliable analysis of somatic structural variants and copy number aberrations using long-read sequencing." Nature Methods.

**核心概念**: ML-based somatic/artifact 分類，在 single-haplotype resolution 下運作。99 tumors benchmarked, 13-82x higher specificity。

### 2.5 Phasing 基礎設施

| 工具 | 特點 | 與本假說的關係 |
|------|------|---------------|
| **WhatsHap** | Read-backed phasing, 支援 PacBio/ONT | smrest 使用; 可直接用於 tumor-only phasing |
| **LongPhase** | SNP/SV/甲基化共同 phasing | 我們已在用；提供 haplotagging 基礎 |
| **LongPhase-TO** | Tumor-only phasing + LOH | 我們已在用；提供 HP:i tag + PS tag |
| **HapCUT2** | 多平台 phasing | 不直接處理 somatic; cancer genome 會干擾 |

---

## 3. 與本研究假說的直接對照分析

### 3.1 假說回顧

> 如果一條 read 同時覆蓋 (1) 一個 PON 標記的 germline 位點（germline anchor）和 (2) 一個未被 PON 標記的候選位點，那麼兩者的 co-occurrence/linkage 可作為「候選位點也是 germline」的證據。

### 3.2 文獻中最接近的方法

| 方法 | 相似度 | 差異 |
|------|--------|------|
| **SomaticHaplotype (linked alleles)** | 最高 | 用已知 germline haplotype 推斷同一分子上的未知 variant 的 haplotype；但目的是 phasing 而非 germline identity |
| **smrest** | 高 | 用 haplotype-level allele mixture 判斷 somatic vs germline；但不使用 PON anchors |
| **ClairS ancestral haplotype filter** | 中高 | 判斷 somatic candidate 是否有 single haplotype origin；但需 matched normal |
| **Octopus co-phasing** | 中 | 同時 phase germline + somatic；但主要用於 short reads |
| **Lancet haplotype structure** | 中 | Barcode-aware filtering；但需 linked-reads + matched normal |

### 3.3 文獻空白（Gap）

**本假說的核心創新點在於**：

1. **逆向 guilt by association**: 現有方法大多是「利用 germline haplotype 來 phase somatic variants」（正向），我們的假說是「利用 known germline 來 tag unknown variants as germline」（逆向）。
2. **Tumor-only + PON anchor**: 現有 haplotype-based methods（ClairS, LongPhase-S, Octopus）幾乎都需要 matched normal。在 tumor-only 場景下利用 PON-flagged variants 作為 germline anchor 的方法尚未見文獻。
3. **Read-level co-occurrence 而非 haplotype reconstruction**: 我們不需要完整的 haplotype reconstruction，只需要在 individual read 上觀察兩個 variant 是否 co-occur。這比全域 phasing 更直接且計算更簡單。

---

## 4. 技術可行性評估

### 4.1 理論基礎

| 原則 | 支持度 | 說明 |
|------|--------|------|
| Germline variants 在同一 haplotype 上 co-segregate | 強 | 這是 LD/phasing 的基本原理，long reads 直接觀察 |
| Somatic variants 不應與 germline variants 有 LD | 中-強 | Somatic mutations 獨立於 germline haplotype 產生，但仍發生在某一 parental haplotype 上 |
| PON-flagged germline variants 是可靠 anchor | 中 | PON sensitivity 有限（~99.5% 但非 100%），anchor 本身的準確性取決於 PON 品質 |

### 4.2 潛在問題與挑戰

| 問題 | 嚴重度 | 說明 |
|------|--------|------|
| **Somatic variant 也出現在同一 haplotype 上** | 高 | Somatic mutation 必然發生在某一 parental haplotype 上，因此也會與該 haplotype 上的 germline variants co-occur。單純的 co-occurrence 無法區分「同一 haplotype 上的 germline pair」和「同一 haplotype 上的 germline + somatic pair」。 |
| **需要 higher-order logic** | 高 | 單純的「同一 read 上有 PON germline」不夠；需要進一步利用「PON germline 在 HP1 上，候選位點也在 HP1 上，且 HP2 上也有候選位點的 ALT」（= heterozygous germline）vs「PON germline 在 HP1 上，候選位點只在 HP1 上有 ALT」（= could be somatic）。這基本上就回到了 smrest/ClairS 的 haplotype-level analysis。 |
| **LOH 區域干擾** | 中 | LOH 區域中 germline heterozygous 變成 homozygous，haplotype 資訊失效 |
| **Phase switch errors** | 中 | ONT ~1-5% switch error rate 會引入 noise |
| **Read length vs variant spacing** | 低-中 | ONT N50 ~10-40kb，兩個 variant 需在同一 read 覆蓋範圍內 |
| **Coverage asymmetry** | 低 | 某些位點覆蓋不足可能導致 linkage evidence 不夠 |

### 4.3 關鍵洞察：為何單純 co-occurrence 不夠

**核心問題**: 一個 somatic variant 出現在（例如）HP1 haplotype 上。同一條 HP1 read 上也會有 germline variants（包括 PON-flagged 的）。因此：

- PON-flagged germline (on HP1) + 真正 germline variant (on HP1) -> co-occur on same read -> 正確推斷
- PON-flagged germline (on HP1) + somatic variant (on HP1) -> co-occur on same read -> **錯誤推斷為 germline**

**解法**: 不能只看「是否 co-occur」，而要看「是否在 **兩個** haplotype 上都有 ALT」。真正的 germline heterozygous variant 在 HP1 和 HP2 上都有 ALT reads（各約 50%），而 somatic variant 只在一個 haplotype 上有 ALT。這就是 smrest/ClairS 的 ancestral haplotype filter 的核心邏輯。

### 4.4 可行的修正方案

原始假說（pure co-occurrence）需要修正為：

**修正假說 A: Haplotype-Consistent Co-Occurrence**
> 如果一個候選位點的 ALT alleles 分布在兩個 haplotype 上（HP1 和 HP2），且在兩個 haplotype 上都與 PON-flagged germline variants co-occur，那麼該候選位點很可能是 germline（heterozygous）。

**修正假說 B: Multi-Anchor Density**
> 如果一個候選位點周圍（同一 read 覆蓋範圍內）有異常多的 PON-flagged germline variants 與其 co-occur（相較 genome-wide 平均），這可能暗示該區域的 PON coverage 不足（private germline 較多的區域）。

**修正假說 C: Haplotype Imbalance as Somatic Evidence**
> 反向使用：如果候選位點的 ALT alleles 只出現在一個 haplotype 上（且 AF 符合 purity 預期），這是 somatic 的正向證據。這就是 smrest/ClairS 已實作的方法。

---

## 5. 資料來源評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| Simpson 2024 (smrest) | bioRxiv preprint | 中-高 | 理論框架嚴謹，prototype 階段，高突變 cell line 驗證 |
| Zheng et al. 2023 (ClairS) | bioRxiv preprint | 高 | HKU-BAL 團隊，廣泛 benchmarked |
| Chen, Zheng et al. 2025 (ClairS-TO) | Nature Communications | 高 | 已正式發表，多中心驗證 |
| Ho et al. 2025 (LongPhase-S) | bioRxiv preprint | 高 | CCU Lab，與本專案直接相關 |
| Cooke et al. 2021 (Octopus) | Nature Biotechnology | 高 | 頂刊正式發表 |
| Khandkar et al. 2024 (SomaticHaplotype) | Genome Biology | 高 | 正式發表，23 samples 驗證 |
| Wala et al. 2021 (Lancet) | Bioinformatics | 高 | NYGC 團隊 |
| Sun et al. 2018 (SGZ) | PLoS Computational Biology | 中-高 | Foundation Medicine, 已被臨床使用 |
| Kolmogorov et al. 2025 (Severus) | Nature Biotechnology | 高 | SV-focused 但概念可借鑑 |
| Espejo Valle-Inclan et al. 2025 (SAVANA) | Nature Methods | 高 | 99 tumor-normal pairs 驗證 |
| O'Neill et al. 2024 (Long-read POG) | Cell Genomics | 高 | 189 patient cohort |
| LRSomatic 2026 | bioRxiv preprint | 中 | 極新，尚無外部驗證 |

---

## 6. 結論與建議行動

### 6.1 文獻結論

1. **Haplotype-based germline/somatic discrimination 是已被多方驗證的有效方法**，但所有已發表方法的核心邏輯是「somatic variant 在單一 haplotype 上」而非「read-level co-occurrence with known germline」。

2. **原始假說（pure co-occurrence = germline evidence）有根本性缺陷**：somatic variants 也會與 germline variants co-occur（因為它們必然發生在某一 parental haplotype 上）。因此，單純的 co-occurrence 無法區分 germline 和 somatic。

3. **修正後的假說等價於 smrest/ClairS 的 haplotype-level analysis**：判斷 candidate 是否在兩個 haplotype 上都有 ALT（= germline）vs 只在一個 haplotype 上有 ALT（= somatic）。

4. **Tumor-only 場景下的 haplotype-based filtering 是一個公認的未解難題**。smrest 是目前唯一直接針對此場景的 prototype，但僅在高突變 cell lines 上驗證。

5. **文獻空白**: 沒有找到任何已發表方法直接利用 PON-flagged variants 作為 germline anchor 進行 read-level linkage analysis。

### 6.2 建議行動

| 優先級 | 行動 | 理由 |
|--------|------|------|
| **P0** | 放棄 pure co-occurrence 假說 | 理論上有根本缺陷（somatic 也 co-occur with germline） |
| **P1** | 評估 smrest-like haplotype imbalance approach | 在 TO 場景下，利用 LongPhase-TO 已產生的 HP tags，計算每個 candidate 的 haplotype balance ratio（HP1 ALT / HP2 ALT） |
| **P1** | 評估 ClairS-TO MultiHap flag 的 residual value | 目前 ClairS-TO 已有 MultiHap flag，確認這是否已捕獲 haplotype imbalance 資訊 |
| **P2** | 設計 "germline density" metric | 計算 candidate 周圍固定窗口內 PON-flagged variants 的密度，作為「PON-poor region」指標 |
| **P2** | 測試 haplotype-level AF consistency | 在 HP1 和 HP2 上分別計算 candidate 的 AF，如果兩者都 ~0.5（即兩個 haplotype 上 ALT 都存在），這是 germline 的強證據 |
| **P3** | 與 LongPhase-TO 團隊討論是否可增加 haplotype-based germline filter | 作為 tool-level 整合的可能性 |

### 6.3 與現有研究的定位

本調查確認了**haplotype-based approach 在 TO 場景下的理論有效性但實作空白**。結合 MEMORY 中記錄的「TO Germline FP 鑑別 NO-GO」（60+ 特徵 AUC < 0.64），haplotype imbalance 可能是**尚未被系統性測試的新維度**：

- 現有 60+ 特徵全部是 **site-level** 特徵（AF, GQ, depth, methylation, etc.）
- Haplotype balance ratio 是 **read-level aggregated** 特徵，概念上不在之前的搜索空間內
- 但需注意：ClairS-TO 的 MultiHap flag 已部分捕獲此資訊，需先確認其效果

---

## 7. 參考文獻完整列表

1. Simpson, J.T. (2024). Detecting Somatic Mutations Without Matched Normal Samples Using Long Reads. bioRxiv. DOI: 10.1101/2024.02.26.582089
   - URL: https://www.biorxiv.org/content/10.1101/2024.02.26.582089v1

2. Zheng, Z. et al. (2023). ClairS: a deep-learning method for long-read somatic small variant calling. bioRxiv. DOI: 10.1101/2023.08.17.553778
   - URL: https://www.biorxiv.org/content/10.1101/2023.08.17.553778v1.full

3. Chen, L., Zheng, Z. et al. (2025). ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling. Nature Communications, 16:9630.
   - URL: https://www.nature.com/articles/s41467-025-64547-z

4. Ho, M.-E. et al. (2025). LongPhase-S: purity estimation and variant recalibration with somatic haplotagging for long-read sequencing. bioRxiv. DOI: 10.1101/2025.11.20.689492
   - URL: https://www.biorxiv.org/content/10.1101/2025.11.20.689492v1.full

5. Cooke, D. et al. (2021). A unified haplotype-based method for accurate and comprehensive variant calling. Nature Biotechnology, 39:1284-1292.
   - URL: https://www.nature.com/articles/s41587-021-00861-3

6. Khandkar, I. et al. (2024). Somatic mutation phasing and haplotype extension using linked-reads in multiple myeloma. Genome Biology, 25:218.
   - URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC11326269/

7. Wala, J. et al. (2021). Somatic variant analysis of linked-reads sequencing data with Lancet. Bioinformatics, 37(9):1310-1316.
   - URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC8487631/

8. Sun, J.X. et al. (2018). A computational approach to distinguish somatic vs. germline origin of genomic alterations from deep sequencing of cancer specimens without a matched normal. PLoS Computational Biology, 14(2):e1005965.
   - URL: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005965

9. Kolmogorov, M. et al. (2025). Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing. Nature Biotechnology.
   - URL: https://www.nature.com/articles/s41587-025-02618-8

10. Espejo Valle-Inclan, J. et al. (2025). SAVANA: reliable analysis of somatic structural variants and copy number aberrations using long-read sequencing. Nature Methods.
    - URL: https://www.nature.com/articles/s41592-025-02708-0

11. O'Neill, K. et al. (2024). Long-read sequencing of an advanced cancer cohort resolves rearrangements, unravels haplotypes, and reveals methylation landscapes. Cell Genomics, 4(11).
    - URL: https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00293-3

12. LRSomatic (2026). A highly scalable and robust pipeline for somatic variant calling in long-read sequencing data. bioRxiv. DOI: 10.64898/2026.02.26.707772
    - URL: https://www.biorxiv.org/content/10.64898/2026.02.26.707772v1

13. Tan, K.T. et al. (2021). Haplotype-resolved germline and somatic alterations in renal medullary carcinomas. Genome Medicine, 13:114.
    - URL: https://link.springer.com/article/10.1186/s13073-021-00929-4

14. Shafin, K. et al. (2021). Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads. Nature Methods, 18:1322-1332.
    - URL: https://www.nature.com/articles/s41592-021-01299-w

15. gnomAD Variant Co-Occurrence (Phasing) Information (2021).
    - URL: https://gnomad.broadinstitute.org/news/2021-07-variant-co-occurrence-phasing-information-in-gnomad/

---

**文檔版本**: v1.0
**最後更新**: 2026-04-01
