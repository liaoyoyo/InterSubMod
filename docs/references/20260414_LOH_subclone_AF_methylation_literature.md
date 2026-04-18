# LOH 區域中間 AF 與亞克隆結構 — AF + 甲基化文獻調查

<!--
建立時間: 2026-04-14 00:00
目標: 調查 LOH 區域中 intermediate allele frequency 是否指示亞克隆結構，以及甲基化能否提供額外證據
處理範圍: 學術文獻搜尋，涵蓋 deletion LOH、cnLOH、multi-site linkage、ASM、epi-allele diversity、長讀長工具
關聯檔案:
  - docs/references/20260331_LOH區域甲基化模式與事件順序文獻調查_01.md
  - docs/references/20260401_ASM偵測方法與預期比例文獻調查_01.md
  - docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
-->

## 1. 搜尋概述

- **核心問題**: 在腫瘤 LOH 區域中，intermediate AF 是否指示 subclonal structure？甲基化模式能否提供額外證據？
- **搜尋時間**: 2026-04-14
- **資料來源數**: 50+ 篇論文 / 預印本 / 工具文件 / 綜述
- **關鍵詞**: subclonal LOH, deletion LOH CN=1, cnLOH CN=2, allele frequency tumor purity, cancer cell fraction, allele-specific methylation, epi-allele diversity, long-read subclone, nanopore methylation heterogeneity, TITAN, Battenberg, CAMDAC, PyClone, EVOFLUx, MethPhaser, NANOME

---

## 2. 核心 VAF 數學框架

### 2.1 通用 VAF 公式

在腫瘤樣本中，一個體細胞突變的觀測 VAF 取決於四個因素（Dentro et al., Cold Spring Harbor Perspectives in Medicine, 2017; Tarabichi et al., Nature Methods, 2021）：

```
VAF = (CCF * m * rho) / (rho * CN_tumor + (1 - rho) * CN_normal)
```

其中：
- **rho (ρ)**: 腫瘤純度（tumor purity），腫瘤細胞佔總細胞的比例
- **CCF (phi)**: 癌細胞分數（cancer cell fraction），帶有該突變的腫瘤細胞比例
- **m**: 突變多重性（multiplicity），突變所在的拷貝數
- **CN_tumor**: 該位點在腫瘤細胞中的總拷貝數
- **CN_normal**: 該位點在正常細胞中的總拷貝數（通常 = 2）

此公式可以反推 CCF：

```
CCF = VAF * (rho * CN_tumor + (1 - rho) * CN_normal) / (m * rho)
```

對於 **clonal mutation（CCF = 1）** 的期望 VAF：

```
VAF_expected = (m * rho) / (rho * CN_tumor + (1 - rho) * 2)
```

---

## 3. Deletion LOH (CN=1): 期望 AF 與偏離的意義

### 3.1 理論期望

在 deletion LOH 區域（一個 allele 被刪除，保留一個拷貝，CN_tumor = 1）：

**Germline heterozygous variant（保留的 allele 上有 ALT）：**

```
VAF_expected = (1 * rho + 1 * (1 - rho)) / (rho * 1 + (1 - rho) * 2)
         = 1 / (2 - rho)
```

- rho = 1.0 (100% purity) => VAF = 1.0
- rho = 0.8 => VAF = 1 / 1.2 = 0.833
- rho = 0.5 => VAF = 1 / 1.5 = 0.667

**Germline heterozygous variant（被刪除的 allele 上有 ALT）：**

```
VAF_expected = (0 * rho + 1 * (1 - rho)) / (rho * 1 + (1 - rho) * 2)
         = (1 - rho) / (2 - rho)
```

- rho = 1.0 => VAF = 0.0（完全消失）
- rho = 0.8 => VAF = 0.2 / 1.2 = 0.167
- rho = 0.5 => VAF = 0.5 / 1.5 = 0.333

**Somatic mutation（在保留的 allele 上，clonal, m=1, CCF=1）：**

```
VAF_expected = (1 * rho) / (rho * 1 + (1 - rho) * 2)
         = rho / (2 - rho)
```

- rho = 1.0 => VAF = 1.0
- rho = 0.8 => VAF = 0.8 / 1.2 = 0.667
- rho = 0.5 => VAF = 0.5 / 1.5 = 0.333

### 3.2 Intermediate AF 的兩種解釋

當觀測到 AF 不在上述期望值時，有兩種主要解釋：

#### (A) Subclonal LOH（部分腫瘤細胞有 deletion）

如果只有一部分腫瘤細胞（fraction = s）發生了 deletion：

```
VAF = ALT_reads / total_reads

對於 germline het（ALT 在被刪除的 allele 上）：
  - 正常細胞: 1 ALT / 2 total (佔比 1-rho)
  - 腫瘤無 deletion: 1 ALT / 2 total (佔比 rho * (1-s))
  - 腫瘤有 deletion: 0 ALT / 1 total (佔比 rho * s)
  
  VAF = [(1-rho)*1 + rho*(1-s)*1 + rho*s*0] / [(1-rho)*2 + rho*(1-s)*2 + rho*s*1]
      = [1 - rho*s] / [2 - rho*s]
```

這會產生介於 「完全 LOH」 和 「無 LOH」 之間的 intermediate VAF。

TITAN (Ha et al., Genome Research, 2014) 的核心模型正是基於此原理：它使用 Hidden Markov Model 將觀測到的 B-allele frequency 分解為三群細胞的混合：正常細胞 (proportion n)、無事件的腫瘤細胞 (proportion (1-n)*s_z)、有事件的腫瘤細胞 (proportion (1-n)*(1-s_z))。不同的 s_z 值對應不同的 subclonal cellular prevalence。

**文獻支持**: Nik-Zainal et al. (Cell, 2012) 使用 Battenberg 演算法首次在乳腺癌 WGS 資料中系統性地發現 subclonal copy number events，包括 subclonal LOH。Battenberg 基於 ASCAT (Van Loo et al., PNAS, 2010) 的原理，但額外加入了 haplotype phasing 和 subclonal CNA 模型。

#### (B) Normal cell contamination（正常細胞稀釋）

純度不完美導致 ALT allele 被正常細胞的 REF reads 稀釋。這與上述公式中 rho < 1.0 時的預期一致。

### 3.3 如何區分 Subclonal LOH vs Purity Effect

| 區分方法 | 原理 | 工具 |
|---------|------|------|
| **多位點一致性** | 若區域內所有 germline het 的 BAF 偏移程度一致 → 更可能是 purity effect 或 clonal LOH；若偏移程度不一致 → 可能有 subclonal 結構 | TITAN, Battenberg |
| **已知 purity 校正** | 先用其他方法估計 purity（如 ABSOLUTE, FACETS），再檢查校正後 BAF 是否仍偏離預期 | ABSOLUTE, FACETS, All-FIT |
| **深度比 (Log-R ratio)** | Deletion LOH 會同時影響 BAF 和 read depth；subclonal deletion 的 depth 下降幅度較小 | ASCAT, Battenberg |
| **相鄰區域比較** | LOH boundary 的 sharpness 可區分 clonal（sharp boundary）vs subclonal（fuzzy boundary） | TITAN |
| **Long-read phasing** | 長讀長可直接 phase variants，觀察是否一個 haplotype 的 reads 一致性減少 | LongPhase-TO |

**關鍵文獻**: Sun et al. (JCO Precision Oncology, 2022) 證實在 tumor-only 定序中，germline heterozygous variant 在 LOH 區域的 VAF 偏移至 60-80%，導致 TO caller 誤判為 somatic。LOHGIC (Khiabanian et al., Journal of Molecular Diagnostics, 2018) 開發了 Bayesian model 來區分 germline-with-LOH vs somatic，使用 Akaike weights 跨越 purity 和 VAF 的信賴區間。

---

## 4. Copy-Neutral LOH (cnLOH, CN=2): AF 與亞克隆

### 4.1 理論期望

cnLOH 中一個 allele 被刪除後另一個被複製，總拷貝數維持 2。

**Germline heterozygous variant：**
- 被保留（複製）的 allele 有 ALT => 腫瘤中 ALT copy = 2, total = 2

```
VAF = (2 * rho + 1 * (1-rho)) / (2 * rho + 2 * (1-rho))
    = (rho + 1) / 2
```

- rho = 1.0 => VAF = 1.0
- rho = 0.8 => VAF = 0.9
- rho = 0.5 => VAF = 0.75

- 被刪除的 allele 有 ALT => 腫瘤中 ALT copy = 0, total = 2

```
VAF = (0 * rho + 1 * (1-rho)) / (2 * rho + 2 * (1-rho))
    = (1 - rho) / 2
```

- rho = 1.0 => VAF = 0.0
- rho = 0.8 => VAF = 0.1
- rho = 0.5 => VAF = 0.25

### 4.2 cnLOH 後的 Somatic Mutation

這是本研究的關鍵場景。在 cnLOH 之後，若在兩個相同拷貝之一上發生 somatic mutation：

```
腫瘤細胞中: ALT 在 1 of 2 copies => m = 1, CN_tumor = 2
VAF = (1 * rho) / (2 * rho + 2 * (1-rho))
    = rho / 2
```

- rho = 1.0 => **VAF = 0.5**
- rho = 0.8 => VAF = 0.4
- rho = 0.5 => VAF = 0.25

**重要結論**: 在 cnLOH 區域中，clonal somatic mutation 的 expected VAF = rho/2，在 100% purity 下為 0.5。這看起來像一般 diploid 區域的 heterozygous variant，但實際上是 cnLOH 區域 post-duplication 的 somatic event。

**文獻支持**: Gerstung et al. (Nature, 2020) 的 PCAWG 分析框架明確使用此原理：duplicated mutations（gain 前發生）出現在兩個 copies 上，VAF 高；non-duplicated mutations（gain 後發生）僅在一個 copy 上，VAF 低（約 clonal VAF 的一半）。CNAqc (Caravagna et al., Genome Biology, 2024) 使用此 VAF peak 結構來驗證 copy number calls 的正確性。

### 4.3 cnLOH 區域 AF 不等於 0, 0.5, 或 1 的意義

若觀測到的 AF 偏離上述三個期望值（0, rho/2, 1），有以下可能：

| 觀測 AF | 可能解釋 |
|---------|---------|
| 0 < AF < rho/2 | **Subclonal somatic mutation**（CCF < 1），只有部分腫瘤細胞帶有此突變 |
| rho/2 < AF < 1 | **Pre-duplication mutation（m=2）但 subclonal**，或 **post-duplication mutation + subclonal LOH** |
| AF 略偏離 rho/2 | 測序噪音、mapping bias、或 **partial cnLOH**（subclonal cnLOH，只有部分細胞經歷 cnLOH） |
| AF 不等於整數 CN 的任何期望值 | **多個 subclone 混合**，或 **tumor heterogeneity** |

**文獻支持**: Tarabichi et al. (Nature Methods, 2021) 的 PCAWG11 subclonal reconstruction benchmark 強調：在 cnLOH 區域中，mutation multiplicity 的推斷是 subclonal reconstruction 的最大挑戰之一。CNVkit documentation (Talevich et al.) 明確指出：「variants cannot be distinguished between different somatic mutational models (heterozygous, under LOH, or under copy-neutral LOH) without knowing specimen's tumor content.」

### 4.4 Subclonal cnLOH 的 AF 模型

若只有 fraction s 的腫瘤細胞經歷 cnLOH（其餘仍為 normal diploid het）：

```
對於 germline het（ALT 在保留 + 複製的 allele 上）：

正常細胞: 1 ALT / 2 total (佔比 1-rho)
腫瘤無 cnLOH: 1 ALT / 2 total (佔比 rho * (1-s))
腫瘤有 cnLOH: 2 ALT / 2 total (佔比 rho * s)

VAF = [(1-rho) + rho*(1-s) + rho*s*2] / [(1-rho)*2 + rho*(1-s)*2 + rho*s*2]
    = [1 + rho*s] / 2
```

- 當 s = 1 (clonal cnLOH): VAF = (1 + rho) / 2（與 4.1 一致）
- 當 s = 0 (no cnLOH): VAF = 0.5（正常 het）
- 當 s = 0.5, rho = 0.8: VAF = (1 + 0.4) / 2 = 0.7

**文獻支持**: Shen & Seshan (Nucleic Acids Research, 2016) 在 FACETS 中明確模型化了 subclonal copy number events，通過 mixture model 將 BAF 和 log-R ratio 分解為 clonal 和 subclonal 成分。

---

## 5. Multi-site Linkage: 多位點一致性作為 Subclonal LOH 證據

### 5.1 核心原理

如果一個 LOH 區域內的多個相鄰 germline heterozygous variants **全部**展現一致的 intermediate BAF 偏移，這是比單個位點更強的 subclonal LOH 證據。

**數學論證**: 對於一個 LOH 事件影響了 n 個相鄰 het SNPs：
- 若是 clonal LOH + purity effect：所有 n 個 SNPs 的 BAF 偏移應高度一致（取決於 local purity 和 CN state）
- 若是 subclonal LOH：BAF 偏移程度取決於 cellular prevalence s，但同樣在所有 n 個 SNPs 上一致
- 若是隨機噪音或 artifact：各 SNP 的 BAF 偏移應不一致

### 5.2 工具如何利用 multi-site linkage

| 工具 | Multi-site 使用方式 | 來源 |
|------|-------------------|------|
| **TITAN** | HMM segmentation，相鄰 SNPs 的 BAF 被 jointly modeled，segment 內所有 SNPs 共享同一 cellular prevalence | Ha et al., Genome Research, 2014 |
| **Battenberg** | 使用 PELT (Pruned Exact Linear Time) changepoint detection，在 haplotyped BAF track 上做 segmentation | Nik-Zainal et al., Cell, 2012 |
| **FACETS** | Joint segmentation of log-R ratio and BAF，使用 multi-resolution approach | Shen & Seshan, Nucleic Acids Research, 2016 |
| **ASCAT** | Allele-specific PCAF (piecewise constant fitting) smooth BAF and log-R ratio jointly | Van Loo et al., PNAS, 2010 |
| **CalicoST** | 空間轉錄體中，聚合同一 haplotype 上相鄰位點的 allele counts 以增強信號 | Chen et al., Nature Methods, 2024 |

### 5.3 Long-read 的獨特優勢

長讀長定序可以直接在單一 read 上觀察多個 het SNPs 的 allele state，提供 **physical linkage** 而非統計推斷：

- 若一條 read 覆蓋 3 個 het SNPs 且全部展現同一 haplotype => 強烈支持 LOH
- 若混合 reads 中部分展現 LOH pattern、部分展現 het pattern => 直接證據支持 subclonal LOH

**文獻支持**: Viswanathan et al. (Cell Genomics, 2024) 使用 ONT 長讀長定序 189 個晚期癌症腫瘤，展示了長讀長 phasing 如何揭示 allele-specific copy number 和 LOH 事件，特別是在 BRCA1/RAD51C 等基因中發現了 haplotype-specific 甲基化與 LOH 的共存。

---

## 6. 甲基化 + 亞克隆：表觀遺傳學在 LOH 亞克隆偵測中的角色

### 6.1 Allele-Specific Methylation (ASM) 在 LOH 區域

#### 正常 ASM 在 LOH 後的表現

**CAMDAC (Tarazona et al., Nature Genetics, 2025)** 是此領域的里程碑工作：
- 建立了 copy number-aware methylation deconvolution 模型
- 在 LOH 區域中，所有攜帶特定 variant allele 的 reads 可確定來自腫瘤細胞
- 利用此特性達到 Pearson r = 0.97 的驗證相關性
- **核心公式**: 將 bulk tumor methylation rate 分解為 tumor 和 normal 成分，考慮 copy number 和 purity

**預期表現**:
| LOH 狀態 | ASM 表現 | 原因 |
|---------|---------|------|
| **Clonal LOH** | ASM 消失（甲基化均質化） | 只剩一個 allele，allelic 差異消除 |
| **Subclonal LOH** | ASM 部分保留 | 有 LOH 的細胞只有一個 allele，無 LOH 的細胞保有兩個 allele；混合後 ASM 信號減弱但不消失 |
| **No LOH** | ASM 正常存在 | 兩個 allele 的甲基化模式不同 |

**文獻支持**: Do et al. (Genome Biology, 2020) 報告癌症中 ASM 比正常組織高 5-9 倍，主要源於 allele-specific loss of methylation (72-76%)。Shoemaker et al. (Genome Research, 2010) 發現 23-37% 的 germline het SNPs 展現 ASM。

#### ASM 作為 Subclonal LOH 的 orthogonal 證據

**如果 BAF 顯示 intermediate LOH，且同一區域的 ASM 也顯示部分保留 => 兩條獨立證據鏈互相支持 subclonal LOH 假說。**

這是因為：
1. BAF 偏移反映 DNA copy number 層面的 allelic imbalance
2. ASM 反映 epigenetic 層面的 allelic difference
3. 兩者應在 clonal LOH 時同時消失，在 subclonal LOH 時同時部分保留

### 6.2 Epi-allele Diversity 作為 Subclone Marker

#### 基本概念

Epi-allele 指同一基因組位置上具有不同甲基化模式的 DNA 分子。在長讀長定序中，每條 read 是一個獨立的 epi-allele observation。

**Epipolymorphism 定義** (Landan et al., Genome Research, 2012):
```
Epipolymorphism = 1 - sum(fi^2)
```
其中 fi 是第 i 種 epi-allele pattern 的頻率。

#### Epi-allele diversity 與 subclone 的關係

| 場景 | Epi-allele diversity | 原因 |
|------|---------------------|------|
| **Clonal, single subclone** | 低（均質） | 所有腫瘤細胞共享相似甲基化模式 |
| **Multiple subclones** | 高（異質） | 不同 subclone 可能有不同甲基化模式 |
| **Subclonal LOH** | 中等 | LOH subclone 有一種模式，non-LOH cells 有另一種 |
| **Stochastic drift (PMD)** | 高但非結構化 | 隨機波動，不形成 discrete clusters |

**文獻支持**:
- Brocks et al. (Nature Reviews Cancer, 2024) 綜述了 epigenomic heterogeneity 在腫瘤演化中的角色，指出 enhancer 區域的甲基化異質性是 subclonal diversification 的主要來源
- Johnstone et al. (Nature Communications, 2019) 發現 PMDs (Partially Methylated Domains) 是甲基化 variation 的主要來源，覆蓋基因組的 24-63%
- **關鍵區分**: promoter CpG islands 的甲基化在 subclonal diversification 中相對穩定，而 enhancer 和 open sea 區域波動更大

#### 長讀長的獨特優勢

傳統短讀長只能觀察少量 CpG 的共甲基化模式；長讀長可以：

1. **觀察 extended methylation haplotype blocks (MHBs)**: Guo et al. (Nature Genetics, 2017) 鑑定了 ~148K MHBs，長讀長可以完整捕捉這些 block 的甲基化模式
2. **同時觀察 genetic + epigenetic variants**: 一條 read 上同時看到 SNV genotype 和 methylation state，直接 link genotype to epitype
3. **區分 heterogeneous vs homogeneous hypomethylation**: 長讀長可以判斷一個區域的低甲基化是均質的（所有分子一致）還是異質的（不同分子有不同模式）

**文獻支持**: Gershman et al. (Nature Methods, 2022) 展示了 ONT 直接測序同時偵測 5mC 和 5hmC 的能力。Viswanathan et al. (Cell Genomics, 2024) 在 189 個晚期癌症中使用 ONT 揭示了 haplotype-specific methylation landscapes。

### 6.3 Long-read Phasing + Methylation 用於 Subclone Detection

#### 方法論框架

```
Step 1: Variant calling + Phasing
  - 使用 germline het SNPs 建立 haplotype phase blocks
  - Long-read phasing (LongPhase, WhatsHap) 可產生 Mb-scale phase blocks

Step 2: Haplotagging
  - 將每條 read 標記為 HP1 或 HP2

Step 3: Allele-specific methylation analysis
  - 比較 HP1 reads 和 HP2 reads 的甲基化模式
  - 在 LOH 區域: 兩個 HP 的 read 數量不平衡

Step 4: Subclone inference
  - 若 LOH 是 clonal: HP1 和 HP2 的甲基化應相似（都來自同一 allele）
  - 若 LOH 是 subclonal: HP1 和 HP2 的甲基化可能不同（部分來自 retained allele，部分來自 non-LOH cells 的兩個 alleles）
```

#### 已發表的工具與方法

| 工具/方法 | 功能 | 平台 | 來源 |
|---------|------|------|------|
| **CAMDAC** | Copy number-aware methylation deconvolution | RRBS/WGBS (短讀長) | Tarazona et al., Nature Genetics, 2025 |
| **MethPhaser** | 利用甲基化信號延伸 SNV-based phasing | ONT + PacBio | Youn et al., Nature Communications, 2024 |
| **NANOME** | Haplotype-aware allele-specific consensus methylation detection | ONT | Li et al., 2025 |
| **EVOFLUx** | 利用 fluctuating CpGs 作為分子鐘推論腫瘤演化 | Array/WGBS | Gabbutt et al., Nature, 2025 |
| **nanoNOMe** | Nanopore + NOMe-seq 結合甲基化和染色質可及性 | ONT | Lee et al., 2020 |
| **MethHaplo** | 結合 ASM 和 SNP 做 haplotype region identification | 通用 | Xu et al., BMC Bioinformatics, 2020 |
| **LongPhase-TO** | Tumor-only phasing + LOH detection + purity estimation | ONT | Chen (CCU Thesis), 2025 |

**特別值得注意**: 目前尚無專門針對「在 LOH 區域中結合 AF + methylation 偵測 subclone」的工具。CAMDAC 最接近此目標，但它使用短讀長 bisulfite sequencing。將 CAMDAC 的概念框架應用於 ONT 長讀長資料是一個尚未充分探索的方向。

---

## 7. 既有工具/方法綜覽

### 7.1 AF-based Subclonal LOH Detection Tools

| 工具 | 原理 | Subclonal 支援 | LOH 支援 | 長讀長支援 | 來源 |
|------|------|---------------|---------|----------|------|
| **TITAN** | HMM on BAF + log-R | 是（cellular prevalence clusters） | 是 | 否 | Ha et al., Genome Research, 2014 |
| **Battenberg** | ASCAT + changepoint + subclonal mixture | 是 | 是 | 否 | Nik-Zainal et al., Cell, 2012 |
| **FACETS** | Joint BAF + log-R segmentation | 是（cellular fraction） | 是 | 否 | Shen & Seshan, Nucleic Acids Research, 2016 |
| **ABSOLUTE** | BAF + log-R + mutation VAF | 是 | 是 | 否 | Carter et al., Nature Biotechnology, 2012 |
| **ASCAT** | Allele-specific segmentation | 否（clonal only） | 是 | 否 | Van Loo et al., PNAS, 2010 |
| **PyClone** | Bayesian clustering of CCF | 是（需 CNA input） | 間接 | 否 | Roth et al., Nature Methods, 2014 |
| **LOHGIC** | Bayesian germline vs somatic in LOH | 否 | 是 | 否 | Khiabanian et al., JMD, 2018 |
| **All-FIT** | VAF-based purity + CCF estimation | 是 | 間接 | 否 | Loh et al., Bioinformatics, 2020 |
| **CNAqc** | Clonal/subclonal CNA validation | 是 | 是 | 否 | Caravagna et al., Genome Biology, 2024 |
| **LongPhase-TO** | Graph-based somatic calling + LOH | 否（clonal） | 是 | **是 (ONT)** | Chen (CCU), 2025 |

### 7.2 Methylation-based Subclone Detection Tools

| 工具 | 原理 | Subclone 支援 | 長讀長支援 | 來源 |
|------|------|---------------|----------|------|
| **CAMDAC** | CN-aware methylation deconvolution | 是（methylation-based subclone phylogeny） | 否 | Tarazona et al., Nature Genetics, 2025 |
| **EVOFLUx** | Fluctuating CpG 分子鐘 | 是（growth rate, epimutation rate） | 否 | Gabbutt et al., Nature, 2025 |
| **MethPhaser** | Methylation-based phasing | 間接 | **是 (ONT)** | Youn et al., Nature Communications, 2024 |
| **NANOME** | Haplotype-aware ASM detection | 間接 | **是 (ONT)** | Li et al., 2025 |
| **InterSubMod** | Read-level methylation distance + clustering | 是（ISM 特徵） | **是 (ONT)** | 本專案 |

### 7.3 Combined AF + Methylation Approaches (Gap in Literature)

**目前文獻中的缺口**: 將 AF-based subclonal LOH detection 與 methylation-based subclone evidence 整合的工具極少。最接近的是：

1. **CAMDAC**: 同時使用 copy number (from WGS) 和 methylation (from RRBS/WGBS)，但它們來自不同的定序平台
2. **Viswanathan et al. (Cell Genomics, 2024)**: 在 ONT 上同時分析 structural variants + methylation landscapes，但未專門針對 subclonal LOH + methylation 的整合分析
3. **本專案 (InterSubMod)**: 在 ONT 上分析 read-level methylation patterns，具有潛力將 AF 和 methylation 整合用於 subclonal LOH detection

---

## 8. 核心發現整合

### 8.1 共識觀點

1. **Intermediate AF 在 LOH 區域確實可以指示 subclonal LOH** — 多個工具（TITAN, Battenberg, FACETS）的核心假設即是此原理（來源: Ha et al. 2014, Nik-Zainal et al. 2012, Shen & Seshan 2016）

2. **但 intermediate AF 也可能來自 normal contamination** — 必須先估計 purity 才能區分（來源: Carter et al. 2012, Khiabanian et al. 2018, All-FIT）

3. **cnLOH 後的 somatic mutation 期望 AF = rho/2** — 這是 PCAWG timing 分析的基礎（來源: Gerstung et al. 2020）

4. **Multi-site 一致性是更強的 subclonal LOH 證據** — 所有 CNA 工具都使用 segmentation 而非單點分析（來源: 所有上述工具）

5. **ASM 在 clonal LOH 後消失** — CAMDAC 利用此特性做驗證（來源: Tarazona et al. 2025）

6. **長讀長在此領域有獨特優勢但尚未充分開發** — 同時提供 phasing + methylation + variant calling（來源: Viswanathan et al. 2024, MethPhaser 2024）

### 8.2 衝突觀點

| 觀點 A | 觀點 B | 可能原因 |
|--------|--------|----------|
| Epi-allele diversity 可用於 subclone detection (Brocks 2024) | Epi-allele diversity 主要反映 stochastic drift 而非 subclonal structure (本專案 O11 結論) | Context 不同：promoter CpG islands vs PMD/open sea；工具解析度不同 |
| ASM 可區分 germline vs somatic (理論預測) | ASM 在實際癌症中高度異質，難以作為分類器 (Do et al. 2020, 本專案結論) | Cancer-driven ASM 是 stochastic 的，與 germline-driven ASM 混淆 |
| AF alone 足以偵測 subclonal LOH (TITAN, Battenberg) | 需要 methylation 等 orthogonal 證據 (CAMDAC 方向) | AF-based 方法在低 cellularity 和 noisy data 中可能不足 |

---

## 9. 資料來源評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| Gerstung et al. 2020, Nature (PCAWG) | 期刊論文 | 高 | 2,658 tumors, VAF-based timing |
| Tarazona et al. 2025, Nature Genetics (CAMDAC) | 期刊論文 | 高 | TRACERx 122 samples |
| Ha et al. 2014, Genome Research (TITAN) | 期刊論文 | 高 | 標準工具 |
| Nik-Zainal et al. 2012, Cell (Battenberg) | 期刊論文 | 高 | 乳腺癌 WGS 標竿 |
| Gabbutt et al. 2025, Nature (EVOFLUx) | 期刊論文 | 高 | 1,976 lymphoid cancers |
| Viswanathan et al. 2024, Cell Genomics | 期刊論文 | 高 | 189 tumors, ONT |
| Dentro et al. 2017, CSH Perspectives | 綜述 | 高 | Subclonal reconstruction principles |
| Brocks et al. 2024, Nature Reviews Cancer | 綜述 | 高 | Epigenomic heterogeneity |
| Khiabanian et al. 2018, JMD (LOHGIC) | 期刊論文 | 中-高 | 臨床驗證 |
| Sun et al. 2022, JCO Precision Oncology | 期刊論文 | 中-高 | TO LOH germline FP |
| CNVkit documentation | 工具文件 | 中 | 實用但非 peer-reviewed |
| Biostars 討論 | 社群 | 低-中 | 實務經驗但需驗證 |

---

## 10. 與 InterSubMod 的關聯和建議行動

### 10.1 直接相關

本專案已確認的結論（來自 docs/reports/research_landscape/）：
- Self-phasing circular dependency 影響 TO scenario 的 HP 分配
- LOH 區域的甲基化趨向均質（單 allele pattern）
- ASM 在 FP vs TP 區分中效果有限（AUC < 0.64）

### 10.2 新的研究方向可能性

1. **Subclonal LOH Detection via Long-Read**: InterSubMod 可以直接觀察 read-level 的 haplotype + methylation，理論上可以偵測 subclonal LOH — 但需要先解決 self-phasing 問題（使用 normal BAM 或 PON-only phasing）

2. **VAF + Methylation Integration**: 在已知 LOH 區域中，結合 variant AF 和 read-level methylation pattern 可能提供更精確的 subclonal structure 估計。這是文獻中的明確 gap。

3. **cnLOH Post-Duplication Mutation Detection**: InterSubMod 可以利用 long-read phasing 觀察 cnLOH 區域中 AF=0.5 的 somatic mutations 的甲基化 context，這是短讀長工具無法做到的。

### 10.3 注意事項

- 本專案 O11 已確認 epipolymorphism 的 AUC 受 n_reads confound 影響 — 任何基於 epi-allele diversity 的分析必須控制 read depth
- 本專案已確認 HP-free features 在 TO scenario 中 AUC < 0.56 — subclonal LOH detection 可能需要 normal BAM 才能有效
- CAMDAC 的 copy number-aware deconvolution 是最可能 translatable 到長讀長的方法框架

---

## 附錄: 關鍵公式速查表

### A1. Deletion LOH (CN_tumor = 1)

| 場景 | VAF 公式 | VAF (rho=1.0) | VAF (rho=0.8) | VAF (rho=0.5) |
|------|---------|---------------|---------------|---------------|
| Germline, ALT retained | 1/(2-rho) | 1.00 | 0.83 | 0.67 |
| Germline, ALT deleted | (1-rho)/(2-rho) | 0.00 | 0.17 | 0.33 |
| Somatic, clonal, m=1 | rho/(2-rho) | 1.00 | 0.67 | 0.33 |
| Somatic, subclonal, CCF=0.5 | 0.5*rho/(2-rho) | 0.50 | 0.33 | 0.17 |

### A2. Copy-Neutral LOH (CN_tumor = 2)

| 場景 | VAF 公式 | VAF (rho=1.0) | VAF (rho=0.8) | VAF (rho=0.5) |
|------|---------|---------------|---------------|---------------|
| Germline, ALT duplicated (m=2) | (rho+1)/2 | 1.00 | 0.90 | 0.75 |
| Germline, ALT deleted | (1-rho)/2 | 0.00 | 0.10 | 0.25 |
| Somatic, pre-cnLOH (m=2) | rho | 1.00 | 0.80 | 0.50 |
| Somatic, post-cnLOH (m=1) | rho/2 | 0.50 | 0.40 | 0.25 |
| Somatic, post-cnLOH, subclonal CCF=0.5 | rho/4 | 0.25 | 0.20 | 0.13 |

### A3. Normal Diploid (CN_tumor = 2, no LOH)

| 場景 | VAF 公式 | VAF (rho=1.0) | VAF (rho=0.8) | VAF (rho=0.5) |
|------|---------|---------------|---------------|---------------|
| Germline het | 0.5 | 0.50 | 0.50 | 0.50 |
| Somatic, clonal | rho/2 | 0.50 | 0.40 | 0.25 |
| Somatic, subclonal CCF=0.5 | rho/4 | 0.25 | 0.20 | 0.13 |

---

**文檔版本**: v1.0
**建立者**: AI Agent (Researcher)
**最後更新**: 2026-04-14
