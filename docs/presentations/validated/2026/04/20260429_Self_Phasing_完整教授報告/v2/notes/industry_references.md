# 業界對照文獻清單（PPT v2 · S26 用）

> 2026-04-29 WebSearch 取得；用於 Section 5「為何改動有效 + 業界對照」slide 內容支撐。

## 1. 直接相關（同一研究組／同一工具家族）

### 1.1 LongPhase-S（2025 bioRxiv，CCU-Bioinformatics-Lab）★★★

- **連結**: <https://www.biorxiv.org/content/10.1101/2025.11.20.689492v1.full>
- **GitHub**: <https://github.com/CCU-Bioinformatics-Lab/longphase-s>
- **作者所屬**: 中正大學生資實驗室（與 LongPhase 原論文 Lin J-H 同實驗室）
- **核心做法**: "anchoring somatic reads to their parental germline haplotypes"（將 somatic reads 錨定到 germline haplotype）→ 解開 germline 與 somatic haplotype
- **重要對照**: 
  - LongPhase-S 是 **tumor-normal paired** 模式（有 normal）
  - 我們使用的 **longphase-to / longphase-to-mod** 是 tumor-only 模式（無 normal）
  - 兩者來自同一實驗室；本實驗室是 longphase-to 的下游 fork
- **F1 改善**: ClairS / DeepSomatic 上 SNV +4.5% / +1.2%、indel +7.1% / +0.5%
- **與我們的關係**: ✅ **同實驗室相鄰工作**（CCU-Bioinformatics-Lab）— LongPhase-S 提供 paired 場景的「somatic 錨在 germline scaffold」實作；本工作是 **TO + PON 條件下的本地實作與修補**，並非業界共識（WhatsHap、HapCUT2 等業界主流工具實際上不直接處理 tumor-only 模式）

### 1.2 LongPhase-TO（GitHub repo，無正式論文）

- **連結**: <https://github.com/CCU-Bioinformatics-Lab/longphase-to>
- **狀態**: 本實驗室 fork 的上游源頭；本 PPT 描述的所有 4-commit 修補（V2b/V3F/INDEL guard/V5）都建立於此 repo

## 2. 平行 caller / phasing 工具

### 2.1 DeepSomatic（Nature Biotechnology 2025）★★

- **連結**: <https://www.nature.com/articles/s41587-025-02839-x>
- **bioRxiv**: <https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1>
- **重點**: 深度學習 somatic SNV/indel caller，支援 short-read + long-read，**含 tumor-only 模式**
- **與我們的關係**: 平行於 ClairS-TO 的 caller；可作為 future cross-caller validation 候選

### 2.2 WhatsHap 2.8（2024-2025 持續更新）

- **連結**: <https://whatshap.readthedocs.io/en/latest/>
- **重點**: 業界主流 germline phasing 工具，支援 Illumina / PacBio / Nanopore；不直接處理 tumor-only self-phasing
- **與我們的關係**: ❌ 不適用 tumor-only；**證明業界缺乏 tumor-only 標準解，本實作填補此 gap**

### 2.3 PEPPER-Margin-DeepVariant（2021 Nature Methods）

- **連結**: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8571015/>
- **重點**: Haplotype-aware nanopore variant calling
- **與我們的關係**: germline 為主，tumor-only 應用有限

### 2.4 MethPhaser（Nature Comm 2024）

- **連結**: <https://www.nature.com/articles/s41467-024-49588-0>
- **重點**: 用甲基化模式做 long-read phasing
- **與我們的關係**: 互補方向（甲基化幫助 phasing）；本實作目前是 phasing 修對後再算甲基化，未來可探索 cross-direction

## 3. 癌症 long-read phasing 早期文獻

### 3.1 Phasing analysis of lung cancer genomes（Nature Comm 2022）

- **連結**: <https://www.nature.com/articles/s41467-022-31133-6>
- **PMC**: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9203510/>
- **重點**: 肺癌長讀長 phasing 分析，揭示腫瘤異質性與 LOH region 的 phase block 限制
- **與我們的關係**: 早期癌症 long-read phasing 開創性工作；提供「LOH region 限制 phasing」的先前依據

### 3.2 Unraveling cancer through long-read sequencing（Genome Research 2025）

- **連結**: <https://genome.cshlp.org/content/early/2025/03/19/gr.280041.124.full.pdf>
- **重點**: 2025 review 文章；強調 long-read 癌症基因組的隱藏複雜性
- **與我們的關係**: 提供 cancer long-read 領域 2025 整體 review 背景

### 3.3 Somatic mutation phasing in multiple myeloma（PubMed 2024）

- **連結**: <https://pubmed.ncbi.nlm.nih.gov/39149342/>
- **重點**: Linked-reads 在多發性骨髓瘤 somatic phasing 應用
- **與我們的關係**: 短讀／linked-reads 對照，補強「somatic phasing 跨樣本一致性」需求

## 4. 業界共識總結（給 S26 slide 用）

| 工具 | 模式 | Tumor-only 處理 | 與本實作的關係 |
|------|------|---------------|---------------|
| LongPhase（Lin 2022 *Bioinformatics*）| Germline | N/A | 上游基底 |
| **LongPhase-S（bioRxiv 2025）** | **Tumor-Normal Paired** | N/A | **Anchoring 概念來源** |
| **longphase-to / longphase-to-mod** | **Tumor-Only** | **PON 替代 normal + 4-commit 修補（本工作）** | **本 PPT 焦點** |
| WhatsHap 2.8 | Germline | 不支援 | 業界 gap |
| DeepSomatic（Nature Biotech 2025）| Caller（含 TO 模式）| caller 層 | 未來 cross-caller |
| MethPhaser（Nature Comm 2024）| Methylation-based phasing | N/A | 互補方向 |

## 5. 一句話總結（直接寫進 S2 / S21 speaker note）

> 「LongPhase-S（2025 bioRxiv，paired 模式）提供 **somatic-to-germline anchoring** 的 paired 場景相鄰參考；本工作是 **TO + PON 條件下的本地實作與修補**。當 normal 不可得時，本實作以 PON 替代 normal 作 anchor，修補了 LongPhase-TO 在 `--pon-only-phasing` 啟用後集中暴露的 3 層 tag-side bug，使 read-level concordance 從 82.2% 提升至 90.5%（+8.3 pp，clean PS blocks 全基因組）。WhatsHap、HapCUT2 等業界 germline phasing 主流工具不直接處理 tumor-only 模式；本實作不宣稱『業界共識』，而是同實驗室相鄰工作的 TO 變體。」

**口徑校準**（避免 over-claim）：
- ❌ 不寫：「業界共識」「業界標準替代」
- ✅ 改寫：「同實驗室相鄰工作」「TO + PON 本地實作」「填補 tumor-only 在公開工具中的 gap」

## 6. Sources

- [LongPhase-S bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.11.20.689492v1.full)
- [LongPhase-S GitHub](https://github.com/CCU-Bioinformatics-Lab/longphase-s)
- [LongPhase-TO GitHub](https://github.com/CCU-Bioinformatics-Lab/longphase-to)
- [DeepSomatic Nature Biotech 2025](https://www.nature.com/articles/s41587-025-02839-x)
- [DeepSomatic bioRxiv](https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1)
- [Lung cancer phasing Nature Comm 2022](https://www.nature.com/articles/s41467-022-31133-6)
- [Cancer long-read review Genome Research 2025](https://genome.cshlp.org/content/early/2025/03/19/gr.280041.124.full.pdf)
- [PEPPER-Margin-DeepVariant PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8571015/)
- [WhatsHap 2.8 docs](https://whatshap.readthedocs.io/en/latest/)
- [MethPhaser Nature Comm 2024](https://www.nature.com/articles/s41467-024-49588-0)
- [Multiple myeloma phasing PubMed 2024](https://pubmed.ncbi.nlm.nih.gov/39149342/)
