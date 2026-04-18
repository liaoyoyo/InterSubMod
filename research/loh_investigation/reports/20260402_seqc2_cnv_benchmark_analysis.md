<!--
建立時間: 2026-04-02 17:00
目標: SEQC2 CNV Benchmark 論文與數據集整理，用於驗證 LongPhase-TO LOH.bed 與 ISM LOH 判定
處理範圍: PMC11188507 論文摘要 + Zenodo v4 數據集分析 + 與我們 LOH 研究的對應
關聯檔案:
  - research/loh_investigation/data/seqc2_cnv_benchmark_v4/
  - docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md
-->

# SEQC2 CNV Benchmark — 論文與數據集整理

## 1. 論文資訊

| 項目 | 內容 |
|------|------|
| 標題 | Evaluation of somatic copy number variation detection by NGS technologies and bioinformatics tools on a hyper-diploid cancer genome |
| 期刊 | Genome Biology, 25:163 |
| 日期 | 2024-06-20 |
| DOI | 10.1186/s13059-024-03294-8 |
| PMC | PMC11188507 |
| 作者 | Daniall Masood, Luyao Ren, Cu Nguyen, ... Wenming Xiao, Daoud M Meerzaman |

## 2. 核心內容

### 2.1 研究目標

使用 HCC1395（三陰性乳癌細胞株）及其配對正常 HCC1395BL，系統性評估 6 個 CNV calling 工具的準確性、敏感度與再現性。結合三種正交驗證技術（microarray、Bionano optical mapping）建立 high-confidence CNV benchmark。

### 2.2 樣本

| 樣本 | 描述 | 特性 |
|------|------|------|
| **HCC1395** | Triple-negative breast cancer cell line (ATCC) | **Hyper-diploid, ploidy = 2.85** |
| **HCC1395BL** | Matched normal B-lymphoblastoid | Diploid |

**關鍵：這是我們 InterSubMod 使用的同一個 HCC1395 細胞株。**

### 2.3 定序平台

**本論文僅使用 Illumina 短讀長定序**（HiSeq + NovaSeq），跨 7 個中心、21 WGS + 12 WES replicates。

- 未使用 ONT 或 PacBio 長讀長
- 正交驗證：Affymetrix CytoScan HD、Illumina CytoSNP-850K、Bionano Saphyr

### 2.4 評估的 CNV Callers

| Caller | 版本 | Ploidy-aware | 表現 |
|--------|------|-------------|------|
| ascatNgs | v4.2.1 | Yes | 最佳之一 |
| CNVkit | v0.9.1 | 需手動校正 | 最佳之一（需 --center-at -0.51） |
| DRAGEN | v4.0.x | Yes | 最佳之一 |
| FACETS | v0.6.0 | Yes | 中等，LOH 一致性高 |
| HATCHet | v1.0.4 | Yes | 再現性最差 |
| Control-FREEC | v11.6 | Yes | 不穩定 |

### 2.5 Truth Set 建構方法

三階段流程：

1. **Data Partitioning**：BedTools multiinter 從 126 WGS call sets 建立不重疊區間
2. **Confidence Scoring**：跨 5 組中心 × 6 callers 的再現性評分（0-3 per group → 累積 → strong/medium/weak/neutral）
3. **正交驗證**：Strong NGS calls + 有 ≥2 正交支持的 medium/weak calls → high-confidence benchmark

### 2.6 Key Performance Metrics

- **F1 Score**：best callers (ascatNgs, CNVkit, DRAGEN) ~0.85-0.95 for gains
- **Ploidy 是最大影響因素**：7 個 WGS runs 誤判 ploidy=5（vs 真實 2.85）→ 大量假 gain/loss
- **Tumor purity < 20%**：所有 callers 大幅失效
- **FFPE vs Fresh**：precision 顯著下降（p=1.49e-07）

---

## 3. SEQC2 CNV Benchmark 數據集（Zenodo v4）

### 3.1 基本資訊

| 項目 | 內容 |
|------|------|
| DOI | 10.5281/zenodo.14619054 |
| 版本 | v4 (2024-06-20) |
| Genome Build | **GRCh38 (hg38)** — 與我們使用的 reference 一致 |
| License | CC BY 4.0 |

### 3.2 下載檔案

| 檔案 | 大小 | 內容 |
|------|------|------|
| `ngs_benchmark_cnv_gain_loss_loh.bed` | 19.4 KB | **核心：688 regions (340 gain + 315 loh + 33 loss)** |
| `ngs_benchmark_cnv_gain_cn.bed` | 39.3 KB | Gain regions with copy number |
| `ngs_benchmark_cnv_loss_cn.bed` | 1.6 KB | Loss regions with copy number |
| `cnv_benchmark_calls.vcf` | 176.1 KB | VCF 格式 benchmark calls (1452 entries) |
| `exclusion.bed` | 14.8 KB | 排除區域（633 regions） |
| `cnv_gain_cn_median.txt` | 26.8 KB | Gain regions median CN |
| `cnv_loss_cn_median.txt` | 1.4 KB | Loss regions median CN |

### 3.3 LOH Benchmark 統計

| CNV Type | Regions | Total bp | Coverage |
|----------|---------|----------|----------|
| **LOH** | **315** | **1,490.4 Mb** | **49.2% genome** |
| Gain | 340 | 1,525.5 Mb | 50.3% |
| Loss | 33 | 87.9 Mb | 2.9% |

**LOH 佔近一半基因組** — 這證實 HCC1395 是一個 LOH 非常廣泛的細胞株。

### 3.4 LOH Per Chromosome

| Chr | Regions | Mb | 佔該 chr 比例 |
|-----|---------|-----|--------------|
| chr3 | 32 | 157.7 | ~80% |
| chr8 | 26 | 139.3 | ~96% |
| chr5 | 21 | 130.2 | ~72% |
| chr2 | 39 | 122.7 | ~51% |
| chr11 | 7 | 114.6 | ~85% |
| chr1 | 37 | 112.9 | ~45% |
| chr12 | 5 | 102.6 | ~77% |
| chr10 | 31 | 104.7 | ~78% |

### 3.5 VCF 中的 LOH 類型

| 類型 | VCF entries | 說明 |
|------|------------|------|
| GAINLOH | 249 | Copy gain + LOH（CN=3 且雜合度消失） |
| LOH | 196 | Copy-neutral LOH（CN=2 但雜合度消失） |
| GAIN | 957 | Copy gain（有或無 LOH） |
| LOSS | 50 | Copy loss |

---

## 4. 與我們 LOH 研究的對應

### 4.1 LongPhase-TO LOH.bed vs SEQC2 Benchmark

| 指標 | LongPhase-TO LOH.bed | SEQC2 Benchmark |
|------|---------------------|-----------------|
| HCC1395 LOH 覆蓋 | **53.9%** (1,632 Mb) | **49.2%** (1,490 Mb) |
| LOH regions | 1,094 | 315 |
| Region 平均大小 | 1.49 Mb | 4.73 Mb |
| Genome build | GRCh38 | GRCh38 |

**LongPhase-TO 的 53.9% 與 SEQC2 的 49.2% 相近** — 但 LongPhase-TO 有更多碎片化的小 regions（1094 vs 315），可能是 self-phasing 造成邊界不穩定。

### 4.2 意義

1. **SEQC2 確認 HCC1395 確實有 ~49% 基因組是 LOH** — LongPhase-TO 的高覆蓋率不完全是 artifact，細胞株本身就有大量結構性 LOH
2. **但 LongPhase-TO 多出 4%**（53.9% vs 49.2%）— 這 4% 差距可能是 self-phasing 造成的額外假 LOH
3. **SEQC2 LOH benchmark 可作為 ground truth** 驗證 LongPhase-TO LOH.bed 的準確性
4. **Chr8 LOH 驗證**：SEQC2 確認 chr8 有 139.3 Mb LOH（~96% 覆蓋），與我們之前觀察到的 chr8 LOH hotspot 一致

### 4.3 注意事項

- **SEQC2 使用 Illumina 短讀長**建立 benchmark，不是 ONT 長讀長
- **Ploidy 2.85** — HCC1395 是 hyper-diploid，CN=3 的 gain 實際上接近 "normal" for this cell line
- **LOH 定義差異**：SEQC2 的 LOH 是 copy-number aware（CN=2 loss of het 或 CN=3 gain with LOH），LongPhase-TO 的是 phased genotype ratio based
- **GAINLOH (249 entries) 佔多數** — 大部分 LOH 是伴隨 copy gain 的，不是 copy-neutral

### 4.4 後續可用分析

1. **BedTools intersect**: LongPhase-TO LOH.bed vs SEQC2 LOH benchmark → 量化 overlap
2. **ISM Potential_LOH 驗證**: ISM 判 LOH 的 regions 是否落在 SEQC2 benchmark 內
3. **False LOH 區域識別**: LongPhase-TO 判 LOH 但 SEQC2 判非 LOH 的區域 → 這些可能是 self-phasing artifacts
4. **Per-chromosome LOH 校準**: SEQC2 提供逐染色體的 LOH 邊界，可精確對比

---

## 5. 數據存放位置

```
research/loh_investigation/data/seqc2_cnv_benchmark_v4/
├── ngs_benchmark_cnv_gain_loss_loh.bed   # 核心：688 regions (gain/loss/loh)
├── ngs_benchmark_cnv_gain_cn.bed         # Gain + CN
├── ngs_benchmark_cnv_loss_cn.bed         # Loss + CN
├── cnv_benchmark_calls.vcf               # VCF 格式完整 benchmark
├── exclusion.bed                         # 排除區域
├── cnv_gain_cn_median.txt                # Gain median CN
└── cnv_loss_cn_median.txt                # Loss median CN
```

## 6. 原始數據來源

| 資源 | URL |
|------|-----|
| 論文 | https://pmc.ncbi.nlm.nih.gov/articles/PMC11188507/ |
| Zenodo v4 | https://zenodo.org/records/14619054 |
| BAM files | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/ |
| Illumina array | NCBI GEO (ref 26 in paper) |
