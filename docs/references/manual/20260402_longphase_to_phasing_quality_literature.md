<!--
建立時間: 2026-04-02 12:00
目標: LongPhase-TO phasing 品質問題與改善方案文獻調查
處理範圍: LongPhase 演算法機制、TO vs Paired phasing 差異、self-phasing circular dependency 根因、相關工具改善方案、phase block 品質評估
關聯檔案:
  - docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md
  - /big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-to.md
  - /big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-s.md
  - /big8_disk/liaoyoyo2001/Knowledge/06_workflows/phasing-workflow.md
-->

# LongPhase-TO Phasing 品質問題與改善方案：文獻調查報告

## 1. 搜尋概述

- **關鍵詞**: LongPhase, tumor-only phasing, self-phasing, circular dependency, somatic-aware phasing, phase block quality, switch error, WhatsHap, HapCUT2, MethPhaser, Octopus, SAVANA, Severus, Wakhan, PhaseME
- **搜尋時間**: 2026-04-02
- **資料來源數**: 18（6 本地知識庫 + 2 本地原始碼 + 5 論文 PDF + 5 網路搜尋）
- **涵蓋工具**: LongPhase / LongPhase-S / LongPhase-TO / WhatsHap / HapCUT2 / MethPhaser / Octopus / SAVANA / Severus / Wakhan / PhaseME

---

## 2. LongPhase 演算法核心機制

### 2.1 兩階段 Phasing（LongPhase 原版，Germline）

根據 Lin et al. (Bioinformatics, 2022) 原始論文：

**Stage 1: Initial Phasing（局部圖）**
- SNPs 和 SVs 作為有向無環圖（DAG）的頂點，按基因組座標排序
- 每個 heterozygous SNP/SV 建立兩個頂點（major allele 和 minor allele）
- 邊的權重 = 覆蓋兩個 variant 的 reads 數量
- 求解 longest pairs of disjoint paths → 初步 phase blocks
- 時間複雜度 O(V)，因為 DAG 的特殊結構允許 greedy 演算法

**Stage 2: Block Re-phasing（長距離連接）**
- 初始 phased blocks 成為新圖的頂點
- 長 reads 的遠距離 linkage 作為邊（跨 block 的 reads）
- 再次求解 longest pairs of disjoint paths → 合併成更大的 blocks
- 此階段去除因假 SNP/SV call、spurious alignment、sequencing error 造成的斷點

**關鍵參數**:
- `--connectAdjacent` (default 35): 每個 variant 與相鄰 N 個 variants 建立邊（germline 密度高）
- `--distance` (default 300,000 bp): 兩個 variant 間最大距離
- `--snpConfidence` (default 0.75): SNP allele 分配到 haplotype 的信心門檻
- `--readConfidence` (default 0.65): Read 分配到 haplotype 的信心門檻

**效能**: ONT 30x WGS ~2 分鐘，N50 = 25 Mbp（ONT ultra-long reads），switch error rate ~0.11-0.16%

### 2.2 LongPhase-TO 的 Phasing 流程（原始碼確認）

透過閱讀 LongPhase-TO 原始碼 (`PhasingProcess.cpp`, `PhasingGraph.cpp`) 確認實際流程：

```
Step 1: 載入 SNP VCF（包含 caller 的所有 variants：PASS + NonSomatic + LowQual + ...）
Step 2: 載入 PON，呼叫 snpFile.setGermline() 標記 germline variants
Step 3: addEdge() — 將 ALL reads 的 ALL variants 加入 phasing graph
         *** 此時 somatic 和 germline variants 全部參與圖建構 ***
Step 4: somaticCalling() 或 tagSomatic() — 在圖建構之後標記 somatic variants
Step 5: phasingProcess() → edgeConnectResult() + readCorrection()
         使用已包含 somatic variants 的圖進行 haplotype assignment
Step 6: （高純度樣本 >0.95）convertNonGermlineToSomatic() 後重跑 phasingProcess()
```

**關鍵發現**：`addEdge()` 在 `somaticCalling()` 之前執行，意味著 somatic variants 在被識別之前就已經作為 phasing anchor 參與了圖的構建。這是 self-phasing circular dependency 的根因。

### 2.3 LongPhase-TO 特有的 Somatic 參數

| 參數 | 預設值 | 意義 |
|------|--------|------|
| `--somaticConnectAdjacent` | 6 | Somatic variant 連接相鄰 SNP 數量（germline 為 35） |
| `--disable-calling` | false | 停用 LongPhase 內建 somatic calling |
| `--disable-pon-tag` | false | 停用從 VCF FILTER 讀取 PON tag |

**`somaticConnectAdjacent=6` 的意義**：限制 somatic variant 在圖中的連接範圍（僅連接最近 6 個 variants），相比 germline 的 35 個更保守，但仍允許 somatic variant 參與 phasing scaffold。

---

## 3. TO vs Paired Phasing 的根本差異

### 3.1 Phasing Scaffold 組成

| 方面 | LongPhase-S (Paired) | LongPhase-TO |
|------|---------------------|-------------|
| **Scaffold 來源** | Normal sample 的 germline SNPs | Tumor VCF 所有 variants（PON-filtered） |
| **Somatic 參與 scaffold** | 否（somatic 不進入 phase graph） | 是（89.9% TP 有 PS assignment） |
| **Germline 來源** | 直接從 normal BAM 呼叫（高品質） | PON database 間接推斷（有遺漏） |
| **Variant 密度** | 高（~3M germline het SNPs） | 低（主要依賴 NonSomatic 標記的 variants） |
| **Phase block 特徵** | 較多、較短（密集 SNP → 更多 switch point） | 較少、超長（稀疏 → mega-blocks） |

### 3.2 Self-Phasing Circular Dependency（已驗證）

根據我們的因果鏈驗證報告（20260402_longphase_to_vs_s_causal_chain_report_01.md）：

```
LongPhase-TO 的循環依賴：

  Somatic variant 進入 VCF
       ↓
  LongPhase-TO addEdge(): 所有 variants（含 somatic）建構 phasing graph
       ↓  ← 循環依賴點
  Somatic variant 作為 anchor 將其自身的 ALT reads 定向到同一個 HP
       ↓
  Haplotag: Reads 根據 phased VCF 分配 HP tag
       ↓
  結果: Somatic ALT reads 系統性偏向一個 HP → 虛假 HP imbalance → 虛假 LOH
```

**量化證據**:
- 31.2% 的 TO LOH 由 self-phasing 造成
- 移除 self-phasing 後 62% LOH 消失（AF 0.1-0.8 範圍內接近 100%）
- 同位點 TO vs Paired HP_Ratio Cohen's d = -1.20（巨大效應）
- 67.6% phase blocks 被 somatic variant 汙染

### 3.3 Phase Block 特徵比較（HCC1395 實測）

| 指標 | TO | Paired | 比值 |
|------|-----|--------|------|
| Block 數 | 1,594 | 4,615 | 0.35x |
| Median block | 293 Kbp | 222 Kbp | 1.3x |
| N50 | 11.9 Mbp | 1.2 Mbp | **10.2x** |
| Max block | 77.7 Mbp | 7.9 Mbp | 9.8x |
| HP assign rate | 0.924 | 0.853 | 1.08x |
| LOH.bed 覆蓋率 | 13-61% | N/A | 遠超預期 |

**反直覺現象**：TO 的 phase block 更大，是因為 phasing variants 較少導致 switch error 機會減少，LongPhase 將稀疏 variants 串成 mega-blocks。更大的 blocks 強制更多 reads 進入 HP assignment，但 HP 品質差（因 somatic 汙染）。

---

## 4. 相關工具與改善方案

### 4.1 WhatsHap — Read-Based Phasing 基準工具

**論文**: Martin et al. (bioRxiv, 2016); Patterson et al. (J Comput Biol, 2015)

**核心演算法**:
- 將 phasing 建模為 Minimum Error Correction (MEC) 問題
- 使用動態規劃 (DP) 求解，精確求解 NP-hard 問題
- 透過 downsampling reads（預設 15x）加速

**Phasing 品質評估工具**:
- `whatshap compare`: 計算 switch error rate、Hamming distance、phase block N50
- `whatshap stats`: 報告 phase block 統計數據
- 這是目前最廣泛使用的 phasing benchmark 工具

**與本研究的關聯**:
- WhatsHap 本身是 germline phasing 工具，不處理 somatic variants
- 其 `compare` 模組可用於量化 LongPhase-TO 的 switch error rate
- **局限**: 在 tumor-only 模式下沒有 truth phasing 可供比對

**Switch error rate 基準值**:
- WhatsHap germline: ~0.02% at full coverage
- LongPhase germline: ~0.055%（ONT）
- **TO somatic phasing: 無公開基準值**（目前文獻缺口）

### 4.2 HapCUT2 — 多平台 Haplotype Assembly

**論文**: Edge et al. (Genome Research, 2017)

**核心特點**:
- 支援多種數據類型（WGS, linked-reads, Hi-C, SMRT）
- Max-likelihood 最佳化
- 迭代切割和合併 haplotype fragments

**與本研究的關聯**:
- HapCUT2 同樣是 germline-focused，不特別處理 tumor-only 情境
- 可作為替代 phasing 引擎的 benchmark 對照
- 不提供 somatic-aware phasing 功能

### 4.3 MethPhaser — 甲基化輔助 Phasing

**論文**: Nature Communications, 2024 (DOI: 10.1038/s41467-024-49588-0)

**核心方法**:
1. 在已有 SNV phase blocks 內識別 haplotype-specific methylation patterns（Wilcoxon rank sum test, p<0.05）
2. 利用甲基化信號將未 phased reads 分配到 haplotype
3. 跨越 homozygous regions 橋接斷開的 phase blocks

**N50 改善**: ONT R9 60X 數據上提升 1.6-2.5 倍，accuracy 83.4-98.7%

**與本研究的關聯**:
- **潛在改善方向 A**: 在 TO 模式下，可用甲基化信號作為額外的 phasing 證據，減少對 somatic variants 的依賴
- **限制**: 仍需 pre-phased SNVs 作為 seed（不能完全獨立於 SNP phasing）
- **Cancer 應用**: 作者承認尚未針對 cancer genomes 最佳化，amplified regions 的不同甲基化信號可能造成問題
- **實用性評估**: 中等。可以改善 phase block 的連續性，但不直接解決 self-phasing 問題

### 4.4 Octopus — Haplotype-Aware Variant Calling with Phase Quality

**論文**: Cooke et al. (Nature Biotechnology, 2021; DOI: 10.1038/s41587-021-00861-3)

**核心特點**:
- Bayesian haplotype-aware framework
- 支援 germline、somatic paired、somatic tumor-only 三種模式
- **同時呼叫 variant 和 phase**，而非先 call 再 phase
- 提供 phase quality score（Phred-scaled）

**Phase Quality Score**:
- QUAL: variant 存在的後驗機率
- PP (INFO field): variant 被正確分類的後驗機率
- **Somatic-germline phasing accuracy**: Phase quality >= 10 時 98% 正確，>= 20 時 99% 正確

**Tumor-Only 模式**:
- Sensitivity 0.45 vs Paired 0.84（顯著下降）
- 大部分 FN 實際上被 call 出來但被錯誤分類為 germline
- Quality score calibration 仍然良好

**與本研究的關聯**:
- **潛在改善方向 B**: Octopus 的 joint calling + phasing 架構避免了 self-phasing 問題，因為 variant calling 和 phasing 同時進行而非循序
- **限制**: Octopus 是 short-read 為主的工具，long-read 支援有限
- **Phase quality score 概念可借鑑**: 為每個 phased variant 提供可信度分數，下游可據此過濾

### 4.5 SAVANA — Haplotype-Resolved SV/CNA 分析

**論文**: Nature Methods, 2025

**核心方法**:
- Haplotype-aware SV 偵測 + CNA profiling + purity/ploidy estimation
- 支援 tumor-only 模式
- 對 phased reads 的 SV 進行 haplotype assignment

**與本研究的關聯**:
- 展示了 haplotype-resolved cancer analysis 的最新進展
- 在 99 個 tumor-normal pairs 上驗證，sensitivity 和 specificity 均優於競爭工具
- **啟示**: haplotype-aware 分析對 cancer genomics 至關重要，phasing 品質直接影響下游分析

### 4.6 Severus — Breakpoint Graph SV Calling

**論文**: Nature Biotechnology, 2025

**核心方法**:
- 從 split alignments 呼叫 haplotype-aware junctions
- **Phasing 策略**: 使用 normal sample 的 SNV calling + phasing → haplotag tumor and normal
- 支援 unbalanced cancer karyotypes

**與本研究的關聯**:
- Severus 明確使用 **normal sample 的 germline phasing** 作為基礎，避免了 tumor variants 汙染 phasing scaffold
- 這是 tumor-normal 模式的最佳實踐，但不適用於 tumor-only

### 4.7 Wakhan — Phase Block 延伸與 Switch Error 修正

**論文**: Kolmogorov Lab (與 Severus 同組)

**核心方法**:
- 利用 **copy number 差異** 延伸 phase blocks 和修正 switch errors
- 在 haplotype 有不同 copy number 的區域，read depth 差異提供額外的 phasing 信號
- 包含 `hapcorrect()` 功能，專門修正 phasing errors

**支援模式**:
- Tumor-Normal: tumor BAM + normal phased VCF
- **Tumor-Only**: tumor BAM + tumor phased VCF
- Unphased: 低雜合度物種

**與本研究的關聯**:
- **潛在改善方向 C**: Wakhan 的 copy-number-based phase correction 可作為 post-hoc 修正手段
- 在 TO 模式下，Wakhan 可利用 CNA 信號修正 LongPhase-TO 的 switch errors
- **限制**: 依賴 CNA 信號存在，在 copy-neutral 區域可能無效
- **實用性評估**: 高。直接適用於 tumor-only，且修正已有的 phasing 結果

### 4.8 PhaseME — Phasing 品質評估與改善

**論文**: GigaScience, 2020 (DOI: 10.1093/gigascience/giaa078)

**核心功能**:
- Phase quality ratio（每個 block 內 SNV 受 population linkage 支持的比例）
- Switch error 偵測與修正
- 區分 short switch (2-20 bp) 和 long switch (>=21 bp)

**效能**: 平均降低 Hamming error rate 22.4%，PacBio HiFi long switch 改善 54.6%

**與本研究的關聯**:
- **潛在改善方向 D**: 可用於評估 LongPhase-TO 的 phasing 品質，並修正明顯的 switch errors
- **限制**: 針對 germline phasing 開發，未驗證 cancer/somatic 情境
- 不直接解決 self-phasing 問題，但可偵測因 self-phasing 造成的異常 switch patterns

### 4.9 Lung Cancer Phasing 分析

**論文**: Sato et al. (Nature Communications, 2022; DOI: 10.1038/s41467-022-31133-6)

**研究設計**:
- 20 名日本非小細胞肺癌患者
- Short + long read WGS 結合
- WhatsHap (v1.0) 進行 phasing
- N50 = 834 kbp，concordance rate >99%

**關鍵發現**:
- Cancer genomes 中 mutations 在兩個 haplotypes 間不均勻分佈
- Chromothripsis 事件僅發生在單一 chromosome，造成 haplotype-biased mutation 分佈
- EGFR 突變陽性肺腺癌的特徵性事件

**與本研究的關聯**:
- 展示了 cancer genome phasing 的生物學意義
- 使用 germline phasing（WhatsHap）作為 scaffold，somatic mutations 被動分配到 haplotype
- **未使用 somatic variants 作為 phasing anchor** — 與 LongPhase-TO 的做法形成對比

### 4.10 Long-Read Cancer Cohort（POG Study）

**論文**: O'Neill et al. (Cell Genomics, 2024)

**研究設計**:
- 189 個 tumor + 41 個 matched normal，ONT PromethION
- Long-range phasing 發現 allelically differentially methylated regions (aDMRs)
- 包括 RET、CDKN2A、BRCA1、RAD51C 等 cancer genes 的 aDMRs

**與本研究的關聯**:
- 展示 haplotype-resolved methylation 在 cancer genomics 的臨床價值
- 使用 matched normal 進行 phasing，保證 phasing 品質
- aDMR 分析依賴高品質 phasing，self-phasing 可能導致虛假 aDMR 結果

### 4.11 Somatic Mutation Phasing in Multiple Myeloma

**論文**: PMC11326269 (2024)

**核心方法**:
- 使用 linked-reads 的 barcode 信息 phase somatic mutations
- Linked alleles method: 要求 >=91% linked alleles 來自同一 haplotype（precision 0.997, recall 0.936）
- Cross-sample haplotype extension: 利用同一患者不同樣本的共享 germline 信息擴展 phase blocks（4.6 倍延伸）

**與本研究的關聯**:
- **Key insight**: Somatic mutations 是 **被動分配** 到 germline haplotypes 的，而非作為 anchor
- 要求高 concordance threshold (91%) 確保 phasing 正確性
- **啟示**: 即使在 linked-reads 中，somatic mutation phasing 也是 secondary 步驟

### 4.12 LongPhase-S — Paired Mode 的正確做法

**論文**: Ho et al. (bioRxiv, 2025; DOI: 10.1101/2025.11.20.689492)

**核心方法**:
1. Normal sample 的 germline SNPs 經 LongPhase 進行 phasing → 純 germline scaffold
2. Somatic haplotagging: tumor reads 根據 germline scaffold 分配 HP（HP1-1, HP2-1, HP3）
3. Purity estimation: 利用 germline haplotype imbalance ratio (GHIR)
4. Variant recalibration: 結合 purity 和 haplotype consistency 進行 somatic variant 過濾

**與 TO 的關鍵差異**:
- **Somatic variants 完全不參與 phasing scaffold**
- Germline phasing 品質有保證（來自 normal sample）
- HP assignment 不受 somatic variant 影響
- Variant recalibration 可利用 haplotype consistency 作為過濾依據

---

## 5. 衝突觀點分析

| 觀點 A | 觀點 B | 可能原因 |
|--------|--------|----------|
| TO mega-blocks（N50=11.9 Mbp）看起來更好 | Paired blocks（N50=1.2 Mbp）實際品質更高 | TO 的大 blocks 是因為稀疏 variants 缺少 break point，不代表 phasing 正確 |
| MethPhaser 可改善 TO phasing | MethPhaser 需要 pre-phased SNVs 作為 seed | 如果 seed phasing 已有 self-phasing 汙染，甲基化擴展可能放大錯誤 |
| Octopus joint calling+phasing 可避免 self-phasing | Octopus 主要支援 short-reads | Long-read 和 short-read phasing 的演算法基礎不同 |
| TO phasing 因為 somatic anchor 造成虛假 LOH | 部分 LOH 可能是真實的（AF >0.9 區域） | Self-phasing 主要影響 AF 0.1-0.8 範圍；AF >0.9 的 LOH 多為真實 structural LOH |

---

## 6. 資料來源評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| LongPhase 原始碼 (PhasingProcess.cpp, PhasingGraph.cpp) | 一手原始碼 | **最高** | 直接確認 addEdge() 在 somaticCalling() 之前執行 |
| LongPhase-TO vs -S 因果鏈報告 | 一手實驗數據 | **最高** | 748K regions, 288K same-locus pairs, 7 samples |
| Lin et al. 2022 (Bioinformatics) | 已發表論文 | 高 | LongPhase 原始演算法 |
| Ho et al. 2025 (bioRxiv) | 預印本 | 高 | LongPhase-S，審稿中但數據完整 |
| Cooke et al. 2021 (Nature Biotech) | 已發表論文 | 高 | Octopus phase quality score |
| Cell Genomics 2024 (POG Study) | 已發表論文 | 高 | 189 tumor long-read cancer cohort |
| Sato et al. 2022 (Nature Comms) | 已發表論文 | 高 | Lung cancer phasing |
| SAVANA 2025 (Nature Methods) | 已發表論文 | 高 | Haplotype-resolved SV |
| Severus 2025 (Nature Biotech) | 已發表論文 | 高 | Breakpoint graph SV |
| MethPhaser 2024 (Nature Comms) | 已發表論文 | 高 | Methylation-based phasing |
| PhaseME 2020 (GigaScience) | 已發表論文 | 中 | Germline-only 驗證 |
| Wakhan (GitHub) | 工具文件 | 中 | 尚未正式發表 |
| LongPhase-TO GitHub README | 工具文件 | 中 | 演算法細節不足 |

---

## 7. 可行改善方案評估

### 方案 A: Pre-filtering Somatic Variants from Phasing Input（最直接）

**概念**: 在 LongPhase-TO phasing 之前，從 VCF 移除已知/疑似 somatic variants，僅保留高信度 germline variants 作為 phasing scaffold。

**實施方式**:
1. 使用 PON 嚴格過濾：僅保留 NonSomatic 標記的 variants
2. 額外使用 AF 過濾：排除 AF < 0.3 或 AF > 0.7 的 variants（可能是 somatic 或 LOH）
3. 排除 PASS variants（ClairS-TO 的 somatic candidates）
4. 用清理後的 VCF 重跑 LongPhase-TO phase

**優點**: 直接從根因解決 self-phasing；不需修改 LongPhase-TO 原始碼
**缺點**: germline variant 密度可能不足（PON 過濾後剩餘的 variants 品質不均）；可能產生更短的 phase blocks
**可行性**: **高** — 可立即用現有工具實施
**預期效果**: 消除 AF 0.1-0.8 範圍的虛假 LOH；phase blocks 可能變短但品質提升

### 方案 B: Two-Pass Phasing（先 Phase 再 Tag）

**概念**: 第一輪用全部 variants 建構 phasing scaffold；第二輪排除 somatic variants 重新 phase。

**實施方式**:
1. Round 1: LongPhase-TO phase（現有流程）→ 標記 somatic variants
2. Round 2: 移除 Round 1 標記的 somatic variants，僅用 germline variants 重新 phase
3. Round 3: 用 Round 2 的 phased VCF 進行 haplotag

**優點**: 利用 Round 1 的 somatic calling 結果改善 Round 2 的 scaffold
**缺點**: 需跑兩輪 phasing（時間加倍）；Round 1 的 somatic calling 品質可能受 self-phasing 影響
**可行性**: **中** — 需要修改 pipeline 但不需改原始碼
**備註**: LongPhase-TO 原始碼顯示已有類似邏輯（高純度 >0.95 時 `convertNonGermlineToSomatic()` 後重跑 `phasingProcess()`），但方向相反——它是將更多 variants 標為 somatic 後重跑，而非排除 somatic

### 方案 C: Post-Hoc Phase Correction with Wakhan

**概念**: 使用 Wakhan 的 copy-number-based phase correction 修正 LongPhase-TO 的 phasing errors。

**實施方式**:
1. 正常執行 LongPhase-TO phase + haplotag
2. 使用 Wakhan `hapcorrect()` 修正 phasing errors
3. 輸出 rephased VCF 供下游分析

**優點**: 不需修改 LongPhase-TO；利用獨立信號（CNA）修正 phasing
**缺點**: 僅在有 CNA 的區域有效；copy-neutral 區域的 self-phasing 無法修正
**可行性**: **中** — Wakhan 支援 tumor-only 模式，但效果受 CNA 分佈限制

### 方案 D: MethPhaser-Enhanced Phase Block Extension

**概念**: 用 MethPhaser 擴展以 germline-only variants 建構的（較短）phase blocks。

**實施方式**:
1. 方案 A 的 germline-only phasing → 較短但品質好的 phase blocks
2. MethPhaser 利用甲基化信號延伸和橋接 phase blocks

**優點**: 結合兩種獨立信號（germline SNP + methylation）
**缺點**: MethPhaser 尚未針對 cancer genomes 最佳化；需額外處理步驟
**可行性**: **低-中** — 概念可行但需要驗證甲基化信號在 tumor 中的可靠性

### 方案 E: Somatic Variant 的 Phase Confidence Scoring

**概念**: 借鑑 Octopus 的 phase quality score 概念，為 LongPhase-TO 的每個 phased somatic variant 添加可信度分數。

**實施方式**:
1. 計算每個 somatic variant 的 self-phasing score: 排除該 variant 自身 reads 後的 HP ratio
2. 若排除後 HP ratio 大幅改變（例如從 0.95 降到 0.55）→ 標記為 self-phased
3. 提供 corrected HP ratio 供下游使用

**優點**: 不需重跑 phasing；可識別 self-phasing 的嚴重程度
**缺點**: 需要開發新的分析模組；計算量可能較大
**可行性**: **中-高** — 概念清晰，可作為 ISM 的下游改善
**備註**: 這與因果鏈驗證報告中的分析方法一致

### 方案 F: LongPhase-TO 原始碼修改（上游修正）

**概念**: 修改 LongPhase-TO 原始碼，在 `addEdge()` 階段排除 somatic variants。

**實施方式**:
1. 在 `PhasingProcess.cpp` 中，將 `somaticCalling()` / `tagSomatic()` 移到 `addEdge()` 之前
2. 在 `addEdge()` 中跳過已標記為 somatic 的 variants
3. Somatic variants 僅在 haplotag 階段被動分配 HP

**優點**: 從根本解決 self-phasing
**缺點**: 需要修改上游工具原始碼；可能需要重新驗證整個 pipeline
**可行性**: **低**（短期）/ **高**（長期）— 需要與 LongPhase 開發團隊合作

---

## 8. 建議行動（按優先順序）

### 立即可行（1-2 天）

1. **方案 A 驗證**: 在 1 個樣本上測試 germline-only phasing scaffold：
   - 從 ClairS-TO VCF 僅提取 NonSomatic variants（排除 PASS + LowQual + MultiHap）
   - 用清理後 VCF 跑 LongPhase-TO phase + haplotag
   - 比較 HP ratio 分佈、LOH.bed 覆蓋率、ISM 結果

2. **方案 E 原型**: 在現有 ISM 分析中計算 self-phasing score：
   - 對每個 somatic variant，計算「排除自身 ALT reads 後的 HP ratio」
   - 與原始 HP ratio 比較，量化 self-phasing 效應

### 短期（1-2 週）

3. **方案 B 實施**: 若方案 A 效果好，設計 two-pass pipeline
4. **WhatsHap compare 分析**: 用 WhatsHap 量化 LongPhase-TO vs LongPhase-S 的 phasing 一致性
5. **Wakhan 評估**: 在 1 個樣本上測試 Wakhan post-hoc correction 的效果

### 中期（1-2 月）

6. **方案 F 討論**: 與 LongPhase 開發團隊討論原始碼修改的可行性
7. **MethPhaser 整合**: 評估甲基化輔助 phasing 在 cancer genomes 的效果
8. **Benchmarking**: 在所有 7 個樣本上系統性比較各改善方案

---

## 9. 關鍵文獻參考

1. Lin et al. "LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants." *Bioinformatics*, 2022. DOI: [10.1093/bioinformatics/btac058](https://doi.org/10.1093/bioinformatics/btac058)
2. Ho et al. "LongPhase-S: purity estimation and variant recalibration with somatic haplotyping for long-read sequencing." *bioRxiv*, 2025. DOI: [10.1101/2025.11.20.689492](https://doi.org/10.1101/2025.11.20.689492)
3. Cooke et al. "A unified haplotype-based method for accurate and comprehensive variant calling." *Nature Biotechnology*, 2021. DOI: [10.1038/s41587-021-00861-3](https://doi.org/10.1038/s41587-021-00861-3)
4. MethPhaser. "Methylation-based long-read haplotype phasing of human genomes." *Nature Communications*, 2024. DOI: [10.1038/s41467-024-49588-0](https://doi.org/10.1038/s41467-024-49588-0)
5. PhaseME. "Automatic rapid assessment of phasing quality and phasing improvement." *GigaScience*, 2020. DOI: [10.1093/gigascience/giaa078](https://doi.org/10.1093/gigascience/giaa078)
6. Sato et al. "Phasing analysis of lung cancer genomes using a long read sequencer." *Nature Communications*, 2022. DOI: [10.1038/s41467-022-31133-6](https://doi.org/10.1038/s41467-022-31133-6)
7. O'Neill et al. "Long-read sequencing of an advanced cancer cohort resolves rearrangements, unravels haplotypes, and reveals methylation landscapes." *Cell Genomics*, 2024. DOI: [10.1016/j.xgen.2024.100653](https://doi.org/10.1016/j.xgen.2024.100653)
8. SAVANA. "Reliable analysis of somatic structural variants and copy number aberrations using long-read sequencing." *Nature Methods*, 2025. DOI: [10.1038/s41592-025-02708-0](https://doi.org/10.1038/s41592-025-02708-0)
9. Severus. "Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing." *Nature Biotechnology*, 2025. DOI: [10.1038/s41587-025-02618-8](https://doi.org/10.1038/s41587-025-02618-8)
10. Chen et al. "ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling." *Nature Communications*, 2025. DOI: [10.1038/s41467-025-64547-z](https://doi.org/10.1038/s41467-025-64547-z)
11. Martin et al. "WhatsHap: fast and accurate read-based phasing." *bioRxiv*, 2016. DOI: [10.1101/085050](https://doi.org/10.1101/085050)
12. Rare variant phasing using paired tumor:normal sequence data. *BMC Bioinformatics*, 2019. DOI: [10.1186/s12859-019-2753-1](https://doi.org/10.1186/s12859-019-2753-1)
13. Somatic mutation phasing and haplotype extension using linked-reads in multiple myeloma. 2024. PMC11326269.
