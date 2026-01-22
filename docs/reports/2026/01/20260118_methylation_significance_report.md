<!--
建立時間: 2026-01-18
範圍: 甲基關聯顯著性方法與限制 + 文獻回顧 + HCC1395 輸出分析 + 後續建議
資料來源:
  - output/bip8_disk_output/20260113_all-with-w5000_3/analysis/report.md
  - output/bip8_disk_output/20260113_all-with-w5000_3/filtered_snv_*/*.txt
  - data/README.md
  - src/core/*, include/core/*
  - PubMed/DOI (見參考資料)
-->
# 甲基關聯顯著性分析：方法、限制與後續方向 (HCC1395 ONT)

報告日期: 2026-01-18  
分析輸出: `output/bip8_disk_output/20260113_all-with-w5000_3`  
樣本: HCC1395 ONT tumor/normal (tumor 75x / normal 45x, purity 1.0)

---

## 1. 研究目標與主軸

**核心假設**: 若 sSNV 為真，則其周邊 read-level 甲基化模式可能與 ALT/REF 或 HP (haplotype) 產生可偵測的差異，進而作為 TP/FP 的輔助驗證訊號。

**研究目的**:
1. 以甲基化訊號輔助辨識 TP/FP (不作為唯一過濾條件)。
2. 在 read 層級探索亞克隆 (subclone) 的甲基化結構。
3. 整合 genotype (SNV)、haplotype (HP)、epigenotype (methylation) 的多視角分析。

---

## 2. 現有實作方式與方法細節 (InterSubMod)

### 2.1 輸入資料與前處理
- **BAM**: ONT tumor/normal，含 MM/ML (5mC) 與 HP/PS 標籤。
- **VCF**: ClairS pileup + longphase-s 強篩選後的 sSNV (TP/FP 以 SEQC2 驗證集標記)。
- **視窗**: 本次輸出為 ±5000 bp (`--window-size 5000`)。
- **過濾條件**: `MAPQ >= 20`、`read length >= 1000`，需有 MM/ML。

**甲基化解析**:
- 解析 `MM/ML` 的 `C+m?` 5mC，並處理反股 read 的序列方向。
- 只保留 CpG context 的 methylation calls。

### 2.2 Read×CpG 矩陣與距離矩陣
- 產生 `Read × CpG` 的甲基化矩陣 (`-1` 表示缺失)。
- 距離矩陣以 read 對為單位計算，只在共同 CpG 位點上計算。
- **本次使用距離**: Bernoulli (預設二值化，`C_min=3`)。

### 2.3 聚類與顯著性分析
**Cluster-First (Phase 1-3)**
1. **GlobalTest**: Fisher-Freeman-Halton + Cramer's V。
2. **Gating**: `p <= 0.1` 且 `V >= 0.1` 才進入後續測試。
3. **LocalTest**: one-vs-rest。
4. **StructureTest**: PERMANOVA + dispersion。

**Label-First (Phase 4)**
- 以 HP/Allele 作為分組，計算 Delta = d_between - d_within，並做 permutation test。
- HP 也提供多群 PERMANOVA (HP1/HP2/HP0/HP3...)。

**Phase 5: 雙向驗證分類**
- **Strong**: cluster 與 label 同時顯著。
- **Subclone**: cluster 顯著但 label 不顯著。
- **Weak**: label 顯著但 cluster 不顯著。
- **Noise**: 兩者皆不顯著。

### 2.4 輸出指標
- `significance_summary.csv`: 每個 SNV 區域的顯著性指標 (P, V, Delta, 類別)。
- `significance_statistics.txt`: 匯總、染色體統計。
- `analysis/report.md`: TP/FP 比較與 ROC。

---

## 3. 讀段層級相關研究回顧 (SNV × 甲基化)

### 3.1 Read-level / haplotype-aware / allele-specific 方法

**重點**: 關注「read-level / haplotype-aware / allele-specific」的甲基化分析方法，並整理輸入資料、方法、輸出與限制。

| 研究 | 輸入資料與特性 | 方法重點 | 輸出/結論 | 主要限制 |
|---|---|---|---|---|
| **NanoMethPhase** (Genome Biol 2021, DOI: 10.1186/s13059-021-02283-5) | ONT long-read + 5mC calls + phased SNPs | 將 reads 依 heterozygous SNP 分相，計算 haplotype-specific methylation | 提供 megabase-scale haplotype methylation | 需足夠 SNP 密度與深度，phasing 誤差會放大 |
| **MethPhaser** (Nat Commun 2024, DOI: 10.1038/s41467-024-49588-0) | ONT long-read methylation | 利用甲基化訊號延伸 SNV-based phasing | 擴展 haplotype block、支援 ASM 分析 | 依賴 methylation 變異度，深度不足時效果受限 |
| **NANOME** (bioRxiv 2025, DOI: 10.1101/2025.06.29.662079) | ONT reads + methylation callers | pipeline 化 allele-specific methylation | 產出 haplotype-aware consensus methylation | pipeline 複雜、對輸入品質與 caller 選擇敏感 |
| **Long-range phasing of regulatory elements** (Nat Genet 2022, DOI: 10.1038/s41588-022-01188-8) | targeted ONT reads | 在單分子層級同時偵測 chromatin accessibility + methylation | 證實 allele-specific 調控元素可被長讀段分相 | 主要為 targeted assay，外推至 WGS 需注意 |
| **Fine mapping regulatory variants** (HGG Adv 2025, DOI: 10.1016/j.xhgg.2025.100532) | ONT long reads (native CpG) | 連結變異與甲基化模式、定位調控變異 | 建立變異-甲基化關聯的 read-level 框架 | 需要高 coverage、結果受 mapping bias 影響 |
| **Cancer long-read cohort** (Cell Genomics 2024, DOI: 10.1016/j.xgen.2024.100674) | 189 tumor + 41 normal ONT | WGS long-read + methylation landscape | 證實腫瘤中 aDMR/ASM 存在、與癌症基因相關 | 腫瘤異質性與純度影響顯著 |
| **Medulloblastoma long-read** (Cell Genomics 2023, DOI: 10.1016/j.xgen.2023.100281) | tumor paired ONT | haplotype-resolved epigenetic landscape | 展示 read-level epigenetic signatures | 單案例研究，泛化需多樣本 |

### 3.2 特定腫瘤類型的 long-read 甲基化研究

| 腫瘤/樣本 | 研究 (DOI) | 輸入資料與方法 | 主要輸出 | 限制/提醒 |
|---|---|---|---|---|
| **癌症 cfDNA** | Genome Med 2023 (10.1186/s13073-023-01178-3) | nanopore 單分子 cfDNA methylome | 以 read-level methylation profile 做癌症相關分類/檢測 | cfDNA 片段化嚴重，需高 read 數 |
| **白血病** | Commun Biol 2023 (10.1038/s42003-023-04756-8) | nanopore WGS methylome maps | 發現 CpG-poor 區域異常甲基化與耐藥相關 | cohort/疾病特異，難直接外推 |
| **急性白血病 (臨床分類)** | Nat Genet 2025 (10.1038/s41588-025-02321-z) | nanopore methylation profiling + 分類模型 | 2 小時內達到臨床分類等級的準確度 | 以分類為主，非 SNV 驗證 |
| **兒童腦瘤 (術中樣本)** | J Neurooncol 2024 (10.1007/s11060-024-04702-6) | 低深度 nanopore WGS + methylation classification | 快速分類兒童腦瘤亞型 | 低深度下 read-level 穩定性有限 |
| **膀胱癌 (尿液 DNA)** | Clin Epigenetics 2025 (10.1186/s13148-025-01946-5) | 尿液 DNA long-read WGS + 直接甲基化偵測 | genome-wide methylation 變化可用於偵測 | 尿液 DNA 品質/含量變異大 |
| **卵巢癌 (HGSOC)** | Sci Rep 2025 (10.1038/s41598-025-21907-5) | ONT LRS 探索重複區的 genomic/epigenomic 變化 | 揭露重複序列區域的甲基化異常 | 受重複區 mapping 與樣本數限制 |
| **腎臟腫瘤 (配對樣本)** | NAR Genom Bioinform 2025 (10.1093/nargab/lqae190) | nanopore + optical mapping | 同時刻畫結構變異與表觀變化 | 個案型研究，統計力有限 |

### 3.3 Read-level 方法堆疊與評估基礎

- **甲基化 calling**: 以 nanopore raw signal 或 basecaller/modbase 產生 `MM/ML` 標籤；效能依 caller 與化學版本而異。
- **讀段分相**: 常見做法是先用 SNP/phase set (WhatsHap/longphase-s) 分相，再建立 haplotype-specific methylation。
- **工具評估**: 系統性 benchmark 顯示 caller 之間差異顯著，且 coverage、CpG 密度與 read quality 會明顯影響偵測率與 false positives。
- **大型方法資源**: 近年已可在全基因組層級建立 haplotype-resolved methylation（支援後續 ASM 或 read-level clustering）。

**方法型參考**:
- Systematic benchmarking of tools for CpG methylation detection from nanopore sequencing (Nat Commun 2021, DOI: 10.1038/s41467-021-23778-6)
- DNA methylation-calling tools for Oxford Nanopore sequencing: a survey and human epigenome-wide evaluation (Genome Biol 2021, DOI: 10.1186/s13059-021-02510-z)
- Scalable Nanopore sequencing of human genomes provides a comprehensive view of haplotype-resolved variation and methylation (Nat Methods 2023, DOI: 10.1038/s41592-023-01993-x)

**正向觀點**:
- 長讀段可在單分子層級同時觀察 SNV/phase/methylation，技術可行。
- 研究已證明 haplotype-specific methylation (ASM) 廣泛存在。

**反向/限制觀點**:
- 需足夠深度與高品質 phasing；read-level 噪音會稀釋效應。
- 不是所有 SNV 都會造成 methylation 變化，效應量偏小。
- 甲基化訊號常受腫瘤純度、CNV、LOH 影響。

---

## 4. HCC1395 輸出結果分析 (20260113_all-with-w5000_3)

### 4.1 基本統計 (TP vs FP)

| 指標 | TP | FP |
|---|---:|---:|
| Regions (成功分析) | 30,476 | 4,822 |
| 顯著數 | 1,860 | 101 |
| 顯著率 | 6.10% | 2.09% |
| Reads 平均 | 71.4 | 59.8 |
| CpG 平均 | 97.4 | 98.0 |
| P-value 平均 | 0.527 | 0.429 |
| Cramer's V 平均 | 0.052 | 0.017 |
| Label Delta 平均 | 0.039 | 0.077 |

**觀察**:
- TP 顯著率高於 FP，但幅度有限 (6.1% vs 2.1%)。
- FP 的平均 Label Delta 反而較高，顯示 label delta 更可能反映噪音或偏差。

### 4.2 分類分布
- **VerificationClass**:
  - TP: Strong 24.5%, Weak 43.5%, Noise 27.9%, Subclone 4.0%
  - FP: Strong 37.1%, Weak 31.4%, Noise 27.7%, Subclone 3.8%
- **DominantLabel**:
  - TP: hp 68.6%, allele 16.9%, none 14.5%
  - FP: hp 69.2%, allele 26.2%, none 4.6%

**異常**: FP 的 Strong 比例高於 TP，與原始假設相反。

### 4.3 ROC 鑑別力
| 指標 | AUC |
|---|---:|
| Cramer's V | 0.519 |
| Heuristic Score | 0.444 |
| -log10(P) | 0.439 |
| Label Delta | 0.409 |

**結論**: 指標鑑別力接近隨機 (AUC ≈ 0.5)。

### 4.4 染色體層級異常
- TP 顯著率最高: chr16 (17.3%), chr19 (10.3%), chr20 (9.7%)。
- FP 顯著率異常: chr9 (17.5%), chr14 (10.8%)。

**解讀**: 某些染色體可能存在結構性偏差或特定調控區域，需要進一步排查。

---

## 5. 方法有效性與問題診斷

### 5.1 是否足夠完整有效?
- **完整性**: 有完整輸出 (全區域顯著性、分類、ROC)，但缺少多重比較校正與生物註解 (CpG density/功能區域)。
- **有效性**: TP/FP 鑑別力弱 (AUC ≈ 0.5)，直接作為過濾會嚴重惡化 F1。

### 5.2 與原始假設的出入
- 原假設: TP 會有顯著甲基化關聯。
- 觀察: 大多數 TP 無顯著訊號；FP 的 Strong 比例更高。

**可能原因**:
1. **生物學效應量小**: 多數 SNV 不位於甲基化調控區域。
2. **資料稀疏與缺失**: CpG 覆蓋不足造成距離矩陣不穩定。
3. **phasing/allele 標籤不穩**: HP0/HP3 或 ALT/REF 判定誤差。
4. **統計門檻設定**: gating (p<=0.1, V>=0.1) 可能過嚴或不符合實際分布。
5. **只用 tumor 1.0**: 缺少 tumor-normal 差異對照。

### 5.3 特殊發現
- **Label Delta 在 FP 較高**: 代表標籤與距離結構可能在 FP 受到技術偏差影響。
- **chr9/chr14 FP 顯著率異常**: 可能是重複區或 mapping 偏差，需要位點級排查。

---

## 6. 後續分析與實作建議

### 6.1 其他可能有用的資料 (5.1)
- **SNV 層級**: VAF, QUAL, DP, AD, strand bias (FAU/RAU)。
- **結構變異**: CNV/LOH、copy-number-adjusted VAF。
- **樣本資訊**: tumor purity (0.2-1.0)、phase block 長度、HP0/HP3 比例。
- **區域註解**: CpG density、promoter/enhancer、PMD/LMR。

### 6.2 需要再驗證的部分 (5.2)
- **Strong 分類異常**: FP > TP 的原因 (邏輯或標籤問題)。
- **phasing 品質**: HP1/HP2 的一致性與錯配率。
- **距離矩陣稀疏性**: 有效距離比例 vs 顯著率。
- **chr9/chr14 異常**: 檢查 mapping bias、重複序列。

### 6.3 重點觀察分析的實作方式 (5.3)
- **多特徵整合模型**:
  - 以 `VCF + methylation` 特徵建立 logistic/RF/XGBoost。
  - 使用 LOCO (leave-one-chromosome-out) 交叉驗證。
- **差異甲基化對照**:
  - 以 tumor vs normal 的 read-level methylation 差異作為特徵 (delta_meth)。
- **門檻搜尋**:
  - Grid search: `p, V, Delta, min_reads, min_cpgs`。
- **亞克隆分析**:
  - 聚合 Subclone 類別，檢查是否對應 VAF/CNV 群。

### 6.4 多樣本與純度梯度 (5.4)
- **HCC1395**: 三組 ONT T/N，純度 1.0/0.8/0.6/0.4/0.2。
- **其他樣本**: COLO829、HCC1937、HCC1954、H2009、H1437。
- **建議流程**:
  - 在每個樣本計算 TP/FP 顯著率、AUC、chromosome-level 偏差。
  - 評估「顯著率 vs purity」曲線是否呈正相關。

### 6.5 影響因素與是否加入篩選標準的判定流程

#### 6.5.1 可能影響甲基顯著關聯的因素
- **生物學層面**: ASM 效應量小、調控區域分佈不均、腫瘤異質性/亞克隆、tumor purity、CNV/LOH 造成 VAF/甲基背景漂移。
- **標籤品質**: HP 標籤錯配/相位區塊過短、ALT/REF 標記錯誤或偏差。
- **覆蓋與稀疏**: reads 數不足、CpG 密度低、共同覆蓋不足導致距離矩陣不完整。
- **測序/比對偏差**: flow cell/chemistry 差異、basecaller/modbase 模型差異、mapping bias、strand bias。
- **方法設定**: distance metric、二值化閾值 (0.2/0.8)、min_common_coverage、window size、聚類方法/參數、gating 閾值、無效距離策略。
- **統計層面**: 多重比較未校正、p-value 分布偏態、effect size 被稀釋。

#### 6.5.2 需要完成的觀察分析，才能決定是否加入特定篩選標準
1) **數據可用性檢核**: 統計每個區域的 reads/CpG/有效距離比例分布，確認「顯著」是否只出現在高覆蓋區。
2) **標籤可靠性檢核**: HP0/HP3 比例、phase block 長度、ALT/REF read 支援一致性，排查標籤造成的假訊號。
3) **混雜因子分層**: 以 VAF、DP、QUAL、CNV/LOH 分層，觀察顯著率是否只是品質或結構偏差的副產物。
4) **距離度量敏感性**: 比較 NHD/L1/L2/Bernoulli 等距離對顯著率與 AUC 的影響，確認結果穩定性。
5) **染色體/區域偏差**: 針對 chr9/chr14 等異常區域進行重複序列與 mapping bias 檢查。
6) **多樣本/純度驗證**: 不同樣本/純度下是否仍維持相同趨勢，避免單一樣本過度擬合。
7) **效能評估**: 用 TP/FP 交換率與 F1 變化確認任何新門檻是否有實際提升。

#### 6.5.3 加入篩選標準的判定原則 (建議門檻)
- **效果可量化**: 在 TP/FP 上達到明確提升 (例如 F1 +0.01 以上或 FP/TP 交換率優於基準 1:0.69)。
- **穩定可重現**: 於不同樣本/純度/染色體分層下仍成立。
- **不依賴單一 confounder**: 檢查是否只是被 DP/QUAL/coverage 驅動，而非真正的 methylation signal。
- **門檻可解釋**: 與生物意義一致 (如 CpG density/ASM 區域偏好)。

#### 6.5.4 整合分析與其他路徑建議
- **多特徵整合**: 將 VAF/QUAL/DP/strand bias + methylation 指標 (V, P, Delta, reads, CpGs) 以 logistic/RF 模型整合，避免單一門檻過度簡化。
- **區域層級替代**: 以 aDMR/ASM 區域為單位驗證，而非逐 SNV。
- **tumor-normal 差異化**: 直接使用 tumor vs normal 的 read-level methylation 差異 (delta_meth) 作為特徵。
- **亞克隆/相位導向**: 將 Subclone 類別與 phase block、VAF 分布連結，形成亞克隆甲基化訊號的獨立驗證路徑。

---

## 7. 結論

- 目前甲基化顯著性 **不能作為單獨 TP/FP 過濾器**。
- 主要價值在 **輔助特徵**、**亞克隆探索** 與 **特殊區域深挖**。
- 後續需結合 VAF/QUAL/CNV/LOH 等變異特徵，才能提升整體辨識力。

---

## 參考資料

- NanoMethPhase: https://pubmed.ncbi.nlm.nih.gov/33618748/ (DOI: 10.1186/s13059-021-02283-5)
- MethPhaser: https://pubmed.ncbi.nlm.nih.gov/38909018/ (DOI: 10.1038/s41467-024-49588-0)
- NANOME: https://pubmed.ncbi.nlm.nih.gov/40631091/ (DOI: 10.1101/2025.06.29.662079)
- Long-range phasing (Nat Genet 2022): https://pubmed.ncbi.nlm.nih.gov/36195755/ (DOI: 10.1038/s41588-022-01188-8)
- Fine mapping regulatory variants (HGG Adv 2025): https://pubmed.ncbi.nlm.nih.gov/41109954/ (DOI: 10.1016/j.xhgg.2025.100532)
- Cancer long-read cohort (Cell Genomics 2024): https://pubmed.ncbi.nlm.nih.gov/39406235/ (DOI: 10.1016/j.xgen.2024.100674)
- Medulloblastoma long-read (Cell Genomics 2023): https://pubmed.ncbi.nlm.nih.gov/37082141/ (DOI: 10.1016/j.xgen.2023.100281)
- cfDNA single-molecule methylation (Genome Med 2023): https://pubmed.ncbi.nlm.nih.gov/37138315/ (DOI: 10.1186/s13073-023-01178-3)
- Leukemia nanopore methylome maps (Commun Biol 2023): https://pubmed.ncbi.nlm.nih.gov/37031307/ (DOI: 10.1038/s42003-023-04756-8)
- Acute leukemia nanopore classification (Nat Genet 2025): https://pubmed.ncbi.nlm.nih.gov/40983754/ (DOI: 10.1038/s41588-025-02321-z)
- Pediatric brain tumors, low-pass nanopore (J Neurooncol 2024): https://pubmed.ncbi.nlm.nih.gov/38769169/ (DOI: 10.1007/s11060-024-04702-6)
- Bladder cancer urinary DNA methylation (Clin Epigenetics 2025): https://pubmed.ncbi.nlm.nih.gov/40790606/ (DOI: 10.1186/s13148-025-01946-5)
- Ovarian cancer repeats epigenomics (Sci Rep 2025): https://pubmed.ncbi.nlm.nih.gov/41168225/ (DOI: 10.1038/s41598-025-21907-5)
- Kidney tumor matched sample (NAR Genom Bioinform 2025): https://pubmed.ncbi.nlm.nih.gov/39781516/ (DOI: 10.1093/nargab/lqae190)
- Benchmarking nanopore methylation callers (Nat Commun 2021): https://pubmed.ncbi.nlm.nih.gov/34103501/ (DOI: 10.1038/s41467-021-23778-6)
- Methylation-calling tools survey (Genome Biol 2021): https://pubmed.ncbi.nlm.nih.gov/34663425/ (DOI: 10.1186/s13059-021-02510-z)
- Haplotype-resolved methylation at scale (Nat Methods 2023): https://pubmed.ncbi.nlm.nih.gov/37710018/ (DOI: 10.1038/s41592-023-01993-x)

- 專案內部文件: `data/README.md`, `docs/plans/methylation_significance/20260116_甲基化顯著性分析研究報告.md`
