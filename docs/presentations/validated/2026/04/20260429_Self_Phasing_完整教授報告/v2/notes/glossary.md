# 術語表（v2 PPT 配套）

> 30+ 術語對照表，按出現順序排列；每個術語標明定義 + 在哪張 slide 首次出現。

## A. 工具與系統

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **ClairS-TO** | Tumor-Only somatic variant caller（外部 binary）；輸出 VCF | S2 |
| **LongPhase**（Lin 2022 *Bioinformatics*）| Long-read phasing tool 上游基底；germline 模式 | S2 |
| **LongPhase-S**（bioRxiv 2025）| 同實驗室 paired 模式版本；anchoring 概念來源 | S2 |
| **longphase-to-mod** | LongPhase 本地 fork @ `/big7_disk/liaoyoyo2001/longphase-to-mod/`（**獨立 git repo**）；4-commit 修補在此 | S2 / S3 |
| **InterSubMod (ISM)** | read-level epigenetic characterization（本 repo）；下游消費者，**無 C++ 改動** | S3 |
| **WhatsHap** | 業界主流 germline phasing 工具；不支援 tumor-only | S21 |
| **DeepSomatic**（Nature Biotech 2025）| 平行 caller，含 tumor-only 模式 | S21 |
| **MethPhaser**（Nature Comm 2024）| 用甲基化模式做 long-read phasing | industry_references |

## B. HP tag 與 phasing 概念

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **HP tag**（HP:i） | BAM read-level haplotype tag；整數值 0/1/2/11/21/33 | S4 |
| **HP:i:1, HP:i:2** | germline HP1 / HP2 | S4 |
| **HP:i:11, HP:i:21** | somatic on HP1 / HP2（somatic-linked）| S4 |
| **HP:i:33** | somatic ambiguous（兩 hap 都不確定）| S4 |
| **HP:i:0** | unphased | S4 |
| **HAPLOTYPE1_1 (enum=3)** | longphase C++ 內部 enum；對應 HP:i:11 | S13 |
| **enum vs integer literal mismatch** | 型別語意失配：caller 端 `if(hpResult != HAPLOTYPE1_1)` 用 enum=3，但 hpResult 已是 HP tag integer 11 | S13 |
| **PS** (phase set) | 同一 phase block 的 reads 共享的 ID；VCF/BAM 都可有 | S4 (note) |
| **phased VCF** | 含 GT (genotype with `\|`) + PS 的 VCF | S2 (note) |
| **GT2 / GT3** | LongPhase phased VCF 的 sub-genotype 欄位 | source 03 |
| **phasing graph** | 把 reads 上的 variants 串成同一條 haplotype 的圖結構；edge weight = 共現 | S6 |
| **phase block** | 同一 PS 內的連續 phased region；長度受 read length 與 het 密度限制 | S6 (note) |

## C. PON / Tumor-Only 概念

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **PON** (Panel of Normals) | 群體 germline 變異資料庫（1000g、CoLoRSdb、dbSNP、gnomAD）| S12 |
| **`--pon-only-phasing`** | LongPhase-TO flag；true 時 phasing graph 只用 PON-confirmed germline anchor | S12 |
| **somatic anchor** | 將 somatic variant 當作 phasing graph anchor（baseline 的問題）| S12 |
| **germline anchor** | 將 PON-confirmed germline 當作 phasing graph anchor（V2b 起）| S12 |
| **TO mode (Tumor-Only)** | 無 matched normal sample 的 sequencing 模式 | S2 |
| **Paired mode** | 含 matched normal sample 的對照模式 | S7 |
| **tumor-normal paired** | LongPhase-S 設計模式 | S2 |

## D. Self-phasing 機制術語

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **Self-phasing** | somatic variants 互相 phasing → reads 整體偏向一條 hap，造成 tag bias 17.3:1 | S5 |
| **Self-phasing scaffold**（Phase 層）| phasing graph 用 somatic 當 anchor 形成 self-phasing block | S13 |
| **getVote priority bug**（Tag 層）| `getVote()` 迴圈順序 `{HP1_1,HP2_1}` 在 `{HP1,HP2}` 之前；任一 somatic vote 即 break | S13 |
| **HP1 family** | HP1 + HP1_1 reads（包含 germline HP1 與 somatic on HP1）| S5 |
| **HP2 family** | HP2 + HP2_1 reads | S5 |
| **HP_Ratio** | HP1 family / (HP1 family + HP2 family)；理論 0.5、baseline 0.94 | S6 |
| **HP_Ratio_norm** | HP_Ratio normalized by sample baseline | S9 (note) |
| **AMB%** | Ambiguous reads (HP:i:33) 比例；V5 從 17.5% 降至 8.0% | S17 |
| **somatic bias 17.3:1** | HCC1395 baseline HP1 reads / HP2 reads = 614K/35.5K | S5 |
| **94.6% on HP1** | baseline 全基因組 somatic reads 集中於 HP1 的比例 | S5 |

## E. ISM 特徵術語

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **4-bucket 分群** | ISM 核心：HP1 / HP1-1 / HP2 / HP2-1 四群 | S9 |
| **NGroups (HPFineNGroups)** | 4-bucket 中 occupancy 數 (1/2/3/4) | S22 (note) |
| **HPSig** | HP-derived significance metric | S10 |
| **HPMergedDelta** | HP merged methylation delta | S10 |
| **AlleleDelta** | Allele-based methylation delta | S10 |
| **PairwiseMeanDist** | Pairwise read distance（HP-free, 不受 self-phasing 影響）| S10 |
| **Potential_LOH** | ISM 內基於 HP_Ratio 的 LOH 推測；62% 為 self-phasing artifact | S10 |
| **CramersV** | Cramér's V 效應量；HP/Allele 取最大 | S10 |
| **GlobalP** | Global p-value；HP/Allele 取最小 | S10 |

## F. 結構生物學術語

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **Germline het** | Germline heterozygous variant；應隨機分散於 hap | S6 |
| **Somatic SNV** | Somatic single-nucleotide variant；tumor-acquired | S2 |
| **LOH** (Loss of Heterozygosity) | 雜合變純合區；本工作區分兩層（見下）| S4 |
| **ISM HP_Ratio LOH** | BAM HP tag → HP_Ratio < 0.1 or > 0.9 推測的 LOH；62% 為 artifact | S8 |
| **LOH.bed region-level LOH** | LongPhase region detection 從 VCF AD 產生；不受 self-phasing 影響（Jaccard=1.0）| S8 |
| **cnLOH** (copy-neutral LOH) | 雙親同源無 het 區；V5 仍未解 | S20 |
| **subclone** | tumor sub-population；HPFineNGroups 試圖偵測 | S24 |

## G. 統計與評估指標

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **Concordance** | per-read 比對 V5 / Baseline 與 paired GT 一致率 | S1 |
| **Clean PS blocks** | germline accuracy ≥ 70% 的 PS blocks（高品質 phase block）| S1, S17 |
| **Problem PS blocks** | 與 clean PS 對應的低品質 PS blocks | S20 |
| **F1, Precision, Recall** | 標準二元分類指標 | S20 |
| **AUC** | Area Under Curve；feature discriminative power | S10 (note) |
| **Cohen's d** | 效應量；本工作 −1.20（巨大）| §12 |
| **Wilcoxon signed-rank** | 配對非參數檢定；Thread D NG=2 用 p=0.0156 | S24 |
| **Jaccard** | 集合相似度；LOH.bed Jaccard=1.0 | S8 |

## H. 內部專案命名

| 術語 | 定義 | 首次出現 slide |
|------|------|:-:|
| **Thread D** | LOH-constrained phasing signatures 主軸（TO 層論文主軸候選）| S22 |
| **HPFineNGroups marker** | 早期視為 subclone marker；2026-04-23 降級需重驗 | S22 |
| **LOH-constrained phasing discovery** | 2026-04-22 新發現：NG=2 in Inner 93-99% same-hap（6/6 樣本）| S24 |
| **Phase 1A** | paired-pure delta F1=+0.0112 已鎖定（不受 self-phasing 影響）| - |
| **Phase 2A** | Normal Methylation Reference baseline（依 V5 重跑）| S22 |
| **五大研究目標** | InterSubMod 願景文件正式定錨；目標 1-5 詳見 §13 background | S23 |
| **V2b / V3-Fixed / V3F / V5** | 4-commit 漸進演進命名 | S17 |
| **V5max1, V5max2, V5max3** | V3F→V5 重分配最戲劇化的 3 個位點（chr19:4639528 等）| S18 |
| **SP1, SP2, SP3** | Self-phasing extreme 3 個位點（HP2:HP1=113:0 等）| S19 |

## I. 樣本相關

| 術語 | 定義 |
|------|------|
| **HCC1395** | 主驗證樣本（5kHz）；V3F → V5 audit 主用；其他 6 樣本為後續 |
| **HCC1395_DORADO** | DORADO basecaller 版本（DORADO replicate）|
| **HCC1937, HCC1954, H1437, H2009, COLO829** | 其他 6 個樣本 |
| **SEQC2** | FDA / SEQC2 cancer benchmark；HCC1395 truth set 來源 |

## J. 開放議題識別碼

| ID | 含義 |
|----|------|
| R1-R9 | 已知限制（cnLOH、AF>0.9 邊界、cherry-picked、Confidence、未 commit、7 樣本、metric 略遜、problem PS、Paired ground truth 自身用 LongPhase）|
| F1-F8 | 後續行動（commit V5、vote log、7 樣本擴展、master×flag、manifest、cnLOH 評估、trio-phased、100 隨機位點 cross-validate）|

詳見 `00_background_context.md` §14, §15。
