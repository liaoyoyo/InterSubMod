<!--
建立: 2026-06-19
報告類型: subclone 分析+判讀 完整文獻框架 — clone/subclone 特徵 + confounder + DNA 確認法 + 甲基虛實 + SNV-only/LOH-CNV 困難 + 缺口映射
任務類型: D handoff（論文 Ch1-Ch2「why-hard / 為何需要正交訊號」漏斗底座）
狀態: 文獻映射（L3）— evidence_quote 為 agent WebFetch（ASCAT/Gerstung PDF 圖像親讀=較可靠）；投稿前 PI 須從 PMC OA 全文逐字核對 + 過 /citation-verification
🔴 用戶硬要求: **只映射研究缺口，不下任何「ISM 的貢獻/創新」結論**；ISM 對應欄僅客觀標「對應哪面向 / 已處理 or 未處理」。
data_sources:
  - workflow wil25g7bu（5 researcher agents gap search + synthesis）
  - 既有: docs/method_comparison/20260618_distance_matrix_cluster_validation_methods_for_methylation_01.md / 20260615_subclone_external_validation_and_chr2_18M_AI_explainer_01.md / docs/methodology/20260618_ism_method_soundness_validation_01.md
  - ISM 內部: memory project_apriori_subclone_classification_model / project_tumor_only_axis_negative_subclone_classification / project_ism_complete_tpfpfn_existence_cis
provenance_note: 不產新數字；每筆附 citation。困難類 claim 附該文 limitations/intro 原句。甲基 claim 沿用既有虛實判準。
-->

# Subclone 分析+判讀 完整文獻框架（clone/subclone 特徵 · confounder · DNA 確認 · 甲基虛實 · SNV-only/LOH-CNV 困難 · 缺口映射）

> **用途**：論文 Ch1-Ch2「為何 subclone 重建難、為何需要正交訊號」的一站式文獻漏斗。**本檔只映射研究空白，不下 ISM 貢獻結論**（用戶要求）。
> **一句話**：subclone 重建在純 DNA(SNV±LOH/CNV)下有**系統性的不可識別性**（multiplicity、purity/ploidy、neutral-vs-selection、低-AF floor）→ 加 CNV **不消除**困難；甲基提供與 SNV/CNV 正交的 subclone/lineage 訊號**在 single-cell 是 REAL**，但 **bulk read-level 直接補 SNV-only 不足 = 文獻未建立**。

## §1 哪些特徵 = clone vs subclone（特徵 ledger）
| 特徵 | clone(truncal/MRCA) | subclone(branch) | 如何量化 | ISM 對應（客觀）|
|---|---|---|---|---|
| **SNV CCF** | CCF≈1.0（全細胞）| 中/低 CCF（子集）| VAF→CCF=f·m/(ρ·N_T+2(1−ρ)) | ISM 走 read-level haplotag 群集，非 VAF/CCF cluster（未做 CCF 推斷）|
| **mutation multiplicity** | — | 決定 CCF 的關鍵隱變數 | 須估非觀測（見 §5/GAP-A）| ISM 未處理 multiplicity |
| **LOH-phasing** | clonal LOH | subclonal LOH | BAF/logR + haplotag | ISM 用 LOH BED + HP-ratio（context）|
| **CNA** | clonal CNA | subclonal CNA（低 fraction 難）| allele-specific CN | ISM CN 為 context 非 driver（未完整解）|
| **mutational signature** | 早期/clonal signature | 晚期/subclonal signature | SBS time-resolved | ISM 未用 |
| **甲基（germline-HP 層）** | — | — | read×CpG 距離 | ISM **強**（HP-AUC normal median 0.788）|
| **甲基（within-HP somatic 層）** | — | 候選 subclone 訊號（弱/未定）| HP-axis Δβ | ISM **undetermined**（G-B gate 未跑）|

## §2 影響判讀的 confounder / noise（ledger + ISM 已處理/未處理）
| confounder | 為何影響判讀 | ISM 現況（客觀）|
|---|---|---|
| **CN-dosage** | VAF 受 copy number 扭曲 | ISM HP-axis 控 CN；CN 本身未完整解 |
| **purity/ploidy 耦合** | 同 BAF/logR 多解（見 §6）| ISM 未做 purity/ploidy 推斷 |
| **multiplicity 歧義** | clonal 被誤判 subclonal（§5）| 未處理 |
| **coverage / read-sampling floor** | 低-CCF 訊號被噪聲淹沒 | ISM min_reads=5 gate；sens~4% 受同 floor |
| **caller VAF<2% FP floor** | 上游變異集低 VAF 不可靠 | ISM 用已 call 變異，未改 caller floor |
| **over-dispersion（計數）** | binomial 殘餘離散→假群 | ISM per-CpG Fisher 仍 binomial-family（KB 必修）|
| **over-clustering / k 過估** | 無結構切成假群 | ISM read-read 層有 PERMDISP/check_dispersion；計數層未處理 |
| **homopolymer（ONT）** | 偽變異/錯 call | ISM per-position 手動審（chr2:18M pos4）|
| **dispersion（假結構）** | PERMANOVA 顯著≠型態異 | ✅ ISM check_dispersion 已實作 |
| **技術雜訊偽冒異質性** | tumor-only ~69% FP；discordant 34-80% 為雜訊 | ISM tumor-only 非監督 = NEGATIVE（已關閉）|

## §3 他人如何用 DNA 確認 clone/subclone + 移除雜訊
- **確認**：VAF→（CN+purity 校正）→CCF → mutation clustering（PyClone-VI/DPClust/SciClone，beta-binomial）→ phylogeny（PhyloWGS/Canopy）→ benchmark（Salcedo DREAM SMC-Het 7-task）。
- **移除雜訊（文獻共識組合）**：beta-binomial over-dispersion 建模、K 上界 heuristic（PyClone-VI；TRACERx K=40/PCAWG K=20）、neutral-tail 移除（MOBSTER）、提高深度（>60×, NRPCC>10）、low-mappability mask、matched-normal、**多區域/多 caller/正交驗證**（單管線不可自證）。
- **驗證判讀**：simulated 已知 phylogeny（51 腫瘤）+ multi-region 交叉 + bootstrap。

## §4 甲基虛實 + 助繪演化（接 20260618/20260615 ledger）
- **REAL（single-cell）**：MethylTree（lineage，~100% 是 supervised Q metric 非訊號率）、Sgootr（距離式甲基樹，距離式重建**非新穎**）、Gaiti（SF3B1→epimutation clade P=7.4e-9）、Epiclomal（epiclone transcend CN-clone）、**EPI-Clone Nature 2025**（epimutation 作 single-CpG digital barcode 補 SNV 稀疏，>230k cells）。
- **REAL（bulk-clock）**：evoflux fCpG 時鐘推演化時間（自承 1,610/1,976 癌無 subclonal/independent clone）。
- **🔴 INFLATED/未建立**：bulk read-level 甲基直接重建 SNV-subclone；甲基當 variant filter（DEAD）。

## §5 只用 SNV 的困難（problem-statement + 原句）
1. **非唯一性-tree**：少樣本多 lineage → 多 phylogeny 同等相容。*"multiple phylogenies are usually consistent with a given set of lineage CPs"*（Tarabichi 2021, Nat Methods, PMC7867630）
2. **非唯一性-cluster/clone 數**：CP 相近 → VAF 重疊被併群。*"clustering algorithms will merge these lineages without further information"*（同上）
3. **非唯一性-演化模型**：1/f 線性 neutral vs selection 不可識別。*"models with selection can also lead to a linear relationship between M and 1/f"*（Balaparya & De 2018, Nat Genet, DOI 10.1038/s41588-018-0217-6）
4. **multiplicity 歧義**：VAF=purity×prevalence×multiplicity 複合；CCF 公式 m 必估非觀測。*"allelic prevalence... is a compound measure"*（Roth 2014 PyClone, PMC4864026）
5. **低-AF floor①（read-sampling）**：*"As NRPCC increases, the signal of true CCF peaks becomes clearer relative to read-sampling noise"*（Tarabichi 2021）
6. **低-AF floor②（caller）**：*"detection of variants at VAFs <2% is affected by a high risk of a false positive result, regardless of the coverage depth"*（Petrackova 2019, Front Oncol, PMID 31552176）
7. **低-AF 偵測偏差定量**：observable subclone 內平均 14% SNV 低於偵測限；CCF<30% subclone 平均漏 21%（Tarabichi 2021）
8. **breadth 梯度**：*"WES detects ~50-fold fewer SNVs... hindering subclonal lineage detection; Targeted sequencing... rarely supports meaningful subclonal reconstruction"*（Tarabichi 2021）⚠ ~50-fold 投稿前核

## §6 加 LOH+CNV 仍困難嗎？→ 仍困難（problem-statement + 原句）
1. **purity-ploidy 不可識別**：*"it is always possible to posit a larger copy number and smaller CP that explain the BAF and logR equivalently well"*（Tarabichi 2021）
2. **WGD 等價解**：任一 CNA 重建 ≡「copy number 加倍 + purity 降低」（Tarabichi 2021）
3. **CN 已知仍不足定 multiplicity**：*"Knowledge of copy numbers, even allele-specific... is insufficient"* → CCF 無法僅從 DNA 辨識（Satas DeCiFer 2021）
4. **subclonal CNA/LOH 低 fraction 漏檢**：*"diminished statistical signals... false negative detection"*（Ha TITAN 2014）
5. **aneuploidy+admixture 雙重未知**：*"have proven difficult, because tumors... deviate from diploid and contain... multiple cell populations"*（Van Loo ASCAT 2010, PNAS）
6. **CN 校正污染下游 clustering**：未完全解 subclonal CN 歧義 → 偽群；建議只用 normal-diploid 區 SNV（Tarabichi 2021）
7. **SNV-tree vs CNA-tree 整合不一致**：CNA 系統發生更難且為 SNV 解讀 confounder（Fu/TUSV-ext 2022）
8. **大 cohort 亦僅 partial**：*"partial reconstruction... some clonal variants could not be timed"*（PCAWG/Gerstung 2020, Nature）

## §7 甲基能否「有效處理」？（正交性虛實 + regime 邊界）
- **正交性 REAL（但 100% single-cell）**：甲基提供與 SNV/CNV **正交**的 subclone/lineage 訊號——Gaiti（甲基不依賴 SNV 即重現遺傳分群，P=7.4e-9，但 self-caveat=共祖結構非 driver DMR）、Epiclomal（捕捉 CN 以外結構）、EPI-Clone（SNV 太稀疏，甲基 barcode 補上）。→ 「甲基在 SNV/CNV 偵測不到處仍有訊號」**概念上 REAL**。
- **🔴 bulk read-level 邊界（關鍵）**：bulk DNA 甲基平均遮蔽 subpopulation；long-read 甲基至今多 bulk 解析度；最新 long-read subclone 研究靠 physical 預分選 → **「bulk read-level 甲基直接 subclone 重建 / 補 SNV-only 不足」在文獻=未建立；bulk regime 若直接 claim=INFLATED**。
- **ISM 對應（客觀）**：ISM 在 bulk ONT read-level regime；甲基-subclone 正交性**是否實證=未**（G-B gate 未跑）；characterization 非 reconstruction。**（不下貢獻結論）**

## §8 研究空白客觀映射（GAP-A~L；不下 ISM 貢獻結論）
| GAP | 內容 | 主要源 |
|---|---|---|
| A | multiplicity/copy-number-per-mutation 非可識別 | Roth 2014; Tarabichi 2021; Satas DeCiFer 2021 |
| B | tree/cluster/演化 三層非唯一性 | Tarabichi 2021; Balaparya 2018 |
| C | 低-AF 雙重 floor（read-sampling + caller VAF<2%）| Tarabichi 2021; Petrackova 2019 |
| D | 定序 breadth 解析力梯度（WGS>WES>panel）| Tarabichi 2021 |
| E | CNA/LOH 加入仍困難（purity-ploidy/WGD/低 fraction/整數 CN 亟需）| Tarabichi 2021; TITAN; ASCAT |
| F | SNV-tree↔CNA-tree 整合不一致；大 cohort 僅 partial | Fu 2022; Gerstung 2020 |
| G | VAF noise model 選擇=群數偏誤根因（→beta-binomial）| Tarabichi 2021; Roth 2014 |
| H | over-clustering/k 過估（CC 切假群）| Şenbabaoğlu 2014; Salcedo 2020 |
| I | 雜訊移除已成共識組合 | PyClone-VI; MOBSTER |
| J | 技術雜訊偽冒異質性（tumor-only ~69% FP）| Shi 2018 |
| K | **甲基正交 subclone 訊號 REAL 但 100% single-cell；bulk read-level 未建立** | Gaiti; Epiclomal; EPI-Clone |
| L | （庫缺口）8 篇 why-hard 一手源缺卡（見下節） |  |

> **客觀研究空白總述**：純 DNA(SNV±LOH/CNV) 的 subclone 重建有系統性不可識別性 + 雙重偵測 floor + 技術雜訊偽異質性（GAP-A~J）；甲基正交訊號在 single-cell REAL、bulk read-level 未建立（GAP-K）。**這些是客觀存在的空白；是否、如何由任何特定方法（含 ISM）填補，本檔不做結論。**

## §9 引用紀律 + 庫缺口
- 全部 evidence_quote = agent WebFetch **L3**（ASCAT/Gerstung PDF 圖像親讀較可靠）；投稿前 PI 須 PMC OA 全文逐字核（尤 Tarabichi ~50-fold、N_T 定義）+ 補 Williams 2016 / DeCiFer / EPI-Clone 全文。過 `/citation-verification` 入 .bib。
- **新增 external_validation 卡（8 篇 why-hard 一手源，本輪確認 ABSENT）**：dentro-pcawg11-2021 · satas-decifer-2021 · senbabaoglu-2014 · petrackova-2019 · shi-wes-itgh-2018 · balaparya-de-2018 · ha-titan-2014 · vanloo-ascat-2010。
