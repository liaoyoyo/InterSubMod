<!--
建立時間: 2026-06-09 (workflow 完成 2026-06-10)
狀態: in_progress (method-comparison study, part 2/6 — 外部工具 survey)
報告類型: external_tool_survey
受眾: PI · 外部協作者 · 論文 Related Work 骨架
framework: MECE 分類(依「分析軸」) + 四層分層揭露
data_sources:
  - _assets/workflow_raw_result.json (workflow wf_9e64e51d-03d, 24 agents, 2.3M tok, 2026-06-10 完成, 全程 WebFetch 實證)
  - _assets/survey_digest.md (萃取版)
provenance_note: 每筆工具方法/統計/口徑皆由 subagent WebSearch+WebFetch 自一手來源(官方文件/論文/Bioconductor/GitHub)驗證, 附真實 PMID/DOI/URL; 標 ✗UNVERIFIED 者為無法 fetch 確認。本檔為文獻 paper-derived 整理, 非本專案新跑分析。
-->
<!-- provenance-verified: 外部工具事實由 workflow wf_9e64e51d-03d 之 WebFetch 一手來源查證, 識別碼存 _assets/workflow_raw_result.json; 撰寫與查證在不同 batch(§13.0)。 -->

# 02 — 外部工具 survey（ISM 相關甲基化方法地景）

> **這份是什麼**：方法比較研究的**第 2 部分** —— 把網路上與 ISM 任務相關的甲基化軟體/論文，**依「分析軸」分類盤點**（82 個工具/變體，全程 WebFetch 實證來源）。重點工具：**ONT modkit、二代短讀癌症 DMR caller、甲基位點關聯方法**。
> **基準對照**：我們的方法定義見 `01_ism_method_spec_from_source.md`。

---

## L0 — 一句話地景（最關鍵的發現）

外部甲基化方法可乾淨地分成 **6 條軸**。ISM 主要站在**「軸 C：read-to-read 距離 + clustering」**，而**絕大多數主流工具（含 modkit、DSS、methylKit）站在「軸 A：per-position 甲基率差」**。兩者**問的問題本質不同**：軸 A 問「這裡兩組的平均甲基率差多少」，ISM（軸 C）問「reads 是否在結構上分成不同亞群/單倍型」。**沒有任何一個工具佔據 ISM 的完整組合**（軸 C 距離矩陣 + PERMANOVA 顯著性 + normal-anchored somatic cis-test + 5mC/5hmC 分軌 + LOH/CN 耦合），但**單獨拆開的每個元件都已有人做**。

---

## L1 — 6 條分析軸 + 各軸代表工具（一表定位全局）

| 軸 | 問的問題 | 代表工具 | ISM 在這條軸？ |
|----|---------|---------|---------------|
| **A. per-position 甲基率差**（DMR caller）| 「兩組在這個位點/區域平均甲基率差多少」| **modkit dmr** · DSS · methylKit · metilene · BSmooth/bsseq · dmrseq · DMRcate · comb-p · MOABS · RnBeads · eDMR | 部分（CORE 4 Fisher / Δβ 屬此，但非主引擎）|
| **B. within-read CpG-CpG 共甲基/LD**（沿基因組、同分子）| 「同一條 read 上 CpGᵢ 與 CpGⱼ 是否連動」| **MHB/MHL (Guo 2017)** · mHapTk · Metheor LPMD · DAMEfinder tuple · CAMDA | ❌（CORE 4 NME/epipoly 是 disorder 版，非 LD；可學習）|
| **C. between-read read-read 距離 + clustering**（跨分子）| 「reads 是否分成結構上分離的亞群」| **ISM** · cvlr · **ASMS (Raineri 2024)** · modbamtools(HDBSCAN) · MeConcord · (qFDRP 同 kernel 但塌成 scalar) | ✅ **主引擎（CORE 1-3）** |
| **D. cohort 跨樣本相關網路**（跨個體）| 「跨樣本間哪些區域甲基共變/驅動表現」| CoMeBack · Comethyl(WGCNA) · ELMER · MethylMix | ❌（單樣本設計，不適用）|
| **E. disorder / WSH scalar**（失序量）| 「甲基有多隨機/失序」| PDR · epipolymorphism · entropy · FDRP · **qFDRP** · MHL · informME · CpelNano · methclone · epihet · epiCHAOS | 部分（CORE 4 內建作對照）|
| **F. 長讀 phasing / ASM**（單倍型）| 「用甲基 phase，或 phase 後算 ASM」| MethPhaser · NanoMethPhase · NANOME · HapBridge · longphase modcall · modkit pileup --phased | 反向（ISM **消費** HP tag 去量 ASM 結構，不靠甲基 phase）|

> 🔑 **三個最像 ISM 的工具**（投稿時必正面對照）：
> 1. **cvlr**（軸 C，最近）：ONT read-level、不需 phasing、Bernoulli-mixture EM 把 reads 分群 → 但**無 read-read 距離矩陣、無顯著性檢定、germline**。
> 2. **ASMS (Raineri 2024 bioRxiv)**：ONT modBAM、不需 phasing、EM 2-component 把 reads 依甲基 pattern 分群找 ASM → **單一最像 ISM 的已發表工具**（⚠ digest 標 UNVERIFIED，須補查）。
> 3. **qFDRP (Scherer 2020)**：用**和 ISM NHD 相同**的 normalized-Hamming read-pair kernel → 但塌成 per-CpG scalar、丟掉幾何、無 label/PERMANOVA。

---

## L2 — 各軸工具細節（method / granularity / 平台 / 引用，全 WebFetch 實證）

### 軸 A — per-position 甲基率差（DMR / DMC caller）

| 工具 | 平台 | 方法核心 | 統計 | 引用 |
|------|------|---------|------|------|
| **modkit dmr** ⭐ (ONT 官方, #1 對照) | ONT 長讀 | bedMethyl(pileup 聚合)上比兩組；**Bayesian conjugate-prior marginal likelihood ratio**（非 frequentist LRT）；Beta-Binomial(二態)/Dirichlet-Multinomial(多態)+Jeffreys prior；single-site 加 MAP-based effect-size p；region 報 **Cohen's h**；可選 2-state HMM 分段 | Bayes-factor 風格 LR score + MAP p + Cohen's h | nanoporetech.github.io/modkit/dmr_scoring_details.html（無正式論文）|
| **DSS** (黃金標準短讀) | 短讀 WGBS/RRBS | Beta-Binomial 計數模型 + **lognormal dispersion 經驗貝氏 shrinkage**（跨 CpG 借強度）+ per-CpG **Wald**；無 replicate 用 500bp 滑窗 pseudo-replicate；2016 版 arcsine-link beta-binomial regression | beta-binomial + EB shrinkage + Wald | Feng/Conneely/Wu 2014 NAR PMID 24561809; Park/Wu 2016 Bioinf |
| **methylKit** | 短讀 (+array) | 多樣本→logistic regression(Chi²/F)；單樣本/組→**Fisher exact**；可選 overdispersion 'MN' | logistic / Fisher；SLIM q-value | Akalin 2012 GB PMID 23034086 |
| **metilene** | 短讀/任意 level | **circular binary segmentation** on 組間均差 → de-novo DMR → **2D-KS** + Mann-Whitney | CBS + 2D-KS | Jühling 2016 GR PMID 26631489 |
| **BSmooth/bsseq** | 短讀 (低覆蓋友善) | locfit **局部似然平滑** → signal-to-noise **t-stat**（平滑均差/平滑 SD）；需 ≥2 replicate | 平滑 t-statistic | Hansen 2012 GB PMID 23034175 |
| **dmrseq** | 短讀 | 兩階段：候選 DMR → **GLS + nested AR(1)** 相關誤差（arcsine prop）→ **region permutation** | GLS-AR(1) + 標籤置換 | Korthauer 2018 Biostatistics PMID 29481604 |
| **DMRcate** | array (+WGBS) | limma moderated-t² → **Gaussian kernel** 沿基因組平滑 → Satterthwaite scaled χ² | limma + kernel smoothing | Peters 2015 Epigen Chrom PMID 25972926 |
| **comb-p** | 任意 p-value | **Stouffer-Liptak-Kechris** 用空間自相關加權合鄰近 p；Sidak region 校正 | 空間自相關 p-combiner | Pedersen 2012 Bioinf PMID 22954632 |
| **MOABS** | 短讀 | Beta-Binomial 階層 + 全基因組 EB prior → **CDIF**（credible 甲基差，非 p）| credible-interval effect | Sun 2014 GB PMID 24565500 |
| **RnBeads 2.0** | array+WGBS | limma 階層線性模型(M-value)+ **max-of-3-ranks**（均差/比值/p）；region = Fisher 合 p | limma + rank fusion | Müller 2019 GB; Assenov 2014 NatMeth |
| **eDMR** | 短讀 | 雙峰分佈定分段距離 + weighted Z-test 合 p（methylKit 上游）| weighted Z (Stouffer) | Li 2013 BMC Bioinf PMID 23735126 |

> ⚠ **重要事實（modkit dmr 的 haplotype 限制）**：**modkit dmr 內部沒有 `--haplotype` flag**（跨 intro_dmr/advanced_usage/GitHub book 原始碼確認）。要做單倍型 DMR = 先 `modkit pileup --phased`（或 `--partition-tag HP`）拆 hp1/hp2 bedMethyl，再 `dmr pair -a hp1 -b hp2`。**haplotype 是資料切分，不是模型項** — 這正是 ISM 把 HP1 vs HP1-1 做 somatic-controlled 軸的差別所在。

### 軸 B — within-read CpG-CpG 共甲基 / LD（同分子、沿基因組）

| 工具 | 平台 | 方法核心 | 引用 |
|------|------|---------|------|
| **MHB / MHL (Guo 2017)** ⭐ | 短讀 BS | **Methylation Haplotype Block** = 相鄰 CpG 間甲基 **LD r²** 超閾值合併（r²=D²/[pA(1−pA)pB(1−pB)]，2×2 跨 reads；seed ≥2 CpG r²>0.5）；**MHL** = 長度加權「完全甲基化子串」比例。147,888 blocks；cfDNA tissue-of-origin | Guo 2017 Nat Genet PMID 28263317 |
| **mHapTk** | 短讀 | MHB/MHL/R²(signed)/PDR/entropy 的 Python 實作；R²>=0.5 p<=0.05 | Ding 2022 Bioinf PMID 36179079 |
| **Metheor LPMD** | 短讀 BS (Rust) | **LPMD** = 同 read 上固定距離 d 的 CpG-pair 不一致比例（去除 PDR read-length bias，預設 2-16bp）；one-sweep；快 300× | Lee 2023 PLoS CB PMID 36940213 |
| **DAMEfinder tuple** | 短讀 BS | tuple 模式 ASM score = CpG-pair **log-odds (MM·UU)/(MU·UM)** × Beta 權重（同 read 共甲基連動）；150bp cap | Orjuela 2020 Epigen Chrom PMID 32487212 |
| **CAMDA** | 短讀 BS | 同 read 內甲基+未甲基**並存比例**(concurrence) | Shi 2021 Nat Commun |

### 軸 C — between-read read-read 距離 + clustering（**ISM 主場**）

| 工具 | 平台 | 方法核心 | vs ISM 一句話 | 引用 |
|------|------|---------|--------------|------|
| **ISM (我們)** | ONT 長讀 癌症 | read×CpG 矩陣 → read-read 距離(NHD/L1/Bernoulli/Pearson/Jaccard) → UPGMA+silhouette → **PERMANOVA** + normal-anchored cis-test | — | `01_*.md` |
| **cvlr** ⭐最近 | ONT 長讀 germline | **Bernoulli-mixture EM** 把 reads 分 k 群（soft posterior，無需 phasing）| 模型式 vs ISM 距離式；**無距離矩陣、無顯著檢定、germline、固定 k 無 BIC** | Raineri 2023 Bioinf Adv PMID 36726731 |
| **ASMS** ⭐最像 (⚠UNVERIFIED) | ONT modBAM | EM 2-component 依甲基 pattern 分 reads 找 ASM，**不需 phasing** | 單一最像 ISM 的已發表工具；ISM 多了距離矩陣+PERMANOVA+normal-anchor | Raineri 2024 bioRxiv 2024.12.18.629129 |
| **modbamtools cluster** | ONT 長讀 | **HDBSCAN** 依甲基 state 分 reads，回傳 cluster 數 | 回傳數量非距離矩陣，無顯著性 | Razaghi 2022 bioRxiv |
| **MeConcord** | 短讀 BS | read-read + CpG Hamming concordance | pairwise Hamming 但無 label PERMANOVA | Bioinformatics 2022 |
| **(qFDRP)** | 短讀 BS | **與 ISM NHD 相同的 normalized-Hamming read-pair kernel** → 但平均成 per-CpG scalar | 同 kernel、丟幾何、無 label | Scherer 2020 NAR PMID 32338758 |

### 軸 D — cohort 跨樣本相關 / 表現整合（單樣本 ISM 不適用，僅對照）

| 工具 | 方法核心 | 引用 |
|------|---------|------|
| **CoMeBack** | array 跨個體相鄰 probe 相關 → co-methylated regions（背景 CpG 密度過濾）| Gatev 2020 Bioinf |
| **Comethyl** | **WGCNA** region 甲基跨樣本相關 → module + eigennode-trait | Mordaunt 2022 Brief Bioinf PMID 35037016 |
| **ELMER** | 遠端 DMC × 鄰近基因表現反相關 → master TF 網路（需 RNA-seq）| Silva 2019 Bioinf PMID 30364927 |
| **MethylMix** | per-gene **beta mixture(BIC)** patient 分群 + 甲基→表現 | Gevaert 2015 Bioinf |

### 軸 E — disorder / WSH scalar（ISM CORE 4 內建作對照）

| 工具/度量 | 量 | 平台 | 引用 |
|----------|----|------|------|
| **PDR** (Landau 2014) | within-read 不一致 read 比例；**明文排除 ASM 是來源**（disorder=stochastic）| 短讀 癌症 | Landau 2014 Cancer Cell PMID 25490447 |
| **epipolymorphism** (Landan 2012) | 4-CpG 16-epiallele Simpson diversity 1−Σpₖ² | 短讀 | Landan 2012 Nat Genet PMID 23064413 |
| **methylation entropy** (Xie/Landan) | 16-epiallele Shannon entropy | 短讀 | Xie 2011 NAR |
| **FDRP / qFDRP** (Scherer 2020) | read-pair (binary/normalized-Hamming) per-CpG 失序；**qFDRP 最 coverage-robust** | 短讀 | Scherer 2020 NAR PMID 32338758 |
| **informME (NME/dNME)** | 1D **Ising/CPEL** 區域聯合分佈熵；低覆蓋(5×)強 | 短讀 WGBS | Jenkinson 2018 BMC Bioinf PMID 29514626 |
| **CpelNano** ⭐同平台 | informME 的 **ONT 長讀**版（Ising-HMM）；breast cancer 差異 landscape | ONT 長讀 | Abante 2021 Sci Rep PMID 34732768 |
| **methclone / epihet** | 4-CpG epiallele 組合熵 shift / ITH | 短讀 癌症 | Li 2014 GB; Chen 2021 Sci Rep |
| **epiCHAOS** | single-cell chance-corrected Jaccard pairwise | scATAC/scWGBS | Kelly 2024 GB PMID 39623476 |

### 軸 F — 長讀 phasing / ASM（單倍型；ISM 反向消費 HP）

| 工具 | 方向 | germline/cancer | 方法 | 引用 |
|------|------|----------------|------|------|
| **MethPhaser** | 甲基 → phase | germline 株 | SNV-first，homozygous gap 用甲基 Wilcoxon 投票延伸 phase | Fu 2024 Nat Commun PMID 38909018 |
| **NanoMethPhase** ⭐ HP 軸最近 | phase → ASM | germline imprinting | phased VCF→haplotag→**DSS** beta-binomial HP1 vs HP2 DMR | Akbari 2021 GB PMID 33618748 |
| **NANOME** | consensus caller + ASM | germline | XGBoost 共識 caller + Clair3 phase + methylKit χ² ASM | Liu 2025 |
| **HapBridge** (⚠UNVERIFIED) | 甲基 → 修 switch error + bridge | germline | HSM 位點 inter-site 一致性 + 迭代 tag unphased reads | Chen 2025 bioRxiv |
| **longphase modcall** | 甲基 co-phase | germline | binarize(mod>0.8/unmod<0.2，**與 ISM 同閾值**)→site→graph co-phase SNP/indel/SV | Lin 2022 Bioinf |

### 提取 / 視覺化 utility（非分析引擎，但 pipeline 常用）

| 工具 | 角色 | 引用 |
|------|------|------|
| **modkit pileup** | modBAM→bedMethyl（**--phased** 拆 HP）；ISM 用 read×CpG 矩陣取代此聚合 | ONT 官方 |
| **modbam2bed** | htslib pileup→per-site %（一次一 hap）| epi2me-labs |
| **methylartist** | 單分子 locus 圖 + 輸出 DSS input；**ISM 沒有的視覺化可借** | Cheetham 2022 Bioinf |
| **NanoMethViz** | R spaghetti/aggregate 圖 + 接 bsseq/DSS/edgeR | Su 2021 PLoS CB |
| **pb-CpG-tools** | PacBio HiFi per-site（5mC only，自動 hap1/hap2）| PacBio |

### 癌症專用（array/短讀；多為 cohort，對照用）

| 工具 | 功能 | 引用 |
|------|------|------|
| **conumee** | array 總強度→**CNV 分段**（丟甲基訊號；ISM 反而把 CN 當 covariate 耦合）| Hovestadt/DKFZ |
| **InfiniumPurify** | iDMC beta KDE → **tumor purity** 標量 | Zheng 2017 GB |
| **CancerLocator** | cfDNA Beta-deconvolution → tumor fraction + tissue-of-origin | Kang 2017 GB |

### Δβ / effect-size 方法學（與我們 Δβ 直接對照 — 見 04/05）

| 方法 | 要點 | 引用 |
|------|------|------|
| **beta vs M-value** (Du 2010) | raw beta 在 0/1 端 **heteroscedastic**；建議 M-value 檢定、report beta | Du 2010 BMC Bioinf PMID 21118553 |
| **raw Δβ + |Δβ|≥0.1** | 社群通用 effect-size 閾值（methylKit 25%）| 慣例 |
| **Xie 2019** | 檢定尺度(M/logit)與 report 尺度(beta)須一致；naive beta 在極端偏 | Xie 2019 Bioinf |
| **arcsine vs logit** (Warton-Hui) | proportion 變異穩定化；ISM 用 rank-based Wilcoxon 是合法第三選 | Warton-Hui 2011 Ecology |
| **CPEL/Ising ASM** (Abante 2020) | read-level 多-CpG pattern 比兩 allele（MML/NME/JSD）；**ISM NME/epipoly 幾乎 1:1 對應** | Abante 2020 Nat Commun |
| **Cohen's d (paired d_z)** | ISM 用 Cohen's d 配 Wilcoxon 符合最佳實踐；應明說是 paired d_z | Lakens 2013 |

---

## L3 — 完整來源
- 全 82 工具結構化資料（含每筆 method/stat/hap/cancer/url/citation/verified）：`_assets/workflow_raw_result.json`
- 萃取版：`_assets/survey_digest.md`
- 查證方式：workflow `wf_9e64e51d-03d`（24 agents：10 survey + 8 deep-read + 5 cross + 1 synth-critic；每筆 WebSearch→WebFetch 一手來源；無法 fetch 標 ✗UNVERIFIED）

## ⚠ 已標記的不確定 / 須補查
- **ASMS (Raineri 2024)** = 最像 ISM 的已發表工具，但 digest 標 UNVERIFIED（bioRxiv，須親讀確認 EM/CLI）— **投稿前必補查**。
- **cvlr GitHub repo 404**（github.com/EmanueleRaineri/cvlr，作者 0 public repos @2026-06-10）— CLI/實作語言未能一手確認，僅論文 PMC 文字。
- **HapBridge** preprint 精確統計檢定未能 fetch。

## 與後續文件
→ `03` 方法對照矩陣（6 核心逐格）｜`04` 我們結果×外部數據交叉｜`05` 改進/學習｜`06` Phase B benchmark
