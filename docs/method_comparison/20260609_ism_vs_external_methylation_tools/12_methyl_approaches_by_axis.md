<!--
建立時間: 2026-06-11
狀態: in_progress (重新確認 — 業界甲基數據處理 4 軸地景 + 難點 + ISM 對照)
報告類型: methyl_approaches_landscape_by_axis
受眾: PI · 理解業界方向 · 論文 Related Work · ISM 定位
framework: MECE 4 軸分類(用戶框架) + 每軸(canonical/tools/ISM 差別/難點) + 共同難點 + Verdict
data_sources:
  - _assets/axis_reconfirm_raw.json + axis_reconfirm_digest.md (workflow wf_22ba1b9b-db4, 7 agents, WebFetch 一手實證, 2026-06-11)
  - 01_ism_method_spec_from_source.md (ISM 源碼)
provenance_note: 外部工具/canonical pipeline/難點皆 WebFetch 一手來源(附 PMID/DOI); ISM 側引源碼。撰寫與查證不同 batch(§13.0)。標 UNVERIF 者為無法 fetch。
-->
<!-- provenance-verified: 外部事實引 workflow wf_22ba1b9b-db4 WebFetch 一手(_assets/axis_reconfirm_raw.json); ISM 引 01 源碼; 撰寫與查證分批。 -->

# 12 — 業界甲基數據處理「4 軸地景」+ 難點 + ISM 對照

> **這份回答你**：業界處理甲基數據有哪些方向 —— 依你關心的 **4 條軸**組織（① tag 分群→比 tag 甲基平均差 ② read 層甲基距離 ③ 位點層甲基差 ④ 甲基判別方式），每軸列 **canonical 做法 + 工具 + ISM 在哪/差別 + 難點**，最後 **8 個大家共同的難點** + ISM 獨有/欠缺。**全程 WebFetch 一手實證。**

---

## L0 — 一張圖看懂（最關鍵 4 句）

1. **軸1（tag 分群→比平均）是業界主流**，但有個**天花板**：它比的是 **HP1 vs HP2 = germline 兩條親本軸**（混 imprinting + 隨機單等位 + somatic），**不是 somatic 軸** —— 顯著 ≠ tumor-specific。
2. **ISM 的結構性差異 = somatic-controlled HP1-vs-HP1-1 軸 + read-level PERMANOVA 結構檢定 + normal-anchored cis-test**（量測而非只過濾 confound）—— 這組合無人佔（軸2 最完整實例）。
3. **但 ISM 也有明確欠缺**：無 spatial region-FDR、無 beta-binomial、無 phasing-free 模式、**無 copy-number/purity deconvolution（CAMDAC 有）**、無 cohort/tissue-of-origin。
4. **8 個共同難點所有人都卡**：覆蓋 / 異質性(epipolymorphism) / germline-somatic confound / copy-number / 單樣本 power / 純度 / phasing+LOH / 5mC-5hmC —— **ISM 對其中 2 個（somatic 軸 + cis-test）有獨特解法，其餘共卡。**

---

## L1 — 4 軸總表

| 軸 | canonical 做法 | 代表工具 | ISM 在哪 | 最大難點 |
|----|---------------|---------|---------|---------|
| **① tag 分群→比平均** | phase SNV→haplotag→拆 per-HP 甲基→比 HP1 vs HP2 mean(DMR)| NanoMethPhase+DSS · **modkit pileup --phased+dmr** · pycoMeth · O'Neill-POG · NANOME | **下游消費 HP tag；但用 HP1-vs-HP1-1 somatic 軸**（非 HP1-vs-HP2）| HP1-vs-HP2 是 germline 軸；LOH 毀 het anchor |
| **② read 層甲基距離** | read×CpG 矩陣→read-read 距離 or model-EM→分群 | **ISM** · cvlr · ASMS · modbamtools · qFDRP/FDRP · MeConcord | **主場、最完整實例**（唯一有距離矩陣+clustering+PERMANOVA+confound 控）| 多數無顯著性檢定；O(N²) |
| **③ 位點層甲基差** | per-CpG count 模型檢定+MT 校正+spatial 聚合→DMR | **DSS** · methylKit · metilene · BSmooth · dmrseq · DMRcate · comb-p · MOABS · modkit dmr | **次要 characterization 層（Fisher+Δβ，是此軸子集且較弱）** | over-dispersion / genome-wide MT / spatial 相關 |
| **④ 甲基判別** | 依目標分流：allele/tissue/tumor-fraction/expression | DAMEfinder · CGmapTools · amrfinder · **MHB/MHL** · CancerLocator · MethylMix · ELMER | **問了無人問的 target：ASM 能否判別真/假變異？→ 誠實 NEGATIVE** | germline-somatic confound；能否判別是 open question |

---

## L2 — 逐軸細節

### ① tag 分群 → 比 tag 甲基平均差（業界主流；canonical 5 階段）
**canonical pipeline**（O'Neill-POG / NanoMethPhase）：
1. **SNV calling**（Clair3）→ germline het SNV
2. **phasing**（WhatsHap / longphase）→ phased VCF（PS phase-set）
3. **haplotag**（whatshap/longphase haplotag）→ 每 read 貼 HP:Z:1/2（無 het SNV 的 read untagged）
4. **拆甲基**（NanoMethPhase / modkit pileup --partition-tag HP）→ per-CpG β 分 HP1/HP2
5. **比 HP1 vs HP2**：per-CpG 檢定 + region 聚合 → aDMR；效果 = Δβ=mean(HP1)−mean(HP2)，gate `|Δβ|≥0.1`

| 工具 | 比 HP 平均的統計 | 引用 |
|------|----------------|------|
| **NanoMethPhase + DSS** | DSS beta-binomial **dispersion shrinkage + smoothing + Wald**（callDML/callDMR）| Akbari 2021 GB PMID 33618748 |
| **modkit pileup --phased + dmr** | Bayesian **marginal-LR + MAP p + Cohen's h**（Beta-Binomial/Dirichlet-Mult）| ONT 官方 docs |
| **pycoMeth Meth_Comp** | **Mann-Whitney U**（2 群）/ Kruskal-Wallis + BH | Snajder 2023 GB |
| **O'Neill-POG**（cancer cohort）| DSS Wald + `|Δβ|≥0.1`；**somatic 特異靠事後排除 imprinting/normal** | O'Neill 2024 Cell Genomics PMID 39406235 |
| **NANOME** | XGBoost 共識 caller + longphase + 下游 ASM | bioRxiv 2025 PMID 40631091 |

**🔑 ISM 在這軸的差別**：ISM **下游消費** 這條 pipeline 的 HP tag（用 LOH-constrained longphase-TO/-S 變體），但：
- **單位不同**：軸1 把每 HP 塌成 per-CpG **mean β** 比 HP1-mean vs HP2-mean（盲於 read-to-read 結構）；ISM 保留 read×CpG 矩陣 → read-read 距離 → **PERMANOVA 結構檢定**。
- **軸不同（最重要）**：軸1 = **HP1 vs HP2（germline 兩親本軸）**；ISM = **HP1 vs HP1-1（germline-HP1 vs 其上的 somatic 重建）= somatic-controlled 軸**。
**難點**：① **phasing 依賴**（無 het SNV→無法分群）② **LOH 毀 het anchor**（癌症最關鍵的刪除區 CDKN2A 反而最難 phase）③ **HP1-vs-HP2 是 germline 軸**（顯著 ≠ somatic，需事後排 imprinting）④ 覆蓋減半（分 HP 後每堆只剩一半深度）⑤ 純度稀釋 Δβ ⑥ 單樣本 power（per-CpG 把 read counts 當 replicate = pseudo-replication）⑦ **mean 塌掉 read-level 異質性** ⑧ copy-number 偏 β ⑨ 5mC/5hmC 混 ⑩ region 邊界靠 smoothing/HMM 參數（任意）。

### ② read 層甲基距離（ISM 主場）
**canonical**：read×CpG 矩陣（0/1/missing）→ (a) **顯式 read-read 距離**（NHD over 共覆蓋 CpG）→ 階層/HDBSCAN 分群，或 (b) **model-EM 跳過距離矩陣**（multivariate Bernoulli mixture，cvlr/ASMS）。

| 工具 | 方法 | 有顯著性檢定？| 引用 |
|------|------|:---:|------|
| **ISM** | NHD 距離矩陣→UPGMA+silhouette→**PERMANOVA pseudo-F + dispersion + LabelTest Δ** + normal-anchored cis | ✅ PERMANOVA | 源碼 |
| **cvlr** | Bernoulli-mixture **EM**（無距離矩陣，K 固定）| ⚠ perm+Fisher | Raineri 2023 Bioinf Adv PMID 36698762 |
| **ASMS** | Bernoulli-mixture EM **固定 2 群**（無距離矩陣，phasing-free）| ⚠ cluster-by-snp p | Raineri 2024 bioRxiv（preprint）|
| **modbamtools cluster** | **HDBSCAN**（無矩陣、**無檢定**，探索/視覺化）| ❌ | Razaghi 2022 bioRxiv |
| **qFDRP/FDRP** | read-**pair** normalized-Hamming → **per-CpG scalar**（無分群、無檢定）| ❌ | Scherer 2020 NAR PMID 32103242 |
| **MeConcord** | read-read Hamming（矩陣乘法，不存矩陣）→ binomial p | ⚠ binomial | Wang 2022 Bioinf PMID 35758820 |
| **Metheor** ⚠ | **within-read**（LPMD 同 read 上 CpG-CpG），**非 read-to-read** | ❌ | Lee 2023 PLoS CB PMID 36976832 |

**🔑 ISM 差別**：**唯一同時具備**（1）顯式 read-read NHD 距離矩陣為主引擎（cvlr/ASMS 用 EM 跳過矩陣；Metheor 是 within-read；qFDRP/MeConcord 做 pairwise 但不存/不分群矩陣）（2）UPGMA+silhouette 切 subclonal 群（3）**真正的群間結構檢定（PERMANOVA pseudo-F + permutation）**（4）cancer 特異 confound 控（normal-anchored cis）。**多數對手停在 per-CpG/per-bin scalar，從不檢定「這裡有沒有真讀層結構」。**
**難點**：① 多數工具**無原生顯著性檢定**（只 cvlr 有）② 覆蓋（FDRP 需 ~40x 上限、<10x 不穩；EM 需每群夠 read）③ 距離度量選擇敏感（missing/重疊長度/CpG 密度/read 長都偏距離）④ **compute O(N²)~O(N³)**（read-pair WSH 還靠 subsample → 不可重現）⑤ **距離/分群的生物意義模糊**（allele？subclone？copy？epipolymorphic 噪音？）⑥ germline-somatic confound（7 工具無一原生分）⑦ 單樣本無 replicate ⑧ bisulfite-era 出身（FDRP/MeConcord/Metheor 為短讀設計，套長讀是 off-design）。

### ③ 位點層甲基差（成熟工具軍火庫；ISM 是子集）
**canonical**：per-CpG count 矩陣 →（可選 smoothing 借鄰近 CpG）→ count 模型檢定 → genome-wide MT 校正 → spatial 聚合 → ranked DMR。

| 工具 | 統計核心 | 引用 |
|------|---------|------|
| **DSS** | beta-binomial **dispersion shrinkage + Wald**（黃金標準）| Feng/Wu 2014 NAR; Park/Wu 2016 |
| **methylKit** | logistic / Fisher + overdispersion 'MN' | Akalin 2012 GB |
| **metilene** | CBS + 2D-KS + Mann-Whitney | Jühling 2016 GR PMID 26631489 |
| **BSmooth** | 局部似然 smoothing + signal-to-noise t | Hansen 2012 GB PMID 23034175 |
| **dmrseq** | **GLS + nested AR(1)** + region permutation | Korthauer 2018 Biostatistics |
| **DMRcate** | limma moderated-t² + Gaussian kernel | Peters 2015 Epigen Chrom |
| **comb-p** | **Stouffer-Liptak** ACF-加權 + Sidak | Pedersen 2012 Bioinf PMID 22954632 |
| **MOABS** | beta-binomial 階層 + **CDIF**（credible Δ）| Sun 2014 GB |
| **modkit dmr** | Bayesian marginal-LR + Cohen's h | ONT 官方 |

**🔑 ISM 差別（誠實）**：ISM 的位點層（per-CpG Fisher+BH-FDR / Cramér's V+Cochran gate / Δβ Wilcoxon+Cohen's d）**方法學上是此軍火庫的子集，且較弱** —— ≈ methylKit-Fisher + effect-size 層（Δβ/Cohen's d 類比 MOABS CDIF / modkit Cohen's h）。**ISM 不跟這軸比 machinery**（無 AR(1)/ACF/CBS/kernel smoothing、無原則性 region-FDR、無 beta-binomial over-dispersion）—— 它把這層當**保守的 characterization/confirmation 層**（Cramér's V 被 Cochran gate 主動歸零、抑制稀疏假陽性），真正的引擎在軸2。
**難點**：① **over-dispersion**（純 Fisher/binomial 偏樂觀，每個可信工具都模 dispersion）② **genome-wide MT**（~28M CpG 相關檢定，BH 假設違反）③ **spatial 自相關**（鄰 CpG 共甲基，naive region 雙重計數）④ 覆蓋低深 ⑤ beta 0/1 端 heteroscedastic ⑥ replicate 需求 vs 單樣本 ⑦ germline-somatic/allelic confound（無 off-the-shelf 控）⑧ **copy-number/purity/ploidy**（β 是跨 copy/subclone/normal 混合，CNV 純算術造假 DMR）⑨ 異質性塌掉 ⑩ phasing 依賴（somatic 軸框架下）。

### ④ 甲基判別方式（4 個判別目標）
業界把「判別」分 4 個 target：

| target | 工具 | 判別什麼 / 統計 | 引用 |
|--------|------|----------------|------|
| **allele（ASM）** | DAMEfinder（tuple log-odds + SNP）· CGmapTools asr/ass（Fisher/t）· amrfinder（1-vs-2-allele LR, Fang 2012）| 兩等位甲基差 | Orjuela 2020; Guo 2018; Fang 2012 PMID 22523239 |
| **tissue-of-origin** | **MHB/MHL（Guo 2017）**· CancerLocator | cfDNA 組織來源（甲基 haplotype block 去卷積）| Guo 2017 Nat Genet PMID 28263317 |
| **tumor-fraction** | CancerLocator（Beta 混合 MLE）| 腫瘤負荷 θ | Kang 2017 GB PMID 28335812 |
| **expression-driver** | MethylMix（DM-value）· ELMER（enhancer-gene+TF）| 甲基→表現驅動基因 | Gevaert 2015; Silva 2019 |

**🔑 ISM 差別**：ISM 問了一個**無人直接測的 target —— 甲基/ASM 能否判別「真 vs 假 somatic 變異」（variant-QC filter）**，並給出**誠實 NEGATIVE**（四道：LOSO circularity / AUC=0.505 / strong-ASM 在 FP enriched ~5× OR=0.194 p=1.8e-28 / 5th-rank vestigial）。**這個 NEGATIVE 本身就是貢獻**（沒人測過、且 ASM 真實但 anti-discriminative）。
**難點**：① **germline-somatic confound**（SNP-anchored caller 分不出 germline imprinting vs somatic）② **異質性 vs 真 ASM**（DAMEfinder 明說「很難從異質性區分真 ASM」）③ copy-number/LOH（cnLOH 毀 1:1 等位比）④ 單樣本 vs cohort（MethylMix 需 ~100 樣本）⑤ 覆蓋/phasing（每等位需夠 read；LOH/SNP-desert 失效）⑥ bisulfite/短讀出身（無 5hmC、co-CpG 少）⑦ **能否判別 = open question**（ISM 的 NEGATIVE 是目前唯一直接答案）⑧ 純度地板。

---

## L2.5 — 🔴 8 個大家共同的難點（所有軸都卡，附一手來源）

| 難點 | 是什麼 | 影響哪些軸 | 業界怎麼處理 | 源 |
|------|--------|-----------|-------------|----|
| **1 覆蓋地板** | read-level/ASM 度量有硬覆蓋下限（FDRP/qFDRP ~10x/CpG + ≥4 共覆蓋 CpG；分 HP 後深度減半）| 全 | 設 per-metric floor；cap reads；目標 ~30x 讓每 HP ~15x | Scherer 2020 NAR |
| **2 異質性/epipolymorphism** | 腫瘤是 subclone 混合 → 50% 甲基可能是「一半細胞全甲基」非 hemi-methyl；**mean 法塌掉** | 全（尤其 ①③）| 用 read-level（非 pooled）統計 PDR/qFDRP；≥4 CpG/read | Landau 2014 Cancer Cell PMID 25490447 |
| **3 germline-somatic confound** | ASM 三來源（imprinting ~150 域 / germline cis hap-ASM / somatic）難分；**HP1-vs-HP2、REF-vs-ALT 本質是 germline 軸** | ①②③④ | **強制 matched-normal 扣基線** + 丟 imprinted 區 + allele-switch test | Do 2020 GB PMID 32594908 |
| **4 copy-number（CNA）** | CN 改變等位比 → 純算術造出 apparent ASM（4:1 amp → 75:25，零表觀改變）| 全 | 算 CNA-predicted 甲基並回歸扣除（Martin-Trujillo）；**CAMDAC 聯合模 purity+allele-CN 反推純 tumor** | Martin-Trujillo 2017 Nat Commun PMID 28883545 |
| **5 單樣本 power** | 一個腫瘤無 biological replicate → 只能估 sampling 變異非 biological | 全 | beta-binomial 用鄰 CpG 當 pseudo-replicate（DSS-single）；HP-shuffle 置換 null；matched-normal | Wu 2015 NAR（DSS-single）|
| **6 純度** | β_obs = p·tumor + (1−p)·normal → 低純度把訊號拉向 50/50，Δβ 掉到 <0.1 門檻下 | 全 | 估純度 + purity+CN 混合模型反卷（**CAMDAC**）；matched-normal anchor | CAMDAC bioRxiv 2020 |
| **7 phasing + LOH** | ASM 需 phase；het SNV desert / LOH 失 anchor（~95% autosomal CpG 可 phase，缺口在最關鍵的 LOH 區）| ①②③④ | 長讀跨更多 SNV；**甲基-based phase-block bridging（MethPhaser）** | MethPhaser Fu 2024 Nat Commun PMID 38909018 |
| **8 5mC/5hmC 混** | bisulfite 分不出（全部 5mC+5hmC 混）；5hmC 在癌常相反意義 | 全 | ONT basecaller 分軌（MM/ML 兩機率）；oxBS/TAB-seq | Bohrer 2024 BMB Rep PMID 38449301 |

> 🔑 **ISM 對其中 2 個有獨特解法**：難點 3（germline-somatic）→ **somatic HP1-vs-HP1-1 軸 + normal-anchored cis-test（量測 d_cis vs d_drift）**；難點 2（異質性）→ **read-level PERMANOVA 結構（不塌 mean）**。**其餘 6 個（覆蓋/CN/純度/單樣本/phasing/5hmC）ISM 與大家共卡** —— 尤其**難點 4/6 copy-number+purity，CAMDAC 有反卷積、ISM 只有 cis-test 部分控**。

---

## L3 — ISM 獨有 vs 欠缺（誠實總結）

### ISM 獨有（窄但真）
1. **somatic-controlled HP1-vs-HP1-1 軸**（非業界 germline HP1-vs-HP2 / REF-vs-ALT）—— 唯一結構性新穎。
2. **normal-anchored cis-test（d_cis vs d_drift）—— 量測而非只過濾** somatic-cis vs germline-drift/copy confound（操作化「異質性 vs 真 ASM」）。
3. **read-level 多變量 PERMANOVA 結構檢定**為 headline（多數對手停在 scalar）。
4. **誠實 hard NEGATIVE**（甲基非 variant filter）本身是貢獻。
5. 保守 gated 位點層（Cramér's V 被 Cochran gate 歸零，抑制稀疏假陽性）。

### ISM 欠缺（業界比我們強的）
1. **無原則性 spatial region 聚合 / region-FDR**（dmrseq AR(1) / comb-p / DMRcate / DSS smoothing 全遠超）。
2. **無更強 per-site 統計**（Fisher+Δβ 是 DSS/MOABS/methylKit/modkit 的子集且較弱；無 beta-binomial）。
3. **無 phasing-free 模式**（cvlr/ASMS/modbamtools/qFDRP/MeConcord 都不需 HP tag）。
4. **無 copy-number/purity deconvolution**（**CAMDAC** 反推純 tumor allele-specific 甲基；CancerLocator 聯合估 purity）。
5. **無 cohort / tissue-of-origin / expression-driver**（MethylMix/ELMER/MHB/MHL/CancerLocator）。
6. **無 scale 優勢**（O(N²)~O(N³) + read-pair subsample 不可重現問題）。

## Provenance
- 外部：workflow **wf_22ba1b9b-db4**（7 agents：4 軸 + cross-difficulties + canonical-pipeline + synth；WebSearch→WebFetch 一手；2026-06-11）。原始 `_assets/axis_reconfirm_raw.json` + digest。
- ISM：`01_ism_method_spec_from_source.md`（源碼 file:line）。
- 新工具浮現：**CAMDAC**（purity+CN deconvolution，bioRxiv 2020.11.03.366252）、**MethPhaser**（甲基-bridging phasing）—— 投稿前值得親讀。
- ⚠ ASMS（preprint）/ NANOME（preprint）/ CAMDAC 識別碼部分未完全 fetch；引用前再查。
