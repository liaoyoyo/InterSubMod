<!--
建立時間: 2026-06-10
狀態: in_progress (method-comparison study, part 3/6 — 方法對照矩陣)
報告類型: method_comparison_matrix
受眾: PI · 論文 Related Work/Methods 對照 · reviewer 防守
framework: 對照矩陣(capability × tool) + 逐工具 point-by-point
data_sources:
  - 01_ism_method_spec_from_source.md (ISM 側, 源碼確認)
  - _assets/survey_digest.md (deep-read vs_ism_detailed, WebFetch 實證)
provenance_note: ISM 側 = 源碼 file:line; 外部側 = WebFetch 一手來源。✓/✗/部分 為能力是否存在的事實判斷, 非優劣評分。
-->
<!-- provenance-verified: ISM 能力引自 01 源碼確認; 外部能力引自 workflow deep-read WebFetch。 -->

# 03 — 方法對照矩陣（ISM 6 核心 × 外部工具，逐格比）

> **這份是什麼**：第 3 部分 —— 把 ISM 的能力**逐項**對到最相關的外部工具，標「有/無/部分」+ **細節差異**。回答用戶「**是否有細節不同**」。
> ✓=有此能力 ／ ⚠=部分或需外部前置 ／ ✗=無。

---

## L0 — 一句話結論

**沒有任何單一外部工具同時具備 ISM 的「read-read 距離矩陣 + PERMANOVA 結構顯著性 + normal-anchored somatic cis-test + 5mC/5hmC 分軌 + LOH/CN 耦合」五項組合**；但**每一項單獨拆開都已有更成熟的對手**（cvlr 做 read 分群、DSS 做嚴謹率差統計、MHB 做 within-read LD、qFDRP 用同一 Hamming kernel、CpelNano 做 ONT 失序模型）。ISM 的可防守新穎性 = **這個特定組合 + 「結構非失序」的 supervised PERMANOVA + somatic cis 控制**，而非任何單一步驟。

---

## L1 — 能力對照矩陣（ISM 10 項 × 11 個最相關工具）

| ISM 能力 ＼ 工具 | **ISM** | modkit dmr | cvlr | ASMS | pycoMeth | NanoMethPhase | DAMEfinder | DSS | MHB/MHL | Metheor | qFDRP/WSH |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **平台 = ONT 長讀** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗短讀 | ✗短讀 | ✗短讀 | ✗短讀 | ✗短讀 |
| **read×CpG 矩陣保留**（不塌成 per-site %）| ✓ | ✗聚合 | ✓ | ✓ | ⚠塌 scalar | ⚠塌 freq | ⚠tuple | ✗計數 | ⚠per-read | ✓read-sweep | ⚠read-pair |
| **C1 read-read 距離矩陣（顯式 N×N）** | ✓ | ✗ | ✗(模型式) | ✗(模型式) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ⚠同 kernel,塌 scalar |
| **C2 階層 clustering 分亞群** | ✓UPGMA+silhouette | ✗ | ✓EM(固定k) | ✓EM(2-comp) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| **C3 結構顯著性檢定（PERMANOVA）** | ✓999perm | ✗ | ✗無檢定 | ⚠1v2 LR | ⚠MW/KW(per-interval) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| **per-CpG 率差 / Δβ** | ✓Fisher+Wilcoxon | ✓Bayes LR+MAP | ✗ | ✗ | ✓MW/KW | ✓(DSS) | ✓ASMsnp | ✓✓Wald+shrinkage | ✗ | ✗ | ✗ |
| **disorder 度量（NME/epipoly）** | ✓(CORE4 對照) | ✓entropy | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ⚠MHL | ✓✓全套 | ✓✓全套 |
| **within-read CpG-CpG LD/tuple** | ✗(可學) | ✗ | ✗ | ✗ | ✗ | ✗ | ✓✓tuple log-odds | ✗ | ✓✓LD r² | ✓LPMD | ✗ |
| **5mC / 5hmC 分軌** | ⚠有資料但曾 max-collapse | ⚠adjust-mods 可分 | ✗ | ✗ | ✗ | ✗5mC only | ✗ | ✗ | ✗ | ✗ | ✗ |
| **somatic-controlled（HP1 vs HP1-1）** | ✓✓ | ✗對稱 A/B | ✗ | ✗ | ✗ | ✗germline HP1/HP2 | ✗ | ✗ | ✗ | ✗ | ✗ |
| **normal-anchored cis-test（d_cis vs d_drift）** | ✓✓**唯一** | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| **LOH/CN 耦合** | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| **顯著性 + effect-size 雙報** | ✓Cohen d | ✓Cohen h | ✗ | ✗ | ⚠ | ✓(DSS) | ✓ | ✓ | ✗ | ✗ | ✗ |
| **變異穩定化 / dispersion shrinkage** | ✗(可學) | ✓Jeffreys+Cohen h | ✗ | ✗ | ✗ | ✓✓ | ✓limma eBayes | ✓✓lognormal shrink | ✗ | ✗ | ✗ |

> 讀法：✓✓ = 該工具此能力**比 ISM 更成熟**（可學習對象）；✓ = 有；⚠ = 部分/需前置/塌成標量；✗ = 無。

---

## L2 — 逐工具 point-by-point（細節差異，最相關 8 個）

### 1. modkit dmr（ONT 官方，最強 production 長讀基準）
- **同**：都讀 ONT MM/ML、都做 5mC/5hmC、binarize 概念近（modkit pileup 閾值 vs ISM >0.8/<0.2）。
- **異（關鍵）**：
  - modkit dmr 在 **per-position 聚合計數**上比兩組 → 量「率差(β)」；ISM 保留 read×CpG 量「結構」。modkit dmr ≈ ISM 的 **Δβ/CORE4** 層做成 Bayesian，**完全沒有** ISM 的 CORE 1-3（距離/clustering/PERMANOVA）。
  - **modkit dmr 無內部 `--haplotype`** → 單倍型 DMR 須先 `pileup --phased` 拆 hp1/hp2 當「兩個 sample」餵入 → **失去分子連結、把 haplotype 當靜態樣本**；ISM 的 HP1 vs HP1-1 是 somatic-controlled。
  - modkit **無** normal-anchored cis-test、無 LOH/CN 耦合、無 read clustering。
- **modkit 比 ISM 強之處**：Beta-Binomial marginal-LR + Jeffreys prior（小計數穩健）、**Cohen's h**（飽和覆蓋 robust effect size）、HMM 分段、速度（Rust）。→ **ISM 可學**（見 05）。

### 2. cvlr（軸 C 最近鄰）
- **同**：ONT、read-level、**不需 phasing**、把 reads 依甲基 pattern 分群。
- **異**：cvlr = **模型式**（multivariate Bernoulli mixture + EM，soft posterior）；ISM = **距離式**（顯式距離矩陣 + UPGMA）。cvlr **無 read-read 距離矩陣、無顯著性檢定、固定 k 無 BIC、germline-demonstrated（無 cancer/somatic/normal-anchor）**。
- **cvlr 比 ISM 強之處**：soft 機率成員（surface 邊界模糊 read，正是 cancer subclone 誤判處）、缺值 -1 的 EM imputation。→ ISM 可學 soft clustering + BIC 選 k。

### 3. ASMS（Raineri 2024，單一最像 ISM，⚠UNVERIFIED）
- ONT modBAM + 不需 phasing + EM 2-component 依甲基分 reads 找 ASM。**最該正面對照**。ISM 多了：顯式多-metric 距離矩陣、PERMANOVA 顯著性、normal-anchored somatic cis、>2 群 silhouette。**投稿前必親讀確認**。

### 4. pycoMeth（haplotype-aware 但 region-level）
- **同**：ONT、read-level 儲存(MetH5)、haplotype-aware（WhatsHap read-group）、Bayesian HMM 分段。
- **異**：差異檢定 **塌到 region/interval**（llr_diff 池化 LLR 做 MW/KW；bs_diff per-read β scalar）→ **無 read-read 距離矩陣、無 read clustering、germline-centric、無 somatic normal-anchor、無 5mC/5hmC 分軌**。
- **pycoMeth 比 ISM 強之處**：**soft-call emission**（保留 p(methylated) 而非硬二值）、joint multi-group consensus 分段。→ ISM 可學 soft 機率建距離矩陣。

### 5. NanoMethPhase（HP 軸最直接對照，long-read）
- **同**：ONT、haplotag、比 HP1 vs HP2 甲基。
- **異**：需 **germline phased SNV**、測 **HP1 vs HP2（germline imprinting）**；差異檢定 **外包給 DSS**（beta-binomial）；**5mC only**。ISM 加 **somatic HP1 vs HP1-1 + normal-anchored cis + LOH（germline-het 消失處）+ 5hmC**。
- 這是「別人用 ONT haplotype 甲基做什麼」的標準引用。

### 6. DAMEfinder（ISM 源碼明文引用的對象）
- ISM `PerCpgAsm.cpp` 把 per-CpG Fisher 標註出處為 DAMEfinder。
- **異**：DAMEfinder 招牌是 **tuple（CpG-pair）log-odds 共甲基連動**（ISM 的 per-CpG Fisher 反而**丟掉** cross-CpG pattern，那部分 ISM 是放在 NME/epipoly）；短讀 BS、germline/imprinting、150bp tuple cap、SNP 模式 re-derive allele（不如 ISM 既有 phased 長讀 HP tag）。
- 🔴 **發現的引用錯誤**（已查證+修正）：ISM 源碼 **`PerCpgAsm.cpp:6`** 標 **"De Waele 2020"** 是錯（`.hpp:8` 寫 "(DAMEfinder, pycoMeth)" 無誤）；**正確 = Orjuela et al. 2020 (Epigenetics & Chromatin 13:25)**（一手 PMC7268773 確認）。根源 = `docs/references/20260415_..._survey.md` 捏造作者「De Waele L, Fourne L, Lent J」（見 05 #9）。

### 7. DSS（短讀黃金標準 DMR 統計，方法藍圖）
- **同**：都要算「兩組 per-CpG 甲基差」。
- **異（最大方法 gap）**：ISM per-CpG 用 **Fisher exact**（假設每 read 獨立 Bernoulli，φ=0）；DSS 用 **Beta-Binomial + lognormal dispersion shrinkage + Wald**，把「同單倍型/clone 的 reads 非獨立」的 over-dispersion 顯式吸收。**ONT 長讀 somatic 中 reads 高度非獨立 → ISM Fisher 低估變異 → p 偏小 → 顯著 CpG 偏多**。
- → **ISM 最高 ROI 的可學項**：per-CpG 換 beta-binomial + shrinkage，給 Δβ 一個 model-based SE/CI（見 05）。

### 8. MHB/MHL（Guo 2017）— 軸的關鍵釐清
- **載重式釐清**：MHB/MHL 在 **within-read、沿基因組**軸（CpGᵢ vs CpGⱼ 在同分子上的 LD r²，reads 只是估 pair 連動的觀測，**reads 彼此不互比**）；ISM 在 **between-read**軸（read_a vs read_b 的距離，CpG 軸被 collapse）。**兩軸正交** —— 這是回答用戶「甲基位點間關聯方法 vs 我們」的核心：MHB 是「位點×位點」，ISM 是「分子×分子」。
- → ISM 可借 MHB 的 LD-r² greedy block 來**資料驅動定 analysis window**，再在 block 內跑 read-read。

### 補：Metheor LPMD / qFDRP — 與 ISM 距離的精確區別
- **Metheor LPMD** = within-read、固定距離 CpG-pair 不一致（**位置×位置，同分子**）；ISM read-read = between-molecule。**必在論文明確區分，否則 reviewer 易混**。
- **qFDRP** = 用**和 ISM NHD 完全相同**的 normalized-Hamming read-pair kernel，但**平均成 per-CpG scalar 丟掉幾何、無 label**。→ ISM 可把 qFDRP 當 coverage-robust 的「失序對照」與 PERMANOVA「結構」做 2×2（見 04/05）。

---

## L3 — 來源
- ISM 能力：`01_ism_method_spec_from_source.md`（源碼 file:line）
- 外部能力：`_assets/survey_digest.md` 的 DEEP-READ 段（每筆 WebFetch verified_claims + 一手 URL）
- 完整：`_assets/workflow_raw_result.json`

→ 下一份 `04` 把這些方法差異**接到實際數據**：我們的 Δβ/ASM 結果 vs 外部論文觀察與結論。
