<!--
建立時間: 2026-06-10
狀態: in_progress (method-comparison study, part 4/6 — 結果×文獻交叉分析)
報告類型: cross_analysis_internal_vs_external
受眾: PI · 論文 Discussion · reviewer 防守
framework: 5 維對照(external↔ours↔verdict↔improve) + alignment 分級 (CONFIRM/CALIBER/DIFFERENT-VIEW)
data_sources:
  - _assets/survey_digest.md (CROSS 段, 每維 WebFetch 實證來源)
  - knowledge/11_external_literature/02,05,06,07 (內部 validated 數字)
provenance_note: 「我們的數字」皆引自既有 validated KB(已標來源,非本檔新跑); 外部數字皆 WebFetch 一手實證附 PMID/DOI。alignment verdict 由 cross-analysis agent 產出, 已對齊既有 KB doc 結論。
-->
<!-- provenance-verified: 內部數字引 validated KB(11_external_literature/02,05,06,07); 外部數字引 workflow CROSS 段 WebFetch 一手來源; 撰寫與分析不同 batch(§13.0)。 -->

# 04 — 交叉分析：我們的結果 × 外部論文數據（觀察與方法結論）

> **這份是什麼**：第 4 部分 —— 把**我們的實際研究結果（Δβ、ASM 量值/方向、filter-NEGATIVE、phasing）**逐一對到**外部論文的數據觀察與方法結論**，判定 **一致 / 口徑差 / 不同觀點**，並指出可改進處。回答用戶「**比對外部論文數據的觀察與方法結論**」。
> 每維 alignment：**CONFIRM**（一致）／**CALIBER**（口徑差非矛盾）／**DIFFERENT-VIEW**（不同方法觀點）。

---

## L0 — 五維一句話

| 維度 | alignment | 一句話 |
|------|-----------|--------|
| **1. Δβ 方法學** | 🟡 PARTIAL / SOUND-FOR-PURPOSE | ISM 的 raw 均-Δβ + Wilcoxon 描述層對齊文獻，但**缺 variance-stabilization** —— magnitude claim 需 M-value/beta-binomial 交叉檢核 |
| **2. ASM 量值/方向** | 🟡 CALIBER（非矛盾）| 「稀有 0.34%、無方向、BRCA2 hypo」與 Do2020「+5×、hypo-dominant」是**不同分母**（within-sample somatic-controlled read-率 vs cross-individual SNP-anchored DMR-率），非衝突 |
| **3. structure vs disorder** | 🟢 STRONGLY ALIGNED | **Landau 2014 明文排除 ASM 是 disorder 來源** = 我們「結構非失序」的最強外部背書；O11 崩潰 = 公認 WSH coverage artifact |
| **4. locus-association 軸** | 🟢 COMPLEMENTARY + 真空缺 | within-read LD(MHB) / between-read 距離(ISM) / cohort 網路 三軸正交；ISM 的組合無人佔 |
| **5. filter-NEG + phasing 白地** | 🟢 BOTH SUPPORTED | 無 caller 用甲基當 standalone FP filter；LOSO circularity 是公認 ML 陷阱；germline methyl-phaser 明確留 cancer/LOH 為 future |

---

## L1 — 逐維交叉（external ↔ ours ↔ verdict）

### 維 1：Δβ / 甲基差 effect-size 方法學 🟡 PARTIAL
- **外部（一手實證）**：領域已收斂到**變異穩定、dispersion-aware、模型式**差異估計，而非 raw per-CpG 均差。
  - DSS：beta-binomial + lognormal dispersion **shrinkage** + Wald（Feng/Wu 2014, PMID 24561809）
  - Du 2010（PMID 21118553）：raw beta 在 [0,0.2]∪[0.8,1] **嚴重 heteroscedastic** → 建議 M-value 檢定、report beta
  - Xie 2019：檢定尺度與 report 尺度須一致；naive beta 在極端偏
  - modkit dmr：Cohen's h（arcsine 穩定）作飽和-robust effect size
- **我們**：Δβ = per-CpG 均甲基差（HP1 vs HP1-1 / ALT vs REF），按 CpG 配對 → **Wilcoxon signed-rank + Cohen's d**，raw 機率二值化(>0.8/<0.2)，per-CpG Fisher+BH-FDR；5mC/5hmC **max-collapse**。
- **VERDICT = SOUND-FOR-PURPOSE，但 magnitude 需交叉檢核**：核心架構（read-read 距離 + PERMANOVA + normal-anchor）健全且新穎；但 effect-size **magnitude** 有已知有界偏差：(a) raw Δβ 在 0/1 端 heteroscedastic（BRCA2=−0.122、TBC1D16=+0.142 都在飽和區附近）；(b) Fisher 假設 read 獨立，低估變異；(c) 5mC+5hmC max-collapse 會稀釋（已知 MSA dup-bug 同源）。Wilcoxon（rank-based）部分繞過 (a)，是合法選擇。
- → **改進見 05 #1/#2/#5**（M-value 交叉、beta-binomial、Cohen's h、停 max-collapse）。

### 維 2：ASM 量值與方向 🟡 CALIBER（口徑差，非矛盾）
- **外部（一手實證）**：
  - **Do 2020**（PMC7322865）：15,112 recurrent ASM DMR = 0.7% informative regions；cancer ASM **5–9×** > normal；**hypo-dominant**。← 但這是 **cross-individual、SNP-anchored、多樣本 DMR-率**。
  - **Martin-Trujillo 2017**（PMC5589900）：CN-normal 樣本中 imprinted allelic 甲基異常 **82–92% 由 CN 解釋**。← 直接背書我們 HP-axis held-CN 設計擊敗的 confound。
  - **O'Neill/POG 2024**（PMC11605692）：4.46M aDMR，mean **23.6K/tumor**（~5× normal），**79% 落 CNV/LOH**，**recurrent** RET(26)/CDKN2A(9)。
  - Schalkwyk 2010：germline >35,000 cis-ASM（baseline floor）。
- **我們**（HCC1395 單樣本 ⭐3 單 pipeline）：strong-ASM **0.34%**、**無方向**(hypo44/hyper56)；BRCA2 HP-axis Δβ=−0.122（~80% copy + 邊際真 cis d_within=−0.023 p=0.022）；ASM×CN ρ=−0.055（**非 CN-driven**）；cross-sample 6/6 excess；somatic private **0/38**。
- **VERDICT = DIFFERENT CALIBER（主），一處低-power 方向張力**：
  - 我們 0.34% 是 **stringent within-sample somatic-controlled HP-axis read-率**；文獻 5–9× / 4.46M 是 **cross-individual / haplotype-anchored DMR 計數**（含 germline cis-ASM + CN/LOH-driven allelic 甲基）→ **不同分母，非矛盾**。
  - 「無方向」是 somatic-controlled 軸的**預期簽名**（germline allelic bias 被扣除後 hypo/hyper 應趨平衡）；Do2020 的 hypo-dominance 在 cross-individual 軸（含未扣的 germline）。
  - **Martin-Trujillo 是最強外部對齊** —— 他們警告的 CN-主導正是我們 HP-axis 設計 neutralize 的、ρ=−0.055 確認非 CN-driven。
  - **誠實限制**：recurrence 是決定性 gap —— POG 在 189 cohort 偵 recurrent driver aDMR，我們單 pipeline 6 樣本偵不到（private 0/38 = underpowered 非 biological privacy）。
- 🔴 **引用修正**：digest 指出「Sakamoto 2022 = patient-private hap-DMR」不精確（Sakamoto 是肺癌 phasing 論文）；乾淨 anchor 應改 **O'Neill/POG 2024**。

### 維 3：structure-not-disorder 框架 🟢 STRONGLY ALIGNED
- **外部（一手實證）**：
  - **Landau 2014 PDR**（PMC4302418）：**明文測試並排除 ASM 是 discordant reads 的來源**（「discordant reads 涉及兩個 parental allele」、>2 patterns/locus 排除 1-2 ASM patterns），disorder = 「stochastic genome-wide process」。← **單一最強背書**：disorder 工具構造上量不到 haplotype-linked 結構。
  - **Scherer 2020 WSH**：epipoly/entropy/FDRP/qFDRP 皆 **inter-molecule**；epipoly/entropy 需 >10×、coverage-fragile；qFDRP 最 robust。
- **我們**：read×CpG → read-read 距離 → **PERMANOVA**（檢定 **pre-labeled** HP/somatic 群是否結構分離）；O11 epipolymorphism AUC **0.845→0.530**（n_reads 校正後崩）。
- **VERDICT = SUPPORTED，一處框架修正**：ISM PERMANOVA 確實量「**supervised、label-tested 的 between-molecule 結構**」，與 disorder 文獻的「unsupervised、explicitly-stochastic、non-allelic 失序」不同量。O11 崩潰 = **Scherer 證實的 WSH coverage artifact**（真方法學 win）。
  - 🔴 **框架修正**：**不要**講「within-read(disorder) vs between-read(ISM)」—— Scherer 證 epipoly/FDRP/qFDRP 都是 inter-molecule。正確框架 = 「**label-tested 結構 vs label-free 熵**」。
- → **改進見 05 #?**（加 qFDRP 作 coverage-robust disorder 對照，與 PERMANOVA 做 2×2）。

### 維 4：methylation-locus-association 軸 🟢 COMPLEMENTARY + 真空缺
- **外部（一手實證）**：三軸正交，**無一佔 ISM 的格**：
  - (a) within-read CpG-CpG LD：MHB/MHL（Guo 2017，147,888 blocks）、Metheor LPMD、DAMEfinder tuple
  - (b) between-read read-read 距離：cvlr（EM、無顯著性）、ISM
  - (c) cohort CpG-CpG 相關網路：CoMeBack、Comethyl（WGCNA）
- **我們**：站 (b)，但配置無外部匹配 = 多-metric 距離矩陣 + **PERMANOVA 顯著性** + normal-anchored cis + 5mC/5hmC + LOH/CN。
- **VERDICT = DEFENSIBLE NICHE，窄但真**：單獨元件都非新（read 分群=cvlr；within-read LD=MHB；haplotype DMR=NanoMethPhase；disorder=Metheor）；**真正未填的組合** = read-read 距離(距離式非模型式) + PERMANOVA 結構顯著性 + normal-anchored somatic cis + 5mC/5hmC + cancer/LOH。
- → 論文須 **foreground PERMANOVA-on-read-read-distance 顯著性檢定**（cvlr 缺的）作最可防守單一新穎點。

### 維 5：filter-NEGATIVE + phasing 白地 🟢 BOTH SUPPORTED
- **外部（一手實證）**：
  - **無任何已發表 somatic caller 用甲基當 standalone FP filter**：ClairS-TO（9 hard-filter + 4 PoN + verdict）、ClairS+LongPhase-S（haplotype 一致性）、DeepSomatic（pileup-CNN）—— 全不用甲基。
  - **ML/sample-level circularity 是公認陷阱**：Soneson 2014（PMID 24967636）、Bernett 2024 Nat Methods、Kapoor-Narayanan leakage survey（feature-selection leakage 平均膨脹 AUROC ~0.18）。
  - germline methyl-phaser（MethPhaser/NanoMethPhase/HapBridge）**明確留 cancer/LOH 為 future**；cancer LOH-phaser（Wakhan/HiCancer）用 CNA/SV/Hi-C **不用甲基**。
- **我們**：filter NEGATIVE 四道（LOSO 100% circularity / AUC=0.505 / strong-ASM anti-discriminative OR=0.194 / 5th-rank vestigial）；LOH-constrained phasing NG=2 inner 93–99% 6/7、n=7 W=28 **p=0.0078**（Grade B+）。
- **VERDICT = 兩者皆被一手文獻支持**：filter-NEGATIVE 是「未被嘗試之想法上的乾淨負結果」（非與既有矛盾）；phasing 白地是「germline 工具明確止步、cancer LOH-phaser 不用甲基」的精確未佔格。
- → **改進見 05**：把 LOH-phasing 正面定位成 gap-filler（不是新方法類）。

---

## L2 — 「見樹也見林」：我們結果在外部地景的四層定位

| 層 | 內容 |
|----|------|
| **林（aggregate）** | 我們的全部 ASM 結論 = 單樣本 HCC1395 ⭐3 單 pipeline；外部 cohort（Do2020 n=16 癌、POG n=189）在更高 caliber 確認**現象真實**但我們**無 recurrence power** |
| **樹1（canonical）** | BRCA2 promoter：我們 −0.122（~80% copy）；POG 報 BRCA1/RAD51C/MLH1 是 **promoter hyper**(HRD/Lynch) 非 aDMR → 我們 BRCA2 hypo **勿當 TSG-silencing 證據** |
| **樹2（extreme/outlier）** | chr17/TBC1D16 = 我們唯一乾淨 cis exemplar（d_within=0.142 p=0.001）；對應 Do2020/Schalkwyk 的 cis-driven by sequence |
| **樹3（well-explained）** | structure-vs-disorder：Landau 2014 **明文排除 ASM 是 disorder 來源** → 機制示範我們 PERMANOVA 量的是不同東西 |

---

## L3 — 來源
- 內部數字：`knowledge/11_external_literature/02`(O11/0.34%)、`05`(filter-NEG)、`06`(phasing)、`07`(ASM/BRCA2/CN) — 皆 validated。
- 外部數字：`_assets/survey_digest.md` CROSS 段（每維 sources 附 PMID/DOI/PMC + WebFetch verified）。
- ⚠ 須補：Sakamoto 2022 改引 O'Neill/POG 2024；ASMS（Raineri 2024）親讀。

→ 下一份 `05` 把這些差異轉成**可改進/可學習的具體建議（排序）**。
