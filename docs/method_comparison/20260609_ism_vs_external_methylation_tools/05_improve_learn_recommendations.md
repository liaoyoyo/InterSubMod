<!--
建立時間: 2026-06-10
狀態: in_progress (method-comparison study, part 5/6 — 改進/學習建議)
報告類型: improvement_recommendations
受眾: PI · 開發決策 · 論文 Methods 強化
framework: 排序建議(ROI×effort) + KEEP/LEARN/FIX 三類
data_sources:
  - _assets/survey_digest.md (SYNTHESIS top_improvements + novelty_defensible)
  - 03_method_comparison_matrix.md, 04_cross_analysis_data_vs_literature.md
provenance_note: 每條建議對應一個 WebFetch-verified 外部方法 + 一個源碼-verified 的 ISM gap(file:line)。effort 為粗估。
-->
<!-- provenance-verified: 建議來源工具皆 WebFetch 實證; ISM gap 皆源碼確認(01); 撰寫與分析不同 batch。 -->

# 05 — 我們的方法可改進 / 可學習處（排序建議）

> **這份是什麼**：第 5 部分 —— 回答用戶「**確認我們的方法是否有可以改進或學習的**」。每條 = 一個外部 verified 方法 → 對應一個 ISM 源碼層 gap → 具體動作 + effort + ROI。分 **FIX（修正錯誤）/ LEARN（學習採納）/ KEEP（守住的新穎性）** 三類。

---

## L0 — 三個最高 ROI（先做這些）

1. **🔴 FIX：per-CpG Fisher → beta-binomial + dispersion shrinkage**（最高 ROI）。ISM 現用 Fisher exact 假設 read 獨立；但同單倍型/clone 的長讀**非獨立** → 低估變異 → 顯著 CpG 偏多。學 **DSS**。
2. **🟡 LEARN：read-read 距離建在 soft 機率上**（非硬 >0.8/<0.2 二值）。學 **cvlr / pycoMeth**，降低覆蓋噪音、傳遞 ML 不確定度。
3. **🔴 FIX：停止 5mC/5hmC max-collapse**（`cis_asm_core.py:31-32`）。改 **Dirichlet-Multinomial 分軌**（modkit 思路）；此 max-collapse 與已知 MSA dup-bug 同源，會稀釋訊號。

---

## L1 — 完整建議排序表

| # | 類 | 建議 | 來源工具 | ISM gap (源碼) | effort | ROI |
|---|----|------|---------|---------------|--------|-----|
| 1 | LEARN | per-CpG **beta-binomial + lognormal dispersion shrinkage**（給 Δβ model-based SE/CI）| DSS (Feng/Wu 2014) | `PerCpgAsm.cpp` Fisher 無 overdispersion | High | ★★★ |
| 2 | LEARN | read-read 距離用 **soft 機率**（保留 p(meth)）| cvlr / pycoMeth BernoulliPosterior | `DistanceMatrix.hpp:27-28` 硬二值 | Medium | ★★★ |
| 3 | FIX | **停 5mC/5hmC max-collapse** → Dirichlet-Multinomial 分軌 | modkit (Dirichlet-Mult); Du 2010 | `cis_asm_core.py:31-32` max-collapse | Low | ★★★ |
| 4 | LEARN | **Cohen's h（arcsine）** 作飽和-robust effect size，配 M-value 交叉檢核 | modkit dmr; Du 2010; Xie 2019 | Δβ 只有 raw + Cohen d | Low | ★★ |
| 5 | LEARN | **CpG-pair（tuple）共甲基 log-odds**（HP-stratified），補 per-CpG Fisher 漏掉的 on-read 連動（長讀可超 150bp） | DAMEfinder tuple | 無 tuple/co-methylation code | Medium | ★★ |
| 6 | LEARN | **modkit dmr 當快速 pre-filter** 篩候選 region，再跑昂貴 O(N³) UPGMA+999perm | modkit dmr; MHL screen | CORE 1-3 重 | Low-Med | ★★ |
| 7 | LEARN | **qFDRP 作 coverage-robust disorder 對照**，與 PERMANOVA 做 2×2（操作化「結構非失序」）| Scherer 2020 qFDRP | CORE4 用 coverage-fragile epipoly | Low | ★★ |
| 8 | LEARN | **imprinting/ICR 正控**（在已知 ASM 真值處驗證引擎偵得到）| pycoMeth/cvlr/NanoMethPhase germline benchmark | ISM 只有 HP-shuffle null，無 known-truth 正控 | Low-Med | ★★ |
| 9 | FIX | **修引用**：`PerCpgAsm.cpp:6` 「De Waele 2020」→ **Orjuela 2020**（DAMEfinder 正確作者；一手 PMC7268773 確認）| — | 源碼註解誤引 | Low | ★ |
| 10 | FIX | 論文 **Sakamoto 2022 → O'Neill/POG 2024**（private/recurrent hap-DMR 乾淨 anchor）| — | KB 引用不精確 | Low | ★ |

---

## L2 — 逐條細節（動作 + 為何 + 適配性）

### #1 LEARN — beta-binomial + dispersion shrinkage（最高 ROI）
- **動作**：per-CpG 把 `Y ~ beta-binomial(coverage, π_group, φ)`，φ 用 lognormal prior 跨 CpG shrinkage（DSS namesake）。輸出 Δβ 的 model-based SE/CI + Wald p。
- **為何**：ONT 長讀同 haplotype/clone reads 非獨立 → true over-dispersion φ>0 → Fisher（φ=0）低估變異、p 偏小 → 顯著位點膨脹。這是 ISM 在低 HP 覆蓋的核心不穩定來源。
- **適配**：完美適配單樣本（DSS 無 replicate 模式用鄰近 CpG 當 pseudo-replicate）。**注意**：O11 epipolymorphism NEGATIVE 曾是 n_reads confound → 任何新統計仍須過 coverage-校正關。
- 來源：Feng/Conneely/Wu 2014 NAR PMID 24561809；Wu 2015 NAR PMID 26184873。

### #2 LEARN — soft-probability 距離矩陣
- **動作**：`DistanceMatrix` 加 soft 選項：用 p(methylated)（從 ML/expected methylation）算距離，而非先二值化。
- **為何**：硬二值丟掉 ML 信心；低覆蓋（ISM 罰到 MAX_DIST=1.0 處）噪音大。cvlr Bernoulli-mixture + pycoMeth BernoulliPosterior 都保留 soft。
- **適配**：高（ISM 已存 raw 機率）。可作 NHD 的 soft 變體 + Bernoulli 距離天然吃機率。

### #3 FIX — 停 5mC/5hmC max-collapse
- **動作**：`collapse_modtype`（`cis_asm_core.py:31-32`）的 max-collapse 改成保留兩 mod-type；用 Dirichlet-Multinomial（h/m/C）或分軌平行算。
- **為何**：max-collapse 把 5mC+5hmC 兩列取最大 → 5hmC marginal 訊號被稀釋；與已知 MSA Level1 dup-bug 同源（memory 已記）。modkit 原生分 mod-code。
- **適配**：高，且修掉一個已知 artifact。**5mC/5hmC 分軌做 ASM 是文獻幾乎沒人做的真新穎軸**（見 KEEP）。

### #4 LEARN — Cohen's h + M-value 交叉
- **動作**：headline loci（BRCA2 −0.122 / TBC1D16 +0.142）加算 **Δ(M-value)** 確認 sign+顯著性在 variance-stabilization 後存活；report 仍用 Δβ（Du 建議），但 inference 用 M-value 版；effect size 加 **Cohen's h**（arcsine，飽和-robust）。
- **為何**：raw Δβ 在 0/1 端 heteroscedastic，而我們 headline 正在飽和區。
- 來源：Du 2010 PMID 21118553；Xie 2019 Bioinf；Cohen's h（modkit dmr 原生）。

### #5 LEARN — CpG-pair tuple 共甲基 log-odds
- **動作**：加 HP-stratified 的 tuple log-odds（MM·UU/MU·UM），長讀可在任意距離（不受短讀 150bp cap）算。
- **為何**：DAMEfinder 核心洞見 = on-read MM/UU vs MU/UM 連動帶 allele-specificity，single-CpG 測不到；ISM per-CpG Fisher 丟掉這資訊（只在 NME/epipoly 用 unsupervised 版）。ISM 有理想基質（完整 read×CpG + 長讀）。
- 來源：Orjuela 2020 Epigen Chrom PMID 32487212。

### #6 LEARN — modkit dmr 當 pre-filter
- **動作**：全基因組先跑 modkit dmr（快 Rust，Beta-Binomial LR）或 MHL coordination screen 篩候選 region，再對通過者跑昂貴的 UPGMA(O(N³))+PERMANOVA(999perm)。
- **為何**：降 compute 不改結論層。modkit dmr 重疊 ISM 的 β-rate/Δβ 層，正好當快速 triage。
- **適配**：高（modkit 單一 Rust binary 易裝）— 也直接接到 06 benchmark。

### #7 LEARN — qFDRP 作 coverage-robust disorder 對照
- **動作**：CORE4 對照槽用 **qFDRP**（Scherer 推薦、最 coverage-robust，且用**和 ISM NHD 相同的 normalized-Hamming kernel**）取代 coverage-fragile 的 epipolymorphism；做 2×2：高 qFDRP + 非顯著 PERMANOVA = stochastic erosion；高 qFDRP + 顯著 PERMANOVA = 結構化 allelic ASM。
- **為何**：直接操作化「結構非失序」論點，且修掉 O11 暴露的 epipoly coverage artifact。
- 來源：Scherer 2020 NAR PMID 32338758。

### #8 LEARN — imprinting/ICR 正控
- **動作**：在 canonical imprinted ICR（已知 ASM 真值）跑 ISM，證明引擎在真有 ASM 時偵得到（cvlr/pycoMeth/NanoMethPhase 都這樣 validate）。
- **為何**：ISM 目前只有 HP-shuffle null，**無 known-truth 正控** —— reviewer 會問。便宜、高 credibility。

### #9 / #10 FIX — 引用修正
- **#9**（2026-06-10 一手查證後校正）：**只有 `src/core/PerCpgAsm.cpp:6`** 註解寫「DAMEfinder (**De Waele 2020**)」是錯的；`PerCpgAsm.hpp:8` 寫「(DAMEfinder, pycoMeth)」**無誤**（survey agent 誤報 .hpp:8 也錯，已 grep 推翻 — §13.7 實例）。正確作者 = **Orjuela S, Machlab D, Menigatti M, Marra G, Robinson MD. Epigenetics & Chromatin 2020;13:25. DOI 10.1186/s13072-020-00346-8**（一手 PMC7268773 確認）。
  - **根源**：`InterSubMod/docs/references/20260415_ASM_subclone_methods_literature_survey.md:44,519` 把作者捏造成「De Waele L, Fourne L, Lent J, et al.」（同 DOI）—— C++ 註解由此沿用。已一併訂正該 doc。
  - （亦：task prompt 一度把 DAMEfinder 期刊寫成 NAR Genom Bioinform，也錯 —— 是 Epigenetics & Chromatin。）
- **#10**：KB/論文若引「Sakamoto 2022 = patient-private hap-DMR」→ 不精確（Sakamoto 是肺癌 phasing 論文）；private/recurrent hap-DMR 的乾淨 cancer anchor 改 **O'Neill/POG 2024 (Cell Genomics, PMC11605692)**。

---

## L3 — KEEP：守住的可防守新穎性（不要被「對手都做過」誤導）

> synth critic 的誠實結論：ISM 單一步驟都非新，但**這個組合 + 兩個能力是無人佔的**。

1. **PERMANOVA-on-read-read-distance 的「結構非失序」supervised 顯著性檢定** —— 整個 disorder 家族（PDR/epipoly/entropy/FDRP/qFDRP/NME，ISM CORE4 已內建）回答「分子內多隨機」；ISM CORE 1-3 回答「分子是否分成結構分離的、label-tested 的亞群」。**cvlr（最近鄰）有分群但無顯著性檢定** → 這是最該 foreground 的單一新穎點。
2. **normal-anchored somatic cis-control（d_cis vs d_drift）** —— **唯一**：所有對手（DAMEfinder/NanoMethPhase/pycoMeth/AMRfinder/modkit dmr）都是對稱 A-vs-B 或 trio，**無法分 somatic-specific allelic 甲基 vs germline baseline**。
3. **5mC/5hmC 分軌做 ASM** —— 文獻幾乎沒人做（前提是先修 #3 max-collapse）。

> ⚠ 但 synth 也誠實指出：**「組合無人佔」≠「該組合更好/有用」**。沒有 worked discordant case 證明「modkit 說無率差、ISM PERMANOVA 說有結構」時 ISM 的 call 是真訊號還是噪音。→ **這正是 06 Phase B benchmark 要回答的**。

→ 下一份 `06`：把 #1-#8 與「驗證新穎性是否 matter」轉成**可下載執行的 Phase B benchmark 方案**。
