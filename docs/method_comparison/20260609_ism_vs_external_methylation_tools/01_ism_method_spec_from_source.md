<!--
建立時間: 2026-06-09
狀態: in_progress (method-comparison study, part 1/6 — ISM-side authoritative basis)
報告類型: method_specification_from_source
受眾: 比較分析自身基準 · PI · 外部協作者 · 未來的自己
framework: Diátaxis(reference) + 四層分層揭露(L0→L3)
data_sources:
  - src/core/*.cpp, include/core/*.hpp (C++ 方法本體, read 2026-06-09)
  - research/tsg_promoter_asm_reviewer/scripts/03,34,37*.py + pipeline/lib/cis_asm_core.py (Δβ/cis-test/copy-partition, read 2026-06-09)
  - knowledge/11_external_literature/03_ism_method_comparison.md (validated 內部定位)
  - knowledge/05_tools/intersubmod.md, methyl-somatic-analysis.md (KB tool docs)
provenance_note: 方法描述(演算法/參數/閾值) = 直接讀 C++/Python 源碼, 附 file:line; 引用的「validated 結果數字」= 取自既有 validated KB 文件(已標來源), 非本文件新跑分析。本檔 Write 與任何分析不同 batch(§13.0)。
-->
<!-- provenance-verified: 方法參數(NHD/C_min=5/UPGMA/999perm/THR=0.5/MIN_N=3) 直接引自 src/core 與 tsg scripts 源碼行號; 結果數字(Δβ=-0.122 等)引自 knowledge/11_external_literature/07 validated 文件。撰寫與源碼閱讀分離。 -->

# 01 — ISM (InterSubMod) 方法基準：直接從源碼確認

> **這份是什麼**：方法比較研究的**第 1 部分** —— 把「我們自己的方法（ISM）到底做了什麼」精確寫到**演算法 / 統計 / 閾值 / file:line** 層級，作為後續對照外部工具（modkit、pycoMeth、DSS、DAMEfinder、MHB…）的**權威基準**。所有方法描述皆來自 2026-06-09 直接讀源碼；引用的結果數字皆標 validated 來源。
> **導航**：本研究 6 文件總入口 → `InterSubMod/docs/method_comparison/20260609_ism_vs_external_methylation_tools/00_INDEX.md`

---

## L0 — 一句話定位（給最忙的人）

**ISM = 在 ONT 長讀癌症資料上，以「單 read × CpG 甲基矩陣 → read-to-read 距離矩陣 → 階層 subclonal clustering → 對 haplotype/somatic 標籤做 PERMANOVA 結構檢定 + normal-anchored cis-test」量化「haplotype-linked 甲基結構（between-molecule structure）」的 C++17 引擎。** 它測的是**結構/連結**（誰跟誰同一甲基亞群），不是 within-read 的**失序**（epipolymorphism/PDR）；它鎖定**癌症 somatic（tumor-normal paired）**，不是 germline imprinting。

---

## L1 — 6 個方法核心（一表看完，附源碼位置）

| # | 核心 | 做什麼 | 關鍵參數（源碼預設）| 源碼 |
|---|------|--------|------|------|
| **1** | **read-read 距離矩陣** | 對 region 內每對 read 算甲基模式距離 | metric=**NHD**(預設)；另 L1/Bernoulli/Pearson/Jaccard；`min_common_coverage=5`(C_min)；NaN→`MAX_DIST=1.0`；binary 閾值 高>0.8/低<0.2；可 strand-aware | `include/core/DistanceMatrix.hpp:21-36,44-57` |
| **2** | **階層 subclonal clustering** | 距離矩陣 → 樹 → 切成 read 亞群 | **UPGMA**(預設，分子鐘假設)；另 Ward/Single/Complete；切樹 by distance / k / **silhouette 自動選 k=2..10** | `include/core/HierarchicalClustering.hpp:23-37,176-207` |
| **3** | **PERMANOVA 結構檢定**（"結構"引擎）| 檢定 HP/somatic 群在距離空間是否**結構分離** | pseudo-F=(SS_between/(k−1))/(SS_within/(N−k))；**999 permutations**；min reads=5；另 dispersion 檢定(ANOVA-F on dist-to-centroid) | `include/core/StructureTest.hpp:24-36`；`src/core/StructureTest.cpp:141-205` |
| **4** | **per-CpG ASM + disorder 度量** | 同時算 3 個「文獻認證」度量供對照 | (a) per-CpG **Fisher exact + BH-FDR**；(b) **NME**(normalized methylation entropy)；(c) **epipolymorphism**(4-CpG 滑窗) | `include/core/PerCpgAsm.hpp:3-13,36-98` |
| **5** | **Cramér's V + 可靠性閘** | 全域 HP×甲基關聯強度 | Cramér's V；**Cochran 最小期望格≥5 reliability flag**；`min_cramers_v` 閾值 | `src/core/GlobalTest.cpp:128-141` |
| **6** | **NormalBaseline + normal-anchored cis-test** | 扣除 germline 基線、分辨真 cis-driven ASM vs drift | residual=raw−per-CpG normal mean；cis 三角 d_somatic/d_cis/d_drift（見 L2-F）| `include/core/NormalBaseline.hpp:38-67`；tsg `scripts/34_*.py:108-127` |

> 🔑 **程式碼本身就 cite 了競品**：`PerCpgAsm.hpp:7-13` 註解明寫三度量家族出處 —— Fisher per-CpG「(DAMEfinder, pycoMeth)」、NME「(CPEL, Jenkinson 2020)」、epipolymorphism「(methclone 2014, Metheor 2023)」。也就是說 ISM 的 disorder 度量是**刻意實作既有文獻度量供 contrast**，真正的差異化在核心 1/2/3/6（結構軸）。

---

## L2 — 逐核心細節（含公式與閾值）

### A. 資料表示 — MethylationMatrix（reads × CpG）
- 從 BAM `MM/ML` tags 取 5mC（以及 5hmC）modification 機率，建 region-level 矩陣（行=read、列=CpG，cell=raw mod 機率）。
- 二值化閾值：`binary_threshold_high=0.8`（甲基化）/ `binary_threshold_low=0.2`（未甲基化），介於兩者間為 ambiguous（`DistanceMatrix.hpp:27-28`）。
- HP 標籤分群（`LabelTest.cpp:244-305`）：
  - **HP1-family** = `1` + `1-1`（HP1 germline + HP1 上 somatic 重建）→ group 0
  - **HP2-family** = `2` + `2-1` → group 1
  - **fine 4-group**：HP1(0)/HP1-1(1)/HP2(2)/HP2-1(3)；HP3/HP0/unphased 排除(−1)

### B. read-read 距離（CORE 1）
- 預設 **NHD**（normalized Hamming distance on binary methylation patterns）；可選 L1、Bernoulli、Pearson（mean-centered）、Jaccard。
- **每對 read 需共同覆蓋 ≥5 個 CpG（C_min=5）** 才算有效；不足則距離記 `MAX_DIST=1.0`（懲罰式，而非丟棄）。
- 可分 forward/reverse strand 各自算（`DistanceCalculator::compute_strand_specific`）。
- **這是 ISM 最 distinctive 的一步：顯式建構 N×N read↔read 距離矩陣**（between-molecule），下游 clustering 與 PERMANOVA 全建立於此。

### C. 階層 clustering（CORE 2）
- `build_upgma`（預設）O(N³)；另 Ward / Single / Complete linkage（strategy pattern，`HierarchicalClustering.hpp:151-156`）。
- `TreeCutter`：依距離閾值切、依固定 k 切、或 `find_optimal_clusters` 用 **silhouette score 在 k=2..10 自動選最佳群數**（`HierarchicalClustering.hpp:198-206`）。

### D. PERMANOVA 結構檢定（CORE 3 — "結構不是失序" 的核心）
- SS 計算（`StructureTest.cpp:compute_ss` / `LabelTest.cpp:520-551`）：
  - `SS_total = (Σ_{i<j} d_ij²) / N`
  - `SS_within = Σ_group (Σ_{i<j ∈ group} d_ij²) / n_k`
  - `SS_between = SS_total − SS_within`
- pseudo-F = `(SS_between/(k−1)) / (SS_within/(N−k))`（`StructureTest.cpp:141-150`）。
- **999 次 label permutation** 建 null，p = (#{perm_F ≥ obs_F} + 1)/(999+1)（`StructureTest.cpp:186-205`）。
- 另有 **dispersion 同質性檢定**（各群到 centroid 平均距離的 one-way ANOVA-F），用來警示「群差是位置差還是離散度差」。

### E. LabelTest Δ = d_between − d_within（CORE 3 衍生，⚠ 與 Δβ 不同量）
- 對二元標籤（HP-merged / ALLELE / sample tumor-vs-normal）計算：
  - `within_mean` = 同群 read 對的平均距離；`between_mean` = 異群平均距離（`LabelTest.cpp:340-379`）
  - **Δ = between_mean − within_mean**（`LabelTest.cpp:213`）
- **只在 Δ>0 才跑 permutation**（Δ≤0 代表群內比群間還遠 = 反生物方向，直接 p=1）（`LabelTest.cpp:156-167`）。
- fine 4-group 還會輸出 6 條 pairwise 群間平均距離（`d_hp1_hp1s`、`d_hp1_hp2`…，`LabelTest.cpp:743-748`），以及 HP3/HP0 unassigned reads 的 **affinity score**（靠近 HP1 還是 HP2，含 permutation p，`LabelTest.cpp:751-898`）。
- 🔴 **務必區分**：這個 **Δ 是「距離」差（在 read-read 距離矩陣上）**，不是甲基率差。下面 F 的 **Δβ 是「甲基 β（rate）」差**。兩者是不同量、不同單位。

### F. Δβ（delta-beta）— 在 Python/tsg 層（CORE 6 配套）
> 主要實作在 `research/tsg_promoter_asm_reviewer/scripts/03_step4_ism_methylation_diff.py` + `pipeline/lib/cis_asm_core.py`，不是 C++ binary。
- **per-CpG β**：某 HP group 在某 CpG 上「甲基化 read 比例」。每 read 先去重取 max meth_call；ML 閾值 `>=200(~0.78)=甲基 / <=50(~0.20)=未甲基`（script 03:30-31）或 `THR=0.5`（cis_asm_core）。
- **Δβ**：把兩 group 的 per-CpG β **按 CpG 配對**，算：
  - `mean_delta = mean(β_som − β_germ)`（script 03:161-163）
  - `max_abs_delta`、**Cohen's d**（pooled SD）、**Wilcoxon signed-rank p**（配對）+ Mann-Whitney p（非配對 sanity）
  - 需每 group ≥5 CpG、每 CpG 每 group ≥3 reads（script 03:144,150）
- **兩條軸**：
  - **HP-axis** = HP1 vs **HP1-1**（同 germline 母單倍型，germline-tag vs somatic-重建-tag）→ somatic-controlled
  - **ALLELE-axis** = ALT vs REF
- **5mC/5hmC max-collapse**：`collapse_modtype` 把每 read/CpG 的 5mC、5hmC 兩列**取最大值合成 "any-modification"**（`cis_asm_core.py:31-32`）。

### G. normal-anchored cis-test（CORE 6 — 分辨真 cis vs drift）
（`scripts/34_control_loci_cohesion_cistest.py:108-152`，HP1 軸 per-CpG 配對）
- 三組 per-CpG β：**A=normal HP1**（突變前基線）、**B=tumor HP1**（germline-單倍型 reads，無 somatic allele）、**C=tumor HP1-1**（somatic-allele reads）
- 三個 delta：
  - `d_somatic = C − B`（tumor 內 ASM）
  - `d_cis = C − A`（somatic vs 正常基線）
  - `d_drift = B − A`（germline 單倍型自身漂移）
- **cis-candidate 判準**：`p_cis<0.05 AND |d_cis| > 1.8×|d_drift| + 0.02`（script 34:149）→ 真 cis-driven；否則 drift / weak。
- 設計上 held-constant CN/ploidy/alignment（同單倍型內比較）。

### H. copy-partition 分解（CORE 6 — 拆掉 copy confound）
（`scripts/37_copy_partition_confirm.py:33-45`）
- 四個量：
  - `d_HP` = HP1 vs HP1-1（原始；混入 copy）
  - `d_within` = HP1-1 內 alt vs ref（**真 somatic-allele cis，copy-controlled**）
  - `d_allele` = 全 ALT vs 全 REF（ALLELE 軸；也混 copy）
  - `d_copy` = HP1/ref vs HP1-1/ref（同 allele、純 copy/tag 差）
- **塌陷判準**：`|d_within| < 0.5|d_HP|` → cis 效應遠小於 HP 軸 = copy-partition confound 主導（script 37:64）。
- BRCA2 加 **per-read residual permutation**：HP1-1 內 alt-ref 觀測值 vs 2000 次隨機 split null（script 37:114-122）。

### I. Cramér's V 可靠性閘（CORE 5）
- `GlobalTest.cpp:128`：`cramers_v` + `cramers_v_reliable`（Cochran：列聯表最小期望格≥5 才可靠）；`min_cramers_v` 門檻 gate（`GlobalTest.cpp:141`）。稀疏表時 reliable=false（避免假高 V）。

---

## L3 — 引用的「已驗證結果」（標來源，供對照用，非本文件新跑）

> ⚠ 下列數字皆來自既有 **validated KB / ledger**，本文件只引用以利對照外部文獻；tier 與限制照原文件。

| 結果 | 值 | tier | 來源（validated）|
|------|----|------|------|
| strong-ASM 全基因組稀有率 | 0.34%（hypo 44% / hyper 56%，無方向偏好）| ⭐3 單樣本 | `knowledge/11_external_literature/02,07` |
| BRCA2 promoter HP-axis Δβ | −0.122 (n=197)；~80% copy + 邊際真 cis d_within=−0.023 (perm p=0.022) | ⭐3 | `knowledge/11_external_literature/07` |
| 唯一乾淨 cis exemplar | chr17/TBC1D16 d_within=0.142 (perm p=0.001) | ⭐3 | `knowledge/11_external_literature/07` |
| ASM × CN partial ρ | −0.055（**非** copy-driven）| ⭐3 pilot | `knowledge/11_external_literature/07` |
| 跨樣本 excess-over-null | 6/6 正（mean 0.168, 3 癌種）；somatic private 0/38 | ⭐3 | `knowledge/11_external_literature/07,08` |
| methyl→TP/FP filter | **NEGATIVE 四道**：LOSO 100% circularity / AUC=0.505 / strong-ASM anti-discriminative OR=0.194 | ⭐2 L4 DEAD | `knowledge/11_external_literature/05` |
| LOH-constrained phasing | NG=2 inner 93–99% same-HP 6/7；n=7 Wilcoxon W=28 p=0.0078 | 🟡 Grade B+ ⭐3 | `knowledge/11_external_literature/06` |
| O11 epipolymorphism | AUC 0.845→0.530（coverage 校正後 = artifact）| NEGATIVE | `knowledge/11_external_literature/02` |
| 5mC vs 5hmC | 5mC-driven；5hmC marginal (Δ 0.03–0.07) | ⭐3 | `knowledge/11_external_literature/07` |

---

## 與後續文件的關係
- **02** 外部工具 survey（modkit / 短讀 cancer DMR / co-methylation / ASM caller…）— web 實證
- **03** 方法對照矩陣（本文件 6 核心 × 各外部工具，逐格比）
- **04** 我們結果 × 外部論文數據觀察/結論 交叉比對
- **05** 可改進 / 可學習建議
- **06** Phase B 實機 benchmark 方案（下載運行）

## Provenance
- 方法描述：2026-06-09 直接讀 `src/core/{DistanceMatrix,HierarchicalClustering,StructureTest,LabelTest,GlobalTest,PerCpgAsm,NormalBaseline}` + `include/core/*.hpp` + `research/tsg_promoter_asm_reviewer/scripts/{03,34,37}*.py` + `pipeline/lib/cis_asm_core.py`，附 file:line。
- 結果數字：引自 `knowledge/11_external_literature/` validated 文件（已各標來源），非本檔新分析。
- 本檔為 method-comparison study part 1/6（in_progress）。
