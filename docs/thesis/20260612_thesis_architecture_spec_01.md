<!--
建立時間: 2026-06-12
報告類型: 論文架構 spec（thesis architecture design / 碩士論文藍圖）
任務類型: D handoff — 中文碩士論文撰寫架構，供逐章撰寫依循
狀態: design spec / 待用戶 review（已採方案 1：單一整合 characterization）；非正文
data_sources: docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md, docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md
參考: github.com/ythuang0522/Slide2Thesis（PPT→thesis 7 階段流程）；MethPhaser(Nat Commun 2024)/NanoMethPhase(Genome Biol 2021)/Tarabichi(Nat Methods 2021) 敘事模板；Swales CARS；Weber et al. 2019 benchmarking guidelines
-->

# 碩士論文架構 Spec — Subclonal reconstruction using somatic haplotagging and methylation profiles (ONT)

> **本文件職責**：把「研究成果 → 碩士論文（中文）」的**架構、格式、敘述慣例、素材對應、待補清單**定義清楚，作為逐章撰寫的藍圖。
> **不是正文**。正文逐章撰寫於 `docs/thesis/chapters/`（待 spec 核准後開工）。
> **決策（2026-06-12，用戶確認）**：
> - 骨架 = **方案 1**（單一整合 characterization 論文）；語言＝中文、層級＝碩士。
> - 標題定位 = **A 折衷（2026-06-12 用戶指定標題後改定）**：保留 "Subclonal reconstruction" 標題，但「reconstruction」歸 somatic haplotagging 骨幹之能力、「methylation profiles」為被 characterize 之有界附加層。（覆蓋稍早 D2=B；B 之紅線全數續守，因紅線約束的是甲基化主張而非骨幹重建。）
> - 撰寫順序 = **A 先寫不依賴 G-A 的部分**（Fig1 schematic + Table1 + Ch3 Methods + Results 單樣本 R1/R2/R3/R5 立即動工；G-A 跨樣本平行/稍後）。
> - **論文標題（用戶指定，2026-06-12）**：**Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing**（中文版可譯「以奈米孔定序之 somatic haplotagging 與甲基化圖譜重建腫瘤亞克隆結構」）。
> - 🔴 **此標題之 integrity 條件**：摘要首段必須立即點明分工——重建由 haplotagging 骨幹完成、甲基化之貢獻被誠實 characterize 且有界（germline-haplotype 強、subclone 僅存在性、非 variant filter）；否則純 reconstruction 標題會被誤讀為「甲基化重建了亞克隆」。
> **真值 SoT**：所有數字以 `VERIFIED_RESULTS.md`（V1-V12）為唯一來源；撰寫每章前先 Read 回真值再寫（§13 鐵則）。

---

## §0 三個必先確認的定位決策（請 review）

### D1. 論文類型 = characterization / 方法貢獻型（非 tool paper）— 已定
判準（PLOS Comp Biol 編輯口徑）：「把計算方法以外全拿掉，還剩夠多新穎性嗎？」我們**沒有 TP/FP-filter 的 win 可賣**（該方向 concluded DEAD），故**不走 tool paper**。貢獻 = 「**量化 per-read 甲基化在 somatic-haplotagging 骨幹上、不同解析度能加什麼**」這個科學洞察本身。領域對誠實有界結果加分（MethPhaser/EVOFLUx/Do2020 皆明標自身 bound）。

### D2. 🔴 標題與主張的「誠實張力」— 需你定奪
現行 foundation 標題寫 "**subclonal reconstruction**"，但已驗證紅線是「甲基訊號 = germline-haplotype 層級，subclone 層級僅存在性弱訊號、不可 read 救援」。純 "subclonal reconstruction" 標題會讓人**以為甲基化重建了亞克隆**，與真值衝突。三個標題定位候選：

| 候選 | 標題方向 | 主張 | 風險 |
|---|---|---|---|
| **A（折衷，建議）** | 「以 somatic haplotagging 重建亞克隆，並 characterize 甲基化的有界貢獻」 | 骨幹做重建（真）＋甲基化是被 characterize 的附加層（誠實）| 低；最對齊 MethPhaser「extend not replace」 |
| B（保守誠實）| 「Per-read 甲基化對 somatic-haplotagging 腫瘤分析的多解析度 characterization」 | 完全以 characterization 為主張，不提 reconstruction | 較不吸睛，但最防 overclaim |
| C（保留現標題）| 維持 "Subclonal reconstruction using ... methylation profiles" | 強調 reconstruction | 🔴 高 overclaim 風險，須在摘要首段立刻 bound，否則誤導 |

> **✅ 最終定 = A（用戶指定 "Subclonal reconstruction" 標題後，2026-06-12 改定；覆蓋稍早 B）**：標題保留 reconstruction，但**「reconstruction」明確歸 somatic haplotagging 骨幹之能力**（生物學成立、不 overclaim），**「methylation profiles」為被 characterize 之有界附加層**。論文貢獻 = 特徵化甲基化在骨幹上之加值與解析度上限。動詞紀律：對**甲基化**仍 **characterize / quantify the bounded contribution / extend**、**禁 enable / improve / reconstruct（指甲基化時）**；reconstruction 一詞僅用於描述 haplotagging 骨幹。摘要首段即點明此分工（見上 integrity 條件）。

### D3. 證據 tier 現況（影響「能宣稱多強」）— 需知悉
| 研究線 | 樣本範圍 | tier | 缺口 |
|---|---|---|---|
| 甲基救相位 / longphase-S tag 輔助（V1-V12）| **單樣本 HCC1395** | ⭐3 | 🔴 **G-A：V10/V11c 跨 6 樣本未跑**（內容缺口，非格式）|
| ASM characterization（ZAR1L/BRCA2、cross-sample ASM）| 6 樣本 × 3 癌種 | ⭐3 | 同位點 private（0/38），靠 excess-over-null 6/6 復現 |
| ISM 完整 TP/FP/FN（存在性+cis）| 多樣本 | ⭐4 | ASM 廣存但非 usable filter |

> **意涵**：論文若以 V1-V12 為主軸，**現狀只能宣稱單樣本 ⭐3**。要在論文裡宣稱「跨癌種復現」必須先跑 G-A（見 §6 待補 P0）。這是寫作前要決定的：**(a) 先跑 G-A 再寫（衝 ⭐4）** 還是 **(b) 先以單樣本 ⭐3 寫初稿、G-A 結果後補**。

---

## §1 誠實紅線（全文一致，不可違反 — 來自 VERIFIED_RESULTS + foundation §4）

撰寫任何章節時，下列每條都必須在對應位置誠實呈現，**不可在 polish/expand 階段被沖淡**：

1. **甲基 = germline-haplotype 層級**（V10 決定性）：分不同 haplotype 強、within-haplotype 弱。
2. **絕對 AUC 有方法樂觀成分**（V3/V4/V10）：normal ~0.98 全基因組不可能全是真 locus-specific ASM。**報告一律引相對 null 效應**（real 超 shuffle p95 median +0.260；V6 held-out 88.5% vs null 52.4%），**不引絕對 0.87–0.98**。
3. **T1 headline = V6 全基因組外推 0.885（null 0.524）**，**不是** 0.926（後者是事後篩可分位點的條件值，不可當 headline）。
4. **T2 OVERSTATED**：只證「有 germline 真值的 1-1 vs 2-1 可分（AUC 0.90）」；「用甲基歸 H3」**未驗證**（H3 無真值）。
5. **T3 = 存在性窄翻案、可用性 NEGATIVE**：亞群 ALT vs 母本 REF 甲基「可分」（farCpG AUC 0.85/0.87 vs 噪音 ~0.50），但 ambiguous read 偏向**母本**（8/8 frac <0.5）→「可分」≠「可救 ambiguous read」。勿無條件説「T3 翻案」。
6. **絕非 variant TP/FP filter**（concluded DEAD）；是 subclonal characterization。
7. **CpG-SNP pseudo-ASM 尚未完全排除**（V3/V4 最大嫌疑）→ 須在 limitations 誠實標。
8. **single-pipeline 自我參照**：被矯正的 HP tag 本身是 longphase-S 產物（T2/T3 self-reference）→ tier 硬上限 ⭐2-3，需正交對照（G-E）。
9. **救援適用性窄**（V12）：unphase 88% 有甲基，但僅 ~6% 在有 HP1/HP2 錨點區可嘗試，94% 無本地參考無法救。88.5% 是「可嘗試池內」正確率。
10. **cohesion ≠ cis**（V11c）：farCpG ±100bp 太窄，subclone-program vs 突變 cis 足跡未分。

---

## §2 章節架構（中文碩士論文 6 章 + 前後件）

採 Slide2Thesis 6-章 taxonomy ＝ 中文碩士論文標準骨架。每章標「purpose / 敘述慣例 / 素材對應（含真值 V#）」。

### 前件
- **封面 / 書名頁 / 口試委員審定書**（系所模板）
- **中文摘要**（Background-Methods-Results-Conclusions；ABT 鉸鏈；headline 數字＋誠實天花板都要進）
- **英文 Abstract**（250–300 字，同結構）
- **誌謝**
- **目錄 / 圖目錄 / 表目錄**

### 第一章　緒論（Introduction）— CARS 漏斗
- **Purpose**：從廣（腫瘤異質性、ONT 長讀甲基化、somatic haplotagging）→ gap（甲基化對 somatic 變異骨幹的「加值」未知、且其解析度上限未被量化）→ 佔據 gap（研究目的＋預告有界發現）。
- **敘述慣例**：5–6 段，逐步收窄；末段明列研究目的（Q：甲基化在 germline-haplotype / subclone / variant-filter 三個解析度各能加什麼？）＋論文架構預告。
- **素材**：landscape `01_TO_FP問題全貌` / `02_Self_Phasing根因`（phasing 問題背景）；foundation §1 一句話定位；method_comparison `12_methyl_approaches_by_axis`（業界 4 軸）。
- **Figure**：（可選）Fig 1.1 研究問題示意（與 Fig 1 schematic 共用）。

### 第二章　文獻探討（Related Works）
- **Purpose**：建立讀者理解所需背景 + 定位本研究於業界地景。
- **次節建議**：
  1. 腫瘤亞克隆重建的基因骨幹（Tarabichi Nat Methods 2021 workflow：mutation→CNA→VAF→CCF→clustering→phylogeny；PyClone-VI）。
  2. ONT 原生甲基化與 haplotype-resolved methylation（NanoMethPhase、MethPhaser「甲基延伸 SNV phasing」）。
  3. 癌症中的 allele-specific methylation（Do2020、Martin-Trujillo；ASM×LOH）。
  4. 甲基追蹤腫瘤演化的可能與邊界（EVOFLUx Nature 2025 — aspirational sibling，但須對比：它宣稱甲基加**時間/演化**解析；我們誠實發現甲基**不**加可用 subclone 判別，只存在性）。
  5. 甲基分群/距離方法地景（modkit/DSS per-position 率差 vs cvlr/ASMS/qFDRP read-level；ISM 站「read-read 距離+clustering」軸）。
- **素材**：method_comparison `02_external_tools_survey`（82 工具）、`12_methyl_approaches_by_axis`、`03_method_comparison_matrix`；landscape `00_INDEX`。
  - ⭐ **外部驗證庫（2026-06-13，撰 Ch2 主用）**：`InterSubMod/docs/method_comparison/20260613_external_validation_library_index_01.md` → 指向 repo 外 `/big7_disk/liaoyoyo2001/external_validation/`（**49 源親讀 CONTEXT 卡**：28 已驗 + 21 P0/P1 audit-driven 補卡，含 Ch2 逐軸對應表）。
  - 🔴 **文獻缺口稽核（framing 凍結前必讀）**：`InterSubMod/docs/reports/research_landscape/20260612_external_validation_literature_gap_audit_01.md` — framing 裁決：①「subclonal reconstruction」用詞偏強，正文先用「**characterize methylation-associated subclonal structure conditioned on somatic haplotags**」（除非 G-D 重建 demo + G-E 正交完成）②距離式/read-level 甲基非天然新穎（Sgootr/Gaiti 為前例）③6 cell line=跨模型 reproducibility **非病人 cohort**（Methods 記 subline/passage/basecaller）。
  - **Ch2 引用守則**：① ISM vs cvlr/ASMS/MethylBERT 差異**禁**寫「平台（二代）」或「缺顯著性檢定」（皆錯，會被打臉）→ 改述「無監督距離矩陣結構 PERMANOVA + normal-baseline cis-test + somatic-subclone 目標」；② 與 49 源**無同 regime 直接衝突，但有 framing tensions**（非「0 conflict」）；③ LongHap 2026=germline-only 不威脅白地、反佐證 R2。每源 citation/DOI 見其 CONTEXT 卡 §7。

### 第三章　材料與方法（Materials & Methods）— 最詳盡，可重現
- **Purpose**：讓獨立讀者能複現。
- **次節**：
  1. **樣本與資料**：6 cell line × 3 癌種（HCC1395/1937/1954 乳腺、H1437/H2009 肺、COLO829 黑色素瘤）；tumor+matched normal；ONT R10 5mCG。🔴 **缺 Table 1（樣本統計：purity/depth/ploidy/mutations/CNV）— 待補**。
  2. **somatic haplotagging 骨幹**：longphase-S tag cascade（`HaplotagStrategy.cpp:452` judgeSomaticReadHap + `SomaticHaplotagProcess.cpp:461` inheritHaplotype，門檻 0.6）；6-state 語意（unTag/H1/H2/H1-1/H2-1/H3）；unphase pool 45.84%（V1）。
  3. **per-read 甲基化萃取**：BAM MM/ML tag → 5mCG 機率矩陣（`MethylationParser`；strand-flip 處理；二值化 0.2/0.8）。🔴 **MM/ML 解析零單元測試 — code audit 待補 golden test**。
  4. **ISM 結構分析 6 核心**（`01_ism_method_spec_from_source` 有 file:line）：NHD/L1/Bernoulli/Jaccard 距離（`DistanceMatrix`）、UPGMA/Ward 聚類+silhouette（`HierarchicalClustering`）、PERMANOVA 999 排列（`StructureTest`）、per-CpG Fisher+NME+epipolymorphism（`PerCpgAsm`）、Cramér's V+Cochran gate（`GlobalTest`）、normal-baseline cis-test（`NormalBaseline`）。
  5. **甲基救相位 / tag 輔助方法**：anchor AUC（leave-one-out 質心）、shuffle-label null、held-out PS-block 外推（V6 設計）、allele-based REF/ALT AUC（V10 不靠 HP tag）。
  6. **null / 對照設計**：germline-het null（V3）、shuffle-label（V4）、matched-normal copy-clean 對照（V10）、SEQC2 CN/LOH confound 控制（V5）。
  7. **統計**：效應量＋CI（非僅 p）；多重比較校正；獨立 test set（避免純 LOO — OUP 退稿）。
  8. **軟體版本/參數 + Data/Code Availability**。🔴 **待補**。
- **Figure**：**Fig 1 方法總覽 schematic（待製，最高優先）** — tumor reads→longphase-S haplotag(HP1/HP2+子單倍型)→疊 per-read 5mCG→問「不同解析度加什麼」（MethPhaser Fig1 + Tarabichi Fig1 hybrid）。
- **可重現素材**：20+ Python 腳本（`a0a_hp_distribution.sh`/`extract_per_read_methyl.py`/`germline_het_null.py`/`shuffle_control.py`/`extrapolation_validation.py`/`allele_asm_auc.py`/`rigor_t1-3.py`/`t3_local_allele.py`/`unphase_inventory.py`…）；C++ 23 模組；method_comparison `08_external_program_audit`（modkit/DSS 互驗 r=0.98）。

### 第四章　結果（Results）— 降強度階梯（一圖一節，先報不解釋）
按「解析度 / 證據強度」由強到弱排，最強領頭、NEGATIVE 主動擁有：

| 節 | 訊息（會說話的標題）| 真值（SoT）| Figure 候選（實檔）|
|---|---|---|---|
| **R1** 骨幹與救援池 | somatic haplotagging 把 tumor read 分到 germline 單倍型，留下 45.8% unphase 救援池 | V1（30.15M read，unphase 45.84%）| （表）HP 分布表 |
| **R2** 甲基化在 germline-haplotype 層級**強**分離 | 甲基化攜帶真 haplotype-linked 訊號，且非 copy artifact | V10（normal 0.979≥tumor 0.866，6/6 chr；depth-matched 否證；GNAS Δ=0.49 正控）；V3/V4（95% het 區 real 超 shuffle p95，+0.260）| `fig_allele_asm_tumor_vs_normal.png`；`brca2_core_heatmap.png`（+gnas/h19 對照）|
| **R3** 甲基化能把 unphase read 救回真 HP | held-out 外推救援 88.5%（null 52.4%，77% 窗顯著）| V6（acc median 0.8852，AUC 0.958）；V12（適用性窄：可嘗試池 ~6%）| `fig_extrap_hist.png` / `fig_extrap_perchr.png`；`showcase_chr20_rescue.png` |
| **R4** 甲基化判 germline 分支可行、但 subclone 層僅**存在性** | T1 救 unphase ✅；T2 分 1-1/2-1 ✅但不可外推 H3；T3 亞群可分但不可救 ambiguous | V11/V11b/V11c（T1 reframe→0.885；T2 OVERSTATED；T3 farCpG 0.85/0.87 vs 0.50 但 ambiguous 8/8<0.5）| `fig_rigor_corrected.png`；`fig_assist_3targets.png` |
| **R5** CN/LOH 不是 confound；LOH 區反而分不開 | CN vs AUC 無相關（rho=0.035）；LOH 顯著低（p=0.0019）→ 副產品：低 AUC 可偵 LOH | V5 | `fig_seqc2_cn_vs_auc.png` / `fig_seqc2_auc_by_status.png` |
| **R6** 跨樣本復現（⭐4 證據）| ASM 同位點 private（0/38）但 6/6 excess-over-null>0 跨 3 癌種 | cross-sample ASM 線（記憶 SoT）🔴 **V10/V11c 跨樣本=G-A 待跑** | 跨樣本彙整圖 🔴 **待製** |
| **R7** owned NEGATIVE：非 usable variant filter | ASM 廣存但 TP≈FP 不判別、cis FP 最高 | ISM 完整 TP/FP/FN 線 | NEGATIVE 摘要表 🔴 **待製** |

### 第五章　討論（Discussion）— 先最強、再 bound、再限制
- **開頭**：單一最強 takeaway（甲基化在 germline-haplotype 層級強分離，跨樣本/癌種穩健）。
- **再 bound**：subclone 層僅存在性、不可 read 救援（V11c）；機制＝甲基是 germline-haplotype 層級，tumor-acquired ASM 不足製造強 subclone-特異 epiallele。
- **對齊文獻**（method_comparison `04`）：Δβ 方法、ASM 方向口徑差、structure vs disorder（Landau 2014 支持）、V9 aDMR×LOH「數字符合 ≠ 機制驗證」的誠實案例。
- **專段 bounded-claims（仿 MethPhaser limitations 段）**：單樣本/單 pipeline 自我參照、CpG-SNP pseudo-ASM 未排除、絕對 AUC 樂觀、cohesion≠cis、3 癌種 scope。
- **implications / future work**：G-B/G-C 對照、G-D 重建 demo、G-E 正交 pipeline、beta-binomial 改進。

### 第六章　結論與未來工作（Conclusions & Future Work）
- 2–4 句 take-home：甲基化**加了什麼**（germline-haplotype 解析的 read 救援/分支判別）＋**沒加什麼**（subclone 判別、variant filter）。

### 後件
- **參考文獻**（🔴 待建 .bib，過 `/citation-verification`）
- **附錄**：Supplementary Methods（PERMANOVA/normal-cis-test pseudo-code）、Supplementary Tables（816 aDMR 詳單、cross-sample）、13 個 standalone HTML explainer 連結、benchmark 細節（Phase B）。

---

## §3 圖表計畫（descending-strength；Fig 1 = schematic）

| 編號 | 類型 | 內容 | 狀態 |
|---|---|---|---|
| **Fig 1** | 方法 schematic | haplotag 骨幹 + 5mCG 疊加 + 三解析度問題 | 🔴 **待製**（最高優先；可用 `/methods-example` 或 image-gen）|
| Fig 2 | 結果 | V10 normal>tumor 不是 copy（決定性）| ✅ `fig_allele_asm_tumor_vs_normal.png` |
| Fig 3 | 結果 | V6 救援外推 88.5% | ✅ `fig_extrap_hist.png` / `fig_extrap_perchr.png` |
| Fig 4 | 結果 | V11 三 target（含 T3 NEGATIVE 面板）| ✅ `fig_rigor_corrected.png` / `fig_assist_3targets.png` |
| Fig 5 | 結果 | V5 CN 非 confound + LOH 反而低 | ✅ `fig_seqc2_*.png` |
| Fig 6 | 結果 | 跨樣本復現（excess-over-null）| 🔴 **待製**（依 G-A）|
| Table 1 | 表 | 6 樣本統計（purity/depth/ploidy/CNV）| 🔴 **待製** |
| Table 2 | 表 | ISM vs 外部工具對照矩陣 | ✅ `03_method_comparison_matrix` 可轉 |
| Table 3 | 表 | NEGATIVE 結果摘要（T3/T2/filter）| 🔴 **待製** |
| Supp Figs | 附錄 | IGV showcase（chr8/15/20）、cross-sample landscape | ✅ 既有 PNG/HTML |

> 圖檔路徑：`docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/*.png`。撰寫時用 pandoc-crossref `![cap](file){#fig:id}` + `@fig:id`，不手寫圖號。

---

## §4 摘要結構（中英對照 — ABT + Background-Methods-Results-Conclusions）

骨架（填數字前先 Read VERIFIED_RESULTS）：
1. **Background/territory**：腫瘤亞克隆重建靠 somatic 變異骨幹；ONT 同時給 per-read 甲基化。
2. **However（gap）**：甲基化對此骨幹的「加值」與其**解析度上限**未被量化。
3. **Methods**：6 cell line × 3 癌種；longphase-S haplotag + per-read 5mCG；anchor AUC / held-out 外推 / matched-normal 對照 / PERMANOVA。
4. **Results（含 bound）**：甲基化在 germline-haplotype 層級強分離（非 copy；救援外推 {{V6}}）；subclone 層僅存在性、不可 read 救援；非 variant filter。
5. **Conclusions**：甲基化為 somatic-haplotagging 骨幹加上 germline-haplotype 解析的補強，但 subclone 判別有解析度天花板。

> 動詞紀律：**quantify the bounded contribution / extend / characterize**，禁 enable/improve/rescue-all。

---

## §5 Slide2Thesis 轉換流程 + §13 反捏造防線

**借用其三遍轉換法**（PPT/報告 → 學術散文）：
1. **generate**：把投影片 bullet / 報告重點 → 連貫學術散文（正式研究生語氣、markdown 小標、**不寫章節號**、剝除「投影片 N 説」痕跡）。
2. **expand**：拿原始報告 diff，補回漏掉的 Methods/Results 細節（**只補真實有的，不 fill gaps with 幻覺**）。
3. **polish**：去投影片語氣 + 統一術語。**polish 不可沖淡 §1 紅線**。

**🔴 必加防線（Slide2Thesis 沒有）**：
- **數字**：每個進正文的數字先 Read VERIFIED_RESULTS 回真值，**產數字與寫正文不同 batch**（§13.0）；frontmatter 標 `data_sources:`；validated 階段過 `number_provenance_check.sh`。
- **引用**：禁 LLM 自動挖 PubMed；每條過 `/citation-verification`（WebSearch+Scholar）才入 .bib。
- **NEGATIVE**：T3/T2/filter 的 bound 手動疊上，確保 polish 後仍在。
- **圖**：per-read 甲基圖需手動 figure-matching（非 auto-crop）。

---

## §6 待補 checklist（prioritized）

### P0 — 內容缺口（影響能宣稱多強，非格式）
- [ ] **G-A 跨樣本**：V10（matched-normal not-copy）+ V11c（local-allele 存在性）跑遍 6 樣本（normal 甲基 5/6 ready，COLO829 缺）→ R2/R4/R6 升 ⭐4。**決定 (a) 先跑再寫 vs (b) 單樣本 ⭐3 先寫後補**。
- [ ] **G-B**：within-haplotype somatic-vs-baseline 正確 null（非 germline-het null）→ reopen T3 前置。

### P1 — 高 ROI 格式/素材缺口
- [ ] **Fig 1 方法 schematic**（最高優先；`/methods-example`）。
- [ ] **Table 1 樣本統計**（purity/depth/ploidy/mutations/CNV）。
- [ ] **Data Availability + Code Repository**（GitHub + commit hash；BAM/VCF/甲基矩陣存放）。
- [ ] **參考文獻 .bib**（過 citation-verification）。

### P2 — 補強
- [ ] 跨樣本彙整圖（Fig 6）、NEGATIVE 摘要表（Table 3）。
- [ ] Phase B benchmark（modkit/cvlr vs ISM；待核准）→ Results/Supplementary。
- [ ] Supplementary Methods pseudo-code。
- [ ] code audit 3 必修（MM/ML 解析測試 / Fisher over-dispersion / fill_report 洞）— 影響 Methods 可信度。

---

## §7 撰寫順序建議（spec 核准後）

1. **定 D2 標題 / D3 先跑 G-A 或先寫**（本 spec review 一併決定）。
2. **Fig 1 schematic + Table 1**（圖表先行，Results 才有依附）。
3. **第三章 Methods**（素材最齊、最不依賴 G-A）→ 先寫可立刻產出。
4. **第四章 Results**（一節一圖，先 R1/R2/R3/R5 單樣本已驗證者；R4/R6 視 G-A）。
5. **第一章 緒論 + 第二章 文獻**（Results 定型後回填 gap 敘事）。
6. **第五章 討論 + 第六章 結論**（最後，含 bounded-claims 專段）。
7. **摘要（中英）+ 參考文獻 + 附錄**（收尾）。

> 對齊 Slide2Thesis：每章 generate→expand→polish 三遍；每章撰寫前 Read VERIFIED_RESULTS 對應 V#。

---

## §8 風險 / 必守紅線總表（一頁速查）

| # | 紅線 | 出現章節 |
|---|---|---|
| 1 | 甲基=germline-haplotype 層級，subclone 弱 | Abstract/R2/R4/Discussion |
| 2 | 引相對 null 效應，非絕對 AUC | R2/R3/Methods |
| 3 | T1 headline=0.885 非 0.926 | R3/R4 |
| 4 | T2 不可外推 H3 | R4/Discussion |
| 5 | T3 存在性窄翻、可用性 NEGATIVE | R4/Discussion |
| 6 | 非 variant filter | Abstract/R7/Discussion |
| 7 | CpG-SNP pseudo-ASM 未排除 | Methods/Discussion limitations |
| 8 | single-pipeline 自我參照 ⭐3 上限 | Discussion limitations |
| 9 | 救援適用性窄（~6% 可嘗試）| R3/Discussion |
| 10 | cohesion≠cis | Discussion limitations |
| — | 標題不得讓人以為甲基化重建亞克隆（見 D2）| Title/Abstract |
