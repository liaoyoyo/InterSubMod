<!--
建立時間: 2026-06-13
狀態: reference (外部研究 vs ISM 五桶對照 + 論文敘述評估; 已過 4-lens 對抗稽核並套用修正)
報告類型: external_vs_ism_comparison_for_paper
受眾: 廖子游 · PI · 碩論 Ch2(文獻探討) + Discussion 撰寫
framework: 五桶分類(優勢/相似/不同/矛盾/佐證) + Assertion-Evidence + tier 分層
data_sources:
  - /big7_disk/liaoyoyo2001/external_validation/_landscape/05_context_cards_index.md (52源逐源關係標籤)
  - /big7_disk/liaoyoyo2001/external_validation/_landscape/06_source_read_corrections.md (源碼親讀校正 + 3投稿口徑)
  - InterSubMod/docs/reports/research_landscape/20260612_external_validation_literature_gap_audit_01.md (framing-tension reframe)
  - InterSubMod/docs/method_comparison/20260609_ism_vs_external_methylation_tools/01_ism_method_spec_from_source.md (ISM 6 cores + file:line + L3 結果表)
  - workflow wf_b40031c1-517 (11 agents: 6軸複驗 + 合成 + 4 lens 對抗稽核)
provenance_note: 對照結論逐源溯 external_validation CONTEXT 卡(一手親讀 repo 源碼/全文)。ISM 結果數字引自 method spec L3 表(各標 validated 來源)。本報告為「組織既有 verified 內容」, 不產新數字; 與任何分析不同 batch(§13.0)。對抗稽核 3-lens NEEDS_WORK 之修正已逐項套用(見 §對抗稽核)。
-->
<!-- provenance-verified: 數字(0.34%/Δβ=-0.122/d_within=-0.023/AUC=0.505/82-92%/18.2%)皆引自上列 validated 來源, 撰寫與分析分離; reviewer-trap 口徑對齊 06_source_read_corrections F-1/F-1b。 -->

# 外部研究 vs ISM — 五桶對照與論文敘述評估

> **這份是什麼**：把 52 個外部驗證源（`external_validation/`，一手親讀 repo 源碼/全文）逐源獨立複驗後，歸成你問的五桶 —— **優勢/強項、相似、不同、矛盾、可佐證** —— 並評估論文敘述方式。已過 4-lens 對抗稽核並套用修正。供逐節確認。
> **方法**：11-agent workflow（6 軸各 1 agent 重讀實際 CONTEXT 卡複驗關係 → 合成五桶 → 4 lens 對抗檢查 reviewer-trap / 優勢防守性 / 矛盾誠實 / 佐證有效）。主 agent 再用自讀的 ISM spec + landscape 05/06 獨立複核 load-bearing 數字。

---

## L0 — 一眼結論

**五桶分布（52 源逐源獨立判定）**：

| 桶 | 數 | 一句 |
|---|:--:|---|
| 🟢 **可佐證**（CORROBORATES）| **25** | 外部數據**反向支撐** ISM 內部結論（A-framing 骨幹分工 / R2 甲基=germline-haplotype 層級 / cancer 甲基-phasing 白地 / 弱-subclone 紅線 / 甲基非 filter）|
| 🔴 **矛盾/限制**（TENSION）| **12（+2 補入=14）** | 不是實驗結果衝突，是 **framing tension** —— 限制 ISM「甲基新穎」與「reconstruction」用詞的可宣稱範圍 |
| ⚪ **不同**（DIFFERENT）| **11** | 不同 regime/目標（井水不犯河水，不威脅紅線）|
| 🟡 **相似**（SIMILAR）| **4** | 與 ISM 同軸最近的 read-clustering 對手（cvlr/ASMS/qFDRP/DAMEfinder）|

> **一句總評**：外部文獻**整體支撐** ISM 的 A-framing（reconstruction 歸 somatic-haplotagging 骨幹、甲基被 characterize 有界），但對「甲基新穎」施加**嚴格邊界**。ISM 可防守優勢 = **3 件特定組合**，且對抗稽核確認此三件**僅是 characterization 層能力、不是 variant utility 優勢**（TP/FP 判別 AUC≈0.50、filter-DEAD）。
>
> **gap audit 鐵則（必守）**：**無同 regime 直接實驗衝突，但有多條高影響 framing tension** —— 絕不可裸寫「0 conflict」。

---

## L1 — 三個必守 reviewer-trap（投稿前先記住，否則被打臉）

> 這三條來自源碼親讀（`06_source_read_corrections` F-1/F-1b），是**最容易自爆**的點。對抗稽核發現連我們自己的優勢敘述都差點犯 T2 的鏡像版（見 §對抗稽核）。

| Trap | ❌ 禁寫 | ✅ 為何錯 + 正確措辭 |
|---|---|---|
| **T1 平台** | 「對手用二代定序」 | cvlr/ASMS/MethylBERT/ClairS-TO/DeepSomatic/Severus/PoreMeth2/EVOFLUx/MethPhaser/LongHap **全 ONT/PacBio-capable**（MethylBERT `cli.py:94 -m dorado` 直讀 MM/ML）。差異放「方法/目標/檢定對象」，不放平台世代。對 bisulfite/array 源（DAMEfinder/qFDRP/Do2020）只標「cross-platform 數字不可並列」。|
| **T2 顯著性** | 「對手缺顯著性檢定」 | cvlr（median\|Δmeth\| vs 1000 random + Fisher）/ASMS（`pval.rs` 5000-perm + `cluster_by_snp` 1000-perm+BH）/DAMEfinder（bumphunter perm）/EVOFLUx（貝氏 nested sampling）**都有**。正確：ISM 差異是**不同檢定對象**（PERMANOVA 測 read×read 距離矩陣*結構分離* vs 對手測 cluster 間甲基*差量*）。⚠ **qFDRP 例外**：它確實**無**檢定（純摘要統計）—— 但對 qFDRP 的差異走 structure-vs-disorder，**不可**用「有檢定 vs 沒檢定」（否則 reviewer 拿 qFDRP 來打你「你說對手皆有」）。|
| **T3 用詞/新穎** | 「用甲基重建 subclone」/「距離式甲基重建是新穎 method-class」/「0 conflict」/「甲基不能定義 subclone」 | Sgootr（距離式甲基 lineage tree, 真腫瘤單細胞）/Gaiti2019/Epiclomal/MethylTree/colorectal-2026 已在**單細胞或 bulk regime 證甲基能重建 lineage**；PoreMeth2 已佔「長讀甲基演化」；ClairS-TO Verdict 已做 variant-AF subclonal 分類。把這些降為 **prior art 並 cite**；ISM 新穎收斂三件；reconstruction 用詞**只掛 somatic-haplotagging 骨幹**，甲基面一律用 **characterization**。|

---

## L2 — 五桶逐桶（已套用對抗稽核修正）

### 桶① 優勢/強項（ISM 可防守的三件）

> 🔴 **對抗稽核最重要的修正**：這三件**只是 characterization 層能力，不是 variant utility 優勢**。TP/FP 判別已證 NEGATIVE/anti（AUC≈0.5049、strong-ASM OR=0.194 FP enriched 5×、filter-DEAD）。論文中「advantage」一詞**絕不可漂移成「改善 variant calling / reconstruction 效能」**。每件都標「characterization capability, orthogonal to variant TP/FP utility」。

| # | 可防守優勢 | 源碼 | ⚠ 修正後的誠實邊界（必帶）|
|---|---|---|---|
| **A1** | **無監督 read×read 距離矩陣的*結構*檢定**：顯式建 N×N 距離矩陣 → 升格為分群引擎 → PERMANOVA pseudo-F（999 perm）測 read 在距離空間是否結構分離（非 cluster 甲基差量、非 within-read disorder）| `StructureTest.cpp:141-205` | (a) NHD per-pair kernel 與 **qFDRP 數學同式** → Methods 須 cite qFDRP（Scherer 2020 NAR）為先例，新增值在 downstream（保留全矩陣→分群→結構檢定）非 kernel。(b) **每個元件都是現成的**（NHD=qFDRP、階層分群=教科書、PERMANOVA=Anderson 2001 生態學標準）→ 防守靠的是「**組合** + somatic-haplotag conditioning + normal anchor + somatic-subclone 目標」，非任一單元件。(c) 🔴 **structure-vs-magnitude 須實證**：reviewer 會問「給我一個 pseudo-F 與 median\|Δmeth\| 不一致的位點」。若在乾淨兩-allele 群上兩者高度共線，這區分形式成立但實務空洞 → **投稿前須 demo 至少一個 divergence 位點**，否則框成「在保留矩陣上測不同統計量」而非「判別優勢」。|
| **A2** | **normal-baseline somatic cis-test**：以 matched-normal per-CpG 基線扣除（residual=raw−normal mean）+ 三角 d_somatic/d_cis/d_drift 分辨真 cis vs drift。**所有 germline 工具（cvlr/ASMS/MethPhaser/NanoMethPhase/qFDRP）與 modkit 皆無 normal anchor** —— 這是 ISM **真正獨有**的一件 | residual 在 `NormalBaseline.hpp:38-67`；三角在 `scripts/34_*.py:108-152` | (a) 🔴 **歸屬**：C++ engine **只做 residual**；d_cis/d_drift 三角 + cis-candidate 判準在 **Python scripts/34**（分析層，hardcoded path / 無 regression test）→ 論文**不可**把 cis-test 説成 compiled-engine feature。(b) 🔴 **可靠性降級**：旗艦 BRCA2 表面 d_cis=−0.142 但 ~80% 是 subclone/copy confound，真 focal within-somatic 只 **d_within=−0.023（perm p=0.024，marginal，% split 不 robust）**；乾淨 exemplar 是 **chr17/TBC1D16 d_within=0.142（p=0.001）**。→ 可防守的是「**分解層本身**」（暴露 raw 差有多少是 copy/subclone 背景），**不是**單位點穩健 cis 定量。|
| **A3** | **somatic-subclone characterization 目標**（tumor-normal paired, read-level）：對手是 germline ASM/imprinting 發現（cvlr/ASMS/MethPhaser）或 supervised tumour-vs-normal deconvolution（MethylBERT 需 fine-tuned model + DMR 標記）| `01_ism_method_spec` | (a) 用「somatic-subclone **structure characterization**」非「subclone reconstruction」。(b) MethylBERT 分的是 tumour-vs-normal（較大訊號）非 subclone-vs-subclone → 不構成「subclone 弱」紅線反例。(c) 🟡 unsupervised 標籤**只適用甲基分群步**：整條 pipeline 仍 conditioned on 上游 HP tags（longphase-S）+ matched-normal mask → 引用時帶此限定，MethylBERT 對照**以「目標差」（結構 characterization vs fraction 估計）領頭**，而非 supervised/unsupervised 軸（後者可被模糊）。|

### 桶② 相似（同軸最近對手 + ISM 邊際新增值）

| 外部 | 重疊 | ISM 邊際新增值 |
|---|---|---|
| **cvlr** (Raineri 2023) | 同 ONT read×CpG 二元甲基矩陣底材、同找 read 甲基子群、同 no-phasing 精神。ASMS 直系前身 | cvlr=model-based EM soft k 群、測 cluster 間 median\|Δmeth\|；ISM 把距離矩陣升格分群引擎 + PERMANOVA 測*結構* + normal-baseline cis-test（cvlr 無）+ somatic 框架（cvlr 純 germline/imprinting）。⚠ cvlr 有 randomization 檢定，禁寫缺 |
| **ASMS** (Raineri 2024) | **單一最像**：同底材、同 no-phasing | ASMS 硬編碼 K=2 soft EM（無法 >2 subpop），距離僅 silhouette 副產品（`utils.rs:183`）；ISM 全距離矩陣 + UPGMA/Ward 可 >2 群 + 升格分群引擎。⚠ ASMS 有 2 permutation 檢定，禁寫缺 |
| **qFDRP** (Scherer 2020) | **同一 NHD kernel**（per-pair normalized-Hamming 完全同式）→ ISM pairwise 距離的 NAR 同儕審查先例 | qFDRP 把距離**塌成 per-CpG scalar**（只量 disorder magnitude）；ISM 保留全矩陣→分群→PERMANOVA 測 structure。⚠ qFDRP **無顯著性檢定**（摘要統計）；bisulfite 短讀（cross-platform 數字不可並列，但非 T1，ISM 從未宣稱它用二代）|
| **DAMEfinder** (Orjuela 2020) | 軸3 方法級最像：allele/共甲基結構直覺（tuple）；ISM `PerCpgAsm.cpp:6` 實際引用它 | ISM 保留 read×CpG 維度做無監督距離（DAMEfinder 先聚合成 per-tuple/region 率，丟 read 身份）+ normal-anchored cis（DAMEfinder germline het-SNP 主、無 normal baseline）。⚠ DAMEfinder 有 bumphunter perm。🔴 `PerCpgAsm.cpp:6` 註解寫「Per-CpG Fisher: DAMEfinder」但 DAMEfinder 原生**非 Fisher**（SNP-mode 率差/tuple-mode odds-ratio）→ 方法用詞須改「per-CpG ASM 概念源自 DAMEfinder，我方以 Fisher exact 實作」（同行的 pycoMeth 統計也須先查再引）|

### 桶③ 不同（不同 regime/目標，井水不犯河水）

| 外部 | 維度差 | 註 |
|---|---|---|
| **MethylBERT** (Jeong 2025) | supervised deconvolution（估 tumour fraction、需 fine-tuned model + DMR 標記）vs ISM 無監督結構 characterization | ⚠ 碼層支援 ONT（T1 禁寫只能二代）；分 tumour-vs-normal 非 subclone（不構成紅線反例）|
| **Dorado + Modkit** | 上游 provenance（raw→modBAM；pileup=per-position 跨 read 聚合率，塌 read 身份）vs ISM 直讀 read-level MM/ML | Methods 義務：必標 dorado model 版本字串（v0.7 vs v0.8 輸出可變）|
| **DeepSomatic / DeepVariant / Severus** | variant/SV calling vs ISM 甲基結構檢定 | 全 ONT-capable（T1 禁二代）；全甲基零命中（佐證 filter-DEAD）。🔴 landscape 05 寫 deepsomatic「CASTLE 含 COLO829」**過度宣稱** —— repo 只有 HCC1395，CASTLE/COLO829 在獨立 landing |
| **LongPhase-base / WhatsHap** | germline 共相位上游 vs ISM somatic haplotagging | 🔴 LongPhase base modcall 內建 germline-haplotype 級 ASM-as-phasing-marker（`ModCallParsingBam.cpp:633-704`）→「甲基有 allelic 結構/可輔助 phasing」**非 ISM 首創**；引用須明引 DOI 10.1093/bioinformatics/btac058 |
| **軸1 基因型 SR**（PyClone/PyClone-VI/PhyloWGS/Canopy/SciClone/DPClust/CliP）+ Smallwood/Capper/Fu | bulk 短讀 DNA-seq VAF/CCF 重建骨幹（源碼層 grep 零甲基）vs ISM read-level 甲基 | 不同 regime，但 9 源同時**佐證**「reconstruction=基因骨幹責任」（見桶⑤）|
| **methClone** (Li 2014) | 縱向兩階段 epiallele 組成熵差 ΔS vs ISM 橫向單時點距離結構檢定 | 共用 read 內 phased-CpG substrate 但統計目標正交 |

### 桶④ 矛盾/限制（framing tension —— 對抗稽核補入 2 個漏歸的高 tier 源）

> **這是你最在意的「是否矛盾」**。重申：**沒有實驗結果衝突**，全部是**用詞/新穎邊界**的 framing tension。每條附「限制什麼 + 論文如何處理」。
> 🔴 **對抗稽核發現原合成漏了 MethylTree + EVOFLUx**（兩個 ⚡-marked、Nature/Nat Methods tier A 源被只放進「佐證」），已補入為**雙角色**。

| 嚴重度 | 源 | 限制 ISM 什麼 | 論文如何處理 |
|:--:|---|---|---|
| 🔴 HIGH | **Sgootr** (2023) | 「距離式甲基 lineage tree 重建」first principled method（真腫瘤單細胞）→ 縮限「距離式甲基重建是新穎 method-class」| 降為 prior art 並 cite；新穎收斂三件（Sgootr 無 phasing/無 normal-anchor/終態=tree 非結構檢定/somatic 軌放棄）。⚠ Sgootr 仍 cited_secondary（repo 未下載、PDF 未讀）→ **投稿前須升 primary** 確認它無 normal-anchor + 無 read-level 結構檢定（A2/A3 正靠這缺口）|
| 🔴 HIGH | **Gaiti2019** (CLL) | 單細胞 regime 直接用「reconstruction」宣稱甲基能重建 lineage + genetic subclone（SF3B1）準確映射純 epimutation clade（P=7.4e-9）| 明界定 regime 差（單細胞+監督 ML phylogeny vs ISM bulk 無監督距離）；次要**佐證** R2（提供 ISM 缺的單細胞正交 ground truth）|
| 🔴 HIGH | **Epiclomal** (2020) | **最尖對撞**：sc-WGBS 甲基可解析次選殖模式（epiclones usable，甚至 transcend CN）→ 挑戰「subclone 層甲基弱、不可 read 救援」| scope 化為「bulk long-read 單 region read×CpG、為 variant 判別的 regime 下」subclone 甲基不 usable。⚠ 正文效能數字（V-measure/CN-transcend）未讀 PDF，引用前 fetch PMC7546467 |
| 🔴 HIGH | **PoreMeth2** (2025) | 標題「decoding the evolution of methylome with nanopore」+ AML 演化 → 壓縮「長讀甲基演化」泛詞（**同平台**，T1 不可當差異）| 🔴 必寫 Discussion：它是 per-CpG Δβ+ΔS 雙變量分段（bulk，無 phasing/無 read×read 距離/無 PERMANOVA）→ 正交。🔴 landscape 05 標 COMPLEMENTARY **低估**，建議加 ⚡ |
| 🔴 HIGH | **CAMDAC/TRACERx** (2025) | (1) 59 病人真 cohort → 限制「6 cell line」泛化；(2) 證 region-level deconvolution「即足以」對齊甲基與基因組演化 → 施壓「read-level 必要性」| Discussion 明示 read-level 防守點（單分子 allele 連鎖 + 無監督結構 + somatic read-level cis）。🔴 它的「allele-specific」**實指 copy number（ASCAT.m）非甲基** → 不可混引/借它背書 ISM「ASM」|
| 🔴 HIGH | **HapBridge + NANOME** (2025) | 「甲基=robust linkage 可修 germline phasing switch error」「haplotype-aware allele-specific 甲基化」已成熟/head-to-head → ISM **絕不可**宣稱「甲基輔助 phasing / haplotype-aware ASM」本身新穎 | 主動劃界「此 white-space 已被佔」；新穎收斂三件。⚠ HapBridge 無檢定=啟發式閾值（T2 禁用反推優勢）。🔴🔴 landscape 05 兩源仍標 cited/R- 但卡已升 source_verified |
| 🔴 HIGH（補入）| **MethylTree** (Nat Methods 2025) | 雙角色：單細胞~100% 甲基譜系重建可用 → 同 Gaiti/Epiclomal 的「reconstruction 用詞 / subclone 層甲基弱」壓力 | regime-bound 處理 + **必帶 Q-metric 口徑**（~100%=有監督 leaf-order 正確率，baseline Q≈0.05-0.18，**非訊號存在率**）。同時保留弱-subclone **佐證**角色（其成功校準 ISM 天花板是 modality/coverage-bounded）|
| 🔴 HIGH（補入）| **EVOFLUx** (Nature 2025) | 雙角色：標題級「甲基化 tracks cancer EVOLUTION」最強已發表 claim → 壓縮「長讀甲基演化/reconstruction」用詞 | bulk-population fCpG clock 層（無 read×read 距離/PERMANOVA/normal-cis）。同時**佐證**弱-subclone（自承 1,610/1,976 effectively-neutral、subclone 偵測限強選擇大 subclone）|
| 🟡 MED | **colorectal-growth** (2026) | 正面「甲基→腫瘤生長史 reconstruction」（fCpG 分子鐘 + ABC-SMC）| Discussion 尺度區隔（它走 model-based parametric，ISM 刻意不採）。摘要首段點明分工。⚠ cited_secondary（403）|
| 🟡 MED | **nanomethphase** (2021) | 摘要近似點名「甲基-ASM 可 leverage 於 LOH 癌症」→ 縮限「LOH×甲基新穎」措辭 | 新穎收緊三件。⚠ LOH 句僅摘要層未逐字核 → 對抗引用前升 V2。🔴 landscape 標 COMPLEMENTARY，主桶應 TENSION |
| 🟡 MED | **mcf7-subline** (2026) | 把 genomic+epigenomic divergence 連到 subline 分化 → 若被援引「甲基差異即可標記 lineage」逼近紅線 | 嚴守「subline（已分群固定）≠subclone（待重建）、allele 層≠subclone 層」。反證 ISM「6 cell line=跨模型 reproducibility 非泛化」。⚠ tier C 大量 UNVERIFIED（403）|
| 🟡 MED | **BitPhylogeny** (2015) | 把「reconstruction」定義為含 clone 數+演化樹 → 限制標題用 'reconstruction' 卻不出樹 | Discussion 明界定 ISM reconstruction=read→subclone 指派/分群結構，不出 lineage tree。⚠ 須承認它有 colon cancer 甲基建樹例（per-site 二元 genotype 群體層，非 read-level）—— 防守點是 granularity + 不建樹，非「甲基無關」|
| 🟡 MED（升）| **Lee2025** melanoma | 用 SNV「reconstruct」演化樹 + bulk DMR 標記 aggressive subclone → 用詞撞車 | Related Works 明區分 reconstruction 來源（ISM=somatic haplotag vs Lee=單細胞 ground-truth+SNV）。其 r=0.842 甲基鏡像 SNV + 樹由 SNV 建 → 強**佐證** A-framing |
| 🟡 MED | **ClairS-TO** (2025) | Verdict 模組已做 variant-AF subclonal 分類 → 縮限「subclonal 用詞原創性」+ 是 single-pipeline 自參照 ⭐3 上限根因 | Discussion 聲明 ISM 消費（非超越）caller 的 subclonal call；禁裸寫「首觸 subclonal」。🔴 landscape 05 未標此 T3 tension |
| ⚪ LOW | **Smallwood2014**（補入）| 單細胞甲基異質性**金標準** → underwrite Epiclomal/MethylTree 同類挑戰 | regime-bound：scBS-seq 異質性不轉移到 bulk long-read 單 region variant 判別；Discussion regime-bounding 引它為上游金標準 |
| ⚪ LOW | **LRSomatic** (2026) | abstract 宣稱「Fiber-seq haplotype-specific chromatin accessibility」全文未讀（403）→ 若含 read-level epigenetic clustering 會縮限 A1/A3 | 升 verification 時第一個查。目前證據未顯示 subclone/lineage 重建，暫不縮限但不可定論 |

### 桶⑤ 可佐證（外部數據反向支撐 ISM 內部結論）

> **這是 ISM 最厚的外部支撐**。每條附 strength + 對抗稽核校正。

| ISM claim | strength | 主要佐證源（已校正）|
|---|:--:|---|
| **A-framing：reconstruction=基因/somatic-haplotag 骨幹，甲基非驅動** | **strong** | 軸1 九源源碼層 grep **零甲基**（PhyloWGS/PyClone(-VI)/Canopy/SciClone/DPClust/CliP）+ Tarabichi 領域 review + longphase-S 骨幹本體 somatic 路徑零甲基 + Lee2025（r=0.842 甲基鏡像 SNV、樹由 SNV 建）|
| **R2：甲基=germline-haplotype 層級且為主導可回收訊號** | **strong** | cvlr（8 imprinted gene Δm P<10⁻³；🔴**H19 例外** phasing P=1，7/8 顯著 →「甲基子群≠必然=haplotype」**反而強化 ISM 誠實**）+ ASMS（Table1 Fisher 16/20 ICR）+ O'Neill（4.46M aDMR, tier A）+ longhap（Fig3C 強不對稱 DMS、甲基橋接 18.2% phase block；🔴標 **bioRxiv preprint tier C**）+ MethPhaser/sakamoto/HapBridge/NANOME/gabbutt-fcpg + Landau（structure≠disorder）|
| **白地：cancer-subclone 甲基-phasing 真實白地** | **strong** | keystone=MethPhaser（tier A Nat Commun，**逐字**把 cancer 列 future work，獨立第三方）+ Wakhan（源碼+全文零甲基雙證）+ HiCancer（Hi-C 跨模態第二例）+ longhap（最強 germline phaser 仍 germline-only）。⚠ longphase-S 須誠實標**非獨立第三方**（同實驗室上游）|
| **弱-subclone 紅線：subclone 層甲基弱-存在性、不可 read 救援、⭐3 封頂** | **strong** | EVOFLUx（Nature；1,610/1,976 effectively-neutral）+ MethylTree（Nat Methods；~100% 校準天花板=modality-bounded，🔴帶 Q-metric 口徑）+ Wakhan（polyclonal 需多樣本）+ Salcedo benchmark（演算法吃 19-35% 變異、無單一最佳）+ PyClone README（單樣本性能差）+ longphase-S（HP3 F1=0.746 量化最弱邊界）|
| **甲基非 variant TP/FP filter（filter-DEAD 正交）** | moderate-strong | 軸6 caller 全甲基零命中（DeepSomatic source-verified repo grep=0 / ClairS-TO / Severus / DeepVariant）|
| **CN-confound 必控（HP-axis held-CN）+ normal-baseline 必要** | moderate | Martin-Trujillo（印記 37 DMR：🔴**CN-independent 真表觀變化只 8.1-17.6%/癌種**，即 CN 解釋互補的 ~82-92%；分母限 imprinted DMR）+ Canopy/SciClone/PyClone-VI（顯式建模 CN）+ cancer-genome-standards（matched-normal depth 是 somatic 準確關鍵，設計理據層佐證非同 metric）|
| **癌症 ASM 真實存在且升高（ISM 6/6 excess-over-null 生物學先驗）** | moderate | Do2020（cancer ASM 5-9×、hypo-dominant）+ MCF7 + DAMEfinder。🔴 **Do2020 5-9× 是 cross-individual bisulfite 規模 vs ISM 0.34% within-sample somatic-anchored**，不同分母/平台/單位**絕不可並列**；BRCA2 hypo≠canonical Do 口徑（3-不對齊）|

---

## L2.5 — 對抗稽核結果（4 lens，已套用修正）

| lens | verdict | 關鍵發現（已修正進上表）|
|---|:--:|---|
| **reviewer-trap** | NEEDS_WORK | 🔴 原合成把 qFDRP 算進「對手皆有顯著性檢定」（事實錯：qFDRP 無檢定）→ 已改精確口徑（桶①A1 caveat + 桶② qFDRP）。其餘 T1/T2/T3 處理逐條 PASS |
| **advantage-defensibility** | NEEDS_WORK | 🔴 A2 cis-test 過度宣稱（BRCA2 塌陷、chr17 才乾淨）→ 降級。🔴 A2 歸屬錯誤（三角在 Python 非 C++）→ 已分。🔴 三件皆 characterization-only 非 utility（filter-DEAD）→ 已標。A1 須 demo structure-vs-magnitude divergence。Sgootr 等支撐卡 cited_secondary 須升 primary |
| **tension-honesty** | NEEDS_WORK | 🔴 MethylTree + EVOFLUx（2/10 ⚡源）漏歸 tension → 已補入雙角色。Smallwood 升 tension。Lee2025 severity low→med |
| **corroboration-validity** | **PASS** | strength 標籤校準良好未灌水。微修：H19 7/8、Martin-Trujillo 8.1-17.6%、longhap preprint tag、landau「explicitly excluded」載重句投稿前須核 axis4 卡 |

---

## L3 — 論文敘述方式評估

### 建議整體定位（A-framing 定稿句）

> ISM 是在 **ONT 長讀癌症 tumor-normal paired** 資料上，以 **somatic-haplotagging 為 subclonal reconstruction 骨幹**、並對 read×CpG 甲基矩陣做**有界 characterization** 的整合引擎。甲基面的可防守貢獻**嚴格限於三件正交於既有工具的特定組合**：①把 read-to-read 距離矩陣升格為分群引擎並用 PERMANOVA pseudo-F 檢定 read 在距離空間的*結構分離*；②以 matched-normal 為錨點的 somatic cis/drift 分解層；③somatic-subclone structure characterization 目標。甲基定位為 **germline-haplotype 層級強、within-haplotype subclone 層級弱-存在性、不可 read 救援、非 variant TP/FP filter（DEAD）**。

### DO（該做）
- Related Works 把已被佔的地盤全列 **prior art 並 cite**：甲基-輔助 phasing（MethPhaser/HapBridge/NANOME/LongPhase modcall）、距離式甲基重建（Sgootr）、單細胞甲基 lineage（Gaiti/Epiclomal/MethylTree）、haplotype-aware ASM（NanoMethPhase）、甲基演化（EVOFLUx/PoreMeth2/colorectal）。
- 甲基面一律 **'characterization'/'structure test'**；reconstruction 用詞**只掛 somatic-haplotagging 骨幹**，明界定 = read→subclone 指派非 lineage tree。
- 對手對照差異放「**不同檢定對象**（距離結構 pseudo-F vs cluster 甲基差量/disorder magnitude）」+「normal-anchor + somatic 框架」，**絕不**放平台世代或「缺顯著性」。
- Methods cite qFDRP 為 NHD kernel 先例、cite DAMEfinder/cvlr/ASMS 為 read-level ASM 概念先例、cite PERMANOVA（Anderson 2001）、reconstruction 骨幹明引 LongPhase DOI。
- 三件優勢一律標「**characterization-layer capability, orthogonal to variant TP/FP utility（filter-DEAD, AUC≈0.50）**」。
- 誠實標 cell-line scope（6 cell line=跨模型 reproducibility 非 cohort）、single-pipeline ⭐3 封頂。
- 修 `PerCpgAsm.cpp:6` 方法用詞（DAMEfinder 非 Fisher）。

### DON'T（不該做）
- 不裸寫「用甲基重建 subclone」/「距離式甲基重建是新穎 method-class」/「甲基有 allelic 結構/可輔助 phasing 本身新穎」。
- 不寫「對手用二代定序」/「對手缺顯著性檢定」/「0 真 CONFLICT」/「甲基不能定義 subclone」。
- 不把 CAMDAC「allele-specific（=copy number）」當同類 ASM 混引；不把 MCF7 subline 差異當 subclone 重建。
- 不把 ISM 0.34% / 6 cell line 與 cohort-scale 數字（Do2020 5-9×、O'Neill 4.46M、CAMDAC 59 病人）並列。
- 不超出 pre-reg 宣稱 lineage tree / 祖裔推斷 / 臨床級甲基分類（Capper scope 不可借）。
- 三件「advantage」不可漂移成 utility/效能 claim。

### 第二章（文獻探討）6 節骨架（對應 6 軸，每節「領域定位→與 ISM regime/目標差→ISM 邊際位置」）
- **§2.1 亞克隆重建骨幹**（軸1）— 確立 reconstruction=基因骨幹責任、甲基非驅動、single-pipeline ⭐3 紀律；BitPhylogeny 界定 reconstruction 用詞範圍。
- **§2.2 長讀定相/haplotag**（軸2）— 骨幹分工（longphase-S）、甲基-輔助 phasing 已被佔（MethPhaser/NANOME/HapBridge）、cancer-subclone 甲基白地真實（Wakhan/HiCancer 零甲基雙證）。
- **§2.3 癌症 ASM**（軸3）— ASM 廣存須控 CN/purity（Martin-Trujillo/CAMDAC 入場券）、用詞精準切割（allele-specific CN≠甲基）。
- **§2.4 甲基追演化/lineage**（軸4）— **最密集 T3 戰場**：Sgootr/Gaiti/Epiclomal/MethylTree/EVOFLUx 降 prior art、新穎收斂三件、反向校準弱-subclone。
- **§2.5 read-level 甲基分群/距離**（軸5）— **三件可防守差異的核心論證節**：cvlr/ASMS 同軸最近對手 + qFDRP NHD kernel 先例 + MethylBERT supervised 對照。
- **§2.6 somatic callers**（軸6）— 定位 ISM 消費 caller、甲基非 filter、第二 pipeline（LRSomatic）破自參照。

### 誠實天花板（5 條，必寫 Limitation）
1. 單樣本 single-pipeline 自參照 → 信心封頂 **⭐3**（TP/FP 真值取自 tumor-only caller 非獨立 orthogonal ground truth）；6 cell line 資產可往 ⭐4 但須跨樣本 cross-validation。
2. 6 cell line = **跨模型 reproducibility，非病人 cohort 泛化**；purity≈100% 規避 deconvolution 是 scope 不同非優越。
3. **無 lineage tree 輸出** —— reconstruction 用詞限 read→subclone 指派/分群結構。
4. subclone 層甲基僅**弱-存在性**（modality/coverage-bounded，非方法失敗也非普世定律：MethylTree 證需單細胞+全基因組+ground-truth 高資訊 regime 才達~100%），不可 read 救援、非 variant filter（DEAD）。
5. 大量軸3/軸4 **C-tier 源為 cited_secondary**，關鍵數字（r=0.842 / Epiclomal 效能 / LRSomatic Fiber-seq 段落 / landau exclusion 句）未親讀全文 PDF，對外正式引用/對抗前須升 V2 親讀。

---

## §趨同驗證 — 外部驗證 × ISM 實驗 × 一致性 × 獨特性

> **這節回答**：每個 ISM 目標，問 ①有人獨立驗證過嗎（有結論?）→ ②我們也實驗了嗎 → ③結果一致嗎 → ④我們方法更獨特嗎。
> **內部數字皆 2026-06-13/14 grep-verified**（`knowledge/11_external_literature/` + ISM spec L3 + methyl_phasing VERIFIED_RESULTS），非記憶。

| ISM 目標 | ①有人驗證？（who + 結論 + 獨立性）| ②我們實驗？（result + tier + 來源）| ③一致？ | ④更獨特？ |
|---|---|---|:--:|---|
| **G1 LOH-constrained phasing**（subclone 骨幹的甲基-assist 面）| ✅ Wakhan/HiCancer/Sakamoto 用 LOH/imbalance phasing，≥3 癌症論文背書「LOH-as-phasing 可行」（用 CNA/SV/Hi-C **非甲基**）| ✅ NG=2 inner same-hap **93–99% 6/6**（HCC1395 93.2%…DORADO 99.0%）；7/7 **W=28 p=0.0078**（median gap 0.345）；`LabelTest.cpp:265-302`；paired 負控 B3 NS（⭐3 B+）| 🟢 原理 CONFIRMED | ✅ **甲基貢獻 distinct（他們未測甲基）**。⚠ HD-1 循環對照未跑 |
| **G2 ASM 在癌症存在且升高** | ✅ Do2020（5-9×、hypo-dominant）/Martin-Trujillo/cvlr/ASMS(16/20 ICR) → ASM 真實廣存 | ✅ **6/6 excess-over-null mean 0.168**（HCC1395 0.2413…COLO829 0.101，3 癌種）；0.34% within-sample；BRCA2/ZAR1L ⭐3 | 🟡 **caliber 差**：方向一致；Do2020 5-9× cross-individual bisulfite vs 0.34% within-sample somatic **不同分母不可並列**；BRCA2 hypo≠canonical hyper | ✅ within-sample somatic-controlled（HP-axis held-CN）|
| **G3 甲基=germline-haplotype 層級（R2）**| ✅ **強、多平台**：cvlr(P<10⁻³)/ASMS(16/20 ICR)/O'Neill(4.46M aDMR)/longhap(Fig3C) | ✅ BRCA2 HP-axis **Δβ=−0.122** | 🟢 strong（H19 P=1 例外，**7/8**：甲基子群≠必然=haplotype，**反強化誠實**）| ⚪ **部分**：無監督 read-level+normal-anchored，但「甲基=haplotype 層」**非我們獨有**（cvlr/ASMS 先做）|
| **G4 cancer-subclone 甲基-phasing 白地** | ✅（白地真實）MethPhaser **逐字**列 cancer=future；Wakhan/HiCancer 零甲基；longhap germline-only | ✅ methyl-rescue pilot：T1 V6 全基因組 **0.885**（0.926 僅可分子集、LOH 未測）；T2 **0.90 僅在有真值 1-1/2-1**（H3 未證、15-18% 可指派 **[OVERSTATED]**）；**T3 DEAD**（0.50, 亞群−噪音 p=0.92）| 🟢 白地確認（沒人做）| ✅ 我們是**填白地者**。⚠ 較難版本 T3 不 work；longphase-S **非獨立**（同實驗室上游）；無 LOH 內 phase-truth |
| **G5 弱-subclone 紅線**（subclone 層甲基弱）| ✅（**反向佐證**）EVOFLUx(1610/1976 effectively-neutral)/MethylTree(~100% 但 Q=有監督, baseline 0.05-0.18, modality-bounded)/Wakhan(polyclonal 需多樣本)| ✅ subclone 甲基 undetermined/弱；HPFineNGroups pipeline-dependent；methyl filter ΔF1 marginal | 🟢 外部說 subclone 罕見/需高資訊 regime，我們 bulk read-level 找到弱訊號 | ⚪ 誠實標 bulk long-read regime 天花板（別人不在此 regime）|
| **G6 甲基非 variant TP/FP filter（filter-DEAD）** ⭐最強趨同+獨特 | ✅（**confirmed-by-absence + ML 文獻**）無 production caller 用甲基當 filter（ClairS-TO/VarNet/AIVariant 用 alignment/PoN/CN）；Kapoor leakage taxonomy(L3.2/L1.3)/Soneson **直接命名**我們 LOSO 循環 | ✅ **詳盡**：LOSO 100% circularity（in-dist +0.0224 vs held-out −0.0001、HCC1954 −0.377）；**\|Δβ\| AUC=0.5049**(perm p=0.496)；strong-ASM FP 富集 5× **OR=0.194 p=1.8e-28**；TP 3.95%>FP 1.07% sens~4% | 🟢 外部 CONFIRM 負結果 | ✅ **我們唯一實測**（多數沒測）；strong-ASM-in-FP anti-discriminative = **NOVEL 外部未 addressed**。⚠ filter-negative 單 pipeline（全 ClairS-TO）|
| **G7 reconstruction=基因骨幹、甲基非驅動（A-framing）**| ✅ 軸1 九源全零甲基重建；Lee2025 樹由 SNV 建（r=0.842 甲基鏡像）| ✅ A-framing；reconstruction 走 somatic haplotagging | 🟢 | ⚪ 在基因骨幹上**加**甲基 characterization（別人完全不做甲基），但甲基 characterization-only |

### 🔴 核心洞察：趨同驗證 與 獨特性 往往互斥

| | 外部強驗證 | 我們獨特 | 判定 |
|---|:--:|:--:|---|
| ASM 存在(G2)、甲基=haplotype(G3) | ✅ 強 | ❌ prior art | 被驗證 = 但**不獨特**（cvlr/ASMS/Do2020 先做）|
| normal-baseline cis-test、距離矩陣結構**組合** | ❌ 無人獨立測 | ✅ 真獨特 | 獨特 = 但**外部尚未獨立驗證** |
| **filter-DEAD(G6)** | ✅ 外部 confirm 負 | ✅ 唯一實測+novel anti-signal | 🟢 **兩者皆有 → 最強材料** |
| **白地 phasing(G4/G1)** | ✅ niche 真實 | ✅ 我們填 | 🟢 兩者皆有，⚠ 但自參照/HD-1/無 phase-truth 限制 |

**一句裁決**：真正獨特的只有 3 件（normal-baseline somatic cis-test〔唯一無人有〕、read×read 距離矩陣結構檢定的**組合**、somatic-subclone 目標）—— 但這 3 件**正好是外部還沒獨立驗證的**；反之被外部強佐證的（ASM 存在、甲基=haplotype）**恰恰非我們獨有**。**唯二「既被外部佐證又是我們獨特貢獻」= ① filter-DEAD 負結果 ② 甲基-assist LOH phasing 白地**。

---

## §ISM 的獨特改進 — 可防守貢獻的敘述

> 把「我們相對既有方法的獨特推進」分三層敘述。每件標 tier + **characterization-only**（非 variant utility）+ 誠實 caveat + 外部對照（何以獨特）。

### 第 1 層 — 唯一無人有的方法元件：normal-baseline somatic cis-test

**是什麼**：以 **matched-normal per-CpG 基線**扣除（residual = raw − normal mean，C++ `NormalBaseline.hpp:38-67`），再用三角分解 `d_somatic=C−B / d_cis=C−A / d_drift=B−A`（Python `scripts/34_*.py:108-152`）+ cis-candidate 判準（`p_cis<0.05 AND |d_cis|>1.8×|d_drift|+0.02`）分辨**真 cis-driven ASM vs germline 漂移**。

**何以獨特**：**所有 survey 的 germline 甲基工具（cvlr/ASMS/MethPhaser/NanoMethPhase/qFDRP）與 modkit 皆無 normal anchor** —— 這是 ISM **真正無人佔的一件**。matched-normal 充足 depth 是 somatic 準確度關鍵，從 variant-calling 側獨立印證（cancer-genome-standards-2026，設計理據層）。

**誠實 caveat**：(a) 🔴 歸屬：C++ engine **只做 residual**；三角在 **Python 分析層**（hardcoded path/無 regression test）→ 論文禁稱 compiled-engine feature。(b) 🔴 可靠性 = **分解層本身**，非穩健單位點 cis 定量：旗艦 BRCA2 表面 d_cis=−0.142 約 80% 是 copy/subclone confound，真 focal **d_within=−0.023（perm p=0.024 marginal、% 不 robust）**；乾淨 exemplar 是 **chr17/TBC1D16 d_within=0.142（p=0.001）**。

### 第 2 層 — 獨特的方法組合（非單一元件，是 composition）

**是什麼**：把 read×read 距離矩陣（預設 NHD）**升格為分群引擎**（UPGMA/Ward 階層分群）+ 對 HP/somatic 標籤做 **PERMANOVA pseudo-F 結構檢定**（999 perm，`StructureTest.cpp:141-205`），檢定 read 在距離空間是否**結構分離**（非 cluster 甲基差量、非 within-read disorder）。

**何以獨特**：cvlr/ASMS 測 cluster 間 median|Δmeth|；qFDRP 把同一 NHD kernel **塌成 per-CpG scalar**（只量 disorder magnitude）。ISM **保留全距離矩陣 → 分群 → 測 structure**，在 **somatic-haplotag conditioning + normal-anchor + somatic-subclone 目標**下，是無人佔的組合。

**誠實 caveat**：🔴 **每個元件都是現成的**（NHD=qFDRP kernel；階層分群=教科書；PERMANOVA=Anderson 2001 生態學標準；距離式甲基分群=Sgootr 已有）→ 防守靠的是**組合 + 條件化**，**非任一單元件新穎**。⚠ structure-vs-magnitude 區分**須實證**：投稿前須 demo ≥1 個 pseudo-F 與 median|Δmeth| 不一致的位點，否則只能框成「在保留矩陣上測不同統計量」而非「判別優勢」。Methods 須 cite qFDRP(NHD 先例)+PERMANOVA(Anderson 2001)+DAMEfinder/cvlr/ASMS(ASM 概念先例)。

### 第 3 層 — 獨特的科學貢獻（findings，非方法新穎）

| 貢獻 | 是什麼 | 何以獨特 | tier / caveat |
|---|---|---|---|
| **filter-DEAD 嚴格證偽 + anti-discriminative novel** | 甲基當 variant TP/FP filter：LOSO 100% circularity / AUC=0.5049 / strong-ASM 在 FP 富集 5×(OR=0.194) | **我們唯一實際嚴格測過**（多數工具不試）；strong-ASM-in-FP anti-discriminative 是**外部完全未 addressed 的 novel 觀察**；外部 confirm-by-absence + Kapoor/Soneson 背書方法學 | ⭐2 L4 DEAD（locked，勿 reopen）。⚠ 單 pipeline（全 ClairS-TO，無第二 caller 重 call）|
| **within-sample somatic-controlled ASM 量化** | HP-axis held-const-CN 下量 somatic ASM（0.34% 稀有、無方向偏好；6/6 excess-over-null）；ASM×CN ρ=−0.055 證非 copy-driven | 對手是 cross-individual cohort（Do2020）或 germline（cvlr/ASMS）—— **無人在 within-sample somatic-controlled 口徑量過** | ⭐3 單樣本主導。⚠ 與 cohort-scale 數字不可並列；private 0/38 underpowered |
| **cancer 甲基-assist LOH phasing 白地填補** | 在 germline-het-absent 的 LOH 區用甲基 tag unphased read（LOH-constrained phasing Grade A）| 文獻 niche 真實（MethPhaser 逐字留 cancer future、germline 工具全不碰）；LOH-as-phasing 被 ≥3 癌症論文背書但走非甲基訊號 → 甲基貢獻 distinct | ⭐3 Grade B+。⚠ 無 vs MethPhaser/HapBridge switch-error/N50 head-to-head；無 LOH 內 orthogonal phase-truth |

### 🔴 誠實分級：哪些是「真改進」，哪些是「prior art 上的整合」

- **真獨特改進**（無人佔）：normal-baseline somatic cis-test（方法元件）+ filter-DEAD 嚴格證偽與 anti-discriminative 觀察（findings）+ within-sample somatic-controlled ASM 口徑（findings）。
- **整合式貢獻**（元件皆 prior art，價值在組合+條件化+把甲基 characterization 接上 somatic-haplotag 骨幹）：距離矩陣結構檢定 composition、A-framing 骨幹分工。
- **被外部佐證但非我們獨有**（只能當 corroboration 不能當 novelty）：ASM 存在、甲基=germline-haplotype 層級。
- **鐵則**：以上全部是 **characterization 層貢獻，永不可寫成改善 variant calling/reconstruction 效能的 utility 優勢**（TP/FP 判別 filter-DEAD AUC≈0.50）。

### 方法精進方向（LEARN，從對照中浮現的可再改進處）

1. **Fisher → beta-binomial + dispersion shrinkage**（修 read 非獨立膨脹；有 pyclone/evoflux/cvlr 發表級實作可對接，最高 ROI）。
2. **停 5mC/5hmC max-collapse**（`cis_asm_core.py:31-32`，與 MSA dup-bug 同源）。
3. **客觀選 k**：cvlr 用 BIC/AIC（ISM 目前靠 PERMANOVA 非 model selection）。
4. **破自參照**：第二 caller（DeepSomatic/Mutect2）重 call 驗 filter-DEAD；LOH 內 orthogonal phase-truth（trio/Hi-C/strand-seq）驗 phasing → ⭐3→⭐4 路徑。
5. **demo structure-vs-magnitude divergence 位點**（坐實第 2 層的「不同檢定對象」非空洞）。
6. **Sgootr 升 primary**（確認其無 normal-anchor + 無 read-level 結構檢定，鞏固第 1/2 層獨特性的 prior-art 缺口）。
7. **修 `PerCpgAsm.cpp:6` 方法用詞**（DAMEfinder 原生非 Fisher）。

---

## 附：external_validation 庫維護發現（其他 session owner，僅標不改）

逐源複驗時 agent 發現的 landscape SoT 漂移（影響可尋性，非結論）：
- **verification drift**：`bitphylogeny` / `hapbridge` / `nanome` 卡已升 source_verified（repo cloned 2026-06-13）但 `_landscape/05` 索引列仍標 cited/R-。
- **過度宣稱**：`05` 寫 `deepsomatic`「CASTLE 含 COLO829」→ repo 實際只有 HCC1395（COLO829/CASTLE 在獨立 landing），據此抓 COLO829 真值會撲空（對齊 `06` F-3）。
- **tension 低估**：`nanomethphase`（LOH×甲基）、`clairs-to`（subclonal 用詞）、`poremeth2`（同平台甲基演化）、`hapbridge/nanome` 在 `05` 標 COMPLEMENTARY/DIFFERENT-VIEW，主桶應含 framing-tension（建議加 ⚡）。

---

## Provenance
- 逐源關係 + 五桶：workflow `wf_b40031c1-517`（6 軸 agent 親讀 `external_validation/axis*/<slug>/CONTEXT.md` + landscape 05/06 + ISM spec）。
- ISM 結果數字（0.34% / Δβ=−0.122 / d_within=−0.023 perm p=0.024 / chr17 d_within=0.142 p=0.001 / AUC=0.505 / 6/6 excess）：引自 `01_ism_method_spec_from_source.md` L3 表（各標 `knowledge/11_external_literature/` validated 來源）。
- reviewer-trap 口徑：`06_source_read_corrections.md` F-1/F-1b（源碼親讀）。
- framing-tension reframe：`20260612_external_validation_literature_gap_audit_01.md` §6.5。
- 對抗稽核：4 lens（reviewer-trap / advantage-defensibility / tension-honesty / corroboration-validity），3 NEEDS_WORK + 1 PASS，修正逐項套用。
- 主 agent 獨立複核：BRCA2/chr17/0.34%/NormalBaseline.hpp-vs-scripts/34 歸屬皆對照自讀的 ISM spec 確認（§13.7）。
- **§趨同驗證 + §ISM 獨特改進（2026-06-14 新增）**：內部數字 2026-06-13/14 grep-verified —— 6/6 excess mean 0.168 / LOH same-hap 93–99% 6/6 W=28 p=0.0078 / AUC=0.5049 / OR=0.194 p=1.8e-28 / ASM×CN ρ=−0.055 引自 `knowledge/11_external_literature/{05,06,07,08,10}`；methyl-rescue **T1 0.885（非 0.926，後者僅可分子集）/ T2 0.90 僅在有真值[OVERSTATED] / T3 DEAD** 引自 `docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/`（plot_rigor_corrected.py + VERIFIED_RESULTS）—— 🔴 此處與 memory `project_methyl_phasing_assist_line` **本體一致**（本體已誠實標 T1 0.885/T2 OVERSTATED/T3 DEAD）；已順手校正 MEMORY.md **索引行**的樂觀縮寫（原寫 T1✅0.926/T2✅0.90）。
