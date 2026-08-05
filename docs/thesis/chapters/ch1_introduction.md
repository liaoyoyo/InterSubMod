<!--
建立時間: 2026-06-15（v3 對齊 06-14 Intro 主軸定稿 v1：ISM 首個整合 C1/C2/C3 + 旗艦 BRCA2/chr2 + 甲基↔突變因果歸文獻 + characterize-conditioned 措辭）
報告類型: 碩論正文草稿 — 第一章 緒論
狀態: draft v3（取代 v2）；對齊 docs/paper_focus/02_paper_framework/20260614_Intro主軸定稿_thesis_statement_01.md；citation {{cite:...}} 待 /citation-verification
data_sources:
  - docs/paper_focus/02_paper_framework/20260614_Intro主軸定稿_thesis_statement_01.md (用戶逐句確認主軸句 + C1/C2/C3 + 旗艦 + citation 清單)
  - docs/reports/research_landscape/20260612_external_validation_literature_gap_audit_01.md (framing 凍結裁決)
provenance_note: 主軸句逐句對應 06-14 定稿（用戶確認）；BRCA2 Δβ=−0.122/d_within=−0.023、chr2 F=10.6/V=0.67/V=0.97 grep-verified；甲基→突變因果機制 = 既有文獻背景 {{cite:SBS1/5mC-deamination 待 citation-verify，庫未全文驗勿憑記憶引}}。措辭守 framing 凍結裁決。
-->

# 第一章　緒論

## 1.1 研究背景與動機

腫瘤於演化過程中分化出多個帶有不同體細胞變異之亞克隆（subclone），形成腫瘤內異質性（intratumor heterogeneity），與惡性進展、轉移及抗藥性密切相關。重建腫瘤之亞克隆結構與演化（subclonal reconstruction）為癌症基因體學之核心課題，其主流做法以體細胞變異為骨幹——經變異呼叫、拷貝數重建、由變異等位基因頻率推得癌細胞分率、分群並建立譜系{{cite:Tarabichi 2021 Nat Methods 待 citation-verification}}{{cite:PyClone-VI 2020 待 citation-verification}}；正式之重建需 clone 數、譜系樹、拷貝數／純度去捲積，或正交之 ground truth{{cite:BitPhylogeny 2015 待 citation-verification}}。

Oxford Nanopore（ONT）長讀定序能於**同一條 read** 同時取得體細胞變異、germline 單倍型與原生 5-methylcytosine（5mCG）甲基化。在此之上，somatic haplotagging（如 longphase-S）可將腫瘤 read 依其證據指派至 germline 單倍型（HP1/HP2）與體細胞子單倍型（HP1-1/HP2-1）{{cite:LongPhase Lin 2022 待 citation-verification}}；而以甲基化延伸相位之方法亦已成熟，惟多侷限於 **germline 層級**且該領域已相當擁擠{{cite:MethPhaser 2024 / NanoMethPhase 2021 / HapBridge / NANOME 待 citation-verification}}。

## 1.2 研究缺口

在 germline 甲基-相位已擁擠、而 bulk 之亞克隆重建工具多不處理 read 層級甲基之背景下，存在一個明確而未被覆蓋之利基（niche）：**在僅有單樣本的癌症資料上，以 somatic haplotag 為條件（conditioned on somatic haplotags），於 read 層級對「甲基與體細胞變異之聯合結構」進行 normal-anchored 之 characterization**。具體而言：（一）read 層級之甲基-變異聯合結構未在 somatic haplotag 框架下被量化；（二）甲基差異究屬真 focal cis 抑或拷貝數／subclone 背景，未以 matched-normal 為錨點於 read 層級分解；（三）「以甲基重建 subclone」之過寬主張缺乏 clone 數／譜系樹／正交 ground truth 之支撐，故本研究於完成大規模重建 demo 與正交驗證前，正文一律採較保守之措辭——**「characterize methylation-associated subclonal structure conditioned on somatic haplotags」**，而非「reconstruct subclones using methylation」。

## 1.3 研究目的、貢獻與方法定位

本研究提出並評估 **InterSubMod（ISM）**——**首個**在 **ONT 長讀 tumor/normal paired** 資料上、整合下列三者於 read 層級之軟體：**① read 層次甲基結構檢定（read×read 距離 → 分群 → PERMANOVA）+ ② matched-normal-anchored somatic cis-test（殘差 + d_cis／d_drift／d_within 三角分解）+ ③ somatic-subclone 目標**。需明確界定：此三項之**個別**元件皆有前例（如 qFDRP 之距離核、Do2020 之 de-novo ASM 子分析亦用 within-patient matched-normal、Sgootr 之距離式甲基重建），**本研究之新穎性在於三者於 tumor/normal read-level 之整合**，而非距離式或 read-level 甲基本身{{cite:qFDRP 2020 / Do2020 / Sgootr 2023 待 citation-verification}}。

ISM 以 somatic haplotagging 之 HP tag 與 somatic 變異標記為框架，於每個 somatic 位點鄰域量化等位基因特異性甲基化（ASM）結構，輸出三類資訊：（1）與 somatic 變異高度關聯、於 read×read 距離空間呈顯著結構分離之甲基位點（PERMANOVA）；（2）tumor 內 haplotype 間／tumor-normal 間有顯著差異之甲基位點，並以 matched-normal 為錨點**分解**此差異屬真 focal cis 抑或 copy／subclone 背景（**BRCA2 promoter 為旗艦示範**）；（3）甲基差異與哪一軸（HP／ALLELE／copy）相關之分解。進一步，以「多個 somatic 變異 + 鄰域甲基結構」之組合對單一 region 之 read 做 subclone 結構 characterization，並以 **chr2:18M 案例示範**甲基如何與突變共同標示同一子克隆事件、據以推導時序假說。

需特別界定甲基與突變之關係：甲基化造成體細胞突變之因果機制（如 5mCG 去胺化導致 C>T、對應 COSMIC SBS1 特徵）為**既有文獻背景**{{cite:Pfeifer 2006 / COSMIC SBS1 / Alexandrov 2020 待 citation-verification（庫內未全文驗，勿憑記憶引）}}；本研究**不宣稱**建立因果或先後順序，ISM 之貢獻為量化 somatic 變異與甲基位點之**關聯強度**。又因甲基→突變與突變→局部 cis 甲基{{cite:Do2020 de-novo ASM / Zhang2019 待 citation-verification}}雙向機制並存，不可天真歸因單一方向。

## 1.4 誠實天花板與論文架構

本研究須明標下列邊界：（一）**case-demonstration**——chr2／BRCA2 為案例示範主軸概念，**非已驗證之自動化 subclone reconstruction**；不輸出譜系樹或演化參數推論，時序為數據支撐之假說而非正交確認；自動化分析流程與大規模驗證為 ISM 之延伸（後續工作）。（二）**單樣本／單一定相流程，證據層級 ⭐3**——TP/FP 真值取自 tumor-only caller（非獨立 ground truth）；六株 cell line 為**跨模型 reproducibility 而非病人 cohort 泛化**。（三）**甲基為 characterization 而非 utility**——不改善變異呼叫（甲基做變異真偽過濾之 AUC≈0.50，已證偽）。（四）**方向誠實**——BRCA2 之甲基差為 hypomethylation，與 canonical TSG promoter hyper-silencing 方向不同，不寫「沉默 BRCA2」；其 focal cis 因果於分解層為 marginal。

本論文其餘章節安排如下：第二章以五個面向（novelty／corroboration／white-space／prior-art／framing-tension）回顧相關研究並定位 ISM；第三章説明資料、somatic haplotagging 骨幹、ISM 之 read-level 方法（結構檢定／cis-test／subclone characterization）與對照設計；第四章先建立 de-confounded 之全基因組脈絡（4.1–4.4），再以兩旗艦案例示範主軸（BRCA2 之偵測與分解、chr2 之甲基-突變共同標示與時序假說），並説明所依附之相位骨幹與跨樣本復現；第五章討論機制詮釋、與文獻之對照及限制；第六章總結並提出後續工作。
