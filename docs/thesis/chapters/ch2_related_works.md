<!--
建立時間: 2026-06-15（v2 對齊 framing 凍結裁決 + 06-14 定稿 5 桶定位：novelty=組合 / 距離-read-level 甲基非新穎劃界 / framing tensions 非 0 conflict）
報告類型: 碩論正文草稿 — 第二章 文獻探討
狀態: draft v2（取代 v1）；citation {{cite:...}} 待 /citation-verification（庫內多數已 verified，見 06-14 定稿 L5）
data_sources:
  - docs/reports/research_landscape/20260612_external_validation_literature_gap_audit_01.md (framing 凍結裁決 + P0/P1 來源)
  - docs/paper_focus/02_paper_framework/20260614_Intro主軸定稿_thesis_statement_01.md (5 桶 + citation 清單)
  - docs/method_comparison/20260613_external_validation_library_index_01.md (49 源 CONTEXT 卡)
provenance_note: 每源 citation/DOI 見 external_validation CONTEXT 卡 §7 + 06-14 定稿 L5（庫內 verified）。守裁決：①ISM vs cvlr/ASMS/MethylBERT 禁寫「平台/缺顯著性檢定」②「無同 regime 衝突但有 framing tensions」非「0 conflict」③距離/read-level 甲基非新穎（Sgootr/Gaiti 前例）④6 cell line=跨模型 reproducibility 非 cohort。
-->

# 第二章　文獻探討

本章以五個面向定位 ISM 於業界地景：何謂正式之亞克隆重建及其對 framing 之約束（2.1）、ONT 長讀與 germline 甲基-相位之擁擠與開放白地（2.2）、癌症 ASM／cis／拷貝數混淆與 normal-anchored cis-test 之前例（2.3）、距離式與單細胞甲基重建之前例（劃界本研究之新穎性，2.4）、以及甲基分群／距離方法地景與 ISM 之可防守差異（2.5）。

## 2.1 腫瘤亞克隆重建之骨幹與「正式重建」之定義

亞克隆重建之主流以體細胞變異為骨幹（變異呼叫 → 拷貝數 → 癌細胞分率 → 分群 → 譜系）{{cite:Tarabichi 2021 / PyClone-VI 2020 待 citation-verification}}。**何謂「正式之重建」於文獻中有明確標準**：需推估 clone 數與組成、建立譜系樹或演化模型、處理拷貝數／純度去捲積，並具正交或單細胞 ground truth{{cite:BitPhylogeny 2015 待 citation-verification}}；單細胞甲基可直接重建 lineage tree{{cite:Gaiti 2019 / Sgootr 2023 待 citation-verification}}，bulk 甲基亦能推演化參數{{cite:EVOFLUx 2025 / 2026 CRC growth-history 待 citation-verification}}，而 CAMDAC 等則正式處理 CN／purity 之去捲積{{cite:CAMDAC/TRACERx 2025 待 citation-verification}}。

此一文獻地景**不直接否證本研究之結果，而是約束過寬之動詞**：在未完成大規模重建 demo 與正交驗證（本研究之 G-D／G-E 後續工作）前，本研究之正文一律採 **「characterize methylation-associated subclonal structure conditioned on somatic haplotags」**，而非「reconstruct subclones using methylation」；chr2 等案例為**案例示範與時序假說**，不輸出譜系樹。

## 2.2 ONT 長讀與 germline 甲基-相位：擁擠領域與開放白地

以甲基化延伸相位之方法已相當成熟，但**主要侷限於 germline 層級**：NanoMethPhase 於無親代基因型下偵測 ASM{{cite:NanoMethPhase 2021 待 citation-verification}}，MethPhaser／HapBridge／NANOME 以甲基補足 SNV 相位連續性{{cite:MethPhaser 2024 / HapBridge / NANOME 待 citation-verification}}；最新之 LongHap 亦為 germline-only{{cite:LongHap 2026 待 citation-verification}}。此一擁擠之 germline 地景反而**佐證**了一個開放白地：**tumor／LOH 之 methyl-assisted phasing 尚無完成同軸之直接論文**。惟此白地要成立為論文主張，須與 MethPhaser／HapBridge／LongHap 在相同 tumor 資料 head-to-head 比較 switch error／block N50／read-tag recovery、於 LOH 內具正交相位真值、並以多樣本與第二 pipeline 驗證——本研究將此列為後續工作而非已成立之主張。

## 2.3 癌症 ASM、cis 與拷貝數混淆：corroboration 與 cis-test 之前例

癌症中 ASM 升高已有文獻支持{{cite:Do2020 待 citation-verification（caliber 口徑差須標）}}，並與拷貝數變異之混淆密切相關{{cite:Martin-Trujillo 2017 待 citation-verification}}；read-level 之差異甲基偵測亦有 DAMEfinder 等前例{{cite:DAMEfinder 2020 待 citation-verification}}。需特別劃界 ISM 之 normal-anchored cis-test 之**增量**：**並非「獨家使用 matched-normal 錨點」**——Do2020 之 de-novo ASM 子分析亦使用 within-patient matched-normal，且其「somatic 局部 cis 甲基稀少／數值貢獻不大」之觀察**正佐證**本研究之 1/816 稀有{{cite:Do2020 de-novo ASM / Zhang2019 SV 待 citation-verification}}；ISM 之增量在於 **read-level within-sample 距離結構 + LOH／HP 軸 + 三者整合**。normal-anchor 所扣除之 germline-cis 基線有 meQTL atlas 為背景{{cite:GoDMC Min2021 / Onuchic2018 待 citation-verification}}；而 CpG-SNP 造成之假性 ASM 則為須誠實標註之 limitation{{cite:Shoemaker2010 待 citation-verification}}。

## 2.4 距離式與單細胞甲基重建之前例：劃界本研究之新穎性

本研究**不主張**距離式或 read-level 甲基本身為新穎：以距離／單細胞甲基重建 lineage 之方法已有 Sgootr（距離式甲基 lineage tree，真腫瘤單細胞）、Gaiti（單細胞甲基 lineage）、Epiclomal、MethylTree 等前例{{cite:Sgootr 2023 / Gaiti 2019 / Epiclomal 2020 / MethylTree 2025 待 citation-verification}}；ONT epiallele 多樣性與甲基頻率 DMR 亦有 PoreMeth2{{cite:PoreMeth2 2025 待 citation-verification}}。因此本研究之**可防守新穎性在於組合**：somatic haplotag-conditioned 之 read 結構 + matched-normal-anchored somatic cis-test + LOH／CN-aware 詮釋 + 有界之 5mC／5hmC 分析——於 tumor/normal read-level 之整合，而非任一單件。

## 2.5 甲基分群／距離方法地景與 ISM 之 5 桶定位

read-level ASM 與甲基分群之最相近工具包括 qFDRP（NHD 核之先例）、cvlr、ASMS、MethylBERT 等{{cite:qFDRP 2020 / cvlr 2023 / ASMS 2024 待 citation-verification}}。ISM 與之差異**須以方法本身描述**：ISM 為「**無監督之 read×read 距離矩陣 → 結構檢定（PERMANOVA）+ normal-baseline cis-test + somatic-subclone 目標**」，而非以「定序平台世代」或「缺顯著性檢定」等不準確之描述劃界（此類描述會被合理質疑）。與 read 層級異質性（epipolymorphism／methylation entropy）之度量{{cite:Landan 2012 / methclone 2014 待 citation-verification}}亦須區隔：本研究之先前分析（O11）證實此類熵度量之表觀判別力為覆蓋（read 數）混淆之假象（校正後 AUC 自 0.845 降至 0.530），ISM 不依賴熵而以可靠度閘控分群 + normal-anchored cis-test 控制此混淆。

綜合上述，本研究與外部文獻之關係可分為五類：**(1) 新穎（組合）**——C1／C2／C3 於 tumor/normal read-level 之整合；**(2) 佐證**——ASM 存在（Do2020）、甲基為 germline-haplotype 層級（cvlr／ASMS／O'Neill／LongHap）、弱-subclone 紅線（EVOFLUx／MethylTree）、甲基非變異過濾器（Kapoor／Soneson）；**(3) 開放白地**——tumor／LOH methyl-assisted phasing（待 head-to-head）；**(4) 前例（劃界）**——距離式／read-level 甲基重建（Sgootr／Gaiti）、somatic 局部 cis（Do2020 de-novo ASM／Zhang2019）；**(5) framing tensions**——「正式重建」之定義、距離法非新穎、6 cell line 非 cohort。需誠實陳述：本研究**尚未發現與其有界內部結果在同資料 regime 下直接相反之實驗結論，但存在多個高可信度之 framing tensions**，限制 reconstruction、novelty 與 generalization 之可宣稱範圍——此較「零衝突」更符合完整文獻地景。
