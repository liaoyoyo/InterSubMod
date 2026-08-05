<!--
建立時間: 2026-06-15（v3 對齊 06-14 Intro 主軸定稿：ISM 首個整合 C1/C2/C3 + 旗艦 BRCA2/chr2 + characterize-conditioned 措辭；取代 v2 de-confound-headline 版）
報告類型: 碩論摘要（中文摘要 + English Abstract）
狀態: draft v3；待 review
data_sources: docs/paper_focus/02_paper_framework/20260614_Intro主軸定稿_thesis_statement_01.md, docs/paper_focus/00_共識證據台帳_20260612_01.md, flagship_chr2_18086020_subclone HTML
provenance_note: ISM 定位/旗艦引 06-14 定稿（用戶確認）；catalog 332705/12868/3.87%/1-816、BRCA2 Δβ=−0.122/d_within=−0.023、chr2 F=10.6/V=0.67/V=0.97 grep-verified；甲基→突變因果歸文獻 {{cite 待補}}。措辭守 framing 凍結（characterize-conditioned 非 reconstruct）。
-->

# 摘要

腫瘤亞克隆結構與演化之重建傳統上以體細胞變異為骨幹。Oxford Nanopore（ONT）長讀定序能於同一條 read 同時取得體細胞變異、germline 單倍型與原生 5-methylcytosine（5mCG）甲基化，惟以甲基化延伸相位之方法多侷限於 germline 層級，而既有之亞克隆重建工具多不處理 read 層級之甲基。在僅有單樣本之癌症資料上、以 somatic haplotag 為條件、於 read 層級對「甲基與體細胞變異之聯合結構」進行 normal-anchored 之特徵化，仍是一個未被覆蓋之利基。

本研究提出並評估 **InterSubMod（ISM）**——首個於 ONT 長讀 tumor/normal paired 資料上、整合「① read 層次甲基結構檢定（read×read 距離→分群→PERMANOVA）+ ② matched-normal-anchored somatic cis-test（殘差 + d_cis/d_drift/d_within 分解）+ ③ somatic-subclone 目標」之 read-level 軟體。其新穎性在於三者於 tumor/normal read-level 之**整合**（個別元件皆有前例），用以量化 somatic 變異與甲基位點之**關聯強度**。

於六株癌細胞株（三癌別）之 332,705 位點目錄，可靠之甲基 read 分群僅 12,868（3.87%）且壓倒性跟隨 germline 而非 somatic、乾淨之 somatic cis 於可測位點中少見（816 中 1，與既有文獻「somatic 局部 cis 甲基稀少」一致），且甲基化不能做變異真偽過濾器。**正因如此**，本研究以兩旗艦案例示範主軸：（一）**BRCA2 promoter**——HP-axis 甲基差 Δβ=−0.122 與 tumor-normal 差並存，以 matched-normal 錨點**分解**顯示此差主要 track somatic subclone tag、focal cis 效應 marginal（d_within=−0.023，perm p=0.024）；ISM 不只標記「此處有與突變共定位之甲基差」，更**量化**其屬 subclone 結構抑或直接 focal cis。（二）**chr2:18M**——以 longphase-S 系譜（Normal=G → HP1 母本 19 reads 全 G → HP1-1 子克隆 36 reads 含 28 A）為框架，於標示子克隆之位點甲基分群顯著（PERMANOVA F=10.6，分群↔單倍型 Cramér's V=0.67），於另一 somatic 等位位點甲基分群↔等位 V=0.97；兩種互補模式並存，示範甲基如何與突變**共同標示**同一子克隆事件，據以推導 Normal→HP1→HP1-1 之時序假說。

本研究之主要邊界與後續工作：上述為**案例示範**而非已驗證之自動化重建（不輸出譜系樹，時序為數據支撐之假說，自動化與大規模驗證為延伸）；機制建於單樣本、單一定相流程（證據層級 ⭐3，TP/FP 取自 tumor-only caller，六 cell line 為跨模型 reproducibility 非病人 cohort）；甲基與突變之因果機制屬既有文獻背景，ISM 量化關聯強度而不宣稱因果或順序；within-haplotype 之 subclone 甲基 somatic-特異性尚待正確對照確認。

**關鍵詞**：奈米孔長讀定序、等位基因特異性甲基化、somatic haplotagging、normal-anchored cis-test、腫瘤亞克隆結構、單樣本

---

# Abstract

Reconstructing the subclonal architecture and evolution of tumors has conventionally relied on a somatic-variant backbone. Oxford Nanopore (ONT) long-read sequencing captures somatic mutations, germline haplotypes and native 5-methylcytosine (5mCG) on the same read, yet methylation-based phasing is largely germline-level and existing subclonal-reconstruction tools rarely handle read-level methylation. Characterizing the joint structure of methylation and somatic variants at read level—conditioned on somatic haplotags, normal-anchored, in single-sample cancer data—remains an uncovered niche.

We present and evaluate **InterSubMod (ISM)**, the first read-level software on ONT tumor/normal paired data to integrate (i) a read-level methylation structure test (read×read distance → clustering → PERMANOVA), (ii) a matched-normal-anchored somatic cis-test (residual + d_cis/d_drift/d_within decomposition), and (iii) a somatic-subclone target. Its novelty lies in the **integration** of these three at tumor/normal read level (each component has prior art), quantifying the **association strength** between somatic variants and methylation loci.

In a catalogue of 332,705 loci across six cancer cell lines (three cancer types), reliable methylation read-clustering exists at only 12,868 loci (3.87%) and is overwhelmingly germline- rather than somatic-allelic; clean somatic cis is rare among testable loci (1/816, consistent with literature on the scarcity of local somatic cis-methylation); and methylation does not discriminate true from false somatic variants. **Against this backdrop**, two flagship cases demonstrate the main thread. (i) The **BRCA2 promoter**: an HP-axis methylation difference (Δβ=−0.122) coexists with a tumor–normal difference; normal-anchored **decomposition** shows it mainly tracks the somatic subclone tag, with a marginal focal-cis effect (d_within=−0.023, perm p=0.024)—ISM not only flags a methylation difference co-localized with a mutation but **quantifies** whether it reflects subclone structure or direct focal cis. (ii) **chr2:18M**: within a longphase-S lineage (Normal=G → HP1 parent, 19 all-G reads → HP1-1 subclone, 36 reads incl. 28 A), methylation clustering at the subclone-marking locus is significant (PERMANOVA F=10.6; clustering↔haplotype Cramér's V=0.67), while at a somatic-allele locus clustering↔allele V=0.97; these complementary patterns demonstrate how methylation and mutation jointly mark the same subclonal event, supporting a Normal→HP1→HP1-1 temporal hypothesis.

Principal boundaries and future work: these are **case demonstrations**, not validated automated reconstruction (no lineage tree is output; the temporal order is a data-supported hypothesis; automation and large-scale validation are extensions); the mechanism rests on a single cell line and single pipeline (evidence tier ⭐3; TP/FP from a tumor-only caller; six cell lines = cross-model reproducibility, not patient-cohort generalization); methylation–mutation causality is existing literature background, with ISM quantifying association rather than claiming causality or order; and the somatic-specificity of within-haplotype subclone methylation awaits the correct control.

**Keywords**: nanopore long-read sequencing; allele-specific methylation; somatic haplotagging; normal-anchored cis-test; tumor subclonal structure; single-sample
