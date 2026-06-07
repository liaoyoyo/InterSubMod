<!--
建立時間: 2026-06-08
狀態: in_progress (G6 論文 methods-NEGATIVE 骨幹草稿 — 獨立於 HD-1，今天就防彈)
報告類型: paper_section_draft
受眾: 論文撰寫 (本人/PI) — 可直接改寫為投稿 Methods/Results-Negative
framework: SRE-blameless-postmortem × Assertion-Evidence (每段 = 主張 + 接地證據 + 誠實 caveat)
data_sources:
  - knowledge/11_external_literature/10_paper_readiness_convergence.md
  - research/methyl_augmented_filter_phase2/cycle4/loso_validation/
  - research/tsg_promoter_asm_reviewer/genome_survey_v2/{fast_cnv_validation,condition_fp_consensus}.json
  - docs/experiments/in_progress/2026/05/20260531_ISM_aux_tag_observation_funnel_01.md
  - research/autoresearch/evidence_ledger.jsonl
provenance_note: 每個數字取自收斂稽核 wf_9e169112-573（auditor 實讀上列 landed 檔）；草稿，投稿前各數字再 grep 原檔 + tsg 定稿核對。
-->
<!-- provenance-verified: 全部 metric 取自 knowledge/11_external_literature/10_paper_readiness_convergence.md（wf_9e169112-573 auditor 實讀 landed json/md）；本草稿 Write 與產數字分析不同 batch（§13.0）。投稿前 /citation-verification + 各數字 grep 原檔。 -->

# G6 論文 — Methods-NEGATIVE 骨幹草稿

> **這份是什麼**：論文最可防守的質量 = 三~四道防彈 NEGATIVE。本草稿把它們寫成可直接改寫投稿的段落（每段 = 主張 + 接地證據 + 誠實 caveat + 外部背書）。**獨立於 HD-1**（不論 phasing 脊柱怎麼決，這些今天就能寫）。
> **核心 framing（FT-1，每段守）**：這些是「甲基**不能**當 filter / cis-dosage **不是** driver」的**乾淨負結果 + 機制**，不是 filter。characterization-set（628/15391）絕不寫成 filter（15391 TP/FP=1.0）。
> **status**：DRAFT。tsg 專案仍 active（T2/T3 數字可能再動）；投稿前各數字 grep 原檔 + 以 tsg 定稿為準。

---

## §N1 — 甲基化不能判別 somatic 變異真/假（三獨立角度）

**主張**：以甲基化/表觀特徵篩 somatic variant TP/FP，在三個獨立角度上失敗；這不是方法瑕疵，而是公認 ML 陷阱 + 訊號本質的中性。

**證據**：
- **(a) LOSO sample-level 循環性**：HCC1395 in-distribution ΔF1 **+0.02236** 在留一驗證時塌成 **−0.00012**（drop +0.02248 = 100% sample-level circularity）；5-樣本 mean **−0.00004**，Wilcoxon p=0.125（0/5 正向），HCC1954 transfer **−0.377**。所有 best-τ 退化為 0.10（= keep everything）。[src: `methyl_augmented_filter_phase2/cycle4/loso_validation/loso_findings.md` + `loso_summary.json`]
- **(b) 直接判別力中性**：連續 |Δβ| per-position **AUC=0.505**，落在 2,000× label-shuffle null 內；即使最乾淨的「credible」regime AUC 0.570 仍在 null。[src: `20260531_ISM_aux_tag_observation_funnel_01.md`]
- **(c) powered 6-樣本 ISM native run（~280k loci）**：ASM 確實存在且 somatic-enriched（TP significant **3.95%** > FP **1.07%**，subhap-matched 3.77% vs 1.09% ≈3.5×），但絕對率僅 ~4% = 低 sensitivity；COLO829 TP **0.0089** ≈ FP **0.0103**（判別力消失）；cis HP_Residual 非判別（FP 最高 0.3652）。[src: ledger `20260604_ISM_complete_TPFPFN`]

**外部背書（L3）**：Kapoor & Narayanan 2023 (*Patterns*, PMID 37720327) 的 leakage taxonomy **L3.2（樣本非獨立）+ L1.3（select-then-train）** 正命名我們的 LOSO 循環；Soneson 2014 (*PLoS One*, DOI 10.1371/journal.pone.0100335) 模擬顯示 confounded-CV 的「in-dist 完美 / external 隨機」標誌——我們 in-dist +0.0224 / held-out −0.0001 正是此標誌。**無任何 production somatic caller（ClairS-TO / VarNet / AIVariant / cfDNA-ML）用 CpG 甲基當 standalone FP 判別器**——它們用 alignment/PoN/purity/CN。

**caveat**：filter-negative 全 single-pipeline（ClairS-TO，無第二 caller 重驗）；ASM 機制工作單樣本 HCC1395-主導；未測甲基當*邊際* DL 特徵（AIVariant 唯一報有益的配置）。**provenance 修**：narrative「0.5049」→ source 真值 **0.505**（去掉假第 4 位小數）。

> **FT-1 守則**：本段是「甲基非 filter」的證明；prose 絕不出現 filter 動詞描述 628/15391。

---

## §N2 — Copy-dosage 不是甲基分群的 driver（乾淨 falsification）

**主張**：腫瘤甲基分群的 magnitude 不由 copy-number dosage 驅動——這是全程式最乾淨的單一結果，且直接 back-stop 文獻的 CN-confound 警告。

**證據**：CN-state class-contrast |Δβ| 跨 neutral/gain/cnLOH/loss 平坦（0.070–0.082）；neutral-nonLOH vs gain-nonLOH **MW p=0.6183（不可區分）**；signed **ρ(|Δβ|,CN)=−0.083（p=2.6e-9, 反向）→ dosage-artifact 模型 REFUTED**；cnLOH（copy-neutral）分群 ≈ neutral（0.082 vs 0.081）= 改 copy 結構不改分群 magnitude。[src: `tsg.../genome_survey_v2/fast_cnv_validation.json` @06-07]

**外部背書（L3）**：Martin-Trujillo 2017 (*Nat Commun*, PMID 28883545)「CN 解釋腫瘤 imprinted allelic 甲基 82–92%」——我們的 somatic-controlled HP-axis 設計 held-constant CN，dosage falsification 顯示我們**非** CN-dosage-driven，正是該文警告的 confound 被我們控掉。

**caveat（關鍵·誠實）**：copy-**partition**（非 dosage）仍在最強訊號上 confound HP-axis magnitude——HP-axis cis 在 amplified 區 conflate copy 與 allele。816 HP-axis loci 經 Bonferroni+copy-test 後**唯一乾淨 cis = chr17:79,991,120 (TBC1D16)**（within d=0.142 > HP-axis 0.122, perm p=0.001）；**BRCA2 (chr13) + chr5 = copy-artifact**（BRCA2 HP-axis −0.122 ≈ d_copy −0.11 + 邊際 d_within −0.023, perm p=0.022, n=197）。**6 個最強訊號 untestable**（pure-ALT tag, CGI-desert；copy-test 只能算 leaky-tag 7/60）。partition 未 anchor SEQC2 CN ground-truth。⇒ **勿寫「BRCA2 真 cis-driven」；BRCA2 = 小 copy-confounded 例，chr17/TBC1D16 = 唯一乾淨 cis exemplar（稀有但真實，單 locus 單樣本 n=16/14）**。

---

## §N3 — 高-Δβ loci 因 regression-to-extreme 聚集於 FP（機制解明）

**主張**：「strong-ASM 在 FP 富集」不是甲基的判別力，而是一個可機制解釋的 artifact——把負結果轉成*已解釋*的負結果（論文加分）。

**證據**：HCC1395 條件 Fisher：**非低覆蓋**（OR=0.87，非因素），而是 **no-clustering（OR=8.63）+ LOH（OR=4.09）+ spurious small-subhaplotype（OR=5.84）**；combo 高-Δβ+no-clustering+LOH 給 TP **0.86%** vs FP **7.97%**（≈9×）。機制：LOH→single-haplotype→低覆蓋→極端 baseline 是 mechanistically correlated 的 ONE regime，非 3 個獨立維度。**ISM CramersV gate 正確排除這些**（需 discrete clustering，Δβ-only 抓不到）→ 兩方法互補（BRCA2 是 clustering-only, Δβ~0.007）。[src: `tsg.../condition_fp_consensus.json` + ledger `20260608_..._capstone`]

**caveat**：單樣本 HCC1395（唯一 FP 夠多 627–4,842 能 power Fisher）；bootstrap stability 未跑；FDR = BH-on-min-p（anti-conservative）；跨樣本未驗。⇒ 寫「single-sample HCC1395 worked example，bootstrap pending」。先前「FP |Δβ| 較大」claim **已 RETRACTED**（abs_dbeta MW p=0.137 NS）；strong-ASM 5× FP **OR=0.194 anti-discriminative**（61% 來自單一 chr8 hotspot，drop-chr8 → OR 5.16→2.84）。

---

## §N4 — somatic 分群是 germline-allelic 非 somatic-specific

**主張**：read-level 甲基分群（subclone-甲基訊號）不是 somatic 驅動，是 germline allelic 背景。

**證據**：validated blind-ARI ruler（imprinting 正控 GNAS/RB1 = **1.000**，normal DMR median 0.758 → ruler 證實 sound）：somatic TP median ARI **0.135** < germline-het null **0.177**，MW p=**0.922**（無 somatic-specific enrichment）。[src: master_draft §2.2 line 111; ledger 2026-05-30/31]

**caveat**：single-pipeline；HP-axis OR=1.79「somatic-specific」是 threshold-leniency artifact（threshold 收緊到 Bonferroni 時 OR 1.79→0.81 反轉）。

---

## §N5 — 跨樣本 locus/gene 復發 absent / underpowered

**主張**：somatic ASM 現象**現象層**跨樣本復現，但**locus/gene 層**復發 absent，且這個 absent 是 underpowered（誠實三層區分）。

**證據（三層）**：
- **phenomenon-level YES**：6/6 樣本 excess-over-null **+0.101–0.241**（mean 0.168, CV 0.26, 3 癌種）。[src: `tsg.../cn_confound/cross_sample/*_gwasm.json`]
- **locus/gene-level NO**：同-locus somatic 復發 **0/38**（HCC1395 33/38 somatic, 28 somatic-ASM, 其他 5 樣本 0/38）；credible-gene 復發 **0/94**（Jaccard=0）。[src: `cross_sample_synthesis.json`]
- **CN-truth generalization 結構不可達**：只 HCC1395 有 SEQC2 ground-truth。

**caveat（必隨行·FT-3）**：0/38 是 **UNDERPOWERED**（E[隨機 overlap]=0.16，觀察 0 與 chance 不可分）→ 立「**復發稀有**」**非「loci 是 private**」。6/6 是 **allelic-methylation-difference 復現**，**非 cis-ASM 復現**（FT-4；只 chr17 乾淨 cis）。single-pipeline ClairS-TO → ⭐3。raw median |Δβ|<0.10（訊號僅 null 減後浮現），COLO829（melanoma 腿）最薄（+0.101, 0 credible）。**外部對照**：O'Neill/POG 2024（189 腫瘤）偵 RET/CDKN2A 復發 aDMR——我們 6 樣本單 pipeline 偵不到（這句使 0/38 誠實而非 null finding）。

---

## §N6 — 補充負結果（同方向 corroboration，強化底座）

> critic 補：filter 不是死三道，是死**四**道——以下三個獨立死法強化 §N1。

- **Zone-Aware QS-simulation F1 NEGATIVE**：5 zone characterization H1/H3 POSITIVE，但 QS-simulation filter max **ΔF1=+0.001**（TO QS AUC=0.497 隨機）。[src: INDEX:225-226]
- **CNV zone-aware filter（2026-04 CLOSED）**：per-sample zone-specific mean AUC ≤0.641，cross-sample 不一致，Simpson's-Paradox-negated。[src: INDEX:198-199,318]——這是 §N2 dosage-refutation 的 4月先祖（可 pre-empt reviewer「試過 CN-zone filter 嗎」）。
- **V6 LR baseline-equivalence**：4-way head-to-head 顯 baseline LR ΔF1 **+0.02302 ≥ V6 LR +0.02236**——整個 +0.023 uplift 是 **BAM-independent**（caller_af/LOH/Cov），「NG 只是 10 特徵之一，機制非 V6-specific」。[src: INDEX:263]——bounds T1 tooling claim + reinforce §N1。

---

## 撰寫狀態 + Provenance

- **就緒度**：§N1/§N2/§N4/§N5 = 今天可寫（grounded）；§N3 加「single-sample HCC1395 + bootstrap pending」caveat；§N6 = corroboration。**全部獨立於 HD-1**（phasing 脊柱怎麼決都不影響這些）。
- **每段已內建 FT-1~5 + OR-1~6 防護**（見 `knowledge/11_external_literature/10_paper_readiness_convergence.md §7`）。
- 數字取自收斂稽核 **wf_9e169112-573**（auditor 實讀 landed json/md）。**投稿前**：(1) 各數字 grep 原檔 + 以 tsg 定稿為準；(2) `/citation-verification` 核外部 PMID/DOI；(3) §N2 BRCA2 口徑以 06-08 amendment 為準。
- 撰寫此草稿的 Write 與產數字分析**不同 batch**（§13.0）。
