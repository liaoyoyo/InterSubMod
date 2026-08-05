---
id: ism-kb-11-external-literature-paper-readiness-convergence
name: "論文就緒收斂 — 終極目標·可寫清單·矛盾·GO-NO-GO"
description: "跨 7 條研究線獨立 fresh-context 稽核（實讀 ledger/json/report）的收斂：終極目標 + 論文敘述主軸、哪些 usable/failed/needs-thought、現在就能寫的章節、4 個方向逐一裁決、2 個真矛盾、6 個 HD 決策、決定性 R-SELFREF GO/NO-GO。經 workflow wf_9e169112-573。"
status: active
last_verified: 2026-06-08
content_nature: paper-derived
doc_type: explanation
verified_scope: "7-thread independent audit (general-purpose agents read landed ledger/json/reports) in wf_9e169112-573 (2026-06-08); all numbers grounded in cited files"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-asm-cis-cancer-impact
tags: [paper-readiness, convergence, ultimate-goal, thesis, usable-failed, contradictions, go-no-go, r-selfref, audit]
canonical_paths: [11_external_literature/10_paper_readiness_convergence.md]
alias_paths: []
---

<!-- data_sources: research/autoresearch/evidence_ledger.jsonl, research/tsg_promoter_asm_reviewer/genome_survey_v2/*.json, research/methyl_augmented_filter_phase2/, research/loh_subclone_af_paired/, docs/reports/validated/2026/06/.../master_draft.md, docs/experiments/INDEX.md -->
<!-- provenance-verified: 全部 metric 由 workflow wf_9e169112-573 的 7 個獨立 auditor 實讀上列 landed 檔取得（fresh-context，不信既有結論、自行重判）；本 .md Write 與 auditor 的檔案讀取在不同 batch（§13.0）。內部 L1-L2 / 外部 L3。 -->

# 論文就緒收斂 — 終極目標·可寫清單·矛盾·GO/NO-GO

- **一句結論**：**一篇能過 review 的論文今天就存在** —— read-level LOH/haplotype + 甲基化的 **characterization + tooling** 論文，主體是**三~四道防彈的 NEGATIVE** + 一個方向一致的跨樣本 phasing 訊號(Grade **B+** 非 A) + ASM 當誠實的 copy-confounded 支撐。**唯一決定性卡關 = 一個 GO/NO-GO 研究決定（HD-1，見 §5），由你決定。**
- **這份怎麼來的**：7 條研究線各派 1 個獨立 fresh-context auditor，**實讀** ledger/json/report 真值、不信既有結論自行重判、抓 overclaim/矛盾、分類 usable/failed/needs-thought/paper-ready + 跨線矛盾獵手 + 收斂 + critic（workflow wf_9e169112-573，10 agents，2026-06-08）。**這同時就是在跟其他 session 的落地產出對帳。**

> ⚠ **本文件性質**：paper-readiness 決策層。內部數字為 auditor 實讀 json 的真值（grounded），非本輪新跑分析。引用前看每項的 caveat。

---

## 1. 終極目標 + 論文敘述主軸（收斂後）

**終極目標（最可防守）**：在癌症長讀（ONT 5mC/5hmC）資料中，以**單 read 解析度** characterize **somatic 突變與 LOH 如何重塑 haplotype/甲基地景**，並交付利用它的 tooling（LOH-constrained phasing signature）。**這是 read-level epigenetic & haplotype CHARACTERIZATION + TOOLING 論文，明確不是 variant-calling filter。**

**一段話主軸（所有存活證據都支持的）**：
> *在 LOH 區，腫瘤物理上只保留單一 haplotype，故 somatic SNV 形成同-haplotype 子家族，跨樣本可偵測（NG=2 Inner same-HP1 > Outer cross-het，7/7 方向一致，附乾淨 paired 負控）。Allele-specific 甲基化（ASM）在這些腫瘤中確實存在、且現象層跨 6 樣本/3 癌種復現，但它小、非方向、與 copy-dosage 無關、且關鍵地**不能判別真/假變異**；甲基化因此**不能當 variant filter**（三~四種獨立方式證死）。論文貢獻 = LOH-driven haplotype 結構的 read-level characterization + 一份誠實、機制解明的「甲基能做/不能做什麼」目錄。*

**對既有 G6 三軸 framing 的必要修正（因 06-07/06-08 收緊）**：三軸結構（主=LOH phasing / 支撐=ASM characterization / methods=誠實負結果）**存活且正確**，但 4 個標籤投稿前必改：
1. 🔴 **「Grade A」→「Grade B+ characterization」** —— 決定性 R-SELFREF circularity 對照（全 7 樣本 flag-on，~25-50hr C++）**經磁碟確認未跑**；headline 統計是方向一致性非 effect 大小。
2. **ASM 支撐軸是 characterization + supplement-null 軸，非 positive-discovery 軸** —— 所有 positive-discrimination framing 在數據中已撤回。
3. 🔴 **BRCA2 不再是 cis 錨點 = subclone/copy-confounded（06-07）** —— HP1-1 是 longphase-S 的 somatic subclone tag〔germline-H1+somatic-ALT〕**非 copy、非 CN-dosage**〔dosage 已 REFUTED〕；唯一乾淨 cis exemplar 是 **chr17:79991120 (TBC1D16)**；BRCA2 hero 敘事退役（% split 不 robust，勿寫精確 %）。
4. **methyl-assist phasing (umtag) 是 Discussion/future-work + 1 張 methods-grade copy-control 圖（V10）**，非 result section。

---

## 2. 論文就緒主表（7 線）

| 線 | 論文角色 | USABLE 現在可寫 | FAILED / 反例 | 要再想清楚 | 就緒度 |
|---|---|---|---|---|---|
| **T1 LOH-phasing** | **脊柱（主 +）** | NG=2 Inner same-HP1 > Outer **7/7** 方向一致 (W=28, p=0.0078125, median gap 0.345, CI[0.104,0.459])；Inner×NG=2 ≥93% same-hap 6/6；C++ 機制 (LabelTest.cpp:265-302)；**乾淨 paired 負控** (B3 median≈0, NS) | methyl→phasing-derived filter DEAD；HPFineNGroups ⭐4→⭐3；V5 Layer-1.5 繼承 priority-bug | 🔴 **R-SELFREF circularity 對照未跑**（Inner 用 somatic-attributed HP1-1 定義 → by-construction 部分循環）；HCC1954/H2009 fragile；n=7 = bio-n=6 | **NEEDS-WORK（Grade B+ 非 A）** |
| **T2 ASM 存在** | 支撐（characterization）| ASM 真但 ~4% sens、非 filter；copy-dosage REFUTED (MW p=0.6183, ρ=−0.083)；high-Δβ→FP 機制 (no-clust OR=8.6/LOH OR=4.1)；6/6 跨樣本 excess | somatic_specificity NEG；B2 clustering NEG (germline-allelic)；FP strong-ASM = regression-to-extreme ONE regime；「628/15391 當 filter」矛盾來源 (TP/FP=1.0) | 313 vs 172 Bonferroni 軸；COLO829 才能 ⭐4；「pervasive」overclaim (median \|Δβ\|<0.10) | **READY-AS-NEGATIVE / characterization** |
| **T3 cis-test / subclone-copy-partition** | 支撐（per-locus + falsification）| copy-**dosage** 乾淨 falsified；normal-anchored cis-test 方法；subclone/copy-partition 分解（HP1-1=somatic subclone tag 非 copy）；**chr17/TBC1D16 = 唯一乾淨 cis** (perm p=0.001, within>HP) | 🔴 **BRCA2 + chr5 = subclone/copy-confounded**（cis 錨退役；非 CN-dosage）；chr18 mechanical CpG-create；「60 T3」塌成 1 | 6 strongest UNTESTABLE (pureALT CGI-desert)；chr17 單 locus/單樣本；未 anchor SEQC2 CN；% split 不 robust | **READY-AS-NEGATIVE / methods** |
| **T4 methyl→TP/FP filter** | methods（誠實負結果）| LOSO 100% circularity (in-dist +0.0224 → held-out −0.0001)；\|Δβ\| **AUC=0.505 (null)**；ISM vestigial (drop-5-methyl ΔF1 +0.00065)；FP-enrichment 機制解明 | 整個 production filter DEAD；「FP \|Δβ\| 較大」RETRACTED (p=0.137 NS)；strong-ASM 5× = anti-discriminative (OR=0.194) | 所有機制單樣本 HCC1395；bootstrap stability 未跑；COLO829 TP≈FP | **READY-AS-NEGATIVE（強）** |
| **T5 跨樣本 ASM** | 支撐（generality）| 6/6 excess-over-null +0.101–0.241 (3 癌種)；discovery 方法 transfers 6/6；38/38 keypos ≥10× | 同-locus somatic 復發 **0/38**；credible-gene 復發 0 (Jaccard=0)；「4 imprinting」= sign-concordant only | 0/38 UNDERPOWERED (E[overlap]=0.16)；single-pipeline 天花板；COLO829 最薄腿 | **READY-AS-NEGATIVE / characterization** |
| **T6 subclone (AF→NGroups)** | 支撐（1 正向描述）| AF→NGroups **7/7** (Δng +0.637–0.910)，survives NumReads 分層；dosage refuted；LOH dual-hap diagnostic partition (72/18/9) | B2 somatic clustering NEG；HPFineNGroups-as-filter ⭐4→⭐3；cis line 推翻 (1/816 clean) | 🔴 **AF→NGroups: 甲基 vs phasing 歸因 UNRESOLVED**；「~64% partition」= 7/60 convenience subset | **READY-AS-NEGATIVE (+1 正向 IF 歸因解決)** |
| **T7 methyl-assist phasing (umtag)** | future-work + 1 methods 圖 | opportunity space (45.84% unphased)；held-out extrap median 0.8852 (null 0.5236)；**V10「非 copy」控** (normal 0.979 > tumor 0.866, depth-matched) | aDMR×LOH NOT enriched (Fisher p=0.382)；germline-het null 非判別 (0.974)；2 次 fabrication（已修） | held-out ≠ 真-unphase 救援；無 switch-error/N50 yardstick；無 orthogonal phase-truth；單樣本；**ledger+INDEX 未註冊** | **NOT-READY（prototype / Discussion-only）** |

---

## 3. 現在就能寫進論文的（具體章節）

### POSITIVE 結果（今天可寫，caveat 必隨行）

- **§主 — LOH-constrained phasing signature（唯一跨樣本 positive）**：NG=2 Inner same-HP1 TP-rate > Outer cross-het，7/7，W=28, p=0.0078125, median gap 0.345, CI[0.104,0.459]；Inner×NG=2 ≥93% same-hap 6/6；C++ 機制 (LabelTest.cpp:265-302)。**必隨行 caveat**：(a) frame 為 **Grade B+ characterization 非「Grade A」**；(b) p=0.0078 量**方向一致性**（W=28 = n=7 最大 signed-rank，全正即被 forced），magnitude 跨樣本變 11×（+0.05~+0.52），勿當大 effect；(c)「n=7」生物上是 **n=6**（HCC1395 算兩次：5kHz + DORADO）；(d) single-pipeline ClairS-TO → ⭐3；(e) **R-SELFREF circularity caveat 首段明說**（Inner same_HP1 用 somatic-attributed HP1-1 bucket）。
- **§支撐 — paired 負控**：paired gap median≈−0.0002, W=10 p=0.578 NS；TO-vs-paired p=0.0625。乾淨、可寫，證訊號 TO-specific。
- **§支撐 — AF→NGroups 關聯（7/7）→ 🔴 HD-4 已判定為 phasing 非甲基**：Δng +0.637–0.910, MW p≈0, survives NumReads 分層，**但 06-08 HD-4 證實 NGroups 是 HP-tag sub-family count（非甲基矩陣）→ AF→NGroups 是 phasing/haplotype-occupancy signature 非甲基發現**（見 §5 HD-4）。**改述為 phasing signature**（反而支撐 phasing/LOH 主軸，非甲基支撐）；**勿當「subclone 帶不同甲基」之正向甲基結果**。
- **§支撐 — chr17/TBC1D16 唯一乾淨 cis exemplar**：within-somatic d_within=0.142 > HP-axis 0.122（唯一這種 locus），perm p=0.001, CGI 差異 CpG (p=1.8e-7)。**caveat**：單 locus、單樣本 (n=16/14)、nominal p、無 replication →「稀有但真實 exemplar」非可推廣。
- **§支撐 — 跨樣本現象復現**：6/6 excess-over-null +0.101–0.241, mean 0.168, 3 癌種。**caveat**：這是 *allelic-methylation-difference* 復現，**不是 cis-ASM 復現**；只能 phenomenon-level 陳述。

### NEGATIVE 當貢獻（methods 主體 — 最強、可自信寫）

> 🔴 **filter 不是死三次，是死四次**（critic 補）：T4 三角度 + Zone-Aware QS F1 NEG + CNV-zone filter（4月關閉）+ V6-LR baseline-equivalence（uplift 是 BAM-independent）。

- **§Methods-Neg 1 — 甲基不能 filter 變異（三獨立角度）**：(a) LOSO 100% sample-level circularity (in-dist +0.02236 → held-out −0.00012；HCC1954 transfer −0.377；5-sample mean −0.00004, Wilcoxon p=0.125)；(b) \|Δβ\| **AUC=0.505**（2,000× label-shuffle null 內；最乾淨 credible regime 0.570 仍在 null）；(c) powered 6-sample ISM：ASM 真 (TP 3.95% > FP 1.07%, ~3.5× subhap-matched) 但 sensitivity ~4%, COLO829 TP≈FP。**provenance 修正**：「0.5049」是 narrative-only，source funnel 檔是 **AUC=0.505**——投稿前去掉假的第 4 位小數。
- **§Methods-Neg 2 — copy-dosage 非 driver**：MW neutral-vs-gain p=0.6183 不可區分；signed ρ(\|Δβ\|,CN)=−0.083 反向 → dosage-artifact 模型 REFUTED。**全程式最乾淨的單一結果**，直接 back-stop Martin-Trujillo CN-confound 警告。**caveat**：copy-*partition*（非 dosage）仍在最強訊號上 confound HP-axis；未 anchor SEQC2 CN。
- **§Methods-Neg 3 — high-Δβ loci 因 regression-to-extreme 聚在 FP（機制解明）**：非低覆蓋 (OR=0.87) 而是 no-clustering (OR=8.6)+LOH (OR=4.1)+small-subhap (OR=5.8)；combo TP 0.86% vs FP 7.97% (~9×)。把負結果轉成 *已解釋的 artifact*。**caveat**：單樣本 HCC1395（唯一 FP 夠多能 power Fisher）；bootstrap stability 未跑。
- **§Methods-Neg 4 — somatic clustering/subclone-甲基是 germline-allelic 非 somatic-specific**：blind-ARI ruler（imprinting 正控 GNAS/RB1=1.000 validated）：somatic TP median 0.135 < germline-het null 0.177, MW p=0.922。
- **§Methods-Neg 5 — 跨樣本 locus/gene 復發 absent/underpowered**：同-locus 0/38, credible-gene 0 (Jaccard=0)。**必隨行 caveat**：0/38 是 **UNDERPOWERED** (E[overlap]=0.16) → 立「復發稀有」**非「loci 是 private」**。

---

## 4. 你的 4 方向 + 脊柱 — 逐一誠實裁決（能否完成、當作什麼）

- **phasing 脊柱 (LOH-constrained)**：**能完成** —— 當 **characterization + tooling positive（Grade B+, ⭐3 single-pipeline）**；「Grade A」需未跑的 R-SELFREF C++ 對照，否則寫時 circularity caveat 首段 + p=0.0078 frame 為方向一致性非 effect 大小。
- **找甲基差異 + 驗證（ASM 存在/characterization）**：**能完成** —— 當 **characterization + 防彈 supplement-null**（ASM 真、copy-dosage-independent、非方向、非判別；chr17/TBC1D16 唯一乾淨 cis；BRCA2 = 小 copy-confounded 例）；**不是** positive discovery（單樣本 ⭐3，⭐4 需 COLO829）。
- **methyl → TP/FP filter**：**能完成** —— 當**強、機制解明的 NEGATIVE**（LOSO 100% circularity + \|Δβ\| AUC=0.505 null + ~4% sens + regression-to-extreme）；全程式最可防守的負結果。**DEAD as filter — 勿重開。**
- **subclone (HPFineNGroups / AF→NGroups)**：**能完成但非甲基正向** —— HD-4(06-08) 判定 AF→NGroups(7/7) 是 **phasing/haplotype-occupancy signature 非甲基發現**（NGroups=HP-tag count）→ 它**支撐 phasing 主軸不是甲基支撐**；甲基側 subclone = **乾淨負結果包**（dosage refuted、B2 clustering 非 specific、HPFineNGroups-as-filter 降級）；**不是 filter、不是乾淨 somatic-cis story**（1/816 clean、單樣本硬 confound 天花板）。**T6 不再有「正向甲基發現」**。
- **methyl-assist phasing / 救 unphase read (umtag)**：**不能當 result section** —— **future-work / Discussion novelty + 1 張 methods-grade copy-control 圖 (V10)**；無 scale 級真-unphase 救援、無 switch-error/N50 yardstick、無 orthogonal phase-truth、單樣本、ledger+INDEX 未註冊。白地真實且外部已驗，但交付是 roadmap 非 result。

---

## 5. 投稿前必先想清楚（排序）+ 2 個真矛盾

> 只有 **2 項是真·數據對數據矛盾**（須改文件）；其餘是看似矛盾的 scope/caliber 差異，明説即可。

🔴 **HD-1（決定性 — 卡論文最強 claim）：跑 R-SELFREF 還是把 phasing 降為 characterization？** 全 7 樣本 flag-on circularity 對照**經磁碟確認未跑**。因 Inner same-HP1 *用* somatic-attributed HP1-1 bucket 定義，headline positive **by construction 部分循環**。critic 加碼：**「加 caveat 然後保留 Grade-A-adjacent framing」不是 methods reviewer 會接受的選項**——只有 (1) **投稿前跑 ~25-50hr flag-on 7-樣本對照**（唯一讓脊柱防彈），或 (2) **把 phasing 從 positive-spine 降為 characterization observation，讓三~四個 NEGATIVE 扛論文**（接受它變 methods/negative 論文）。**這是 GO/NO-GO 研究決定（Hard-Gate-adjacent）→ 你決定，AI 不代決。**

🔴 **HD-3（真矛盾 #1 — ✅ 2026-06-09 R1/R2/R3 統一口徑落地）：把 BRCA2 改成 06-07 狀態。** `master_draft.md` line 122 原寫「BRCA2 真 cis-driven (d_cis=−0.142/−0.152)…cis-test 雙軸通過」；`07_asm_cis_cancer_impact.md` 原寫「cis-driven ASM 真…60 T3+4 survivor」。兩者被 06-07 重分析（BRCA2 = subclone/copy-confounded；只 chr17 clean）矛盾。✅ **honest 口徑（已落地）**：BRCA2 不是純 artifact，HP-axis Δβ=−0.122 是 **somatic-subclone vs germline 甲基差**（**HP1-1 = longphase-S 的 somatic subclone tag〔germline-H1+somatic-ALT, HaplotagStrategy.cpp:505-516〕，非 copy；非 CN-dosage**），focal cis 殘餘 d_within=−0.023 邊際（perm p=0.024）；單樣本與 subclone 不可分離（CAMDAC 同此）。**勿寫精確 %（split 不 robust）、勿寫「copy number」、勿從一個 overclaim 換成相反 overclaim。** ⚠ **Martin-Trujillo (CN-DOSAGE) 不 corroborate** —— copy-DOSAGE 已決定性 REFUTED（MW p=0.6183, signed ρ=−0.083 反向）；BRCA2 confound 是 subclone 非 CN-dosage，引 CN-dosage 文獻當 corroborate 是 **category error**（其正確角色 = 我們控掉的 confound 警告）。

🔴 **HD-2（真矛盾 #2 內含 — Grade A vs B+）：active-cycle title / master §0 寫「Grade A」，但同文件 body + Grade A doc §4/§17 寫 Grade B+（req #2 未達）。** 收斂裁決：**B+ 非 A**。投稿前統一。

**HD-4 ✅ RESOLVED（2026-06-08 背景分析 ledger `20260608_HD4_af_ngroups_attribution`）= PHASING-DRIVEN，非甲基發現**：C++ 確認 NGroups = `LabelTest::hp_to_fine_labels` 純從 HP tag 數 {HP1,HP1-1,HP2,HP2-1} sub-family（capped 4；甲基矩陣只進 HPFineF PERMANOVA，**不進 count**，LabelTest.cpp:265-305）。LOH 區 NGroups 1→2 = somatic ALT sub-family 出現 = **AF 算術函數**。baseline Spearman r=0.656；控甲基特徵幾乎不動（0.655, wrong-direction collider）；NGroups 99.5%≤2 capped count；within-class 連續 AF 崩塌（Extreme r=0.13, Near-half r=−0.01）；mechanism chain AF→hp_balance 0.164→NGroups 0.190。⇒ **AF→NGroups 報告為 phasing/haplotype-occupancy signature（genetic, AF-mechanical），不當『intermediate-AF subclone 帶不同甲基』之正向甲基發現**。⚠ 完整 partial-correlation knockout 待**一次 ISM re-run** 匯出 `fine_group_counts{HP1,HP1-1,HP2,HP2-1}/variant` + `n_clusters`(甲基-only count) 當乾淨負控（definitional+mechanism-chain+分布四證已一致指 phasing）。[src: research/loh_subclone_af_paired/data/hd4_attribution.json]

**HD-5：umtag (T7) 正式狀態 + 註冊** —— commit fe18c06 over-state（「能真的救 unphase」）、master line 162 under-state（「DESIGN-ONLY 無實測」）、ledger+INDEX 缺席。挑定 canonical（建議 Discussion/future-work）並註冊使三源停止衝突。

**HD-6（資料衛生）：15391 vs 15557 對賬**（branch-sum vs total_selected）+ 永遠標 CHARACTERIZATION set（TP/FP=1.0），**絕不**寫 filter。

> **每個 section 都要防的主 overclaim（FT-1）：把 characterization set（628 / 15391 / strong-ASM）當 filter。** 每條線內部都已撤回，但會悄悄爬回 summary slide。15391 union 的 TP/FP=1.0 —— 可證它不是 filter。

---

## 6. 線的關聯 + 依賴 + 撰寫順序

**依賴**：
- **T1 (phasing)** 脊柱 #1 —— 依賴 V3F/V5/V6 self-phasing C++ fork + priority-bug audit（same-hap bucket 只因該 haplotagging 才存在；R-SELFREF 對照就是 C++ flag-on 重跑 → T1 最終 grade 押在它不擁有的 infra）。
- **T2/T3/T5/T6** 共用 `genome_survey_v2` + copy-partition + CN-confound machinery → **一根支撐柱四個面**（存在 / per-locus cis / 跨樣本 / subclone）。共用 substrate 故彼此不矛盾；聯合產物 =「ASM 真但非-cis-判別、copy-dosage-independent、單樣本天花板」。
- **T4 (filter-dead)** 是「試圖把 T2 的 ASM 變 filter」的副產物 → 它 FEED T2 使其成「乾淨 characterization 非 filter」。最強、最 reviewer-proof 的負結果。
- **T6 的 HPFineNGroups 降級 = T1 獲得 phasing signature 的同一事件**（曾是「subclone 甲基 marker」變成 LOH-constrained phasing signature）。
- **T7 (umtag)** future-work 腿 + 供 V10 copy-control 機制圖（為何甲基能在 LOH 內 phase）。

**撰寫關鍵路徑**：①**先定 HD-1**（R-SELFREF vs reframe，gates 脊柱 grade）→ ②**先寫三~四個 NEGATIVE methods 段**（filter-dead / dosage-refuted / regression-to-extreme）—— 今天就防彈、構成論文最可防守的質量 → ③在 step① 選的 grade 寫 phasing 脊柱 → ④寫 ASM characterization + chr17 exemplar + 跨樣本現象（BRCA2 降為「小-but-真 copy-confounded 例」）→ ⑤改 stale 文件（HD-3）。
**Optional/future-work（不在關鍵路徑）**：umtag 救援、chr17 COLO829 deep-dive 升 ⭐4、6 個 untestable cis、AF→NGroups 甲基-vs-phasing 分解。

---

## 7. 殘存 overclaim 風險 + framing traps（每節防）

- **OR-1**：「n=7 p=0.0078」當強證據 → 是方向一致性 + by-construction 部分循環，n=6 生物獨立。
- **OR-2**：「628 / 15391」當 validated/filter → 單樣本 method-agreement characterization tally；15391 TP/FP=1.0。
- **OR-3**：「~64% partition」當 genome-wide → 7/60 convenience subset、threshold-sensitive（median-of-ratios 實 47%）。
- **OR-4**：「ASM pervasive / 跨癌復現」→ raw median \|Δβ\|<0.10、訊號僅 null 減後浮現、COLO829 最薄、single-pipeline。改「detectable-above-null but small」。
- **OR-6**：chr17/TBC1D16 當「真 cis locus」→ 單 locus 單樣本 n=16/14 nominal p、無第二樣本 → 「稀有但真實 exemplar」。
- **FT-1 characterization-set vs filter**（最高頻）｜**FT-2 cis-anchor vs subclone/copy-confound**（BRCA2 = subclone/copy 主導，HP1-1=subclone tag 非 copy，% 不 robust；且 hypo≠canonical hyper，勿當 TSG-silencing 證據）｜**FT-3 private vs underpowered**（0/38）｜**FT-4 phenomenon-recurrence vs cis-recurrence**（6/6 是 allelic-meth-diff 復現非 cis-ASM 復現）｜**FT-5 held-out simulation vs true-unphase rescue**（umtag 88.5% 是 held-out 模擬非真救援）。

---

## Provenance

- **稽核**：workflow **wf_9e169112-573**（2026-06-08；10 agents = 7 獨立 thread auditor (general-purpose，實讀 ledger/json/report) + 矛盾獵手 + 收斂 + critic）。每個數字 grounded 在 cited 檔（json/tsv/md），auditor fresh-context 不信既有結論自行重判。
- **內部數字來源**：`evidence_ledger.jsonl` + `tsg_promoter_asm_reviewer/genome_survey_v2/*.json`（其他 session 06-07 仍在寫）+ `methyl_augmented_filter_phase2/` + `loh_subclone_af_paired/` + validated master_draft + INDEX。內部 L1-L2 / 外部文獻 L3（見 01-09）。
- 撰寫此 .md 的 Write 與 auditor 檔案讀取**不同 batch**（§13.0）。
- ⚠ tsg 專案仍 active，T2/T3 數字可能再動；投稿前以該專案定稿 + HD-1~6 決議為準。

**相關**：[00 索引](00_index.md) ｜ [08 跨問題鏈](08_cross_question_correlation.md) ｜ [07 ASM-cis（待 HD-3 改）](07_asm_cis_cancer_impact.md) ｜ 視覺版 `10_paper_readiness_convergence.standalone.html`
