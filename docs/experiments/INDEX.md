<!--
建立時間: 2026-03-05 10:00
更新時間: 2026-05-30 (D4 reconcile: 補 5/21 PI V6 signoff / 5/24 HKU handoff / 5/29 ZAR1L ASM；前 header 2026-04-11 為 drift 誤值，內文實含至 5/20)
目標: 實驗研究歷史主索引（精簡版），供 AI 快速掌握研究脈絡
處理範圍: InterSubMod 專案 2025-11 至今的所有研究方向
關聯檔案:
  - docs/experiments/INDEX_DETAIL_ARCHIVE.md（完整詳情）
  - docs/CURRENT_FOCUS.md
  - docs/reports/research_landscape/00_INDEX.md
-->

# InterSubMod 實驗研究索引

> **AI 使用提示**：本索引為精簡版（Layer 1 總覽）。每個條目 1-2 行摘要 + 狀態。
> 完整詳情見 [INDEX_DETAIL_ARCHIVE.md](INDEX_DETAIL_ARCHIVE.md)。
> Live 主線請看 [CURRENT_FOCUS.md](../CURRENT_FOCUS.md)。

## 狀態圖例

- ✅ 成功完成 | ❌ 失敗/無效 | ⏳ 進行中 | 🔄 值得再探索

---

## 最新（2026-06）

### chr2:18M × SEQC2 外部 subclone 比對工具 + demo — 🟡 2026-06-15（F demo / PARTIAL 6 位點）
可重用比對工具 `seqc2_subclone_concordance.py`（loci TSV → SEQC2 sSNV 共識 TVAF/CCF + CNV/LOH + HC 可評估性，3 層 concordance；可套全基因組/COLO829）。**chr2:18M 6 位點 demo 真值**（script 直讀 SEQC2 truth，exit 0）：6/6 落 SEQC2 **LOH** 區（L1 100% 一致）；**5/6** 在 SEQC2 sSNV truth（4 HighConf+1 MedConf）；CCF 梯度 **0.10→0.49→0.77-0.81**（18,099,697 CCF≈0.10=**minor subclone**，吻合 Fang 2021 VAF peak 0.15/0.08）→ 外部 TVAF 獨立佐證此區有亞群 fraction 差異；🔴 樹 pivot **18,086,020 落 HC 空隙=truth-unevaluable**（外部無法確認）。外部僅作 context/sanity（PhyloWGS 單 pipeline，DREAM 19-35% 變異，非金標準）；甲基層無外部檔可比。[demo 報告](in_progress/2026/06/20260615_chr2_18M_seqc2_concordance_demo_01.md) · script+data 在 `in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/`

### 論文就緒收斂稽核 + 外部文獻地景 + T7 umtag 註冊 — ✅ 2026-06-08（meta：跨 7 線獨立稽核）
**(a) 7 角度外部文獻地景**（`knowledge/11_external_literature/` 00-09 + HTML，workflow wf_37b2cc97-663 實際 WebFetch ~38 篇）：我們做法/困難/觀察 × 外部 CONFIRM/REFUTE/DIFFERENT-VIEW + 跨問題依賴鏈。**(b) 論文就緒收斂稽核**（doc 10 + HTML，workflow wf_9e169112-573 7 線獨立 fresh-context auditor 實讀 ledger/json）：論文主軸 = read-level LOH/haplotype+甲基 **characterization+tooling**；主體三~四道防彈 NEGATIVE + phasing 脊柱(**Grade B+ 非 A**) + ASM copy-confounded 支撐。**2 真矛盾**：🔴 BRCA2 cis→**subclone/copy-confounded**(HP1-1=somatic subclone tag 非 copy、非 CN-dosage；% 不 robust；chr17/TBC1D16 唯一乾淨 cis；06-09 R1/R2/R3 已 amend lit-07+master_draft+全相關文件) / Grade A→B+。**決定性 HD-1（用戶決定·hold）**：phasing by-construction 循環 → 跑 R-SELFREF(~25-50hr C++) or 降為 characterization。**(c) T7 umtag 註冊（HD-5）**：methyl-assist phasing / 救 unphase read prototype 之前未進 SoT，本輪 append ledger `20260608_methyl_assisted_phasing_umtag_registration` = FUTURE-WORK/Discussion-only（held-out 0.8852 是模擬非真救援；白地真實但 NOT-READY）。**(d) NEGATIVE backbone 草稿**：[in_progress/2026/06/20260608_G6_methods_negative_backbone_draft_01.md](in_progress/2026/06/20260608_G6_methods_negative_backbone_draft_01.md)。現況地圖（給其他 AI）：[concepts/2026/06/20260608_研究現況地圖_整體目標與流程_給其他AI_01.md](../concepts/2026/06/20260608_研究現況地圖_整體目標與流程_給其他AI_01.md)。⚠ tsg 專案仍 active，T2/T3 數字可能再動。

### 甲基「可分群」非 coverage/CN 假象 + LOH-抑制機制 — ✅ 2026-06-10（confound 驗證 characterization，§N2b）

clusterability confound 驗證（§Methods-Neg §N2b 素材）：甲基 read 可分群性**非 coverage 也非 CN-dosage 假象**——coverage 18 sample×class cells ρ 全 |ρ|<0.25（power floor<50× 後 plateau，多變量控制後 OR=0.91 轉負）；CN 強相關但反 dosage 方向（clustered-rate neutral 0.222>gain 0.179≫**loh 0.022**，LOH logistic OR=0.07，gain 固定覆蓋亦不增）；真 driver=**n_CpG**（OR=1.20，cov+nCpG+CN 僅解釋 8%）。**LOH-抑制分群機制 6/6 樣本 OR<1（0.01–0.107, 3 癌種）= 正向 mechanistic 證據**（read-clustering 需 ≥2 haplotype 等位甲基差）→ 分群是真等位生物學非技術假象。⚠ NGroups(HPFineNGroups)與 coverage ρ=0.51 但屬 HP-tag occupancy 非甲基(HD-4)。FP 不顯示 LOH 抑制(loh n=70 低 power, characterization-only)。[報告](in_progress/2026/06/20260610_methylation_clusterability_vs_coverage_CN_characterization_01.md) · scripts 101/102 · figs `04_figures/clusterability_vs_cov_cn{,_extended}.png` · ledger 97/98 · 單樣本 CN truth ⭐3

### ASM 甲基分群 characterization：regime-A credible 子集 + LOH 成因診斷 — ⏳ 2026-06-02（HCC1395 單樣本 Tier-A）
G1 + G2 雙目標 consolidation（既有 scripts 15-35 約 70% 已做，本輪補決定性 gate + 全新 LOH 診斷）。**G1**：somatic-controlled HP-axis 的 credible regime（高cov+中baseline+nonLOH）子集**有真 somatic 甲基 cluster** — 23/62 過 length-placebo cut、**15/23 撐過完整 collider battery**（M4c 96% + 強化 random-anchor M8 70%），regime-A blind-ARI **顯著高於 germline-het baseline**（0.229 vs 0.038, Cliff δ=0.37, p=2.3e-4, M3 rarefied 0.308 vs 0.080）。**refine（非推翻）5/31**：全基因組仍 artifact-dominated，乾淨子集有真訊號。⚠ **強化 collider 刷掉 SOX2/HOTTIP/SDHAF1**（genomic-context collider）；BRCA2 在 regime-A 內但未過 placebo cut（borderline）；**單樣本 Tier-A 天花板，M9 cross-sample 不可能 → 不可 generalize**。**G2**：LOH 表觀雙 haplotype（longphase-S 不驗證 LOH 而生）診斷分類 = **72% self-phasing artifact / 18% candidate subclone / 9% CN-regression**（驗證用戶盲點假設；診斷非強制歸類，強制歸類=V5 Layer 1.5 已失敗）。filter-NEGATIVE 未觸碰（characterization-only）。[報告](in_progress/2026/06/20260602_ASM_cluster_characterization_regimeA_LOH_diagnostic_01.md) · 4 figs (fig1-4) · ledger `ASM-cluster-characterization-20260602`

---

## 基礎設施（2025-11 ~ 2025-12）

### 01. 甲基化解析（MM/ML 標籤）— ✅
ONT BAM MM/ML 標籤解碼 + CIGAR 座標校正正確。[詳情](INDEX_DETAIL_ARCHIVE.md#01-甲基化解析mml-標籤)

### 02. CIGAR 座標映射 — ✅
反向股 MM delta encoding 解碼正確。[詳情](INDEX_DETAIL_ARCHIVE.md#02-cigar-座標映射)

### 03. 距離計算與聚類分析 — ✅
NHD/L1/L2/CORR/BERNOULLI/JACCARD 六種距離 + UPGMA 聚類。[詳情](INDEX_DETAIL_ARCHIVE.md#03-距離計算與聚類分析)

### 04. OpenMP 平行化 — ✅
Per-region 平行化，單 Region < 300ms。[詳情](INDEX_DETAIL_ARCHIVE.md#04-openmp-平行化)

### 05. 統計顯著性：Fisher / χ² — ✅
Fisher-Freeman-Halton RxC + Monte Carlo early stopping。[詳情](INDEX_DETAIL_ARCHIVE.md#05-統計顯著性fisher--χ²)

### 06. 統計顯著性：PERMANOVA — ✅
NHD 距離矩陣 + 999 permutations + Cramér's V 效應量。[詳情](INDEX_DETAIL_ARCHIVE.md#06-統計顯著性permanova)

### 07. Bernoulli 距離度量 — ✅
Bernoulli probability 距離補充 NHD，對稀疏 CpG 更穩定。[詳情](INDEX_DETAIL_ARCHIVE.md#07-bernoulli-距離度量)

---

## 特徵分析與過濾策略（2025-12 ~ 2026-02）

### 08. TP/FP 特徵富集分析 — ✅
CramersV>0.15 + PairwiseMedianDist>0.10 初步有效。[詳情](INDEX_DETAIL_ARCHIVE.md#08-tpfp-特徵富集分析)

### 09. 雙向驗證（Label-First vs Cluster-First）— ✅
兩路線結論一致性 > 80%，可互補但非替代。[詳情](INDEX_DETAIL_ARCHIVE.md#09-雙向驗證label-first-vs-cluster-first)

### 10. F1 最佳化 — ✅
最佳組合 F1 +0.001 級別；甲基是 support 非 primary。[詳情](INDEX_DETAIL_ARCHIVE.md#10-f1-最佳化)

### 11. Subsample 混樣甲基化偏差 — ❌
purity-aware 受 tumor-normal 組織差異混淆。[詳情](INDEX_DETAIL_ARCHIVE.md#11-subsample-混樣甲基化偏差)

### 12. Purity-Aware 過濾策略 — ❌
subsample 無法可靠模擬純度效應。[詳情](INDEX_DETAIL_ARCHIVE.md#12-purity-aware-過濾策略)

---

## 純樣本基線研究（2026-03-07 ~ 2026-03-12）

### 13. 多樣本甲基化距離基準線 — ✅ 2026-03
7 個 pure paired 基線建立。[詳情](INDEX_DETAIL_ARCHIVE.md#13-多樣本甲基化距離基準線2026-03)

### 14. no-filter 模式 Segfault 修正 — ✅ 2026-03
NumReads=0 guard + empty matrix check。[詳情](INDEX_DETAIL_ARCHIVE.md#14-no-filter-模式-segmentation-fault-修正2026-03)

### 15. 全樣本放寬標準拆解與 F1 提升驗證 — ✅ 2026-03
GQ>6 跨 7 樣本 F1 中位升幅 +0.002。[詳情](INDEX_DETAIL_ARCHIVE.md#15-全樣本放寬標準拆解與-f1-提升驗證2026-03)

### 16-23. 純樣本 round 自動化與特徵分析 — ✅ 2026-03
TO pilot 完成（F1=0.7130, +0.0003）；低 VAF+高 AlleleDelta 僅 HCC1395 有效，不可全域化。[詳情](INDEX_DETAIL_ARCHIVE.md#16-純樣本-round-自動化實際執行2026-03-07)

### 24-28. ClairS 邊緣 TP rescue 與甲基輔助 — ✅ 2026-03
Caller-first（GQ≥15/20）仍最佳；甲基 rescue 不超過 caller-only。[詳情](INDEX_DETAIL_ARCHIVE.md#24-clairs-邊緣-tp-rescue-與甲基輔助評估2026-03-08)

### 29-33. GQ 與甲基 rescue 系統性驗證 — ✅ 2026-03
GQ 最穩定；AlleleDelta 適合 artifact triage 非 rescue 主規則。[詳情](INDEX_DETAIL_ARCHIVE.md#29-gq-與甲基-rescue-的系統性驗證2026-03-11)

### 34-37. TO 雙模型、scope control、Pool B 收尾 — ✅ 2026-03
pileup-route only 確認；same-scope control 18/18 identical；Pool B FN caller-first 最佳。[詳情](INDEX_DETAIL_ARCHIVE.md#34-to-雙模型可得性snapshot-scope-與-pool-b-fn-收尾盤點2026-03-12)

### 38-39. TO FP 來源分解與 residual FP — ✅ 2026-03
PoN 移除率 99.48%；persistent FP 只有 ~80 個/樣本；ISM 無法再改善。[詳情](INDEX_DETAIL_ARCHIVE.md#38-to-fp-來源分解與-paired-對照分析2026-03-22)

---

## 方法學審查與研究方向定調（2026-03-24）

### 40. 研究方法與文獻突破方向全域分析 — ✅
E→A+D→B→C 優先序確認；ISM 獨特定位確認。[詳情](INDEX_DETAIL_ARCHIVE.md#40-研究方法與相關文獻突破方向全域分析2026-03-24)

---

## Phase 1A ML Read Classification（2026-03-25 ~ 2026-03-28）

### 41-43. Phase 1 研究啟動、schema、manifest — ✅
141,014 regions manifest；discovery/external split 70K/70K。[詳情](INDEX_DETAIL_ARCHIVE.md#41-phase-1-ml-read-classification-研究啟動與資料缺口盤點2026-03-25)

### 44. Phase 1A benchmark round 1 — ✅
context-rich 模型穩定；pure methyl 不足；to-pure 是主要困難。[詳情](INDEX_DETAIL_ARCHIVE.md#44-phase-1a-read-classifier-benchmark-round-12026-03-25)

### 45. Phase 1A round 2 multi-bio validation — ✅
**paired-pure delta F1=+0.0112, CI=[+0.0044,+0.0188]** 已鎖定。[詳情](INDEX_DETAIL_ARCHIVE.md#45-phase-1a-incremental-test-與-multi-bio-validation-round-22026-03-25)

### 46-47. LOH Round 1-4 cross-sample audit — ✅ 2026-03-27
LOH evidence panel 完成；推薦 feature `tier_aplus_loh`（enrichment=2.02×）。[詳情](INDEX_DETAIL_ARCHIVE.md#46-loh-round-1-cross-sample-audit-啟動與決策紀錄2026-03-27)

---

## 系統性觀察（2026-03-31 ~ 2026-04-02）

### O1-O10. 系統性大規模觀察 — ✅
82 圖表 × 9 主題。**TO 所有 AUC<0.58**；LOH penalty 是 QS 根因；Paired/TO 根本不同。[詳情](INDEX_DETAIL_ARCHIVE.md#29-系統性大規模數據觀察-o1-o102026-03-31)

### O11. Within-Group Methylation Heterogeneity — ❌ NEGATIVE
epipolymorphism AUC=0.845 完全是 n_reads confound（殘差化後 0.530）。[詳情](INDEX_DETAIL_ARCHIVE.md#30-o11-within-group-methylation-heterogeneity2026-03-31)

### O12. TO LOH Methylation Scenario Discrimination — ❌ NEGATIVE
AlleleDelta 是 AF confound；CramersV L2 AUC=0.80 是 collider bias。[詳情](INDEX_DETAIL_ARCHIVE.md#31-o12-to-loh-methylation-scenario-discrimination2026-03-31)

### O13. Cross-Region Methylation Correlation — ❌ NEGATIVE
FP>TP correlation 是 shared read count confound；分層+residualize+matching 後消失。[詳情](INDEX_DETAIL_ARCHIVE.md#32-o13-cross-region-methylation-correlation2026-03-31--2026-04-01)

### ASM 定量. SNV-甲基化關聯性 — ✅ POSITIVE
32-66% SNV 位點有 ASM；FP>TP 但重疊大；ISM PERMANOVA 唯一捕捉 entropy imbalance。[詳情](INDEX_DETAIL_ARCHIVE.md#33-snv-甲基化關聯性定量分析2026-04-01)

### G1-G7. TO Germline FP 鑑別 — ❌ NO-GO
48 圖表 × 60+ 特徵全 AUC<0.64；**TP loss ≤2% 下 FP removal=0%**。[詳情](INDEX_DETAIL_ARCHIVE.md#34-to-germline-fp-鑑別研究2026-04-01)

### Read-Level. Haplotype & Methylation 鑑別 — ⚠️ CONDITIONAL NO-GO
LOSO AUC=0.721（首次>0.70）但安全約束 FP removal=0%；根因=高純度細胞株。[詳情](INDEX_DETAIL_ARCHIVE.md#35-read-level-haplotype--methylation-germline-fp-鑑別2026-04-02)

### O14. LOH Haplotag Bias — ✅
99.5% HP_Ratio 極端偏向；TO vs Paired 方向一致性 52.9%（≈隨機）。[詳情](INDEX_DETAIL_ARCHIVE.md#36-o14-loh-haplotag-bias-觀察2026-04-02)

### O15. LOH 區域 ISM 指標全面量化 — ✅
32 圖 + 2 報告。LOH 內甲基化全面失效（7/7 samples AUC~0.50）；只有 caller 特徵保留區分力。[詳情](INDEX_DETAIL_ARCHIVE.md#37-o15-loh-區域-ism-指標全面量化2026-04-02)

---

## 因果鏈驗證與修正（2026-04-02 ~ 2026-04-06）

### Self-Phasing 循環依賴因果鏈 — ✅ CONFIRMED
完整五步因果鏈驗證：62% LOH 消失（d=-1.20）、somatic bias 17.3:1→消除、31.2% self-phasing LOH。23 TSV + 7 PNG。[詳情](../reports/validated/2026/04/) [肉眼觀察](in_progress/2026/04/20260403_multilayer_HP_before_after_肉眼觀察.md)

### VCF 來源錯誤矯正與 V3-Fixed 重分析 — ✅ 2026-04-04
pileup symlink 指向 ClairS paired(ssrs) 非 ClairS-TO；TP=30490/FP=4842 全為 paired 結果；正確 TO=28509/11606；V3F ISM F1=0.0124 非 0.0945。多層 HP 驗證 + ClairS-TO VCF 矯正完成。[VCF 矯正](in_progress/2026/04/20260404_VCF來源錯誤矯正報告_01.md) [ClairSTO 重分析](in_progress/2026/04/20260404_ClairSTO_VCF矯正全面重分析報告_01.md) [LongPhase v2b 驗證](in_progress/2026/04/20260404_LongPhase_TO_v2b_多層HP_ISM驗證完整報告_01.md)

### PON-Only Phasing 驗證 — ✅
LOH.bed Jaccard=1.0（不變）；somatic bias 消除；N50 +99.7%；phased rate +23.6pp。待 haplotag+ISM 全量重跑。

### TO Feature Deep Study Q1-Q6 — ✅ 2026-04-04~05
HP2FamilyN AUC 0.72 是循環 artifact（校正後 0.54）；RF ceiling 0.69-0.77；baseline>v2b in LOH；甲基化特徵全無效。[Q1-Q6 深度研究](in_progress/2026/04/20260405_TO_ClairSTO特徵區分力深度研究_01.md) [LOH 內外觀察](in_progress/2026/04/20260405_TO_LOH內外特徵觀察報告_01.md) [HPFineP+QS](in_progress/2026/04/20260404_HPFineP_QS整合完整研究報告_01.md)

### QS Mode-Aware 修正 — ✅
TO 模式下停用 LOH penalty 與 verify bonus；C++ 已實作並部署。

### LOH 雙定義交叉分析 Wave 1+2+3 — ✅ 確定性結論
W1+2: HP Imbalance 是 LOH.bed 超集 (Sens 99.7%); LOH.bed vs SEQC2 Jaccard=0.928; HP PERMANOVA LOH 內不可用 (valid 5%); **LOH 不可作為 variant filter** (10/10 FAIL); AlleleDelta AUC=0.556 真實但微弱。W3: **Non-LOH 區分力同樣有限** (max AUC=0.643 read count proxy); **多特徵組合不可行** (Voting AUC=0.577<0.58); **cnLOH PairwiseMeanDist 0.587 是 Simpson's Paradox** (per-sample mean=0.50); **CramersV 被 NumReads confound** (0.511→0.464); AlleleDelta 是 LOH 唯一 confound-free 信號但不足。40 位點肉眼檢視文件+120 圖表已生成。[完整報告](../reports/validated/2026/04/20260406_LOH雙定義交叉分析報告/00_INDEX.md) [AlleleDelta 深度](in_progress/2026/04/20260404_LOH_AlleleDelta_深度驗證報告_01.md) [Strong/Weak AUC](in_progress/2026/04/20260404_LOH_Strong_Weak_7feature_AUC驗證報告_01.md)

### R1-R5. 特徵設計局限與改進方向研究 — ✅ 2026-04-07
CramersV 93% 為零 = 2×2 缺陷（R1）；Excess groups 無新信號（R2）；結構清楚子集 AUC 下降 = identifiability 非設計問題（R3）；N≥4+NR≥80 TP 89.1%（R4）；PwDist 正交但太弱（R5）。**HPFineNGroups 確認為 somatic heterogeneity 標記**。正負效果+驗證標記已寫入報告。

### 肉眼檢視推理鏈與 TP/FP 可區分性 — ✅ 2026-04-06~07
40 站點（4 類 × TP/FP × 5）系統性肉眼比對 + 全量 7 樣本定量分析。LOH 結構反向確認；Non-LOH ALT-HP alignment 群體不可分；F1 影響分析。完整報告含 R1-R5。[完整報告](../reports/validated/2026/04/20260406_肉眼檢視推理鏈與TP_FP可區分性分析_01.md)

### Option C: LabelTest 雙路測試（HP-free vs HP-dependent）— ❌ 2026-04-07
架構調查確認 cluster_labels 已是 HP-free。ClusterPermanovaF AUC=0.512（隨機）；HP-free 5 features combo AUC=0.564 vs HP-dependent 0.598 vs 全部 0.607。HP-free 僅增加 +0.009。**純甲基化 clustering 無區分力，identifiability problem 再次確認。C++ 修改取消。**

### O9: FN 特徵觀察 — ❌ 2026-04-08
7 samples × 2 modes (Paired+TO)，122,790 FN regions 完整 ISM 執行。HP-free 甲基化特徵全 AUC<0.53（random）；最強信號 LabelAllelePermanovaF=0.664 是 AF 代理非甲基化；TO QS AUC=0.338（FN>TP 反轉，Cohen's d=-0.671）。**NO-GO: ISM 無法 rescue FN，甲基化空間 FN≡TP。**

### TO-pure 獨立建模 — ❌ 2026-04-08
7 samples LOSO 四組模型：HP-free AUC=0.53（random）、All-ISM AUC=0.60-0.64、Caller-only AUC=0.63、ISM+Caller AUC=0.66。ISM 僅在 Caller 基礎上增加 +0.003~+0.030。**Caller_af (0.654) 獨自超越全部 ISM 特徵，TO 模式 ISM 近乎無效。**

### Fine-Pairwise Distance 分析 — ❌ 2026-04-08
HP 四群組（HP1/HP1-1/HP2/HP2-1）6 組 pairwise 距離，748,391 regions。**Paired 全 AUC<0.50（反轉 — germline ASM > somatic ASM）**；TO 最高 0.579（<0.58 門檻）。LOH 層 Paired 極端反轉（0.132）。Group count 確認 self-phasing 指紋。多 agent 驗證通過。**ISM 甲基化特徵空間已耗盡。**[報告](finalized/2026/04/20260408_fine_pairwise_distance_analysis_01.md)

### Beyond-AUC 7 方法綜合驗證 — ❌ EXHAUSTED CONFIRMED 2026-04-09
7 種互補統計方法（PR-AUC, Residualization, Subgroup Discovery, Calibration, Bootstrap, GradientBoosting, Distribution/Tail）× 多層分層（LOH×AF×Sample）× 25 特徵 × 748K regions 全面驗證。**Pure methylation 特徵全部 AUC ≤ 0.58**（within-group residualized 後更低）；pooled OLS residualization 被揭露為 data snooping（confounder-only AUC 0.66-0.85 完全解釋 residualized 提升）；PR-AUC 與 ROC-AUC 排序一致（Richardson 2024 確認 AUC 穩健）；BSS 全為負值。唯一正面：HPFineNGroups TO non-LOH low AF AUC=0.72（7/7 samples consistent）確認為 somatic heterogeneity marker。45+ 篇論文支持：germline ASM >> somatic passenger SNV 效應（3-6×）。**ISM 甲基化特徵空間正式關閉。**[報告](finalized/2026/04/20260409_beyond_auc_comprehensive_validation_01.md)

### SEQC2 CNV 分層觀察 — ❌ CNV zone-aware filter 關閉 2026-04-10
SEQC2 正交驗證 CNV truth set（6 callers × 21 replicates × 3 technologies）分層分析。**Phase 1 (HCC1395)**：Gain+LOH 集中 57.5% FP（FP rate=10.2%, 全域 2.6×）；zone-specific AUC 達 0.782（AlleleDelta）；Coverage_Multiple vs SEQC2 CN r=0.831 代理可信。**Phase 2 (7 樣本跨樣本驗證)**：FP rate 模式跨樣本不一致（CN_HighGain > CN_Normal 僅 4/7）；per-sample zone-specific mean AUC ≤ 0.641 未突破上限；Simpson's Paradox 否定（Quality_Score 反向 -0.042）。**Phase 3 (根因)**：Gain+LOH 內 CN=3 FP rate 12.9% 最高（FP rate 隨 CN 增加下降）— 中等 CN gain + LOH 的 allele 不平衡環境是樣本特異性根因。**結論：CNV 不是特徵空間耗盡的根因；zone 排除策略 trade-off 全不可行；CNV zone-aware filter 正式關閉。**15 張圖 + 5 TSV。[報告](in_progress/2026/04/20260409_SEQC2_CNV分層觀察_01.md)

### Coverage_Multiple GC 校正與甲基化-CN 驗證 — ❌ 全 NO-GO 2026-04-11
TO Pipeline（ClairS-TO → LongPhase-TO → ISM, TP=28,383, FP=11,830）三方向驗證。**GC-Content 校正**：delta-r=-0.0002（門檻 ≥0.03）；ONT 5kHz GC bias 極小（98.7% regions 變化<5%）；TP/FP AUC 不變（0.5095→0.5097）。**甲基化-CN 相關**：所有 HP-free 特徵 residualized |r|<0.07（CN-blind）；HPFineNGroups 68% 是 NumReads confound（raw 0.495→resid 0.160）；CramersV resid r=-0.726 是零值 artifact。**KDE 已實作確認正確**：CN 分類準確度 6.2%→43.8%（`--expected-coverage` CLI 可用）。**結論：GC 校正不需實作；甲基化無法驗證 CN；Coverage_Multiple 現有精度已足夠。**6 張圖 + 1 TSV。[腳本](../../scripts/analysis/gc_correction_to_validation.py)

### PON 跨樣本 FP 移除率驗證 — ✅ 穩定度升級 2026-04-11
7 樣本 raw VCF-level PON rate 95.16%-98.81%（mean 97.77%, SD 1.12%）。Refined FP-level 全 > 98%。H2009 最低（95.16%）因高突變負荷（168K TP）非 PON 失效。HCC1954 refined 98.44% 最低但仍遠超門檻。**結論穩定度 3→4：「PON 是最強 filter」跨樣本確認。**[報告](in_progress/2026/04/20260411_PON跨樣本移除率驗證_01.md)

### H2009 Phase 1A 負向根因診斷 — ✅ 根因確認 2026-04-11
H2009 paired FP rate 僅 0.06%（86/132,994）— 7 樣本中最低之一。改進天花板 +0.0004 F1。76.7% FP 在 LOH 區域、89.5% Noise class。ISM 特徵區分力反而優於平均（QS gap +25 vs others +2.9）。**根因=caller 已近乎完美，ISM 只能誤傷 TP，非甲基化方法學問題。** H2009 可在論文中定位為「caller 已極佳時 ISM 無增量價值」的示範。[報告](in_progress/2026/04/20260411_H2009負向根因診斷_01.md)

### 文獻假說交叉驗證（L1-L4）— ❌ 大部分 NEGATIVE 2026-04-12
60+ 篇文獻 4 大假說 vs 340K 區域 × 7 樣本實證驗證。**L1 方向性 ASM (epiTRACERx)**: TP signed_delta=-0.0003（p=0.854, 完全隨機），FP 有方向性 (p=1.16e-14)。**L2 PMD 分層**: CpG variance PMD>Non-PMD confirmed (p<1e-100) 但 TP/FP 判別 AUC 不超過 0.622。**L3 Normal baseline**: 7 Normal BAM 均含 MM/ML 標籤，**Phase 2A 可行**。**L4 fCpG (EVOFLUx)**: TP/FP CpG variance 完全相同 (p=0.77)。文獻正面結果源於任務差異（tumor vs normal 細胞區分 ≠ TP vs FP variant 判別）。[綜合報告](../../research/literature_validation/reports/20260412_文獻驗證綜合報告_01.md) [L1 詳細](../../research/literature_validation/reports/20260412_L1_directional_ASM_report.md)

### Phase B/C/D Dual-BAM 架構驗證 — ✅ PASS 2026-04-13
HCC1395 paired (31K TP + 1.3K FP) 全基因體驗證。Phase B: Sample ASM 97.3% 顯著、Normal Baseline 100% 有效。Phase C: LOH BED 94.1% concordance with ISM hp_ratio。Phase D: 4-group subclone 分層（Normal Diploid 17.5%/Epi. Heterogeneity 12.9%/LOH 2.6%/Tumor-Specific 67.0%）。173/173 unit tests pass。[報告](validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md)

### LOH Subclone AF × Methylation 雙重證據鏈 — ✅ POSITIVE 2026-04-14 + Paired Extension 2026-04-15
LOH 區域 intermediate AF (0.1-0.4/0.6-0.9) 的 variants 具有顯著更高的甲基化多樣性。**TO mode**: Intermediate vs Extreme NGroups: 1.796 vs 1.091（+0.705, **7/7 p<10^-39**），NR 控制後持續（r=0.48-0.71）；Segment ρ=0.270（6/7 positive）。**Paired mode 延伸驗證**: ΔNG=+0.787（7/7 p<10^-65），效應量更強（median |r|=**0.755** vs TO 0.630），segment rho=0.382；**4/4 假說全部 POSITIVE，7/7 效應方向跨模式一致**。構成跨模式確認的 subclonal LOH 雙重證據鏈（genetic AF + epigenetic ASM）。[報告](validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md)

### Per-CpG ASM + Epiallele 異質性指標 — ⚠️ CONDITIONAL POSITIVE 2026-04-15
30+ 方法文獻調查 → 6 指標家族 24 metrics Python PoC（Fisher, NME, Epipolymorphism, PDR, Shannon, VEF）。**Filter value=0**（全 AUC 0.38-0.54，與 Beyond-AUC 一致）。**Characterization value 確認**：PERMANOVA concordance 84%；LOH 降 fisher_frac_sig 12.5pp；NGroups=1→0.001/NGroups=3→0.270。PDR/VEF 飽和（0.917/0.918）。推薦精簡 C++ 整合 10 欄位（Fisher 4 + NME 3 + Epipoly 3）。[報告](in_progress/2026/04/20260415_PerCpG_ASM_Epiallele_Metrics_01.md) [文獻](../../references/20260415_ASM_subclone_methods_literature_survey.md)

### LOH/CN/AF 結論驗證 — ✅ 5/6 confirmed, 1 needs correction 2026-04-17
6 個關鍵結論系統化驗證：Q1 cnLOH Simpson's Paradox（確認），Q2 HPFineNGroups CN confound 68%（確認），Q3 Coverage_Multiple CN proxy r=0.997（確認），Q4 LOH Subclone DNGroups 7/7 positive（確認），Q5 Zone exclusion abs ratio 0.082（確認），**Q6 LOH Tier direction reversal（需更正：post-HP-fix ALL tiers TP-enriched，0.43×/2.018× 為 artifact）**。[報告](../../research/loh_cn_af_verification/20260417_LOH_CN_AF_結論驗證報告_01.md) [總整理](../reports/research_landscape/07_LOH_CN_AF_研究總整理.md)

### Zone-Aware Confidence Framework — ⚠️ Characterization POSITIVE, F1 NEGATIVE 2026-04-17
**H1/H3 驗證**：5 zone（Z1-Z5）全量驗證。H1 Z1 方向確認但覆蓋率不足 → Z1b（NGroups≥2）放寬後 TO 4.6% 覆蓋率、TP rate 0.965。H3 Paired 7/7 ≥ 89.1% 確認；TO 不成立（mean 0.716）但 6/7 significant。**Z1 放寬分析**：7 變體 Pareto 最佳為 Z1b。**QS 模擬 ❌ NEGATIVE**：5 delta configs × 21 thresholds × 7 samples，max delta F1=+0.001。根因：TO QS AUC=0.497 隨機，zone 調整無法修正。**結論：Zone-Aware 價值僅在 characterization annotation，不在 F1 改進。** [H1/H3 報告](../../research/zone_aware_validation/20260417_H1_H3_Zone_Validation_Report_01.md) [QS 模擬](../../research/zone_aware_validation/20260417_QS_Simulation_Report_01.md) [Framework](../concepts/2026/04/20260417_Zone_Aware_Confidence_Framework_01.md)

### P3 Window Aggregation Pilot — ❌ NEGATIVE 2026-04-17
1 Mb window aggregation pilot 驗證 gene-level 方向。Naive 結果 4/7 樣本看似突破（H2009 ΔAUC=+0.342），但 mid-TP-rate window confound check 後反轉：H2009 Δ=-0.346，8 個 TO 測試中僅 1 個（H1437 AlleleDelta, Δ_mid=+0.097）保留明確 gain。**結論：window aggregation AUC 幾乎完全由 TP/FP 基因組空間 auto-correlation 驅動。** Phase 2B gene-level 若要推進必須以 shuffle-within-chr null model + mid-TP-rate Δ>+0.03 雙門檻驗收。新增 pitfall 規則於 memory。[報告](in_progress/2026/04/20260417_P3_window_aggregation_pilot_NEGATIVE_01.md) [Pilot 數據](../../research/P3_window_aggregation_pilot/)

### P4 Second-Hit Order Pilot — ⚠️ CONDITIONAL NEGATIVE 2026-04-17
ISM 目標 3（二次打擊順序）在 region-level AF × methylation 聯合分佈是否可區分。Intermediate AF (0.3-0.5) LOH vs NonLOH：**Bimodality Coefficient mean |Δ|=0.043**（遠低於 0.15 門檻），AlleleDelta 分佈不顯現 biallelic/monoallelic 分離。HPFineNGroups Δ=-0.637 (7/6 樣本一致) 但為 LOH 機制副作用（單 haplotype）非 order 訊號。**結論：單 region 特徵 summary 無法推斷 two-hit order；需 per-read epigenotype（ISM 目標 1）。目標 3 依賴目標 1，應暫緩等目標 1 進度。** [報告](in_progress/2026/04/20260417_P4_second_hit_order_pilot_CONDITIONAL_NEG_01.md) [Pilot 數據](../../research/P4_second_hit_order_pilot/)

### F HPFineNGroups Deepening Pilot — ✅ POSITIVE (refined) 2026-04-18
Part B.1 對 HPFineNGroups POSITIVE 結論 4 項質疑（residualized AUC、飽和、7/7 強度、LOH/AF 混淆）三階段驗證。**三大新發現**：(1) **NGroups 非單調**（NG=2 TP rate 0.643<NG=1 0.763）根因為 germline AF confound（NG=2 AF mean 0.471 最接近 0.5），(2) **HCC1954 失效根因 = FP 在 AF≥0.4 極端富集**（AF<0.2 TP=0.874 正常→AF 0.8-1.0 TP=0.022），(3) **Paired 99.85% 非真 gain**（baseline 98.96%）。**新 canonical filter: NG=4+AF<0.4+NR≥80 NonLOH TP rate=0.9281（舊 0.8912 +3.7pp）**；HCC1954 +21.0pp、HCC1937 +15.4pp 挽救；5/7 樣本 ≥0.85（舊 4/7）；AF<0.2 則 6/7≥0.87。**Confound checks PASS**：chr-shuffle null Z=43.5（套用 P3 教訓）、Coverage_Multiple 跨 CN tiers 0.90-0.94 穩定。**COLO829 out-of-scope**（ONT_R10 無 methylation basecall）。**結論穩定度補充結論 16 升級 ⭐3 → ⭐4**。[報告](in_progress/2026/04/20260418_F_HPFineNGroups_deepening_POSITIVE_01.md) [Pilot 數據](../../research/F_hpfinengroups_deepening/)

### CN KDE 雙 Pass 校準 — ✅ POSITIVE (HCC1395 pilot) 2026-04-20
commits `374fad4` + `12d9b3e`：`--expected-coverage` 從 hardcoded 75.0 default 升級為 per-sample KDE 2× diploid peak。HCC1395 pilot KDE estimate = **53.0×**（SEQC2 neutral median 54× ground truth，bias −1.9%），舊 75× bias +39%。CovM median 0.880 → 1.245（×1.415 恰為 75/53 等比例拉伸）；Coverage_Category 重分類 Normal −5,718 / CNV_Gain +2,956 / High_Copy +2,710。KDE fallback paths 5 種明確記錄 + audit TSV column `Diploid_Coverage_Used`；unit tests 202/202 通過。影響邊界：影響跨樣本 CovM 絕對值比較、CovM percentile filter；**不影響** per-sample 內部排序、LOH.bed（Jaccard=1.0）、percentile-based filter（Variant A 證明 scale-invariant）。**結論**：`--expected-coverage 75.0` 共用為 stale binary artifact；端到端驗收 PASS bias <2%；master dataset 7×2 modes 待全量重跑。[驗收報告](in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md)

### `--germline-hp-only` flag Phase 1 — ⚠️ CONDITIONAL NEGATIVE on filter 2026-04-21
commit `775036c` 實作 ReadParser `--germline-hp-only` flag（PON 當 germline 錨點，somatic HP demote）。HCC1395 TO 全量 40,115 sites × 4 runs 驗證。**H-C1 機制 ✓**：Audit 獨立性 NHP_Somatic11/21/33 全基因 sum 在 TP+FP 合併下兩 flag identical，排除 flag 誤 remove reads 替代解釋。**H-C2 Filter Gate FAIL**：18 個 HP-related 特徵無一達 +0.02 Gate（最大 HP1FamilyN +0.0084），4 個 HP-derived（HPFineNGroups/NHP3/HPMergedDelta/HPFine_NGroups_CF）反向下降 ≤−0.025；AlleleDelta (HP-independent 對照) AUC=0.6294 不動符合機制預期。**H-C3 Paired 一致性**：本週未完成對照（需 Phase 2B）。**flag=on 下 NG≥3 regions=0**（TP 22,732+FP 8,148 → 0，97.3% 歸 NG=2）→ 揭露 HPFineN marker 訊號源至少部分依賴 somatic HP tag。**結論**：機制保留為 characterization 工具（default=off），filter 方向懸掛 Phase 2B 判定；不進 Phase 2 全量重跑。Plan R3 「HP_Ratio 歸零」條件未觸發，bug 不在上游 hp_tag 解析。[Phase 1 報告](in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md) [Marker re-audit](in_progress/2026/04/20260421_HPFineNGroups_Marker_Reaudit_01.md)

### LOH-constrained phasing — n=7 跨樣本 Wilcoxon ✅ POSITIVE（Grade A 要求 #1 達成）2026-05-30
論文主軸 Grade B→A 升級。NG=2 Inner same_HP1 vs Outer cross_het TP-rate gap 全 **7/7 樣本一致正向**，Wilcoxon **W=28 p=0.0078**（exact, = 1/2⁷ 理論最小；從 n=6 W=21 p=0.0156 強化）。加入 COLO829（gap +0.104，Inner 0.724/Outer 0.620）— 親自驗證 `20260423_colo829_to_pilot` 真實 KDE-corrected ISM 存在（先前誤判 blocked 源於 2026-03-17 dry-run stub）。**便宜 Python aggregation**（COLO829 ISM 已存在，非 100GB full pilot）。**provenance fix**：原 B1 Wilcoxon 腳本從未 commit → 寫 `obs18b_wilcoxon_ng2_gap.py`（重現 B1 reproduces=True）。Grade A 要求 #2（R-SELFREF 全 7-sample flag-on 負控）仍部分（X3 HCC1395 + V6C 4-sample chr19 機制已證，全版需 ~25-50hr C++ 重跑）→ 整體 **Grade B+ 接近 A**。[報告](in_progress/2026/05/20260530_LOH_phasing_n7_cross_sample_Grade_A.md) [data](../../research/tpfp_loh_af_kde_discrimination/data/obs18b_wilcoxon_ng2_gap_n7.json) [ledger 20260530_loh_phasing_n7_cross_sample_wilcoxon]

### ISM ZAR1L/BRCA2 ASM — ⭐3 方向 POSITIVE + magnitude RESOLVED + B-discrimination anti-discriminative 2026-05-29（5/30 全 TP 重跑定案）
既有 artifact 驗證 (tsg_promoter_asm_reviewer) + Phase 2 全基因組擴展。BRCA2/ZAR1L bidirectional promoter (chr13:32,315,128 G>A, TVAF 0.189) ASM **方向真實穩健** (hypo, 4 切法一致, 抗 Bonferroni)。✅ **magnitude RESOLVED (5/30 全 39,447 TP max-collapse 修正口徑重跑, 已驗證 v2 TSV)**：**BRCA2 HP-axis Δβ=−0.122 (n=197, p=6.62e-10) / ALLELE −0.099 (n=198)** — 與 script03 5mC-only −0.122 收斂；buggy −0.054 證實為 MSA Level1 對 5mCG_5hmCG BAM 每 (read,CpG) 發 5mC+5hmC 兩列各獨立計數把 beta 砍半的 artifact (非 CpG-set)。修法已落地 Python `18_dual_axis_pivot.py` MAX-collapse (C++ 根因在外部 repo 未動)。**全基因組 A-characterization**: 訊號漏斗 51,091→Bonferroni 313→strong-ASM (Bonf AND |Δβ|≥0.1) **172/51,171 (0.34%, ~26 獨立位點)**；hypo 44%/hyper 56% 無偏好 (binomial p=0.147)；dual-axis concordance 77%。**B-discrimination NEGATIVE 且 anti-discriminative**: strong-ASM 在 FP enriched 5× (OR=0.194, p=1.8e-28) → ASM 真實但不能判別 TP/FP。**單樣本 HCC1395** (⭐4 需 COLO829)。[報告](in_progress/2026/05/20260529_ISM_ZAR1L_BRCA2_ASM_verification_01.md)（5/30 RESOLVED banner）[synthesis](../../research/tsg_promoter_asm_reviewer/output/phase2_dataverify_synthesis.md) [ledger 20260529 + 20260530_correction + 20260530_full_tp_characterization]

### ISM 輔助 tag — 甲基 × haplotype × somatic 觀察推論漏斗 ⭐3 NEGATIVE 整合（5/29→5/31）2026-05-31
ZAR1L/BRCA2 ASM 調查全鏈收斂為**觀察推論漏斗**：39,447 SEQC2 高信度 somatic TP → 30,511 進 dual-axis ASM（51,171 records，篩除 8,936=23% 無可測軸）→ Wilcoxon p<0.05 11,830 → Bonferroni 313 → strong-ASM (|Δβ|≥0.1) **172 (0.34%)**。**兩方法兩軸全 NEGATIVE**：B1 平均（門檻 artifact，OR 1.79→0.81 反轉）；B2 read-level 聚類（TP median ARI 0.135 << validated ruler imprinted 0.758，GNAS/RB1=1.0；TP vs het MW p=0.922）。深查 H2/H3/H5：穩定但非方向（sign-test p=0.943）、不能預測 TP（AUC 0.505 落 null）、LOH-ASM 非 CNV-dose（ρ=−0.087）。TO 對比 FP/TP=3.11×（強化 TO germline-FP NO-GO）。per-locus 歸因 Q1/Q2/Q3：聚類=germline allelic Layer A；唯一 open lead = chr16:34707014（allele-only，須 BAM）。**單樣本**。**[5/31 紅隊 critique]** 6-agent 對抗式紅隊 + 主 agent 親自重算 5 個硬發現：核心 NEGATIVE 更穩（FN 比 TP 更聚類、null 0.55%>TP 0.34%、somatic<<imprinted 尺），但 5 個正面子 claim 要改口徑（172→FDR 3,018 17.5×、anti-discriminative 5× 中 61% 來自 chr8 hotspot、HP-axis 43% 在 LOH、ahl VAF confound ρ=0.35 未排除、5 候選僅 chr16 robust）；3 條真 open（chr16 BAM / H5 深覆蓋 AUC 0.605 / 5mC-5hmC 未分離）。[觀察文件 md](in_progress/2026/05/20260531_ISM_aux_tag_observation_funnel_01.md) [觀察 HTML](in_progress/2026/05/20260531_ISM_aux_tag_observation_funnel_01.standalone.html) [紅隊 critique HTML](in_progress/2026/05/20260531_ISM_aux_tag_critique_and_improvement_01.standalone.html) [圖片互動審查 HTML](in_progress/2026/06/20260601_ISM_aux_tag_figure_review_interactive_01.standalone.html) [attribution](../../research/tsg_promoter_asm_reviewer/output/clustering_source_attribution.md) [ledger 56-65]

**[6/1 credibility patch + 圖片互動審查]** ①+② 數據化（`scripts/32_credibility_patch.py`）：P1 172(Bonf)→3018(FDR) 17.5× 但 null≥TP 兩校正都成立（correction-robust）；P2 OR 5.16→2.84 drop-chr8（chr8=61% FP-strong）；P3 36%(18,447) records n_cpg≤20 power-censored；P4 bootstrap(TP−null)=−0.11% CI[−0.50,+0.21] 跨 0（somatic 與 baseline 無法區分，L1）。誠實修正「null>TP」→「無法區分」+ 新 surface 高覆蓋 crossover residual（cov≥60 TP>null，與 H5 AUC 0.605 同向，underpowered）。圖片互動審查 master（11 圖 + 8 判斷 localStorage 三態審查 + L3 校對記錄）。**判斷：ASM = G6 phasing 論文有效且足夠的 supplement null，STOP；5mC-5hmC/chr16/COLO829 off-goal。** [ledger 20260531_asm_credibility_patch_correction_robust]

**[6/1 BRCA2 cis 深掘 + 方法驗證]** 釐清 **MSA(MethylSomaticAnalysis)=外部抽取工具≠ISM binary**（Level-3 統計不可用改 Level-1+Python）。對 BRCA2 + 9 對照位點重跑小範圍 MSA（tumor+normal, 8s）做 3 分析：(1) per-tag cohesion BRCA2 HP1-1 somatic 最內聚(silhouette 0.267)；(2) **normal-anchored cis-test**（normal HP1/tumor HP1/tumor HP1-1）—phasing 一致性已驗證—**BRCA2 兩軸(HP+phasing-free ALLELE)都通過 cis-cand**(d_cis=−0.142/−0.152 p<.0001, d_drift≈−0.02)，**germline-het 即使 cohesion 0.598 也只 drift，高 ARI 0.926 TP 也 drift → cohesion/ARI≠cis**；(3) 5mC/5hmC from BAM：BRCA2 5hmC 可忽略(~0.03)、ASM 純 5mC、max-collapse 有效。BRCA2=nonLOH+CN-gain。**結論：cis-test 是判別真 cis-driven 的方法；BRCA2 是稀少真案例，與全域 NEGATIVE 並存**。[BRCA2 cis deepdive HTML](in_progress/2026/06/20260601_ISM_BRCA2_cis_deepdive_methods_01.standalone.html) [ledger 20260601_brca2_cis_test_normal_anchored]

### ClairS-TO HKU Luo Lab handoff — bimodal HP3 finding ⭐4 PI-scope（單樣本 HCC1395）2026-05-24
對 HKU LongPhase-S reviewer response (3 reviewer 23 條) 對齊：5 強對應 + 6 部分 + 12 skip。**核心 finding — HP3 是 bimodal**（全 24 chr scan）：Type A 20/23 well-phased chr HP3 frac 0.22-0.38% / TP rate 78.6-95.4% (mean 90.4%, ~6,900× random 富集) = somatic-evidence bucket；Type B 3 outlier chr6(HLA 9.41%/1.37%)/chr16(segdup 8.58%/3.53%)/chrX(21.6%/0% truth=0) = phasing-failure fallback。Pooled 1.76%/14.3% (outliers 主導 85.4%)。**單樣本 [F] L1**；cross-sample + allele-level concordance 仍 pending。T1/T2/T10/T11 + A2/A6/A7/A8 8 模組全 provenance。[PI 報告](../reports/pi_reports/2026/05/20260524_ClairS_TO_HKU_methyl_subclone_handoff_01.md) [coverage matrix](../concepts/2026/05/20260524_reviewer_coverage_matrix_01.md) [T10 24-chr findings](../../research/hku_collaboration/findings_5_24/T10_HP3_24chr_findings.md)

### PI V6 4-Goal sign-off — ⭐⭐⭐⭐ L2 meta-rollup 2026-05-21
V6 production tag 對 ISM 五目標的 PI 簽核 meta-rollup（非單一 cycle 升級）：Goal 2 marker filter + Goal 5 hard-threshold F1 = V6 嚴格 dominate baseline（GO with calibrated caveats）；Goal 5 LR-F1 BAM-independent（不可 over-claim V6 改進）；Goal 3+4 ⭐0 NEGATIVE。[PI signoff draft](../reports/validated/2026/05/20260521_PI_V6_signoff_email_draft_5goal.md) [ledger 20260522_PI_V6_4goal_signoff]

### V6 vs baseline HCC1395 ISM TP/FP head-to-head — 🟢 V6 strictly dominates baseline 2026-05-19 (Day 2-4 complete: F1 analysis + counter-example + 5-Goal validation + V6 production tag finalized)
**Day 2 V6 production tag supporting experiment**。對 HCC1395 全基因組 30,490 TP + 4,842 FP 跑完整 **4-way ISM head-to-head**（baseline / V3F / V5 / V6 × 2 flag × 2 label = 16 ISM runs，新加 baseline × 4 共 14m38s wall clock）。**V6 vs baseline 主結果**：marker count **+52.4%** (15,738 → 23,980 NG≥3 regions)、marker rate **+1.26pp** (0.8967 → 0.9093 TP purity)、hp=3 ambiguous reads **×13.2** (10,440 → 138,317；TO 模式下 hp=33 bucket 全為 0)、HP1:HP2 all-reads ratio 輕微改善 (1.696 → 1.609)。**V5 over-promote 證實**：V5 marker rate 0.8937 < baseline 0.8967 (-0.3 pp，V5 Layer 1.5 增加 marker 候選但 TP 純度沒保持)。**V3F vs V6 trade-off**：V3F marker rate 0.9175 最高 vs V6 0.9093；V6 marker count +9% over V3F — V6 = V3F-style 保守 + V5 phasing 紅利 hybrid 直接驗證。**LOH zone**：V6 LOH-positive marker rate 0.9801 全 4 BAM 最高 (vs baseline 0.9710)。**Audit 修正 (5/19)**：(1) NG_on=2 公式從 retention rate 修為 V6 doc §5.3 canonical TP purity 公式 — 修正後 V3F=0.8579 / V5=0.8285 / V6=0.8285 與 V6 doc 逐位元相同 ✓；新測 baseline=0.8570；(2) `hp33_reads` 欄改名 `hp_ambig_reads` + 拆解 hp3/hp33（HCC1395 TO 100% 在 hp=3）。**Counter-example hunt (5/19)**：5 維度反例搜尋 — Q1 marker membership: +8,043 V6-only TP / -350 baseline-only TP / +560 V6-only FP → **net +7,693 TP : +549 FP = 14:1 ratio**；Q2 8/22 chrs baseline > V6（最大 chr12 -1.58pp 為 cnLOH chr 殘餘 priority bug）；Q3 AF<0.1 bin baseline > V6 -3.04pp（marker filter 不適用 sub-clonal）；Q4 LOH=NA / Cov=NA zone -1.37pp（框架外灰區）；Q5 350 lost TP **全在 LOH=NA**、**97% V6_NG_off=2 boundary loss**、V3F 仍保留 93% — V6 hp=3 保守化的設計取捨。**無 critical 反例**，主結論不被推翻。**F1 analysis (5/19 加, §11 + F6-F9)**: caller-alone F1=0.7165 sanity reproduce ✓；hard NG threshold T=1 V5/V6=0.7169 僅 marginal +0.0004 略勝 caller (Cohen ribbon 內), T≥2 全 BAM &lt; caller — **hard threshold 非可行路徑**；ΔF1 V6-baseline 跨 T=1-5 全部正 (+0.001 至 +0.197) — **V6 strictly dominates baseline at F1 level**；LOH-positive zone V6 T=3 F1=0.3159 &gt; caller-zone 0.2798 (+0.036)；要打敗 canonical caller F1 需 multi-axis LR (Phase 2 Cycle 1 已證 +0.02236, V6 only, single-sample, cross-sample 受 caller-F1-headroom 限制)。**LR-applied-to-baseline (5/19 §12 + F10)**: 把 Cycle 1 10-feature LR 套用 baseline_off_NG 替換 V6_off_NG (保留其他 9 features) — **驚人發現：baseline LR ΔF1 = +0.02302 反而略勝 V6 LR +0.02236 (+0.00066)，multi-seed 兩者 std=0.00005 同穩定**。LR 自動 shift τ 0.39→0.48 補償 baseline NG 分布差。意味 Cycle 1 +0.023 ΔF1 主要 leverage BAM-independent features (caller_af/LOH/Cov/methylation)，**NG 只是 10-feature 之一，機制 NOT V6-specific**。對 V6 production tag：升級價值在 filter-level + hard-threshold F1 (V6 嚴格 dominate baseline)；不是 LR-based F1 ceiling (LR 不挑 BAM)。**Day 3 5-Goal validation (5/20 §13 + F11-F13)**: (Goal 5 LR ablation methylation) V6 no-meth ΔF1=+0.02170 / baseline no-meth +0.02253 — **methylation contribution 僅 +0.0005-0.0007 遠低 Cohen ribbon +0.005**，H_ABL_1 confirmed ISM 對 F1 沒實質貢獻；(Goal 1 per-CpG × HP × ALT) Cramer's V: baseline 0.1068 ≈ V5 0.1069 (V5 Layer 1.5 沒實質修補 HP × ALT 關聯), V3F 0.0675 最低, V6 0.0899 中等 — V6 比 V3F worse 27% imbalance; (Goal 2 HPFineN cross-sample V6 × 5 samples) avg Spearman ρ=0.845 PASS, HCC1937 outlier (BRCA1 sample); (Goal 3 二次打擊) ⭐0 NEGATIVE blocked by Goal 1 redesign; (Goal 4 TO normal) ⭐0 NEGATIVE caller_af 主導。**V6 production tag final verdict**: 🟢 GO with calibrated caveats — 升級價值在 Goal 2 marker filter + Goal 5 hard threshold F1；**不應 over-claim** V6 改進 Goal 5 LR-F1 (BAM-independent) 或 Goal 3+4。**Day 4 (5/20) finalized**: V6 patch commit `8a90532` + git tag `v6-prod-20260520` (local on longphase-to-mod `fix/pon-only-phasing` branch, +31/-19 lines); manifest.yaml v6_binary_commit + v6_git_tag 填入；PI sign-off email draft `InterSubMod/docs/reports/validated/2026/05/20260521_PI_V6_signoff_email_draft_5goal.md` ready Day 5 review；remote push 待用戶授權。Caveat：ALT-only ratio 未測（PI 4-29 17.3:1 口徑，未來加 alt_support filter 2hr 重跑）。**Verdict**: 🟢 V6 production tag 證據加強 — 4-way 直接 quantify「V6 BAM 經 ISM 後嚴格優於 baseline」，反例 transparency + F1 analysis 加強可信度。[主報告](in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md) [step1 4-way master](../../research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_master_four_way.tsv) [summary TSV](../../research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_baseline_vs_v6_summary.tsv) [F1 圖](../../research/v6_bam_tpfp_hp_loh_cn/figures/baseline_vs_v6/F1_global_summary_4panel.png)

### Phase 2 LOSO Sample-Level CV — Cycle 1 +0.02236 ≡ Sample-Level Circularity Bias — Filter Direction FAILED 2026-05-20 ⭐
**用戶 2026-05-20 直接質疑**「LR filter 用 HCC1395 數據訓練,又用 HCC1395 數據驗證,這敘述合理嗎」觸發 LOSO (Leave-One-Sample-Out) cross-validation。**結果**: HCC1395 in-distribution 5-fold OOF +0.02236 vs LOSO held-out (4-sample combined train) **-0.00012** → drop +0.02248 = **100% effect size 來自 sample-level circularity bias**。5 樣本 LOSO ΔF1 全部 ≈ 0 (HCC1395 -0.00012 / HCC1937 +0.00000 / HCC1954 -0.00008 / H1437 -0.00001 / H2009 -0.00001), Wilcoxon p=0.125, **0/5 positive, 4/5 negative, 1/5 zero → DIRECTION_NEGATIVE**。**Best τ 全部退化為 0.10** (=keep everything) → 4-sample-trained LR 對任何 held-out sample 找不到 useful threshold。**核心 reframe**: Cycle 1 +0.02236 是 HCC1395 self-overfit 非真正 filter effect；ISM vestigial 結論 scope 變窄 (LR 本身在 sample-level 無 signal, ISM contribution 變 moot)；Cross-sample transfer NEGATIVE 從 ⭐⭐⭐ 升 ⭐⭐⭐⭐ (LOSO 強化)。**Tier 重評**: Cycle 1 +0.02236 從 ⭐⭐⭐⭐ L2 降為 ⭐⭐ L4；LOSO sample-level negative 為 new ⭐⭐⭐⭐ L2 結論；Caller-F1-headroom mechanism ⭐⭐⭐ L3 → ⭐⭐ L4 (filter 整體失效讓 mechanism moot)。**Production filter direction FAILED**; Paper §3 完全撤回 "ISM-augmented filter" 宣稱，改 "ISM characterization study + LR sample-level negative finding"; Cycle 4 Trial A/B/C prior PASS 顯著下降 (75%→15% / 30%→10% / 35%→10%); 建議 pivot phase_block_3d 或 thread_d 而非繼續 cycle 4 LR variants。**ISM characterization (v0.3 cycle ⭐3 不變)** — filter direction 失敗 不影響 mechanistic understanding 結論。[LOSO findings](../../research/methyl_augmented_filter_phase2/cycle4/loso_validation/loso_findings.md) [Script](../../research/methyl_augmented_filter_phase2/cycle4/loso_validation/scripts/run_loso_cv.py) [Data TSV](../../research/methyl_augmented_filter_phase2/cycle4/loso_validation/data/loso_cv_results.tsv) [Figure](../../research/methyl_augmented_filter_phase2/cycle4/loso_validation/figures/loso_5sample_dF1.png) [Updated PI Trust HTML](../../research/methyl_augmented_filter_phase2/phase2_pi_verification/phase2_pi_trust_framework.standalone.html) (Section 1.5 added)

### Phase 2 Cycle 3 Step 1.5 — Ablation: ISM Vestigial + caller_af Dominant — Filter Reframe 2026-05-20
**Pre-reg 4 H ablation** 對 cycle 1+2 既有 5 樣本 master TSV 跑 4 variants (full / no-caller-af / no-methyl / no-both) × 2 modes = 40 combinations (~25 sec, cycle 2 master TSV reuse 無 BAM rerun)。**H_M1a FAIL**: drop ISM 5 methyl features → HCC1395 refit ΔF1 +0.02236 → +0.02171 (drop +0.00065 = 3% of uplift)；**4/5 樣本 no-methyl ΔF1 ≥ full** (HCC1937 +0.00037 / HCC1954 +0.00084 / H1437 / H2009 ≈0)；5-sample mean full +0.00619 vs no-methyl +0.00630 (**no-methyl 略勝**) — ISM 為 vestigial covariate 確認，與 cycle 1 step 5c "methylation 訊號實為 caller_af proxy" 線索一致 + V6 vs baseline 5/20 §13 H_ABL_1 結論同方向。**H_A1 MARGINAL**: caller_af coef shrinkage → HCC1954 transfer -0.37744 → -0.12629 improvement +0.25114 (67% disaster mitigated)；caller_af 為主 confound 但餘 33% 未解釋。**機制解讀**: HPFineF coef +0.75 (cycle 1 LR rank 5) 為 caller_af L2 ridge-split 而非 incremental signal；L2 regularization 把相關特徵 split — drop 5 methyl 後 caller_af + LOH + Coverage + NG 重新分配權重 capture 同訊息。**Cycle 3 Reframe Action**: (1) Filter 命名從 "methyl-augmented multi-axis" 改為 "**caller-F1-headroom-gated 4-feature filter** (caller_af + LOH_inner + Coverage_Multiple + NG)"；(2) Paper §3 撤回 "methylation-augmented filter" 宣稱，ΔF1 +0.02236 改歸 4-feature multi-axis；(3) ISM characterization (v0.3 cycle ⭐3) 保留作 mechanistic understanding 但不在 filter design；(4) Cycle 3 Step 2 panel survey 仍 GO 驗 4-feature filter cross-sample；(5) 大規模 ablation (M2/M3/M0/LOFO) 降級 paper supplementary，僅 A2 caller_af KS test 建議跑 (~30 min) 補強 H_A1 殘餘 33% 機制。**Pre-reg prediction 命中** (plan v2.1 事前估 55% A1 PASS + M1a FAIL → 實際 A1 MARGINAL + M1a FAIL)。[Min ablation findings](../../research/methyl_augmented_filter_phase2/cycle3/ablation/cycle3_step1_5_min_ablation_findings.md) [Script](../../research/methyl_augmented_filter_phase2/cycle3/ablation/scripts/run_ablation_variants.py) [Data TSV](../../research/methyl_augmented_filter_phase2/cycle3/ablation/data/cycle3_step1_5_min_ablation.tsv) [Figure](../../research/methyl_augmented_filter_phase2/cycle3/ablation/figures/cycle3_step1_5_min_ablation.png)

### Phase 2 Cycle 2 — Cross-Sample n=5 + V3F/V5/V6 Cross-Binary — ⭐3 保持 + caveat 2026-05-19
**Main session Coordinator + 4 subagent (B1'/B2'+B3'+B4'/C1)** 對 cycle 1 filter 做 cross-sample n=5 Wilcoxon (4 新樣本 V6 ISM rerun 2 hr) + HCC1395 V3F/V5/V6 cross-binary sanity。**H_C1_5 DIRECTION_NEGATIVE**：transfer τ=0.39 fixed 1+/4-（HCC1395 +0.02232 / HCC1954 -0.37744 災難 / 其他 -0.03~-0.07）Wilcoxon p=0.1875 ΔRecall p=0.0625；re-fit per-sample 3+/0-/2≈0 MIXED 僅 HCC1937 +0.00761 有意義 uplift。**H_C1_6 PASS**：V3F +0.02295 / V5 +0.02222 / V6 +0.02232，max_var transfer 0.00073 / re-fit 0.00055（6.8-9.1× under 0.005 threshold）；V6 re-fit bit-exact reproduce cycle 1 +0.022358 drift 0。**Mechanism**：caller-F1-headroom-bounded — 3/4 新樣本 caller F1 > 0.83 + FP density < 4% 不留 filtering room（best re-fit τ→0.10 keep everything）；HCC1937 (F1=0.37, FP density 16%) only 有意義 uplift；HCC1954 transfer 災難 = caller_af coef +3.44 overfit HCC1395 AF 分布非 methylation 失敗。**Decision** (user 5/19)：(1) Cycle 1 ⭐3 保持 + caveat = "caller-F1-headroom-bounded, not cross-sample transferable as production rule"；(2) PI HTML preview defer 到 cycle 3 結果；(3) 啟動 Cycle 3 Caller-F1-headroom-aware redesign (H_C3_1/2/3)。**Implication**：V6 production deployment zero F1 regression（H_C1_6 強化）；filter 需 per-sample re-fit + caller F1 < 0.80 + FP density > 0.10 gate。[Cycle 2 synthesis](../../research/methyl_augmented_filter_phase2/cycle2/cycle2_findings.md) [B3'+B4' findings](../../research/methyl_augmented_filter_phase2/cycle2/cycle2_step_b3_b4_findings.md) [C1 sanity](../../research/methyl_augmented_filter_phase2/cycle2/cycle2_step_c1_h_c1_6_sanity.md)

### Phase 2 Cycle 1 — Global FP Filter (HCC1395) — ⭐3 strong + caveat (cross-sample NEGATIVE 5/19) 2026-05-18
**Multi-agent fan-out (A1/A2/B1/BAM-search + Coordinator)** 從 v1.0 cell-gated pivot 為 global FP exploration + multi-axis LR (no cell gate)。HCC1395 paired-pileup 35,332 rows × 10-feature LR (L2 C=1.0, drop NumReads VIF=217)。**ΔF1 = +0.02236 @ τ*=0.39 (9.24× v1.0 +0.00242 baseline, 2.24× +0.01 Cohen ribbon)**。**3/4 H PASS**: H_C1_1 top 10 cells 94.22% all FP / H_C1_2 +0.02236 > +0.00242 / H_C1_3 ≥ +0.01 / H_C1_4 high-AF FAIL (-0.00011)。**核心發現**: (1) **Step 5c lost TP 17/21 (81%) rescued by cycle 1** — v1.0 cell-gated 過嚴的證據；(2) **Methylation 是 5th-rank covariate** (HPFineF coef +0.75)，主導為 caller_af (+3.44) > LOH (+1.46) > Cov (+1.27) > NG (+1.07) → filter naming 改 "multi-axis filter incl. methylation"；(3) **Multi-seed std=5e-5** (20× below threshold, intra-sample stable)；(4) **NaN MNAR confirmed** (60.4% NME/Epipoly NaN U-shape across AF bins, impute Strategy B 正確)；(5) **VIF=217 NumReads × Cov collinearity** drop NumReads_master 解決 (max VIF 1.83)。**Track B (cross-sample) DEFERRED to cycle 2**: V3F/V5 4 樣本 BAM 物理不存在 (deep search confirmed only HCC1395)；用戶決定 single-sample 收尾。**Tier ⭐3 strong** (HCC1395-internal validated)；⭐4 升級需 cycle 2 V6 4 樣本 ISM rerun + Wilcoxon n=5。[主報告](in_progress/2026/05/20260518_Phase2_Cycle1_Global_FP_Filter_01.md) [Coordinator synthesis](../../research/methyl_augmented_filter_phase2/cycle1/cycle1_findings.md) [Step 0 audit](../../research/methyl_augmented_filter_phase2/cycle1/cycle1_step0_global_fp_audit.md) [Step 1 filter design](../../research/methyl_augmented_filter_phase2/cycle1/cycle1_step1_filter_design.md) [Track A findings](../../research/methyl_augmented_filter_phase2/cycle1/cycle1_track_a_findings.md) [Deployable filter](../../research/methyl_augmented_filter_phase2/cycle1/cycle1_track_a_filter.json)

### Methylation-Augmented FP Filter Pilot (HCC1395) — ⭐3 PARTIAL POSITIVE marginal 2026-05-18
**Multi-agent fan-out (F2/G/H/I + Coordinator)** 從 v0.3 characterization (⭐3) 跨越到 filter F1 evaluation。重跑 phaseC ISM 12 runs with significance (80 min, 移除 `--no-distance-matrix` flag)、build augmented master TSV (35,332 × 202 cols, 13 methylation features × V3F+V5+V6 × off/on)、跑 138 augmented LR + LRT、FP-rich τ sweep、ΔF1 vs caller F1=0.7166。**ALL 5 H POSITIVE**：H1 LRT 16/30 cells q<0.05 (top `Outer|other|cov_normal` p=1.8e-58) / H2 ΔF1 **+0.00242 marginal (< +0.005 Cohen ribbon)** / H3 FP_remov 98.3% > TP_loss 35.0% (filter_signal +0.633) / H4 13 mechanism candidates + 14 PubMed refs (cis-mQTL/cancer ASM/allele-imbalance/repeat/replication timing) / H5 V5 ≈ V6 Δβ=1.87e-5 (median, V6 重用 V5 phased VCF)。**Step 5c TP rescue NEGATIVE**：95.2% lost TP 是 low-AF subclone (caller_af<0.3) + 52.4% 與 FP IQR 重疊；methylation 訊號實為 caller_af proxy；所有 rescue rule NET NEGATIVE vs Step 3 baseline。**Step 5d robustness GREEN with caveats**：ΔF1 stable across 5 seeds (std 2e-5) / 4 unique LRT cells (2 universal cross BAM/flag) / NaN median n_eff/n_total 0.503 borderline YELLOW (NME/Epipoly drop 85-91% Inner cells)。**ΔF1 marginal 主因**: 35% TP loss 不可接受；single-sample HCC1395 未驗 generalizability。**Next**: phase 2 project methyl_augmented_filter_phase2 + cycle 1 pivot 為 global FP exploration + heterogeneous threshold filter (目標 ΔF1 ≥ +0.01)。[主報告](in_progress/2026/05/20260518_V6_Methyl_Filter_Pilot_01.md) [step5_findings](../../research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5_findings.md) [step5c TP rescue](../../research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5c_tp_rescue_analysis.md) [step5d audit](../../research/v6_bam_tpfp_hp_loh_cn/step5_methyl_filter_pilot/step5d_robustness_audit.md) [phase 2 project](../../research/methyl_augmented_filter_phase2/)

### V3F/V5/V6 ISM 三向 × LOH × HP × CN — HCC1395 pilot + 4 樣本擴展 TP/FP characterization — ⭐3 PARTIAL POSITIVE 2026-05-15
**Multi-agent fan-out (A/B/C/D/E + Coordinator)** 完整後分析（無 BAM/ISM 重跑），HCC1395 paired-pileup 30,490 TP + 4,842 FP 跑 3 軸 50-cell grid + 7 道 confound guard + 4 樣本 (H1437/H2009/HCC1954/HCC1937, COLO829 deferred) Wilcoxon n=5。**Paradigm reframe**: 2/3 預設「FP-rich zone」實際 TP-pure — Z-OCH Outer cross_het (FP rate 0.017 << global 0.137, Fisher p=3.8e-62 for **TP**-enrichment) 與 Z-GL Inner gain+LOH 都是 somatic-evidence signatures 非 FP markers。**H4 POSITIVE**：chr8 hotspot LR deviance decomp caller_af 0.393 > **CN 0.211** > HP 0.063 > LOH 0.038，(LOH+CN)-HP = +0.186 (3.7× 0.05 threshold)，CN+AF 主導 chr8 FP enrichment 2.31×。**Z-AUTO** KDE top 5% FP density 捕獲 25.1% 全 FP (5.02×) 但 ~70% 在 known zones 外 (Jaccard 0.184) → **novel FP mechanism lacks LOH/HP/CN signature**。Z-CHR8 + Z-AUTO 共捕獲 ~37% FP，**剩 ~63% 不被此 framework 解釋**。**V5 over-promote** Inner LOH NG=2 region +60% (V5=8,136 vs V3F=5,064 vs V6=5,353) 但 TP rate 沒升 — Layer 1.5 marker over-promotion 直接證據；V6 修補回 V3F 水準。**Cross-sample**: 唯一 signature candidate `Outer|other|cov_high_gain` 5/5 同方向 (Wilcoxon p=0.0625, Δ=+0.0069 受 caller 飽和 ≥0.998 限制)。**Prior art** TumorLens/ROCIT/SGZ/Wakhan/SAVANA 全無同口徑 (Agent D)。**Limitations**: master-join FP loss 90% / ceiling effect / Z-CHR8 sample-specific / HP1/HP2 asymmetry 為設計選擇非 artifact。[主報告](in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md) [研究目錄](../../research/v6_bam_tpfp_hp_loh_cn/) [Step 1](../../research/v6_bam_tpfp_hp_loh_cn/step1_v3f_v5_v6_three_way/step1_findings.md) [Step 2](../../research/v6_bam_tpfp_hp_loh_cn/step2_three_axis_grid/step2_findings.md) [Step 3](../../research/v6_bam_tpfp_hp_loh_cn/step3_fp_zone_zoom/step3_synthesis.md) [Step 4](../../research/v6_bam_tpfp_hp_loh_cn/step4_cross_sample_extension/step4_findings.md) [Prior art](../../research/v6_bam_tpfp_hp_loh_cn/02_prior_art_notes.md)

### `[RETRACTED 2026-04-26]` LOH × AF × CN Biology-Informed Stratified Filter — ⚠ 跨樣本 filter 用途撤回 / characterization 保留 (原 2026-04-22 POSITIVE)
> **撤回理由**：X6 caller_af 跨樣本驗證（2026-04-24）顯示 S3 TP≥0.85 達成比例 1/6（僅 H2009 飽和 artifact）、Wilcoxon S3 > baseline p=1（系統性低於 baseline）、原 v2 HCC1395 TO S3 95.5% 為 stale-binary CN_tier artifact 無法複現。詳見[撤回宣告](../reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md)。HCC1395 single-sample case study + per-sample characterization 保留。

HCC1395 TO 40,115 sites × 新 KDE CN tier × LOH_Subtype (5) × AF_class (3) × cn_tier_F (5) = 32 有效 cells 切分。**雙極 TP/FP 分佈確認**。S1-S7 scheme 比較：**S3 Diploid Het (LOH=None ∩ Near-half ∩ CN∈T1/T2) TP rate 95.5%**（n=380, FP reduction 99.85%, TP:FP fold 8.69×）⭐ / S1 LOH_Strong+Extreme 90.1% / S2 Subclonal+Intermediate 87.4% / **S5 combo (S1∨S2∨S3 \ S4) TP 91.8% (FP reduction 99.37%)** / **S4 ambiguous (LOH=None ∩ Extreme AF) baseline 71.1%**（含 75% TP + 76% FP 混雜，需二級判別）。NG≥3 加入對 S1/S5 邊際貢獻 <1pp。paired_full 7 樣本 per-sample scheme 穩定性驗證（S3 每樣本 ≥99% TP，S1 HCC1937 n=0）。**結論**：~~biology-informed stratified filter 框架 POSITIVE；S3 最純單模組~~ `[RETRACTED]` cross-sample filter 用途撤回（X6 跨樣本反證）；S4 需二級判別觀察保留。Archive TO 6 樣本重跑（下週）補齊跨樣本 TO 驗證。[報告 v2](in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md)

### LOH-constrained phasing 機制揭露 — ⭐ POSITIVE (TO-layer, H-D1/D2/D3) 2026-04-22
**本週最大發現 · TO-層論文主軸 pivot**。用戶 2026-04-22 晚間提問「NG=2 與甲基有關係嗎」觸發 C++ 原始碼回查。`src/core/LabelTest.cpp:265-302` `hp_to_fine_labels()` 證實 HPFineNGroups = `{HP1, HP1-1, HP2, HP2-1}` 4-bucket occupancy count，**與 methylation 無直接計算關係**；設計者 `Stats.hpp:324` 註 "Same haplotype, germline vs somatic"。obs18 跨 6 TO 樣本驗證：**Inner × NG=2 全部 ≥93% same-haplotype**（HCC1395 93.2% / HCC1395_DORADO 99.0% / HCC1937 98.8% / HCC1954 96.5% / H2009 98.3% / H1437 97.0%，median 97%，6/6 一致）；TP gap (Inner same_HP1 − Outer cross_het) median **+0.37**（range +0.05~+0.52，6/6 正向，0 反向）。**機制**：LOH 區物理限制 → somatic SNV 必產 same-hap 分裂；非 LOH 區 germline het 與真 somatic het 產生相同 cross-het phasing，TO caller 無法區分，是 germline-leak FP 物理根源。**論文主軸**：從 "Haplotype-loss-dependent methylation bimodality" 正式 pivot 為 **"LOH-constrained phasing signatures distinguish somatic from germline-like variants in tumor-only sequencing"**（TO 層）；paired AF×NGroups POSITIVE 保留不撤回 + 加註。**H-D4** (flag=on 應消除 Inner/Outer gap) 待下週 P0 paired 對照完整驗證；Thread C flag=on NG≥3=0 直接驗證 H-D1（NG 非 methylation）。⭐5 TO-only on H-D1/D2/D3 三重證實。[09 章](../../research/tpfp_loh_af_kde_discrimination/09_TO_sample_af_lohside_ng.md) [obs18 data](../../research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv)

### Paired Mode HP1-1/HP2-1 Priority Bug Audit — ✅ NEGATIVE（整體）+ V5 Layer 1.5 設計缺陷揭露（局部）2026-05-09 / V6 eval 2026-05-10
用戶 5/9 兩階段提問：(I)「paired tag 是否也會出現 HP1-1/HP2-1 同位點共現偏向 HP1-1（類似 TO 17.3:1 priority bug）」(II)「同位點 germline 缺席但 somatic 共存的 reads 是否會被 longphase-to 錯標到 HP1 系列」。**核心釐清**：paired 用 `longphase-s somatic_haplotag` 是不同 binary、不同 codebase（vs TO mode `longphase-to`），HP tag 編碼 `HP:Z:` 字串（1/2/1-1/2-1/3）vs TO 的 `HP:i:` 整數。**Step A** chr19 整體 HP1-1:HP2-1 = 12,401:14,504 = 1:1.169（近 1:1，對比 TO baseline 17.3:1）。**Step C** 1Mb window som_ratio mean **0.462** / median 0.494 / stdev 0.332，跨 chr 變化反映**真實 sub-clone signal**（LOH / sub-clone heterogeneity）而非 systematic bias。**Step D（5/9 補充）** germline-absent cross-ref（5,789 chr19 events）：baseline hp=11:hp=21 = 3,312:791 = **4.19:1 偏 HP1**（priority bug 次峰）；V3F 全標 hp=33（保守不選邊）；**V5 與 baseline 4.19:1 偏移完全相同 — Layer 1.5 設計繼承 priority bug 而非修補**；V3F 在該區域比 V5 更穩健。**Step B** (binary patch) skip。後續 F-paired-D1（全基因組擴展）/ D2（axis-aligned）/ D3（V5 改回 V3F 的 ISM 影響量化）。[整體 audit](../../research/paired_priority_bug_audit/00_audit_report.md) [Step D 設計缺陷](../../research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md)

### PI Report 4-29 Errata Companion — ✅ 2026-05-09 / amended 2026-05-10
F1 PI 報告 errata patch（5/8 整合報告 §9.2 4 條訊息變更後續 + 5/10 加 E5）。**主結論不撤回**，5 處表述修訂：(E1) chr19 SP1/2/3 從「主要 hotspot」→「可重現案例」（全基因組占比 2.16%, rank 19）/ (E2) V5 working tree → 已 commit（5 commits, d0bcd8c+938f0df 4-30）/ (E3) priority bug 機制證據升級至 read-level 34,855 全基因組鐵證 / (E4) V5 §6 數值為 Pass 1 only，主要功勞 V3F + Layer 1.5 / **(E5, 5/10 加) V5 Layer 1.5 在 germline-absent 區域與 baseline 4.19:1 偏 HP1 完全相同 — 是 priority bug 的 feature 化非修補；V3F 標 hp=33 反而更穩健**。erratum 獨立 companion + 原 PI 報告頂部 banner（5 條 errata 摘要）。[erratum 檔](../reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md) [原 PI 報告（含 banner）](../reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md)

### Self-Phasing 完整觀察整合報告（observation-first synthesis） — ✅ 2026-05-08 / amended 2026-05-10 (§8.6 + F6+F7) / V6 attribution errata 2026-05-13
單一 PI / lab meeting 入口（10 段 + 3 附錄 + 7 figures，~1,211 行）。**整合 7 處觀察**：PI 報告 17.3:1 + chr19 SP1/2/3、V5 audit Pass 1 only 揭露、T1.2 chr19 752 victims、T1.2-F1 全基因組 34,855 victims、**5/9 paired Step A+C 整體 NEGATIVE**、**5/9 Step D V5 Layer 1.5 設計缺陷**、5 commits 修補演進。**機制因果確立**（17.3:1 → 100% V3F/V5 修正、6 cell × 4 metric caller F1 三版完全相同 0.7166 / 0.6273）。**§8.4 V5 vs V3F zero-sum 已釐清** = `d0bcd8c` Pass 2 reclassify 104,457 germline het。**§8.5 20+9 指標 no regression** + caller F1 三版相同 + LOH bed 三版相同。**§8.5.3 三路徑算法不依賴 purity**。**§8.6 (5/10 整合)** Paired Mode Cross-Reference Audit 8 子段：paired 整體無 priority bug + germline-absent 區域 V5 與 baseline 4.19:1 偏移完全相同（V5 Layer 1.5 設計缺陷）。**HEAD `938f0df` = 最新有效演算法**。對 PI 報告 4-29 5 處 errata 已 patch（E1-E4 5/9 + E5 5/10）。**5/13 V6 attribution errata** — 3 處模糊補強（17.3:1 修對主力 V3F / +13.3 pp = V5 V6 未實測 / 34,855 = V3F 主力 V5/V6 logic 繼承），主結論不撤回。9 commits 鏈完整封存。[報告（含 5/13 banner）](../reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md) [PI 4-29 errata](../reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md) [V6 attribution errata](../reports/validated/2026/05/20260513_V6_Attribution_Errata_01.md)

### V6 Attribution Errata Companion — ✅ 2026-05-13
5/8 整合報告補強 erratum：3 處 V6 attribution 模糊修正（不撤回主結論）。源於 PPT v1.6/v1.7 誠信檢查 (b/c/d 三輪 audit)。**E1** 17.3:1 → ~1:1 修對主力為 **V3F** (commit 41ff147)；V5 加 Layer 1.5 退步至 2.00:1；V6 改善至 1.84:1 vs baseline 但比 V3F 退步（trade-off marker eng）。**E2** +13.3 pp paired GT @ 0.93 = **V5 (V3F + Layer 1.5) vs baseline** 達成（PI 4-29 line 824）；V6 重用 V5 phased VCF 預期保留但**未重跑 15-site Clean PS metric**（V6 audit 0 hits）。**E3** 34,855 victims 100% 修對 = **V3F (41ff147 tagging fix) 主力**；V5/V6 此 germline-existent 子集 logic 繼承 V3F（V6 唯一改動 Layer 1.5 revert 對此子集不適用）；V6 audit 沒重跑但 logic V6=V3F=V5 valid。對應 PPT 改動已 apply 至 preview/slide_08/09d/12/17 + REHEARSAL_CHEATSHEET。[erratum 檔](../reports/validated/2026/05/20260513_V6_Attribution_Errata_01.md) [PPT preview](../presentations/validated/2026/05/self_phasing_synthesis_PI/preview/) [REHEARSAL_CHEATSHEET](../presentations/validated/2026/05/self_phasing_synthesis_PI/REHEARSAL_CHEATSHEET.md)

### longphase-to-mod 4-Commit 演進與 V5 Layer 1.5：Self-Phasing 因果鏈整理 — ✅ 2026-04-29 (post fact-check)
`/structured-tech-report` skill 首個示範案例 + **反例教材**（初版 V5 歸屬／函數名／數字錯誤經用戶質疑後全面 fact-check 重寫）。整合 v5_audit_suite 6 子報告 + 4 MEMORY 條目 + PI 報告 4 全基因組數字。**核心釐清**：(a) 修補在 `/big7_disk/liaoyoyo2001/longphase-to-mod/` 獨立 fork 內，**InterSubMod 本 repo 無 C++ 改動**（ISM 為下游消費者）；(b) self-phasing 由 4 commit 漸進修復：`8b8c1fd` (V2b PON-only) → `41ff147` (V3F getVote 兩層+enum→int) → `380e8d2` (INDEL guard) → V5 working tree (Layer 1.5+SNP alt guard, **uncommitted**)；(c) 介面契約零變動 (HaplotagProcess.h:66-68)、+68/-36 行集中 3 函式。**HCC1395 5kHz 結果**：4 項 sanity check 15/15 PASS、SEQC2 F1 0.7153→0.7154、AMB% 17.5→8.0%、HP:i:33 −54%、clean PS（11/15 sites）+13.3 pp、aggregate +6.65 pp。**全基因組（PI 報告 4）**：clean PS V5=90.5%/BL=82.2% (+8.3 pp)。**Open issues**：V5 working tree uncommitted、Confidence threshold 0.6 未直接驗證、7 樣本擴展未做、cnLOH 雙親同源區。skill 因此事件補強 checklist B9（fact-check）+ B10（歸屬正確性）。[報告](../reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md) [skill](../../.claude/skills/structured-tech-report/SKILL.md) [反例教材](../../.claude/skills/structured-tech-report/examples/longphase_TO_vs_V5_example.md)

### 週報 0423 · NG=2 LOH-constrained phasing 發現與 TO-層論文主軸 pivot — ✅ 2026-04-23 ｜ Thread B 段落 `[RETRACTED 2026-04-26]`（[撤回宣告](../reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md)）
涵蓋 2026-04-16 ~ 04-23（8 天）。4 Thread 整合敘事：A (CN KDE 方法學) → B (LOH×AF×CN filter framework) → C (`--germline-hp-only` flag Phase 1) → D (NG=2 LOH-constrained phasing 機制揭露 · 論文 pivot)。Wave 1 驗證 CONDITIONAL PASS，H-D4 宣告已修正為「⏸ 待 paired P0」、H-D1/D2/D3 ⭐5 TO-only 保留；HCC1954 Inner TP 0.43 加限定（物理必然指 phasing pattern 非 TP rate）。結論表 6 項變動：CL-LCP-001 ⭐5 TO-only / CL-S3-001 ⭐4 / CL-CN-KDE-001 ⭐4（新增）· CL-016 ⭐4→⭐3 降級 · CL-LAF-001 保留加註 · CL-HP-ONLY-001 ⭐3 懸掛。[週報](../reports/validated/2026/04/20260423_研究週報_20260416_20260423_NG2_LOH_constrained_phasing與TO_pivot_01.md) [PPT 資料夾](../presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/)

---

## 附錄：待驗證方向（尚未正式啟動）

| 方向 | 期望收益 | 難度 | 依據 |
|---|---|---|---|
| 5hmC 雙通道距離矩陣 | 可能更有亞克隆特異性 | 高 | ONT 5kHz 同時提供 C+h |
| ~~跨 Region 亞克隆一致性~~ | ~~已驗證 NEGATIVE~~ | — | O13/O13v2: shared read count confound |
| ~~機器學習組合特徵分類器~~ | ~~整合 15 個特徵的 ensemble~~ | — | Wave 3 J13: AUC=0.577<0.58; Beyond-AUC 2026-04-09 確認特徵空間耗盡 |
| ~~CNV zone-aware filter~~ | ~~利用 SEQC2 CNV 分層過濾~~ | — | 2026-04-10: 跨樣本不一致；zone-specific mean AUC≤0.641；trade-off 全不可行 |
| ~~GC 校正與甲基化-CN 驗證~~ | ~~GC bias 校正 + CN 驗證~~ | — | 2026-04-11: delta-r=-0.0002（門檻≥0.03）；甲基化 CN-blind；不需實作 |
| ~~PMD/ChromHMM Gating 啟用~~ | ~~降低甲基化背景噪聲~~ | — | 2026-04-12: PMD CpG variance TP=FP; AUC≤0.622; fCpG p=0.77 |
| paired-only multi-bio 全量擴充 | 確認 +0.0112 F1 穩定性 | 中 | round 2 sample637 驗證完成 |
