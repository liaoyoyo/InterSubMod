<!--
建立時間: 2026-03-05 10:00
更新時間: 2026-04-07 01:00
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
完整五步因果鏈驗證：62% LOH 消失（d=-1.20）、somatic bias 17.3:1→消除、31.2% self-phasing LOH。23 TSV + 7 PNG。[詳情](../reports/validated/2026/04/)

### PON-Only Phasing 驗證 — ✅
LOH.bed Jaccard=1.0（不變）；somatic bias 消除；N50 +99.7%；phased rate +23.6pp。待 haplotag+ISM 全量重跑。

### TO Feature Deep Study Q1-Q6 — ✅
HP2FamilyN AUC 0.72 是循環 artifact（校正後 0.54）；RF ceiling 0.69-0.77；baseline>v2b in LOH；甲基化特徵全無效。

### QS Mode-Aware 修正 — ✅
TO 模式下停用 LOH penalty 與 verify bonus；C++ 已實作並部署。

### LOH 雙定義交叉分析 Wave 1+2+3 — ✅ 確定性結論
W1+2: HP Imbalance 是 LOH.bed 超集 (Sens 99.7%); LOH.bed vs SEQC2 Jaccard=0.928; HP PERMANOVA LOH 內不可用 (valid 5%); **LOH 不可作為 variant filter** (10/10 FAIL); AlleleDelta AUC=0.556 真實但微弱。W3: **Non-LOH 區分力同樣有限** (max AUC=0.643 read count proxy); **多特徵組合不可行** (Voting AUC=0.577<0.58); **cnLOH PairwiseMeanDist 0.587 是 Simpson's Paradox** (per-sample mean=0.50); **CramersV 被 NumReads confound** (0.511→0.464); AlleleDelta 是 LOH 唯一 confound-free 信號但不足。40 位點肉眼檢視文件+120 圖表已生成。[完整報告](../reports/validated/2026/04/20260406_LOH雙定義交叉分析報告/00_INDEX.md)

### R1-R5. 特徵設計局限與改進方向研究 — ✅ 2026-04-07
CramersV 93% 為零 = 2×2 缺陷（R1）；Excess groups 無新信號（R2）；結構清楚子集 AUC 下降 = identifiability 非設計問題（R3）；N≥4+NR≥80 TP 89.1%（R4）；PwDist 正交但太弱（R5）。**HPFineNGroups 確認為 somatic heterogeneity 標記**。正負效果+驗證標記已寫入報告。

### 肉眼檢視推理鏈與 TP/FP 可區分性 — ✅ 2026-04-06~07
40 站點（4 類 × TP/FP × 5）系統性肉眼比對 + 全量 7 樣本定量分析。LOH 結構反向確認；Non-LOH ALT-HP alignment 群體不可分；F1 影響分析。完整報告含 R1-R5。[完整報告](../reports/validated/2026/04/20260406_肉眼檢視推理鏈與TP_FP可區分性分析_01.md)

---

## 附錄：待驗證方向（尚未正式啟動）

| 方向 | 期望收益 | 難度 | 依據 |
|---|---|---|---|
| 5hmC 雙通道距離矩陣 | 可能更有亞克隆特異性 | 高 | ONT 5kHz 同時提供 C+h |
| ~~跨 Region 亞克隆一致性~~ | ~~已驗證 NEGATIVE~~ | — | O13/O13v2: shared read count confound |
| ~~機器學習組合特徵分類器~~ | ~~整合 15 個特徵的 ensemble~~ | — | Wave 3 J13: 多特徵組合 AUC=0.577<0.58, LOH/non-LOH 均無突破 |
| PMD/ChromHMM Gating 啟用 | 降低甲基化背景噪聲 | 中 | 架構已有 `is_pmd` 欄位 |
| paired-only multi-bio 全量擴充 | 確認 +0.0112 F1 穩定性 | 中 | round 2 sample637 驗證完成 |
