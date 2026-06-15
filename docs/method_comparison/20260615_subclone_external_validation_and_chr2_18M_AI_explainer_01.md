<!--
建立: 2026-06-15
報告類型: 外部驗證綜整 + 給「說明 AI」的 chr2:18M 位點解釋指南（synthesis + AI explainer）
任務類型: D handoff（讓其他 AI session 正確解釋此位點，不過度宣稱）
狀態: pointer/synthesis — 數字皆本 session grep/script 親驗或一手論文 WebFetch；外部庫在 repo 外
data_sources:
  - InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_seqc2_concordance_demo_01.md（demo 真值）
  - InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/chr2_18M_seqc2_concordance.tsv（script 輸出）
  - InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_independent_verdict_02.md（內部獨立複核）
  - /big8_disk/data/HCC1395/SEQC2/*（LOH BED + sSNV TVAF truth，磁碟親驗）
  - external_validation/（repo 外；Fang2021/PhyloWGS/DREAM 等卡）+ memory project_external_validation_library
provenance_note: 本檔不產新數字；每個數字標來源。CCF=TVAF×CN/(purity×mult) 近似值。
-->

# HCC1395 chr2:18M 亞克隆 — 外部驗證綜整 + 給「說明 AI」的解釋指南

> **這份給誰**：要對外/對內**解釋這個旗艦位點**的 AI session 或人。先讀 §5（可/不可宣稱）+ §7（紅旗），再按 §6 分眾範本說明。**目的＝讓解釋既有力又不過度宣稱。**
> **一句話**：在 SEQC2 已知含多 subclone 的 HCC1395、於外部確認有 somatic LOH+SNV 的 chr2:18M 區，我們用 somatic-haplotag + 甲基提出一個**有界的「區域性亞克隆狀態」假說**；外部資料**佐證此區有亞群 fraction 差異**，但**沒有任何外部研究專門驗證「這裡是一個亞克隆」**，且演化樹的關鍵 pivot SNV 外部無法評估。

---

## 1. 位點身份（事實層）
| 項目 | 值 | 來源 |
|---|---|---|
| 區域 | chr2:18,066,480–18,110,828（≈18.07–18.11 Mb, GRCh38）| 內部分析區 |
| 樣本 | HCC1395（乳腺 TNBC, ductal carcinoma）/ 配對 normal HCC1395BL | SEQC2 標準樣本 |
| cytoband / 基因 | 2p24.2；**KCNS3 intron**（電壓門控鉀通道亞基）| Ensembl |
| 基因癌症角色 | **非 cancer driver**（不在 COSMIC Cancer Gene Census）→ 此區為**資料驅動**選點，非已知 driver 生物學預測 | COSMIC/HPA |
| 樹 pivot SNV | chr2:18,086,020 G>A | 內部 |

## 2. 證據三層（內部 + 外部 + demo 真值）

**內部（2026-06-15 獨立複核 verdict_02）**：可防守 = **≥3 個 regional operational subclonal states**（alpha `(3)=A`；alpha-1 `+(5)=C`；beta-like `(1)(2)`+`(4)(6)`），`(3)→(5)` 有 nested 支持；LOH confirmed；甲基 coherence 跨 basecaller（HKU+DORADO，**12 個 allele-CpG association 中 10 個兩 basecaller 皆過 BH-FDR、另 2 個〔CpG 3.4/3.5〕方向複製但 DORADO 未達 FDR**）。**不可支持**＝5 個 biological subclone / 完整 clone identity/CCF/phylogeny / 完整演化順序 / 甲基造成突變。

**外部（本 session）**：
- SEQC2 Fang 2021（PMC8532138, **Fig 4**）把 HCC1395 **sample-level** characterize 為**多 subclone**（Fig4=S1–S10；**count=10 為 WebFetch 二手、一手待覆核**；PhyloWGS+superFreq+subHMM+10x 單細胞 CNV）——但**全基因組層，無 chr2:18M 標註**；其 subclone driver 在 **17q/6p/16q/X，不在 chr2**。
- 此區落 SEQC2 confirmed **LOH**（`chr2:16,146,119–22,100,000=loh`，磁碟親驗）。
- **5/6** 候選 SNV 在 SEQC2 sSNV truth；pivot `18,086,020` 落 **HC 空隙=truth-unevaluable**。

**demo 真值（`seqc2_subclone_concordance.py`，script 直讀，exit 0）**：
| id | 位點 | SEQC2 狀態 | TVAF | CCF≈ | 分類 |
|----|------|-----------|:--:|:--:|------|
| 1 | 18,068,480 C>G | HighConf | 0.403 | 0.814 | major |
| 2 | 18,072,546 G>C | HighConf | 0.389 | 0.786 | major |
| 3 | **18,086,020 G>A** | **out-of-HC 不可評** | . | . | n/a |
| 4 | 18,096,269 C>G | MedConf | 0.242 | 0.489 | major |
| 5 | **18,099,697 G>C** | HighConf† | **0.048** | **0.097** | **subclonal-minor** |
| 6 | 18,108,828 C>G | HighConf | 0.382 | 0.772 | major |

> † 18,099,697 為低-VAF 位點：SEQC2 共識 = HighConf（TVAF=0.048），但 per-aligner（bwa/bowtie/novo）分類為 LikelyFalsePositive，經 300× 深度重分類為 HighConf（`FLAGS=LowConf_to_HighConf_by_300X`）——低 VAF 個別 caller 不確定、深度+共識救回，與「真低比例 subclonal SNV」一致（非削弱證據，但解釋時須一併說明）。

→ **L1 LOH 100% 一致**；**L2 CCF 梯度 0.10→0.49→0.77-0.81**（18,099,697 CCF≈0.10 = subclonal-fraction，與 Fang 報告該區存在 subclonal-fraction 峰一致〔0.15/0.08 為 Fang 二手、待一手覆核〕）→ **外部 TVAF 獨立佐證此區有亞群 fraction 差異**。

## 3. 六樣本 clone/subclone 外部解答對照（一好一中四缺）
| 樣本 | 有 subclone 解答？ | 來源 + 方法/軟體 | 可當正交真值？ |
|---|---|---|---|
| HCC1395 | ✅ 專門研究 | Fang 2021：PhyloWGS+superFreq+subHMM+10x scCNV → 10 subclone | 🟡 PARTIAL（CNV/LOH-level）|
| COLO829 | ✅ 專門研究 | Velazquez-Villarreal 2020：10x scCNV+CellRanger DNA+DAPC → 4 subclone（A–D，chr8 3-vs-4 copy）| 🟡 PARTIAL（CN-level）|
| NCI-H1437 | 🔴 無 | 僅 SKY 核型 + 變異/CNV truth | 🟠 sanity-check |
| NCI-H2009 | 🔴 無（最缺）| 連 SV truth 都無 | ❌ NO |
| HCC1937 | 🔴 無 | Daemen 2015 變異 + SKY；BRCA1/TP53/PTEN | 🟠 sanity-check |
| HCC1954 | 🔴 無 | bulk SV/CNV truth + 單細胞 CNV 驗證對象 | 🟠 boundary |
> ⚠ **provenance（投稿前必補一手）**：§3 外部 subclone count（HCC1395=10 / COLO829=4）+ §4 SRA accession（SRP199641/PRJNA504037）+ Fang VAF peak（0.15/0.08）皆 **WebFetch 二手**，external_validation 庫尚無對應 verified 卡（PMC8532138/PMC7316972 僅出現於本文件）→ 引用前過 `/citation-verification` 補一手。
> 跨樣本最常用 subclonal 分析法：**單細胞 CNV（10x+CellRanger DNA）**、**PhyloWGS**、**superFreq**、subHMM/DNACopy；核型靠 SKY+FISH。標準工具群（PyClone-VI/DPClust/SciClone/Battenberg/FACETS）在 4 缺樣本上**方法存在但無 line-specific 解答**＝若要自證的 Methods 選項。

## 4. SEQC2 方法可信度 + 限制
- **可信度兩層**：「HCC1395 有 subclone」= 高可信（四法交叉 + 真正 orthogonal 的單細胞 CNV）；「精確 10-subclone 樹 + CCF」= 較不確定（單一 PhyloWGS pipeline）。
- **DREAM benchmark（Salcedo/Boutros, Nat Biotechnol）**：即使純遺傳 SR，**演算法選擇就吃掉 19–35% 變異、無單一最佳法** → SEQC2 樹**非金標準**。
- **可下載驗證**：原始資料 NCBI SRA `SRP199641`/BioProject `PRJNA504037`（+ Zhao 2021 Sci Data）；scRNA 矩陣在 Figshare；程式碼在 sites.google.com/view/seqc2。**但 subclone 成品（樹/每群 CCF/每細胞指派）無乾淨下載檔**，須抽 Fig4/supplement 或重跑。

## 5. 🟢 可宣稱 / 🔴 不可宣稱（解釋此位點的核心守則）

**🟢 可宣稱（有證據）**
- 「chr2:18M 區落在 SEQC2 confirmed LOH 內，含 SEQC2 truth somatic SNV（5/6），其中一個 SEQC2 HighConf SNV 的 TVAF≈0.05（CCF≈0.10）= **subclonal-fraction SNV（minor subclonal VAF）**，外部佐證此區存在低比例細胞亞群」（注意：是「低比例亞群的 SNV」，**非**「一個被命名的 subclone」）。
- 「HCC1395 為 SEQC2 已 characterize 的 **multi-subclone** reference cell line（sample-level）」。
- 「本研究在此**雙重外部基底**（sample-level 多 subclone + locus-level somatic LOH+SNV）上，個案示範『同一長讀聯合觀測 somatic allele + haplotag + native 5mC』，提出**有界的 regional subclonal-state 假說**，並有 cross-basecaller（HKU/DORADO）甲基 coherence」。
- 建議定型用語：**"a regional, LOH-constrained, somatic-haplotag-conditioned subclonal structure with cross-basecaller methylation coherence."**

**🔴 不可宣稱（會被打臉）**
- ✗ 「確認 5 個 subclone / 完整演化樹 / 確認演化順序」（內部 verdict 已降為 ≥3 operational states）。
- ✗ 「甲基造成突變」「甲基獨立重建 subclone」（甲基是 corroborate/characterize，T3 為負）。
- ✗ 「SEQC2 驗證了 chr2:18M 是一個 subclone」（外部只到 sample-level + somatic-event level，**無 locus-level subclone 對應**）。
- ✗ 「branch SNV 18,086,020 被外部/SEQC2 confirm」（落 HC 空隙=unevaluable，**樹 pivot 只有內部證據**）。
- ✗ 「DORADO 是獨立病人/biological replicate」（同 cell line 的技術/資料重現）。
- ✗ 把 SEQC2 樹/TVAF 當**金標準逐 subclone 對拍**（PhyloWGS 單 pipeline，DREAM 19–35%；只能 concordance/sanity）。
- ✗ 「1:1 CCF 一致」（passage drift + 平台差 + 軸不同〔他們 SNV-VAF/CNV，我們 haplotag+甲基〕）。
- ✗ 把 LongPhase-S 當「獨立外部驗證」（自家 CCU 工具，只證 method 可行）。

## 6. 分眾解釋範本（說明 AI 可直接改寫）
- **一句話（非專家）**：在乳腺癌細胞株 HCC1395 一個「確認失去一條染色體拷貝（LOH）」的 chr2:18M 區，我們用奈米孔長讀同時看突變、單倍型與甲基，提出一個**有界的「區域性亞克隆狀態」假說**；外部資料確認此區確有不同比例的細胞亞群（含一個約佔 10% 的小亞群），但**沒有外部研究專門驗證「這裡是一個亞克隆」**。
- **給 reviewer**：強調「characterization + tooling，非 reconstruction」；明列三層證據與其邊界（LOH 外部 confirmed、somatic SNV 5/6 confirmed、subclonal fraction 梯度外部佐證；但 pivot SNV unevaluable、樹為 hypothetical、無外部 locus-level subclone、單 cell line 技術重現）。
- **給 PI**：「這是 proof-of-concept 旗艦個案——展示 ISM 能在 LongPhase-S fine tag 拆不開時補上 regional molecular-state characterization，並用 matched normal 揭露既存 ASM confound；外部（SEQC2）只能在 CNV/LOH 層 sanity-check，故結論封頂 ⭐3、需 6 樣本/normal-baseline cross-validation 才升 ⭐4。」

## 7. 說明 AI 紅旗清單（出現以下措辭＝停下修正）
1. 「已重建/確認出 N 個 subclone 的演化樹」→ 改「有界的 ≥3 regional operational subclonal states + 假想樹」。
2. 「SEQC2 驗證了這個 subclone」→ 改「SEQC2 驗證 sample-level 多 subclone + 此區 somatic 事件；locus-level subclone 無外部對應」。
3. 「甲基揭示/驅動了 subclone」→ 改「甲基 corroborate/characterize lineage coherence；不驅動、不獨立重建」。
4. 提到 18,086,020 卻沒標「out-of-HC、外部不可評」→ 必補。
5. 「DORADO 獨立驗證」→ 改「同 cell line 跨 basecaller 技術重現」。
6. 把任何 CCF/VAF 數字寫成「與 SEQC2 一致到小數」→ 改「fraction 帶一致（concordance），非 1:1」。

## 8. 指引（檔案/script/citation）
- **重算/比對**：`InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/seqc2_subclone_concordance.py`（可套全基因組/COLO829）+ demo 報告 `…/20260615_chr2_18M_seqc2_concordance_demo_01.md`。
- **內部複核**：`InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_independent_verdict_02.md`。
- **外部庫**（repo 外）：`/big7_disk/liaoyoyo2001/external_validation/`（Fang2021/PhyloWGS/DREAM/Velazquez-Villarreal 等卡）；橋接索引 `InterSubMod/docs/method_comparison/20260613_external_validation_library_index_01.md`；memory `project_external_validation_library`。
- **citation**：HCC1395 subclone=Fang 2021 Nat Biotechnol(PMC8532138, Fig4)；COLO829=Velazquez-Villarreal 2020 Comm Biol(PMC7316972)；方法可信度=PhyloWGS(Deshwar 2015)+DREAM(Salcedo/Boutros)；資料=SRP199641/PRJNA504037。
- ⚠ 一手待覆核：Fang Fig4 的 S1–S10 完整 CCF（正文只給 S2=60%，須抽 supplement）。
