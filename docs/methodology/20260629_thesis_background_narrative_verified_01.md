<!-- provenance-verified: 數據錨點為 6 路只讀查證 agent 於 2026-06-29 grep 讀回之真值（sm_summary.json / single_snv_accounting.json / sm_locus_master_summary.json / region_shape_distribution.json 三檔交叉一致 + BRCA2 報告）；文獻引用分「庫內親讀已驗」與「WebSearch 候選須過 /citation-verification」兩級標明，未杜撰任何 DOI/作者/年份。 -->
---
title: 論文背景敘述（verified）— Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing
date: 2026-06-29
type: methodology / thesis-background
status: draft（背景段定稿前的 verified 整合；引用待 /citation-verification）
task_type: D-handoff-prep（碩論背景段 + 口試講稿）
tier: ⭐3（單樣本 characterization；對外口徑同 SoT）
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_summary.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/single_snv_accounting.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/sm_locus_master_summary.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/region_shape_distribution.json, docs/experiments/in_progress/2026/05/20260529_ISM_ZAR1L_BRCA2_ASM_verification_01.md
related:
  - docs/methodology/20260627_subclone_unified_verified_narrative_01.md
  - docs/methodology/20260628_subclone_unified_narrative_reverify_addendum_01.md
  - memory: reference_hp_tag_definition_and_subclone_caveat
---

# 論文背景敘述（verified）

> **一句話（BLUF）**：本研究在 **ONT 長讀單樣本 tumor/normal** 上，用 **somatic sSNV 的單分子共現**重建腫瘤的 subclone 結構（**唯一非循環骨幹**），用 **germline haplotagging** 當乾淨的等位鑑別軸，把 **甲基化**定位為被 characterize 的**有界輔助層**。背景敘述的脊椎是「**三條軸不可混**：clone/subclone（譜系）≠ haplotype（等位）≠ methylation（表觀）」。

> **本文件用途**：① 自己鞏固理解；② 內部單一真值（SoT）；③ 教授口頭報告（週報/口試）。
> **狀態**：背景段 verified 整合稿。數據已 grep 真值；文獻引用分兩級（庫內親讀 vs 須過 `/citation-verification`）。**未經 citation-verification 的候選書目不可直接進 .bib**。

---

## §1 背景敘述（概念階梯，含引用標記）

### §1.1 癌症 = DNA（與其修飾）發生改變的結果 — 三種失活路徑

癌症本質是**體細胞基因組（與其表觀修飾）累積改變**而失去正常調控。以一個抑癌基因為例，同一個功能可被**三種獨立路徑**關掉——這正好對應本研究要同時看的三種資訊：

| 失活路徑 | 教學例子 | 對應本研究的軸 | 引用狀態 |
|---|---|---|---|
| **germline 天生突變**（遺傳第一擊）| **BRCA2** 遺傳性致病變異 | haplotype phase | 〔需引用：Wooster 1995 BRCA2；須 /citation-verification〕 |
| **somatic 後天突變**（第二擊）| **BRCA2** 腫瘤內體細胞變異 / LOH | somatic sSNV | 〔需引用：Knudson 1971 two-hit；須 /citation-verification〕 |
| **啟動子甲基化沉默**（表觀第二擊）| **BRCA1**（同 HRR 路徑姊妹基因）| 甲基（methylation profiles）| 已驗+候選（見下）|

> 🔴 **科學精準度（查證 VERIFIED）**：**BRCA2 啟動子甲基化在腫瘤中罕見/有爭議**〔候選：PMID 19340607「BRCA2 promoter methylation has not been reported in sporadic tumors」；須 /citation-verification〕，**故甲基沉默那條腿改用 BRCA1**（與 BRCA2 同屬 HRR/乳卵巢路徑，是最自然的姊妹基因，且貼合本研究 HCC1395 乳癌樣本）。
> - BRCA1 啟動子高甲基化致散發性乳/卵巢癌沉默：〔候選一手源：JNCI 92(7):564, PMID 10749912（經典）；乳癌 ~10–15%（PMC10550062）；卵巢 35.2%（PMC3961372）；**須 /citation-verification 鎖癌種與比例**〕。
> - BRCA1 cis-等位特異性甲基化沉默先例：**Evans et al. 2018, Am J Hum Genet 103:213-220（PMID 30075112）**【KB 已驗】——此源也直接呼應本研究 ASM 框架。
> - 🔴 **避免與本研究結果混淆**：本研究的 BRCA2 是 **等位特異甲基化（ASM, HP-axis Δβ=−0.122, hypo 方向）**，**不是**啟動子甲基化「沉默」；背景段的 BRCA2（germline/somatic 兩擊）與結果段的 BRCA2（ASM）是兩件事，敘述須分開。

**Knudson 兩擊與「四路徑共識」措辭（查證 PARTIAL）**：抑癌基因常需兩擊失活；失活可由 germline 突變 / somatic 突變 / LOH / 啟動子甲基化沉默等**多種機制**達成。⚠ 我們的驗證庫**無**把這四條列為單一「共識路徑集」的綜述，故**措辭採「多種機制可致失活」**而非宣稱「四路徑共識」，各機制各掛來源。

---

### §1.2 從正常細胞到 clone / subclone（譜系軸）

- 腫瘤是**體細胞 clonal evolution** 的產物：源自單一祖先細胞，累積 somatic 突變後分化。〔需引用：Nowell 1976「The clonal evolution of tumor cell populations」, Science；須 /citation-verification〕
- **clone** = 同祖先、共享 genotype 的細胞群；**subclone** = clone 內再獲突變的**子集**，關係有兩型：
  - **nested（祖先–後代）**：後代帶祖先全部突變 + 新增突變。
  - **sibling（並列）**：兩支各自獲得不同突變、互斥。
  〔定義來源：BitPhylogeny 2015（PMID 25786108）等系統發生樹範式「樹隱含 nested + branching」【庫內已驗 tier C】；正式定義句建議補 clonal-evolution 綜述（McGranahan & Swanton 2017 Cell / Greaves & Maley 2012），須 /citation-verification〕

> **🔴 取樣時機點（決定能 claim 什麼）**：ONT bulk tumor/normal = **單一時間點的混合快照**。
> - 看到的是當下所有 subclone 的**混合**；
> - **看不到時間順序**——祖先/後代是從**巢狀子集關係推論**（CCF/pigeonhole 約束），非時間軸觀察〔Tarabichi 2021「illusion of clonality」，PMID 33398189【庫內全文已驗】〕；
> - 因此單 bulk 重建 = **characterization（⭐3）**。

---

### §1.3 三軸不可混 + HP scaffold（等位軸）— **全篇最 load-bearing**

有三種「看起來都是 2+ 群」的分群，意義完全不同，**絕不能混**：

1. **兩個 haplotype**（等位軸）：germline 雜合把 read 分成 HP1/HP2 = 哪一條父母染色體。
2. **兩個甲基 cluster**（表觀軸）：可能是 ASM。
3. **兩個 subclone**（譜系軸）：才是重建目標。

**HP tag 精確定義（源碼 L1 已驗，canonical 見 memory `reference_hp_tag_definition_and_subclone_caveat`）**：

- **`HP:Z:1` / `2`（HP1/HP2）** = 由 **matched normal 的 germline 雜合 SNP** 相位定出（LongPhase-S 命令 `-s <germline_phased_vcf> -b <normal_bam>`，`scripts/pipeline/steps/01_longphase_s.sh:155-168`）→ **完全不碰 tumor somatic、無循環依賴 = 乾淨的等位鑑別軸**。
- **`1-1` / `2-1` / `3`** = 在 germline 骨幹**之上**再用 somatic 變異定的子相位（`1-1`=錨 HP1、`3`=germline 無支撐 HP3）。
- 🔴 **HP1-1 ≠ 已確認 subclone**：它是 somatic-phasing 衍生 tag，把它當 somatic 結構的獨立證據 = **循環依賴**（故 `--germline-hp-only` 把 `1-1/2-1/3` 降為 unphased，只信 `1/2`）。
- 5 個 tag 類別（+unphased）：`1 / 2 / 1-1 / 2-1 / 3`〔LongPhase-S 卡 §4 / README:152-155 親讀；「9 態」為某查證口誤，採 5 類〕。

> **一句話**：HP1/HP2 是乾淨等位軸；**subclone 重建的非循環骨幹 = 同一 germline-HP 上 ≥2 個 somatic sSNV 的單分子共現**（互斥=sibling / 嵌套=nested），不是 HP1-1 tag 本身。

---

### §1.4 為何 subclone 難 + 確認黃金標準

- 傳統 bulk short-read 靠 **VAF 聚類**重建，有 non-identifiability / multiplicity 問題〔Tarabichi 2021；Dentro/PCAWG-11 2021（Cell 184:2239）：即使 2,658 WGS 單切片仍 partial 重建【庫內已驗】〕。
- **確認黃金標準 = single-cell / multi-region**（提供 single-bulk 缺乏的確認性證據）〔Tarabichi 2021 PMID 33398189「tree inference should generally be restricted to multi-sample studies」【全文已驗】；Gaiti 2019 Nature 569:576, PMID 31092926 single-cell lineage【全文已驗】〕。
- 🔴 **「黃金標準」措辭軟化**：無單一外部來源逐字稱 "gold standard"，建議寫「提供 single-bulk 缺乏的**確認性/正交證據**」。
- 🔴 **「Tarabichi LEARN framework」是誤稱**（"LEARN" 是我方 context 卡章節標籤）→ 正名 **Tarabichi 2021, "A practical guide to cancer subclonal reconstruction"**。

---

### §1.5 ONT 平台定位（三優勢 → 為何選 ONT）

| 維度 | 短讀（二代）| ONT 長讀 | 引用 |
|---|---|---|---|
| 讀長 | ~100–150 bp，跨不過長距離雜合位點 | 數十 kb，可把遠處變異連到同一 haplotype | 〔讀長為定序常識；Smallwood 2014 100bp 短讀實例【tier A 已驗】〕 |
| 甲基資訊 | **需另做實驗**（bisulfite / EM-seq）| **native 5mC，basecaller(Dorado) 從電流訊號一次輸出** | Smallwood 2014（bisulfite, PMID 25042786【A】）；Dorado/Modkit【已驗】；EM-seq〔需引用 NEBNext〕 |
| 單分子共現 | 只能看 ≤讀長內共現 | **同一條分子上多個 somatic 變異的共現** = 骨幹使能者 | LongPhase-S（Ho 2025）【已驗 B】；Foltz 2024（PMID 39149342，linked-read genetic-only 先驗）【已驗 B】 |

**ONT 三優勢**：① 長程 phasing/haplotagging；② native 甲基（無需額外實驗）；③ **單分子多變異共現**（重建骨幹的命脈，**勿漏**）。

> 🔴 **定位紅線（查證強制）**：
> 1. 共現措辭要精確為「**短讀無法觀測跨長距離的單分子多變異共現**」（短讀可看讀長內共現；linked-read 即為克服此限而生）。
> 2. **禁**寫「對手用短讀 / 缺顯著性檢定」——cvlr / ASMS / MethylBERT 皆 **ONT-capable 且都有 randomization 檢定**，會被 reviewer 打臉。Claim「短讀取甲基需 bisulfite」是**平台事實**，與「對手平台」是兩回事，措辭勿混。
> 3. 優勢定位 = **實用/低成本 + 一次資料三軸**，**不宣稱解析度勝過 single-cell**（single-cell 仍是確認黃金標準）。

---

## §2 研究範圍 / 目標 / 非目標

| 維度 | 定義 |
|---|---|
| **Scope** | HCC1395 單樣本、**單一時間點**、tumor + matched normal、ONT 長讀；全基因組 |
| **Goal** | 從 somatic 共現骨幹**重建 subclone 結構**（區域 partition + nested/sibling）；haplotag **鑑別** allelic vs clonal；甲基**輔助 characterize** |
| **Non-goal（誠實邊界）** | ❌ 不產完整 CCF 系統發生樹；❌ 不觀察時間順序（靠巢狀推論）；❌ 不達 single-cell 解析；❌ 甲基非 subclone 驅動/確認器（subclone-specificity **UNDETERMINED**）|
| **Tier** | ⭐3（單樣本 characterization；升 ⭐4 需 ≥5/7 樣本，COLO829/matched-normal 為共同 blocker）|

---

## §3 數據錨點（grep-confirmed，可直接寫進論文）

> 全部由 2026-06-29 只讀查證 agent grep 讀回，三檔交叉一致（`sm_summary.json` / `single_snv_accounting.json` / `sm_locus_master_summary.json`）。

| 指標 | 值 | 來源檔 | 注意 |
|---|---|---|---|
| union somatic sSNV | **35,332**（TP 30,490 + FP 4,842，sum✓）| sm_summary.json `universe` | — |
| read-span 區域 | **7,143** | region_shape_distribution.json `total_regions` | — |
| 有確認結構區（對外口徑）| **4,678**（65%，含 858 單-lineage）| region_shape_distribution.json `structured_plus_colinked` | 🔴 非 3,820；引用標明含單 lineage |
| full_tree | **677**（9.5%）| region_shape_distribution.json | — |
| CN-gain somatic-weighted | **52.8%**（12,569）| sm_locus_master_summary.json `cn_somatic_pct.gain` | 🔴 禁與 68.9%（linked 子集）/ 77.8%（segment, UNVERIFIABLE）並列 |
| BRCA2 ASM Δβ（HP-axis）| **−0.122**（n=197 CpG, p=6.62e-10）| 20260529_ISM_ZAR1L_BRCA2 報告:14,80 | 🔴 是 ASM 非 silencing；舊值 −0.05~−0.07 **已撤回勿抄** |
| Tier | **⭐3** 單樣本 | sm_summary.json `scope` + 兩 SoT | — |

> 🔴 **跨 build 禁並列**：本 build FP=4,842 / TP=30,490；舊 build = 34,736 / TP 30,077 / FP 3,643——不可混用。

---

## §4 說明解釋方式（圖例注意 + 敘述強調）

**🖼 圖例最大地雷**：任何把 reads 分左右兩欄/兩色的圖（HP heatmap、甲基 cluster），**圖名與圖例不可讓讀者推斷「左=subclone A、右=subclone B」**——必明標那欄是「等位 HP」還是「甲基狀態」。

**🔊 必強調（不講難理解）**：
1. 三軸不可混（clone ≠ haplotype ≠ methylation）。
2. HP1/HP2 = germline 乾淨等位軸；**HP1-1 ≠ subclone**。
3. 骨幹 = somatic sSNV 單分子共現（互斥=sibling / 嵌套=nested）。
4. 單樣本 = characterization（⭐3），不是 confirmation。
5. 甲基 = 有界輔助，subclone-specificity UNDETERMINED。

---

## §4.5 Confound 教學框：「HP1-1 內甲基次結構」5 假說裁決（核心 confound SoT）

> 此框是「為何甲基只能 bounded-auxiliary」的可對外完整論證。源於 BRCA2 案例 + 群內次結構窮盡測試（數字皆 grep 真值）。

### 前提定義
- **HP1-1 ≈ clone**（帶某 somatic 變異 X 的 read 群之 read-level 代理；遺傳定義、非循環）。
- 🔴 **HP1-1 ≠ 已確認 subclone**（單一 marker、read 尺度；subclone 需樹中位置 + 第 2 個遺傳錨）。源碼 canonical 見 memory `reference_hp_tag_definition_and_subclone_caveat`。

### 現象（grep 真值）
依 carrier(REF/ALT) / HP / HP-fine 軸，ISM **有時**看到明顯甲基分群：clear≥2 群 = 9.6%（2,941 位點），對齊 **germline/carrier 85% / HP1-HP2 15% / 隨機 0%**〔`project_subcluster_cluster_count_determination`〕。→ 非亂切，但「確認 ≥2 群 ≠ 確認 subclone」。

### 5 假說 discrimination table（能否用現有資料排除/確認）
| 機制 | 為何會產生群內甲基多群 | 現有資料能測 | 單樣本能排除/確認 |
|---|---|---|---|
| **1 subclone** | 子群再獲 somatic 變異 + 表觀分化 | ✅ 第 2 sSNV 共現（骨幹 census 7,143 區 / full_tree 677）| 有遺傳錨→確認；無→不能只靠甲基 |
| **2 CNV** | 部分細胞 CN 改變 → read 混合/劑量移位 | ✅ CN 註解（cn_somatic；CN-gain 52.8%）| ✅ 可標 gain/LOH/neutral |
| **3 cis-ASM** | 序列在 cis 驅動甲基 → 依等位分流 | ⚠️ matched-normal cis-control（已做 06-28）| CROSS-HP 35.4% 可分 / **SAME-HP 59% 結構性不可分**（軸正交 corr −0.026）|
| **4 epi-drift / intrinsic epi-heterogeneity** | 基因型相同細胞隨分裂累積隨機 epimutation | ❌ 需 single-cell lineage | 🔴 **單樣本不可消的殘餘 null** |
| **5 artifact** | 低覆蓋 / 飽和 / 鏈偏 / basecall 噪聲 | ✅ 覆蓋/飽和/鏈 QC（飽和 96.5% 已知假象）| ✅ 可篩 |
| **＋ dispersion 假象** | 群差在**變異度**非**均值** | 🔴 PERMANOVA 不分 | 「85%」分離 **~72% 是 dispersion**，扣 HP clean-location 僅 23.8% |

**裁決**：subclone / CNV / artifact **可確認或排除**；cis-ASM **部分**（CROSS-HP 可、SAME-HP 不可）；**epi-drift 單樣本結構性不可消** → **這就是 ⭐3 天花板的精確根源**（非覆蓋/方法不足，是 epi-drift 與 SAME-HP subclone 在單一 bulk 不可分離，需 single-cell）。

### 🔴 不對稱可用性（最值得記）
- ✅ **負向（可信，但要排飽和）**：HP1-1 內甲基 unimodal（實測 95.8%）→ 「無證據顯示有次群」。⚠ unimodal 含「真平坦」與「**飽和無法解析**」兩種（飽和 96.5% 是已知假象）→ 不等於「證明無次群」。補：「切不出 ≠ 沒訊號」，a-priori Δβ 可救回 ~21% 真差。
- ❌ **正向（double-dip）**：HP1-1 內甲基多群 → **不能直接說有 subclone**（5 機制皆可能、需獨立遺傳/CN 錨）。
- 22 個 substructure 候選（0.64%）**無遺傳佐證 = L3 旗標非確認**（多數可能 intrinsic epi-heterogeneity）。

### epi-drift null 的 transferability 註記
epi-drift 作為 null **合理**（與 dispersion 主導一致）；它是 **Gaiti 2019 single-cell（MscRRBS 2,435 cells）** 譜系追蹤的基礎，但那**需 single-cell**，**不直接移植到我們的 bulk read-level**——在我們資料裡 epi-drift 只是**擾亂項，不是可用譜系訊號**。

### 甲基的真實輔助價值與邊界（= bounded-auxiliary）
- ✅ **有價值處**：① 補 **sSNV-poor read** 的覆蓋（每 read 多 CpG 特徵，可 cluster 無 somatic marker 的 read）；② 比 **segment-scale CNV 更細的 read-level 解析度**（CN 段內 flag 次結構）。
- 🔴 **冗餘處**：若是 **cis-ASM**，甲基與遺傳軸共分離 = **冗餘非獨立**（baseline-dependence）；獨立價值集中在「遺傳/CN 不足以定義」的區。
- 🔑 **角色**：在遺傳/CN 框架**之上**的**註解/特徵化層**（標記、描述、L3 旗標、負向篩選），**不獨立確認 subclone**。

### 可辯護定稿口徑（對外/口試）
> 甲基在 clone（HP1-1）內：**unimodal → 保守負篩（排飽和）**；**multi-group → L3 旗標**（subclone/CNV/cis-ASM/epi-drift/假影/dispersion 皆可能）。甲基的真實價值 = **比 CNV 更細的 read-level 註解層**，對 sSNV-poor / CN-only 區最有用，但 **single-sample 下「實際意義」UNDETERMINED、不獨立確認 subclone**。reviewer 問「為何不用甲基定 subclone」→ 引此表。

---

## §5 口試講稿骨架（9 步）

1. **臨床鉤子**（短）：癌=clonal evolution → ITH → 抗藥/免疫逃逸 subclone 驅動復發〔需引用 McGranahan & Swanton 2017 / Gerlinger 2012〕→ 為何要重建 subclone。
2. **三機制**：同一抑癌基因三種失活（BRCA2 germline/somatic + BRCA1 甲基）→ 帶出三軸動機。
3. **概念階梯**：somatic → clone → subclone（nested/sibling）。
4. **三軸不可混 + HP scaffold**：HP1/HP2 乾淨等位軸、HP1-1≠subclone。
5. **方法 gap**：VAF 聚類困難 + 甲基 read-level 未建立。
6. **ONT 三優勢**：phasing + native 甲基 + 單分子共現。
7. **本研究骨幹**：somatic 共現重建 + 甲基有界輔助。
8. **誠實邊界**：⭐3 characterization、單樣本快照、subclone-specificity UNDETERMINED。
9. **收束**：實用/低成本一次資料三軸的定位。

---

## §6 引用台帳

### 6a. 庫內親讀已驗（可引，標 tier）
- **Tarabichi et al. 2021**, A practical guide to cancer subclonal reconstruction, Nat Methods 18(2):144-155, PMID 33398189【B, 全文】
- **Gaiti et al. 2019**, Nature 569:576-580, PMID 31092926【B, 全文】（single-cell lineage）
- **Smallwood et al. 2014**, Nat Methods 11(8):817-820, PMID 25042786【A, 全文】（scBS-seq / bisulfite）
- **LongPhase-S (Ho et al. 2025)**, bioRxiv DOI 10.1101/2025.11.20.689492【B, 源碼+全文】（somatic haplotagging 骨幹）
- **Foltz et al. 2024**, PMID 39149342【B, 全文】（單分子共現→subclone genetic 先驗；linked-read 非 ONT）
- **Evans et al. 2018**, Am J Hum Genet 103:213-220, PMID 30075112【KB 已驗】（BRCA1 cis-ASM silencing）
- **Dentro/PCAWG-11 2021**, Cell 184(8):2239-2254【C】；**BitPhylogeny 2015**, PMID 25786108【C】
- ONT **Dorado / Modkit**（工具，引 repo URL + 版本 + dorado model id）

### 6b. WebSearch 候選，**須過 /citation-verification 才可進 .bib**
- BRCA1 啟動子甲基化：JNCI 92(7):564, PMID 10749912；PMC10550062（10–15% 乳）；PMC3961372（35.2% 卵巢）
- BRCA2 罕見：PMID 19340607
- MLH1：PMC3612054（PLOS One）；PMC2375277, PMID 11433526
- CDKN2A：PMC4919535；PMC3618325（NSCLC）；PMC4383544（HNSCC）

### 6c. 作者須補的基礎/臨床經典（庫內無，須 /citation-verification）
- Knudson 1971（two-hit）；Wooster 1995（BRCA2）；Miki 1994（BRCA1）
- Nowell 1976（clonal evolution）；Greaves & Maley 2012；McGranahan & Swanton 2017（ITH 臨床）；Gerlinger 2012（NEJM）
- Jamal-Hanjani 2017（TRACERx 基因組 subclone 重建，庫內僅 CAMDAC≠tree）
- NEBNext EM-seq（短讀甲基替代技術）

---

## §7 6 路查證裁決表（traceability）

| 群 | 主張 | 裁決 | 關鍵依據 |
|---|---|---|---|
| 1 | 兩擊/四路徑/BRCA 遺傳性 | NO_SOURCE/PARTIAL | 庫為方法導向、無奠基書目；措辭軟化 + 作者補經典 |
| 2 | 甲基沉默例子 | BRCA2 罕見=VERIFIED；MLH1/CDKN2A 鐵證；BRCA1 中等 | 用戶拍板 **BRCA2+BRCA1 並用** |
| 3 | clonal evolution / ITH 臨床 | PARTIAL | ITH 存在有 Dentro；臨床句須補綜述 |
| 4 | 確認黃金標準 / 單樣本快照 | SUPPORTED | Tarabichi 2021 + Gaiti 2019；正名 practical guide |
| 5 | ONT vs 短讀 | SUPPORTED（4b 措辭收斂）| Smallwood/Dorado/LongPhase-S/Foltz；3 紅線 |
| 6 | 我們的數據真值 | VERIFIED（grep-confirmed）| 三檔交叉一致 |

> **方法論誠信**：本文件數據為 grep 真值；文獻分兩級；§6b/§6c 候選**未過 /citation-verification 前不可進論文**。HP 定義 canonical SoT = memory `reference_hp_tag_definition_and_subclone_caveat`。
