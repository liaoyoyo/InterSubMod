<!--
建立時間: 2026-06-20
報告類型: methodology — 體細胞位點甲基化可觀測量「組合 → 生物學詮釋」窮舉對應
任務類型: 研究判讀框架（延伸 20260617 situation inventory，加 LOH + phase-set lens）
build_branch: research/subclonal-reconstruction-202606 (HEAD c660179)
status: data-grounded（內部數字全溯源 situation inventory；外部 citation = L3 agent-fetched，投稿前過 /citation-verification）
data_sources:
  - docs/methodology/20260617_structure_label_situation_inventory_01.md (內部數字 SoT；verify_breakdown.py 25/25 PASS)
  - docs/methodology/_assets/20260616_readset_provenance/{label_breakdown,level2_calibration,tn_dependency,signal_tally,tumor_signal_tally,decomp_tn_v2}.json
  - workflow wr4u3752w (6-agent：4 文獻驗證 + 合成 + 對抗稽核 PASS) — 外部 citation 來源
provenance_note: 內部 metric 數字皆引自 20260617 situation inventory（§0-§7）；外部文獻 claim 為 agent WebFetch（L3），各帶 PMID/DOI，正式引用前須過 /citation-verification 從 PMC 全文逐字核。
adversarial_audit: workflow wr4u3752w evaluator = PASS（0 fabrication / 0 overclaim / 0 redline 違反 / 0 缺漏組合；4 non-blocking 修正已折入）。
-->

# 體細胞位點 × 甲基化可觀測量：組合窮舉 → 生物學詮釋對應

> **這是什麼**：延伸 `20260617_structure_label_situation_inventory`（三狀況 × 標籤 3 軸盤點），改用 **「單一位點的可觀測量組合」** 視角，回答：給定一個 somatic 位點，看到某個組合代表什麼？哪些與 clone/subclone 有關、哪些是混合多種可能、哪些不太可能、哪些需要更多資料；以及**同一 phase set 內 vs 跨 phase set** 的適用邊界。
> **紅線（全程遵守）**：(a) 無監督聚類+PERMANOVA = double-dip NEGATIVE，subclone 只走 a-priori conditioned 軸；(b) 確認 ≥2 群 ≠ 確認 subclone；(c) within_clean ≠ subclone（Jaccard 0.123）；(d) clean somatic cis 罕見（1/816）；(e) 甲基化虛實 = single-cell REAL / bulk read-level NOT-established / filter DEAD；(f) ISM = regional characterization，非 genome-wide reconstruction；**(g) 只 MAP GAPS，不下 ISM 貢獻宣稱**。

---

## §0 一頁總結（先看這裡）

**問題**：一個 somatic 位點，看到（CN context／HP-tag 組成／ALT 分布／甲基幾群／對齊哪條標籤軸／normal 歸因）的某組合 → 這代表什麼？

**三句話結論**：

1. **沒有任何單一組合能在單一 bulk 樣本上「確認 subclone」。** 最像 subclone 的組合（LOH × 半帶 ALT × k≥2 × 對齊 ALT）**全部 undecidable** — 因為四個競爭詮釋在讀數上退化：**I2 LOH-unmask（既存 ASM 顯形）／ I3 cis-ASM 單事件／ I4 真 subclone／ I7 over-cluster 假象**。
2. **「對齊一條獨立標籤軸」是唯一非循環的確認器**（CramérV≥0.7 + Fisher，避開 double-dip），**但「確認 ≥2 群 ≠ 確認 subclone」**：內部 clear≥2 之 **85% 對齊 germline/carrier**，多數其實是 cis-ASM；外部文獻佐證 — tumor 內 imprinted-DMR 甲基變化 **82–91% 由 CN/LOH（含 copy-neutral LOH）解釋，非 epimutation**（Martin-Trujillo 2017, PMID 28883545）。
3. **此窮舉只適用「同一 phase set（haplotype block）內」**。跨 phase set 在 **bulk 無 read/cell 連結**（MethPhaser PMID 38909018 明文：兩 block 間無法定 SNV 屬哪條 haplotype）→ 要串成 genome-wide lineage 須外部 **sample-level anchor**（CCF / CN-LOH segment / phylogeny）或 **single-cell** → 這是 reconstruction GAP，**不宣稱 ISM 解決**。

**四個 resolve 工具（缺一就 undecidable）**：① **normal cis-control**（砍 I2、分 germline vs somatic）② **multi-sSNV CCF gradient（chr2:18M 法）**（分 I3 單事件 vs I4 lineage、解半帶 multiplicity）③ **COLO829 跨樣本**（excess-over-null 復發、補單樣本封頂）④ **single-cell ground truth**（I4 終判，bulk 結構性不可獨立確證）。

**內部數字錨**（單一 HCC1395 paired run，30,490 TP SNV；全溯源 situation inventory §0-§7）：clear≥2 僅 **2,941（9.6%）**；其對齊 germline/carrier **85%** / HP1-HP2 **15%** / random **0%**；三態 subclone 決策 cant-ask **39.5%** / confirmed-1-group **45.4%** / multigroup-candidate **15.1%**；嚴格 within-HP PATTERN 僅 **354（1.16%）**。

---

> **框架前提（讀前必看）**：
> - **觀察單位** = per somatic locus，±5000bp 窗口，read×CpG 甲基矩陣（raw ML/255）→ **Bernoulli** read-read 距離（`delta=p(1−q)+(1−p)q`，權重 `w=2|p−0.5|`，`C_min=3`）→ **UPGMA** average-linkage → k 由 silhouette 在 `k∈[2,min(6,n/2)]` 取最大；unbalanced-bump k+1。**完全無監督**（標籤不進入聚類）。
> - **甲基化虛實（紅線 e）**：single-cell lineage methylation = **REAL**（Gaiti/Epiclomal/MethylTree/EPI-Clone）；**bulk read-level subclone = NOT-established**；methylation-as-filter = **DEAD**。
> - **scope（紅線 f）**：「subclonal reconstruction」對 ISM bulk long-read 是 **OVERCLAIM**；ISM 做的是 **REGIONAL CHARACTERIZATION（phase-block 尺度）**，非 genome-wide reconstruction。
> - 本窮舉**只適用於 PER phase set**（單一 haplotype block）；跨 phase set 整合是另一個更難的問題（見 §5，標為 GAP）。
> - **本文只 MAP GAPS（紅線 g）**：陳述什麼未知 / 需更多資料；貢獻判斷留給用戶，不在此草擬。

---

## §1 可觀測量維度表（D1–D7）

| 維度 | 取值 | 觀測的是什麼 | ISM 如何量測 |
|---|---|---|---|
| **D1** CN / genotype context | CN-neutral het-retained（2 haplotypes）/ LOH（1 haplotype 存活）/ CN-gain·amplification / CN-loss-non-LOH | 該位點窗口的拷貝數與雜合狀態；決定還剩幾個 germline allele 可作 anchor | ISM 自身**不算 CN**（無 copy-number/multiplicity model）；CN/LOH 屬**外部輸入**（需 CN/LOH segment）。內部以 germline-het SNP 雜合度 + HP-tag 構成可間接旁證 LOH（單 haplotype → 只見 1/1-1，無 2/2-1） |
| **D2** HP-tag composition（出現哪些） | {HP1, HP2, HP1-1, HP2-1, HP1-2, HP2-2, HP3, unphase} 的子集 | longphase-S `somatic_haplotag`：`1`/`2`=germline phase（normal-derived）；`1-1`/`2-1`=**第一條**支持 somatic ALT 且可溯源到 germline HP1/HP2 的 somatic haplotype；suffix `-2`（即 `1-2`/`2-2`）=**第二條** somatic haplotype = **MULTI-SUBCLONE marker**；`3`=unassigned；`unphase`=no phase | 從 BAM HP-tag 直接讀取（`HaplotagType.h:327-337` + ReadParser 源碼驗證）。⚠ HCC1395 區域（如 chr2:18M）只見 `1/2/1-1/2-1/3`，**從未見 `-2`**（HCC1395 該主題缺第二 somatic haplotype 證據） |
| **D3** ALT（somatic）分布 | ALT clonal on one somatic-haplotag（carrier reads 全帶）/ ALT subclonal within a haplotag（**只一部分**帶）/ ALT balanced·absent | somatic ALT 在 carrier reads 內的滲透率 → 是否「全帶」vs「半帶」 | 由 read 層 ALT 計數 / HP-family 內 ALT-vs-REF 拆分（`AltSubcloneDbeta`，**TUMOR-ONLY** `!is_tumor continue`）。⚠ 「半帶」formally **NON-IDENTIFIABLE** between {multiplicity=1 clonal / 真 subclone CCF<1 / ADO·error}（Tarabichi 2021, DeCiFer 2021）——ISM **無 multiplicity model 無法升格為 CCF** |
| **D4** 甲基 cluster 數 k | 1（no split）/ 2 / 3+ | UPGMA 聚類在該窗口分出幾群 read | silhouette 在 `k∈[2,min(6,n/2)]` 取 max；unbalanced-bump k+1。⚠ silhouette + k+1 bump 在 imbalance/outlier 下**已知 over-cluster**（gap statistic 文獻 Tibshirani 2001）；read-internal CpG 相關（DNMT1 processivity ~36bp, PMID 36285438）亦灌 k（紅線 c） |
| **D5** structure test | PERMANOVA sig + location（真）/ PERMANOVA sig 但 PERMDISP dispersion-driven（artifact）/ not sig | 群間差異是「位置（中心）不同」還是「離散度（spread）不同」 | StructureTest.hpp 內 `check_dispersion`（PERMDISP）已實作；ISM 內部以 PERMDISP 攔 PERMANOVA 假象（**Dispersion 1,700 被攔**屬 REMOVE）。⚠ 純無監督「挑最分得開的群再測」=**double-dip / 循環**（紅線 a；Kriegeskorte 2009；noise approx structure 83–100% under any null） |
| **D6** alignment（甲基 cluster 對齊哪條 a-priori 軸）+ CramerV 強度 | allele ALT/REF / HP-family / HP-fine carrier-vs-germline / none（cluster⊥label） | 用**獨立**的 a-priori 標籤檢驗聚類是否「跟著標籤走」（避開 double-dip） | cluster×label 列聯表 GlobalP（**PAIRED tumor+normal**，GlobalTest 無 is_tumor 過濾）；clear≥2 用 **CramerV≥0.7 + Fisher sig**。⚠ 對齊是 paired 訊號（混 normal reads），Level-2 tumor-only 對齊保持 51.6%，median V 0.665→0.319（源 doc §0/§6.4） |
| **D7** normal attribution | germline-cis（在 normal 即有）/ residual（tumor-specific）/ unattributed / SampleASM-only | 甲基差異是否 normal 也有（=germline，非 somatic）vs 相減後 tumor-specific | `GermlineAsmDbeta`=**NORMAL-only**（定義 germline）；`SomaticResidual/HP_Residual`=**NEEDS BOTH**（tumor−normal）；`SampleASM`=tumor-vs-normal 整體差（**飽和命中 96.53%（29,433），無區辨力，勿引用為訊號**，源 doc §0.5）。⚠ germline 與 residual **意義相反，禁合併** |

---

## §2 詮釋類別表（I1–I8）

| 代碼 | 名稱 | 定義 | ISM anchor | 類別 |
|---|---|---|---|---|
| **I1** | germline ASM / imprinting | cis、在 normal 即存在、在**兩條 haplotype 之間**的等位特異甲基 | `GermlineAsmDbeta`（NORMAL-only）；D7=germline-cis；需兩條 allele 對齊。LOH 下 **I1 不可能**（只剩 1 haplotype） | **unlikely**（指「在 LOH 區作為 subclone 解釋」時不可能；非 LOH 區則為常見良性背景） |
| **I2** | LOH-unmasked pre-existing/imprinted ASM | LOH-special：存活 haplotype 自身**遺傳來的** ASM 因對側 allele 消失而**顯形**——**非新 somatic 事件** | D1=LOH + D7 需 **normal cis-control** 才能區分（normal 該 allele 即帶此甲基）；甲基 magnitude 隨 LOH allelic-imbalance scale（Chase 2015, r=0.76 signature） | **mixed**（會偽裝成 within-HP 多群；需 normal cis-control 才能與 I3/I4 分開） |
| **I3** | cis-ASM somatic single event | ALT 在**局部**改變甲基；**單一事件**，非 lineage | `SubcloneDbeta`/`AltSubcloneDbeta`（TUMOR-ONLY）對齊 ALT carrier-vs-non-carrier；clean somatic cis **罕見 1/816**（chr17/TBC1D16，紅線 d） | **mixed**（與 I4 在讀數上退化；需 multi-sSNV linkage + CCF gradient 才能與 lineage 分開） |
| **I4** | somatic subclone | 多個 cell lineage；甲基與 **SUBCLONAL ALT** 共分離；tumor-specific | a-priori conditioned 軸（haplotag carrier vs germline，扣 normal baseline）；`1-2`/`2-2` 第二 somatic haplotype 為 marker（HCC1395 缺）。⚠ bulk read-level subclone NOT-established（紅線 e） | **subclone_related**（但 bulk 單樣本**無法獨立確證**；需外部 anchor / single-cell） |
| **I5** | somatic single-clone characterization | 一條 somatic haplotype、tumor-specific 甲基、**非 multi-subclone** | 對齊 carrier-vs-germline 但 k=2 且無第二 somatic haplotype（無 `-2`）；屬 regional characterization（紅線 f） | **subclone_related**（屬「單 clone 刻畫」非多 subclone 重建；是 I4 的保守替代，常為更可能解） |
| **I6** | CN / dosage artifact | 甲基「差異」由拷貝數/劑量改變製造，非真表觀事件 | D1=CN-gain/loss；keep 分類中 CN=1,486；HP-axis held-const CN 通道（confound pilot 排除） | **unlikely**（作為 subclone 解釋時；但作為 confound channel 真實存在，須排除） |
| **I7** | over-clustering / epiallele / stochastic | k 被 read-internal 相關灌大；無真 lineage | Bernoulli 距離 + UPGMA + silhouette 在相關 CpG 下 over-cluster；within_clean≠subclone（Jaccard 0.123, 55% epiallele，紅線 c）；DNMT1 processivity ~36bp（PMID 36285438）；stochastic disorder/PDR（Landau 2014） | **unlikely**（作為真結構時）／ 但**作為 null 必須先排除**——是 k≥2 的預設解釋直到證偽 |
| **I8** | technical（coverage imbalance / strand bias / mapping） | 覆蓋不均、strand bias、mapping artifact 製造假群 | D5 PERMDISP 攔 dispersion-driven；reads.tsv 側欄（HP/ALT/T-N/Strand）肉眼檢視；homopolymer（chr2:18M pos4 20bp poly-T）製造假 read 群 | **unlikely**（作為生物訊號時）／ 屬技術 null，須先排除 |

> **類別語義**：`subclone_related` = 該詮釋若成立則與 clone/subclone 有關；`mixed` = 在讀數上與 subclone 退化、需額外資料才能分開；`unlikely` = 作為 subclone 解釋時不太可能（但作為 confound/null 真實，須排除）；`needs_more_data` = 單一 paired 樣本不可判定。

---

## §3 組合 → 詮釋對應（最具資訊量的 ~20 組）

> 組織順序：**D1 context → D3 ALT 分布 → D4 k → D6 alignment**。每列給 plausible 詮釋集 + 類別 + resolve 所需額外資料。
> ⚠ 全部**只適用 PER phase set**（§5）。內部佐證數字見 §0 / 源 doc。

### A. CN-neutral het-retained（D1，兩條 haplotype 都在 → I1 可能）

| # | D3 ALT | D4 k | D6 alignment | plausible 詮釋集 | 類別 | resolve 需要 |
|---|---|---|---|---|---|---|
| A1 | balanced/absent | 1 | — | 平坦無結構（§0.5 定義A 無結構 5,725=18.78%，皆≥20 reads=真平坦） | unlikely(任何 subclone) | — |
| A2 | balanced/absent | 2 | **HP-family**（兩 haplotype 間） | **I1** germline ASM/imprinting | unlikely(subclone)；屬良性 | normal cis-control 確認在 normal 即有；intersect Rosenski 2025 460 parental-ASM + 34,426 SD-ASM 座標 |
| A3 | balanced/absent | 2 | none（cluster⊥label） | **I7** over-cluster / epiallele（狀況③ 510=1.67%） | unlikely | correlation-preserving null（M3C 原理）；Bayesian epiallele count（Barrett 2017） |
| A4 | clonal on one HP-tag | 2 | **HP-fine carrier-vs-germline** | **I3** cis-ASM single event / **I5** single-clone char. / I7 | mixed | 排 I7：epiallele null；分 I3 vs I5/I4：multi-sSNV linkage + CCF gradient（chr2:18M 法） |
| A5 | **subclonal within HP-tag（半帶）** | 2 | **allele ALT/REF** | **I4** subclone / **I3** cis-ASM / multiplicity=1 假象 | mixed→subclone候選 | CCF gradient（VAF+CN+purity）；multi-sSNV co-seg；single-cell |
| A6 | subclonal（半帶） | 3+ | allele ALT/REF + 部分 none | **I4** / **I7**（k 灌大） / I3 | mixed | epiallele null 先壓 I7；存活群再對齊 ALT 才談 I4 |
| A7 | clonal | 2 | **HP-family**（cluster 跟著 germline HP 而非 ALT） | **I1** germline ASM（與 somatic 無關）/ I6 | unlikely(subclone) | normal cis-control；CN 檢查排 I6 |

### B. LOH（D1，只剩 1 haplotype → **I1 不可能**，候選收斂 I2/I3/I4/I7）

| # | D3 ALT | D4 k | D6 alignment | plausible 詮釋集 | 類別 | resolve 需要 |
|---|---|---|---|---|---|---|
| B1 | balanced/absent | 1 | — | 平坦（LOH 區無結構） | unlikely | — |
| B2 | balanced/absent | 2 | none / 唯一 HP-fine | **I2** LOH-unmask（最 parsimonious）/ **I7** stochastic | mixed | **normal cis-control**（Martin-Trujillo 82–91% 由 CN/LOH 解釋）；甲基 magnitude 是否隨 LOH BAF scale（Chase r=0.76） |
| B3 | balanced/absent | **3+** | none | **I7** over-cluster（**HIGH artifact risk**）/ I2 | unlikely(subclone) | epiallele/PDR null（Landau）；Bayesian epiallele（Barrett）；單 haplotype 無 allele-anchor 對齊→確認受限 |
| B4 | clonal on one HP-tag | 2 | **HP-fine carrier-vs-germline** | **I3** cis-ASM single / **I2** unmask / **I5** single-clone | mixed | normal cis-control 分 I2；multi-sSNV+CCF 分 I3 vs I5/I4 |
| B5 | **subclonal（半帶）** | 2 | **carrier-vs-non-carrier** | **I4** subclone / **I3** cis / **I2** unmask / multiplicity=1 / ADO | **mixed→subclone候選** | (1)normal cis-control 排 I2；(2)multi-sSNV linkage+CCF gradient 排 I3、定 I4；(3)single-cell 終判 |
| **B6** | **subclonal（半帶），HP-tags 僅 HP1+HP1-1** | **3** | **同時對齊 ALT 與 HP labels** | **見 §4 WORKED EXAMPLE**（I2/I3/I4/I7 並列） | **mixed→subclone候選** | §4 decisive data |
| B7 | clonal | 2 | **allele ALT/REF**（存活 allele 內部仍見 ALT/REF 拆分） | I3 / I2 / I7 | mixed | normal cis-control；注意 LOH 已移除對側 allele anchor（最乾淨 discriminator 消失） |
| B8 | 任意 | 2/3+ | sig 但 **PERMDISP dispersion-driven** | **I8** technical / I7 | unlikely | check_dispersion 已攔（Dispersion 1,700 REMOVE）；coverage/strand 檢查 |

### C. CN-gain / amplification（D1）

| # | D3 ALT | D4 k | D6 alignment | plausible 詮釋集 | 類別 | resolve 需要 |
|---|---|---|---|---|---|---|
| Cg1 | subclonal（半帶） | 2/3+ | allele ALT/REF | **I6** dosage artifact / **I4** subclone / I3 | mixed | CN-aware multiplicity（半帶可由 gain 造成非 subclone）；HP-axis held-const CN；CCF |
| Cg2 | clonal | 2 | HP-family | I6 / I1 / I3 | mixed | CN segment + normal cis-control |

### D. CN-loss-non-LOH（D1，仍有兩 allele 但劑量降）

| # | D3 ALT | D4 k | D6 alignment | plausible 詮釋集 | 類別 | resolve 需要 |
|---|---|---|---|---|---|---|
| Dl1 | balanced | 2 | HP-family | I1 / I6 | unlikely(subclone) | normal cis-control + CN |
| Dl2 | subclonal（半帶） | 2 | carrier-vs-non | I4 / I3 / I6 | mixed | CCF（CN 校正後）；multi-sSNV |

### E. 通用 needs_more_data 觸發（任何 D1）

| # | 觸發條件 | 類別 | 為何不可判 |
|---|---|---|---|
| E1 | D4 k≥2 但 D6=none（cluster⊥label） | needs_more_data | 510(1.67%) 狀況③；無 a-priori 對齊→無法分 I7 vs 真結構 |
| E2 | D5 not sig 但 mean-β 雙峰（LEVEL low-confidence，主導 7,910） | needs_more_data | LEVEL 非 silhouette-verified（PATTERN 嚴格僅 354=1.16%） |
| E3 | D7=unattributed（① 4,698=15.4%） | needs_more_data | tumor 偵測到結構但 germline-ASM 與 residual 軸皆不顯著→未定 germline/somatic |
| E4 | D3 半帶但 carrier reads<3（carrier-limited 414=1.36%，源 doc §0.6） | needs_more_data | 低 VAF→somatic-ALT reads<3，subclone 軸測不了 |

> **跨表規律**：(1) D6=**HP-family 對齊** → 偏 I1（germline，非 subclone）；(2) D6=**carrier-vs-non/ALT-REF 對齊 + D3 半帶** → subclone 候選（但需排 I2/I3/multiplicity）；(3) D4=**3+ 在單 haplotype 且 D6=none** → I7 預設直到對齊證偽（紅線 a/c）；(4) **clear≥2 對齊 germline/carrier 85%** → 多數「確認≥2 群」其實是 cis-ASM 不是 subclone（紅線 b）。

---

## §4 WORKED EXAMPLE（用戶案例）

> **觀測**：LOH 區；出現的 HP-tags **只有 HP1 + HP1-1**；HP1-1 內**只有一半**帶 ALT；甲基聚成**剛好 3 群**，與 ALT 及 HP 標籤對齊。

### 4.1 逐維度判讀（D1–D7）

| 維度 | 本案取值 | 判讀 |
|---|---|---|
| **D1** | **LOH**（只 HP1 系存活，無 HP2/HP2-1） | **關鍵**：只剩 1 條 germline haplotype → **I1（兩 haplotype 間 germline ASM）IMPOSSIBLE**。候選收斂到 **I2 / I3 / I4 / I7** |
| **D2** | {HP1, HP1-1}；**無 HP2/HP2-1，無 1-2/2-2** | 與 LOH 自洽（單系）；**無 `-2` suffix → 無第二 somatic haplotype marker** → multi-subclone 的最直接證據缺席 |
| **D3** | **HP1-1 內半帶 ALT** | 半帶 = subclonal-ALT / multiplicity signal。**formally NON-IDENTIFIABLE** between {multiplicity=1 clonal-on-one-copy / 真 subclone CCF<1 arose-after / ADO·error}（Tarabichi 2021, DeCiFer 2021）。ISM **無 multiplicity model → 不能升格為 CCF 或確認 subclone** |
| **D4** | **k=3** | 單一 haplotype 內 3 群 = **over-clustering HIGH risk**（silhouette + k+1 bump；Şenbabaoğlu chance-partition；read-internal CpG 相關）。**UNLESS** 切法對齊 ALT carrier-vs-non-carrier，則 I7 降權、I3/I4 升權 |
| **D5** | 需檢查 | 必跑 PERMDISP：若 dispersion-driven → 退回 I7/I8。本案需確認是 location（真）非 spread |
| **D6** | **3 群對齊 ALT + HP labels** | 對齊 = 避開 double-dip 的正確檢驗（紅線 a）。但 **clear≥2 對齊 germline/carrier 85% → 對齊≠subclone**（紅線 b），多數是 cis-ASM |
| **D7** | **LOH → cant-ask（無 control，39.5% 類）** | 單 haplotype + LOH 移除對側 allele anchor → normal-anchored cis-control 此處塌縮為「cant-ask 12,039(39.5%)」。**這是分 I2 與 I3/I4 的唯一鑰匙，卻在 LOH 最難取得** |

### 4.2 全部 plausible 詮釋（**不收斂成單一**），按可能性排序

| 排名 | 詮釋 | 為何可能 | 會 confirm 它的證據 | 會 refute 它的證據 |
|---|---|---|---|---|
| **1** | **I7 over-clustering / epiallele**（k=3 在單 haplotype 的預設 null） | silhouette+k+1 在 imbalance 下 over-cluster；read-internal CpG 相關（DNMT1 ~36bp）；單 clone 即可有 stochastic disorder（Landau 高 PDR）。Jaccard 0.123 / 55% epiallele | correlation-preserving null 下 k=3 **不**勝過 stochastic baseline；Bayesian epiallele（Barrett）估出 ~1 dominant epiallele+noise | k=3 切法**穩定對齊 ALT carrier-vs-non**，且勝過 epiallele-aware null → I7 被排除 |
| **2** | **I2 LOH-unmasked pre-existing/imprinted ASM**（非新 somatic 事件） | LOH 區最 parsimonious（Martin-Trujillo 82–91% 由 CN/LOH 解釋非 epimutation）；存活 haplotype 帶遺傳來的 ASM 顯形 | **normal cis-control**：normal 該 allele 即帶此甲基 pattern；甲基 magnitude 隨 LOH BAF scale（Chase r=0.76）；位點落 Rosenski 460 parental-ASM / 11p15·MEG3-DLK1·RB1 blocklist | normal 完全無此甲基 pattern，且 pattern 與 **somatic ALT carrier** 共分離（非與 germline genotype） |
| **3** | **I3 cis-ASM somatic single event**（一個事件，非 lineage） | ALT 局部改甲基；對齊 ALT carrier 可解 k 之一分界 | normal 無此甲基（排 I2）+ 甲基與 ALT carrier 共分離；但 **無 CCF gradient、無 multi-sSNV linkage** | 多 sSNV 顯示 CCF 階梯/巢狀 lineage 結構 → 升 I4 |
| **4** | **I4 somatic subclone**（多 lineage 共分離） | 半帶 ALT + 3 群對齊 ALT/HP → 表面像 subclone | **multi-sSNV linkage + CCF gradient（chr2:18M 法）** 顯示 lineage 巢狀；single-cell ground truth | 半帶可由 multiplicity=1 或 ADO 解釋；無 `-2` 第二 somatic haplotype；single sample bulk **無法獨立確證**（紅線 e/f） |
| (附) | **I8 technical** | LOH 區 coverage 不均 / homopolymer（chr2:18M pos4 20bp poly-T 製造假群）/ strand bias | PERMDISP=dispersion-driven；reads.tsv 側欄 HP/ALT/T-N/Strand 肉眼見不均 | location-driven 且側欄均衡 → 排除 |

### 4.3 關鍵推理鏈（用戶指定要點，逐條落實）

1. **LOH ⇒ I1 不可能**：只剩 1 haplotype，兩-haplotype-間 germline ASM 無從談起 → 候選收斂 **I2/I3/I4/I7**。
2. **「HP1-1 半帶 ALT」= subclonal-ALT / multiplicity 訊號**：可能是 **I4**，但**同樣可能**是 phasing-assignment 誤判 / ADO / multiplicity=1（NON-IDENTIFIABLE，Tarabichi/DeCiFer）。**不可單憑半帶宣稱 subclone**。
3. **單 haplotype 內 3 群 = over-clustering 風險（I7）**——**UNLESS** 切法對齊 ALT carrier-vs-non-carrier，則轉為 **I3 或 I4**（對齊避開 double-dip，紅線 a）。
4. **分 I2（LOH-unmask）需 normal cis-control**：normal 該存活 allele 是否本來就帶此甲基（Martin-Trujillo/Chase 簽名）。
5. **分 I3（單 cis 事件）vs I4（lineage）需 multi-sSNV linkage + CCF gradient（chr2:18M 法）**。

### 4.4 decisive 額外資料（單點突破優先序）

1. **normal cis-control（同位點 normal reads）** — 先砍 I2：若 normal 即帶此甲基 → I2，全案降為「LOH-unmask 良性顯形」非 somatic。**LOH 下這正是 cant-ask 39.5% 的核心困難**。
2. **epiallele-aware / correlation-preserving null + Bayesian epiallele count** — 砍 I7：k=3 是否勝過 stochastic baseline。
3. **multi-sSNV linkage + CCF gradient（chr2:18M 法）** — 分 I3 vs I4：是否有 lineage 巢狀 CCF 階梯。
4. **single-cell ground truth** — 終判 I4（bulk 單樣本對 I4 結構性不可獨立確證，紅線 e/f）。

> **本案誠實結論（紅線 g，不下貢獻宣稱）**：在**單一 paired 樣本**上，此組合**無法判定** I2/I3/I4/I7 中何者為真——最 parsimonious 起點是 **I7（over-cluster）∨ I2（LOH-unmask）**，**subclone（I4）是 last resort**，需上述 1→4 序貫排除。ISM 此處能做的是 **regional characterization**（描述 epiallelic/ASM 結構），**非 subclone 確認**。

---

## §5 PHASE-SET 分析（窮舉的適用邊界）

### 5.1 within a SINGLE phase set — 可解讀

- 一個 phase set = 一個 haplotype block，reads 由 germline SNP 連結成 read-connectivity connected component（VCF/WhatsHap 定義：PS = connected components；標籤 arbitrary）。
- ONT phase block 尺度 = **kb–Mb**（tool N50 ~13–25 Mb；per-block median ~80–590 kb）。
- ISM 觀察窗 = **±5000bp（10 kb）** → **比 phase block 小 1–2 個數量級** → **每個 ISM 窗口必落在單一 block 內部**，其 HP1/HP2 / 1-1/2-1 標籤在窗內**局部自洽且可比**。
- 因此：**within-block 的 HP-tag + within-block 甲基共變是有意義的 → within-block subclone 候選可解讀**（§3/§4 的窮舉適用於此層）。

### 5.2 ACROSS phase sets — 無連結（GAP）

- **MethPhaser（PMID 38909018）**：「between two phase blocks, the assignment of SNV to haplotype 1 or 2 cannot be determined as no connection information is available」。
- 故 **block A 的 HP1 不保證是 block B 的 HP1 同一條物理 homolog**；在 **bulk** 中 reads **無法溯源到單一細胞**（無 cell-level linkage）。
- **per-block call 不能直接串接成 genome-wide lineage**——除非有**外部 SAMPLE-level anchor**：CCF clustering（VAF+CN+purity）、CN/LOH segments、phylogeny（Tarabichi 2021），或 allele-specific-methylation block-bridging（MethPhaser/HapBridge，但那是 germline-ASM 用途，**非** subclone lineage）。
- genome-wide 串接的 established 解只在 **single-cell**（MethylTree ~100% / EPI-Clone / Gaiti / Epiclomal）。

### 5.3 明確聲明

> **本窮舉（§3/§4）只適用 PER phase set。** Cross-phase-set 整合是**另一個、更難的問題** = **genetic-anchor / reconstruction gap**。**這是 GAP，不宣稱 ISM 解決它**（紅線 f/g）。ISM 做 phase-block-scale REGIONAL characterization；要變 genome-wide lineage statement 需外部 anchor（CCF / CN-LOH / phylogeny）或 single-cell——皆 ISM 此 pipeline 不具備。

---

## §6 「needs more data」彙總（單一 paired 樣本上 undecidable 的組合）

| undecidable 組合 | 為何不可判（單樣本） | 會 resolve 它的資料 |
|---|---|---|
| LOH 區任何 within-HP 多群（§4 B6 / B2 / B3） | LOH 移除對側 allele anchor → cis-control 塌縮為 cant-ask（**39.5%, 12,039**）；I2/I3/I4/I7 退化 | **normal cis-control**（同位點 normal reads）+ Rosenski/imprinted blocklist intersect |
| 半帶 ALT → subclone vs multiplicity=1 vs ADO（A5/B5/B6/Cg1/Dl2） | read fraction 對 multiplicity formally NON-IDENTIFIABLE（Tarabichi/DeCiFer）；ISM 無 CN/multiplicity model | **multi-sSNV CCF gradient（chr2:18M 法）** + CN segment 校正 |
| k≥2 cluster⊥label（狀況③ **510, 1.67%**） | 無 a-priori 對齊→無法分 I7 vs 真結構 | correlation-preserving null + Bayesian epiallele count |
| D7=unattributed（① **4,698, 15.4%**） | tumor 見結構但 germline/residual 軸皆不顯著 → germline/somatic 未定 | **COLO829 cross-sample** 復發檢查 + normal 深覆蓋 |
| 「對齊 germline/carrier」是否 subclone（clear≥2 之 **85%**） | 對齊≠subclone，多為 cis-ASM（紅線 b）；clean somatic cis 罕見（**1/816**） | normal-anchored cis-control + multi-sSNV linkage |
| I4 subclone 之終局確認（任何上述升到 I4 候選者） | bulk read-level subclone NOT-established（紅線 e）；single sample 結構性不可獨立確證（紅線 f） | **single-cell ground truth**（MethylTree/EPI-Clone class） |

**四大 resolve 工具一覽**：
1. **normal cis-control** → 砍 I2（LOH-unmask）、分 germline vs somatic。
2. **COLO829 cross-sample** → 跨樣本復發（excess-over-null）、補單樣本封頂。
3. **multi-sSNV CCF gradient（chr2:18M 法）** → 分 I3 單事件 vs I4 lineage、解半帶 multiplicity。
4. **single-cell ground truth** → I4 終判（bulk 不可獨立確證）。

> **總結（紅線 g）**：在單一 paired 樣本上，最具 subclone 嫌疑的組合（LOH × 半帶 ALT × k≥2 × 對齊）**全部** undecidable；ISM 能交付的是**有 documented confounds 的 regional characterization**，subclone/reconstruction 的貢獻判斷**留給用戶**。

---

## §7 外部文獻錨（L3 — agent WebFetch，投稿前過 /citation-verification）

> ⚠ 以下為 workflow wr4u3752w 之 4 文獻驗證 agent 一手 fetch（WebSearch/WebFetch = **L3**）；各帶 PMID/DOI。**正式引用前須過 `/citation-verification` 從 PMC OA 全文逐字核**（含 author list、年份、頁碼）。部分 DOI 來自 search-result 抽取（如 MethylBERT bioRxiv↔Nat Commun 雙版、Wolf processivity 標題作者）需優先重核。

### 7.1 LOH-unmask ASM（I2 的文獻底座 — 本輪最關鍵新增）
- **Martin-Trujillo A, et al.** Copy number rather than epigenetic alterations are the major dictator of imprinted methylation in tumors. *Nat Commun* 2017;8:467. **PMID 28883545**. → tumor imprinted-DMR 甲基變化 **82–91% 由 CN/LOH 解釋**，非 epimutation。【I2 最強錨】
- **Chase A, et al.** Profound parental bias associated with chromosome 14 acquired uniparental disomy. *Leukemia* 2015;29(10):2069-74. **PMID 26114957**. → 14q MEG3-DLK1 甲基隨 aUPD 程度 scale，**r=0.76**（unmasking 簽名）。
- **Rosenski J, et al.** Atlas of imprinted and allele-specific DNA methylation in the human body. *Nat Commun* 2025;16:2455. **PMID 40069157**. → **460 parental-ASM + 34,426 SD-ASM** 座標可作 intersect blocklist；SD-ASM = genotype-driven cis-ASM 正式名。
- **Steenman MJC, et al.** *Nat Genet* 1994;7:433-9. **PMID 7920659**（+ Wilms UPD survey PMID 18464243）→ 11p15/IGF2-H19：UPD 與 epimutation 在甲基讀數上退化。
- **Kanber D, et al.** The human retinoblastoma gene is imprinted. *PLoS Genet* 2009;5(12):e1000790. **PMID 20041224**（+ Greger PMID 8163810）→ RB1 intron-2 imprinted DMR。
- **O'Keefe C, McDevitt MA, Maciejewski JP.** Copy neutral LOH: a novel chromosomal lesion. *Blood* 2010;115(14):2731-9. **PMID 20107232**. → cnLOH 為通用 homozygosing lesion，帶任何既存 ASM forward。
- **Gigante S, et al.** Using long-read sequencing to detect imprinted DNA methylation. *NAR* 2019;47(8):e46. **PMID 30793194**（+ haplotype-dependent ASM *Nat Commun* 2020;11:5238, PMC7567826）→ ASM 偵測需兩 allele + allele-of-origin；LOH 移除對側 allele = 最乾淨 discriminator 消失。

### 7.2 within-haplotype bulk subclone = NOT-established（I4 / 半帶 ALT）
- **Tarabichi M, et al.** A practical guide to cancer subclonal reconstruction. *Nat Methods* 2021;18:144-55. **PMID 33398189**（+ Dentro principles PMC5538405）→ multiplicity 歧義、CCF 框架。
- **Satas G, et al.** DeCiFering the elusive cancer cell fraction. *Cell Syst* 2021;12(10):1004-18. **PMID 34416171**. → 半帶 ALT multiplicity NON-IDENTIFIABLE。
- **Bonfiglio S, et al.** cvlr（*Bioinform Adv* 2023, DOI 10.1093/bioadv/vbac101）+ **ASMS**（bioRxiv 2024, DOI 10.1101/2024.12.18.629129）→ bulk ONT read-clustering 最近似物；**ASMS 有 permutation test**（修正舊 memory「缺顯著性」口徑）；但無 ASM-vs-subclone discriminator。
- **Jeong Y, et al.** MethylBERT. *Nat Commun* 2025;16:788. **PMID 39824848**. → supervised tumor/normal deconvolution（對齊「subclone 只走 a-priori 軸」）。
- **Wolf G, et al.** DNMT1 processivity（PMC9597173, **PMID 36285438**）+ Busto-Moner（PLoS Comput Biol 2020）→ DNMT1 ~36bp 連續甲基化 = read-internal CpG 相關機制根因（紅線 c）。
- single-cell REAL：**Gaiti** *Nature* 2019 PMID 31748746 / **Epiclomal** PMID 33001985 / **MethylTree** PMID 39820752 / **EPI-Clone** *Nature* 2025。physical 預分選：**Liu/Goretsky** bioRxiv 2025 PMID 40950124。

### 7.3 phase-block 跨 block 無連結（§5 GAP）
- **Tan et al.** MethPhaser. *Nat Commun* 2024;15:5327. **PMID 38909018**. → 「between two phase blocks … no connection information is available」直接引句。
- **Akbari V, et al.** NanoMethPhase. *Genome Biol* 2021;22:68. **PMID 33618748**. → megabase-scale phasing；block N50 尺度。
- **Cornish et al.** Determinants of phasing accuracy. bioRxiv 2026 + **Zhou et al.** HapBridge bioRxiv 2025（block-bridging = germline-ASM 用途非 subclone）。
- VCF spec / WhatsHap：PS = read-connectivity connected components（arbitrary labels）。

### 7.4 over-clustering / k 灌大 + 如何確認切群為真（I7 / D4）
- **Şenbabaoğlu Y, et al.** Critical limitations of consensus clustering. *Sci Rep* 2014;4:6207. **PMID 25158761**. → cluster-less 資料被切成「看似穩定」群；valid null 須保 feature 相關。
- **John CR, et al.** M3C: Monte Carlo reference-based consensus clustering. *Sci Rep* 2020;10:1816. **PMID 32020004**. → 程序偏高 K、從不測 K=1；correlation-preserving null 原理。
- **Landau DA, et al.** Locally disordered methylation … intratumor methylome variation. *Cancer Cell* 2014;26(6):813-25. **PMID 25490447**. → 單 clone 內即有 stochastic disorder（高 PDR）= 多 read 群零 subclone。
- **Scherer M, et al.** within-sample heterogeneity scores. *NAR* 2020;48(8):e46. **PMID 32103242**（+ Landan epipolymorphism *Nat Genet* 2012, PMID 23064413）。
- **Kimes PK, et al.** SigClust2/SHC. *Biometrics* 2017;73(3):811-21. **PMID 28099990**（+ Liu 2008 SigClust）→ 1-vs-≥2 正規檢定框架，**但 Gaussian null 不合 Bernoulli/epiallele 甲基**。
- **Barrett JE, et al.** Bayesian epiallele detection. *BMC Bioinformatics* 2017;18:354. **PMID 28743252**. → model-based epiallele 計數（正解方向，但 bisulfite multi-region，未驗 single-run ONT 5mC）。
- **Kriegeskorte N, et al.** Circular analysis / double dipping. *Nat Neurosci* 2009;12:535-40. **PMID 19396166**. → 無監督「挑最分得開再測」= 循環（紅線 a 理論根據）。

---

*內部數字 SoT：`InterSubMod/docs/methodology/20260617_structure_label_situation_inventory_01.md`（verify_breakdown.py 25/25 PASS）。本 doc 為其組合視角延伸。*
*關聯 memory：`project_somatic_locus_methylation_combination_interpretation`（本 doc 摘要）· `project_cluster_label_alignment_readset_paired`（situation inventory）· `project_subclone_snv_difficulty_methylation_framework`（SNV-only/LOH-CNV 困難）· `project_chr2_18m_subclone_locus_verification`（multi-sSNV CCF 法）· `project_apriori_subclone_classification_model`。*
