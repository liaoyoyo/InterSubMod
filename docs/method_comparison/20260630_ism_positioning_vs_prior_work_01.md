<!--
建立時間: 2026-06-30
報告類型: method_comparison — ISM 定位「與過去研究的差別」(positioning vs prior work；論文 Related Works/Discussion 防守底稿)
任務類型: D handoff — 論文定位口徑單一真值
build_branch: research/subclonal-reconstruction-202606
status: 定位定稿（用戶 2026-06-30 crystallize 設定+方法+3 精準點）；文獻裁決=74-源庫+7-競品全文核 + agent a702cc1c 確認「ISM 白地為空」（HCC1395/COLO829 真值皆短讀/單細胞非 ONT；單分子共現骨幹=最清楚無佔者）；🔴 投稿前 /citation-verification 優先序：PCDH(403)>McAuley(403)>Fang2021 subclone 數(Fig4 二手)
data_sources:
  - memory project_subclone_confirmation_construction_ont_transferability（白地裁決 + 移植）
  - memory project_clone_subclone_landscape_and_ism_feasibility（6 學派 + 競品）
  - memory project_subclone_snv_linkage_verification_pipeline（somatic 共現=唯一非循環骨幹）
  - external_validation 74 源（競品卡：oneill/tumorlens/rocit/pcdh/mcauley/foltz/lee/wakhan/severus）
provenance_note: 文獻「無 ONT 單樣本 subclone 重建」裁決源自 74-源庫 + 7-競品全文核（verified）；HCC1395-specific fresh check = agent a702cc1c（進行中）。3 精準點 = 用戶定稿，與既有紅線一致。
-->

# ISM 定位：與過去研究的差別（positioning vs prior work）

> **這份是什麼**：論文 Related Works / Discussion 的**定位防守底稿**。把 ISM 的精確設定、方法、與過去研究的 3 個對比軸、3 個決定能否扛 reviewer 的精準點、ONT 三優勢，整理成單一口徑真值。
> **可信度**：文獻裁決 = 74-源庫 + 7-競品全文核（verified）；外部 citation = **L3**（投稿前 /citation-verification）。

---

## §0 ISM 的精確設定 + 方法（用戶 2026-06-30 定稿）

- **設定**：**單樣本（HCC1395 primary）、單一時間點、tumor + matched normal、ONT 長讀**。
- **方法**：從 **somatic 共現骨幹**重建 subclone 結構（**區域 partition + nested/sibling 關係**）；用 **haplotag 鑑別 allelic vs clonal**；**甲基輔助 characterize**（bounded-auxiliary）。
- **對齊內部真值**：somatic sSNV 共現 = **唯一非循環骨幹**（memory `project_subclone_snv_linkage_verification_pipeline`）；HP = **鑑別器非確認器**；甲基用途已窮盡 = **bounded-auxiliary**（負篩 + L3 旗標，非確認）。

---

## §1 文獻裁決：有沒有人做過 ONT 單樣本 subclone 重建？

**裁決：NO —— 無已發表方法在「單一 bulk ONT 樣本（tumor+normal、單時間點）」上、用「單分子共現骨幹」重建 clone/subclone 結構。** ISM 的精確交集是空的（但鄰居在 2024-26 快速逼近，故 framing 須誠實窄）。

**最近鄰 + 各自為何不是同一件事**：
| 鄰居 | 平台 | 做什麼 | 為何不是 ISM 的事 |
|---|---|---|---|
| **O'Neill 2024**（Cell Genomics）| bulk ONT | cancer cohort landscape（SV + phasing + 甲基地景）| **term-search subclone=0**，不建 subclone 樹；甲基=germline-aDMR |
| **TumorLens 2026**（medRxiv）| bulk ONT | 統一多組學**偵測** + HLA-ASM immune-escape | 偵測+annotation，**非 subclone 重建** |
| **Wakhan / Severus 2025** | bulk ONT | CNA/LOH/SV + phasing | 遺傳骨幹，**無 subclone 樹** |
| **PCDH 2026**（Hackett/Blundell）| **EM-seq+ONT** targeted-capture | PCDH 特製 barcode locus 甲基-clock lineage | **甲基當 DRIVER / 特製 locus / 無 normal-cis**（反 ISM 紅線）|
| **Lee/Liu-Goretsky 2025** | ONT | 23 melanoma subline 多組學演化 | **single-cell-DERIVED（物理分選）**，非單 bulk |
| **Foltz 2024**（SomaticHaplotype）| **10X linked-read（非 ONT）** | 成對 somatic mutation 共享 barcode → same-haplotype vs subclone | **非 ONT / 零甲基 / 無樹**（= ISM 共現骨幹的 genetic 半邊先驗）|
| **McAuley 2025** | ONT+Illumina | 單病人 concordance + DMR-vs-normal | **非 subclone 重建** |

**HCC1395 / COLO829 的 subclone 外部真值（agent a702cc1c 全文核；皆非 ONT、皆非單分子共現）**：
- **HCC1395 = Fang 2021 SEQC2**（Nat Biotechnol，PMC8532138）：SuperFreq + PhyloWGS + subHMM 於**短讀 WES/WGS VAF + 10x scCNV** → ~10 subclone（⚠ 數字二手待核 Fig 4）。🔴 **PacBio 50X 長讀有產，但只用於 variant-detection benchmark、未用於 subclone 重建** → 連 HCC1395 的長讀資料都沒人拿來建 subclone。chr2:18M pivot（18,086,020）落 SEQC2 high-conf 空隙 = truth-unevaluable。
- **COLO829 = Velazquez 2020**（Comm Biol，PMC7316972）：1475-cell **single-cell DNA CN**，4 subclone A-D（47.6/8.5/3.1/40.8%），**CN-level 非 SNV-level**。
- 其餘 4 株（H1437/H2009/HCC1937/HCC1954）無 line-specific subclone 解（只 SKY/SV/CNV）。🔴 HCC1937 常被誤搜成 HCC1395。
- → **ISM 是首個在 ONT 單樣本上用單分子共現骨幹做 regional subclone 重建。**

**🔴🔴 最空的格子 = 單分子共現骨幹（sub-question 3，最清楚無佔者）**：唯二先例 = **Foltz 2024**（linked-read 非 ONT、零甲基、descriptive count 無樹）+ **TreeClone 2017**（短讀 paired-end co-map 的 mutation-pair **統計 VAF 模型**，非單分子 read-level、非長讀）。**native ONT 單分子多 sSNV 共現作 subclone 骨幹 = 無佔者** → ISM 增量 = 原生 ONT + 甲基輔助軸 + 結構檢定（PERMANOVA/PERMDISP）+ normal-anchored cis。
> 2026-06 web sweep 補核：新名字皆 off-niche（LongSom=長讀單細胞 RNA / Canopy2=bulk+scRNA / TreeClone=短讀 / EVOFLUx=bulk fCpG clock）。**無新佔者**。

---

## §2 三個對比軸（✅ 用戶確認）

| 對比軸 | 過去研究 | ISM |
|---|---|---|
| **樣本數 / 解析度** | multi-region（TRACERx）/ single-cell（**確認黃金標準**）| **單樣本** bulk |
| **定序平台** | 短讀 **~100–150 bp**（150 常見）| **ONT 長讀**（kb–Mb）|
| **甲基取得** | 需**另做實驗**（bisulfite / EM-seq）| **native 5mC**（basecaller 直接從電流訊號出）= **一次資料一次分析** |

---

## §3 ONT 三優勢（🔴 精準點 2：單分子共現是命脈，別漏）

ISM 的方法成立**完全靠 ONT 長讀的三個優勢，三個都要列**：

1. **長程 phasing / haplotagging** —— 把 read 指派到 germline haplotype（鑑別 allelic vs clonal）。
2. **native 5mC** —— basecaller 直接出甲基，無需 bisulfite/EM-seq 額外實驗。
3. 🔴 **單分子多 somatic sSNV 共現** —— 長讀讓**多個 somatic sSNV 出現在同一條分子上 → read-level 直接共現觀測**。**這是重建骨幹（four-gamete / 共現 census）能成立的唯一前提，短讀做不到。** ← **本方法的命脈使能者，論文必強調，別只講 phasing+甲基。**

---

## §4 三個精準點（決定這段能不能扛 reviewer）

**🔴 精準點 1（最關鍵的誠實邊界）**：multi-region / single-cell 是 subclone 確認的**黃金標準**（也是 ISM 自己 ⭐3 天花板的原因，Tarabichi LEARN）。
- ❌ **不可**寫「我們比單細胞好 / 不需要單細胞」。
- ✅ **正確口徑**：「single-cell 是確認黃金標準但**昂貴、臨床不易普及**；我們證明 **ONT 單樣本能達到 characterization 級重建**，並**誠實標出哪些是單樣本結構性無法解決、需 single-cell 的部分**（如 subclone-specificity）。」
- → **優勢 = 實用性 / 成本 / 可近性，不是解析度。** 這樣定位才同時有貢獻又站得住。

**🔴 精準點 2（ONT 優勢別漏單分子共現）**：見 §3 —— 三優勢全列，單分子共現是骨幹使能者。

**✅ 精準點 3（小修，已正確）**：短讀 ~100–150 bp（150 常見）；「需另做實驗才知甲基」= bisulfite / EM-seq；ONT「一次資料一次分析」賣點成立。

---

## §5 ISM 創新口徑（sanctioned）+ 紅線

**創新（一句）**：**首個在 ONT 單樣本（tumor+normal）上、以單分子 somatic 共現為非循環骨幹做 regional subclone 重建（partition + nested/sibling），並用 haplotag 鑑別 allelic-vs-clonal、甲基 bounded-auxiliary characterize。**

**🔴 紅線（投稿必守）**：
1. **reconstruction = regional partition / nested-sibling，非 genome-wide clonal phylogeny / CCF tree**（內文第一次出現即界定；禁無限定用「reconstruction」）。
2. **甲基 = bounded-auxiliary characterize / corroborate，非 driver / 非 confirmer**（甲基用途已窮盡，負篩+L3 旗標可用、確認不可）。
3. **單樣本 = characterization / identification，非 confirmation**（confirmation 需 single-cell/multi-region；⭐3 封頂）。
4. **優勢 = 實用性/成本/可近性，非解析度**；禁「勝過/不需 single-cell」。
5. **co-occurrence → lineage 是 necessary-not-sufficient**（互斥→不同 subclone 可靠；共現→同譜系須加 VAF/CCF；Foltz NRAS 反例）。
6. **禁「對手用短讀」或「對手缺顯著性檢定」**（cvlr/ASMS/MethylBERT 都 ONT-capable + 有 randomization）。
7. **haplotag 非循環性**靠「phasing 用零甲基」（甲基若進 phasing→循環）；Methods 須明寫。

---

## §6 多區域整合成「一棵樹」有人做乾淨嗎？→ 沒有；報 ensemble 是 field-standard（2026-06-30 agent aa2c509e）

**裁決**：多區域是黃金標準，但**不產一棵定論樹** —— 全部工具輸出 **posterior / ensemble / 枚舉解集**（PhyloWGS MCMC posterior · Canopy 所有 config+confidence · PairTree posterior+consensus graph · Orchard posterior 抽樣 · CITUP/SPRUCE 枚舉所有相容樹 · LICHeE 評估所有最佳分樹 · BitPhylogeny full-Bayesian posterior）。非唯一性是**內在的**（multiplicity/CCF 由 DNA 不可識別，DeCiFer）；更多樣本**縮小**ambiguity 但不保證歸零（SPRUCE 逐字）。
- 🔴🔴 **killer citation = MACH2（2024，bioRxiv 2024.11.19.624301）**：領域特地發 MACHINA 後繼工具，**唯一目的=枚舉所有等可能解、取代回單一棵樹**，因「只回一個解漏掉資料同等支持的其他解」→ **領域自己宣告：報完整等可能解集才是正解、回單樹才是 bug**。
- 🔴 **MACHINA**（Nat Genet 2018）：「parsimony 不足以區分…需 biologically-motivated 額外準則」→ **背書 ISM「排序絕不用甲基」**（不可硬塞無依據 tie-breaker）。
- → **ISM 的「定不出來即答案 + 完整候選樹集 + 等機率 badge + 不用甲基排序」= field-endorsed 嚴謹輸出**（4 軸強 match）。🔴 **唯一邊界**：領域是用**多區域黃金標準**撞這牆；ISM 是**單樣本**（Tarabichi「勉強才該建樹」regime）→ 優勢在**輸出紀律**非**解析度平起平坐**；⭐3 不變；**禁**「領域也失敗故 ISM 等同」。
- **投稿口徑**：「報 non-identifiability / 等可能 ensemble 而非硬湊一棵樹，是領域背書標準做法（MACH2 2024 即為此而生）；ISM 誠實輸出形式與黃金標準一致，非缺陷。」

## §7 相鄰 phase set 用甲基/VAF 串聯成更長？→ germline-only prior art，不延長 somatic 天花板（agent afdcaf94，源碼核）

**裁決**：甲基串聯相鄰 phase set **有人做且是 prior art，但全 GERMLINE-only** —— 延長的是 germline-haplotype 連續性，**非 cell/subclone 連結**。
- **甲基-bridging（germline）**：MethPhaser 2024（PMID 38909018，N50 +78-151%）→ HapBridge 2025 → **LongHap 2026**（5mC+SNV 聯合，補 18.2% 純序列接不起的 block）。皆 germline，cancer=future work。
- **VAF**：read 層**不能**串（族群統計量無共享分子；MethPhaser「兩 block 間無連結資訊」）；sample 層 CCF/CN 是**統計縫合**（歸 subclone 族群）非分子連結。
- **somatic cross-block bridging = 白地未解**。
- 🔴 **crux**：甲基串起兩 germline block 給「同一條 germline haplotype」連續性，**≠ 哪個 cell/subclone**（germline-haplotype 連續 ≠ somatic-cell-lineage 連續）。
- **🔴🔴 對 ISM（載重）**：(1) 甲基-bridging 是 prior art → **不可宣稱首創**；白地精準框 somatic-subclone。(2) germline-only → 不延長 ISM somatic 天花板。(3) **non-circularity 保住**：ISM HP-tag 來自 LongPhase-S **零甲基路徑**（源碼 grep=0）→ PERMANOVA(甲基~HP) 非循環；**MethPhaser/HapBridge/LongHap 的存在正好證明「ISM 若用甲基建 HP-tag 就循環」→ phasing 零甲基是必要設計非疏忽**。

---

## §8 主流 subclone 重建方法 vs ISM 優缺點對照表（2026-07-01）

> 每方法皆已 74-源庫/競品全文核（verified）；外部 = L3 投稿前 /citation-verification。「確認級別」= confirmation(對 ground truth) vs characterization(對 null/label)。

| 方法/學派（代表）| 樣本·平台 | 重建骨幹 | 甲基角色 | 輸出 | 確認級別 | ✅ 優點 | ⚠ 缺點/限制 |
|---|---|---|---|---|---|---|---|
| **SNV-CCF clustering**(PyClone-VI/DPClust/SciClone/MOBSTER) | bulk 短讀 WGS/WES | VAF→CCF | 無 | CCF clonal/subclonal 群 | characterization | cohort-scale；成熟標準；beta-binomial 處理 over-dispersion | multiplicity 非識別→CCF elusive；只產群非樹；零 read-level/零甲基 |
| **Phylogeny tree**(PhyloWGS/Canopy/CITUP/PairTree) | bulk 短讀，**偏多樣本** | CCF + infinite-sites/pigeonhole | 無 | **posterior/ensemble over trees** | characterization | 建完整 lineage 樹；報不確定性(PairTree/MACH2) | 樹非唯一；單樣本嚴重欠定；零甲基/零 read-level |
| **Allele-CNA**(ASCAT/Battenberg/TITAN) | bulk 短讀 | BAF+logR→CN | 無 | allele-specific CN(+CN tree) | characterization | purity/ploidy/LOH；ISM 上游 context | purity-ploidy 不可識別；WGD 等價；低 fraction CNA 漏檢 |
| **多區域**(TRACERx+CONIPHER/MACHINA→MACH2) | **多空間 biopsy** 短讀 | 跨區 CCF+pigeonhole | 無 | **枚舉等可能解集**(MACH2) | **接近 confirmation**(空間復現) | 黃金標準之一；空間約束更緊 | 需多 biopsy(臨床難)；仍報 ensemble 非單樹；零甲基 |
| **single-cell DNA**(10x scCNV/SCITE) | **單細胞** | per-cell CNV/SNV | 無 | genome-wide lineage tree | **confirmation(黃金標準)** | 真 genome-wide lineage；終極 arbiter | 昂貴/臨床不普及；dropout/低覆蓋；零甲基 |
| **single-cell 甲基**(Epiclomal/MethylTree/Sgootr/Gaiti) | **單細胞** sc-WGBS | per-cell epimutation barcode | **DRIVER(真 lineage)** | genome-wide epiclone tree | **confirmation** | 甲基 lineage ~100%(REAL)；正交 SNV | 需單細胞解析度+label；昂貴；**不同 regime** |
| **read-level 甲基**(cvlr/ASMS/CpelNano/qFDRP) | bulk ONT | 甲基 read 分群 | 主訊號 | region ASM/cluster | characterization | ONT-capable+有 randomization 檢定；read-level | **無 normal-cis→ASM≠subclone**；target=germline ASM；無 subclone target |
| **supervised read-origin**(ROCIT/MethylBERT) | bulk(ROCIT=**PacBio** demo) | 甲基(supervised) | DRIVER | tumour/非tumour 二分 | — | read-level；改善 variant calling | supervised/binary/無 normal；非 subclone；非無監督 |
| **long-read 多組學偵測**(O'Neill/TumorLens/Wakhan) | bulk ONT | SNV/SV/CNA(+甲基地景) | aDMR/window DMR | 偵測+annotation | — | 統一多組學一次資料；read-level phased aDMR | **不建 subclone 樹**(O'Neill term-search=0)；甲基=germline-aDMR/immune-escape |
| **linked-read 共現**(Foltz SomaticHaplotype) | bulk **10X linked-read** | **成對 somatic 共享 barcode** | 無 | same-hap vs separate-subclone | characterization | **共現邏輯先驗**(= ISM genetic 半邊) | 非原生 ONT；零甲基；descriptive 無結構檢定；無樹 |
| **🟣 ISM**(本研究) | **單樣本 bulk ONT**(tumor+normal) | **單分子 somatic 共現**(four-gamete)+haplotag | **bounded-auxiliary corroborate** | regional partition+nested/sibling+**等可能候選集** | characterization(⭐3) | **見下方專欄** | **見下方專欄** |

### 🟣 ISM 優缺點（專欄）

**✅ 優點**
1. **單樣本可行 + 成本/可近性高**（不需多區域 biopsy 或單細胞；臨床友善）。
2. **native 5mC 一次資料一次分析**（無需 bisulfite/EM-seq 額外實驗）。
3. **單分子多 sSNV 共現 = 直接局部譜系證據**（short-read 做不到；= 重建骨幹命脈）。
4. **normal-anchored somatic cis-test**（germline 工具皆無此錨點）。
5. **無監督 read×read 距離矩陣 PERMANOVA 結構檢定**（非循環，HP 標籤獨立於甲基）。
6. **誠實 ensemble 輸出（定不出來即答案）= field-endorsed**（MACH2 2026 同向）。
7. **佔尚空交集**：unsupervised 結構檢定 + normal-cis + somatic-subclone target on bulk ONT = 無人佔。

**⚠ 缺點/限制**
1. **單樣本 = characterization 非 confirmation**（⭐3；confirmation 需 single-cell/multi-region，黃金標準 ISM 無）。
2. **跨 phase-set 無連結 → regional partition 非 genome-wide tree**（物理天花板）。
3. **定不出單一樹**（候選等可能；誠實但非「一個乾淨答案」）。
4. **甲基 bounded-auxiliary 非 driver**；**G-B（subclone 甲基 somatic-specific）UNDETERMINED**（624 needs_methyl 乾淨可用≈0）。
5. **single-pipeline**（跨樣本一致=內部一致非獨立平台；需正交真值）。
6. **繼承 multiplicity 非識別**（VAF→CCF 不算，刻意繞但也不解）。

### 🎯 ISM 在各軸「贏 / 平 / 讓」
- **贏（實用性軸）**：單樣本 + 一次資料 + read-level 直接共現 → 比 single-cell/多區域**成本低可近**；比短讀 SR**多 read-level 分子證據 + 原生甲基**。
- **平（誠實輸出軸）**：報等可能 ensemble 與 MACH2/PairTree **同為 field-standard**。
- **讓（解析度/確認軸）**：genome-wide lineage 樹**讓給 single-cell/多區域**；正式 CCF-tree **讓給 SNV-CCF/phylogeny 學派**。**禁宣稱在這些軸領先**。

---

> **🔴 §6/§7 投稿前 /citation-verification 補核（L3）**：MACH2 bioRxiv 2024.11.19.624301（最載重，「枚舉等解」逐字）· SPRUCE「ambiguity decreases with samples」· PairTree/Orchard/CITUP/LICHeE 輸出型態 · HapBridge perf 數字（abstract，PDF 403）· LongHap bioRxiv 2026.03.11.710820 · CONIPHER 是否內文報單樹 vs 建模不確定性（勿誤述）。MethPhaser PMID 38909018 + Foltz PMID 39149342 + LongPhase-S 零甲基 = L1 源碼/全文核可直引。

*關聯 memory：`project_subclone_confirmation_construction_ont_transferability`（白地+移植）· `project_subclone_snv_linkage_verification_pipeline`（共現骨幹）· `project_methylation_use_exhausted_bounded_auxiliary`（甲基窮盡）· `project_candidate_tree_ranking_impossible`（定不出來即答案）· `project_clone_subclone_landscape_and_ism_feasibility`（競品）。詞典紅線 SoT：`20260621_baseline_discipline_conceptual_definitions_01.md`。*
