<!--
建立時間: 2026-06-21
報告類型: method_comparison — clone/subclone 重建全景知識紀錄 + ISM 可行性裁決（durable reference）
任務類型: D handoff — 跨 session/AI 可參考、對比、驗證、確認的外部知識留底
build_branch: research/subclonal-reconstruction-202606 (HEAD c660179)
status: synthesis（landscape + verdict）；外部 citation = L3 agent-fetched，投稿前過 /citation-verification；ISM 內部數字標 tier+來源
data_sources:
  - workflow ww3nitd8n（5-agent：3 survey 先讀 69-源庫再補 2024-2026 缺口 + 合成 + 對抗稽核 PASS）
  - /big7_disk/liaoyoyo2001/external_validation/_landscape/{05_context_cards_index,08_paper_literature_map}.md（69 源庫）
  - InterSubMod/docs/method_comparison/20260619_subclone_analysis_interpretation_full_framework_01.md（why-hard GAP-A~L）
  - InterSubMod/docs/method_comparison/20260620_somatic_locus_methylation_combination_enumeration_01.md（組合窮舉）
  - InterSubMod/docs/methodology/20260617_structure_label_situation_inventory_01.md（內部數字 SoT）
  - memories: project_{subclonal_reconstruction_paper_focus, thesis_writing_architecture, ism_complete_tpfpfn_existence_cis, cross_sample_asm_reproducibility, tumor_only_axis_negative_subclone_classification, chr2_18m_subclone_locus_verification, apriori_subclone_classification_model}
provenance_note: ISM 內部 metric 全溯源上列 memories/docs（verified）；外部 claim 為 agent WebFetch（L3，各帶 PMID）；2025-26 競品（ROCIT/TumorLens/MethPhaser-Cancer/DPClust-on-LR）為 L3-preprint-abstract 未全文核。
adversarial_audit: workflow ww3nitd8n evaluator = PASS（0 fabrication / 0 overclaim / 0 redline；honesty check 確認裁決清楚分層未閃避；4 non-blocking 修正已折入）。
-->

# Clone/Subclone 重建：全景知識紀錄 + ISM 可行性裁決

> **這份是什麼**：把「腫瘤 clone/subclone 重建」的**全流派 / 方法 / 論文**、**難點與共識**、**未知與分歧**、**可學習方向**、**盲點檢核**整理成一份**可長期參考的知識紀錄**，並回答「ISM 能否完成論文目標」的**分級裁決**。設計給後續 session / 其他 AI **對比、驗證、確認**用。
> **可信度標註**：外部 = **L3**（agent WebFetch，帶 PMID，投稿前過 `/citation-verification`）；ISM 內部 = tier ⭐ + L1-L5 + 來源檔。
> **配套教學指南**（給其他教學 AI）：`InterSubMod/docs/method_comparison/20260621_clone_subclone_teaching_guide_for_ai_01.md`。

---

## §0 一頁總結（重點解釋敘述 — Verdict-Pyramid）

**核心裁決（先講結論）**：論文目標 **「Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing」可完成 —— 但僅在嚴格 frame 下**，且**條件性受 HD-1 未解所限**：

1. **可完成（as framed）**：作為 **proof-of-concept / case study** 的「**區域內（phase-block 尺度）、read 層次、多證據（haplotag + multi-sSNV + methylation）subclone characterization / verification**」，達 **⭐3 單樣本**（chr2:18M 旗艦已 PRIMARY-verified），G-A 解後 6 樣本 → **⭐4**。**重建骨幹 = somatic haplotagging（genetic：longphase-S + LOH + CN）**；**methylation = characterization + honest-negative + corroborate（與多點 somatic SNV 共用，非獨立 driver）**。
2. **不可宣稱（overclaim 紅線）**：methylation 驅動的 genome-wide subclonal reconstruction / 完整 clonal phylogeny-CCF tree / tumor-only 無監督 subclone 偵測 / methylation 當獨立 discriminator 或 filter —— 現有證據**明確不支持**。標題的 "reconstruction" **必須在內文界定為 regional LOH-constrained same-haplotype 子家族 partition**，否則撞學派 1-2 的正式定義被 reviewer reject。
3. **「達到別人做不到」的誠實措辭**：ISM **不佔據任何沒人有的能力**；它佔據**先驗元件的新穎組合於一個尚空的交集** = 「first **INTEGRATION** of（unsupervised read×read 距離矩陣 PERMANOVA 結構檢定 + normal-baseline somatic cis-test + somatic-subclone target）at tumor/normal read level on **bulk ONT**」。single-cell 在 lineage 甲基**更強**、SNV-CCF 在 genome-wide 重建**更完整**；ISM 的位置 = **regional + read-level + multi-evidence + bulk-long-read + 甲基當佐證**。
4. **決定強度的 gate（解前不可 overclaim）**：**HD-1**（phasing spine by-construction 循環，未解前骨幹宣稱受限）+ **G-B**（subclone 甲基 somatic-specific vs germline-allelic = UNDETERMINED，受 LOH-unmask confound 壓制）+ **G-A**（COLO829 缺 + zhenyu112 SPOF → 封頂 ⭐3）+ **V-1**（乾淨 cis 1/816 under-tested）+ 大規模驗證 = future work。

> **一句話**：論文做得成，但它是「**在 genetic 骨幹上、用 read 層次多證據、對特定區域的 subclone 做有 documented confounds 的個案刻畫與驗證**」，不是「**用甲基重建全基因組 subclone 樹、勝過既有方法**」。把標題的每個詞放回這個界定，論文就站得住。

---

## §1 難點（why-hard — 為何 subclone 重建本質困難）

> 每條附一手 problem-statement 來源；對映 why-hard 缺口圖 GAP-A~L（doc `20260619_subclone_analysis_interpretation_full_framework_01.md`）。

| # | 難點 | 一手陳述 | 對 ISM 的意涵 |
|---|---|---|---|
| **H1** | **SNV-only 不可識別性（三層）** | tree / cluster / evolution-model 非唯一（少樣本多 phylogeny 同等相容）；Tarabichi 2021(PMID 33398189)、Balaparya&De 2018 | ISM 不走 VAF→CCF，故不直接撞；但也**不能宣稱解決** |
| **H2** | **multiplicity 歧義** | VAF=purity×prevalence×multiplicity 複合量，**即使 allele-specific CN 已知仍不足定 m**；DeCiFer Satas 2021(PMID 34416171) | 「半帶 ALT」formally non-identifiable（multiplicity=1 / subclone CCF<1 / ADO）；ISM 無 multiplicity model 不能升 CCF |
| **H3** | **加 LOH+CNV 仍難** | purity-ploidy 不可識別 / WGD 等價（ASCAT Van Loo 2010 PMID 20837533）；subclonal-CNA 低 fraction 漏檢（TITAN Ha 2014 PMID 25060187）；SNV-tree↔CNA-tree 不一致；cohort 僅 partial（PCAWG-11 Dentro 2021） | CN 是 context 非 driver；**加 CN 不消除困難** |
| **H4** | **低-AF 雙重 floor** | read-sampling（NRPCC）+ caller floor（VAF<2% FP 高不論 coverage）；Petrackova 2019 | 低 VAF carrier reads<3 → subclone 軸測不了（ISM carrier-limited 414） |
| **H5** | **over-clustering（無結構被切成假群）** | consensus/最佳-k 把無結構資料切成「看似穩定」群（Şenbabaoğlu 2014 PMID 25158761）；M3C 偏高 K 不測 K=1（PMID 32020004）；epiallele read-內相關（DNMT1 ~36bp PMID 36285438）；單 clone 高 PDR（Landau 2014 PMID 25490447） | ISM tumor-only 無監督 NEGATIVE 的外部根據；「切 3 群」需先排 over-cluster |
| **H6** | **bulk 觀測單元 = READ（非 CELL）** | bulk long-read per-read 局部單分子、**無跨 locus 細胞 linkage、無跨 phase-set linkage**（MethPhaser PMID 38909018「兩 block 間無連結資訊」） | ISM 無法在單 bulk 樣本繼承 single-cell genome-wide lineage → **物理上限於 regional** |
| **H7** | **LOH-unmask 既存 ASM（甲基特有 confound）** | tumor imprinted-DMR 甲基變化 **82–91% 由 CN/LOH 解釋非 epimutation**（Martin-Trujillo 2017 PMID 28883545）；甲基隨 aUPD scale r=0.76（Chase 2015 PMID 26114957） | LOH 區的「等位甲基」最可能是既存 ASM 顯形，**非新 somatic**；G-B 未定的決定性壓制因素 |

**漏斗式總結**：純 DNA（SNV ± LOH/CNV）subclone 重建**系統性欠定**（H1-H4）→ 動機找 orthogonal 訊號（甲基）；但 bulk 甲基受 H5-H7 三重夾擊 → 甲基的 subclone 力在 bulk read-level **未建立**。

---

## §2 共識（領域已定論 — 可直接引用不需重證）

| # | 共識 | 來源 | 強度 |
|---|---|---|---|
| C1 | 判 somatic **必須** matched-normal | 全領域（caller/SR 骨幹）| 鐵共識 |
| C2 | **single-cell 甲基 = lineage 訊號 REAL**（近-100%，有 label）| Gaiti 2019(PMID 31092926)/Epiclomal(PMID 32966276)/MethylTree(PMID 39820752)/EPI-Clone(PMID 40399669) | 強 |
| C3 | **bulk read-level 甲基 subclone = NOT-established** | 領域空白 + ISM G-B undetermined；最新長讀 subclone 靠物理預分選（Liu-Goretsky 2025 PMID 40950124）| 強（空白本身是共識）|
| C4 | **正式 reconstruction = clone# + CCF + lineage/CN tree + purity deconv** | 學派 1-3 定義（PyClone/PhyloWGS/ASCAT）| 定義性 |
| C5 | 純 DNA subclone 重建**系統性欠定** | Tarabichi 2021 / DeCiFer 2021 / Şenbabaoğlu 2014 | 強 |
| C6 | **distance-based 甲基 reconstruction 本身非新穎** | Sgootr 2023（真腫瘤單細胞 distance-based methylation lineage tree）| 強（framing 限制）|
| C7 | 🔴 **cvlr/ASMS/CpelNano/MethylBERT 全 ONT-capable 且全有 randomization/permutation 檢定** | 各源碼/論文親讀 | 鐵共識（**innovation 口徑紅線**）|
| C8 | methylation-as-variant-filter = **DEAD**（ISM 內部 4 路證偽）| ISM；對齊外部 ASM≠subclone | 內部定論 |

---

## §3 未知與分歧（哪裡還沒答案 / 哪裡有爭議 / 我們曾誤判被更正）

| # | 未知/分歧 | 狀態 | 如何收斂 |
|---|---|---|---|
| U1 | **G-B：subclone-discriminating 甲基 = somatic-specific 還是 germline-allelic？** | **UNDETERMINED**（最關鍵開放問題）| within-haplotype somatic-vs-baseline control；受 H7 LOH-unmask 壓制 |
| U2 | **bulk 區域 read-level 多證據能否在單樣本「確認」subclone？** | 結構性存疑（非僅資料量問題）| 組合窮舉證多數 undecidable；需 normal cis-control + multi-sSNV CCF + single-cell 終判 |
| U3 | **乾淨 somatic-cis 真實稀有度** | 僅測 HCC1395 **1/816** loci（under-tested）| full-locus × 6-sample 才能定量 |
| U4 | 🔴 **我們曾誤判被更正：ASMS「缺顯著性檢定」口徑錯** | 已更正 — **ASMS 實有 5000-perm 檢定** | 對外**禁**用「對手缺檢定」；C7 紅線 |
| U5 | **2025-26 競品正快速填補白地** | L3-preprint，未全文核 | ROCIT/TumorLens/MethPhaser-Cancer/DPClust-on-LR 須全文核 + head-to-head（見 §11）|
| U6 | **HD-1：phasing spine by-construction 循環** | 未解 → 骨幹宣稱受限 | R-SELFREF（~25-50hr C++）或降 characterization |
| U7 | chr2:18M 18,086,020 timing | hypothesis（不在 somatic VCF；ASM 無法形式排除）| 需 cis-control |

> **分歧的本質**：與 69 源 **0 真 CONFLICT**，但有 **framing tensions**（Sgootr/Gaiti/BitPhylogeny 限制「reconstruction 用詞」與「距離式甲基新穎」）→ 強制窄 framing，非數據衝突。

---

## §4 可學習方向（有明確下一步參考的外部方法）

| 借鏡來源 | 可學什麼 | 對 ISM 適配度 |
|---|---|---|
| **PyClone-VI**（PMID PMC7730797）| over-complete mixture **auto-k** + **beta-binomial over-dispersion** | 高 — ISM per-CpG Fisher 仍 binomial-family（over-dispersion = KB 必修）；dendrogram k 仍手動 |
| **M3C**（PMID 32020004）| **correlation-preserving null**（取代 bare best-silhouette→PERMANOVA）| 高 — 直擊 ISM over-cluster 根因；但需 binary/Bernoulli 版本 |
| **Barrett 2017 Bayesian epiallele**（PMID 28743252）| **model-based epiallele 計數**（「幾群」正解，含 noise term + 複雜度懲罰）| 中-高 — 正解方向但需 read×CpG，ONT 單 locus 大改 |
| **DeCiFer**（PMID 34416171）| DCF / multiplicity model（把 fractional-ALT 升 CCF）| 中 — 可能超出 ISM regional-characterization 定位（紅線 f）|
| **Rosenski 2025 atlas**（PMID 40069157）| **460 parental-ASM + 34,426 SD-ASM 座標 intersect blocklist** | 高 — 標既存 ASM（分 I2 LOH-unmask vs somatic）|
| **SciClone**（PMID 25102416）| **CN-neutral 紀律** | 高 — = ISM「HP-axis 固定 CN」鏡像，正當化設計 |
| **Canopy**（PMID 27573852）| **多 config + confidence** 誠實-不確定性呈現 | 中 — 報告多可能解非單一，符合 ISM 誠實框架 |
| **single-cell ground truth**（MethylTree/EPI-Clone）| I4 subclone 終判 | 高 — bulk 不可獨立確證的唯一出路 |

---

## §5 Blind-spot 檢核表（容易沒考慮到 — 每次分析前對照）

- [ ] **LOH-unmask 既存 ASM**（82–91%）：LOH 區「等位甲基」先當既存 ASM 顯形，非新 somatic
- [ ] **epiallele over-cluster**：read 內相關把 1 clone 切成 k=3-4（切群前先過 correlation-preserving null）
- [ ] **double-dip**：無監督「挑最分得開再測」= 循環（noise≈structure 83-100%）
- [ ] **homopolymer 假群**：chr2 pos4 20bp poly-T 製造假 read 群（看 reads.tsv 側欄）
- [ ] **SampleASM 飽和假象**：96.5% 命中 → **勿引用為「幾乎都有甲基訊號」**
- [ ] **germline vs residual 意義相反**：都用 normal 但結論相反（germline=非 somatic / residual=tumor-specific），禁合併
- [ ] **cross-phase-set 無連結**：per-block call **不能**串成 genome-wide lineage（需外部 anchor）
- [ ] **confirmed ≥2 群 ≠ subclone**：85% 對齊 germline/carrier = cis-ASM
- [ ] **half-ALT non-identifiable**：半帶不能單獨宣稱 subclone（multiplicity）
- [ ] 🔴 **reviewer-reject 雷**：禁說「對手用 short-read」或「對手缺顯著性檢定」（C7）

---

## §6 後續對比驗證指引（讓未來 session 能 compare-verify-confirm）

> 每條 claim → 證據層 L1-L5 + 來源 + **如何重新驗證** + 與 ISM 現況的對照。

| Claim | 層級 | 來源 | 如何重驗 |
|---|---|---|---|
| ASM 存在非 filter（TP 3.95%/FP 1.07%/sens~4%）| L1-L2 ⭐2-4 | memory `ism_complete_tpfpfn_existence_cis` | grep 原始 TP/FP TSV 重算 |
| 跨樣本 6/6 excess-over-null | L2 ⭐3 | memory `cross_sample_asm_reproducibility` | 6 樣本 excess 重跑（非 raw rate）|
| tumor-only 無監督 NEGATIVE（noise≈structure 83-100%）| Python=C++ verified | memory `tumor_only_axis_negative` | gap-null Python 重算對 C++ |
| chr2:18M 骨幹（0 violations/LOH 1.4%）| L1 ⭐3 PRIMARY | memory `chr2_18m_subclone_locus_verification` | `research/flagship_chr2_18086020_*/` ism_out 重讀 |
| a-priori 軸 non-circular ~7.6% | ⭐3 | memory `apriori_subclone_classification_model` | legacy-VC/gap-null 重算 |
| LOH-unmask 82-91% | L3 | Martin-Trujillo PMID 28883545 | PMC 全文核（過 /citation-verification）|
| 2025-26 競品 ROCIT/TumorLens/… | **L3-preprint 未全文核** | 見 §11 | **全文核 + head-to-head**（投稿前必做）|

---

## §7 全景 taxonomy — 6 學派

> 前 4 學派 = **遺傳骨幹**（aggregate VAF/CN，零甲基、非 read-level）= ISM 紅線「reconstruction 由 genetic 驅動」所指，與 ISM **互補（0 真 CONFLICT）**。學派 5 = methylation-lineage 唯一真實 regime（single-cell）。學派 6 = ISM **最近鄰與真實邊界**。

### 學派 1 — SNV / VAF-CCF clustering
- **核心**：SNV 的 VAF 經 CN/purity 校正轉 CCF，按 CCF clustering 成 clonal/subclonal 群。
- **代表**：PyClone/PyClone-VI（CCF grid=100，default binomial）· DPClust（PCAWG 標準）· SciClone（CN-neutral 過濾）· MOBSTER（Beta+Pareto 中性尾）。
- **論文**：Roth 2014 PMID 24633410 · Gillis&Roth 2020 PMC7730797 · Miller 2014 PMID 25102416 · Caravagna 2020 MOBSTER（待核）。
- **regime**：bulk WGS/WES；genome-wide scope 但 per-mutation；零甲基、非 read-level。
- **CANNOT**：multiplicity 不可識別；noise-model = over/under-cluster 根因；只產 cluster 非 tree。
- **vs ISM**：互補骨幹；SciClone CN-neutral = ISM HP-axis 固定 CN 鏡像；PyClone-VI auto-k+beta-binomial 可借。

### 學派 2 — Phylogeny / clonal-tree
- **核心**：建祖裔 lineage **TREE**（infinite-sites/pigeonhole 約束）。
- **代表**：PhyloWGS（TSSB+MCMC）· Canopy（SNA-in-CNA 統計定相）· CITUP · LICHeE · B-SCITE。
- **論文**：Deshwar 2015 PMID 25786235 · Jiang 2016 PMID 27573852 · Malikic 2015 CITUP PMID 25568283 ·（LICHeE/B-SCITE 待核）。
- **regime**：bulk，偏好 multi-sample；輸出 lineage tree；零甲基、非 read-level（Canopy「phase」是 site-level 統計後驗非物理 haplotagging）。
- **CANNOT**：tree 非唯一（少樣本尤甚）；單 bulk 嚴重欠約束。
- **vs ISM**：**界定 scope 護欄** — 本學派定義「reconstruction」，ISM **絕不可宣稱 genome-wide phylogeny**；ISM 餵骨幹不建樹。

### 學派 3 — Allele-specific CNA（+CNA-phylogeny）
- **核心**：BAF+logR 聯估 purity/ploidy + per-segment 整數 CN，按 cellular fraction 分 clonal/subclonal CN。
- **代表**：ASCAT · Battenberg（PCAWG，餵 DPClust）· TITAN · FACETS · Sequenza · MEDICC2（CN phylogeny）。
- **論文**：Van Loo 2010 PMID 20837533 · Ha 2014 TITAN PMID 25060187 ·（餘待核）。
- **regime**：bulk；操作於連續 CN segment（最接近「regional」但 region=CN 段非 haplotype read-block）；零甲基。
- **CANNOT**：purity-ploidy 不可識別；WGD 等價；低 fraction subclonal-CNA 漏檢；CN-known 仍不足定 multiplicity。
- **vs ISM**：互補上游 context（ISM 消費 LOH BED+HP-ratio CN）。🔴 **Martin-Trujillo 82-91% imprinted-DMR by CN/LOH = G-B confound 根源**。

### 學派 4 — Multiplicity / timing
- **核心**：正面處理 multiplicity + 把 mutation 置於相對演化時間。
- **代表**：DeCiFer（DCF 取代 CCF，允許 mutation loss）· MutationTimeR（定年）· Dentro/PCAWG-11。
- **論文**：Satas 2021 PMID 34416171 · Gerstung 2020（PCAWG timing）· Dentro 2021。
- **regime**：bulk cohort；零甲基（timing 靠 mutation 累積非 epigenetic clock）。
- **CANNOT**：CN-known 仍不唯一定 m（DeCiFer 核心：CCF from DNA alone is elusive）；單 slice 仍只部分重建。
- **vs ISM**：**最乾淨 why-hard 錨**（GAP-A/E/F）= ISM Ch1-Ch2 漏斗動機；正當化找 orthogonal 訊號，但誠實 caveat：該訊號只 single-cell 真實。

### 學派 5 — Single-cell（DNA + methylation lineage）
- **核心**：以**細胞**為觀測單元 — genome-wide epimutation barcode，各 locus 在每 cell 內**物理連鎖**（bulk long-read 所缺）。
- **代表**：(甲基) Gaiti · Epiclomal · MethylTree（~100% supervised）· EPI-Clone · Sgootr（真腫瘤單細胞 distance-based）；(DNA) 10x scCNV；(物理預分選) Liu-Goretsky 2025 / Lee 2025。
- **論文**：Gaiti 2019 PMID 31092926 · Epiclomal PMID 32966276 · MethylTree PMID 39820752 · EPI-Clone PMID 40399669 · Sgootr 2023 · Liu-Goretsky PMID 40950124。
- **CAN**：甲基近-100% 重建 lineage/epiclone；Lee 2025 出完整多組學 subclone phylogeny。
- **CANNOT**：需 single-cell 解析度+label；Liu-Goretsky/Lee 靠**物理預分選**取得 linkage → 不解 ISM 的 bulk 單樣本問題。
- **vs ISM**：**framing tension（不同 regime 非競爭）**。Sgootr 證距離式甲基重建**非新穎**；single-cell 在 lineage 甲基**遠勝**。**單元根本差異** = single-cell 觀測 CELL（跨 locus linkage 在）vs bulk long-read 觀測 READ（無跨 locus/phase-set linkage）。

### 學派 6 — Long-read phasing & read-level methylation clustering（ISM 最近鄰）
- **核心**：(6a) 從 long read 解 haplotype、刻畫 per-haplotype 甲基/ASM 於 germline-haplotype 層；(6b) bulk long-read 不靠 phasing 按甲基把 READ 分群。
- **代表**：(6a) longphase-S somatic_haplotagging（ISM 骨幹，HP3 F1~0.746 由 somatic-allele 非甲基）· MethPhaser+MethPhaser-Cancer · NanoMethPhase · Wakhan · Severus；(6b) cvlr（EM Bernoulli+1000-random）· ASMS（5000-perm）· qFDRP（NHD 距離）· CpelNano（Ising/NME 熵+perm，ONT-native）· MethylBERT（supervised 反卷積）。**2025-26 新**（全文核見 §12）：TumorLens（統一 bulk-ONT 多組學偵測，ASM 限 74 HLA/APM 為 immune-escape，**非 subclone**）· ROCIT（read-level 甲基 tumour-origin 分類，**supervised/binary/無 normal**；PacBio 實證 ONT aspirational）。🔴 **更正**：O'Neill 2024 **不做 DPClust/subclone**（term-search=0），其甲基是 read-level haplotype-phased aDMR（germline-allele 軸，**非**「分離 DMR track」）；「DPClust + 平行 DMR 的 field-default pipeline」**作為單一已發表論文並不存在**（McAuley 2025 rectal 是單病人 concordance+DMR-vs-normal，才是「分離 DMR track」的例）。
- **論文**：LongPhase 2022 · LongPhase-S bioRxiv 2025 · MethPhaser PMID 38909018 · cvlr PMID 36726731 · ASMS bioRxiv 2024 · CpelNano · MethylBERT 2025 · TumorLens/ROCIT/O'Neill（**L3-preprint，投稿前全文核**）。
- **CANNOT**：6a 多數不做 subclone（bulk 無跨 phase-set linkage）；6b **不分離 ASM 與 somatic subclone**（無 matched-normal 錨點）；ROCIT binary/supervised/無 normal；DPClust-on-LR 甲基是分離 track 非 read-level 整合。
- **vs ISM**：**ISM 真實邊界 + 唯一可防守 delta**。🔴 **守則**：cvlr/ASMS/CpelNano/MethylBERT 全 ONT-capable 且全有 randomization/permutation 檢定 → **禁**說「對手 short-read」或「缺檢定」。

### 🎯 Where ISM Sits
ISM 佔據**先驗元件的新穎組合於尚空交集**：
> **「first INTEGRATION of（unsupervised read×read 距離矩陣 STRUCTURE 檢定 [PERMANOVA] + normal-baseline somatic CIS-test + somatic-subclone TARGET）at tumor/normal read level on bulk ONT」** — 每單元件皆 prior art，**整合**是 claim。
- vs 學派 1-4：ISM 不做 VAF→CCF/tree/CN-deconv/multiplicity，**消費**它們當 context；reconstruction 紅線歸 genetic。
- vs 學派 5：single-cell lineage 甲基**更好**；單元差異使 ISM 限 regional+bulk+characterization。
- vs 學派 6：三軸同時的交集目前無單一方法佔滿，但 2025-26 正快速填補 → **窄誠實 framing 強制**。

---

## §8 ISM 進度 × 全景對映

### A. 已 DONE（落學派 6；骨幹消費學派 1/3）
| 發現 | 對映 | tier | 關鍵數字 |
|---|---|---|---|
| ASM 存在非 filter | 6b | ⭐2-4 L1-L2 | TP 3.95% / FP 1.07% / sens ~4% / COLO TP≈FP |
| 跨樣本 reproducibility | 6b+4 | ⭐3 封頂 | 0/38 private；**6/6 excess-over-null>0** |
| BRCA2/ZAR1L ASM 真實非判別 | 6b+3 | ⭐3 | B-discrim NEGATIVE+anti（strong-ASM FP enriched 5×）；Δβ −0.122 |
| Situation inventory | 6+5 | L1-L2 | 30,490 TP；結構 81.2%；tumor-only 38.6%；clear≥2 9.6%，其中 85% 對齊 germline/carrier |
| a-priori 軸 C′ | 6 | ⭐3 verified | non-circular ~7.6%（≈unsup-clean ~2×）；within_clean≠subclone Jaccard 0.123/55% epiallele |
| chr2:18M 旗艦 | 6 旗艦 | ⭐3 PRIMARY | 0 violations；LOH 1.4%；「5 群」→~3 lineage（pos4 20bp poly-T）；18086020 不在 VCF |
| all-7 split accounting | 6+5 | verified | Jaccard 0.091-0.161 全低（3 癌種）|
| 5-state classification | 6 | verified | state-3 label-first ≠ tumor-only NEGATIVE；~52% 真關聯但 cis-ASM 非 subclone |
| combination enumeration | 1-4 落地 | PASS | subclone-like 組合單樣本全 undecidable；LOH-unmask 82-91% |
| methyl-phasing-assist V1-V12 | 6a | verified | T1 0.885(~6% unphase)；T2 OVERSTATED；T3 DEAD；甲基=germline-haplotype 層 |

### B. 已 NEGATIVE / DEAD（不可再開）
- **tumor-only 無監督 clustering+PERMANOVA = NEGATIVE**（double-dip；noise≈structure 83-100%）→ 不可宣稱 tumor-only 無監督 subclone 偵測。
- **methylation-as-filter = DEAD**（4 路）。
- **cluster-count 無乾淨無監督解**；「clear≥2」只靠 a-priori 對齊確認；確認 ≥2 ≠ subclone。

### C. OPEN GATES
HD-1（循環，未解骨幹受限）· G-B（甲基 somatic vs germline UNDETERMINED）· G-A（COLO829 缺 + SPOF，封頂 ⭐3）· V-1（cis 1/816 under-tested）· 大規模驗證=future。

---

## §9 可行性裁決（分級）

**OVERALL**：可行 —— **僅在嚴格 frame 下，且 HD-1 未解前 reconstruction-spine 宣稱受限**。論文目標可作 **proof-of-concept / case study** 完成，⭐3 single-sample（G-A 解後 ⭐4）。reconstruction 骨幹歸 somatic haplotagging（genetic），methylation 限 characterization+corroborate。任何把 methylation 升 genome-wide reconstruction 驅動 / 宣稱完整 phylogeny-CCF tree / tumor-only 無監督偵測 = OVERCLAIM。標題 "reconstruction" 須內文界定為 regional LOH-constrained same-haplotype partition。

**✅ 可完成（as framed）**：
1. Regional/read-level/multi-evidence subclone **characterization/verification**（chr2:18M 旗艦 PRIMARY-verified）
2. ASM 存在+somatic-enrichment 刻畫 + **誠實 non-filter**（honest-negative 本身可發表）
3. 跨 6 樣本/3 癌種 reproducibility ⭐3（G-A 後 ⭐4）
4. Sanctioned **整合 innovation**（三件式）
5. methyl-phasing-assist T1（0.885）supporting
6. a-priori 軸 C′ = 唯一合法 non-circular 管道展示（優勢 = **可解釋非偵測力**）

**❌ 不可完成（overclaim）**：methylation-driven genome-wide reconstruction surpassing existing · 完整 phylogeny/CCF tree · tumor-only 無監督偵測 · methylation 獨立 discriminator/filter · 「對手 short-read/缺檢定」· 「distance-based 甲基重建新穎」（Sgootr 已建）· 18086020 timing 當定論。

**Gate 條件**：HD-1 / G-B / G-A / V-1 / 大規模驗證=future（見 §8C）。

**「ISM 能做別人不能」的誠實答案**：不是無人有的能力，是**先驗元件於尚空交集的新穎組合 + 資料-regime niche**（regional+read-level+multi-evidence+bulk-long-read+甲基當佐證）。**別人更強處**：single-cell 甲基做 genome-wide lineage（觀測 CELL，真重建）；SNV-CCF+phylogeny+allele-CNA 做正式 reconstruction；DeCiFer/MutationTimeR 解 multiplicity/timing；Lee 2025 出完整長讀多組學 phylogeny（靠物理分選）；ROCIT 已做 read-level 甲基 tumour-origin genome-wide（supervised/無 normal）。ISM 在這些軸**不領先、不宜宣稱領先**。

---

## §10 ISM vs 最近鄰對照表

> 🔴 差異**不是**「有無檢定」或「short vs long read」，而是**檢定對象 + normal 錨點 + target**。PERMANOVA perm 次數待 C++ 源碼核。

| 方法 | regime | read-level? | methyl? | normal-anchored cis? | unsupervised structure test? | subclone target? | genome-wide recon? | validation |
|---|---|---|---|---|---|---|---|---|
| **ISM** | **bulk ONT** | **✅ read×read 距離矩陣** | ✅ corroborate（germline-hap 層，非 driver）| **✅ normal-baseline cis-test** | **✅ PERMANOVA+PERMDISP existence-gate** | **✅ a-priori-conditioned ~7.6%** | ❌ regional（非 tree）| ⭐3 single（HCC1395）；6-sample→⭐4 待 G-A；HD-1/G-B/V-1 open |
| cvlr | bulk ONT | ✅ EM Bernoulli | ✅ | ❌ | ◐ 甲基**差量** vs 1000-random（非結構）| ❌ germline ASM | ❌ | published 2023 |
| ASMS | bulk ONT | ✅ no-phasing | ✅ | ❌ | ◐ **5000-perm** 甲基差量 | ❌ germline ASM | ❌ | preprint 2024 |
| MethylBERT | bulk ONT-capable | ✅ transformer | ✅ | ❌ | ❌ **supervised** deconv | ◐ tumour/normal deconv | ❌ | published 2025 |
| single-cell 甲基（Epiclomal/MethylTree/Sgootr）| **single-cell** | per-CELL | ✅ lineage signal | ❌ | ◐ sc regime | **✅ 真重建** | **✅ sc lineage tree** | published；**不同 regime 非競爭** |
| SNV-CCF（PyClone-VI/DPClust）| bulk WGS/WES | ❌ aggregate VAF | ❌ | ❌ | ❌ CCF cluster | ✅ CCF 群 | **✅ per-mutation deconv** | cohort-validated；**互補骨幹** |

**ISM 唯一在「unsupervised 距離矩陣結構檢定 + normal-anchored cis + somatic-subclone target」三軸同時 ✅**（= sanctioned 整合），但每軸單元件皆 prior art。`genome-wide reconstruction` 欄 ISM 刻意 ❌ = others do better 之處。

---

## §11 外部文獻錨（L3 — 投稿前過 /citation-verification）

> ⚠ 皆 agent WebFetch（L3，帶 PMID）。**2025-26 競品（下方標 🆕）= L3-preprint-abstract（bioRxiv/medRxiv 403-blocked，未全文核）→ 投稿前必全文核 + head-to-head，並建 CONTEXT 卡入庫**。

**遺傳骨幹（學派 1-4）**：Roth 2014 PyClone PMID 24633410 · Gillis&Roth 2020 PyClone-VI PMC7730797 · Miller 2014 SciClone PMID 25102416 · Caravagna 2020 MOBSTER · Deshwar 2015 PhyloWGS PMID 25786235 · Jiang 2016 Canopy PMID 27573852 · Van Loo 2010 ASCAT PMID 20837533 · Ha 2014 TITAN PMID 25060187 · Satas 2021 DeCiFer PMID 34416171 · Dentro 2021 PCAWG-11 · Tarabichi 2021 guide PMID 33398189 · Şenbabaoğlu 2014 PMID 25158761.

**single-cell / 甲基-lineage（學派 5）**：Gaiti 2019 PMID 31092926 · Epiclomal PMID 32966276 · MethylTree PMID 39820752 · EPI-Clone PMID 40399669 · Sgootr 2023 · Liu-Goretsky 2025 PMID 40950124.

**long-read phasing + read-level 甲基（學派 6）**：MethPhaser 2024 PMID 38909018 · cvlr 2023 PMID 36726731 · ASMS bioRxiv 2024 · CpelNano · MethylBERT 2025 · O'Neill 2024 PMC11605692.

**🆕 2025-26 競品（2026-06-21 已全文核，見 §12）**：**ROCIT**（Baker/Spellman bioRxiv 2026.03，PMC12991090，**fulltext-verified**，Tier A — read-level 甲基 tumour-origin 分類，supervised/binary/**無 normal**；🔴 **PacBio HiFi 實證 ONT 僅 aspirational**）· **TumorLens**（Paulin/Sedlazeck medRxiv 2026.03，PMC13015642，**fulltext-verified**，Tier B — 統一 bulk-ONT 多組學偵測，ASM 限 74 HLA/APM 基因為 immune-escape）· **MethPhaser-Cancer**（🔴 **非論文 = Yilei Fu 2023 Rice 博論章節**，無 DOI/code，purity+k=2 read 二分非 subclone，Tier C L3）· **O'Neill 2024**（Cell Genomics，PMC11605692，Tier A peer-reviewed，**非 subclone 論文**：term-search DPClust/subclone=0）+ **McAuley 2025 rectal**（bioRxiv，單病人，Tier C abstract-only）。

**LOH-unmask / over-cluster confound（甲基特有）**：Martin-Trujillo 2017 PMID 28883545 · Chase 2015 PMID 26114957 · Rosenski 2025 PMID 40069157 · M3C 2020 PMID 32020004 · Barrett 2017 PMID 28743252 · Landau 2014 PMID 25490447 · Wolf DNMT1 processivity PMID 36285438.

---

## §12 競品全文核結果（2026-06-21 — 4 背景 agent，已建 CONTEXT 卡）

> 對 §11 標 🆕 的 4 競品做一手全文核（Europe PMC / PMC / Semantic Scholar 繞過 bioRxiv 403）。**頭條：四者 whitespace_threat 全 = LOW → ISM 三軸整合白地經實測仍站得住。** 已建/更新 5 張 external_validation 卡（庫 69→72）。

| 競品 | 全文核狀態 | 對 ISM 三軸 | whitespace | tier | 卡 |
|---|---|---|---|---|---|
| **ROCIT**（Baker/Spellman 2026，PMC12991090）| **fulltext-verified** | supervised transformer / binary tumour-vs-非tumour / **無 normal** / 甲基=driver / 非 subclone → **0/3 重疊** | LOW（最尖反襯）| **A** | 新 axis5 `rocit-tumor-read-classifier-2026` |
| **TumorLens**（Paulin/Sedlazeck 2026，PMC13015642）| **fulltext-verified** | 多組學偵測 + 74-基因 ASM(immune-escape)；10kb window 非 read-level；非 subclone → 0/3 | LOW-MED（同生態系 must-cite）| **B** | 新 axis2 `tumorlens-multiomic-longread-2026` |
| **MethPhaser-Cancer** | dissertation-only（L3 second-hand）| purity + k=2 read 二分；停 phasing 層；非 subclone | LOW | **C** | addendum 於 `methphaser-germline-methyl-phaser-2024` §8 |
| **O'Neill 2024**（Cell Genomics，PMC11605692）| peer-reviewed fulltext | germline-aDMR read-level phased（**非 subclone**，term-search=0）| LOW(niche)/MED(broad) | **A** | 更新既有 `oneill-...-2024`（+term-search 佐證）|
| **McAuley 2025 rectal** | **abstract-only**（全文 403）| 單病人 concordance + DMR-vs-normal（「分離 DMR track」例）| LOW | **C** | 新 axis2 `mcauley-rectal-cancer-longread-2025` |

**🔴 三個一手更正（誠實核對的價值）**：
1. **MethPhaser-Cancer 非已發表論文** = Yilei Fu 2023 Rice 博論章節（無 DOI/code/僅模擬）；做 purity + k=2 read 二分非 subclone；先前「2025-26 延伸論文窄化白地」= Europe PMC 查無、疑 search-summary 幻覺，已排除。
2. **ROCIT = PacBio HiFi 實證（6 樣本），ONT 僅 aspirational** —— 引用勿稱「bulk ONT」，須寫「PacBio-demonstrated / ONT-planned」；且 ROCIT 無 matched-normal（in-silico 從 clonal SNV/LOH 標籤）。
3. **O'Neill 2024 不做 DPClust/subclone**（全文 term-search=0）——我先前「DPClust-on-LR field default」框架**錯**；「DPClust + 平行 DMR」的 field-default pipeline **作為單一已發表論文並不存在**；O'Neill 的甲基是 read-level haplotype-phased aDMR（germline 軸）非「分離 DMR track」。

**白地裁決（更新）= PARTLY-FILLED**：ISM 三軸交集（unsupervised read×read 結構檢定 + normal-baseline somatic cis-test + somatic-subclone target on bulk ONT）**經 4 競品實測無單一方法佔滿**；最尖 ROCIT 三軸 0/3。屬「先驗元件的新穎整合」，**窄誠實 framing 仍強制**（ROCIT/TumorLens 是 must-cite must-distinguish 鄰居）。⚠ McAuley 全文不可得（403）、MethPhaser-Cancer 博論 PDF 未讀 → 投稿前須補。

---

*配套教學指南：`InterSubMod/docs/method_comparison/20260621_clone_subclone_teaching_guide_for_ai_01.md`。*
*內部數字 SoT：situation inventory `20260617` + memories（見 frontmatter）。外部庫：`/big7_disk/liaoyoyo2001/external_validation/`（72 源，2026-06-21 +ROCIT/TumorLens/McAuley）。*
*關聯 memory：`project_clone_subclone_landscape_and_ism_feasibility`（本 doc 摘要）。*
