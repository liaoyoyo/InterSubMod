<!--
建立時間: 2026-06-21
報告類型: methodology — baseline-discipline 概念定義詞典 + 字詞關係圖（治本：防「中間產物當結果 / baseline 藏起來 / 同詞兩義」概念模糊）
任務類型: D handoff — 全專案 + 接手 AI 引用的單一概念真值；所有 doc 用詞應對齊本檔
build_branch: research/subclonal-reconstruction-202606
status: 內部定義已立 + 外部驗證完成（3 agent a1609556/a9ea12bb/ae55213f）；§2.5 = verdict（1 誤用 SampleASM + 4 NEEDS-QUALIFIER + 1 關係更正 co-occurrence→lineage；其餘 MATCH）；外部 citation=L3，投稿前 /citation-verification
data_sources:
  - 本 session 已立紅線 + memories（apriori/tumor_only/cluster_label/situation_inventory/confirmation_construction）
  - 外部驗證庫 74 源（agent 進行中確認名詞定義+關係 vs 領域）
provenance_note: 內部定義溯源 memories/docs（verified）；§2「外部標準定義」欄 = agent 文獻核（L3，進行中），回來填 + 標 MATCH/DIVERGE。
trigger: 2026-06-21 用戶 catch「甲基共分離 ≠ 佐證（cis vs lineage baseline 不同）」→ 一般化為 baseline-dependence 原理。
-->

# Baseline-discipline 概念定義詞典 + 字詞關係圖

> **這份是什麼**：全專案 + 接手 AI 的**單一概念真值**。把容易把「中間產物當結果」「藏起 baseline」「同詞兩義」的詞，逐一定義成 **測什麼 × 對照誰(baseline) × 因此意義 × 不能說**，並畫出**字詞之間的關係**（§3）。**所有 doc / 回覆用到這些詞時應對齊本檔。**
> **緣起**：2026-06-21 用戶 catch「甲基隨 SNV 共分離不能直接當佐證 —— cis 基準 vs lineage 基準不同」→ 一般化。

---

## §0 統一原理（baseline-dependence 鐵則）

> **任何「測量出來的量」（甲基分群 / Δβ / ASM / 對齊 / 結構 / k 群 / co-occurrence）都不是結果，是中間產物。它的科學意義 = f(量, 你拿它對照的 baseline)。當一個詞把 baseline 藏起來，就會把中間產物誤當答案。**

**操作守則**：寫任何此類 claim，必明寫三件 —— **(1) 測什麼、(2) 對照誰(baseline/null/label/control)、(3) 因此意義 + 不能說什麼**。缺 (2) = 概念模糊，禁出。

> 🔴 **外部驗證（2026-06-21，stats agent）**：「測量意義 = f(量, baseline)」是 **sound 原理但非單一可引用的定理** → **投稿時禁當一個 citation**；要拆成 ① 聚類顯著性的 null-reference 原理（M3C / Şenbabaoğlu correlation-preserving null）+ ② confound 的 control/counterfactual design。這是唯一「true 但 mis-attributable」的點。

---

## §1 🔴🔴 紅線（傳播到全專案 + 接手 AI 必守）

1. **甲基分群/共分離 ≠ 結果/佐證**：意義依 baseline 而定。甲基隨 **somatic-變異-標籤**（ALT/carrier）分開，在 **cis 基準**下 = **cis-ASM（同一突變事件、非獨立證據）**；須先過 **normal cis-control + 遠端獨立位點**扣 cis，**只有殘餘共分離才算獨立 lineage 佐證**。
2. **對齊 a-priori 標籤 ≠ subclone**（85% 對齊 germline/carrier = cis-ASM）。
3. **有結構（PERMANOVA）≠ subclone**（存在性 ≠ 身分）。
4. **切 k 群依 null 而定**（i.i.d. null 過度偵測；須 correlation-preserving / a-priori 對照）。
5. **單樣本能「驗證」不能「確認」**（confirm 需 single-cell/multi-region/longitudinal ground truth）。
6. **co-occurrence ≠ confirmed subclone**（需 normal 排 germline/artifact）。
7. **同詞兩義必界定**：「reconstruction」「ASM」「Δβ」「驗證」單詞出現必標 baseline/軸。

---

## §2 詞典（每詞：測什麼 × 對照誰 × 意義 × 不能說 × 外部定義[pending]）

> 🔴 「外部標準定義」欄 = 文獻驗證進行中（agent），回來補 + 標 MATCH/DIVERGE。

### 甲基類
| 詞 | 測什麼 | 對照誰(baseline) | 因此意義 | 🔴 不能說 | 外部定義 |
|---|---|---|---|---|---|
| **甲基分群**（read×read 距離 cluster）| read 間甲基距離→群 | ①i.i.d.null ②corr-null ③a-priori 標籤 ④normal | 依 baseline | 分群=subclone | {{pending}} |
| **甲基共分離 / SubcloneDbeta** | 甲基隨 carrier-vs-germline tag（或 SNV 組合）Δβ | **cis（normal+遠端）vs lineage** | 扣 cis 前=cis-vs-lineage 未分(=G-B)；扣 cis 後殘餘=lineage 佐證 | 共分離=獨立佐證/subclone | {{pending}} |
| **ASM**（allele-specific methylation）| 兩 allele/haplotype 間甲基差 | **哪兩組**：HP-axis(HP1/HP2) / ALLELE-axis(ALT/REF) / SampleASM(tumor/normal) | 依軸；germline-het ASM ≠ somatic ASM | 「ASM」不標軸；SampleASM(飽和96.5%)當訊號 | {{pending}} |
| **Δβ**（beta 差）| 兩組 mean β(=M/(M+U)) 相減 | **哪兩組**：SubcloneDbeta / GermlineAsmDbeta / SomaticResidual | 依組（意義天差地別）| Δβ 不標哪兩組 | {{pending}} |
| **cis-ASM / SD-ASM** | 變異→局部甲基 | 是否變異**直接**造成（vs 譜系獨立）| 同一事件的下游 | cis-ASM=獨立 lineage 證據 | {{pending}} |
| **epiallele / PDR** | read 內 CpG 共甲基 pattern | — | 單 clone 內亦有（= over-cluster 根源）| epiallele 多群=subclone | {{pending}} |

### subclone / reconstruction 類
| 詞 | 測什麼/指什麼 | baseline/前提 | 意義 | 🔴 不能說 | 外部定義 |
|---|---|---|---|---|---|
| **reconstruction** | ISM=regional LOH-constrained same-haplotype partition | 領域定義=clone#+CCF+lineage tree+CN/purity deconv | 同詞兩義 | ISM 做領域定義的 full reconstruction | {{pending}} |
| **clone vs subclone** | clone=共享突變集群體；subclone=子集/後代 | CCF~1(truncal) vs 低(branch) | — | clone=subclone | {{pending}} |
| **CCF** | f·m/(ρN+2(1−ρ)) | **需 multiplicity m + purity + CN** | 細胞分率 | VAF=CCF（少 m）| {{pending}} |
| **multiplicity** | 帶突變的拷貝數 | VAF 單獨不可定 | — | 半帶 ALT=subclone | {{pending}} |
| **co-occurrence**（same-molecule vs same-haplotype）| read 上 SNV 組合 | same-molecule(強) vs haplotag 延伸(弱,依 phasing) | 局部譜系直接證據 | co-occurrence=confirmed subclone | {{pending}} |
| **HP-tag 1-1/2-1/1-2/2-2** | somatic_haplotag **演算法指派** | baseline=phasing 正確 | 1-1/2-1=第一 somatic haplotype；-2=第二(multi-subclone marker) | tag=已確認 subclone | {{pending}} |

### 統計 / baseline 類
| 詞 | 測什麼 | baseline | 意義 | 🔴 不能說 | 外部定義 |
|---|---|---|---|---|---|
| **結構 / PERMANOVA** | 群間 vs 群內距離 + null | null + PERMDISP(location vs dispersion) | 有結構存在 | 有結構=subclone | {{pending}} |
| **對齊 / CramérV** | 切群×a-priori 標籤列聯 | **哪標籤**：random/HP-family/carrier-germline | 對齊獨立軸=非循環確認≥2群 | 對齊=subclone(85% cis-ASM) | {{pending}} |
| **double-dipping** | 在資料上挑群再測同資料 | — | 循環推論 | unsupervised 切群+測=判別器 | {{pending}} |
| **excess-over-null** | 跨樣本復發率 − null | null | 跨樣本現象 | raw rate（必看 excess）| {{pending}} |
| **驗證 vs 確認** | validate(對 null/label) vs confirm(對 ground truth) | single-cell/multi-region/longitudinal | 單樣本=validate=characterization | 單樣本「確認」subclone | {{pending}} |

---

## §2.5 外部驗證結果（3 agent，2026-06-21；外部=L3，投稿前過 /citation-verification）

> **結論：詞彙絕大多數 field-correct。1 個誤用（SampleASM）+ 4 個 NEEDS-QUALIFIER + 1 個關係更正。** 填 §2 各詞「外部定義」欄 = 下表。

| 詞 | verdict | field 重點（citation）| 🔴 action |
|---|---|---|---|
| **ASM** | NEEDS-QUALIFIER | field ASM **永遠 axis-relative**（Shoemaker2010 PMID20418490 Cat I/II/III）| 永遠標軸；🔴 **「SampleASM」=misnomer（tumor-normal 是跨樣本 DMR 非 allelic）→ 改「tumor–normal differential methylation」**；「ASM」只留 HP-axis/allele-axis |
| **cis-ASM/SD-ASM** | MATCH | 「cis」正確（Min2021 34493871/Onuchic2018 30150324/Do2020 32594908）| 明寫「somatic cis」（SD-ASM 通常指 germline）|
| **甲基共分離（TERM6 最大暴露）** | 🔴 RISK | field **用 CpG-SNP/haplotype cis-control 分 cis vs lineage**（Shoemaker Cat I/II；Onuchic-Do haplotype-ASM PMID33077733）| 改「consistent with cis-ASM **或** lineage，單區證據不能分」=G-B；**扣 cis 後**才可說「刻畫 subclone」|
| **epiallele/PDR** | MATCH | Landan2012/Landau2014 | 定位 **structure-not-disorder**（對比 PDR）|
| **Δβ/β** | MATCH 定義/RISK 應用 | **ONT β=read-fraction(threshold/depth-dependent)≠array intensity ratio**；M-value 較佳 | Methods 講清估計量；Δβ 顯著性走 **beta-binomial**（=FIX-1）|
| **甲基分群（read×read）** | MATCH 框架/RISK | cvlr+ASMS **都有 permutation 檢定** | 🔴 **禁「只有 ISM 有顯著性檢定」**；delta=檢定**對象**(距離矩陣結構 vs Δmeth 量級)+normal-anchor+somatic target；配 PERMDISP |
| **reconstruction** | NEEDS-QUALIFIER | Tarabichi 5-step 定義（PMC7867630）| 「reconstruction」**綁 somatic-haplotagging**；甲基永遠用「characterize」動詞；禁「we reconstruct subclones using methylation」|
| **clone/subclone** | MATCH（scope caveat）| Tarabichi：subclone=descendant CCF<1 的 cell-population | ISM 是 regional read-partition→寫「**consistent with** subclone」非「=subclone」（除非 CCF 錨）|
| **CCF** | MATCH（abstention）| =f·m/(ρN+2(1−ρ))，需 m+purity+CN（DeCiFer 34416171）| ISM 不算 CCF→**禁把 VAF 當 CCF**、禁 prevalence 排序 |
| **multiplicity** | MATCH | VAF 單獨不可定 m（需 allele-specific CN）| ISM held-CN/LOH **刻意繞 multiplicity=strength**（cite DeCiFer/Tarabichi）|
| **lineage/phylogeny** | MATCH | **lineage(descent path)≠cluster(CCF 相似群)** | 🔴 ISM 有 **cluster/partition 非 lineage**；禁說 phylogeny/tree/lineage tree |
| **co-occurrence/linkage** | MATCH | Foltz PMID39149342 genetic 先驗 | 增量=加甲基軸**非 phasing 本身**（longphase-S 已做）|
| **same-molecule vs same-haplotype** | NEEDS-QUALIFIER | molecule(直接)>haplotype(phasing-依賴)；**同 haplotype ≠ 同 subclone** | ISM 多在 haplotype 層→勿當 same-molecule/same-cell |
| **PERMANOVA** | MATCH | Anderson2001；**非循環**因 HP 標籤獨立 | 配 PERMDISP；非 Euclidean 距離合法（McArdle&Anderson）|
| **PERMDISP** | MATCH | Anderson2006 location≠dispersion | dispersion 警告**須 surface**（LOH 不平衡易假陽）|
| **double-dipping** | MATCH | Kriegeskorte 19396166/Gao-Bien-Witten | tumor-only 無監督=double-dip（已 NEGATIVE 正確）；a-priori 軸=remedy |
| **baseline/null** | MATCH 原理/非單一定理 | M3C/Şenbabaoğlu correlation-null | §0 拆引（見上）|
| **validate vs confirm** | MATCH | Tarabichi 單樣本 weakly-informative | project convention over real hierarchy→**Methods 明寫定義** |
| **HP-tag** | MATCH | 演算法指派（phasing 依賴）| 非循環性**完全靠零甲基 phasing 路徑**（載重前提）|

🔴 **兩個載重前提（Methods 必寫 + 引 longphase-S 源碼/commit）**：(a) **HP 標籤非循環性 = 完全靠「phasing 用零甲基」**（甲基若進 phasing 當 tie-breaker → PERMANOVA 變循環）；(b) **over-clustering 的 null 必須 correlation-preserving 非 i.i.d.**。
⚠ 註：subclone agent 回報「axis*/CONTEXT.md 卡不存在」= 其 ripgrep 失效之 artifact；卡實存（本 session 親編 74 張，grep 實測 74）。

---

## §3 字詞關係圖 + 外部驗證 verdict（user 要的「字詞的關係與意義」）

### (a) 升級鏈（necessary-not-sufficient，每步是下一步的**必要非充分**條件）
```
有 reads → 切得出群 → 有結構(PERMANOVA sig + 非 PERMDISP)
        → 對齊獨立 a-priori 標籤(CramérV≥0.7)
        → 扣 cis 後仍共分離(normal cis-control)
        → [仍只是 subclone 候選 / characterization]
        → 需 ground truth(single-cell/multi-region) 才 confirm
```
> 每個箭頭都是「過了才有資格談下一步」，**任一步不可跳級當結論**。多數位點停在前段（85% 對齊但 = cis-ASM）。

### (b) 包含關係（⊃）
- **reconstruction ⊃ {CCF clustering + phylogeny + CN/purity deconv}**；ISM 只做 regional partition = **真子集**。
- **clone ⊃ subclone**（subclone 是子集/後代）。
- **ASM ⊃ {germline-ASM, cis-ASM/SD-ASM, somatic-ASM}**；「ASM」單詞涵蓋多源，必標。
- **confirm ⊃ validate**（confirm 需 ground truth，validate 只需 null/label；validate 是 confirm 的弱版）。

### (c) 對比對（易混，必分清）
- **cis-ASM ↔ lineage**（同事件 vs 獨立；扣 cis 才分）— **本詞典緣起**。
- **驗證 ↔ 確認**（null/label vs ground truth）。
- **same-molecule ↔ same-haplotype**（read 直接 vs phasing 延伸）。
- **location ↔ dispersion**（PERMANOVA vs PERMDISP）。
- **germline-ASM ↔ somatic-ASM**（差異在 normal vs tumor-specific）。
- **切群/對齊/結構 ↔ subclone**（中間產物 vs 身分）。

### (d) 前提關係（需要 X 才能談 Y）
- **CCF 需 multiplicity**（無 m 不能算 CCF）。
- **判 somatic 需 normal**（無 normal 分不出 germline）。
- **lineage 佐證需扣 cis**（不扣 cis = 可能同事件重複計）。
- **confirm 需 ground truth**（single-cell/multi-region/longitudinal）。

### (e) 因果關係
- **epiallele / read-內相關 → over-clustering**（相關灌大 k）。
- **somatic SNV → cis-ASM**（局部、同分子；= 須扣的 baseline）。
- **LOH → 既存 ASM unmask**（Martin-Trujillo 82-91%；非新 somatic）。

### (f) 🔴 外部驗證 verdict（2026-06-21，3 agent）
- ✅ **CONFIRMED（field-standard）**：升級鏈 R1（含 PERMDISP 作 location≠dispersion 子步）· 對齊避 double-dip R2（須標籤獨立於甲基）· epiallele→over-cluster R3 · validate⊊confirm R4（Tarabichi 單樣本「weakly informative of the tree」逐字）· **reconstruction⊃{CCF+tree+CN deconv}**（ISM 真子集）· **CCF 需 multiplicity** · **clone⊃subclone** · **same-molecule > same-haplotype** · identification⊊confirmation。
- 🔴 **CORRECTED — co-occurrence → lineage（我原本寫錯）**：「同-molecule/同-haplotype co-occurrence → 同 lineage」是 **necessary-not-sufficient**。**互斥 → 不同 subclone 可靠**；但**共現 → 同譜系須再加 cellular fraction（VAF/CCF）**才能定 lineage（Foltz NRAS PMID39149342：G13R+Q61K 同 haplotype 卻不同 subclone，靠 VAF 35.7%/22.2% 區分）。已修 `20260621_singlemolecule_multisnv_cooccurrence_task_spec_01.md` §1。
- ⚠ **兩個 Methods 必寫載重前提**：(a) R5「baseline-dependence」非單一定理（拆引 §0）；(b) HP 標籤非循環性完全靠「phasing 用零甲基」+ over-cluster null 須 correlation-preserving。

---

## §4 傳播指引 + 使用方式
- **所有 doc / 回覆**用到 §2 的詞 → 必對齊本檔定義 + 標 baseline；缺 baseline 的句子視為概念模糊。
- **接手 AI** 先讀本檔 §0+§1，再動手。
- **🔴 回頭標註**：既有「甲基有用 A 刻畫 subclone（SubcloneDbeta）」這句 = **未扣 cis → cis-vs-lineage 未分（=G-B）**，引用時必加「待 G-B 扣 cis」前提（見 §1.1 / 詞典「甲基共分離」列）。
- **外部驗證回來後** → 填 §2「外部定義」欄 + §3 標 MATCH/DIVERGE + 修任何用詞與領域衝突處。

---

*關聯 memory：`feedback_baseline_dependence_not_result`（紅線）· `project_apriori_subclone_classification_model`（SubcloneDbeta cis caveat）· `project_subclone_confirmation_construction_ont_transferability`（驗證 vs 確認）。*
