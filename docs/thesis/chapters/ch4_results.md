<!--
建立時間: 2026-06-15；**2026-07-06 主軸重置 v4**：正向核心改為 sSNV 共現骨幹(新 4.1–4.3)、甲基降襯托(原 4.1–4.7 順移 4.4–4.10)。canonical 對回 master_spec §6(mtime 07-05)。演算法定位 pending /citation-verification。
報告類型: 碩論正文草稿 — 第四章 結果
狀態: draft v3；待 review
data_sources:
  - docs/paper_focus/00_共識證據台帳_20260612_01.md (§1 catalog/filter/copy/cross-sample/cis grep-confirmed)
  - docs/paper_focus/02_paper_framework/20260614_Intro主軸定稿_thesis_statement_01.md (旗艦 BRCA2/chr2 數字)
  - research/flagship_chr2_18086020_20260612/ism_out/anchor_18086020/chr2/chr2_18086020/chr2_18071020_18101020/{metadata.txt,clustering/significance.json,reads/reads.tsv} (chr2 PRIMARY ISM 輸出，2026-06-15 本 session 逐項 grep-verified)
  - research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/V6_on_fp/filtered_snv_fp/chr2/chr2_18096269/chr2_18091269_18101269/clustering/significance.json (18096269 對照位點 PRIMARY，2026-06-15 grep-verified)
provenance_note: 4.1-4.4/4.7 引共識台帳 §1（🟢）；4.6 chr2 兩位點 2026-06-15 本 session 皆從 PRIMARY ISM 輸出獨立驗證 — 18086020(57×203/32+25−、HP1 19 全REF、HP1-1 36=28ALT/8REF、F=10.6120/V=0.6689 p=0.027/k=2、allele V=0.2699 p=0.097 NS)、18096269(32 reads、global_alt V=0.9671 p≈0、global_hp V=0、k=5、cluster F=39.9956)，全吻合；4.5 BRCA2 Δβ=−0.122/d_within=−0.023 p=0.024 引 06-14 定稿（待 V-4 primary 重驗）。已套用 4 校正。誠實護欄：case-demo、時序=假說、hypo≠silencing、association 非 causality、⭐3。
-->

# 第四章　結果

本章以體細胞 sSNV 單分子共現重建之區域級克隆骨幹為正向核心（4.1–4.3），再以甲基化之能力邊界與去混淆脈絡為襯托（4.4–4.7），最後以兩案例示範甲基與遺傳訊號之交會（4.8–4.9）及所依附之相位骨幹與跨樣本復現（4.10）。正向段顯示：全基因體 35,332 個 somatic sSNV 中 61.0% 於 read-span 內具可建樹共讀連鎖（4.1），據此於 3,885 個具 genotype 向量之區域量化可決定性——45.4% 達可唯一確定之克隆拓撲、為軟上界，逐層加嚴後穩固核心收斂至 44 區（4.2）；三項演算法成果進一步從既有共現資料萃取拓撲資訊（4.3）。襯托段則界定甲基化之邊界：體細胞 ASM 雖存在且略偏 somatic（4.4），但甲基化無法做變異真偽過濾器（4.5）、read 層級可靠之甲基分群壓倒性跟隨 germline 而非 somatic（4.6）、且此分群強度非源於拷貝數（4.7）。

本章紅線：機制建於單樣本 HCC1395 與單一定相流程（⭐3，TP/FP 取自 tumor-only caller）；旗艦案例為**案例示範（case-demonstration）**而非已驗證之自動化重建，不輸出譜系樹，時序為數據支撐之假說；甲基與突變之關係表述為**關聯／共定位**而非因果。各數值真值來源標於各節（共識台帳 §1 grep-confirmed；旗艦數字引 06-14 定稿與 flagship 補跑）。

---

## 4.1　sSNV 宇宙帳本、連鎖普查與區域拓樸型態

### 4.1.1　全 sSNV 宇宙與 per-site 連鎖普查

以 HCC1395（SEQC2 truth set）ONT 長讀資料為主分析對象，全基因體 somatic sSNV 宇宙為
**35,332 個位點**（其中 SEQC2 標記 TP 30,490、FP 4,842；TP/FP 標記僅供事後觀察，
**不進入任何前處理或分桶**）{{cite:SEQC2_truth}}。

對每一個 sSNV，我們並行普查其在 read-span（≤ 50 kb）內的**共讀連鎖屬性**，得到三桶
（此為 *per-sSNV census*，非過濾步驟——三桶本身不決定後續哪些位點進入重建）：

| 連鎖桶 | 位點數 | 佔比 | 意義 |
|---|--:|--:|---|
| **linked（可建樹）** | 21,554 | 61.0 % | read-span 內有 same-PS partner 且有共讀連鎖 |
| **underpowered** | 5,458 | 15.4 % | 有 partner 但共讀不足；加深覆蓋可救，且 census VAF→CCF 仍可刻畫 |
| **isolated** | 8,320 | 23.5 % | read-span 內無 partner；caller VAF 仍可刻畫 |

三桶合計 35,332（校驗 ✓）。**須強調**：單位點**並非「全然無法處理」**——underpowered 與
isolated 皆可經 VAF→CCF 刻畫，真正拓樸-dead 僅限 isolated 中無 same-PS partner 者。

> ⚠ **資料細節（已閉環）**：HCC1395 tumor BAM 有 ~18 %（1.18×）read 具兩份 primary alignment。
> 因骨幹之 co-read 計數為 **QNAME-keyed 自動去重**，三桶分布**完全不受此重複灌水影響**，
> 無需 dedup 重跑（唯一受影響量=絕對 depth，該量不用於分桶）{{cite:internal_S1}}。

### 4.1.2　區域切分與拓樸型態分佈

重建區域由 S1 階段以**相鄰 gap > 50 kb** 切割 linkage group 而成，此步將**全部 35,332 個位點**
分組為區域（50 kb 為 ONT read-span 的物理上界），**非僅取 linked 桶**。全基因體共得 **7,143 個
含 ≥ 2 sSNV 的區域**，其 pairwise 拓樸型態分佈如下（07-04 canonical）：

| 拓樸型態 | 區數 | 說明 |
|---|--:|---|
| single | 1,995 | 單一克隆 genotype |
| branched（直系＋姊妹） | 1,110 | 分支結構，祖先已觀測 |
| branched（需隱藏祖先） | 30 | 分支需推斷未觀測祖先 clone（A_det 23 + A_amb 7） |
| linear | 741 | 線性巢狀 |
| germline | 371 | 僅 germline 向量、無 somatic 結構 |
| no-vector | 2,887 | 無 read 全覆蓋所有 sSNV，拿不到聯合 genotype 向量 |
| incompatible | 9 | 四配子違反（重分類後之殘餘，見 §4.3.3） |

合計 7,143（校驗 ✓）。**須注意三點**：(i) 「需隱藏祖先」30 區**不併入** branched 1,110；
(ii) germline 371 區為獨立型態、不可省略；(iii) no-vector 2,887 區並非無樹，可經 pairwise
傳遞式串接救回（§4.3.1）。

---

## 4.2　determinacy census 與 robustness ladder

### 4.2.1　determinacy census

在具備 genotype 向量的 **3,885 個區域**（determinacy 分母）上，逐區判定其克隆拓樸的可決定性：

| determinacy 類別 | 區數 | 佔 3,885 |
|---|--:|--:|
| **A_determined**（唯一確定樹） | **1,764** | **45.4 %** |
| A_ambiguous | 69 | 1.8 % |
| B_pairwise（僅 pairwise 可定） | 943 | 24.3 % |
| C_underdetermined | 544 | 14.0 % |
| recurrence_required | 70 | 1.8 % |
| incompatible | 18 | 0.5 % |
| （未歸類殘差） | 477 | 12.3 % |

**誠實標註**：上表六個具名類別合計 **3,408** 區；其餘 **477 區於來源 canonical 表中未歸類**，
本文據實列為「未歸類殘差」，**不另行命名或填補**（六類 3,408 + 未歸類 477 = 3,885，校驗 ✓）。

> ⚠ **gap#1 subcube 分母擴充（2026-07-04 canonical）**：以 partial-read subcube 救回
> 「無 full-coverage 但 ≥2 位點」的區域後，determinacy region universe 由 **3,885 擴至 6,288**
> （HCC1395 救回 **2,403 區**，透明標 E_subcube_recovered，不僭稱 A_determined）。
> 🔴 **A_determined 1,764 完全不變**，determinacy rate 仍以 full-cov census 分母 3,885 報（45.4%）。
> （僅 HCC1395；六樣本 subcube 未重跑。）

### 4.2.2　robustness ladder：45.4 % 是軟上界

A_determined 的 45.4 % **是可決定性的軟上界（soft upper bound）**，而非穩固克隆核心的比例。
逐層加嚴條件後，真正穩固的核心大幅收斂：

| 層 | 累加條件 | 區數 | 佔 3,885 |
|---|---|--:|--:|
| **L1** | A_determined（單分子向量可定） | 1,764 | 45.4 % |
| **L2** | ＋ TP-backed（fp = 0, tp > 0） | 910 | 23.4 % |
| **L3** | ＋ 多位點（n_sSNV ≥ 3） | 239 | 6.2 % |
| **L4** | ＋ 非 CN-gain（真正穩固核心） | **44** | **1.1 %** |

> 🔴 **最大未控 confound = CN-gain multiplicity**。determinacy rate 在 CN-gain 區最高
> （gain 53.9 % vs neutral 44.5 % / loh 36 % / loss 10.8 %），其機制為拷貝數增益造成的
> multiplicity 可能製造**假共現**、進而製造假 determined。因此任何 tree-level 主張前均應
> mask 或控制 CN-gain；L4 的 44 區即為施加此控制後之殘存穩固核心。
> **另須誠實標**：L4 核心中 80.1 % 僅含 2 位點、且救回的隱藏祖先分支含**推斷（非觀測）祖先 clone**。

> **可辨識度的物理根因**：拓樸可辨識比例 ~11 %（穿越充分 12.3 % / 拓樸可辨識 10.9 %）；
> A_determined 率隨 read-span 增大而崩潰（> 50 kb 僅 ~11 %）——這正是本方法「**區域級
> （≤ read-span）而非 genome-wide 克隆樹**」定位的物理根因，亦是 ⭐3 邊界的來源。

---

## 4.3　三項演算法正向成果

本節報告三項可主張為方法貢獻的演算法成果（定位 `pending /citation-verification`）。
三者皆**不新增乾淨 subclone 生物證據**，而是**更正確地從既有共現資料中萃取拓樸資訊**。

### 4.3.1　pairwise 傳遞式串接：救回空向量區

no-vector 區（無 read 全覆蓋所有 sSNV）可經 **A–B 連、B–C 連 ⟹ A→B→C** 的傳遞式串接重建，
**不需單分子整跨**（理論根據：二元字元 perfect-phylogeny 之 pairwise 相容即足以定樹）。
coverage-tier 層級共救回 **2,140 區**（tree_shape：linear_nested 760 / sibling_only 546 /
full_tree 496 / co_linked_lineage 338）；partial-read soft-likelihood 量化顯示，空向量 2,887 區中
**可救 80.2 %（2,316 區）**、其中 chain ≥ 3 loci 者 58.6 %（1,692 區），真單位點不可救僅 19.8 %（571 區）。

> ⚠ **勿混母體**：coverage-tier 救回的 **2,140** 區 ≠ determinacy 桶內的 **B_pairwise 943**。
> 前者是全 7,143 分層中 no-vector 靠 pairwise 成樹的數；後者是 3,885 有向量區**內**的 determinacy 標籤。

### 4.3.2　gained-pair 定序：消除假等機率

原枚舉器（`enumerate_candidate_trees.py`）將「缺中間群」的 37 區一律標為**等機率候選**，
但其僅依全跨 read 的 population 向量判斷，**未查缺失中間群突變的 pairwise 2×2**（該 2×2 已含定序資訊）。
補上 pairwise 三層判定後（3,885 有向量區）：clean 單突變邊 3,845 / BLOCK 15 / **ORDERED 24** /
CONFLICT 1 / **AMBIG_NOCOREAD（真等機率）= 0**。此結果已落地 production：等機率候選 **37 → 0**、
每區候選數 total 84 → 37（每區收斂單一解），且 24 個定序方向與 VAF 交叉核對 **24/24 一致**。

> ⚠ **與 gap#2 identifiability 不衝突**：此處「真等機率 = 0」指 gained-pair 定序層
> 消除「缺中間群」假等機率（37 → 0）。另一層 gap#2 對 69 個 A_ambiguous 全枚舉 → **44 收斂 determined
> + 25 真歧義**（candidate_trees 上層，不改 determinacy census 的 A_determined 1,764）。兩者為不同分析層，
> 「25 真歧義」仍存在、非被此步消除。

### 4.3.3　incompatible 重分類與隱藏節點發現

由 Model A 形式化反查，發現 `solve_topology` 的 `_laminar()` 檢定**過嚴**：它檢查 genotype
ROWS 的 altset 巢狀/互斥，但 perfect-phylogeny 的正確檢定應為**字元（欄/位置）兩兩四型相容**
（Gusfield {{cite:Gusfield1991}}）。姊妹 clone 共享祖先突變、各有私有突變，其 altset 天生
非 row-laminar，即使祖先已觀測仍被誤殺。以 pairwise 四型正確拆解 118 個 incompatible 區：

| 類別 | 區數 | 判準 | 動作 |
|---|--:|---|---|
| ✅ 真可救 | **30** | nic = 0 且無環 | 建隱藏祖先深分支樹（A_det 23 / A_amb 7） |
| 🟡 非-gain 真四型 | **70** | nic > 0 且 CN ≠ gain 且無環 | 重標 `recurrence_required`（Model A 候選） |
| 🔴 真 likely-artifact | **18** | 有環，或（nic > 0 且 CN = gain） | 維持 incompatible |

HCC1395 headline：incompatible **118 → 18**（校驗 30 + 70 + 18 = 118 ✓）。

> ⚠ **v1→v2 誠實更正**：早期曾宣稱「117 可救 / 118→1」，此為 over-claim（跑在有損 population 向量上）；
> 以 pairwise co-read 四型正確檢定後，真可救僅 **30** 區。此修正只把**假 incompatible 正名**，
> **不新增乾淨 subclone 證據**。

### 4.3.4　Model A m-通道拆分 recurrence_required（V4, 2026-07-04）

70 個 recurrence_required 區的核心難題是：**多重度假象（multiplicity artifact）與真 recurrence
在觀測上不可分**，須有**獨立於共現的多重度通道 m_f(x)** 才能拆分。此通道以外部正交拷貝數
（SEQC2 整數 CN ＋ SAVANA WGS allele-specific CN，purity ρ = 0.96）落地，**非由共現自身導出**：

```
拆分前：118 = 30 救回 + 70 recurrence_required(待拆) + 18 incompatible
拆分後：118 = 30 救回 + 24 候選真 recurrence(m=1) + 64 廣義 artifact(m>1，含原 18)
        其中 24 = 9 SEQC2 neutral/loss + 15 SAVANA LOH copy-neutral
校驗：30 + 24 + 64 = 118 ✓；70 = 24 留 + 46 棄 ✓
```

<!-- 合併者注意（不入正文）：master_spec §6 目前把 recurrence 70 回填至 SEQC2 整數 CN 階段(artifact 11+candidate 9+LOH_unresolved 50);本節引 V4 SAVANA allele-specific 最終拆分(24=9 SEQC2+15 SAVANA LOH,46 棄)。同步 master_spec §6 至 V4 最終後刪此註。 -->

> **誠實邊界**：SAVANA 對 50 個 LOH 區的貢獻是把「total = 2 backbone」從**假設**變為**正交量測證實**
> （49/50 confirm copy-neutral），**不翻轉計數**。24 個候選真 recurrence 為 **⭐3 單樣本弱主張**，
> 升級需 orthogonal 非甲基佐證（多樣本 recurrence 一致性），非本單樣本可定案。

---

## 小結（§4.1–4.3）

正向骨幹在 HCC1395 上重建腫瘤區域級克隆結構：35,332 sSNV 中 61.0 % linked 可建樹；
7,143 區的拓樸型態完整列舉；3,885 有向量區中 45.4 % A_determined（**軟上界**），
施加 TP-backed／多位點／非-CN-gain 三層控制後，真正穩固核心為 **L4 的 44 區（1.1 %）**。
三項演算法成果（pairwise 傳遞救回、gained-pair 定序 37→0、incompatible 重分類 118→18）
更正確地萃取既有共現資料的拓樸資訊；Model A m-通道以外部正交 CN 拆 70 recurrence_required
為 24 候選真 recurrence。全段 ⭐3 單樣本封頂，regional（≤ read-span）非 genome-wide，
分子共現 ≠ single-cell 生物確認。

## 4.4 體細胞 ASM 存在且略偏 somatic，但靈敏度低

TP 位點之顯著 ASM 比例（3.95%）高於 FP（1.07%），子單倍型配對後仍維持（3.77% 對 1.09%，約 3.5 倍），惟絕對靈敏度約 4%（多數 TP 無可偵測 ASM）【🟢 ledger:85】（見 fig2）。

## 4.5 甲基化無法做為體細胞變異真偽之過濾器（四道）

四道獨立證據一致為陰性【🟢，見 fig4】：(一) LOSO 自 +0.02236 塌至 −0.00012、五樣本平均 −0.00004（Wilcoxon p=0.125）；(二) |Δβ| 區分 TP/FP 之 AUC=0.505（落虛無）；(三) 整體靈敏度約 4%，COLO829 TP≈FP（0.0089≈0.0103）；(四) 分層後 FP 反而富集（OR 8.63／4.09／5.84，regression-to-extreme 假象）。此一決策相關之負結果界定甲基化之邊界。

## 4.6 ⭐ 規則：甲基分群跟隨 germline 而非 somatic，乾淨 somatic cis 少見

跨六樣本 332,705 位點之七類標籤目錄中，可靠之甲基 read 分群（TAG-C）12,868（3.87%）經 cis-test 顯示壓倒性為 germline-allelic 背景；latent（TAG-E）28,254 富集與背景相當；無訊號（TAG-G）291,518（87.6%）。乾淨之 somatic cis 於 HCC1395 之 816 可測 HP-axis 位點中僅 1（chr17/TBC1D16，約 0.12%）【🟢 catalog 6/6 驗，見 fig3 主圖】。此「稀少」與既有文獻「somatic 局部 cis 甲基稀少／數值貢獻不大」一致{{cite:Do2020 de-novo ASM / Zhang2019 待 citation-verification}}。

> ⚠ **範圍與口徑**：「少見」之分母為 816 可測位點，**非** 332,705（後者未全 cis-tested，不可將未測當已驗陰性，見驗證 handoff V-1）；可靠分群計數依口徑（TAG-C 12,868／reliable TP+FP 11,689／all 12,876）須標明。

## 4.7 拷貝數非分群來源；以 matched-normal 錨點分解 cis 與 copy

拷貝數整數與 |Δβ| 無顯著關聯（signed ρ=−0.0829），gain 與 neutral 之分群強度無顯著差異（Mann-Whitney p=0.6183），否證 copy-dosage 驅動【🟢，見 fig6】。ISM 以 matched-normal 為錨點計算殘差（raw − normal mean）並以三角量 d_cis／d_drift／d_within 分解差異之來源（此為第三章所述之 normal-anchored cis-test，C2）。此分解使後續旗艦案例得以判別「與突變共定位之甲基差」究屬真 focal cis 抑或 copy／subclone 背景。

## 4.8 ⭐ 旗艦一：BRCA2 promoter — ISM 偵測並分解甲基-變異共定位（C1+C2）

**為何檢視**：BRCA2 為已知癌症基因；於 HCC1395（乳腺）之 BRCA2 promoter，可同時觀察到 read 層級之 haplotype 間甲基差與 tumor-normal 差，且共定位於一 somatic 變異——是檢驗「ISM 能否不只偵測、更分解此類差異」之理想旗艦。

於 BRCA2 promoter，HP-axis（HP1 germline-tag vs HP1-1 somatic-subclone-tag）之甲基率差 Δβ=−0.122，且相對 matched-normal 基線有差（normal-baseline residual 非零），共定位於 BRCA2（TSG）promoter【實際數據，⭐3】。以 normal-anchored 分解（C2）進一步判別其來源：此甲基差**主要 track somatic subclone tag（HP1-1）**，而**直接 focal cis 效應 marginal**（d_within=−0.023，perm p=0.024，惟百分比不 robust）【🟡 分解層】。

此案例之意義在於：ISM **不只標記**「BRCA2 promoter 有與突變共定位之甲基差」，**更量化**此差有多少屬 subclone 結構、多少屬直接 focal cis——這正是本研究欲提供之核心資訊（C2 之價值）。

> ⚠ **誠實邊界**：BRCA2 之甲基差方向為 **hypomethylation**，與 canonical TSG promoter hyper-silencing 方向不同，故**不寫**「promoter 高甲基沉默 BRCA2」；「與突變之關係」為**共定位／候選關聯**而非已證之 cis 因果（focal cis 於分解層 marginal）；單樣本 ⭐3。

## 4.9 ⭐ 旗艦二：chr2:18M — 甲基與多點突變共同標示同一子克隆事件（C3）

**為何檢視**：本節示範主軸——以「多個 somatic 變異 + 鄰域甲基結構」對單一 region 之 read 做 subclone 結構 characterization，並示範甲基如何與突變共同標示同一子克隆事件。

於 chr2:18,071,020–18,101,020（30 kb，HCC1395，ISM 補跑 57 reads〔32 正股／25 反股〕× 203 CpG），錨點 chr2:18,086,020（G→A）。somatic haplotagging 之系譜為：Normal=G → **HP1 母本**（19 reads 全為 REF/G，germline）→ **HP1-1 子克隆**（36 reads＝28 ALT/A／8 REF/G，somatic A 等位富集），另 HP2 為 2 reads（REF）【primary ISM 輸出 `reads.tsv`／`significance.json`／`metadata.txt`，build 069cadb，逐 read 與逐項 grep-verified】。於同一區段之兩個互補位點觀察到兩種模式：

- **18,086,020（標示子克隆之 het 位點）**：非監督甲基分群顯著（cluster PERMANOVA F=10.61，p=0.01，最佳分群數 k=2，無離散假象），且分群與**單倍型／子克隆譜系**高度對應（Cramér's V=0.669，p=0.027；HP-family V=0.694，p=0.003），而與**該位點 het 等位本身**之對應較弱（V=0.270，p=0.097，不顯著）——即甲基分群跟隨的是 haplotype／子克隆譜系，而非單純的 het 等位。
- **18,096,269（區段內 somatic 等位位點，獨立 anchor 窗 18,091,269–18,101,269，32 reads）**：甲基分群與 somatic 等位本身高度對應（Cramér's V=0.967，p≈0，cluster PERMANOVA F=39.996，k=5），而與 HP 無對應（V=0）——甲基跟隨 somatic 等位本身【primary ISM 輸出 `V6_on_fp/.../chr2_18096269/.../significance.json`，2026-06-15 本 session grep-verified】。

兩種互補模式並存——其一甲基跟隨子克隆譜系、其二甲基跟隨 somatic 等位——構成「甲基與多點 somatic 變異於 read 層級共同標示同一子克隆事件」之機制基礎；據此巢狀等位結構與甲基之共變，推導 Normal→HP1→HP1-1 之**時序假說**。

> ⚠ **誠實邊界**：18,086,020 **不在 somatic VCF**（為標示子克隆之 het 位點，其甲基↔HP 對應含 germline-allelic 成分）；subclone 甲基之 somatic-專一性待 within-haplotype 對照（G-B 未跑）；**時序為數據支撐之推論而非正交確認**，**不輸出譜系樹**；自動化 subclone 分析流程與大規模驗證為 ISM 延伸（後續工作，G-D/G-E）；單樣本單 pipeline ⭐3。

## 4.10 相位骨幹與跨樣本復現

旗艦案例所依附之子單倍型骨幹由 LOH-constrained 相位提供：NG=2 時內側 same-HP1 高於外側，方向 7/7（Wilcoxon W=28，p=0.0078125），same-haplotype 占用率（occupancy，非 TP-rate）六樣本皆 ≥93%（0.932／0.990／0.988／0.965／0.983／0.970）；AF 與分群數之關聯源於相位非甲基（HD-4，partial r=0.656）【🟢，見 fig5】。⚠ 此骨幹之正面性部分受 by-construction 循環影響（R-SELFREF 未跑，handoff V-8），單一 pipeline ⭐3。

跨六樣本，相對虛無模型之超出量皆為正（6/6，+0.101 至 +0.241，平均 0.168），跨三癌別復現現象層訊號；惟同一位點之 somatic ASM 為樣本私有（0/38，E[overlap]=0.16，underpowered）。此「現象復現、位點私有」之模式説明旗艦案例為樣本特異，其跨樣本系統性須以標準化判準確立（handoff V-2/V-5）。
