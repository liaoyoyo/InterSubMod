<!--
建立時間: 2026-06-20
類型: 可行性研究 / 方法設計 (A-pilot feasibility)
狀態: in_progress
主軸: Subclonal reconstruction (somatic haplotag + methylation, ONT)
產出來源: Dynamic Workflow wf_5497277c-620 (10 agents, 943K tokens, 174 tool calls)
data_sources: /big8_disk/data/HCC1395/SEQC2/CNV/, InterSubMod/src/core/SignificanceAnalyzer.cpp, InterSubMod/src/core/RegionProcessor.cpp, InterSubMod/docs/CURRENT_FOCUS.md, InterSubMod/docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md, InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_seqc2_concordance_demo_01.md
-->

# ONT CNV/SV → ISM 甲基 read-clustering → SEQC2 truth：subclone 驗證可行性裁決

> **版本**: v1.0 (2026-06-20) · **任務分類**: A pilot（可行性研究）· **scope**: HCC1395-focused（SEQC2 truth 樣本）
> **證據 tier**: **L1**=源碼/磁碟數據第一手重現 ｜ **L2**=文獻/README 二手 ｜ **L3**=推論

---

## 0. 一句話結論（BLUF）

工具選擇明確（**SAVANA 主力 + Wakhan haplotype-specific 互補**），但你提的核心方法路徑——「**用 CN 分群 ↔ 甲基 read-clustering 對齊 → 驗 subclone**」——在現有條件下 **REFUTED as a standalone validator**：因為 CN 軸與甲基軸經 **LOH-unmask-ASM（W2/I2）隱性耦合**，「對齊」可以由**單一個 CN 事件人為製造**，所以對齊≠subclone≠獨立確認。

**路徑「可救」但天花板明確**——唯一能切開混淆的三件事：
1. **CN-neutral（雙 allele 存活）區**仍對齊 somatic label（唯一非循環的肯定式設計）
2. **SV 斷點軸**（split-read 物理證據，是三軸中**唯一**對甲基 β 沒有因果耦合的標籤）
3. **多 sSNV CCF 梯度**（chr2:18M pilot 原型）+ **normal cis-control**

**整體 verdict：PROBE**（值得做，但對外**不可**宣稱「甲基驗證 subclone」——甲基在此框架是 characterize〔有界〕，reconstruction 由 somatic haplotag 骨幹負責）。

---

## 1. ONT CNV 工具推薦（對抗驗證後修正排序）

> ⚠ **對「哪個最好」的天真預期已被對抗驗證修正**：B2 發現 Wakhan 的兩個首選理由（integer-CN、benchmark）都建立在 **medRxiv preprint 自評 + 全文 403 無法一手核實**之上。以**可驗證性 + 成熟度**為主軸，務實排序如下。

| 工具 | 定位 | 關鍵證據 | tier |
|---|---|---|---|
| **SAVANA**（cortes-ciriano-lab）| 🥇 **務實 PRIMARY** | Nature Methods 2025 **已正式發表**；23 releases、v1.3.7 (2026-03)；README **一手明寫** `total and minor ABSOLUTE copy number`(=integer) + purity/ploidy；ONT 原生(`--ont`)；tumor-normal + tumor-only(`savana to`)；99 tumor-normal pairs peer-reviewed benchmark；同時做 somatic SV + CNA | L2 |
| **Wakhan**（KolmogorovLab）| 🥈 **haplotype-specific COMPLEMENT** | README 一手 `haplotype-specific copy number`、HP1/HP2 分離 = **SAVANA 沒有的不可替代能力**；吃 phased VCF 與既有 longphase haplotag 流程**天然銜接**；`--copynumbers-subclonal-enable` 給 subclonal/fractional CN + CCF（僅 segment 層）| L2（preprint 2025-12，未同儕審查）|

**B2 對抗驗證五判準結果**：
- **(a) 跑 ONT BAM** — 兩者 CONFIRMED。
- **(b) per-segment integer CN** — 🔴 **SAVANA 的 integer 證據較硬**（README 一手 "absolute CN"）；Wakhan 的 integer 宣稱依賴 403 全文 = L2 二手。
- **(c) haplotype/allele-specific CN** — **Wakhan 唯一不可替代**；SAVANA 只有 allele-specific total+minor，無 HP 拆分。
- **(d) 維護** — SAVANA 更成熟（peer-reviewed + 23 releases）；Wakhan v0.4.3 (2026-05-28) 也活躍。
- **(e) benchmark** — 🔴 **REFUTED 強讀法**：Wakhan「6 樣本 benchmark」是**作者自建 CASTLE panel 的 self-benchmark**（preprint + AACR abstract），**非獨立第三方評估**，不可寫成「已獨立驗證於我方 5/6 樣本」。

**不推薦**（L2）：araCNA（short-read SNP-based 非 long-read 原生）/ HapCNV（single-cell/haploid 用途不符）/ Spectre（germline-only、無 somatic purity/haplotype CN，僅 germline 大段 CNV sanity check）/ CNVpytor / QDNAseq / Delly / HiFiCNV / Ploidetect。

**Severus**（KolmogorovLab）= Wakhan 的姊妹 SV caller（見 §6），三者構成 `Severus(SV) → Wakhan(haplotype-CN)` pipeline。

> **命名核對未解（L3，落地前必做）**：Wakhan 論文用 `H1937/H1954/H2009/H1437`，與本專案 `HCC1937/HCC1954` 命名不同（既有 HCC1937≠HCC1395 混淆警示）→ 須逐一比對 SEQC2/ATCC ID 確認同細胞系。

---

## 2. per-locus CN 能力（你問的「能否定位每個位置的 CN 數量」）

| 能力 | 可行？ | 工具 | tier |
|---|---|---|---|
| **per-segment integer CN** | ✅ | SAVANA(absolute) / Wakhan | L2 |
| **allele-specific CN**（major/minor）| ✅ | SAVANA(total+minor absolute) / Wakhan | L2 |
| **haplotype-specific CN**（HP1/HP2 分離）| ✅ **唯 Wakhan** | Wakhan（用 tumor allelic imbalance 延伸 phasing 到整染色體、糾正 phase-switch）| L2 |
| **subclonal/fractional CN + CCF** | ✅ 但**僅 segment 層** | Wakhan `--copynumbers-subclonal-enable` | L2 |
| **per-base integer CN** | 🔴 **不可** | CN 是**分段常數**，解析度受 bin（SAVANA 10kb / Spectre 1kb）+ SV breakpoint boundary 限制 | L2 |
| **per-read CN tag / read→subclone 指派** | 🔴 **無任何工具** | Wakhan/SAVANA **全程 segment-level**；README 確認無 per-read CN 輸出 | L2 |

**這正補上內部缺口**：現有 `Coverage_Multiple`（`RegionProcessor.cpp:63-162`）只給**連續 depth ratio**（vs SEQC2 CN r=0.831/0.827），**無 integer/allele CN**；且 Gain recall 僅 14.6%（CN=3 中位 CovM 0.85，系統性右偏壓縮，全域 r 隱藏分層偏誤）。SAVANA/Wakhan 補上 integer + allele/haplotype CN。

---

## 3. segment-CN ↔ read-level subclone 的「鴻溝」（決定能否用 CN 把 read 分群）

把 read 分到 subclone 需要的是 **read-level allele assignment**（哪條 read 帶哪個 allele/multiplicity/屬哪個 CCF 群），而非 segment-level CN。兩者之間有一道結構性鴻溝（L2）：

- **multiplicity 歧義**：`u_i = CCF_i × m_i` 不可單獨解（Tarabichi 2021）——知道 segment CN **不足以**定 SNV 的 multiplicity。
- **purity-ploidy 不可識別**：homozygous deletion@30% purity ≡ heterozygous deletion@60% purity（PyLOH）——須引入 **allele-specific 資訊（het-SNP BAF / LOH）**才能打破（這正支撐既有 LOH.bed + haplotag 方向）。
- **業界框架硬邊界**：subclonal reconstruction 工具**只把 mutations/segments 指派到 subclone，不指派 individual reads**（Principles 2017）；reads 本身不被指派。

> **InterSubMod 的真正增量定位（L3）**：填補「segment-level CN ↔ read-level subclone 標籤」的鴻溝。ISM 相對 bulk SNV 方法的增量 = **多 sSNV + 甲基 epiallele 共現於同一條長 read/分子**，提供 segment CN 給不了的額外約束（chr2:18M 即此原型）。**對外不可宣稱「把 read 分到 subclonal CN」為已解** = ONT 領域開放缺口。

---

## 4. 三軸獨立性 + 核心裁決（B4 對抗 reviewer = REFUTED as standalone validator）

### 三軸獨立性前置判斷（全案地基，L3）

| 軸 | 訊號來源 | 對 subclone 角色 | 與甲基軸獨立性 |
|---|---|---|---|
| **CN 軸**（Wakhan/SAVANA）| read depth + het-SNP BAF | segment-level 結構標籤 | 🔴 **非獨立** — CN=1(LOH) 直接改寫甲基 allelic 結構；透過「哪條 allele 還在」共同上游耦合 |
| **SV 軸**（Severus/SAVANA breakpoint）| split/supplementary alignment | read 二分標籤 | 🟢 **真獨立** — read 是否物理跨 junction 與其 CpG β 值無因果耦合 |
| **甲基軸**（ISM read×read 距離 + UPGMA + PERMANOVA）| per-read 5mC (MM/ML) | 待驗證 clustering | — |

### 為何「CN 分群 ↔ 甲基分群對齊」是循環的（L1 源碼自證 + 三條獨立證據）

**頭號威脅 = W2/I2 LOH-unmask-ASM**：
- **機制（L2）**：tumor 內 imprinted-DMR 甲基變化 **82–91% 由 CN/LOH 解釋**（含 copy-neutral LOH），非 epimutation（Martin-Trujillo PMID28883545；Chase r=0.76）。LOH 區只剩單 allele → 甲基「看似 monoallelic」是 **CN artifact 不是新 epimutation**。
- **為何致命（L1 源碼）**：親查 `src/core/SignificanceAnalyzer.cpp:307-340`，`verification_class="Subclone"` 判準 = `!label_sig && cluster_significant`（有 cluster 結構但與 somatic label 不相關）—— **此標籤本身不含任何 CN/cis 排除**；`potential_loh`（hp_ratio<0.1 或 >0.9）在 `Strong` 分支**只被標記不被降級**。把一條 CN 軸塞回去當 label，「LOH 分群 ↔ 甲基分群對齊」**機械必然為真**（I2）→ **零 subclone 資訊**。

**三條獨立證據鎖死「成功對齊反而是壞消息」（L1）**：
1. 跨 7 樣本 **Jaccard 0.091–0.161** —— 無監督切群與標籤對齊 PERMANOVA **大致不相交**（跨 3 癌種穩健）。
2. clear≥2 中 **85% 對齊 germline/carrier = cis-ASM**（非 subclone）。
3. within_clean ≠ subclone（Jaccard 0.123，55% epiallele）。

**W1 援引紀律（framing 風險，L1）**：region-level 甲基 `|r|<0.07`（CN-blind）是 **feature×CN residualized 相關**，**不是** read×read 距離分群對齊；兩者**不可互相援引**。但這對提案是**雙刃**：若 read-level 甲基對 CN 無感 → CN 無法驅動甲基分群（無相關基礎）；若有感並對齊 → 落入「甲基只是 CN 鏡子」。

---

## 5. SEQC2 truth 能驗什麼 / 不能驗什麼（B3 全 CONFIRMED，第一手 awk 重現 = L1）

> SEQC2 CNV truth **已在本機** `/big8_disk/data/HCC1395/SEQC2/CNV/`，數字皆 B3 第一手 `awk` 重算。

| 用途 | 可行？ | 第一手依據（tier）|
|---|---|---|
| L_CN **gain/loss/loh 方向**驗證（segment 層）| ✅ **GO** | `ngs_benchmark_cnvs_gain_loss_loh.bed` = **660 段（307 gain / 33 loss / 320 loh）**；bp 總和 gain 1503.9Mb / loss 87.9Mb / loh 1490.5Mb（L1）|
| **chr8 LOH 96% / chr2:18M LOH 區段**外部成立 | ✅ **GO** | chr8 LOH = **139.3Mb = 96.0%**（疊加 93.2Mb GAIN = amplified-on-LOH = W2 高發區）；chr2:18,086,020 落單一 LOH 段 16,146,119–22,100,000（L1）|
| **per-base integer CN** 驗證 | 🔴 **NO-GO** | 來源論文（Genome Biol 2024, PMC11188507）**逐字明文**：`"our confidence was only limited to qualification of CNV calls... should not be used as 'gold-standard' for benchmarking"` + `"break points... imprecise"`；median-CN BED 含**半整數**（3.5/4.5…11.5）（L1+L2）|
| **allele-specific CN** truth | 🔴 **NO-GO** | `Additional_file_5` VCF 1451 records，**FORMAT=GT 全為 `./.`**（1451/1451）；grep BAF/MCN/major/minor/ASCN **全空**（L1）|
| **read-level** subclone 驗證 | 🔴 **NO-GO** | 無 read/cell linkage + IMPRECISE breakpoint（L1+L2）|
| SV 正交確認（chr8 易位）| ⚠ **PROBE** | Talsania 2022（1777 SV / 5 平台 consensus）存在但**本機無 SV 檔**（find 空），需下載 `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/`（L2 medium）|

**關鍵副產品（L1）**：chr2:18,086,020 pivot 同時落 **High-Confidence-region GAP**（18,085,984 與 18,086,058 之間）→ SNV-truth **unevaluable**，與既有 chr2:18M 結論精確吻合。CNV truth 來源是 **Genome Biol 2024 PMC11188507**（6 caller × 21 WGS replicate × 3 orthogonal），**非** Fang2021(SNV/indel) / Talsania2022(SV)——三者 truth 分屬不同論文，引用勿混。

> **未解口徑歧異（L1，引用前必釘）**：LOH.bed Jaccard —— C13 卡片 + 20260409 報告記 **0.847**（⭐3 單樣本），prompt/memory 引用 **0.928**。屬 region-level vs bp-level 重疊基準的口徑歧異。

---

## 6. SV 斷點的角色 = joint resolve 的核心增量（A3，L2）

> **SV 軸是三軸中唯一對甲基 β 沒有因果耦合的標籤** = 不可替代價值。

**具體用法（GO）**：
- **SAVANA `{S}_sv_breakpoints_read_support.tsv`**（第二欄 tumour-supporting read IDs，第三欄 normal）+ **Severus `--output-read-ids`**（haplotype-specific，HP/PHASESET_ID）。
- 形式化：跨斷點 read → **「攜帶 SV / 不攜帶」二分群** → 當 a-priori label 餵 ISM 甲基距離 → PERMANOVA(`StructureTest`) + PERMDISP(`check_dispersion`) + CramérV/Fisher 對齊檢定（沿用既有 cluster×label 框架）。

**三大優勢**：
1. **非循環、不依賴甲基 β**（split-read 來自 CIGAR/supplementary，與 β 正交）。
2. **排 C1 的 LOH-unmask 假象**（不被 LOH 製造的 monoallelic 假象污染）。
3. **補 MethPhaser cross-block GAP**：MethPhaser（Nat Commun 2024, PMC11193733）明言它**缺乏物理跨越兩 block 的 read 連結確認力**；Severus/SAVANA 的 SV adjacency edge 提供——一條跨 junction split-read 同時 anchor 兩個遠端 loci。

**⚠ power 限制（C5）**：somatic SV 數**遠少於 SNV**，單 region 可能無 informative SV → SV-label 是「**有則為強 anchor**」非全基因組普遍可用，與 SNV/甲基**互補非取代**；且 low-CN allele 的 SV-supporting read 變少 → 用 SV anchor 控 C1 時須**同帶 Wakhan allele-CN** 一起判讀。

---

## 7. Confound 處置表（裁決成敗核心）

| # | Confound | 為何致命 | 處置 | 殘餘風險 |
|---|---|---|---|---|
| **C1** | W2/I2 LOH-unmask-ASM | CN 分群與甲基分群因同一 CN 事件對齊 → 對齊≠subclone | (a) normal cis-control；(b) CN-held-constant（只在 CN 相同 read 子集內分群）；(c) multi-sSNV CCF 獨立定 lineage | 🔴🔴 **最高**。LOH 區 normal cis-control 塌縮 can't-ask（39.5%/12039 位點），分 I2 鑰匙恰在最難取得處 |
| **C2** | W1 region-level CN-blind 是否外推 read-level | 若 read-level 也無訊號 → 路徑死 | **明確區分層次**：W1 是 feature 相關非 read×read 分群對齊；read-level 須獨立實測（chr2:18M ⭐3 已證特定位點**有**結構）| 🟡 framing 風險 |
| **C3** | 循環性（CN 用 depth+BAF，甲基用 5mC）| 若隱性耦合，對齊是同義反覆 | depth/BAF 與 5mC 表面獨立但透過「哪條 allele 存活」共同上游耦合 → **不是純循環，是 confound 耦合**；SV 軸是唯一無此耦合 | 🟠 中 |
| **C4** | alignment≠subclone（cis-ASM）| 「確認≥2 群」「對齊獨立軸」都≠confirmed subclone | 對齊後一律標 "aligned-not-attributed"；subclone 歸因須額外 multi-sSNV CCF / single-cell / COLO829 | 🟠 中 |
| **C5** | SV-read 計數被 CN 偏掉 | SV anchor 在 LOH 區可能假陰 | 用 SV anchor 控 C1 時同帶 Wakhan allele-CN | 🟡 中低 |
| **C6** | purity-ploidy 不可識別 + multiplicity 歧義 | segment CN 無法定 read CCF | allele-specific(LOH.bed+BAF) 打破歧義；read-level 多 sSNV 共現補不足 | 🟡 結構限制 |

---

## 8. HCC1395 Pilot 計畫（含 CN-neutral control = B4 指定主實驗）

> **核心設計原則（B4）**：把「甲基對齊 CN 軸」當 **NULL/壞消息假設**，不當 subclone 證據。只有 **CN-neutral 區仍對齊 somatic label** 才是非循環候選。

### Pilot-1：chr2:18M（已 ⭐3，三軸交叉首選，含 CN-neutral 對照位點）
- **Region**：chr2:16,146,119–22,100,000（SEQC2 單一 LOH 段，6 sSNV）。
- **資料流**：Wakhan 取此段 allele-CN → Severus/SAVANA 找此窗 somatic SV breakpoint → ISM 甲基 clustering → 三軸對齊。
- **預期結果（既有 ⭐3，L1）**：骨幹獲證（α/β allele 互斥 0 違反、LOH HP1 1.4%、normal=REF）；**~3 lineage**（非「5 群」；pos4 20bp poly-T homopolymer 假象已知）；pos **18,086,020 落 SNV-HC 空隙 → unevaluable**；多 sSNV CCF 梯度 **0.10→0.81** → 甲基 cluster 對到 CCF 階序（**非循環的 lineage 確認，不靠 CN 對齊**）。
- 🔴 因落 LOH 段，甲基 monoallelic **必先掛 C1 警示**，lineage 歸因走 multi-sSNV CCF。

### Pilot-2：chr8 全 LOH 高密度 40kb 子窗（C1 confound 壓力測試，**非** CN-neutral 對照）
- **Region**：chr8 子窗 **133.68M / 144.84M / 132.52M**。
- **目的**：C1 對抗測試——若甲基在此「分群」，先驗證是否 CN artifact。
- **預期結果（L1+L3）**：normal cis-control 塌縮 can't-ask → I2 無法形式排除 → 天花板 = undecidable；**正確輸出 = 標 "LOH-confounded, not subclone"，不升級**。
- ⚠ chr8 全 LOH 96% 是 **7.4× FP 雙面刃**，整段不可用作 CN-neutral 對照。

### 關鍵區分實驗（主實驗）
分層成 **(i) CN-neutral het-retained 區（2 haplotypes 存活）** vs **(ii) LOH 區**。
- 在 (i)：CN 無變異 → 無法用 CN 解釋甲基分群 → 若甲基 read×read clustering **仍與 somatic ALT-carrier 標籤對齊**（CramérV≥0.7 + Fisher + PERMDISP 排 dispersion 假象）→ 這才是「甲基帶獨立 subclone/cis 訊號」的**非循環證據**。
- 在 (ii)：對齊**一律先歸 I2**，不計入 subclone。
- 進一步分 cis-ASM(I3) vs subclone(I4)：CN-neutral 區對齊後，仍需 multi-sSNV CCF gradient（chr2:18M 法）證甲基 pattern 跟著 ALT-carrier 的 CCF 梯度（subclone）而非單一 cis 位點。

---

## 9. 逐子問題 GO/PROBE/NO-GO 總表

| 子問題 | 裁決 | 理由 |
|---|---|---|
| 最佳 ONT CNV 工具 | ✅ **GO** | SAVANA + Wakhan |
| 能否每位置整數 CN | ✅ per-segment / 🔴 per-base | 分段常數 |
| 能否用 CN 直接把 read 分群 | 🔴 **NO-GO 直接** / PROBE（haplotype-CN 經 HP tag 間接）| multiplicity 歧義 + purity-ploidy 不可識別 |
| 與既有 label 關係 | ✅ **GO（互補）** | L_CN HP1/HP2 正交補強 somatic haplotag；tumor-only 無監督軸**已 NEGATIVE 勿再開** |
| CN-derived 分群 ↔ 甲基對齊驗 subclone | 🔴 **REFUTED standalone / PROBE under strict gate** | C1 LOH-unmask 循環耦合 |
| CN+SV+甲基 joint resolve | ⚠ **PROBE（SV 軸為核心增量）** | 各軸切不同混淆，非投票 |
| SV 斷點當標籤 | ✅ **GO（有 informative SV 時為強 anchor）** | 唯一非循環 anchor |
| SEQC2 驗證 | ✅ segment 方向 / 🔴 per-base/allele/read-level | 定性 truth + imprecise breakpoint |

---

## 10. 對外撰寫紅線 + 落地前必補清單

**三條對外撰寫紅線**：
1. **「對齊」≠「subclone 確認」**——只有 SV 軸對齊 + multi-sSNV CCF 是非循環。
2. **「把 read 分到 subclonal CN」是 ONT 開放缺口**，非已解；甲基在此框架是 characterize（有界），reconstruction 由 somatic haplotag 骨幹負責。
3. **單樣本/單 pipeline 封頂 ⭐3**；cross-platform 復現需 COLO829 + 外部工具 + SEQC2。

**落地前必補（L1/L3）**：
1. **Wakhan/SAVANA 在 HCC1395 + COLO829 先跑**，與內部 Coverage_Multiple(r~0.83) + LOH.bed(Jaccard 待釘 0.847 vs 0.928) 對齊驗證後，再寫入 KB `05_tools`（本地 KB **未索引** Wakhan/SAVANA）。
2. 🔴 **驗證對齊必在 `feat/summary-nreadsvalid` branch 跑**（含 `strength_tumor` 去 double-dip + `SubcloneDbeta`/`compute_group_dbeta_test`，commit 2d5692f 等）——當前 worktree `research/subclonal-reconstruction-202606`（HEAD f8212f6）grep 零命中，跑出的「Subclone」是**舊 4 類裁決**（有 double-dip 缺陷）。
3. **Wakhan 命名 H1937/H1954 vs HCC1937/HCC1954 逐一核 ATCC ID**。
4. **SV 正交確認需從 ftp-trace 下載 Talsania2022 1777 SV**（本機無）。
5. **pre-registration 寫死**：「若甲基高度對齊 CN/LOH → 預設裁決 = I2 LOH-unmask / CN-mirror，非 subclone」，否則重蹈 tumor-only 軸 NEGATIVE 的 double-dip 覆轍。

---

## 11. Provenance / 數據來源 / tier

- **產出**：Dynamic Workflow `wf_5497277c-620`（10 agents：5 調研 + 1 方法設計 + 3 對抗驗證 + 1 綜整；943K subagent tokens；174 tool calls；698s）。完整 raw 輸出在 session task 檔 `w9bkhvtp7.output`。
- **L1（第一手重現）**：SEQC2 CNV truth 660 段/307·33·320、chr8 96%/139.3Mb、GT=./.（1451/1451）、SignificanceAnalyzer.cpp:307-340 裁決邏輯、Coverage_Multiple C++ 位置、跨 7 樣本 Jaccard 0.091-0.161、85% cis-ASM、binary 在 feat/summary-nreadsvalid。
- **L2（文獻/README 二手）**：Wakhan/SAVANA/Severus 能力（GitHub README + 搜尋摘要；medRxiv/PMC 全文 403/reCAPTCHA 未一手核實）、Tarabichi/Principles multiplicity、CAMDAC、MethPhaser、Martin-Trujillo 82-91%。
- **L3（推論）**：三軸獨立性判斷、ISM 增量定位、命名核對。
- **未解項**：LOH.bed Jaccard 0.847 vs 0.928 口徑、r=0.827 來源、Saraband（查無此工具，疑名稱有誤需用戶補來源）、SV truth 1777 細分（未一手）。

---

## 12. 下一步（需用戶 go-ahead；含 compute + 工具安裝）

1. **裝 SAVANA + Wakhan + Severus**（conda），在 HCC1395 tagged BAM（`big7_disk_output/canonical/HCC1395/paired_full/`）跑 tumor-normal，產 segment integer/allele/haplotype CN + SV read_support.tsv。
2. 與 SEQC2 `gain_loss_loh.bed` 做 **segment 方向 concordance**（合法）+ 與內部 Coverage_Multiple/LOH.bed 對齊。
3. 切 `feat/summary-nreadsvalid` worktree，在 **chr2:18M + chr8 40kb 子窗 + CN-neutral 對照區**跑「關鍵區分實驗」（§8）。
4. 把 SV-read 二分群當 a-priori label 餵 ISM 甲基距離（§6）。
