<!--
建立時間: 2026-06-09
狀態: focus (ISM 工具方法細節 + 應用意義 — 餵 Methods §4.2 + benchmark 任務 #23 + tooling 定位)
報告類型: paper_focus_tool_method_significance
受眾: 廖子游 · PI · 執行 benchmark 的其他 session
provenance_note: 方法描述為本工具設計(🔵)；與外部工具差異一律標「待 benchmark/KB 確認」(P-14 outside-claim 勿臆測)；內部數字沿用已驗證集合。
-->
<!-- provenance-verified: ISM 方法描述為自有設計(🔵)；外部工具能力宣稱標「待 #23 benchmark + KB 確認」未臆測；目標 venue Genome Biology。 -->

# ISM 工具方法細節 + 應用意義

> **L0 一眼結論**：ISM 的價值**不在單一新演算法**，而在**一條 read-level 的「去 confound 鏡頭」** —— 把 somatic 子單倍型(HP1-1) 解析度 + reliability gating + **normal-anchored cis-test** 串起來，回答「這個 allelic 甲基差異是**真 cis**、還是 copy / germline-allelic / drift 假象」。**應用意義 = 誠實的 de-confounding + characterization，不是 discovery/filter。**
>
> **L1 重點邏輯**：
> ① **方法 = 7 個組件的整合**（read×CpG → 距離 → 分群 → reliability gate → HP/LOH 比對 → normal-anchored cis-test → catalog）。
> ② **差異化在「組合 + somatic-sub-haplotype + normal-anchoring」**，個別演算法多為既有 → GB 定位為「useful pipeline / de-confounding lens」非「new algorithm」（誠實，避免 novelty overclaim）。
> ③ **應用意義 = de-confound**：文獻（Martin-Trujillo）證 apparent ASM 多被 copy confound，多數工具不去 confound → ISM 的 normal-anchored cis-test 是真實 niche。
> ④ **誠實限制**：normal-anchored 只在 **有 matched normal**（目前只 HCC1395）才完整；single-pipeline ⭐3。

---

## A. 方法細節（7 組件 → 對映 Methods §4.2）

| # | 組件 | 做什麼 | 輸入 → 輸出 | 關鍵參數/檔 |
|---|------|--------|-----------|------------|
| 1 | **read×CpG 甲基矩陣** | 每 region 建 reads × CpG sites 的甲基狀態矩陣（5mC/5hmC 分軌）| BAM MM/ML → matrix | MatrixBuilder；⚠ 5mC/5hmC dup-bug 待修 |
| 2 | **read-read 距離矩陣** | reads 間甲基 profile 距離（6 度量）| matrix → distance | NHD/L1/L2/Pearson/Bernoulli/Jaccard；DistanceMatrix |
| 3 | **階層分群** | 依甲基相似度把 reads 分群 | distance → clusters | hierarchical |
| 4 | **CramersV reliability gate** | 判分群是否統計可靠（非稀疏表 artifact）+ PERMANOVA 顯著性 | clusters → reliable/latent/none | Cochran 最小期望格≥5；RegionProcessor.cpp:1592；487 latent(ISM-2 待 audit) |
| 5 | **HP/LOH tag 比對** | 把 haplotype（HP1/HP2/HP1-1/HP2-1）+ LOH 疊到甲基分群上 | clusters + tags → 對齊關係 | somatic-sub-haplotype 解析度 |
| 6 | **normal-anchored cis-test** ⭐ | normal HP1 / tumor HP1 / tumor HP1-1 三路比 → 分真 cis vs drift vs copy | tumor+normal → cis_status | **核心去 confound 步驟** |
| 7 | **characterization 輸出** | per-locus 統計 → catalog（16 欄 7 TAG）| 上述 → catalog.tsv + 分佈圖 | ARI vs imprinting ruler；Δβ(ISM-3 待併) |

> **方法層級的洞**：甲基分群是 **germline-haplotype 層級**（分不同 haplotype 強、同 haplotype 內 somatic 亞群弱）—— 這解釋為何 ISM 能 characterize haplotype 結構，但「同 haplotype 內 subclone 甲基」訊號弱（D6 NEG）。

---

## B. 與現有工具的差異（🔴 全部「待 #23 benchmark + KB 確認」— 不臆測）

> **P-14**：以下是**待驗證的差異化假設**，不是已證宣稱。#23 benchmark 先查 `/big8 Knowledge/05_tools` + `mcp__knowledge__*` 確認各工具實際能力，再 head-to-head。

| 工具 | 一般定位（待 KB 確認）| ISM 假設的差異（待 benchmark 證） |
|------|---------------------|-------------------------------|
| **modkit** | per-position 甲基 calling / pileup / DMR | ISM 是 **read-level 分群**非 per-position 聚合（待證）|
| **NanoMethPhase** | haplotype-resolved 甲基（HP1 vs HP2 ASM）| ISM 多 **somatic-sub-haplotype(HP1-1)** + read-read 距離 + normal-anchored 三路 cis（待證）|
| **nanomethviz** | 甲基視覺化 | ISM 是分析+統計 gate 非純視覺（待證）|
| **(epipolymorphism 類)** | read-level 甲基熵/異質性 | ISM 多 **HP/LOH 對齊 + cis de-confound**（待證）|

→ **GB 賣點（若 benchmark 成立）**：「read-level read-read 距離 + 階層分群 + reliability gate + somatic-sub-haplotype 解析 + normal-anchored cis de-confound 的**整合**，是現有單一工具做不到的。」

---

## C. 應用意義（能回答什麼 + 實際用途 + 誠實限制）

**ISM 能回答、且有人在乎的問題**：
1. **「這個 allelic 甲基差異是真 cis 還是假象？」** ← 最有價值。normal-anchored cis-test 把 copy / germline-allelic / drift 從 cis 分開（chr17 過、BRCA2 被判 copy）。
2. **「某位點 reads 是否分甲基亞群、且與 haplotype/somatic 結構對齊？」** ← read-level characterization。
3. **「哪些位點甲基分群統計可靠 vs 稀疏表 artifact？」** ← CramersV/PERMANOVA gate（487 latent 發現）。

**實際用途**：
- 建**位點甲基分群 catalog**（交付物）。
- 當**de-confounding lens**：文獻（Martin-Trujillo, CN 解釋 allelic 甲基 82–92%）證 apparent ASM 多被 copy confound，**多數工具不去 confound** → ISM 補這個。
- read-level 腫瘤甲基異質性 characterization。

**🔴 誠實限制（GB reviewer 會抓，先寫）**：
- **不做 variant filter**（已證死四道）—— 是 characterization/de-confound 工具，非 discovery/filter。
- **normal-anchored 只在有 matched normal 才完整** → 目前只 HCC1395（其他 5 樣本無法跑完整 cis-test）= 應用範圍窄。
- **single-pipeline ClairS-TO ⭐3**；5mC/5hmC 分軌待 dup-bug 修。
- **演算法多為既有組裝**（見 D）→ 定位 useful pipeline 非 new algorithm。

---

## D. 演算法 novelty 自評（誠實 — 避免 GB novelty overclaim）

| 組件 | novelty |
|------|---------|
| read×CpG 矩陣 / read-read 距離 / 階層分群 | 既有（epipolymorphism/methylation entropy 類已有）|
| CramersV/PERMANOVA gate | 標準統計，**應用於 reliability gating 的紀律**是貢獻 |
| somatic-sub-haplotype(HP1-1) 解析 | 較少見（依賴 longphase-S somatic tag）|
| **normal-anchored cis-test（de-confound）** | **最接近原創貢獻** —— somatic-controlled HP-axis 設計 held-constant copy/ploidy/alignment |

→ **誠實定位**：ISM = 「**整合既有 read-level 方法 + 一個 somatic-controlled de-confounding cis-test**」的 pipeline。GB 接受 useful pipeline，但 **Methods/Discussion 不可宣稱「新演算法」**，賣點放**整合 + de-confound 紀律 + 誠實 characterization**。

---

## E. 對映
- **Methods §4.2**（`論文架構_正式學術版_Slide2Thesis格式.md`）：本檔 A 節 = 方法描述來源。
- **任務 #23 benchmark**：B 節差異假設 = benchmark 要證的；先查 KB。
- **tooling 定位**：C 應用意義 + D novelty 自評 = Discussion §3.2/§3.3 + tooling claim 的誠實邊界。
- **catalog（#15）**：組件 7 = catalog 輸出。
