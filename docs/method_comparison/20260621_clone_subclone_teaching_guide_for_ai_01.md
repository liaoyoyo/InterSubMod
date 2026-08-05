<!--
建立時間: 2026-06-21
報告類型: 教學指南（meta-teaching guide）— 給「其他 AI」用來教學人類用戶理解 clone/subclone 重建 × ISM 主軸
任務類型: D handoff — 跨 AI 教學交接
build_branch: research/subclonal-reconstruction-202606
status: teaching scaffold（搭配知識紀錄主 doc 使用）
data_sources:
  - InterSubMod/docs/method_comparison/20260621_clone_subclone_reconstruction_landscape_and_ism_feasibility_01.md（內容真值來源）
  - InterSubMod/docs/method_comparison/20260620_somatic_locus_methylation_combination_enumeration_01.md
  - InterSubMod/docs/methodology/20260617_structure_label_situation_inventory_01.md
provenance_note: 本檔為教學「如何教」的 scaffold，不重述數字真值；所有事實以上列主 doc 為準。教學 AI 引用數字時須回主 doc grep。
-->

# 教學指南（給其他 AI）：如何教用戶理解 clone/subclone 重建 × ISM 主軸

> **給你（教學 AI）的話**：這份不是內容本身，是**「怎麼教」的腳本**。內容真值在主 doc（見 frontmatter）。你的任務 = 用下面的順序、類比、checkpoint、紅線，把用戶從「知道有 clone/subclone」帶到「理解 ISM 能做什麼、不能做什麼、為什麼」。**先讀主 doc 再教**；每個數字都回主 doc grep，不憑記憶。

---

## 0. 學習者背景（先校準）

- **誰**：碩士生，正在寫論文「Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing」。技術夠（懂 phasing/methylation/somatic calling），要的是**把全領域與自己方法的關係、邊界、可行性想清楚**。
- **已知**：ISM pipeline 細節、自己跑出的 NEGATIVE/POSITIVE 結果、external_validation 69 源庫。
- **想要**：(a) 不踩沒考慮到的盲點 (b) 知道哪些是共識可引用、哪些未知 (c) 確認論文目標站得住 (d) 能對外防守（reviewer）。
- **語言**：繁體中文；技術詞保留英文。
- **🔴 用戶最在意「誠實」**：寧可說「未定/做不到」，不要 overclaim。教學時**主動標不確定性**，不要為了讓他開心而誇大。

---

## 1. 概念依賴順序（先教 A 才教得懂 B）

```
L0 clone vs subclone 定義（CCF~1 truncal vs 中低 CCF branch）
   ↓
L1 怎麼量：VAF → CCF = f·multiplicity/(purity·ploidy)  ← multiplicity 是關鍵未知
   ↓
L2 為何只用 SNV 難：non-identifiability（tree/cluster/evolution）+ multiplicity 歧義
   ↓
L3 加 LOH+CNV 還是難：purity-ploidy 不可識別 / CN-known 仍不足定 m
   ↓
L4 甲基的「虛實」：single-cell lineage REAL / bulk read-level 未建立 / filter DEAD
   ↓
L5 觀測單元根本差異：single-cell=每細胞全基因組條碼 vs bulk long-read=單分子局部、無跨 locus 細胞連結
   ↓
L6 ISM 的 niche：在 genetic 骨幹上、read 層次、區域內、多證據、甲基當佐證（proof-of-concept）
   ↓
L7 可行性裁決：可完成 as framed（⭐3）/ 不可宣稱 genome-wide 甲基驅動重建 / gate HD-1·G-B
```

**不要跳級**：用戶若還沒接受 L5（觀測單元差異），就講 L6（ISM niche）會聽成「ISM 也能像 single-cell 一樣重建」→ overclaim 誤解。**L5 是整個理解的樞紐**。

---

## 2. 必須主動破除的迷思（教學 AI 要預先攔截）

| 迷思 | 為何錯 | 正確說法 |
|---|---|---|
| 「甲基切出 ≥2 群 = 找到 subclone」 | read 內甲基相關（epiallele）把 1 clone 切成 3-4 群；85% 清楚切群對齊 germline/carrier = cis-ASM | 確認 ≥2 群 ≠ 確認 subclone；需 a-priori 對齊 + normal cis-control |
| 「HP1-1 半帶 ALT = subclone」 | 半帶對 multiplicity formally non-identifiable（=1 / CCF<1 / ADO）| 半帶是訊號**之一**，不能單獨宣稱 subclone |
| 「甲基能驅動重建」 | bulk read-level 甲基 subclone 未建立（G-B undetermined）；甲基=germline-haplotype 層 | 重建由 **somatic haplotagging（genetic）**驅動；甲基 = corroborate |
| 「tumor-only 無監督聚類能偵測 subclone」 | double-dip 循環；noise≈structure 83-100% | 只有 a-priori-conditioned 軸合法 |
| 「ISM 比 cvlr/ASMS 強因為對手用二代定序/缺檢定」 | 🔴 cvlr/ASMS/MethylBERT 全 ONT-capable 且全有 randomization 檢定 | 差異 = **檢定對象 + normal 錨點 + target 的整合**，不是定序平台或有無檢定 |
| 「ISM 做 subclonal reconstruction（= 像 PhyloWGS 建樹）」 | 那是學派 1-2 的正式定義，ISM 單 bulk 樣本欠定 | ISM 的 "reconstruction" = **regional LOH-constrained same-haplotype partition**，非 tree |
| 「ISM 填補了沒人做過的空白」 | 每個元件都 prior art；2025-26 競品正填補 | **先驗元件的新穎整合**於尚空交集，非無人嘗試的真空 |

---

## 3. 有效的類比（用戶用這些「秒懂」）

- **觀測單元（最重要）**：single-cell = 「每個細胞拿到一張**完整基因組的條碼卡**，跨位點天然連在同一張卡上」；bulk long-read = 「一地散落的**單張碎紙片（read）**，每張只記局部、**無法知道哪幾張來自同一個細胞**」。→ 所以 single-cell 能拼全基因組家譜，bulk 只能在「碎片夠長能跨幾個位點」的**區域內**講故事。
- **phase set / haplotype block**：像「**分段沒接起來的拼圖**」——每段內部拼得對（block 內 HP1/HP2 一致），但**段與段之間哪邊是哪邊接不起來**（block A 的 HP1 ≠ block B 的 HP1）。
- **LOH-unmask**：像「**雙聲道音樂關掉一邊聲道**」——你聽到的「特別」不是新錄的音軌，是本來就有的另一聲道**因為對側被關掉而顯出來**（既存 ASM，非新 somatic）。
- **over-clustering**：像「**把一團均勻的霧硬分成三塊**」——分得出界線不代表真有三塊。
- **multiplicity 歧義**：像「**看到一半的人舉手，分不出是『一種人各舉一隻手』還是『兩種人，一種全舉一種全不舉』**」。
- **ISM 的 niche**：像「**不是發明新樂器，是第一個把三件現成樂器（結構檢定 + normal 對照 cis-test + subclone 目標）合奏在 bulk ONT 的 read 層次上**」。

---

## 4. 建議教學順序 + Socratic checkpoint（每階段一個驗證理解的問題）

| 階段 | 教什麼 | Checkpoint 問題（用戶答對才往下）|
|---|---|---|
| 1 | clone/subclone + CCF + multiplicity | 「為什麼看到 VAF=0.25 不能直接說是 25% 細胞帶的 subclone？」（答：要除 multiplicity 與 purity/CN）|
| 2 | SNV-only 為何難 | 「給同一組 VAF，為什麼可能有多棵不同的演化樹都符合？」（答：non-identifiability）|
| 3 | 加 LOH+CNV 仍難 | 「知道了 allele-specific CN，multiplicity 就定了嗎？」（答：不，DeCiFer 證仍不足）|
| 4 | 甲基虛實 | 「為什麼 single-cell 甲基能重建家譜，bulk 卻不行？」（答：觀測單元 = cell vs read）|
| 5 | ISM niche | 「ISM 的『重建』和 PhyloWGS 的『重建』是同一件事嗎？」（答：否；regional partition vs genome-wide tree）|
| 6 | 可行性裁決 | 「論文目標能完成嗎？在什麼條件下？」（答：能，as proof-of-concept，骨幹=haplotag、甲基=corroborate、⭐3、gate HD-1/G-B）|

---

## 5. 🔴 教學紅線（不可教成定論的東西）

1. **不可**把「ISM 用甲基重建 subclone」教成 ISM 的成果 → 甲基是 corroborate，骨幹是 genetic。
2. **不可**把 chr2:18M 個案教成「已證實的 subclone 演化樹」→ 是 hypothetical / consistent-with-data，⭐3，18086020 timing=假說。
3. **不可**教「tumor-only 無監督能偵測 subclone」→ 已 NEGATIVE。
4. **不可**教「確認 ≥2 群 = subclone」→ 85% 是 cis-ASM。
5. **不可**教「ISM 勝過 single-cell / SNV-CCF」→ 不同軸，別人在各自軸更強。
6. **必標**：single-sample / single-pipeline / by-construction circularity（HD-1）/ G-B undetermined。
7. **外部數字必標 L3**（投稿前過 /citation-verification）；ISM 數字回主 doc grep，不憑記憶。

---

## 6. 教完的「成功標準」（用戶理解到位的訊號）

用戶能用自己的話講出：
- (a) 為何純 DNA subclone 重建欠定（multiplicity + non-identifiability）；
- (b) 為何甲基在 single-cell 真實、bulk read-level 未建立（觀測單元差異）；
- (c) ISM 的 "reconstruction" 是 regional partition 不是 genome-wide tree；
- (d) ISM 的可防守創新 = **三件式整合 on bulk ONT**，不是「別人做不到的能力」；
- (e) 論文**可完成 as proof-of-concept**，但受 HD-1/G-B gate 限制、甲基是 corroborate。

---

## 7. 權威來源指標（教學 AI 要 grep 的地方）

| 要教的點 | 去哪查真值 |
|---|---|
| 全景 / 裁決 / 難點共識 | 主 doc `20260621_clone_subclone_reconstruction_landscape_and_ism_feasibility_01.md`（§1-§11）|
| 組合判讀（LOH/half-ALT/k群）| `20260620_somatic_locus_methylation_combination_enumeration_01.md` |
| 內部數字 SoT | `20260617_structure_label_situation_inventory_01.md`（verify 25/25 PASS）|
| why-hard 缺口 GAP-A~L | `20260619_subclone_analysis_interpretation_full_framework_01.md` |
| 外部 69 源庫 | `/big7_disk/liaoyoyo2001/external_validation/_landscape/{05,08}.md` |
| memory 快速 recall | `project_clone_subclone_landscape_and_ism_feasibility`、`project_somatic_locus_methylation_combination_interpretation` |

> **最後叮嚀（給教學 AI）**：用戶要的是**有方向、不漏盲點、能對比驗證**。教學時多用「這點是共識可引用 / 這點還未定 / 這點我們曾誤判已更正」的標籤，讓他建立**可信度地圖**而非單純記結論。
