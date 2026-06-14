<!--
建立時間: 2026-06-13
報告類型: external_validation 集中庫 → InterSubMod 橋接索引（discoverability bridge）
任務類型: D handoff — 讓 ISM 開發 session 與論文撰寫 session 都能找到 + 使用 repo 外的外部驗證庫
狀態: pointer/index（指向 repo 外實體庫；本檔不複製內容，只導航 + 摘要 paper-critical 結論）
data_sources:
  - /big7_disk/liaoyoyo2001/external_validation/REGISTRY.md
  - /big7_disk/liaoyoyo2001/external_validation/_landscape/05_context_cards_index.md (59 源索引)
  - /big7_disk/liaoyoyo2001/external_validation/_landscape/08_paper_literature_map.md (整體論文文獻地圖)
  - /big7_disk/liaoyoyo2001/external_validation/_landscape/06_source_read_corrections.md (源碼校正+佐證)
provenance_note: 摘要結論皆引自 external_validation 各源 CONTEXT 卡(一手親讀 repo 源碼/全文 PDF)；數字溯源見各卡 §7。
-->

# 外部驗證庫橋接索引 — external_validation（repo 外實體庫）

> **這是什麼**：把 ISM 論文相關的**外部論文/程式碼實體下載 + 親讀驗證庫**接進 InterSubMod。庫在 **repo 外**（git 不追蹤），本檔是 repo 內的導航 + paper-critical 結論摘要。
> **庫位置**：`/big7_disk/liaoyoyo2001/external_validation/`
> **入口**：`REGISTRY.md`（主索引）→ `_landscape/05_context_cards_index.md`（59 源卡索引）+ `_landscape/08_paper_literature_map.md`（整體論文文獻地圖, Ch1-Ch5 一站式導航）+ `_landscape/06_source_read_corrections.md`（源碼校正與佐證）+ `_landscape/07_cis_basis_synthesis.md`（cis-test 科學依據）。
> **單源細節**：`axis{1-6}_*/＜slug＞/CONTEXT.md`（每源完整驗證紀錄：可信度/scope/主軸/論點/方法對照/共識⇄矛盾/操作差別/結論差別/親讀 code 對照）。
> **memory**：`project_external_validation_library`（其他 session 啟動可 recall）。

## 庫狀態（2026-06-14）
- **59 源**；tier **A×13 · B×31 · C×15**；39 repo clone + 36 PDF（2026-06-14 磁碟實點）。**57/59 已親讀**（primary/source/fulltext）；2 cited_secondary（deamination + mcf7-subline，皆無可下載 PDF/preprint）。
- **2026-06-15 cis-basis 全升級**：5 張 cis-basis 卡親讀 PDF 升 fulltext_verified（Onuchic2018/Min-GoDMC2021/Turcan2018/Shoemaker2010/Do2020，4 篇 C→B + Do2020 C→B）。fulltext 親讀同步抓出並修正 3 處原 web-二手瑕疵：① Onuchic 卡捏造的方法名「NPD」（PDF 無此詞 → 改 Fisher's exact + Shannon entropy）+ 誤稱「無 repo」（實有 github.com/BRL-BCM/allelic_epigenome）；② Min/GoDMC 卡 heritability 55%/12% 分層數字 PDF 查無（web-二手）→ 移除；③ Turcan 卡 citation 錯（卡寫 2012 Nature，磁碟 PDF 實為 2018 Nat Genet follow-up）→ primary_citation 改 2018、2012 降 related。每張卡均經對抗複查逐頁碼核對（含主 agent 圖面親驗 Turcan Fig 1c）。
- **整體準確性對抗稽核（20 卡，wmnyd7ytw）= CLEAN**（0 MAJOR，15 MINOR「更誠實」修正已折入）；末次補卡 = CpelNano（axis5, B — ISM structure-vs-disorder 最近 ONT 對手，最後一個真 P0 gap）+ methylation-to-mutation-deamination（axis3, C — 甲基→突變 SBS1 反方向）。
- 既有 **82-工具文獻 survey** 在 `InterSubMod/docs/method_comparison/20260609_ism_vs_external_methylation_tools/`（modkit/DSS/NanoMethPhase/MHB/Metheor… 文獻層，未建獨立卡）。
- 🔴 **文獻缺口稽核（2026-06-13，論文 framing 凍結前必讀）**：`InterSubMod/docs/reports/research_landscape/20260612_external_validation_literature_gap_audit_01.md` — 裁決：①「subclonal reconstruction」用詞偏強（正式 reconstruction 須 clone 數+lineage tree+CN/purity deconv+正交 ground truth；Sgootr/Gaiti/CAMDAC/Lee/BitPhylogeny 為壓力鏈）②「距離式/read-level 甲基」非天然新穎 ③6 cell line=跨模型 reproducibility 非病人 cohort ④「0 conflict」改「無同 regime 衝突但有 framing tensions」。據此 2026-06-13 補 P0+P1 ~21 源（庫→49）→ 2026-06-14 再補 cis-basis + CpelNano + deamination 收口至 **59 源、稽核 CLEAN**。

## 🔴 paper-critical 結論（論文 Related Works / Discussion / 投稿防守必讀）

### 1. 與我方 0 真 CONFLICT
最尖對比物（Epiclomal/MethylTree/EVOFLUx）皆 regime 差（單細胞 sc-WGBS / bulk-clock vs 我方 bulk 長讀 read-level）；EVOFLUx 自承 1,610/1,976 effectively-neutral **反而印證我方弱-subclone 紅線**。

### 2. 🔴🔴 ISM 創新點口徑（投稿必守，避 3 個會被打臉的錯誤差異）
與 read-level 甲基工具（cvlr / ASMS / MethylBERT）的差異 **不是** 下列任一，**禁止**這樣寫：
- ❌「對手用二代定序」— cvlr/ASMS/MethylBERT **都 ONT-capable**（MethylBERT 源碼 `cli.py:94 -m dorado` 讀 MM/ML）。
- ❌「對手缺顯著性檢定」— cvlr（Table1 median|Δmeth| vs 1000 random）/ ASMS（pval.rs 5000-perm）**都有** randomization 檢定。

**✅ 真正創新點（三件）**：① **無監督 read×read 距離矩陣的*結構*檢定（PERMANOVA 999-perm）** ② **normal-baseline somatic cis-test**（germline 工具皆無 matched-normal 錨點）③ **somatic-subclone characterization 目標**（它們是 germline ASM/imprinting 發現，或 supervised tumour-vs-normal deconvolution）。

### 3. 白地（cancer 甲基-phasing 開放）source+全文雙重證實
cancer phaser（Wakhan/HiCancer/LongPhase-S somatic 路徑）+ SR 骨幹（dpclust/sciclone/phylowgs/canopy/clip）+ somatic caller（DeepSomatic grep=0）**跨工具一致零甲基**；MethPhaser 明文留 cancer 為 future；**LongHap（2026）= germline-only、不威脅白地**（Fig3C 強不對稱甲基反佐證 R2）。

### 4. 多源獨立佐證 R2（甲基=germline-haplotype 層級強）
ASMS Table1（16/20 ICR 顯著 IGF1R p=2e-25）· cvlr Table1（8 imprinted gene Δm P<10⁻³）· LongHap Fig3C（germline-haplotype 強不對稱甲基）· Sakamoto/O'Neill（SNP-phase→per-hap 甲基 aDMR）。

### 5. 🔬 ISM normal-anchored somatic cis-test 的科學依據（2026-06-14；詳 `external_validation/_landscape/07_cis_basis_synthesis.md`）
**三前提皆有文獻基礎，分層 germline/somatic × cis/trans**：
- **前提① 變異→鄰近甲基(cis)真實普遍** = CONFIRM（**germline-cis 鐵錨**）：Min2021/GoDMC cis-mQTL 248,607≫trans 23,117（45.2% 位點帶遺傳關聯）+ Onuchic2018 sequence-dependent ASM → 也是 R2 機制背書；**但 germline-cis 正是我方 normal-anchor 要扣的基線**。
- **前提② somatic 變異→cis 甲基** = 有前例但**稀有**：Do2020「rare somatic mutations…track with de novo ASM」「does not make a large numerical contribution」+ Zhang2019 SV-breakpoint cis → **我方乾淨 somatic cis 稀有（1/816, chr17/TBC1D16）與文獻一致，非弱點**。
- **前提③ CpG-SNP confound** = CONFIRM：Shoemaker2010「38–88% germline ASM 由 CpG-SNP」= 我方 limitation 錨。
- **甲基↔somatic 突變總答**：YES，但主通道 = modifier 突變（IDH/TET2/DNMT3A, Turcan2012）的 **trans/genome-wide** 重塑，**非** 我方 cis-test 量的 locus-local cis。
- **Discussion 可用句（07 §L3）**：「Cis genetic effects on methylation are pervasive at germline level (GoDMC); our normal-anchored somatic cis-test subtracts this baseline to isolate the somatic-controlled increment; consistent with somatic-locus cis being rare (Do2020), clean somatic cis is infrequent (1/816).」

## 論文第二章（文獻探討）逐軸對應 → 外部驗證卡
| Ch2 次節 | 對應軸 / 卡（slug）|
|---|---|
| 2.1 亞克隆重建骨幹 | axis1: tarabichi-guide · salcedo-DREAM · pyclone(-vi) · phylowgs · dpclust · sciclone · canopy · clip |
| 2.2 ONT haplotype-resolved methylation / phasing | axis2: methphaser · longphase-s · sakamoto · oneill · longhap · wakhan · hicancer |
| 2.3 癌症 ASM | KB（Do2020/Martin-Trujillo）+ axis5: cvlr · asms · qfdrp |
| 2.4 甲基追演化 / lineage | axis4: evoflux · gabbutt-fcpg · epiclomal · methyltree · capper2018 · smallwood · fu-review |
| 2.5 甲基分群/距離地景 | axis5 + 82-工具 survey（method_comparison）|
| (Methods 兩-caller / TP-FP) | axis6: deepsomatic-castle（CASTLE 含 HCC1395+COLO829 外部真值）|

## 待補（非阻塞）
- 5 篇原源 PDF（dpclust/roth-pyclone/tarabichi/clip/epiclomal）+ MethylBERT/DeepSomatic 的 PMC PDF — 見 `external_validation/_provenance/manual_download_needed.md`（用戶補中）。
- 補後該源 CONTEXT §7 升級 verification/tier。
