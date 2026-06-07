# 基因組簡報繪圖慣例（從 2 份口試簡報萃取 + 通例驗證）

> 來源：DECK1 何明恩「Identification of Somatic Haplotypes in Tumor/Normal Paired Long-Read Sequencing」(12 張) + DECK2 陳振裕「LongPhase-TO / Joint Somatic Phasing + Purity + LOH」(36 張)，同 PI Yao-Ting Huang。萃取 workflow `wf_b6aa4c96-be9`（7 agents）。
> **Tier**：L1 = deck 實見畫法 / L3 = web 領域通例。**§13 誠實**：deck 內數字（F1/precision/purity/N50…）= **該 deck 自身 claim，非本專案真值**；本檔只取**畫法/物件語意/排版慣例**。
> 原始抽取：`InterSubMod/docs/references/manual/20260607_pptx_convention_study/extract/{deck1,deck}/{text.md,slides.json,media/}`。

## 1. 名詞 Glossary（白話）

| 名詞 | 白話 | 色 |
|---|---|---|
| Germline variant | 生殖細胞突變、可遺傳；het germline = phasing 錨點 | **藍** |
| Somatic variant | 體細胞突變、不遺傳、可致癌 | **橘/紅** |
| Phasing | 把 het 變異依單套染色體分到 HP1/HP2 | — |
| Haplotagging | 把每條 read 依攜帶的 phased germline 變異 tag 到 HP1/HP2 | — |
| Haplotype HP1/HP2 | germline 兩條母/父本鏈（DECK1 標 Maternal/Paternal）| 藍系 |
| Somatic haplotype HP1-1/HP2-1 | germline 鏈上因 somatic 變異再細分的子單倍型（樹狀命名 HPx-1）| 橘系 |
| Tumor-Normal paired / Tumor-Only(TO) | 配對 normal vs 只有 tumor | normal 藍人偶 / tumor 橘人偶 |
| LOH | 雜合性缺失，chromosome-scale | LOH interval 紫/洋紅 |
| Haplotagged read pileup (IGV) | read 依 HP 分群堆疊 + SNV tick | lane fill 編 HP |
| PoN (Panel of Normals) | 正常族群變異庫；候選 − PoN = filtered | — |
| ASM / cis（本專案核心）| 等位特異甲基化 / 同鏈作用 — **兩 deck 皆無此圖**，走 L3 | 甲基 ramp |

## 2. 色彩約定（最 load-bearing；grep 驗證 L1）

跨兩 deck 的 SVG 主色（我親自 grep 抽出 SVG 確認次數）：

| 語意 | hex | grep 次數 | 用途 |
|---|---|---|---|
| germline / normal / reference / HP1-HP2 | `#2E75B6` / `#2F5597`(深) / `#8FAADC`(淺) | 40 / 24 / — | **冷=遺傳** |
| somatic / tumor / alt / HPx-1 | `#C55A11` / `#ED7D31` | 16 / 4 | **暖=後天** |
| unphased / 中性 | `#7F7F7F` / `#BFBFBF` | 12 / 8 | 灰 |

**鐵則**：冷色=遺傳(germline/normal)、暖色=後天(somatic/tumor)，跨 icon/label/bar/pileup 一致。用 CSS `<style>` class（`.grp-germline{fill:#2E75B6}` / `.grp-somatic{fill:#C55A11}`）綁定 → 改一處全域重繪。**category-by-FILL 不換形狀**（變異 marker 一律小圓，類別只靠 fill）。

> ⚠ **甲基化例外（§6 關鍵）**：deck **無甲基化圖**，故 P1/P5/P6 的甲基 cell **不套 deck 字母/藍橘編色**，改 L3 通例：5mC **red>0.5 / blue<0.5 ramp**（IGV/PacBio）或 filled/open lollipop（NanoMethViz）。**藍橘 = haplotype/variant 軸；red/blue ramp = 甲基化軸**，兩套不混。

## 3. 排版 + 字型模板（16:9）

Canvas：13.33×7.5 in = `viewBox="0 0 1280 720"`（px = inch × 96）。

| role | pt | px@1280 | 用途 |
|---|---|---|---|
| deck_title_cover | 38–48 | 50–64 | 封面/封底 hero |
| slide_title | 24（D1）/ 44→18 autofit（D2）| 24–28 左對齊 | 每頁標題（左上 band，**不置中**）|
| section_subhead | 18 | 18–20 | 標題下主題行 |
| body | 14（D1）/ 18（D2）**bold** | 14–18 | 內文/區塊標籤 |
| caption_label | 12–14 | 12–14 | legend 標籤/小註 |
| diagram_microlabel | 9–10.5 | 9–11 | 圖內標/軸刻度/nucleotide |
| footer_pagenum | 18 | 12–14 | 右下 `< N >` |

字級 6 階遞減：標題 44/24 → 區塊 18 → 內文 14-18 → 圖內 12 → 軸 9-10 → micro 5-10.5。

**layout pattern**：① TITLE-BAND（左對齊 bold，y 0..0.6in）② PERSISTENT LEGEND（右上角固定 4 色 key，每張 method slide 重複）③ SECTION-SUBHEAD ④ FIGURE-CENTRIC（大圖佔 40-73%、文字 4-16%）⑤ LEFT-RIGHT TWO-PANEL（compare）⑥ RESULT GRID（等大 2×2 panel ~504×302px）⑦ COVER+CLOSING bookends ⑧ PIPELINE CHEVRON ROW（圓角 stage box h≈30px + 右箭頭）。

## 4. 切格與組合 Guide

**頁層 zone（水平 band）**：Title `y0..0.6in` → Subhead `~0.6..1.1` → Content `~1.1..6.9`（圖/欄）→ Footer `6.9..7.5`（頁碼右下）；Persistent legend 覆蓋 content 右上。

**圖層 segmentation**：method 圖三帶 = 上(輸入軌)–中(pileup/count)–下(Result/interval)；read 矩陣三欄 = 左(read 標籤)–中(序列/甲基格)–右(W 視窗/✕判定)。**基本切格單位 = haplotype 群組**（read 水平堆疊成群、群左貼 HP 標籤、群間細分隔線分層：germline 上 / somatic 下）。

**組合語法（atomic 物件 → method pipeline）**：
1. **stage-flow connector**：圓角階段框 + 單箭頭 + `✕`(排除)/`∈`(歸屬)/算式符號串 primitive。
2. **就地算式可視化**：步驟畫成算式（`candidate − PoN = filtered`）= deck 教學式核心。
3. **brace-aggregator**：多縮圖 → 大花括號 → arrow → equation → 純量輸出（purity 式）。
4. **progressive-build / (cont.)**：複雜圖跨多張漸進揭露（Parsing→Peak→Calling→Result）。
5. **small-multiples**：固定一軸=SAMPLE、一軸=METRIC → 跨樣本一致性一眼讀；legend 抽出共用。
6. **method-vs-baseline 同色 solid/dashed**：同 hue solid=method、dashed=baseline。

對齊 skill case template：每 case = `title+subhead`(zone) → 依目標選 primitive → stage-flow/brace 連接 → persistent legend → footer；真值走 data.json(verified+src)、示意走 schematic（§13-A 物理分離）。

## 5. 通例符合度 Verdict（L3 驗證）

**兩 deck 皆 CONFORM 基因組通例**：IGV haplotype-split pileup（read 依 HP 分色 band、grey=unphased、彩色 mismatch tick）= whatshap/EPI2ME/JBrowse 通例；藍/橘 HP 色法直覺；benchmark grouped bar = 標準 tool-comparison；LOH cytoband ideogram + array 對照軌 = CNVkit/array 通例；faceted small-multiples / y=x calibration / solid-dashed = variant-calling 論文通例；typography 階層乾淨。

**1 個 CAUTION**：benchmark bar 的 **zoomed 非零軸**（如 0.84-0.91）放大微小 delta = 已知感知陷阱 → skill 支援但**必標 caveat**；且這些是 deck claim 非本專案真值。

## 6. Master technique（所有 primitive 共用）

1. **單一藍/橘 accent-pair 當全域 group code**（CSS class 綁 fill）。
2. **category-by-FILL 不換形狀**。
3. **relocatable `<g transform>`**：每 haplotype block 包成可定位/複製的 group。
4. **vector-first for schematics, raster only for real data**：method 示意全手建 SVG；真實 data（IGV/matplotlib）才 generate+embed（對齊 §13-A）。
5. **schematic factory / 小倍數**：一個 panel `<symbol>` reuse N 次（換 title+data），讀者學一次 legend。

**全部可 HTML SVG 重現**：read=`<polygon>`五邊形箭頭、tick=`<rect>`/`<circle>`、legend=`<circle>`+`<text>`、流程=`<path>`箭頭、IGV pileup=多`<rect>`+彩色 tick、grouped bar=`<rect>`+末端`<text>`、ideogram=帶狀`<rect>`+centromere、graph=`<ellipse>`node+`<line>`edge。deck 的 IGV/array 為 raster 截圖 → **須 re-synthesize 成 SVG**。

---

## 7. 頂尖期刊整合（Cell/Nature/Science，2026-06-07 wf_76c92c4c-0cf；全 L3）

> 全期刊規格 = **L3**（作者指南，非本專案真值）。完整貯備計畫 + 改進 backlog 見 `stockpile_plan.md`。

### Top 共通特性（三家都有 → 優先）
1. **固定投稿欄寬**：單欄 ~89mm / 欄半 ~120mm / 雙欄 ~180mm；最大高 ~170mm（Nature 最嚴）。
2. **body 字 5-7pt（絕對下限 5pt）**；但 **px≠pt**：slide 模式維持 ≥8px 下限，journal 模式才在 mm 換算層查 5-7pt。
3. **panel label = bold 8-9pt upright，never italic**（Nature 小寫 a/b/c、Cell/Sci 大寫 A/B/C）。← 最大結構缺口（待 A1）。
4. **line/stroke 下限 0.5pt**；相鄰灰填亮度差 ≥20%。
5. **色盲安全（Wong/Okabe-Ito 金標）+ 絕不單靠顏色**（配形狀/位置/直接標籤；避 red-green；對比 ≥4.5:1）。
6. **sans-serif 全圖一致**（Arial/Helvetica）；只變數/希臘字斜體。
7. **vector 優先（PDF/EPS）、text 保 live 不 outline** ← 我們 SVG-first 正中要求 ✓。
8. **RGB 投稿**（期刊自轉 CMYK）；印刷飽和橘紅偏移 → 優先 vermillion `#D55E00`。
9. **no-decoration / data-ink**（無 drop-shadow/3D/裝飾格線）← 待 lint C12 enforce。
10. **每圖 ≤6-8 categorical 色**（4-6 更安全）← 我們 4 色 legend 在限內 ✓。
11. **一圖一概念 + 明確 reading direction**（線性 L→R / 階層 T→B / 並列 / 環狀）。

### Wong 色盲安全 palette（journal_a11y 平行 token；slide 仍用 deck 藍橘）
| 語意 | deck（pi_slide，現用）| Wong（journal_a11y）|
|---|---|---|
| germline / HP1 | `#2F5597` / `#2E75B6` | `#0072B2`（藍）|
| somatic / HP1-1 | `#C55A11` | `#E69F00`（橘）/ 印刷 `#D55E00`（vermillion）|
| 其他 Wong | — | 黑 #000 / sky #56B4E9 / green #009E73 / yellow #F0E442 / purple #CC79A7 |

> 藍↔橘**本就是最安全 CVD 對**；我們現用 hex 已 CVD-safe，journal 投稿時換 Wong canonical hex。甲基 red/blue ramp 與 haplotype 藍橘軸**正交**，且甲基已 filled/open 形狀冗餘（符「不單靠顏色」）。

### Cell graphical abstract（選配 preset）
方形 1200×1200px @300dpi、Arial 12-16pt（比 in-figure 大，因縮圖）、單一 take-home、**無 data items（純 schematic）**、明確 reading flow。Nature/Science 不要求。

### 改進落點（待逐批；見 stockpile_plan.md backlog）
renderer A1-A7（panel 編號 / export_target / type_scale / color mode / stroke floor / font-family / reading_flow）；lint C8-C12（CVD / grayscale / panel-label / single-take-home / no-decoration）。
