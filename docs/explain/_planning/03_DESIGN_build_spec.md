---
title: 解釋中心（docs/explain/）建置 spec — 設計定案與第一批範圍
doc_type: design_spec
created: 2026-06-12
status: approved_design
audience: 實作者（自己 / AI）
supersedes: none
references:
  - InterSubMod/docs/explain/_planning/00_盤點_各任務敘述與邏輯.md
  - InterSubMod/docs/explain/_planning/01_敘述節點與模板_taxonomy.md
  - InterSubMod/docs/explain/_planning/02_名詞表_glossary.md
note: 本 spec 為 brainstorming 產出的設計定案；數據鐵則對齊 CLAUDE.md §13.0。
---

# 解釋中心 build spec

> **這份做什麼**：把與用戶確認的設計決策 + 第一批建置範圍 + 數據誠信前置條件 + 迭代協議，固化成單一可執行依據。
> **怎麼讀**：L0 看決策與範圍 → §1-§4 看建置細節 → §5 看迭代回饋協議。

## L0 — 設計定案（用戶 2026-06-12 確認）

| 維度 | 定案 |
|---|---|
| **目的** | 服務新主軸論文與後續研究的「方法/概念解釋」，讓有實驗室背景但非長讀計算基因體背景的讀者（主＝用戶自己，副＝可萃取給 PI）全面理解 |
| **讀者分層** | L0 BLUF（含 `pi_hook` 給 PI 一句話）→ L1 重點邏輯/直覺 → L2 嚴謹細節（`<details>` 折疊）→ L3 溯源 |
| **與舊檔關係** | **並存**：`docs/explain/` = 去時間戳 canonical 教學層（持續更新）；`docs/concepts/*.standalone.html` = 歷史快照（不動） |
| **檔案組織** | INDEX 總覽 + 每概念一個 standalone HTML 子頁 |
| **複用機制** | 14 種敘述節點積木（01_taxonomy §2）+ 固定 9 段頁模板（01 §3）+ 統一名詞表 SoT（02_glossary）|
| **SVG 策略** | 第 0 步**一次批量生 6 個共用 primitive**；既有 5 張 `fig*.svg` 直接複用 |
| **涵蓋範圍** | 當前主軸 + 支撐線優先；新概念即時觸發式加入 |
| **第一批** | INDEX 樞紐 + 背景名詞地基頁 + ISM 五階段頁 |

## §1 目錄結構

```
InterSubMod/docs/explain/
├── index.standalone.html              ← 樞紐：N-DASH 儀表板 + 名詞 anchor 註冊表 + 全站守則(4 鐵則)
├── _assets/
│   ├── explain.css                    ← 共用 CSS（自 CS學生分層版 HTML 抽出 --sp-*/.card/.badge/.funnel/details）
│   └── svg/                           ← 共用 SVG primitive（第 0 步批量生）
│       ├── svg-readxcpg-grid.svg
│       ├── svg-dist-matrix.svg
│       ├── svg-upgma-tree.svg
│       ├── svg-hp-haplotype.svg
│       ├── svg-cis-3way.svg
│       └── svg-2x2-cochran.svg
├── _planning/                         ← 盤點底稿（已產出 00/01/02 + 本 spec 03）
├── 01_background-glossary.standalone.html   ← 名詞地基（所有頁 term 連回此頁 anchor）
└── 02_ism-core.standalone.html              ← ISM 五階段方法本體
```

> 既有 SVG 可直接嵌入：`fig1_ism_5stage.svg`（ISM 頁段 3）、`fig2_single_locus_pipeline.svg`、`fig5_axisC_vs_axisA.svg`（均在 `docs/paper_focus/04_figures/`）。

## §2 第一批三頁的節點組成

| 頁 | 主用節點（01_taxonomy §2）| 關鍵 SVG | 核心內容 |
|---|---|---|---|
| **index** | `N-DASH`（10 線就緒度儀表板）· 名詞 anchor 註冊表 · `N-REDFLAG`（4 鐵則守則）| 無 | 全站導覽 + 樞紐節點順序（methyl-phasing-assist 入度 7 最高）|
| **01_background-glossary** | `N-DEF`×~23（名詞卡）· `N-ROLE`(fig1_5stage) · `N-EVO`(TO vs paired) | svg-hp-haplotype | 全站 term 定義 SoT + 樣本/pipeline 速查 |
| **02_ism-core** | `N-ROLE`(fig1_5stage / fig2_single-locus) · `N-DEF`(Δ-vs-Δβ/NHD/PERMANOVA) · `N-EVO` · `N-CMP` | svg-readxcpg-grid · svg-dist-matrix · svg-upgma-tree | ISM **作用位置**（上下游 I/O）+ 五階段 + old-vs-new + 數據比較 |

## §3 🔴 數據誠信前置條件（建頁前必先解決，守 §13.0）

盤點抓出 4 處同主題不同口徑數字 + 2 個源碼參數歧異。**建到用該數字的頁之前**，先 grep 源碼/原檔定案單一口徑，不並列、不憑記憶：

| 爭議 | 來源歧異 | 解法 | 影響頁 |
|---|---|---|---|
| significance_summary.csv 欄數 | 59（knowledge/05_data_formats）vs 117（method_comparison/01_spec）| grep `RegionResult.hpp` 實際欄位計數 | 名詞地基 + ISM |
| PERMANOVA permutation 預設 | 99 vs 999 | grep `StructureTest.cpp` 生產預設值 | ISM |
| BRCA2 Δβ | −0.054（buggy 雙列砍半）vs −0.122（修正 max-collapse）| 用修正值 −0.122 並標版本 | （第一批外，但名詞卡 example 可能引）|
| phasing n | 6（bio-n）vs 7（含 basecall 變體）| 兩者都對，標明「樣本 instance 7 / 生物 n=6」| index 儀表板 |

**鐵則**：頁內每個數字必能在某 source_path grep 到；grep 不到 → 留 `{{待查}}` 徽章不補值。撰寫頁的 `Write` 與產數字的 `Bash`/grep **不同 batch**。

## §4 建置順序（第一批）

```
第 0 步（一次性基建）
  0a. 解決 §3 數字歧異（grep 源碼定案）
  0b. 抽 explain.css（自 CS學生分層版 HTML）
  0c. 批量生 6 個 svg-* primitive（/methods-example，data_ref 注入；schematic 示意可無真值，但若含數字必走 verified）
第 1 步：01_background-glossary.standalone.html（建全站 term anchor）
第 2 步：02_ism-core.standalone.html（連回名詞 anchor）
第 3 步：index.standalone.html（最後做，因需連到已存在的子頁）
第 4 步：人眼/瀏覽器檢視 → 交付用戶看風格
```

## §5 迭代回饋協議（用戶持續改進「敘述方式 + 圖例」）

1. 第一批產出後，用戶對**敘述方式**（太深/太淺/順序/類比）與**圖例**（SVG 細節/配色/標註）提要求。
2. 我針對性修，並把「好的敘述/圖例範式」沉澱回 `01_taxonomy`（node 模板）與 `_assets/svg/`（primitive），使後續頁自動一致。
3. 新概念/新名詞出現 → 先用本 spec 的 `N-DEF` 三件套（嚴謹定義+直覺類比+例子）；不確定先與用戶討論確認理解，再寫完整頁。
4. 風格穩定後，再依 01_taxonomy §5 建置順序推進第 3-7 頁。

## §6 邊界與不做（YAGNI）

- 不動既有 `docs/concepts/` HTML（並存決策）。
- 不在第一批做全部 7 頁（風格未確認前不批量）。
- 不 commit 到當前 branch（`docs/method-comparison-...`）—— explain 是新主題，commit 需先決定 branch（git governance）；第一批先留工作樹，交付後再議。
- 名詞表只定義概念，不放專案特定數字（數字留各研究線頁引用並附 source）。
