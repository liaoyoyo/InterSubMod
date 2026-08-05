# 可選 section 模組 + genome-scale 模式

`tools/build_workstation.py` 的核心 = header + banner + changelog + **逐項卡片 + 判讀 UI + 匯出**。其餘觀察/分析模組透過 `spec.sections[]`（escape-hatch 自訂 HTML）加入，或在 genome-scale 改用 script-102 式密集陣列。

## 透過 `sections[]` 加的可選模組（自訂 HTML 字串）

每個 = `{title, html}`，html 可含 inline SVG / table / 自寫 `<script>` 互動。常見模組（參考 HCC1395 ASM script 102 的實作）：

| 模組 | 用途 | 102 對應 |
|---|---|---|
| **篩選漏斗 funnel** | caller→phasing→ISM→gate 真實數字流 | §① coverage `FN` |
| **即時門檻試算** | slider(CV/Δβ/reads)+mode 即時看通過 TP/FP 與富集 | §⑤ `live()` |
| **圖庫 filter/sort/分頁** | 大量項目按 tier/驗證/判讀/CN 篩、10 鍵排序、24/頁 | §⑥ `filtered()` |
| **即時 charts** | 由當前資料畫分類/Tier/驗證/CN 分布長條 | §④ `bars()` |
| **confound 關連表** | 真值(CN)×分類 Cramér's V + 每類比例 | §④b `renderSqvs()` |
| **redesign meta-loop** | 統計人工判讀(誤收/漏收)→建議標準調整 | §③ `redesign()` |
| **圖例與資料字典** | 兩張圖怎麼讀 + 每欄位直覺/判讀 | §⓪ |

把這些當 `sections[].html` 注入，或在 genome-scale 直接沿用 102。

## genome-scale（>~1000 項）：密集陣列模式

本 generator 的 `items[]` 是 per-item 物件陣列，富但冗長；**> 約一千項**請改 **script-102 式 compact column array**：

- 一個 `const D=[[...22 cols...], ...]` 扁平陣列（每列一項，欄位用 `const C={ck:0,ps:1,...}` 索引）→ 檔案小、載入快。
- 圖一律外部 PNG（`figs/{key}_meth.png` / `_dist.png`），**gitignore**，由 render 腳本重生（HCC1395 = script 101 從 ISM per-region 輸出 render 60,700 張）。
- 卡片/modal/filter/judge/export JS 與本 generator 同核心（localStorage、JSON/CSV、3 態 + reason）。
- 範例：`InterSubMod/research/tsg_promoter_asm_reviewer/scripts/102_build_workstation.py`（30,350 位點），其 §⓪b 修正過程紀錄即本 skill 的 changelog 模組落地。

**判準（圖策略）**：≤ 數百項 + 示意/可 data-bind 的圖 → 本 generator + `figures.mode="svg"`（單檔可攜）。> 數百項 / 真實資料點陣（熱圖、距離矩陣 PNG）→ 102 式 compact array + 外部 PNG + 重生指令。

## ISM 觀察圖（genomics 預設兩張）

- `_meth.png` = read×CpG 甲基熱圖：列依 HP 子單倍型分群再依 mean β 排序；HP + ALT/REF 側欄；變異 CpG 橘虛線；RdBu_r 0–1。驅動檔 `methylation/methylation.csv` + `reads/reads.tsv`。
- `_dist.png` = read×read 距離矩陣：依分群樹 `clustering/leaf_order.txt` 重排；雙軸 HP 側欄（看分群是否對齊 HP = CramersV 是否為真）；magma_r。驅動檔 `distance/<METRIC>/matrix.csv` + `leaf_order.txt` + `reads.tsv`。
- render 參考 `research/tsg_promoter_asm_reviewer/scripts/{85,101}_render_*.py`。
- ⚠ ISM per-region 真實輸出是**巢狀**（`methylation/`、`distance/<METRIC>/`、`clustering/`、`reads/`），比 `.claude/rules/output-structure.md` 的扁平描述更細。

---

## quick-compute vs SVG vs PNG 決策表（每個分析輸出 → 該畫哪種）

> **核心啟發**：**精確/多欄/逐項 → 表（quick-compute）；比例/分佈/流程/比較 → inline-SVG；真實點陣（熱圖/散點）→ 外部 PNG**。任何「要不要畫圖、畫哪種」先查此表，不要每次重想。

| 資料性質 | 呈現 | 為何 | 落地 |
|---|---|---|---|
| 每項精確值（座標 / n / 群數 / V / cn_value） | **quick-compute 數字 / kv 表** | 要可讀精確值、可複製、可排序 | card kv + 展開 `<table>` |
| 總量 / 計數 / TP-FP / 各態 N | **quick-compute 數字卡** | 單一純量，畫圖反而噪音 | `.statgrid` 數字卡（big-number + 小標）|
| **類別比例**（CN gain/loss/loh/neutral、structure 5 態、tier 3 態） | **SVG 水平 stacked 比例條** | 一眼看佔比，不需精確值 | `svg_stacked`（見下目錄）|
| **兩組成比例**（tumor vs normal reads、TP vs FP 某指標） | **SVG 比例條 / 比較條** | 二元/少組對比 | `svg_stacked` 或 `svg_compare` |
| **數值分佈**（coarse/fine/geometry 切群數 1,2,3… 計數） | **SVG 直方** | 看形狀（單峰/長尾/離散） | `svg_hist`（離散 bar）|
| **階段流失流程**（全位點→有結構→對齊→候選） | **SVG funnel** | 看每關卡保留率 | `svg_funnel`（見可選模組表 §①）|
| 每項真實點陣（甲基 read×CpG 熱圖、read×read 距離矩陣） | **外部 PNG**（path-based lazy） | inline 會爆 HTML 體積、瀏覽器卡 | `figs/{key}.png` + gitignore + 重生指令（見密集陣列模式）|

**門檻提醒**：SVG 比例/分佈圖屬「**aggregate 層**」（整份一張，數字來自 stats JSON）；PNG 屬「**per-item 層**」（每項一張，數量大走 path）。**不要把 per-item 真實熱圖塞 inline，也不要把 aggregate 比例硬畫成 PNG**。

## SVG-proportion section module 目錄（aggregate 比例/分佈/流程）

> 全部 = 從 stats JSON 注入數字的 **Python 端 SVG 生成器**（§13-A：generator 不手打數字，數字缺 → 不 render）。viewBox + `role="img"` + `<title>`；色用**規約色**（見「aggregate-chart 配色」段）。

| 模組 | 函式簽名（參考） | 注入契約（從 stats JSON） | 複用 component |
|---|---|---|---|
| **水平 stacked 比例條** | `svg_stacked(segs=[(label,count,color)…], w, h)` | 每段 `(label, S[key][cat], color)`；條內標 `label NN%` | `html-report-build/components/svg_compare_bar.html` |
| **離散直方** | `svg_hist(dist={k:count}, …)` | `S["coarse_dist"]` 等 `{群數: 位點數}` | 手刻（FULL dashboard 12-bar pattern）|
| **流失 funnel** | `svg_funnel(stages=[{stage,n}…])` | `S["funnel"]`（每關 `{stage, n}`，標保留率） | `verify-workstation/references` §① funnel |
| **TP-vs-FP 比較條** | `svg_compare(a_pct, b_pct, labels)` | `S["confident_multi_tp_pct"]` / `…_fp_pct` | `svg_compare_bar.html`（比較模式）|

**注入契約鐵則**：每張 SVG 的每個數字 = `S[...]` 直取；**sum-check**（各段加總 = N）在 aggregate 腳本算好寫進 stats，generator 只讀。worked 實例 = `docs/methodology/_assets/20260618_subcluster_pilot/scripts/{dashboard_aggregate.py（算 stats）, build_phylo_cpp_dashboard_v3.py（注入 SVG）}`。

**aggregate-chart 配色（與規約 1 的關係）**：規約 1（`ism_heatmap_std.py`）定義的是 **per-read 熱圖側欄**色（HP/ALT/T-N/cluster/strand）。aggregate 比例圖是**獨立、帶 legend label** 的圖，**collision 風險不適用**（無相鄰 HP 側欄）。但仍守兩條：① **T/N 比例條必用 canonical** `tumor=#f97316 / normal=#22c55e`；② 其餘 aggregate 軸（CN-state / structure-state / tier）用**有 label 的區別色**、避免整圖只有一個色相。理想上集中到 `ism_heatmap_std` 的 `AGG_*` dict（目前 dashboard generator 因 branch 未合併暫硬編，屬已知小債）。

## 5 層漸進揭露 dashboard 結構（genome-scale 觀察面板範本）

> 參考 MAIN repo `…/obs_ws/wg/20260618_HCC1395_observation_dashboard_FULL.standalone.html` 抽出的層次。**aggregate（看林）在上、per-item（看樹）在下**，全用 `<details>` 漸進揭露 + sticky 控制列。

| 層 | 內容 | 元件 |
|---|---|---|
| **L0** Header + sticky 控制列 | 標題 + 位點搜尋 + 全 filter + sort + 升降序 | sticky bar |
| **L1** 摘要 + 紅線 | `<details open>` scope / ⭐ tier / 2-3 條 caveat（confound / 反判別 / 單樣本）| note callout |
| **L2** 整體統計（aggregate SVG）| `.statgrid` — 每指標一張 stacked/hist/funnel（**SVG-proportion 目錄**）| svg_stacked/hist/funnel |
| **L3** 分類交叉表 | `<details>` cross-tab（類別×類別）+ 小 stacked-bar | quick-compute 表 + svg_stacked |
| **L4** per-locus 卡（分頁 grid）| 預設 badge（set/CN/T-N/切群/LOH）+ 對齊**只敘述** + `<details>` 展開（歸類/切群數/T-N）+ 大圖（path PNG）+ 判讀 UI | compact array + figName JS + localStorage |
| **L5** 參考/方法 | `<details>` display-standard + 方法 + §13 provenance footer | — |

**genome-scale 落地**：L4 用**密集陣列模式**（`const D=[[…cols…]]` + `const C={…}`，圖名 JS 決定式無 FIGMAP）。worked 實例 = `build_phylo_cpp_dashboard_v3.py`（34,736 位點、22 欄、2.1MB、800/800 抽樣圖解析、9 張 aggregate SVG）。
