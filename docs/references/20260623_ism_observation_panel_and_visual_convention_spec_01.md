<!--
建立時間: 2026-06-23
類型: 規約 SoT (convention registry + spec) — 觀察面板 + 視覺規約
狀態: active SoT
緣由: 同一組視覺規約(色盤/側欄/圖片模式/面板)反覆漂移、被重複重修(geom→phylo_cpp_render→dashboard legend 同一個 bug 修三次)。
      根因 = 規約「有紀錄但散落 + 無強制 + 寫時不 surface」。本檔 = 單一合併紀錄 + 關聯 + 定義 + 合規清單。
authority: 本檔不取代以下權威源, 而是「索引+定義+合規」彙整 →
  - 色盤/側欄/colormap/render 函式 = `InterSubMod/scripts/ism_heatmap_std.py`(程式 SoT)
  - 面板 8 步 + 圖策略 = `/verify-workstation` SKILL.md + references/section_modules.md
  - memory: feedback_ism_case_heatmap_standard_sidebars
-->

# ISM 觀察面板 + 視覺規約 — 單一 SoT（registry + spec）

> **一句話**：任何 ISM 觀察/判讀面板（dashboard / workstation）的視覺規約 = 4 條，全部已有權威源；本檔把它們**合併成一份可查、可複用、可稽核**的紀錄，並列**合規狀態**與**複用路徑**。
>
> 🔴 **鐵則**：renderer **一律 import `ism_heatmap_std`、呼叫其 render 函式；禁止本地硬編色盤/側欄/colormap**。違反 = 漂移（已發生多次）。

---

## 規約 1 — 色盤（palette）

**權威源** = `InterSubMod/scripts/ism_heatmap_std.py`（程式 SoT）。**禁止任何 renderer 自定義這些色**。

| 軸 | canonical 值 | 語意/為何 |
|---|---|---|
| **HP1 族（藍）** | `1`=#60a5fa(淺) / `1-1`=#1e3a8a(深) | germline 單倍型骨幹；色相=單倍型、**明度=germline(淺) vs somatic 子標籤(深)** |
| **HP2 族（紫）** | `2`=#c084fc(淺) / `2-1`=#6b21a8(深) | 同上，HP2 |
| HP3 / unphase | `3`=#14b8a6(青) / `0`=#9ca3af(灰) | |
| **ALT 軸** | ALT=#dc2626(紅) / REF=#fbbf24(琥珀) | 帶變異與否 |
| **T/N 軸** | tumor=#f97316(橘) / normal=#22c55e(綠) | |
| **Strand** | +=#334155 / −=#94a3b8(板岩) | 非語意色相 |
| **群 id（cluster）** | `CLUSTER_COL`=[#db2777,#0d9488,#0891b2,#65a30d,#7c3aed,#475569] | **teal-pink，刻意不撞任何語意軸**（無橘/綠/藍/紫/紅/琥珀）|
| SNV marker | #facc15(黃線, 在 CpG 之間內插) | |
| meth colormap | RdBu_r（紅=甲基/藍=未甲基, NaN 灰）| dist=magma（暗=近/亮=遠）|

**最常見違規（本 session 修三次的 bug）**：HP1 寫成 #dc2626(紅) → **撞 ALT 軸的紅**；群色用紅/藍/綠 → 撞 ALT/T-N。

## 規約 2 — 側欄順序

固定（左→右，熱圖在右）：**`[cluster?] HP | ALT | T/N | Strand`**，全部取自 `reads.tsv` + `phylo_groups.tsv`（hp/alt_support/is_tumor/strand）。`sidebar_specs()` 已產好；勿自排。

## 規約 3 — 圖片顯示模式（path vs inline）⭐ 本次確認

**判準 = 圖的數量與性質**（門檻「~數百」為約定非硬數字）：

| 情境 | 模式 | 機制 | 實例（grounded）|
|---|---|---|---|
| **全集 / 大量檢視**：> ~數百項、或**真實資料點陣圖**（甲基熱圖、距離矩陣 PNG）| **外部路徑 PNG** | `figs/{key}.png` 相對路徑 + `loading="lazy"` + **figs 目錄 gitignore** + render 腳本**重生指令** + HTML 內嵌資料陣列(只放路徑) | geom 400(59KB) / FULL phylo 10,707(2MB) / HCC1395 ASM 60,700(script 102) |
| **快速 / 少量檢視**：≤ ~數百項、**示意圖 / 可 data-bind 的圖** | **嵌入 inline** | inline **SVG**(手刻示意，首選) 或 base64 PNG(少量真實圖)；**單檔可攜、無 sibling 依賴** | chr2:18M ~10 項(65KB, 3 inline SVG) |

**鐵則**：大量真實點陣**永遠走 path**（inline 會讓 HTML 爆到數十 MB、瀏覽器卡）；少量示意**走 inline**（單檔好攜帶、不怕 sibling 圖遺失）。**path 模式必附重生指令**（圖 gitignore，否則他人開啟空白）。

## 規約 4 — 檢視面板「如何定義與處理」（panel spec）

面板 = `/verify-workstation` 8 步（W0–W7）的產物。**定義一個面板 = 給一份 data JSON + 一個 generator**，不是手刻 HTML：

1. **W0 數字先驗證落檔**（§13.0）：分析跑完 → 數字落 `.json` → Read 回真值。產數字的 Bash 與寫面板的 generator **不同批**。
2. **W1 卡片 schema**：每項 `{id,title,badges[],metrics[{k,v,src}],figures,modal_stats}`。
3. **W2 圖綁定**：依**規約 3** 選 path / inline。
4. **W3 分類在資料層算**（generator 不重算、人不手填；判準可由 src grep）。
5. **W4 判讀 UI**：每項 verdict 按鈕（同意/存疑/否定 + reason）→ localStorage（key=`meta.lskey`）+ JSON/CSV 匯出。
6. **W5 provenance + changelog**：每 metric 帶 `src`（⌖）；修正過程紀錄 data-driven。
7. **W6 §13-A refuse-on-missing render**：`tools/build_workstation.py <spec.json>`，缺必填欄 → exit 3 拒繪（杜絕捏造/漂移）。數字由 JSON 注入，**不手打**。
8. **W7 visual QC + export**：playwright 截圖檢無溢出；確認 ⌖ 來源、判讀、匯出可動。
9. **色盤/側欄**：面板內任何圖例 SVG / 嵌圖 → 用**規約 1+2**（embedded legend 也禁硬編，FULL dashboard 圖例就犯過）。

**複用路徑（避免重刻）**：`ism_heatmap_std.mpl_dual_panel()` / `sidebar_specs()` / `grouped_legend()` 已提供整套渲染；通用面板 generator = `tools/build_workstation.py`。**新 renderer 先呼叫這些，不要新寫色盤或側欄。**

## 規約 5 — 資料盤點 → 顯示映射（每個分析輸出該畫哪種圖）⭐ 本次新增

> **目的**：消除「這個數據要不要畫圖、畫哪種」每次重想。**先盤點每個分析輸出的性質 → 查下表 → 直接知道呈現方式。** 對齊 `/verify-workstation references/section_modules.md` 的決策表（skill 端 SoT，本檔同步索引）。

**核心啟發**：**精確/多欄/逐項 → 表；比例/分佈/流程/比較 → inline-SVG；真實點陣（熱圖/散點）→ 外部 PNG。**

| 分析輸出性質 | 呈現 | 圖型 | 數字來源 |
|---|---|---|---|
| 每項精確值（座標 / n / 群數 / V / cn_value） | quick-compute | card kv + 展開 `<table>` | per-item 陣列 |
| 總量 / TP-FP / 各態 N | quick-compute | `.statgrid` 數字卡 | stats JSON |
| 類別比例（CN 4 態 / structure 5 態 / tier 3 態） | **SVG** | 水平 stacked 比例條 | stats JSON（每段 count）|
| 兩組成比例（tumor vs normal、TP vs FP 指標） | **SVG** | 比例條 / 比較條 | stats JSON |
| 數值分佈（coarse/fine/geometry 切群數計數） | **SVG** | 離散直方 | stats JSON（`{群數:位點數}`）|
| 階段流失（全位點→有結構→對齊→候選） | **SVG** | funnel | stats JSON（`[{stage,n}]`）|
| 每項真實點陣（甲基熱圖 / 距離矩陣） | **外部 PNG** | path-based lazy（規約 3） | per-region 巢狀輸出 |

**aggregate（看林）vs per-item（看樹）兩層**：SVG 比例/分佈/流程屬 aggregate（整份一張、數字來自 stats JSON、§13-A 注入不手打）；真實熱圖屬 per-item（每項一張、量大走 path）。**禁把 per-item 真實熱圖塞 inline、禁把 aggregate 比例硬畫 PNG。**

**aggregate-chart 配色**：規約 1 的色盤是 **per-read 熱圖側欄**用；aggregate 比例圖帶 legend label、無相鄰 HP 側欄 → collision 不適用，但仍守：① **T/N 比例條必用 canonical** `#f97316`/`#22c55e`；② 其餘軸（CN/structure/tier）用有 label 區別色。理想集中到 `ism_heatmap_std.AGG_*`（目前 dashboard generator 因 branch 未合併暫硬編 = 已知小債）。

**worked 實例（全鏈）**：`docs/methodology/_assets/20260618_subcluster_pilot/scripts/` → `dashboard_aggregate.py`（算 stats JSON + sum-check）→ `build_phylo_cpp_dashboard_v3.py`（SVG 注入 + 5 層 collapsible，34,736 位點/22 欄/9 張 aggregate SVG）→ `obs_ws/cpp_wg/20260623_phylo_cpp_observation_dashboard_v3.html`。

---

## 合規狀態（2026-06-23 稽核，pilot scripts = 104 個 .py）

| 類別 | 數量 | 處置 |
|---|---:|---|
| import 標準模組但**仍硬編本地色盤**（最糟）| **6**：phylo_render_v1v2 / phylo_undersplit_diagnostic / phylo_v3_render_all / phylo_v2_render_all / phylo_edge_render / phylo_v4_edge_render | 改呼叫 H.* |
| 定義本地色盤（候選違規，需逐檔確認是否畫標準側欄）| **15**（含上 6）| 逐檔確認 → 改 import |
| 乾淨 import | ~19 | 合規 |
| HTML 嵌死色圖例 | **2**：`20260620_decisionflow_5state_..._01.standalone.html` + `obs_ws/cpp_wg/...FULL.html` | regen（須先重畫圖）|
| **已修（本 session）** | `geom_render.py` ✅ / `phylo_cpp_render.py` ✅ / `build_phylo_cpp_full_dashboard.py` 圖例 ✅ | 圖待重跑重畫 |
| **強制層（hook/audit）** | **已建** ✅ | `scripts/palette_drift_check.py`（AST 稽核）+ `scripts/hooks/palette_drift_advisor.sh`（PostToolUse advisory）|

## 機械層（已落地 — 為何規約有紀錄卻仍漂移 → 補強制層）

紀錄自 2026-06-18 已存在（module + memory），但先前**無 enforcement + 寫時不 surface** → 靠 copy-paste 在 ad-hoc renderer 間繁殖。**2026-06-23 已建兩件**（對齊本專案 §13 教訓「純文字規則失效，需機械防線」）：
1. **`scripts/palette_drift_check.py`** — AST 稽核：import `ism_heatmap_std` 取 canonical，偵測本地 HP/ALT **hex** 色 dict（要求 hex 值以避開 LABMAP label-remap 假陽性）+ HTML 嵌死圖例漂移。本 session 實測 6 WRONG → 修後 0。
2. **`scripts/hooks/palette_drift_advisor.sh`** — PostToolUse Edit|Write advisory（exit 0）：grep 本地語意色 dict（quote-agnostic），排除 `ism_heatmap_std.py` + `palette_drift_check.py` 本身。比照 `skill_registry_sync` / `number_provenance` 的 advisory 模式。

---
> 關聯：`ism_heatmap_std.py`（規約1-2程式SoT）、`/verify-workstation`（規約4面板SoT）、memory `feedback_ism_case_heatmap_standard_sidebars`、`20260623_weak_structure_classification_and_geometry_divergence_01.md`（消費此規約的最新面板）。
