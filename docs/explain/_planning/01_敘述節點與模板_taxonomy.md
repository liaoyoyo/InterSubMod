---
doc_type: planning / system-design
title: 敘述節點型錄 + 概念子頁模板 taxonomy（解釋系統地基）
created: 2026-06-12
branch: docs/method-comparison-ism-external-202606
scope: 為 docs/explain/ 下每個概念子頁提供「可複用敘述節點 + 固定頁模板 + 一致性規則」
inputs: 10 條研究線盤點 JSON（reusable_nodes / high_cog_load / glossary / key_data / old_vs_new）
status: design-spec（尚未開始實作子頁；本檔=構造規範，非含真值報告）
provenance_note: |
  本檔是「解釋系統的構造規範」，不含任何 metric 真值——所有示意數字均為佔位符 {{...}}
  或標注「資料來自盤點 JSON / 待 grep 原檔」。實作子頁時每個數字必走 §13.0 序列（分析→寫檔→Read→才寫）。
related_layered_disclosure: feedback_reports_need_layered_disclosure（L0→L1→L2→L3）
related_examples:
  - InterSubMod/docs/concepts/2026/06/20260609_G6_研究構想三層架構_漏斗研究卡證據鏈_01.standalone.html（既有 L0-L3 + 三色 badge 範本）
  - InterSubMod/docs/paper_focus/04_figures/fig1_ism_5stage.svg（既有 5-stage SVG primitive）
  - InterSubMod/docs/paper_focus/04_figures/fig5_axisC_vs_axisA.svg（既有 axis 對照 SVG primitive）
---

# 敘述節點型錄 + 概念子頁模板 taxonomy

> **這份是什麼**：把 10 條研究線盤點裡所有 `reusable_nodes` 去重、抽象成一套**可複用的敘述積木（node catalog）**，再定一個**固定的概念子頁骨架（page template）**，讓後續每個 `docs/explain/<concept>.standalone.html` 都用「組合節點」產生，而非每頁重寫。目標：一致、好懂、低認知負擔、防捏造。
>
> **怎麼用**（給實作子頁的 AI / 人）：① 查 §2 node catalog 挑需要的節點 → ② 照 §3 page template 把節點按固定順序排 → ③ 套 §4 一致性規則（分層/標記/term/provenance）→ ④ 照 §5 建置順序決定先做哪頁。

---

## L0 一眼結論（先讀這段就夠導航）

- **14 種敘述節點**（§2）：8 種「常見必備型」（名詞定義 / 方法作用位置 / old-vs-new / NEGATIVE 卡 / 數據比較表 / confound 警告 / claim card / SVG primitive 清單）+ 6 種「本專案高頻特化型」（機制天花板卡 / memory-stale 修正卡 / 軸對照矩陣 / 漏斗圖 / tier 儀表板 / overclaim 紅旗表）。
- **1 個概念子頁模板**（§3）：固定 9 段骨架 = `frontmatter → L0 BLUF → 名詞地基 → 方法作用位置 → 主敘述（為何/發現/結論）→ old-vs-new → 高認知負擔拆解 → 數據/NEGATIVE/confound 證據區 → 關聯導覽 + provenance footer`。L0 平鋪、L2-L3 用 `<details>` 收。
- **一致性規則**（§4）：4 層揭露（L0-L3）+ 三色×符號標記（🔴⭐◽ 重要性 / 🟢🟡🔵 證據 tier）+ term 首現必定義（連到名詞節點 anchor）+ provenance footer 帶 commit hash。
- **建置順序**（§5）：先做支撐「當前主軸」最直接、且節點複用率最高的 7 頁——背景名詞表 → ISM 五階段 → HP-axis vs ALLELE-axis → 軸 C vs 軸 A → 四道 NEGATIVE → LOH-phasing 脊柱 → normal-anchored cis-test。
- **🔴 鐵則**：本型錄定義的是**構造**；任何節點一旦填入 metric，必走 §13.0（先有 verified 數字才撰寫；產數字與寫頁不同 batch）。型錄裡所有數字都是佔位符。

---

## §1 設計原則（為什麼是「節點 + 模板」而非每頁重寫）

| 原則 | 說明 | 來源依據 |
|------|------|---------|
| **複用優先** | 10 條線的 `reusable_nodes` 高度重疊（NEGATIVE 卡出現 6 次、confound 警告 6 次、HP-axis/ALLELE-axis 定義出現 4 次、5-stage SVG 出現 3 次）→ 抽成單一 node，多頁引用同一 anchor，改一次同步全部。 | 盤點 JSON 去重統計（§2 表） |
| **分層揭露** | 同一頁服務「掃一眼就走」與「要鑽到源碼行號」兩種讀者 → L0 平鋪、L2-L3 折疊。 | `feedback_reports_need_layered_disclosure`；既有 G6 三層 HTML |
| **構造防捏造** | 數據比較表 / claim card / key_data 一律「從盤點 JSON 的 `value` + `source_path` 注入」，缺 source 直接留 `{{待查}}` 不杜撰。 | CLAUDE.md §13-A（template+data 注入）、§13.0 |
| **誠實邊界先行** | 本專案核心張力是「存在 ≠ 可判別 / 可說明 ≠ 可解釋有證據」→ NEGATIVE 卡與 confound 警告是一等公民節點，不是附註。 | main-axis / negatives / asm-* 三線一致 |
| **既有資產不重造** | `fig1_ism_5stage.svg`、`fig5_axisC_vs_axisA.svg`、`fig2_single_locus_pipeline.svg` 已存在 → SVG primitive 節點直接 `<img>` 嵌入或 inline，不重畫。 | `docs/paper_focus/04_figures/` 實測 |

---

## §2 敘述節點型錄（Node Catalog）

### 2.0 去重彙總（10 線 reusable_nodes → 14 種 node type）

| node id | 中文名 | 出現線數 | 是否含 SVG | 對應盤點 reusable_nodes |
|---------|--------|:---:|:---:|---|
| `N-DEF` | 名詞定義節點 | 9 | 選配 | 各線 glossary + 「名詞定義節點」 |
| `N-ROLE` | 方法作用位置節點 | 6 | **是** | 「SVG 示意（pipeline 位置圖）」「pipeline 流程圖」 |
| `N-EVO` | 方法演進 old-vs-new 對照節點 | 9 | 選配 | 各線 old_vs_new + 「方法演進對照節點」 |
| `N-NEG` | NEGATIVE 結論卡 | 8 | 否 | 「NEGATIVE 結論卡」×8 線 |
| `N-CMP` | 數據比較表節點（含 tier badge） | 9 | 否 | 「數據比較表」「能力對照矩陣」 |
| `N-CONF` | confound 警告節點 | 7 | 選配 | 「confound 警告」×7 線 |
| `N-CLAIM` | 證據分級 claim card（L1-L5 badge） | 5 | 否 | main-axis「tier 標記節點」+ scientific-rigor 對齊 |
| `N-SVG` | SVG 示意 primitive 清單 | 10 | **是** | 各線「SVG 示意」 |
| `N-CEIL` | 機制天花板卡（特化） | 4 | 選配 | methyl-phasing「機制天花板節點」 |
| `N-STALE` | memory-stale 修正卡（特化） | 2 | 否 | ism-evolution「memory stale 修正卡」「機制更正卡」 |
| `N-AXIS` | 軸分類矩陣節點（特化） | 3 | **是** | ism-vs-external「6軸分類矩陣」「能力對照矩陣」 |
| `N-FUNNEL` | 漏斗/收斂圖節點（特化） | 3 | **是** | methyl-phasing「unphase 漏斗」negatives「多道 NEGATIVE 收斂圖」 |
| `N-DASH` | 就緒度 / tier 儀表板節點（特化） | 2 | 選配 | main-axis「tier 標記節點」「研究線就緒度儀表板」 |
| `N-REDFLAG` | overclaim 紅旗表節點（特化） | 2 | 否 | main-axis「overclaim 紅旗列表」 |

> **複用率最高 3 名**：`N-DEF`(9)、`N-EVO`(9)、`N-CMP`(9) → 這三種模板最該先寫好、最穩定。

---

### 2.1 常見必備型（8 種）

#### `N-DEF` — 名詞定義節點

- **用途**：term 首次出現時，給「嚴謹定義 + 直覺類比 + 具體例子」三件套，並成為全站 anchor 供其他頁 `#term-xxx` 連回。
- **L0/L1/L2 內容結構**：
  - **L0**（卡片標題行，永遠可見）：`<術語> = <一句白話>`（≤20 字）。
  - **L1**（直覺類比，預設可見）：盤點 glossary 的 `intuition` 欄。
  - **L2**（`<details>` 折疊）：盤點 glossary 的 `rigorous_def`（含 file:line / 公式 / 判準）+ 1 個 worked example（取 high_cog_load 的 `suggested_example` 文字部分）。
- **何時用**：每頁的「名詞地基段」；任何 term 在該頁首現處。**同一 term 全站只定義一次**（在其 home 概念頁），其他頁用 inline `<abbr title>` + 連回 anchor。
- **是否含 SVG**：選配——`HP1-1`、`Δ vs Δβ`、`HP-axis vs ALLELE-axis` 等空間性術語建議配 `N-SVG`（見 2.8）。
- **資料來源**：直接注入各線 `glossary[].rigorous_def / intuition`，無需重寫。
- **候選 term 全集（去重，建議 anchor 命名）**：`#term-nhd` `#term-hp11` `#term-permanova` `#term-delta-vs-dbeta`（合卡）`#term-cis-test` `#term-cramersv-gate` `#term-epipoly` `#term-hp-axis` `#term-allele-axis` `#term-over-dispersion` `#term-beta-binomial` `#term-max-collapse` `#term-loh` `#term-cnloh` `#term-asm` `#term-excess-over-null` `#term-loso` `#term-circularity` `#term-vestigial` `#term-5mc-5hmc` `#term-unphase` `#term-hp3` `#term-anchor-auc`。

#### `N-ROLE` — 方法作用位置節點

- **用途**：用 pipeline 示意圖回答「這個方法/工具坐在整條流程的哪裡，吃什麼、吐什麼」。
- **L0/L1/L2**：
  - **L0**：一句話 = `<上游> → 【本方法】 → <下游>`（取自各線 `method_role` 首句）。
  - **L1**：輸入清單 + 輸出清單（bullet，4-5 項）+ 嵌入 SVG（`fig1_ism_5stage.svg` 或 `fig2_single_locus_pipeline.svg`）。
  - **L2**：各階段模組名 + 源碼參數（如 `MAPQ≥20·len≥1000`、`C_min=5`、`±1000 bp`），與上下游邊界（哪些在 C++ binary 內、哪些在 Python 後處理）。
- **何時用**：ISM 本體頁、ism-evolution、ism-vs-external、background-glossary 的 pipeline 段。
- **是否含 SVG**：**是**——優先嵌既有 `fig1_ism_5stage.svg`（5 階段）；單位點細節用 `fig2_single_locus_pipeline.svg`。
- **⚠ 注意**：盤點裡 PERMANOVA permutation 數 99 vs 999 在不同檔不一致（ism-core uncertainties）→ 此節點標的源碼參數必逐一 grep 確認，不可照抄。

#### `N-EVO` — 方法演進 old-vs-new 對照節點

- **用途**：把「當時做法 vs 新做法」並排，加「差異理由 + 數據比較」，講清楚每一步為何改。
- **L0/L1/L2**：
  - **L0**：`<aspect>：old → new`（一行）。
  - **L1**：兩欄並排卡（左 old / 右 new_）+ 一句 `diff`（為何改）。
  - **L2**（`<details>`）：`data_compare`（含 tier badge，走 `N-CMP` 格式）+ `source`（file:line）。
- **何時用**：ism-evolution（12 個 aspect）、ism-vs-external（4 個）、main-axis（5 個）、loh-phasing（4 個）、negatives（3 個）、background-glossary（2 個 pipeline 對照）。
- **是否含 SVG**：選配——軸別/分群策略類建議配並排 SVG（如 fig5）。
- **資料來源**：各線 `old_vs_new[]`（已有 old/new_/diff/data_compare/source 五欄，直接注入）。
- **🔴 注意**：盤點裡有多個 `data_compare: "{{待查}}"` 或 `ILLUSTRATIVE`（如 ism-evolution Fisher 合成示意 ρ=0.3 → 53-68%）→ 必照原樣標 `{{待查}}` / `ILLUSTRATIVE 非真資料`，不可補成真值。

#### `N-NEG` — NEGATIVE 結論卡

- **用途**：把一個「已證偽 / DEAD」方向講成正資產——結論 + 判定方法 + 為何是論文價值 + 防重開邊界。
- **L0/L1/L2**：
  - **L0**：紅框卡標題 `<方向> = DEAD / NEGATIVE`（如「甲基→TP/FP filter = DEAD」）。
  - **L1**：三件套 = ①核心判定數字（如 LOSO held-out `{{val}}`、AUC `{{val}}`）②判定方法一句（LOSO / shuffle null / OR enrichment）③**為何是貢獻**（reviewer 難攻破 / 界定方法邊界）。
  - **L2**：完整證據鏈（多道 NEGATIVE 列點）+ 失敗機制（如「甲基=caller_af 代理」「regression-to-extreme」）+ **reopen 條件**（C1 新數據 / C2 新方法 / C3 新前置，對齊 `feedback_productive_failure_reopen_threshold`）。
- **何時用**：negatives 線（四道收斂）、main-axis、ism-vs-external、asm-zar1l-brca2（B-discrimination）、asm-cis-confound、methyl-phasing（T2 OVERSTATED / T3 NEGATIVE）、background-glossary（TO pipeline DEAD）。
- **是否含 SVG**：否（必要時配 `N-FUNNEL` 多道收斂圖，見 2.12）。
- **🔴 邊界紀律**：每張 NEGATIVE 卡必明標「這個 NEGATIVE 是哪一層」——存在性 vs 判別力（asm 線最易混）；filter-direction vs characterization-value。

#### `N-CMP` — 數據比較表節點（含 source / tier badge）

- **用途**：任何多值對照（樣本×指標、old×new、TP×FP×FN）的標準表格，**每格數字後面掛 tier badge + source anchor**。
- **L0/L1/L2**：
  - **L0**：表標題 + 一句 takeaway。
  - **L1**：表格本體（行=樣本/類別，列=指標），每個 metric cell = `值 + 🟢/🟡/🔵 + (來源 hover)`。
  - **L2**：`<details>` 展開完整 source path（file:line）與口徑 caveat（單樣本？single-pipeline？subhap-matched？）。
- **何時用**：幾乎每頁的證據區；尤以 asm-cis-confound（6 樣本 TP/FP/FN 矩陣）、loh-phasing（7 樣本 gap 表）、negatives（in-dist vs LOSO）。
- **tier badge 對應**（與盤點 `key_data[].tier` 對齊）：
  | badge | 盤點 tier | 意義 |
  |---|---|---|
  | 🟢 P | `P` | 原始可信，原檔可 grep 對賬 |
  | 🟡 P-caveat / S | `P-caveat` `S` | 方法存疑 / 二次紀錄；引用前查口徑 |
  | 🔵 framing | （無 tier 的定位句） | 立場敘述，非硬數字 |
  | 🔴 F | （已捏造，不應出現） | 一旦發現即下架 |
- **🔴 構造防捏造**：表格一律 `value` + `source_path` 從盤點 `key_data[]` 注入；無 source 的格留 `{{待查}}`。

#### `N-CONF` — confound 警告節點

- **用途**：標出「這個數字/結論有一個混淆，沒控制就會誤讀」，並指出正確控制法。
- **L0/L1/L2**：
  - **L0**：黃框 `⚠ confound：<一句>`（如「ALLELE-axis 被 germline-het baseline 污染」）。
  - **L1**：混淆是什麼 + 後果（誤讀方向）+ 正確控制（如「跑 germline-het negative control」「只信 HP-axis」「引 excess-over-null 不引 raw rate」「引相對 null delta 不引絕對 AUC」）。
  - **L2**：對照數字（confounded vs controlled，如 `ALLELE TP 11.1% < null 15.2%` vs `HP-axis TP 7.2% > null 4.1% OR=1.79`）。
- **何時用**：asm-* 三線、methyl-phasing（絕對 AUC 膨脹 null 0.974）、main-axis（BRCA2 退役）、background-glossary（HP-axis held-constant CN）。
- **canonical confound 清單（全站共用）**：①ALLELE-axis baseline allelic confound ②normal-anchored cis-test 只在有 matched normal 才完整（限 HCC1395）③絕對 AUC 方法樂觀（用相對 null）④single-pipeline 共用偏差（tier 封頂 ⭐3）⑤LOH→低覆蓋→極端 baseline 三角 ⑥spatial autocorrelation（coverage 校正後 O11 AUC 0.845→0.530）⑦by-construction circularity（Inner same-HP1 用 HP1-1 定義）。

#### `N-CLAIM` — 證據分級 claim card（L1-L5 badge）

- **用途**：把一個論文 claim 標上證據等級（L1 實測 / L2 / L3 推測 / L4 降級 / L5 ceiling），對齊 `scientific-rigor` evidence ladder + tier ⭐。
- **L0/L1/L2**：
  - **L0**：claim 一句 + `[L?]` badge + `⭐?` tier。
  - **L1**：支持證據一句 + 限制一句（單樣本？n=?）。
  - **L2**：DAG / pre-reg 對照（重大 claim 才有）+ ledger entry 引用。
- **何時用**：main-axis（R1-R6 results）、需要區分「可說明 vs 可解釋有證據」的任何頁（沿用既有 G6 HTML 的三色證據階梯）。
- **是否含 SVG**：否。
- **與 `N-NEG` 區別**：claim card 是「正向 claim 標級」，NEGATIVE 卡是「DEAD 方向講價值」；二者可同頁並存。

#### `N-SVG` — SVG 示意 primitive 清單（共用圖庫）

- **用途**：登記全站可複用的 SVG 視覺積木，避免每頁重畫。詳見 §2.8。

---

### 2.8 `N-SVG` SVG primitive 清單（共用圖庫）

> 命名規範：`svg-<concept>`；既有檔直接引用，新需求才用 `/methods-example` 或 `/image-gen` 生成，且**含嵌入數字者必走 data_ref 注入（缺 verified 真值 refuse）**。

| primitive id | 描述 | 狀態 | 來源 / 生成方式 | 被哪些頁用 |
|---|---|---|---|---|
| `svg-5stage` | ISM 五階段 pipeline（色塊+模組名+參數） | ✅ 既有 | `04_figures/fig1_ism_5stage.svg` | ISM 本體、ism-evolution、background |
| `svg-single-locus` | 單位點 read 擷取→矩陣→距離→分群 流程 | ✅ 既有 | `04_figures/fig2_single_locus_pipeline.svg` | ISM 本體、ism-vs-external |
| `svg-axisC-vs-axisA` | 軸 A 率差條圖 ‖ 軸 C read×CpG 格陣 並排 | ✅ 既有 | `04_figures/fig5_axisC_vs_axisA.svg` | ism-vs-external、main-axis |
| `svg-significance` | 顯著性檢定套件示意 | ✅ 既有 | `04_figures/fig6_significance_schematic.svg` | ISM 本體 |
| `svg-brca2-dbeta` | BRCA2 per-CpG Δβ 條形 | ✅ 既有 | `04_figures/fig_brca2_dbeta.svg` | asm-zar1l-brca2 |
| `svg-readxcpg-grid` | N reads × M CpG 二元格（甲基=實/未=空/ambiguous=灰） | ⬜ 待生成 | `/methods-example` case；data_ref 注入 | high_cog_load「為何要距離矩陣」×多線 |
| `svg-dist-matrix` | N×N read-read 距離熱圖（顏色=NHD） | ⬜ 待生成 | `/methods-example` | Δ vs Δβ 卡、ISM 本體 |
| `svg-upgma-tree` | UPGMA 樹 + HP label 疊色 | ⬜ 待生成 | `/methods-example` | ISM 本體、軸 C 解釋 |
| `svg-hp-haplotype` | 兩條染色體 HP1/HP2 + somatic ALT 分裂出 HP1-1 | ⬜ 待生成 | `/methods-example` | HP1-1 定義、HP-axis vs ALLELE-axis、LOH-phasing |
| `svg-cis-3way` | normal-HP1 / tumor-HP1 / tumor-HP1-1 三盒 + d_cis/d_drift/d_somatic 箭頭 | ⬜ 待生成 | `/methods-example` | cis-test 卡、ism-evolution 4-way |
| `svg-circularity` | HP1-1 定義 → Inner bucket 篩選 閉環 + R-SELFREF 打破 | ⬜ 待生成 | `/image-gen` concept | loh-phasing HD-1、main-axis |
| `svg-unphase-funnel` | 13.8M unphase → 有甲基 → 有錨點 → 可救 ~6% 漏斗 | ⬜ 待生成 | `/image-gen` flow | methyl-phasing chicken-egg |
| `svg-6axis-landscape` | 業界甲基 6 軸地景 + ISM 站位 | ⬜ 待生成 | `/image-gen` concept | ism-vs-external（亦見 `N-AXIS`） |
| `svg-2x2-cochran` | Cramér's V 稀疏表 vs 密集表 2×2 期望格 | ⬜ 待生成 | `/methods-example` | CramérV gate 卡、latent 結構 |
| `svg-roc-operating` | ROC 曲線 + 「TP loss ≤2%」操作點 → FPR≈1 | ⬜ 待生成 | `/methods-example` | negatives「安全約束下 FP removal=0%」 |

---

### 2.9 本專案高頻特化型（6 種）

#### `N-CEIL` — 機制天花板卡

- **用途**：用一句「機制本質」統一解釋一組看似矛盾的結果（為何 A 成立 B 失敗）。
- **L0/L1/L2**：L0=機制一句（如「甲基 = germline-haplotype 層級訊號：分不同 haplotype 強、within-haplotype 弱」）；L1=這句如何解釋 T1/T2 成立、T3 失敗；L2=決定性證據（如 V10 matched normal 0.979≥tumor 0.866 6/6 chr）。
- **何時用**：methyl-phasing（T1/T2/T3 統一框架）、asm-cis-confound（存在 vs 判別）。是「把散落結論收斂成一句」的高價值節點。

#### `N-STALE` — memory-stale 修正卡

- **用途**：列「舊筆記宣稱 vs 稽核後真相」，防舊 memory 誤導。
- **L0/L1/L2**：L0=「3 條 stale 已修」；L1=表格（舊宣稱 | 真相 | commit）；L2=各條根因。
- **何時用**：ism-evolution（De Waele→Orjuela 已修 891e04b / KDE 已 wired / CramérV gate 合理）、loh-phasing（HPFineNGroups 從「甲基雙峰」更正為「phasing bucket count」）。

#### `N-AXIS` — 軸分類矩陣節點

- **用途**：MECE 軸分類（業界甲基 6 軸）+ 能力對照矩陣（ISM 10 能力 × 11 工具 ✓/⚠/✗）。
- **L0/L1/L2**：L0=「ISM 站軸 C」；L1=6 軸地景圖（`svg-6axis-landscape`）+ 每軸代表工具；L2=10×11 能力矩陣 + ISM 獨有 5 / 欠缺 6。
- **何時用**：ism-vs-external（主節點）。
- **是否含 SVG**：是。

#### `N-FUNNEL` — 漏斗 / 收斂圖節點

- **用途**：把「逐層篩選」或「多道證據收斂到一結論」視覺化。
- **何時用**：methyl-phasing（unphase 漏斗 `svg-unphase-funnel`）、negatives（O1-O13→G1-G7→...→DEAD 收斂時間軸）、ism-vs-external（研究構想漏斗，沿用 G6 HTML 的 `.funnel` CSS）。

#### `N-DASH` — 就緒度 / tier 儀表板節點

- **用途**：一表掃完所有研究線的狀態（phase / tier ⭐ / USABLE-FAILED-NEEDS-WORK / 就緒度）。
- **何時用**：main-axis（T1-T7 儀表板）、explain 首頁 INDEX（兼當全站導覽）。

#### `N-REDFLAG` — overclaim 紅旗表節點

- **用途**：OR-1~OR-6 + FT-1~FT-5 速查，每條一句「識別模式 + 正確寫法」。
- **何時用**：main-axis；亦可放 explain 首頁當「讀本站時的防誤讀守則」。

---

## §3 概念子頁模板（Page Template）

> 每個 `docs/explain/<concept>.standalone.html` 一律照此 9 段骨架。L0 平鋪、L2-L3 用 `<details open?>` 收。節點以 `[N-XXX]` 標出該段用哪種積木。

```
┌─ 0. FRONTMATTER（HTML <head> + meta）
│   doc_type / title / created / commit / source_paths[] / provenance_note
│   <style> 引共用 CSS（沿用 G6 HTML：--sp-* 變數 / .card / .badge / .funnel / details）
│
├─ 1. L0 BLUF 卡（永遠可見，≤6 行）          ［N-CLAIM 縮版 + takeaway］
│   · 一句話結論（one_line）
│   · pi_hook（給 PI 的一句）
│   · 3-5 個 headline 數字（每個帶 🟢/🟡 badge + anchor）
│   · 「怎麼讀這頁」一行（指向各 L 段）
│
├─ 2. 名詞地基段（term 首現都在這定義）       ［N-DEF × n］
│   · 該頁所有新 term 的定義卡（L1 類比可見 / L2 嚴謹折疊）
│   · 已在他頁定義的 term → inline <abbr> + 連回 anchor，不重複
│
├─ 3. 方法作用位置段（方法/工具類頁才有）     ［N-ROLE + svg-5stage/single-locus］
│   · 上游 → 本方法 → 下游 一句 + 輸入/輸出 bullet + 嵌 SVG
│   （純結論頁如 negatives 可省此段，改放「研究設計位置」一句）
│
├─ 4. 主敘述段：為何做 / 發現什麼 / 結論       ［散文 + 內嵌 N-CMP/N-NEG anchor］
│   · 三幕：Why（動機 + gap）→ What（發現，分點）→ So-what（結論 + 誠實限制）
│   · 取自盤點 narrative_logic，但每個數字連到 §8 證據區 anchor
│
├─ 5. old-vs-new 演進段（有演進的頁才有）      ［N-EVO × n + 選配 svg-axisC-vs-axisA］
│   · 並排對照卡，L2 折疊 data_compare + source
│
├─ 6. 高認知負擔拆解段（每頁 2-5 個難點）      ［N-DEF 加強版 + N-SVG + N-CEIL］
│   · 每個 high_cog_load: 為何難（why_hard）+ 類比/SVG（suggested_example）
│   · 能用一句機制統一的 → 升級成 N-CEIL 機制天花板卡
│
├─ 7. 證據區（數據 / NEGATIVE / confound）     ［N-CMP + N-NEG + N-CONF + N-CLAIM］
│   · 7a 數據比較表（key_data 注入，tier badge）
│   · 7b NEGATIVE 卡（DEAD 方向 + reopen 條件）
│   · 7c confound 警告（混淆 + 控制法）
│   · 7d stale 修正卡（若該線有 memory drift）  ［N-STALE］
│
├─ 8. 關聯導覽段                               ［related_concepts 連結 + N-DASH 縮版］
│   · related_concepts[] → 連到其他 explain 頁
│   · 「這頁在全站的位置」一句（指向 INDEX 儀表板）
│
└─ 9. PROVENANCE FOOTER（永遠可見）
    · source_paths[] 全列（絕對路徑）
    · uncertainties[]（{{待查}} 清單，誠實揭露）
    · commit hash + 生成日期 + 「Write 與產數字 batch 分離」聲明
    · <!-- provenance-verified: ... --> 註解（若引他處 validated 數字）
```

**模板變體**（不同 thread 類型微調，但骨架不變）：
- **方法/工具頁**（ism-core, ism-evolution, ism-vs-external, background-glossary）：段 3 必有、段 5-6 重。
- **結論線頁**（negatives, asm-zar1l-brca2, asm-cis-confound）：段 3 改一句設計位置；段 7 重（NEGATIVE/confound 是主角）。
- **混合脊柱頁**（main-axis, loh-phasing, methyl-phasing）：段 4 主敘述 + 段 7 全開（claim + NEGATIVE + confound + circularity）。
- **名詞地基頁**（background-glossary）：段 2 是整頁主體，其餘段縮為「指向各概念頁」。

---

## §4 一致性規則（全站強制）

### 4.1 分層揭露 L0-L3

| 層 | 內容 | HTML 呈現 |
|---|---|---|
| **L0** | 一眼結論 / headline 數字 / takeaway | 永遠可見（`.bluf` / `.card`，置頂） |
| **L1** | 重點邏輯 / 直覺類比 / 並排對照 / 表格本體 | 預設可見（不折疊） |
| **L2** | 細節 / 嚴謹定義 / data_compare / 機制 | `<details>` 折疊（標題give一句） |
| **L3** | 溯源 / file:line / ledger entry / 公式推導 | `<details>` 折疊 or footer |
- **讀到夠就停**：每段 L0 必能獨立讀懂，不需展開 L2/L3。
- 對應 `feedback_reports_need_layered_disclosure`；多檔資料夾的 `INDEX` 本身 = L0+L1。

### 4.2 標記系統（雙軌：重要性 + 證據 tier）

- **重要性軌**（行內前綴）：🔴 必看/警告 · ⭐ 關鍵發現 · ◽ 次要/背景。
- **證據 tier 軌**（數字後綴 badge）：🟢 P（原檔對賬）· 🟡 P-caveat/S（口徑/二次紀錄，引用前查）· 🔵 framing（立場非硬數字）· 🔴 F（已捏造，禁出現）。
- 兩軌可疊加：`🔴 ⭐ ASM 真實但不可判別 TP/FP 🟢`。
- badge 必須對齊盤點 `key_data[].tier`，不可自行升級。

### 4.3 term 首次出現必定義

- 每個 term 在**其 home 概念頁**用 `N-DEF` 完整定義一次（建 anchor `#term-xxx`）。
- 其他頁首現 → `<abbr title="一句">term</abbr>` + 連回 home anchor，不重複完整定義。
- home 頁分配：NHD/HP1-1/PERMANOVA/Δ-vs-Δβ/cis-test/CramérV/epipoly → `ism-core`；HP-axis/ALLELE-axis/excess-over-null/ASM → `asm` 系；LOSO/circularity/vestigial → `negatives`；LOH/cnLOH/5mC-5hmC/unphase/HP3/anchor-AUC → `background-glossary` 或 `methyl-phasing`。
- §2.1 `N-DEF` 列了 23 個 anchor 名作為全站 term 註冊表。

### 4.4 provenance footer（每頁強制）

- 列 `source_paths[]`（絕對路徑，`/big7_disk/...`）。
- 列 `uncertainties[]` 為 `{{待查}}` 清單——誠實揭露未對賬項。
- 帶 commit hash + 生成日期。
- 聲明「本頁 Write 與產數字分析在不同 batch」（§13.0）。
- 引他處 validated 數字 → 加 `<!-- provenance-verified: 來源 -->`。
- **🔴 鐵則**：footer 列出的每個數字，必能在某 source_path grep 到；grep 不到 = 留 `{{待查}}`。

---

## §5 建議建置順序（前 7 頁）

> 排序準則：①支撐「當前主軸（Subclonal Reconstruction 論文）」最直接 ②節點複用率高（先建好被多頁引用的 anchor/SVG）③背景先於應用（名詞地基先行）。

| # | 概念頁 | thread 來源 | 主用節點 | 為何先做 |
|---|--------|------------|---------|---------|
| **1** | **背景名詞地基**（`background-glossary`） | background-glossary | `N-DEF`×~23 · `N-ROLE`(svg-5stage) · `N-EVO`(TO vs paired) | 🔴 **地基**：建全站 term anchor 註冊表 + 樣本/pipeline 速查，其餘 6 頁都連回這裡 |
| **2** | **ISM 五階段方法本體**（`ism-core`） | ism-core | `N-ROLE`(svg-5stage/single-locus) · `N-DEF`(Δ-vs-Δβ/NHD/PERMANOVA) · `N-SVG`(readxcpg-grid/dist-matrix) | 主軸引擎；建好 `svg-readxcpg-grid`/`svg-dist-matrix` 供多頁複用 |
| **3** | **HP-axis vs ALLELE-axis**（軸別 confound） | asm-cis-confound + asm-zar1l-brca2 | `N-DEF`(hp-axis/allele-axis) · `N-CONF`(baseline confound) · `N-SVG`(svg-hp-haplotype) | 🔴 全專案最易誤讀的混淆；`svg-hp-haplotype` 被 cis-test/LOH 頁複用 |
| **4** | **軸 C vs 軸 A**（ISM vs 業界） | ism-vs-external | `N-AXIS`(6軸+能力矩陣) · `N-EVO` · `svg-axisC-vs-axisA`(既有) | 論文 Related Work 防守；複用既有 SVG，產出快 |
| **5** | **四道 NEGATIVE**（filter DEAD） | negatives + main-axis | `N-NEG`×4 · `N-FUNNEL`(收斂) · `N-CMP`(in-dist vs LOSO) · `N-DEF`(loso/circularity) | 🔴 論文最強武器（防彈負結果）；建好 NEGATIVE 卡模板供 asm 頁複用 |
| **6** | **LOH-phasing 脊柱 + 循環性**（HD-1） | loh-phasing + main-axis | `N-CLAIM`(Grade B+) · `N-CONF`(by-construction circularity) · `N-STALE`(HPFineNGroups 更正) · `svg-circularity` | 主軸 positive headline；`svg-circularity` 是 main-axis 共用 |
| **7** | **normal-anchored cis-test**（去 confound 核心） | ism-core + asm-cis-confound + methyl-phasing | `N-DEF`(cis-test) · `N-CONF`(限 HCC1395) · `N-SVG`(svg-cis-3way) · `N-CEIL` | ISM 唯一可防守原創貢獻；三向比較需專屬 SVG |

**第 8 頁起（次優先）**：ASM 跨樣本復現（asm-zar1l-brca2 + cross-sample）· 甲基救 unphase（methyl-phasing 全線，含 T1/T2/T3 + chicken-egg）· ism-evolution 統計口徑（Fisher/beta-binomial/5mC-5hmC，偏方法學讀者）。

**第 0 步（建頁前一次性）**：建 `docs/explain/index.standalone.html`（全站 INDEX = `N-DASH` 儀表板 + term anchor 註冊表 + `N-REDFLAG` 讀本站守則）+ 抽共用 CSS 到 `docs/explain/_assets/explain.css`（沿用 G6 HTML 的 `--sp-*`/`.card`/`.badge`/`.funnel`/`details` 樣式）+ 用 `/methods-example` 一次生成待補 SVG primitive（svg-readxcpg-grid / dist-matrix / upgma-tree / hp-haplotype / cis-3way / 2x2-cochran），因為這些被多頁複用，先做攤平成本。

---

## §6 開放問題（需用戶 / 後續決策）

1. **explain 頁與既有 concepts/ HTML 的關係**：`docs/concepts/2026/06/` 已有 5+ 份 standalone HTML（G6 三層、研究史講解、CS 學生分層版）。explain 系列是**取代**它們、還是**並存（explain=穩定教學版 / concepts=時間戳快照）**？建議並存，explain 為去時間戳的 canonical 教學層，但需用戶確認。
2. **SVG 生成批次核准**：§5 第 0 步要一次生成 6 個 `svg-*` primitive（走 `/methods-example` data_ref 注入）。是否核准批次生成，或逐頁按需生成？
3. **PERMANOVA permutation 數歧異**：盤點 uncertainties 指出 99 vs 999 在不同檔不一致 → `N-ROLE`/`N-DEF` 標此參數前需 grep 源碼定案，這筆是否納入第 2 頁建置前置？
4. **term home 頁衝突**：HP-axis/ALLELE-axis 在 asm 三線都出現，home 頁暫定第 3 頁；但若第 1 頁背景地基也要收這兩個 term，需定「背景頁給 L1 直覺 + 連到第 3 頁 L2 嚴謹」的分工，避免雙定義 drift。
5. **tier badge 與 scientific-rigor L1-L5 的對應**：`N-CMP` 用 🟢P/🟡S（provenance tier），`N-CLAIM` 用 L1-L5（evidence ladder）——兩套分級系統是否要在 INDEX 給一張對照表，避免讀者混淆？

---

## §7 v2.1 風格定案（2026-06-12 用戶迭代後鎖定）⭐

> 第一批（INDEX + 名詞地基 + ISM）已建並經 2 輪用戶回饋，以下為**鎖定的全站風格**。新頁一律遵守；改樣式改 `_assets/explain.css` 一處即全站套用。

### 開放問題的決議

| # | 問題 | 決議（用戶確認）|
|---|------|------|
| 1 | explain vs concepts 關係 | ✅ **並存**：explain = 去時間戳 canonical 教學層 / concepts = 歷史快照（不動）|
| 2 | SVG 生成方式 | ✅ 風格探針階段**手繪 inline schematic**（快、可即時迭代圖例）；風格鎖定後再批量生 `_assets/svg/` primitive |
| 3 | PERMANOVA 99 vs 999 | ✅ **定案**：生產 99（`RegionProcessor.cpp:1573` 覆寫）/ 庫預設 999（`StructureTest.hpp:26`）。significance.csv = **59 欄**（非 117，實測）。見 `04_number_resolution.md` |
| 4 | term home 頁 | ✅ 名詞地基頁 = 全站 term SoT（28 核心詞 + anchor）；他頁首現連回，不重複定義 |
| 5 | 兩套分級 | 🟢P/🟡S（provenance）與 ⭐tier（成熟度）並用，badge 對齊盤點 `key_data[].tier`；L1-L5 留給未來 claim card |

### 用戶回饋落地的 4 條強制風格（v2 → v2.1）

1. **🔴 圖優先**：任何「重要敘述」都要配解釋圖（不只難點）。ISM 頁 3→7 圖。schematic SVG 用暖色盤、左右對照、箭頭標流向。
2. **資訊密集改表格**：發現清單、輸入/輸出、新舊對照、數據——能表格化就表格化（含 tier 欄、重要列 `td.imp` 粗體、斑馬紋）。
3. **重點提取盒 `.takeaway`**：每段最重要一句提取成醒目盒（accent/pass/dead 三色 + `.lab` 標籤）。頻率：每節 0-1 個，勿濫用。
4. **視覺層級**：h1 2.1rem / h3 1.62rem / `h4.mini`（橘色 ▎標）/ 行內 `.key`（黃 highlight）/ KPI 數字 1.6rem。

### 新增節點型 `N-WORKED`（統計算式 worked example）⭐ v2.1

- **用途**：方法學的<u>統計詞</u>（NHD/PERMANOVA/Cramér's V/Fisher/BH-FDR/AUC/excess-over-null…）必附「怎麼從數字算出來」的 worked example。源自用戶回饋「方法學的詞語的統計方式要補例子」。
- **結構**（`.worked` class）：`算式（附源碼 file:line）` → `.calc 示意數字代入` → `.res 結果` → 一句 caveat。
- **🔴 反捏造規則**：worked example 的**輸入數字為示意**（必標「示意」），但**算式與常數取自源碼**（NHD `DistanceMatrix.cpp:23` / pseudo-F `StructureTest.cpp:150` / Cramér's V `MathUtils.hpp:85`）。arithmetic 必須正確（如 Cramér's V [[30,8],[6,25]] → χ²≈24.3 → V≈0.59 已驗算）。
- **已套用**：名詞地基頁 §4（7 個 worked example）+ ISM 頁 §1/§5（3 個）。

### 共用 CSS 架構

- `_assets/explain.css` = canonical 樣式（改一處全站套用）。`index` + `01_background-glossary` 用 `<link>`；`02_ism-core` 暫 inline（視覺一致，待收斂為 link）。
- 新增 class：`.takeaway(.dead/.pass)` · `.worked(.lab/.calc/.res)` · `.key` · `h4.mini` · `.grid/.lcard`（INDEX 用）· `td.imp`。

### 第一批建置狀態（2026-06-12）

- ✅ `index.standalone.html`（10 研究線卡 + 閱讀順序 + 建置進度）
- ✅ `01_background-glossary.standalone.html`（28 詞 + 7 worked example）
- ✅ `02_ism-core.standalone.html`（9 段 + 7 圖 + 3 worked example）
- ✅ `_assets/explain.css`
- ✅ `03_methylation-read-filter.standalone.html`（2026-06-13；甲基讀取/篩選細節，3 SVG + 1 雙軌對照圖 + 3 worked example + 14 術語卡 + N-PIPE 模板首用；經 3-agent 對抗驗證 PASS/PASS/修正後）
- ⏳ 次批：ism-vs-external / 四道 NEGATIVE / LOH-phasing 脊柱 / normal-anchored cis-test（依 §5 順序）

---

## §8 新節點型 `N-PIPE`（演算法管線解説）⭐ 2026-06-13（03 頁沉澱）

> **緣起**：03 頁要把「ISM 怎麼讀甲基、怎麼篩」講到教授能聽懂的細節。這類「**程式管線 / 資料流**」解説與研究線敘述（N-EVO/N-CLAIM）不同 — 它的權威來自**源碼行號**而非實驗 tier。沉澱成可複用節點。

### N-PIPE 節點規格

- **用途**：解説一段確定性處理管線（解析 / 轉換 / 過濾），如 MM/ML→矩陣、距離計算、二值化。讀者要的是「每一步做什麼 + 門檻是多少 + 在哪行」。
- **結構（每個 PIPE 段）**：
  1. **流程 SVG**（左→右階段框 + 箭頭；每框一動作；底部 `svgnote` 標源碼行號）
  2. **逐步 worked example**（具體輸入 → 每步輸出；輸入標「示意」，常數/公式標 file:line）
  3. **I/O 表 或 門檻表**（含「來源/行號」欄，每列一個可 grep 的事實）
  4. **高認知負荷折疊**（`details.card.hcl`）放最易錯的一段機制（如反向股反掃）
- **🔴 反捏造鐵則（N-PIPE 專屬，延續 N-WORKED）**：
  - **每個門檻/常數/公式都必須親讀源碼**並在表格「來源」欄標 `檔案:行號`；不可憑記憶或舊 memory（舊 memory 的 file:line 可能 stale，如本輪「C_min=5」實為被覆寫的 struct 預設，生產值=3）。
  - **「宣告但未接線」要誠實標**：若 Config 有宣告但 grep 不到使用點（如 `min_site_coverage=5`），標「宣告未接線」而非假設有效（§13.0 範例）。
  - 數字歧異（如 99 vs 999、59 vs 117 欄）→ 採驗證值並在 footer 註明來源差異。

### N-FUNNEL 子型（過濾漏斗）

- **用途**：N-PIPE 中「逐關淘汰」的步驟（read QC、位點過濾）用漏斗 SVG（梯形遞縮）+ 編號守門清單（`.gate` ol，G1/G2…每關標源碼行）。
- **已套用**：03 頁 §④（七道守門漏斗）。

### 03 頁新增 CSS（本頁 inline；穩定後併入 explain.css）

- `.track2/.track(.proc/.app)`：兩種「篩選」左右對照軌道（處理層 info 藍 / 應用層 dead 紅）。
- `.gate`（編號守門 ol，G1/G2 方塊標號 + `.src` 行號）。
- **教學要點**：03 頁核心是「**一個詞兩個意思**」框架 — 「篩選」= 處理層（data cleaning，照做）vs 應用層（variant filter，DEAD）。此「同詞歧義拆解」框架可複用於其他易混詞（如「ASM 存在 vs 判別」、「phasing 訊號 vs 循環」）。

### N-PIPE 適用頁（未來）

- ism-evolution（統計口徑稽核：Fisher→beta-binomial、5mC/5hmC）← 強 N-PIPE 候選
- normal-anchored cis-test（d_cis / d_copy / d_drift 三向計算管線）
- methyl-phasing-assist（甲基模式→haplotype 投票管線）
