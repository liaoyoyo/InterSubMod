---
name: methods-example
description: |
  方法解釋圖「生成 + 驗證 + 圖例細節迭代修正」的整體 skill。通用、可複用的「物件庫(primitive) + 案例模板(case) + data_ref 注入 + 圖例細節 verify loop」架構，產論文/報告級方法解釋圖（read-level haplotype 示意 / 新舊方法互補對照 / 統計流程示意）。圖 = 宣告式 figure_spec.json（依目標組哪些小物件）+ data.json（verified 真值附 src ＋ schematic 示意，物理分離）→ render_figure_spec.py 注入（缺 verified 真值直接 refuse ＝ §13-A 構造防捏造）。
  USE WHEN：「方法解釋圖」「method figure」「解釋這個分析方法」「read×CpG/haplotype 示意圖」「新舊方法對照圖」「論文 methods 圖」「把方法畫給 PI/教授看」「Δβ/cis-test/分群怎麼運作的圖」、為「分析方法」（非資料結果）產生可重複的解釋例子或圖、要持續檢視/修正既有方法圖的圖例細節時。
  SKIP WHEN：真實資料統計圖（用 InterSubMod/scripts/lib/plot_setup.py matplotlib + 人眼驗證）、generic 概念/流程/icon（用 image-gen）、PNG 品檢（用 image-vision-check）、純文字方法敘述無圖（用 narrative-frame）、PPTX 內已 review 過的圖微調、無嵌入數字的純裝飾 icon。
---

# methods-example — 方法解釋圖 生成+驗證+迭代

把「怎麼產生並驗證一個論文/報告級的方法解釋圖」做成**可重複、可複用、物件組合式**的流程。核心信念：**一張方法圖 = 小物件(primitive) 依解釋目標(case) 的組合，真值只能從 data_ref 注入（無法手打捏造），圖例細節靠 verify loop 持續修正。**

## Phase & Chain Position

- Stand-alone skill；retrospective 報告/簡報階段觸發（P5-P6 或對 PI 解釋方法時）。
- **Upstream**：`narrative-frame`（N1-N4 決定要講哪個方法、用哪個教學框架 Feynman/Diátaxis）→ 本 skill 把敘事落成圖。分析真值來自既有研究 JSON（如 `research/**/genome_survey_v2/*.json`）。
- **Downstream**：產出 SVG 可嵌入 `html-report-build`（standalone/report）、`pptx-build`、validated/PI 報告。
- **與 image-gen/image-vision-check 互補不重疊**：image-gen = generic 概念/流程/icon → PNG；本 skill = 領域特定方法 explainer（手控 SVG + 內嵌已驗證 metric + 示意分子）。PNG 佐證子圖才回 image-vision-check 6 維；SVG explainer 走本 skill 的 12-criteria + SciFig R1-R6。

## Dependencies

| 方向 | 對象 | 用途 |
|---|---|---|
| **Uses** | `tools/render_figure_spec.py` | 讀 spec+data → 注入真值 → 出 SVG；缺 verified 真值 refuse；12 primitive（hap_split_track / igv_pileup / cpg_beta_matrix / grouped_bar / facet_grid / stacked_bar / loh_ideogram / readtrack_legend / read_cpg_matrix / beta_bar / dbeta_step / normal_cis_triplet）|
| **Uses** | `tools/lint_figure.py` | pixel-free 設計自檢（字型/溢出/認知負擔/label-flip/色約定/示意）→ render→lint→fix 自我迭代 |
| **Uses** | `tools/build_gallery.py` | 把 examples/*/method_explainer.svg 併成單頁 `_gallery.html` 供一次眼驗 |
| **Uses** | `templates/figure_spec.schema.json` | spec 格式（case + primitives + shared） |
| **Reads** | `references/detail_levels.md` | 分級顯示（L0 必要/L1 輔助/L2 on-demand/L3 隱藏）+ detail 參數 + 字型對齊目的 |
| **Reads** | `references/object_library.md` | 物件 primitive 定義（P1-P7 + 如何擴新物件） |
| **Reads** | `references/case_templates.md` | 案例模板（C1-C3 + goal→objects 組合法） |
| **Reads** | `references/explainer_checklist.md` | 12-criteria + 圖例細節迭代 verify + SciFig R1-R6 |
| **Reads** | `references/genomics_figure_conventions.md` | 名詞 glossary + 色彩約定(藍 germline/橘 somatic；甲基用 ramp) + 字型/排版 16:9 模板 + 切格組合 + 新 primitive roadmap（從 2 份 PI 簡報萃取，wf_b6aa4c96-be9）|
| **Reads** | 研究 JSON（`research/**/*.json` 等） | verified 真值來源（每數字附 `src`） |
| **Used-by** | `html-report-build` / `pptx-build` | 消費產出的 SVG |
| **Writes** | `<case_dir>/method_explainer.svg`(+.html) | 渲染輸出（非 SoT，可重生） |

## 架構（你要的「可複用高架構 + 物件定義 + 依目標組合」）

```
figure_spec.json                      data.json
  case_template_id  ─────┐             ├─ verified.<k> = {value, src}   ← 真值, 附可 grep 來源
  primitives:[{type,params,*_ref}] ◀───┤
  shared:{color, glyph, coord}         └─ schematic.* (synthetic:true)  ← 示意, 打 watermark
          │                                          │
          ▼  render_figure_spec.py（物件 dispatch）   ▼
   每 primitive 從 data_ref 拉值 ── verified 缺 key/null → REFUSE(exit 3) ＝ §13-A 構造防捏造
          │
          ▼
   method_explainer.svg  +  provenance 稽核表（metric → value → src）
```

- **物件 primitive（小物件）**：P1 read×CpG 矩陣 / P2 haplotype 分軌 / P3 β bar / P4 Δβ step(mean 主鍵+max|Δ| 副鍵) / P5 normal-cis 三條 / P6 cluster↔hp + Cramér's V box(舊法) / P7 交叉型盲點(synthetic)。定義見 `references/object_library.md`，並含**如何為新領域/新方法擴一個新 primitive**。
- **案例模板 case（依目標組合）**：C1 方向型 anchor 對照 / C2 新舊方法並列 / C3 盲點合成模擬。挑哪些物件、怎麼組，由「解釋目標」決定 → 見 `references/case_templates.md` 的 goal→objects 對照。
- **真值物理隔離**：所有要顯示的聚合數字走 `verified.<k>`（附 src），renderer 缺值即 refuse → **「每數字可 grep」從紀律升級為架構性質**；示意分子走 `schematic.*` + synthetic watermark，與真值不混。

## 6 步流程（LOCK 先於 GENERATE，對齊 §13.0 物理隔離）

1. **LOCK-AND-GATHER** — 鎖定要講的方法 + 一個真實 anchor（位點/案例）+ 解釋目標；**先把要嵌的每個數字落檔→Read 回確認非 error/未完成**，寫進 `data.json` 的 `verified.<k>`（附 `src`）。示意分子寫 `schematic.*`(synthetic)。**不可手打預期值**。
2. **GENERATE** — 依目標選 case + primitives 寫 `figure_spec.json` → 跑 `render_figure_spec.py`（缺 verified 真值會 refuse）。
3. **VERIFY-12** — 對 `references/explainer_checklist.md` 12 條（BIO/ENG/CLR/HON/CMP）逐項判 PASS/PARTIAL/FAIL + 跑 `number_provenance_check.sh`（HTML 嵌入時）。
4. **FIX**（minimally）— 補缺口；⚠ **「不在白名單」≠「刪已驗證真值」**（§13 要 traceable 不是 delete）；分析 Bash 與 spec/data Write 不同 batch。
5. **RE-VERIFY** — 複跑 12 條 + SciFig R1-R6；先前 FAIL 已解或誠實降為示意。
6. **FINALIZE** — caption 明標示意 vs 真值；轉 validated/PI 前確認 frontmatter `data_sources:` 齊（否則 `number_provenance_check.sh` 對 validated 路徑 exit 2）。

> **圖例細節迭代 verify（持續檢視→驗證→修正）**：每輪改完重跑 renderer（diff SVG），對照 checklist 的「圖例/legend 細節」段（label-flip 註明、單位、色標/形狀冗餘一致、示意標註、anchor 跨圖一致、符號方向不畫反）。圖例層問題（非結構）就地修 data/spec 再重生，不必重畫。

## 通用化（不只這張 BRCA2 圖）

- **新位點/新樣本**：只換 `data.json` 的 verified 值（附新 src）+ schematic pattern，spec/renderer 不動。
- **新解釋目標**：在 `case_templates.md` 加 case（選既有 primitives 組合）。
- **新物件型別**：在 `object_library.md` 定義新 primitive（params + data_ref 契約 + SVG 片段），在 renderer `RENDERERS` 註冊一個 render 函式。P6/P7（C2/C3 用）即此擴充路徑。

## Failure Mode & Diagnostics

| 症狀 | 可能原因 | 修法 |
|---|---|---|
| `[REFUSE] value is null` exit 3 | verified 真值缺/未 LOCK-AND-GATHER | 先落檔真值 + Read 回 + 填 `verified.<k>`；**勿手打** |
| `missing data path` | spec 的 `*_ref` 打錯 dotted path | 對 `data.json` 鍵名核對 ref |
| `primitive type ... 尚未實作` | 用了 P6/P7 但 renderer 未註冊 | 依 object_library.md 擴一個 render 函式 |
| 圖例 label-flip 漏標 | HP1-1≡HP2-1 未註明（常見領域捏造風險） | data.json 補 `label_flip_note` + spec annotations 引用 |
| 數字與報告其他處不一致 | 不同口徑並列（cis-test 0.25 vs deepdive 5mC-only 0.627） | 同一圖只用同口徑；src 標明出處 |
| SVG 看起來擠/重疊 | primitive 高度估算偏差 | 調 render 函式回傳 height 或 shared.width；無 rasterizer 時用瀏覽器開 .html 目視 |

## Worked example（已驗證 pilot）

`examples/c1_brca2_dbeta/`：BRCA2 chr13:32,315,128 Δβ somatic-cis。
- `data.json`（11 verified 真值，src 指向 `control_cohesion_cistest.json` / `brca2_pertag_cohesion.json` / `copy_partition_confirm.json`）+ `figure_spec.json`（C1，P1+P3+P4+P5）。
- 跑：`python3 tools/render_figure_spec.py examples/c1_brca2_dbeta/figure_spec.json --html --audit`
- 輸出 `method_explainer.svg`(+.html)；audit 表列每 metric→src。REFUSE 行為已測（beta value=null → exit 3）。

## 限制

- renderer 目前實作 P1-P5（C1 覆蓋）；P6/P7（C2 新舊對照 / C3 盲點）為待擴介面。
- 無內建 SVG→PNG rasterizer；像素級 vision 檢查需瀏覽器開 .html 或外部轉換器（rsvg-convert/cairosvg）。
- 程式產 SVG 的視覺精緻度 < 手工微調；目標是「data-binding 反捏造 + 可重複」，極致排版仍可事後手調 SVG（但會失去 data-binding，建議只調 shared token）。
