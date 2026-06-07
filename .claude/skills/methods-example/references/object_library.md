# 物件庫 — primitive taxonomy（可複用小物件）

> 一張方法解釋圖 = 這些小物件依案例組合。每個 primitive 有：**用途 / 參數 / data_ref 契約（verified vs schematic）/ SVG 輸出 / 擴充點**。
> 對齊業界先例：grammar-of-graphics 圖層（geom）、Vega-Lite spec-data 分離、NanoMethViz（read×CpG 為一級物件）、gggenomes（multi-track）、design-token（共用色標）。

## data_ref 契約（所有物件共用）

- 參數凡 `*_ref` 結尾 = `data.json` 的 dotted path。
- 指向 `verified.<k>` → 必須有非 null `value`（附 `src`），否則 renderer **refuse**（§13-A）。
- 指向 `schematic.* / synthetic` → 示意，允許無 src，但 renderer 打 watermark。
- 顯示用數字一律 verified；示意分子一律 schematic。**兩者不可互換**。

---

## P1 read×CpG matrix（read-level 甲基示意軌）
- **用途**：讓人看出 β＝某列某組 reads 的甲基比例；rows=reads（依 HP 分組），cols=CpG 座標，cell=甲基狀態。
- **參數**：`groups_ref`(schematic.groups: 每組 label/sublabel/color_key/reads) · `cpg_labels_ref` · `snv_after_ref`(SNV 插在第幾個 CpG 後) · `synthetic_note_ref`。
- **data_ref**：cell 狀態為 **schematic（synthetic）**；read 數 n 可引 verified（如 `n_reads_hp1`）。
- **SVG**：CpG tick 軸 + SNV 三角標 + 每 read 一條 line + 實心紅圓(甲基)/空心藍框圓(未甲基)（形狀冗餘 → 色盲友善）。
- **擴充**：多樣本 → 加 facet 軸；真實資料版 → 改吃 methylation matrix（標「真實資料」非示意）。

## P2 haplotype 分軌（normalHP1 / tumorHP1 / tumorHP1-1）
- **用途**：把 reads 依 haplotype/allele 分軌；**HP1-1 = HP1 germline 下 somatic-allele 子標，非第三 haplotype**；圖側若用 HP2-1 命名須標 label-flip。
- **參數**：`track_label` · `read_subset` · `color_key`。（目前折入 P1 的 groups 分組實作；獨立化時拆出。）
- **data_ref**：軌標籤靜態；read 子集示意。
- **擴充**：>3 軌時加 y 位 token；命名差異一律走 `label_flip_note`。

## P3 β bar（per-group 群率）
- **用途**：某 haplotype 群在某 CpG 集合的甲基比例 β（0-1，第 1 層平均）。
- **參數**：`bars:[{label, value_ref(→verified.beta_*), color_key}]`。
- **data_ref**：value **verified**（缺即 refuse）。
- **SVG**：水平 bar（寬 ∝ β）+ 數值標。

## P4 Δβ step（mean 主鍵 + max|Δ| 副鍵 + Wilcoxon）
- **用途**：Δβ = 平均(per-CpG Δ)（第 2 層平均）；同時保留 max|Δ| 副鍵抓交叉/位置型盲點。
- **參數**：`mean_ref`(→verified.dbeta_mean) · `n_ref`(n_shared_cpg) · `p_ref`(p_somatic) · `maxabs_ref`(verified.max_abs_delta，可 null)。
- **data_ref**：mean/n/p **verified**；max|Δ| 若 verified 無值 → 渲染成「方法說明」不顯示假數字。
- **SVG**：主鍵 Δβ + Wilcoxon p + 副鍵 max|Δ| note。**主+副鍵必並列**（spec ENG-1）。

## P5 normal-cis triplet（d_cis vs d_drift vs d_within）
- **用途**：加 normal 基線，拆「真 cis」與「drift」+ copy-partition 後純 allele cis（誠實標小）。
- **參數**：`normal_ref / tumorHP1_ref / tumorHP11_ref`(三條 β) · `dcis_ref/pcis_ref · ddrift_ref/pdrift_ref · dwithin_ref/pwithin_ref`。
- **data_ref**：全 **verified**。
- **SVG**：三條 β bar（normal/tumorHP1/tumorHP1-1）+ d_drift(小,NS)/d_cis(大)/d_within(小但真) 三行標註。

## P6 cluster↔hp + Cramér's V box（舊 ISM，C2 用）— *待擴*
- **用途**：10kb region + k=4 無監督分群 + cluster↔hp/alt 列聯 + Cramér's V + gating（舊法給「存在性 YES」但混淆/無方向/無 normal）。
- **參數**：`contingency_ref` · `v_hp_ref` · `v_alt_ref`(須含 0.51) · `n_reads_ref` · `n_cpg_ref`(須含 229) · `gating_ref`。
- **data_ref**：全 **verified**（⚠ BRCA2 舊 ISM 值在離線 `bip8_output_archive significance.json`；轉 validated 前須落檔到可 grep 源）。
- **擴充點**：在 `render_figure_spec.py` 加 `p_cluster_v_box()` 並註冊到 `RENDERERS`。

## P7 交叉型盲點情境（合成模擬，C3 用）— *待擴*
- **用途**：前後反向 pattern → Δβ mean≈0/Wilcoxon NS（主鍵漏）但 max|Δ| 大/ARI 抓得到（副鍵救）。
- **參數**：`pattern`(合成) · `mean≈0` · `p` · `maxabs` · `ari` — **全標 synthetic**，獨立 namespace + watermark，**不准混真值**。
- **擴充點**：加 `p_blindspot_scenario()`；合成數字不進 verified（無 src）。

---

## 新 primitive roadmap（從簡報慣例研究 2026-06-07，workflow wf_b6aa4c96-be9）

> 來源：2 份 PI 實驗室口試簡報萃取（詳見 `references/genomics_figure_conventions.md`）。**狀態（2026-06-07）：U1-U7 已實作**（`hap_split_track` / `igv_pileup` / `cpg_beta_matrix` / `grouped_bar` / `facet_grid` / `stacked_bar` / `loh_ideogram` / `readtrack_legend` 已註冊於 `tools/render_figure_spec.py`，對應 `examples/u1-u7/`，lint 全 PASS）。**仍待實作**：P-CALIB / P8 axis-dot / P9 allele-graph。色彩約定見 conventions §2（藍=germline/橘=somatic；**甲基化軸仍用 red/blue ramp 不套藍橘**）。分級顯示見 `references/detail_levels.md`。

| ID | primitive | 代表 | data_ref 契約 | 優先 |
|----|-----------|------|--------------|------|
| **P-READTRACK-LEGEND** | 右上角固定 4 色 legend（germline-HP1/HP2 / somatic / seq-error）| 全圖共用色 key | 靜態色 token（shared.color_scale）| 高（所有 P1-P7 共用）|
| **P-IGV-PILEUP** | IGV haplotype-lane read pileup（真實版 P2/P5）| N lane × M read `<rect>`，lane fill 編 HP，SNV 彩色 tick(A綠/C藍/G棕/T紅)，頂部 truth 軌 | reads/HP 來自真實 BAM（verified）| 高 |
| **P-FACET-GRID** | faceted small-multiples（metric × sample）| N×M panel，邊緣才畫軸，solid=method/dashed=baseline | 各 cell value verified + src | 高（跨 7 樣本 TP/FP）|
| **P-GROUPED-BAR** | 跨樣本 N-method 群組長條 | 每樣本 N bar + 末端 4dp + 旋轉 x 標 | verified + src；⚠ zoomed 非零軸**必標 caveat** | 高 |
| **P-CALIB** | predicted-vs-true scatter + y=x | 方形 panel + 灰虛 45° + 彩色 marker | verified pairs + src | 中 |
| **P-LOH-ideogram** | 單染色體 Giemsa 帶狀 + LOH interval 軌 | 帶狀 `<rect>` + centromere pinch + 對齊軌疊紫 interval | interval 座標 verified | 中（LOH 主軸）|
| **P8 axis-dot-pileup** | clipping count 沿基因座 + 視窗 bracket | y 軸 + 散點 `<ellipse>` + 逐窗算式 | count verified | 低（CNV/clipping）|
| **P9 allele-path-graph** | rr/ra/ar/aa node-edge 圖 | `<ellipse>` node(R/A) + `<line>` edge + 四色 | DECK2 caller 專屬 | 低（本專案複用低）|

**layout 件（非物件）**：`P-LEGEND-STRIP`（多 panel 共用 detached legend）、`stage-flow-connector`（圓角框+箭頭+算式串 pipeline）、`brace-aggregator`（花括號收束 → 純量）。

**強化既有**（同來源）：P1 加 sliding-window `W` bracket + `✕` 判定；P2 → 兩層階層（germline 上/somatic 下 + 分隔線）+ read 改五邊形箭頭 + ATGC 子格 + HP label；P3 → 雙色分組 + 末端數值標籤 + zoomed-axis（含 caveat）。

---

## 如何定義一個新 primitive（通用化路徑）

1. 在本檔加一段：用途 / 參數（哪些 `*_ref` 指 verified、哪些示意）/ SVG 概念。
2. 在 `tools/render_figure_spec.py` 寫 `def p_<name>(data, prim, sh, y) -> (svg, height)`，用 `req_val()`(verified, 缺即 refuse) / `opt_val()`(可空) / `resolve()`(schematic)。
3. 註冊到 `RENDERERS` dict + 在 `figure_spec.schema.json` 的 `type` enum 加值。
4. 在 `case_templates.md` 標哪些 case 會用它。
