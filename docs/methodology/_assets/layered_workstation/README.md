<!--
建立時間: 2026-07-07
更新時間: 2026-07-14
目標: 說明 canonical v5 layered reconstruction 全基因工作站的證據邊界、資料綁定與重生方式
處理範圍: InterSubMod/docs/methodology/_assets/layered_workstation/
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
  - InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md
-->

# layered_workstation/ — chr1–22 分層候選拓撲工作站（每 dataset 一檔）

> **Canonical v5**：唯一主骨幹是 LongPhase-S recalibrated `FILTER=PASS`。HP1／HP2 mutation-bearing units 分開重建；H3／H4／none／reference-only 只作 auxiliary/control。
> **可攜版**：7 個 dataset HTML + `index.html`；可離線開啟且零外部依賴。driver 會核對 current summary、`_SUCCESS`、verifier、code provenance 與所有輸入 SHA-256，任何 drift 都 fail closed。

## 先理解證據邊界

- **回答什麼**：chr1–22 的哪些 regional mutation-state candidate sets 完整；完整區的 exact candidate 組合數 `C` 與 unlabeled topology-shape 組合數 `Topo` 各是多少。
- **證據是什麼**：L0 HP family partition + L1 read-state constraints；candidate enumeration 不排序。Observed state、inferred hidden state、primary HP lane 與 auxiliary/control 分層顯示。
- **不能宣稱什麼**：候選樹不是已確認的 cell clone、真實 ancestry、細胞比例或全腫瘤 phylogeny。Read ALT fraction 不是 CCF。
- **Sidecar 邊界**：CN 只作 post-tree context；PS 只作 phase-block QC；L3 是 `not_evaluated / bounded_auxiliary`，只允許 negative screen / residual flag，禁止 tree ranking 或 lineage/clone confirmation。
- **PS 查詢限制**：canonical v5 有 cohort／dataset mixed-PS aggregate，但目前 region payload 沒有逐 region PS membership；因此頁面明示限制、不提供 PS filter，也不建立任何 PS-derived edge。若未來要逐 region 查詢，必須先擴充 upstream schema，而不能由 HTML 猜測。
- **Backbone sensitivity**：目前 verdict 是 `backbone_sensitive`；ClairS `FILTER=PASS` 只作 sensitivity，不可升格或混入 canonical 分母。全 cohort minimum retained-position Jaccard=`0.577`、primary-unit-key Jaccard=`0.474`、shared topology digest concordance=`0.936`，因此每頁都必須保留醒目的 sensitivity warning。

## 檔案
| 檔 | 內容 | 大小 |
|---|---|--:|
| **`index.html`** | Cohort command center：authority、W 漏斗、C/Topo composition、7 dataset 全基因入口、claim boundary 與收合 provenance。**先開這個。** | 由生成檔案決定 |
| 各 dataset `.html` | chr1–22 overview → region facets → verdict → observed evidence → primary HP candidate network → sidecars/raw drawer。 | 主頁逐列顯示實際大小 |

每個 dataset 頁把完整 detail 依 chromosome 嵌成 22 個 JSON chunks；初始只解析 lightweight all-genome index，選到 region 才解析該 chromosome。網路圖在手機維持最小可讀寬度並提供具名局部捲動區，不再把文字壓成不可讀尺寸。

- 座標搜尋支援單點與區間 overlap（例如 `chr8:34220481`、`chr8:34200000-34300000`）；染色體、C/Topo、read evidence、independent facet、查詢與 selected region 都寫進 URL hash，可直接複製目前檢視。
- Network 箭頭只表 mutation-state transition constraint 的方向；實線／虛線／點線分別表示 complete exact set 的 forced、candidate-variable 與未評估 edge，不等同已確認的細胞 ancestry。
- `.json` provenance 連結只在預設收合的「方法、驗證雜湊與原始資料連結」區出現，不插入主要閱讀數字。

## 資料模型（分子/分母定義 SoT）
`InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md`

- `W_tree = W_primary + no_primary_lineage`；`W_primary = complete + incomplete`。
- 完整區：`C_region = ∏ n_trees`；`Topo_region = ∏ n_distinct_shapes_exact`，且必須滿足 `C ≥ Topo ≥ 1`。
- Incomplete 不以 0 代替，C／Topo 顯示 unavailable；`display_trees_complete=false` 時不得從 stored prefix 宣稱 forced edge。
- `partial-only` 表示 overlapping partial reads 約束候選集合；不可解讀為 no-data 或 no-clone。
- 7 datasets 來自 6 biological samples；HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 dataset。

## 重生
```bash
cd InterSubMod
# Rebuild 7 hash-bound dataset pages and the landing page.
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py

# Rebuild only the landing page; fails if any existing dataset page is stale.
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py --index-only
```
- Current renderer：`build_layered_workstation_v5.py`；`build_layered_workstation.py` 是相容入口並委派給同一 fail-closed renderer。
- 中間格式：`build_region_view.py`（layered_reconstruction → region_view）
- 改 builder 後必做 Python/JavaScript syntax、canonical reconciliation、Chromium runtime、keyboard、零網路請求，以及 320 / 390 / 1440px overflow 與視覺截圖檢查。
- 唯一 authority 是 `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`；driver 不含可回退到歷史資料的 sample path。

## 與舊 topology_workstation/ 關係
舊資料夾是歷史 pooled／三軸展示；本資料夾是 current canonical v5 的 per-primary-HP candidate-set 工作站。兩者可並存供方法史對照，但舊工作站不可當 canonical evidence source。
