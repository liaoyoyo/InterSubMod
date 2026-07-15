<!--
建立時間: 2026-07-07
更新時間: 2026-07-15
目標: 說明 canonical v5 layered reconstruction 全基因工作站的證據邊界、資料綁定與重生方式
處理範圍: InterSubMod/docs/methodology/_assets/layered_workstation/
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology/current_v5_read_af_topology.index.json
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
  - InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md
-->

# layered_workstation/ — chr1–22 分層候選拓撲工作站（每 dataset 一檔）

> **Canonical v5**：唯一主骨幹是 LongPhase-S recalibrated `FILTER=PASS`。HP1／HP2 mutation-bearing units 分開重建；H3／H4／none／reference-only 只作 auxiliary/control。
> **可攜版**：7 個 dataset HTML + `index.html`；可離線開啟且零外部依賴。driver 會核對 current summary、`_SUCCESS`、verifier、code provenance 與所有輸入 SHA-256，任何 drift 都 fail closed。

## 先理解證據邊界

- **回答什麼**：chr1–22 的哪些 regional mutation-state candidate sets 完整；完整區的 exact candidate 組合數 `C` 與 unlabeled topology-shape 組合數 `Topo` 各是多少。
- **證據是什麼**：L0 HP family partition + L1 read-state constraints。Structural exact universe 由 current-v5 solver 以 `tree_cap=0` exhaustive 重算並對帳；HTML structural browser 仍可能只展開 stored display subset。Read-AF 對完整候選集合計分，頁面另顯示 ranking preview 與 co-top shape representatives。
- **不能宣稱什麼**：候選樹不是已確認的 cell clone、真實 ancestry、細胞比例或全腫瘤 phylogeny。Read ALT fraction 不是 CCF。
- **順位邊界**：排序所用 read-AF 是同一 primary HP family 內的 `ALT/(REF+ALT)`，不是 VCF caller AF、tumor-wide VAF 或 CCF。`score = Σ(AF_ancestor − AF_newly-acquired)` 只用於同一 HP unit 內排序，不得跨 unit／region 比較；`read-AF 唯一第一順位` 不等於「最可能真樹」，score 也不是 probability、posterior 或 calibrated confidence。`0/0` coverage 顯示 `N/A`，不會被當成 0%。
- **形態邊界**：單支、直系、旁系、直系＋旁系是 primary-HP mutation-state graph 的相容型態。直系以 depth≥2、旁系以 outdegree≥2 判定；ROOT 與 hidden nodes 都計入，HP 間只做 region-level OR、不建立跨 HP 親緣，也不等於 clone census。
- **Sidecar 邊界**：CN 只作 post-tree context；PS 只作 phase-block QC；L3 是 `not_evaluated / bounded_auxiliary`，只允許 negative screen / residual flag，禁止 tree ranking 或 lineage/clone confirmation。
- **PS 查詢限制**：canonical v5 有 cohort／dataset mixed-PS aggregate，但目前 region payload 沒有逐 region PS membership；因此頁面明示限制、不提供 PS filter，也不建立任何 PS-derived edge。若未來要逐 region 查詢，必須先擴充 upstream schema，而不能由 HTML 猜測。
- **Backbone sensitivity**：目前 verdict 是 `backbone_sensitive`；ClairS `FILTER=PASS` 只作 sensitivity，不可升格或混入 canonical 分母。全 cohort minimum retained-position Jaccard=`0.577`、primary-unit-key Jaccard=`0.474`、shared topology digest concordance=`0.936`，因此每頁都必須保留醒目的 sensitivity warning。

## 檔案
| 檔 | 內容 | 大小 |
|---|---|--:|
| **`index.html`** | Cohort command center：authority、W 漏斗、C/Topo composition、7 dataset 全基因入口、claim boundary 與收合 provenance。**先開這個。** | 由生成檔案決定 |
| 各 dataset `.html` | sample-wide 七組重點觀察 → GRCh38 座標比例分布 → chromosome grid → region facets → 四維 verdict → observed evidence → read-AF 第一順位 network + structural candidate browser → sidecars/raw drawer。 | 主頁逐列顯示實際大小 |

每個 dataset 頁把完整 detail 依 chromosome 嵌成 22 個 JSON chunks；初始只解析 lightweight all-genome index，選到 region 才解析該 chromosome。網路圖在手機維持最小可讀寬度、初始置中並提供具名局部捲動區；單一第一順位 shape 在桌面自動撐滿，兩棵以上才分欄。

- **GRCh38 全基因分布**：22 條染色體依 NCBI GRC 的 GRCh38 bp 長度縮放，每個 W_tree region 以 midpoint 落點；可切換 determinacy、read-AF selection、clone/subclone-compatible morphology、read evidence、primary HP、region size 與 CN region-sidecar 七種著色。CN 不是連續 segment track。
- **圖例多選**：點 A 顯示 A；再點 B 顯示 A∪B；再點已選類別會取消；沒有任何類別被選取就顯示全部。重點同一 mode 不會清空選取，真正切換 mode 才會清空，避免把不同 ontology 的 key 混在一起。
- **Sample-wide 七圖**：Topo count、C exact-candidate count、structural determinacy、read-AF selection、clone/subclone-compatible morphology、primary HP × H3 auxiliary、region retained-sSNV size。前五者以 `W_primary` 為分母；HP×H3 與 region size 以 `W_tree` 為分母。
- **舊名詞邊界**：舊 `single/linear/branched/star`、舊 observed ALT 群數 `c`、A/B/C/E determinacy 不直接搬入 current v5。新版只沿用可追溯的形態算法，對 current-v5 exhaustive candidates 重算單支／直系／旁系／直系＋旁系／未解；read-AF top representative 即使位於舊 stored top-32 外，也由 versioned sidecar 保留 edges，不由展示 prefix 猜測。
- **MAX_SNV=8**：region-size 圖的 8-site bin 可能同時包含 natural 8 與 cap-compressed groups；current region-view 無法再拆開，頁面不捏造兩者數量。
- 座標搜尋支援單點與區間 overlap（例如 `chr8:34220481`、`chr8:34200000-34300000`）；染色體、C/Topo、read-AF selection、morphology、read evidence、independent facet、查詢與 selected region 都寫進 URL hash，可直接複製目前檢視。
- Network 箭頭只表 mutation-state transition constraint 的方向；實線／虛線／點線分別表示 complete exact set 的 forced、candidate-variable 與未評估 edge，不等同已確認的細胞 ancestry。
- `.json` provenance 連結只在 index 的「機器證據與原始 JSON」或樣本頁的「方法、驗證雜湊與原始資料連結」預設收合區出現，不插入主要閱讀數字。

## 資料模型（分子/分母定義 SoT）
`InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md`

- `W_tree = W_primary + no_primary_lineage`；`W_primary = complete + incomplete`。
- 完整區：`C_region = ∏ n_trees`；`Topo_region = ∏ n_distinct_shapes_exact`，且必須滿足 `C ≥ Topo ≥ 1`。
- Incomplete 不以 0 代替，C／Topo 顯示 unavailable；`display_trees_complete=false` 時不得從 stored prefix 宣稱 forced edge。
- Read-AF selection 六類在 `W_primary` 內互斥守恆；morphology 五類也在 `W_primary` 內互斥守恆。無 primary region 只在 W_tree ideogram 顯示 `N/A`，不混進上述分母。
- `partial-only` 表示 overlapping partial reads 約束候選集合；不可解讀為 no-data 或 no-clone。
- 7 datasets 來自 6 biological samples；HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 dataset。

## 重生
```bash
cd InterSubMod
# Rebuild current-v5 exhaustive read-AF ranking + morphology sidecars.
python3 research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py \
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 \
  --input-manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json \
  --current-summary research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json \
  --method-script-dir docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts \
  --output-dir research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology

# Rebuild 7 hash-bound dataset pages and the landing page.
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py

# Rebuild only the landing page; fails if any existing dataset page is stale.
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py --index-only
```
- Current renderer：`build_layered_workstation_v5.py`；`build_layered_workstation.py` 是相容入口並委派給同一 fail-closed renderer。
- Freshness：sample page 除 summary / region-view / read-AF sidecar / sample marker 外，還必須符合 renderer SHA-256 與 `layered-workstation-v5-grch38-topology-multiselect-2` UI contract；`--index-only` 不會接受舊 renderer 產出的頁面。
- 中間格式：`build_region_view.py`（layered_reconstruction → region_view）
- 改 builder 後必做 Python/JavaScript syntax、canonical reconciliation、Chromium runtime、keyboard、零網路請求，以及 320 / 390 / 1440px overflow 與視覺截圖檢查。
- 持續 smoke：`python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_workstation_ui.py`。
- 本輪完整視覺稽核：`python3 research/20260715_layered_workstation_genome_topology_multiselect/scripts/audit_layered_workstation_playwright.py`（輸出 7 datasets × desktop/mobile 的 metrics、截圖與 JSON receipt）。
- 稽核解讀與證據鏈：`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/20260715_GRCh38拓撲形態多選工作站完整檢視_01.md`。
- Structural／W 指標 authority 是 `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`；read-AF selection 與 morphology authority 是與同一 run、summary SHA、region-view SHA 綁定的 current-v5 sidecar/index。Driver 不允許回退到歷史 ranking 資料。

## 與舊 topology_workstation/ 關係
舊資料夾是歷史 pooled／三軸展示；本資料夾是 current canonical v5 的 per-primary-HP candidate-set 工作站。兩者可並存供方法史對照，但舊工作站不可當 canonical evidence source。
