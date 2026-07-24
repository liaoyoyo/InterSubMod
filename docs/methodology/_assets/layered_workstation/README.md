<!--
建立時間: 2026-07-07
更新時間: 2026-07-24
目標: 說明 exact-PS×HP layered reconstruction 全基因工作站的證據邊界、資料綁定與重生方式
處理範圍: InterSubMod/docs/methodology/_assets/layered_workstation/
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology/current_v5_read_af_topology.index.json
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
  - InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md
-->

# layered_workstation/ — exact-PS chr1–22 分層拓撲工作站

> **Primary authority（2026-07-24）**：exact nonmissing PS × primary HP1/HP2 ×
> threshold-3 strict endpoint read-linkage × bounded block。預設 builder 不再以
> `adjacent_gap<=50000` 傳遞合併。
>
> **可攜版**：7 個 dataset HTML + `index.html`；可離線開啟且零外部依賴。
> driver 會核對 cohort/pipeline/MLHP/topology/census receipt、SHA-256、逐列 join
> 與分母守恆，任何 drift 都 fail closed。

## 新版數據範圍與正確分母

- 98,955 final exact-PS topology groups。
- 85,941 mutation-bearing groups。
- 71,955 `ranked_complete` groups；只在此分母計算 determinacy/signature。
- 10,717 `ABSTAIN_RESOURCE_LIMIT`、3,224 zero denominator、45
  recurrence-required，均保留為獨立狀態。
- 680,527 棵 AF-global-best trees 已作 exact signature census：
  - `UNIQUE_TREE` 39,648（55.1011%）；
  - `TIED_SAME_TOPOLOGY` 23,858（33.1568%）；
  - `TIED_CROSS_TOPOLOGY` 8,449（11.7421%）；
  - 單一 root-preserving unlabeled exact shape 63,506（88.2579%）。
- 7 technical datasets = 6 biological samples；HCC1395 與
  HCC1395_DORADO 是同一 biological sample 的 technical comparison。

**Claim ceiling**：read-AF 只排序 recurrence-allowed minimum mathematical
graph models；不是 caller VAF、CCF、calibrated posterior、cellular clone 或真實
ancestry。7/7 technical PASS 不等於 topology-complete。

## 檔案

| 檔 | 內容 | 大小 |
|---|---|--:|
| **`index.html`** | Exact-PS cohort command center：authority、分母漏斗、determinacy、coarse geometry、7 dataset 入口、HCC1395 技術比較與收合 provenance。**先開這個。** | 約 22 KB |
| 各 dataset `.html` | 完整 compact final-group index、GRCh38 canvas、多選/搜尋、exact PS/HP boundary、signature census、代表最佳樹與 read-ALT 表。 | 約 9–89 MB |

## 新版互動與資料模型

- **GRCh38 全基因分布**：22 條染色體按 GRCh38 bp 長度縮放；全部
  98,955 final groups 依 midpoint 落點。
- **圖例多選**：點 A 選 A，再點 B 為 A∪B，再點已選類別取消；沒有選擇
  時顯示全部。切換 ontology 會清空前一 ontology 的 selection。
- **著色維度**：Determinacy、拓撲形態、HP family、判定狀態、active k。
- **座標搜尋**：支援單點與區間 overlap，例如
  `chr10:87818272-87928739`。搜尋只回傳 final groups；source singleton
  abstain 不會假裝成 topology=0。
- **Region grain**：
  `exact_ps_hp_component_bounded_block`；同一 exact PS、同一 primary HP，
  且 threshold-3 strict endpoint read-linkage 支持。PS 是 boundary，
  不是 ancestry edge。
- **Determinacy**：`UNIQUE_TREE`、`TIED_SAME_TOPOLOGY`、
  `TIED_CROSS_TOPOLOGY` 只以 71,955 ranked units 為分母。
- **形態**：Single-only、Direct-only、Sister-only、Sister+direct 是
  rooted mutation-state graph geometry；不是 clone 類別。
- **Network**：只畫 frozen topology row 的 deterministic global-best
  representative。Census 沒有逐棵 edge list；tied rows 以 signature/count
  顯示，頁面不捏造完整 edge union。
- **Provenance**：所有 `.json`/receipt 路徑只在預設收合的證據區出現，
  不插入主要數字敘事。

## Legacy canonical-v5（僅歷史重現）

舊頁面使用 `adjacent_gap<=50000` 傳遞合併，以及
`W_tree/W_primary/C/Topo` 等不同 grain 與分母。它只能由 `--legacy-v5`
重現，不能與新版 98,955/71,955 分母、88.2579% exact-shape 結果混用。

## 重生
```bash
cd InterSubMod
# 預設：重建 exact-PS authority 的 7 個樣本頁與 cohort index。
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py

# 等價的 direct primary renderer。
python3 \
  research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py

# 只重建 index；7 個既有頁必須已綁定同一 census receipt。
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py --index-only

# 資料與跨面板 regression。
python3 -m unittest \
  research/20260724_exact_ps_cpp_topology_signature_census/tests/test_exact_ps_layered_workstation.py

# Chromium 桌機 / 手機真實互動與視覺稽核。
python3 \
  research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py

# 僅需重現舊 canonical-v5 / 50-kb 頁面時才使用。
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py --legacy-v5
```
- Primary renderer：
  `research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py`。
- Legacy renderer：`build_layered_workstation_v5.py`；只由 `--legacy-v5` 觸發。
- Freshness：sample page 必須符合 `intersubmod-authority =
  20260724-exact-ps-hp-strict-read-linkage`、同一 census receipt SHA 與
  `layered-workstation-exact-ps-v1` UI contract；`--index-only` 不接受舊頁。
- Primary authority 是 2026-07-24 topology cohort receipt/all7 summary、
  per-sample MLHP/topology outputs，以及 signature `receipt.v2.json/summary.json`。
- 改 builder 後必做 syntax、跨面板 regression、Chromium runtime、
  多選/再次取消、座標搜尋、零網路請求、390/1440px overflow 與視覺截圖。
- 本輪 audit receipt：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit/receipt.json`。

## 與舊 topology_workstation/ 關係
舊資料夾是歷史 pooled／三軸展示；本資料夾是 exact-PS × primary-HP
regional mathematical topology 工作站。兩者可並存供方法史對照，但舊工作站
不可當新版 canonical evidence source。
