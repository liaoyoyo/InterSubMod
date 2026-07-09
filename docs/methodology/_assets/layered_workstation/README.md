# layered_workstation/ — 分層 per-HP-家族樹枚舉重建工作站（每樣本一檔）

> **新分層模型**（2026-07-07；取代舊 `topology_workstation/` 的 pooled 三軸模型）：樹 **per germline-HP-家族分開建**（家族優先於算法，修 allelic/clonal 混淆），每 lineage 枚舉**所有等機率最小樹**。
> **可攜版**：每樣本一檔（3–31MB）+ `index.html` 主頁，任何電腦可離線開、零外部依賴。

## 檔案
| 檔 | 內容 | 大小 |
|---|---|--:|
| **`index.html`** | 🏠 主頁：7 樣本跨樣本表（somatic sSNV / multi-HP% / region-determinacy / avg樹）+ 開啟連結。**先開這個。** | ~10 KB |
| `HCC1395.html`（SEQC2 主樣本）/ `COLO829` / `H1437` / `H2009` / `HCC1395_DORADO` / `HCC1937` / `HCC1954` | 各樣本工作站 | 3–31 MB |

**每樣本檔內**（2026-07-09 改 topology 式兩欄版面，UI 更清楚）：① dashboard 4 表（U1-U6 層級 / HP-multiplicity / region-determinacy / 比例字典）② **兩欄 region 瀏覽器**：左側緊湊 region 列表（篩選 region-determinacy·HP-multiplicity·chrom·**排序**·搜尋）→ 點一個區 → 右側詳情面板 ③ 右側逐 germline-HP 家族色框樹（hp1 藍 / hp2 粉 / hp3 青 / none 灰）④ **樹切換器**（◀▶ + thumbnail 跳轉，一次一棵大圖，標「N 棵等機率樹 = M 種形狀」）⑤ L0→L3 判斷軌跡 ⑥ V1-V7 驗證。

> **UI = topology_workstation 的清楚兩欄格式（左列表+右詳情+篩選/排序）＋ 本 layered 資料模型（per-HP-家族枚舉樹 + 左右切換）**（2026-07-09 合併;使用者要求）。

## 資料模型（分子/分母定義 SoT）
`InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md`
- 🔴 **主分母 = region**；region 確定 ⟺ 其**所有** germline(1/2) lineage 都唯一樹（多-HP 需雙確定）。
- 🔴 分母鐵則：region% ≠ lineage-unit% ≠ 舊 pooled%（單位不同不可比）；跨樣本只比比例（絕對受深度混淆）。
- 🔴 somatic sSNV = census somatic==True（HCC1395=23,810，非總 census 35,332）。
- 6 樣本 cn=unknown → L2 CN 通道不可用（只 HCC1395 有 SEQC2 整數 CN）。

## 重生
```bash
cd .../20260627_subclone_4axis_teaching/scripts
# 全 7 樣本(需各 layered_region_view_{S}.json 已在)
python3 build_layered_per_sample.py
# 6 樣本新資料重跑(讀 tagged BAM,~2hr)
bash run_layered_6samples.sh      # → mlhp + layered + region-view
```
- 單樣本 HTML：`build_layered_workstation.py`（SM_RV=region_view.json SM_OUT=out.html）
- 中間格式：`build_region_view.py`（layered_reconstruction → region_view）
- ⚠ 改 build JS 後必 `node --check` + runtime smoke（曾有 runtime bug 讓白畫面）
- HCC1395 資料=`20260618_subcluster_pilot/`；6 樣本=`multisample_subclone/{S}/`

## 與舊 topology_workstation/ 關係
舊資料夾 = pooled 單樹 + 三軸（c/branched/confirmed）模型；本資料夾 = 新 per-HP-家族枚舉模型。**兩者並存**（舊未刪，供對照）；新模型修正舊 pooled 的 allelic/clonal 混淆（多-HP 56-91% 區被舊法混家族）。
