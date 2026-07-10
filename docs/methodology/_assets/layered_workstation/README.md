<!--
建立時間: 2026-07-07
更新時間: 2026-07-10
目標: 說明 canonical layered reconstruction 可攜工作站的證據邊界與重生方式
處理範圍: InterSubMod/docs/methodology/_assets/layered_workstation/
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md
-->

# layered_workstation/ — 分層 per-HP-家族樹枚舉重建工作站（每樣本一檔）

> **Current / canonical 分層模型**：樹 **per germline-HP-家族分開建**（家族優先於算法，避免 pooled allelic/clonal 混淆），每個 lineage 枚舉**所有等機率最小樹**。
> **可攜版**：每樣本一個 HTML + `index.html` 主頁；可離線開啟且零外部依賴。實際檔案大小、census 數字、資料版本與產生時間都由 builder 寫入主頁，不在 README 手動維護。

## 先理解證據邊界

- **回答什麼**：在 eligible multilocus region 內，read-level sSNV 共現資料允許哪些 per-HP-family 最小狀態樹，以及結果是唯一、歧義或 capped。
- **證據是什麼**：L0 germline HP family + L1 read-level sSNV 共現；L2 CN 只在資料可用時作事後裁決；V1–V7 顯示內部一致性。
- **不能宣稱什麼**：工作站不能單獨確認全腫瘤 clone phylogeny、細胞比例或生物 subclone。**L3 methylation 目前 pending，不參與樹排序或選擇。**

## 檔案
| 檔 | 內容 | 大小 |
|---|---|--:|
| **`index.html`** | Canonical landing：證據邊界、資料 provenance、L0–L3 狀態、跨樣本 census 與開啟連結。**先開這個。** | 由生成檔案決定 |
| 各樣本 `.html` | Dashboard + region browser + per-HP-family 枚舉樹 + V1–V7。 | 主頁逐列顯示實際大小 |

**每樣本檔內**：① dashboard 4 表（U1–U6 層級 / HP-multiplicity / region-determinacy / 比例字典）② **兩欄 region 瀏覽器**：左側 region 列表（determinacy、HP-multiplicity、chrom、排序、搜尋）→ 右側詳情 ③ 逐 germline-HP family 色框樹（HP1 藍 / HP2 粉 / HP3 青 / none 灰）④ **樹切換器**（前後切換 + thumbnail，一次一棵大圖，標示 exact trees 與 unique shapes）⑤ L0→L3 判斷軌跡 ⑥ V1–V7 驗證。

> **UI = topology_workstation 的兩欄瀏覽模式（左列表 + 右詳情 + 篩選 / 排序）＋ layered 資料模型（per-HP-family 枚舉樹 + 解集合切換）**。

## 資料模型（分子/分母定義 SoT）
`InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md`

- **主分母 = region**；region 確定 ⟺ 其**所有 census lineage** 都唯一樹。`hp_multiplicity=2` 是 span-two family census，可能含 root-only control，**不可直接解讀為 both-mutation-bearing**。
- **分母鐵則**：region% ≠ lineage-unit% ≠ 舊 pooled%（單位不同不可比）；跨樣本絕對數受資料量與 eligibility 影響，只作 inventory。
- Somatic sSNV、region、determinacy、recurrence 與驗證狀態全部由各樣本 `census` 物件生成。
- L2 / L3 的 available、unavailable 或 pending 狀態必須依資料顯示；不可把缺資料當 neutral，也不可用 L3 排樹。

## 重生
```bash
cd .../20260627_subclone_4axis_teaching/scripts
# Rebuild sample HTML files and the landing page.
python3 build_layered_per_sample.py

# Safely refresh only the landing page from existing census inputs.
python3 build_layered_per_sample.py --index-only
```
- 單樣本 HTML：`build_layered_workstation.py`（SM_RV=region_view.json SM_OUT=out.html）
- 中間格式：`build_region_view.py`（layered_reconstruction → region_view）
- 改 builder 後必做 Python compile、browser runtime smoke、390px / 1440px overflow 與截圖檢查。
- 輸入位置由 `build_layered_per_sample.py` 的 `SAMPLES` 設定；主頁顯示最新輸入 snapshot 時間。

## 與舊 topology_workstation/ 關係
舊資料夾 = pooled 單樹 + 三軸（c/branched/confirmed）模型；本資料夾 = current / canonical per-HP-family 枚舉模型。兩者可並存供方法史對照，但舊工作站不可當 canonical evidence source。
