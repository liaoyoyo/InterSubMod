# topology_workstation/ — [ARCHIVED] 舊克隆樹拓樸歷史快照

> # ⛔ DEPRECATED（2026-07-10 正式停用）
> **此工作站建於舊 `is_somatic` 骨幹**（normal VAF<5% 粗重檢，已證誤殺 429 SEQC2-TP 真 somatic），
> 已由 **ClairS v0.4.0 paired PASS** 新骨幹的**分層樹枚舉工作站**取代。
> **canonical = `../layered_workstation/`**（driver `build_layered_per_sample.py` / 單樣本 `build_layered_workstation.py`）。
> layered 版特點：per germline-HP 家族分建樹（修 allelic/clonal 混淆）+ 移植本站全部可觀察面板
> （位點證據改用原始 col_coverage / locus / pairwise 共現 / 前後相鄰區 / 每家族觀測 read）
> + 新增**位點層實測 vs 組合層(phasing)推斷**的正確區分（`n_full_cov_reads` co-phase 缺口）。
> **本資料夾數字為舊骨幹，僅供歷史對照，勿引用。** 不再接受研究更新；只有在修復封存介面或驗證可重現性時，才可由原始資料重生。
>
> ---
>
> **（以下為停用前原始 README，歷史保留）**

## 封存使用規則

- `index.html` 是 archive landing page；舊樣本表預設收合，主要操作會導向 `../layered_workstation/index.html`。
- 所有 7 個 direct sample HTML 都必須在首屏顯示 `ARCHIVED / DEPRECATED`、舊 `is_somatic` backbone、不可引用與 canonical layered 連結，不能只依賴 index 警示。
- 樣本頁只保留歷史篩選與查閱；人工 review、瀏覽器儲存與判讀匯出已停用，避免舊決策跨樣本或跨版本污染。
- 允許用途：重現歷史畫面、理解方法沿革、比對介面演進。禁止用途：validation evidence、目前跨樣本結論、paper 引用與新人工判讀。

> **這是「可攜版」工作站**：把原本 63.6MB 的單檔拆成「每樣本一檔(4–21MB) + 主頁 index.html」，
> 讓**任何電腦**（低 RAM / 一般瀏覽器）都開得動、也不會被 email/Slack 的 25–50MB 上限截斷。
> 產生日：2026-07-06。

## 檔案內容
| 檔 | 內容 | 大小 |
|---|---|--:|
| **`index.html`** | 🏠 **主頁**：7 樣本跨樣本摘要表(總區/c分布/incompatible%/branched%/confirmed%/檔案大小)+ 每樣本「開啟」連結 + 欄位意義 + 相關文件。**先開這個。** | ~5 KB |
| `HCC1395.html` | 凍結主樣本(SEQC2 truth) | ~12 MB |
| `COLO829.html` / `H1437.html` / `H2009.html` / `HCC1937.html` / `HCC1954.html` / `HCC1395_DORADO.html` | 各樣本完整工作站(單分頁) | 4–21 MB |

**用法**：整個資料夾複製到任何電腦 → 開 `index.html` → 點樣本「開啟」。每檔獨立、離線可用、零外部依賴。

## 為何要拆（原單檔的問題）
原 `docs/methodology/_assets/20260629_multisample_topology_workstation.standalone.html` 已達 **63.6MB，其中 99.96% 是 7 樣本 inline JSON**。換到低 RAM 電腦，瀏覽器一次 parse 63MB inline `<script>`（需 300–600MB+ RAM）會 hang/白畫面/崩潰；且超過多數傳輸管道上限會被截斷成壞檔。拆每樣本後最大也只 ~21MB，安全。

## 歷史重現（非研究更新）

**產生程式（唯一 SoT，不變）**：`../20260627_subclone_4axis_teaching/scripts/build_topology_workstation.py`（單樣本邏輯共用）
**拆檔 driver（本資料夾的產生器）**：`../20260627_subclone_4axis_teaching/scripts/build_per_sample_workstation.py`

**重生本資料夾**（僅限 archive reproduction / 封存介面修補）：
```bash
cd .../20260627_subclone_4axis_teaching/scripts
python3 build_per_sample_workstation.py   # → 覆寫本資料夾全部 8 檔
```
- driver 複用 build_topology_workstation.py 的 `SM_SAMPLES`(子集)+`SM_OUT`(輸出路徑) env，對 7 樣本各跑一次 → `{sample}.html`，再產 `index.html`。
- 資料來源同單檔：HCC1395=凍結 `data/`；6 樣本=`/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{S}/`。改資料就重跑上游 pipeline 後再跑本 driver。
- **只可修封存顯示/安全性**：改 `build_topology_workstation.py`（單樣本+單檔共用同一份 JS）→ 跑 driver 重生本資料夾。新研究功能一律改 layered generators，不得加回人工 review。

> ⚠ **改 build script 後必 `node --check` 或 playwright 驗**：2026-07-06 曾因(a) 並行 session 加的 m-通道三元運算子缺 `)`、(b) 三軸面板 `joint[c]` 未防 n_clusters=0 區 → 整站白畫面。兩 bug 已修（見驗證記錄 log）。

## 相關文件
- 三軸分析 + 定義：`InterSubMod/docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md`
- 對抗驗證記錄：`InterSubMod/docs/methodology/20260702_topology_three_axis_verification_log_01.md`
- 建置地圖：`InterSubMod/docs/methodology/20260701_multisample_topology_workstation_build_map_01.md`
