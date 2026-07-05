# topology_workstation/ — 多樣本克隆樹拓樸工作站（每樣本一檔）

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

## 🔧 日後改動與資料 → 都在新位置處理

**產生程式（唯一 SoT，不變）**：`../20260627_subclone_4axis_teaching/scripts/build_topology_workstation.py`（單樣本邏輯共用）
**拆檔 driver（本資料夾的產生器）**：`../20260627_subclone_4axis_teaching/scripts/build_per_sample_workstation.py`

**重生本資料夾**（資料更新 / 改了 build script 後）：
```bash
cd .../20260627_subclone_4axis_teaching/scripts
python3 build_per_sample_workstation.py   # → 覆寫本資料夾全部 8 檔
```
- driver 複用 build_topology_workstation.py 的 `SM_SAMPLES`(子集)+`SM_OUT`(輸出路徑) env，對 7 樣本各跑一次 → `{sample}.html`，再產 `index.html`。**不改 build_topology_workstation.py。**
- 資料來源同單檔：HCC1395=凍結 `data/`；6 樣本=`/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{S}/`。改資料就重跑上游 pipeline 後再跑本 driver。
- **改顯示/功能**：改 `build_topology_workstation.py`（單樣本+單檔共用同一份 JS）→ 兩者一起生效 → 跑 driver 重生本資料夾 + 跑 `build_topology_workstation.py` 重生單檔。

> ⚠ **改 build script 後必 `node --check` 或 playwright 驗**：2026-07-06 曾因(a) 並行 session 加的 m-通道三元運算子缺 `)`、(b) 三軸面板 `joint[c]` 未防 n_clusters=0 區 → 整站白畫面。兩 bug 已修（見驗證記錄 log）。

## 相關文件
- 三軸分析 + 定義：`InterSubMod/docs/methodology/20260702_topology_cluster_shape_three_axis_analysis_01.md`
- 對抗驗證記錄：`InterSubMod/docs/methodology/20260702_topology_three_axis_verification_log_01.md`
- 建置地圖：`InterSubMod/docs/methodology/20260701_multisample_topology_workstation_build_map_01.md`
