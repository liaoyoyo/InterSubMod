# CLAUDE.md - Claude Code 專案指南

## 繼續研究前的必讀清單（每次對話開始時強制執行）

**每次開始研究/分析任務前，必須依序閱讀以下文件，不得省略：**

1. **`docs/CURRENT_FOCUS.md`** — 當前進行中的事項、阻塞點與風險
2. **`docs/experiments/INDEX.md`** — 過去所有研究方向的成功/失敗結論與建議後續
3. **`docs/README.md`** — 如需了解文件導航與查閱路徑

**目的**：避免重複已失敗的方向、對齊當前最優先目標、了解哪些結論已驗證、哪些尚未解決。

**觸發條件**：開始任何研究分析、實驗設計、程式改進、或延續前次工作時，此步驟為必要前置。

---

## 專案概述

**InterSubMod (Inter-Subclonal Methylation Analysis)** 是一個生物資訊工具，專門用於透過長讀取 (Long-read) 測序數據偵測腫瘤樣本中的亞克隆結構 (Subclonal Structure)。本專案整合甲基化模式 (Methylation Patterns)、體細胞變異 (Somatic SNVs) 以及單倍體型 (Haplotypes) 來分析表觀遺傳異質性 (Epigenetic Heterogeneity)。

### 核心技術特點

- **高效能 C++17 核心**: 結合 OpenMP 平行運算，單 Region 平均處理時間 < 300ms
- **精確甲基化解析**: 支援 BAM 格式 MM/ML 標籤，精確定位 CpG 位點甲基化狀態
- **多元距離度量**: L1 / L2 / NHD / Bernoulli / Jaccard 等多種距離算法
- **統計顯著性分析**: PERMANOVA、卡方檢驗、Cramér's V 效應量計算
- **自動化視覺化**: Python 工具生成距離熱圖 (Distance Heatmap) 與分群熱圖 (Cluster Heatmap)

---

## 建置命令

```bash
# 配置並建置 (從專案根目錄)
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Debug 模式建置
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# 執行所有測試
cd build && ctest --output-on-failure

# 執行特定測試
./build/tests/test_<name>
```

## 程式碼規範

- **C++17** 標準
- 遵循 `.clang-format` (Google 風格、120 字元行寬、4 空格縮排)
- 提交前格式化：`clang-format -i <file>`
- **預設回應與註解語言**：繁體中文
- **程式碼註解語言**：英文

---

## 關鍵依賴

| 依賴 | 用途 |
| :--- | :--- |
| HTSlib | BAM 檔案處理 |
| OpenMP | 平行運算 |
| Eigen3 | 線性代數 |
| GoogleTest | 單元測試 |
| jemalloc 5.3.0 | 記憶體分配 |
| Python3 + Matplotlib/Seaborn/Scipy/Pandas | 視覺化 |


```bash
# 執行所有測試
cd build && ctest
```

---

## 實驗室知識庫 (Knowledge Base)

**路徑**：`/big8_disk/liaoyoyo2001/Knowledge/`

當對話涉及以下主題時，**必須**先查閱知識庫對應文件確認細節，再進行回答或操作：

| 主題 | 查閱路徑 | 觸發關鍵字 |
|------|---------|-----------|
| 資料總覽與路徑 | `01_data_overview/` | 資料位置、目錄結構、儲存空間 |
| 癌症樣本資訊 | `02_samples/` | HCC1395, COLO829, H1437, H2009, HG002, purity, subsample |
| 檔案格式規格 | `03_file_formats/` | VCF, BAM, MM/ML, FILTER, phased VCF, modcall, HP tag |
| 資料庫與參考集 | `04_databases/` | PON, gnomAD, dbSNP, CoLoRSdb, SEQC2, truth set, reference genome |
| 工具使用與參數 | `05_tools/` | LongPhase, ClairS, ClairS-TO, DeepSomatic, InterSubMod |
| 分析流程 | `06_workflows/` | somatic calling, phasing, haplotagging, methylation analysis, benchmark |
| 腳本操作說明 | `07_scripts/` | auto_run.sh, benchmark script, 自動化腳本 |
| 論文與參考資料 | `08_references/` + `paper/` | paper, 論文, server paths |

### 查閱深度指引

| 情境 | 查閱深度 | 動作 |
|------|---------|------|
| 快速確認（路徑、名稱） | 淺層 | 讀 `README.md` 速查表 |
| 格式或參數細節 | 中層 | 讀對應子目錄的特定文件 |
| 完整流程或工具操作 | 深層 | 讀 workflow + tool 文件，交叉驗證 |
| 工具原始碼邏輯 | 最深層 | 讀 `codebase/` 目錄下的原始碼 |

### 查閱原則

- **不要憑記憶回答可以查證的事實**：檔案路徑、工具參數、VCF 欄位定義等務必查閱確認
- **引用來源**：回答時標註「根據 Knowledge/03_file_formats/vcf_clairs_to.md」
- **發現過時資訊時主動提醒使用者**

---

## 常用工作流程

### 1. 完整 VCF 分析 (預設執行命令)

```bash
./scripts/run_vcf_all_snv.sh --mode all-with-w5000
```

### 2. 批次分析 (TP/FP 比較)

```bash
./scripts/run_batch_vcf_analysis.sh
```

### 3. 單點快速驗證

```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification
```

### 4. 新增分析模組

1. 在 `include/core/` 建立標頭檔
2. 在 `src/core/` 建立實作檔
3. 更新 `CMakeLists.txt` (如需要)
4. 在 `tests/` 撰寫單元測試

---

## 輸出檔案結構

每個 Region 目錄包含：

- `metadata.txt`: 區域基本資訊與統計數據
- `reads.tsv`: Read 詳細資訊
- `methylation.csv`: 甲基化矩陣 (Read × CpG)
- `distance_matrix_[METRIC].csv`: Read-Read 距離矩陣
- `significance_summary.csv`: 顯著性分析結果彙總
- `*.png`: 視覺化圖表

---

## 開發重點

1. **甲基化解析精確度**: MM/ML 標籤解碼與 CIGAR 座標校正的正確性
2. **距離計算效能**: OpenMP 平行化與稀疏矩陣最佳化
3. **統計假設檢驗**: PERMANOVA p-value 與 Cramér's V 效應量計算
4. **視覺化品質**: 熱圖標註 (HP tag / Allele) 與分群樹狀圖
5. **批次分析流程**: TP/FP 位點比較與 F1 分數優化

---

## 文件資源

- `README_PROJECT_SUMMARY.md`: 專案完整技術總結
- `QUICKSTART.md`: 快速入門指南
- `docs/`: 完整技術文件 (API、架構、開發、報告)

---

## 文檔管理規範

### 文檔目錄結構

```
docs/
├── architecture/        # 專案主軸架構說明
├── concepts/            # 構想紀錄（扁平結構）
├── plans/               # 計劃文件（YYYY/MM 分層）
├── solutions/           # 問題解決報告（任務目標/YYYY/MM 分層）
├── experiments/         # 實驗紀錄
│   ├── outputs/         # 測試輸出（任務目標/YYYY/MM 分層）
│   └── parameters/      # 參數研究
├── ai_sessions/         # AI 對話紀錄（YYYY/MM 分層）
├── data_specs/          # 數據規格
├── references/          # 參考資料
├── archive/             # 歸檔
│   └── deep/            # 深度歸檔（需查詢歷史）
└── trash/               # 暫存待刪除
```

### 檔案命名格式

```
{YYYYMMDD}_{中文說明目標}_{流水號}.md
```

範例：`20260111_文檔庫重整計劃_01.md`

### 檔案元數據

每個 Markdown 檔案開頭需包含：

```markdown
<!--
建立時間: YYYY-MM-DD HH:MM
目標: [本檔案的目標或用途]
處理範圍: [涵蓋的工作範圍]
關聯檔案:
  - [相關檔案路徑 1]
  - [相關檔案路徑 2]
-->
```

### AI 對話紀錄撰寫

每次 AI 對話完成重要任務後，撰寫執行報告：

1. **報告位置**：`docs/provenance/ai_sessions/YYYY/MM/`
2. **檔案格式**：`{YYYYMMDD}_{對話主題}_{流水號}.md`
3. **必要內容**：
   - 對話目標
   - 關鍵決策
   - 修改的檔案清單
   - 後續行動建議

---

## 程式碼改動檢查清單

修改程式碼後，必須完成以下步驟：

```bash
# 1. 編譯程式碼
cd build && make -j$(nproc)

# 2. 執行測試
./scripts/run_batch_vcf_analysis.sh

# 3. 確認測試結果合理

# 4. 更新 Docker 配置（如需要）
```

---

## Hooks 配置

專案使用 Claude Code Hooks 自動化檢查流程，配置於 `.claude/settings.local.json`。

### 現有 Hooks

| Hook 類型 | 觸發條件 | 動作 |
|-----------|----------|------|
| PostToolUse | 編輯 `.cpp`/`.hpp`/`.h` 檔案 | 提醒編譯和測試 |
| PostToolUse | 執行 `git commit` | 提醒確認測試和文檔 |
| Stop | 會話結束 | 提醒撰寫執行報告 |

### 注意事項

- Hooks 會在對應事件觸發時自動執行
- 務必遵循 Hooks 提示完成必要步驟
- 如需新增或修改 Hooks，編輯 `.claude/settings.local.json` 的 `hooks` 區段
