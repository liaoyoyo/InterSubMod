# CLAUDE.md - Claude Code 專案指南

## 繼續研究前的必讀清單（每次對話開始時強制執行）

**每次開始研究/分析任務前，必須依序閱讀以下文件，不得省略：**

1. **`docs/CURRENT_FOCUS.md`** — 當前進行中的事項、阻塞點與風險
2. **`docs/experiments/INDEX.md`** — 過去所有研究方向的成功/失敗結論與建議後續
3. **`docs/README.md`** — 如需了解文件導航與查閱路徑
4. **`docs/concepts/2026/04/20260409_研究構想總索引_01.md`** — 研究大圖景、發展樹、理論基礎、論文規劃
5. **`docs/reports/research_landscape/00_INDEX.md`** — 完整研究推論鏈、14 結論穩定性、8 證據鏈（需深度理解時）

**目的**：避免重複已失敗的方向、對齊當前最優先目標、了解哪些結論已驗證、哪些尚未解決。

**觸發條件**：開始任何研究分析、實驗設計、程式改進、或延續前次工作時，此步驟為必要前置。

### 重點資訊速查

| 我想知道... | 去哪裡找 |
|-------------|---------|
| TO FP 問題全貌 | `docs/reports/research_landscape/01_TO_FP問題全貌.md` |
| Self-phasing 根因與影響 | `docs/reports/research_landscape/02_Self_Phasing根因.md` |
| 哪些 ISM 特徵可信 | `docs/reports/research_landscape/03_ISM分析價值界定.md` |
| 暫停判定與修正後預期 | `docs/reports/research_landscape/04_暫停判定與重評估.md` |
| 完整證據鏈推論 | `docs/reports/research_landscape/05_證據鏈總覽.md` |
| 結論穩定性評分 | `docs/reports/research_landscape/06_結論穩定性審查.md` |

---

## Context 壓縮保留指令

壓縮（compact）時**必須保留**：架構決策及理由、未解決 blockers、涉及檔案路徑清單、用戶限制條件、假說 ID 與驗證層級、未完成步驟清單、所有 CLAUDE.md 規則。

**可安全壓縮**：中間計算過程、已通過的測試完整輸出、冗長程式碼輸出、重複的工具呼叫結果。

---

## 模型執行特性與 Prompt 策略（Opus 4.7，2026-04-17 起適用）

**Opus 4.7 的關鍵行為差異直接影響本專案流程**。所有 skills 與 hooks 設計均假設以下特性：

### 執行行為（無法以 prompt 反轉）

| 特性 | 實際行為 | 對本專案的意涵 |
|------|---------|---------------|
| **Literal 指令遵循** | 不會悄悄泛化指令、不推斷未明講需求 | 模糊指令 = 模型直接按字面做，責任在 prompt 完整度 |
| **預設少 tool calls** | 優先用 reasoning 解決，而非反覆讀檔/搜尋 | 需要廣泛掃描時，明確寫「讀遍 src/core/*.cpp」 |
| **預設少 subagent** | 單回合完成優先，除非明確 fan-out | 跨檔案/跨樣本平行任務須明寫「spawn N agents」 |
| **回應長度動態化** | 隨任務複雜度調整（不再固定冗長） | 簡單問題不會得到過度解釋；複雜分析仍詳盡 |
| **主動 progress update** | 長 trace 中會自發回報進度 | 不需再加「每 5 步報告一次」類 scaffolding |
| **Tokenizer 改版** | token 用量為 4.6 的 1.0–1.35× | 長任務 `max_tokens` 與 compaction 閾值需給足 headroom |
| **Thinking 預設不輸出** | 需設 `display: "summarized"` 才回傳推理內容 | Hooks 若依賴 thinking 內容需調整 |

### Prompt 策略（本專案要求）

1. **First-turn completeness**：首 turn 給完整規格（意圖、約束、檔案路徑含行號、驗收標準、specialist 分派），避免多 turn 迭代補資訊
2. **正向範例優於否定列表**：寫「用 within-group OLS 這樣做：...」優於「不要用 pooled OLS」
3. **Subagent 明確觸發**：需要平行時寫「spawn parallel-benchmark for HCC1395_5kHz, COLO829, H2009」；否則預期模型單回合完成
4. **Effort 建議**：
   - `xhigh`（預設）— 大多數 agentic 研究任務
   - `high` — 已有詳細 plan 的執行階段
   - `max` — 僅真正困難的問題（會 overthink）
5. **Task budgets（beta，可選）**：長迴圈（如 HPFineNGroups 全樣本）可設 `task_budget` 讓模型自我節流

### 可移除的過度 scaffolding（4.6 遺留）

以下語句在 4.7 下多為冗餘，審查既有文件/prompt 時可刪除：
- 「double-check 結果合理性」「確認看起來正確」
- 「每 N 步給 interim status」「持續回報進度」
- 「預設 fan-out 給多個 subagent」（除非真需要）
- 「先渲染 PPTX 再自我檢查 layout」（4.7 slide/chart 能力已內建自檢）
- 「每個決策都 FYI 告知」→ 改為**僅低影響+高信心時用一行告知**，其他情境依暫停判斷矩陣

### 不可放寬的硬性規則（與模型無關）

以下為本專案結構性要求，不因模型升級而鬆動：
- C++ 修改的 6 步驟 PDD 協議（`/cpp-change`）
- C++ commit 前必編譯（Hard Gate hook）
- 研究方向 NO-GO 判定需用戶確認（Hard Gate）
- 刪除/搬移檔案需用戶確認（Hard Gate）
- Evidence ledger 每輪必記錄

---

## 實作前行為準則（Think Before Code）

### 假設陳述規則

**開始任何實作、分析、或研究決策前，必須先列出關鍵假設。**

> **⚠ Opus 4.7 literal 特性**：模型不會推斷未明講需求、不會悄悄泛化指令。模糊輸入將被直接按字面執行。本節規則比 4.6 時代更為關鍵 — 假設陳述是避免「按字面做錯事」的唯一屏障。

- 明確寫出你正在做的假設。不確定的標記為「⚠ 待確認」，不要默默選擇。
- 存在多種解讀時，列出所有合理解讀並說明你傾向哪一種及理由 — 不要自行挑選後沉默執行。
- 如果存在更簡單的替代方案，主動提出。必要時推回過度複雜的需求。
- **遇到模糊不清的指令 → 查「暫停判斷矩陣」**（見下節）：
  - 不可逆操作 → 🔴 立即暫停（與影響/信心無關）
  - 高影響 OR 低信心 → 🔴 立即暫停報告
  - 低影響 + 高信心 + 可逆 → 🟢 **一行告知後繼續**，不打斷工作流
  - 其餘中間地帶 → 🟡/🟠 列假設後繼續，節點給進度
- **首 turn 規格完整度檢核**：接到任務後，先列出「我理解的任務規格」（意圖／約束／檔案路徑／驗收標準／scope 邊界）。**規格缺項 ≥2 且涉及高影響決策**時必須回問；否則列出預設選擇與理由後繼續。

### 暫停判斷矩陣（2D + 可逆性 override）

| 影響度 ＼ 信心度 | 高信心（有先例/單一解讀） | 中信心（2-3 種解讀，偏向明確） | 低信心（無先例/多重解讀/規格缺項） |
|----------------|------------------------|----------------------------|-----------------------------|
| **低**（可逆/局部/<10min） | 🟢 一行告知後繼續 | 🟡 列假設後繼續 | 🟠 節點暫停報告 |
| **中**（多檔案/需重跑/10min-1h） | 🟡 列假設後繼續 | 🟠 節點暫停報告 | 🔴 立即暫停 |
| **高**（跨模組/影響結論/>1h） | 🟠 節點暫停報告 | 🔴 立即暫停 | 🔴 立即暫停 |

**可逆性 override**：任何不可逆操作（刪檔／C++ commit／研究方向 NO-GO／覆寫 evidence_ledger／啟動 >10min 計算）永遠 🔴 立即暫停，不論矩陣落點。對應現有 `/confirmation-protocol` 的 Hard Gate。

**一行告知格式**：`[決策]（影響: 低/中/高, 信心: 低/中/高, 理由: 一句）`
**範例**：`以 within-group OLS 殘差化（影響: 低, 信心: 高, 理由: O12 已驗證 pooled OLS 的 collider bias）`

**完整維度定義、場景對照、決策流程** → 詳見 `/confirmation-protocol` skill「自主推進 vs 暫停的動態判斷」章節。

### 高影響場景清單（預設落在矩陣的中/高影響區）

下列場景**典型落在高影響區**，需有明確高信心理由才能自主推進；否則預設暫停：

| 場景 | 典型影響評估 | 為什麼危險 | 正確做法 |
|------|------------|-----------|---------|
| 研究重點排序 | 高（>1h 投入） | 選錯優先序浪費數天工作 | 列候選方向 + 預期收益/風險，等用戶決定 |
| 假說選擇 | 高（影響結論） | 隱含假設可能導致結論偏誤 | 寫出假說前提條件和可能 confound |
| 統計方法選擇 | 中-高（影響結論） | 不同方法可能得相反結論 | 說明為何選此方法、替代方案、預期差異 |
| 「改進」/「優化」模糊指令 | 高（方向未定） | 改進方向無數，選錯浪費時間 | 要求用戶定義成功標準，或提 2-3 方向供選 |
| 多檔案/多步驟重構 | 中-高（影響範圍廣） | 影響範圍不明確 | 列受影響檔案 + 預期改動，確認後再動手 |

### 目標驅動驗證格式（Step → Verify）

**多步驟任務必須轉化為可驗證的執行計劃。** 每一步都要有明確的驗證標準：

```
1. [步驟描述] → 驗證: [具體可觀察的檢查方式]
2. [步驟描述] → 驗證: [具體可觀察的檢查方式]
3. [步驟描述] → 驗證: [具體可觀察的檢查方式]
```

**範例**：

```
任務：加入 Normal BAM read 過濾邏輯

1. 讀取現有 ReadParser 了解 BAM 讀取流程
   → 驗證: 能指出 normal read 進入的函數位置和資料結構

2. 新增 normal BAM 路徑參數與讀取邏輯
   → 驗證: make -j$(nproc) 編譯通過，無 warning

3. 實作 normal read 過濾條件
   → 驗證: 單元測試覆蓋 normal-only / tumor-only / mixed 三種情境

4. 全量測試確認無回歸
   → 驗證: run_batch_vcf_analysis.sh 結果與修改前 F1 差異 < 0.01
```

**弱驗證標準（禁止）**：「讓它能動」、「確認結果合理」、「看起來正確」、「double-check 一下」
**強驗證標準（要求）**：具體數值、特定檔案輸出（含路徑）、可重現的測試命令、命令退出碼、預期輸出片段

> **Opus 4.7 備註**：模型自我驗證能力已提升，不需再加「確認看起來正確」類軟性 scaffolding。驗證標準必須是**外部可觀察**的（檔案存在、數值範圍、命令成功），不是「模型主觀判斷合理」。

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

## 建置與程式碼規範

建置命令、C++ 改動檢查清單、新增模組流程 → 自動載入於 `.claude/rules/cpp-build.md`（觸發：編輯 src/include/tests/CMakeLists.txt 時）。

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

---

## 實驗室知識庫 (Knowledge Base)

**路徑**：`/big8_disk/liaoyoyo2001/knowledge/`（MCP server 已接入）。主題表與查閱深度詳見 `AGENTS.md`「實驗室知識庫」區段。
- **不要憑記憶回答可查證的事實** — 查閱後引用來源
- **涉及 OLS/VCF 來源/特徵設計時** — 先查 `/known-pitfalls` skill

---

## 常用工作流程

VCF 分析、批次測試、快速驗證命令 → 自動載入於 `.claude/rules/workflow-commands.md`（觸發：編輯 scripts/ 時）。

---

## 輸出檔案結構

Region 目錄結構說明 → 自動載入於 `.claude/rules/output-structure.md`（觸發：操作 output/ 時）。

### Output 頂層結構速覽

| 目錄 | 分類 | 說明 |
|------|------|------|
| `canonical/` | Canonical | 7 樣本 × 3 模式 ISM baseline（19 runs） |
| `synthesis/` | Synthesis | 觀察、研究 rounds、批次、concluded、manifests |
| `bip8_output_archive/` | Archive | bip8 歷史（39 dirs） |
| `big8_output_archive/` | Archive | big8 歷史（5 dirs） |
| `multilayer_hp_benchmark/` | Standalone | HP 比較 benchmark |

**完整索引**: `output/OUTPUT_INDEX.md`
**信任度規範**: `docs/data_specs/20260414_output資料信任度與生命週期_01.md`

---

## 開發重點

1. **Phase 2 Normal Methylation Reference (方向 A+D)**: Normal BAM 整合、Sample ASM、LOH BED 標註、跨區域亞克隆分析 — 當前最高優先
2. **甲基化解析精確度**: MM/ML 標籤解碼與 CIGAR 座標校正的正確性
3. **距離計算效能**: OpenMP 平行化與稀疏矩陣最佳化
4. **統計假設檢驗**: PERMANOVA p-value 與 Cramér's V 效應量計算
5. **視覺化品質**: 熱圖標註 (HP tag / Allele) 與分群樹狀圖

> **定位說明（2026-04 確認）**：ISM 的核心價值在 **read-level epigenetic characterization**，而非 variant filter。Phase 1A F1 優化已鎖定（paired-pure delta F1=+0.0112），TO 模式甲基化增益為負。當前重心是 Phase 2 生物學特徵化方向。

---

## 文件資源

- `README_PROJECT_SUMMARY.md`: 專案完整技術總結
- `QUICKSTART.md`: 快速入門指南
- `docs/`: 完整技術文件 (API、架構、開發、報告)

---

## 文檔管理規範

檔案命名、元數據、圖片、目錄結構等規範詳見 `/doc-standards` skill。建立新 .md 檔案時自動觸發。

---

## Hooks 配置

專案使用 Claude Code Hooks 自動化檢查流程，配置於 `.claude/settings.local.json`。

### 現有 Hooks

| Hook 類型 | 觸發條件 | 動作 | 模式 |
|-----------|----------|------|------|
| UserPromptSubmit | 用戶提問 | 知識庫相關主題提醒 | 提醒 |
| UserPromptSubmit | 用戶提問 | **已關閉研究方向警告** | 警告 |
| PreToolUse | `git commit` | **未編譯 C++ 時阻擋 commit** | **阻擋** (`exit 2`) |
| PreToolUse | 危險命令 | rm -rf / reset --hard 警告 | 提醒 |
| PostToolUse | 編輯 `.cpp`/`.hpp`/`.h` | **建立待編譯標記** | 標記 |
| PostToolUse | `make` 成功 | **清除待編譯標記** | 清除 |
| PostToolUse | Edit/Write 任意檔案 | **觸發路由提示**（依檔案類型建議對應流程） | 提醒 |
| PostToolUse | `git commit` | 提醒確認測試和文檔 | 提醒 |
| SubagentStop | 子代理完成 | 提醒檢查結果 | 提醒 |
| Stop | 會話結束 | 提醒撰寫執行報告 | 提醒 |

### 注意事項

- Hooks 會在對應事件觸發時自動執行
- 務必遵循 Hooks 提示完成必要步驟
- 如需新增或修改 Hooks，編輯 `.claude/settings.local.json` 的 `hooks` 區段

---

## 確認時機協議

兩種模式：**互動模式**(預設，觸發詞：`互動模式`) / **全自動模式**(觸發詞：`全自動`、`auto`)。
暫停級別：**Hard Gate**(必停) > **Gate**(互動停/自動過) > **Review**(互動停/自動過) > **FYI**(顯示/靜默)。
Hard Gate 場景：刪除檔案、C++ commit、研究方向 NO-GO 判定。
**Opus 4.7 subagent 觸發規則**：預設不 spawn subagent；跨樣本平行、>3 查詢探索、PR 審查等需明確觸發語。
完整規則、場景分類表、subagent 觸發條件詳見 `/confirmation-protocol` skill。
