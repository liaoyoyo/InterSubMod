# CLAUDE.md - Claude Code 專案指南

## 研究上下文載入策略（觸發式 on-demand）

每次對話開始時，CLAUDE.md 自身 + `docs/CURRENT_FOCUS.md` 為**唯一 always-loaded** 上下文（~3k tokens）。

**研究分析啟動時**：呼叫 `/research-context-loader` skill，依問題深度分 3 tier 載入：
- Tier 1（always）：CURRENT_FOCUS.md
- Tier 2（light, ~3-5k tokens）：experiments INDEX + landscape 00_INDEX + concepts 索引
- Tier 3（deep, ~5-15k tokens）：依問題主題載入特定 landscape 檔案

完整 landscape 速查表移到 skill 內 — 詳見 `.claude/skills/research-context-loader/SKILL.md`。

**何時不需要 loader**：純 code edits（make/test/commit）、單檔 doc 寫作、簡單問答 — 預設不載入 landscape。

---

## Agent 上下文控制面（2026-04-27）

本專案同時有 repo 規範、Claude 執行規則、研究壓縮上下文、memory 與 AutoResearch queue。為避免 agent 啟動時混用過時來源，分工如下：

| 入口 | 權威範圍 | 使用時機 |
|------|----------|----------|
| `AGENTS.md` | repo 內硬規則、Knowledge Base 查閱義務、output/runbook/meeting 分流 | 判斷可讀寫範圍、新檔落點、不可刪檔規則 |
| `.claude/CLAUDE.md` | Claude Code 行為、確認矩陣、hooks、壓縮保留、C++ gate | 約束 agent 怎麼做事 |
| `docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` | 研究壓縮上下文、重要數據、任務順序、待決策矩陣 | 每次研究/分析的快速定向 |
| `docs/CURRENT_FOCUS.md` | live 主軸、阻塞、最新撤回或主軸切換 | 每次對話開始確認現況 |
| `research/autoresearch/research_direction.md` | AutoResearch 候選 queue | 只作候選，不作自動執行觸發 |

啟動研究任務時先回答 5 問：是否涉及 Thread D、是否碰到 Thread B 撤回範圍、資料是否 KDE-corrected、是否需要 VCF caller AF、是否觸及長計算/C++/搬移/NO-GO gate。

---

## Context 壓縮保留指令

壓縮（compact）時**必須保留**：架構決策及理由、未解決 blockers、涉及檔案路徑清單、用戶限制條件、假說 ID 與驗證層級、未完成步驟清單、所有 CLAUDE.md 規則。

若要跨 session 交接，先以 `docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` 作為壓縮上下文骨架，再補充本輪新增的 artifacts、verdict 與未完成 gate。

**可安全壓縮**：中間計算過程、已通過的測試完整輸出、冗長程式碼輸出、重複的工具呼叫結果。

---

## 模型執行特性與 Prompt 策略

詳見 `.claude/rules/opus47-behavior.md`（編輯 SKILL.md / *.json 時自動觸發載入）。

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

**可逆性 override**：任何不可逆操作（刪檔／C++ commit／研究方向 NO-GO／覆寫 evidence_ledger）永遠 🔴 立即暫停，不論矩陣落點。對應現有 `/confirmation-protocol` 的 Hard Gate。

**長時間計算（>10 min）例外**：用戶在**當輪對話明示**啟動長計算（例：「跑全量 benchmark」「平行 7 樣本」）時，**不再二次確認**，直接執行並一行告知。若用戶未明示而模型自行判斷需長計算 → 仍 🔴 暫停。

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

### 外部工具 / 論文 claim 查詢順序（2026-05-13 新增）★

> 觸發：討論外部工具 (longphase / clairs-to / claude code) F1/AUC/metric claim、論文 claim 對照、業界標準口徑（bcftools isec / hap.py / som.py）、tool README/論文行為敘述時 — **不可從本專案 report 推論外部行為**

**強制查詢順序**：
1. **先查 KB**：`mcp__knowledge__knowledge_search "<tool> F1"` 或 `Read /big8_disk/liaoyoyo2001/Knowledge/05_tools/<tool>.md`
2. **引用 paper §N 或 KB 段落**（明示出處）
3. **才下結論** — 並對照本專案 report 是否口徑一致
4. 若 KB + 論文與本專案 report **看似衝突** → 通常是**不同口徑**（如 caller-level F1 vs post-filter F1），KB + 論文優先；除非本專案 report 明寫「我們的修補 vs 論文此處」

**反例（不該做）**：
- ❌ 「§8.6.2 寫 F1 三版相同 → 結論 F1 不變」（沒查論文 §4.3 的不同 metric 口徑）
- ❌ 「我記得 longphase 不改 FILTER」（憑記憶未引用）

**正例（該做）**：
- ✅ 「先 Read `Knowledge/05_tools/longphase-to.md` 確認論文 §4.3 F1 口徑為 V_H/V_L post-filter」→ 引用後再對照本專案 §8.6.2 caller-level F1

**對應**：`/known-pitfalls` P-14 + memory `feedback_outside_claim_must_query_kb.md` + UserPromptSubmit hook `[Knowledge Base]` 提示

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
