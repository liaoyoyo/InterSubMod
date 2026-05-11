# AGENT.md — AI Agent 使用 InterSubMod KB 的指引

**目標讀者**：Claude Code、其他 AI agent、使用 KB 回答問題的 LLM
**前置**：先讀 [README.md](README.md) 的 5 秒決策樹
**核心原則**：KB 是**導航層**，`docs/reports/` 與 source code 才是權威層

---

## ⚡ 行為協議快速入口（必讀）

執行任何任務前，依序確認（**09 → 10 → 07** 為正確依賴方向；先學政策再學方法）：

1. **[09 確認時機協議](00_governance/09_confirmation_protocol.md)**（政策）— Hard Gate / Gate / Review / FYI 4 級、執行模式、25 場景、2D 矩陣 **SoT**
2. **[10 實作前準則](00_governance/10_think_before_code.md)**（方法）— 假設陳述、Step→Verify、首 turn 規格檢核（暫停矩陣連結 09）
3. **[07 新資訊補充協議](00_governance/07_new_info_protocol.md)**（流程）— 新增結論/參數/特徵時遵循 SoT 規則

---

## 🧭 查詢決策 SOP（4 步驟）

收到使用者請求後，**嚴格按順序**執行：

### Step 1: 分類問題類型
對照 [00_governance/05_query_protocol.md](00_governance/05_query_protocol.md) 的 5 種情境：
- **A**：參數/欄位格式 → `04_parameters/` or `05_data_formats/`
- **B**：pipeline 怎麼跑 → `03_pipelines/` or `06_workflows/`
- **C**：研究結論/歷史 → `09_conclusions/`
- **D**：樣本/truth set → `02_samples/` or `08_truth_and_benchmark/`
- **E**：現在在做什麼 → `10_research_status/`

### Step 2: 進對應目錄的 `00_index.md`
不要直接猜檔名，先讀 index。每個 index 會列出：
- 該目錄所有文件的「一句結論」
- 何時該進哪份文件的決策樹

### Step 3: 讀具體文件，看 frontmatter
進入文件後**優先看 frontmatter**：
- `last_verified`：若 > 90 天，跑文件內「可直接執行命令」驗證
- `source_type`：`runtime-fact` 最新可信 / `frozen-decision` 不會變 / `historical-note` 已過期須謹慎
- `related_ids`：若發現要跨主題，從這裡跳

### Step 4: 引用來源
回答使用者時**必須**：
- 引用 KB 文件路徑：`knowledge/03_pipelines/01_paired_full.md`
- 若連結到 docs/ 權威報告：同時給 docs 路徑
- 若引用 source code：給 `file:line`

---

## 🚨 重要規則

### ✅ 該做的
1. **相信 canonical_paths**：frontmatter 的 `canonical_paths` 是單一權威路徑
2. **快照文件有時效**：`10_research_status/` 是快照，必要時追溯到 `docs/CURRENT_FOCUS.md` 最新版
3. **NEGATIVE 結論要尊重**：`09_conclusions/03_concluded_negative.md` 列出已證偽方向，**別建議重做**除非有新證據
4. **用可執行命令驗證**：每份 reference/howto 文件開頭有「可直接執行命令」，用它確認現況

### ❌ 不該做的
1. **別重寫 research_landscape**：`09_conclusions/` 只做索引，結論細節在 `docs/reports/research_landscape/`
2. **別把 KB 當 memory**：KB 是公開共享的；個人偏好/學到的教訓寫在 `MEMORY.md`
3. **別繞過 index 直接讀細節文件**：使用者請求模糊時，從 index 往下鑽
4. **別創造不存在的 file:line 引用**：不確定 source code 位置時，用 Grep/Read 確認

---

## 📋 典型對話範例

### 範例 1：使用者問「--min-mapq 預設是多少？」
```
→ Step 1: 情境 A（參數查詢）
→ Step 2: 讀 04_parameters/00_index.md → 指向 01_cli_arguments.md
→ Step 3: 讀 01_cli_arguments.md 找表格「1.3 讀取過濾參數」
→ 答案：20（取值 [0, 60]，出處 include/utils/ArgParser.hpp:67-68）
→ 引用：knowledge/04_parameters/01_cli_arguments.md + 原始碼 file:line
```

### 範例 2：使用者問「HPFineNGroups 結論如何？」
```
→ Step 1: 情境 C（研究結論）
→ Step 2: 讀 09_conclusions/00_index.md → 狀態為 positive characterization
→ Step 3: 讀 07_derived_features/01_hpfinengroups.md 看 feature 細節
→ 跨跳：docs/reports/research_landscape/09_Part_B.md（權威研究報告）
→ 答案：N≥4+NR≥80 canonical filter, TP rate 92.81%；characterization only 不可做 variant filter
→ 警告：2026-04-21 flag=on 下 N≥3 完全消失，可能 somatic HP tag artifact
```

### 範例 3：使用者問「現在最高優先在做什麼？」
```
→ Step 1: 情境 E（當前狀態）
→ Step 2: 讀 10_research_status/00_index.md → 指向 01_current_focus_snapshot.md
→ Step 3: 讀 01_current_focus_snapshot.md，警告「此為快照，最新狀態 docs/CURRENT_FOCUS.md」
→ 若快照 > 1 週，用 Read 讀 docs/CURRENT_FOCUS.md 核對
→ 答案：Phase 2 方向 A+D（Normal Methylation Reference + CN/Purity-aware correction），阻塞於 haplotag + ISM 全量重跑
```

---

## 🔧 KB 維運腳本（出錯時跑）

| 症狀 | 跑什麼 |
|------|--------|
| frontmatter 欄位缺失 | `python3 knowledge/scripts/validate_frontmatter.py knowledge/` |
| related_ids 雙向不對稱 | `python3 knowledge/scripts/check_related_ids_symmetry.py knowledge/` |
| 連結到檔案不存在 | `python3 knowledge/scripts/check_canonical_paths.py knowledge/` |
| 需要批次更新驗證日期 | `python3 knowledge/scripts/refresh_last_verified.py knowledge/ --file <path>` |

所有腳本在 [scripts/](scripts/) 目錄，執行需 python3 + pyyaml。

---

## 🚀 Opus 4.7 Prompt 策略（本專案要求）

對照 CLAUDE.md「模型執行特性」章節，AI agent 執行任務時遵循：

### First-turn completeness
收到任務**首回合**即給完整規格（意圖、約束、檔案路徑含行號、驗收標準、specialist 分派），避免多 turn 迭代補資訊。若規格缺項 ≥2 且涉及高影響決策 → 必須回問。

### Literal 指令特性（避免按字面做錯事）
Opus 4.7 **不會**推斷未明講需求、不會悄悄泛化指令。防禦做法：
- 開工前先列「我理解的任務規格」（意圖 / 約束 / 路徑 / 驗收 / scope 邊界）
- 模糊指令時列假設後繼續（低影響）或暫停（中高影響）

### Subagent 明確觸發
**預設不 spawn subagent**；需要時明確寫「spawn parallel-benchmark for HCC1395, COLO829, H2009」（樣本名與 `02_samples/` 對齊）。
- 跨樣本並行（同時 ≥3 樣本在背景平行）→ `parallel-benchmark`（時間上並行；序列跑不觸發）
- >3 查詢探索 → `Explore` agent
- 結論可信度驗證 → `general-purpose` 紅隊
- PR 審查 → `pr-review-toolkit:*`

### 避免過度 scaffolding（4.6 遺留）
以下語句在 4.7 下多為冗餘，審查 prompt 時可刪：
- 「double-check 結果合理性」「確認看起來正確」
- 「每 N 步給 interim status」
- 「預設 fan-out 給多個 subagent」

### Effort 設定
- `xhigh`（預設）— 大多數 agentic 研究任務
- `high` — 已有詳細 plan 的執行階段
- `max` — 僅真正困難的問題

---

## 🔬 實驗室 MCP 知識庫對照

本 KB（`/big7_disk/.../knowledge/`）聚焦 **ISM 專案**；另有跨專案實驗室 KB：

**路徑**：`/big8_disk/liaoyoyo2001/knowledge/`（MCP server 已接入）

**5 個 MCP 工具**：
| 工具 | 用途 |
|------|------|
| `knowledge_search` | 全文搜尋（OR+AND bonus、`-排除詞`、中英混）|
| `knowledge_list` | 列所有文件（可按 tag/status 篩）|
| `knowledge_get_doc` | 取文件全文或章節 |
| `knowledge_resolve_path` | 路徑/別名解析 |
| `knowledge_backlinks` | 查誰引用該文件 |

**用途邊界**：
- 跨專案事實（樣本資料路徑、工具版本）→ **MCP KB**
- ISM 專案規格（pipeline、參數、結論）→ **本 KB**
- 兩者 frontmatter 格式一致，便於 future merge

**AI agent 原則**：
- 查 ISM 相關先用本 KB（`knowledge/`）
- 查 `/big8_disk/data/` 路徑或跨工具 workflow 再跳 MCP
- 不要憑記憶回答可查證事實 — 查閱後引用來源

---

## 🔗 外部權威文件（KB 不複製，只連結）

| 主題 | 權威路徑（專案根相對） |
|------|------------------------|
| 當前進度 | `docs/CURRENT_FOCUS.md` |
| 研究景觀總索引 | `docs/reports/research_landscape/00_INDEX.md` |
| 實驗歷史索引 | `docs/experiments/INDEX.md` |
| 研究構想總索引 | `docs/concepts/2026/04/20260409_研究構想總索引_01.md` |
| Truth set / F1 協議 | `docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md` |
| significance_summary 字典 | `docs/data_specs/20260411_significance_summary欄位字典_01.md` |
| master_dataset 字典 | `docs/data_specs/20260411_master_dataset欄位字典_01.md` |
| 專案規則 | `.claude/CLAUDE.md` |
| 假說佇列 | `research/autoresearch/hypothesis_queue.json` |
| Evidence ledger | `research/autoresearch/evidence_ledger.jsonl` |
| AI memory index | `/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md` |

---

## ⚠️ 常見陷阱（避免誤答）

1. **誤把 MEMORY.md 當 KB**：MEMORY 是個人偏好與教訓，KB 是公開規格；兩者不重疊
2. **誤以為 TO mode 等於 paired mode 的子集**：TO 有 self-phasing bias，結果不可直接比較
3. **誤以為 HPFineNGroups 可做 variant filter**：結論是 characterization only（會產生負 F1 增益）
4. **誤用舊 canonical run**：2026-04-04 前的 TO run 有 VCF 來源錯誤（見 MEMORY: project_vcf_source_error_correction）
5. **誤以為 docs/CURRENT_FOCUS.md 與 10_research_status/ 完全同步**：KB 快照更新較慢，時效敏感查詢以 docs/ 為準

---

**本文件版本**：v0.1
**最後更新**：2026-04-22
