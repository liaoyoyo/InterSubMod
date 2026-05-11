---
description: Opus 4.7 model behavior characteristics and prompt strategy. Auto-loaded when editing skills, rules, or prompt configs.
globs:
  - .claude/skills/**/SKILL.md
  - .claude/skills/**/*.json
  - .claude/rules/**/*.md
  - .claude/skills/**/playbook.md
  - .claude/skills/**/prompts/*.md
---

# Opus 4.7 Model Execution Characteristics

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
