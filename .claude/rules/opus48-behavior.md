---
description: Opus 4.8 model behavior characteristics and prompt strategy. Auto-loaded when editing skills, rules, or prompt configs.
globs:
  - .claude/skills/**/SKILL.md
  - .claude/skills/**/*.json
  - .claude/rules/**/*.md
  - .claude/skills/**/playbook.md
  - .claude/skills/**/prompts/*.md
---

# Opus 4.8 Model Execution Characteristics

## 模型執行特性與 Prompt 策略（Opus 4.8，2026-05 起適用；取代 opus47-behavior.md）

> **版本前提**：本專案目前 model = `claude-opus-4-8`（pin `claude-opus-4-8[1m]`，1M context），Claude Code ≥ v2.1.154。所有 skills 與 hooks 設計均假設以下特性。
> **官方來源**：code.claude.com/docs/en/model-config + platform.claude.com/docs/en/build-with-claude/effort。

### 執行行為（無法以 prompt 反轉 — 與 4.7 一致、4.8 強化）

| 特性 | 實際行為 | 對本專案的意涵 |
|------|---------|---------------|
| **Literal 指令遵循（4.8 更甚）** | 不會悄悄泛化指令、不推斷未明講需求 | 模糊指令 = 模型直接按字面做，責任在 prompt 完整度 |
| **預設少 tool calls** | 優先用 reasoning 解決，而非反覆讀檔/搜尋 | 需要廣泛掃描時，明確寫「讀遍 src/core/*.cpp」 |
| **預設少 subagent** | 單回合完成優先，除非明確 fan-out | 跨檔案/跨樣本平行任務須明寫「spawn N agents」或用 Dynamic Workflow（見 CLAUDE.md §8）|
| **回應長度動態化** | 隨任務複雜度調整 | 簡單問題不過度解釋；複雜分析仍詳盡 |
| **主動 progress update** | 長 trace 中自發回報進度 | 不需再加「每 5 步報告一次」類 scaffolding |
| **Adaptive thinking（4.8 變更）** | 預設 **OFF**；取代 4.7 的 `budget_tokens` / `display:"summarized"` 機制 | Hooks 若依賴 thinking 內容需改用其他訊號；需深推理時用 prompt 內 `ultrathink` 關鍵字 |
| **Tokenizer** | 長任務需給足 `max_tokens` headroom | xhigh/max 起手 64k max_tokens（官方建議）；compaction 閾值留 headroom |

### Effort 設定（4.8 核心 — 比任何前代 Opus 更重要）

> 官方：Opus 4.8 支援 5 級 `low / medium / high / xhigh / max`（4.6/Sonnet 4.6 無 xhigh）；**模型預設 high**；agentic/coding **官方建議從 xhigh 起跳**。本專案 settings 已 pin `effortLevel: xhigh`（project 層，防 /model 默默降級）。

| Effort | 適用 task type | 說明 |
|--------|---------------|------|
| `xhigh`（本專案預設）| B validation / C production / 長 agentic 研究 / 跨樣本 benchmark | 長時程（>30min）+ 百萬級 token budget 的 long-horizon work |
| `high` | A pilot 已有詳細 plan 的執行階段 / 一般 intelligence-sensitive | 模型原生預設 |
| `medium` / `low` | F demo / trivial lookup | 僅在 eval 證實低 effort 仍保品質時降 |
| `max` | 真正 frontier 問題 | 多數 workload 只增成本不增品質、結構化輸出易 overthink — **不設長期預設** |
| `ultracode`（session-only）| 用戶當輪明示的重 fan-out | xhigh + 自動 workflow 編排；非 settings effortLevel 值 |

**策略**：遇淺推理 → **raise effort**，而非 prompt around it。需單輪深推理 → prompt 內 `ultrathink`（不改 session effort）。

### Prompt 策略（本專案要求 — 與 4.7 一致）

1. **First-turn completeness**：首 turn 給完整規格（意圖、約束、檔案路徑含行號、驗收標準、specialist 分派）。
2. **正向範例優於否定列表**：寫「用 within-group OLS 這樣做：...」優於「不要用 pooled OLS」。
3. **Subagent / Workflow 明確觸發**：平行寫「spawn parallel-benchmark for ...」或 prompt 含 `workflow` keyword；否則預期單回合完成。
4. **Dynamic Workflow 路由**：大規模 fan-out 無 Hard Gate → workflow；含 Hard Gate / 需中間決策 → 主 agent 編排（見 CLAUDE.md §8）。

### 可移除的過度 scaffolding（4.6 遺留，4.8 下冗餘）

- 「double-check 結果合理性」「確認看起來正確」
- 「每 N 步給 interim status」「持續回報進度」
- 「預設 fan-out 給多個 subagent」（除非真需要）
- 「先渲染 PPTX 再自我檢查 layout」（內建自檢）
- 「每個決策都 FYI 告知」→ 改為**僅低影響+高信心時用一行告知**，其他依暫停判斷矩陣

### 不可放寬的硬性規則（與模型無關）

- C++ 修改的 6 步驟 PDD 協議（`/cpp-change`）
- C++ commit 前必編譯（Hard Gate hook）
- 研究方向 NO-GO 判定需用戶確認（Hard Gate）
- 刪除/搬移檔案需用戶確認（Hard Gate）
- Evidence ledger 每輪必記錄

---
> **2026-05-30 變更**：本檔取代 `opus47-behavior.md`（後者已轉為 deprecation stub，物理移除待 Hard Gate ack）。CLAUDE.md §2 + §5 已同步指向本檔。
