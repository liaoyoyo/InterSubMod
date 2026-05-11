---
id: ism-kb-00-governance-index
name: "Governance 目錄索引"
description: "本 KB 自身的規範文件索引：frontmatter schema、命名、交叉引用、freshness、查詢 SOP、更新流程。寫入新文件前必讀。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "governance directory structure"
related_ids:
  - ism-kb-00-governance-new-info-protocol
  - ism-kb-00-governance-hooks-and-automation
  - ism-kb-00-governance-confirmation-protocol
  - ism-kb-00-governance-think-before-code
tags: [governance, index, standards, meta]
canonical_paths: [00_governance/00_index.md]
alias_paths: []
---

# Governance 目錄索引

- 一句結論：寫入或修改 KB 任何文件前，先確認本目錄的規範
- 適用對象：KB 文件作者、維運者、新增 agent 時
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls -la /big7_disk/liaoyoyo2001/InterSubMod/knowledge/00_governance/
  ```

---

## 文件列表

| 檔案 | 主題 | 何時讀 |
|------|------|--------|
| [01_frontmatter_schema.md](01_frontmatter_schema.md) | YAML frontmatter 各欄位定義與取值 | 寫新文件或修改 metadata 前 |
| [02_naming_conventions.md](02_naming_conventions.md) | 檔名、id、canonical_path 命名規則 | 新增文件、重新命名時 |
| [03_cross_reference_rules.md](03_cross_reference_rules.md) | related_ids 雙向、連結 docs/research 標準 | 建立跨文件連結時 |
| [04_freshness_policy.md](04_freshness_policy.md) | last_verified 管理、失效處理 | 維護老化文件前 |
| [05_query_protocol.md](05_query_protocol.md) | ★ 5 種典型查詢情境 SOP | AI agent 首次使用 KB 或模糊查詢時 |
| [06_update_workflow.md](06_update_workflow.md) | 新增/修改/淘汰文件流程 | 執行任何文件變動前 |
| [07_new_info_protocol.md](07_new_info_protocol.md) | ★ 新資訊補充與驗證協議（5 類別 × 3 層驗證 × SoT 規則） | 新增結論/feature/數值時 |
| [08_hooks_and_automation.md](08_hooks_and_automation.md) | 10+ Hooks 配置與自動化（含 C++ commit Hard Gate） | 了解被擋 commit 原因、新增 hook 前 |
| [09_confirmation_protocol.md](09_confirmation_protocol.md) | ★ 確認時機協議速查（Hard Gate / Gate / Review / FYI + 模式切換） | AI 決策時、模式切換時 |
| [10_think_before_code.md](10_think_before_code.md) | ★ 實作前準則（假設陳述 + 暫停矩陣 + 目標驅動驗證） | 接到任何實作任務前 |

---

## 核心規範速查

- **每份文件必有 frontmatter** — 見 `01_frontmatter_schema.md`
- **每份文件必有「一句結論 + 適用對象 + 可直接執行命令」** — 文件開頭三欄
- **related_ids 必須雙向** — 見 `03_cross_reference_rules.md`
- **last_verified 每次實質更新必改** — 見 `04_freshness_policy.md`
- **命名格式 `NN_keyword_keyword.md`** — 見 `02_naming_conventions.md`

---

## 與既有專案規範的關係

| 規範 | 此 KB 處理 | 專案既有規範 | 衝突時依哪個 |
|------|-----------|--------------|--------------|
| 文件命名 | `NN_keyword.md`（下劃線 + 小寫） | `/doc-standards` skill（繁中檔名） | **KB 內用本規範；docs/ 用專案規範** |
| Frontmatter | 對齊 `/big8_disk/knowledge/` | docs/ 通常無 frontmatter | KB 獨立 |
| 語言 | 繁中為主，術語/code 保英文 | 同專案規範 | 一致 |

---

**相關**：
- 上層：[../README.md](../README.md)、[../AGENT.md](../AGENT.md)
