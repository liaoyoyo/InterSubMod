# Research AI Agent Patterns for Scientific Workflow (2025-2026)

<!--
建立時間: 2026-04-14 16:00
目標: 收集 AI agent 架構模式，用於改進 InterSubMod 的研究自動化流程
處理範圍: 科學研究 AI 編排、證據鏈管理、反回歸模式、知識圖譜、自動報告、實驗筆記
關聯檔案:
  - .claude/skills/research-loop/SKILL.md
  - AGENTS.md
-->

## 1. 搜尋概述

- **關鍵詞**: scientific AI agent, research workflow orchestration, evidence provenance, anti-regression, knowledge graph memory, automated reporting, lab notebook
- **搜尋時間**: 2026-04-14
- **資料來源**: arXiv (cs.AI/cs.SE), Nature, Frontiers, GitHub, Google Research, bioRxiv
- **涵蓋時間範圍**: 2025-01 至 2026-04

## 2. 科學工作流 AI 編排

### 2.1 主要框架

| 系統 | 架構 | 核心特點 | 適用性 |
|------|------|---------|--------|
| **Google AI Co-Scientist** | 6 專化 agent (Generation/Reflection/Ranking/Evolution/Proximity/Meta-review) | Elo 錦標賽假說排名；generate-debate-evolve 循環 | 假說產生與排名機制可借鑑 |
| **Karpathy AutoResearch** | 單 agent + program.md + 固定時間預算 | 630 行 Python；700 實驗/2天；只保留改進結果 | **最直接可用** -- 與 headless-research agent 架構相似 |
| **Agent Laboratory** | 3 階段 (Literature Review / Experimentation / Report) | 支援 autonomous + co-pilot 模式；84% 成本降低 | 模式切換概念已在 InterSubMod 確認協議中實現 |
| **SelfAI** | Multi-agent: User Agent + Cognitive Agent | 軌跡推理；自適應停止；效率-多樣性平衡 | **停止機制與反冗餘最值得採用** |
| **InternAgent 1.5** | 統一框架跨 Algorithm + Empirical Discovery | 12 科學任務端到端自動化 | 規模過大，但任務分解模式可參考 |
| **AI Scientist v2** | Agentic tree search | 工作坊等級論文自動產出 | 樹狀搜索比線性循環更能避免死路 |

### 2.2 對 InterSubMod 的啟發

**AutoResearch 模式最接近現有架構**：program.md = 我們的 SKILL.md，固定預算 = 我們的 5 分鐘 region 分析，keep/revert = 我們的 hypothesis queue accept/reject。差異在於 AutoResearch 只追蹤單一指標 (val_bpb)，而我們需要多指標 (AUC, F1, effect size)。

**SelfAI 的軌跡推理**值得整合：不只記錄「這個假說失敗了」，而是記錄「在什麼條件下、用什麼方法、得到什麼結果」的完整軌跡，讓 Cognitive Agent 推斷下一步。

**Google Co-Scientist 的 Elo 排名**可用於假說優先序管理，比目前的手動排序更客觀。

## 3. 證據鏈管理

### 3.1 架構比較

| 方案 | 代表 | 優點 | 缺點 | 適用場景 |
|------|------|------|------|---------|
| **JSONL Ledger** | 現有 evidence_ledger | 簡單、git 友善、可 grep | 無關聯查詢、無時態 | 小規模 (<500 條目) |
| **Temporal KG** | Graphiti/Zep | 雙時態建模 (系統時間+事實有效期)；完整溯源 | 需 Neo4j；部署複雜 | 大規模多 session |
| **Multi-Graph** | MAGMA | 4 層正交圖 (語義/時態/因果/實體)；95% token 節省 | 2026-01 論文，尚無穩定工具 | 理論最優但不成熟 |
| **Agentic Memory** | A-MEM (NeurIPS 2025) | Zettelkasten 式自組織；記憶自動演化 | 依賴 LLM 品質 | 中規模，需自適應 |

### 3.2 建議方向

**短期 (現在可做)**：保留 JSONL ledger，但增加結構化欄位：

```jsonl
{"id":"E042","hypothesis":"H12","result":"NEGATIVE","conditions":{"mode":"TO","samples":7,"metric":"AUC"},"conclusion_stability":0.95,"timestamp":"2026-04-14T10:00:00Z","valid_until":null,"superseded_by":null}
```

**中期**：評估 Graphiti 作為跨 session 記憶層。其 episode-to-fact 溯源與 bi-temporal 模型直接解決「這個結論是什麼時候、基於什麼證據得出的」問題。

## 4. 反回歸模式

### 4.1 現有技術

目前學術界**尚無成熟的「反回歸」專用框架**。相關機制散布在不同系統中：

| 機制 | 來源 | 原理 | InterSubMod 對應 |
|------|------|------|-----------------|
| **Keep/Revert Gate** | AutoResearch | 只保留改進結果，自動回滾 | 部分實現 (hypothesis queue) |
| **Adaptive Stopping** | SelfAI | 偵測搜索路徑無產出時主動終止 | **缺少** -- 目前依賴人工 NO-GO |
| **Utility-based Deletion** | Memory Management (arXiv 2505.16067) | 基於效用分數刪除低品質記憶，防止錯誤傳播 | **缺少** |
| **Experience-Following 去偏** | 同上 | 高相似輸入導致重複輸出；需 selective storage | known-pitfalls 文件部分覆蓋 |
| **Trajectory Reasoning** | SelfAI | 累積實驗軌跡推斷下一步，避免重複 | **缺少** |
| **Elo Tournament** | Google Co-Scientist | 假說排名自動更新，低分假說自然沉底 | 可整合到 hypothesis queue |

### 4.2 建議：Confidence Decay + Tombstone 機制

```
假說狀態機:
  ACTIVE → TESTED(result) → {POSITIVE | NEGATIVE | INCONCLUSIVE}
  NEGATIVE → TOMBSTONE(reason, conditions, expiry=180d)
  TOMBSTONE → 自動阻擋相似假說 (cosine similarity > 0.85)
  TOMBSTONE → REVIVABLE (如果新證據改變前提條件)
```

- **Tombstone 記錄**：失敗假說帶有失敗條件簽名，新假說提交前自動比對
- **Confidence Decay**：超過 90 天未被引用的 POSITIVE 結論自動降級為 NEEDS_REVALIDATION
- **條件式復活**：如果前提條件改變 (如新資料、新方法)，允許 TOMBSTONE 復活但需明確標註差異

## 5. 研究知識圖譜

### 5.1 工具景觀

| 工具 | 類型 | 記憶模型 | 整合難度 |
|------|------|---------|---------|
| **Graphiti** (getzep/graphiti) | Temporal KG engine | Entity+Relation+Episode；bi-temporal | 中 (需 Neo4j) |
| **A-MEM** (NeurIPS 2025) | Zettelkasten-style | Note 自動生成 context/keywords/tags/links | 低 (純 LLM) |
| **MAGMA** | Multi-graph | 4 正交圖層 + intent-aware traversal | 高 (研究原型) |
| **Cognee** | Hybrid KG+Vector | 圖結構 + 向量嵌入統一記憶 | 中 |

### 5.2 建議：輕量 Zettelkasten + JSONL

對 InterSubMod 規模 (約 30 已結案研究方向、14 穩定結論、8 證據鏈)，**A-MEM 的 Zettelkasten 方法最適合**：

- 每個結論 = 一張卡片，含 context、keywords、links to evidence
- 卡片間自動建立連結 (因果、矛盾、支持)
- 不需外部資料庫，可用 JSONL/Markdown 實現
- 與現有 `docs/reports/research_landscape/` 結構高度相容

## 6. 自動研究報告

### 6.1 引用驗證

| 工具 | 功能 | 開源 | 適用性 |
|------|------|------|--------|
| **SemanticCite** | 4 階段 pipeline：預處理→混合檢索→神經重排→LLM 分析；4 類分類 (Supported/Partially/Unsupported/Uncertain) | 是 | 論文寫作階段 |
| **sciwrite-lint** | 22 項完整性檢查：統計一致性、引用驗證、圖文對齊、撤稿偵測 | 是 | **高度推薦** -- 可整合到報告產生流程 |
| **Scite.ai** | 16 億+ 引用索引；Smart Citation 分類 | 商業 | 文獻搜尋階段 |
| **Statcheck** | 自動掃描 p-value 與統計量一致性 | 是 | 統計報告驗證 |

### 6.2 建議整合

**sciwrite-lint** (arXiv 2604.08501, 2026-04-09) 最值得整合：
- 22 項檢查中至少 10 項直接適用於 ISM 研究報告
- 統計數字 vs 表格一致性檢查可防止報告中的數值錯誤
- 使用 Qwen3 8B 本地運行，不需外部 API
- pip 可安裝，可加入 CI/CD

## 7. 實驗筆記模式

### 7.1 AI 時代的 ELN 趨勢

2025-2026 ELN 趨勢聚焦於：
- **自動溯源**：每個操作自動記錄時間戳、軟體版本、參數
- **可執行筆記**：筆記本身可重現 (bundled code + data + environment)
- **AI 整合**：語音命令、自動建議、趨勢分析

### 7.2 InterSubMod 現有 vs 建議改進

| 面向 | 現有 | 建議改進 |
|------|------|---------|
| 實驗紀錄 | `docs/experiments/` Markdown | 增加機器可讀 frontmatter (YAML) |
| 溯源追蹤 | `docs/provenance/ai_sessions/` | 增加 session-to-conclusion 自動連結 |
| 可重現性 | 手動記錄命令 | 每個實驗紀錄嵌入 `reproduce:` 區塊含完整命令 |
| 版本鎖定 | git commit hash | 增加 `environment:` 區塊記錄 tool versions |

## 8. 生物資訊專用 AI Agent

### 8.1 領域系統

| 系統 | 發表 | 架構 | 特點 |
|------|------|------|------|
| **BioAgents** | Nature Sci Reports 2025 | 小模型 + RAG + 自評估 | 發現反覆精煉有遞減回報 |
| **PromptBio** | bioRxiv 2025 | DataAgent + OmicsAgent + AnalysisAgent + QAgent | stepwise human-in-the-loop |
| **AutoBA** | PMC 2024 | LLM-driven pipeline 生成 | 最少人工輸入 |
| **CellAtria** | Nature 2025 | AI agent + 自動化 pipeline | 對話驅動分析 |

### 8.2 對 InterSubMod 的啟示

BioAgents 的發現 -- **反覆精煉有遞減回報** -- 與我們的經驗一致 (O1-O13 系統性觀察中多輪迭代收益遞減)。這支持 SelfAI 的自適應停止機制：偵測到改進幅度低於閾值時主動終止。

## 9. 優先行動建議

| 優先序 | 行動 | 預期收益 | 實作複雜度 |
|--------|------|---------|-----------|
| 1 | **Tombstone 機制**：JSONL 記錄失敗假說簽名，新假說自動比對 | 防止重複已結案的 29 個方向 | 低 |
| 2 | **JSONL Ledger 結構化升級**：增加 conditions/stability/superseded_by 欄位 | 證據鏈可查詢性 | 低 |
| 3 | **Adaptive Stopping**：在 research-loop SKILL 中增加改進幅度閾值判斷 | 避免無效迭代 | 低 |
| 4 | **sciwrite-lint 整合**：報告產生後自動執行統計一致性檢查 | 報告品質保證 | 中 |
| 5 | **Zettelkasten 結論卡片**：為 14 穩定結論建立結構化卡片 + 自動連結 | 知識可查詢性 | 中 |
| 6 | **Graphiti 評估**：PoC 測試跨 session 記憶持久化 | 長期記憶管理 | 高 |

## 10. 資料來源評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| Google AI Co-Scientist (arXiv 2502.18864) | 預印本+部落格 | 高 | 有實驗驗證 (肝纖維化、AML) |
| Karpathy AutoResearch (GitHub) | 開源工具 | 高 | 21K stars；Shopify 驗證 |
| SelfAI (arXiv 2512.00403) | 預印本 | 中-高 | 跨領域驗證；v2 版本 2026-02 |
| MAGMA (arXiv 2601.03236) | 預印本 | 中 | 2026-01；尚無廣泛採用 |
| A-MEM (NeurIPS 2025) | 頂會論文 | 高 | 同行評審通過 |
| Graphiti/Zep (arXiv 2501.13956) | 預印本+產品 | 高 | 生產級部署；DMR benchmark 領先 |
| sciwrite-lint (arXiv 2604.08501) | 預印本 | 中 | 極新 (2026-04-09)；開源可驗證 |
| BioAgents (Nature Sci Reports) | 期刊論文 | 高 | 同行評審 |
| Agent Laboratory (arXiv 2501.04227) | 預印本 | 中-高 | 有 AgentRxiv 後續 |
| Memory Management Study (arXiv 2505.16067) | 預印本 | 中-高 | 2025 empirical study |
| SemanticCite (arXiv 2511.16198) | 預印本 | 中-高 | 開源 + 1000 引用資料集 |
| InternAgent 1.5 (GitHub) | 技術報告 | 中 | 2026-02；規模大但通用 |
