<!--
建立時間: 2026-01-12 00:00
更新時間: 2026-03-05 10:00
狀態: validated
資料來源:
  - docs/standards/20260228_文件命名與狀態管理規範_01.md
  - docs/standards/20260228_output軟連結與版本控管規範_01.md
  - scripts/analysis/check_ai_agent_readiness.sh
-->

# 當前目標

## 1. 目前狀態

1. docs 重整已完成核心落地：
   - 命名：`YYYYMMDD_主題_流水號.md`
   - 報告分層：`reports/validated|finalized`
   - 實驗分層：`experiments/in_progress|validated|finalized`
2. `output/` 入口已固定為軟連結：
   - `output -> /big8_disk/liaoyoyo2001/InterSubMod_runs/output`
3. Knowledge MCP 已接入：
   - `.mcp.json` 指向 `/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py`

## 2. AI Agent 主要入口

1. docs 導航：`docs/README.md`
2. 研究歷史索引：`docs/experiments/INDEX.md`（已試驗方向、成功/失敗結論、建議後續）
3. Agent 手冊：`docs/references/manual/20260301_AI_Agent_快速操作手冊_01.md`
4. 健康檢查：`scripts/analysis/check_ai_agent_readiness.sh`
5. 文件規範：`docs/standards/README.md`

## 3. 當前進行中

1. 降低舊語彙殘留（例如舊輸出路徑與舊文件分層名稱）。
2. 收斂腳本路徑設定（優先使用可配置 `OUTPUT_ROOT` 或 `output/` 入口）。
3. 清理不必要空目錄，降低 Agent 探索噪音。

## 4. 阻塞與風險

1. `archive/deep/` 為 immutable 快照，保留歷史失效連結（不回寫）。
2. 活躍腳本仍有部分硬編碼歷史路徑，可能影響跨專案重用。
