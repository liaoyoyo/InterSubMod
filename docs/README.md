<!--
建立時間: 2026-01-11 10:00
更新時間: 2026-03-05 10:00
目標: 提供 docs/ 目錄的最新結構、命名規範與工作流程，並提供 AI 漸進查閱指引
處理範圍: docs/ 全目錄（archive/deep 歷史快照除外）
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
  - docs/standards/20260228_文件命名與狀態管理規範_01.md
  - docs/standards/20260228_output軟連結與版本控管規範_01.md
  - docs/standards/20260303_文件盤點分類與歸檔流程規範_01.md
-->

# InterSubMod 文檔庫

## 快速導航

- [當前目標](CURRENT_FOCUS.md)
- [AI Agent 快速操作手冊](references/manual/20260301_AI_Agent_快速操作手冊_01.md)
- [專案檔案清冊與查找流程手冊](references/manual/20260303_專案檔案清冊與查找流程手冊_01.md)
- [2026-03-03 全量盤點報告](reports/finalized/2026/03/20260303_docs與專案檔案全量盤點報告_01.md)
- [2026-03-04 scripts/tools/output 重整與遷移報告](reports/finalized/2026/03/20260304_scripts_tools_output流程重整與空間遷移報告_01.md)
- `standards/`：規範與治理文件
- `plans/`：計畫
- `reports/`：報告（validated/finalized）
- `experiments/`：實驗（in_progress/validated/finalized）
- `archive/`：歸檔（含 `deep/` 歷史快照）

## AI Agent 建議起手式

```bash
scripts/analysis/check_ai_agent_readiness.sh
```

## AI 漸進查閱指引

本文件庫採用 4 層漸進披露架構，避免 AI 直接陷入大量細節文件：

| 層 | 文件 | 目的 |
|---|---|---|
| Layer 0 (地圖) | `docs/README.md`（本頁）| 全局導航，方向感 |
| Layer 1 (當下) | `docs/CURRENT_FOCUS.md` | 現在在做什麼、阻塞點 |
| Layer 2 (歷史) | `docs/experiments/INDEX.md` | 已試驗方向，成功/失敗總覽 |
| Layer 3 (細節) | 各 experiments/reports/solutions 文件 | 完整數據與推導 |

### 建議查閱流程

**新任務起手式：**
1. 讀 `CURRENT_FOCUS.md` → 確認當前優先事項
2. 讀 `experiments/INDEX.md` → 避免重複已失敗的方向
3. 讀相關 `architecture/*.md` → 了解系統設計約束
4. 執行 `scripts/analysis/check_ai_agent_readiness.sh` → 確認環境狀態

**特定問題查閱：**
- 找過去的解決方案 → `solutions/{topic}/`
- 找已驗證的實驗結果 → `experiments/validated/` 或 `experiments/finalized/`
- 找系統架構決策 → `architecture/` 或 `decisions/`
- 找研究歷史與方向優先級 → `experiments/INDEX.md`

## 目錄結構（現行）

| 目錄 | 用途 |
|---|---|
| `architecture/` | 長期架構設計文件 |
| `concepts/` | 構想與設計草稿 |
| `plans/YYYY/MM/` | 執行計畫與里程碑 |
| `reports/validated/YYYY/MM/` | 已驗證分析報告 |
| `reports/finalized/YYYY/MM/` | 最終結論與決策報告 |
| `experiments/in_progress/YYYY/MM/` | 實驗草稿 |
| `experiments/validated/YYYY/MM/` | 可重現驗證結果 |
| `experiments/finalized/YYYY/MM/` | 實驗最終結論 |
| `solutions/{topic}/YYYY/MM/` | 問題解法與修復紀錄 |
| `research/{topic}/YYYY/MM/` | 專題研究文件 |
| `references/` | 參考文件與手冊 |
| `decisions/YYYY/MM/` | 重整與治理決策紀錄 |
| `archive/YYYY/MM/` | 一般歸檔 |
| `archive/deep/` | 歷史快照（保留原貌，不回溯改名） |

## 研究歷史與實驗索引

→ **[實驗總索引](experiments/INDEX.md)**：所有已嘗試方向的成功/失敗/建議後續

主要研究主題（依時間軸）：
- 甲基化解析與 CIGAR 座標映射（2025-11）✅ 已完成
- 距離計算、聚類分析、Bernoulli 度量（2025-11 ~ 2025-12）✅ 已完成
- 統計顯著性分析（Fisher / PERMANOVA / Cramér's V）（2025-12 ~ 2026-01）✅ 已驗證
- TP/FP 特徵富集分析與 F1 最佳化（2026-01）✅ F1=0.8481
- Subsample 混樣甲基化偏差分析（2026-02 ~ 2026-03）⏳ 進行中
- Purity-Aware 策略驗證（2026-02 ~ 2026-03）🔄 值得再探索

## 命名規範

### 標準格式

```text
YYYYMMDD_主題_流水號.md
```

### 例外

1. 固定名稱：`README.md`、`CURRENT_FOCUS.md`、`INDEX.md`
2. 長期架構文件：`snake_case.md`
3. `archive/deep/` 歷史快照：維持原檔名

## 文件狀態規範

### reports/

1. `validated`：可重跑驗證，可供內部引用
2. `finalized`：最終對外口徑

### experiments/

1. `in_progress`：草稿/探索中
2. `validated`：驗證完成
3. `finalized`：整體結論定稿

## 建議工作流程

1. 新增文件：先決定狀態層（in_progress/validated/finalized）
2. 套用命名規範並填寫 metadata
3. 更新關聯報告索引
4. 需要歷史保留時移至 `archive/YYYY/MM/`

## 如何新增文件

1. 確定文件類型（實驗、報告、解決方案、計畫、AI 對話紀錄）
2. 選擇狀態層（in_progress / validated / finalized）
3. 套用命名：`YYYYMMDD_主題_流水號.md`
4. 填寫 metadata 區塊（見下方範本）
5. 若為實驗，同步更新 `experiments/INDEX.md`
6. 若解決了重要問題，在 `solutions/` 下補充紀錄

**Metadata 範本：**
```markdown
<!--
建立時間: YYYY-MM-DD HH:MM
目標: [本檔案的目標或用途]
處理範圍: [涵蓋的工作範圍]
關聯檔案:
  - [相關檔案路徑]
-->
```

## 2026-03-03 盤點輸出

1. 全專案檔案清冊：`reports/validated/2026/03/assets/20260303_repository_full_file_inventory_01.tsv`
2. 全專案目錄清冊：`reports/validated/2026/03/assets/20260303_repository_directory_inventory_01.txt`
3. docs 子樹清冊：`reports/validated/2026/03/assets/20260303_docs_file_inventory_01.tsv`
4. 本次歸檔待審區：`archive/2026/03/20260303_ai_sessions_raw_artifacts_pending_review_01/`

## 相容性說明

- 2026-02-28 後已啟用新分層與命名規範。
- 舊路徑可能已重整，請優先從 `reports/finalized` 與 `reports/validated` 查找。
