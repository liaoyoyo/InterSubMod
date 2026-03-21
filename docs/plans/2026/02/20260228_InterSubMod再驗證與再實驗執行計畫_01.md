<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# InterSubMod 再驗證與再實驗執行計畫

- 文件日期：2026-02-28
- 計畫期間：建議 6 週（可壓縮為 4 週）
- 計畫目的：把目前「已觀察到的現象」轉成「可重現、可比較、可決策」的研究證據

## 1. 實作總目標

1. 建立統一 benchmark 口徑（跨樣本/跨純度/跨區域）。
2. 將現行顯著性硬閾值策略升級為多特徵融合策略。
3. 補齊 purity-aware methylation 驗證資料缺口。
4. 形成可提交與可審查的一致化報告組。

## 2. Workstreams 與任務清單

### WS-A：基準口徑統一（P0）

| 任務ID | 任務 | 輸入 | 產出 | 驗收條件 |
|---|---|---|---|---|
| A1 | 固定 benchmark 規範（BED/PASS/SNV-only） | SEQC2 + caller VCF | `docs/references/manual/` 規範文檔 | 所有後續報告引用同一規範 |
| A2 | 跨樣本重跑最小集合（HCC1395 + 2 樣本） | TP/FP/FN VCF + BAM | 統一格式 metrics 表 | 每樣本都可重現 TP/FP/FN 指標 |
| A3 | 建立結果彙總表 | A2 輸出 | `docs/experiments/outputs/...` 匯總 CSV | 支援橫向比較與排序 |

### WS-B：特徵與判定策略（P0）

| 任務ID | 任務 | 輸入 | 產出 | 驗收條件 |
|---|---|---|---|---|
| B1 | 把 `Significant` 由 hard filter 改為 feature | `significance_summary.csv` | 新版 feature schema | 可同時保留 recall/precision 調整空間 |
| B2 | 加入 HP balance、coverage、region 註解 | BAM/VCF + annotation | 增強版 feature matrix | 至少可重現 1 次 end-to-end |
| B3 | 比較 rule-based vs 輕量模型 | B2 特徵 | 比較報告（F1、PR-AUC、校準） | 新策略在 >=2 樣本優於 baseline |

### WS-C：Subclone 與高風險區驗證（P0/P1）

| 任務ID | 任務 | 產出 | 驗收條件 |
|---|---|---|---|
| C1 | `Subclone` 召回提升試驗（分層閾值） | Subclone 專項報告 | 召回提升且 precision 可接受 |
| C2 | chr9/高風險熱區系統檢查 | 問題區域清單 | 可以解釋主要高訊號 FP 來源 |
| C3 | 黑/灰名單策略驗證 | blacklist 版結果比較 | 高風險區 FP 明顯下降 |

### WS-D：資料與流程可維運性（P1）

| 任務ID | 任務 | 產出 | 驗收條件 |
|---|---|---|---|
| D1 | 補齊或替代 `data/bam/test.bam` 測試 fixture | 測試資料+說明 | BamReader 測試不再長期 skip |
| D2 | 修正文檔舊命令與舊輸出結構 | README/Quickstart 同步 | 新人可依文檔一次跑通 |
| D3 | 統一報告模板（含來源欄位） | `docs/reports` 模板 | 每份報告都有來源與命令 |

### WS-E：MCP 與知識庫接入治理（P1）

| 任務ID | 任務 | 產出 | 驗收條件 |
|---|---|---|---|
| E1 | 在本專案與其他專案統一 `.mcp.json` 範本 | 配置模板文件 | `knowledge` server 可在各專案列出資源 |
| E2 | 建立「回覆前先查 Knowledge」檢查清單 | AGENTS/流程文件 | 對路徑/格式/參數問題先查證再回答 |
| E3 | 建立 MCP smoke test 命令 | 驗證腳本或命令清單 | 每次重開後可快速確認連線狀態 |

## 3. 建議時程（6 週）

| 週次 | 主要工作 |
|---|---|
| 第 1 週 | WS-A（A1/A2）+ WS-D（D2） |
| 第 2 週 | WS-B（B1/B2） |
| 第 3 週 | WS-B（B3）+ WS-C（C1） |
| 第 4 週 | WS-C（C2/C3） |
| 第 5 週 | WS-D（D1）+ WS-E（E1/E3） |
| 第 6 週 | 匯總報告、決策 Gate 評估、下一輪計畫 |

## 4. 實驗設計（最小可行版）

### Experiment-01：Baseline 重現
- 目的：建立所有後續比較基準
- 方法：固定 BED + PASS + SNV-only，重跑 3 樣本
- 成功條件：每樣本輸出完整 TP/FP/FN + feature 表

### Experiment-02：Feature 擴充
- 目的：測試多特徵是否穩定優於單一顯著性
- 方法：以 `Significant`、`HPMergedDelta`、coverage、HP balance、region 註解建模
- 成功條件：至少 2 樣本 F1 與 PR-AUC 同步提升

### Experiment-03：Subclone 專項
- 目的：提升 `VerificationClass=Subclone` 的實用召回
- 方法：建立 Subclone 專屬閾值/分類器與誤差分析
- 成功條件：Subclone 召回提升且整體 precision 不顯著下降

### Experiment-04：高風險區治理
- 目的：降低高訊號 FP
- 方法：黑/灰名單 + 區域註解分層報告
- 成功條件：高風險區 FP 下降，且主要 TP 保留

## 5. KPI 與驗收指標

| 指標類型 | 指標 |
|---|---|
| 模型效能 | F1、Precision、Recall、PR-AUC、Calibration |
| 亞群能力 | Subclone 召回/precision |
| 泛化能力 | 跨樣本表現變異係數（CV） |
| 可重現性 | 同命令重跑結果一致性 |
| 維運性 | 測試 skip 下降、文檔可執行率 |

## 6. 里程碑輸出文件（本輪應新增）

1. `docs/experiments/outputs/...`：跨樣本統一 benchmark 結果表
2. `docs/reports/2026/..`：策略比較與決策報告
3. `docs/references/manual/...`：唯一正式流程命令文檔
4. `docs/plans/2026/...`：下一輪計畫與 Gate 結果

## 7. MCP 跨專案接入範本（可直接用）

```json
{
  "mcpServers": {
    "knowledge": {
      "type": "stdio",
      "command": "python3",
      "args": ["/big8_disk/liaoyoyo2001/Knowledge/scripts/mcp/knowledge_server.py"]
    }
  }
}
```

### 重開後最小驗證清單
1. 檢查 `.mcp.json` 是否有 `knowledge` 條目。
2. 執行 MCP 資源列舉（應可見 `knowledge://stats`）。
3. 隨機讀一份 `knowledge://doc/...` 確認可取回內容。
4. 針對路徑/格式/參數問題，先查對應 Knowledge 文件再回覆。

## 8. 依據來源

- 根據 `Knowledge/05_tools/InterSubMod.md`
- 根據 `Knowledge/06_workflows/methylation_analysis.md`
- 根據 `Knowledge/04_databases/seqc2_truth_set.md`
- 根據 `Knowledge/06_workflows/benchmark_workflow.md`
- 專案輸出：`output/20260118_vcf_all_w5000_t120/*`、`output/f1_evaluation_20260119/*`

