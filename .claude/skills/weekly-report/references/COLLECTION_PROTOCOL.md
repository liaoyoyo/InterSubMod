# Phase 0-1：收集與協作規劃

## Phase 0：自動收集本週進展（AI 獨立完成）

1. 讀取研究狀態：
   ```
   docs/CURRENT_FOCUS.md
   docs/experiments/INDEX.md
   docs/reports/research_landscape/00_INDEX.md
   ```

2. 掃描本週新增檔案：
   ```bash
   git log --since="last Monday" --oneline --name-only
   ```
   重點目錄：`docs/experiments/in_progress/`, `docs/reports/validated/`, `research/*/`, `scripts/analysis/`

3. 掃描 Memory 中近期 project/feedback 記錄

4. Knowledge Base 查閱（為 Layer 0.2 準備）：
   a. 辨識 Layer 2 主題
   b. 對照「主題→知識庫對照表」（見 LAYER_STRUCTURE.md）
   c. 使用 `knowledge_search` / `knowledge_get_doc` 讀取概念定義
   d. 產出「背景知識計畫」
   e. 讀取上週週報 Layer 3 + Layer 4

5. 輸出：結構化進展清單

   | 日期 | 主題 | 類型 | 狀態 | 關鍵數字 | 來源檔案 |
   |------|------|------|------|---------|---------|

---

## Phase 1：協作規劃（AI 整理，用戶決策）

### 1.1 本週進展優先排序草稿

依價值排序規則排列：

| 優先序 | 類型 | 說明 |
|--------|------|------|
| 1 | 方向性翻轉 / 異常發現 | 預期與結果相反 |
| 2 | 跨樣本一致性共識 | 7/7 樣本一致 |
| 3 | 方法學修正 / Bug Fix | 影響全域數據 |
| 4 | Reproducibility 多線收斂 | 多方法同結論 |
| 5 | NEGATIVE 正式關閉 | 假說否決 |
| 6 | POSITIVE 量化確認 | 假說支持 |
| 7 | 流程改進 / 治理 | 文件、工具 |

**請用戶確認**排序和取捨。

### 1.2 敘事主軸建議

| 類型 | 適用場景 | 結構 |
|------|---------|------|
| 偵探故事 | 重大發現/Bug Fix | 發現→懷疑→驗證→根因→修正 |
| 里程碑收斂 | 假說驗證密集 | N 假說→逐一驗證→收斂 |
| 系統性掃描 | 觀察系列 (O1-O15) | 維度1→維度2→交叉分析 |
| 方向定錨 | 策略/治理 | 現況→選項→評估→路線圖 |

**請用戶選擇**。

### 1.3 背景知識與前情提要計畫

**Part A** — Layer 0.2 背景知識：概念群組（最多 3 組）+ KB 來源 + PI 已知/需提醒
**Part B** — Layer 0.3 前情提要：上週哪些結論與本週直接相關
**Part C** — Layer 2 定義區塊需求

### 1.4 遺漏檢查

1. 上週待辦追蹤
2. 未記錄決策
3. 數據完整性（圖表已產出但未寫入報告？）

### 1.5 PPTX 規劃

- 預估頁數（12-18）和時長（15-25 min）
- 圖表來源：已有現成 vs 需新製
- 腳本選擇

---

## 過去週報範本

| 週報 | 特色 | 路徑 |
|------|------|------|
| 0310 | 三層方法分工 | `docs/reports/validated/2026/03/20260310_研究主線週報_20260305_20260310_01.md` |
| 0330 | HP Bug Fix + LOH | `docs/reports/validated/2026/03/20260330_研究週報_20260325_20260330_HP_bug_fix與LOH_evidence_panel_01.md` |
| LOH 簡報 | 偵探故事線 PPTX | `docs/presentations/validated/2026/04/20260401_LOH_weekly_report_draft/` |

## 研究脈絡速查

| 資訊 | 來源 |
|------|------|
| 當前狀態 + 阻塞 | `docs/CURRENT_FOCUS.md` |
| 實驗歷史索引 | `docs/experiments/INDEX.md` |
| 完整推論鏈 | `docs/reports/research_landscape/00_INDEX.md` |
| TO FP 問題全貌 | `docs/reports/research_landscape/01_TO_FP問題全貌.md` |
| ISM 分析價值界定 | `docs/reports/research_landscape/03_ISM分析價值界定.md` |
| 結論穩定性 | `docs/reports/research_landscape/06_結論穩定性審查.md` |
