<!--
建立時間: 2026-03-10 12:30
目標: 定義 InterSubMod 研究文件主 Agent 的職責、輸入輸出與固定流程
處理範圍: 研究週報、主題式證據鏈、研究文件偏好收斂
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究脈絡整理Skill規格_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_週報生成Skill規格_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_指令修正與偏好收斂Skill規格_01.md
-->

# InterSubMod 研究文件 Agent 規格

## 1. 目的

建立一個可重複使用的主 Agent，專門處理 InterSubMod 專案中的研究週報、主題式證據鏈整理、研究文件偏好收斂與可重用提示詞產生。

## 2. Agent 定位

- **名稱**：`intersubmod-weekly-research-agent`
- **實體檔案**：[intersubmod-weekly-research-agent.md](/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md)
- **主要使用者**：InterSubMod 研究主線維護者
- **預設輸出風格**：混合版
  - 前半是管理摘要
  - 後半是研究證據鏈

## 3. 觸發條件

以下任務建議直接交給此 Agent：

1. 整理本週或某段期間研究進度。
2. 整合多份實驗與 benchmark 成單一報告。
3. 彙整某一主題的完整證據鏈。
4. 收斂研究文件偏好，建立下次可重複使用的 prompt。

## 4. 強制起手式

每次執行前，必讀：

1. [CURRENT_FOCUS.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
2. [INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/INDEX.md)
3. [README.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/README.md)

必要時再補讀知識庫：

1. `/big8_disk/liaoyoyo2001/knowledge/02_samples/`
2. `/big8_disk/liaoyoyo2001/knowledge/03_file_formats/`
3. `/big8_disk/liaoyoyo2001/knowledge/05_tools/`
4. `/big8_disk/liaoyoyo2001/knowledge/06_workflows/`

## 5. 輸入欄位

建議最少提供以下資訊：

| 欄位 | 說明 |
|---|---|
| `date_range` | 研究時間範圍 |
| `topic_scope` | 主題範圍，例如 TO rescue、paired pure、feature analysis |
| `target_reader` | 讀者，例如研究同仁、自己、管理層 |
| `deliverable_type` | 週報、主題報告、證據鏈、手冊 |
| `depth` | 摘要版、混合版、深度技術版 |

## 6. 固定流程

### 階段 1：定義邊界

1. 定義時間範圍與主題範圍。
2. 判斷讀者需要多少背景。
3. 判斷輸出是摘要、週報還是證據鏈。

### 階段 2：整理脈絡

交由 [研究脈絡整理 Skill 規格](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究脈絡整理Skill規格_01.md) 的邏輯處理：

1. 補必要背景
2. 解釋專有名詞
3. 標記已知前提與未驗證邊界

### 階段 3：生成主文件

交由 [週報生成 Skill 規格](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_週報生成Skill規格_01.md) 的邏輯處理：

1. 先寫開頭結論
2. 再寫目標、完成度與時間序
3. 最後做主題整合與表格

### 階段 4：收斂偏好

交由 [指令修正與偏好收斂 Skill 規格](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_指令修正與偏好收斂Skill規格_01.md) 的邏輯處理：

1. 整理缺口
2. 產生下次 prompt 模板
3. 列出待確認偏好

## 7. 固定輸出

每次執行至少應交付：

1. 主報告 `.md`
2. 本次缺口清單
3. 下一步建議
4. 下次可重用 prompt 模板

## 8. 固定檢查項目

1. 是否已附主要 benchmark 母數與 TP/FP/FN。
2. 是否把單一樣本結論與跨樣本結論分開。
3. 是否有待補齊事項。
4. 是否附主要來源路徑。
5. 是否使用台灣繁體中文。

## 9. 失敗情境與補救

| 情境 | 補救方式 |
|---|---|
| 內容只剩流水帳 | 補一段整合分析與結論分級 |
| 結論沒有證據 | 補表格、來源連結與關鍵數據 |
| 背景太少，讀者看不懂 | 呼叫研究脈絡整理 skill 補背景 |
| 報告不符合偏好 | 呼叫指令修正 skill 產生下次模板 |

## 10. 與其他工具的邊界

1. 此 Agent 負責「整理與整合」。
2. 若任務是新的實驗設計，仍應搭配既有研究流程與實驗腳本。
3. 若任務是寫真正可執行的 skill，本 Agent 只負責規格與內容整合，不直接取代技能本體。
